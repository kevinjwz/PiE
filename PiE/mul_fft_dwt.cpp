#include<complex>
#include<math.h>
#include<assert.h>
#include "bigint.h"
#include "micro_helpers.inl"
#include "array_helpers.h"
#include "array_basic_arith.h"
#include "fft_rad2_rad4.h"
#include "fft_constants.h"
#include "dcomplex.inl"

using namespace std;

template <class T>
static void free(vector<T>& v)
{
	vector<T> empty;
	v.swap(empty);
}

/*
transform
[re0, im0, re1, im1, ...]
to
[re0, re1, re2, re3, ...
im0, im1, im2, im3, ...]
*/
static void complex_h2v(dcomplex result[], const dcomplex x[], size_t n)
{
	double* r_re = (double*)result;
	double* r_im = (double*)result + n;
	for (size_t i = 0; i < n; i += 1) {
		r_re[i] = x[i].real();
		r_im[i] = x[i].imag();
	}
}



static double estimate_error(size_t k, size_t elem_bits)
{
	double base = exp2(elem_bits);
	double kd = (double)k;
	//return (6 * 1e-16) * kd*kd * base*base * log2(kd);
	return (3 * 1e-16) * kd * base*base;

}

static size_t search_bits(size_t min, size_t max, size_t n)
{
	if (min > max)
		return 0;
	size_t mid = (min + max) / 2;
	size_t k = clp2((n * 32 + mid - 1) / mid);
	double error = estimate_error(k, mid);
	if (error < 0.5) {
		size_t right = search_bits(mid + 1, max, n);
		return right == 0 ? mid : right;
	}
	else {
		return search_bits(min, mid - 1, n);
	}
}

static void select_sizes(size_t& k, size_t& elem_bits, size_t n)
{
	// no overflow during computing k
	assert(n <= SIZE_MAX / 32);
	// binary search
	elem_bits = search_bits(1, 32, n);
	assert(elem_bits > 0);	// if this fails, n is too large 
	k = clp2((n * 32 + elem_bits - 1) / elem_bits);
}

static void extract_elems(vector<dcomplex>& elems, const bigint& src, size_t elem_bits)
{
	assert(elem_bits <= 32);
	uint32_t mask = elem_bits == 32 ? -1 : ((1 << elem_bits) - 1);
	size_t word_i = 0;
	size_t bit_offset = 0;
	size_t i = 0;
	while (word_i < src.size()) {
		if (bit_offset + elem_bits <= 32 || word_i == src.size() - 1) {
			elems[i] = (src[word_i] >> bit_offset) & mask;
		}
		else {
			uint32_t lo = src[word_i] >> bit_offset;
			uint32_t hi = src[word_i + 1] << (32 - bit_offset);
			elems[i] = (lo | hi) & mask;
		}
		bit_offset += elem_bits;
		if (bit_offset >= 32) {
			word_i += 1;
			bit_offset %= 32;
		}
		i += 1;
	}
	for (; i < elems.size(); ++i)
		elems[i] = 0;
}

static void combine_elems(bigint& dst, const vector<dcomplex>& elems, size_t elem_bits)
{
	assert(elem_bits <= 32);
	size_t word_i = 0;
	size_t bit_offset = 0;
	size_t i = 0;
	set(&dst[0], 0, dst.size());
	while (word_i < dst.size()) {
		uint32_t elem = (uint32_t)elems[i].real();
		if (bit_offset + elem_bits <= 32 || word_i == dst.size() - 1) {
			dst[word_i] |= elem << bit_offset;
		}
		else {
			dst[word_i] |= elem << bit_offset;
			dst[word_i + 1] |= elem >> (32 - bit_offset);
		}
		bit_offset += elem_bits;
		if (bit_offset >= 32) {
			word_i += 1;
			bit_offset %= 32;
		}
		i += 1;
	}
}

static void adjust_carry(vector<dcomplex>& elems, size_t elem_bits)
{
	assert(elem_bits <= 32);
	uint64_t carry = 0;
	for (size_t i = 0; i < elems.size(); ++i) {
		// r <= k*(B-1)^2
		uint64_t r = (uint64_t)max(0.0, round(elems[i].real()));
		r += carry;
		assert(r >= carry);
		carry = r >> elem_bits;
		r &= (1 << elem_bits) - 1;
		elems[i].real((double)r);
	}
}

void mul_fft_fp(bigint& prod, const bigint& a, const bigint& b)
{
	size_t prod_n = a.size() + b.size();
	size_t k, elem_bits;
	select_sizes(k, elem_bits, prod_n);
	
	vector<dcomplex> fa(k);
	extract_elems(fa, a, elem_bits);
	fft_mm_radix_4<false>(&fa[0], k);

	vector<dcomplex> fb(k);
	extract_elems(fb, b, elem_bits);
	fft_mm_radix_4<false>(&fb[0], k);

	for (size_t i = 0; i < k; ++i)
		fa[i] *= fb[i];
	free(fb);
	
	fft_mm_radix_4<true>(&fa[0], k);

	adjust_carry(fa, elem_bits);
	prod.clear();
	prod.resize(prod_n);
	combine_elems(prod, fa, elem_bits);
	prod.resize(prod_n - high_zeros(prod));
}

/*
size[elems]=h
size[src]=n
*/
static void extract_elems_dwt(dcomplex elems[], size_t h, const uint32_t src[], size_t n, size_t elem_bits)
{
	assert(elem_bits <= 32);
	uint32_t mask = elem_bits == 32 ? -1 : ((1 << elem_bits) - 1);
	size_t word_i = 0;
	size_t bit_offset = 0;
	size_t i = 0;
	double(*cmplx_arr)[2] = (double(*)[2])&elems[0];
	size_t re_im = 0;	// 0: re; 1:im

	for (re_im = 0; re_im < 2; ++re_im)
		for (i = 0; i < h; ++i) {
			if (word_i >= n)
				goto fill_zeros;
			if (bit_offset + elem_bits <= 32 || word_i == n - 1) {
				cmplx_arr[i][re_im] = (src[word_i] >> bit_offset) & mask;
			}
			else {
				uint32_t lo = src[word_i] >> bit_offset;
				uint32_t hi = src[word_i + 1] << (32 - bit_offset);
				cmplx_arr[i][re_im] = (lo | hi) & mask;
			}
			bit_offset += elem_bits;
			if (bit_offset >= 32) {
				word_i += 1;
				bit_offset %= 32;
			}
		}
fill_zeros:
	if (re_im == 0) {
		for (; i < h; ++i)
			cmplx_arr[i][0] = 0;
		for (i = 0; i < h; ++i)
			cmplx_arr[i][1] = 0;
	}
	else if (re_im == 1) {
		for (; i < h; ++i)
			cmplx_arr[i][1] = 0;
	}
}

/*
for all e in elems, e in [0, k*(B-1)^2]
*/
static void adjust_carry(double elems[], size_t n, size_t elem_bits)
{
	assert(elem_bits <= 32);
	uint64_t carry = 0;
	for (size_t i = 0; i < n; ++i) {
		// r <= k*(B-1)^2 < k*B^2 < 1.7e15
		uint64_t r = max(0LL, round_to_int64(elems[i]));
		r += carry;	// r <= k*(B-1)^2 + k*B < k*B^2
		assert(r >= carry);	// no overflow
		carry = r >> elem_bits;	// carry <= k*B
		r &= (1 << elem_bits) - 1;
		elems[i] = (double)r;
	}
}

/*
size(dst)=n
size(elems)=k
*/
static void combine_elems(uint32_t dst[], size_t n, const double elems[], [[maybe_unused]]size_t k, size_t elem_bits)
{
	assert(elem_bits <= 32);
	size_t word_i = 0;
	size_t bit_offset = 0;
	size_t i = 0;
	set(dst, 0, n);
	while (word_i < n) {
		uint32_t elem = (uint32_t)elems[i];
		if (bit_offset + elem_bits <= 32 || word_i == n - 1) {
			dst[word_i] |= elem << bit_offset;
		}
		else {
			dst[word_i] |= elem << bit_offset;
			dst[word_i + 1] |= elem >> (32 - bit_offset);
		}
		bit_offset += elem_bits;
		if (bit_offset >= 32) {
			word_i += 1;
			bit_offset %= 32;
		}
		i += 1;
	}
}

static void weight(vector<dcomplex>& x, const vector<dcomplex>& w)
{
	size_t h = x.size();
	if (h <= 1)
		return;
	for (size_t i = 0; i < h / 2; ++i) {
		x[i] *= w[i];
	}
	x[h / 2] *= {SQRT2 / 2, SQRT2 / 2};
	for (size_t i = h / 2 + 1; i < h; ++i) {
		x[i] *= reflect(w[h - i]);
	}
}

static void unweight(vector<dcomplex>& x, const vector<dcomplex>& w)
{
	size_t h = x.size();
	if (h <= 1)
		return;
	for (size_t i = 0; i < h / 2; ++i) {
		x[i] *= conj(w[i]);
	}
	x[h / 2] *= {SQRT2 / 2, -SQRT2 / 2};
	for (size_t i = h / 2 + 1; i < h; ++i) {
		x[i] *= conj(reflect(w[h - i]));
	}
}

/*
size[prod] = an+bn
*/

void mul_dwt_raw(uint32_t prod[], const uint32_t a[], size_t an, const uint32_t b[], size_t bn)
{
	size_t prod_n = an + bn;
	assert(prod_n >= an);
	if (prod_n == 0)
		return;

	size_t k, elem_bits;
	select_sizes(k, elem_bits, prod_n);
	size_t h = k / 2;	// prod_n>=1, k>=2, h>=1 
	vector<dcomplex> w(h / 2);
	for (size_t i = 0; i < h / 2; ++i) {
		double angle = (2 * PI) * ((double)i / (4 * h));
		w[i] = { cos(angle), sin(angle) };
	}

	vector<dcomplex> fa(h);
	extract_elems_dwt(&fa[0], h, a, an, elem_bits);
	weight(fa, w);
	fft_mm_radix_4<false>(&fa[0], h);

	vector<dcomplex> fb(h);
	extract_elems_dwt(&fb[0], h, b, bn, elem_bits);
	weight(fb, w);
	fft_mm_radix_4<false>(&fb[0], h);

	for (size_t i = 0; i < h; ++i)
		fa[i] *= fb[i];
	
	fft_mm_radix_4<true>(&fa[0], h);
	unweight(fa, w);
	free(w);

	complex_h2v(&fb[0], &fa[0], h);
	free(fa);
	double* d_prod = (double*)&fb[0];

	adjust_carry(d_prod, k, elem_bits);
	combine_elems(prod, prod_n, d_prod, k, elem_bits);
}

void mul_dwt(bigint& prod, const bigint_src& a, const bigint_src& b)
{
	assert(!is_from(a, prod) && !is_from(b, prod));
	size_t an = a.size(), bn = b.size();
	if (an == 0 || bn == 0) {
		prod.clear();
		return;
	}
	size_t prod_n = an + bn;
	assert(prod_n >= an);
	prod.clear();
	prod.resize(prod_n);
	mul_dwt_raw(&prod[0], &a[0], an, &b[0], bn);
	if (prod[prod_n - 1] == 0)
		prod.pop_back();
}

/*
for all e in elems, e in [-(k-1)*(B-1)^2, k*(B-1)^2]
return: carry (may < 0)
*/
static int64_t adjust_carry_mod(double elems[], size_t n, size_t elem_bits)
{
	assert(elem_bits <= 32);
	int64_t carry = 0;
	for (size_t i = 0; i < n; ++i) {
		// r in [-(k-1)*(B-1)^2, k*(B-1)^2] in [-1.7e15, 1.7e15]
		int64_t r = round_to_int64(elems[i]);
		r += carry;				// r in [-(k-1)*B^2, k*B^2] in [-1.7e15, 1.7e15]		r >= -(k-1)*(B-1)^2-(k-1)*B >= -(k-1)*B^2
		carry = r >> elem_bits;	// carry = floor(r/B), carry in [-(k-1)*B, k*B]			carry >= floor((-(k-1)*B^2)/B) = -(k-1)*B
		r &= (1 << elem_bits) - 1;
		elems[i] = (double)r;
	}
	return carry;
}

static int64_t mod_0x100000001(int64_t x)
{
	int64_t xl = x & 0xffffffff;
	int64_t xh = x >> 32;
	int64_t M = 0x100000001;
	x = xl - xh;	// x in [-2^31+1, 3*2^31-1]
	if (x < 0) x += M;
	if (x >= M) x -= M;
	return x;
}

/*
prod is a ... ok
prod is b ... ok
prod is a and b ... ok
*/
void mul_dwt_mod_pow2_plus_1(uint32_t prod[], uint32_t& prod_hi, size_t n, const uint32_t a[], size_t an, const uint32_t b[], size_t bn)
{
	// k >= n*32/b;  k is power of 2
	// k*b = n*32;  n*32 mod k = 0
	assert(n >= an && n >= bn);
	if (n == 0)
		return;
	size_t k, elem_bits;
	select_sizes(k, elem_bits, n);	// k*elem_bits >= n*32
	size_t logk = trailing_zeros(k);
	assert(trailing_zeros(n) + 5 >= logk);	// n*32 mod k = 0

	// elem_bits = n*32/k: choose a smaller elem_bits, so that k*elem_bits = n*32
	if (logk >= 5) {
		elem_bits = n >> (logk - 5);
	}
	else {
		elem_bits = n << (5 - logk);
	}

	size_t h = k / 2;	// n>=1, k>=2, h>=1 
	vector<dcomplex> w(h/2);
	for (size_t i = 0; i < h / 2; ++i) {
		double angle = (2 * PI) * ((double)i / (4 * h));
		w[i] = { cos(angle), sin(angle) };
	}

	vector<dcomplex> fa(h);
	extract_elems_dwt(&fa[0], h, a, an, elem_bits);
	weight(fa, w);
	fft_mm_radix_4<false>(&fa[0], h);

	vector<dcomplex> fb(h);
	extract_elems_dwt(&fb[0], h, b, bn, elem_bits);
	weight(fb, w);
	fft_mm_radix_4<false>(&fb[0], h);

	for (size_t i = 0; i < h; ++i)
		fa[i] *= fb[i];

	fft_mm_radix_4<true>(&fa[0], h);
	unweight(fa, w);
	free(w);

	complex_h2v(&fb[0], &fa[0], h);
	free(fa);
	double* d_prod = (double*)&fb[0];

	int64_t carry = adjust_carry_mod(d_prod, k, elem_bits);
	// carry in [-1.7e15, 1.7e15]
	combine_elems(prod, n, d_prod, k, elem_bits);
	if (carry == 0 || carry == 1) {
		prod_hi = (uint32_t)carry;
		return;
	}
	// carry > 1 or carry < 0 
	carry *= -1;
	// carry > 0 or carry < -1
	// add carry to prod
	if (n == 1) {
		// evaluate carry mod 2^32+1
		carry = mod_0x100000001(carry);
		// carry in [0, 2^32]
		uint64_t sum = (uint64_t)prod[0] + (uint64_t)carry;
		prod[0] = (uint32_t)sum;
		prod_hi = sum >> 32;
		return;
	}
	// n>1
	uint32_t carry_little_endian[2];
	if (carry > 0) {
		carry_little_endian[0] = (uint32_t)carry;
		carry_little_endian[1] = (uint32_t)(carry >> 32);
		prod_hi = acc_short(prod, n, carry_little_endian, 2);
	}
	else {
		// add 2^(n*32)+1, subtract -carry, equivalent to add 2^(n*32), subtract -carry-1
		// carry < -1
		carry = -carry - 1;
		// carry > 0
		carry_little_endian[0] = (uint32_t)carry;
		carry_little_endian[1] = (uint32_t)(carry >> 32);
		prod_hi = 1 - sub_short(prod, n, carry_little_endian, 2);
	}
}