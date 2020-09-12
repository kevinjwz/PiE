#include <stdint.h>
#include <assert.h>
#include <algorithm>
#include "bigint.h"
#include "micro_helpers.inl"
#include "array_basic_arith.h"
#include "array_helpers.h"
#include "mul_karatsuba.h"
#include "mul_fft_dwt.h"
/*
c = a + ~b
return: carry
c is a ... ok
c is b ... ok
c is a and b ... ok
*/
static uint32_t add_not_n(uint32_t c[], const uint32_t a[], const uint32_t b[], size_t n)
{
	uint64_t carry = 0;
	for (size_t i = 0; i < n; ++i) {
		uint64_t sum = (uint64_t)a[i] + (uint64_t)~b[i] + carry;
		c[i] = (uint32_t)sum;
		carry = sum >> 32;
	}
	return (uint32_t)carry;
}


/*
size(c)=size(a)=size(b)=n
a_hi, b_hi in {0,1}
c_hi in {0,1}

c is a ... ok
c is b ... ok
c is a and b ... ok
*/
static void add_mod_pow2_plus_1(uint32_t c[], uint32_t& c_hi,
								const uint32_t a[], uint32_t a_hi,
								const uint32_t b[], uint32_t b_hi,
								size_t n)
{
	assert(a_hi <= 1 && b_hi <= 1);
	uint32_t ch = add_n(c, a, b, n);
	uint32_t hi_sum = ch + a_hi + b_hi;
	if (hi_sum == 0) {
		c_hi = 0;
	}
	else {
		uint32_t borrow = sub_1(c, hi_sum - 1, n);
		c_hi = 1 - borrow;
	}
}

/*
size(c)=size(a)=size(b)=n
a_hi, b_hi in {0,1}
c_hi in {0,1}

c is a ... ok
c is b ... ok
c is a and b ... ok
*/
static void sub_mod_pow2_plus_1(uint32_t c[], uint32_t& c_hi,
								const uint32_t a[], uint32_t a_hi,
								const uint32_t b[], uint32_t b_hi,
								size_t n)
{
	assert(a_hi <= 1 && b_hi <= 1);
	uint32_t ch = add_not_n(c, a, b, n);
	c_hi = add_1(c, 2 + b_hi - a_hi - ch, n);
}

static void neg_mod_pow2_plus_1(uint32_t a[], uint32_t& a_hi, size_t n)
{
	for (size_t i = 0; i < n; ++i)
		a[i] = ~a[i];
	a_hi = add_1(a, a_hi + 2, n);
}


/*
[src_hi:src] => [hi:lo]
size(lo)=size(hi)=size(src)=n
*/
static void shift_left(uint32_t lo[], uint32_t hi[], const uint32_t src[], uint32_t src_hi, size_t n, size_t shift_words, size_t shift_bits)
{
	assert(shift_words < n);
	assert(shift_bits < 32);
	assert(lo != src && hi != src && lo != hi);
	size_t w = shift_words;
	size_t b = shift_bits;

	if (b == 0) {
		set(lo, 0, w);
		copy(lo + w, src, n - w);

		copy(hi, src + n - w, w);
		hi[w] = src_hi;
		set(hi + w + 1, 0, n - w - 1);
	}
	else {
		set(lo, 0, w);
		copy_shift_left(lo + w, src, n - w, b);

		copy_shift_left(hi, src + n - w, w, b);
		hi[0] |= src[n - w - 1] >> (32 - b);
		hi[w] = (src_hi << b) | (src[n - 1] >> (32 - b));
		set(hi + w + 1, 0, n - w - 1);
	}
}

/*
c = a * 2^k
k in [0, n*32 -1]
*/
static void mul_pow2_mod_pow2_plus_1(uint32_t c[], uint32_t& c_hi,
									 const uint32_t a[], uint32_t a_hi,
									 size_t k, size_t n)
{

	assert(k < n * 32);
	assert(a_hi <= 1);
	size_t w = k / 32;	// shift words: w < n
	size_t b = k % 32;	// shift bits:  b < 32
	uint32_t* a_shift = new uint32_t[2 * n];
	shift_left(a_shift, a_shift + n, a, a_hi, n, w, b);
	sub_mod_pow2_plus_1(c, c_hi, a_shift, 0, a_shift + n, 0, n);
	delete[] a_shift;
}

/*
c = a / 2^k
k in [1, n*32]
*/
static void div_pow2_mod_pow2_plus_1(uint32_t c[], uint32_t& c_hi,
									 const uint32_t a[], uint32_t a_hi,
									 size_t k, size_t n)
{
	assert(k - 1 < n * 32);
	assert(a_hi <= 1);
	size_t w = (n * 32 - k) / 32;
	size_t b = (n * 32 - k) % 32;
	uint32_t* a_shift = new uint32_t[2 * n];
	shift_left(a_shift, a_shift + n, a, a_hi, n, w, b);
	sub_mod_pow2_plus_1(c, c_hi, a_shift + n, 0, a_shift, 0, n);
	delete[] a_shift;
}

/*
n is power of 2
*/
static void mul_pow2_mod_pow2_plus_1_ext(uint32_t c[], uint32_t& c_hi,
										 const uint32_t a[], uint32_t a_hi,
										 int64_t k, size_t n)
{
	/*
	N = n*32
	2^2N = 1
	m = k mod 2N
	a * 2^k = a * 2^m

	m in [0, N-1]:  simple
	m in [N, 2N-1]: a*2^m = a/2^(-m) = a/2^(2N-m)  2N-m in [1, N]
	*/
	assert(is_pow2(n));
	assert(a_hi <= 1);
	size_t m = k & (n * 32 * 2 - 1);
	if (m < n * 32) {
		mul_pow2_mod_pow2_plus_1(c, c_hi, a, a_hi, m, n);
	}
	else {
		div_pow2_mod_pow2_plus_1(c, c_hi, a, a_hi, n * 32 * 2 - m, n);
	}

}

/*
	semi-normalized --> fully normalized
	semi-normalized: ah in {0, 1}
	fully normalized: ah=0 | ah=1, a=0
*/
static void normalize(uint32_t a[], uint32_t& ah, size_t n)
{
	assert(ah <= 1);
	if (ah == 1) {
		uint32_t borrow = sub_1(a, 1, n);
		if (borrow == 1) {
			// original a = 0
			set(a, 0, n);
		}
		else {
			// original a > 0
			ah = 0;
		}
	}
}

/*
mod 2^(n*32)+1
semi-normalized
*/
static void mod_pow2_plus_1(uint32_t a[], uint32_t& ah, size_t n)
{
	if (ah > 0) {
		uint32_t borrow = sub_1(a, ah - 1, n);
		ah = 1 - borrow;
	}
}

/*
size(x)=nx*k
size(xh)=k
size(a)=n*k
nx > n
*/
static void extend(uint32_t x[], uint32_t xh[],
				   const uint32_t a[], uint32_t ah,
				   size_t nx, size_t n, size_t k)
{
	for (size_t i = 0; i < k - 1; ++i) {
		copy(x + i * nx, a + i * n, n);
		set(x + i * nx + n, 0, nx - n);
	}
	copy(x + (k - 1) * nx, a + (k - 1) * n, n);
	x[(k - 1) * nx + n] = ah;
	set(x + (k - 1) * nx + n + 1, 0, nx - n - 1);

	set(xh, 0, k);
}

/*
size(x) = nx * k
size(xh) = k

mul 1, p, p^2, ..., p^(k-1)
*/
static void weight(uint32_t x[], uint32_t xh[], size_t nx, size_t k, int64_t logp)
{
	for (size_t i = 1; i < k; ++i) {
		mul_pow2_mod_pow2_plus_1_ext(x + i * nx, xh[i], x + i * nx, xh[i], logp*i, nx);
	}
}

/*
size(x) = nx * k
element-size is nx
*/
static void bitrev_perm_inplace(uint32_t x[], size_t nx, size_t k)
{
	size_t logk = trailing_zeros(k);
	size_t rev_i = 0;
	for (size_t i = 0; i < k; ++i) {
		if (i < rev_i) {
			std::swap_ranges(x + i * nx, x + (i + 1)*nx, x + rev_i * nx);
		}
		rev_i = rev_inc(rev_i, logk);
	}
}

/*
size(r0) = size(r1) = n
*/
//static void fft_radix_2_kernel(uint32_t r0[], uint32_t& r0h, uint32_t r1[], uint32_t& r1h, size_t n, int64_t logwi)
//{
//	uint32_t* r0_copy = new uint32_t[n];
//	copy(r0_copy, r0, n);
//	uint32_t r0_copy_h = r0h;
//	
//	mul_pow2_mod_pow2_plus_1_ext(r1, r1h, r1, r1h, logwi, n);
//	add_mod_pow2_plus_1(r0, r0h, r0_copy, r0_copy_h, r1, r1h, n);
//	sub_mod_pow2_plus_1(r1, r1h, r0_copy, r0_copy_h, r1, r1h, n);
//	delete[] r0_copy;
//}

/*
size(r0) = size(r1) = n
*/
static void fft_radix_2_kernel(uint32_t r0[], uint32_t& r0h, uint32_t r1[], uint32_t& r1h, size_t n, int64_t logwi)
{
	/*
		r0 = in0 + in1*w
		r1 = in0 - in1*w
	*/
	uint32_t* prod = new uint32_t[n];
	uint32_t prodh;

	mul_pow2_mod_pow2_plus_1_ext(prod, prodh, r1, r1h, logwi, n);
	copy(r1, r0, n);
	r1h = r0h;
	add_mod_pow2_plus_1(r0, r0h, r0, r0h, prod, prodh, n);
	sub_mod_pow2_plus_1(r1, r1h, r1, r1h, prod, prodh, n);
	delete[] prod;
}

/*
size(x) = nx * k
size(xh) = k
w^k = 1
*/
static void fft(uint32_t x[], uint32_t xh[], size_t nx, size_t k, int64_t logw)
{
	assert(is_pow2(k));
	bitrev_perm_inplace(x, nx, k);
	bitrev_perm_inplace(xh, 1, k);

	int64_t expon = logw * k / 2;
	for (size_t n_part = 2; n_part <= k; n_part *= 2, expon /= 2) {
		size_t n_half = n_part / 2;
		for (size_t part_i = 0; part_i != k; part_i += n_part) {
			for (size_t half_i = 0; half_i < n_half; ++half_i) {
				fft_radix_2_kernel(x + (part_i + half_i)*nx, xh[part_i + half_i],
								   x + (part_i + half_i + n_half)*nx, xh[part_i + half_i + n_half],
								   nx, expon * half_i);
			}
		}
	}
}

/*
size(x) = nx * k
size(xh) = k
w^k = 1
*/
static void ifft(uint32_t x[], uint32_t xh[], size_t nx, size_t k, int64_t logw)
{
	fft(x, xh, nx, k, -logw);
	size_t logk = trailing_zeros(k);
	for (size_t i = 0; i < k; ++i) {
		uint32_t* elem = x + i * nx;
		uint32_t& elem_hi = xh[i];
		div_pow2_mod_pow2_plus_1(elem, elem_hi, elem, elem_hi, logk, nx);
	}
}

/*
size(src) = m*k
bits <= 32
*/
static void extract_lowbits(uint32_t lowbits[], size_t lowbits_size, const uint32_t src[], size_t m, size_t k, size_t bits, size_t unit_bits)
{
	assert(bits <= 32);
	set(lowbits, 0, lowbits_size);

	uint32_t mask = bits == 32 ? -1 : ((1 << bits) - 1);
	size_t dst_i = 0;
	size_t bit_offset = 0;
	for (size_t src_i = 0; src_i < m*k; src_i += m) {
		uint32_t val = src[src_i] & mask;
		lowbits[dst_i] |= val << bit_offset;
		if (bit_offset + bits > 32) {
			lowbits[dst_i + 1] |= val >> (32 - bit_offset);
		}
		bit_offset += unit_bits;
		if (bit_offset >= 32) {
			dst_i += bit_offset / 32;
			bit_offset %= 32;
		}
	}
}

/*
size(src) = m*k
bits <= 32
*/
static uint32_t read_bits(const uint32_t lowbits[], size_t i, size_t bits, size_t unit_bits)
{
	assert(bits <= 32);
	size_t word_i = i * unit_bits / 32;
	size_t bit_offset = (i * unit_bits) % 32;
	uint32_t mask = bits == 32 ? -1 : ((1 << bits) - 1);
	if (bit_offset + bits <= 32) {
		return (lowbits[word_i] >> bit_offset) & mask;
	}
	else {
		uint32_t lo = lowbits[word_i] >> bit_offset;
		uint32_t hi = lowbits[word_i + 1] << (32 - bit_offset);
		return (lo | hi) & mask;
	}
}



static const size_t DWT_THRESHOLD = 16;

/*
a *= b (mod 2^(n*32)+1)
input: fully normalized
output: semi-normalized
*/
static void base_mul_dwt(uint32_t a[], uint32_t& a_hi,
						 const uint32_t b[], uint32_t b_hi,
						 size_t n)
{
	if (b_hi == 1) {
		// b=-1
		neg_mod_pow2_plus_1(a, a_hi, n);
	}
	else if (a_hi == 1) {
		// a=-1
		copy(a, b, n);
		a_hi = b_hi;
		neg_mod_pow2_plus_1(a, a_hi, n);
	}
	else {
		// a, b in [0, 2^(n*32)-1]
		mul_dwt_mod_pow2_plus_1(a, a_hi, n, a, n, b, n);
	}
}

/*
a *= b (mod 2^(n*32)+1)
input: fully normalized
output: semi-normalized
*/
static void base_mul_krtb(uint32_t a[], uint32_t& a_hi,
						  const uint32_t b[], uint32_t b_hi,
						  size_t n)
{
	if (b_hi == 1) {
		// b=-1
		neg_mod_pow2_plus_1(a, a_hi, n);
	}
	else if (a_hi == 1) {
		// a=-1
		copy(a, b, n);
		a_hi = b_hi;
		neg_mod_pow2_plus_1(a, a_hi, n);
	}
	else {
		// a, b in [0, 2^(n*32)-1]
		uint32_t* prod = new uint32_t[2 * n];
		mul_karatsuba_p2_destroy_a(prod, a, b, n);
		copy(a, prod, n);
		sub_mod_pow2_plus_1(a, a_hi, a, a_hi, prod + n, 0, n);
		delete[] prod;
	}
}

typedef void(*func_mul_mod_pow2_plus_1)(uint32_t a[], uint32_t& a_hi,
										const uint32_t b[], uint32_t b_hi,
										size_t n);

/*
a *= b (mod 2^(n*32)+1)
input: fully normalized
output: semi-normalized
*/
template <size_t base_threshold, func_mul_mod_pow2_plus_1 base_mul>
static void mul_ssa_mod_pow2_plus_1(uint32_t a[], uint32_t& a_hi,
									const uint32_t b[], uint32_t b_hi,
									size_t n)
{
	assert(is_pow2(n));
	assert(a_hi <= 1 && b_hi <= 1);

	if (n <= base_threshold) {
		base_mul(a, a_hi, b, b_hi, n);
		return;
	}
	assert(n >= 16);	// n < 16 ==> m == 0

	// calculate sizes
	size_t nbits = n * 32;
	assert(nbits > 0);
	size_t log2_nbits = trailing_zeros(nbits);
	size_t log2_k = log2_nbits / 2;
	size_t log2_mbits = log2_nbits - log2_k;
	size_t log2_nxbits = log2_mbits + 1;

	size_t k = 1ULL << log2_k;					// fft parts
	size_t m = (1ULL << log2_mbits) / 32;	// words per part
	assert(m > 0);
	size_t nx = m * 2;

	// if large n makes k > 2^31, REWRITE this SSA program.
	assert(k <= 1ull << 31);

	int64_t log2_p = 1ULL << (log2_nxbits - log2_k);	// nxbits / k
	int64_t log2_w = log2_p * 2;

	// calculate ci mod 2^(2*m*32)+1

	uint32_t* ax = new uint32_t[nx*k];
	uint32_t* axh = new uint32_t[k];
	extend(ax, axh, a, a_hi, nx, m, k);
	weight(ax, axh, nx, k, log2_p);
	fft(ax, axh, nx, k, log2_w);

	uint32_t* bx = new uint32_t[nx*k];
	uint32_t* bxh = new uint32_t[k];
	extend(bx, bxh, b, b_hi, nx, m, k);
	weight(bx, bxh, nx, k, log2_p);
	fft(bx, bxh, nx, k, log2_w);

	for (size_t i = 0; i < k; ++i) {
		uint32_t* a_part = ax + i * nx;
		uint32_t& a_part_hi = axh[i];
		uint32_t* b_part = bx + i * nx;
		uint32_t& b_part_hi = bxh[i];
		normalize(a_part, a_part_hi, nx);
		normalize(b_part, b_part_hi, nx);
		mul_ssa_mod_pow2_plus_1<base_threshold, base_mul>(a_part, a_part_hi, b_part, b_part_hi, nx);
	}
	delete[] bx;
	delete[] bxh;

	ifft(ax, axh, nx, k, log2_w);
	weight(ax, axh, nx, k, -log2_p);

	// calculate ci mod k
	size_t lowbits_size = (k * 3 * log2_k + 31) / 32;	// >= k * 3logk bits
	size_t unit_bits = 3 * log2_k;	// k < 2^64, log2_k < 64, unit_bits < 192
	uint32_t* a_lowbits = new uint32_t[lowbits_size];
	extract_lowbits(a_lowbits, lowbits_size, a, m, k, log2_k, unit_bits);

	uint32_t* b_lowbits = new uint32_t[lowbits_size];
	extract_lowbits(b_lowbits, lowbits_size, b, m, k, log2_k, unit_bits);
	uint32_t* prod = new uint32_t[2 * lowbits_size];
	mul_karatsuba(prod, a_lowbits, b_lowbits, lowbits_size);
	delete[] a_lowbits;
	delete[] b_lowbits;

	// calculate ci mod k*(B^2+1), and sum of ci
	set(a, 0, n);
	uint32_t sum_hi = 0;
	uint32_t borrow = 0;
	std::vector<size_t> sub_indices;
	sub_indices.reserve(k);
	for (size_t i = 0; i < k; ++i) {
		uint32_t* a_part = ax + i * nx;
		uint32_t& a_part_hi = axh[i];
		uint32_t ci_mod_k = read_bits(prod, i, log2_k, unit_bits);
		uint32_t ciplusk_mod_k = read_bits(prod, i + k, log2_k, unit_bits);
		ci_mod_k = ci_mod_k - ciplusk_mod_k;
		normalize(a_part, a_part_hi, nx);
		uint32_t diff_mod_k = (ci_mod_k - a_part[0]) & (k - 1);
		a_part_hi += add_1(a_part, diff_mod_k, nx);
		a_part_hi += diff_mod_k;
		if (a_part_hi >= i + 1) {
			// delay subtracting k*(B^2+1)*B^i
			// to avoid frequent sign reversals leading to O(n^2) complexity
			sub_indices.push_back(i);
			//printf("k=%zd i=%zd\n", k, i);
		}
		if (i <= k - 3) {
			// sum <= B^k
			sum_hi += acc_short(a + i * m, n - i * m, a_part, nx);
			sum_hi += add_1(a + i * m + nx, a_part_hi, n - i * m - nx);
		}
		else if (i == k - 2) {
			// sum <= (k+1)*B^k
			sum_hi += acc_short(a + i * m, n - i * m, a_part, nx); // n-(k-2)*m = 2*m,  nx = 2*m
			sum_hi += a_part_hi;
		}
		else {
			// i == k-1
			// sum <= (k+2)*B^k
			sum_hi += acc_short(a + i * m, n - i * m, a_part, m); // n - (k-1) * m = m
			borrow += sub_short(a, n, a_part + m, m);
			borrow += sub_1(a + m, a_part_hi, n - m);
		}
	}
	delete[] ax;
	delete[] axh;
	delete[] prod;

	// delayed subtraction
	for (size_t i : sub_indices) {
		// max(i) = k-2
		if (i < k - 2) {
			// subtract k*B^i + k*B^(i+2)
			// borrow <= 1
			borrow += sub_1(a + i * m, (uint32_t)k, n - i * m);
			borrow += sub_1(a + (i + 2) * m, (uint32_t)k, n - (i + 2) * m);
		}
		else {
			// i == k-2
			borrow += sub_1(a + i * m, (uint32_t)k, n - i * m);
			sum_hi += add_1(a, (uint32_t)k, n);
		}
	}
	// sum_hi <= k+3
	// borrow <= 3

	if (sum_hi >= borrow) {
		a_hi = sum_hi - borrow;
		mod_pow2_plus_1(a, a_hi, n);
	}
	else {
		// borrow > sum_hi
		a_hi = add_1(a, borrow - sum_hi, n);
	}
}

/*
	n >= src.size()
*/
static void zero_extend(bigint& dst, const bigint_src& src, size_t n)
{
	dst.clear();
	dst.reserve(n);
	dst.insert(dst.end(), src.begin(), src.end());
	dst.insert(dst.end(), n - src.size(), 0);
}

static const size_t KRTB_THRESHOLD = 1024;

template <size_t base_threshold, func_mul_mod_pow2_plus_1 base_mul>
void mul_ssa(bigint& c, const bigint_src& a, const bigint_src& b)
{
	assert(!is_from(a, c) && !is_from(b, c));
	size_t an = a.size(), bn = b.size();
	if (an == 0 || bn == 0) {
		c.clear();
		return;
	}
	assert(an + bn >= an);	// no overflow
	size_t n = clp2(an + bn);
	assert(n > 0);			// clp2(an+bn)<=SIZE_MAX
	n = std::max<size_t>(n, 16);

	zero_extend(c, a, n);
	bigint b_copy;
	zero_extend(b_copy, b, n);

	uint32_t c_hi = 0;
	mul_ssa_mod_pow2_plus_1<base_threshold, base_mul>(&c[0], c_hi, &b_copy[0], 0, n);
	normalize(&c[0], c_hi, n);
	// c_hi must be 0
	// c[an+bn .. n-1] must be 0
	if (c[an + bn - 1] == 0)
		c.resize(an + bn - 1);
	else
		c.resize(an + bn);
}

void mul_ssa_fast(bigint& c, const bigint_src& a, const bigint_src& b)
{
	mul_ssa< (1 << 16), base_mul_dwt>(c, a, b);
}

void mul_ssa_pure_int(bigint& c, const bigint_src& a, const bigint_src& b)
{
	mul_ssa< (1 << 12), base_mul_krtb>(c, a, b);
}

// tests
/*
#include <random>
static std::mt19937 rng;

bool test_fft(size_t nx, size_t k)
{
	bigint a;
	rand_bigint(a, nx*k);
	bigint ah(k);
	for (size_t i = 0; i < k; ++i) {
		ah[i] = rand() % 2;
		normalize(&a[i * nx], ah[i], nx);
	}

	bigint a_copy = a;
	bigint ah_copy = ah;
	int64_t logw = nx * 32 * 2 / k;
	fft(&a[0], &ah[0], nx, k, logw);
	ifft(&a[0], &ah[0], nx, k, logw);

	for (size_t i = 0; i < k; ++i) {
		normalize(&a[i * nx], ah[i], nx);
	}

	return a == a_copy && ah == ah_copy;
}

bool mod_equal(bigint a, bigint b)
{
	size_t n = a.size() - 1;
	assert(a.size() == b.size());
	assert(a[n] <= 1 && b[n] <= 1);
	normalize(&a[0], a[n], n);
	normalize(&b[0], b[n], n);
	return a == b;
}

void mul_pow2_krtb(uint32_t c[], uint32_t& c_hi,
				   const uint32_t a[], uint32_t a_hi,
				   size_t k, size_t n)
{
	assert(k < n * 32);
	assert(a_hi <= 1);
	size_t w = k / 32;	// shift words: w < n
	size_t b = k % 32;	// shift bits:  b < 32

	bigint pow2(n);
	pow2[w] = 1 << b;

	bigint a_copy(a, a + n);
	normalize(&a_copy[0], a_hi, n);
	if (a_hi == 1) {
		copy(c, &pow2[0], n);
		c_hi = 0;
		neg_mod_pow2_plus_1(c, c_hi, n);
	}
	else {
		bigint prod(2 * n);
		mul_karatsuba(&prod[0], &a_copy[0], &pow2[0], n);
		sub_mod_pow2_plus_1(c, c_hi, &prod[0], 0, &prod[n], 0, n);
	}

}

void div_pow2_krtb(uint32_t c[], uint32_t& c_hi,
				   const uint32_t a[], uint32_t a_hi,
				   size_t k, size_t n)
{
	assert(k - 1 < n * 32);
	assert(a_hi <= 1);
	mul_pow2_krtb(c, c_hi, a, a_hi, n * 32 - k, n);
	neg_mod_pow2_plus_1(c, c_hi, n);
}

// PASS!
bool test_mul_pow2(size_t n)
{
	bigint a(n + 1);
	rand_bigint(a, n + 1);
	a[n] &= 1;

	size_t k = rng() % (n * 32);
	bigint prod(n + 1), prod_gt(n + 1);
	mul_pow2_mod_pow2_plus_1(&prod[0], prod[n], &a[0], a[n], k, n);
	mul_pow2_krtb(&prod_gt[0], prod_gt[n], &a[0], a[n], k, n);
	return mod_equal(prod, prod_gt);
}

// PASS!
bool test_div_pow2(size_t n)
{
	bigint a(n + 1);
	rand_bigint(a, n + 1);
	a[n] &= 1;

	size_t k = rng() % (n * 32) + 1;
	bigint prod(n + 1), prod_gt(n + 1);
	div_pow2_mod_pow2_plus_1(&prod[0], prod[n], &a[0], a[n], k, n);
	div_pow2_krtb(&prod_gt[0], prod_gt[n], &a[0], a[n], k, n);

	return mod_equal(prod, prod_gt);
}

bool test_mul_div_pow2(size_t n)
{
	bigint a(n + 1);
	rand_bigint(a, n + 1);
	a[n] &= 1;

	bigint a_copy = a;

	size_t k = rng() % (n * 32 - 1) + 1;
	div_pow2_mod_pow2_plus_1(&a[0], a[n], &a[0], a[n], k, n);
	mul_pow2_mod_pow2_plus_1(&a[0], a[n], &a[0], a[n], k, n);

	return mod_equal(a, a_copy);
}

bool test_mul_pow2_ext(size_t n)
{
	bigint a(n + 1);
	rand_bigint(a, n + 1);
	a[n] &= 1;

	bigint a_copy = a;

	size_t k = rng() % (n * 32 - 1) + 1;
	div_pow2_mod_pow2_plus_1(&a[0], a[n], &a[0], a[n], k, n);
	mul_pow2_mod_pow2_plus_1(&a[0], a[n], &a[0], a[n], k, n);

	return mod_equal(a, a_copy);
}

int mainee()
{
	//assert(test_fft(1, 2));
	while (true) {
		//assert(test_mul_div_pow2(4096));
		assert(test_fft(32, 32));
		printf("PASS!\n");

	}
}
*/