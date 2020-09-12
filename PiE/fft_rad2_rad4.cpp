#include<complex>
#include<math.h>
#include<assert.h>
#include "micro_helpers.inl"
#include "dcomplex.inl"
#include "dcomplex_sse.inl"
#include "fft_constants.h"

typedef std::complex<double> dcomplex;


/*
size(w) = 2n
w[0] = undefined
w[1] = w_p^0
w[2..3] = w_2p^0, w_2p^1
w[4..7] = w_4p^0, w_4p^1, w_4p^2, w_4p^3
...
w[n..2n-1] = w_np^0, w_np^1, w_np^2, ...,  w_np^(n-1)
*/
template<size_t p, int dir>
static void compute_w_mipmap(dcomplex w[], size_t n)
{
	assert(is_pow2(n));	// n is power of 2
	w[1] = 1;
	for (size_t n_layer = 2; n_layer <= n; n_layer *= 2) {
		dcomplex* layer = &w[n_layer];
		for (size_t i = 0; i < n_layer; ++i) {
			if (i & 1) {
				double angle = dir * (2 * PI) * i / (n_layer * p);
				layer[i] = { cos(angle), sin(angle) };
			}
			else {
				layer[i] = w[n_layer / 2 + i / 2];
			}
		}
	}
}

static dcomplex w_mipmap_lookup(dcomplex w[], size_t i, size_t n)
{
	return w[n + i];
}


/*
size(w) = n
w[0] = 1
w[1] = w_2p
w[2..3] = w_4p, w_4p^3
w[4..7] = w_8p, w_8p^3, w_8p^5, w_8p^7
...
w[n/2..n-1] = w_np, w_np^3, w_np^5, ...,  w_np^(n-1)
*/

template<size_t p, int dir>
static void compute_w_mipmap_reuse(dcomplex w[], size_t n)
{
	assert(is_pow2(n));	// n is power of 2
	w[0] = 1;
	for (size_t n_layer = 1; n_layer <= n / 2; n_layer *= 2) {
		dcomplex* layer = &w[n_layer];
		for (size_t i = 0; i < n_layer; ++i) {
			double angle = dir * (2 * PI) * (i * 2 + 1) / (n_layer * 2 * p);
			layer[i] = { cos(angle), sin(angle) };
		}
	}
}

[[maybe_unused]]
static dcomplex w_mipmap_reuse_lookup(dcomplex w[], size_t i, size_t n)
{
	size_t index = n + i;
	index >>= 1 + trailing_zeros(index);
	return w[index];
}


/*
i in [0, 4n-1]
*/
template<bool inv, dcomplex lookup(dcomplex w[], size_t i, size_t n)>
static dcomplex lookup_ext4(dcomplex w[], size_t i, size_t n)
{
	/*
	if (i < n) {
	return w_mipmap_lookup(w, i, n);
	}
	else if (i < 2 * n) {
	return dir * mul_i(w_mipmap_lookup(w, i - n, n));
	}
	else if (i < 3 * n) {
	return -w_mipmap_lookup(w, i - 2 * n, n);
	}
	else {
	return -dir * mul_i(w_mipmap_lookup(w, i - 3 * n, n));
	}
	*/
	constexpr double dir = inv ? 1 : -1;
	if (i < 2 * n) {
		if (i < n)
			return lookup(w, i, n);
		else
			return dir * mul_i(lookup(w, i - n, n));
	}
	else {
		if (i < 3 * n)
			return -lookup(w, i - 2 * n, n);
		else
			return -dir * mul_i(lookup(w, i - 3 * n, n));
	}

}

template<bool inv, dcomplex lookup(dcomplex w[], size_t i, size_t n)>
static __m128d lookup_ext4_sse(dcomplex w[], size_t i, size_t n)
{
	constexpr double dir = inv ? 1 : -1;
	if (i < 2 * n) {
		if (i < n) {
			return dcomplex_to_m128d(lookup(w, i, n));
		}
		else {
			__m128d m = dcomplex_to_m128d(lookup(w, i - n, n));
			return dir > 0 ? mul_i_sse(m) : mul_neg_i_sse(m);
		}
	}
	else {
		if (i < 3 * n) {
			return -dcomplex_to_m128d(lookup(w, i - 2 * n, n));
		}
		else {
			__m128d m = dcomplex_to_m128d(lookup(w, i - 3 * n, n));
			return dir > 0 ? mul_neg_i_sse(m) : mul_i_sse(m);
		}
	}
}



/*
i in [0, 8n-1]
*/
template<bool inv, dcomplex lookup(dcomplex w[], size_t i, size_t n)>
static dcomplex lookup_ext8(dcomplex w[], size_t i, size_t n)
{
	constexpr double dir = inv ? 1 : -1;
	size_t q = i >> trailing_zeros(n);
	size_t s = (i&(2 * n - 1)) == n ? 8 : q;
	size_t mask = n - 1;
	switch (s) {
		case 0:return lookup(w, i, n);
		case 1:return dir * reflect(lookup(w, (-i)&mask, n));
		case 2:return dir * mul_i(lookup(w, i&mask, n));
		case 3:return mul_i(reflect(lookup(w, (-i)&mask, n)));
		case 4:return -lookup(w, i&mask, n);
		case 5:return -dir * reflect(lookup(w, (-i)&mask, n));
		case 6:return -dir * mul_i(lookup(w, i&mask, n));
		case 7:return -mul_i(reflect(lookup(w, (-i)&mask, n)));
		default:
		{
			// i is odd multiple of n
			constexpr dcomplex w8 = { SQRT2 / 2, dir * SQRT2 / 2 };
			constexpr dcomplex w8_pow_3 = { -SQRT2 / 2, dir * SQRT2 / 2 };
			constexpr dcomplex w8_pow_5 = { -SQRT2 / 2, -dir * SQRT2 / 2 };
			constexpr dcomplex w8_pow_7 = { SQRT2 / 2, -dir * SQRT2 / 2 };
			static constexpr dcomplex w_diag[4] = { w8, w8_pow_3, w8_pow_5, w8_pow_7 };
			return w_diag[q / 2];
		}
	}
}


static void load_blocks(dcomplex buffer[], const dcomplex x[], size_t mid_bits, unsigned mid_len, unsigned pfx_len)
{
	size_t pfx = 0;
	size_t pfx_rev = 0;
	size_t pfx_lo_bit = 1ull << (mid_len + pfx_len);
	size_t pfx_end = 1ull << (mid_len + pfx_len * 2);
	while (pfx < pfx_end) {
		for (size_t sfx = 0; sfx < (1ull << pfx_len); ++sfx) {
			buffer[pfx_rev + sfx] = x[pfx + mid_bits + sfx];
		}
		pfx += pfx_lo_bit;
		pfx_rev = rev_inc(pfx_rev, pfx_len * 2);
	}
}

static void store_blocks_mul_coef(dcomplex x[], const dcomplex buffer[], size_t mid_bits_rev, unsigned mid_len, unsigned pfx_len, double coef)
{
	size_t sfx = 0;
	size_t sfx_rev = 0;
	unsigned idx_len = mid_len + pfx_len * 2;
	size_t sfx_end = 1ull << pfx_len;
	while (sfx < sfx_end) {
		size_t pfx_rev = 0;
		size_t pfx_rev_shift = 0;
		while (pfx_rev < (1ull << pfx_len)) {
			x[sfx_rev + mid_bits_rev + pfx_rev] = buffer[pfx_rev_shift + sfx] * coef;
			pfx_rev += 1;
			pfx_rev_shift += 1ull << pfx_len;
		}
		sfx += 1;
		sfx_rev = rev_inc(sfx_rev, idx_len);
	}
}

/*
cache optimal bit-reversal permutation
*/
static const unsigned COBRA_Q = 7;
alignas(64) static dcomplex cobra_buffer1[1 << (2 * COBRA_Q)];
alignas(64) static dcomplex cobra_buffer2[1 << (2 * COBRA_Q)];

static void cobra_in_place(dcomplex x[], size_t n, double coef)
{
	assert(is_pow2(n));

	unsigned idx_len = trailing_zeros(n);
	unsigned pfx_len = std::min(COBRA_Q, idx_len / 2);
	unsigned mid_len = idx_len - pfx_len * 2;
	size_t mid_bits_end = 1ull << (mid_len + pfx_len);
	size_t mid_bits = 0;
	size_t mid_bits_rev = 0;
	while (mid_bits < mid_bits_end) {
		if (mid_bits < mid_bits_rev) {
			load_blocks(cobra_buffer1, x, mid_bits, mid_len, pfx_len);
			load_blocks(cobra_buffer2, x, mid_bits_rev, mid_len, pfx_len);
			store_blocks_mul_coef(x, cobra_buffer2, mid_bits, mid_len, pfx_len, coef);
			store_blocks_mul_coef(x, cobra_buffer1, mid_bits_rev, mid_len, pfx_len, coef);
		}
		else if (mid_bits == mid_bits_rev) {
			load_blocks(cobra_buffer1, x, mid_bits, mid_len, pfx_len);
			store_blocks_mul_coef(x, cobra_buffer1, mid_bits, mid_len, pfx_len, coef);
		}
		mid_bits += 1ull << pfx_len;
		mid_bits_rev = rev_inc(mid_bits_rev, mid_len + pfx_len);
	}
}

static void fft_radix_2_kernel(dcomplex& r0, dcomplex& r1, dcomplex wi)
{
	dcomplex in0 = r0, in1 = r1;
	r0 = in0 + in1 * wi;
	r1 = in0 - in1 * wi;
}

template<bool inverse>
void fft_mm_radix_2(dcomplex x[], size_t n)
{
	assert(is_pow2(n));	// n is power of 2
	if (n == 1)
		return;

	double coef = inverse ? 1.0 / n : 1;
	double dir = inverse ? 1 : -1;

	cobra_in_place(x, n, coef);

	if (n >= 2) {
		for (dcomplex* part = x; part != x + n; part += 2) {
			fft_radix_2_kernel(part[0], part[1], 1.0);
		}
	}
	if (n < 4)
		return;

	dcomplex* w_mm = new dcomplex[n / 2];
	compute_w_mipmap<4, inverse ? 1 : -1>(w_mm, n / 4);

	for (size_t n_part = 4; n_part <= n; n_part *= 2) {
		size_t n_half = n_part / 2;
		size_t n_qt = n_part / 4;
		for (dcomplex* part = x; part != x + n; part += n_part) {
			size_t half_i = 0;
			for (; half_i < n_qt; ++half_i) {
				dcomplex w = w_mipmap_lookup(w_mm, half_i, n_qt);
				fft_radix_2_kernel(part[half_i], part[half_i + n_half], w);
			}
			for (; half_i < n_half; ++half_i) {
				dcomplex w = dir * mul_i(w_mipmap_lookup(w_mm, half_i - n_qt, n_qt));
				fft_radix_2_kernel(part[half_i], part[half_i + n_half], w);
			}
		}
	}
	delete[] w_mm;
}

template
void fft_mm_radix_2<false>(dcomplex x[], size_t n);

template
void fft_mm_radix_2<true>(dcomplex x[], size_t n);


template<bool inverse>
static void fft_radix_4_kernel(dcomplex& r0, dcomplex& r1, dcomplex& r2, dcomplex& r3, dcomplex w, dcomplex w2, dcomplex w3)
{
	dcomplex in0 = r0;
	// swap in1 and in2
	dcomplex in2 = r1;
	dcomplex in1 = r2;
	dcomplex in3 = r3;

	//dcomplex w2 = w * w;
	//dcomplex w3 = w * w2;

	dcomplex in1_mul_wi = in1 * w;
	dcomplex in2_mul_w2i = in2 * w2;
	dcomplex in3_mul_w3i = in3 * w3;

	double dir = inverse ? 1 : -1;

	r0 = (in0 + in2_mul_w2i) + (in1_mul_wi + in3_mul_w3i);

	// in0 + in1 * wi *(dir * 1.0i) + in2 * w2i*(-1.0) + in3 * w3i*(-dir*1.0i);
	r1 = (in0 - in2_mul_w2i) + dir * mul_i(in1_mul_wi - in3_mul_w3i);

	// in0 + in1 * wi *(-1.0) + in2 * w2i + in3 * w3i*(-1.0);
	r2 = (in0 + in2_mul_w2i) - (in1_mul_wi + in3_mul_w3i);

	// in0 + in1 * wi *(-dir * 1.0i) + in2 * w2i*(-1.0) + in3 * w3i*(dir * 1.0i);
	r3 = (in0 - in2_mul_w2i) - dir * mul_i(in1_mul_wi - in3_mul_w3i);
}

template<bool inverse>
static void fft_radix_4_kernel_sse(dcomplex& r0, dcomplex& r1, dcomplex& r2, dcomplex& r3, __m128d w, __m128d w2, __m128d w3)
{
	__m128d in0 = dcomplex_to_m128d(r0);
	// swap in1 and in2
	__m128d in2 = dcomplex_to_m128d(r1);
	__m128d in1 = dcomplex_to_m128d(r2);
	__m128d in3 = dcomplex_to_m128d(r3);

	//dcomplex w2 = w * w;
	//dcomplex w3 = w * w2;

	__m128d in1_mul_wi = mul_dcomplex(in1, w);
	__m128d in2_mul_w2i = mul_dcomplex(in2, w2);
	__m128d in3_mul_w3i = mul_dcomplex(in3, w3);

	__m128d r0m, r1m, r2m, r3m;
	double dir = inverse ? 1 : -1;
	if (dir > 0) {
		r0m = (in0 + in2_mul_w2i) + (in1_mul_wi + in3_mul_w3i);
		r1m = (in0 - in2_mul_w2i) + mul_i_sse(in1_mul_wi - in3_mul_w3i);
		r2m = (in0 + in2_mul_w2i) - (in1_mul_wi + in3_mul_w3i);
		r3m = (in0 - in2_mul_w2i) - mul_i_sse(in1_mul_wi - in3_mul_w3i);
	}
	else {
		r0m = (in0 + in2_mul_w2i) + (in1_mul_wi + in3_mul_w3i);
		r1m = (in0 - in2_mul_w2i) - mul_i_sse(in1_mul_wi - in3_mul_w3i);
		r2m = (in0 + in2_mul_w2i) - (in1_mul_wi + in3_mul_w3i);
		r3m = (in0 - in2_mul_w2i) + mul_i_sse(in1_mul_wi - in3_mul_w3i);
	}
	r0 = m128d_to_dcomplex(r0m);
	r1 = m128d_to_dcomplex(r1m);
	r2 = m128d_to_dcomplex(r2m);
	r3 = m128d_to_dcomplex(r3m);
}


/*
fastest
*/
template<bool inverse>
void fft_mm_radix_4(dcomplex x[], size_t n)
{
	assert(is_pow2(n));	// n is power of 2
	if (n == 1)
		return;

	double coef = inverse ? 1.0 / n : 1.0;

	cobra_in_place(x, n, coef);

	size_t n_part = 4;
	size_t bits = trailing_zeros(n);
	if (bits & 1) {
		// n is not power of 4
		// do simple radix-2 first
		for (dcomplex* part = x; part != x + n; part += 2) {
			fft_radix_2_kernel(part[0], part[1], 1.0);
		}
		n_part *= 2;
	}
	else {
		// n is power of 4
		// do simple radix-4 first
		for (dcomplex* part = x; part != x + n; part += 4) {
			fft_radix_4_kernel<inverse>(part[0], part[1], part[2], part[3], 1.0, 1.0, 1.0);
		}
		n_part *= 4;
	}
	if (n_part > n)
		return;
	// n > 4

	dcomplex* w_mm = new dcomplex[n / 2];
	compute_w_mipmap<4, inverse ? 1 : -1>(w_mm, n / 4);

	for (; n_part <= n; n_part *= 4) {
		size_t n_qt = n_part / 4;
		for (dcomplex* part = x; part != x + n; part += n_part) {
			size_t qt_i = 0;
			for (; qt_i < n_qt; qt_i += 1) {
				__m128d w = dcomplex_to_m128d(w_mipmap_lookup(w_mm, qt_i, n_qt));
				__m128d w2 = lookup_ext4_sse<inverse, w_mipmap_lookup>(w_mm, 2 * qt_i, n_qt);
				__m128d w3 = lookup_ext4_sse<inverse, w_mipmap_lookup>(w_mm, 3 * qt_i, n_qt);
				fft_radix_4_kernel_sse<inverse>(part[qt_i], part[qt_i + n_qt], part[qt_i + n_qt * 2], part[qt_i + n_qt * 3], w, w2, w3);
			}
		}
	}

	delete[] w_mm;
}


template
void fft_mm_radix_4<false>(dcomplex x[], size_t n);

template
void fft_mm_radix_4<true>(dcomplex x[], size_t n);