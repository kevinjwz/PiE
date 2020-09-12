#pragma once
#include<complex>
#include<intrin.h>

typedef std::complex<double> dcomplex;

inline __m128d dcomplex_to_m128d(dcomplex z)
{
	return _mm_loadu_pd((double*)&z);
}

inline dcomplex m128d_to_dcomplex(__m128d m)
{
	dcomplex z;
	_mm_storeu_pd((double*)&z, m);
	return z;
}

inline __m128d operator+(__m128d x, __m128d y)
{
	return _mm_add_pd(x, y);
}

inline __m128d operator-(__m128d x, __m128d y)
{
	return _mm_sub_pd(x, y);
}

inline const __m128d SIGN_BITS = { -0.0,-0.0 };
inline __m128d operator-(__m128d x)
{
	return _mm_xor_pd(x, SIGN_BITS);
}

inline __m128d operator*(__m128d x, __m128d y)
{
	return _mm_mul_pd(x, y);
}

inline __m128d mul_dcomplex(__m128d x, __m128d y)
{
	// x = [a, b], y = [c, d]
	// x*y = [a*c-b*d, a*d+b*c]
	__m128d a_a = _mm_movedup_pd(x);
	__m128d b_b = _mm_unpackhi_pd(x, x);
	__m128d d_c = _mm_shuffle_pd(y, y, 0b01);
	__m128d c_d = y;
	__m128d ac_ad = _mm_mul_pd(a_a, c_d);
	__m128d bd_bc = _mm_mul_pd(b_b, d_c);
	return _mm_addsub_pd(ac_ad, bd_bc);
}

/*
z0 = x0*y0
z1 = x1*y1
*/
inline void mul_dcomplex_v(__m128d& z0_z1_re, __m128d& z0_z1_im, __m128d x0_x1_re, __m128d x0_x1_im, __m128d y0_y1_re, __m128d y0_y1_im)
{
	z0_z1_re = x0_x1_re * y0_y1_re - x0_x1_im * y0_y1_im;
	z0_z1_im = x0_x1_re * y0_y1_im + x0_x1_im * y0_y1_re;
}

inline __m128d mul_i_sse(__m128d x)
{
	// x = [a, b]
	// x*i = [-b, a]
	__m128d zero = _mm_setzero_pd();
	__m128d b_a = _mm_shuffle_pd(x, x, 0b01);
	return _mm_addsub_pd(zero, b_a);
}

inline __m128d mul_neg_i_sse(__m128d x)
{
	// x = [a, b]
	// x*i = [b, -a]
	__m128d zero = _mm_setzero_pd();
	__m128d a_b = x;
	__m128d nega_b = _mm_addsub_pd(zero, a_b);
	return _mm_shuffle_pd(nega_b, nega_b, 0b01);
}