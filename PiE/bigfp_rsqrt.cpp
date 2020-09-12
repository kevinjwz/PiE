#include<assert.h>
#include<math.h>
#include"bigint.h"
#include"bigfp.h"
#include"int128.h"
#include"micro_helpers.inl"

static bool is_ge_1_lt_B2(const bigfp_src& a)
{
	return a.sign == 1 && (a.exp == 1 || a.exp == 2);
}


void rsqrt_newton_iterate(bigfp& r, const bigfp_src& r0, const bigfp_src& s, size_t p)
{
	bigfp r0_sq;
	mul(r0_sq, r0, r0);

	bigfp s_mul_r0_sq;
	mul(s_mul_r0_sq, s, r0_sq);
	free(r0_sq);
	round(s_mul_r0_sq, -(int128)p);


	bigfp one(1);
	bigfp e0;
	sub(e0, one, s_mul_r0_sq);
	free(s_mul_r0_sq);
	mul_pow2(e0, -1);
	round(e0, -(int128)p);

	bigfp r0_mul_e0;
	mul(r0_mul_e0, r0, e0);
	free(e0);
	round(r0_mul_e0, -(int128)p);

	add(r, r0, r0_mul_e0);
}


void rsqrt_newton_base(bigfp& r, const bigfp_src& s, size_t p)
{
	if (p == 1) {
		double sd = to_double(s);
		from_double(r, 1 / sqrt(sd));
		round(r, -2);
		return;
	}

	// p=2 or 3

	// s = s_norm * 2^(2*e)
	// s_norm in [1, 4)
	size_t lz = leading_zeros(*s.data.rbegin());	// [0, 31]
	int e = (31 - (int)lz) / 2;	// [0, 15]

	double sd_norm = ldexp(to_double(s), -e * 2);
	// sd_norm in [1, 4]
	// |sd_norm - s_norm| <= 0.5*2^(-51)

	double sqrt_norm = sqrt(sd_norm);
	// sqrt_norm in [1, 2]
	// |sqrt_norm - s_norm^0.5| <= 0.25*2^(-51) + 0.5*2^(-52) = 0.5*2^(-51)

	bigfp r_norm;
	from_double(r_norm, 1 / sqrt_norm);
	// r_norm in [0.5, 1]
	// |r_norm - s_norm^-0.5| <= 0.5*2^(-51) + 0.5*2^(-53) <= 0.5*2^(-50)
	// precision needed: 0.5*2^(p*32 - e)

	size_t prec_bits_needed = p * 32 - (size_t)e;		// [49, 96]
	if (prec_bits_needed > 50) {
		bigfp s_round;
		round(s_round, s, -(int128)(p + 1));
		bigfp sr_norm;
		mul_pow2(sr_norm, s_round, -e * 2);
		round(sr_norm, -(int128)(p + 1));
		// sr_norm in [1, 4]
		// |sr_norm - s_norm| <= 0.5*2^(-p*32-31)
		rsqrt_newton_iterate(r_norm, bigfp(r_norm), sr_norm, p + 1);
		// |r_norm - s_norm^-0.5| <= 0.5*2^(-p*32)
	}
	mul_pow2(r, r_norm, -e);
	round(r, -(int128)(p + 1));
}


void rsqrt_newton_rec(bigfp& r, const bigfp_src& s, size_t p)
{
	assert(p <= SIZE_MAX - 1);
	if (p <= 3) {
		rsqrt_newton_base(r, s, p);
		return;
	}

	bigfp r0;
	size_t p0 = (p + 3) / 2;	// p0 = ceil((p+2)/2)
	rsqrt_newton_rec(r0, s, p0);

	bigfp s_round;
	round(s_round, s, -(int128)(p + 1));
	rsqrt_newton_iterate(r, r0, s_round, p + 1);
}


void rsqrt_newton(bigfp& r, const bigfp_src& s, size_t p)
{
	assert(p > 0);
	assert(s.sign == 1);

	/*
	normalize s:
	s = s_norm * B^(2*e)
	s_norm in [1, B^2)
	*/

	bigfp_src s_norm(s);
	s_norm.exp = s.exp % 2 == 0 ? 2 : 1;
	int64_t e = s.exp % 2 == 0 ?
		s.exp / 2 - 1 :
		floor_div(s.exp, 2);	// e = (s.exp - s_norm.exp)/2


	assert(is_ge_1_lt_B2(s_norm));
	rsqrt_newton_rec(r, s_norm, p);
	// r_norm in [1/B, 1]
	// r = r_norm / B^e
	r.exp -= e;
}


double validate_rsqrt(const bigfp_src& r, const bigfp_src& s, size_t p)
{
	bigfp r_sq;
	mul(r_sq, r, r);
	bigfp r_sq_mul_s;
	mul(r_sq_mul_s, r_sq, s);
	bigfp one(1);
	bigfp diff;
	sub(diff, r_sq_mul_s, one);

	int64_t e = s.exp % 2 == 0 ? s.exp / 2 - 1 : floor_div(s.exp, 2);
	bigfp r_mul_s;
	mul(r_mul_s, r, s);
	bigfp upper_bound;
	mul_pow2(upper_bound, r_mul_s, -((int64_t)p + e) * 32);

	diff.sign = abs(diff.sign);
	//return is_less_or_equal(diff, upper_bound);
	return div_apxm(diff, upper_bound);
}

void sqrt(bigfp& r, const bigfp_src& s, size_t p)
{
	if (s.sign == 0) {
		r.sign = 0;
		r.data.clear();
		return;
	}

	bigfp rsqrt;
	rsqrt_newton(rsqrt, s, p + 2);
	bigfp s_round;
	round(s_round, s, (int128)s.exp - 2 - p);
	mul(r, s_round, rsqrt);
	round(r, (int128)r.exp - p - 1);
}

void sqrt(bigfp& r, size_t p)
{
	sqrt(r, bigfp(r), p);
}


double validate_sqrt(const bigfp_src& r, const bigfp_src& s, size_t p)
{
	bigfp r_sq;
	mul(r_sq, r, r);
	bigfp diff;
	sub(diff, r_sq, s);

	int64_t e = s.exp % 2 == 0 ? s.exp / 2 - 1 : floor_div(s.exp, 2);
	bigfp upper_bound;
	mul_pow2(upper_bound, r, -((int64_t)p - e - 1) * 32);

	diff.sign = abs(diff.sign);
	//return is_less_or_equal(diff, upper_bound);
	return div_apxm(diff, upper_bound);
}