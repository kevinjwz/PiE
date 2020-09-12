#include<assert.h>
#include<math.h>
#include<algorithm>
#include"bigint.h"
#include"bigfp.h"
#include"int128.h"
#include"micro_helpers.inl"

static bool is_ge_1_lt_B(const bigfp_src& a)
{
	return a.sign == 1 && a.exp == 1;
}

static bool is_eq_B(const bigfp_src& a)
{
	return a.sign == 1 && a.exp == 2 && a.data.size() == 1 && a.data[0] == 1;
}

static bool is_ge_1_le_B(const bigfp_src& a)
{
	return is_ge_1_lt_B(a) || is_eq_B(a);
}


void rcp_newton_iterate(bigfp& r, const bigfp_src& r0, const bigfp_src& d, size_t p)
{
	assert(is_ge_1_le_B(d));

	bigfp d_mul_r0;
	mul(d_mul_r0, d, r0);
	round(d_mul_r0, -(int128)p);

	bigfp one(1);
	bigfp e0;
	sub(e0, one, d_mul_r0);
	free(d_mul_r0);

	bigfp r0_mul_e0;
	mul(r0_mul_e0, r0, e0);
	free(e0);
	round(r0_mul_e0, -(int128)p);

	add(r, r0, r0_mul_e0);
}

/*
d in [1, B)
precision of result:
p fraction words + 1 preserved word
( |r-1/d|<=0.5*B^(-p) )
*/
void rcp_newton_base(bigfp& r, const bigfp_src& d, size_t p)
{
	assert(p == 1 || p == 2);
	assert(is_ge_1_lt_B(d));

	if (p == 1) {
		double dd = to_double(d);
		from_double(r, 1 / dd);
		// precision of dd: 53 bits
		// precision of r:  52 bits
		round(r, -2);
		//assert(validate_rcp(r, d, 1));

		return;
	}
	if (p == 2) {
		size_t lz = leading_zeros(*d.data.rbegin());	// [0, 31]
		int e = 31 - (int)lz;	// [0, 31]

		double dd_norm = ldexp(to_double(d), -e);
		// dd_norm in [1, 2]
		// |dd_norm - d/2^e| <= 0.5*2^(-52)
		bigfp r_norm;
		from_double(r_norm, 1 / dd_norm);
		// r_norm in [0.5, 1]
		// |r_norm - 1/dd_norm| <= 0.5*2^(-53)
		// |r_norm - 1/d*2^e| <= 0.5*2^(-51)
		// need precision: 0.5*2^(-64+e)

		if (51 < 64 - e) {
			bigfp d_round;
			round(d_round, d, -3);	// [1, B]
			bigfp d_norm;
			mul_pow2(d_norm, d_round, -e);
			round(d_norm, -3);
			// d_norm in [1, 2]
			// |d_norm - d/2^e| <= 0.5*2^(-95)
			rcp_newton_iterate(r_norm, bigfp(r_norm), d_norm, 3);
			// |r_norm - 1/d*2^e| <= 0.5*2^(-93)
		}

		mul_pow2(r, r_norm, -e);
		round(r, -3);
		//assert(validate_rcp(r, d, 2));
	}
}


/*
d in [1, B)
precision of result:
p fraction words + 1 preserved word
( |r-1/d|<=0.5*B^(-p) )
*/
void rcp_newton_rec(bigfp& r, const bigfp_src& d, size_t p)
{
	assert(p <= SIZE_MAX - 1);
	if (p <= 2) {
		rcp_newton_base(r, d, p);
		return;
	}

	bigfp r0;
	size_t p0 = p / 2 + 1;	// p0 = ceil((p+1)/2)
	rcp_newton_rec(r0, d, p0);

	bigfp d_round;
	round(d_round, d, -(int128)(p + 1));
	rcp_newton_iterate(r, r0, d_round, p + 1);
	//assert(validate_rcp(r, d, p));
}

/*
reciprocal with relative precision: p precise words + 1 preserved word
( |r-1/d| <= 0.5*B^(-p) * B^(1-d.exp) )
*/
void rcp_newton(bigfp& r, const bigfp_src& d, size_t p)
{
	assert(p > 0);
	assert(d.sign != 0);

	/*
		normalize d:
		d = d_norm * B^(d.exp-1)
		d_norm in [1, B)

		|r_norm - 1/d_norm| <= 0.5*B^(-p)
		r = r_norm * B^(1-d.exp)
		d = d_norm * B^(d.exp-1)
		|r - 1/d| <= 0.5*B^(-p+1-d.exp)
	*/


	bigfp_src d_norm(d);
	d_norm.sign = 1;
	d_norm.exp = 1;

	assert(is_ge_1_lt_B(d_norm));
	rcp_newton_rec(r, d_norm, p);
	// r in (1/B, 1]
	int128 r_exp = (int128)r.exp + 1 - d.exp;
	assert(r_exp >= INT64_MIN && r_exp <= INT64_MAX);
	r.exp = (int64_t)r_exp;
	r.sign = d.sign;
}


double div_apxm(const bigfp_src& a, const bigfp_src& b)
{
	bigfp q, rcp, a_round;
	rcp_newton(rcp, b, 3);
	round(a_round, a, (int128)a.exp - 3);
	mul(q, a_round, rcp);
	return to_double(q);
}

double validate_rcp(const bigfp_src& r, const bigfp_src& d, size_t p)
{
	bigfp r_mul_d;
	mul(r_mul_d, r, d);
	bigfp d_mul_h;
	mul_pow2(d_mul_h, d, -((int64_t)p + d.exp - 1) * 32 - 1);
	bigfp r_mul_d_sub_1;
	bigfp one(1);
	sub(r_mul_d_sub_1, r_mul_d, one);

	d_mul_h.sign = abs(d_mul_h.sign);
	r_mul_d_sub_1.sign = abs(r_mul_d_sub_1.sign);
	return div_apxm(r_mul_d_sub_1, d_mul_h);
}



/*
division with relative precision: p precise words + 1 preserved word
( |q-dvd/dvs| <= 0.5*B^(-p) * B^q.exp )
*/
void div(bigfp& q, const bigfp_src& dividend, const bigfp_src& divisor, size_t p)
{
	assert(p + 2 > p);
	bigfp rcp;
	rcp_newton(rcp, divisor, p + 2);
	bigfp dividend_round;
	round(dividend_round, dividend, (int128)dividend.exp - p - 2);
	mul(q, dividend_round, rcp);
	round(q, (int128)q.exp - p - 1);
}


double validate_div(const bigfp_src& q, const bigfp_src& dividend, const bigfp_src& divisor, size_t p)
{
	bigfp_src dvd_n(dividend), dvs_n(divisor);
	dvd_n.sign = 1;
	dvd_n.exp = 1;
	dvs_n.sign = 1;
	dvs_n.exp = 1;

	int64_t q_prec = is_less_or_equal(dvs_n, dvd_n) ? p - 1 : p;

	bigfp q_mul_dvs;
	mul(q_mul_dvs, q, divisor);
	bigfp diff;
	sub(diff, q_mul_dvs, dividend);

	bigfp upper_bound;
	mul_pow2(upper_bound, divisor, -(q_prec + divisor.exp - dividend.exp) * 32 - 1);

	diff.sign = abs(diff.sign);
	upper_bound.sign = abs(upper_bound.sign);
	return div_apxm(diff, upper_bound);
}


/*
reciprocal with absolute precision:
|r-1/d| <= 0.5*B^p
*/
void rcp_abs_prec(bigfp& r, const bigfp_src& d, int64_t p)
{
	assert(d.sign != 0);

	/*
	normalize d:
	d = d_norm * B^e
	d_norm in [1, B)
	r_norm = 1/d_norm
	r_norm in [1/B, 1]
	r = r_norm * B^(-e)

	prec(r_norm) * B^(-e) = 0.5*B^p
	prec(r_norm) = 0.5*B^(p+e)
	*/


	bigfp_src d_norm(d);
	d_norm.sign = 1;
	d_norm.exp = 1;

	int128 e = (int128)d.exp - 1;	// [INT64_MIN-1, INT64_MAX-1]
	int128 prec_norm = p + e;		// [2*INT64_MIN-1, 2*INT64_MAX-1]
	prec_norm = std::min<int128>(-1, prec_norm);	// [2*INT64_MIN-1, -1]
	rcp_newton_rec(r, d_norm, checked_cast<size_t>(-prec_norm));
	r.exp = checked_cast<int64_t>(r.exp - e);
	r.sign = d.sign;
}


double validate_rcp_abs_prec(const bigfp_src& r, const bigfp_src& d, int64_t p)
{
	bigfp r_mul_d;
	mul(r_mul_d, r, d);
	bigfp d_mul_h;
	mul_pow2(d_mul_h, d, p * 32 - 1);
	bigfp r_mul_d_sub_1;
	bigfp one(1);
	sub(r_mul_d_sub_1, r_mul_d, one);

	d_mul_h.sign = abs(d_mul_h.sign);
	r_mul_d_sub_1.sign = abs(r_mul_d_sub_1.sign);
	return div_apxm(r_mul_d_sub_1, d_mul_h);
}

/*
division with absolute precision:
|q-dvd/dvs| <= 0.5*B^p
*/
void div_abs_prec(bigfp& q, const bigfp_src& dividend, const bigfp_src& divisor, int64_t p)
{
	assert(divisor.sign != 0);
	if (dividend.sign == 0) {
		q.sign = 0;
		q.data = {};
		return;
	}
	/*
	normalize dvd, dvs:
	dvs = dvs_norm * B^e1
	dvs_norm in [1, B)
	dvd = dvd_norm * B^e2
	dvd_norm in [1/B, 1)

	r_norm = 1/dvs_norm
	q_norm = dvd_norm * r_norm
	q = q_norm * B^(e2-e1)

	prec(q_norm) * B^(e2-e1) = 0.5*B^p
	prec(q_norm) = 0.5*B^(p+e1-e2)
	prec(r_norm) = 0.5*B^(p+e1-e2-1)
	prec(dvd_norm) = 0.5*B^(p+e1-e2-1)
	*/

	bigfp_src dvs_norm(divisor);
	dvs_norm.sign = 1;
	dvs_norm.exp = 1;

	bigfp_src dvd_norm(dividend);
	dvd_norm.sign = 1;
	dvd_norm.exp = 0;

	int128 e1 = (int128)divisor.exp - 1;
	int128 e2 = (int128)dividend.exp;
	int128 prec_norm = p + e1 - e2 - 1;
	prec_norm = std::min<int128>(-1, prec_norm);

	bigfp r_norm;
	rcp_newton_rec(r_norm, dvs_norm, checked_cast<size_t>(-prec_norm));
	bigfp dvd_round;
	round(dvd_round, dvd_norm, prec_norm);
	mul(q, dvd_round, r_norm);

	q.exp = checked_cast<int64_t>(q.exp + e2 - e1);
	q.sign = dividend.sign * divisor.sign;
	round(q, p - 1);
}

double validate_div_abs_prec(const bigfp_src& q, const bigfp_src& dividend, const bigfp_src& divisor, int64_t p)
{
	bigfp q_mul_dvs;
	mul(q_mul_dvs, q, divisor);
	bigfp diff;
	sub(diff, q_mul_dvs, dividend);

	bigfp upper_bound;
	mul_pow2(upper_bound, divisor, p * 32 - 1);

	diff.sign = abs(diff.sign);
	upper_bound.sign = abs(upper_bound.sign);
	return div_apxm(diff, upper_bound);
}