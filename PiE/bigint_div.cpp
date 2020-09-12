#include <stdint.h>
#include <assert.h>
#include <algorithm>
#include "bigint.h"
#include "bigfp.h"
#include "int128.h"

/*
a>0, b>0
|a-b|<=1
*/
static bool floor_equal(const bigfp_src& a, const bigfp_src& b)
{
	if (a.exp != b.exp)
		return false;
	if (a.exp <= 0)	// a<1, b<1
		return true;

	/*
	a = aaaaaa.aaaaaaaaaa
	b = bbbbbb.bbbbbbb
			 ^
			 equal?
	*/
	uint64_t int_n = (uint64_t)a.exp;	// integer part length
	size_t an = a.data.size(), bn = b.data.size();
	uint32_t ai_lo = int_n > an ? 0 : a.data[an - (size_t)int_n];
	uint32_t bi_lo = int_n > bn ? 0 : b.data[bn - (size_t)int_n];
	return ai_lo == bi_lo;
}

/*
	dvd/dvs > B^(-1)
	|q_est - dvd/dvs| <= 0.5*B^(-1)
*/
static void exact_quot_rem(bigint& q, bigint& r, const bigfp& q_est, const bigint& dvd, const bigint& dvs)
{
	bigfp radius(0.5 / (1ULL << 32));	// 0.5*B^(-1)
	bigfp upb, lob;
	add(upb, q_est, radius);
	sub(lob, q_est, radius);
	// lob = q_est - 0.5/B >= dvd/dvs - 1/B > 0
	assert(lob.sign > 0);
	if (floor_equal(lob, upb)) {
		to_bigint(q, q_est);
		bigint q_mul_dvs;
		mul_auto(q_mul_dvs, q, dvs);
		int sign = sub_sign_mag(r, dvd, q_mul_dvs);
		assert(sign >= 0);
	}
	else {
		to_bigint(q, upb);
		bigint q_mul_dvs;
		mul_auto(q_mul_dvs, q, dvs);
		int sign = sub_sign_mag(r, dvd, q_mul_dvs);
		if (sign < 0) {
			sign = sub_1_sign_mag(q, 1);
			assert(sign >= 0);
			// r == -(r_gt - dvs) = dvs - r_gt
			// r_gt = dvs-r = -(r-dvs)
			sign = -sub_sign_mag(r, dvs);
			assert(sign >= 0);
		}
	}
}


void div_rem(bigint& q, bigint& r, const bigint& dvd, const bigint& dvs)
{
	size_t dvd_n = dvd.size(), dvs_n = dvs.size();
	assert(dvs_n > 0);
	if (dvd_n < dvs_n) {
		q = {};
		r = dvd;
		return;
	}
	// dvd_n >= dvs_n
	// dvd/dvs >= B^(-1)
	bigfp q_est;
	div_abs_prec(q_est, dvd, dvs, -1);
	// |q_est - dvd/dvs| <= 0.5*B^(-1)
	exact_quot_rem(q, r, q_est, dvd, dvs);
}



void precompute_rcp(bigfp& rcp, const bigint& dvs, size_t dvd_n_max)
{
	assert(is_non_neg_integer(dvs));

	int128 e2 = dvd_n_max;
	/*
		prec(r_norm) = 0.5*B^(p+e1-e2-1)
		prec(r) = prec(r_norm)/B^e1 = 0.5*B^(p-e2-1)
	*/
	rcp_abs_prec(rcp, dvs, checked_cast<int64_t>(-e2 - 2));
}



void div_rem_with_rcp(bigint& q, bigint& r, const bigint& dvd, const bigint& dvs, const bigfp& rcp)
{
	/*
		prec(dvd_norm) = 0.5*B^(p+e1-e2-1)
		prec(dvd) = prec(dvd_norm)*B^e2 = 0.5*B^(p+e1-1)
	*/
	int128 e1 = (int128)dvs.size() - 1;

	bigfp dvd_round(dvd);
	round(dvd_round, e1 - 2);
	bigfp q_est;
	mul(q_est, dvd_round, rcp);
	round(q_est, -2);
	exact_quot_rem(q, r, q_est, dvd, dvs);
}
