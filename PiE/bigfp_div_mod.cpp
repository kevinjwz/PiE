#include <assert.h>
#include "bigfp.h"
#include "int128.h"
#include "timer.h"

static bigfp radius(0.5 / (1ULL << 32));	// 0.5*B^(-1)

/*
	dvd.sign == dvs.sign, dvd!=0
	|q_est - dvd/dvs| <= 0.5*B^(-1)
*/
static void exact_quot_rem_same_sign(bigfp& q, bigfp& r, const bigfp_src& q_est, const bigfp_src& dvd, const bigfp_src& dvs)
{
	// dvd/dvs > 0
	add(q, q_est, radius);	// upperbound of dvd/dvs
	trunc(q, 0);			// floor(q_est+r) in [floor(dvd/dvs), floor(dvd/dvs+2r)] <= [floor(dvd/dvs), floor(dvd/dvs)+1]
	bigfp q_mul_dvs;
	mul(q_mul_dvs, q, dvs);
	sub(r, dvd, q_mul_dvs);
	if (r.sign + dvd.sign == 0) {	// r.sign = -dvd.sign
		sub(q, bigfp(q), bigfp_1);
		add(r, bigfp(r), dvs);
		assert(r.sign + dvd.sign != 0);	// r==0, or r.sign == dvd.sign
	}
}


/*
	dvd != 0
	|q_est - dvd/dvs| <= 0.5*B^(-1)
*/
static void exact_quot_rem(bigfp& q, bigfp& r, const bigfp& q_est, const bigfp_src& dvd, const bigfp_src& dvs)
{
	if (dvd.sign == dvs.sign) {
		exact_quot_rem_same_sign(q, r, q_est, dvd, dvs);
	}
	else {
		/*	dvd % dvs = dvd % (-dvs)
			dvd \ dvs = (dvd - dvd % dvs) / dvs
					  = (dvd - dvd % (-dvs)) / dvs
					  = -(dvd - dvd % (-dvs)) / (-dvs)
					  = -(dvd \ (-dvs))
		*/
		bigfp_src neg_dvs(dvs);
		neg_dvs.sign *= -1;
		bigfp_src neg_q_est(q_est);
		neg_q_est.sign *= -1;
		exact_quot_rem_same_sign(q, r, neg_q_est, dvd, neg_dvs);
		q.sign *= -1;
	}
}

void div_rem(bigfp& q, bigfp& r, const bigfp_src& dvd, const bigfp_src& dvs)
{
	assert(dvs.sign != 0);
	if (dvd.sign == 0) {
		q = bigfp_0;
		r = bigfp_0;
		return;
	}
	bigfp q_est;
	div_abs_prec(q_est, dvd, dvs, -1);	// |q_est - dvd/dvs| <= 0.5*B^(-1)
	exact_quot_rem(q, r, q_est, dvd, dvs);
}

void precompute_rcp(bigfp& rcp, const bigfp_src& dvs, int64_t dvd_exp_max)
{
	assert(dvs.sign != 0);

	int128 e2 = dvd_exp_max;
	/*
		prec(r_norm) = 0.5*B^(p+e1-e2-1)
		prec(r) = prec(r_norm)/B^e1 = 0.5*B^(p-e2-1)
	*/
	rcp_abs_prec(rcp, dvs, checked_cast<int64_t>(-e2 - 2));
}

void div_rem_with_rcp(bigfp& q, bigfp& r, const bigfp_src& dvd, const bigfp_src& dvs, const bigfp& rcp)
{
	/*
		prec(dvd_norm) = 0.5*B^(p+e1-e2-1)
		prec(dvd) = prec(dvd_norm)*B^e2 = 0.5*B^(p+e1-1)
	*/
	int128 e1 = (int128)dvs.exp - 1;

	bigfp dvd_round;
	round(dvd_round, dvd, e1 - 2);
	bigfp q_est;
	mul(q_est, dvd_round, rcp);
	round(q_est, -2);
	exact_quot_rem(q, r, q_est, dvd, dvs);
}

