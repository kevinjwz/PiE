#include <random>
#include <assert.h>
#include <string>
#include <iostream>
#include "bigint.h"
#include "bigfp.h"
#include "timer.h"
using namespace std;

static std::mt19937 rng;

int maini1()
{
	int test_i = 0;
	while (true) {
		printf("%d\t", test_i);
		size_t n1 = rng() % 4, n2 = rng() % 4;
		bigint a, b, a_sub_b;
		int a_sub_b_sign;
		rand_bigint(a, n1);
		rand_bigint(b, n2);
		// assume add is correct
		a_sub_b = a;
		a_sub_b_sign = sub_sign_mag(a_sub_b, b);
		//a_sub_b_sign = sub_sign_mag(a_sub_b,a, b);
		if (a_sub_b_sign == 0)
			assert(a_sub_b.size() == 0);
		if (a_sub_b_sign >= 0) {
			bigint a_sub_b_add_b;
			add(a_sub_b_add_b, a_sub_b, b);
			assert(a_sub_b_add_b == a);
		}
		else {
			// a_sub_b == b-a
			bigint b_sub_a_add_a;
			add(b_sub_a_add_a, a_sub_b, a);
			assert(b_sub_a_add_a == b);
		}
		printf("PASS!\n");
		test_i += 1;
	}
}

bool validate_divmod(const bigint& q, const bigint& r, const bigint& dvd, const bigint& dvs)
{
	if (compare(r, dvs) >= 0)
		return false;

	bigint dvs_mul_q_add_r;
	mul_auto(dvs_mul_q_add_r, dvs, q);
	add(dvs_mul_q_add_r, r);
	if (dvs_mul_q_add_r != dvd)
		return false;
	
	return true;
}

int maini2()
{
	int test_i = 0;
	while (true) {
		printf("%d\t", test_i);
		size_t dvs_n = rng() % 50 + 1, q_n = rng() % 50 + 1, diff_n = rng()%20;
		bigint dvd, dvs, q, r, diff;
		int a_sub_b_sign;
		rand_bigint(q, q_n);
		rand_bigint(dvs, dvs_n);
		rand_bigint(diff, diff_n);
		mul_auto(dvd, dvs, q);
		if (rng() % 2)
			add(dvd, diff);
		else
			sub_sign_mag(dvd, diff);

		div_rem(q, r, dvd, dvs);
		assert(validate_divmod(q, r, dvd, dvs));
		printf("PASS!\n");
		test_i += 1;
	}
}

