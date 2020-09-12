#include<stdio.h>
#include<random>
#include<omp.h>
#include<iostream>
#include<time.h>
#include <assert.h>
#include <functional>
#include "bigint.h"
#include "timer.h"
#include "mul_karatsuba.h"
#include "mul_ssa.h"
#include "mul_fft_dwt.h"
#include "parse_unparse.h"

int main1()
{
	while (true) {
		std::cout << "================" << std::endl;
		std::string astr, bstr, cstr;
		std::cin >> astr >> bstr;
		bigint a, b, c;
		parse_naive(a, astr);
		parse_naive(b, bstr);
		int sign;
		//mul_naive(c, a, b);
		mul_dwt(c, a, b);
		unparse_naive(cstr, c);
		std::cout << cstr << std::endl;

	}
}

static std::mt19937 rng;

void rand_bigint(bigint& a, size_t len)
{
	if (len == 0) {
		a.resize(0);
		return;
	}

	a.resize(len);

	for (size_t i = 0; i != len - 1; ++i) {
		uint32_t r = rng();
		a[i] = r;
	}
	while (true) {
		uint32_t r = rng();
		if (r != 0) {
			a[len - 1] = r;
			return;
		}
	}
}

void max_bigint(bigint& a, size_t len)
{
	a.resize(len);

	for (size_t i = 0; i != len; ++i) {
		a[i] = 0xffffffff;
	}
}

void rand_large_bigint(bigint& a, size_t len)
{
	a.resize(len);

	for (size_t i = 0; i != len; ++i) {
		a[i] = 0xffff0000 | rng();
	}
}

void bigint_ones(bigint& a, size_t len)
{
	a.resize(len);

	for (size_t i = 0; i != len; ++i) {
		a[i] = 1;
	}
}

template <class f_mul, class f_init>
bool test_ub(f_mul mul, f_init init, size_t len1, size_t len2)
{
	bigint a, b, prod, prod_gt;
	init(a, len1);
	init(b, len2);
	mul(a, b, prod);
	mul_ssa_pure_int(prod_gt, a, b);
	bool r = prod == prod_gt;
	return r;
}

template <class f_mul, class f_init>
void run_ub(f_mul mul, f_init init, size_t len1, size_t len2)
{
	bigint a, b, prod;
	init(a, len1);
	init(b, len2);
	mul(a, b, prod);
}

template <class func>
void test_forever(func mul, size_t len_min, size_t len_max)
{
	while (true) {
		size_t len = len_min + rng() % (len_max - len_min + 1);
		printf("================\n");
		printf("len = %d\n", len);

		assert(test_ub(mul, rand_bigint, len, len));
		printf("rand:  PASS!\n");

		assert(test_ub(mul, rand_large_bigint, len, len));
		printf("rlg:   PASS!\n");

		//assert(test_ub(mul, max_bigint, len, len));
		//printf("max:   PASS!\n");
	}
}

template <class func>
void test_forever_ub(func mul, size_t len_min, size_t len_max)
{
	while (true) {
		size_t len1 = len_min + rng() % (len_max - len_min + 1);
		size_t len2 = len_min + rng() % (len_max - len_min + 1);
		printf("================\n");
		printf("len1 = %zd  len2 = %zd\n", len1, len2);

		assert(test_ub(mul, rand_bigint, len1, len2));
		printf("rand:  PASS!\n");

		assert(test_ub(mul, rand_large_bigint, len1, len2));
		printf("rlg:   PASS!\n");

		assert(test_ub(mul, max_bigint, len1, len2));
		printf("max:   PASS!\n");
	}
}

void warm_up(uint32_t n)
{
	while (n--) {
		run_ub(mul_naive, rand_bigint, 100, 100);
		run_ub(mul_karatsuba_unbalance, rand_bigint, 1000, 100);
		run_ub(mul_dwt, rand_bigint, 1024, 1024);
		run_ub(mul_ssa_fast, rand_bigint, 65536, 65536);
	}
}

void fastest_mul(size_t an_min, size_t an_max, size_t bn_min, size_t bn_max)
{
	rng.seed((uint32_t)time(nullptr));
	// warm up

	bigint a, b, r;
	size_t n = 0;
	while (true) {
		size_t an = rng() % (an_max - an_min + 1) + an_min;
		size_t bn = rng() % (bn_max - bn_min + 1) + bn_min;

		rand_bigint(a, an);
		rand_bigint(b, bn);
		r.resize(an + bn);

		double time_list[4];
		if (an*bn >= 1024 * 1024 && std::min(an,bn)>100)
			time_list[0] = 86400;
		else
			time_list[0] = timer("", [&] {mul_naive(r, a, b); });

		if (an*bn >= 16384 * 16384 && std::min(an,bn)>300)
			time_list[1] = 86400;
		else
			time_list[1] = timer("", [&] {mul_karatsuba_unbalance(r, a, b); });

		time_list[2] = timer("", [&] {mul_dwt(r, a, b); });
		time_list[3] = timer("", [&] {mul_ssa_fast(r, a, b); });

		const char* name_list[4] = { "naive" , "krtb" ,"dwt", "ssa" };
		size_t min_index = std::min_element(time_list, time_list + 4) - time_list;
		printf("%zu %zu %s\n", an, bn, name_list[min_index]);
		n += 1;
		fprintf(stderr, "\r                \r");
		fprintf(stderr, "%zu", n);
	}
}

int main2(int argc, char** argv)
{
	//warm_up(100);

	//size_t an_min, an_max, bn_min, bn_max;
	//sscanf(argv[1], "%zu", &an_min);
	//sscanf(argv[2], "%zu", &an_max);
	//sscanf(argv[3], "%zu", &bn_min);
	//sscanf(argv[4], "%zu", &bn_max);
	//fastest_mul(an_min, an_max, bn_min, bn_max);
	bigint a, b, c(1<<16);
	//rand_bigint(a, 1 << 15);
	//rand_bigint(b, 1 << 15);
	//int n = 5;
	//while (n--) {
	//	mul_naive(c, a, b);
	//	printf("========\n");
	//}
	//return 0;

	for (size_t lenb = 16; lenb <= 1<<20; lenb *= 2) {
		//size_t len = lenb + 1;
		size_t len1 = lenb, len2 = lenb;
		printf("================\n");
		//printf("len = %zd\n", len);
		printf("len1 = %zd  len2 = %zd\n", len1, len2);
		c.resize(len1 + len2);
		rand_bigint(a, len1);
		rand_bigint(b, len2);
		//bigint_ones(a, len1);
		//bigint_ones(b, len2);
		//timer("naive:               %f sec\n", [&] {mul_naive(c, a, b); });
		//timer("krtb_proto:          %f sec\n", [&] {mul_karatsuba_prototype(c, a, b); });
		//timer("krtb_ub:             %f sec\n", [&] {mul_karatsuba_unbalance(c, a, b); });

		timer("ssa_i:               %f sec\n", [&] {mul_ssa_pure_int(c, a, b); });
		timer("ssa_f:               %f sec\n", [&] {mul_ssa_fast(c, a, b); });
		//timer("fft:                 %f sec\n", [&] {mul_fft_fp(c, a, b); });
		//timer("dwt:                 %f sec\n", [&] {mul_dwt(c, a, b); });
	}
	return 0;
}


int main3()
{
	using namespace std;
	//size_t len = 1<<20;
	//while (1) {
	//	bigint a(len), b(len), c(2*len);
	//	a[len - 1] = 1;
	//	b[len/2] = 1;
	//	mul_ssa(a, b, c);
	//	printf("========\n");
	//}
	while (1) {
		size_t e = rand() % 5 + 12;
		size_t prod_len = 1ULL << e;
		size_t len1 = rng() % (prod_len - 1) + 1; // 1~prod_len-1
		size_t len2 = prod_len - len1;
		//size_t e1=10,e2=10;
		printf("e = %2zu, len1 = %10zu, len2 = %10zu\t", e, len1, len2);
		assert(test_ub(mul_dwt, rand_bigint, len1, len2));
		printf("rand:PASS! ");
		assert(test_ub(mul_dwt, rand_large_bigint, len1, len2));
		printf("rlg:PASS!\n");
	}
	return 0;
}


