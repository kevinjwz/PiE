#pragma once
#include<vector>
#include<string>
#include<stdint.h>
#include"readonly_span.h"

typedef std::vector<uint32_t> bigint;
typedef readonly_span<uint32_t> bigint_src;

inline bool is_alias(const bigint& a, const bigint& b)
{
	return &a == &b;
}

#define noinline __declspec(noinline)

void rand_bigint(bigint& a, size_t len);

void add(bigint& c, const bigint& a, const bigint& b);
void add(bigint& a, const bigint& b);
void add_1(bigint& a, uint32_t b);

int sub_sign_mag(bigint& c, const bigint& a, const bigint& b);
int sub_sign_mag(bigint& a, const bigint& b);
int sub_1_sign_mag(bigint& a, uint32_t b);

int compare(const bigint& a, const bigint& b);

void mul_naive(bigint& c, const bigint_src& a, const bigint_src& b);

void div_rem_1_naive(bigint& quotient, uint32_t& remainder, const bigint& dividend, uint32_t divisor);

void mul_auto(bigint& c, const bigint_src& a, const bigint_src& b);

void div_rem(bigint& quotient, bigint& remainder, const bigint& dividend, const bigint& divisor);

struct bigfp;
void precompute_rcp(bigfp& rcp, const bigint& dvs, size_t dvd_n_max);
void div_rem_with_rcp(bigint& q, bigint& r, const bigint& dvd, const bigint& dvs, const bigfp& rcp);