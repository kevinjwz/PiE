#include <assert.h>
#include <math.h>
#include <algorithm>
#include <random>
#include <iostream>
#include <limits>
#include "bigint.h"
#include "bigfp.h"
#include "micro_helpers.inl"
#include "bigfp_tmpl_tools.inl"
#include "timer.h"

using namespace std;

void parse_naive(bigint& a, const std::string& str)
{
	bigint prod, ten = { 10 }, add_tail = { 0 };
	a.clear();
	size_t n = (size_t)ceil(str.size()*(log2(10) / 32));
	a.reserve(n);
	prod.reserve(n);
	for (size_t si = 0; si < str.size(); ++si) {
		mul_naive(prod, a, ten);
		add_tail[0] = str[si] - '0';
		add(a, prod, add_tail);
	}
}

void unparse_naive(std::string& str, const bigint& a)
{
	if (a.size() == 0) {
		str = "0";
		return;
	}

	bigint copy = a;
	str.clear();
	while (copy.size() != 0) {
		uint32_t r;
		div_rem_1_naive(copy, r, copy, 10);
		str.append(1, (char)(r + '0'));
	}
	std::reverse(str.begin(), str.end());
}

/*
n > 0
*/
static void parse_base(bigfp& a, const char str[], size_t n)
{
	uint64_t au = str[0] - '0';
	for (size_t i = 1; i < n; ++i) {
		au = au * 10 + str[i] - '0';
	}
	from_uint(a, au);
}

// 10^n-1 <= UINT64_MAX
static const size_t PARSE_BASE_N = numeric_limits<uint64_t>::digits10;

/*
n>0
size(pow10)>=ceil_log2(n)
*/
static void parse_rec(bigfp& a, const char str[], size_t n, const std::vector<bigfp>& pow10)
{
	if (n <= PARSE_BASE_N) {
		parse_base(a, str, n);
		return;
	}
	// n >= 2
	// pow10_i >= 0
	// lo_str_n >= 1
	// lo_str_n = clp2(n)/2 < n
	size_t pow10_i = ceil_log2(n) - 1;
	size_t lo_str_n = 1ull << pow10_i;
	size_t hi_str_n = n - lo_str_n;
	const char* lo_str = str + hi_str_n;
	const char* hi_str = str;

	bigfp lo, hi;
	parse_rec(hi, hi_str, hi_str_n, pow10);
	bigfp hi_mul_pow10;
	mul(hi_mul_pow10, hi, pow10[pow10_i]);

	parse_rec(lo, lo_str, lo_str_n, pow10);
	add(a, hi_mul_pow10, lo);
}

/*
str must be trimmed (no leading/trailing spaces)
str is natural writting order
0        n-1
1234567890
^        ^
|        least significant
most significant
*/
void parse_fast(bigfp& a, const char str[], size_t n)
{
	if (n == 0) {
		a = bigfp_0;
		return;
	}
	if (n <= PARSE_BASE_N) {
		parse_base(a, str, n);
		return;
	}
	// n >= 2
	// i_max >= 0
	size_t i_max = ceil_log2(n) - 1;
	// 2^i_max in [n/2, n-1]
	// pow10[i] = 10^(2^i): 10, 100, 10000, 100000000, ...
	std::vector<bigfp> pow10(i_max + 1);
	from_uint(pow10[0], 10u);
	for (size_t i = 1; i <= i_max; ++i)
		mul(pow10[i], pow10[i - 1], pow10[i - 1]);
	parse_rec(a, str, n, pow10);
}



static void unparse_base(char str[], size_t n, const bigfp_src& a)
{
	uint64_t au = to_uint_checked<uint64_t>(a);
	for (size_t i = n; i != 0;) {
		i -= 1;
		str[i] = au % 10 + '0';
		au /= 10;
	}
	assert(au == 0);
}

static const size_t UNPARSE_BASE_N = numeric_limits<uint64_t>::digits10;

void unparse_rec(char str[], size_t n, const bigfp_src& a, const std::vector<bigfp>& pow10, const std::vector<bigfp>& pow10_rcp)
{
	if (n <= UNPARSE_BASE_N) {
		unparse_base(str, n, a);
		return;
	}

	// n >= 2
	// lo_str_n in [n/2, n-1]
	// hi_str_n in [1, n/2]
	size_t pow10_i = ceil_log2(n) - 1;
	size_t lo_str_n = 1ull << pow10_i;
	size_t hi_str_n = n - lo_str_n;
	char* lo_str = str + hi_str_n;
	char* hi_str = str;

	bigfp q, r;
	div_rem_with_rcp(q, r, a, pow10[pow10_i], pow10_rcp[pow10_i]);
	unparse_rec(hi_str, hi_str_n, q, pow10, pow10_rcp);
	unparse_rec(lo_str, lo_str_n, r, pow10, pow10_rcp);
}

/*
a must be non-negative integer
a < 10^n  (representable in n digits)
*/
void unparse_fast(char str[], size_t n, const bigfp_src& a)
{
	assert(is_non_neg_integer(a));
	if (n <= UNPARSE_BASE_N) {
		unparse_base(str, n, a);
		return;
	}
	size_t i_max = ceil_log2(n) - 1;
	// 2^i_max in [n/2, n-1]
	// pow10[i] = 10^(2^i): 10, 100, 10000, 100000000, ...
	std::vector<bigfp> pow10(i_max + 1);
	std::vector<bigfp> pow10_rcp(i_max + 1);
	from_uint(pow10[0], 10u);
	for (size_t i = 1; i <= i_max; ++i)
		mul(pow10[i], pow10[i - 1], pow10[i - 1]);
	for (size_t i = 0; i <= i_max; ++i)
		precompute_rcp(pow10_rcp[i], pow10[i], pow10[i].exp * 2);
	unparse_rec(str, n, a, pow10, pow10_rcp);
}

size_t estimate_unparse_length(size_t word_n)
{
	size_t an = word_n;
	if (an == 0)
		return 1;
	/*
		10^n-1 >= a
		n >= log10(a+1)
		n >= log10(B^an) = an*log10(B)
	*/
	double log10_B = log10(1ull << 32);
	/*
		log10(B) ~ 9.6
		(1)001.xxxx....xxxx
		       | 49 bits  |
		|log10_B - log10(B)| <= 0.5*2^(-49) = 2^(-50)
		upper bound of log10(B) = log10_B + 2^(-50)
	*/
	bigfp log10_B_bf(log10_B);
	bigfp log10_B_err(1.0 / (1ull << 50));
	bigfp log10_B_up;
	add(log10_B_up, log10_B_bf, log10_B_err);
	
	bigfp an_bf;
	from_int(an_bf, an);

	bigfp prod;
	mul(prod, an_bf, log10_B_up);

	size_t n = to_uint_checked<size_t>(prod);
	assert(prod.exp >= 0);
	if (prod.data.size() > (uint64_t)prod.exp) {
		n += 1;
		assert(n != 0);	// no overflow
	}
	return n;
}

size_t estimate_unparse_length(const bigfp_src& src)
{
	assert(is_non_neg_integer(src));
	if (src.sign == 0)
		return 1;
	size_t int_word_n = checked_cast<size_t>(src.exp);
	return estimate_unparse_length(int_word_n);
}

static std::mt19937 rng;

int mainp3()
{
	int test_i = 0;
	while (true) {
		printf("%d\t", test_i);
		size_t n = rng() % 10000;
		bigint a;
		rand_bigint(a, n);
		size_t str_n = estimate_unparse_length(n);
		string str1(str_n, 0), str2(str_n, 0);
		unparse_fast(str1.data(), str_n, a);
		unparse_fast(str2.data(), str_n, a);
		assert(str1 == str2);
		printf("PASS!\n");
		test_i += 1;
	}
}

int mainp4()
{
	for (size_t lenb = 16; lenb <= 1 << 22; lenb *= 2) {
		size_t str_n = lenb;
		//str_n = (1<<20)*3;
		printf("================\n");
		printf("n = %zu\n", str_n);
		string str(str_n, 0);
		for (char& c : str)
			c = rng() % 10 + '0';
		bigfp f;
		f.data.resize(str_n / 9);
		timer("parse:    %f\n", [&] {parse_fast(f, str.data(), str_n); });
		timer("unparse:  %f\n", [&] {unparse_fast(str.data(), str_n, f); });

	}
	return 0;
}

int mainp5()
{
	while (true) {
		string str;
		cin >> str;
		bigfp a;
		parse_fast(a, str.data(), str.size());
		bigint ai;
		to_bigint(ai, a);
		size_t parse_n_est = estimate_unparse_length(a.exp);
		string str_fast(parse_n_est, 0), gt;
		unparse_fast(str_fast.data(), parse_n_est, a);
		unparse_naive(gt, ai);

		cout << str_fast << endl;
		cout << gt << endl;

	}
}
