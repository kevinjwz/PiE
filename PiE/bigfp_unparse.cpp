#include<assert.h>
#include<random>
#include<string>
#include<algorithm>
#include<vector>
#include"bigfp.h"
#include"timer.h"
#include"parse_unparse.h"
#include"micro_helpers.inl"

using namespace std;

void bigfp_pow10(bigfp& r, uint64_t exp)
{
	// using simplified 2^k method, k=3
	// pow10_small[i] = 10^i
	bigfp pow10_small[8] = { bigfp(1) };
	bigfp ten(10);
	for (size_t i = 1; i < 8; ++i) {
		mul(pow10_small[i], pow10_small[i - 1], ten);
	}
	if (exp < 8) {
		r = pow10_small[(size_t)exp];
		return;
	}

	unsigned bits = 64 - leading_zeros(exp);
	unsigned shift_amt = (bits - 1) / 3 * 3;	// largest multiple of 3 < bits

	bigfp prod(1), temp;
	for (; shift_amt > 0; shift_amt -= 3) {
		size_t i = (exp >> shift_amt) & 7;
		mul(temp, prod, pow10_small[i]);
		swap(prod, temp);
		unsigned square_times = 3;
		while (square_times--) {
			mul(temp, prod, prod);
			swap(prod, temp);
		}
	}
	free(temp);
	mul(r, prod, pow10_small[exp & 7]);
}

/*
frac in [0, 1)
no decimal point
precision = 0.5*10^(-n)
return: carry
(if frac is rounded to 1.0, str="0000...", return 1;
else return 0)
*/
int unparse_fraction(char str[], size_t n, const bigfp_src& frac)
{
	assert(frac.sign >= 0 && frac.exp <= 0);
	bigfp pow10;
	bigfp_pow10(pow10, checked_cast<uint64_t>(n));
	bigfp frac_mul_pow10;
	mul(frac_mul_pow10, frac, pow10);
	round(frac_mul_pow10, 0);
	
	// round to 10^n ??
	if (is_equal(frac_mul_pow10, pow10)) {
		fill_n(str, n, '0');
		return 1;
	}
	unparse_fast(str, n, frac_mul_pow10);
	return 0;
}

void unparse_integer_fraction(vector<char>& int_str, vector<char>& frac_str, const bigfp_src& src, size_t frac_digit_n)
{
	assert(src.sign >= 0);
	int_str.clear();
	frac_str.clear();
	frac_str.resize(frac_digit_n);
	bigfp_src frac = fractional_part(src);
	int carry = unparse_fraction(frac_str.data(), frac_digit_n, frac);
	bigfp int_part_add_1;
	if (carry != 0) {
		add(int_part_add_1, integral_part(src), bigfp_1);
	}
	bigfp_src int_part = carry != 0 ? bigfp_src(int_part_add_1) : integral_part(src);
	if (int_part.sign == 0) {
		int_str.push_back('0');
		return;
	}

	size_t int_digit_n_est = estimate_unparse_length(int_part);
	int_str.resize(int_digit_n_est);
	unparse_fast(int_str.data(), int_digit_n_est, int_part);
	auto nz_it = find_if(int_str.begin(), int_str.end(), [](char c) { return c != '0'; });
	int_str.erase(int_str.begin(), nz_it);
}


int mainu1()
{
	while (true) {
		uint64_t exp;
		scanf("%llu", &exp);
		bigfp pow10;
		bigfp_pow10(pow10, exp);
		bigint bi;
		to_bigint(bi, pow10);
		string str;
		unparse_naive(str, bi);
		printf("%s\n", str.data());
	}
}

int mainu2()
{
	bigfp pow10;
	for (size_t lenb = 16; lenb <= 1 << 22; lenb *= 2) {
		printf("================\n");
		printf("n = %zu\n", lenb);
		timer("pow10:    %f\n", [&] {bigfp_pow10(pow10, lenb); });
		printf("n = %zu\n", lenb * 2 - 1);
		timer("pow10:    %f\n", [&] {bigfp_pow10(pow10, lenb * 2 - 1); });
	}
	return 0;
}

int mainu3()
{
	for (size_t lenb = 16; lenb <= 1 << 22; lenb *= 2) {
		size_t len = lenb;
		//str_n = (1<<20)*3;
		printf("================\n");
		printf("n = %zu\n", len);
		bigint dvd, dvs;
		rand_bigint(dvd, len);
		rand_bigint(dvs, len / 2);
		bigint q(len / 2 + 1), r(len / 2 + 1);
		bigfp qf, rf;
		qf.data.reserve(len / 2 + 1);
		rf.data.reserve(len / 2 + 1);
		bigfp rcp;
		precompute_rcp(rcp, dvs, len);
		timer("bigint:    %f sec\n", [&] {div_rem_with_rcp(q, r, dvd, dvs, rcp); });
		timer("bigfp:     %f sec\n", [&] {div_rem_with_rcp(qf, rf, dvd, dvs, rcp); });
	}
	return 0;
}