#include<utility>
#include<assert.h>
#include<random>
#include<string>
#include<vector>
#include<sstream>
#include<iostream>
#include<fstream>
#include"bigfp.h"
#include"timer.h"
#include"micro_helpers.inl"
#include"bigfp_tmpl_tools.inl"
#include"parse_unparse.h"

using namespace std;

void pi_agm(bigfp& pi, size_t p)
{
	bigfp ai(1);
	bigfp bi;
	sqrt(bi, bigfp(0.5), p + 1);
	bigfp ti(0.25);
	bigfp ci1, ci1_sq;
	bigfp ai1, bi1, ti1;
	int i = 0;
	bigfp ci1_upperbound((double)UINT32_MAX - 48);
	mul_pow2(ci1_upperbound, -((int64_t)p + 1) * 32 - 4);
	bigfp ci1_upperbound_dec(40);
	mul_pow2(ci1_upperbound_dec, -((int64_t)p + 1) * 32 - 4);
	while (true) {
		//printf("i = %d\n", i);
		add(ai1, ai, bi);
		mul_pow2(ai1, -1);
		round(ai1, -(int128)p - 1);

		sub(ci1, ai, ai1);
		if (is_less_or_equal(ci1, ci1_upperbound)) {
			break;
		}
		sub(ci1_upperbound, bigfp(ci1_upperbound), ci1_upperbound_dec);
		assert(ci1_upperbound.sign > 0);

		mul(ci1_sq, ci1, ci1);
		mul_pow2(ci1_sq, i);
		round(ci1_sq, -(int128)p - 1);
		sub(ti1, ti, ci1_sq);

		mul(bi1, ai, bi);
		sqrt(bi1, p + 1);

		swap(ai, ai1);
		swap(bi, bi1);
		swap(ti, ti1);
		i += 1;
	}
	bigfp ai_add_bi;
	add(ai_add_bi, ai, bi);
	bigfp numerator;
	mul(numerator, ai_add_bi, ai_add_bi);
	mul_pow2(ti, 2);
	div(pi, numerator, ti, p + 2);
}

/*
return: p
0.5*B^(-p) <= 0.45*10^(-dn)
*/
size_t digits_to_prec(size_t dn)
{
	/*	p >= 1/32 * (dn*log2(10) - log2(0.9))
		p >= 1/32 * ((dn+1)*log2(10) - log2(9))
	*/
	double log2_9 = log2(9);
	/*	log2(9) ~ 3.17
		(1)1.xxxx....xxxx
			 | 51 bits  |
		lower bound of log2(9) = log2_9 - 2^(-52)
	*/
	bigfp log2_9_lo;
	sub(log2_9_lo, bigfp(log2_9), bigfp(1.0 / (1ull << 52)));

	double log2_10 = log2(10);
	/*	log2(10) ~ 3.32
		(1)1.xxxx....xxxx
			 | 51 bits  |
		upper bound of log2(10) = log2_10 + 2^(-52)
	*/
	bigfp log2_10_up;
	add(log2_10_up, bigfp(log2_10), bigfp(1.0 / (1ull << 52)));

	bigfp dn_bf, dn_add_1;
	from_uint(dn_bf, dn);
	add(dn_add_1, dn_bf, bigfp_1);

	bigfp prod;
	mul(prod, dn_add_1, log2_10_up);
	bigfp diff;
	sub(diff, prod, log2_9_lo);
	mul_pow2(diff, -5);	// multiply 1/32

	size_t p = to_uint_checked<size_t>(diff);
	if (safe_compare(diff.data.size(), diff.exp) > 0) {
		p += 1;
		assert(p != 0);
	}
	return p;
}
static std::mt19937 rng;

int mainpi()
{
	bigfp pi_test, pi_gt, diff;
	bigfp upperbound(0.5);
	pi_agm(pi_gt, 60000);
	while (true) {
		size_t p = rng() % 50000 + 100;
		printf("p = %zu\t", p);
		pi_agm(pi_test, p);
		sub(diff, pi_test, pi_gt);
		diff.sign = abs(diff.sign);
		upperbound.exp = -(int64_t)p;
		assert(is_less_or_equal(diff, upperbound));
		printf("PASS!\n");
	}


}

int main(int argc, char** argv)
{
	assert(argc == 2);
	istringstream dn_str(argv[1]);
	size_t dn;
	dn_str >> dn;
	if (dn_str.fail() || !dn_str.eof()) {
		cerr << "argument error!\n";
	}

	bigfp pi;

	timer("compute time: %f sec\n", [&] {pi_agm(pi, digits_to_prec(dn)); });
	
	vector<char> i_str, f_str;
	timer("unparse time: %f sec\n", [&] {unparse_integer_fraction(i_str, f_str, pi, dn + 1); });

	ostringstream filename("pi_", ios::app);
	filename << dn << ".txt";
	
	ofstream file(filename.str(), ios::trunc);
	timer("output time:  %f sec\n", [&] {print_number_formatted(file, i_str, f_str, 1); });
}