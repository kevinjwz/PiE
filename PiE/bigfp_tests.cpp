#include<assert.h>
#include<math.h>
#include<map>
#include<functional>
#include<random>
#include"bigint.h"
#include"bigfp.h"
#include"int128.h"
#include"timer.h"

using namespace std;




map<char, function<double(double, double)>> double_arith = {
	{ '+', plus<double>()},
	{ '-', minus<double>() },
	{ '*', multiplies<double>() },
	{ '/', divides<double>() }
};

map<char, function<void(bigfp&, const bigfp_src&, const bigfp_src&)>> bigfp_arith = {
	{ '+', (void(*)(bigfp&, const bigfp_src&, const bigfp_src&))add },
	{ '-', sub },
	{ '*', mul },

};

static mt19937_64 rng64;
static mt19937 rng32;

uint32_t rand32_nonzero()
{
	uint32_t r = 0;
	while (r == 0)
		r = rng32();
	return r;
}

int64_t rand_between(int64_t min, int64_t max)
{
	assert(max >= min);
	if (min == INT64_MIN && max == INT64_MAX) {
		return (int64_t)rng64();
	}
	uint64_t range = (uint64_t)max - (uint64_t)min + 1;	// [1, UINT64_MAX]
	return (int64_t)(rng64() % range + (uint64_t)min);
}

double rand_double(int e_min, int e_max)
{
	uint64_t r = rng64();
	uint64_t frac_bits = r % (1ull << 52) + (1ull << 52);
	double frac = ldexp(frac_bits, -53);
	int e = (int)(r >> 53) % (e_max - e_min + 1) + e_min;
	int sign = (r&(1ull << 52)) != 0 ? -1 : 1;
	return sign * ldexp(frac, e);
}

void rand_bigfp(bigfp& f, size_t p, int64_t e_min, int64_t e_max)
{
	assert(p >= 1);

	f.sign = rng32() & 1 ? -1 : 1;
	f.exp = rand_between(e_min, e_max);

	f.data.resize(p);
	for (size_t i = 1; i < p - 1; ++i)
		f.data[i] = rng32();

	size_t loz = rng32() % 32;
	size_t hiz = rng32() % 32;

	f.data[0] = (rng32() | 1) << loz;
	f.data[p - 1] = (rng32() >> hiz) | (1 << (31 - hiz));
}

int mainf()
{
	while (true) {
		printf("================\n");
		double d1, d2, gt;
		char op[2];
		//int n = scanf("%lf %1[-+*] %lf", &d1, op, &d2);
		//assert(n == 3);		
		d1 = rand_double(-64, 64);
		d2 = rand_double(-64, 64);
		//d2 = (float)d1;
		//op[0] = (rand() & 1) ? '+' : '-';
		op[0] = '*';
		printf("%a %c %a\n", d1, op[0], d2);

		gt = double_arith[op[0]](d1, d2);
		bigfp bf1, bf2, bf3;
		from_double(bf1, d1);
		from_double(bf2, d2);
		bigfp_arith[op[0]](bf3, bf1, bf2);
		double from_bf = to_double(bf3);
		printf("gt = %a\n", gt);
		printf("bf = %a\n", from_bf);
		if (gt != from_bf)
			getchar();
	}
}



int mainf2()
{
	int n = 1;
	while (true) {
		printf("%d ", n);
		size_t len = abs(rand_double(6,6));
		int e = rand_between(-100, 100);
		bigfp x, test, gt;
		rand_bigfp(x, len, -100, 100);
		int64_t lo_exp = x.exp - len;
		int64_t hi_exp = x.exp - 1;
		int64_t p = rand_between(lo_exp - 2, hi_exp + 2);
		
		mul_pow2(test, x, p);
		gt = x;
		mul_pow2(gt, p);

		if (compare(test, gt)!=0) {
			getchar();
		}
		printf("PASS!\n");

		n += 1;
	}
	while (true) {
		double s;
		scanf("%lf", &s);
		bigfp bs, br, b1;
		from_double(bs, s);
		print(bs);
		sqrt(br, bs, 10);
		print(br);
		bigfp br2;
		mul(br2, br, br);
		print(br2);
	}
}

int mainf3()
{
	bigfp a, b, c;
	size_t p = 1 << 20;
	rand_bigfp(a, p, 1, 1);
	while (true) {
		rsqrt_newton(c, a, p);
		printf("========\n");
	}
	bigfp half(0.5);
	for (size_t p = 15; p <= 1 << 25; p *= 2) {
		printf("================\n");
		printf("p = %zd\n", p);
		rand_bigfp(a, p, 1, 1);
		c.data.reserve(p+1);
		//rand_bigfp(b, p, 1, 1);

		timer("mul:        %f sec\n", [&] {mul_pow2(c, a, 100); });
		timer("mul_copy:   %f sec\n", [&] {mul_pow2(c, bigfp(a), 100); });
		timer("mul_self:   %f sec\n", [&] {mul_pow2(a, 100); });
		//timer("mul:      %f sec\n", [&] {mul(c, a, b); });
		//timer("mul1:     %f sec\n", [&] {mul(c, a, half); });
		//timer("rcp:      %f sec\n", [&] {rcp_newton(c, a, p); });
		//timer("div:      %f sec\n", [&] {div(c, a, b, p); });
		//timer("rsqrt:    %f sec\n", [&] {rsqrt_newton(c, a, p); });
		//timer("sqrt:     %f sec\n", [&] {sqrt(c, a, p); });
	}
	return 0;
}