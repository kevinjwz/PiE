#include<assert.h>
#include<random>
#include "dcomplex.inl"
#include "fft_rad2_rad4.h"
#include "fft_constants.h"
#include "timer.h"

template<bool inverse>
static void dft(dcomplex result[], const dcomplex x[], size_t n)
{
	dcomplex* w = new dcomplex[n];
	double coef = inverse ? 1.0 / n : 1;
	double dir = inverse ? 1 : -1;
	for (size_t wi = 0; wi < n; ++wi) {
		double angle = dir * (2 * PI) * ((double)wi / n);
		w[wi] = { cos(angle), sin(angle) };
	}
	for (size_t ri = 0; ri < n; ++ri) {
		dcomplex sum = 0;
		size_t wi = 0;
		for (size_t xi = 0; xi < n; ++xi) {
			sum += x[xi] * w[wi];
			wi += ri;
			if (wi >= n || wi < ri)
				wi -= n;
		}
		result[ri] = sum * coef;
	}
	delete[] w;
}



static std::mt19937 rng;
void rand_complex_vector(std::vector<dcomplex>& vec, size_t n)
{
	vec.clear();
	vec.resize(n);
	for (size_t i = 0; i < n; ++i) {
		vec[i] = { (double)rng(), (double)rng() };
	}
}



int mainx()
{
	while (true) {
		printf("================\n");
		size_t n;
		scanf_s("%zd", &n);
		if (n == 0)
			break;
		assert((n&(n - 1)) == 0);	// n is power of 2
		std::vector<dcomplex> x, fx1, fx2;
		rand_complex_vector(x, n);
		//x.resize(n);
		//for (size_t i = 0; i < n; ++i) {
		//	double re;
		//	scanf_s("%lf", &re);
		//	x[i] = { re, 0 };
		//}
		fx1.resize(n);
		fx2.resize(n);
		dft<false>(&fx1[0], &x[0], n);
		fx2 = x;
		fft_mm_radix_4<false>(&fx2[0], n);
		printf("%24s%24s\n", "dft", "fft");
		for (size_t i = 0; i < n; ++i) {
			printf("%24.16g%24.16g\n", fx1[i].real(), fx2[i].real());
		}
	}
	return  0;
}



int mainxx()
{
	//size_t n = 1 << 24;
	//std::vector<dcomplex> x(n), fx(n);
	//rand_complex_vector(x, n);
	//while (true) {
	//	printf("================\n");
	//	fft_mm_radix_4<false>(&fx[0], &x[0], n);
	//}
	for (size_t n = 256; n <= 500000000; n *= 2) {
		printf("================\n");
		printf("n = %zd\n", n);
		std::vector<dcomplex> x, fx;
		rand_complex_vector(x, n);
		fx.resize(n);
		printf("[forward]\n");
		fx = x;
		timer("fft_mm2:    %f sec\n", [&] {fft_mm_radix_2<false>(&fx[0], n); });
		fx = x;
		timer("fft_mm4:    %f sec\n", [&] {fft_mm_radix_4<false>(&fx[0], n); });
		printf("[inverse]\n");
		x = fx;
		timer("fft_mm2:    %f sec\n", [&] {fft_mm_radix_2<true>(&x[0], n); });
		x = fx;
		timer("fft_mm4:    %f sec\n", [&] {fft_mm_radix_4<true>(&x[0], n); });
	}
	return 0;
}

void compute_error(double& avg_err, double& max_err, const std::vector<dcomplex>& v, const std::vector<dcomplex> gt)
{
	using namespace std;
	assert(v.size() == gt.size());
	double sum = 0;
	double max_ = 0;
	for (size_t i = 0; i < v.size(); ++i) {
		dcomplex diff = v[i] - gt[i];
		double err = sqrt(norm(diff) / norm(gt[i]));
		sum += err;
		max_ = max(max_, err);
	}
	avg_err = sum / v.size();
	max_err = max_;
}

//template<class func>
//void forward_inverse(std::vector<dcomplex>& r, const std::vector<dcomplex>& x, func fft_func)
//{
//	std::vector<dcomplex> fx(x.size());
//	fft_func
//}

int mainxxx()
{
	for (size_t n = 1; n <= 100000000; n *= 2) {
		printf("================\n");
		printf("n = %zd\n", n);
		std::vector<dcomplex> x, fx, ffx;
		rand_complex_vector(x, n);
		fx.resize(n);
		ffx.resize(n);
		printf("            %12s%12s\n", "avg", "max");
		double avg, max;

		std::vector<dcomplex> fx_gt(n);
		//dft<false>(&fx_gt[0], &x[0], n);

#define test_fwd_inv(fft) \
		ffx = x;\
		fft<false>(&ffx[0], n);\
		fft<true>(&ffx[0], n);\
		compute_error(avg, max, ffx, x);

#define test_fwd(fft) \
		fx = x;\
		fft<false>(&fx[0], n);\
		compute_error(avg, max, fx, fx_gt);



		//dft<true>(&ffx[0], &fx[0], n);
		//compute_error(avg, max, ffx, x);
		//printf("dft:        %12.2g%12.2g\n", avg, max);

		//test_fwd(dft)
		//	printf("dft:        %12.2g%12.2g\n", avg, max);

		//	printf("fft_nonrec: %12.2g%12.2g\n", avg, max);

		test_fwd_inv(fft_mm_radix_2)
			printf("fft_mm2:   %12.2g%12.2g\n", avg, max);

		test_fwd_inv(fft_mm_radix_4)
			printf("fft_mm4:    %12.2g%12.2g\n", avg, max);

	}
	return 0;
}

