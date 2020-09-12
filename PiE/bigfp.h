#pragma once
#include <vector>
#include <stdint.h>
#include <utility>
#include "int128.h"
#include "bigint.h"
#include "readonly_span.h"
#include "array_helpers.h"

struct bigfp;
struct bigfp_src;


void add(bigfp& c, const bigfp_src& a, const bigfp_src& b);

void sub(bigfp& c, const bigfp_src& a, const bigfp_src& b);

int compare(const bigfp_src& a, const bigfp_src& b);
bool is_less_or_equal(const bigfp_src& a, const bigfp_src& b);
bool is_equal(const bigfp_src& a, const bigfp_src& b);

void mul(bigfp& c, const bigfp_src& a, const bigfp_src& b);
void mul_pow2(bigfp& c, const bigfp_src& a, int64_t e);
void mul_pow2(bigfp& c, int64_t e);

void rcp_newton(bigfp& r, const bigfp_src& d, size_t p);
double div_apxm(const bigfp_src& a, const bigfp_src& b);
double validate_rcp(const bigfp_src& r, const bigfp_src& d, size_t p);

void div(bigfp& q, const bigfp_src& dividend, const bigfp_src& divisor, size_t p);
double validate_div(const bigfp_src& q, const bigfp_src& dividend, const bigfp_src& divisor, size_t p);

void rcp_abs_prec(bigfp& r, const bigfp_src& d, int64_t p);
double validate_rcp_abs_prec(const bigfp_src& r, const bigfp_src& d, int64_t p);

void div_abs_prec(bigfp& q, const bigfp_src& dividend, const bigfp_src& divisor, int64_t p);
double validate_div_abs_prec(const bigfp_src& q, const bigfp_src& dividend, const bigfp_src& divisor, int64_t p);


void rsqrt_newton(bigfp& r, const bigfp_src& s, size_t p);
double validate_rsqrt(const bigfp_src& r, const bigfp_src& s, size_t p);

void sqrt(bigfp& r, const bigfp_src& s, size_t p);
void sqrt(bigfp& r, size_t p);
double validate_sqrt(const bigfp_src& r, const bigfp_src& s, size_t p);


void div_rem(bigfp& quotient, bigfp& remainder, const bigfp_src& dividend, const bigfp_src& divisor);
void precompute_rcp(bigfp& rcp, const bigfp_src& dvs, int64_t dvd_n_max);
void div_rem_with_rcp(bigfp& q, bigfp& r, const bigfp_src& dvd, const bigfp_src& dvs, const bigfp& rcp);

void round(bigfp& r, const bigfp_src& a, int128 prec);
void round(bigfp& a, int128 prec);
void trunc(bigfp& r, const bigfp_src& a, int128 prec);
void trunc(bigfp& a, int128 prec);
void to_bigint(bigint& i, const bigfp_src& a);
bigfp_src integral_part(const bigfp_src& a);
bigfp_src fractional_part(const bigfp_src& a);
void free(bigfp& a);

void from_double(bigfp& a, double d);

double to_double(const bigfp_src& a);

void print(const bigfp_src& a);


/*
data = [w0, w1, w2, ..., w(n-1)]
B = 2^32
|value| = 0.w(n-1) w(n-2) ... w1 w0 * B^exp
w(n-1) != 0
w0 != 0
represent zero:
sign = 0, size(data)=0
*/
struct bigfp
{
	std::vector<uint32_t> data;
	int64_t exp;
	int sign;

	explicit bigfp(int i = 0)
	{
		if (i == 0) {
			sign = 0;
			return;
		}
		sign = i > 0 ? 1 : -1;
		data = { i > 0 ? (uint32_t)i : -(uint32_t)i };
		exp = 1;
	}

	explicit bigfp(double d)
	{
		from_double(*this, d);
	}

	/*
	explicit bigfp(const bigint& i)
	{
		if (i.size() == 0) {
			sign = 0;
		}
		else {
			sign = 1;
			exp = i.size();
			size_t loz = low_zeros(i);
			data.insert(data.begin(), i.begin() + loz, i.end());
		}
	}
	*/

	explicit bigfp(const bigfp_src& src);

	void operator=(const bigfp_src& src);
};

inline void swap(bigfp& a, bigfp& b)
{
	std::swap(a.data, b.data);
	std::swap(a.exp, b.exp);
	std::swap(a.sign, b.sign);
}


/*
struct bigfp_src
{
	const std::vector<uint32_t>& data;
	int64_t exp;
	int sign;

	bigfp_src(const bigfp& f) :data(f.data), exp(f.exp), sign(f.sign) {}

	// copyable, not movable, not assignable
	bigfp_src(const bigfp_src& f) = default;
	bigfp_src(bigfp_src&& f) = delete;
	bigfp_src& operator=(const bigfp_src&) = delete;
	bigfp_src& operator=(bigfp_src&&) = delete;
};
*/

struct bigfp_src
{
	const readonly_span<uint32_t> data;
	int64_t exp;
	int sign;

	bigfp_src(const bigfp& f) :data(f.data), exp(f.exp), sign(f.sign) {}

	bigfp_src(const bigint& i) :
		bigfp_src(i.data() + low_zeros(i),
				  i.data() + i.size(),
				  i.size(),
				  i.empty() ? 0 : 1)
	{
	}

	bigfp_src(const readonly_span<uint32_t>& span, int64_t exp, int sign) :
		data(span), exp(exp), sign(sign)
	{
	}

	// copyable, not assignable
	bigfp_src(const bigfp_src&) = default;
	bigfp_src& operator=(const bigfp_src&) = delete;
	bigfp_src& operator=(bigfp_src&&) = delete;

private:
	bigfp_src(const uint32_t* begin, const uint32_t* end, int64_t exp, int sign) :
		data(begin, end - begin), exp(exp), sign(sign)
	{
	}
};

inline bigfp::bigfp(const bigfp_src& src) :data(src.data.begin(), src.data.end()), exp(src.exp), sign(src.sign) {}

inline void bigfp::operator=(const bigfp_src& src)
{
	data.clear();
	data.insert(data.end(), src.data.begin(), src.data.end());
	exp = src.exp;
	sign = src.sign;
}

/*
detect if a references b
*/
inline bool is_from(const bigfp_src& a, const bigint& b)
{
	return is_from(a.data, b);
}

inline bool is_from(const bigfp_src& a, const bigfp& b)
{
	return is_from(a.data, b.data);
}

inline bool is_non_neg_integer(const bigfp_src& a)
{
	return a.sign == 0 || a.sign > 0 && safe_compare(a.exp, a.data.size()) >= 0;
}

inline const bigfp bigfp_0(0);
inline const bigfp bigfp_1(1);


