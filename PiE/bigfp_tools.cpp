#include<assert.h>
#include<math.h>
#include<algorithm>
#include"bigint.h"
#include"bigfp.h"
#include"int128.h"
#include"micro_helpers.inl"
#include"array_helpers.h"
#include"array_basic_arith.h"

namespace // classes only visible in this file
{
	struct self_modifier
	{
		bigfp& a;
		self_modifier(bigfp& a) :a(a) {}
		void gen_same() const {}
		void gen_zero() const
		{
			a = bigfp_0;
		};
		void gen_pow2(int64_t exp) const
		{
			// a = a.sign * (2^32)^exp
			a.data = { 1 };
			a.exp = checked_cast<int64_t>((int128)exp + 1);
		};
		void do_trunc(size_t trunc_len) const
		{
			a.data.erase(a.data.begin(), a.data.begin() + trunc_len);
		};
		void do_trunc_add_1(size_t trunc_len) const
		{
			do_trunc(trunc_len);
			a.data[0] += 1;
		};
	};

	struct copy_modifier
	{
		bigfp& r;
		const bigfp_src& a;
		copy_modifier(bigfp& r, const bigfp_src& a) :r(r), a(a) {}
		void gen_same() const
		{
			r = a;
		}
		void gen_zero() const
		{
			r = bigfp_0;
		};
		void gen_pow2(int64_t exp) const
		{
			// r = a.sign * (2^32)^exp
			r.sign = a.sign;
			r.data = { 1 };
			r.exp = checked_cast<int64_t>((int128)exp + 1);
		};
		void do_trunc(size_t trunc_len) const
		{
			r.sign = a.sign;
			r.exp = a.exp;
			r.data.clear();
			r.data.insert(r.data.begin(), a.data.begin() + trunc_len, a.data.end());
		};
		void do_trunc_add_1(size_t trunc_len) const
		{
			do_trunc(trunc_len);
			r.data[0] += 1;
		};
	};
}

static const uint32_t HALF = 1u << 31;

static int cmp_half(const uint32_t a[], size_t n)
{
	if (a[n - 1] < HALF)
		return -1;
	if (a[n - 1] > HALF)
		return 1;
	/*
		a[0]~a[n-2] cannot be all zeros
	*/
	return n > 1 ? 1 : 0;
}

/*
round a to multiple of B^prec
*/
template < class modifier_t >
static void round_skeleton(const bigfp_src& a, int128 prec, const modifier_t& modifier)
{
	if (a.sign == 0) {
		modifier.gen_same();
		return;
	}

	int128 lo_exp = (int128)a.exp - a.data.size();
	int128 hi_exp = (int128)a.exp - 1;
	if (prec <= lo_exp) {
		modifier.gen_same();
		return;
	}
	if (prec >= hi_exp + 2) {
		modifier.gen_zero();
		return;
	}

	if (prec == hi_exp + 1) {
		int cmp_r = cmp_half(&a.data[0], a.data.size());
		if (cmp_r <= 0) {
			modifier.gen_zero();
		}
		else {
			modifier.gen_pow2(a.exp);
		}
		return;
	}
	/* prec in [lo_exp+1, hi_exp]
		|  head    || tail |
		aaaaaaaaaaaaaaaaaaaa
		^          ^       ^
		hi_exp     prec    lo_exp
	*/
	size_t tail_n = (size_t)(prec - lo_exp);	// [1, size-1]
	size_t head_n = a.data.size() - tail_n;		// [1, size-1]
	int cmp_r = cmp_half(&a.data[0], tail_n);
	if (cmp_r < 0 || cmp_r == 0 && (a.data[tail_n] & 1) == 0) {
		/*
			aaaaaaaa0000aaaaaaaa
			^          ^       ^
			hi_exp     prec    lo_exp
		*/
		size_t head_loz = low_zeros(&a.data[tail_n], head_n);
		modifier.do_trunc(tail_n + head_loz);
	}
	else {
		/*
			aaaaaaaaffffaaaaaaaa
			^          ^       ^
			hi_exp     prec    lo_exp
		*/
		size_t head_lof = count_val_low(&a.data[tail_n], head_n, UINT32_MAX);
		if (head_lof == head_n) {
			modifier.gen_pow2(a.exp);
			return;
		}
		modifier.do_trunc_add_1(tail_n + head_lof);
	}
}

void round(bigfp& r, const bigfp_src& a, int128 prec)
{
	assert(!is_from(a, r));	
	round_skeleton(a, prec, copy_modifier(r, a));
}

void round(bigfp& a, int128 prec)
{
	round_skeleton(a, prec, self_modifier(a));
}


template < class modifier_t >
static void trunc_skeleton(const bigfp_src& a, int128 prec, const modifier_t& modifier)
{
	if (a.sign == 0) {
		modifier.gen_same();
		return;
	}
	int128 lo_exp = (int128)a.exp - a.data.size();
	int128 hi_exp = (int128)a.exp - 1;

	if (prec <= lo_exp) {
		modifier.gen_same();
		return;
	}
	if (prec >= hi_exp + 1) {
		modifier.gen_zero();
		return;
	}
	/* prec in [lo_exp+1, hi_exp]
		|  head    || tail |
		aaaaaaaa0000aaaaaaaa
		^          ^       ^
		hi_exp     prec    lo_exp
	*/
	size_t tail_n = (size_t)(prec - lo_exp);
	size_t head_n = a.data.size() - tail_n;
	size_t head_loz = low_zeros(&a.data[tail_n], head_n);
	modifier.do_trunc(tail_n + head_loz);
}

void trunc(bigfp& r, const bigfp_src& a, int128 prec)
{
	assert(!is_from(a, r));
	trunc_skeleton(a, prec, copy_modifier(r,a));
}

void trunc(bigfp& a, int128 prec)
{
	trunc_skeleton(a, prec, self_modifier(a));
}

void free(bigfp& a)
{
	std::vector<uint32_t> empty;
	a.data.swap(empty);
}

void from_double(bigfp& a, double d)
{
	assert(isfinite(d));

	if (fabs(d) == 0) {
		a.sign = 0;
		a.data = {};
		return;
	}

	a.sign = d > 0 ? 1 : -1;
	int exp2;
	double frac = frexp(fabs(d), &exp2);
	// frac in [0.5, 1)
	a.exp = ((int64_t)exp2 + 31) >> 5;	// ceil(exp2/32)
	int shift_r = (int)(a.exp * 32 - exp2);	// [0, 31]
	frac = ldexp(frac, -shift_r);

	a.data.clear();
	while (frac > 0) {
		double i;
		frac = modf(ldexp(frac, 32), &i);
		a.data.insert(a.data.begin(), (uint32_t)i);
	}
}

double to_double(const bigfp_src& a)
{
	if (a.sign == 0)
		return 0;

	size_t n = a.data.size();
	size_t lz = leading_zeros(*a.data.rbegin());	// [0, 31]
	double frac;
	if (53 + lz < 64) {
		uint32_t hi_2[2] = {
			n >= 2 ? a.data[n - 2] : 0,
			a.data[n - 1]
		};

		size_t tail_bits = 64 - (53 + lz);	// [1, 11]
		uint32_t tail_mask = (1 << tail_bits) - 1;

		frac = ldexp(hi_2[1], -32) + ldexp(hi_2[0] & ~tail_mask, -64);

		uint32_t tail_half = 1 << (tail_bits - 1);
		uint32_t tail = hi_2[0] & tail_mask;
		uint32_t add_tail = (tail >= tail_half) * 2 + ((tail > tail_half) | (n > 2));
		/*
			bool tail_gt_half = tail > tail_half || tail == tail_half && n > 2;
			bool tail_eq_half = tail == tail_half && n <= 2;
			bool tail_lt_half = tail < tail_half;

			add_tail =
				11: tail_gt_half
				10: tail_eq_half
				01, 00: tail_lt_half
		*/
		frac += ldexp(add_tail, -64 + (int)tail_bits - 2);
	}
	else {
		uint32_t hi_3[3] = {
			n >= 3 ? a.data[n - 3] : 0,
			n >= 2 ? a.data[n - 2] : 0,
			a.data[n - 1]
		};
		size_t tail_bits = 96 - (53 + lz);	// [12, 32]
		uint32_t tail_mask = (1 << tail_bits) - 1;

		frac = ldexp(hi_3[2], -32) + ldexp(hi_3[1], -64) + ldexp(hi_3[0] & ~tail_mask, -96);

		uint32_t tail_half = 1 << (tail_bits - 1);
		uint32_t tail = hi_3[0] & tail_mask;
		uint32_t add_tail = (tail >= tail_half) * 2 + ((tail > tail_half) | (n > 3));

		frac += ldexp(add_tail, -96 + (int)tail_bits - 2);
	}

	int128 exp2 = (int128)a.exp * 32;
	exp2 = std::min<int128>(INT_MAX, exp2);
	exp2 = std::max<int128>(INT_MIN, exp2);
	return a.sign * ldexp(frac, (int)exp2);
}

void print(const bigfp_src& a)
{
	if (a.sign < 0)
		printf("-");
	printf("0.");
	for (size_t i = a.data.size() - 1; i != SIZE_MAX; --i) {
		printf("%08x ", a.data[i]);
	}
	printf("B%+lld\n", a.exp);
}

/*
ignore sign, truncate fraction part
dst = floor(|src|)
*/
void to_bigint(bigint& dst, const bigfp_src& src)
{
	dst.clear();
	if (src.sign == 0 || src.exp<=0)
		return;
	
	/*
		a = aaaaaa.aaaaaaaaaa
	*/
	uint64_t int_n = (uint64_t)src.exp;	// integer part length
	assert(int_n <= SIZE_MAX);
	size_t int_n_s = (size_t)int_n;
	size_t n = src.data.size();
	if (int_n_s <= n) {
		dst.insert(dst.end(), src.data.begin() + (n - int_n_s), src.data.end());
	}
	else {
		size_t zero_n = int_n_s - n;
		dst.insert(dst.end(), zero_n, 0);
		dst.insert(dst.end(), src.data.begin(), src.data.end());
	}
}

bigfp_src integral_part(const bigfp_src& a)
{
	if (a.sign == 0 || a.exp <= 0) {
		return bigfp_0;
	}
	/*
	aaaaaa.aaaaaaaaaa
	aaa000.aaaaaaaaaa
	aaa000.
	*/

	size_t n = a.data.size();
	if (safe_compare(a.exp, n) >= 0) {
		return a;
	}
	else { // a.exp < n  =>  a.exp <= SIZE_MAX
		size_t int_n = (size_t)a.exp;
		size_t frac_n = n - int_n;
		size_t tz = low_zeros(&a.data[frac_n], int_n);
		return bigfp_src(
			readonly_span<uint32_t>(&a.data[frac_n + tz], int_n - tz),
			a.exp, a.sign);
	}
}

bigfp_src fractional_part(const bigfp_src& a)
{
	if (a.sign == 0)
		return bigfp_0;
	if (a.exp <= 0)
		return a;

	/*
	aaaaaa.aaaaaaaaaa
	aaaaaa.0000aaaaaa
	aaa000.
	*/

	size_t n = a.data.size();
	if (safe_compare(a.exp, n) >= 0) {
		return bigfp_0;
	}
	else { // a.exp < n  =>  a.exp <= SIZE_MAX
		size_t int_n = (size_t)a.exp;
		size_t frac_n = a.data.size() - int_n;
		size_t lz = high_zeros(&a.data[0], frac_n);
		int64_t exp = checked_negate<int64_t>(lz);
		return bigfp_src(
			readonly_span<uint32_t>(&a.data[0], frac_n - lz),
			exp, a.sign);
	}
}