#pragma once
#include <type_traits>
#include <limits>
#include <assert.h>
#include "bigfp.h"
#include "bigint_tmpl_tools.inl"
#include "micro_helpers.inl"

template < class uint_t,
	std::enable_if_t<std::is_unsigned_v<uint_t>, int> dummy = 0
>
inline void from_uint(bigfp& a, uint_t u)
{
	if (u == 0) {
		a.data.clear();
		a.sign = 0;
		return;
	}
	a.sign = 1;

	unsigned bits = std::numeric_limits<uint_t>::digits;

	if (bits <= 32) {
		from_uint(a.data, u);
		a.exp = a.data.size();
	}
	else {
		a.exp = 0;
		while (true) {
			uint32_t lo32 = (uint32_t)u;
			if (lo32 != 0)
				break;
			a.exp += 1;
			u >>= 32;
		}
		from_uint(a.data, u);
		a.exp += a.data.size();
	}
}

template < class int_t,
	std::enable_if_t<std::is_integral_v<int_t>, int> dummy = 0
>
inline void from_int(bigfp& a, int_t i)
{
	using uint_t = std::make_unsigned_t<int_t>;
	if (std::is_unsigned_v<int_t> || i >= 0) {
		from_uint(a, (uint_t)i);
	}
	else {
		uint_t u = -(uint_t)i;	// u = 2^b-(2^b+i) = -i
		from_uint(a, u);
		a.sign = -1;
	}
}

template < class uint_t,
	std::enable_if_t<std::is_unsigned_v<uint_t>, int> dummy = 0
>
inline uint_t to_uint_checked(const bigfp_src& a)
{
	assert(a.sign >= 0);
	if (a.sign == 0 || a.exp <= 0)
		return 0;

	unsigned bits = std::numeric_limits<uint_t>::digits;
	if (bits <= 32) {
		assert(a.exp == 1);
		uint_t max = std::numeric_limits<uint_t>::max();
		uint32_t int_part = *a.data.rbegin();
		assert(int_part <= max);
		return (uint_t)int_part;
	}
	else {
		/*
			aaaaaaaa.aaaaaaaa
			^      ^
			|      e_min=0
			int_n-1

			or
			aaaaa000.
			^   ^
			|   e_min=3
			int_n-1
		*/
		size_t n = a.data.size();	// n>=1
		assert(a.exp <= SIZE_MAX);
		size_t int_n = (size_t)a.exp;	// int_n>=1
		size_t e_min = n >= int_n ? 0 : int_n - n;
		size_t e_to_i = n - int_n;

		uint_t u = a.data[n - 1];
		// int_n-2 downto 0
		for (size_t e = int_n - 1; e != 0;) {
			e -= 1;
			u = shift_left_checked(u, 32);
			if (e >= e_min)
				u += a.data[e + e_to_i];
		}
		return u;
	}
}

template < class int_t,
	std::enable_if_t<std::is_signed_v<int_t>, int> dummy = 0
>
inline int_t to_int_checked(const bigfp_src& a)
{
	using uint_t = std::make_unsigned_t<int_t>;

	if (a.sign >= 0) {
		uint_t u = to_uint_checked<uint_t>(a);
		return checked_cast<int_t>(u);
	}
	else {
		bigfp_src abs_a(a);
		abs_a.sign = -a.sign;
		uint_t u = to_uint_checked<uint_t>(abs_a);
		return checked_negate<int_t>(u);
	}
}