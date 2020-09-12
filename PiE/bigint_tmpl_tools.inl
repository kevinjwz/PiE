#pragma once
#include <type_traits>
#include <limits>
#include <assert.h>
#include "bigint.h"
#include "micro_helpers.inl"

template < class uint_t,
	std::enable_if_t<std::is_unsigned_v<uint_t>, int> dummy = 0
>
inline void from_uint(bigint& a, uint_t u)
{
	a.clear();
	if (u == 0)
		return;

	unsigned bits = std::numeric_limits<uint_t>::digits;
	a.reserve((bits + 31) / 32);
	if (bits <= 32) {
		a.push_back((uint32_t)u);
	}
	else {
		while (u > 0) {
			a.push_back((uint32_t)u);
			u >>= 32;
		}
	}
}

template < class uint_t,
	std::enable_if_t<std::is_unsigned_v<uint_t>, int> dummy = 0
>
inline uint_t to_uint_checked(const bigint& a)
{
	unsigned bits = std::numeric_limits<uint_t>::digits;
	size_t n = a.size();
	if (n == 0)
		return 0;
	if (bits <= 32) {
		assert(n == 1);
		uint_t max = std::numeric_limits<uint_t>::max();
		assert(a[0] <= max);
		return (uint_t)a[0];
	}
	else {
		// n>=1
		uint_t u = a[n - 1];
		// n-2 downto 0
		for (size_t i = n - 1; i > 0;) {
			i -= 1;
			u = shift_left_checked(u, 32);
			u += a[i];
		}
		return u;
	}
}