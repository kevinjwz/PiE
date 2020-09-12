#pragma once
#include <stdint.h>
#include <intrin.h>
#include <limits>
#include <type_traits>
#include "micro_helpers.inl"


struct int128
{
	// value = hi * 2^64 + lo
	uint64_t lo;
	int64_t hi;

	int128() {}

	int128(int64_t x)
	{
		lo = (uint64_t)x;
		hi = x >= 0 ? 0 : -1;
	}

	int128(uint64_t x)
	{
		lo = x;
		hi = 0;
	}

	int128(int32_t x) :int128((int64_t)x) {}

	int128(uint32_t x) :int128((uint64_t)x) {}

	template < typename int_t,
				std::enable_if_t<std::is_integral_v<int_t>, int> dummy = 0 >
	explicit operator int_t() const
	{
		if (std::is_signed_v<int_t>) {
			using unsigned_t = std::make_unsigned_t<int_t>;
			return unsigned_to_signed((unsigned_t)lo);
		}
		else {
			return (int_t)lo;
		}

	}
};

template < typename int_t,
	std::enable_if_t<std::is_integral_v<int_t>, int> dummy = 0 >
inline int_t checked_cast(const int128& a)
{
	using limits = typename std::numeric_limits<int_t>;

	int_t min = limits::min();
	int_t max = limits::max();
	assert(a >= min && a <= max);
	return (int_t)a;
}


inline int128 operator+(const int128& a, const int128& b)
{
	int128 r;
	r.lo = a.lo + b.lo;
	int64_t carry = r.lo < a.lo;
	r.hi = a.hi + b.hi + carry;
	return r;
}

inline int128 operator-(const int128& a, const int128& b)
{
	int128 r;
	r.lo = a.lo - b.lo;
	int64_t borrow = a.lo < b.lo;
	r.hi = a.hi - b.hi - borrow;
	return r;
}

inline int128 operator*(const int128& a, const int128& b)
{
	int128 r;
	r.lo = _umul128(a.lo, b.lo, (uint64_t*)&r.hi);
	r.hi += a.hi*b.lo + a.lo*b.hi;
	return r;
}

inline int128 operator-(const int128& a)
{
	return 0 - a;
}

inline bool operator<(const int128& a, const int128& b)
{
	bool hi_le = a.hi <= b.hi;
	bool hi_lt = a.hi < b.hi;
	return a.lo<b.lo ? hi_le : hi_lt;
}

inline bool operator>(const int128& a, const int128& b)
{
	return b < a;
}

inline bool operator<=(const int128& a, const int128& b)
{
	return !(a > b);
}

inline bool operator>=(const int128& a, const int128& b)
{
	return !(a < b);
}

inline bool operator==(const int128& a, const int128& b)
{
	bool hi_eq = a.hi == b.hi;
	bool lo_eq = a.lo == b.lo;

	return hi_eq & lo_eq;
}

inline bool operator!=(const int128& a, const int128& b)
{
	return !(a == b);
}
