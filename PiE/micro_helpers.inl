#pragma once
#include <intrin.h>
#include <algorithm>
#include <string.h>


template <class uint_t,
	std::enable_if_t<std::is_unsigned_v<uint_t>, int> dummy = 0
>
inline unsigned leading_zeros(uint_t x)
{
	constexpr unsigned bits = std::numeric_limits<uint_t>::digits;
	static_assert(bits <= 64);
	if (bits <= 32) {
		unsigned lz32 = leading_zeros((uint32_t)x);
		return lz32 - (32 - bits);
	}
	else {
		unsigned lz64 = leading_zeros((uint64_t)x);
		return lz64 - (64 - bits);
	}
}


template <>
inline unsigned leading_zeros(uint64_t x)
{
	unsigned long index;
	bool not_zero = _BitScanReverse64(&index, x);
	if (!not_zero)
		return 64;
	else
		return 63 - index;
}

template <>
inline unsigned leading_zeros(uint32_t x)
{
	unsigned long index;
	bool not_zero = _BitScanReverse(&index, x);
	if (!not_zero)
		return 32;
	else
		return 31 - index;
}

template <class uint_t,
	std::enable_if_t<std::is_unsigned_v<uint_t>, int> dummy = 0
>
inline unsigned trailing_zeros(uint_t x)
{
	constexpr unsigned bits = std::numeric_limits<uint_t>::digits;
	static_assert(bits <= 64);
	if (bits <= 32) {
		return trailing_zeros((uint32_t)x);
	}
	else {
		return trailing_zeros((uint64_t)x);
	}
}


template <>
inline unsigned trailing_zeros(uint64_t x)
{
	unsigned long index;
	bool not_zero = _BitScanForward64(&index, x);
	if (!not_zero)
		return 64;
	else
		return index;
}

template <>
inline unsigned trailing_zeros(uint32_t x)
{
	unsigned long index;
	bool not_zero = _BitScanForward(&index, x);
	if (!not_zero)
		return 32;
	else
		return index;
}

/*
return 0 for n=0
*/
inline size_t ceil_log2(size_t n)
{
	unsigned long index;
	bool not_zero = _BitScanReverse64(&index, n - 1);
	if (n == 0 || !not_zero)
		return 0;
	else
		return index + 1;
}

/*
smallest power of 2 >= n
return 0 for n > 2^63
*/
inline size_t clp2(size_t n)
{
	size_t shift = ceil_log2(n);
	size_t size_t_bits = CHAR_BIT * sizeof(size_t);
	if (shift >= size_t_bits)
		return 0;
	return 1ULL << shift;
}

inline bool is_pow2(size_t n)
{
	return (n&(n - 1)) == 0 && n > 0;
}

inline int64_t round_to_int64(double x)
{
	x = x + 0.5 - (x < 0);
	return (int64_t)x;
}

/*
	floor(dividend/divisor)
*/
inline int64_t floor_div(int64_t dividend, int64_t divisor)
{
	int64_t q = dividend / divisor;
	int64_t r = dividend % divisor;

	if (divisor > 0)
		return r >= 0 ? q : q - 1;
	else
		return r <= 0 ? q : q - 1;
}

/*
	dividend - floor(dividend/divisor)*divisor
*/
inline int64_t floor_rem(int64_t dividend, int64_t divisor)
{
	int64_t r = dividend % divisor;
	if (divisor > 0)
		return r >= 0 ? r : r + divisor;
	else
		return r <= 0 ? r : r + divisor;
}

/*
reverse the low bits, add 1, reverse again
*/

inline size_t rev_inc(size_t x, size_t bits)
{
	size_t hi_bits = 64 - bits;
	x <<= hi_bits;
	unsigned long zero_index;
	unsigned char all_1 = !_BitScanReverse64(&zero_index, ~x);
	if (all_1)
		return 0;
	else {
		x ^= ((size_t)(-1) << zero_index);
		x >>= hi_bits;
		return x;
	}
}

template < class to_int_t, class from_int_t,
	std::enable_if_t<std::is_integral_v<to_int_t> && std::is_integral_v<from_int_t>, int> dummy = 0
>
inline to_int_t checked_cast(from_int_t x)
{
	using to_limits = std::numeric_limits<to_int_t>;

	to_int_t min = to_limits::min();
	to_int_t max = to_limits::max();

	if (std::is_signed_v<from_int_t> && std::is_unsigned_v<to_int_t>) {
		using from_uint_t = std::make_unsigned_t<from_int_t>;
		assert(x >= 0 && (from_uint_t)x <= max);
	}
	else if (std::is_unsigned_v<from_int_t> && std::is_signed_v<to_int_t>) {
		using to_uint_t = std::make_unsigned_t<to_int_t>;
		assert(x <= (to_uint_t)max);
	}
	else { // same signedness
		assert(x >= min && x <= max);
	}
	return (to_int_t)x;
}

template < class int_t1, class int_t2,
	std::enable_if_t<std::is_integral_v<int_t1> && std::is_integral_v<int_t2>, int> dummy = 0
>
inline int safe_compare(int_t1 a, int_t2 b)
{
	if (std::is_signed_v<int_t1> && std::is_unsigned_v<int_t2>) {
		using uint_t1 = std::make_unsigned_t<int_t1>;
		if (a < 0)
			return -1;
		return safe_compare((uint_t1)a, b);;
	}
	else if (std::is_unsigned_v<int_t1> && std::is_signed_v<int_t2>) {
		using uint_t2 = std::make_unsigned_t<int_t2>;
		if (b < 0)
			return 1;
		return safe_compare(a ,(uint_t2)b);
	}
	else { // same signedness
		if (a < b)
			return -1;
		return a > b ? 1 : 0;
	}
}


/*
signed equivalent modulo 2^b  (b = bits of unsigned_t)
because unsigned to signed cast is UB
*/
template < typename unsigned_t,
	std::enable_if_t<std::is_unsigned_v<unsigned_t>, int> dummy = 0
>
inline std::make_signed_t<unsigned_t> unsigned_to_signed(unsigned_t u)
{
	using signed_t = std::make_signed_t<unsigned_t>;
	using signed_limits = std::numeric_limits<signed_t>;

	if (u <= (unsigned_t)signed_limits::max()) {
		return (signed_t)u;
	}
	else { // u >= 0x80000000

		// alternative way to compute -(~u+1)
		// avoids underflow when u = 0x80000000

		signed_t inv = ~u;   // [0, 0x7fffffff]
		return -inv - 1;
		// in a non-2's complement system ( -INT_MAX = INT_MIN )
		// when u = 0x80000000
		// it underflows anyway
	}
}

template < typename to_int_t, typename from_uint_t,
	std::enable_if_t<std::is_unsigned_v<from_uint_t> && std::is_signed_v<to_int_t>, int> dummy = 0
>
inline to_int_t checked_negate(from_uint_t u)
{
	using to_uint_t = std::make_unsigned_t<to_int_t>;
	to_int_t to_int_min = std::numeric_limits<to_int_t>::min();

	// -u >= to_int_min
	// u <= -(2^b + to_int_min) + 2^b	b = bits of to_uint_t
	assert(u <= -(to_uint_t)to_int_min);
	// u < 2^b
	// -u mod 2^b = 
	//   u = 0: 0       = -(to_uint_t)u
	//   u > 0: 2^b - u = -(to_uint_t)u
	return unsigned_to_signed(-(to_uint_t)u);
}


template < class uint_t,
	std::enable_if_t<std::is_unsigned_v<uint_t>, int> dummy = 0
>
inline uint_t shift_left_checked(uint_t v, unsigned shift_amt)
{
	unsigned bits = std::numeric_limits<uint_t>::digits;
	assert(shift_amt < bits);
	uint_t shifted = v << shift_amt;
	assert((shifted >> shift_amt) == v);
	return shifted;
}