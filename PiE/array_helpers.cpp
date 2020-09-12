#include <stdint.h>
#include <algorithm>
#include <vector>
#include <string.h>
#include <assert.h>

void copy(uint32_t dst[], const uint32_t src[], size_t n)
{
	if (n == 0)
		return;
	memcpy(dst, src, n * sizeof(uint32_t));
}

void set(uint32_t dst[], uint32_t val, size_t n)
{
	std::fill_n(dst, n, val);
}

/*
	0				n-1
	vvvvvvvvxxxxxxxxx
	| len  |
*/
size_t count_val_low(const uint32_t arr[], size_t n, uint32_t val)
{
	for (size_t i = 0; i < n; ++i) {
		if (arr[i] != val)
			return i;
	}
	return n;
}

/*
	0				n-1
	xxxxxxxxxvvvvvvvv
			 | len  |
*/
size_t count_val_high(const uint32_t arr[], size_t n, uint32_t val)
{
	if (n == 0)
		return 0;
	size_t i = n;
	do {
		i -= 1;
		if (arr[i] != val) {
			return n - 1 - i;
		}
	} while (i != 0);
	return n;
}


size_t low_zeros(const uint32_t arr[], size_t n)
{
	return count_val_low(arr, n, 0);
}


size_t high_zeros(const uint32_t arr[], size_t n)
{
	return count_val_high(arr, n, 0);
}

/*
copy src to dst, then left-shift whole dst (shift in 0 bits)
shift_bits < 32
*/
void copy_shift_left(uint32_t dst[], const uint32_t src[], size_t n, size_t shift_bits)
{
	assert(shift_bits < 32);
	if (shift_bits == 0) {
		copy(dst, src, n);
		return;
	}
	uint32_t prev = 0;
	for (size_t i = 0; i < n; ++i) {
		uint32_t cur = src[i];
		dst[i] = (cur << shift_bits) | (prev >> (32 - shift_bits));
		prev = cur;
	}
}

/*
copy src to dst, then right-shift whole dst (shift in 0 bits)
shift_bits < 32
*/
void copy_shift_right(uint32_t dst[], const uint32_t src[], size_t n, size_t shift_bits)
{
	assert(shift_bits < 32);
	if (shift_bits == 0) {
		copy(dst, src, n);
		return;
	}
	uint32_t prev = 0;
	for (size_t i = n - 1; i < n; --i) {
		uint32_t cur = src[i];
		dst[i] = (cur >> shift_bits) | (prev << (32 - shift_bits));
		prev = cur;
	}
}

/*
left-shift whole arr (shift in 0 bits)
shift_bits < 32
*/
void shift_left(uint32_t arr[], size_t n, size_t shift_bits)
{
	if (shift_bits == 0)
		return;
	copy_shift_left(arr, arr, n, shift_bits);
}

/*
right-shift whole arr (shift in 0 bits)
shift_bits < 32
*/
void shift_right(uint32_t arr[], size_t n, size_t shift_bits)
{
	if (shift_bits == 0)
		return;
	copy_shift_right(arr, arr, n, shift_bits);
}

size_t low_zeros(const std::vector<uint32_t>& v)
{
	if (v.size() == 0)
		return 0;
	return count_val_low(&v[0], v.size(), 0);
}

size_t high_zeros(const std::vector<uint32_t>& v)
{
	if (v.size() == 0)
		return 0;
	return count_val_high(&v[0], v.size(), 0);
}
