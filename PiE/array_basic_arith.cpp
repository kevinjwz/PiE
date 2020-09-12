#include <stdint.h>
#include <assert.h>
#include "array_helpers.h"

/*
return: carry
c is a ... ok
c is b ... ok
c is a and b ... ok
*/
uint32_t add_n(uint32_t c[], const uint32_t a[], const uint32_t b[], size_t n)
{
	uint64_t carry = 0;
	for (size_t i = 0; i < n; ++i) {
		uint64_t sum = (uint64_t)a[i] + b[i] + carry;
		c[i] = (uint32_t)sum;
		carry = sum >> 32;
	}
	return (uint32_t)carry;
}

/*
c = a - b
return: borrow
c is a ... ok
c is b ... ok
c is a and b ... ok
*/
uint32_t sub_n(uint32_t c[], const uint32_t a[], const uint32_t b[], size_t n)
{
	uint64_t carry = 1;
	for (size_t i = 0; i < n; ++i) {
		uint64_t sum = (uint64_t)a[i] + (uint64_t)~b[i] + carry;
		c[i] = (uint32_t)sum;
		carry = sum >> 32;
	}
	return 1 - (uint32_t)carry;
}

/*
n must >= 1
return: carry
*/
uint32_t add_1(uint32_t a[], uint32_t b, size_t n)
{
	assert(n >= 1);
	uint32_t a0 = a[0];
	a[0] = a0 + b;
	if (a[0] >= b)
		return 0;

	for (size_t i = 1; i < n; ++i) {
		a[i] += 1;
		if (a[i] != 0) {
			return 0;
		}
	}
	return 1;
}


/*
n must >= 1

return:
borrow: 1
no borrow: 0
*/
uint32_t sub_1(uint32_t a[], uint32_t b, size_t n)
{
	assert(n >= 1);
	uint32_t a0 = a[0];
	a[0] = a0 - b;
	if (a0 >= b)
		return 0;

	for (size_t i = 1; i < n; ++i) {
		if (a[i] != 0) {
			a[i] -= 1;
			return 0;
		}
		else {
			a[i] = (uint32_t)-1;
		}
	}
	return 1;
}


// TODO: optimize copy_add_1, copy_sub_1
/*
c = a + b
n must >= 1
return: carry
*/
uint32_t copy_add_1(uint32_t c[], const uint32_t a[], uint32_t b, size_t n)
{
	assert(n >= 1);
	c[0] = a[0] + b;
	uint64_t carry = c[0] < b;
	for (size_t i = 1; i < n; ++i) {
		uint64_t sum = (uint64_t)a[i] + carry;
		c[i] = (uint32_t)sum;
		carry = sum >> 32;
	}
	return (uint32_t)carry;
}

/*
c = a - b
n must >= 1
return: borrow
*/
uint32_t copy_sub_1(uint32_t c[], const uint32_t a[], uint32_t b, size_t n)
{
	assert(n >= 1);
	c[0] = a[0] - b;
	uint64_t carry = a[0] >= b;
	for (size_t i = 1; i < n; ++i) {
		uint64_t sum = (uint64_t)a[i] + 0xffffffff + carry;
		c[i] = (uint32_t)sum;
		carry = sum >> 32;
	}
	return 1 - (uint32_t)carry;
}

/*
return: carry
bn <= an
b does NOT overlap with a
or b>=a (address)
*/
uint32_t acc_short(uint32_t a[], size_t an, const uint32_t b[], size_t bn)
{
	assert(bn <= an);
	uint32_t carry = add_n(a, a, b, bn);
	if (an > bn)
		return add_1(a + bn, carry, an - bn);
	else
		return carry;
}

/*
size(c) = an
bn <= an
return: carry
*/
uint32_t add_short(uint32_t c[], const uint32_t a[], size_t an, const uint32_t b[], size_t bn)
{
	assert(bn <= an);
	uint32_t carry = add_n(c, a, b, bn);
	if (an > bn)
		return copy_add_1(c + bn, a + bn, carry, an - bn);
	else
		return carry;
}


/*
return: borrow
bn <= an
b does NOT overlap with a
or b>=a (address)
*/
uint32_t sub_short(uint32_t a[], size_t an, const uint32_t b[], size_t bn)
{
	assert(bn <= an);
	uint32_t borrow = sub_n(a, a, b, bn);
	if (an > bn)
		return sub_1(a + bn, borrow, an - bn);
	else
		return borrow;
}

/*
size(c) = an
bn <= an
return: borrow
*/
uint32_t sub_short(uint32_t c[], const uint32_t a[], size_t an, const uint32_t b[], size_t bn)
{
	assert(bn <= an);
	uint32_t borrow = sub_n(c, a, b, bn);
	if (an > bn)
		return copy_sub_1(c + bn, a + bn, borrow, an - bn);
	else
		return borrow;
}


/*
cmp =
	a>b: 1
	a=b: 0
	a<b: -1

	n-1                  0
a =   12345678901234567890
b =   12345666666666666666
			^
			diffidx
*/

void cmp_diff_index_max(int& cmp, size_t& diffidx, const uint32_t a[], const uint32_t b[], size_t n)
{
	for (size_t i = n - 1; i != -1; --i) {
		if (a[i] > b[i]) {
			cmp = 1;
			diffidx = i;
			return;
		}
		else if (a[i] < b[i]) {
			cmp = -1;
			diffidx = i;
			return;
		}
	}
	cmp = 0;
}


/*
c = |a - b|
return: sign
a>b: sign = 1
a=b: sign = 0
a<b: sign = -1
*/
int sub_sign_mag(uint32_t c[], const uint32_t a[], const uint32_t b[], size_t n)
{
	int cmp;
	size_t idxmax;
	cmp_diff_index_max(cmp, idxmax, a, b, n);
	if (cmp == 0) {
		set(c, 0, n);
		return 0;
	}

	const uint32_t* max = cmp > 0 ? a : b;
	const uint32_t* min = cmp > 0 ? b : a;
	size_t sub_len = idxmax + 1;

	uint32_t borrow = sub_n(c, max, min, sub_len);
	assert(borrow == 0);
	set(c + sub_len, 0, n - sub_len);
	return cmp;
}

/*
c = a * b
size(c)=size(a)=n
return: carry
*/
static uint32_t mul_1(uint32_t c[], const uint32_t a[], uint32_t b, size_t n)
{
	// carry < 2^32
	uint64_t carry = 0;
	for (size_t i = 0; i < n; ++i) {
		// guaranteed: prod < 2^64
		uint64_t prod = (uint64_t)a[i] * b + carry;
		c[i] = (uint32_t)prod;
		carry = prod >> 32;
	}
	return (uint32_t)carry;
}

/*
c += a * b
size(c)=size(a)=n
return: carry
*/
static uint32_t mul_1_acc(uint32_t c[], const uint32_t a[], uint32_t b, size_t n)
{
	// carry < 2^32
	uint64_t carry = 0;
	for (size_t i = 0; i < n; ++i) {
		// guaranteed: prod_acc < 2^64
		uint64_t prod_acc = (uint64_t)a[i] * b + c[i] + carry;
		c[i] = (uint32_t)prod_acc;
		carry = prod_acc >> 32;
	}
	return (uint32_t)carry;
}

/*
c = a * b
size(a) = an
size(b) = bn
size(c) = an+bn
*/

void mul_naive_ub(uint32_t c[], const uint32_t a[], size_t an, const uint32_t b[], size_t bn)
{
	if (bn == 0) {
		set(c, 0, an);
		return;
	}
	uint32_t carry = mul_1(c, a, b[0], an);
	c[an] = carry;
	for (size_t bi = 1; bi < bn; ++bi) {
		carry = mul_1_acc(c + bi, a, b[bi], an);
		c[bi + an] = carry;
	}
}