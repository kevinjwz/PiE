#include <assert.h>
#include <algorithm>
#include "bigint.h"
#include "array_helpers.h"
#include "array_basic_arith.h"


/*
start >= a.size():      part = 0
start + len > a.size(): part = a[start .. (a.size()-1)]
*/

static void cut(const bigint& a, size_t start, size_t len, bigint& part)
{
	if (start >= a.size()) {
		part.resize(0);
		return;
	}
	
	len = std::min(len, a.size() - start);
	size_t end = start + len - 1;
	size_t ai;
	// end downto start
	for (ai = end; ai != start - 1; --ai) {
		if (a[ai] != 0)
			break;
	}
	part.resize(ai - start + 1);
	for (; ai != start - 1; --ai) {
		part[ai - start] = a[ai];
	}
}

/*
a = a + (b<<(shift_amt*32))

*/

static void shift_acc(bigint& a, const bigint& b, size_t shift_amt)
{
	size_t an = a.size(), bn = b.size();
	if (bn == 0)
		return;
	
	assert(bn + shift_amt >= bn);
	if (bn + shift_amt <= an) {
		assert(an + 1 != 0);
		a.resize(an + 1);
		uint32_t carry = acc_short(&a[shift_amt], an - shift_amt, &b[0], bn);
		if (carry == 0)
			a.resize(an);
		else
			a[an] = carry;
	}
	else if (shift_amt < an) {
		// bn + shift_amt > an
		assert(bn + shift_amt + 1 != 0);
		a.resize(bn + shift_amt + 1);
		size_t add_len1 = an - shift_amt;
		size_t add_len2 = bn - add_len1;	// add_len2 >= 1
		uint32_t carry = add_n(&a[shift_amt], &a[shift_amt], &b[0], add_len1);
		carry = copy_add_1(&a[an], &b[add_len1], carry, add_len2);
		if (carry == 0)
			a.resize(bn + shift_amt);
		else
			a[bn + shift_amt] = carry;
	}
	else {
		// shift_amt >= an
		a.resize(bn + shift_amt);
		set(&a[an], 0, shift_amt - an);
		copy(&a[shift_amt], &b[0], bn);
	}
}


/*
a = a - (b<<(shift_amt*32))
result must >= 0
*/
static void shift_sub(bigint& a, const bigint& b, size_t shift_amt)
{
	size_t an = a.size(), bn = b.size();
	if (bn == 0)
		return;

	assert(bn + shift_amt >= bn); // no overflow
	assert(bn + shift_amt <= an);
	
	uint32_t borrow = sub_short(&a[shift_amt], an - shift_amt, &b[0], bn);
	assert(borrow == 0); // result >= 0	
	a.resize(an - high_zeros(a));
}

static const size_t MUL_NAIVE_THRESHOLD = 64;

void mul_karatsuba_prototype(bigint& c, const bigint& a, const bigint& b)
{
	size_t an = a.size(), bn = b.size();
	if (an == 0 || bn == 0) {
		c.resize(0);
		return;
	}
	if (an <= MUL_NAIVE_THRESHOLD) {
		mul_naive(c, b, a);
		return;
	}
	if (bn <= MUL_NAIVE_THRESHOLD) {
		mul_naive(c, a, b);
		return;
	}

	size_t n = std::max(an, bn);
	size_t m = (n + 1) / 2;
	bigint a0, a1, b0, b1;
	bigint a0b0, a1b1;
	
	cut(a, 0, m, a0);	
	cut(b, 0, m, b0);		
	mul_karatsuba_prototype(a0b0, a0, b0);
	
	cut(a, m, n - m, a1);
	cut(b, m, n - m, b1);
	mul_karatsuba_prototype(a1b1, a1, b1);
	
	c = a0b0;
	shift_acc(c, a1b1, 2 * m);
	shift_acc(c, a0b0, m);
	shift_acc(c, a1b1, m);
	
	bigint a1_minus_a0, b1_minus_b0;
	int a1_minus_a0_sign = sub_sign_mag(a1_minus_a0, a1, a0);
	int b1_minus_b0_sign = sub_sign_mag(b1_minus_b0, b1, b0);
	bigint diff_prod;
	mul_karatsuba_prototype(diff_prod, a1_minus_a0, b1_minus_b0);
	if (a1_minus_a0_sign + b1_minus_b0_sign == 0) {
		shift_acc(c, diff_prod, m);
	}
	else {
		shift_sub(c, diff_prod, m);
	}
}