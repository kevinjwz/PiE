//#define NDEBUG
#include <stdint.h>
#include <algorithm>
#include <assert.h>
#include "bigint.h"
#include "micro_helpers.inl"
#include "array_basic_arith.h"
#include "array_helpers.h"

static const size_t MUL_NAIVE_THRESHOLD = 16;



/*
a = |a - b|
return: sign
a>b: sign = 1
a=b: sign = 0
a<b: sign = -1

size(c) = size(a) = an
size(b) = bn
an must >= bn
*/

static int sub_sign_mag_ub(uint32_t c[], const uint32_t a[], size_t an, const uint32_t b[], size_t bn)
{
	size_t a_ext_n = an - bn;
	if (a_ext_n > high_zeros(a + bn, a_ext_n)) {
		uint32_t borrow = sub_short(c, a, an, b, bn);
		assert(borrow == 0);
		return 1;
	}
	else {
		int sign = sub_sign_mag(c, a, b, bn);
		set(c + bn, 0, a_ext_n);
		return sign;
	}
}

/*
a = a + (b<<(shift_amt*32))
an must >= bn + shift_amt
when b overlays with a,  b >= a + shift_amt
*/

static void shift_add(uint32_t a[], size_t an, const uint32_t b[], size_t bn, size_t shift_amt)
{
	//assert(bn + shift_amt >= bn);
	assert(an >= bn + shift_amt);
	//if (b < a + an || a < b + bn)
	//	assert(b >= a + shift_amt);
	[[maybe_unused]]
	uint32_t carry = acc_short(a + shift_amt, an - shift_amt, b, bn);
	// here, carry = 1 (overflow) is allowed
}

/*
a = a - (b<<(shift_amt*32))
an must >= bn + shift_amt
when b overlays with a,  b >= a + shift_amt
*/

static void shift_sub(uint32_t a[], size_t an, const uint32_t b[], size_t bn, size_t shift_amt)
{
	assert(an >= bn + shift_amt);
	[[maybe_unused]]
	uint32_t borrow = sub_short(a + shift_amt, an - shift_amt, b, bn);
	// here, borrow = 1 (underflow) is allowed
}


/*
size(a) = n
size(b) = n
size(prod) = 2n
n must > 0
n+1 elements in a will be destroyed !!!
*/

static void mul_krtb_rec(uint32_t prod[], uint32_t a[], const uint32_t b[], uint32_t space[], size_t n)
{
	/*
	when n is odd, ceil(n/2) must >= 3, so n must >= 4
	if n = 1 or 3, the intermediate value of prod may become too big
	*/
	if (n <= MUL_NAIVE_THRESHOLD) {
		mul_naive_ub(prod, a, n, b, n);
		return;
	}

	size_t n0 = (n + 1) / 2;
	size_t n1 = n - n0;
	uint32_t* a0 = a;
	uint32_t* a1 = a + n0;
	const uint32_t* b0 = b;
	const uint32_t* b1 = b + n0;

	uint32_t* a0_minus_a1 = space;		// size: n0
	uint32_t* a0b0 = prod;			// size: 2 * n0;
	uint32_t* a1b1 = prod + 2 * n0;	// size: 2 * n1;

	int a0_minus_a1_sign = sub_sign_mag_ub(a0_minus_a1, a0, n0, a1, n1);

	mul_krtb_rec(a1b1, a1, b1, space + n0, n1);
	mul_krtb_rec(a0b0, a0, b0, space + n0, n0);

	uint32_t* a0b0_copy = space + n0;	// size: 2*n0
	copy(a0b0_copy, a0b0, 2 * n0);

	shift_add(prod, 2 * n, a1b1, 2 * n1, n0);
	shift_add(prod, 2 * n, a0b0_copy, 2 * n0, n0);

	uint32_t* b0_minus_b1 = space + n0 + 1;		// size: n0, +1 for overwriting
	int b0_minus_b1_sign = sub_sign_mag_ub(b0_minus_b1, b0, n0, b1, n1);
	uint32_t* diff_prod = a;	// size: 2 * n0
	mul_krtb_rec(diff_prod, a0_minus_a1, b0_minus_b1, space + 2 * n0 + 1, n0);
	if (a0_minus_a1_sign + b0_minus_b1_sign == 0) {
		shift_add(prod, 2 * n, diff_prod, 2 * n0, n0);
	}
	else {
		shift_sub(prod, 2 * n, diff_prod, 2 * n0, n0);
	}
}

/*
size(a) = size(b) = n
size(c) = 2n
for SSA computing "ci mod k" 
*/
void mul_karatsuba(uint32_t c[], const uint32_t a[], const uint32_t b[], size_t n)
{
	uint32_t* _a = new uint32_t[n + 1];
	uint32_t* _space = new uint32_t[2 * n + 3 * ceil_log2(n)];

	copy(_a, &a[0], n);

	mul_krtb_rec(&c[0], _a, b, _space, n);
	delete[] _a;
	delete[] _space;
}


/*
size(a) = an
size(b) = bn
size(prod) = an+bn

an must >= bn

an+1 elements in a will be destroyed !!!

*/

static void mul_krtb_rec_ub(uint32_t prod[], uint32_t a[], size_t an, const uint32_t b[], size_t bn, uint32_t space[])
{
	/*
	when n is odd, ceil(n/2) must >= 3, so n must >= 4
	if n = 1 or 3, the intermediate value of prod may become too big
	*/
	if (bn <= MUL_NAIVE_THRESHOLD) {
		mul_naive_ub(prod, a, an, b, bn);
		return;
	}

	if (an == bn) {
		mul_krtb_rec(prod, a, b, space, an);
		return;
	}

	size_t n0 = (an + 1) / 2;
	size_t n1 = an - n0;
	uint32_t* a0 = a;
	uint32_t* a1 = a + n0;

	if (bn <= n1) {
		uint32_t* a0b = space;	// size: n0+bn
		uint32_t* a1b = prod + n0;  // size: n1+bn

		mul_krtb_rec_ub(a1b, a1, n1, b, bn, space);
		mul_krtb_rec_ub(a0b, a0, n0, b, bn, space + n0 + bn);

		// copy low n0 limbs of a0b, then add high bn limbs of a0b
		copy(prod, a0b, n0);
		shift_add(prod, an + bn, a0b + n0, bn, n0);
		return;
	}
	if (bn == n0) {
		// n1 = n0-1 && bn == n0

		uint32_t* a0b = prod;	// size: n0+bn
		uint32_t* a1b = space;	// size: n0+bn

		uint32_t a1_0_save = a1[0];
		mul_krtb_rec_ub(a0b, a0, n0, b, bn, space);

		// size(a1_copy)=n0, because an must >= bn
		uint32_t* a1_copy = a0;	// size: n0
		a1_copy[0] = a1_0_save;
		copy(a1_copy + 1, a1 + 1, n1 - 1);
		a1_copy[n1] = 0;

		mul_krtb_rec_ub(a1b, a1_copy, n0, b, bn, space + n0 + bn);

		// copy high n1 limbs of a1b, then add low bn limbs of a1b
		copy(prod + n0 + bn, a1b + bn, n1);
		shift_add(prod, an + bn, a1b, bn, n0);
		return;
	}
	// bn>n0
	size_t bn1 = bn - n0;

	const uint32_t* b0 = b;
	const uint32_t* b1 = b + n0;

	uint32_t* a0_minus_a1 = space;		// size: n0
	uint32_t* a0b0 = prod;			// size: 2 * n0;
	uint32_t* a1b1 = prod + 2 * n0;	// size: n1 + bn1;

	int a0_minus_a1_sign = sub_sign_mag_ub(a0_minus_a1, a0, n0, a1, n1);

	mul_krtb_rec_ub(a1b1, a1, n1, b1, bn1, space + n0);
	mul_krtb_rec(a0b0, a0, b0, space + n0, n0);

	uint32_t* a0b0_copy = space + n0;	// size: 2*n0
	copy(a0b0_copy, a0b0, 2 * n0);

	shift_add(prod, an + bn, a1b1, n1 + bn1, n0);
	shift_add(prod, an + bn, a0b0_copy, 2 * n0, n0);	// an+bn >= 3*n0

	uint32_t* b0_minus_b1 = space + n0 + 1;		// size: n0, +1 for overwriting
	int b0_minus_b1_sign = sub_sign_mag_ub(b0_minus_b1, b0, n0, b1, bn1);
	uint32_t* diff_prod = a;	// size: 2 * n0
	mul_krtb_rec(diff_prod, a0_minus_a1, b0_minus_b1, space + 2 * n0 + 1, n0);
	if (a0_minus_a1_sign + b0_minus_b1_sign == 0) {
		shift_add(prod, an + bn, diff_prod, 2 * n0, n0);
	}
	else {
		shift_sub(prod, an + bn, diff_prod, 2 * n0, n0);
	}
}

/*
size(prod) = an+bn
an must >= bn
*/
static void mul_krtb_ub_raw(uint32_t c[], const uint32_t a[], size_t an, const uint32_t b[], size_t bn)
{
	assert(an >= bn);
	uint32_t* a_copy = new uint32_t[an + 1];
	uint32_t* space = new uint32_t[2 * an + 3 * ceil_log2(an)];
	copy(a_copy, a, an);
	mul_krtb_rec_ub(c, a_copy, an, b, bn, space);
	delete[] a_copy;
	delete[] space;
}

void mul_karatsuba_unbalance(bigint& c, const bigint_src& a, const bigint_src& b)
{
	assert(!is_from(a, c) && !is_from(b, c));	
	size_t an = a.size(), bn = b.size();

	const bigint_src& l = an >= bn ? a : b;
	const bigint_src& s = an >= bn ? b : a;
	size_t ln = an >= bn ? an : bn;
	size_t sn = an >= bn ? bn : an;

	if (sn == 0) {
		c.clear();
		return;
	}

	size_t cn = ln + sn;
	assert(cn >= ln);	// assert no overflow
	c.clear();
	c.resize(cn);
	mul_krtb_ub_raw(&c[0], &l[0], ln, &s[0], sn);
	if (c[cn - 1] == 0)
		c.pop_back();
}