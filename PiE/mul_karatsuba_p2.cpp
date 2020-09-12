//#define NDEBUG
#include <stdint.h>
#include <assert.h>
#include "bigint.h"
#include "array_basic_arith.h"
#include "array_helpers.h"
#include "micro_helpers.inl"

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

	uint32_t carry = acc_short(a + shift_amt, an - shift_amt, b, bn);
	assert(carry == 0); // no overflow
}

/*
a = a - (b<<(shift_amt*32))
an must >= bn + shift_amt
when b overlays with a,  b >= a + shift_amt
*/

static void shift_sub(uint32_t a[], size_t an, const uint32_t b[], size_t bn, size_t shift_amt)
{
	assert(an >= bn + shift_amt);
	uint32_t borrow = sub_short(a + shift_amt, an - shift_amt, b, bn);
	assert(borrow == 0); // no underflow
}

static const size_t MUL_NAIVE_THRESHOLD = 16;

/*
size(a) = n
size(b) = n
size(prod) = 2n
n must be power of 2
a will be destroyed
*/

static void mul_krtb_rec_p2(uint32_t prod[], uint32_t a[], const uint32_t b[], uint32_t space[], size_t n)
{
	if (n <= MUL_NAIVE_THRESHOLD) {
		mul_naive_ub(prod, a, n, b, n);
		return;
	}

	size_t n0 = n / 2;
	//size_t n1 = n - n0;
	uint32_t* a0 = a;
	uint32_t* a1 = a + n0;
	const uint32_t* b0 = b;
	const uint32_t* b1 = b + n0;

	uint32_t* a0_minus_a1 = space;		// size: n0
	uint32_t* a0b0 = prod;			// size: 2 * n0;
	uint32_t* a1b1 = prod + 2 * n0;	// size: 2 * n0;

	int a0_minus_a1_sign = sub_sign_mag(a0_minus_a1, a0, a1, n0);
	mul_krtb_rec_p2(a0b0, a0, b0, space + n0, n0);
	mul_krtb_rec_p2(a1b1, a1, b1, space + n0, n0);

	uint32_t* a0b0_copy = space + n0;	// size: 2*n0
	copy(a0b0_copy, a0b0, 2 * n0);

	shift_add(prod, 2 * n, a1b1, 2 * n0, n0);
	shift_add(prod, 2 * n, a0b0_copy, 2 * n0, n0);

	uint32_t* b0_minus_b1 = space + n0;		// size: n0
	int b0_minus_b1_sign = sub_sign_mag(b0_minus_b1, b0, b1, n0);
	uint32_t* diff_prod = a;	// size: 2 * n0
	mul_krtb_rec_p2(diff_prod, a0_minus_a1, b0_minus_b1, space + 2 * n0, n0);
	if (a0_minus_a1_sign + b0_minus_b1_sign == 0) {
		shift_add(prod, 2 * n, diff_prod, 2 * n0, n0);
	}
	else {
		// before subtraction, prod < 2^(2*n) is guaranteed
		shift_sub(prod, 2 * n, diff_prod, 2 * n0, n0);
	}
}


/*
size(a) = size(b) = n
size(prod) = 2n
faster version for SSA basecase
*/
void mul_karatsuba_p2_destroy_a(uint32_t prod[], uint32_t a[], const uint32_t b[], size_t n)
{
	assert(is_pow2(n));
	if (n <= MUL_NAIVE_THRESHOLD) {
		mul_naive_ub(prod, a, n, b, n);
		return;
	}
	uint32_t* _space = new uint32_t[2 * n];

	mul_krtb_rec_p2(prod, a, b, _space, n);
	delete[] _space;
}

