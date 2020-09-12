#include<stdint.h>
#include<assert.h>

#include "bigint.h"
#include "array_basic_arith.h"
#include "array_helpers.h"


void add(bigint& c, const bigint& a, const bigint& b)
{
	assert(!is_alias(a, c) && !is_alias(b, c));
	c.clear();
	size_t an = a.size();
	size_t bn = b.size();

	const bigint& l = an >= bn ? a : b;
	const bigint& s = an >= bn ? b : a;
	size_t ln = an >= bn ? an : bn;
	size_t sn = an >= bn ? bn : an;

	if (sn == 0) {
		c = l;
		return;
	}
	// ln >= sn > 0
	c.resize(ln);
	uint32_t carry = add_short(&c[0], &l[0], ln, &s[0], sn);
	if (carry > 0)
		c.push_back(carry);
}

/*
	b is a ... ok
*/
void add(bigint& a, const bigint& b)
{
	if (b.size() == 0)
		return;
	size_t an = a.size();
	size_t bn = b.size();

	uint32_t carry;
	if (an >= bn) {
		// ok if b is a (just as add_n(a, a, a, an))
		carry = acc_short(&a[0], an, &b[0], bn);
	}
	else {	// an < bn
		a.resize(bn);
		carry = add_n(&a[0], &a[0], &b[0], an);
		carry = copy_add_1(&a[an], &b[an], carry, bn - an);
	}
	if (carry > 0)
		a.push_back(carry);
}

/*
	c = |a-b|
	return: sign(a-b)
*/
int sub_sign_mag(bigint& c, const bigint& a, const bigint& b)
{
	assert(!is_alias(a, c) && !is_alias(b, c));
	c.clear();
	size_t an = a.size();
	size_t bn = b.size();
	if (bn == 0) {
		c = a;
		return an == 0 ? 0 : 1;
	}
	if (an == 0) {
		// bn>0
		c = b;
		return -1;
	}

	if (an == bn) {
		int cmp;
		size_t idxmax;
		cmp_diff_index_max(cmp, idxmax, &a[0], &b[0], an);
		if (cmp == 0) {
			return 0;
		}

		const bigint& max = cmp > 0 ? a : b;
		const bigint& min = cmp > 0 ? b : a;
		size_t sub_len = idxmax + 1;
		c.resize(sub_len);
		uint32_t borrow = sub_n(&c[0], &max[0], &min[0], sub_len);
		assert(borrow == 0);
		c.resize(sub_len - high_zeros(c));
		return cmp;
	}
	else {
		const bigint& l = an > bn ? a : b;
		const bigint& s = an > bn ? b : a;
		size_t ln = an > bn ? an : bn;
		size_t sn = an > bn ? bn : an;
		c.resize(ln);
		uint32_t borrow = sub_short(&c[0], &l[0], ln, &s[0], sn);
		assert(borrow == 0);
		c.resize(ln - high_zeros(c));
		return an > bn ? 1 : -1;
	}
}

/*
	a = |a-b|
	return: sign(a-b)
	b is a ... ok
*/
int sub_sign_mag(bigint& a, const bigint& b)
{
	size_t an = a.size();
	size_t bn = b.size();
	if (bn == 0) {
		return an == 0 ? 0 : 1;
	}
	if (an == 0) {
		// bn>0
		a = b;
		return -1;
	}

	if (an == bn) {
		int cmp;
		size_t idxmax;
		cmp_diff_index_max(cmp, idxmax, &a[0], &b[0], an);
		if (cmp == 0) {
			a.clear();
			return 0;
		}
		// is_alias(a, b) == false
		const bigint& max = cmp > 0 ? a : b;
		const bigint& min = cmp > 0 ? b : a;
		size_t sub_len = idxmax + 1;	// sub_len<=an
		uint32_t borrow = sub_n(&a[0], &max[0], &min[0], sub_len);
		assert(borrow == 0);
		a.resize(sub_len - high_zeros(&a[0], sub_len));
		return cmp;
	}
	else {
		if (an > bn) {
			uint32_t borrow = sub_short(&a[0], an, &b[0], bn);
			assert(borrow == 0);
			a.resize(an - high_zeros(a));
			return 1;
		}
		else {	// bn>an
			a.resize(bn);
			uint32_t borrow = sub_n(&a[0], &b[0], &a[0], an);
			borrow = copy_sub_1(&a[an], &b[an], borrow, bn - an);
			assert(borrow == 0);
			a.resize(bn - high_zeros(a));
			return -1;
		}

	}
}

/*
	a += b
*/
void add_1(bigint& a, uint32_t b)
{
	if (b == 0)
		return;
	if (a.size() == 0) {
		a = { b };
		return;
	}
	uint32_t carry = add_1(&a[0], b, a.size());
	if (carry > 0)
		a.push_back(carry);
}

/*
	a = |a-b|
	return: sign(a-b)
*/
int sub_1_sign_mag(bigint& a, uint32_t b)
{
	if (b == 0)
		return a.size() == 0 ? 0 : 1;
	if (a.size() == 0) {
		a = { b };
		return -1;
	}
	uint32_t borrow = sub_1(&a[0], b, a.size());
	if (borrow > 0) {
		// a.size() == 1
		// a[0] == 2^32 + a0 - b
		a[0] = -a[0];
		return -1;
	}
	if (*a.rbegin() == 0)
		a.pop_back();
	return a.size() == 0 ? 0 : 1;
}

int compare(const bigint& a, const bigint& b)
{
	size_t an = a.size(), bn = b.size();
	if (an < bn)
		return -1;
	if (an > bn)
		return 1;
	if (an == 0)
		return 0;
	int cmp;
	size_t dummy;
	cmp_diff_index_max(cmp, dummy, &a[0], &b[0], an);
	return cmp;
}

/*
result c is trimmed
*/
void mul_naive(bigint& c, const bigint_src& a, const bigint_src& b)
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
	// better performance if ln>=sn
	// longer inner loop, shorter outer loop
	mul_naive_ub(&c[0], &l[0], ln, &s[0], sn);
	if (c[cn - 1] == 0)
		c.pop_back();
}


/*
dividend and quotient can be the same bigint
*/
void div_rem_1_naive(bigint& quotient, uint32_t& remainder, const bigint& dividend, uint32_t divisor)
{
	assert(divisor != 0);
	size_t dn = dividend.size();
	quotient.clear();
	quotient.resize(dn);

	uint64_t r = 0;
	// dn-1 downto 0 (poor unsigned size_t)
	for (size_t i = dn - 1; i < dn; --i) {
		uint64_t littled = (r << 32) + dividend[i];
		uint64_t littleq = littled / divisor;  // guaranteed: littleq < 2^32
		r = littled % divisor;
		quotient[i] = (uint32_t)littleq;
	}
	remainder = (uint32_t)r;
	quotient.resize(dn - high_zeros(quotient));
}