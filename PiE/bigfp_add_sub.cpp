#include <assert.h>
#include <algorithm>
#include "bigint.h"
#include "bigfp.h"
#include "int128.h"
#include "array_basic_arith.h"

static size_t trim_get_hi_zeros(std::vector<uint32_t>& data)
{
	size_t hi_z = high_zeros(data);
	data.resize(data.size() - hi_z);

	size_t lo_z = low_zeros(data);
	data.erase(data.begin(), data.begin() + lo_z);

	return hi_z;
}

/*
ignore sign
a.exp >= b.exp
*/
static void add_abs(bigfp& c, const bigfp_src& a, const bigfp_src& b)
{
	/*
		aaaaaaaaaaaaaa
			bbbb

		aaaaaa		==>		aaaaaa00
			bbbb			    bbbb

		aaa			==>		aaa00000
			bbbb			    bbbb

	*/
	assert(a.exp >= b.exp);
	int64_t e1 = a.exp, e2 = b.exp;
	size_t n1 = a.data.size(), n2 = b.data.size();

	int128 indent_n = (int128)e1 - e2;
	int128 zeropad_n = indent_n + n2 - n1;
	size_t cn;
	c.data.clear();
	if (zeropad_n < 0) {
		int128 cnx = (int128)n1 + 1;
		assert(cnx <= SIZE_MAX);

		cn = (size_t)cnx;
		size_t add_n = (size_t)(indent_n + n2);

		c.data.resize(cn);
		copy(&c.data[0], &a.data[0], n1);

		c.data[cn - 1] = acc_short(&c.data[n1 - add_n], add_n, &b.data[0], n2);
		c.exp = a.exp + 1;
		// overflow is ok, if c.data[cn - 1]==0, c.exp will be decreased 
	}
	else if (zeropad_n < n2) {
		int128 cnx = n1 + zeropad_n + 1;
		assert(cnx <= SIZE_MAX);

		cn = (size_t)cnx;
		size_t zn = (size_t)zeropad_n;

		c.data.resize(cn);
		set(&c.data[0], 0, zn);
		copy(&c.data[zn], &a.data[0], n1);

		c.data[cn - 1] = acc_short(&c.data[0], cn - 1, &b.data[0], n2);
		c.exp = a.exp + 1;
		// overflow is ok
	}
	else {
		// zeropad_n >= n2
		int128 cnx = n1 + zeropad_n;
		assert(cnx <= SIZE_MAX);

		cn = (size_t)cnx;
		size_t zn = (size_t)zeropad_n;

		c.data.resize(cn);
		copy(&c.data[0], &b.data[0], n2);
		set(&c.data[n2], 0, zn - n2);
		copy(&c.data[zn], &a.data[0], n1);
		c.exp = a.exp;
	}
	size_t hi_z = trim_get_hi_zeros(c.data);
	assert(hi_z <= 1);
	c.exp -= hi_z;
	assert(c.exp >= a.exp);
}

/*
ignore sign
a.exp > b.exp
*/
static void sub_abs(bigfp& c, const bigfp_src& a, const bigfp_src& b)
{
	/*
		aaaaaaaaaaaaaa
			bbbb

		aaaaaa		==>		aaaaaa00
			bbbb			    bbbb
	*/
	assert(a.exp > b.exp);
	int64_t e1 = a.exp, e2 = b.exp;
	size_t n1 = a.data.size(), n2 = b.data.size();


	int128 indent_n = (int128)e1 - e2;

	int128 zeropad_n = indent_n + n2 - n1;
	size_t cn;
	c.data.clear();
	if (zeropad_n < 0) {
		cn = n1;
		size_t sub_n = (size_t)(indent_n + n2);

		c.data.resize(cn);
		copy(&c.data[0], &a.data[0], n1);

		sub_short(&c.data[n1 - sub_n], sub_n, &b.data[0], n2);
	}
	else {
		int128 cnx = n1 + zeropad_n;
		assert(cnx <= SIZE_MAX);

		cn = (size_t)cnx;
		size_t zn = (size_t)zeropad_n;

		c.data.resize(cn);
		set(&c.data[0], 0, zn);
		copy(&c.data[zn], &a.data[0], n1);

		sub_short(&c.data[0], cn, &b.data[0], n2);
	}
	size_t hi_z = trim_get_hi_zeros(c.data);
	assert(hi_z < cn);
	int128 cexpx = (int128)a.exp - hi_z;
	assert(cexpx >= INT64_MIN);
	c.exp = (int64_t)cexpx;
}

/*
ignore sign
a.exp = b.exp
size(a.data) >= size(b.data)
return: sign of result
*/
static int sub_abs_exp_eq(bigfp& c, const bigfp_src& a, const bigfp_src& b)
{
	/*
		aaaaaaaaaaaaaaa
		bbbbbbbb
			^
			diff
	*/

	int64_t e1 = a.exp, e2 = b.exp;
	size_t n1 = a.data.size(), n2 = b.data.size();
	assert(e1 == e2);
	assert(n1 >= n2);

	int cmp_r;
	size_t diff_i;
	cmp_diff_index_max(cmp_r, diff_i, &a.data[n1 - n2], &b.data[0], n2);
	size_t cn;
	int128 cexpx;
	int sign;
	c.data.clear();
	if (cmp_r == 0) {
		if (n1 == n2) {
			return 0;
		}
		cn = n1 - n2;
		c.data.resize(cn);
		copy(&c.data[0], &a.data[0], cn);
		cexpx = (int128)e1 - n2;
		sign = 1;
	}
	else if (cmp_r > 0) {
		size_t sub_n = diff_i + 1;
		cn = sub_n + n1 - n2;
		c.data.resize(cn);
		copy(&c.data[0], &a.data[0], cn);
		sub_short(&c.data[n1 - n2], sub_n, &b.data[0], sub_n);
		cexpx = (int128)e1 - (n2 - sub_n);
		sign = 1;
	}
	else {
		// cmp_r < 0
		size_t zn = n1 - n2;
		cn = diff_i + 1 + zn;
		c.data.resize(cn);
		set(&c.data[0], 0, zn);
		copy(&c.data[zn], &b.data[0], diff_i + 1);
		sub_short(&c.data[0], cn, &a.data[0], cn);
		cexpx = (int128)e1 - (n2 - diff_i - 1);
		sign = -1;
	}
	size_t hi_z = trim_get_hi_zeros(c.data);
	assert(hi_z < cn);
	cexpx = cexpx - hi_z;
	assert(cexpx >= INT64_MIN);
	c.exp = (int64_t)cexpx;
	return sign;
}

/*
using b_sign, not b.sign
*/
void add(bigfp& c, const bigfp_src& a, const bigfp_src& b)
{
	assert(!is_from(a, c) && !is_from(b, c));
	if (a.sign == 0) {
		c = b;
		return;
	}
	if (b.sign == 0) {
		c = a;
		return;
	}

	if (a.sign == b.sign) {
		const bigfp_src& first = a.exp >= b.exp ? a : b;
		const bigfp_src& second = a.exp >= b.exp ? b : a;
		add_abs(c, first, second);
		c.sign = a.sign;
	}
	else if (a.exp != b.exp) {
		// a.sign != b.sign
		const bigfp_src& first = a.exp > b.exp ? a : b;
		const bigfp_src& second = a.exp > b.exp ? b : a;
		sub_abs(c, first, second);
		c.sign = a.exp > b.exp ? a.sign : b.sign;
	}
	else {
		// a.sign != b.sign, a.exp == b.exp
		bool a_longer = a.data.size() >= b.data.size();
		const bigfp_src& first = a_longer ? a : b;
		const bigfp_src& second = a_longer ? b : a;
		int sign = sub_abs_exp_eq(c, first, second);
		c.sign = sign * first.sign;
	}
}

void sub(bigfp& c, const bigfp_src& a, const bigfp_src& b)
{
	bigfp_src neg_b = b;
	neg_b.sign *= -1;
	add(c, a, neg_b);
}

int compare(const bigfp_src& a, const bigfp_src& b)
{
	// only deal with a.sign >= 0
	if (a.sign < 0) {
		bigfp_src neg_a = a, neg_b = b;
		neg_a.sign *= -1;
		neg_b.sign *= -1;
		return -compare(neg_a, neg_b);
		// WRONG: return compare(neg_b, neg_a);
		// guess why?
	}

	if (a.sign != b.sign)
		return a.sign > b.sign ? 1 : -1;
	if (a.sign == 0)
		return 0;
	// a.sign == b.sign == 1
	if (a.exp != b.exp)
		return a.exp > b.exp ? 1 : -1;

	size_t an = a.data.size(), bn = b.data.size();
	int cmp;
	size_t dummy;
	size_t cmp_len = std::min(an, bn);
	cmp_diff_index_max(cmp, dummy, &a.data[0], &b.data[0], cmp_len);
	if (cmp != 0)
		return cmp;
	if (an != bn)
		return an > bn ? 1 : -1;
	return 0;
}

bool is_less_or_equal(const bigfp_src& a, const bigfp_src& b)
{
	return compare(a, b) <= 0;
}


bool is_equal(const bigfp_src& a, const bigfp_src& b)
{
	if (a.sign == 0)
		return a.sign == b.sign;
	return a.sign == b.sign && a.exp == b.exp
		&& std::equal(a.data.begin(), a.data.end(),
					  b.data.begin(), b.data.end());
}