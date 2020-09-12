#include<assert.h>
#include"bigint.h"
#include"bigfp.h"
#include"int128.h"
#include"micro_helpers.inl"

void mul(bigfp& c, const bigfp_src& a, const bigfp_src& b)
{
	if (a.sign == 0 || b.sign == 0) {
		c.sign = 0;
		c.data.clear();
		return;
	}

	c.sign = a.sign*b.sign;
	mul_auto(c.data, a.data, b.data);
	int128 cexpx = (int128)a.exp + b.exp;
	if (c.data.size() < a.data.size() + b.data.size()) {
		cexpx = cexpx - 1;
	}
	size_t lo_z = low_zeros(c.data);
	c.data.erase(c.data.begin(), c.data.begin() + lo_z);
	c.exp = checked_cast<int64_t>(cexpx);
}

void mul_pow2(bigfp& c, const bigfp_src& a, int64_t e)
{
	if (a.sign == 0 || e == 0) {
		c = a;
		return;
	}
	
	// e = E*32 + shift, shift in [0, 31]	
	int64_t E = floor_div(e, 32);	// floor(e/32)
	size_t shift = (size_t)floor_rem(e, 32);

	/*
		0. 00001000 * B^(E+1)
		=  00001000 * B^E	
	*/
	bigfp pow2;
	pow2.sign = 1;
	pow2.exp = E + 1;
	pow2.data = { 1u << shift };

	mul(c, a, pow2);
	// this method is fast enough
	// no need to use shift tricks (maybe)
}

void mul_pow2(bigfp& a, int64_t e)
{
	if (a.sign == 0 || e == 0)
		return;

	// e = E*32 + shift, shift in [0, 31]	
	int64_t E = floor_div(e, 32);	// floor(e/32)
	size_t shift = (size_t)floor_rem(e, 32);

	int128 exp = (int128)a.exp + E;

	if (shift > 0) {
		size_t n = a.data.size();
		size_t lo_0bits = trailing_zeros(a.data[0]);
		/*
			0. 0000aaaa aaaaaaaa aaaa0000
			or
			0. 000aaa00
		*/
		if (32 - lo_0bits > shift) {
			// after shift, the low word != 0
			uint32_t hi = a.data[n - 1];
			size_t hi_0bits = leading_zeros(hi);
			shift_left(&a.data[0], n, shift);
			if (hi_0bits < shift) {
				a.data.push_back(hi >> (32 - shift));
				exp = exp + 1;
			}
		}
		else {
			/* after shift, the low word == 0
				0. xxxx1234 56789abc defg0000 = A
			 xxxx. 12345678 9abcdefg 00000000 = L    after left-shift (imaginary)
				0. 0000xxxx 12345678 9abcdefg = R1   after right-shift
				         0. 12345678 9abcdefg = R2   after remove hi
				if xxxx = 0,  L = R2
				if xxxx != 0, L = R1*2^32
			*/
			shift_right(&a.data[0], n, 32 - shift);
			if (a.data[n - 1] == 0) {
				a.data.pop_back();
			}
			else {
				exp = exp + 1;
			}
		}
	}
	a.exp = checked_cast<int64_t>(exp);
}