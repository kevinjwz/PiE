#include<stdint.h>
#include<assert.h>
#include<intrin.h>

uint32_t add_sat(uint32_t a, uint32_t b)
{
	uint32_t s = a + b;
	return s >= a ? s : UINT32_MAX;
}

uint32_t add_4(uint32_t s[4], const uint32_t a[4], const uint32_t b[4], uint32_t c_in)
{
	assert(c_in <= 1);
	bool p[4], g[4];
	for (int i = 0; i < 4; ++i) {
		s[i] = a[i] + b[i];
		uint32_t sum_sat = add_sat(a[i], b[i]);
		p[i] = sum_sat == UINT32_MAX;
		g[i] = sum_sat != s[i];
	}

	bool c_add[4], c_out;
	c_add[0] = c_in;
	c_add[1] = c_add[0] & p[0] | g[0];
	c_add[2] = c_add[1] & p[1] | g[1];
	c_add[3] = c_add[2] & p[2] | g[2];
	c_out    = c_add[3] & p[3] | g[3];

	for (int i = 0; i < 4; ++i) {
		s[i] += c_add[i];
	}
	return c_out;
}