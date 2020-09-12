#include <algorithm>
#include "bigint.h"
#include "mul_karatsuba.h"
#include "mul_ssa.h"
#include "mul_fft_dwt.h"


void mul_auto(bigint& c, const bigint_src& a, const bigint_src& b)
{
	size_t an = a.size(), bn = b.size();
	size_t shorter_n = std::min(an, bn);

	if (shorter_n <= 16) {
		mul_naive(c, a, b);
		return;
	}

	if (an + bn <= 800 || shorter_n <= 225) {
		mul_karatsuba_unbalance(c, a, b);
		return;
	}

	// when size is large, floating-point dwt eats too much memory!
	if (an + bn <= 1 << 16) {
		mul_dwt(c, a, b);
		return;
	}
	mul_ssa_fast(c, a, b);
}