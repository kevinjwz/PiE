#pragma once
#include<stdint.h>
#include "bigint.h"

void mul_fft_fp(bigint& c, const bigint& a, const bigint& b);

void mul_dwt(bigint& prod, const bigint_src& a, const bigint_src& b);

void mul_dwt_mod_pow2_plus_1(uint32_t prod[], uint32_t& prod_hi, size_t n,
							 const uint32_t a[], size_t an, const uint32_t b[], size_t bn);
