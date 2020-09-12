#pragma once
#include <stdint.h>
#include "bigint.h"

void mul_karatsuba_prototype(bigint& c, const bigint& a, const bigint& b);

void mul_karatsuba(uint32_t c[], const uint32_t a[], const uint32_t b[], size_t n);

void mul_karatsuba_unbalance(bigint& c, const bigint_src& a, const bigint_src& b);

void mul_karatsuba_p2_destroy_a(uint32_t prod[], uint32_t a[], const uint32_t b[], size_t n);
