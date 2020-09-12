#pragma once
#include <stdint.h>
#include "bigint.h"

void mul_ssa_fast(bigint& c, const bigint_src& a, const bigint_src& b);

void mul_ssa_pure_int(bigint& c, const bigint_src& a, const bigint_src& b);