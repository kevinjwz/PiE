#pragma once
#include <stdint.h>

uint32_t add_n(uint32_t c[], const uint32_t a[], const uint32_t b[], size_t n);

uint32_t sub_n(uint32_t c[], const uint32_t a[], const uint32_t b[], size_t n);

uint32_t add_1(uint32_t a[], uint32_t b, size_t n);
uint32_t copy_add_1(uint32_t c[], const uint32_t a[], uint32_t b, size_t n);

uint32_t sub_1(uint32_t a[], uint32_t b, size_t n);
uint32_t copy_sub_1(uint32_t c[], const uint32_t a[], uint32_t b, size_t n);

uint32_t acc_short(uint32_t a[], size_t n, const uint32_t b[], size_t m);
uint32_t add_short(uint32_t c[], const uint32_t a[], size_t an, const uint32_t b[], size_t bn);

uint32_t sub_short(uint32_t a[], size_t n, const uint32_t b[], size_t m);
uint32_t sub_short(uint32_t c[], const uint32_t a[], size_t an, const uint32_t b[], size_t bn);

void cmp_diff_index_max(int& cmp, size_t& diffidx, const uint32_t a[], const uint32_t b[], size_t n);
int sub_sign_mag(uint32_t c[], const uint32_t a[], const uint32_t b[], size_t n);

void mul_naive_ub(uint32_t c[], const uint32_t a[], size_t an, const uint32_t b[], size_t bn);