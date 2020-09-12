#pragma once
#include <stdint.h>
#include <vector>

void copy(uint32_t dst[], const uint32_t src[], size_t n);

void set(uint32_t dst[], uint32_t val, size_t n);

size_t count_val_low(const uint32_t arr[], size_t n, uint32_t val);

size_t count_val_high(const uint32_t arr[], size_t n, uint32_t val);

size_t low_zeros(const uint32_t arr[], size_t n);

size_t high_zeros(const uint32_t arr[], size_t n);

size_t low_zeros(const std::vector<uint32_t>& v);

size_t high_zeros(const std::vector<uint32_t>& v);

void copy_shift_left(uint32_t dst[], const uint32_t src[], size_t n, size_t shift_bits);

void copy_shift_right(uint32_t dst[], const uint32_t src[], size_t n, size_t shift_bits);

void shift_left(uint32_t arr[], size_t n, size_t shift_bits);

void shift_right(uint32_t arr[], size_t n, size_t shift_bits);