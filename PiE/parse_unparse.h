#pragma once
#include <vector>
#include <ostream>
#include "bigint.h"
#include "bigfp.h"

void parse_naive(bigint& a, const std::string& str);
void unparse_naive(std::string& str, const bigint& a);

void parse_fast(bigfp& a, const char str[], size_t n);
void unparse_fast(char str[], size_t n, const bigfp_src& a);
size_t estimate_unparse_length(size_t word_n);
size_t estimate_unparse_length(const bigfp_src& src);

int unparse_fraction(char str[], size_t n, const bigfp_src& frac);
void unparse_integer_fraction(std::vector<char>& int_str, std::vector<char>& frac_str, const bigfp_src& src, size_t frac_digit_n);

void print_number_formatted(std::ostream& file, const std::vector<char>& int_str, const std::vector<char>& frac_str, size_t reserved_digit_n);