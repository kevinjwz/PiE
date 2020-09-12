#include <ostream>
#include <vector>
#include <assert.h>

using namespace std;

static void print_1000_digits(ostream& file, const char str[])
{
	file << '\n';
	const char* ptr = str;
	for (int row = 0; row < 20; ++row) {
		for (int col = 0; col < 5; ++col) {
			file << (col == 0 ? '\n' : ' ');
			file.write(ptr, 10);
			ptr += 10;
		}
	}
}

static void print_less_than_1000(ostream& file, const char str[], size_t n)
{
	assert(n > 0 && n < 1000);
	file << '\n';
	const char* ptr = str;
	for (int row = 0; row < 20; ++row) {
		for (int col = 0; col < 5; ++col) {
			file << (col == 0 ? '\n' : ' ');
			if (n <= 10) {
				file.write(ptr, n);
				return;
			}
			file.write(ptr, 10);
			ptr += 10;
			n -= 10;
		}
	}
}

void print_number_formatted(ostream& file, const vector<char>& int_str, const vector<char>& frac_str, size_t reserved_digit_n)
{
	assert(frac_str.size() >= reserved_digit_n);
	if (int_str.empty()) {
		file << '0';
	}
	else {
		file.write(int_str.data(), int_str.size());
	}
	file << '.';

	size_t frac_prec_n = frac_str.size() - reserved_digit_n;
	const char* ptr = frac_str.data();
	while (frac_prec_n >= 1000) {
		print_1000_digits(file, ptr);
		frac_prec_n -= 1000;
		ptr += 1000;
	}
	if (frac_prec_n > 0) {
		print_less_than_1000(file, ptr, frac_prec_n);
		ptr += frac_prec_n;
	}
	if (reserved_digit_n>0) {
		file << '(';
		file.write(ptr, reserved_digit_n);
		file << ')';
	}
	file << '\n';
}