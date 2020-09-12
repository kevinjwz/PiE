#pragma once
#include <iterator>
#include <vector>
#include <functional>

template < typename T>
struct readonly_span
{
private:
	const T* begin_ptr;
	size_t len;
public:
	const T& operator[](size_t index) const { return begin_ptr[index]; }
	size_t size() const { return len; }

	const T* begin() const { return begin_ptr; }
	const T* end() const { return begin_ptr + len; }

	auto rbegin() const { return std::make_reverse_iterator(this->end()); }
	auto rend() const { return std::make_reverse_iterator(this->begin()); }

	readonly_span() :begin_ptr(nullptr), len(0) {}

	readonly_span(const T* begin_ptr, size_t len) :
		begin_ptr(len == 0 ? nullptr : begin_ptr),
		len(len)
	{
	}

	readonly_span(const std::vector<T>& vec) :
		begin_ptr(vec.empty() ? nullptr : &vec[0]),
		len(vec.size())
	{
	}
};

/*
detect if a references b
*/
template < typename T>
inline bool is_from(const readonly_span<T>& a, const std::vector<T>& b)
{
	if (a.size() == 0)
		return false;

	std::less<const uint32_t*> less;

	const uint32_t
		*b_begin = b.data(),
		*b_end = b.data() + b.size();

	// a.begin >= b.begin && a.end <= b.end
	return !less(a.begin(), b_begin) && !less(b_end, a.end());
}