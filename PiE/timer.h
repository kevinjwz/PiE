#pragma once
#include<chrono>
#include<stdio.h>

template<typename FuncT>
inline double timer_result(const char* format, FuncT func)
{
	using clock = std::chrono::high_resolution_clock;
	using double_sec = std::chrono::duration<double>;

	auto begin = clock::now();
	auto result = func();
	auto end = clock::now();
	double sec = double_sec(end - begin).count();
	printf(format, sec, result);
	return sec;
}

template<typename FuncT>
inline double timer(const char* format, FuncT func)
{
	using clock = std::chrono::high_resolution_clock;
	using double_sec = std::chrono::duration<double>;

	auto begin = clock::now();
	func();
	auto end = clock::now();
	double sec = double_sec(end - begin).count();
	printf(format, sec);
	return sec;
}