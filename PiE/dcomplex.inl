#pragma once

#include<complex>
typedef std::complex<double> dcomplex;

/*
x*c*i
*/
inline dcomplex mul_i(dcomplex x)
{
	return { -x.imag(), x.real() };
}

/*
im(x)+re(x)i
*/
inline dcomplex reflect(dcomplex x)
{
	return { x.imag(), x.real() };
}