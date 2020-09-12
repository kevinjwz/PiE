#pragma once

#include<complex>
typedef std::complex<double> dcomplex;

//template < bool inv >
//void fft_mm_radix_2(dcomplex result[], const dcomplex x[], size_t n);
//
//template < bool inv >
//void fft_mm_radix_4(dcomplex result[], const dcomplex x[], size_t n);

template < bool inv >
void fft_mm_radix_2(dcomplex x[], size_t n);

template < bool inv >
void fft_mm_radix_4(dcomplex x[], size_t n);