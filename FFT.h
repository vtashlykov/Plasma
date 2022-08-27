#ifndef FFT_H_
#define FFT_H_

struct Complex;
struct ShortComplex;

/*
    Fast Fourier Transformation: direct (complement= false)
    and complement (complement = true). 'x' is data source.
	'x' contains 2^T items.
	N: N - number of items in array
*/
extern void fft(ShortComplex *x, int T, bool complement);
extern void universal_fft(ShortComplex *x, int N, bool complement);

struct ShortComplex
{
    double re, im;
    inline void operator=(const Complex &y);
};

struct Complex
{
    long double re, im;
    inline void operator= (const Complex &y);
    inline void operator= (const ShortComplex &y);
};


inline void ShortComplex::operator=(const Complex &y)    { re = (double)y.re; im = (double)y.im;  }
inline void Complex::operator= (const Complex &y)       { re = y.re; im = y.im; }
inline void Complex::operator= (const ShortComplex &y)  { re = y.re; im = y.im; }

#endif


