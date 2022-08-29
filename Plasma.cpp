#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <complex>
#include <vector>
#include <map>
// #include <fftw3.h>
#include "Plasma_funcs.h"
#include "FFT.h"

using namespace std;

int main(int argc, char* argv[])
{
    if(argc!=2)
	{
		cout<<"usage: $PATH\\./Plasma Pars.config"<<endl;
		exit(1);
	}

    double f[2*FLENGTH], S[2*FLENGTH], t[LENGTH], R[LENGTH];
    Plasma_pars P;

    Set_pars(argv[1], &P);

    Spectrum(f, S, &P);
	
	ACF(t, R, S);

    ofstream out;
	out.open("Spectrum.dat");
	for(size_t x=0; x<2*FLENGTH; x++)
	{
		out<<f[x]<<"\t";
		out<<S[x]<<"\n";
	}
	out.close();
	out.open("ACF.dat");
	for(size_t x=0; x<LENGTH; x++)
	{
		out<<t[x]<<"\t";
		out<<R[x]<<"\n";
	}
	out.close();
}
