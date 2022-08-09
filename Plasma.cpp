#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <complex>
#include <vector>
#include <map>
#include "Plasma_funcs.h"

using namespace std;

int main(int argc, char* argv[])
{
    if(argc!=2)
	{
		cout<<"usage: $PATH\\./Plasma Pars.config"<<endl;
		exit(1);
	}

    vector<double> f, S;
    Plasma_pars P;

    Set_pars(argv[1], P);

    Spectrum(f, S, P);

    ofstream out;
	out.open("Spectrum.dat");
	for(size_t x=0; x<f.size(); x++)
	{
		out<<f[x]<<"\t";
		out<<S[x]<<"\n";
	}
	out.close();
}
