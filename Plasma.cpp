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

    ifstream inp;
    ofstream out;
    vector<double> S;
    Plasma_pars P;
    map<unsigned, double> C;
    string line;
    int found;

    inp.open(argv[1]);
	getline(inp, line);
	getline(inp, line);
	found=line.find_first_of(" = ");
	P.lambda=stod(line.substr(found+3, line.length()-found-3));
    printf("lambda %f\n", P.lambda);
	getline(inp, line);
    found=line.find_first_of(" = ");
	P.Ne=stod(line.substr(found+3, line.length()-found-3));
	printf("Ne %f\n", P.Ne);
    getline(inp, line);
    found=line.find_first_of(" = ");
	P.Te=stod(line.substr(found+3, line.length()-found-3));
	printf("Te %f\n", P.Te);
    getline(inp, line);
    found=line.find_first_of(" = ");
	P.Ti=stod(line.substr(found+3, line.length()-found-3));
	printf("Ti %f\n", P.Ti);
    getline(inp, line);
    found=line.find_first_of(" = ");
	P.nu_i=stod(line.substr(found+3, line.length()-found-3));
	printf("nu_i %f\n", P.nu_i);
    getline(inp, line);
    found=line.find_first_of(" = ");
	P.nu_e=stod(line.substr(found+3, line.length()-found-3));
	printf("nu_e %f\n", P.nu_e);
    getline(inp, line);
    found=line.find_first_of("[");
    while(found!=-1)
    {
        double M=0.0, Cm=0.0;
        M=stod(line.substr(found+1, 3));
        found=line.find(" = ");
        Cm=stod(line.substr(found+3, 3));
        C[unsigned(M)]=Cm;
        printf("Con[%3.0f] = %3.0f\n", M, Cm);
        getline(inp, line);
        found=line.find_first_of("[");
    }
    P.Con=C;

    S=Spectrum(P);

	out.open("Spectrum.dat");
	for(int f=0; f<2*FLENGTH; f++)
	{
		out<<double(f-FLENGTH)*double(DELTAF)<<"\t";
		out<<S[f]<<"\n";
	}
	out.close();
}
