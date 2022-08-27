#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <complex>
#include <vector>
#include <map>
#include "fftw3.h"
#include "FFT.h"

#define PI 3.14159265
#define FLENGTH 5000
#define DELTAF 10
#define DELTAT 10
#define LENGTH 100

using namespace std;

struct Plasma_pars
{
    double lambda; //Длина волны радара, м
    double Ne; //Концентрация=Ne*10^11, 1/м^3
    double Te; //Температура электронов, K
    double Ti; //Температура ионов, K
    double nu_i; //Частота столкновений ионов, Гц
    double nu_e; //Частота столкновений электронов, Гц. Для бесстолкновительной плазмы обе частоты нужно занулить.
    double Con[103]; //Контейнер для ионного состава, где первый компонент - молярная масса иона,
        //второй - его процентное сожержание в плазме. Например, для однородой O+ плазмы: Con[8]=100
        //Для 30% NO+ и 70% O+: Con[15]=30, Con[8]=70. Размер контейнера определяет кол-во сортов ионов. 
        //Молярные массы сортов частиц, обитающих в ионосфере Земли:
		//O:8
		//H:2
        //He:1
        //N:7
        //N2:14
        //NO:15
        //O2:16
};

complex<double> cintegral(complex<double> a, complex<double> b, unsigned step_count);
extern "C" void Spectrum(double *x, double *y, Plasma_pars *P);
extern "C" void ACF(double *x, double *y, double *S);
extern "C" void Set_pars(char* file, Plasma_pars *P);
