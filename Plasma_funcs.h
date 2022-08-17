#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <complex>
#include <vector>
#include <map>

#define PI 3.14159265
#define FLENGTH 5000
#define DELTAF 10

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
        //N:7
        //N2:14
        //NO:15
        //O2:16
        //O:8
        //H:2
        //He:1
};

complex<double> cintegral(complex<double> a, complex<double> b, unsigned step_count);
void Spectrum(vector<double> &x, vector<double> &y, Plasma_pars P);
void Set_pars(char* file, Plasma_pars &P);