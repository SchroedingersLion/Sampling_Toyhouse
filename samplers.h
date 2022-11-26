#ifndef SAMPLERS_H
#define SAMPLERS_H


#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <chrono>
#include <math.h> 
#include <iomanip>
#include <algorithm>
#include <limits>
#include <numeric>
#include <iterator>
#include "setup_classes.h"
// #include <mpi.h>


constexpr double PI = 3.141592653589793;
constexpr double eps = 1e-12;


class OBABO{

    private:
        const double T;
        const double gamma;
        const double h;

        measurement collect_samples(const int N, const bool tavg, const potential POTCLASS, const int randomseed);

    public:
        // constructors
        OBABO(double T, double gamma, double h): T{T}, gamma{gamma}, h{h} {

        } 

        measurement run_mpi_simulation(const int N, const bool tavg, const potential POTCLASS);  
        // void print_pos_vels();
        void print_sampler_params();




};

















#endif // SAMPLERS_H