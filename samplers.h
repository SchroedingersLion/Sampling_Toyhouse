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



class OBABO{

    private:
        const double T;
        const double gamma;
        const double h;

        measurement collect_samples(const int N, const bool tavg, PROBLEM POTCLASS, const int randomseed);

    public:
        // constructors
        OBABO(double T, double gamma, double h): T{T}, gamma{gamma}, h{h} {
        } 

        measurement run_mpi_simulation(const int N, const bool tavg, PROBLEM POTCLASS);  

        void print_sampler_params();




};

















#endif // SAMPLERS_H