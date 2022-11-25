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
#include <mpi.h>


constexpr double PI = 3.141592653589793;
constexpr double eps = 1e-12;


class OBABO{

    private:
        const double T;
        const double gamma;
        const double h;
        measurement collect_samples(const int N, const bool tavg, const potential POTCLASS);

    public:
        // constructors

        measurement run_mpi_simulation(int argc, char *argv[], const int N, const potential POTCLASS);  
        void print_pos_vels();
        void print_sampler_params();




}

















#endif // SAMPLERS_H