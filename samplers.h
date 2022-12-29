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
#include <mpi.h>


/* The OBABO sampler. It requires the "PROBLEM" class and the "measurement" class in header "setup_classes.h".
   They need to be written by the user and define the sampling problem and what measurements to take. 
   Note that the member functions below require those classes to be structured in a certain way. */

class OBABO{

    private:
        const double T;
        const double gamma;
        const double h;

        measurement collect_samples(const int max_iter, PROBLEM POTCLASS, const int randomseed, const int t_meas);  // draws a single sampling trajectory

    public:
        // constructors
        OBABO(double T, double gamma, double h): T{T}, gamma{gamma}, h{h} {
        } 

        void run_mpi_simulation(int argc, char *argv[], const int max_iter, PROBLEM POTCLASS, const int t_meas, const bool tavg=0, int n_tavg=10, const int n_dist=1);  
        /* sets up mpi environment and calls "collect_samples" on each process within. Also performs averaging. */

        void print_sampler_params();    // print sampler hyperparameters.

};



class SGHMC{

    private:
        const double T;
        const double gamma;
        const double h;

        measurement collect_samples(const int max_iter, PROBLEM POTCLASS, const int randomseed, const int t_meas);

    public:
        SGHMC(double T, double gamma, double h): T{T}, gamma{gamma}, h{h} {
        }

        measurement run_mpi_simulation(const int max_iter, const bool tavg, PROBLEM POTCLASS, const int t_meas);

        void print_sampler_params();

};









#endif // SAMPLERS_H