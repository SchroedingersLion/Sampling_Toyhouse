#include "samplers.h"

int main(int argc, char *argv[]){

    double T = 1;
    double gamma = 10;
    double h = 0.01;

    int N = 10000000;
    bool tavg = false;

    OBABO testsampler(T, gamma, h);

    testsampler.print_sampler_params();

    PROBLEM EMPTY_POT;

    
    measurement EMPTY_meas = testsampler.run_mpi_simulation(N, tavg, EMPTY_POT);



}