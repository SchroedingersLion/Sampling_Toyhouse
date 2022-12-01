#include "samplers.h"

int main(int argc, char *argv[]){

    double T = 1;
    double gamma = 1;
    double h = 0.01;

    int N = 50000;
    bool tavg = false;

    OBABO testsampler(T, gamma, h);


    PROBLEM double_well;
    
    measurement results = testsampler.run_mpi_simulation(N, tavg, double_well);

    results.print_to_csv();

    return 0;


}