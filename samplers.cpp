#include "samplers.h"



void OBABO::print_sampler_params(){
    
    std:: cout << "Print sampler parameters:\n";
    std:: cout << "Temperature = " << T << ",\nFriction = " << gamma << ",\nStepsize = " << h << "." << std:: endl; 

};


measurement OBABO::run_mpi_simulation(const int N, const bool tavg, const potential POTCLASS){
    
    measurement RESULTS = OBABO::collect_samples(N, tavg, POTCLASS);


};



measurement OBABO::collect_samples(const int N, const bool tavg, const potential POTCLASS){

    std:: cout << "Starting to collect samples \n" << std:: endl;
    measurement RESULTS;
    return RESULTS;
};