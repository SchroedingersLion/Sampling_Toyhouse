#include "samplers.h"
#include "setup_classes.h"

int main(int argc, char *argv[]){

    // we want to sample the double well potential using the OBABO sampler


    double T = 1;           // sampler hyperparameters
    double gamma = 1;
    double h = 0.01;

    int N = 50000;          // no. of iteration
    bool tavg = false;      // perform time-average?
    int t_meas = 5;         // take measurement every t_meas iterations (passed to sampler as well as print functions).

    // OBABO testsampler(T, gamma, h);     // construct OBABO object defined in header "samplers.h"
    SGHMC testsampler(T,gamma,h);

    PROBLEM double_well;    /* construct object of the problem class defined by the user in header "setup_classes.h". */
    
    measurement results = testsampler.run_mpi_simulation(N, tavg, double_well, t_meas);  /* run sampler. "measurement" needs to be defined 
                                                                                            by user in header "setup_classes.h". It holds
                                                                                            information of what quantities need to be obtained
                                                                                            by the sampler and how to compute and print them.  */

    results.print_to_csv(t_meas);     // print to file, as specified by the user in the "measurement" class.

    return 0;


}