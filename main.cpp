#include "samplers.h"
#include "setup_classes.h"

int main(int argc, char *argv[]){

    // we want to sample the double well potential using the OBABO sampler


    double T = 1;           // sampler hyperparameters
    double gamma = 1;
    double h = 0.01;

    int N = 10;          // no. of iteration
    bool tavg = false;      // perform time-average?
    int n_tavg = 100;        // t-average over the last n_tavg values.
    int t_meas = 1;         // take measurement every t_meas iterations (passed to sampler as well as print functions).
    int n_dist = 1;         // print and t-average (if activated) only every n_dist-th values.

    OBABO testsampler(T, gamma, h);     // construct OBABO object defined in header "samplers.h"
    // SGHMC testsampler(T,gamma,h);

    std:: string filename = "GM_data_5000.csv";
    const int randomseed = 0;
    const int batchsize = 5;

    PROBLEM gm_mix_esti(filename, batchsize, randomseed);    /* construct object of the problem class defined by the user in header "setup_classes.h". */
    
    testsampler.run_mpi_simulation(argc, argv, N, gm_mix_esti, t_meas, tavg, n_tavg, n_dist);  /* run sampler. "measurement" needs to be defined 
                                                                                                            by user in header "setup_classes.h". It holds
                                                                                                            information of what quantities need to be obtained
                                                                                                            by the sampler and how to compute and print them.  */

    std:: cout<<"reached end"<<std::endl;
    return 0;


}