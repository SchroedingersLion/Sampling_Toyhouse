#ifndef SETUP_CLASSES_H
#define SETUP_CLASSES_H

#include <vector>
#include <cmath>
#include <math.h>
#include <fstream>

constexpr double PI = 3.141592653589793;    // some useful constants
constexpr double eps = 1e-12;

class PROBLEM { 
    
    // THIS PROBLEM IS THE DOUBLE WELL POTENTIAL IN 2 DIMENSIONS
    
    /* In case one wants to change the problem to something else, eg. harmonic oscillator:
       The samplers require a "compute_force" function as well as the vectors "parameters", "velocities", and "forces"
       storing the corresponding quantities. */

    public:

        // parameters, velocities, forces 
        std:: vector <double> parameters {1,0};       // CARE! THE SAMPLERS REQUIRE THIS VECTOR
        std:: vector <double> velocities {0,0};       // CARE! THE SAMPLERS REQUIRE THIS VECTOR
        std:: vector <double> forces {0,0};           // CARE! THE SAMPLERS REQUIRE THIS VECTOR

        // fill force vector
        void compute_force(){                           /* CARE! THE SAMPLERS REQUIRE THIS FUNCTION. 
                                                           It needs to be written by the user. */

                 
            double x = parameters[0];       // current parameters (positions)
            double y = parameters[1];

            double e1, e2, rho;             // some help constants
            double x_mu_diff1 = x-mu1[0]; 
            double x_mu_diff2 = x-mu2[0]; 
            double y_mu_diff1 = y-mu1[1]; 
            double y_mu_diff2 = y-mu2[1];
            
            // the exponentials 
            e1 = exp( -inv_two_det_SIG1 * ( SIG1[2]*x_mu_diff1*x_mu_diff1 - 2*SIG1[1]*x_mu_diff1*y_mu_diff1 + SIG1[0]*y_mu_diff1*y_mu_diff1 ) );
            e2 = exp( -inv_two_det_SIG2 * ( SIG2[2]*x_mu_diff2*x_mu_diff2 - 2*SIG2[1]*x_mu_diff2*y_mu_diff2 + SIG2[0]*y_mu_diff2*y_mu_diff2 ) );

            rho = pref1 * det1 * e1  +  pref2 * det2 * e2;  // the density

            //force part without noise, i.e. double well
            forces[0] = -1/rho * ( pref1 * e1 * (SIG1[2]*x_mu_diff1 - SIG1[1]*y_mu_diff1)   +   pref2 * e2 * (SIG2[2]*x_mu_diff2 - SIG2[1]*y_mu_diff2) );
            forces[1] = -1/rho * ( pref1 * e1 * (-SIG1[1]*x_mu_diff1 + SIG1[0]*y_mu_diff1)  +   pref2 * e2 * (-SIG2[1]*x_mu_diff2 + SIG2[0]*y_mu_diff2) );

            // if(sig!=0){
            //     normal_distribution<> normal{0,sig};	// add noise
            //     F.fx += normal(twister);
            //     F.fy += normal(twister);
            // }

            return;

        };


    private:

        // constants that define the potential
        const std:: vector <double> mu1 {1, 0};		                // (x,y) of the first well 
        const std:: vector <double> mu2 {-1, 0};		            // (x,y) of the second well
        const std:: vector <double> SIG1 {0.6, 0.085, 0.02};	    // elements of the two covariance matrices sig11, sig12, sig22 (note: sig12=sig21)
        const std:: vector <double> SIG2 {0.6, 0.085, 0.02};
        const double phi1 = 0.5;                                    // mixing factor (phi2 = 1-phi1)	

        // constants used by the force computation 
        const double det1 = SIG1[0]*SIG1[2] - SIG1[1]*SIG1[1];      // determinants of cov. matrices
        const double det2 = SIG2[0]*SIG2[2] - SIG2[1]*SIG2[1];      
        const double inv_two_det_SIG1 = 1 / (2*det1) ;              // used in the exponents in the force
        const double inv_two_det_SIG2 = 1 / (2*det2) ;
        const double pref1 = phi1 / ( 2*PI*pow( det1, 1.5) );  	    // prefactors in the force
        const double pref2 = (1-phi1) / ( 2*PI*pow( det2, 1.5) );


};



class measurement{

    /* The measurement class that defines what quantities are collected by the samplers, how they are computed, and
       printed to a file. This class needs to be modified by the user.
       In this example, we only save the first parameter (corresponding to the x-coordinate
       when used in the double well problem above), as well as the kinetic energy. */

    public:

        void take_measurement(std:: vector <double> parameters, std:: vector <double> velocities){  /* CARE!!! The samplers need this function.
                                                                                                       It needs to be written by the user. */
            
            x.push_back(parameters[0]);
            kin_energy.push_back( 0.5* (velocities[0]*velocities[0] + velocities[1]*velocities[1]) );
        
        }

        void print_to_csv(const int t_meas){    /* print results to file. this routine would be called in the main-file. 
                                                   It needs to be written by the user. "t_meas" gives the number of sampler iterations 
                                                   between two measurments (it is passed here to get the correct iteration count in the
                                                   output file.)*/

            std:: ofstream file{"TEST.csv"};
            std:: cout << "Writing to file...\n";

            for ( size_t i = 0; i<x.size(); ++i )
            {
                file << i*t_meas << " " << x.at(i) << " " << kin_energy.at(i) <<  "\n";  
            }
            file.close();

        }

        std:: vector <double> x; 
        std:: vector <double> kin_energy;



};


#endif // SETUP_CLASSES_H