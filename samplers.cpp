#include "samplers.h"



void OBABO::print_sampler_params(){
    
    std:: cout << "Sampler parameters:\n";
    std:: cout << "Temperature = " << T << ",\nFriction = " << gamma << ",\nStepsize = " << h << ".\n" << std:: endl; 

};


measurement OBABO::run_mpi_simulation(const int N, const bool tavg, const PROBLEM POTCLASS){
    
    int seed = 0;
    print_sampler_params();

    measurement RESULTS = OBABO::collect_samples(N, tavg, POTCLASS, seed);

    return RESULTS;

};



measurement OBABO::collect_samples(const int N, const bool tavg, PROBLEM problem, const int randomseed){

    std:: cout << "Starting OBABO simulation...\n" << std:: endl;

    // set integrator constants
    const double a = exp(-1*gamma*h);    
    const double sqrt_a = sqrt(a);
    const double sqrt_aT = sqrt((1-a)*T);
    const double h_half = 0.5*h;   

    size_t No_params = problem.parameters.size();  // number of parameters

    // COMPUTE INITIAL FORCES
    problem.compute_force();

    // COMPUTE INITIAL MEASUREMENTS
    measurement RESULTS;
    RESULTS.take_measurement(problem.parameters, problem.velocities);

    // PREPARE RNG
    std:: mt19937 twister;
    
    std:: seed_seq seq{1,20,3200,403,5*randomseed+1,12000,73667,9474+randomseed,19151-randomseed};
    std:: vector < std::uint32_t > seeds(1);
    seq.generate(seeds.begin(), seeds.end());
    twister.seed(seeds.at(0)); 

	std:: normal_distribution<> normal{0,1};
    double Rn;


	auto t1 = std:: chrono::high_resolution_clock::now();

    // MAIN LOOP
    for ( size_t i = 1;  i <= N;  ++i ) {

        // O + B steps
        for ( size_t j = 0;  j < No_params;  ++j ) {                  
            
            Rn = normal(twister);
            problem.velocities[j] = sqrt_a * problem.velocities[j]  +  sqrt_aT * Rn  +  h_half * problem.forces[j]; 
        
        }

        // A step
        for ( size_t j = 0; j < No_params;  ++j ) {

            problem.parameters[j] += h * problem.velocities[j];
        
        }	
  
        // COMPUTE NEW FORCES
        problem.compute_force();

        // B STEP
        for ( size_t j = 0; j < No_params;  ++j ) {

            problem.velocities[j] += h_half * problem.forces[j];
        
        }   
	
        // O STEP
        for ( size_t j = 0; j < No_params;  ++j ) {
            
		    Rn = normal(twister);
            problem.velocities[j] = sqrt_a * problem.velocities[j]  +  sqrt_aT * Rn;
        
        }   								 

        // TAKE MEASUREMENT
		if( i % 5 == 0 ) {                                                          // turn 5 into variable (needs to be passed 
            RESULTS.take_measurement(problem.parameters, problem.velocities);       // to measurement print function as well).

		}
		
        if( i % int(1e5) == 0 ) std:: cout << "Iteration " << i << " done!\n";
	
	}  // END MAIN LOOP

    // FINALIZE
    auto t2 = std:: chrono:: high_resolution_clock:: now();
	auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
	std:: cout << "Execution took " << ms_int.count() << " seconds!\n ";
        
    return RESULTS;
};








// measurement OBABO_simu(const params param0, const long N, const double h, const double T, const double gamma, const vector <double>& mu1, const vector <double>& mu2, 	
// 					   const vector <double>& SIG1, const vector <double>& SIG2, const double phi1,  const size_t n_meas, const double sig){
    
//    cout<<"Starting OBABO simulation..."<<endl;
// 	auto t1 = chrono::high_resolution_clock::now();

// 	// some constants for the integrator, force, potential...
//    const double a = exp(-1*gamma*h);    // used for integrator
//    const double sqrt_a = sqrt(a);
//    const double sqrt_aT = sqrt((1-a)*T);
//    const double h_half = 0.5*h;      

// 	const double det1 = SIG1[0]*SIG1[2] - SIG1[1]*SIG1[1];
// 	const double det2 = SIG2[0]*SIG2[2] - SIG2[1]*SIG2[1];
// 	const double inv_two_det_SIG1 = 1 / (2*det1) ;    // = 1/2(det(SIGMA1)), used in the exponent of the exponentials in the force
// 	const double inv_two_det_SIG2 = 1 / (2*det2) ;
// 	const double pref1 = phi1 / ( 2*PI*pow( det1, 1.5) );  		 // = prefactors of the exponentials in the force
// 	const double pref2 = (1-phi1) / ( 2*PI*pow( det2, 1.5) );


//    params theta = param0; 
// 	forces force = get_noisy_force_DW(theta, mu1, mu2, SIG1, SIG2, inv_two_det_SIG1, inv_two_det_SIG2, pref1, pref2, sig);  	// HERE FORCES!!
// 	forces full_force = get_noisy_force_DW(theta, mu1, mu2, SIG1, SIG2, inv_two_det_SIG1, inv_two_det_SIG2, pref1, pref2, 0);    // to compute Tconf				
	

// 	// get initial measurements
// 	measurement meas{vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1, 1), 
// 	vector <double> (N/n_meas + 1), vector <double> (N/n_meas + 1)};						// measurements of Tconf, Pacc, <x>, dist. err.	

// 	meas.Tconf[0] = -0.5 * (theta.x*full_force.fx + theta.y*full_force.fy);
// 	meas.x[0] = theta.x;
// 	int k=1;  // gives next index of measurement vectors

// 	// nrbins = nrbins%2==0 ? nrbins : nrbins+1;
// 	// const double deltax = 2*xrange / nrbins;
// 	// const double deltay = 2*yrange / nrbins;
// 	// vector < vector<size_t> > bin_ctr( nrbins+2, vector<size_t>(nrbins+2) );   // counts no of samples in each bin (+2 for outermost bins)
// 	// vector < vector<double> > theo_probs = read_real_probs("true_probs_T0.0005.csv", nrbins);  // reads in theoretical probabilities in bins from file 
// 	// get_bin(bin_ctr, theta, xrange, yrange, deltax, deltay, nrbins);
// 	// meas.dist_err[0] = get_dist_err(bin_ctr, theo_probs, 1, nrbins);

// 	double Rnx, Rny;
// 	normal_distribution<> normal{0,1};
	
// 	// cout << 0 << " "<< theta.y << " "<< theta.py << " "<< force.fy<<endl;
// 	for(size_t i=1; i<=N; ++i){
// 		// cout << "Iteration " << i << "\n";

// 		Rnx = normal(twister);
// 		Rny = normal(twister);
// 		theta.px = sqrt_a*theta.px + sqrt_aT*Rnx + h_half*force.fx;  // O+B step		
// 		theta.py = sqrt_a*theta.py + sqrt_aT*Rny + h_half*force.fy;		
														 
// 		theta.x += h*theta.px;															// A step							
// 		theta.y += h*theta.py;		
  
// 		force = get_noisy_force_DW(theta, mu1, mu2, SIG1, SIG2, inv_two_det_SIG1, inv_two_det_SIG2, pref1, pref2, sig);  // HERE FORCES!!!

// 		theta.px += h_half*force.fx;													// B step
// 		theta.py += h_half*force.fy;		

// 		Rnx = normal(twister);
// 		Rny = normal(twister);								 
// 		theta.px = sqrt_a*theta.px + sqrt_aT*Rnx;							 // O step
// 		theta.py = sqrt_a*theta.py + sqrt_aT*Rny;

// 		// get_bin(bin_ctr, theta, xrange, yrange, deltax, deltay, nrbins);

// 		if(i%n_meas == 0 ) {
// 			full_force = get_noisy_force_DW(theta, mu1, mu2, SIG1, SIG2, inv_two_det_SIG1, inv_two_det_SIG2, pref1, pref2, 0); 	// HERE FORCES!!!!
// 			meas.Tconf[k] = -0.5 * (theta.x*full_force.fx + theta.y*full_force.fy);
// 			meas.x[k] = theta.x;
// 			// meas.dist_err[k] = get_dist_err(bin_ctr, theo_probs, i+1, nrbins);
// 			++k;		
// 		}
// 		if(i%int(1e5)==0) cout<<"Iteration "<<i<<" done!"<<endl;
	
// 		// cout << i << " "<< theta.y << " "<< theta.py << " "<< force.fy<<endl;
	

// 	}

	
// 	auto t2 = chrono::high_resolution_clock::now();
// 	auto ms_int = chrono::duration_cast<chrono::seconds>(t2 - t1);
// 	cout<<"Execution took "<< ms_int.count() << " seconds!"<<endl;
	    
//    return meas;

// }