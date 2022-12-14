#include "samplers.h"

// ##### OBABO METHODS ##### //

void OBABO::print_sampler_params(){
    
    std:: cout << "Sampler parameters:\n";
    std:: cout << "Temperature = " << T << ",\nFriction = " << gamma << ",\nStepsize = " << h << ".\n" << std:: endl; 

};



void OBABO::run_mpi_simulation(int argc, char *argv[], const int max_iter, PROBLEM POTCLASS, const std:: string outputfile, const int t_meas, const bool tavg, int n_tavg, const int n_dist){
    
    MPI_Init(&argc, &argv);				// initialize MPI, use rank as random seed
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nr_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nr_proc);

    const int seed = rank;

    print_sampler_params();
    measurement RESULTS = OBABO::collect_samples(max_iter, POTCLASS, seed, t_meas);    // run sampler

    std:: cout << "Rank " << rank << " reached barrier." << std:: endl;
    MPI_Barrier(comm);

    // sum up results of different processors 
    measurement RESULTS_AVG;
    int row_size;

    for ( size_t i = 0;  i < RESULTS.measured_values.size();  ++i){
        
        row_size = RESULTS.measured_values[i].size();
        if( rank==0 ) RESULTS_AVG.measured_values[i].resize( row_size );     // only on rank 0 to save RAM.

        MPI_Reduce( &RESULTS.measured_values[i][0], &RESULTS_AVG.measured_values[i][0], row_size, MPI_DOUBLE, MPI_SUM, 0, comm);

    }	 

    // compute average
    if( rank==0 ){

        // average over processes
        std:: cout << "Average over processes..." << std:: endl;
        for ( size_t i = 0;  i < RESULTS_AVG.measured_values.size();  ++i){
            row_size = RESULTS_AVG.measured_values[i].size();
            for ( int j = 0;  j < row_size;  ++j ){
                RESULTS_AVG.measured_values[i][j] /= nr_proc;
            }
        }
        
    // time average
        if ( tavg == true ){

            std:: cout << "Time averaging...\n";
            int n_tavg_temp;
            
            for ( size_t row = 0;  row < RESULTS_AVG.measured_values.size();  ++row ){   // t-average one row at a time.
            
                n_tavg_temp = n_tavg;
            
                for ( int i = RESULTS_AVG.measured_values[row].size() - 1;  i >= 0;  i -= n_dist ){
                
                    if ( i <= n_tavg_temp - 1 ) n_tavg_temp = i;
                    for ( int j = i - n_tavg_temp;  j < i;  ++j ){
                        RESULTS_AVG.measured_values[row][i] += RESULTS_AVG.measured_values[row][j];

                    }
                    RESULTS_AVG.measured_values[row][i] /= (n_tavg_temp + 1);

                }
            }
        
        } // end time average

    }

    if ( rank == 0) RESULTS_AVG.print_to_csv(t_meas, n_dist, outputfile);     // print to file, as specified by the user in the "measurement" class.

    MPI_Finalize();

    return;

};



measurement OBABO::collect_samples(const int max_iter, PROBLEM problem, const int randomseed, const int t_meas){

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
    for ( size_t i = 1;  i <= max_iter;  ++i ) {

        // O + B steps
        for ( size_t j = 0;  j < No_params;  ++j ) {                  
            
            Rn = normal(twister);
            problem.velocities[j] = sqrt_a * problem.velocities[j]  +  sqrt_aT * Rn  +  h_half * problem.forces[j]; 
        
        }

        // A step
        for ( size_t j = 0;  j < No_params;  ++j ) {

            problem.parameters[j] += h * problem.velocities[j];
        
        }	
  
        // COMPUTE NEW FORCES
        problem.compute_force();

        // B STEP
        for ( size_t j = 0;  j < No_params;  ++j ) {

            problem.velocities[j] += h_half * problem.forces[j];
        
        }   
	
        // O STEP
        for ( size_t j = 0;  j < No_params;  ++j ) {
            
		    Rn = normal(twister);
            problem.velocities[j] = sqrt_a * problem.velocities[j]  +  sqrt_aT * Rn;
        
        }   								 

        // TAKE MEASUREMENT
		if( i % t_meas == 0 ) {                                                 // take measurement any t_meas steps.        
            RESULTS.take_measurement(problem.parameters, problem.velocities); 
		}
		
        if( i % int(1e5) == 0 ) std:: cout << "Iteration " << i << " done!\n";
	
	}  // END MAIN LOOP

    // FINALIZE
    auto t2 = std:: chrono:: high_resolution_clock:: now();
	auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
	std:: cout << "Execution took " << ms_int.count() << " seconds!\n";
        
    return RESULTS;
};




// ##### SGHMC METHODS ##### //

void SGHMC::print_sampler_params(){
    
    std:: cout << "Sampler parameters:\n";
    std:: cout << "Temperature = " << T << ",\nFriction = " << gamma << ",\nStepsize = " << h << ".\n" << std:: endl; 

};



void SGHMC::run_mpi_simulation(int argc, char *argv[], const int max_iter, PROBLEM POTCLASS, const std:: string outputfile, const int t_meas, const bool tavg, int n_tavg, const int n_dist){
    
    MPI_Init(&argc, &argv);				// initialize MPI, use rank as random seed
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nr_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nr_proc);

    const int seed = rank;

    print_sampler_params();
    measurement RESULTS = SGHMC::collect_samples(max_iter, POTCLASS, seed, t_meas);    // run sampler

    std:: cout << "Rank " << rank << " reached barrier." << std:: endl;
    MPI_Barrier(comm);

    // sum up results of different processors 
    measurement RESULTS_AVG;
    int row_size;

    for ( size_t i = 0;  i < RESULTS.measured_values.size();  ++i){
        
        row_size = RESULTS.measured_values[i].size();
        if( rank==0 ) RESULTS_AVG.measured_values[i].resize( row_size );     // only on rank 0 to save RAM.

        MPI_Reduce( &RESULTS.measured_values[i][0], &RESULTS_AVG.measured_values[i][0], row_size, MPI_DOUBLE, MPI_SUM, 0, comm);

    }	 

    // compute average
    if( rank==0 ){

        // average over processes
        std:: cout << "Average over processes..." << std:: endl;
        for ( size_t i = 0;  i < RESULTS_AVG.measured_values.size();  ++i){
            row_size = RESULTS_AVG.measured_values[i].size();
            for ( int j = 0;  j < row_size;  ++j ){
                RESULTS_AVG.measured_values[i][j] /= nr_proc;
            }
        }
        
    // time average
        if ( tavg == true ){

            std:: cout << "Time averaging...\n";
            int n_tavg_temp;
            
            for ( size_t row = 0;  row < RESULTS_AVG.measured_values.size();  ++row ){   // t-average one row at a time.
            
                n_tavg_temp = n_tavg;
            
                for ( int i = RESULTS_AVG.measured_values[row].size() - 1;  i >= 0;  i -= n_dist ){
                
                    if ( i <= n_tavg_temp - 1 ) n_tavg_temp = i;
                    for ( int j = i - n_tavg_temp;  j < i;  ++j ){
                        RESULTS_AVG.measured_values[row][i] += RESULTS_AVG.measured_values[row][j];

                    }
                    RESULTS_AVG.measured_values[row][i] /= (n_tavg_temp + 1);

                }
            }
        
        } // end time average

    }

    if ( rank == 0) RESULTS_AVG.print_to_csv(t_meas, n_dist, outputfile);     // print to file, as specified by the user in the "measurement" class.

    MPI_Finalize();

    return;

};



measurement SGHMC::collect_samples(const int max_iter, PROBLEM problem, const int randomseed, const int t_meas){

    std:: cout << "Starting SGHMC simulation...\n" << std:: endl;

    // set integrator constants
    const double one_minus_hgamma = 1-h*gamma;    
    const double noise_pref = sqrt(2*h*gamma*T);  

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
    for ( size_t i = 1;  i <= max_iter;  ++i ) {
        
        // UPDATE
        for ( size_t j = 0;  j < No_params;  ++j ) {                  
            
            Rn = normal(twister);
            problem.velocities[j] = one_minus_hgamma * problem.velocities[j]  +  noise_pref * Rn  +  h * problem.forces[j]; // momentum update
            problem.parameters[j] += h * problem.velocities[j];                                                             // parameter update
        
        }
  
        // COMPUTE NEW FORCES
        problem.compute_force();
					 

        // TAKE MEASUREMENT
		if( i % t_meas == 0 ) {                                                 // take measurement any t_meas steps.        
            RESULTS.take_measurement(problem.parameters, problem.velocities); 
		}
		
        if( i % int(1e5) == 0 ) std:: cout << "Iteration " << i << " done!\n";
	
	}  // END MAIN LOOP


    // FINALIZE
    auto t2 = std:: chrono:: high_resolution_clock:: now();
	auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
	std:: cout << "Execution took " << ms_int.count() << " seconds!\n";
        
    return RESULTS;

};





// ##### BBK_AMAGOLD METHODS ##### //

void BBK_AMAGOLD::print_sampler_params(){
    
    std:: cout << "Sampler parameters:\n";
    std:: cout << "Temperature = " << T << ",\nFriction = " << gamma << ",\nStepsize = " << h << ".\n" << std:: endl; 

};



void BBK_AMAGOLD::run_mpi_simulation(int argc, char *argv[], const int max_iter, PROBLEM POTCLASS, const std:: string outputfile, const int t_meas, const bool tavg, int n_tavg, const int n_dist){
    
    MPI_Init(&argc, &argv);				// initialize MPI, use rank as random seed
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nr_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nr_proc);

    const int seed = rank;

    print_sampler_params();
    measurement RESULTS = BBK_AMAGOLD::collect_samples(max_iter, POTCLASS, seed, t_meas);    // run sampler

    std:: cout << "Rank " << rank << " reached barrier." << std:: endl;
    MPI_Barrier(comm);

    // sum up results of different processors 
    measurement RESULTS_AVG;
    int row_size;

    for ( size_t i = 0;  i < RESULTS.measured_values.size();  ++i){
        
        row_size = RESULTS.measured_values[i].size();
        if( rank==0 ) RESULTS_AVG.measured_values[i].resize( row_size );     // only on rank 0 to save RAM.

        MPI_Reduce( &RESULTS.measured_values[i][0], &RESULTS_AVG.measured_values[i][0], row_size, MPI_DOUBLE, MPI_SUM, 0, comm);

    }	 

    // compute average
    if( rank==0 ){

        // average over processes
        std:: cout << "Average over processes..." << std:: endl;
        for ( size_t i = 0;  i < RESULTS_AVG.measured_values.size();  ++i){
            row_size = RESULTS_AVG.measured_values[i].size();
            for ( int j = 0;  j < row_size;  ++j ){
                RESULTS_AVG.measured_values[i][j] /= nr_proc;
            }
        }
        
    // time average
        if ( tavg == true ){

            std:: cout << "Time averaging...\n";
            int n_tavg_temp;
            
            for ( size_t row = 0;  row < RESULTS_AVG.measured_values.size();  ++row ){   // t-average one row at a time.
            
                n_tavg_temp = n_tavg;
            
                for ( int i = RESULTS_AVG.measured_values[row].size() - 1;  i >= 0;  i -= n_dist ){
                
                    if ( i <= n_tavg_temp - 1 ) n_tavg_temp = i;
                    for ( int j = i - n_tavg_temp;  j < i;  ++j ){
                        RESULTS_AVG.measured_values[row][i] += RESULTS_AVG.measured_values[row][j];

                    }
                    RESULTS_AVG.measured_values[row][i] /= (n_tavg_temp + 1);

                }
            }
        
        } // end time average

    }

    if ( rank == 0) RESULTS_AVG.print_to_csv(t_meas, n_dist, outputfile);     // print to file, as specified by the user in the "measurement" class.

    MPI_Finalize();

    return;

};



measurement BBK_AMAGOLD::collect_samples(const int max_iter, PROBLEM problem, const int randomseed, const int t_meas){

    std:: cout << "Starting BBK simulation...\n" << std:: endl;

    // set integrator constants
    const double one_plus_hgamma_half_inv = 1 / (1+0.5*h*gamma);    
    const double one_minus_hgamma_half = 1-0.5*h*gamma;
    const double noise_pref = sqrt(2*h*gamma*T);  

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
    for ( size_t i = 1;  i <= max_iter;  ++i ) {
        
        // UPDATE
        for ( size_t j = 0;  j < No_params;  ++j ) {                  
            
            problem.velocities[j] *= one_minus_hgamma_half;
            Rn = normal(twister);
            problem.velocities[j] = one_plus_hgamma_half_inv  *  ( problem.velocities[j]  +  h * problem.forces[j]  +  noise_pref * Rn );
            problem.parameters[j] += h * problem.velocities[j];                                                             
        
        }
  
        // COMPUTE NEW FORCES
        problem.compute_force();
					 

        // TAKE MEASUREMENT
		if( i % t_meas == 0 ) {                                                 // take measurement any t_meas steps.        
            RESULTS.take_measurement(problem.parameters, problem.velocities); 
		}
		
        if( i % int(1e5) == 0 ) std:: cout << "Iteration " << i << " done!\n";
	
	}  // END MAIN LOOP


    // FINALIZE
    auto t2 = std:: chrono:: high_resolution_clock:: now();
	auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
	std:: cout << "Execution took " << ms_int.count() << " seconds!\n";
        
    return RESULTS;

};