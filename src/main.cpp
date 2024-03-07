#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <filesystem>
#include <chrono>

#include "general.h"
#include "screw.h"
#include "parser.h"


int main(int argc, char *argv[])
{
    /// Start measuring time
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    if (argc != 4)
    {
        std::cerr << "ERROR: Incorrect number of arguments: " + to_string(argc) << std::endl;
        return ARGC_ERROR;
    }

    struct elastic elastic;
    int Lx, Ly, y0, seed;
    double external_stress_init, external_stress_rate, k_spring;
    bool is_edge_original;

    /// Loading defualt parameter values
    load_default_values(&elastic, &k_spring, &Lx, &Ly, &y0, &external_stress_init, &external_stress_rate, &seed, &is_edge_original);

    /// Updating parameters from control file
    int read_status;
    std::string filename = argv[1];
    std::string pinning_filename = argv[2];
    std::string directory = argv[3];
    read_status = parser(filename, &elastic, &k_spring, &Lx, &Ly, &y0, &external_stress_init, &external_stress_rate, &seed, &is_edge_original);
    if (read_status != 0)
    {
    	std::cerr << "ERROR: Inconsistent control file" << std::endl;
    	return BAD_CONTROL_FILE;
    }

    /// Setting factors in the stress field
    double factor_horizontal = ((is_edge_original) ? 1.0-elastic.Poisson : 1.0);
    double factor_vertical = ((is_edge_original) ? 1.0 : 1.0-elastic.Poisson);

    /// Setting random seed
    std::srand(seed);

    /// Initialization external stress option
    double external_stress = external_stress_init;

    /// Ouput directory
    std::filesystem::create_directory(directory);

    /// Creating metadata file
    create_metafile(elastic, k_spring, Lx, Ly, y0, external_stress_init, external_stress_rate, seed, filename, pinning_filename, directory, is_edge_original);
 
    /// Setting pinning stresses
    std::vector<std::vector<bool>> lattice = create_lattice(Lx, Ly, y0);
    std::vector<std::vector<double>> pinning = read_pinning_field(pinning_filename, Lx, Ly);

    /// Setting initial area to scan
    struct boundaries bounds;
    bounds.bottom = 0;
    bounds.top = Ly-1;

    #ifdef SAVE_LATTICE
    std::filesystem::create_directory(directory + "/lattice/");
    std::ofstream outfile ( directory + "/lattice/" + "0.lat");

    for (int j=0; j<Ly; j++)
    {
        for (int i = 0; i < (Lx-1); i++)
        {
            outfile << lattice[i][j] << "\t";
        }
        outfile << lattice[Lx-1][j] << std::endl;
    }

    outfile.close();
    #endif

    #ifdef SAVE_PINNING_FIELD
    std::ofstream pinfile ( directory + "pinning.fld");

    for (int j=0; j<Ly; j++)
    {
        for (int i = 0; i < (Lx-1); i++)
        {
            pinfile << pinning[i][j] << "\t";
        }
        pinfile << pinning[Lx-1][j] << std::endl;
    }

    pinfile.close();
    #endif


    std::vector<struct segment> segments;
    // std::cout << bounds.bottom << "\t" << bounds.top << std::endl;
    segments = find_segments(Lx, Ly, lattice, bounds);
    
    /// Dealing with initially empty case
    if (segments.size() == 0)
    {
    	/// Computation time
    	std::ofstream timefile ( directory + "time");
    
    	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	double elapsed_time = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000;
	timefile << elapsed_time << std::endl;
	std::cout << "Warning: Simulation stopped because the initial lattice was empty." << std::endl;
	
	int minutes = elapsed_time / 60;
	int hours = minutes / 60;
	minutes = minutes - hours * 60;
	int seconds = (int)(elapsed_time)%60;

	std::cout << "Simulation time: " << hours << "h " << minutes << "m " << seconds << "s" <<  std::endl;

	return EMPTY_LATTICE;    
    }

    #ifdef SAVE_SEGMENTS
    std::filesystem::create_directory(directory + "/segment/");
    std::ofstream segfile ( directory + "/segment/" + "0.seg");

    for (unsigned int i = 0; i < segments.size(); i++)
    {
        segfile << segments[i].x << "\t" << segments[i].y << std::endl;
    }

    segfile.close();
    #endif


    /// Flag indicating the status of the simulation
    int current_state;
    current_state = ONGOING;

    std::ofstream stressfile (directory + "stress");


    /// Timestep counter
    int t_step = 0;
    
    #ifdef SHOW_TIME
    /// Number of estimated steps (~upper bound)
    int estimated_steps = Lx * (Ly - 3) + 1;
    #endif

    /// index of randomly picked segment
    int i_rand;
    
    /// new spin state after flipping
    bool new_spin;
    
    /// mean local stress at the picked segment
    double local_stress;
    
    /// plastic strain
    int plastic_strain = 0;

    std::cout << "Simulation started." << std::endl;

    while (current_state == ONGOING )
    {

	t_step += 1;

	segments = find_segments_optimized(Lx, Ly, lattice, segments);

	/// randomly updating lattice
	i_rand = random_segment_index(segments.size());
	local_stress = get_total_stress(segments, i_rand, pinning, elastic, external_stress, Lx, lattice, factor_horizontal, factor_vertical);
	new_spin = stress_to_spin(local_stress);
	current_state = time_step(&lattice, &plastic_strain, Lx, Ly, segments, i_rand, bounds, new_spin);

	/// Elastic stress drop
	change_external_stress(&external_stress, external_stress_init, external_stress_rate, t_step, k_spring, plastic_strain);
	
	/// Appending to stress file
	stressfile << t_step << "\t" << plastic_strain << "\t" <<  external_stress << std::endl;

	#ifndef SHOW_TIME
	std::cout << t_step << std::endl;
	#endif


	#ifdef SHOW_TIME
	// Estimating remaining time
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	double elapsed_time = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000;
	double estimated_time = elapsed_time * ( estimated_steps - t_step ) / t_step;

	if (estimated_time > 0)
	{
	int minutes = estimated_time / 60;
	int hours = minutes / 60;
	minutes = minutes - hours * 60;
	int seconds = (int)(estimated_time)%60;


	std::cout << t_step << " steps done out of estimated " << estimated_steps << "." << "Estimated time to finish: " << hours << "h " << minutes << "m " << seconds << "s" << std::endl;
	}
	else
	{
	std::cout << t_step << " steps done out of estimated " << estimated_steps << "." << "Simulation will finish soon." << std::endl;
	}
	#endif


	#ifdef SAVE_LATTICE
	if ( (t_step % 10) == 0)
	{

	std::ofstream outfile (directory + "/lattice/" + to_string(t_step) + ".lat");

	for (int j=0; j<Ly; j++)
	{
	for (int i = 0; i < (Lx-1); i++)
	{
	    outfile << lattice[i][j] << "\t";
	}
	outfile << lattice[Lx-1][j] << std::endl;
	}

	outfile.close();
	}
	#endif

	#ifdef SAVE_SEGMENTS
        if ( (t_step % 100) == 0)
        {

	std::ofstream segfile ( directory + "/segment/" + to_string(t_step) + ".seg");

	for (unsigned int i = 0; i < segments.size(); i++)
	{
	segfile << segments[i].x << "\t" << segments[i].y << std::endl;
	}

	segfile.close();
	#endif
	}
    }

    stressfile.close();

    /// Computation time
    std::ofstream timefile ( directory + "time");

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double elapsed_time = (double)(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000;
    timefile << elapsed_time << std::endl;

    /// End messages
    switch(current_state) {
  	case ONGOING:
    		std::cout << "Simulation stopped." << std::endl;
		std::cout << "Warning: status flag is 'ongoing', you should not see this message." << std::endl;
    		break;
  	case REACHED_BOTTOM:
    		std::cout << "Simulation stopped." << std::endl;
		std::cout << "Warning: Simulation stopped due to reaching the bottom of the simulation cell. You may consider a larger lattice." << std::endl;
    		break;
  	case REACHED_TOP:
    		std::cout << "Simulation stopped due to reaching the top of the simulation cell." << std::endl;
    		break;
  	case STARVATION:
    		std::cout << "Simulation stopped due to the annihilation of all dislocation segments. This may be something you don't want to see." << std::endl;
    		break;
  	default:
    		std::cout << "Warning: Incorrect status flag: "<< current_state << std::endl;
    }
    
    int minutes = elapsed_time / 60;
    int hours = minutes / 60;
    minutes = minutes - hours * 60;
    int seconds = (int)(elapsed_time)%60;
    
    std::cout << "Simulation time: " << hours << "h " << minutes << "m " << seconds << "s" <<  std::endl;   

    return 0;
}
