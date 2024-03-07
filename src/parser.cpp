#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

#include "general.h"

using namespace std;


void load_default_values(struct elastic *elastic, double* k_spring, int* Lx, int* Ly, int* y0, double* external_stress_init, double* external_stress_rate, int* seed, bool* is_edge_original)
{
	/// ELASTIC PARAMETERS
	(*elastic).b = 1.0;
	(*elastic).mu = 1.0;
	(*elastic).Poisson = 0.35;
	(*elastic).a = 1.0; /// should not be modified to get consistent line tension
	(*k_spring) = 1e-2;

	/// LATTICE PARAMETERS
	(*Lx) = 100;
	(*Ly) = 100;
	(*y0) = 50;
	(*is_edge_original) = false;


	/// LOADING PARAMETERS
	(*external_stress_init) = 0.0;
	(*external_stress_rate) = 0.0;

	/// RANDOM SEED
	(*seed) = 0;

}

int parser(const std::string filename, struct elastic *elastic, double* k_spring, int* Lx, int* Ly, int* y0, double* external_stress_init, double* external_stress_rate, int* seed, bool* is_edge_original)
{

    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        std::cerr << "Error opening control file: " << filename << std::endl;
        return 1;
    }

    int status = 0;

    std::string line;
    std::string token;

    int pos;
    int found_paramter = 0;

    while (std::getline(file, line))
    {
    	/// Ignoring comments starting with '#'
        pos = line.find('#');
        line = line.substr(0, pos);
        line.erase(remove(line.begin(), line.end(), ' '), line.end());
        pos = line.find('=');

        found_paramter = 0;

        if (pos > -1)
        {
            token = line.substr(0, pos);
            line = line.substr(pos, line.size());
            if (line.size() > 1)
            {
                line = line.substr(1, line.size());

                if (token=="Burgers")
                {
                    (*elastic).b = stod(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="mu")
                {
                    (*elastic).mu = stod(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="Poisson")
                {
                    (*elastic).Poisson = stod(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="cellsize")
                {
                    (*elastic).a = stod(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="spring_constant")
                {
                    (*k_spring) = stod(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="Lx")
                {
                    (*Lx) = stoi(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="Ly")
                {
                    (*Ly) = stoi(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="initial_height")
                {
                    (*y0) = stoi(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="stress_init")
                {
                    (*external_stress_init) = stod(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="stress_rate")
                {
                    (*external_stress_rate) = stod(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="seed")
                {
                    (*seed) = stoi(line);
                    found_paramter = 1;
                }
                if (found_paramter == 0 && token=="vertical_burgers")
                {
                	if (stoi(line)==0)
                    		(*is_edge_original) = false;
                    	else
                    		(*is_edge_original) = true;
                    found_paramter = 1;
                }
                if (found_paramter == 0)
                {
                    std::cerr << "ERROR: Wrong parameter: " << token << std::endl;
                    status = 1;
                }
            }
            else
            {
                std::cerr << "ERROR: No parameter given after =" << std::endl;
                status = 1;
            }

        }
        else
        {
            if (line.size()>0)
            {
                std::cerr << "ERROR: Wrong parameter: " << line << std::endl;
                status = 1;
            }
            else
            {
                ; /// Just a comment line
            }
        }
    }

    return status;
     
}

void create_metafile(const struct elastic elastic, const double k_spring, const int Lx, const int Ly, const int y0, const double external_stress_init, const double external_stress_rate, const int seed, const std::string filename, const std::string pinning_filename, const std::string directory, const bool is_edge_original)
{

	std::ofstream metafile ( directory + "/meta");
    	
    	metafile << "ELASTIC PARAMETERS" << std::endl;
    	metafile << "Length of Burgers vector: " << elastic.b << std::endl;
    	metafile << "Shear modulus: " << elastic.mu << std::endl;
    	metafile << "Poisson ratio: " << elastic.Poisson << std::endl;
    	metafile << "Linear size of cells: " << elastic.a << std::endl;
    	metafile << "Spring constant: " << k_spring << std::endl << std::endl;
    
    	metafile << "LATTICE PARAMETERS" << std::endl;
    	metafile << "Lattice width in cells: " << Lx << std::endl;
    	metafile << "Lattice height in cells: " << Ly << std::endl;
    	metafile << "Initial vertical position of dislocation line: " << y0 << std::endl;
    	if (is_edge_original)
    		metafile << "Is the Burgers vector vertical: YES" << std::endl << std::endl;
    	else
    		metafile << "Is the Burgers vector vertical: NO" << std::endl << std::endl;
    
    	metafile << "LOADING PARAMETERS" << std::endl;
    	metafile << "Initial external stress: " << external_stress_init << std::endl;
     	metafile << "External stress rate: " << external_stress_rate << std::endl << std::endl;
    
    	metafile << "RANDOM SEED" << std::endl;
    	metafile << "Random seed: " << seed << std::endl << std::endl;
 
    	metafile << "PATHS" << std::endl;
    	metafile << "Path to control file: " << filename << std::endl;
    	metafile << "Path to pinning field file: " << pinning_filename << std::endl;
    	metafile << "Path to output directory: " << directory << std::endl << std::endl;
    	metafile.close();
}
