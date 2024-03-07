
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>


#include "general.h"

std::vector<std::vector<bool>> create_lattice(const int Lx, const int Ly, const int y0)
{
    /// Creation of lattice with proper initial state
    /// Lx: lattice width
    /// Ly: lattice height
    /// y0: initial height of 'invaded' region
    /// Output: Lx times Ly matrix with true values where invaded

    std::vector<std::vector<bool>> sites(Lx, std::vector<bool>(Ly));

    for(int i = 0; i < Lx; i++)
    {
        for(int j = 0; j < Ly; j++)
        {
            if (j < y0)
            {
                sites[i][j] = true;
            }
            else
            {
                sites[i][j] = false;
            }
        }
    }

    return sites;
}

std::vector<std::vector<double>> read_pinning_field(const std::string filename, const int Lx, const int Ly)
{
    /// Reading pinning stress field from file
    /// Lx: lattice width
    /// Ly: lattice height
    /// Output: Lx times Ly matrix of pinning stresses

    std::vector<std::vector<double>> pinning;//(Lx, std::vector<double>(Ly));

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening pinning stress file: " << filename << std::endl;
        return pinning;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);

        std::vector<double> stresses;
        double stress;
        while (iss >> stress) {
            stresses.push_back(stress);
	}

        pinning.push_back(stresses);
    }

    file.close();

    return pinning;
}

void change_external_stress(double* external_stress, const double external_stress_init, const double external_stress_rate, const int t_step, const double k_spring, const int plastic_strain)
{
    /// Updates external stress based on the applied external stess change and the stress change due to plastic deformation
    /// external_stress: the value of external stress to be updated
    /// external_stress_init: the value of the initial external stress (at the start of the simulation)
    /// external_stress_rate: external stress rate (stress increment per simulation step)
    /// t_step: number of time steps executed so far
    /// k_spring: stress drop per plastic strain unit
    /// plastic_strain: current value of plastic strain in grid cells
    
    (*external_stress) = external_stress_init + external_stress_rate * (double)t_step - k_spring * (double)plastic_strain;

    // external stress can be negative?
    /*
    if ( *external_stress < 0 )
    {
        (*external_stress) = 0.0;
    }
    */
}

struct boundaries find_boundaries(const std::vector<struct segment> segments, const int Ly)
{
    struct boundaries bounds;
    bounds.bottom = Ly-1;
    bounds.top = 0;

    int len = segments.size();
    for (int i=0; i<len; i++)
    {
        if (segments[i].y < bounds.bottom) {bounds.bottom = segments[i].y;}
        if (segments[i].y > bounds.top) {bounds.top = segments[i].y;}
    }

    return bounds;
}

struct boundaries find_area_to_scan(struct boundaries bounds, const int Ly)
{
    if (bounds.bottom < (Ly-1))
    {
        bounds.bottom = floor(bounds.bottom) - 3;
        if ( bounds.bottom < 0)
        {
            bounds.bottom = 0;
        }
    }

    if (bounds.top > 0)
    {
        bounds.top = ceil(bounds.top) + 3;
        if ( bounds.top > (Ly-1))
        {
            bounds.top = (Ly-1);
        }
    }


    /// Error handling
    if ( bounds.bottom < 0 || bounds.bottom > (Ly-1) )
    {
        std::cerr << "ERROR: invalid bottom bound." << std::endl;
    }

    if ( bounds.top < 0 || bounds.top > (Ly-1) )
    {
        std::cerr << "ERROR: invalid top bound." << std::endl;
    }


    return bounds;
}

int random_segment_index(const int n)
{
	/// Generates a random number to pick a segment
	/// n: number of segments
	/// Output: the random integer indexing a segment

	int k = std::rand();
	double j = (double)k / (double)RAND_MAX * (double)n;
	int i = (int)trunc(j);
	if (i == n)
	{
		i = n - 1;
	}
	return i;

}

