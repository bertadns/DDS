#ifndef GENERAL_H
#define GENERAL_H

/// Status flags
#define ONGOING 0
#define REACHED_BOTTOM 1
#define REACHED_TOP 2
#define STARVATION 3

/// Return flags
#define ARGC_ERROR 1
#define EMPTY_LATTICE 2
#define BAD_CONTROL_FILE 3

/// Number of image/ghost dislocations used on each sides
#define N_EWALD 2


//#ifndef SAVE_LATTICE
//#define SAVE_LATTICE
//#endif

#ifndef SAVE_SEGMENTS
#define SAVE_SEGMENTS
#endif

#ifndef SAVE_PINNING_FIELD
#define SAVE_PINNING_FIELD
#endif




#include <iostream>
#include <vector>

using namespace std;

struct elastic{
    /// Elastic and lattice parameters
    /// mu: shear modulus
    /// b: Burgers vector
    /// Poisson: Poisson's ratio
    /// a: lattice constant (of the grid of the simulation!)

    double mu;
    double b;
    double Poisson;
    double a;

};

struct segment{
    /// Strores features of an exisiting segment
    /// x, y: coordinates of the center of gravity of the dislocation segment
    /// is_edge: dislocation character (true: edge; false: screw)
    /// sign: sign of dislocation (positive if the line direction vector points to the positive direction)
    /// stress: the resolved shear stress acting on the dislocation

    double x;
    double y;
    bool is_edge;
    int sign;
    double stress;
};

struct adjacent_pinning{
    /// The pinning stress information on the top and right boundaries of the site
    /// top: pinning stress at the middle of the top boundary
    /// right: pinning stress at the middle of the right boundary

    double top;
    double right;
};


struct boundaries{
    /// Stores the boundaries of the dislocation lines
    /// top: highest vertical segment coordinate
    /// bottom: lowest vertical segment coordinate

    double top;
    double bottom;
};

std::vector<std::vector<bool>> create_lattice(const int Lx, const int Ly, const int y0);

std::vector<std::vector<double>> read_pinning_field(const std::string filename, const int Lx, const int Ly);

void change_external_stress(double* external_stress, const double external_stress_init, const double external_stress_rate, const int t_step, const double k_spring, const int plastic_strain);

struct boundaries find_boundaries(const std::vector<struct segment> segments, const int Ly);

struct boundaries find_area_to_scan(struct boundaries bounds, const int Ly);

int random_segment_index(const int n);

#endif // GENERAL_H
