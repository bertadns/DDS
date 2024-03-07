#include <cstdlib>
#include <vector>
#include <cmath>

#include "general.h"

std::vector<struct segment> find_segments(const int Lx, const int Ly, \
                                          std::vector<std::vector<bool>> lattice, \
                                          struct boundaries bounds)
{
    // Version dealing with initially screw-type line !!!!

    /// Find existing segments
    /// Lx: lattice width
    /// Ly: lattice height
    /// lattice: lattice containing invaded vs not invaded information
    /// bounds: struct containing the vertical bounds of the dislocation
    /// Output list of segment coordinate pairs

    std::vector<struct segment> segs;

    int imin = (int)(bounds.bottom);
    int imax = (int)(bounds.top);

    /// Standard scenario
    for (int x = 0; x < (Lx-1) ; x++)
    {
        for (int y = imin; y < imax; y++)
        {
            /// screw type
            if ( lattice[x][y] != lattice[x][y+1] )
            {
                if (lattice[x][y])
                {
                    segment segment_temp = {(double)x, (double)y+0.5, false, 1};
                    segs.push_back(segment_temp);
                } else
                {
                    segment segment_temp = {(double)x, (double)y+0.5, false, -1};
                    segs.push_back(segment_temp);
                }
            }
            /// edge type
            if ( lattice[x][y] != lattice[x+1][y] )
            {
                if (lattice[x][y])
                {
                    segment segment_temp = {(double)x+0.5, (double)y, true, -1};
                    segs.push_back(segment_temp);
                } else
                {
                    segment segment_temp = {(double)x+0.5, (double)y, true, 1};
                    segs.push_back(segment_temp);
                }
            }
        }
    }

    /// PBC scenario
    for (int y = imin; y < imax; y++)
    {
        /// screw type
        if ( lattice[Lx-1][y] != lattice[Lx-1][y+1] )
        {
            if (lattice[Lx-1][y])
            {
                segment segment_temp = {(double)(Lx-1), (double)y+0.5, false, 1};
                segs.push_back(segment_temp);
            } else
            {
                segment segment_temp = {(double)(Lx-1), (double)y+0.5, false, -1};
                segs.push_back(segment_temp);
            }
        }
        /// edge type
        if ( lattice[Lx-1][y] != lattice[0][y] )
        {
            if (lattice[Lx-1][y])
            {
                segment segment_temp = {(double)(Lx-1)+0.5, (double)y, true, -1};
                segs.push_back(segment_temp);
            } else
            {
                segment segment_temp = {(double)(Lx-1)+0.5, (double)y, true, 1};
                segs.push_back(segment_temp);
            }
        }
    }

    return segs;
}

std::vector<struct segment> find_segments_optimized(const int Lx, const int Ly, \
                                                    std::vector<std::vector<bool>> lattice,\
                                                    const std::vector<struct segment> segments0)
{
    struct boundaries bounds;
    bounds = find_boundaries(segments0, Ly);
    bounds = find_area_to_scan(bounds, Ly);

    std::vector<struct segment> segments;
    segments = find_segments(Lx, Ly, lattice, bounds);

    return segments;
}


double self_stress_screw(const double sign, const struct elastic elastic, const double dx, const double dy, const double factor_horizontal)
{
    /// Computes the self-stress induced by a screw segment
    /// sign: sign depending on the orientation of segments
    /// elastic: the elastic and lattice parameters
    /// dx, dy: relative coordinates of the segments
    /// Output: the signed value of self-stress

    double stress;
    stress = sign * elastic.mu * elastic.b * elastic.a * dy / ( 4 * M_PI * factor_horizontal * pow( dx*dx + dy*dy, 1.5 ) );

    return stress;
}

double self_stress_edge(const double sign, const struct elastic elastic, const double dx, const double dy, const double factor_vertical)
{
    /// Computes the self-stress induced by an edge segment
    /// sign: sign depending on the orientation of segments
    /// elastic: the elastic and lattice parameters
    /// dx, dy: relative coordinates of the segments
    /// Output: the signed value of self-stress

    double stress;
    stress = sign * elastic.mu * elastic.b * elastic.a * dx / ( 4 * M_PI * factor_vertical * pow( dx*dx + dy*dy, 1.5 ) );

    return stress;
}

double ewald_sum( double(*stress_func)(const double sign, const struct elastic elastic, const double dx, const double dy, const double factor), \
                    const double sign, const struct elastic elastic, const double dx, \
                    const double dy, const int Lx, const double xs, const double factor)
{
    /// Performes the Ewald-summation
    /// stress_func: the function pointer to the used self_stress function (screw or edge)
    /// sign: sign depending on the orientation of segments
    /// elastic: the elastic and lattice parameters
    /// dx, dy: relative coordinates of the segments
    /// Lx: the width of the lattice
    /// xs: the x coordinate of the source dislocation
    /// Output: the self-stress obtained by performing the Ewald-summation


    /// Native dislocation
    double stress = stress_func(sign, elastic, dx, dy, factor);

    /// Image dislocations
    for (int n=1; n<N_EWALD; n++)
    {
        stress += stress_func(sign, elastic, dx+Lx*n, dy, factor);
        stress += stress_func(sign, elastic, dx-Lx*n, dy, factor);
    }

    /// Ghost dislocations
    stress += stress_func(sign, elastic, dx+Lx*N_EWALD, dy, factor) * (Lx-xs)/(Lx-1);
    stress += stress_func(sign, elastic, dx-Lx*N_EWALD, dy, factor) * xs/(Lx-1);


    return stress;
}


int interaction_sign(const double sign_s, const double sign_t, const bool is_edge_s, const bool is_edge_t)
{
    /// Computes the orientation-dependent sign of the interaction between two segments
    /// sign_s, sign_t: signs of the source and the target segments, respectively
    /// is_edge_s, is_edge_t: character of the source and the target segments, respectively (edge=true, screw=false)
    /// Output: interaction sign

    int sign;

    /// Computing the sign with respect to directions of axes x, y
    if ( is_edge_s == is_edge_t )
    {
        sign = sign_s * sign_t;
    } else
    {
        sign = -sign_s * sign_t;
    }

    return sign;
}

double get_mean_pinning_stress(const std::vector<struct segment> segs, const int i, const std::vector<std::vector<double>> pinning, \
                          const std::vector<std::vector<bool>> lattice, const int Lx)
{
    /// Finds the pinning stress at the segment location (at the site it would invade)
    /// segs: vector containg the segments
    /// i: the # of the segment at which we compute the resolved shear stress
    /// pinning: matrix containing pinning stress data
    /// lattice: pointer of the lattice
    /// Lx: lattice width
    /// Output: pinning stress value

    double stress;

    /// Segment coordinates
    double x = segs[i].x;
    double y = segs[i].y;

    /// Site coordinates
    double x0, y0;


    if (x == (int)x)
    /// Horizontal segment
    {
        x0 = x;
        y0 = trunc(y);

        stress = ( pinning[x0][y0] + pinning[x0][(int)(y0+1)] ) / 2.0;
    }
    else
    /// Vertical segment
    {
        x0 = trunc(x);
        y0 = y;

        stress = ( pinning[x0][y0] + pinning[(int)((int)(x0+1) % Lx)][y0] ) / 2.0;
    }

    return stress;
}

double PBC(double x, const int Lx)
{
    if (x > (Lx - 1) ) // Lx or Lx-1 ???
    {x = x - Lx;}
    else
    {
        if (x < 0){x = x + Lx;}
    }

    return x;
}



double get_self_stress(const std::vector<struct segment> segs, const int i, \
                       const struct elastic elastic, const int Lx,\
                       const std::vector<std::vector<bool>> lattice, const double factor_horizontal, const double factor_vertical)
{
    /// Performs the sum of self stresses
    /// segs: vector containg the segments
    /// i: the # of the segment at which we compute the resolved shear stress
    /// elastic: the elastic and lattice parameters
    /// Output: the total value of self stress

    double stress = 0.0;

    /// target attributes
    double x = segs[i].x;
    double y = segs[i].y;
    bool is_edge_t = segs[i].is_edge;
    int sign_t = segs[i].sign;

    /// source attributes
    double xs, ys;
    bool is_edge_s;
    int sign_s;

    /// relative attributes
    double dx, dy;

    /// Loop throuh all other pre-existent segments
    int len = segs.size();
    for (int j = 0; j < len; j++)
    {
        if (i != j)
        {
            xs = segs[j].x;
            ys = segs[j].y;
            is_edge_s = segs[j].is_edge;
            sign_s = segs[j].sign;

            dx = x - xs;
            dy = y - ys;

            int sign = interaction_sign(sign_s, sign_t, is_edge_s, is_edge_t);

            if (is_edge_s)
            {
                stress += ewald_sum(self_stress_edge, sign, elastic, dx, dy, Lx, xs, factor_vertical);
            }
            else
            {
                stress += ewald_sum(self_stress_screw, sign, elastic, dx, dy, Lx, xs, factor_horizontal);
            }
        }
    }

    return stress;
}

double get_total_stress(const std::vector<struct segment> segs, const int i, \
                        const std::vector<std::vector<double>> pinning, \
                        const struct elastic elastic, const double external_stress, const int Lx, \
                        const std::vector<std::vector<bool>> lattice, const double factor_horizontal, const double factor_vertical)
{
    /// Sums all stress types for one segment
    /// segs: vector containg the segments
    /// i: the # of the segment at which we compute the resolved shear stress
    /// pinning: matrix containing pinning stress data
    /// elastic: the elastic and lattice parameters
    /// external_stress: value of external stress
    /// lattice: pointer of the lattice
    /// Output: sum of stresses at the segment

    /// PINNING STRESS CONTRIBUTION
    double pinning_stress = get_mean_pinning_stress(segs, i, pinning, lattice, Lx);

    /// SELF STRESS CONTRIBUTION
    double self_stress = get_self_stress(segs, i, elastic, Lx, lattice, factor_horizontal, factor_vertical);

    /// Computing the sign with respect to the extremal direction
    bool is_edge = segs[i].is_edge;
    int sign = segs[i].sign;
    int sign_ext;

    if ( ( (!is_edge) && (sign == 1) ) || ( (is_edge) && (sign==-1) ) )
    {
        sign_ext = 1;
    }
    else
    {
        sign_ext = -1;
    }

    double stress = external_stress + pinning_stress + self_stress * sign_ext;

    return stress;
}

bool  stress_to_spin(const double stress)
{
	/// Converts local total stress value into a flip
	/// stress: total stress
	/// Output: new spin value on the two sides of the segment (true: positive plastic strain, true: negative plastic strain)

	if (stress < 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

int time_step(std::vector<std::vector<bool>>* lattice, int* plastic_strain, const int Lx, const int Ly, \
              const std::vector<struct segment> segs, const int index, struct boundaries bounds, const bool new_spin)
{
    /// Updates lattice and plastic strain and returns information about current system state
    /// lattice: pointer of the lattice
    /// plastic_strain: plastic strain measured in gridcell
    /// Lx, Ly: width and height of lattice
    /// segs: vector of segment information
    /// index: index of segment with extremal stress
    /// bounds: struct containing the vertical bounds of the dislocation
    /// new_spin: new_spin to be flipped to (true: positive plastic strain, false: negative plastic strain)
    /// Output: code of the current state of the system

    /// Stopping scenarios
    if (segs[index].y > (Ly-2) )
    {
        /// The dislocation line reached the limit at the top of the lattice
        return REACHED_TOP;
    }

    if (segs[index].y < 1)
    {
        /// The dislocation line reached the limit at the bottom of the lattice
        return REACHED_BOTTOM;
    }

    /// Moving scenario
    if (segs[index].is_edge)
    /// moving edge segment
    {

        int x1 = (int)(segs[index].x-0.5);
        int x2;
        int y = (int)(segs[index].y);

        if (segs[index].x<(Lx-1))
        {
            /// Standard scenario
            x2 = (int)(segs[index].x+0.5);
        }
        else
        {
            /// PBC scenario
            x2 = 0;
        }
        (*lattice)[x1][y] = new_spin;
        (*lattice)[x2][y] = new_spin;
    }
    else
    /// moving screw segment
    {
        int x = (int)(segs[index].x);
        int y1 = (int)(segs[index].y-0.5);
        int y2 = (int)(segs[index].y+0.5);
        (*lattice)[x][y1] = new_spin;
        (*lattice)[x][y2] = new_spin;
    }
 
    /// Updating plastic strain
    if (new_spin)
    	(*plastic_strain)++;
    else
    	(*plastic_strain)--;


    std::vector<struct segment> segs_new = find_segments(Lx, Ly, *lattice, bounds);

    if (segs_new.size() == 0)
    {
        /// All segments are annilihated
        return STARVATION;
    }


    return ONGOING;
}
