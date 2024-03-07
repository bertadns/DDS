#ifndef SCREW_H
#define SCREW_H

#include "general.h"
#include <cstdlib>
#include <vector>

std::vector<struct segment> find_segments(const int Lx, const int Ly, \
                                          std::vector<std::vector<bool>> lattice, \
                                          struct boundaries bounds);

std::vector<struct segment> find_segments_optimized(const int Lx, const int Ly, \
                                                    std::vector<std::vector<bool>> lattice,\
                                                    const std::vector<struct segment> segments0);


double self_stress_screw(const double sign, const struct elastic elastic, const double dx, \
                         const double dy, const double factor_horizontal);


double self_stress_edge(const double sign, const struct elastic elastic, const double dx, \
                        const double dy, const double factor_vertical);

double ewald_sum( double(*stress_func)(const double sign, const struct elastic elastic, \
                  const double dx, const double dy, const double factor), const double sign, const struct elastic elastic, \
                  const double dx, const double dy, const int Lx, const double xs, const double factor);


int interaction_sign(const double sign_s, const double sign_t, const bool is_edge_s, \
                     const bool is_edge_t);


double get_mean_pinning_stress(const std::vector<struct segment> segs, const int i, \
                          const std::vector<std::vector<double>> pinning, \
                          const std::vector<std::vector<bool>> lattice);

double PBC(double x, const int Lx);

double get_self_stress(const std::vector<struct segment> segs, const int i, \
                       const struct elastic elastic, const int Lx, \
                       const std::vector<std::vector<bool>> lattice, const double factor_horizontal, const double factor_vertical);

double get_total_stress(const std::vector<struct segment> segs, const int i, \
                        const std::vector<std::vector<double>> pinning, \
                        const struct elastic elastic, const double external_stress, const int Lx, \
                        const std::vector<std::vector<bool>> lattice, const double factor_horizontal, const double factor_vertical);

bool  stress_to_spin(const double stress);

int time_step(std::vector<std::vector<bool>>* lattice, int* plastic_strain, const int Lx, const int Ly, \
              const std::vector<struct segment> segs, const int index, struct boundaries bounds, const bool new_spin);


#endif // SCREW_H
