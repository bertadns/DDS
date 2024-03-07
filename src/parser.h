#ifndef PARSER_H
#define PARSER_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

#include "general.h"

void load_default_values(struct elastic *elastic, double* k_spring, int* Lx, int* Ly, int* y0, double* external_stress_init, double* external_stress_rate, int* seed, bool* is_edge_original);

int parser(const std::string filename, struct elastic *elastic, double* k_spring, int* Lx, int* Ly, int* y0, double* external_stress_init, double* external_stress_rate, int* seed, bool* is_edge_original);

void create_metafile(const struct elastic elastic, const double k_spring, const int Lx, const int Ly, const int y0, const double external_stress_init, const double external_stress_rate, const int seed, const std::string filename, const std::string pinning_filename, const std::string directory, const bool is_edge_original);

#endif // PARSER_H
