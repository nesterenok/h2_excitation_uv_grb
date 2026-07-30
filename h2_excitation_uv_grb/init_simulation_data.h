#pragma once
#include<cstring>
#include<vector>
#include "dynamic_array.h"

// sim_path    - path to the folder with simulation data (GRB propagation results)
// output_path - path to the folder where data of the current simulations must be saved,
void init_cloud_parameters(const std::string& sim_path, double& grb_cloud_distance, double& max_hcolumn_density, double& conc_h_tot,
	double& op_ratio_h2, double& h2_fraction, int& jet_type, double& half_opening_angle, double& theta_wing, double& viewing_angle);

void init_layer_coordinates(const std::string& path, std::vector<double>& layer_grid_coord, std::vector<double>& hcolumn_density_coord);

void init_specimen_conc(const std::string& sim_path, std::vector<dynamic_array>& d);

void init_h2_popdensity(const std::string& sim_path, std::vector<dynamic_array>& d, std::vector<int>& qnb_v_arr, std::vector<int>& qnb_j_arr);
