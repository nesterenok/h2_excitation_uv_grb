#pragma once

#include <cstring>
#include <vector>
#include "dynamic_array.h"
#include "spectroscopy.h"


// H2 population densities [cm-3] as a function of time,
void save_h2_populations_evolution(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& h2_popdens_data, double conc_h_tot, double hcolumn_density, int nb_lev_h2);

// only rotational levels of H2, J =< J_max 
void save_h2_populations_evolution_rotational(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& h2_popdens_data, double conc_h_tot, double hcolumn_density, int nb_lev_h2);

//
// line_luminosity_evol_arr - (first) vector index: transition number, (second) array index - time moment number,
void save_h2_line_luminosities(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<transition> & transition_list, 
    const std::vector<dynamic_array>& luminosity_evol_arr, const std::vector<dynamic_array> & luminosity_evol_arr_estimated, 
    double theta_jet_axis, bool is_jet_top_hat, bool rotational_transition_mode);
