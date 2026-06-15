#pragma once

// SUNDIALS CVODE solver headers
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */

#include "spectroscopy.h"


class h2_population_evolution_data {
protected:
    int verbosity, nb_lev_h2;

    const energy_diagram* h2_di;
    const einstein_coeff* h2_einst;

public:
    // function defining the ODE system for the evolution of the parameters
    int f(realtype t, N_Vector y, N_Vector ydot);
    void set_tolerances(N_Vector abs_tol);

    int get_nb_of_h2_lev() const { return nb_lev_h2; }
    int get_level_nb_h2(int v, int j) const;  // returns the nb of the H2 level belonging to the ground electronic state, 
    int get_vibr_nb_h2(int level_nb) const;   // returns the vibration quantum nb of the ground electronic state of H2,

     // returns energy in [eV cm-3] locked in ro-vibrational data:
    double get_energy_h2_rovibr(N_Vector y);

    h2_population_evolution_data(const std::string& data_path, const std::string& output_path, int verb);
    ~h2_population_evolution_data();
};


// saving H2 level quantum numbers and energies of the ground electronic state, for which data are saved,
void save_h2_levels_used(const std::string& output_path, const energy_diagram&);



class func {
protected:
    double a, b;

public:
    void setp(double aa, double bb) { a = aa; b = bb; }
    double operator() (double phi) const;
};


class table_integrals {
protected:
    int    nb_of_phi, nb_of_theta;
    double max_phi, min_theta, max_theta, min_cos_theta, max_cos_theta, dphi, dtheta;
    
    double* phi_arr, * cos_theta_arr, * sin_theta_arr, ** integral_arr;  // 0 <= phi <= pi, 0 < theta < pi,
    func f;

public:
    double get(double phi, double cos_theta) const;
    void setp(double cos_theta_jet_axis, double sin_theta_jet_axis);
    
    table_integrals();
    ~table_integrals();
};
