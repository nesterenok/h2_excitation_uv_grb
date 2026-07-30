#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include<iostream>
#include<fstream>
#include<iomanip>

#include"h2_population_evol_eqs.h"
#include"constants.h"
#include"spectroscopy.h"
#include"parameters.h"
#include"interpolation.h"
#include"integration.h"
#include"utils.h"

using namespace std;
#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "h2_population_evol_eqs.cpp"


h2_population_evolution_data::h2_population_evolution_data(const string& data_path, const string& output_path, int verb) :
	verbosity(verb)
{
	int isotope;

	// new H2 molecule data (Roueff et al. A&A 630, A58, 2019), max number of levels is 302,
	nb_lev_h2 = 302;
	molecule h2_mol("H2", isotope = 1, 2. * ATOMIC_MASS_UNIT);

	h2_di = new h2_diagram_roueff2019(data_path, h2_mol, nb_lev_h2, verbosity);
	h2_einst = new h2_einstein_coeff_roueff2019(data_path, h2_di, verbosity);

	// saving H2 levels of the ground electronic state that are used in simulations,
	save_h2_levels_used(output_path, *h2_di);
}

h2_population_evolution_data::~h2_population_evolution_data()
{
	delete h2_di;
	delete h2_einst;
}

int h2_population_evolution_data::get_level_nb_h2(int v, int j) const
{
	int nb;
	nb = h2_di->get_nb(v, j);
	return nb;
}

int h2_population_evolution_data::get_vibr_nb_h2(int nb) const
{
	if (nb < 0 || nb >= nb_lev_h2)
		return -1;
	return h2_di->lev_array[nb].v;
}

double h2_population_evolution_data::get_energy_h2_rovibr(N_Vector y)
{
	double energy_in_rovibr_levels = 0.;
	realtype* y_data = NV_DATA_S(y);

	// level energy in [cm-1], population density in [cm-3]
	for (int i = 2; i < nb_lev_h2; i++) {
		if (rounding(h2_di->lev_array[i].spin) == 0) {
			energy_in_rovibr_levels += h2_di->lev_array[i].energy * y_data[i];
		}
		else {
			energy_in_rovibr_levels += (h2_di->lev_array[i].energy - h2_di->lev_array[1].energy) * y_data[i];
		}
	}
	return energy_in_rovibr_levels * CM_INVERSE_TO_EV;  // in [eV cm-3]
}


int h2_population_evolution_data::f(realtype t, N_Vector y, N_Vector ydot)
{
	int i, j;
	double x;
	realtype* y_data, * ydot_data;

	y_data = NV_DATA_S(y);
	ydot_data = NV_DATA_S(ydot);

	for (i = 0; i < nb_lev_h2; i++) {
		ydot_data[i] = 0.;
	}

	// Radiative transitions in H2 (only radiative decay)
	for (i = 1; i < nb_lev_h2; i++) {
		for (j = 0; j < i; j++) {
			if (h2_einst->arr[i][j] > 1.e-99)
			{
				x = h2_einst->arr[i][j] * y_data[i];
				ydot_data[j] += x;
				ydot_data[i] -= x;
			}
		}
	}	
	return 0;
}

void h2_population_evolution_data::set_tolerances(N_Vector abs_tol)
{
	int i;
	for (i = 0; i < nb_lev_h2; i++) {
		NV_Ith_S(abs_tol, i) = ABS_POPULATION_ERROR_SOLVER;
	}
}


// Saving the quantum numbers and energies of H2 levels that are used in simulations,
void save_h2_levels_used(const std::string& output_path, const energy_diagram& h2_di)
{
	int i, nb;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_h2_levels.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << "!Energy levels of the ground state of H2 used in simulations;" << endl
		<< "!data: nb, quantum numbers v and j, level energy in [cm-1]" << endl;

	nb = (int)h2_di.lev_array.size();
	output << nb << endl;

	for (i = 0; i < nb; i++) {
		output << left << setw(5) << i + 1 << setw(5) << h2_di.lev_array[i].v << setw(5) << rounding(h2_di.lev_array[i].j)
			<< setw(14) << h2_di.lev_array[i].energy;

		if (i < nb - 1) {
			output << endl;
		}
	}
	output.close();
}


double func::operator()(double phi) const {
	return acos(a + b * cos(phi));
}


table_integrals::table_integrals()
{
	int i;
	nb_of_phi = 501; 
	nb_of_theta = 501;

	max_phi = M_PI;
	min_theta = 0.;
	max_theta = M_PI;

	min_cos_theta = cos(max_theta);
	max_cos_theta = cos(min_theta);
	
	dphi = max_phi / (nb_of_phi - 1);
	dtheta = (max_theta - min_theta)/ (nb_of_theta - 1);

	phi_arr = new double[nb_of_phi];
	cos_theta_arr = new double[nb_of_theta];
	sin_theta_arr = new double[nb_of_theta];

	for (i = 0; i < nb_of_phi; i++) {
		phi_arr[i] = i * dphi;
	}

	for (i = 0; i < nb_of_theta; i++){
		cos_theta_arr[i] = cos(max_theta - i * dtheta);
		sin_theta_arr[i] = sin(max_theta - i * dtheta);
	}

	integral_arr = alloc_2d_array<double>(nb_of_theta, nb_of_phi);
	memset(*integral_arr, 0, nb_of_theta * nb_of_phi * sizeof(integral_arr[0][0]));
}

void table_integrals::setp(double cos_theta_jet_axis, double sin_theta_jet_axis)
{
	int i, j;
	double eps = 1.e-6;

	for (i = 0; i < nb_of_theta; i++) {
		f.setp(cos_theta_jet_axis * cos_theta_arr[i], sin_theta_jet_axis* sin_theta_arr[i]);
		
		for (j = 0; j < nb_of_phi; j++) {
			integral_arr[i][j] = qromb<func>(f, 0., phi_arr[j], eps);
		}
	}
}

table_integrals::~table_integrals() {
	delete[] phi_arr;
	delete[] cos_theta_arr;
	delete[] sin_theta_arr;
	free_2d_array(integral_arr);
}

double table_integrals::get(double phi, double cos_theta) const
{
	int i, j;
	double t, u, a(0.);

	if (cos_theta <= min_cos_theta) {
		i = 0;
		t = 0.;
	}
	else if (cos_theta >= max_cos_theta) {
		i = nb_of_theta - 2;
		t = 1.;
	}
	else {
		locate_index(cos_theta_arr, nb_of_theta, cos_theta, i);
		t = (cos_theta - cos_theta_arr[i]) / (cos_theta_arr[i + 1] - cos_theta_arr[i]);
	}

	if (phi <= 0.) {
		j = 0;
		u = 0.;
	}
	else if (phi >= max_phi) {
		j = nb_of_phi - 2;
		u = 1.;
	}
	else {
		j = (int)(phi / dphi);
		u = (phi - phi_arr[j]) / (phi_arr[j + 1] - phi_arr[j]);
	}

	a = integral_arr[i][j] * (1. - t) * (1. - u) + integral_arr[i][j + 1] * (1. - t) * u
		+ integral_arr[i + 1][j] * t * (1. - u) + integral_arr[i + 1][j + 1] * u * t;
	
	return a;
}
