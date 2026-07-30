#define _CRT_SECURE_NO_WARNINGS

#include<omp.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<sstream>
#include<ctime>
#include<limits>
#include<sstream>
#include<limits>

#include "constants.h"
#include "saving_data.h"
#include "parameters.h"
using namespace std;

#define SAVE_DATA_MAX_J 9   // including this value, J =< J_max
const int h2_ground_level_nb[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 17, 19};


void save_h2_populations_evolution(const string& output_path, const vector<double>& time_moments,
	const vector<dynamic_array>& h2_popdens_data, double conc_h_tot, double hcolumn_density, int nb_lev_h2)
{
	int i, l, nb_of_times;
	double w;

	string fname;
	ofstream output;
	stringstream ss;

	for (i = 0; i < (int)time_moments.size(); i++) {
		if (time_moments[i] > MAX_SAVING_TIME) {
			break;
		}
	}
	nb_of_times = i;

	// suffix of the output file for saving the data on the specific H column density,
	// Here, H column density must be = integer value * 1.e+18
	i = 18;
	l = rounding(hcolumn_density / pow(10., i));

	ss.clear();
	ss.str("");
	ss << "_nh" << l << "e" << i;

	// Saving level populations for all ro-vibrational energy levels,
	fname = output_path + "h2grb_popdens";
	fname += ss.str();
	fname += ".txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Evolution of population densities [cm-3] of H2 energy levels," << endl
		<< "!first row: H nuclei concentration [cm-3], nb of times, number of levels," << endl
		<< "!second row: level nb starting from 1," << endl;

	output << left << "! " << setw(15) << conc_h_tot << setw(7) << nb_of_times << setw(7) << nb_lev_h2 << endl;

	output << left << setw(13) << "!time(s)";
	for (i = 1; i <= nb_lev_h2; i++) {
		output << left << setw(10) << i;
	}
	output << endl;

	for (l = 0; l < nb_of_times; l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		for (i = 0; i < nb_lev_h2; i++) {
			output.precision(3);
			w = h2_popdens_data[l].arr[i];

			if (w < MINIMAL_ABUNDANCE) {
				w = MINIMAL_ABUNDANCE;
			}
			if (w < 1.e-11) {
				output.precision(1);
			}
			output << left << w << " ";
		}
		if (l < nb_of_times - 1) {
			output << endl;
		}
	}
	output.close();
}

//
// Saving level populations for rotational energy levels (0, J), J =< J_max;
void save_h2_populations_evolution_rotational(const string & output_path, const vector<double>&time_moments,
	const vector<dynamic_array>&h2_popdens_data, double conc_h_tot, double hcolumn_density, int nb_lev_h2)
{
	int i, l, nb_of_times;
	double w;

	string fname;
	ofstream output;
	stringstream ss;

	nb_of_times = (int)time_moments.size();

	// suffix of the output file for saving the data on the specific H column density,
	// Here, H column density must be = integer value * 1.e+18
	i = 18;
	l = rounding(hcolumn_density / pow(10., i));

	ss.clear();
	ss.str("");
	ss << "_nh" << l << "e" << i;

	fname = output_path + "rot_h2grb_popdens";
	fname += ss.str();
	fname += ".txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Evolution of population densities [cm-3] of H2 rotational energy levels (0,J)," << endl
		<< "!first row: H nuclei concentration [cm-3], nb of times, number of levels," << endl
		<< "!second row: J of the energy level nb starting from 0," << endl;

	output << left << "! " << setw(15) << conc_h_tot << setw(7) << nb_of_times << setw(7) << h2_ground_level_nb[SAVE_DATA_MAX_J] << endl;

	output << left << setw(13) << "!time(s)";
	for (i = 1; i <= h2_ground_level_nb[SAVE_DATA_MAX_J]; i++) {
		output << left << setw(10) << i;
	}
	output << endl;

	for (l = 0; l < nb_of_times; l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		for (i = 0; i < h2_ground_level_nb[SAVE_DATA_MAX_J]; i++) {
			output.precision(3);
			w = h2_popdens_data[l].arr[i];

			if (w < MINIMAL_ABUNDANCE) {
				w = MINIMAL_ABUNDANCE;
			}
			if (w < 1.e-11) {  // some arbitrary small value
				output.precision(1);
			}
			output << left << w << " ";	
		}
		if (l < nb_of_times - 1) {
			output << endl;
		}
	}
	output.close();
}


void save_h2_line_luminosities(const std::string& output_path, const std::vector<double>& time_moments,
	const vector<transition>& transition_list, 
	const vector<dynamic_array>& luminosity_evol_arr, const vector<dynamic_array>& luminosity_evol_arr_estimated,
	double theta_jet_axis, bool is_jet_top_hat, bool rotational_transition_mode)
{
	int i, j, k, nb_of_lines, nb_of_times;
	double x;

	string fname, str_id;
	ofstream output;
	stringstream ss;

	k = rounding(theta_jet_axis * 100);
	ss.clear();
	ss.str("");

	// jet type
	if (is_jet_top_hat) {
		ss << "th";  // top-hat
	}
	else {
		ss << "sj";  // structural jet
	}

	// viewing angle
	if (k > 0 && k < 10) {
		ss << "00";
	}
	else if (k >= 10 && k < 100) {
		ss << "0";
	}
	ss << k;

	str_id = ss.str();

	nb_of_lines = (int)transition_list.size();
	nb_of_times = (int)time_moments.size();
	
	if (rotational_transition_mode) {
		fname = output_path + "rot_h2_line_list_";
	}
	else {
		fname = output_path + "h2_line_list_";
	}
	fname += str_id;
	fname += ".txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);
	output << left << "!List of H2 line, sorted by wavelength, (v',j') -> (v'', j''); nb is the number of the line in the luminosity data file;" << endl;

	output << left << setw(6) << "!nb" << setw(5) << "v' " << setw(5) << "j' " 
		<< setw(5) << "v'' " << setw(5) << "j'' " << setw(14) << "lambda(um)" << endl;
	
	output << nb_of_lines << endl;

	for (i = 0; i < nb_of_lines; i++) {
		output << left << setw(6) << i << setw(5) << transition_list[i].up_lev.v << setw(5) << (int)transition_list[i].up_lev.j
			<< setw(5) << transition_list[i].low_lev.v << setw(5) << (int)transition_list[i].low_lev.j 
			<< setw(14) << 1.e+4 * SPEED_OF_LIGHT / transition_list[i].freq;
		
		if (i < nb_of_lines - 1) {
			output << endl;
		}
	}
	output.close();

	if (rotational_transition_mode) {
		fname = output_path + "rot_h2_line_luminosities_";
	}
	else {
		fname = output_path + "h2_line_luminosities_";
	}
	fname += str_id;
	fname += ".txt";

	output.open(fname.c_str());
	output << scientific;

	output << left << "!Luminosities of H2 lines," << endl;

	output << left << setw(12) << "!time(s)";
	for (i = 0; i < nb_of_lines; i++) {
		output << left << setw(12) << i;
	}
	output << endl;

	for (j = 0; j < nb_of_times; j++) {
		output << left << setprecision(4) << setw(12) << time_moments[j];

		for (i = 0; i < nb_of_lines; i++) {
			// vector index - transition number
			output << left << setprecision(4) << setw(12) << luminosity_evol_arr[i].arr[j];
		}
		if (j < nb_of_times - 1) {
			output << endl;
		}
	}
	output.close();

	//
	// Is used for analisys,
	// Sorting lines by luminosity,
	vector<double> lum_arr, lum_arr_aux;
	vector<int> nb_arr;

	for (i = 0; i < nb_of_lines; i++) {
		x = 0.;
		for (j = 0; j < nb_of_times; j++) {
			if (x < luminosity_evol_arr[i].arr[j]) {
				x = luminosity_evol_arr[i].arr[j];
			}
		}
		lum_arr.push_back(x);
		nb_arr.push_back(i);

		x = 0.;
		for (j = 0; j < nb_of_times; j++) {
			if (x < luminosity_evol_arr_estimated[i].arr[j]) {
				x = luminosity_evol_arr_estimated[i].arr[j];
			}
		}
		lum_arr_aux.push_back(x);
	}

	for (i = 0; i < nb_of_lines; i++) {
		for (j = i + 1; j < nb_of_lines; j++)
		{
			if (lum_arr[j] > lum_arr[i]) {
				x = lum_arr[i];
				lum_arr[i] = lum_arr[j];
				lum_arr[j] = x;

				x = lum_arr_aux[i];
				lum_arr_aux[i] = lum_arr_aux[j];
				lum_arr_aux[j] = x;

				k = nb_arr[i];
				nb_arr[i] = nb_arr[j];
				nb_arr[j] = k;
			}
		}
	}

	if (rotational_transition_mode) {
		fname = output_path + "rot_h2_line_list_";
	}
	else {
		fname = output_path + "h2_line_list_";
	}
	fname += str_id;
	fname += "_aux.txt";

	output.open(fname.c_str());
	output << scientific;

	output.precision(6);
	output << left << "!List of H2 line, sorted by luminosity, (v',j') -> (v'', j'');" << endl 
		<< "!nb is the number of the line in the luminosity data file;" << endl 
		<< "!maximum luminosity based on accurate and not- accurate dust opacity calculations;" << endl
		<< "!last parameter: -1 means the transition is absent in the initial list of transitions;" << endl;

	output << left << setw(6) << "!nb" << setw(5) << "v' " << setw(5) << "j' "
		<< setw(5) << "v'' " << setw(5) << "j'' " << setw(14) << "lambda(um)" << setw(14) << "max_lum(acc)" << setw(14) << "max_lum(est)" << endl;

	output << "!para-H2" << endl;
	for (i = 0; i < nb_of_lines; i++) {
		if ((int)transition_list[nb_arr[i]].up_lev.j % 2 == 0)
		{
			k = 1;
			if (fabs(lum_arr[i] - lum_arr_aux[i]) < 10. * numeric_limits<double>::epsilon()) {
				k = -1;
			}
			output << left << setw(6) << nb_arr[i]
				<< setw(5) << transition_list[nb_arr[i]].up_lev.v << setw(5) << (int)transition_list[nb_arr[i]].up_lev.j
				<< setw(5) << transition_list[nb_arr[i]].low_lev.v << setw(5) << (int)transition_list[nb_arr[i]].low_lev.j
				<< setw(14) << 1.e+4 * SPEED_OF_LIGHT / transition_list[nb_arr[i]].freq << setw(14) << lum_arr[i] << setw(14) << lum_arr_aux[i]
				<< k << endl;
		}
	}

	output << "!ortho-H2" << endl;
	for (i = 0; i < nb_of_lines; i++) {
		if ((int)transition_list[nb_arr[i]].up_lev.j % 2 == 1)
		{
			k = 1;
			if (fabs(lum_arr[i] - lum_arr_aux[i]) < 10. * numeric_limits<double>::epsilon()) {
				k = -1;
			}
			output << left << setw(6) << nb_arr[i]
				<< setw(5) << transition_list[nb_arr[i]].up_lev.v << setw(5) << (int)transition_list[nb_arr[i]].up_lev.j
				<< setw(5) << transition_list[nb_arr[i]].low_lev.v << setw(5) << (int)transition_list[nb_arr[i]].low_lev.j
				<< setw(14) << 1.e+4 * SPEED_OF_LIGHT / transition_list[nb_arr[i]].freq << setw(14) << lum_arr[i] << setw(14) << lum_arr_aux[i]
				<< k << endl;
		}
	}
	output.close();
}
