#define _CRT_SECURE_NO_WARNINGS

#include <omp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>
#include <ctime>
#include <limits>
#include <sstream>

#include "constants.h"
#include "h2_saving_data_uv_grb.h"
#include "h2_parameters_uv_grb.h"
using namespace std;


void save_h2_populations_evolution(const string& output_path, const vector<double>& time_moments,
	const vector<dynamic_array>& h2_popdens_data, double conc_h_tot, double hcolumn_density, int nb_lev_h2)
{
	int i, l;
	double w;
	
	string fname;
	ofstream output;
	stringstream ss;

	// suffix of the output file for saving the data on the specific H column density,
	// Here, H column density must be = integer value * 1.e+18
	i = 18;
	l = rounding(hcolumn_density / pow(10., i));

	ss.clear();
	ss.str("");
	ss << "_nh" << l << "e" << i;

	fname = output_path + "h2grb_popdens";
	fname += ss.str();
	fname += ".txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Evolution of population densities [cm-3] of H2 energy levels," << endl
		<< "!first row: H nuclei concentration [cm-3], nb of times, number of levels," << endl
		<< "!second row: level nb starting from 1," << endl;

	output << left << "! " << setw(15) << conc_h_tot << setw(7) << (int)h2_popdens_data.size() << setw(7) << nb_lev_h2 << endl;

	output << left << setw(13) << "!time(s)";
	for (i = 1; i <= nb_lev_h2; i++) {
		output << left << setw(10) << i;
	}
	output << endl;

	for (l = 0; l < (int)h2_popdens_data.size(); l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		for (i = 0; i < nb_lev_h2; i++) {
			output.precision(3);
			w = h2_popdens_data[l].arr[i];

			if (w < MINIMAL_ABUNDANCE) {
				w = MINIMAL_ABUNDANCE;
			}
			if (w < ABS_POPULATION_ERROR_SOLVER) {
				output.precision(1);
			}
			output << left << w << " ";
		}
		output << endl;
	}
	output.close();
}


void save_h2_line_luminosities(const std::string& output_path, const std::vector<double>& time_moments,
	const std::vector<transition>& trans_list, const std::vector<dynamic_array>& line_luminosity_evol_arr, 
	double theta_jet_axis, bool is_jet_top_hat)
{
	int i, j, k, nb_of_lines, nb_of_times;
	double x;

	string fname, str_id;
	ofstream output;
	stringstream ss;

	k = rounding(theta_jet_axis * 100);
	ss.clear();
	ss.str("");

	if (is_jet_top_hat) {
		ss << "th";  // top-hat
	}
	else {
		ss << "sj";  // structural jet
	}

	if (k > 0 && k < 10) {
		ss << "00";
	}
	else if (k >= 10 && k < 100) {
		ss << "0";
	}
	ss << k;

	str_id = ss.str();

	nb_of_lines = (int)trans_list.size();
	nb_of_times = (int)time_moments.size();
	
	fname = output_path + "h2_line_list_";
	fname += str_id;
	fname += ".txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);
	output << left << "!List of H2 line, sorted by wavelength, (v',j') -> (v'', j'')" << endl;

	output << left << setw(6) << "!nb" << setw(5) << "v' " << setw(5) << "j' " 
		<< setw(5) << "v'' " << setw(5) << "j'' " << setw(14) << "lambda(um)" << endl;
	
	output << nb_of_lines << endl;

	for (i = 0; i < nb_of_lines; i++) {
		output << left << setw(6) << i << setw(5) << trans_list[i].up_lev.v << setw(5) << (int)trans_list[i].up_lev.j
			<< setw(5) << trans_list[i].low_lev.v << setw(5) << (int)trans_list[i].low_lev.j 
			<< setw(14) << 1.e+4 * SPEED_OF_LIGHT / trans_list[i].freq;
		
		if (i < nb_of_lines - 1) {
			output << endl;
		}
	}
	output.close();

	fname = output_path + "h2_line_luminosities_";
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
			output << left << setprecision(4) << setw(12) << line_luminosity_evol_arr[i].arr[j];
		}
		if (j < nb_of_times - 1) {
			output << endl;
		}
	}
	output.close();

	// Sorting lines by luminosity,
	vector<double> lum_arr;
	vector<int> nb_arr;

	for (i = 0; i < nb_of_lines; i++) {
		x = 0.;
		for (j = 0; j < nb_of_times; j++) {
			if (x < line_luminosity_evol_arr[i].arr[j]) {
				x = line_luminosity_evol_arr[i].arr[j];
			}
		}
		lum_arr.push_back(x);
		nb_arr.push_back(i);
	}

	for (i = 0; i < nb_of_lines; i++) {
		for (j = i + 1; j < nb_of_lines; j++) 
		{
			if (lum_arr[j] > lum_arr[i]) {
				x = lum_arr[i];
				lum_arr[i] = lum_arr[j];
				lum_arr[j] = x;

				k = nb_arr[i];
				nb_arr[i] = nb_arr[j];
				nb_arr[j] = k;
			}
		}
	}

	fname = output_path + "h2_line_list_sorted_";
	fname += str_id;
	fname += ".txt";

	output.open(fname.c_str());
	output << scientific;

	output.precision(6);
	output << left << "!List of H2 line, sorted by luminosity, (v',j') -> (v'', j'')" << endl;

	output << left << setw(6) << "!nb" << setw(5) << "v' " << setw(5) << "j' "
		<< setw(5) << "v'' " << setw(5) << "j'' " << setw(14) << "lambda(um)" << setw(14) << "max_lum(solar)" << endl;

	output << "!para-H2" << endl;
	for (i = 0; i < nb_of_lines; i++) {
		if ((int)trans_list[nb_arr[i]].up_lev.j % 2 == 0)
		{
			output << left << setw(6) << i
				<< setw(5) << trans_list[nb_arr[i]].up_lev.v << setw(5) << (int)trans_list[nb_arr[i]].up_lev.j
				<< setw(5) << trans_list[nb_arr[i]].low_lev.v << setw(5) << (int)trans_list[nb_arr[i]].low_lev.j
				<< setw(14) << 1.e+4 * SPEED_OF_LIGHT / trans_list[nb_arr[i]].freq << setw(14) << lum_arr[i] 
				<< endl;
		}
	}

	output << "!ortho-H2" << endl;
	for (i = 0; i < nb_of_lines; i++) {
		if ((int)trans_list[nb_arr[i]].up_lev.j % 2 == 1)
		{
			output << left << setw(6) << i
				<< setw(5) << trans_list[nb_arr[i]].up_lev.v << setw(5) << (int)trans_list[nb_arr[i]].up_lev.j
				<< setw(5) << trans_list[nb_arr[i]].low_lev.v << setw(5) << (int)trans_list[nb_arr[i]].low_lev.j
				<< setw(14) << 1.e+4 * SPEED_OF_LIGHT / trans_list[nb_arr[i]].freq << setw(14) << lum_arr[i]
				<< endl;
		}
	}
	output.close();
}
