#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include<cstring>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#include"init_simulation_data.h"
#include"dynamic_array.h"
#include"parameters.h"

#define MAX_TEXT_LINE_WIDTH 3500  // check it, maximal size of the comment lines in the files in symbols,
#define SOURCE_NAME "init_simulation_data.cpp"

using namespace std;

// initialization of parameters of the gas-dust cloud
void init_cloud_parameters(const string& sim_path, double& grb_cloud_distance, double& max_hcolumn_density, double& conc_h_tot,
	double& op_ratio_h2, double& h2_fraction, int& jet_type, double& half_opening_angle, double& theta_wing, double& viewing_angle)
{
	int i;
	char text_line[MAX_TEXT_LINE_WIDTH];

	string fname;
	stringstream ss;
	ifstream input;

	fname = sim_path + "grb_cloud_parameters.txt";
	input.open(fname.c_str());

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with data " << fname << endl;
		exit(1);
	}

	while (!input.eof())
	{
		for (i = 0; i < 5; i++) {
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		}

		ss.clear();
		ss.str(text_line);
		ss >> conc_h_tot;  // [cm-3]

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		ss.clear();
		ss.str(text_line);
		ss >> op_ratio_h2;

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		ss.clear();
		ss.str(text_line);
		ss >> h2_fraction;

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		ss.clear();
		ss.str(text_line);
		ss >> grb_cloud_distance;  // [cm], distance from the GRB source to the cloud boundary,

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		ss.clear();
		ss.str(text_line);
		ss >> max_hcolumn_density;  // [cm-2], maximal column density of the cloud,

		for (i = 0; i < 9; i++) {
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		}

		ss.clear();
		ss.str(text_line);
		ss >> jet_type;

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		ss.clear();
		ss.str(text_line);
		ss >> half_opening_angle;  // [rad],

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		ss.clear();
		ss.str(text_line);
		ss >> theta_wing;          // [rad],

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		ss.clear();
		ss.str(text_line);
		ss >> viewing_angle;      // [rad],

		input.close();
	}
}


void init_layer_coordinates(const string& path, vector<double>& layer_grid_coord, vector<double>& hcolumn_density_coord)
{
	int i, l, nb;
	double x1, x2, x3, x4;
	char text_line[MAX_TEXT_LINE_WIDTH];

	string fname;
	ifstream input;

	fname = path + "grb_layer_coordinates.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb;

	layer_grid_coord.clear();
	for (i = 0; i < nb; i++) {
		input >> l >> x1 >> x2 >> x3 >> x4;

		layer_grid_coord.push_back(x1);
		hcolumn_density_coord.push_back(x3);
	}
	layer_grid_coord.push_back(x2);
	hcolumn_density_coord.push_back(x4);

	input.close();
}


void init_specimen_conc(const string& path, vector<dynamic_array>& d)
{
	int i, j, k, nb_l, nb_sp;
	double x, conc_h_tot;
	char ch, text_line[MAX_TEXT_LINE_WIDTH];

	string fname, sn;
	stringstream ss;
	ifstream input;

	vector<string> chemical_species_file;
	dynamic_array specimen_conc(NB_OF_CHEM_SPECIES);

	fname = path + "grb_specimen_conc.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	ss.clear();
	ss.str(text_line);
	ss >> ch >> conc_h_tot >> nb_l >> nb_sp;

	if (nb_sp != NB_OF_CHEM_SPECIES) {
		cout << "Note: nb of species in input file does not coincide with that in the code." << endl;
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	ss.clear();
	ss.str(text_line);

	ss >> sn >> sn;  // "!nb",  "distance(cm)"
	for (i = 0; i < nb_sp; i++) {
		ss >> sn;
		chemical_species_file.push_back(sn);
	}

	for (i = 0; i < nb_l; i++) {
		d.push_back(specimen_conc);

		input >> j >> x;
		for (j = 0; j < nb_sp; j++) {
			input >> x;

			for (k = 0; k < NB_OF_CHEM_SPECIES; k++) {
				if (chemical_species_file[j] == chemical_species[k]) {
					d[i].arr[k] = x;  // vector index - layer nb, array index - specimen nb,
				}
			}
		}
	}
	input.close();
}


// In the file, the levels are denoted by the number,
void init_h2_popdensity(const string& path, vector<dynamic_array>& d, vector<int>& qnb_v_arr, vector<int>& qnb_j_arr)
{
	int i, j, v, nb, nb_l, nb_lev_h2_f;
	double x, conc_h_tot, distance;
	char ch, text_line[MAX_TEXT_LINE_WIDTH];

	string fname, sn;
	stringstream ss;
	ifstream input;

	// Reading the data file with the list of H2 levels of the ground electronic state,
	// these levels were used in GRB propagation simulations,
	qnb_j_arr.clear();
	qnb_v_arr.clear();

	fname = path + "grb_h2_levels.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> nb_lev_h2_f;
	for (i = 0; i < nb_lev_h2_f; i++) {
		input >> nb >> v >> j >> x;

		qnb_v_arr.push_back(v);
		qnb_j_arr.push_back(j);
	}
	input.close();

	// reading the data file with the population densities of H2 molecule,
	fname = path + "grb_h2_popdens.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> ch >> conc_h_tot >> distance >> nb_l >> nb_lev_h2_f;

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);  // reading the end of the previous line
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	dynamic_array population_densities(nb_lev_h2_f);

	for (i = 0; i < nb_l; i++) {
		d.push_back(population_densities);
		input >> j >> x;

		for (j = 0; j < nb_lev_h2_f; j++) {
			input >> x;
			d[i].arr[j] = x;  // vector index - layer nb, array index - molecule level nb,
		}
	}
	input.close();
}


