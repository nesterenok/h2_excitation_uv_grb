// h2_excitation_uv_grb.cpp 
// Important notes:
// The Gauss unit system is adopted elsewhere, exceptions: arguments of the cross sections, 
// concentration [cm-3]


#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#define _CRT_SECURE_NO_WARNINGS

// SUNDIALS CVODE solver headers
// works only with static libraries (shared must not be installed)
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* definition of type realtype          */

#include <omp.h>
#include <cstring>
#include <ctime>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "constants.h"
#include "h2_parameters_uv_grb.h"
#include "dynamic_array.h"
#include "h2_population_evol_eqs.h"
#include "h2_saving_data_uv_grb.h"
#include "integration.h"


#define MAX_TEXT_LINE_WIDTH 3500  // check it, maximal size of the comment lines in the files in symbols,
#define SOURCE_NAME "h2_excitation_uv_grb.cpp"

using namespace std;
const int VERBOSITY = 1;

// sim_path    - path to the folder with simulation data (GRB propagation results)
// output_path - path to the folder where data of the current simulations must be saved,
void init_cloud_parameters(const string& sim_path, double& grb_cloud_distance, double& max_hcolumn_density, double& conc_h_tot, 
	double& op_ratio_h2, double& h2_fraction, int & jet_type, double& half_opening_angle, double & theta_wing, double& viewing_angle);

void init_layer_coordinates(const string& path, vector<double>& layer_grid_coord, vector<double>& hcolumn_density_coord);

void init_specimen_conc(const string& sim_path, vector<dynamic_array>& d);
void init_h2_popdensity(const string& sim_path, vector<dynamic_array>& d, vector<int>& qnb_v_arr, vector<int>& qnb_j_arr);

void init_time_intervals(vector<double>& time_arr);
void simulations_grb_uv(const string& data_path, const string& sim_path, const string& output_path);

// Here, H column density  of the cloud layer centre (from the front cloud boundary) must be = integer value * 1.e+18
void init_h2_popdensity_evol(const string& output_path, const vector<double>& layer_grid_coord, const vector<double>& hcolumn_density_coord,
	vector<double>& time_arr, vector< vector<dynamic_array>>& h2_popdens_evol_arr, double conc_h_tot, double vangle, int nb_lev_h2, int nb_layers);

void calc_line_luminosities(const string& data_path, const string& sim_path, const string& output_path,
	double hcolumn_density_cloud, double theta_jet_axis, bool is_jet_top_hat);


static int f_h2ev(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
	static_cast<h2_population_evolution_data*>(user_data)->f(t, y, ydot);
	return (0);
}



int main()
{
	bool is_jet_top_hat;
	int nb_processors(4), verbosity(1);
	double hcolumn_density_cloud, theta_jet_axis;
	
	string fname;
	string data_path;    // path to the input_data folder with spectroscopic, collisional data and etc.
	string output_path;  // path to the folder where data of the current simulation must be saved,
	string sim_path;     // path to the folder with simulation data (GRB propagation results)

	stringstream ss;
	ifstream input;

#ifdef _OPENMP
	omp_set_num_threads(nb_processors);

#pragma omp parallel 
	{
#pragma omp master 
		{
			cout << "OpenMP is supported" << endl;
			cout << "Nb of threads: " << omp_get_num_threads() << endl;
		}
	}
#else
	cout << "No OpenMP" << endl;
#endif

	data_path = "C:/Users/Александр/Documents/input_data/";
	sim_path = "C:/Users/Александр/Documents/Данные и графики/paper GRB H2 infrared emission ultra-violet/sim_grb_cloud_2026/output_2e2_d100_va020/";
	output_path = "../../../output/nh2e2_va020/";
	//simulations_grb_uv(data_path, sim_path, output_path);


	// the file with viewing angles, for which simulations were performed, is located in the source directory,
	// sim_path is one of the directories with the simulation data of grb-cloud interaction
	// if is_jet_top_hat = true, this data is used in top_hat jet simulations
	sim_path = "C:/Users/Александр/Documents/Данные и графики/paper GRB H2 infrared emission ultra-violet/sim_grb_cloud_2026/output_2e2_d100_va0/";
	output_path = "../../../output/";
	
	hcolumn_density_cloud = 2.e+21;   // [cm-2]
	is_jet_top_hat = true;

	// angle between observation line of sight and jet axis,
	theta_jet_axis = 0.05;   // [radians]
	
	calc_line_luminosities(data_path, sim_path, output_path, hcolumn_density_cloud, theta_jet_axis, is_jet_top_hat);
}


// Note, the time grid is important in integration n(t) dt
void init_time_intervals(vector<double>& time_arr)
{
	int i;
	double x1, x2, x3, x4, x5, x6;
	
	time_arr.clear();
	time_arr.push_back(0.);
	time_arr.push_back(MIN_MODEL_TIME);

	x1 = pow(10., 1. / NB_OF_BINS_PER_ORDER_TIME_1);
	x2 = pow(10., 1. / NB_OF_BINS_PER_ORDER_TIME_2);
	x3 = pow(10., 1. / NB_OF_BINS_PER_ORDER_TIME_3);
	x4 = pow(10., 1. / NB_OF_BINS_PER_ORDER_TIME_4);
	x5 = pow(10., 1. / NB_OF_BINS_PER_ORDER_TIME_5);
	x6 = pow(10., 1. / NB_OF_BINS_PER_ORDER_TIME_6);

	for (i = 0; time_arr.back() < MODEL_TIME_1 * (1. - 1.e-9); i++) {
		time_arr.push_back(time_arr.back() * x1);
	}
	for (i = 0; time_arr.back() < MODEL_TIME_2 * (1. - 1.e-9); i++) {
		time_arr.push_back(time_arr.back() * x2);
	}
	for (i = 0; time_arr.back() < MODEL_TIME_3 * (1. - 1.e-9); i++) {
		time_arr.push_back(time_arr.back() * x3);
	}
	for (i = 0; time_arr.back() < MODEL_TIME_4 * (1. - 1.e-9); i++) {
		time_arr.push_back(time_arr.back() * x4);
	}
	for (i = 0; time_arr.back() < MODEL_TIME_5 * (1. - 1.e-9); i++) {
		time_arr.push_back(time_arr.back() * x5);
	}
	for (i = 0; time_arr.back() < MAX_MODEL_TIME; i++) {
		time_arr.push_back(time_arr.back() * x6);
	}
}


void simulations_grb_uv(const string& data_path, const string& sim_path, const string& output_path)
{
	bool must_be_saved;
	int i, k, lay_nb, nb_of_layers, nb_lev_h2, nb_of_equat, flag, time_nb, step_nb, jet_type;

	double grb_cloud_distance, max_hcolumn_density, hcolumn_density, conc_h_tot, op_ratio_h2, h2_fraction,
		half_opening_angle, theta_wing, viewing_angle;
	
	double rel_tol, model_time_aux, model_time_aux_prev, model_time, model_time_step, model_time_out;

	time_t timer, timer_tot;
	char* ctime_str;
	stringstream ss;

	vector<int>    qnb_v_arr, qnb_j_arr;
	vector<double> layer_grid_coord, layer_centre_distances, hcolumn_density_coord, time_arr;  // in cm

	timer_tot = time(NULL);
	ctime_str = ctime(&timer_tot);
	if (VERBOSITY) {
		cout << ctime_str << endl << "Start of the simulations of the H2 population density evolution" << endl;
	}

	init_time_intervals(time_arr);

	// user data
	h2_population_evolution_data* user_data
		= new h2_population_evolution_data(data_path, output_path, VERBOSITY);

	nb_lev_h2 = user_data->get_nb_of_h2_lev();
	nb_of_equat = nb_lev_h2;

	vector<dynamic_array> h2_popdens_evol, h2_popdens_v_evol, specimen_conc_init_data, h2_popdens_init_data;
	dynamic_array h2_popdens(nb_lev_h2), h2_popdens_v(NB_VIBR_H2);  // arrays are initialized by zeros in the class constructor,

	init_cloud_parameters(sim_path, grb_cloud_distance, max_hcolumn_density, conc_h_tot, op_ratio_h2, h2_fraction,
		jet_type, half_opening_angle, theta_wing, viewing_angle);
	
	init_layer_coordinates(sim_path, layer_grid_coord, hcolumn_density_coord);

	nb_of_layers = (int)layer_grid_coord.size() - 1;
	for (i = 0; i < nb_of_layers; i++) {
		layer_centre_distances.push_back(0.5 * (layer_grid_coord[i] + layer_grid_coord[i + 1]));
	}

	// specimen concentrations are initialized here,
	init_specimen_conc(sim_path, specimen_conc_init_data);

	// population density in the file has dimension [cm-3],
	init_h2_popdensity(sim_path, h2_popdens_init_data, qnb_v_arr, qnb_j_arr);

	// vectors used by solver SUNDIALS CVODE
	N_Vector y, ydot, abs_tol;

	cout << scientific;
	cout.precision(3);

	SUNMatrix A;
	SUNLinearSolver LS;
	void* cvode_mem;

	y = N_VNew_Serial(nb_of_equat);
	ydot = N_VNew_Serial(nb_of_equat);
	abs_tol = N_VNew_Serial(nb_of_equat);

	// Initialization for tolerances:
	rel_tol = REL_ERROR_SOLVER;
	user_data->set_tolerances(abs_tol);

	// Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula and the use of a Newton iteration 
	cvode_mem = CVodeCreate(CV_BDF);

	// Call CVodeInit to initialize the integrator memory and specify the user's right hand side function in y'=f(t,y), 
	// the initial time t0, and the initial dependent variable vector y;
	// y is not defined yet,
	flag = CVodeInit(cvode_mem, f_h2ev, model_time = 0., y);

	// Call CVodeSVtolerances to specify the scalar tolerances:
	flag = CVodeSVtolerances(cvode_mem, rel_tol, abs_tol);

	// The maximal number of steps between simulation stops;
	flag = CVodeSetMaxNumSteps(cvode_mem, 10000);

	// specifies the maximum number of error test failures permitted in attempting one step:
	flag = CVodeSetMaxErrTestFails(cvode_mem, MAX_ERR_TEST_FAILS_SOLVER); // default value is 7; 
	flag = CVodeSetMaxConvFails(cvode_mem, MAX_CONV_FAILS_SOLVER); // default value is 10;

	// Create dense SUNMatrix for use in linear solves 
	A = SUNDenseMatrix(nb_of_equat, nb_of_equat);

	// Create dense SUNLinearSolver object for use by CVode 
	LS = SUNLinSol_Dense(y, A);

	// Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode 
	flag = CVodeSetLinearSolver(cvode_mem, LS, A);

	// The function attaches the user data block to the solver;
	flag = CVodeSetUserData(cvode_mem, user_data);

	// the cycle over layer numbers,
	for (lay_nb = 0; lay_nb < nb_of_layers; lay_nb++)
	{
		if (VERBOSITY) {
			cout << left << "lay_nb: " << setw(7) << lay_nb;
		}
		timer = time(NULL);

		h2_popdens_evol.clear();
		h2_popdens_v_evol.clear();

		// population densities must be in [cm-3],
		for (i = 0; i < h2_popdens_init_data[lay_nb].dim; i++)
		{
			k = user_data->get_level_nb_h2(qnb_v_arr[i], qnb_j_arr[i]);
			if (k >= 0) {
				h2_popdens.arr[k] = h2_popdens_init_data[lay_nb].arr[i];
			}
		}

		// Initialization of initial values of variables,
		for (i = 0; i < nb_of_equat; i++) {
			NV_Ith_S(y, i) = h2_popdens.arr[i];
		}

		h2_popdens_evol.push_back(h2_popdens);

		memset(h2_popdens_v.arr, 0, NB_VIBR_H2 * sizeof(double));
		for (i = 0; i < nb_lev_h2; i++) {
			k = user_data->get_vibr_nb_h2(i);
			if (k >= 0) {
				h2_popdens_v.arr[k] += NV_Ith_S(y, i);
			}
		}
		h2_popdens_v_evol.push_back(h2_popdens_v);

		//
		// Simulations
		//
		model_time_aux = model_time = 0.;

		// restart of the solver with new values of initial conditions,
		flag = CVodeReInit(cvode_mem, model_time, y);

		// update the members of user_data class,
		f_h2ev(model_time, y, ydot, (void*)user_data);

		time_nb = 1;
		// there is integration of the parameters over the time, this time step must be sufficiently small, 
		model_time_step = 1.;

		while (model_time < time_arr.back())
		{
			flag = CV_SUCCESS;
			must_be_saved = false;
			model_time_aux_prev = model_time_aux;

			if (model_time < time_arr[time_nb])
			{
				step_nb = 0;
				model_time_out = model_time + model_time_step;  // ?

				while (step_nb < MAX_NB_STEPS && flag == CV_SUCCESS && model_time < model_time_out)
				{
					// actual model time is stored in cvode_mem, 
					flag = CVode(cvode_mem, model_time_out, y, &model_time, CV_ONE_STEP);  // CV_NORMAL or CV_ONE_STEP	
					step_nb++;
				}
			}

			if (flag == CV_SUCCESS) {
				if (model_time > time_arr[time_nb]) {
					// step back to the time grid value, updating model time,
					flag = CVodeGetDky(cvode_mem, time_arr[time_nb], 0, y);
					model_time_aux = time_arr[time_nb];

					must_be_saved = true;
					time_nb++;

					if (flag != CV_SUCCESS) {
						cout << "Error in CVodeGetDky() " << flag << endl;
					}
				}
				else {
					model_time_aux = model_time;
				}

				// saving data,
				if (must_be_saved) {
					for (i = 0; i < nb_lev_h2; i++) {
						h2_popdens.arr[i] = NV_Ith_S(y, i);
					}
					h2_popdens_evol.push_back(h2_popdens);

					memset(h2_popdens_v.arr, 0, NB_VIBR_H2 * sizeof(double));
					for (i = 0; i < nb_lev_h2; i++) {
						k = user_data->get_vibr_nb_h2(i);
						if (k >= 0) {
							h2_popdens_v.arr[k] += NV_Ith_S(y, i);
						}
					}
					h2_popdens_v_evol.push_back(h2_popdens_v);
				}
			}
			else {
				cout << "Some unexpected error has been occurred: " << flag << endl;
				break;
			}
		}

		if (VERBOSITY) {
			i = (int)(time(NULL) - timer);
			cout << left << "Calc time (s): " << i << endl;
		}

		// H column density from the cloud boundary to the cloud layer centre,
		// Here, H column density must be = integer value * 1.e+18
		hcolumn_density = conc_h_tot * (layer_centre_distances[lay_nb] - grb_cloud_distance);

		save_h2_populations_evolution(output_path, time_arr, h2_popdens_evol, conc_h_tot, hcolumn_density, nb_lev_h2);
	}

	if (VERBOSITY) {
		cout << "The memory is freeing up" << endl;
		cout << left << "Total calc time (min): " << (int)(time(NULL) - timer_tot) / 60 << endl;
	}
	N_VDestroy(y);
	N_VDestroy(ydot);
	N_VDestroy(abs_tol);

	// Free integrator memory
	CVodeFree(&cvode_mem);

	// Free the linear solver memory
	SUNLinSolFree(LS);

	// Free the matrix memory
	SUNMatDestroy(A);
}


// initialization of parameters of the gas-dust cloud
void init_cloud_parameters(const string& sim_path, double& grb_cloud_distance, double & max_hcolumn_density, double& conc_h_tot, 
	double& op_ratio_h2, double & h2_fraction, int& jet_type, double& half_opening_angle, double& theta_wing, double& viewing_angle)
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


void init_layer_coordinates(const string& path, vector<double> & layer_grid_coord, vector<double> & hcolumn_density_coord)
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


void calc_line_luminosities(const string& data_path, const string& sim_path, const string& output_path, 
	double hcolumn_density_cloud, double theta_jet_axis, bool is_jet_top_hat)
{
	bool quit_cycle;
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, k, v, nb, lay, lev, nb_layers, nb_lev_h2, nb_times, nb_observ_times, nb_viewing_angles, isotope, jet_type;
	
	// parameters of the jet and molecular cloud,
	double x, h2_fraction, conc_h_tot, op_ratio_h2, cloud_width, grb_cloud_distance, max_hcolumn_density,
		half_opening_angle, theta_wing, cos_theta_jet_axis, sin_theta_jet_axis;

	double obst_const, observ_time, max_observ_time, max_delta_time, luminosity, vangle, max_cos_theta, cos_theta_observ_time,
		 min_cos_theta_prime, max_cos_theta_prime, cont_abs_const;
	
	vector<double> layer_grid_coord, layer_centre_distances, hcolumn_density_coord, observ_time_arr, time_arr,
		viewing_angle_grid, cos_viewing_angle_grid;
	
	string fname, str_id;
	stringstream ss;
	ifstream input;

	table_integrals table_int;

	time_t timer;
	char* ctime_str;

	timer = time(NULL);
	ctime_str = ctime(&timer);
	if (VERBOSITY) {
		cout << ctime_str << endl << "Start of the simulations of the H2 line luminosities" << endl;
	}

	// new H2 molecule data (Roueff et al. A&A 630, A58, 2019), max number of levels is 302,
	nb_lev_h2 = 302;
	molecule h2_mol("H2", isotope = 1, 2. * ATOMIC_MASS_UNIT);

	h2_diagram_roueff2019 *h2_di 
		= new h2_diagram_roueff2019(data_path, h2_mol, nb_lev_h2, VERBOSITY);
	
	h2_einstein_coeff_roueff2019 *h2_einst 
		= new h2_einstein_coeff_roueff2019(data_path, h2_di, VERBOSITY);

	double* h2_popdens_integr_theta = new double[nb_lev_h2];
	dynamic_array h2_popdens_integr_volume(nb_lev_h2);
	
	// vector index - time moment number, array index - H2 level number,
	vector<dynamic_array> h2_popdens_integr_volume_arr;

	// first vector index - layer number, second vector index - time moment number, third array index - H2 level number,
	vector< vector<dynamic_array> > h2_popdens_evol_arr_1, h2_popdens_evol_arr_2;

	// reading data common to all simulation data, 
	init_cloud_parameters(sim_path, grb_cloud_distance, max_hcolumn_density, conc_h_tot, op_ratio_h2, h2_fraction,
		jet_type, half_opening_angle, theta_wing, vangle);

	init_layer_coordinates(sim_path, layer_grid_coord, hcolumn_density_coord);

	cloud_width = hcolumn_density_cloud / conc_h_tot;
	cont_abs_const = 0.45e-22 * conc_h_tot;  // K band, absorption in continuum, does not depend on molecular fraction

	for (lay = 0; lay < (int)layer_grid_coord.size() - 1; lay++) {
		layer_centre_distances.push_back(0.5 * (layer_grid_coord[lay] + layer_grid_coord[lay + 1]));
	}

	// number of layers that are taken into in luminosity calculations,
	for (lay = 0; (lay < (int)layer_grid_coord.size() - 1) && (layer_grid_coord[lay] - grb_cloud_distance < cloud_width); lay++)
	{;}
	nb_layers = lay;

	cos_theta_jet_axis = cos(theta_jet_axis);
	sin_theta_jet_axis = sin(theta_jet_axis);

	// calculation of integrals
	table_int.setp(cos_theta_jet_axis, sin_theta_jet_axis);

	if (is_jet_top_hat) {
		// only one case with on-axis viewing angle is calculated, given the simulation data path
		// the value of jet_type is ignored,
		nb_viewing_angles = 2;
	
		viewing_angle_grid.push_back(0.);
		viewing_angle_grid.push_back(half_opening_angle);
	}
	else {
		// reading file with viewing angles, for which simulations were performed, is located in the source directory,
		fname = "viewing_angles.txt";
		input.open(fname.c_str(), ios_base::in);

		if (!input.is_open()) {
			cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
			exit(1);
		}

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb_viewing_angles;

		for (i = 0; i < nb_viewing_angles; i++) {
			input >> x;
			viewing_angle_grid.push_back(x);
		}
		input.close();
	}

	for (v = 0; v < nb_viewing_angles; v++) {
		cos_viewing_angle_grid.push_back(cos(viewing_angle_grid[v]));
	}
	
	// the characteristic time of H2 spontaneous de-excitation is of the order of 1.e+7 s,
	max_observ_time = 1.e+7;
	max_observ_time += layer_grid_coord[nb_layers] * (1. - cos(viewing_angle_grid.back() + theta_jet_axis)) / SPEED_OF_LIGHT;  // in [s]
	
	nb_observ_times = (int)(max_observ_time / DELTA_OBSERV_TIME) + 1;
	for (i = 0; i < nb_observ_times; i++) {
		observ_time_arr.push_back(DELTA_OBSERV_TIME * i);
	}

	// the main container with integration data: 
	for (i = 0; i < nb_observ_times; i++) {
		h2_popdens_integr_volume_arr.push_back(h2_popdens_integr_volume);
	}

	// Initialization of data on population densities for viewing angle equal to zero,
	vangle = viewing_angle_grid[0];
	init_h2_popdensity_evol(output_path, layer_grid_coord, hcolumn_density_coord, time_arr, h2_popdens_evol_arr_1, 
		conc_h_tot, vangle, nb_lev_h2, nb_layers);

	nb_times = (int)time_arr.size();

 	for (v = 1; v < nb_viewing_angles; v++) {	
		if (is_jet_top_hat) {
			for (lay = 0; lay < (int)h2_popdens_evol_arr_1.size(); lay++) {
				h2_popdens_evol_arr_2.push_back(h2_popdens_evol_arr_1[lay]);
			}
		}
		else {
			// reading the data on population densities for given viewing angle,
			vangle = viewing_angle_grid[v];

			init_h2_popdensity_evol(output_path, layer_grid_coord, hcolumn_density_coord, time_arr, h2_popdens_evol_arr_2,
				conc_h_tot, vangle, nb_lev_h2, nb_layers);
		}

		max_cos_theta_prime = cos(viewing_angle_grid[v - 1]);   // closer to the axis
		min_cos_theta_prime = cos(viewing_angle_grid[v]);       // farer from the axis

		// for each observational time, the calculation of volume-integrated population densities,
		if (VERBOSITY) {
			cout << "Calculation of volume-integrated population densities for time moments: " << endl;
		}
		for (nb = 0; nb < nb_observ_times; nb++) {	
			observ_time = observ_time_arr[nb];
			
			if (VERBOSITY) {
				cout << left << setw(6) << nb;
				if (nb == nb_observ_times - 1) {
					cout << endl;
				}
			}

			for (lev = 0; lev < nb_lev_h2; lev++) {
				h2_popdens_integr_volume.arr[lev] = 0.;
			}

			// integration by depth (or, equivalently, summation of cloud layers),
			for (lay = 0; lay < nb_layers; lay++)
			{
				memset(h2_popdens_integr_theta, 0, nb_lev_h2 * sizeof(double));
				obst_const = SPEED_OF_LIGHT / layer_centre_distances[lay];

				// upper limit in integration by cos theta, the inner border of the ring
				max_cos_theta = 1.;

				cos_theta_observ_time = 1. - obst_const * observ_time;
				max_delta_time = (max_cos_theta - cos_theta_observ_time) / obst_const;

				if (max_delta_time > time_arr.back()) {
					max_delta_time = time_arr.back();
					max_cos_theta = cos_theta_observ_time + obst_const * max_delta_time;
				}
				
				// integration by theta, or equivalently, time,
#pragma omp parallel private(k, lev, x)
				{
					double x1, x2, tau;
					double cos_theta, cos_theta_1, cos_theta_2, sin_theta, sin_theta_sq, cos_theta_prime_1, cos_theta_prime_2, 
						integral_by_phi, linear_integral_by_phi, max_cos_phi, min_cos_phi, min_phi, max_phi;

					double * arr_h2 = new double[nb_lev_h2];
					memset(arr_h2, 0, nb_lev_h2 * sizeof(double));
					 
#pragma omp for schedule(dynamic, 1)
					for (k = 0; k < nb_times - 1; k++) 
					{
						const dynamic_array & h2_popdens_1 = h2_popdens_evol_arr_1[lay][k];
						const dynamic_array & h2_popdens_2 = h2_popdens_evol_arr_2[lay][k];

						// conversion from time to cos theta
						cos_theta_1 = cos_theta_observ_time + obst_const * time_arr[k];
						cos_theta_2 = cos_theta_observ_time + obst_const * time_arr[k + 1];

						if (cos_theta_1 < 1. && cos_theta_2 > 1.) {
							cos_theta_2 = 1.;
						}
						cos_theta = 0.5 * (cos_theta_1 + cos_theta_2);

						if (cos_theta < max_cos_theta) {
							sin_theta_sq = 1. - cos_theta * cos_theta;
							sin_theta = sqrt(sin_theta_sq);

							// in the frame of reference of the jet axis,
							cos_theta_prime_1 = cos_theta_jet_axis * cos_theta + sin_theta_jet_axis * sin_theta;  // closest point to the axis
							cos_theta_prime_2 = cos_theta_jet_axis * cos_theta - sin_theta_jet_axis * sin_theta;  // more distant point from the axis
							integral_by_phi = linear_integral_by_phi = 0.;

							// the ring of the jet structure is inside the circle of phi integration,
							if (cos_theta_prime_1 > max_cos_theta_prime && cos_theta_prime_2 < min_cos_theta_prime)
							{
								max_cos_phi = (max_cos_theta_prime - cos_theta_jet_axis * cos_theta) / (sin_theta_jet_axis * sin_theta);
								min_phi = acos(max_cos_phi);  // value in the range [0, pi] is returned,

								min_cos_phi = (min_cos_theta_prime - cos_theta_jet_axis * cos_theta) / (sin_theta_jet_axis * sin_theta);
								max_phi = acos(min_cos_phi);  // value in the range [0, pi] is returned,

								integral_by_phi = 2. * (max_phi - min_phi);
								linear_integral_by_phi = 2. * (table_int.get(max_phi, cos_theta) - table_int.get(min_phi, cos_theta));
							}
							// the ring of the jet structure contains the more distant point (from jet axis) of circle of phi integration
							else if (cos_theta_prime_1 > max_cos_theta_prime
								&& cos_theta_prime_2 > min_cos_theta_prime && cos_theta_prime_2 < max_cos_theta_prime)
							{
								max_cos_phi = (max_cos_theta_prime - cos_theta_jet_axis * cos_theta) / (sin_theta_jet_axis * sin_theta);
								min_phi = acos(max_cos_phi);  // value in the range [0, pi] is returned,

								integral_by_phi = 2. * (M_PI - min_phi);
								linear_integral_by_phi = 2. * (table_int.get(M_PI, cos_theta) - table_int.get(min_phi, cos_theta));
							}
							// the ring of the jet structure contains the more closest point of circle of phi integration
							else if (cos_theta_prime_1 < max_cos_theta_prime && cos_theta_prime_1 > min_cos_theta_prime
								&& cos_theta_prime_2 < min_cos_theta_prime)
							{
								min_cos_phi = (min_cos_theta_prime - cos_theta_jet_axis * cos_theta) / (sin_theta_jet_axis * sin_theta);
								max_phi = acos(min_cos_phi);  // value in the range [0, pi] is returned,

								integral_by_phi = 2. * max_phi;
								linear_integral_by_phi = 2. * table_int.get(max_phi, cos_theta);
							}
							// the ring of the jet structure contains the circle of phi integration
							else if (cos_theta_prime_1 < max_cos_theta_prime && cos_theta_prime_2 > min_cos_theta_prime) {
								integral_by_phi = 2. * M_PI;
								linear_integral_by_phi = 2. * table_int.get(M_PI, cos_theta);
							}

							if (integral_by_phi > 0.)
							{
								// calculation of optical depth to the right boundary of the cloud, does not depend on phi
								// it is assumed that absorption is weak, and H2 ionization/dissociation boundary is sharp,
								tau = 0.;
								if (layer_centre_distances[lay] - grb_cloud_distance < cloud_width)
								{
									tau = sqrt((cloud_width + grb_cloud_distance) * (cloud_width + grb_cloud_distance) - layer_centre_distances[lay] * layer_centre_distances[lay] * sin_theta_sq)
										- layer_centre_distances[lay] * cos_theta;

									tau *= cont_abs_const;
								}

								x = (cos_theta_2 - cos_theta_1) * exp(-tau) / (viewing_angle_grid[v] - viewing_angle_grid[v - 1]);
								
								for (lev = 0; lev < nb_lev_h2; lev++) {
									x1 = h2_popdens_1.arr[lev];
									x2 = h2_popdens_2.arr[lev];

									// [cm-3], not divided by statistical weight,
									arr_h2[lev] += (integral_by_phi * (x1 * viewing_angle_grid[v] - x2 * viewing_angle_grid[v - 1])
										+ linear_integral_by_phi * (x2 - x1)) * x;
								}
							}
						}
					}
					
#pragma omp critical
					{
						for (lev = 0; lev < nb_lev_h2; lev++) {
							h2_popdens_integr_theta[lev] += arr_h2[lev];
						}
					}
					delete[] arr_h2;
				}
				
				for (lev = 0; lev < nb_lev_h2; lev++) {
					if (layer_grid_coord[lay + 1] - grb_cloud_distance < cloud_width)
					{
						h2_popdens_integr_volume.arr[lev] += h2_popdens_integr_theta[lev]
							* layer_centre_distances[lay] * layer_centre_distances[lay] * (layer_grid_coord[lay + 1] - layer_grid_coord[lay]);
					}
					else {
						h2_popdens_integr_volume.arr[lev] += h2_popdens_integr_theta[lev]
							* layer_centre_distances[lay] * layer_centre_distances[lay] * (cloud_width + grb_cloud_distance - layer_grid_coord[lay]);
					}
				}
			}
			for (lev = 0; lev < nb_lev_h2; lev++) {
				h2_popdens_integr_volume_arr[nb].arr[lev] += h2_popdens_integr_volume.arr[lev];
			}
		}
	
		if (!is_jet_top_hat) {
			h2_popdens_evol_arr_1.clear();
			
			for (lay = 0; lay < (int)h2_popdens_evol_arr_2.size(); lay++) {
				h2_popdens_evol_arr_1.push_back(h2_popdens_evol_arr_2[lay]);
			}
		}
	}

	vector<transition> trans_list;
	dynamic_array line_luminosity_evol(nb_observ_times);

	// (first) vector index: transition number, (second) array index - time moment number,
	vector<dynamic_array> line_luminosity_evol_arr;

	// Search for high-luminosity lines,
	// loop over H2 ro-vibrational transitions:
	for (i = 1; i < nb_lev_h2; i++) {
		for (j = 0; j < i; j++) 
		{
			if (h2_einst->arr[i][j] > 1.e-99) {
				quit_cycle = false;

				// for each observational time, the calculation of luminosity,
				for (nb = 0; (nb < nb_observ_times) && (!quit_cycle); nb++)
				{
					luminosity = h2_einst->arr[i][j] * (h2_di->lev_array[i].energy - h2_di->lev_array[j].energy) * CM_INVERSE_TO_ERG
						* h2_popdens_integr_volume_arr[nb].arr[i] / SOLAR_LUMINOSITY;

					if (luminosity > MIN_LUMINOSITY) {
						trans_list.push_back(transition(h2_di->lev_array[j], h2_di->lev_array[i]));
						quit_cycle = true;
					}
				}
			}
		}
	}
	sort(trans_list.begin(), trans_list.end());

	for (k = 0; k < (int)trans_list.size(); k++) 
	{
		i = trans_list[k].up_lev.nb;
		j = trans_list[k].low_lev.nb;
		
		for (nb = 0; nb < nb_observ_times; nb++)
		{
			luminosity = h2_einst->arr[i][j] * (h2_di->lev_array[i].energy - h2_di->lev_array[j].energy) * CM_INVERSE_TO_ERG
				* h2_popdens_integr_volume_arr[nb].arr[i] / SOLAR_LUMINOSITY;

			line_luminosity_evol.arr[nb] = luminosity;
		}
		line_luminosity_evol_arr.push_back(line_luminosity_evol);
	}

	save_h2_line_luminosities(output_path, observ_time_arr, trans_list, line_luminosity_evol_arr, theta_jet_axis, is_jet_top_hat);

	if (VERBOSITY) {
		i = (int)(time(NULL) - timer);
		cout << left << "Calc time (min): " << setw(7) << i / 60 << endl;
	}
}


// Here, H column density must be = integer value * 1.e+18
void init_h2_popdensity_evol(const string & output_path, const vector<double> & layer_grid_coord, const vector<double> & hcolumn_density_coord, 
	vector<double> & time_arr, vector< vector<dynamic_array>> & h2_popdens_evol_arr, double conc_h_tot, double vangle, int nb_lev_h2, int nb_layers)
{
	char ch, text_line[MAX_TEXT_LINE_WIDTH];
	int i, k, lev, lay, nb_times, nb_lev_f;
	double x;

	string fname, sub_dir;
	ifstream input;
	stringstream ss;

	dynamic_array h2_popdens(nb_lev_h2);
	vector<dynamic_array> h2_popdens_evol;

	ss.clear();
	ss.str("");

	i = (int)log10(conc_h_tot);  // ? treat correctly values like 1.e+2
	k = rounding(conc_h_tot / pow(10., i));

	ss << "nh" << k << "e" << i;
	ss << "_va";

	k = rounding(vangle * 100);
	if (k > 0 && k < 10) {
		ss << "00";
	}
	else if (k >= 10 && k < 100) {
		ss << "0";
	}
	ss << k;
	sub_dir = ss.str();

	if (VERBOSITY) {
		cout << "Initialization of population density data for viewing angle: " << vangle << endl;
	}

	h2_popdens_evol_arr.clear();
	for (lay = 0; lay < nb_layers; lay++)
	{
		// H column density from the front boundary to the ccloud layer centre,
		x = 0.5 * (hcolumn_density_coord[lay] + hcolumn_density_coord[lay + 1]);
		
		i = 18;
		k = rounding(x / pow(10., i));

		ss.clear();
		ss.str("");
		ss << "_nh" << k << "e" << i;

		fname = output_path + sub_dir;
		fname += "/h2grb_popdens";
		fname += ss.str();
		fname += ".txt";
		input.open(fname.c_str());

		if (!input.is_open()) {
			cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
			exit(1);
		}

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		ss.clear();
		ss.str(text_line);
		ss >> ch >> x >> nb_times >> nb_lev_f;

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		time_arr.clear();
		h2_popdens_evol.clear();

		for (k = 0; k < nb_times; k++) {
			input >> x;
			time_arr.push_back(x);

			for (lev = 0; lev < nb_lev_f; lev++) {
				input >> x;
				h2_popdens.arr[lev] = x;
			}
			h2_popdens_evol.push_back(h2_popdens);
		}
		input.close();

		// averaging of population densities over two time moments
		for (k = 0; k < nb_times - 1; k++) {
			for (lev = 0; lev < nb_lev_h2; lev++) {
				h2_popdens_evol[k].arr[lev] = 0.5 * (h2_popdens_evol[k].arr[lev] + h2_popdens_evol[k + 1].arr[lev]);
			}
		}

		h2_popdens_evol_arr.push_back(h2_popdens_evol);

		if (VERBOSITY) {
			cout << left << setw(6) << lay;
			if (lay == nb_layers - 1) {
				cout << endl;
			}
		}
	}
}






/*
for (lev = 0; lev < nb_lev_h2; lev++) {
						// f(x) = (f1 * x2 - f2 * x1 + x * (f2 - f1)) / (x2 - x1)
						// int f(x) dx = ((f1 * x2 - f2 * x1) * (x'' - x') + 0.5 * (f2 - f1) * (x''2 - x'2)) / (x2 - x1)

						if ((cos_theta_1 < min_cos_theta_jet) && (cos_theta_2 > min_cos_theta_jet) && (cos_theta_2 < max_cos_theta))
						{
							x = ((h2_popdens_evol[k].arr[lev] * cos_theta_2 - h2_popdens_evol[k + 1].arr[lev] * cos_theta_1) * (cos_theta_2 - min_cos_theta_jet)
								+ 0.5 * (h2_popdens_evol[k + 1].arr[lev] - h2_popdens_evol[k].arr[lev]) * (cos_theta_2 * cos_theta_2 - min_cos_theta_jet * min_cos_theta_jet))
								/ (cos_theta_2 - cos_theta_1);
						}
						else if ((cos_theta_1 > min_cos_theta_jet) && (cos_theta_2 < max_cos_theta))
						{
							x = 0.5 * (h2_popdens_evol[k + 1].arr[lev] + h2_popdens_evol[k].arr[lev]) * (cos_theta_2 - cos_theta_1);
						}
						else if ((cos_theta_1 > min_cos_theta_jet) && (cos_theta_1 < max_cos_theta) && (cos_theta_2 > max_cos_theta))
						{
							x = ((h2_popdens_evol[k].arr[lev] * cos_theta_2 - h2_popdens_evol[k + 1].arr[lev] * cos_theta_1) * (max_cos_theta - cos_theta_1)
								+ 0.5 * (h2_popdens_evol[k + 1].arr[lev] - h2_popdens_evol[k].arr[lev]) * (max_cos_theta * max_cos_theta - cos_theta_1 * cos_theta_1))
								/ (cos_theta_2 - cos_theta_1);
						}
						else if ((cos_theta_1 < min_cos_theta_jet) && (cos_theta_2 > max_cos_theta)) {
							x = 0.;  // it is not important?
						}

						// [cm-3], not divided by statistical weight,
						h2_popdens_integr_theta.arr[lev] += x * exp(-tau);
					}

*/