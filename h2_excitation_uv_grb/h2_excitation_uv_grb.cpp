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
#include "parameters.h"
#include "dynamic_array.h"
#include "h2_population_evol_eqs.h"
#include "init_simulation_data.h"
#include "saving_data.h"
#include "integration.h"


#define MAX_TEXT_LINE_WIDTH 3500  // check it, maximal size of the comment lines in the files in symbols,
#define SOURCE_NAME "h2_excitation_uv_grb.cpp"

using namespace std;
const int VERBOSITY = 1;

// sim_path    - path to the folder with simulation data (GRB propagation results)
// output_path - path to the folder where data of the current simulations must be saved,
void init_time_intervals(vector<double>& time_arr);
void simulations_grb_uv(const string& data_path, const string& sim_path, const string& output_path);


void init_viewing_angles(const string& path, const string& fname, vector<double>& viewing_angle_grid, vector<string>& sub_folder_arr);

void init_h2_transition_list(const string& path, const energy_diagram* h2_di, const einstein_coeff* h2_einst,
	vector<transition>& transition_list_file);

// Here, H column density  of the cloud layer centre (from the front cloud boundary) must be = integer value * 1.e+18
void init_h2_popdensity_evol(const string& output_path, const vector<double>& layer_grid_coord, const vector<double>& hcolumn_density_coord,
	vector<double>& time_arr, vector< vector<dynamic_array_template<float>>>& h2_popdens_evol_arr, double conc_h_tot, double vangle, int nb_lev_h2, int nb_layers, 
	bool rotational_transition_mode);

void calc_line_luminosities(const string& data_path, const string& sim_path, const string& output_path, const string &sub_folder,
	const string& fname_viewing_angles, double hcolumn_density_cloud, double theta_jet_axis, bool is_jet_top_hat, bool rotational_transition_mode);


static int f_h2ev(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
	static_cast<h2_population_evolution_data*>(user_data)->f(t, y, ydot);
	return (0);
}


// Indebetouw R.et al., ApJ 619, pp. 931 - 938 (2005)
// AK/AV = 1/8.8 (RV = 3.1, Indebetouw et al., 2005), AV/NH = 5.3e-22 mag cm2 H-1, tau = A/1.086  (Draine, 2011); note, AK/AV = 1/8.5 (Draine, 2011)
// from this follows tau_K = 0.654e-22 * NH,
// wlength in um, return A(lambda)/AK, AK extinction in K band (2.164 um), approximation is valid only at 1.25-7.75 um
double get_extinction_ir(double wlength) {
	if (wlength > 7.75) {
		wlength = 7.75;
	}

	if (wlength > 1.25) {
		// Indebetouw et al. (2005)
		return 4.074 * pow(wlength, -2.22 + 1.21 * log10(wlength));
	}
	else {
		// Draine, Physics of the Interstellar and the Intergalactic medium, chapter 21, p. 239 (2011);
		// at lambda = 1.22 um, A(lambda)/AIC = 0.489
		// at lambda = 0.802 um, A(lambda)/AIC = 1.
		// first, I found A(lambda)/AIC using linear interpolation:
		double x = (0.489 + (1.22 - wlength) * (1. - 0.489) / (1.22 - 0.802));
			
		// A(1.22 um)/AK is multiplied by A(lambda)/AIC / A(1.22 um)/AIC
		x = 4.074 * pow(1.22, -2.22 + 1.21 * log10(1.22)) * x / 0.489;
		return x;
	}
}

// files in simulation folder that are necessary in simulations:
//    grb_cloud_parameters.txt
//    grb_layer_coordinates.txt
//    grb_specimen_conc.txt
//    grb_h2_levels.txt
//    grb_h2_popdens.txt

int main()
{
	bool is_jet_top_hat, rotational_transition_mode;
	int nb_processors(4), verbosity(1);
	double hcolumn_density_cloud, theta_jet_axis;
	
	string fname_viewing_angles;
	string data_path;    // path to the input_data folder with spectroscopic, collisional data and etc.
	string output_path;  // path to the folder where data of the current simulation must be saved, or with auxiliary data
	string sim_path;     // path to the folder with simulation data (GRB propagation results)
	string path, sub_folder;

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
	
	// the file with viewing angles, for which simulations were performed, is located in the simulation data directory,
	// sim_path is the directory with the simulation data of grb-cloud interaction
	sim_path = "C:/Users/Александр/Documents/Данные и графики/paper GRB H2 infrared emission ultra-violet/sim_grb_cloud_lowz_2026/";

	vector<double> viewing_angle_grid;
	vector<string> sub_folder_arr;

	fname_viewing_angles = "viewing_angles_d100.txt";
	init_viewing_angles(sim_path, fname_viewing_angles, viewing_angle_grid, sub_folder_arr);

#pragma omp parallel private(path, output_path)
	{
#pragma omp for schedule(dynamic, 1) nowait
		for (int i = 0; i < (int)viewing_angle_grid.size(); i++)
		{
			path = sim_path + sub_folder_arr[i];
			output_path = "../../output/" + sub_folder_arr[i];
			//simulations_grb_uv(data_path, path, output_path);
		}
	}
	
	sim_path = "C:/Users/Александр/Documents/Данные и графики/paper GRB H2 infrared emission ultra-violet/sim_grb_cloud_lowz_2026/";
	output_path = "../../output/";
	fname_viewing_angles = "viewing_angles_d10.txt";
	sub_folder = "d10_nh2e21/";

	// H column density
	hcolumn_density_cloud = 2.e+21;   // [cm-2]
	
	// if is_jet_top_hat = true, this data is used in top_hat jet simulations
	is_jet_top_hat = false;

	// ro-vibrational and pure rotational transitions are treated separately
	rotational_transition_mode = false;

	// angle between observation line of sight and jet axis,
	theta_jet_axis = 0.;   // [radians]
	
	calc_line_luminosities(data_path, sim_path, output_path, sub_folder, fname_viewing_angles, 
		hcolumn_density_cloud, theta_jet_axis, is_jet_top_hat, rotational_transition_mode);

	rotational_transition_mode = true;
	
	calc_line_luminosities(data_path, sim_path, output_path, sub_folder, fname_viewing_angles,
		hcolumn_density_cloud, theta_jet_axis, is_jet_top_hat, rotational_transition_mode);
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
	
	double rel_tol, model_time, model_time_step, model_time_out;

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
		model_time = 0.;

		// restart of the solver with new values of initial conditions,
		flag = CVodeReInit(cvode_mem, model_time, y);

		// update the members of user_data class,
		f_h2ev(model_time, y, ydot, (void*)user_data);

		time_nb = 1;
		// there is integration of the parameters over the time, this time step must be sufficiently small, 
		model_time_step = 1.;

		while ((model_time < time_arr.back()) || ((model_time > time_arr.back()) && (time_nb < (int)time_arr.size())))
		{
			flag = CV_SUCCESS;
			must_be_saved = false;

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
					
					must_be_saved = true;
					time_nb++;

					if (flag != CV_SUCCESS) {
						cout << "Error in CVodeGetDky() " << flag << endl;
					}
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
		// H column density from the cloud boundary to the cloud layer centre,
		// Here, H column density must be = integer value * 1.e+18
		hcolumn_density = conc_h_tot * (layer_centre_distances[lay_nb] - grb_cloud_distance);

		save_h2_populations_evolution(output_path, time_arr, h2_popdens_evol, conc_h_tot, hcolumn_density, nb_lev_h2);
		save_h2_populations_evolution_rotational(output_path, time_arr, h2_popdens_evol, conc_h_tot, hcolumn_density, nb_lev_h2);
		
		if (VERBOSITY) {
			i = (int)(time(NULL) - timer);
			cout << left << "Calc time (s): " << i << endl;
		}
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


void calc_line_luminosities(const string& data_path, const string& sim_path, const string& output_path, const string& sub_folder, 
	const string& fname_viewing_angles, double hcolumn_density_cloud, double theta_jet_axis, bool is_jet_top_hat, bool rotational_transition_mode)
{
	bool quit_cycle, is_transition_rotational;
	int i, j, k, v, nb, lay, lev, nb_layers, nb_lev_h2, nb_lines_h2, nb_times, nb_observ_times, nb_viewing_angles, isotope, jet_type;
	
	// parameters of the jet and molecular cloud,
	double x, h2_fraction, conc_h_tot, op_ratio_h2, cloud_width, grb_cloud_distance, max_hcolumn_density, half_opening_angle, 
		theta_wing, cos_theta_jet_axis, sin_theta_jet_axis;

	double obst_const, observ_time, max_observ_time, delta_observ_time, max_delta_time, luminosity, vangle, max_cos_theta, cos_theta_observ_time,
		 min_cos_theta_prime, max_cos_theta_prime, abs_const_kband, abs_const_ir, rint_const;
	
	double* abs_const_arr, * h2_lines_integr_theta, * h2_popdens_integr_theta;

	vector<double> layer_grid_coord, layer_centre_distances, hcolumn_density_coord, observ_time_arr, time_arr,
		viewing_angle_grid, cos_viewing_angle_grid;
	
	string fname, path;
	vector<string> sub_folder_arr;
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

	vector<transition> transition_list, transition_list_file;
	vector<transition>::iterator it_tl;

	//
	// reading file with viewing angles, for which simulations were performed,
	init_viewing_angles(sim_path, fname_viewing_angles, viewing_angle_grid, sub_folder_arr);
	nb_viewing_angles = (int) viewing_angle_grid.size();

	// reading file with the strongest transitions
	init_h2_transition_list(sim_path, h2_di, h2_einst, transition_list_file);
	
	nb_lines_h2 = (int) transition_list_file.size();
	sort(transition_list_file.begin(), transition_list_file.end());

	abs_const_arr = new double[nb_lines_h2];
	h2_lines_integr_theta = new double[nb_lines_h2];
	h2_popdens_integr_theta = new double[nb_lev_h2];
	
	dynamic_array h2_popdens_integr_volume(nb_lev_h2);
	dynamic_array h2_lines_integr_volume(nb_lines_h2);
	
	// vector index - obsevration time moment number, array index - H2 level number,
	vector<dynamic_array> h2_popdens_integr_volume_arr;
	// vector index - obsevration time moment number, array index - H2 line number,
	vector<dynamic_array> h2_lines_integr_volume_arr;

	// first vector index - layer number, second vector index - evolution time moment number, third array index - H2 level number,
	vector< vector<dynamic_array_template<float>> > h2_popdens_evol_arr_1, h2_popdens_evol_arr_2;

	//
	// reading data common to all simulation data, 
	path = sim_path + sub_folder_arr[0];

	init_cloud_parameters(path, grb_cloud_distance, max_hcolumn_density, conc_h_tot, op_ratio_h2, h2_fraction,
		jet_type, half_opening_angle, theta_wing, vangle);

	init_layer_coordinates(path, layer_grid_coord, hcolumn_density_coord);

	if (is_jet_top_hat) {
		// only one case with on-axis viewing angle is calculated, given the simulation data path
		// parameter "jet_type" is ignored,
		// "half_opening_angle" must be initialized,
		nb_viewing_angles = 2;

		viewing_angle_grid.clear();
		viewing_angle_grid.push_back(0.);
		viewing_angle_grid.push_back(half_opening_angle);
	}

	cloud_width = hcolumn_density_cloud / conc_h_tot;
	
	// tau = sigma NH is the optical depth for dust extinction, 
	// Bialy S. Communications Physics 3, 32 (2020), the reference to Draine B.T., "Physics of the Interstellar and Intergalactic medium" (2011),:
	//		sigma = 4.5e-23 cm2 is the cross section per hydrogen nucleus (the numeric value is an average over 2-3 um),
	// We use Indebetouw R. et al., ApJ 619, pp. 931-938 (2005),
	// H2 transition J=6->4 has 8.025 um,

	// cross section * n_{h,tot} for K band (2.164 um)
	abs_const_kband = 0.654e-22 * conc_h_tot;  // K band, absorption in continuum, does not depend on molecular fraction
	
	// recalculate to specific wavelength (2.5 um)
	abs_const_ir = abs_const_kband * get_extinction_ir(2.5);

	for (nb = 0; nb < nb_lines_h2; nb++) {
		abs_const_arr[nb] = abs_const_kband * get_extinction_ir(1.e+4 / transition_list_file[nb].energy);  // energy in [cm-1]
	}

	for (lay = 0; lay < (int)layer_grid_coord.size() - 1; lay++) {
		layer_centre_distances.push_back(0.5 * (layer_grid_coord[lay] + layer_grid_coord[lay + 1]));
	}

	// number of layers that are taken into in luminosity calculations,
	for (lay = 0; (lay < (int)layer_grid_coord.size() - 1) && (layer_grid_coord[lay] - grb_cloud_distance < cloud_width); lay++)
	{;}
	nb_layers = lay;

	cos_theta_jet_axis = cos(theta_jet_axis);
	sin_theta_jet_axis = sin(theta_jet_axis);

	// calculation of auxiliary integrals
	table_int.setp(cos_theta_jet_axis, sin_theta_jet_axis);

	for (v = 0; v < nb_viewing_angles; v++) {
		cos_viewing_angle_grid.push_back(cos(viewing_angle_grid[v]));
	}
	
	// the characteristic time of H2 spontaneous de-excitation is of the order of 1.e+7 s,
	max_observ_time = 1.e+7;
	max_observ_time += layer_grid_coord[nb_layers] * (1. - cos(viewing_angle_grid.back() + theta_jet_axis)) / SPEED_OF_LIGHT;  // in [s]
	
	if (max_observ_time > 30. * YEARS_TO_SECONDS) {
		max_observ_time = 30. * YEARS_TO_SECONDS;
	}
	
	if (max_observ_time < 10. * YEARS_TO_SECONDS) {
		delta_observ_time = 5.e+5;  // in [s]
	}
	else if (max_observ_time < 20. * YEARS_TO_SECONDS) {
		delta_observ_time = 10.e+5;  // in [s]
	}
	else {
		delta_observ_time = 15.e+5;  // in [s]
	}

	nb_observ_times = (int)(max_observ_time / delta_observ_time) + 1;
	for (i = 0; i < nb_observ_times; i++) {
		observ_time_arr.push_back(delta_observ_time * i);
	}

	// the main containers with integration data: 
	for (i = 0; i < nb_observ_times; i++) 
	{
		h2_popdens_integr_volume_arr.push_back(h2_popdens_integr_volume);
		h2_lines_integr_volume_arr.push_back(h2_lines_integr_volume);
	}

	// 
	// Initialization of data on population densities for viewing angle equal to zero,
	vangle = viewing_angle_grid[0];
	path = output_path + sub_folder_arr[0];

	init_h2_popdensity_evol(path, layer_grid_coord, hcolumn_density_coord, time_arr, h2_popdens_evol_arr_1, 
		conc_h_tot, vangle, nb_lev_h2, nb_layers, rotational_transition_mode);

	nb_times = (int)time_arr.size();

	//
	// the main loop
 	for (v = 1; v < nb_viewing_angles; v++) {	
		if (is_jet_top_hat) {
			for (lay = 0; lay < (int)h2_popdens_evol_arr_1.size(); lay++) {
				h2_popdens_evol_arr_2.push_back(h2_popdens_evol_arr_1[lay]);
			}
		}
		else {
			// reading the data on population densities for given viewing angle,
			vangle = viewing_angle_grid[v];
			path = output_path + sub_folder_arr[v];

			init_h2_popdensity_evol(path, layer_grid_coord, hcolumn_density_coord, time_arr, h2_popdens_evol_arr_2,
				conc_h_tot, vangle, nb_lev_h2, nb_layers, rotational_transition_mode);
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
			for (i = 0; i < nb_lines_h2; i++) {
				h2_lines_integr_volume.arr[i] = 0.;
			}

			//
			// integration by depth (or, equivalently, summation of cloud layers),
			for (lay = 0; lay < nb_layers; lay++)
			{
				memset(h2_popdens_integr_theta, 0, nb_lev_h2 * sizeof(double));
				memset(h2_lines_integr_theta, 0, nb_lines_h2 * sizeof(double));

				obst_const = SPEED_OF_LIGHT / layer_centre_distances[lay];
				rint_const = layer_centre_distances[lay] * layer_centre_distances[lay] * (layer_grid_coord[lay + 1] - layer_grid_coord[lay]);

				// upper limit in integration by cos theta, the inner border of the ring
				max_cos_theta = 1.;

				cos_theta_observ_time = 1. - obst_const * observ_time;
				max_delta_time = (max_cos_theta - cos_theta_observ_time) / obst_const;

				if (max_delta_time > time_arr.back()) {
					max_delta_time = time_arr.back();
					max_cos_theta = cos_theta_observ_time + obst_const * max_delta_time;
				}
				
				//
				// integration by theta, or equivalently, time,
#pragma omp parallel private(i, k, lev, x)
				{
					double x1, x2, tau;
					double cos_theta, cos_theta_1, cos_theta_2, sin_theta, sin_theta_sq, cos_theta_prime_1, cos_theta_prime_2, 
						integral_by_phi, linear_integral_by_phi, max_cos_phi, min_cos_phi, min_phi, max_phi;

					double * arr_levels = new double[nb_lev_h2];
					memset(arr_levels, 0, nb_lev_h2 * sizeof(double));

					double* arr_lines = new double[nb_lines_h2];
					memset(arr_lines, 0, nb_lines_h2 * sizeof(double));
	 
#pragma omp for schedule(dynamic, 1)
					for (k = 0; k < nb_times - 1; k++) 
					{
						const dynamic_array_template<float> & h2_popdens_1 = h2_popdens_evol_arr_1[lay][k];
						const dynamic_array_template<float> & h2_popdens_2 = h2_popdens_evol_arr_2[lay][k];

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
								}

								x = (cos_theta_2 - cos_theta_1) / (viewing_angle_grid[v] - viewing_angle_grid[v - 1])
									* exp(-tau * abs_const_ir);
									
								for (lev = 0; lev < nb_lev_h2; lev++) {
									x1 = (double) (h2_popdens_1.arr[lev]);
									x2 = (double) (h2_popdens_2.arr[lev]);

									// [cm-3], not divided by statistical weight,
									// average dust absorption for all volume-integrated population densities
									arr_levels[lev] += (integral_by_phi * (x1 * viewing_angle_grid[v] - x2 * viewing_angle_grid[v - 1])
										+ linear_integral_by_phi * (x2 - x1)) * x;
								}
 
								// individual absorption for each H2 line
								x = (cos_theta_2 - cos_theta_1) / (viewing_angle_grid[v] - viewing_angle_grid[v - 1]);

								for (i = 0; i < nb_lines_h2; i++) {
									lev = transition_list_file[i].up_lev.nb;

									x1 = (double)(h2_popdens_1.arr[lev]);
									x2 = (double)(h2_popdens_2.arr[lev]);

									arr_lines[i] += (integral_by_phi * (x1 * viewing_angle_grid[v] - x2 * viewing_angle_grid[v - 1])
										+ linear_integral_by_phi * (x2 - x1)) * x * exp(-tau * abs_const_arr[i]);
								}
							}
						}
					}
					
#pragma omp critical
					{
						for (lev = 0; lev < nb_lev_h2; lev++) {
							h2_popdens_integr_theta[lev] += arr_levels[lev];
						}
						for (i = 0; i < nb_lines_h2; i++) {
							h2_lines_integr_theta[i] += arr_lines[i];
						}
					}
					delete[] arr_levels;
					delete[] arr_lines;
				}
				
				for (lev = 0; lev < nb_lev_h2; lev++) {
					if (layer_grid_coord[lay + 1] - grb_cloud_distance < cloud_width) {
						h2_popdens_integr_volume.arr[lev] += h2_popdens_integr_theta[lev] * rint_const;
					}
					else {
						h2_popdens_integr_volume.arr[lev] += h2_popdens_integr_theta[lev]
							* layer_centre_distances[lay] * layer_centre_distances[lay] * (cloud_width + grb_cloud_distance - layer_grid_coord[lay]);
					}
				}

				for (i = 0; i < nb_lines_h2; i++) {
					if (layer_grid_coord[lay + 1] - grb_cloud_distance < cloud_width) {
						h2_lines_integr_volume.arr[i] += h2_lines_integr_theta[i] * rint_const;
					}
					else {
						h2_lines_integr_volume.arr[i] += h2_lines_integr_theta[i]
							* layer_centre_distances[lay] * layer_centre_distances[lay] * (cloud_width + grb_cloud_distance - layer_grid_coord[lay]);
					}
				}
			}
			for (lev = 0; lev < nb_lev_h2; lev++) {
				h2_popdens_integr_volume_arr[nb].arr[lev] += h2_popdens_integr_volume.arr[lev];
			}
			for (i = 0; i < nb_lines_h2; i++) {
				h2_lines_integr_volume_arr[nb].arr[i] += h2_lines_integr_volume.arr[i];
			}
		}
	
		if (!is_jet_top_hat) {
			h2_popdens_evol_arr_1.clear();
			
			for (lay = 0; lay < (int)h2_popdens_evol_arr_2.size(); lay++) {
				h2_popdens_evol_arr_1.push_back(h2_popdens_evol_arr_2[lay]);
			}
		}
	}

	// (first) vector index: transition number, (second) array index - time moment number,
	vector<dynamic_array> luminosity_evol_arr_estimated, luminosity_evol_arr, fline_luminosity_evol_arr;
	vector<dynamic_array>::iterator it_la;

	dynamic_array luminosity_evol(nb_observ_times);

	for (k = 0; k < (int)transition_list_file.size(); k++) {
		i = transition_list_file[k].up_lev.nb;
		j = transition_list_file[k].low_lev.nb;

		for (nb = 0; nb < nb_observ_times; nb++)
		{
			luminosity = h2_einst->arr[i][j] * (h2_di->lev_array[i].energy - h2_di->lev_array[j].energy) * CM_INVERSE_TO_ERG 
				* h2_lines_integr_volume_arr[nb].arr[k] / SOLAR_LUMINOSITY;

			luminosity_evol.arr[nb] = luminosity;
		}
		fline_luminosity_evol_arr.push_back(luminosity_evol);
	}

	// Search for high-luminosity lines,
	// loop over H2 ro-vibrational transitions:
	for (i = 1; i < nb_lev_h2; i++) {
		for (j = 0; j < i; j++) 
		{
			if (h2_einst->arr[i][j] > 1.e-99) {
				quit_cycle = false;

				is_transition_rotational = false;
				if (h2_di->lev_array[i].v == 0 && h2_di->lev_array[j].v == 0) {
					is_transition_rotational = true;
				}

				// for each observational time, the calculation of luminosity,
				for (nb = 0; (nb < nb_observ_times) && (!quit_cycle); nb++)
				{
					luminosity = h2_einst->arr[i][j] * (h2_di->lev_array[i].energy - h2_di->lev_array[j].energy) * CM_INVERSE_TO_ERG
						* h2_popdens_integr_volume_arr[nb].arr[i] / SOLAR_LUMINOSITY;

					if (luminosity > MIN_LUMINOSITY || (is_transition_rotational && luminosity > MIN_LUMINOSITY_ROT)) {
						transition_list.push_back(transition(h2_di->lev_array[j], h2_di->lev_array[i], h2_einst->arr[i][j]));
						quit_cycle = true;
					}
				}
			}
		}
	}
	sort(transition_list.begin(), transition_list.end());

	for (k = 0; k < (int)transition_list.size(); k++) {
		i = transition_list[k].up_lev.nb;
		j = transition_list[k].low_lev.nb;
		
		for (nb = 0; nb < nb_observ_times; nb++)
		{
			luminosity = h2_einst->arr[i][j] * (h2_di->lev_array[i].energy - h2_di->lev_array[j].energy) * CM_INVERSE_TO_ERG
				* h2_popdens_integr_volume_arr[nb].arr[i] / SOLAR_LUMINOSITY;

			luminosity_evol.arr[nb] = luminosity;
		}
		luminosity_evol_arr_estimated.push_back(luminosity_evol);
	}

	for (k = 0; k < (int)transition_list.size(); k++) {
		for (i = 0; i < (int)transition_list_file.size(); i++) {
			if (transition_list[k] == transition_list_file[i]) {
				break;
			}
		}

		if (i == (int)transition_list_file.size()) {
			luminosity_evol_arr.push_back(luminosity_evol_arr_estimated[k]);
		}
		else {
			luminosity_evol_arr.push_back(fline_luminosity_evol_arr[i]);
		}
	}

	path = output_path + sub_folder;
	save_h2_line_luminosities(path, observ_time_arr, 
		transition_list, luminosity_evol_arr, luminosity_evol_arr_estimated, theta_jet_axis, is_jet_top_hat, rotational_transition_mode);

	if (VERBOSITY) {
		i = (int)(time(NULL) - timer);
		cout << left << "Calc time (min): " << setw(7) << i / 60 << endl;
	}
}


//
void init_viewing_angles(const string& path, const string& fname_viewing_angles, vector<double> & viewing_angle_grid, vector<string> & sub_folder_arr)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, nb_viewing_angles;
	double x;

	string str, fname;
	ifstream input;

	fname = path + fname_viewing_angles;
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_viewing_angles;

	for (i = 0; i < nb_viewing_angles; i++) {
		input >> x >> str;

		viewing_angle_grid.push_back(x);
		sub_folder_arr.push_back(str);
	}
	input.close();
}


//
void init_h2_transition_list(const string& path, const energy_diagram *h2_di, const einstein_coeff* h2_einst,
	vector<transition> & transition_list_file)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int nb, i, j, k, v, nb_lines_h2;
	double x;

	string fname;
	ifstream input;

	fname = path + "h2_line_list.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines_h2;

	for (nb = 0; nb < nb_lines_h2; nb++) {
		input >> k >> v >> j;
		i = h2_di->get_nb(v, j);  // upper energy level

		input >> v >> j >> x;
		k = h2_di->get_nb(v, j);  // lower energy level

		if (i != -1 && k != -1) {
			transition_list_file.push_back(transition(h2_di->lev_array[k], h2_di->lev_array[i], h2_einst->arr[i][k]));
		}
	}
	input.close();
}


// Here, H column density must be = integer value * 1.e+18
void init_h2_popdensity_evol(const string & path, const vector<double> & layer_grid_coord, const vector<double> & hcolumn_density_coord, 
	vector<double> & time_arr, vector< vector<dynamic_array_template<float>>> & h2_popdens_evol_arr, double conc_h_tot, double vangle, int nb_lev_h2, int nb_layers, 
	bool rotational_transition_mode)
{
	char ch, text_line[MAX_TEXT_LINE_WIDTH];
	int i, k, lev, lay, nb_times, nb_lev_f;
	double x;

	string fname;
	ifstream input;
	stringstream ss;
	
	dynamic_array_template<float> h2_popdens(nb_lev_h2);
	vector< dynamic_array_template<float> > h2_popdens_evol;

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

		if (rotational_transition_mode) {
			fname = path + "rot_h2grb_popdens";
		}
		else {
			fname = path + "h2grb_popdens";
		}
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
				h2_popdens.arr[lev] = (float) x;
			}
			h2_popdens_evol.push_back(h2_popdens);
		}
		input.close();

		// averaging of population densities over two time moments
		for (k = 0; k < nb_times - 1; k++) {
			for (lev = 0; lev < nb_lev_h2; lev++) {
				h2_popdens_evol[k].arr[lev] = (float) 0.5 * (h2_popdens_evol[k].arr[lev] + h2_popdens_evol[k + 1].arr[lev]);
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
*/

/*

	for (k = 0; k < (int)transition_list.size(); k++) {
		for (i = 0; i < (int)transition_list_file.size(); i++) {
			if (transition_list[k] == transition_list_file[i]) {
				break;
			}
		}
		// adding the transition to the list:
		if (i == (int)transition_list_file.size()) {
			if (VERBOSITY) {
				cout << left << "Adding transition to the orginal list: "
					<< setw(5) << transition_list[k].up_lev.v << setw(5) << transition_list[k].up_lev.j << " -> "
					<< setw(5) << transition_list[k].low_lev.v << setw(5) << transition_list[k].low_lev.j << endl;
			}

			for (i = 0; i < (int)transition_list_file.size(); i++) {
				if (transition_list[k].energy < transition_list_file[i].energy)
				{
					it_tl = transition_list_file.begin() + i;
					transition_list_file.insert(it_tl, transition_list[k]);

					it_la = fline_luminosity_evol_arr.begin() + i;
					fline_luminosity_evol_arr.insert(it_la, luminosity_evol_arr[k]);
					break;
				}
			}
		}
	}

	for (i = 0; i < (int)transition_list_file.size(); i++) {
		for (k = 0; k < (int)transition_list.size(); k++) {
			if (transition_list[k] == transition_list_file[i]) {
				break;
			}
		}
		// removing i-element
		if (k == (int)transition_list.size()) {
			it_tl = transition_list_file.begin() + i;
			transition_list_file.erase(it_tl);

			it_la = fline_luminosity_evol_arr.begin() + i;
			fline_luminosity_evol_arr.erase(it_la);
			i--;
		}
	}
*/
