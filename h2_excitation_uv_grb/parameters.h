#pragma once

const std::string chemical_species[] = { "e-", "H", "H+", "H2", "H2+", "He", "He+", "He++",
    "CI",  "CII",  "CIII",  "CIV",  "CV",  "CVI",  "CVII",
    "NI",  "NII",  "NIII",  "NIV",  "NV",  "NVI",  "NVII",  "NVIII",
    "OI",  "OII",  "OIII",  "OIV",  "OV",  "OVI",  "OVII",  "OVIII",  "OIX",
    "NeI", "NeII", "NeIII", "NeIV", "NeV", "NeVI", "NeVII", "NeVIII", "NeIX", "NeX", "NeXI",
    "MgI", "MgII", "MgIII", "MgIV", "MgV", "MgVI", "MgVII", "MgVIII", "MgIX", "MgX", "MgXI", "MgXII", "MgXIII",
    "SiI", "SiII", "SiIII", "SiIV", "SiV", "SiVI", "SiVII", "SiVIII", "SiIX", "SiX", "SiXI", "SiXII", "SiXIII", "SiXIV", "SiXV",
    "SI",  "SII",  "SIII",  "SIV",  "SV",  "SVI",  "SVII",  "SVIII",  "SIX",  "SX",  "SXI",  "SXII",  "SXIII",  "SXIV",  "SXV",  "SXVI",  "SXVII",
    "FeI", "FeII", "FeIII", "FeIV", "FeV", "FeVI", "FeVII", "FeVIII", "FeIX", "FeX", "FeXI", "FeXII", "FeXIII", "FeXIV", "FeXV", "FeXVI", "FeXVII",
    "FeXVIII", "FeXIX", "FeXX", "FeXXI", "FeXXII", "FeXXIII", "FeXXIV", "FeXXV", "FeXXVI", "FeXXVII"
};

#define NB_OF_CHEM_SPECIES 8


//-------------------------------------------------
// Numerical parameters
//-------------------------------------------------

#define REL_ERROR_SOLVER 1.e-7
#define ABS_CONCENTRATION_ERROR_SOLVER 1.e-14  // specimen concentration (cm-3)
#define ABS_POPULATION_ERROR_SOLVER 1.e-14     // specimen level population density (cm-3)
#define ABS_ELSPECTRA_ERROR_SOLVER 1.e-14
#define ABS_PARAMETER_ERROR_SOLVER 1.e-11  // e.g., grain radius (cm), temperature (K)
#define MAX_CONV_FAILS_SOLVER 100  // default value is 10;
#define MAX_ERR_TEST_FAILS_SOLVER 14  // default value is 7;


#define MAX_NB_STEPS 15
#define MINIMAL_ABUNDANCE 1.e-99     // for saving in file

// is used in calculations of the time evolution of level populations, the simulations are done up to this time,
// is used in saving the time evolution of populations of energy levels with v = 0,
// decay time of (v,J) = (0,6) is 3.8e+7 s, (0,7) is 1.7e+7 s, (1,0) 1.2e+6 s
#define MAX_MODEL_TIME 1.e+9   

// is used in saving the time evolution of populations of all ro-vibrational levels (v, J),
#define MAX_SAVING_TIME 3.e+7  // in s


#define MIN_MODEL_TIME 1.e+2   // in s, minimal model time at which data are saved,
#define MODEL_TIME_1 1.e+3
#define MODEL_TIME_2 1.e+4
#define MODEL_TIME_3 1.e+5
#define MODEL_TIME_4 1.e+6
#define MODEL_TIME_5 1.e+7

#define NB_OF_BINS_PER_ORDER_TIME_1 16
#define NB_OF_BINS_PER_ORDER_TIME_2 32
#define NB_OF_BINS_PER_ORDER_TIME_3 64
#define NB_OF_BINS_PER_ORDER_TIME_4 128
#define NB_OF_BINS_PER_ORDER_TIME_5 128
#define NB_OF_BINS_PER_ORDER_TIME_6 128     // check the dependence of this parameter on populations of (0,J) energy levels,

// minimal luminosity in the line to save the data (ro-vibrational)
#define MIN_LUMINOSITY 5.       // in solar luminosities
#define MIN_LUMINOSITY_ROT 1.
