#include "init.h"

void initialize( void );

struct parameters { 

	char namelist_file[100];

	int syear, sinterval; // Model start time
	int eyear, einterval; // Model end time
	int use_hindcast; // 0 for normalyear, 1 for hindcast
	int restart_flag; // '1' for restart, '0' for cold start
	int adjoint_integration; // 1 for adjoint
	int tracer_counter; // Use this to keep track of total number of tracers;
	char forcing_path[1000]; // Directory to forcing
	char normalyear_path[1000]; // Directory to forcing
	char hindcast_path[1000]; // Directory to forcing
	char outputfile[100]; // model output name
	char restartfile[100]; // restart file name 
	char new_restartfile[100]; // restart file name
	char timestep[7]; // '5day' or 'month'
	char woa_path[1000]; // Path to WOA files
	char adjoint_initfile[1000]; // Path to adjoint initialization file
	
	int wrint; // Output interval
	int do_averaging; // Whether the output is averaged or not
	int ntime_climatology; // Number of time stamps in climatology forcing
	int ntime_hindcast; // Number of timesteps in hindcast fields

	int conservation_check; // Enables the output of htest and a TTD-like tracer to check conservaiton

	// Tracer flags
	int do_gasex; // Enable gas exchange for all gas tracers
	int do_age;
	int do_cfcs;
	int do_ttd;
	int num_ttd_intervals;
	int do_n2;
	int do_ar;
	int do_oxygen;
	int do_adjttd;
	int do_age_pair;
  	int adjttd_restart;
};

void set_run_parameters( );
