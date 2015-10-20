#include <stdio.h>
#include <math.h>
#include <string.h>
#include "init.h"
#include "netcdf.h"
#include "initialize.h"
#include "alloc.h"
#include "init.h"
#include "io.h"
#include "util.h"
#include "output_variables.h"
#include "read.h"
#include "timekeeper.h"

// Auxiliary variables
int mADJOINT;
// Output arrays
double **mn_adjttd;
// Working arrays
double ***adjttd_init;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];
extern struct vardesc vars[NOVARS];
extern struct parameters run_parameters;
extern struct timekeeper_t timekeeper;

void allocate_adjttd(  ) {
	printf("Allocating adjoint TTD arrays\n");

	// Allocate working and output arrays
	mn_adjttd = alloc2d(NXMEM,NYMEM);
	adjttd_init = alloc3d(NZ,NXMEM,NYMEM);
}

void initialize_adjttd( ) {
	int i,j,k;
//	extern char restart_filename[200];
//	extern double misval;
	char varname[100];

	printf("Adjoint TTD index in main tracer array: %d\n",mADJOINT);
	printf("Setting adjoint TTD variable description...");
	strcpy(varname,"adjttd");
        strcpy(vars[map_variable_to_index(varname)].name,"mn_adjttd");
        strcpy(vars[map_variable_to_index(varname)].longname,"Adjoint TTD from interior impulse");
        vars[map_variable_to_index(varname)].hor_grid='h';
        vars[map_variable_to_index(varname)].z_grid='1';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"none");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;

        strcpy(varname,"adjttd_restart");
        strcpy(vars[map_variable_to_index(varname)].name,"adjttd_restart");
        strcpy(vars[map_variable_to_index(varname)].longname,"Interior values to initialize from interior impulse");
        vars[map_variable_to_index(varname)].hor_grid='h';
        vars[map_variable_to_index(varname)].z_grid='L';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"none");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;
        printf("DONE\n");

	if (run_parameters.restart_flag) {
		printf("Initializing adjoint TTD from restart %s\n",run_parameters.restartfile);
		read_var3d( run_parameters.restartfile, "adjttd_restart", 0, adjttd_init );

	}
	else {
		read_var3d( run_parameters.adjoint_initfile, "adjttd_init", 0, adjttd_init );
		printf("Adjoint TTD initialized from file\n");
	}

	wrap_reentrance_3d(adjttd_init,NZ);
	

	// Copy the initialized tracer value over to main trace array
	for (i=0;i<NXMEM;i++) 
		for (j=0;j<NYMEM;j++)
			for (k=0;k<NZ;k++) {
				if (oceanmask[i][j]) tr[mADJOINT][k][i][j] = adjttd_init[k][i][j];
				else tr[mADJOINT][k][i][j] = 0.0;
			}

}

void step_adjttd(  ){

	// The adjoint TTD from the interior to the surface is the flux of the BIR from the oceanic interior into the mixed layer
	// To do this, we allow tracer transport to happen, make sure that the mixed layer is homogeneous in tracer, store the
	// concentrations in the mixed layer, then set the surface to zero.

	int i,j,k;

	// Inject adjoint tracer at specified locations if necessary
	if (timekeeper.iteration_counter < run_parameters.num_ttd_intervals && !run_parameters.restart_flag)
	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			for(k=0;k<NZ;k++) {
				if (adjttd_init[k][i][j])
					tr[mADJOINT][k][i][j] = 1.0;
			}

	// Storing the concentration from mixed layer into the averaging array
	submit_for_averaging_2d(mn_adjttd,tr[mADJOINT][0]);

	// Setting the mixed layer to 0
	for (k=0;k<NML;k++)
		set_darray2d_zero(tr[mADJOINT][k],NXMEM,NYMEM);

}
