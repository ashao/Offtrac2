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
#include "timekeeper.h"
#include "read.h"

// Auxiliary variables
int mTTD;
// Output arrays
double ***mn_ttd;
// Working arrays
double ***ttd_init;
double *ttd_inventory;
double ***htest;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];
extern struct vardesc vars[NOVARS];
extern struct parameters run_parameters;
extern struct timekeeper_t timekeeper;

void allocate_ttd (  ) {
	printf("Allocating TTD  arrays\n");

	// Set index in tracer array

	// Allocate working and output arrays
	mn_ttd = alloc3d(NZ,NXMEM,NYMEM);
	ttd_init = alloc3d(NZ,NXMEM,NYMEM);
//	htest = alloc3d(NZ,NXMEM,NYMEM);
}

void initialize_ttd( ) {
	int i,j,k;
//	extern char restart_filename[200];
//	extern double misval;
	char varname[100];

	printf("TTD index in main tracer array: %d\n",mTTD);
	printf("Setting TTD variable description...");
	strcpy(varname,"ttd");
        strcpy(vars[map_variable_to_index(varname)].name,"mn_ttd");
        strcpy(vars[map_variable_to_index(varname)].longname,"TTD via boundary propagator");
        vars[map_variable_to_index(varname)].hor_grid='h';
        vars[map_variable_to_index(varname)].z_grid='L';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"none");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;
	printf("DONE\n");

//	printf("Age init example: ");
//	printf("%f\n",tr[mAGE][10][100][100]);
	if (run_parameters.restart_flag) {
		printf("Initializing age from restart %s\n",run_parameters.restartfile);
		read_var3d( run_parameters.restartfile, "mn_ttd", 0, ttd_init );
		// Filter for bad values
		//
/*
		for (i=0;i<NXMEM;i++) 
			for (j=0;j<NYMEM;j++)
				for (k=0;k<NZ;k++)
					if(age_init[k][i][j] > 2050. || age_init[k][i][j] < 0 ) {
						age_init[k][i][j] = 0.0;
					}
		
*/		
	}
	else {
		set_darray3d_zero(ttd_init,NZ,NXMEM,NYMEM);	
		for (i=0;i<NXMEM;i++) 
			for (j=0;j<NYMEM;j++)
				for (k=0;k<NML;k++)
					ttd_init[k][i][j] = 1.0;
	}

	wrap_reentrance_3d(ttd_init,NZ);
	

	// Copy the initialized tracer value over to main trace array
	for (i=0;i<NXMEM;i++) 
		for (j=0;j<NYMEM;j++)
			for (k=0;k<NZ;k++) {
				if (oceanmask[i][j]) tr[mTTD][k][i][j] = ttd_init[k][i][j];
				else tr[mTTD][k][i][j] = 0.0;
			}
	
	copy_darray3d(mn_ttd,tr[mTTD],NZ,NXMEM,NYMEM);
	printf("TTD init example: %f\n",tr[mTTD][15][100][100]);
	free3d(ttd_init,NZ);


}

void step_ttd(  ){

	int i,j;

	if (timekeeper.iteration_counter < run_parameters.num_ttd_intervals) {
		for(i=0;i<NXMEM;i++)
			for(j=0;j<NYMEM;j++)
				if(oceanmask[i][j]) {
					// Set the mixed layer to zero
					tr[mTTD][0][i][j] = 1.0;
					tr[mTTD][1][i][j] = 1.0;
				}
	}
	else {
		for(i=0;i<NXMEM;i++)
			for(j=0;j<NYMEM;j++)
				if(oceanmask[i][j]) {
					// Set the mixed layer to zero
					tr[mTTD][0][i][j] = 0.0;
					tr[mTTD][1][i][j] = 0.0;
				}
	}


}
