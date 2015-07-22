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

// Auxiliary variables
int mAGE;
// Output arrays
double ***mn_age;
// Working arrays
double ***age_init;
double *age_inventory;
double ***htest;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];
extern struct vardesc vars[NOVARS];
extern struct parameters run_parameters;


void allocate_age (  ) {
	printf("Allocating age arrays\n");

	int i, j, k;

	// Set index in tracer array
	mAGE = run_parameters.tracer_counter++;

	// Allocate working and output arrays
	mn_age = alloc3d(NZ,NXMEM,NYMEM);
	age_init = alloc3d(NZ,NXMEM,NYMEM);
	htest = alloc3d(NZ,NXMEM,NYMEM);
}

void initialize_age( ) {
	int i,j,k;
//	extern char restart_filename[200];
//	extern double misval;
	char varname[100];

	printf("Setting age variable description...");
	strcpy(varname,"age");
        strcpy(vars[map_variable_to_index(varname)].name,"mn_age");
        strcpy(vars[map_variable_to_index(varname)].longname,"Ideal age");
        vars[map_variable_to_index(varname)].hor_grid='h';
        vars[map_variable_to_index(varname)].z_grid='L';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"Years");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;
	printf("DONE\n");

//	printf("Age init example: ");
//	printf("%f\n",tr[mAGE][10][100][100]);
	if (run_parameters.restart_flag) {
		printf("Initializing age from restart %s\n",run_parameters.restartfile);
		read_var3d( run_parameters.restartfile, "mn_age", 0, age_init );
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
		set_darray3d_zero(age_init, NZ, NXMEM, NYMEM);
		printf("Ideal age initialized to zero\n");
	}

	wrap_reentrance_3d(age_init,NZ);
	

	// Copy the initialized tracer value over to main trace array
	for (i=0;i<NXMEM;i++) 
		for (j=0;j<NYMEM;j++)
			for (k=0;k<NZ;k++) {
				if (oceanmask[i][j]) tr[mAGE][k][i][j] = age_init[k][i][j];
				else tr[mAGE][k][i][j] = 0.0;
			}
	
	copy_darray3d(mn_age,tr[mAGE],NZ,NXMEM,NYMEM);
	printf("Age init example: %f\n",tr[mAGE][15][100][100]);
	free3d(age_init,NZ);


}

void step_age( double dt ){

	int i,j,k;
	const double numsecsinyear = 365.0*86400.0;

//	calc_inventory(mAGE);

	printf("Age source: %f\n",dt/numsecsinyear);
	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			if(oceanmask[i][j]) {
				// Set the mixed layer to zero
				tr[mAGE][0][i][j] = 0.0;
				tr[mAGE][1][i][j] = 0.0;
				for(k=2;k<NZ;k++)
					tr[mAGE][k][i][j] += dt/numsecsinyear;
			}


/*	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			for(k=0;k<NZ;k++)
				mn_age[k][i][j]+=tr[mAGE][k][i][j];
*/
}
