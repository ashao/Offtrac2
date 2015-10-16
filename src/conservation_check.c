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

// Auxiliary variables
int mTEST;
// Output arrays
double ***mn_test;
// Working arrays
double test_inventory;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];
extern struct vardesc vars[NOVARS];
extern struct parameters run_parameters;


void allocate_test (  ) {
	printf("Allocating age arrays\n");

	int i, j, k;

	// Allocate working and output arrays
	mn_test = alloc3d(NZ,NXMEM,NYMEM);
	mn_htest = alloc3d(NZ,NXMEM,NYMEM);
}

void initialize_test( ) {
	int i,j,k;
//	extern char restart_filename[200];
//	extern double misval;
	char varname[100];

	printf("Setting test variable description...");
	strcpy(varname,"test");
        strcpy(vars[map_variable_to_index(varname)].name,"mn_test");
        strcpy(vars[map_variable_to_index(varname)].longname,"Test tracer");
        vars[map_variable_to_index(varname)].hor_grid='h';
        vars[map_variable_to_index(varname)].z_grid='L';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"none");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;

	strcpy(varname,"test_inventory");
        strcpy(vars[map_variable_to_index(varname)].name,"test_inventory");
        strcpy(vars[map_variable_to_index(varname)].longname,"Test tracer inventory");
        vars[map_variable_to_index(varname)].hor_grid='1';
        vars[map_variable_to_index(varname)].z_grid='n';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"none");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;

	strcpy(varname,"htest");
        strcpy(vars[map_variable_to_index(varname)].name,"htest");
        strcpy(vars[map_variable_to_index(varname)].longname,"hend-hnew");
        vars[map_variable_to_index(varname)].hor_grid='h';
        vars[map_variable_to_index(varname)].z_grid='L';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"meters");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;
	printf("DONE\n");

	
	// Copy the initialized tracer value over to main trace array
	for (i=0;i<NXMEM;i++) 
		for (j=0;j<NYMEM;j++) {
			for (k=0;k<NML;k++)	tr[mTEST][k][i][j] = 1.0;	
			for (k=NML;k<NZ;k++)	tr[mTEST][k][i][j] = 0.0;
		}
}

void step_test( double dt ){

	int i,j,k;
	extern struct timekeeper_t timekeeper;
	if (timekeeper.iteration_counter < 2){
		printf("Injecting tracer into mixed layer\n");
		for(i=0;i<NXMEM;i++)
			for(j=0;j<NYMEM;j++)
				for(k=0;k<NML;k++)
					tr[mTEST][k][i][j] = 1.;
	}
			
//	test_inventory=calc_inventory(tr[mTEST]);
	//printf("Test tracer inventory: %e\n", calc_inventory(tr[mTEST]));

}

