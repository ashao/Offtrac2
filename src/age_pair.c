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

// Auxiliary variables
int mAGE1, mAGE2;
double age_pair_decay;
double age_pair_conversion;
// Output arrays
double ***mn_age_pair1;
double ***mn_age_pair2;

// Working arrays
double ***age_pair1_init;
double ***age_pair2_init;
double *age1_inventory, *age2_inventory;

extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];
extern struct vardesc vars[NOVARS];
extern struct parameters run_parameters;


void allocate_age_pair (  ) {
  printf("Allocating age arrays\n");
  // Allocate working and output arrays
  mn_age_pair1 = alloc3d(NZ,NXMEM,NYMEM);
  mn_age_pair2 = alloc3d(NZ,NXMEM,NYMEM);
  age_pair1_init = alloc3d(NZ,NXMEM,NYMEM);
  age_pair2_init = alloc3d(NZ,NXMEM,NYMEM);
}

void initialize_age_pair( ) {
  int i,j,k;
  //	extern char restart_filename[200];
  //	extern double misval;
  char varname[100];

  printf("Ideal age pairs in main tracer array: %d and %d\n",mAGE1, mAGE2);
  printf("Setting age pair tracers variable description...");
  strcpy(varname,"age_pair1");
  strcpy(vars[map_variable_to_index(varname)].name,"mn_age_pair1");
  strcpy(vars[map_variable_to_index(varname)].longname,"Age pair 1");
  vars[map_variable_to_index(varname)].hor_grid='h';
  vars[map_variable_to_index(varname)].z_grid='L';
  vars[map_variable_to_index(varname)].t_grid='s';
  strcpy(vars[map_variable_to_index(varname)].units,"Years");
  vars[map_variable_to_index(varname)].mem_size='d';
  vars[map_variable_to_index(varname)].mval=MISVAL;

  strcpy(varname,"age_pair2");
  strcpy(vars[map_variable_to_index(varname)].name,"mn_age_pair2");
  strcpy(vars[map_variable_to_index(varname)].longname,"Age pair 2");
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
      read_var3d( run_parameters.restartfile, "mn_age_pair1", 0, age_pair1_init );
      read_var3d( run_parameters.restartfile, "mn_age_pair2", 0, age_pair2_init );
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
      // Age pair 1 is zero, age pair 2 is 1/age_pair_decay
      set_darray3d_zero(age_pair1_init, NZ, NXMEM, NYMEM);
      for (i=0;i<NXMEM;i++)
        for (j=0;j<NYMEM;j++)
          for (k=0;k<NZ;k++) {
              age_pair2_init[k][i][j] = 1.0/age_pair_decay;
      }
      printf("Age pairs initialized to %1.2e and %1.2e\n", 0.0, 1.0/age_pair_decay);
  }

  wrap_reentrance_3d(age_pair1_init,NZ);
  wrap_reentrance_3d(age_pair2_init,NZ);


  // Copy the initialized tracer value over to main trace array
  for (i=0;i<NXMEM;i++)
    for (j=0;j<NYMEM;j++)
      for (k=0;k<NZ;k++) {
	  if (oceanmask[i][j]) {
	      tr[mAGE1][k][i][j] = age_pair1_init[k][i][j];
	      tr[mAGE2][k][i][j] = age_pair2_init[k][i][j];
	  }
	  else {
	      tr[mAGE1][k][i][j] = 0.0;
	      tr[mAGE2][k][i][j] = 0.0;
	  }
      }

  copy_darray3d(mn_age_pair1,tr[mAGE1],NZ,NXMEM,NYMEM);
  copy_darray3d(mn_age_pair2,tr[mAGE2],NZ,NXMEM,NYMEM);
  printf("Age pair 1 init example: %f\n",tr[mAGE1][15][100][100]);
  printf("Age pair 2 init example: %f\n",tr[mAGE2][15][100][100]);
  free3d(age_pair1_init,NZ);
  free3d(age_pair2_init,NZ);


}

void step_age_pair( double dt ){

  int i,j,k;

  for(i=0;i<NXMEM;i++)
    for(j=0;j<NYMEM;j++)
      if(oceanmask[i][j]) {
	  // Set the mixed layer to zero
	  tr[mAGE1][0][i][j] = 0.0;
	  tr[mAGE1][1][i][j] = 0.0;
	  tr[mAGE2][0][i][j] = 0.0;
	  tr[mAGE2][1][i][j] = 0.0;

	  for(k=2;k<NZ;k++) {
	      tr[mAGE1][k][i][j] += dt/age_pair_conversion-age_pair_decay*tr[mAGE1][k][i][j];
	      tr[mAGE2][k][i][j] += dt/age_pair_conversion-age_pair_decay*tr[mAGE2][k][i][j];

	  }

      }


  /*	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			for(k=0;k<NZ;k++)
				mn_age[k][i][j]+=tr[mAGE][k][i][j];
   */
}
