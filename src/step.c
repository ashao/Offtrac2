/********+*********+*********+*********+*********+*********+*********+*
 *      							      *
 *                source/sink terms for offtrac tracers               * 
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "init.h"
#include "metrics.h"
#include "alloc.h"
#include "step.h"
#include "tracadv.h"
#include "util.h"
#include "tracer_utilities.h"
#include "output_variables.h"
#include "timekeeper.h"
#include "initialize.h"

#ifdef AGE
#include "ideal_age.h"
#endif
/*---------------------------------------------------------------------
 *     define variables and subroutines
 *---------------------------------------------------------------------*/

extern double ****tr;
extern double ***hstart, ***h, ***hend;

#ifdef AGE
extern int mAGE;
extern double ***mn_h, ***h;
extern double ***mn_uhtm, ***uhtm;
extern double ***mn_vhtm, ***vhtm;
extern double ***mn_wd, ***wd;
#endif 

extern int oceanmask[NXMEM][NYMEM];
extern struct parameters run_parameters;
extern struct timekeeper_t timekeeper;


/*---------------------------------------------------------------------
 *
 *     BEGIN MAIN EXECUTABLE
 *        Each tracer is stepped forward in 3 consecutive steps:
 *         1) physical transport
 *         2) biological/chemical interior source-sink
 *         3) fluxes across air-sea interface
 *         
 *---------------------------------------------------------------------*/

void step_fields( ) {
	int i, j, k, l, m;
	int ii;

	const double secperyr = 365.25 * 60 * 60 * 24;

	/*-----------------------------------------
	 *
	 *     PARAMETER VALUES
	 *
	 *-----------------------------------------*/

	/*-----------------------------------------
	 *
	 *     CALCULATE PRELIMINARY QUANTITIES
	 *
	 *-----------------------------------------*/

	/*-----------------------------------------
	 *
	 *     PRINT TO SCREEN #1
	 *
	 *-----------------------------------------*/
	/*-----------------------------------------
	 *
	 *     STEP 1: PHYSICAL TRANSPORT
	 *
	 *-----------------------------------------*/

	update_transport_fields( );

	printf("Calculate tracer transport. \n");
	printf("TR[0][15][100][100]: %f\n",tr[0][15][100][100]);
	tracer( 0  ); /* perform transport time step */
	merge_ml_tr();

	/*-----------------------------------------
	 *
	 *     STEP 2: CALCULATE SOURCE/SINK TERMS
	 *
	 *-----------------------------------------*/


	/*-----------------------------------------
	 *
	 *  APPLY BIOTIC SOURCE/SINK TERMS
	 *   step tr array forward using Euler method
	 *   calculate global integral for mass conservation checking
	 *
	 *-----------------------------------------*/


#ifdef AGE
	step_age(timekeeper.dt);
#endif

	merge_ml_tr();





	/*-----------------------------------------
	 *
	 *     OVERWRITE RESULTS IN SPECIAL CASES
	 *
	 *-----------------------------------------*/

	/*-----------------------------------------
	 *
	 *    Sponges and Re-Entrant Boundaries
	 *
	 *-----------------------------------------*/

	/* zonal, meridional re-entrance    */
	for (m = 0; m < NTR; m++) {
		for (k = 0; k < NZ; k++) {
			for (j = 0; j <= NYMEM - 1; j++) {
				tr[m][k][0][j] = tr[m][k][nx - 1][j];
				tr[m][k][1][j] = tr[m][k][nx][j];
				tr[m][k][nx + 1][j] = tr[m][k][2][j];
				tr[m][k][nx + 2][j] = tr[m][k][3][j];
			}
		}
		for (i = 2; i <= nx; i++) {
			ii = 363 - i;
			for (k = 0; k < NZ; k++) {
				tr[m][k][ii][ny + 1] = tr[m][k][i][ny];
				tr[m][k][ii][ny + 2] = tr[m][k][i][ny - 1];
			}
		}
	}


	for (m=0;m<NTR;m++)
		for(i=0;i<NXMEM;i++)
			for(j=0;j<NYMEM;j++)
				if (!oceanmask[i][j])
					for(k=0;k<NZ;k++)
						tr[m][k][i][j]=MISVAL;
					

	submit_for_averaging( mn_h, h );
	submit_for_averaging( mn_uhtm, uhtm );
	submit_for_averaging( mn_vhtm, vhtm );
	submit_for_averaging( mn_wd, wd );

#ifdef AGE
	submit_for_averaging( mn_age, tr[mAGE]) ;
#endif

//	apply_mask(mn_h,oceanmask);
//	apply_mask(mn_uhtm,oceanmask);
//	apply_mask(mn_vhtm,oceanmask);
//	apply_mask(mn_wd,oceanmask);
	

	/*-----------------------------------------
	 *
	 *     PRINT TO SCREEN #3
	 *
	 *-----------------------------------------*/

} /* end time step */

/*---------------------------------------------------------------------
 *
 *     BEGIN SUBROUTINES
 *
 *---------------------------------------------------------------------*/


double tracer_integral(int trnum, double ***hvol) {

	int i, j, k;
	double trintegral;
	trintegral = 0.e0;
	for (k = 0; k < NZ; k++) {
		for (i = 2; i <= NXMEM - 1; i++) {
			for (j = 2; j <= NYMEM - 1; j++) {
				if (hvol[k][i][j] > 1e-6) {
					trintegral += tr[trnum][k][i][j] * hvol[k][i][j];
				}
			}
		}
	}
	return trintegral;

}

void merge_ml_tr( void ) {

	int i, j, k, l;
	double tracsum;

	for(l=0;l<NTR;l++)
		for(i=0;i<NXMEM;i++)
			for(j=0;j<NYMEM;j++) {

				tracsum = 0; // Reset the tracer amount
				for(k=0;k<NML;k++) {
					tracsum+=tr[l][k][i][j]/NML;
				}

				// Set the mixed layer to the average of all the mixed layer 'layers'
				for(k=0;k<NML;k++) {
					tr[l][k][i][j]=tracsum/NML;
				}


			}

}

void update_transport_fields(  ) {

	char file_suffix[100];
	char path[1000];
	int read_index;

	copy_darray3d(hstart,hend,NZ,NXMEM,NYMEM);
	copy_darray3d(h,hstart,NZ,NXMEM,NYMEM);

	if (timekeeper.read_hind_flag) {
		// BUG FIX BECAUSE 1964 is INCOMPLETE REMOVE WHEN HINDCAST RERUN!
		if (timekeeper.current_year == 1964) 
			sprintf(file_suffix,"hind.%d.%s",1965,run_parameters.timestep);
		else 
			sprintf(file_suffix,"hind.%d.%s",timekeeper.current_year,run_parameters.timestep);
		printf("Hindcast Mass transports and isopycnal thickness: Year %d, Interval: %d\n",
				timekeeper.current_year,timekeeper.current_interval);
		read_index = timekeeper.current_interval;
		strcpy(path,run_parameters.hindcast_path);
	}
	else {
		sprintf(file_suffix,"clim.%s",run_parameters.timestep);
		printf("Climatology Mass transports and isopycnal thickness: Climatology index: %d\n",
				timekeeper.climatology_index);
		read_index = timekeeper.climatology_index;
		strcpy(path,run_parameters.normalyear_path);
	}

	read_uvw(read_index,file_suffix,path);
	read_h(read_index,file_suffix,path,hend);;



}

