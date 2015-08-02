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
#include <time.h>
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
#include "read.h"
#include "gas_exchange.h"

#include "ideal_age.h"
#include "cfcs_sf6.h"
#include "n2_module.h"
#include "ar_module.h"
#include "biotic.h"
#include "oxygen.h"
#include "phosphate.h"

#ifdef CONSERVATION_CHECK
#include "conservation_check.h"
extern int mTEST;
extern  double ***mn_test;
extern double test_inventory;
#endif

#include "ttd_bp.h"
/*---------------------------------------------------------------------
 *     define variables and subroutines
 *---------------------------------------------------------------------*/

extern double ****tr;
extern double ***hstart, ***h, ***hend;

extern double ***mn_h, ***h;
extern double ***mn_uhtm, ***uhtm;
extern double ***mn_vhtm, ***vhtm;
extern double ***mn_wd, ***wd;

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
	int i, j, k,  m;
	int ii;

	struct timespec startclock, endclock;
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
	clock_gettime(CLOCK_MONOTONIC, &startclock); 
	printf("Calculate tracer transport. \n");
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


	if (run_parameters.do_age)	step_age(timekeeper.dt);

	if (run_parameters.do_cfcs) {
		surface_cfc11();
		surface_cfc12();
		surface_sf6();
# ifdef NOCONC

		copy_darray3d(pcfc11,tr[mCFC11],NZ,NXMEM,NYMEM);
		copy_darray3d(pcfc12,tr[mCFC12],NZ,NXMEM,NYMEM);
		copy_darray3d(psf6,tr[mSF6],NZ,NXMEM,NYMEM);
# else
		divide_darray3d(pcfc11,tr[mCFC11],cfc11_sol);
		divide_darray3d(pcfc12,tr[mCFC12],cfc12_sol);
		divide_darray3d(psf6,tr[mSF6],sf6_sol);
# endif
	}

#ifdef CONSERVATION_CHECK
	step_test();
	submit_for_averaging(mn_test,tr[mTEST]);
#endif


	if (run_parameters.do_ttd)	step_ttd();
	if (run_parameters.do_n2)	step_n2();
	if (run_parameters.do_ar)	step_ar();

	if (run_parameters.do_oxygen) {

		update_phosphate_fields();
		biotic_sms( 1, timekeeper.dt );
		step_oxygen();
		step_phosphate();

	}

	merge_ml_tr();

        clock_gettime(CLOCK_MONOTONIC, &endclock);
	printf("Current step elapsed time: %fs\n", (double) (endclock.tv_sec-startclock.tv_sec) + (double) (endclock.tv_nsec-startclock.tv_nsec)/1.e9) ;


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
	for (m = 0; m < run_parameters.tracer_counter; m++) {
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


	for (m=0;m<run_parameters.tracer_counter;m++)
		for(i=0;i<NXMEM;i++)
			for(j=0;j<NYMEM;j++)
				if (!oceanmask[i][j])
					for(k=0;k<NZ;k++)
						tr[m][k][i][j]=MISVAL;


	submit_for_averaging( mn_h, h );
	submit_for_averaging( mn_uhtm, uhtm );
	submit_for_averaging( mn_vhtm, vhtm );
	submit_for_averaging( mn_wd, wd );

	if (run_parameters.do_age)	submit_for_averaging( mn_age, tr[mAGE]) ;

	if (run_parameters.do_cfcs) {
		submit_for_averaging( mn_cfc11, tr[mCFC11] );
		submit_for_averaging( mn_pcfc11, pcfc11 );
		submit_for_averaging( mn_cfc12, tr[mCFC12] );
		submit_for_averaging( mn_pcfc12, pcfc12 );
		submit_for_averaging( mn_sf6, tr[mSF6] );
		submit_for_averaging( mn_psf6, psf6 );
	}

	if (run_parameters.do_ttd)	submit_for_averaging( mn_ttd, tr[mTTD] );
	if (run_parameters.do_n2)	{

		submit_for_averaging( mn_n2, tr[mN2] );
		submit_for_averaging( mn_n2sat, n2sat );
	}
	if (run_parameters.do_ar)	{

		submit_for_averaging( mn_ar, tr[mAR] );
		submit_for_averaging( mn_arsat, arsat );
	}

	if (run_parameters.do_oxygen) {

		submit_for_averaging(mn_oxygen,tr[mOXYGEN]);
		submit_for_averaging(mn_jo2,jo2);
		submit_for_averaging(mn_o2sat,o2_sat);
		submit_for_averaging(mn_phos,tr[mPHOSPHATE]);
		submit_for_averaging(mn_dop,tr[mDOP]);
		submit_for_averaging(mn_jpo4,jpo4);

	}


	printf("\n");



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

	for(l=0;l<run_parameters.tracer_counter;l++)
		for(i=0;i<NXMEM;i++)
			for(j=0;j<NYMEM;j++) {

				tracsum = 0; // Reset the tracer amount
				for(k=0;k<NML;k++) {
					tracsum+=tr[l][k][i][j];
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
	read_h(read_index,file_suffix,path,hend);
	z_depth(hend,depth);
	read_temp_and_salt(read_index,file_suffix,path);
	update_gas_exchange_fields();

}

