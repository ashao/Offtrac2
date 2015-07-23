/*
 * sf6.c
 *
 *  Created on: May 2, 2014
 *      Author: ashao
 */
#include <stdio.h>
#include "init.h"
#include "alloc.h"
#include "read.h"
#include "math.h"
#include "cfcs_sf6.h"
#include "tracer_utilities.h"
#include "io.h"
#include "string.h"
#include "initialize.h"
#include "timekeeper.h"
#include "util.h"
#include "output_variables.h"

double ***mn_sf6;
double ***sf6_init;

double **sf6_sat;
double ***sf6_sol;
double **mn_sf6sat;
double **sf6_atmconc;
double ***psf6, ***mn_psf6;

int mSF6;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];

void allocate_sf6 ( ) {

	mn_sf6 = alloc3d(NZ,NXMEM,NYMEM);
	sf6_init = alloc3d(NZ,NXMEM,NYMEM);
	sf6_sat = alloc2d(NXMEM,NYMEM);
	sf6_sol = alloc3d(NZ,NXMEM,NYMEM);
	mn_sf6sat = alloc2d(NXMEM,NYMEM);
	sf6_atmconc = alloc2d(NXMEM,NYMEM);


	mn_psf6 = alloc3d(NZ,NXMEM,NYMEM);
	psf6 = alloc3d(NZ,NXMEM,NYMEM);
}

void initialize_sf6 ( ) {
	int i, j, k;
	char varname[200];
	extern struct parameters run_parameters;	
	extern struct vardesc vars[NOVARS];	

	printf("SF6 index in main tracer array: %d\n",mSF6);
	printf("Setting SF6 variable description...");
	strcpy(varname,"sf6");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_sf6");
	strcpy(vars[map_variable_to_index(varname)].longname,"SF6 volumetric concentration");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"pmol/L");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;

	strcpy(varname,"psf6");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_psf6");
	strcpy(vars[map_variable_to_index(varname)].longname,"Partial pressure of SF6");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"parts-per-trillion");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;
	printf("DONE\n");


	if (run_parameters.restart_flag) {
		printf("Initializing SF6 from restart %s\n",run_parameters.restartfile);
		read_var3d( run_parameters.restartfile, "mn_sf6", 0, sf6_init);
	}
	else	set_darray3d_zero(sf6_init,NZ,NXMEM,NYMEM);

	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			for(k=0;k<NZ;k++)
				tr[mSF6][k][i][j] = sf6_init[k][i][j];

	free3d(sf6_init,NZ);
}

void sf6_saturation(  ) {

	const double solcoeffs[6] = {-80.0343,117.232,29.5817,0.0335183,-0.0373942,0.00774862};
	int i, j, k;
	double work;
	double TempK;
	extern double ***Temptm, ***Salttm;

	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
			{
				if(oceanmask[i][j]) {
					TempK = Temptm[k][i][j]+273.15;
					work = solcoeffs[0] + solcoeffs[1]*(100/TempK) +
							solcoeffs[2]*log( TempK/100 ) +
							Salttm[k][i][j]*(solcoeffs[3]+solcoeffs[4]*(TempK/100)+solcoeffs[5]*pow(TempK/100,2));
					sf6_sol[k][i][j] = exp(work);
					if(k<NML) sf6_sat[i][j] = sf6_sol[k][i][j]*sf6_atmconc[i][j];
				}
				else {
					sf6_sol[k][i][j] = 0.0;
					sf6_sat[i][j] = 0.0;
				}
			}



}

void sf6_find_atmconc(  ) {

	int i,j;
	extern double hlat[NYMEM];
	const double equatorbound[2] = {10,-10}; // Set the latitudes where to start interpolating atmospheric concentrations
	double hemisphere_concentrations[2];
	double currtime;
	extern struct timekeeper_t timekeeper;

	currtime = timekeeper.current_time;

	// Interpolate in time to find the atmospheric concentration
	hemisphere_concentrations[0] = linear_interpolation(
			atmconc[SF6IDX].time, atmconc[SF6IDX].nval, currtime,NUMATMVALS);
	hemisphere_concentrations[1] = linear_interpolation(
			atmconc[SF6IDX].time, atmconc[SF6IDX].sval, currtime,NUMATMVALS);


	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++) {

			if (hlat[j] < equatorbound[0] && hlat[j] > equatorbound[1]) {
				sf6_atmconc[i][j] = linear_interpolation(
						equatorbound,hemisphere_concentrations,hlat[j],2);
			}
			if (hlat[j]>equatorbound[0] ) {
				sf6_atmconc[i][j] = hemisphere_concentrations[0];
			}
			if (hlat[j]<equatorbound[1] ) {
				sf6_atmconc[i][j] = hemisphere_concentrations[1];
			}
		}

}

void surface_sf6( ) {

	int i,j,k;
	extern double ***Salttm, ***Temptm;
	extern struct timekeeper_t timekeeper;
	printf("Setting SF6 surface condition\n");
	// Set oxygen values to saturation at the mixed layer to mimic equilibrium with the atmosphere
	sf6_find_atmconc( );
	sf6_saturation( );

	printf("SF6 Sample Point\n");
	printf("\tTime: %f\n",timekeeper.current_time);
	printf("\tAtmospheric Concentration: %f\n",sf6_atmconc[100][100]);
	printf("\tSalinity: %f\t Temperature: %f\n",Salttm[0][100][100],Temptm[0][100][100]);
	printf("\tSaturation: %f\n\n",sf6_sat[100][100]);

#ifdef NOCONC
	for (k=0;k<NML;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				tr[mSF6][k][i][j] = sf6_atmconc[i][j];

#else

	for (k=0;k<NML;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				tr[mSF6][k][i][j]=sf6_sat[i][j];
#endif
}

