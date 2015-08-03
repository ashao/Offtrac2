/*
 * ar_module.c
 *
 *  Created on: Apr 17, 2013
 *      Author: ashao
 *      Contains all subroutines related to the argon tracer
 */
#include <stdlib.h>
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
#include "gas_tracer_type.h"
#include "ar_module.h"
#include "alloc.h"
#include "init.h"
#include "gsw_src/gswteos-10.h"

// Auxiliary variables
int mAR;
// Output arrays
double ***mn_ar;
double ***mn_arsat;
// Working arrays
double ***ar_init;
double ***arsat;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];
extern struct vardesc vars[NOVARS];
extern struct parameters run_parameters;


struct Gas_const_t ar_props;

void allocate_ar (  ) {
	printf("Allocating ar arrays\n");

	// Allocate working and output arrays
	mn_ar = alloc3d(NZ,NXMEM,NYMEM);
	mn_arsat = alloc3d(NZ,NXMEM,NYMEM);
	arsat = alloc3d(NZ,NXMEM,NYMEM);
	ar_init = alloc3d(NZ,NXMEM,NYMEM);

}

void initialize_ar( ) {
	int i,j,k;
	//	extern char restart_filename[200];
	//	extern double misval;
	char varname[100];


	printf("Setting AR Gas properties\n");
	initialize_ar_properties();
	printf("AR index in main tracer array: %d\n",mAR);
	printf("Setting ar variable description...");
	strcpy(varname,"ar");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_ar");
	strcpy(vars[map_variable_to_index(varname)].longname,"Argon Concentration");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"micromol/m^3");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;

	strcpy(varname,"arsol");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_arsol");
	strcpy(vars[map_variable_to_index(varname)].longname,"Argon Saturation");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"micromol/m^3");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;

	printf("DONE\n");

	if (run_parameters.restart_flag) {
		printf("Initializing ar from restart %s\n",run_parameters.restartfile);
		read_var3d( run_parameters.restartfile, "mn_ar", 0, ar_init );

	}
	else {
		// Initialize AR to saturation everywhere in the ocean
		calc_ar_saturation();
		copy_darray3d( ar_init, arsat, NZ, NXMEM, NYMEM );
		printf("Ideal ar initialized to saturation globally\n");
	}

	wrap_reentrance_3d(ar_init,NZ);


	// Copy the initialized tracer value over to main trace array
	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
			for (k=0;k<NZ;k++) {
				if (oceanmask[i][j]) tr[mAR][k][i][j] = ar_init[k][i][j];
				else tr[mAR][k][i][j] = 0.0;
			}

	copy_darray3d(mn_ar,tr[mAR],NZ,NXMEM,NYMEM);
	printf("AR init example: %f\n",tr[mAR][15][100][100]);
	free3d(ar_init,NZ);


}


void initialize_ar_properties ( )
{

	const int num_Sc_coeffs = 4;
	const int num_sat_coeffs = 7;    ;

	// Store the number of coefficients in the tracer type
	ar_props.num_Sc_coeffs = num_Sc_coeffs;
	ar_props.num_sat_coeffs = num_sat_coeffs;

	// Allocate memory for coefficients
	ar_props.Sc_coeffs = (double *)malloc(num_Sc_coeffs * sizeof(double));
	ar_props.sat_coeffs = (double *)malloc(num_sat_coeffs * sizeof(double));
	ar_props.atmconc = 0.780840;


	// Wanninkhof 1992
	ar_props.Sc_coeffs[0] = 1909.1;
	ar_props.Sc_coeffs[1] = 125.09;
	ar_props.Sc_coeffs[2] = 3.9012;
	ar_props.Sc_coeffs[3] = 0.048953;

	// Hamme and Emerson 2004
	ar_props.sat_coeffs[0] = 2.79150;
	ar_props.sat_coeffs[1] = 3.17609;
	ar_props.sat_coeffs[2] = 4.13116;
	ar_props.sat_coeffs[3] = 4.90379;
	ar_props.sat_coeffs[4] = -6.96233e-3;
	ar_props.sat_coeffs[5] = -7.66670e-3;
	ar_props.sat_coeffs[6] = -1.16888e-2;

}

void calc_ar_saturation( ) {
	// Saturation concentration calculated from HAmme and Emmerson 2004
	// Note that this outputs in volumetric by dividing by potential density
	// (referenced to 2000db as in the model)

	int i,j,k;

	double pden;
	double temp_S, S; // scaled temperature, salinity
	double conc_AR;
	extern double ***Temptm;
	extern double ***Salttm;


	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)	{

				if(oceanmask[i][j])
				{
					temp_S = log((298.15-Temptm[k][i][j])/(273.15+Temptm[k][i][j]));
					S = Salttm[k][i][j];
					conc_AR = ar_props.sat_coeffs[0] +
							ar_props.sat_coeffs[1]*temp_S +
							ar_props.sat_coeffs[2]*pow(temp_S,2.) +
							ar_props.sat_coeffs[3]*pow(temp_S,3. ) +
							S*(ar_props.sat_coeffs[4] +
									ar_props.sat_coeffs[5]*temp_S +
									ar_props.sat_coeffs[6]*pow(temp_S,2.));

					conc_AR = exp(conc_AR);
					pden = gsw_rho_t_exact(S,Temptm[k][i][j],2000.0);
					arsat[k][i][j] = conc_AR / pden;

				}
				else {
					arsat[k][i][j] = 0.0;
				}
			}
}

void step_ar( ) {
	// Surface source of Argon
	int k;


	calc_ar_saturation( );

	if (run_parameters.do_gasex)
		gas_exchange(mAR,ar_props.Sc_coeffs,arsat);
	else
		for (k=0;k<NML;k++)
			copy_darray2d(tr[mAR][k],arsat[k],NXMEM,NYMEM);

}



