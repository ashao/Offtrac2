/*
 * n2_module.c
 *
 *  Created on: Apr 17, 2013
 *      Author: ashao
 *      Contains all subroutines related to the nitrogen tracer
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
#include "n2_module.h"
#include "alloc.h"
#include "init.h"
#include "gsw_src/gswteos-10.h"

// Auxiliary variables
int mN2;
// Output arrays
double ***mn_n2;
double ***mn_n2sol;
// Working arrays
double ***n2_init;
double ***n2sol;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];
extern struct vardesc vars[NOVARS];
extern struct parameters run_parameters;


struct Gas_const_t n2_props;

void allocate_n2 (  ) {
	printf("Allocating n2 arrays\n");

	// Allocate working and output arrays
	mn_n2 = alloc3d(NZ,NXMEM,NYMEM);
	mn_n2sol = alloc3d(NZ,NXMEM,NYMEM);
	n2sol = alloc3d(NZ,NXMEM,NYMEM);
	n2_init = alloc3d(NZ,NXMEM,NYMEM);

}

void initialize_n2( ) {
	int i,j,k;
	//	extern char restart_filename[200];
	//	extern double misval;
	char varname[100];


	printf("Setting N2 Gas properties\n");
	initialize_n2_properties();
	printf("N2 index in main tracer array: %d\n",mN2);
	printf("Setting n2 variable description...");
	strcpy(varname,"n2");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_n2");
	strcpy(vars[map_variable_to_index(varname)].longname,"N2 Concentration");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"mol/m^3");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;

	strcpy(varname,"n2sol");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_n2sol");
	strcpy(vars[map_variable_to_index(varname)].longname,"N2 Saturation");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"mol/m^3");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;

	printf("DONE\n");

	if (run_parameters.restart_flag) {
		printf("Initializing n2 from restart %s\n",run_parameters.restartfile);
		read_var3d( run_parameters.restartfile, "mn_n2", 0, n2_init );

	}
	else {
		// Initialize N2 to saturation everywhere in the ocean
		calc_n2_saturation();
		copy_darray3d( n2_init, n2sol, NZ, NXMEM, NYMEM );
		printf("Ideal n2 initialized to saturation globally\n");
	}

	wrap_reentrance_3d(n2_init,NZ);


	// Copy the initialized tracer value over to main trace array
	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
			for (k=0;k<NZ;k++) {
				if (oceanmask[i][j]) tr[mN2][k][i][j] = n2_init[k][i][j];
				else tr[mN2][k][i][j] = 0.0;
			}

	copy_darray3d(mn_n2,tr[mN2],NZ,NXMEM,NYMEM);
	printf("N2 init example: %f\n",tr[mN2][15][100][100]);
	free3d(n2_init,NZ);


}


void initialize_n2_properties ( )
{

	const int num_Sc_coeffs = 4;
	const int num_sat_coeffs = 7;    ;

	// Store the number of coefficients in the tracer type
	n2_props.num_Sc_coeffs = num_Sc_coeffs;
	n2_props.num_sat_coeffs = num_sat_coeffs;

	// Allocate memory for coefficients
	n2_props.Sc_coeffs = (double *)malloc(num_Sc_coeffs * sizeof(double));
	n2_props.sat_coeffs = (double *)malloc(num_sat_coeffs * sizeof(double));
	n2_props.atmconc = 0.780840;


	// Wanninkhof 1992
	n2_props.Sc_coeffs[0] = 2206.1;
	n2_props.Sc_coeffs[1] = 144.86;
	n2_props.Sc_coeffs[2] = 4.5413;
	n2_props.Sc_coeffs[3] = 0.056988;
	// Hamme and Emerson 2007
	n2_props.sat_coeffs[0] = 6.42931;
	n2_props.sat_coeffs[1] = 2.92704;
	n2_props.sat_coeffs[2] = 4.32531;
	n2_props.sat_coeffs[3] = 4.69149;
	n2_props.sat_coeffs[4] = -7.44129e-3;
	n2_props.sat_coeffs[5] = -8.02566e-3;
	n2_props.sat_coeffs[6] = -1.46775e-2;

}

void calc_n2_saturation( ) {
	// Saturation concentration calculated from HAmme and Emmerson 2004
	// Note that this outputs in volumetric by dividing by potential density
	// (referenced to 2000db as in the model)

	int i,j,k;

	double pden;
	double temp_S, S; // scaled temperature, salinity
	double conc_N2;
	extern double ***Temptm;
	extern double ***Salttm;


	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)	{

				if(oceanmask[i][j])
				{
					temp_S = log((298.15-Temptm[k][i][j])/(273.15+Temptm[k][i][j]));
					S = Salttm[k][i][j];
					conc_N2 = n2_props.sat_coeffs[0] +
							n2_props.sat_coeffs[1]*temp_S +
							n2_props.sat_coeffs[2]*pow(temp_S,2.) +
							n2_props.sat_coeffs[3]*pow(temp_S,3. ) +
							S*(n2_props.sat_coeffs[4] +
									n2_props.sat_coeffs[5]*temp_S +
									n2_props.sat_coeffs[6]*pow(temp_S,2.));

					conc_N2 = exp(conc_N2);
					pden = gsw_rho_t_exact(S,Temptm[k][i][j],2000.0);
					n2sol[k][i][j] = conc_N2 / pden;

				}
				else {
					n2sol[k][i][j] = 0.0;
				}
			}
}

void step_n2( ) {

	int i,j,k;

	calc_n2_saturation( );

	for (k=0;k<NML;k++)
		for (i=0; i<NXMEM; i++)
			for (j=0;j<NYMEM; j++) {
				// Set surface to saturation
				tr[mN2][k][i][j] = n2sol[k][i][j];

			}



}



