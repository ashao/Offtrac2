/*
 * phosphate.c
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */
#include "init.h"
#include "netcdf.h"
#include "phosphate.h"
#include "alloc.h"
#include <stdio.h>

#include "timekeeper.h"
/* PHOSPHATE VARIABLE DECLARATIONS */
// Auxiliary variables
int mPHOSPHATE;
int mDOP;
// Output arrays
double ***mn_phos;
double ***mn_dop;
double ***mn_pobs;
double ***mn_jpo4;
double ***mn_jremin;
double ***mn_jremdop;
double mn_pflux[NXMEM][NYMEM];
double ***mn_jprod;

// Working arrays
double ***phosphate_init;
// double po4[NZ][NXMEM][NYMEM];
double ***jpo4;
double ***po4_star_lay;
double ***jprod;
double ***jremin;
double ***jremdop;
double ***jdop;
double flux_pop[NXMEM][NYMEM];
double ***dop_init;

extern struct timekeeper_t timekeeper;
void allocate_phosphate(  )
{
	int i,j,k;
	// Base the indices of the main tracer array off of the oxygen index
	printf("mDOP: %d mPHOSPHATE: %d\n",mDOP,mPHOSPHATE);
	mn_phos = alloc3d(NZ,NXMEM,NYMEM);
	mn_dop = alloc3d(NZ,NXMEM,NYMEM);
	mn_pobs = alloc3d(NZ,NXMEM,NYMEM);
	mn_jpo4 = alloc3d(NZ,NXMEM,NYMEM);
	mn_jremin = alloc3d(NZ,NXMEM,NYMEM);
	mn_jprod = alloc3d(NZ,NXMEM,NYMEM);

	jpo4 = alloc3d(NZ,NXMEM,NYMEM);
	po4_star_lay = alloc3d(NZ,NXMEM,NYMEM);
	jprod = alloc3d(NZ,NXMEM,NYMEM);
	jremin = alloc3d(NZ,NXMEM,NYMEM);
	jremdop = alloc3d(NZ,NXMEM,NYMEM);
	jdop = alloc3d(NZ,NXMEM,NYMEM);

	phosphate_init = alloc3d(NZ,NXMEM,NYMEM);
	dop_init = alloc3d(NZ,NXMEM,NYMEM);
	po4_star_lay = alloc3d(NZ,NXMEM,NYMEM);
}
void initialize_phosphate( int imon ) {
	int i,j,k;
	extern double ***hend;
	extern double ****tr;
	extern char restart_filename[200];
	char varname[100];

	printf("Phosphate index in main tracer array: %d\n",mPHOSPHATE);
	printf("Setting phosphate variable description...");
	strcpy(varname,"phosphate");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_po4");
	strcpy(vars[map_variable_to_index(varname)].longname,"Phosphate concentration");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"mol/m^3");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;

	strcpy(varname,"dop");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_o2sat");
	strcpy(vars[map_variable_to_index(varname)].longname,"Dissolved organic phosphorous concentration");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"mol/m^3");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;

	strcpy(varname,"jpo4");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_jpo4");
	strcpy(vars[map_variable_to_index(varname)].longname,"Phosphate source/sink term");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"mol/m^3/s");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;


	printf("Initializing phosphate from WOA09\n");
	read_woa_file(imon, hend, phosphate_init, "woa09.phos.nc", "p_an",1e-3);
	printf("Initializing DOP to zero\n");
	set_darray3d_zero(dop_init, NZ, NXMEM, NYMEM);


	if (run_parameters.restart_flag) {
		printf("Initializing phosphate from restart: %s\n",restart_filename);
		read_var3d( restart_filename, "mn_phos", 0, phosphate_init);
		printf("Initializing dop from restart: %s\n",restart_filename);
		read_var3d( restart_filename, "mn_dop", 0, dop_init);
	}


	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for(j=0;j<NYMEM;j++) {
			if (dop_init[k][i][j] > 0.0) tr[mDOP][k][i][j] = dop_init[k][i][j];
			else tr[mDOP][k][i][j] = 0.0;
			if (phosphate_init[k][i][j]>0.0) tr[mPHOSPHATE][k][i][j] = phosphate_init[k][i][j];
			else tr[mPHOSPHATE][k][i][j] = 0.0;

			}
	free3d(dop_init,NZ);
	free3d(phosphate_init,NZ);
	printf("Phosphate Initialization:%e\n",tr[mPHOSPHATE][20][100][100]);
}

void apply_phosphate_jterms( ) {
	int i,j,k;

	extern double ****tr;
	extern int oceanmask[NXMEM][NYMEM];
	extern double ***hend;
	extern double misval;
	// j terms here are calculated from biotic_sms routine in biotic.c
	printf("Applying j terms for phosphate\n");
	printf("Example jpo4/po4: %e/%f\n",jpo4[10][127][127],tr[mPHOSPHATE][10][127][127]);
	for (i = 0; i <= NXMEM - 1; i++) {
		for (j = 0; j <= NYMEM - 1; j++) {
			//BX - reinstated by HF
			if (oceanmask[i][j]) {
				for (k = 0; k < NZ; k++) {
					if (hend[k][i][j]>EPSILON) tr[mPHOSPHATE][k][i][j] += dt * jpo4[k][i][j];
					if (hend[k][i][j]>EPSILON) tr[mDOP][k][i][j] += dt * jdop[k][i][j];
					}
				
			} else {

				for (k = 0; k < NZ; k++) {
				tr[mPHOSPHATE][k][i][j] = misval;
				tr[mDOP][k][i][j] = misval;
				}
			}
		}
	}


}
