/*
 * cfc12.c
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
#include "string.h"
#include "io.h"
#include "initialize.h"
#include "timekeeper.h"

double ***mn_cfc12;
double ***cfc12_init;

double **cfc12_sat;
double ***cfc12_sol;
double **mn_cfc12sat;
double **cfc12_atmconc;
double ***pcfc12, ***mn_pcfc12;

int mCFC12;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];
void allocate_cfc12 ( ) {

	mn_cfc12 = alloc3d(NZ,NXMEM,NYMEM);
	cfc12_init = alloc3d(NZ,NXMEM,NYMEM);
	cfc12_sat = alloc2d(NXMEM,NYMEM);
	cfc12_sol = alloc3d(NZ,NXMEM,NYMEM);
	mn_cfc12sat = alloc2d(NXMEM,NYMEM);
	cfc12_atmconc = alloc2d(NXMEM,NYMEM);

	mn_pcfc12 = alloc3d(NZ,NXMEM,NYMEM);
	pcfc12 = alloc3d(NZ,NXMEM,NYMEM);
}

void initialize_cfc12 ( ) {
	int i, j, k;
	extern struct parameters run_parameters;
	char varname[200];	
	extern struct vardesc vars[NOVARS];	
	mCFC12 = run_parameters.tracer_counter++;
	printf("mCFC12: %d\n",mCFC12);

        printf("Setting CFC-12 variable description...");
        strcpy(varname,"cfc12");
        strcpy(vars[map_variable_to_index(varname)].name,"mn_cfc12");
        strcpy(vars[map_variable_to_index(varname)].longname,"CFC-12 volumetric concentration");
        vars[map_variable_to_index(varname)].hor_grid='h';
        vars[map_variable_to_index(varname)].z_grid='L';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"pmol/L");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;

        strcpy(varname,"pcfc12");
        strcpy(vars[map_variable_to_index(varname)].name,"mn_pcfc12");
        strcpy(vars[map_variable_to_index(varname)].longname,"Partial pressure of CFC-12");
        vars[map_variable_to_index(varname)].hor_grid='h';
        vars[map_variable_to_index(varname)].z_grid='L';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"parts-per-trillion");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;
        printf("DONE\n");



        if (run_parameters.restart_flag) {
                printf("Initializing CFC-1@ from restart %s\n",run_parameters.restartfile);
                read_var3d( run_parameters.restartfile, "mn_cfc12", 0, cfc12_init);
        }
	else	set_darray3d_zero(cfc12_init,NZ,NXMEM,NYMEM);

	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			for(k=0;k<NZ;k++)
				tr[mCFC12][k][i][j] = cfc12_init[k][i][j];

	free3d(cfc12_init,NZ);
}

void cfc12_saturation(  ) {

	const double solcoeffs[7] = {-218.0971,298.9702,113.8049,-1.39165,-0.143566,0.091015,-0.0153924};
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
					work = solcoeffs[0] + solcoeffs[1]*(100/TempK) +solcoeffs[2]*log( TempK/100 ) +
						solcoeffs[3]*pow(TempK/100,2) + Salttm[k][i][j]*(solcoeffs[4]+
						solcoeffs[5]*(TempK/100)+solcoeffs[6]*pow(TempK/100,2));
						cfc12_sol[k][i][j] = exp(work);
	
						if (k<NML) cfc12_sat[i][j] = cfc12_sol[k][i][j]*cfc12_atmconc[i][j];
				}		
				else {
					cfc12_sol[k][i][j] = 0.0;
					cfc12_sat[i][j] = 0.0;
				}
				
			}




}

void cfc12_find_atmconc(  ) {

	int i,j;
	extern double hlat[NYMEM];
	extern double **atmpres;
	const double equatorbound[2] = {10,-10}; // Set the latitudes where to start interpolating atmospheric concentrations
	double currtime;
	double hemisphere_concentrations[2];
	extern struct timekeeper_t timekeeper;
	
	currtime = timekeeper.current_time;

	// Interpolate in time to find the atmospheric concentration
	hemisphere_concentrations[0] = linear_interpolation(
		atmconc[mCFC12].time, atmconc[mCFC12].nval, currtime,NUMATMVALS);
	hemisphere_concentrations[1] = linear_interpolation(
		atmconc[mCFC12].time, atmconc[mCFC12].sval, currtime,NUMATMVALS);
	
for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++) {

			if (hlat[j] < equatorbound[0] && hlat[j] > equatorbound[1]) {
				cfc12_atmconc[i][j] = linear_interpolation(
						equatorbound,hemisphere_concentrations,hlat[j],2);
			}
			if (hlat[j]>equatorbound[0] ) {
				cfc12_atmconc[i][j] = hemisphere_concentrations[0];
			}
			if (hlat[j]<equatorbound[1] ) {
				cfc12_atmconc[i][j] = hemisphere_concentrations[1];
			}
		}

}

void surface_cfc12( ) {

	int i,j,k;
	extern double ***Salttm, ***Temptm;
	extern struct timekeeper_t timekeeper;
	printf("Setting CFC-12 surface condition\n");
	// Set oxygen values to saturation at the mixed layer to mimic equilibrium with the atmosphere
	cfc12_find_atmconc( );
	cfc12_saturation(  );

	printf("CFC12 atmospheric values %d\n",NUMATMVALS);
	i=100;
	printf("%f %f\n",atmconc[mCFC12].time[i],atmconc[mCFC12].nval[i]);

        printf("CFC-12 Sample Point\n");
        printf("\tTime: %f\n",timekeeper.current_time);
        printf("\tAtmospheric Concentration: %f\n",cfc12_atmconc[100][100]);
        printf("\tSalinity: %f\t Temperature: %f\n",Salttm[0][100][100],Temptm[0][100][100]);
        printf("\tSolubility: %f\tSaturation: %f\n\n",cfc12_sol[0][100][100],cfc12_sat[100][100]);




#ifdef NOCONC
        for (k=0;k<NML;k++)
                for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
                        tr[mCFC12][k][i][j] = cfc12_atmconc[i][j];

#else
        for (k=0;k<NML;k++)
                for (i=0;i<NXMEM;i++)
                        for (j=0;j<NYMEM;j++)
                                tr[mCFC12][k][i][j]=cfc12_sat[i][j];

#endif

}

