/*
 * cfc11.c
 *
 *  Created on: May 2, 2014
 *      Author: ashao
 */
#include "init.h"
#include "alloc.h"
#include "read.h"
#include <math.h>
#include "cfcs_sf6.h"
#include "tracer_utilities.h"
#include "timekeeper.h"
#include "initialize.h"
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include "io.h"

double ***mn_cfc11;
double ***cfc11_init;

double **cfc11_sat;
double ***cfc11_sol;
double **cfc11_atmconc;
double **mn_cfc11sat;
double ***pcfc11, ***mn_pcfc11;

int mCFC11;

double nval;
double sval;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];
tracer_boundary atmconc[NUMTRANSIENT];

void allocate_cfc11 ( ) {

	mn_cfc11 = alloc3d(NZ,NXMEM,NYMEM);

	cfc11_init = alloc3d(NZ,NXMEM,NYMEM);

	cfc11_sat = alloc2d(NXMEM,NYMEM);
	mn_cfc11sat = alloc2d(NXMEM,NYMEM);
	cfc11_sol = alloc3d(NZ,NXMEM,NYMEM);

	pcfc11 = alloc3d(NZ,NXMEM,NYMEM);
	mn_pcfc11 = alloc3d(NZ,NXMEM,NYMEM);

	cfc11_atmconc = alloc2d(NXMEM,NYMEM);

}

void read_tracer_boundary ( ) {

	int err, i, j, cdfid, timeid;
	int status, varid;
	char infile[25], inpath[200];
	FILE *file;
	long  start[MAX_NC_VARS];
	long  end[MAX_NC_VARS];
	char varname[20];
	const int numatm = NUMATMVALS;
	extern struct parameters run_parameters;
 
	int tempidx;
	sprintf(infile,"cfc_sf6_bc.nc");
	strcpy(inpath, run_parameters.forcing_path);
	strcat(inpath, infile);
	printf("\nLooking for file '%s'.\n",inpath);
	err = open_input_file(inpath,&file,&cdfid,&timeid);

	start[0] = 0;
	end[0] = NUMATMVALS;

	status = nc_inq_varid(cdfid, "Year", &varid);
	printf("mCFC11: %d, mCFC12: %d, mSF6: %d\n",mCFC11,mCFC12,mSF6);

	for (i=0;i<NUMTRANSIENT;i++) // Put the year into all the	
	    	status = nc_get_vara_double(cdfid, varid, start, end,atmconc[i].time);

	strcpy(varname,"CFC11NH");
	status = nc_inq_varid(cdfid, varname, &varid);
	status = nc_get_vara_double(cdfid, varid, start, end, atmconc[CFC11IDX].nval);

	strcpy(varname,"CFC11SH");
	status = nc_inq_varid(cdfid, varname, &varid);
	status = nc_get_vara_double(cdfid, varid, start, end, atmconc[CFC11IDX].sval);

	strcpy(varname,"CFC12NH");
	status = nc_inq_varid(cdfid, varname, &varid);
	status = nc_get_vara_double(cdfid, varid, start, end, atmconc[CFC12IDX].nval);
	strcpy(varname,"CFC12SH");
	status = nc_inq_varid(cdfid, varname, &varid);
	status = nc_get_vara_double(cdfid, varid, start, end, atmconc[CFC12IDX].sval);

	strcpy(varname,"SF6NH");
	status = nc_inq_varid(cdfid, varname, &varid);
	status = nc_get_vara_double(cdfid, varid, start, end, atmconc[SF6IDX].nval);
	strcpy(varname,"SF6SH");
	status = nc_inq_varid(cdfid, varname, &varid);
	status = nc_get_vara_double(cdfid, varid, start, end, atmconc[SF6IDX].sval);
	
	printf("status=%d\n",status);
	close_file(&cdfid,&file);
}

void initialize_cfc11 ( ) {
	int i, j, k;
	char varname[200];
	extern char restart_filename[200];
	extern struct parameters run_parameters;
	extern struct vardesc vars[NOVARS];


	mCFC11 = run_parameters.tracer_counter++;
        printf("Setting CFC-11 variable description...");
        strcpy(varname,"cfc11");
        strcpy(vars[map_variable_to_index(varname)].name,"mn_cfc11");
        strcpy(vars[map_variable_to_index(varname)].longname,"CFC-11 volumetric concentration");
        vars[map_variable_to_index(varname)].hor_grid='h';
        vars[map_variable_to_index(varname)].z_grid='L';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"pmol/L");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;

        strcpy(varname,"pcfc11");
        strcpy(vars[map_variable_to_index(varname)].name,"mn_pcfc11");
        strcpy(vars[map_variable_to_index(varname)].longname,"Partial pressure of CFC-11");
        vars[map_variable_to_index(varname)].hor_grid='h';
        vars[map_variable_to_index(varname)].z_grid='L';
        vars[map_variable_to_index(varname)].t_grid='s';
        strcpy(vars[map_variable_to_index(varname)].units,"parts-per-trillion");
        vars[map_variable_to_index(varname)].mem_size='d';
        vars[map_variable_to_index(varname)].mval=MISVAL;
        printf("DONE\n");


	if (run_parameters.restart_flag) {
		printf("Initializing CFC-11 from restart %s\n",run_parameters.restartfile);
		read_var3d( run_parameters.restartfile, "mn_cfc11", 0, cfc11_init);
	}
	else	set_darray3d_zero(cfc11_init,NZ,NXMEM,NYMEM);

	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++)
			for(k=0;k<NZ;k++)
				tr[mCFC11][k][i][j] = cfc11_init[k][i][j];

	free3d(cfc11_init,NZ);
}

void cfc11_saturation(  ) {

	const double solcoeffs[7] = {-229.9261,319.6552,119.4471,-1.39165,-0.142382,0.091459,-0.0157274};
	int i, j, k;
	double work;
	double TempK;
	extern double ***Temptm, ***Salttm;
	for (k=0;k<NZ;k++)
	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)	{

			if(oceanmask[i][j])
			{
				TempK = Temptm[k][i][j]+273.15;
				work = solcoeffs[0] + solcoeffs[1]*(100/TempK) +
					solcoeffs[2]*log( TempK/100 ) + solcoeffs[3]*pow(TempK/100,2) +
					Salttm[k][i][j]*(solcoeffs[4]+solcoeffs[5]*(TempK/100)+solcoeffs[6]*pow(TempK/100,2));
				cfc11_sol[k][i][j] = exp(work);
				if (k<NML) cfc11_sat[i][j] = cfc11_sol[k][i][j]*cfc11_atmconc[i][j];
			}
			else {
				cfc11_sol[k][i][j] = 0.0;
				cfc11_sat[i][j] = 0.0;
			}

		}

}

void cfc11_find_atmconc(  ) {

	int i,j;
	extern double hlat[NYMEM];
	const double equatorbound[2] = {10,-10}; // Set the latitudes where to start interpolating atmospheric concentrations
	extern struct timekeeper_t timekeeper;
	double hemisphere_concentrations[2];
	double time[NUMATMVALS], atmval[NUMATMVALS];
	double currtime;

	currtime = timekeeper.current_time;
	
// Interpolate in time to find the atmospheric concentration
	for (i=0;i<NUMATMVALS;i++) {
		time[i] = atmconc[mCFC11].time[i];
		atmval[i] = atmconc[mCFC11].nval[i];

	}
	hemisphere_concentrations[0] = linear_interpolation(
			atmconc[CFC11IDX].time, atmconc[CFC11IDX].nval, currtime, NUMATMVALS);
	hemisphere_concentrations[1] = linear_interpolation(
			atmconc[CFC11IDX].time, atmconc[CFC11IDX].sval, currtime, NUMATMVALS);

	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++) {
			if(hlat[j] < equatorbound[0] && hlat[j]	> equatorbound[1]) {
				cfc11_atmconc[i][j] = linear_interpolation(
						equatorbound,hemisphere_concentrations,hlat[j],2);
			}
			if (hlat[j] > equatorbound[0]) {
				cfc11_atmconc[i][j] = hemisphere_concentrations[0];
			}
			if (hlat[j]<equatorbound[1] ) {
				cfc11_atmconc[i][j] = hemisphere_concentrations[1];
			}
		}

}

void surface_cfc11( ) {

	int i,j,k;
	extern double ***Temptm, ***Salttm;
	extern struct timekeeper_t timekeeper;

	printf("Setting CFC-11 surface condition\n");
	// Set cfc11 values to saturation at the mixed layer to mimic equilibrium with the atmosphere
	//
	cfc11_find_atmconc( );
	cfc11_saturation(  );
	
	printf("CFC-11 Sample Point\n");
	printf("\tTime: %f\n",timekeeper.current_time);
	printf("\tAtmospheric Concentration: %f\n",cfc11_atmconc[100][100]);
	printf("\tSalinity: %f\t Temperature: %f\n",Salttm[0][100][100],Temptm[0][100][100]);
	printf("\tSaturation concentration: %f\n\n",cfc11_sat[100][100]);


#ifdef NOCONC
	for (k=0;k<NML;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				tr[mCFC11][k][i][j] = cfc11_atmconc[i][j];

#else
	for (k=0;k<NML;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++)
				tr[mCFC11][k][i][j]=cfc11_sat[i][j];

	for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
			mn_cfc11sat[i][j]+=cfc11_sat[i][j];
#endif 
}

