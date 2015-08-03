/*
 * gas_exchange.c
 * 17 April 2013
 *
 * MAIN FUNCTION: gas_exchange
 * Calculates the mixed layer concentration of a gas via 1D air-sea gas exchange
 * using the OCMIP-2 Protocol
 *
 *
 *  REFERENCES:
 *
 */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "alloc.h"
#include "netcdf.h"
#include "init.h"
#include "timekeeper.h"
#include "util.h"
#include "gas_exchange.h"
#include "tracer_utilities.h"
#include "io.h"
#include "read.h"

extern struct timekeeper_t timekeeper;
extern int oceanmask[NXMEM][NYMEM];

struct ocmip_t {

	float *time;
	float maxtime;
	float mintime;
	float ***kw;
	float ***atmpres;
	float ***fice;

};

struct ocmip_t ocmip;
double curr_fice[NXMEM][NYMEM];
double curr_kw[NXMEM][NYMEM];
double curr_atmpres[NXMEM][NYMEM];

const int nmonths = 12;
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

void gas_exchange( int tridx, const double Sc_coeffs[4], double **sat) {

	int iter;
	int i,j,k;
	double hml, T;;
	double Sc;
	double  Cout;
	double dt_gas;
	const int gas_nt = 30;

	extern double ***Temptm;;
	extern double ***hend;


	extern double ****tr;

	//	P=101325;
	//	u10 = 5;
	dt_gas = timekeeper.dt / gas_nt;
	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++){

			if (oceanmask[i][j]) {
				Cout = tr[tridx][0][i][j];
				T = Temptm[0][i][j];
				Sc = Sc_coeffs[0]-Sc_coeffs[1]*T+Sc_coeffs[2]*pow(T,2)-Sc_coeffs[3]*pow(T,3);
				hml = 0;

				for(k=0;k<NML;k++)
					hml += hend[k][i][j]; // Calculate depth of mixed layer

				for(iter=0; iter<gas_nt; iter++) // Perform gas exchange
					Cout += (1-curr_fice[i][j])*curr_kw[i][j]*pow(Sc/660.0,-0.5)*
					(sat[i][j]*curr_atmpres[i][j]-Cout)/hml*dt_gas;

				if (curr_atmpres[i][j] < 0.0)
					printf("Cout[%d][%d]: %f Temp: %f\n",i,j,Cout,T);

				for(k=0;k<NML;k++) // Set the mixed layer to the new value
					tr[tridx][k][i][j] = Cout;
			}
		}
}

void read_gas_exchange_fields(char *inpath) {
	// Here all the gas exchange fields are read in for all time so that they can  be linearly interpolated if needed
	int t;


	FILE *file;
	int status;
	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];
	int cdfid, timeid, varid;

	char filepath[400];

	strcpy(filepath,inpath);
	strcat(filepath,"gasx_himgrid.nc");

	ocmip.kw = alloc3d_f(nmonths,NXMEM,NYMEM);
	ocmip.fice = alloc3d_f(nmonths,NXMEM,NYMEM);
	ocmip.atmpres = alloc3d_f(nmonths,NXMEM,NYMEM);
	ocmip.time = (float *)malloc(nmonths * sizeof(float));

	// Use funcitons in read.c to read in the forcing
	for (t=0;t<nmonths;t++) {

		read_var2d_time(filepath,t,"OCMIP_ATMP",ocmip.atmpres[t]);
		read_var2d_time(filepath,t,"OCMIP_FICE",ocmip.fice[t]);
		read_var2d_time(filepath,t,"wanninkhof_xkw",ocmip.kw[t]);

	}

	printf("\tKw: %e Pressure: %e\n",(double) ocmip.kw[0][100][100], (double) ocmip.atmpres[0][100][100]);
	// Use netcdf function to read in time;
	start[0] = 0;
	count[0] = nmonths;

	open_input_file(filepath,&file,&cdfid,&timeid);
	if ((status = nc_inq_varid(cdfid, "MONTH_REG", &varid)))
		ERR(status);

	printf("Reading time variable\n");

	if ((status = nc_get_vara_float(cdfid,varid,start,count,&ocmip.time[0])))
		ERR(status);

	// Convert the time variable to a year fraction from hours elapsed from start of year
	printf("Converted time variable\n");
	for (t=0;t<nmonths;t++)
		ocmip.time[t] = ocmip.time[t]/24./365.;

	ocmip.maxtime = ocmip.time[nmonths-1];
	ocmip.mintime = ocmip.time[0];

	close_file(&cdfid,&file);


}

void update_gas_exchange_fields( ) {
	// Interpolate gas exchange fields to the current timestep
	int i,j;
	float t0, t1; // Interpolation times
	int idx0, idx1; // Indices of the month to interpolate from
	double curr_yearfrac;
	curr_yearfrac = timekeeper.current_time - (double) timekeeper.current_year;
	double tarr[2],datarr[2];
	
	idx0 = 0;
	idx1 = 1;
	if (curr_yearfrac>ocmip.maxtime) {
		idx0 = nmonths-1;
		idx1 = 0;
		t0 = ocmip.time[idx0];
		t1 = ocmip.time[idx1] + 1.0;

	}
	else if (curr_yearfrac<ocmip.mintime) {
		idx0 = 0;
		idx1 = nmonths-1;
		
		t0 = ocmip.time[idx0];
		t1 = ocmip.time[idx1]-1.0;
	}
	else {
		printf("curr_yearfrac: %f ocmip.time[%d]: %f",curr_yearfrac,idx1,ocmip.time[idx1]);
		while (curr_yearfrac>ocmip.time[idx1]) {
			idx0++;
			idx1=idx0+1;
			printf("idx0: %d idx1: %d\n",idx0,idx1);
		}
		t0 = ocmip.time[idx0];
		t1 = ocmip.time[idx1];
		printf("Interpolation Time Interval: %d %d %f %f\n",idx0,idx1,t0,t1);

	}

	tarr[0] = (double) t0;
	tarr[1] = (double) t1;

	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++) {
			
			datarr[0] = (double) ocmip.fice[idx0][i][j];
			datarr[1] = (double) ocmip.fice[idx1][i][j];
			curr_fice[i][j] = linear_interpolation(tarr,datarr,curr_yearfrac,2);
			datarr[0] = (double) ocmip.kw[idx0][i][j];
			datarr[1] = (double) ocmip.kw[idx1][i][j];
			curr_kw[i][j] = linear_interpolation(tarr,datarr,curr_yearfrac,2);
			datarr[0] = (double) ocmip.atmpres[idx0][i][j];
			datarr[1] = (double) ocmip.atmpres[idx1][i][j];
			curr_atmpres[i][j] = linear_interpolation(tarr,datarr,curr_yearfrac,2);
		}
	printf("Gas Exchange @ Year Fraction: %f:\n",curr_yearfrac);
	printf("\tKw: %e Pressure: %e\n",curr_kw[100][100],curr_atmpres[100][100]);

}



