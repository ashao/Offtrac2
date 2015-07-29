/*
 * gas_exchange.c
 * 17 April 2013
 *
 * MAIN FUNCTION: gas_exchange
 * Calculates the mixed layer concentration of a gas via 1D air-sea gas exchange
 * using an RK2 method with bubble injection all according to
 * Stanley et al. 2009
 *
 * Inputs:
 * 	Csurf: Tracer concentration AFTER transport
 * 	h_start: Beginning of month mixed layer thickness
 * 	h_end: End of month mixed layer thickness
 * 	u10_start: Beginning of month 10m windspped
 * 	u10_end: End of month windspeed
 * 	T_start: Beginning of month SST
 * 	T_end: End of month SST
 * 	S_start: Beginning of month SSS
 * 	S_end: End of month SSS
 * 	dt: Timestep to use for gas exchange integration
 * 	nt: number of timesteps
 *
 *  REFERENCES:
 *
 */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "netcdf.h"
#include "init.h"
#include "timekeeper.h"

#include "gas_exchange.h"
extern struct timekeeper_t timekeeper;

struct ocmip_t {

	double *time;
	double maxtime;
	double mintime;
	double ***kw;
	double ***atmpres;
	double ***fice;

};

struct ocmip_t ocmip;
double curr_fice[NXMEM][NYMEM];
double curr_kw[NXMEM][NYMEM];
double curr_atmpres[NXMEM][NYMEM];

const int nmonths = 12;
void gas_exchange( int tridx, double Sc_coeffs[4], double ***sat) {

	int iter;
	int i,j,k;
	double hml, T, S;;
	double Sc;
	double Csat, Cout;
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
			Cout = tr[tridx][0][i][j];
			T = Temptm[0][i][j];
			Sc = Sc_coeffs[0]-Sc_coeffs[1]*T+Sc_coeffs[2]*pow(T,2)-Sc_coeffs[3]*pow(T,3);
			hml = 0;

			for(k=0;k<NML;k++)
				hml += hend[k][i][j]; // Calculate depth of mixed layer

			for(iter=0; iter<gas_nt; iter++) // Perform gas exchange
				Cout += (1-curr_fice[i][j])*curr_kw[i][j]*pow(Sc/660.0,-0.5)*
				(sat[i][j]*curr_atmpres-Cout)/hml*dt_gas;

			for(k=0;k<NML;k++) // Set the mixed layer to the new value
				tr[tridx][k][i][j] = Cout;

		}
}

void read_gas_exchange_fields(char inpath[200]) {
	// Here all the gas exchange fields are read in for all time so that they can  be linearly interpolated if needed
	int i,j,t;


	FILE *file;
	int status;
	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];
	int err, cdfid, timeid, varid;

	char filepath[400];


	strcpy(filepath,inpath);
	strcat(filepath,"gasx_himgrid.nc");

	ocmip.kw = alloc3d(nmonths,NXMEM,NYMEM);
	ocmip.fice = alloc3d(nmonths,NXMEM,NYMEM);
	ocmip.atmpres = alloc3d(nmonths,NXMEM,NYMEM);
	ocmip.time = (double *)malloc(nmonths * sizeof(double));

	// Use funcitons in read.c to read in the forcing
	for (t=0;t<nmonths;t++) {

		read_var2d_time(filepath,t,"OCMIP_ATMP",ocmip.atmpres[t]);
		read_var2d_time(filepath,t,"OCMIP_FICE",ocmip.fice[t]);
		read_var2d_time(filepath,t,"wanninkhof_xkw",ocmip.kw[t]);

	}

	// Use netcdf function to read in time;
	start[0] = 0;
	count[0] = nmonths;

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if ((status = nc_inq_varid(cdfid, "MONTH_REG", &varid)))
		ERR(status);

	if ((status = nc_get_vara_float(cdfid,varid,start,count,ocmip.time)))
		ERR(status);

	// Convert the time variable to a year fraction from hours elapsed from start of year
	for (t=0;t<nmonths;t++)
		ocmip.time[t] = ocmip.time/24./365.;

	ocmip.maxtime = ocmip.time[nmonths-1];
	ocmip.mintime = ocmip.time[0];

	close_file(&cdfid,&file);


}

void update_gas_exchange_fields( ) {
	// Interpolate gas exchange fields to the current timestep
	int i,j;
	double fact0, fact1, factsum; // Interpolation factors
	int idx0, idx1; // Indices of the month to interpolate from
	double curr_yearfrac;
	curr_yearfrac = timekeeper.current_time - (double) timekeeper.current_year;


	idx0 = 0;
	if (curr_yearfrac>ocmip.maxtime) {
		idx0 = nmonths-1;
		idx1 = 0;

		fact0 = ocmip.mintime+1.0-curr_yearfrac;
		fact1 = curr_yearfrac-ocmip.maxtime;
	}
	else if (curr_yearfrac<ocmip.mintime) {
		idx0 = 0;
		idx1 = nmonths-1;
		fact0 = curr_yearfrac - ocmip.maxtime + 1.;
		fact1 = ocmip.mintime-curr_yearfrac;
	}
	else {
		while (curr_yearfrac<ocmip.time[idx0])
			idx0++;

		idx1=idx0+1;
		fact1 = curr_yearfrac - ocmip.time[idx0];
		fact0 = ocmip.time[idx1] - curr_yearfrac;

	}

	factsum = fact0 + fact1;
	fact0 = fact0/factsum;
	fact1 = fact1/factsum;

	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++) {
			curr_fice[i][j] = ocmip.fice[idx0][i][j]*fact0 + ocmip.fice[idx1][i][j]*fact1;
			curr_kw[i][j] = ocmip.kw[idx0][i][j]*fact0 + ocmip.kw[idx1][i][j]*fact1;
			curr_atmpres[i][j] = ocmip.atmpres[idx0][i][j]*fact0 + ocmip.atmpres[idx1][i][j]*fact1;
		}

}



