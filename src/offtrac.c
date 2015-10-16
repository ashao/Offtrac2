/********+*********+*********+*********+*********+*********+*********+*
 *                                                                    *
 *                OFFTRAC - off-line tracer advection code            *
 *                                                                    *
 *                    David Darr - April/May 2002                     *
 *                                                                    *
 *                    OCMIP biogeochemistry added                     *
 *                     cdeutsch 2004                                  *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <netcdf.h>

#include "offtrac.h"
#include "init.h"
#define MAINFILE
#include "metrics.h"
#include "io.h"
#include "util.h"
#include "alloc_arrays.h"
#include "alloc_trac.h"
#include "initialize.h"
#include "iocdf.h"
#include "read.h"
#include "tracer_utilities.h"
#include "timekeeper.h"
#include "output_variables.h"
#include "step.h"
#include "gas_exchange.h"

#include "ideal_age.h"
#include "cfcs_sf6.h"
#include "ttd_bp.h"
#include "n2_module.h"
#include "ar_module.h"
#include "oxygen.h"
#include "phosphate.h"
/*-------------------------------------------------------------------*
 *                                                                   *
 *
 *     define variables and subroutines
 *                                                                   *
 *-------------------------------------------------------------------*/

/*  To add a variable, increase NOVARS (in init.h)  */
/*  then add information to struct vardesc (below) */
//BX  31AUG07  added missing value - last fields
//BX  31OCT07  Changed some variables to double precision (needed for restart)
//BX           All variables were float ('f') before
//BX           This is a bug fix as memory restrictions do not allow to create
//BX           a separate vardesc field for restarts.

//BX-a
//BX  missing values are set here and to the same value in "vars" declaration
//BX  they will be written on their respective mn_ fields with mult_.._mv calls
//HF the parameter and the macro need to be the same!
const double misval = -1.e+6;

// mn_ refers to mean value for output interval, WRINT (see init.h)
// Fields from HIM/GOLD/MOM6 output

double ***h, ***hstart, ***hend;
double ***depth;
double ***htest;
double ***uhtm;
double ***vhtm;
double ***wd;
double ***Temptm, ***Salttm;

double ***mn_h;
double ***mn_uhtm;
double ***mn_vhtm;
double ***mn_wd;

// Metrics related variables
double areagr[NXMEM][NYMEM];
double D[NXMEM][NYMEM];
double umask[NXMEM][NYMEM], vmask[NXMEM][NYMEM];
int oceanmask[NXMEM][NYMEM];
double hlat[NYMEM];
double qlat[NYMEM];
double **geolat;
double **geolon;

// Offtrac related variables
double ****tr; // Main tracer array

struct vardesc vars[NOVARS];

double *var[NOVARS];
long varsize[NOVARS];
int flags[NOVARS];
int rflags[NOVARS];

struct parameters run_parameters;
struct timekeeper_t timekeeper;
/*-------------------------------------------------------------------*
 *
 *     begin main executable
 *
 *-------------------------------------------------------------------*/

int main( int argc, char *argv[] )
{

	int err, i;
	static int m;

	int varmap[NOVARS];
	FILE *fn;

	struct vardesc var_out[NOVARS];
	struct varcdfinfo varinfo[NOVARS];
	int nvar = 0, cdfid, timeid[2];


	double *timeptr;

	//BX-a  for testing only
	int status;
	char message[144];
	//BX-e

	for (i = 0; i <= NOVARS - 1; i++)
		flags[i] = 0;
	for (i = 0; i <= NOVARS - 1; i++)
		rflags[i] = 0;

	strcpy(run_parameters.namelist_file,argv[1]);
	set_run_parameters(  );

	//BX  allocate tracer fields
	err = alloc_arrays();
	err = alloc_trac(run_parameters.tracer_counter);


	/*-----------------------------------
	 *
	 *     allocate and initialize fields
	 *
	 *----------------------------------*/

	read_grid();
	set_metrics();
	printf("areagr[100][100] = %f\n",areagr[100][100]);
	printf("Done setting metrics.\n");
	printf("Reading bathymetry, D.\n");
	read_D();
	initializemasks();

	initialize_timekeeper( );
	read_gas_exchange_fields(run_parameters.forcing_path);
	update_transport_fields( ); // Primarily done to read in the initial layer thicknesses
	printf("Read in fields\n");
	initialize();


	/* Copy the variable descriptions to a list of the actual output variables. */
	for (i = 0; i < NOVARS; i++)
		if (flags[i] > 0)
		{
			var_out[nvar] = vars[i];
			varmap[i] = nvar;
			nvar++;
		}

	// force float precision output with last argument
	printf("Making NETCDF %s file\n", run_parameters.outputfile);
	create_file(run_parameters.outputfile, NETCDF_FILE, var_out, nvar, &fn, &cdfid,
			timeid, varinfo, 1, timekeeper.write_intervals);
	// don't force
	// create_file(output_filename, NETCDF_FILE, var_out, nvar, &fn, &cdfid, timeid, varinfo, 0);
	printf("Closing file \n");
	close_file(&cdfid, &fn);

	/* Allocate the memory for the fields to be calculated.		*/
	printf("allocating fields\n");
	alloc_fields();

	/* initialize tracer pointers					*/
	/*
	for (m = 0; m < NOVARS; m++)
	{
		if (flags[m])
			for (k = 0; k < varsize[m]; k++)
				var[m][k] = 0.0;
	}
	 */


	/*-------------------------------------------------------------------*
	 *
	 *     begin integration loop
	 *
	 *-------------------------------------------------------------------*/

	for (timekeeper.iteration_counter = 0;
			timekeeper.iteration_counter < timekeeper.total_intervals;
			timekeeper.iteration_counter++)
	{

		/*-----------------------------------
		 *
		 *     integrate 1 time step
		 *
		 *----------------------------------*/
		update_timekeeper( );
		//		printf("Iteration counter: %d\n",timekeeper.iteration_counter);
		printf("Iteration: %d Year: %d Interval: %d Ending timestamp: %f\n",
				timekeeper.iteration_counter,timekeeper.current_year,
				timekeeper.current_interval,timekeeper.current_time);
		printf("UHTM: %f VHTM: %f WD: %e\n",uhtm[2][100][100],vhtm[2][100][100],wd[2][102][102]);
		step_fields( );


		/*-------------------------------------------------------------------*
		 *
		 *     end integration loop
		 *
		 *-------------------------------------------------------------------*/

		/*-----------------------------------
		 *
		 *     calculate tracer averages
		 *
		 *----------------------------------*/

		if (timekeeper.averaging_counter == run_parameters.wrint)
		{
			/*-----------------------------------
			 *
			 *     write tracer averages and reset to 0
			 *
			 *----------------------------------*/
			printf("Writing time interval %i variables out to netCDF\n\n", timekeeper.current_interval);
			printf("Accumulated averaging time (days): %f\n",timekeeper.accumulated_time_since_writing/86400);
			status = nc_open(run_parameters.outputfile, NC_WRITE, &cdfid);
			if (status != NC_NOERR)
			{
				strcpy(message, "Opening file");
				strcat(message, run_parameters.outputfile);
				handle_error(message, status);
			}

			timeptr = &timekeeper.current_time;
			err = write_time(cdfid, fn, timeid[0], timekeeper.num_records, timeptr);
			if (err == -1)
				printf("Error writing day.\n");

			for (m = 0; m < NOVARS; m++)
				if (flags[m]==1)
				{
					//					printf("m = %d\n",m);
					err = write_field(cdfid, fn, vars[m],
							varinfo[varmap[m]], timekeeper.num_records, var[m]);
					if (err == -1)
						printf("Error writing %s.\n", vars[m].name);
				}
			close_file(&cdfid, &fn);

			/*    reset all means to zero */
			timekeeper.averaging_counter = 0;
			timekeeper.accumulated_time_since_writing = 0;	
			set_darray3d_zero(mn_uhtm, NZ, NXMEM, NYMEM);
			set_darray3d_zero(mn_vhtm, NZ, NXMEM, NYMEM);
			set_darray3d_zero(mn_wd, NZ+1, NXMEM, NYMEM);
			set_darray3d_zero(mn_h, NZ, NXMEM, NYMEM);

			if (run_parameters.do_age)	set_darray3d_zero(mn_age, NZ, NXMEM, NYMEM);
			if (run_parameters.do_cfcs) {
				set_darray3d_zero(mn_cfc11, NZ, NXMEM, NYMEM);
				set_darray3d_zero(mn_pcfc11, NZ, NXMEM, NYMEM);
				set_darray3d_zero(mn_cfc12, NZ, NXMEM, NYMEM);
				set_darray3d_zero(mn_pcfc12, NZ, NXMEM, NYMEM);
				set_darray3d_zero(mn_sf6, NZ, NXMEM, NYMEM);
				set_darray3d_zero(mn_psf6, NZ, NXMEM, NYMEM);
			}
			if (run_parameters.do_ttd)	set_darray3d_zero(mn_ttd, NZ, NXMEM, NYMEM);

			if (run_parameters.do_n2) {
				set_darray3d_zero(mn_n2, NZ, NXMEM, NYMEM);
				set_darray3d_zero(mn_n2sat, NZ, NXMEM, NYMEM);
			}
			if (run_parameters.do_ar) {
				set_darray3d_zero(mn_ar, NZ, NXMEM, NYMEM);
				set_darray3d_zero(mn_arsat, NZ, NXMEM, NYMEM);
			}
			if (run_parameters.conservation_check) {
				set_darray3d_zero(mn_test, NZ, NXMEM, NYMEM);
				set_darray3d_zero(mn_htest, NZ, NXMEM, NYMEM);
			}
			//			printf("netcdf record = %d\n", timekeeper.num_records + 1);
			timekeeper.num_records++;

		} /*  end if nmn==WRITEINT */

		/*-------------------------------------------------------------------*
		 *
		 *     end integration loop
		 *
		 *-------------------------------------------------------------------*/

	} /* end while */

	//BX-a
	/*-----------------------------------
	 *
	 *     write restart file
	 *
	 *----------------------------------*/

	//   First write the up-to-date tracer values to the mean fields.
	//   In order to save memory the instantaneous values for the
	//   restart will be written on mean fields in the restart file.


	//  Second, create restart file name and open file

	printf("Writing NetCDF restart file '%s'.\n\n", run_parameters.new_restartfile);

	// Force layer thicknesses to be output for the restart
	var[map_variable_to_index("hlay")] = &hend[0][0][0];
	rflags[map_variable_to_index("hlay")] = 1;

	if (run_parameters.do_age) copy_darray3d(mn_age,tr[mAGE],NZ,NXMEM,NYMEM);
	if (run_parameters.do_ttd) copy_darray3d(mn_ttd,tr[mTTD],NZ,NXMEM,NYMEM);
	if (run_parameters.do_cfcs)
		{
			copy_darray3d(mn_cfc11,tr[mCFC11],NZ,NXMEM,NYMEM);
			copy_darray3d(mn_cfc12,tr[mCFC12],NZ,NXMEM,NYMEM);
			copy_darray3d(mn_sf6,tr[mSF6],NZ,NXMEM,NYMEM);
		}
	if (run_parameters.do_n2) copy_darray3d(mn_n2,tr[mN2],NZ,NXMEM,NYMEM);
	if (run_parameters.do_ar) copy_darray3d(mn_ar,tr[mAR],NZ,NXMEM,NYMEM);
	if (run_parameters.do_oxygen) {
		copy_darray3d(mn_oxygen,tr[mOXYGEN],NZ,NXMEM,NYMEM);
		copy_darray3d(mn_phos,tr[mPHOSPHATE],NZ,NXMEM,NYMEM);
		copy_darray3d(mn_dop,tr[mDOP],NZ,NXMEM,NYMEM);
	}

	/* Copy the variable descriptions to a list of the actual restart variables. */
	nvar = 0;
	for (i = 0; i < NOVARS; i++)
		if (rflags[i] > 0)
		{
			var_out[nvar] = vars[i];
			varmap[i] = nvar;
			nvar++;
		}

	// do NOT force float precision output with last argument
	create_file(run_parameters.new_restartfile, NETCDF_FILE, var_out, nvar, &fn, &cdfid,
			timeid, varinfo, 0, 0);


	for (m = 0; m < NOVARS; m++)
		if (rflags[m])
		{
			err = write_field(cdfid, fn, vars[m], varinfo[varmap[m]], 0,
					var[m]);
			if (err == -1)
				printf("Error writing %s.\n", vars[m].name);
		}

	close_file(&cdfid, &fn);

	printf("\n programme termine normallement. \n");
	return (0);
}

/*-------------------------------------------------------------------*
 *
 *     begin subroutines
 *                                                                   *
 *-------------------------------------------------------------------*/

void alloc_fields(void)
{

	int m;

	extern double *var[NOVARS];
	extern long varsize[NOVARS];
	extern int flags[NOVARS];

	for (m = 0; m < NOVARS; m++)
	{

		//HF if ( m>3 )    {
		//HF: added mn_h
		//		if (m >= 3)
		{
			switch (vars[m].z_grid)
			{
			case 'L':
				varsize[m] = NZ * (NXMEM) * (NYMEM);
				break;
			case 'i':
				varsize[m] = (NZ + 1) * (NXMEM) * (NYMEM);
				break;
			case '2':
				varsize[m] = 2 * (NXMEM) * (NYMEM);
				break;
			case '1':
				varsize[m] = (NXMEM) * (NYMEM);
				break;
			case 'n':
				varsize[m] = 1;
				break;
			default:
				printf("Unknown layer axis[%d] %c.\n", m, vars[m].z_grid);
			}

			var[m] = calloc((size_t) varsize[m], sizeof(double));
			if (var[m] == NULL)
			{
				printf("Unable to allocate memory for var[%d].\n", m);
				exit(-21);
			}

			//	temporarily disabling this output
			//      printf("Allocated memory for var[%d],%s.\n\n",m,vars[m].name);

		}
		//		else
		//		{
		varsize[m] = 0;
		//		}
	}

	/*	var[map_variable_to_index("D")] = &D[0][0];
	var[map_variable_to_index("geolat")] = &geolat[0][0];
	var[map_variable_to_index("geolon")] = &geolon[0][0];
	var[map_variable_to_index("wetmask")] = &geolon[0][0];
	 */
	var[map_variable_to_index("uhtm")] = &mn_uhtm[0][0][0];
	var[map_variable_to_index("vhtm")] = &mn_vhtm[0][0][0];
	var[map_variable_to_index("wd")] = &mn_wd[0][0][0];
	var[map_variable_to_index("hlay")] = &mn_h[0][0][0];

	if (run_parameters.do_age)	var[map_variable_to_index("age")] = &mn_age[0][0][0];

	if (run_parameters.do_cfcs) {
		var[map_variable_to_index("cfc11")] = &mn_cfc11[0][0][0];
		var[map_variable_to_index("pcfc11")] = &mn_pcfc11[0][0][0];
		var[map_variable_to_index("cfc12")] = &mn_cfc12[0][0][0];
		var[map_variable_to_index("pcfc12")] = &mn_pcfc12[0][0][0];
		var[map_variable_to_index("sf6")] = &mn_sf6[0][0][0];
		var[map_variable_to_index("psf6")] = &mn_psf6[0][0][0];
	}
	if (run_parameters.do_ttd) var[map_variable_to_index("ttd")] = &mn_ttd[0][0][0];
	if (run_parameters.do_n2) {
		var[map_variable_to_index("n2")] = &mn_n2[0][0][0];
		var[map_variable_to_index("n2sol")] = &mn_n2sat[0][0][0];
	}
	if (run_parameters.do_ar) {
		var[map_variable_to_index("ar")] = &mn_ar[0][0][0];
		var[map_variable_to_index("arsol")] = &mn_arsat[0][0][0];
	}
	if (run_parameters.do_oxygen) {
		var[map_variable_to_index("oxygen")] = &mn_oxygen[0][0][0];
		var[map_variable_to_index("o2sat")] = &mn_o2sat[0][0][0];
		var[map_variable_to_index("jo2")] = &mn_jo2[0][0][0];
		var[map_variable_to_index("po4")] = &mn_phos[0][0][0];
		var[map_variable_to_index("dop")] = &mn_dop[0][0][0];
		var[map_variable_to_index("jpo4")] = &mn_jpo4[0][0][0];
	}

	if (run_parameters.conservation_check) {
		var[map_variable_to_index("test")] = &mn_test[0][0][0];
		var[map_variable_to_index("htest")] = &mn_htest[0][0][0];
	}
	//var[18] = &mn_rml[0][0][0];

	// end ashao

}
