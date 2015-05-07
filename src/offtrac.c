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

#ifdef CFCS
#include "cfcs_sf6.h"
#endif

#ifdef AGE

extern int mAGE;

#endif

/*-------------------------------------------------------------------*
 *                                                                   *
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
double depth[NZ][NXMEM][NYMEM];

double ***h, ***hstart, ***hend;
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
#ifdef AGE

extern double ***mn_age;

#endif
/*-------------------------------------------------------------------*
 *
 *     begin main executable
 *
 *-------------------------------------------------------------------*/

int main( int argc, char *argv[] )
{

	int err, i, j, k;
	double frac;
	static int m;

	int varmap[NOVARS];
	FILE *fn;

	char run_name[200],namelist_file[200];
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
	printf("Namelist file: %s\n",run_parameters.namelist_file);
	set_run_parameters(  );

	//BX  allocate tracer fields
	err = alloc_arrays();
	err = alloc_trac();


	/*-----------------------------------
	 *
	 *     allocate and initialize fields
	 *
	 *----------------------------------*/

	read_grid();
	set_metrics();
	printf("Done setting metrics.\n");
	printf("Reading bathymetry, D.\n");
	read_D();
	initializemasks();

	initialize_timekeeper( );
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
			timeid, varinfo, 1);
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
#ifdef AGE
			set_darray3d_zero(mn_age, NZ, NXMEM, NYMEM);
#endif

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

#ifdef AGE
	copy_darray3d(mn_age,tr[mAGE],NZ,NXMEM,NYMEM);
#endif

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
			timeid, varinfo, 0);


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
	extern double junk[(NZ + 1) * (NXMEM) * (NYMEM)];
	extern long varsize[NOVARS];
	extern int flags[NOVARS];
#ifdef AGE
	extern double ***mn_age;
#endif
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
#ifdef AGE
	var[map_variable_to_index("age")] = &mn_age[0][0][0];
#endif
#ifdef CFCS
	var[map_variable_to_index("cfc11")] = &mn_cfc11[0][0][0];
	var[map_variable_to_index("pcfc11")] = &mn_pcfc11[0][0][0];
	var[map_variable_to_index("cfc12")] = &mn_cfc12[0][0][0];
	var[map_variable_to_index("pcfc12")] = &mn_pcfc12[0][0][0];
	var[map_variable_to_index("sf6")] = &mn_sf6[0][0][0];
	var[map_variable_to_index("psf6")] = &mn_psf6[0][0][0];
#endif
	//var[18] = &mn_rml[0][0][0];

	// end ashao

}
