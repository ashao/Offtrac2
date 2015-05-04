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
#ifdef AGE
#include "ideal_age.h"
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
#define MISVAL -1.e+6

//BX-e
struct vardesc vars[NOVARS] =
{
		{
				"mn_h", "Layer Thickness", 'h', 'L', 's', "meter", 'd', 9.e-10
		}, // 3
		{
				"mn_uhtm",
				"zonal vel - cumul",
				'u',
				'L',
				's',
				"meter second-1",
				'f',
				0
		}, // 4
		{
				"mn_vhtm",
				"merid vel - cumul",
				'v',
				'L',
				's',
				"meter second-1",
				'f',
				0
		}, // 5
		{
				"mn_age",
				"Boundary propagator-like tracer",
				'h',
				'L',
				's',
				"years",
				'd',
				-3.17097930758233E-14
		},// 9
		{
				"h_test",
				"Check for layer convergence",
				'h',
				'L',
				's',
				"mol",
				'd',
				-1e6	
		}, // 10
};

void alloc_fields(void);
void set_metrics(void);

/* Begin edit DT */
void step_fields(int iyear, int itts, int imon, int iterno); // modified ashao
/*  End edit DT  */

// mn_ refers to mean value for output interval, WRINT (see init.h)

double ***h;
double ***hend;
double ***hstart;
double ****tr;
double ***uhtm;
double ***vhtm;
double ***wd;

double ***mn_h;
double ***mn_uhtm;
double ***mn_vhtm;
double ***mn_wd;

double depth[NZ][NXMEM][NYMEM];
double rml[2][NXMEM][NYMEM];

double areagr[NXMEM][NYMEM];
double D[NXMEM][NYMEM];

double umask[NXMEM][NYMEM];
double vmask[NXMEM][NYMEM];
int oceanmask[NXMEM][NYMEM];

double hlat[NYMEM];
double qlat[NYMEM];

double dt;

double *var[NOVARS];
long varsize[NOVARS];
int flags[NOVARS];
int rflags[NOVARS];
char directory[75];
char fname[75];
#ifdef WRTTS
struct varcdfinfo varinfo[NOVARS];
int varmap[NOVARS];
int itts; /* tracer time step counter */
FILE *fn;
char output_filename[200];
double wrts;
int nvar = 0, cdfid, timeid[2];
size_t nrec = 0;
#endif

/*-------------------------------------------------------------------*
 *
 *     begin main executable
 *
 *-------------------------------------------------------------------*/

int main(void)
{
int itts; /* tracer time step counter */
double smon, snxt;
int ismon, isnxt, ihnxt;

size_t nrec = 0;
int err, i, j, k;
double frac;
static int m;
int varmap[NOVARS];

FILE *fn;
char infile[200];
struct vardesc var_out[NOVARS];
struct varcdfinfo varinfo[NOVARS];
int nvar = 0, cdfid, timeid[2];

extern int flags[NOVARS];
extern int rflags[NOVARS];

//BX-a  for testing only
int status;
char message[144];
//BX-e

//BX  allocate tracer fields
err = alloc_arrays();
if (err)
	printf("Error allocating arrays.\n");

err = alloc_trac();
if (err)
	printf("Error allocating tracer field.\n");


/*
	printf("\ninitial month = %g \n", inmon);
	printf("final month = %g \n", inmon + tmon - 1);
	printf("total months = %g \n\n", tmon);
*/
	/*-----------------------------------
	 *
	 *     set output print flags to 0
	 *
	 *     added restart flags as a restart bug fix until
	 *     memory restriction problem is solved 31OCT07 BX
	 *
	 *----------------------------------*/

	for (i = 0; i <= NOVARS - 1; i++)
		flags[i] = 0;
	for (i = 0; i <= NOVARS - 1; i++)
		rflags[i] = 0;

	flags[1] = 0;
	flags[2] = 0; /* u,v */
	rflags[1] = 0;
	rflags[2] = 0; /* u,v */
	flags[0] = 0;
	flags[3] = 0; /* D, h */
	rflags[0] = 0;
	rflags[3] = 0; /* D, h */
	flags[4] = 0;
	flags[5] = 0; /* uhtm, vhtm */
	rflags[4] = 0;
	rflags[5] = 0; /* uhtm, vhtm */
	flags[6] = 0;
	flags[7] = 0;
	flags[8] = 0; /* ea, eb, eaml */
#ifdef ENTRAIN
	rflags[6]=1; rflags[7]=1; rflags[8]=1; /* ea, eb, eaml */
#endif
	flags[18] = 0; /* ML potential density */
	rflags[18] = 0; /* ML potential density */
#ifdef AGE
	flags[9] = 1;
	rflags[9] = 1; /* ideal age tracer*/
	flags[10] = 1; // htest
#endif
#ifdef OXYGEN
	flags[10] = 1; // Oxygen tracer
	flags[11] = 1; // jo2
	flags[12] = 1; // oxygen sat
	flags[13] = 1; // mn_dop
	flags[14] = 1; // mn_phosphate
	flags[15] = 1;
	rflags[10] = 1; // oxygen
	rflags[13] = 1; // dop
	rflags[14] = 1; // phosphate

#endif


/*-----------------------------------
 *
 *     allocate and initialize fields
 *
 *----------------------------------*/

read_grid();
printf("Done reading grid or metric file.\n");

set_metrics();
printf("Done setting metrics.\n");
printf("Reading bathymetry, D.\n");

// initialize D to be all ocean first

printf("calling read_D\n");
read_D();
initializemasks();

/* Copy the variable descriptions to a list of the actual output variables. */
for (i = 0; i < NOVARS; i++)
	if (flags[i] > 0)
	{
		var_out[nvar] = vars[i];
		varmap[i] = nvar;
		nvar++;
	}
// force float precision output with last argument
printf("Making NETCDF %s file\n", output_filename);
create_file(output_filename, NETCDF_FILE, var_out, nvar, &fn, &cdfid,
		timeid, varinfo, 1);
// don't force
// create_file(output_filename, NETCDF_FILE, var_out, nvar, &fn, &cdfid, timeid, varinfo, 0);
printf("Closing file \n");
close_file(&cdfid, &fn);


// Allocate temperature and salinity data arrays to be read in
allocate_ts();

/* Allocate the memory for the fields to be calculated.		*/
printf("allocating fields\n");
alloc_fields();

/* initialize tracer pointers					*/

for (m = 0; m < NOVARS; m++)
{
	if (flags[m])
		for (k = 0; k < varsize[m]; k++)
			var[m][k] = 0.0;
}


// ashao: Read in next month's isopycnal thickness fields as 'hend'
// (will be copied at the beginning of the integration)
// Done this way so that if hindcasts are used the physical fields switch smoothly
// to/from climatological fields
//BX  files are not in regular order (uvw are midmonth and h starts with last month)

// for files in regular order (h before uvw) use code here
//BX
// read_uvw(imon,1);
//BX
// read_h(imon,inxt);

printf("\nSetting up and initializing\n");

/*-----------------------------------
 *
 *     write tracer initial conditions
 *
 *----------------------------------*/

printf("Writing initial condition variables out to netCDF\n\n",cmon);

//  open netcdf file for each writing
status = nc_open(output_filename, NC_WRITE, &cdfid);
if (status != NC_NOERR)
{
	strcpy(message,"Opening file"); strcat(message,output_filename);
	handle_error(message,status);
}

err = write_time(cdfid,fn,timeid[0],nrec, dy);
if (err == -1) printf("Error writing day.\n");

for (m=0;m<NOVARS;m++) if (flags[m])
{
	err = write_field(cdfid, fn, vars[m], varinfo[varmap[m]],
			nrec,var[m]);
	if (err == -1) printf("Error writing %s.\n",vars[m].name);
}
//  close file after each writing
close_file(&cdfid, &fn);
printf("netcdf record = %d\n",nrec++);

/*-------------------------------------------------------------------*
 *
 *     begin integration loop
 *
 *-------------------------------------------------------------------*/

for (cmon = inmon; cmon < inmon + tmon; cmon++)
{
	// Update timekeeper
	// Copy the thickness fields to ensure continuity
	// Read in the current month's transport, thicknesses, and auxiliary fields
	printf("Month- hstart:%d hend:%d\n",ismon,isnxt);
	/*-----------------------------------
	 *
	 *     integrate 1 time step
	 *
	 *----------------------------------*/
	step_fields(iyear, itts, isnxt, iterno); // modified ashao

	/*-----------------------------------
	 *
	 *     write tracer averages and reset to 0
	 *
	 *----------------------------------*/
	printf("Writing month %i variables out to netCDF\n\n", cmon);

	status = nc_open(output_filename, NC_WRITE, &cdfid);
	if (status != NC_NOERR)
	{
		strcpy(message, "Opening file");
		strcat(message, output_filename);
		handle_error(message, status);
	}

	err = write_time(cdfid, fn, timeid[0], nrec, dy);
	if (err == -1)
		printf("Error writing day.\n");

	for (m = 0; m < NOVARS; m++)
		if (flags[m]==1)
		{
			printf("m = %d\n",m);
			err = write_field(cdfid, fn, vars[m],
			varinfo[varmap[m]], nrec, var[m]);
			if (err == -1)
				printf("Error writing %s.\n", vars[m].name);
		}
		close_file(&cdfid, &fn);

		printf("netcdf record = %d\n", nrec + 1);
		nrec++;

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

printf("start of array copying\n");

//  Second, create restart file name and open file

sprintf(restart_filename, "restart.%s.%04d.nc", run_name, cmon);
printf("Writing NetCDF restart file '%s'.\n\n", restart_filename);

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
create_file(restart_filename, NETCDF_FILE, var_out, nvar, &fn, &cdfid,
		timeid, varinfo, 0);

for (m = 0; m < NOVARS; m++)
	if (rflags[m])
	{
		err = write_field(cdfid, &fn, vars[m], varinfo[varmap[m]], 0,
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

#ifdef WRTTS
void write_ts(double wrts)
{
	int i,j,k,m;
	int status;
	char message[144];
	int err;
	double *wrdy;
#ifdef WRTTS
	wrdy = malloc(sizeof(double));
#endif
	/********** set means to zero                                   */
	/*  2d variables */
	if (flags[8]) set_fix_darray2d_zero(mn_eaml);
	/*---------------------------------------
	 *
	 *     write tracer values on output field
	 *
	 *--------------------------------------*/

	printf("write output field\n");

	if (flags[8]) add_fix_darray2d(mn_eaml, eaml);
	//HF outer if, inner for loops, instead of vice versa!
	//HF think about writing appropriate subroutine(s)!!!
	if (flags[1]) copy_darray3d(mn_u, u, NZ, NXMEM, NYMEM);
	/*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_u[k][i][j] = u[k][i][j];
     } */
	if (flags[2]) copy_darray3d(mn_v, v, NZ, NXMEM, NYMEM);
	/* {
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_v[k][i][j] = v[k][i][j];
     }*/
	if (flags[3])
	{
		for (k=0;k<NZ;k++)
			for (i=0;i<NXMEM;i++)
				for (j=0;j<NYMEM;j++)
					mn_h[k][i][j] = h[k][i][j];
	}
	if (flags[4]) copy_darray3d(mn_uhtm, uhtm, NZ, NXMEM, NYMEM);
	/*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_uhtm[k][i][j] = uhtm[k][i][j];
     }*/
	if (flags[5]) copy_darray3d(mn_vhtm, vhtm, NZ, NXMEM, NYMEM);
	/*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_vhtm[k][i][j] = vhtm[k][i][j];
     }*/
	if (flags[6])
	{
		for (k=0;k<NZ;k++)
			for (i=0;i<NXMEM;i++)
				for (j=0;j<NYMEM;j++)
					mn_ea[k][i][j] = ea[k][i][j];
	}
	if (flags[7])
	{
		for (k=0;k<NZ;k++)
			for (i=0;i<NXMEM;i++)
				for (j=0;j<NYMEM;j++)
					mn_eb[k][i][j] = Salttm[k][i][j];
	}
	if (flags[18])
	{
		for (k=0;k<2;k++)
			for (i=0;i<NXMEM;i++)
				for (j=0;j<NYMEM;j++)
					mn_rml[k][i][j] = rml[k][i][j];
	}
	printf("Writing variables for sub timestep %i out to netCDF\n\n",itts);
	//  open netcdf file for each writing
	status = nc_open(output_filename, NC_WRITE, &cdfid);
	if (status != NC_NOERR)
	{
		strcpy(message,"Opening file"); strcat(message,output_filename);
		handle_error(message,status);
	}
	*wrdy = wrts;
	err = write_time(cdfid,fn,timeid[0],nrec, wrdy);
	if (err == -1) printf("Error writing day.\n");

	for (m=0;m<NOVARS;m++) if (flags[m])
	{
		err = write_field(cdfid, fn, vars[m], varinfo[varmap[m]],
				nrec,var[m]);
		if (err == -1) printf("Error writing %s.\n",vars[m].name);
	}
	close_file(&cdfid, &fn);

	printf("netcdf record = %d\n",nrec++);
}
#endif //WRTTS
void alloc_fields(void)
{

	int m;

	extern double *var[NOVARS];
	extern double junk[(NZ + 1) * (NXMEM) * (NYMEM)];
	extern long varsize[NOVARS];
	extern int flags[NOVARS];

	for (m = 0; m < NOVARS; m++)
	{

		//HF if ( m>3 )    {
		//HF: added mn_h
		if (m >= 3)
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
		else
		{
			var[m] = junk;
			varsize[m] = 0;
		}
	}

	var[0] = &D[0][0];
	var[1] = &mn_u[0][0][0];
	var[2] = &mn_v[0][0][0];
	var[3] = &mn_h[0][0][0];
	var[4] = &mn_uhtm[0][0][0];
	var[5] = &mn_vhtm[0][0][0];
	var[6] = &mn_ea[0][0][0];
	var[7] = &mn_eb[0][0][0];
	var[8] = &mn_eaml[0][0];
	//var[18] = &mn_rml[0][0][0];

#ifdef AGE
var[9] = &mn_age[0][0][0];
var[10] = &htest[0][0][0];
#endif

#ifdef OXYGEN
var[10] = &mn_oxygen[0][0][0];
var[11] = &mn_jo2[0][0][0];
var[12] = &mn_o2sat[0][0][0];
var[13] = &mn_dop[0][0][0];
var[14] = &mn_phos[0][0][0];
var[15] = &mn_jpo4[0][0][0];

#endif

// end ashao

}
