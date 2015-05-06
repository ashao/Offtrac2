#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "netcdf.h"
#include "init.h"
#include "io.h"
#include "iocdf.h"
#include "alloc.h"
#include "read.h"
#include "step.h"
#include "util.h"
#include "metrics.h"
#include "initialize.h"
#include "timekeeper.h"

extern double areagr[NXMEM][NYMEM];

extern double dt;

extern double ***uhtm,***vhtm;
extern double ***wd;
extern double D[NXMEM][NYMEM];

extern double qlat[NYMEM],hlat[NYMEM];

extern struct parameters run_parameters;
extern struct timekeeper_t timekeeper;

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

void read_grid()
{
	int err, i, j, cdfid, timeid;
	int status, varid; // HF
	char infile[25], inpath[200];
	FILE *file;
	double latq[NYTOT];
	double lath[NYTOT];
	double Ah[NYTOT][NXTOT];
	//HF
	double dxh_in[NYTOT][NXTOT], dxq_in[NYTOT][NXTOT];
	double dxu_in[NYTOT][NXTOT], dxv_in[NYTOT][NXTOT];
	double dyh_in[NYTOT][NXTOT], dyq_in[NYTOT][NXTOT];
	double dyu_in[NYTOT][NXTOT], dyv_in[NYTOT][NXTOT];

	long  start[MAX_NC_VARS];
	long  end[MAX_NC_VARS];

	sprintf(infile,"metrics.nc");

	strcpy(inpath, run_parameters.forcing_path);
	strcat(inpath, infile);
	printf("\nLooking for file '%s'.\n",inpath);

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
		strcat(inpath, ".cdf");
		err = open_input_file(inpath,&file,&cdfid,&timeid);
		if (err != 0) {
			printf("Unable to find grid file.\n");
			exit(-73);
		}
	}

	end[0] = NYTOT;
	end[1] = NXTOT;

	printf("end %d %d\n",end[0],end[1]);
	//HF
	status = nc_inq_varid(cdfid, "latq", &varid);
	if (status != NC_NOERR) handle_error("inq latq", status);
	status = nc_get_vara_double(cdfid, varid, start, end, latq);
	if (status != NC_NOERR) handle_error("read latq", status);
	//BX-a  ncvarget(cdfid, 2, start, end, latq);

	//CAD should these loop indices be hardwired (to i=0;i<NYTOT) ?
	for (i=0;i<NYTOT;i++) {
		qlat[i+2]=(double)latq[i];
	}

	//HF ncvarget(cdfid, 2, start, end, lath);
	status = nc_inq_varid(cdfid, "lath", &varid);
	if (status != NC_NOERR) handle_error("inq lath", status);
	status = nc_get_vara_double(cdfid, varid, start, end, lath);
	if (status != NC_NOERR) handle_error("read lath", status);

	//CAD should these loop indices be hardwired (to i=0;i<NYTOT) ?
	for (i=0;i<NYTOT;i++) {
		hlat[i+2]=(double)lath[i];
	}
	//BX-e

	//   ncvarget(cdfid,10, start, end, dxv);
	status = nc_inq_varid(cdfid, "dxv", &varid);
	if (status != NC_NOERR) handle_error("inq dxv", status);
	status = nc_get_vara_double(cdfid, varid, start, end, &dxv_in[0][0]);
	if (status != NC_NOERR) handle_error("read dxv", status);

	//   ncvarget(cdfid,12, start, end, dyu);
	status = nc_inq_varid(cdfid, "dyu", &varid);
	if (status != NC_NOERR) handle_error("inq dyu", status);
	status = nc_get_vara_double(cdfid, varid, start, end, &dyu_in[0][0]);
	if (status != NC_NOERR) handle_error("read dyu", status);

	//HF     ncvarget(cdfid,13, start, end, Ah);
	status = nc_inq_varid(cdfid, "Ah", &varid);
	if (status != NC_NOERR) handle_error("inq Ah", status);
	status = nc_get_vara_double(cdfid, varid, start, end, &Ah[0][0]);
	if (status != NC_NOERR) handle_error("read Ah", status);

	//   ncvarget(cdfid,14, start, end, dyv);
	status = nc_inq_varid(cdfid, "dyv", &varid);
	if (status != NC_NOERR) handle_error("inq dyv", status);
	status = nc_get_vara_double(cdfid, varid, start, end, &dyv_in[0][0]);
	if (status != NC_NOERR) handle_error("read dyv", status);

	//   ncvarget(cdfid,15, start, end, dxh);
	status = nc_inq_varid(cdfid, "dxh", &varid);
	if (status != NC_NOERR) handle_error("inq dxh", status);
	status = nc_get_vara_double(cdfid, varid, start, end, &dxh_in[0][0]);
	if (status != NC_NOERR) handle_error("read dxh", status);

	//   ncvarget(cdfid,16, start, end, dyq);
	status = nc_inq_varid(cdfid, "dyq", &varid);
	if (status != NC_NOERR) handle_error("inq dyq", status);
	status = nc_get_vara_double(cdfid, varid, start, end, &dyq_in[0][0]);
	if (status != NC_NOERR) handle_error("read dyq", status);

	//   ncvarget(cdfid,17, start, end, dxq);
	status = nc_inq_varid(cdfid, "dxq", &varid);
	if (status != NC_NOERR) handle_error("inq dxq", status);
	status = nc_get_vara_double(cdfid, varid, start, end, &dxq_in[0][0]);
	if (status != NC_NOERR) handle_error("read dxq", status);

	//   ncvarget(cdfid,18, start, end, dyh);
	status = nc_inq_varid(cdfid, "dyh", &varid);
	if (status != NC_NOERR) handle_error("inq dyh", status);
	status = nc_get_vara_double(cdfid, varid, start, end, &dyh_in[0][0]);
	if (status != NC_NOERR) handle_error("read dyh", status);

	//   ncvarget(cdfid,19, start, end, dxu);
	status = nc_inq_varid(cdfid, "dxu", &varid);
	if (status != NC_NOERR) handle_error("inq dxu", status);
	status = nc_get_vara_double(cdfid, varid, start, end, &dxu_in[0][0]);
	if (status != NC_NOERR) handle_error("read dxu", status);


	for (i=0;i<NXTOT;i++) {
		for (j=0;j<NYTOT;j++) {
			areagr[i+2][j+2]= Ah[j][i];
			DXq(i+2,j+2) = dxq_in[j][i];
			DYq(i+2,j+2) = dyq_in[j][i];
			DXv(i+2,j+2) = dxv_in[j][i];
			DYv(i+2,j+2) = dyv_in[j][i];
			DXh(i+2,j+2) = dxh_in[j][i];
			DYh(i+2,j+2) = dyh_in[j][i];
			DXu(i+2,j+2) = dxu_in[j][i];
			DYu(i+2,j+2) = dyu_in[j][i];
		}
		// HF
		DXv(i+2,1) = DXv(i+2,2);
		DYv(i+2,1) = DYv(i+2,2);
		// HF-e
	}

	close_file(&cdfid,&file);
	printf("\nFinished reading file '%s'.\n",inpath);


	/* zonal re-entrance            */

	for (j=0;j<=NYMEM-1;j++) {
		areagr[nx+1][j] = areagr[2][j];
		areagr[nx+2][j] = areagr[3][j];
		areagr[0][j] = areagr[nx-1][j];
		areagr[1][j] = areagr[nx][j];

		DXq(nx+1,j) = DXq(2,j);
		DXq(nx+2,j) = DXq(3,j);
		DXq(0,j) = DXq(nx-1,j);
		DXq(1,j) = DXq(nx,j);

		DYq(nx+1,j) = DYq(2,j);
		DYq(nx+2,j) = DYq(3,j);
		DYq(0,j) = DYq(nx-1,j);
		DYq(1,j) = DYq(nx,j);

		DXv(nx+1,j) = DXv(2,j);
		DXv(nx+2,j) = DXv(3,j);
		DXv(0,j) = DXv(nx-1,j);
		DXv(1,j) = DXv(nx,j);

		DYv(nx+1,j) = DYv(2,j);
		DYv(nx+2,j) = DYv(3,j);
		DYv(0,j) = DYv(nx-1,j);
		DYv(1,j) = DYv(nx,j);

		DXh(nx+1,j) = DXh(2,j);
		DXh(nx+2,j) = DXh(3,j);
		DXh(0,j) = DXh(nx-1,j);
		DXh(1,j) = DXh(nx,j);

		DYh(nx+1,j) = DYh(2,j);
		DYh(nx+2,j) = DYh(3,j);
		DYh(0,j) = DYh(nx-1,j);
		DYh(1,j) = DYh(nx,j);

		DXu(nx+1,j) = DXu(2,j);
		DXu(nx+2,j) = DXu(3,j);
		DXu(0,j) = DXu(nx-1,j);
		DXu(1,j) = DXu(nx,j);

		DYu(nx+1,j) = DYu(2,j);
		DYu(nx+2,j) = DYu(3,j);
		DYu(0,j) = DYu(nx-1,j);
		DYu(1,j) = DYu(nx,j);
	}

	/* meridional re-entrance            */
#ifdef REENTRANT_Y
	for (i=2;i<=nx;i++) {
		int ii = 363 - i;
		areagr[ii][ny+1] = areagr[i][ny];
		areagr[ii][ny+2] = areagr[i][ny-1];

		DXq(ii,ny+1) = DXq(i,ny);
		DXq(ii,ny+2) = DXq(i,ny-1);

		DYq(ii,ny+1) = DYq(i,ny);
		DYq(ii,ny+2) = DYq(i,ny-1);

		DXv(ii,ny+1) = DXv(i,ny);
		DXv(ii,ny+2) = DXv(i,ny-1);

		DYv(ii,ny+1) = DYv(i,ny);
		DYv(ii,ny+2) = DYv(i,ny-1);

		DXh(ii,ny+1) = DXh(i,ny);
		DXh(ii,ny+2) = DXh(i,ny-1);

		DYh(ii,ny+1) = DYh(i,ny);
		DYh(ii,ny+2) = DYh(i,ny-1);

		DXu(ii,ny+1) = DXu(i,ny);
		DXu(ii,ny+2) = DXu(i,ny-1);

		DYu(ii,ny+1) = DYu(i,ny);
		DYu(ii,ny+2) = DYu(i,ny-1);
	}
#endif
}

void read_D()
{
	int i,ii,j;
	int err, cdfid, timeid;
	char infile[25], inpath[200];
	FILE *file;

	/*  Read the depth and other data in from the saved binary files.     */

	sprintf(infile,"topo.nc");

	strcpy(inpath, run_parameters.forcing_path);
	strcat(inpath, infile);
	printf("Looking for file '%s'.\n",inpath);

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
		strcat(inpath, ".cdf");
		err = open_input_file(inpath,&file,&cdfid,&timeid);
		if (err != 0) {
			printf("Unable to find depth file.\n");
			exit(-73);
		}
	}

	//HF  WARN read_field(cdfid,file,"D", NXTOT,NYTOT,1, 0,0,0, 1,1,0, 1, D[0]);
	read_field(cdfid,file,"D", NXTOT,NYTOT,1, 0,0,0, X1,Y1,0, 1, D[0]);
	close_file(&cdfid,&file);


	//HF10022009	zonal re-entrance
	for (j=0;j<=NYMEM-1;j++) {
		D[nx+1][j] = D[2][j];
		D[nx+2][j] = D[3][j];
		D[0][j] =   D[nx-1][j];
		D[1][j] =   D[nx][j];
	}
#ifdef REENTRANT_Y
	//      meridional re-entrance
	for (i=0;i<=nx+2;i++) {
		ii = 363 - i;
		D[ii][ny+1] = D[i][ny];
		D[ii][ny+2] = D[i][ny-1];
	}
#endif
}

void read_uvw(int imon, char *fieldtype, char *readpath)
{

	int i,j,k;
	int err, cdfid, timeid;
	char infile[25], inpath[1000];
	FILE *file;
	int status;
	int nmonths;
	double fact0,fact1,fact2;

	int uhfileid,vhfileid,wdfileid;

	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];

	float*** tmp3d;
	printf("UVW %s index: %d\n",fieldtype,imon);

	//
	//   read in separate files for U, V, and W
	//
	//    Start with uhtm

	sprintf(infile,"UH.%s.nc",fieldtype);
	strcpy(inpath, readpath);
	strcat(inpath, infile);

//	printf("Looking for file '%s'.\n",inpath);

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
		strcat(inpath, ".cdf");
		err = open_input_file(inpath,&file,&cdfid,&timeid);
		if (err != 0) {
			printf("Unable to find UH file.\n");
			exit(-73);
		}
	}


	if ((status = nc_inq_varid(cdfid, "uh", &uhfileid)))
		ERR(status);
	//    if ((status = nc_inq_varid(cdfid, "UHCLIM", &uhfileid)))
	//if ((status = nc_inq_varid(cdfid, "UH", &uhfileid)))
	//   ERR(status);

	bzero(start, MAX_NC_VARS * sizeof(long));

	count[0] = 1;
	count[1] = NZ;
	count[2] = NYTOT;
	count[3] = NXTOT;


	//    printf("read_clim UH month=  %i\n",imon);

	// allocate NZ+1 so that it can be used for WDCLIM as well
	tmp3d  = alloc3d_f(NZ+1,NYTOT,NXTOT);

	start[0] = imon;
	if ((status = nc_get_vara_float(cdfid,uhfileid,start,count,tmp3d[0][0])))
		ERR(status);
	//printf("read_clim UH month=  %i done\n",imon);
	//    start[0] = iprv;
	//    if ((status = nc_get_vara_float(cdfid,uhclimid,start,count,tmp3dp[0][0])))
	//       ERR(status);
	//printf("read_clim UH month=  %i done\n",iprv);
	//    start[0] = inxt;
	//    if ((status = nc_get_vara_float(cdfid,uhfileid,start,count,tmp3dn[0][0])))
	//       ERR(status);
	//printf("read_clim UH month=  %i done\n",inxt);

	for (k=0;k<NZ;k++)
		for (i=0;i<NXTOT;i++)
			for (j=0;j<NYTOT;j++)
				uhtm[k][i+2][j+2]= tmp3d[k][j][i];
	//uhtm[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
	//    fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

	close_file(&cdfid,&file);
	//printf("READ uhtm=%g;%g,%g,%g\n",uhtm[0][190][26],tmp3d[0][24][188],
	//	   tmp3dp[0][24][188],tmp3dn[0][24][188]);

	//    Next vhtm

	sprintf(infile,"VH.%s.nc",fieldtype);

	strcpy(inpath, readpath);
	strcat(inpath, infile);

//	printf("Looking for file '%s'.\n",inpath);

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
		strcat(inpath, ".cdf");
		err = open_input_file(inpath,&file,&cdfid,&timeid);
		if (err != 0) {
			printf("Unable to find VH file.\n");
			exit(-73);
		}
	}

	if ((status = nc_inq_varid(cdfid, "vh", &vhfileid)))
		ERR(status);
	//    if ((status = nc_inq_varid(cdfid, "VHCLIM", &vhfileid)))
	//if ((status = nc_inq_varid(cdfid, "VH", &vhclimid)))
	//   ERR(status);

	bzero(start, MAX_NC_VARS * sizeof(long));

	count[0] = 1;
	count[1] = NZ;
	count[2] = NYTOT;
	count[3] = NXTOT;

	// printf("read_clim VH month=  %i\n",imon);

	start[0] = imon;
	if ((status = nc_get_vara_float(cdfid,vhfileid,start,count,tmp3d[0][0])))
		ERR(status);
	//    start[0] = iprv;
	//    if ((status = nc_get_vara_float(cdfid,vhfileid,start,count,tmp3dp[0][0])))
	//       ERR(status);
	//    start[0] = inxt;
	//    if ((status = nc_get_vara_float(cdfid,vhfileid,start,count,tmp3dn[0][0])))
	//         ERR(status);

	for (k=0;k<NZ;k++)
		for (i=0;i<NXTOT;i++)
			for (j=0;j<NYTOT;j++)
				vhtm[k][i+2][j+2]= tmp3d[k][j][i];
	//vhtm[k][i+2][j+2]= fact1 * tmp3d[k][j][i] +
	//    fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i];

	close_file(&cdfid,&file);

	// finally wd

	sprintf(infile,"WD.%s.nc",fieldtype);

	strcpy(inpath, readpath);
	strcat(inpath, infile);

//	printf("Looking for file '%s'.\n",inpath);

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
		strcat(inpath, ".cdf");
		err = open_input_file(inpath,&file,&cdfid,&timeid);
		if (err != 0) {
			printf("Unable to find WD file.\n");
			exit(-73);
		}
	}

	if ((status = nc_inq_varid(cdfid, "wd", &wdfileid)))
		ERR(status);

	//    if ((status = nc_inq_varid(cdfid, "WDCLIM", &wdclimid)))
	//if ((status = nc_inq_varid(cdfid, "WD", &wdclimid)))
	//   ERR(status);

	bzero(start, MAX_NC_VARS * sizeof(long));

	count[0] = 1;
	count[1] = NZ;
	count[2] = NYTOT;
	count[3] = NXTOT;

	//    printf("read_clim WD month=  %i\n",imon);

	count[1] = NZ+1;
	start[0] = imon;
	if ((status = nc_get_vara_float(cdfid,wdfileid,start,count,tmp3d[0][0])))
		ERR(status);
	//    start[0] = iprv;
	//    if ((status = nc_get_vara_float(cdfid,wdfileid,start,count,tmp3dp[0][0])))
	//       ERR(status);
	//    start[0] = inxt;
	//    if ((status = nc_get_vara_float(cdfid,wdclimid,start,count,tmp3dn[0][0])))
	//       ERR(status);

	for (k=0;k<NZ+1;k++)
		for (i=0;i<NXTOT;i++)
			for (j=0;j<NYTOT;j++)
				wd[k][i+2][j+2]= timekeeper.dt * tmp3d[k][j][i];
	//wd[k][i+2][j+2]= dt * (fact1 * tmp3d[k][j][i] +
	//    fact0 * tmp3dp[k][j][i] + fact2 * tmp3dn[k][j][i]);

	//BX-a
	// take out PME which is read in onto wd[0][i][j]
	//for (i=0;i<NXTOT;i++)
	//	for (j=0;j<NYTOT;j++)
	//	    wd[k][i+2][j+2]= 0;
	//    printf("ATTENTION: EmP is set to 0. for this run.\n");
	//BX-e


	free3d_f(tmp3d, NZ+1);


	close_file(&cdfid,&file);
}

void read_h(int imon, char *fieldtype, char *readpath, double ***hread)
{
	
	int i,j,k;
	int ii;
	int err, cdfid, timeid;
	char infile[100], inpath[1000];
	FILE *file;
	int status;

	int hfileid;

	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];

	float*** tmp3d;
//	extern double ***hread;

	tmp3d  = alloc3d_f(NZ,NYTOT,NXTOT);
//	printf("Reading hread from month: %d\n",imon);
	sprintf(infile,"H.%s.nc",fieldtype);

	strcpy(inpath, readpath);
	strcat(inpath, infile);

	//  printf("Looking for file '%s'.\n",inpath);

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
		strcat(inpath, ".cdf");
		err = open_input_file(inpath,&file,&cdfid,&timeid);
		if (err != 0) {
			printf("Unable to find H oceanclim file.\n");
			exit(-73);
		}
	}


	//    printf("read_hclim H month=  %i\n",imon);



	if ((status = nc_inq_varid(cdfid, "h", &hfileid)))
		ERR(status);

	//   if ((status = nc_inq_varid(cdfid, "HCLIM", &hclimid)))
	//if ((status = nc_inq_varid(cdfid, "H", &hclimid)))
	//   ERR(status);

	bzero(start, MAX_NC_VARS * sizeof(long));

	count[0] = 1;
	count[1] = NZ;
	count[2] = NYTOT;
	count[3] = NXTOT;

	start[0] = imon;
	// allocate NZ so that it can be used for WDCLIM as well
	if ((status = nc_get_vara_float(cdfid,hfileid,start,count,tmp3d[0][0])))
		ERR(status);
	for (k=0;k<NZ;k++) {
		for (i=0;i<NXTOT;i++) {
			for (j=0;j<NYTOT;j++) {
				hread[k][i+2][j+2]= tmp3d[k][j][i];
			}
		}
	}
	//printf("Read in h=%g, hstart=%g, hread=%g\n",h[1][100][100],hstart[1][100][100],hread[1][100][100]);

//	free3d_f(tmp3d, NZ);


	//	zonal re-entrance
	for (k=0;k<NZ;k++) {
		for (j=0;j<=NYMEM-1;j++) {
			hread[k][nx+1][j] = hread[k][2][j];
			hread[k][nx+2][j] = hread[k][3][j];
			hread[k][0][j] =   hread[k][nx-1][j];
			hread[k][1][j] =   hread[k][nx][j];
		}
	}
#ifdef REENTRANT_Y
	//      meridional re-entrance
	for (i=2;i<=nx;i++) {
		ii = 363 - i;
		for (k=0;k<NZ;k++) {
			hread[k][ii][ny+1] = hread[k][i][ny];
			hread[k][ii][ny+2]   = hread[k][i][ny-1];
		}
	}
#endif

	close_file(&cdfid,&file);

}



/* ashao: Read data routines (UH, VH, WD, T, S) for hindcast runs */
/*
void read_var3d( char *inpath, char *varname, int imon, double ***data)
{

	int i,j,k,ii;
	int err, cdfid, timeid, varid;
	char infile[25];
	FILE *file;
	int status;
	int inxt,iprv;
	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];
	double ***tmp3d;

	start[0] = imon;
	start[1] = 0;
	start[2] = 0;
	start[3] = 0;

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if ((status = nc_inq_varid(cdfid, varname, &varid)))
		bzero(start, MAX_NC_VARS * sizeof(long));


	count[0] = 1;
	count[1] = NZ;
	count[2] = NYTOT;
	count[3] = NXTOT;

	//    for (i=0;i<4;i++)
	//	printf("start[%d]: %d,count[%d]: %d\n",i,start[i],i,count[i]);

	tmp3d  = alloc3d(NZ,NYTOT,NXTOT);
	if ((status = nc_get_vara_double(cdfid,varid,start,count,tmp3d[0][0])))
		ERR(status);
	for (k=0;k<count[1];k++)
		for (i=0;i<NXTOT;i++)
			for (j=0;j<NYTOT;j++)
				data[k][i+2][j+2]= tmp3d[k][j][i];


	for (k=0;k<NZ;k++) {
		for (j=0;j<=NYMEM-1;j++) {
			data[k][nx+1][j] = data[k][2][j];
			data[k][nx+2][j] = data[k][3][j];
			data[k][0][j] =   data[k][nx-1][j];
			data[k][1][j] =   data[k][nx][j];
		}
	}
#ifdef REENTRANT_Y
	//      meridional re-entrance
	for (i=2;i<=nx;i++) {
		ii = 363 - i;
		for (k=0;k<NZ;k++) {
			data[k][ii][ny+1] = data[k][i][ny];
			data[k][ii][ny+2]   = data[k][i][ny-1];
		}
	}
#endif
	free3d(tmp3d,NZ);

	close_file(&cdfid,&file);	
}
*/
void read_var2d( char *inpath, char *varname, double **data)
{

	int i,j;
	int err, cdfid, timeid, varid;
	char infile[100];
	FILE *file;
	int status;
	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];
	float **tmp2d;

	start[0] = 0;
	start[1] = 0;

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if ((status = nc_inq_varid(cdfid, varname, &varid)))
		bzero(start, MAX_NC_VARS * sizeof(long));

	count[0] = NYTOT;
	count[1] = NXTOT;

	tmp2d  = alloc2d_f(NYTOT,NXTOT);
	if ((status = nc_get_vara_float(cdfid,varid,start,count,tmp2d[0])))
		ERR(status);
		for (i=0;i<NXTOT;i++)
			for (j=0;j<NYTOT;j++)
				data[i+2][j+2]= tmp2d[j][i];

	wrap_reentrance_2d(data);

	free2d_f(tmp2d,NZ);
	close_file(&cdfid,&file);

}

void read_var3d(char *readpath, char *varname, int imon, double ***readarray)
{
	
	int i,j,k;
	int ii;
	int err, cdfid, timeid;
	char inpath[1000];
	FILE *file;
	int status;

	int fileid;

	size_t start[MAX_NC_VARS];
	size_t count[MAX_NC_VARS];

	double*** tmp3d;

	tmp3d  = alloc3d(NZ,NYTOT,NXTOT);
	
	strcpy(inpath,readpath);
	printf("Looking for file '%s'.\n",inpath);

	err = open_input_file(inpath,&file,&cdfid,&timeid);
	if (err != 0) {
		strcat(inpath, ".cdf");
		err = open_input_file(inpath,&file,&cdfid,&timeid);
		if (err != 0) {
			printf("Unable to find file %s.\n",inpath);
			exit(-73);
		}
	}

	if ((status = nc_inq_varid(cdfid, varname, &fileid)))
		ERR(status);

	bzero(start, MAX_NC_VARS * sizeof(long));

	count[0] = 1;
	count[1] = NZ;
	count[2] = NYTOT;
	count[3] = NXTOT;

	start[0] = imon;
	if ((status = nc_get_vara_double(cdfid,fileid,start,count,tmp3d[0][0])))
		ERR(status);
	for (k=0;k<NZ;k++) {
		for (i=0;i<NXTOT;i++) {
			for (j=0;j<NYTOT;j++) {
				readarray[k][i+2][j+2]= tmp3d[k][j][i];
			}
		}
	}

	free3d(tmp3d, NZ);

	//	zonal re-entrance
	for (k=0;k<NZ;k++) {
		for (j=0;j<=NYMEM-1;j++) {
			readarray[k][nx+1][j] = readarray[k][2][j];
			readarray[k][nx+2][j] = readarray[k][3][j];
			readarray[k][0][j] =   readarray[k][nx-1][j];
			readarray[k][1][j] =   readarray[k][nx][j];
		}
	}
	//      meridional re-entrance
	for (i=2;i<=nx;i++) {
		ii = 363 - i;
		for (k=0;k<NZ;k++) {
			readarray[k][ii][ny+1] = readarray[k][i][ny];
			readarray[k][ii][ny+2]   = readarray[k][i][ny-1];
		}
	}

	close_file(&cdfid,&file);

}
