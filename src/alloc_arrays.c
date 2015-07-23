#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "init.h"
#include "io.h"
#include "alloc.h"
#include "par_IO.h"

int alloc_error(char* arr_name)
{
	fprintf(stderr,"not enough memory for %s!\n", arr_name);
	quit(2);
	return 2;
}


int alloc_arrays() {
	extern double ***h, ***hstart, ***hend, ***htest;
	extern double ***mn_h;
	extern double ***uhtm, ***mn_uhtm;
	extern double ***vhtm, ***mn_vhtm;
	extern double ***wd, ***mn_wd;
	extern double ***Temptm, ***Salttm;

//	double **areagr;
//	double **D;
//	double **umask, **vmask;

	if (! (h = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("h");
	if (! (hstart = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("hstart");
	if (! (hend = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("hend");
	if (! (htest = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("hend");
	if (! (Temptm = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("Temptm");
	if (! (Salttm = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("Salttm");
	if (! (mn_h = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("mn_h");
	if (! (uhtm = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("uhtm");
	if (! (vhtm = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("vhtm");
	if (! (wd = alloc3d(NZ+1,NXMEM,NYMEM))) alloc_error("wd");
	if (! (mn_uhtm = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("mn_uhtm");
	if (! (mn_vhtm = alloc3d(NZ,NXMEM,NYMEM))) alloc_error("mn_vhtm");
	if (! (mn_wd = alloc3d(NZ+1,NXMEM,NYMEM))) alloc_error("mn_wd");


	return 0;
}

