/*
 * output_variables.c
 *
 *  Created on: Apr 27, 2015
 *      Author: ashao
 */
#include <stdio.h>
#include <stdlib.h>
#include "init.h"
#include "initialize.h"
#include "timekeeper.h"
extern struct parameters run_parameters;
extern struct timekeeper_t timekeeper;
int map_variable_to_index( char *var_name ) {


	/* Maps a variable name to an index number */
	int var_idx=-99;
	// Time-invariant variables
	if (!strcmp(var_name,"depth"))		var_idx = 0;
	if (!strcmp(var_name,"geolat"))		var_idx = 1;
	if (!strcmp(var_name,"geolon"))		var_idx = 2;
	if (!strcmp(var_name,"wetmask"))	var_idx = 3;

	// Time-varying varibles
	if (!strcmp(var_name,"hlay"))		var_idx = 4;
	if (!strcmp(var_name,"uhtm"))		var_idx = 5;
	if (!strcmp(var_name,"vhtm"))		var_idx = 6;
	if (!strcmp(var_name,"wd"))			var_idx = 7;

	// User-defined variables
	if (!strcmp(var_name,"age"))		var_idx = 8;

	if (!strcmp(var_name,"cfc11"))	var_idx = 9;
	if (!strcmp(var_name,"pcfc11"))	var_idx = 10;
	if (!strcmp(var_name,"cfc12"))	var_idx = 11;
	if (!strcmp(var_name,"pcfc12"))	var_idx = 12;
	if (!strcmp(var_name,"sf6"))	var_idx = 13;
	if (!strcmp(var_name,"psf6"))	var_idx = 14;

#ifdef CONSERVATION_CHECK
	if (!strcmp(var_name,"test"))		var_idx = 8;
	if (!strcmp(var_name,"test_inventory"))	var_idx = 9;
	if (!strcmp(var_name,"htest"))		var_idx = 10;
#endif

#ifdef TTD
	if (!strcmp(var_name,"ttd"))		var_idx = 8;

#endif

	if (var_idx < 0)	{
		printf("ERROR: %s is undefined in output_variables.c\n",var_name);
		exit( -2 );
	}
	else	{
//		printf("%s index = %d\n",var_name,var_idx);
		return var_idx;

	}

}

void submit_for_averaging( double ***average_array, double ***variable_snapshot ) {

	int i, j, k;
	for (k=0;k<NZ;k++)
		for (i=0;i<NXMEM;i++)
			for (j=0;j<NYMEM;j++) {
				average_array[k][i][j] += variable_snapshot[k][i][j]*timekeeper.dt;
				if (timekeeper.averaging_counter == run_parameters.wrint)
					average_array[k][i][j] *= 1.0/timekeeper.accumulated_time_since_writing;
			}

}
