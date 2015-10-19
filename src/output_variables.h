/*
 * output_variables.h

 *
 *  Created on: Apr 29, 2015
 *      Author: ashao
 */

#include "init.h"

int map_variable_to_index( char *var_name );
void submit_for_averaging( double ***average_array, double ***variable_snapshot );
void submit_for_averaging_2d( double **average_array, double **variable_snapshot );
