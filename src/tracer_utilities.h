/*
 * tracer_utilities.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */


void z_depth(double ***h, double ***depth);
void allocate_ts( );
void read_temp_and_salt( int imon, char *fieldtype, char *path);
double calc_inventory( double ***tr_array, double ***h_array );
double linear_interpolation(const double xin[], const double yin[], double xi, int numin); 
