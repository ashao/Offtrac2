/*
 * tracer_utilities.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */


void z_depth(double h[NZ][NXMEM][NYMEM], double depth[NZ][NXMEM][NYMEM]);
void allocate_ts( );
void read_temp_and_salt( int imon, char *fieldtype);
double calc_inventory( int tridx );
