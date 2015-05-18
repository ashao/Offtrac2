/*
 * oxygen.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */

#include "init.h"
extern int mTTD;
// Output arrays
extern double ***mn_ttd;
// Working arrays
extern double ***ttd_init;
extern double *ttd_inventory;

/* SUBROUTINE PROTOTYPES */
void allocate_ttd(  );
void initialize_ttd(  );
void step_ttd(  );



