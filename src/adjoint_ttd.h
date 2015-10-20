/*
 * oxygen.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */

#include "init.h"
extern int mADJOINT;
// Output arrays
extern double **mn_adjttd;
// Working arrays
extern double ***adjttd_init;

/* SUBROUTINE PROTOTYPES */
void allocate_adjttd(  );
void initialize_adjttd(  );
void step_adjttd(  );



