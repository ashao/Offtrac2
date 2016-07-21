/*
 * oxygen.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */

#include "init.h"
extern int mAGE1, mAGE2;
extern double age_pair_decay; // Equivalent to lambda
extern double age_pair_conversion; // Conversion from seconds to units of age pair
// Output arrays
extern double ***mn_age_pair1, ***mn_age_pair2;
// Working arrays
extern double ***age_pair1_init;
extern double ***age_pair2_init;
extern double *age1_inventory;
extern double *age2_inventory;


/* SUBROUTINE PROTOTYPES */
void allocate_age_pair(  );
void initialize_age_pair(  );
void step_age_pair( double dt );



