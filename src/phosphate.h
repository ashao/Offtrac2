/*
 * phosphate.h
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */

extern int mPHOSPHATE;
extern int mDOP;
// Output arrays
extern double ***mn_phos;
extern double ***mn_dop;
extern double ***mn_pobs;
extern double ***mn_jpo4;
extern double ***mn_jremin;
extern double ***mn_jremdop;
extern double ***mn_jprod;
//
// // Working arrays
extern double ***phosphate_init;
extern double ***jpo4;
extern double ***po4_star_lay;
extern double ***jprod;
extern double ***jremin;
extern double ***jremdop;
extern double ***jdop;
extern double flux_pop[NXMEM][NYMEM];
//

/* SUBROUTINE prototypes */
void allocate_phosphate( );
void initialize_phosphate(  );
void apply_phosphate_jterms();

