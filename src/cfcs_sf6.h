/*
 * cfcs_sf6.h
 *
 *  Created on: May 3, 2014
 *      Author: ashao
 */


extern double ***mn_cfc11;
extern double ***cfc11_init;
extern double ***cfc11_sol;
extern double ***pcfc11;
extern double ***mn_pcfc11;
extern int mCFC11;

extern double ***mn_cfc12;
extern double ***cfc12_init;
extern double ***cfc12_sol;
extern double ***pcfc12;
extern double ***mn_pcfc12;
extern int mCFC12;

extern double ***mn_sf6;
extern double ***sf6_init;
extern double ***sf6_sol;
extern double ***psf6;
extern double ***mn_psf6;
extern int mSF6;
#define NUMATMVALS 250
typedef struct {

//	int ntime;
	double time[NUMATMVALS];
	double nval[NUMATMVALS];
	double sval[NUMATMVALS];

} tracer_boundary;

#define NUMTRANSIENT 3
extern tracer_boundary atmconc[NUMTRANSIENT];

#define CFC11IDX 0
#define CFC12IDX 1
#define SF6IDX 2

void read_tracer_boundary( void );
void allocate_cfc11( void );
void initialize_cfc11 ( void );
void cfc11_saturation( void );
void cfc11_find_atmconc( void );
void surface_cfc11( void );

void allocate_cfc12( void );
void initialize_cfc12( void );
void cfc12_saturation( void );
void cfc12_find_atmconc( void );
void surface_cfc12( void );

void allocate_sf6( void );
void initialize_sf6( void );
void sf6_saturation( void );
void sf6_find_atmconc( void );
void surface_sf6( void );
