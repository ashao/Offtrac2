/* Header file with function definitions for util.c 
   H. Frenzel, IGPP, UCLA
   December 5, 2007 */

void set_darray2d_zero(double **arr, int NX, int NY);

void set_darray3d_zero(double ***arr, int nz, int NX, int NY);

void set_fix_darray3d_zero(double arr[][NXMEM][NYMEM], int nz);

void set_fix_darray2d_zero(double arr[NXMEM][NYMEM]);

void add_fix_darray2d(double arr1[NXMEM][NYMEM], double arr2[NXMEM][NYMEM]);

void mult_fix_darray2d(double arr[NXMEM][NYMEM], double factor);

void add_darray3d(double ***arr1, double ***arr2, int nz, int NX, int NY);

void mult_darray3d(double ***arr, int nz, int NX, int NY, double factor);

void add_darray2d(double **arr1, double **arr2, int NX, int NY);

void mult_darray2d(double **arr, int NX, int NY, double factor);

void divide_darray3d(double ***arr, double ***num, double ***den);
//BX-a
void mult_fix_darray2d_mv(double arr[NXMEM][NYMEM], double factor, double D[NXMEM][NYMEM], double mv);

void mult_darray3d_mv(double ***arr, int nz, int NX, int NY, double factor, double D[NXMEM][NYMEM], double mv);
//BX-;;e
/*
void copy_fix_darray2d(double arr1[NXMEM][NYMEM], double arr2[NXMEM][NYMEM]);
*/
void copy_darray3d(double ***arr1, double ***arr2, int nz, int NX, int NY);
void copy_darray2d(double **arr1, double **arr2, int NX, int NY);
void copy_darray2d_f(float **arr1, float **arr2, int NX, int NY);

//BX-a
void copy_fix_darray3d(double ***arr1, double arr2[NZ][NXMEM][NYMEM], int nz, int NX, int NY);
//BX-e


int mod(int a, int b);
// ashao
double linear_interp( double x0, double y0,  double x1, double y1, double xstar);
int calc_hindindex(int inmon, int nmonths);

void wrap_reentrance_2d( double **arr );
void wrap_reentrance_3d( double ***arr, int nz );

void apply_mask( double ***arr, int wetmask[NXMEM][NYMEM] );
