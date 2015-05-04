#include <stddef.h>
#include <stdio.h>
#include "init.h"

/* Subroutines in offtrac.c */
void alloc_fields(void);
void set_metrics(void);


/* Subroutines in par_IO.c */
void collect_double_fields(double *field_in, int offx, int offy);
void spread_double_fields(double field_out[NXMEM][NYMEM]);
void spread_double_vector(double vec[], int numpts);
size_t write_layer(double field_in[], FILE *file, int float_vals);
size_t read_layer(double field_out[NXMEM][NYMEM], FILE *file, int float_vals);

/* Subroutine in set_metrics.c */
void set_metrics(void);

/* Subroutine in masks.c */
void initializemasks(void);

/* Subroutines in off_par.c */
void set_up_parallel(void);
void pass_vars(double var[][NXMEM][NYMEM], int nz, int pts_out, int nrows,
               double bounds, int dosync, int safemem, int dirflag);
void pass_var_2d(double var[NXMEM][NYMEM], int nz, int pts_out, int nrows,
                 double bounds, int dosync, int safemem, int dirflag);  
void spread_string(char name[], int length);
void sum_fields(double vec[], int numsum);
void sum_int_fields(int vec[], int numsum);
void check_field(double field_in[NXMEM][NYMEM], char write_string[],
                 int i0, int i1, int j0, int j1);
void quit(int status);

#define TO_NORTH 1
#define TO_SOUTH 2
#define TO_EAST 4
#define TO_WEST 8
