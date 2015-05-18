/*
 * tracer_utilities.c
 *
 *  Created on: Apr 15, 2014
 *      Author: ashao
 */
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "init.h"
#include "tracer_utilities.h"
#include "alloc.h"
#include "read.h"
extern int oceanmask[NXMEM][NYMEM];
extern double D[NXMEM][NYMEM];
extern double ***Temptm,***Salttm;

void allocate_ts( ) {

	Temptm = alloc3d(NZ,NXMEM,NYMEM);
	Salttm = alloc3d(NZ,NXMEM,NYMEM);

}

void read_temp_and_salt( int imon, char *fieldtype, char* path) {
	char filename[50];
	char saltpath[300];
	char temppath[300];

	strcpy(temppath,path);
	sprintf(filename,"TS.%s.nc",fieldtype);
	strcat(temppath,filename);

	printf("Reading temperature and salinity from interval %d\n",imon);
	read_var3d( temppath, "temp", imon, Temptm);
	read_var3d( temppath, "salt", imon, Salttm);

}

void z_depth(double h[NZ][NXMEM][NYMEM], double depth[NZ][NXMEM][NYMEM]) {
	/* compute depth in meters for use in tracer depth dependent functions
	 * Ivan Lima - Nov 2002 */
	int i, j, k;
	double hsum;
	/* initialize variables with zero */

	for (k = 0; k < NZ; k++) {
		for (i = X1; i <= nx; i++) {
			for (j = Y1; j <= ny; j++) {
				depth[k][i][j] = 0.0;
			}
		}
	}
	/* compute depth at layer k as the the depth of the point in the
	 * middle of that layer */
	for (i = X1; i <NXMEM; i++) {
		for (j = Y1; j <NYMEM; j++) {
			//BX - reinstated by HF
			if (oceanmask[i][j]) {
				hsum = h[0][i][j];
				depth[0][i][j] = h[0][i][j] * 0.5;
				for (k = 1; k < NZ; k++) {
					depth[k][i][j] = hsum + h[k][i][j] * 0.5;
					hsum += h[k][i][j];
				}
				//BX - reinstated by HF
			} else {
				for (k = 0; k < NZ; k++) {
					depth[k][i][j] = 0.0;
				}
			}
		}
	}


}

double calc_inventory( double ***array ) {

	extern double areagr[NXMEM][NYMEM];
	extern double ***hend;

	double inventory;
	int i, j, k;

	inventory = 0.0;

	for(i=0;i<NXMEM;i++)
		for(j=0;j<NYMEM;j++) {
			if(D[i][j]>MINIMUM_DEPTH)
				for(k=0;k<NZ;k++) {

					inventory += hend[k][i][j]*areagr[i][j]*array[k][i][j];

				}

		}
	return inventory;
}

double linear_interpolation(const double xin[], const double yin[], double xi, int numin) {

        int i,j,flipidx;
        int intpidx1,intpidx2;
        int decreasing;
        double deltax,dist;
        double x[numin],y[numin];
        double y0,y1,x0,x1,yi;

        // first check to see if xin increases or decreases
        
        if (xin[numin-1]>xin[0]) decreasing = 0;
        if (xin[numin-1]<xin[0]) decreasing = 1; 
        // flip the input vectors if it is decreasing
                if (decreasing) {
                for (i=0;i<numin;i++) {
                        flipidx = numin-1-i;
                        x[i] = xin[flipidx];
                        y[i] = yin[flipidx];
                }
        }
        if (decreasing==0) {
                for (i=0;i<numin;i++) {
                        x[i] = xin[i];
                        y[i] = yin[i];
                }
        }

	deltax = fabs(xi-x[0]);
        intpidx1 = 0;
        intpidx2 = 1;
        for (i=0;i<numin;i++) {

                dist = fabs(xi-x[i]); // Calculate how far away the current x value is from the desired number
                if (dist<deltax) {
                        deltax = dist;
                        intpidx1 = i; // We now know at least one bound of the interpolation interval
                        if ( xi<x[i] ){
                                intpidx2 = i - 1;
                        }
                        else{
                                intpidx2 = i + 1;
                        }
                }
        }
        y0 = y[ intpidx1 ];
        y1 = y[ intpidx2 ];
        x0 = x[ intpidx1 ];
        x1 = x[ intpidx2 ];
        yi = y0 + (y1-y0)*(xi-x0)/(x1-x0);
//        if (x0>1900) printf("x0,x1: %f,%f\ny0,y1: %f,%f\nxi,yi: %f/%f\n",x0,x1,y0,y1,xi,yi);
        return(yi);
}

