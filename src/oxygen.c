#include <stdio.h>
#include <math.h>

#include "oxygen.h"
#include "biotic.h"
#include "init.h"
#include "phosphate.h"
#include "netcdf.h"
#include "alloc.h"
#include "init.h"
#include "util.h"
#include "oxygen.h"
#include "timekeeper.h"

/* OXYGEN VARIABLE DECLARATIONS */
// Auxiliary variables
int mOXYGEN;
// Output arrays
double ***mn_oxygen;
double ***mn_o2sat;
double ***mn_jo2;
double mn_oxyflux[NXMEM][NYMEM];
// Working arrays
double ***oxy_init;
double ***o2_sat;
double oxyflux[NXMEM][NYMEM];
double ***jo2;
extern double ****tr;
extern int oceanmask[NXMEM][NYMEM];

extern struct timekeeper_t timekeeper;

void allocate_oxygen (  ) {
	printf("Allocating oxygen arrays\n");

	int i, j, k;

	// Allocate working and output arrays
	mn_oxygen = alloc3d(NZ,NXMEM,NYMEM);
	jo2 = alloc3d(NZ,NXMEM,NYMEM);
	mn_jo2 = alloc3d(NZ,NXMEM,NYMEM);
	o2_sat = alloc3d(NZ,NXMEM,NYMEM);
	mn_o2sat = alloc3d(NZ,NXMEM,NYMEM);
	oxy_init = alloc3d(NZ,NXMEM,NYMEM);
}

void initialize_oxygen(int imon ) {

	extern double ***hend;
	extern double misval;
	char varname[100];

	printf("Oxygen index in main tracer array: %d\n",mOXYGEN);
	printf("Setting oxygen variable description...");
	strcpy(varname,"oxygen");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_oxygen");
	strcpy(vars[map_variable_to_index(varname)].longname,"Oxygen Concentration");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"mol/m^3");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;

	strcpy(varname,"o2sat");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_o2sat");
	strcpy(vars[map_variable_to_index(varname)].longname,"Oxygen saturation concentration");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"mol/m^3");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;

	strcpy(varname,"jo2");
	strcpy(vars[map_variable_to_index(varname)].name,"mn_jo2");
	strcpy(vars[map_variable_to_index(varname)].longname,"Oxygen source/sink term");
	vars[map_variable_to_index(varname)].hor_grid='h';
	vars[map_variable_to_index(varname)].z_grid='L';
	vars[map_variable_to_index(varname)].t_grid='s';
	strcpy(vars[map_variable_to_index(varname)].units,"mol/m^3/s");
	vars[map_variable_to_index(varname)].mem_size='d';
	vars[map_variable_to_index(varname)].mval=MISVAL;

	// Determine how to set the initial distribution of oxygen

	if (run_parameters.restart_flag) {
		printf("Initializing oxygen from restart %s\n",run_parameters.restartfile);
		read_var3d( run_parameters.restartfile, "mn_oxygen", 0, oxy_init);
	}
	else {
		printf("Initializing oxygen from WOA09\n");
		read_woa_file(imon, hend, oxy_init, "woa09.o2.nc", "LEVO2",1e-3);
	}

	apply_mask(oxy_init,oceanmask);
	// Copy the initialized tracer value over to main trace array
	copy_darray3d(tr[mOXYGEN],oxy_init);

	free3d(oxy_init,NZ);

}

void oxygen_saturation(double ***T, double ***S)
{
	/* compute oxygen saturation concentrations in the mixed layer in ml/l
	 * using Garcia and Gordon 1992 (L&O,1307) and convert to mol/m^3.
	 * Based on Matlab code by Allison Shaw and Steve Emerson.
	 * Ivan Lima (Nov 2002) */

	int i, j, k;
	double TS[NZ][NXMEM][NYMEM], o2s_ml[NZ][NXMEM][NYMEM];
	double work[NZ][NXMEM][NYMEM]; /* temporary work array */
	double molvolo2 = 22.3916; /* molar volumes (liters/mole) */
	double A0, A1, A2, A3, A4, A5, B0, B1, B2, B3, C0; /* constants */

	A0 = 2.00907;
	A1 = 3.22014;
	A2 = 4.05010;
	A3 = 4.94457;
	A4 = -.256847;
	A5 = 3.88767;
	B0 = -.00624523;
	B1 = -.00737614;
	B2 = -.0103410;
	B3 = -.00817083;
	C0 = -4.88682E-07;

	for (k = 0; k < NZ; k++) {
		for (i = X1; i < NXMEM; i++) {
			for (j = Y1; j < NYMEM; j++) {
				// HF this was set to 298.15, but T is in C, not K here!
				// HF the Garcia and Gordon paper uses t=40C as upper bound
				if (T[k][i][j] > 40.0) {
					printf(
							"WARNING: TEMPERATURE OUT OF RANGE in oxygen saturation: %g \n",
							T[k][i][j]);
					printf(
							"  setting Temp to 40 deg C at array position %i,%i,%i \n",
							i, j, k);
					T[k][i][j] = 40.0;
				}
				TS[k][i][j]
				         = log((298.15 - T[k][i][j]) / (273.15 + T[k][i][j]));
				work[k][i][j] = A0 + A1 * TS[k][i][j] + A2 * pow(TS[k][i][j],
						2.) + A3 * pow(TS[k][i][j], 3.) + A4 * pow(TS[k][i][j],
								4.) + A5 * pow(TS[k][i][j], 5.) + S[k][i][j] * (B0 + B1
										* TS[k][i][j] + B2 * pow(TS[k][i][j], 2.) + B3 * pow(
												TS[k][i][j], 3.)) + C0 * pow(S[k][i][j], 2.);
				o2s_ml[k][i][j] = exp(work[k][i][j]);
				o2_sat[k][i][j] = o2s_ml[k][i][j] / molvolo2;
			}
		}
	}
}

void apply_oxygen_jterms( ) {
	int i,j,k;

	extern double hend[NZ][NXMEM][NYMEM];
	extern double misval;
	// j terms here are calculated from biotic_sms routine in biotic.c
	printf("Example jo2/o2: %e/%f\n",jo2[10][127][127],tr[mOXYGEN][10][127][127]);
	for (i = 0; i < NXMEM; i++) {
		for (j = 0; j <NYMEM; j++) {
			if (oceanmask[i][j]) {
				for (k = 0; k < NZ; k++) {
					if (hend[k][i][j]>EPSILON)  tr[mOXYGEN][k][i][j] += timekeeper.dt * jo2[k][i][j];
				}
			} else {
				for (k = 0; k < NZ; k++) {
					tr[mOXYGEN][k][i][j] = misval;
				}
			}
		}
	}


}

void step_oxygen( ) {
	int i,j,k;
	extern double ****tr;

	// Set oxygen values to saturation at the mixed layer to mimic equilibrium with the atmosphere
	extern double ***Temptm;
	extern double ***Salttm;
	const double Sc_coeffs[4] = {1953.4, 128.00, 3.9918, 0.050091};

	oxygen_saturation(Temptm, Salttm, o2_sat);
	apply_oxygen_jterms();

	if (run_parameters.do_gasex)
		gas_exchange(mOXYGEN,Sc_coeffs,o2_sat);
	else
		for (k=0;k<NML;k++)
			copy_darray2d(tr[mOXYGEN][k],o2_sat[0]);
}

