/*
 * gas_tracer_type.h
 *
 *  Created on: Apr 17, 2013
 *      Author: ashao
 *
 *  Defines a structure to be used for gas tracers with allocatable values for
 *  transient tracers
 */

// A variable type for gases with constant atmospheric concentration (e.g. O2,
// N2, Ar)
struct Gas_const_t {

        int num_sat_coeffs;
        int num_Sc_coeffs;
        							 // Number of coeffiicients used in the
        							 // Schmidt number, bunsen solubility, and solubility function
        							 // coefficient equations
        double atmconc; // Atmospheric concentration
        double *Sc_coeffs; // Coefficients for Schmidt number
        double *sat_coeffs;

};
