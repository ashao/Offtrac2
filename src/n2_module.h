extern int mN2;
// Output arrays
extern double ***mn_n2;
extern double ***mn_n2sol;
extern double ***n2sol;
// Working arrays
extern double ***n2_init;

/* SUBROUTINE PROTOTYPES */
void allocate_n2(  );
void initialize_n2(  );
void step_n2(  );
void calc_n2_saturation();
void initialize_n2_properties();

