extern int mAR;
// Output arrays
extern double ***mn_ar;
extern double ***mn_arsol;
extern double ***arsol;
// Working arrays
extern double ***ar_init;

/* SUBROUTINE PROTOTYPES */
void allocate_ar(  );
void initialize_ar(  );
void step_ar(  );
void calc_ar_saturation();
void initialize_ar_properties();

