/********+*********+*********+*********+*********+*********+*********+*
 *         Initialize                                                 *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/


#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "initialize.h"
#include "init.h"
#include "offtrac.h"
#include "alloc.h"
#include "read.h"
#include "output_variables.h"
#include "io.h"

#include "ideal_age.h"
#include "ttd_bp.h"
#include "cfcs_sf6.h"
#include "n2_module.h"
#include "ar_module.h"
#include "oxygen.h"
#include "phosphate.h"
#include "adjoint_ttd.h"
#include "age_pair.h"

#include "conservation_check.h"

extern struct vardesc vars[NOVARS];
extern struct parameters run_parameters;
// end ashao

void initialize( void )
{
  /* Arguments: mon - Time of the start of the run in months.             */
  /* Arguments: inmon - Time of the restart of the run in months.         */

  int i, j, k, m;
  char varname[100];
  extern double ****tr;

  // Setup the variable descriptions

  strcpy(varname,"depth");
  strcpy(vars[map_variable_to_index(varname)].name,"Depth");
  strcpy(vars[map_variable_to_index(varname)].longname,"Depth to bottom");
  vars[map_variable_to_index(varname)].hor_grid='h';
  vars[map_variable_to_index(varname)].z_grid='1';
  vars[map_variable_to_index(varname)].t_grid='s';
  strcpy(vars[map_variable_to_index(varname)].units,"m");
  vars[map_variable_to_index(varname)].mem_size='d';
  vars[map_variable_to_index(varname)].mval=MISVAL;

  strcpy(varname,"geolat");
  strcpy(vars[map_variable_to_index(varname)].name,"geolat");
  strcpy(vars[map_variable_to_index(varname)].longname,"Latitude on tripolar grid");
  vars[map_variable_to_index(varname)].hor_grid='h';
  vars[map_variable_to_index(varname)].z_grid='1';
  vars[map_variable_to_index(varname)].t_grid='s';
  strcpy(vars[map_variable_to_index(varname)].units,"Degrees N");
  vars[map_variable_to_index(varname)].mem_size='d';
  vars[map_variable_to_index(varname)].mval=MISVAL;

  strcpy(varname,"geolon");
  strcpy(vars[map_variable_to_index(varname)].name,"geolon");
  strcpy(vars[map_variable_to_index(varname)].longname,"Longitude on tripolar grid");
  vars[map_variable_to_index(varname)].hor_grid='h';
  vars[map_variable_to_index(varname)].z_grid='1';
  vars[map_variable_to_index(varname)].t_grid='s';
  strcpy(vars[map_variable_to_index(varname)].units,"Degrees E");
  vars[map_variable_to_index(varname)].mem_size='d';
  vars[map_variable_to_index(varname)].mval=MISVAL;

  strcpy(varname,"wetmask");
  strcpy(vars[map_variable_to_index(varname)].name,"wetmask");
  strcpy(vars[map_variable_to_index(varname)].longname,"Wetmask");
  vars[map_variable_to_index(varname)].hor_grid='h';
  vars[map_variable_to_index(varname)].z_grid='1';
  vars[map_variable_to_index(varname)].t_grid='s';
  strcpy(vars[map_variable_to_index(varname)].units,"1 for ocean (0 for land)");
  vars[map_variable_to_index(varname)].mem_size='d';
  vars[map_variable_to_index(varname)].mval=MISVAL;

  strcpy(varname,"hlay");
  strcpy(vars[map_variable_to_index(varname)].name,"mn_h");
  strcpy(vars[map_variable_to_index(varname)].longname,"Isopycnal thickness");
  vars[map_variable_to_index(varname)].hor_grid='h';
  vars[map_variable_to_index(varname)].z_grid='L';
  vars[map_variable_to_index(varname)].t_grid='s';
  strcpy(vars[map_variable_to_index(varname)].units,"m");
  vars[map_variable_to_index(varname)].mem_size='d';
  vars[map_variable_to_index(varname)].mval=MISVAL;

  strcpy(varname,"uhtm");
  strcpy(vars[map_variable_to_index(varname)].name,"mn_uh");
  strcpy(vars[map_variable_to_index(varname)].longname,"Zonal mass transport");
  vars[map_variable_to_index(varname)].hor_grid='h';
  vars[map_variable_to_index(varname)].z_grid='L';
  vars[map_variable_to_index(varname)].t_grid='s';
  strcpy(vars[map_variable_to_index(varname)].units,"kg s^-1");
  vars[map_variable_to_index(varname)].mem_size='d';
  vars[map_variable_to_index(varname)].mval=MISVAL;

  strcpy(varname,"vhtm");
  strcpy(vars[map_variable_to_index(varname)].name,"mn_vh");
  strcpy(vars[map_variable_to_index(varname)].longname,"Meridional mass transport");
  vars[map_variable_to_index(varname)].hor_grid='h';
  vars[map_variable_to_index(varname)].z_grid='L';
  vars[map_variable_to_index(varname)].t_grid='s';
  strcpy(vars[map_variable_to_index(varname)].units,"kg s^-1");
  vars[map_variable_to_index(varname)].mem_size='d';
  vars[map_variable_to_index(varname)].mval=MISVAL;

  strcpy(varname,"wd");
  strcpy(vars[map_variable_to_index(varname)].name,"mn_wd");
  strcpy(vars[map_variable_to_index(varname)].longname,"Vertical velocity");
  vars[map_variable_to_index(varname)].hor_grid='h';
  vars[map_variable_to_index(varname)].z_grid='i';
  vars[map_variable_to_index(varname)].t_grid='s';
  strcpy(vars[map_variable_to_index(varname)].units,"kg s^-1");
  vars[map_variable_to_index(varname)].mem_size='d';
  vars[map_variable_to_index(varname)].mval=MISVAL;


  if (run_parameters.do_age) {
      allocate_age( );
      initialize_age( );
  }

  if (run_parameters.do_cfcs) {
      allocate_cfc11( );
      initialize_cfc11();
      allocate_cfc12( );
      initialize_cfc12();
      allocate_sf6( );
      initialize_sf6();
      read_tracer_boundary( );
  }

  if (run_parameters.conservation_check) {
      allocate_test();
      initialize_test();
  }

  if (run_parameters.do_ttd) {
      allocate_ttd();
      initialize_ttd();
  }

  if (run_parameters.do_n2) {
      allocate_n2();
      initialize_n2();
  }

  if (run_parameters.do_ar) {
      allocate_ar();
      initialize_ar();
  }

  if (run_parameters.do_oxygen) {
      allocate_oxygen( );
      initialize_oxygen( );

      allocate_phosphate( );
      initialize_phosphate();

  }

  if (run_parameters.do_adjttd) {
      allocate_adjttd();
      initialize_adjttd();
  }

  if (run_parameters.do_age_pair) {
      allocate_age_pair();
      initialize_age_pair();
  }


  /* zonal, meridional re-entrance    */
  for (m=0;m<run_parameters.tracer_counter;m++) {
      for (k=0;k<NZ;k++) {
	  for (j=0;j<=NYMEM-1;j++) {
	      tr[m][k][0][j] = tr[m][k][nx-1][j];
	      tr[m][k][1][j] = tr[m][k][nx][j];
	      tr[m][k][nx+1][j] = tr[m][k][2][j];
	      tr[m][k][nx+2][j] = tr[m][k][3][j];
	  }
      }
      for (k=0;k<NZ;k++) {
	  for (i=2;i<=nx;i++) {
	      tr[m][k][363-i][ny+1] = tr[m][k][i][ny];
	      tr[m][k][363-i][ny+2] = tr[m][k][i][ny-1];
	  }
      }
  }
}

void set_run_parameters( void )
{

  FILE *ptr_file;
  char attribute[100], value[100];
  char *line_read = NULL;
  ssize_t read;
  size_t len = 0;
  char fname[100];

  strcpy(fname,run_parameters.namelist_file);

  extern int flags[NOVARS], rflags[NOVARS];

  printf("reading parameters\n");
  ptr_file = fopen(fname,"r");
  if (ptr_file == NULL) {
      printf("Error: Cannot open %s\n",run_parameters.namelist_file);
      exit(1);
  }

  while( (read = getline(&line_read,&len,ptr_file)) != -1 ) {

      sscanf(line_read,"%s %s",attribute,value);

      // Set intervals of time integration
      //		printf("Setting %s to %s\n",attribute,value);
      if (!strcmp(attribute,"syear"))
	run_parameters.syear = atoi(value);
      if (!strcmp(attribute,"sinterval"))
	run_parameters.sinterval = atoi(value);
      if (!strcmp(attribute,"eyear"))
	run_parameters.eyear = atoi(value);
      if (!strcmp(attribute,"einterval"))
	run_parameters.einterval = atoi(value);

      // Define the type of run
      if (!strcmp(attribute,"use_hindcast"))
	run_parameters.use_hindcast = atoi(value);
      if (!strcmp(attribute,"restart_flag"))
	run_parameters.restart_flag = atoi(value);
      if (!strcmp(attribute,"timestep"))
	strcpy(run_parameters.timestep,value);
      if (!strcmp(attribute,"adjoint_integration"))
	run_parameters.adjoint_integration=atoi(value);

      // Set input and output
      if (!strcmp(attribute,"forcing_path"))
	strcpy(run_parameters.forcing_path,value);
      if (!strcmp(attribute,"normalyear_path"))
	strcpy(run_parameters.normalyear_path,value);
      if (!strcmp(attribute,"hindcast_path"))
	strcpy(run_parameters.hindcast_path,value);
      if (!strcmp(attribute,"outputfile")) {
	  sprintf(run_parameters.new_restartfile,"restart.%s.%d.%02d.nc",
		  value, run_parameters.eyear, run_parameters.einterval);
	  sprintf(run_parameters.outputfile,"%s.%d.%02d.nc",value,run_parameters.eyear,run_parameters.einterval);
      }
      if (!strcmp(attribute,"restartfile"))
	strcpy(run_parameters.restartfile,value);
      if (!strcmp(attribute,"wrint"))
	run_parameters.wrint = atoi(value);
      if (!strcmp(attribute,"do_averaging"))
	run_parameters.do_averaging = atoi(value);

      // Set forcing characteristics
      if (!strcmp(attribute,"ntime_climatology"))
	run_parameters.ntime_climatology = atoi(value);
      if (!strcmp(attribute,"ntime_hindcast"))
	run_parameters.ntime_hindcast = atoi(value);


      // Specify whether to save various variables
      if (!strcmp(attribute,"D"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"geolat"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"geolon"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"wetmask"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"hlay"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"uhtm"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"vhtm"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"wd"))
	flags[map_variable_to_index(attribute)] = atoi(value);

      if (!strcmp(attribute,"do_age"))
	run_parameters.do_age = atoi(value);

      if (!strcmp(attribute,"age"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"age_restart"))
	rflags[map_variable_to_index("age")] = atoi(value);

      if (!strcmp(attribute,"do_gasex"))
	run_parameters.do_gasex = atoi(value);
      if (!strcmp(attribute,"do_cfcs")) {
	  run_parameters.do_cfcs = atoi(value);
      }
      if (!strcmp(attribute,"cfc11"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"cfc11_restart"))
	rflags[map_variable_to_index("cfc11")] = atoi(value);
      if (!strcmp(attribute,"pcfc11"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"cfc12"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"cfc12_restart"))
	rflags[map_variable_to_index("cfc12")] = atoi(value);
      if (!strcmp(attribute,"pcfc12"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"sf6"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"sf6_restart"))
	rflags[map_variable_to_index("cfc11")] = atoi(value);
      if (!strcmp(attribute,"psf6"))
	flags[map_variable_to_index(attribute)] = atoi(value);


      if (!strcmp(attribute,"conservation_check")) {
	  run_parameters.conservation_check = atoi(value);
	  flags[map_variable_to_index("test")] = 1;
	  flags[map_variable_to_index("htest")] = 1;
      }

      if (!strcmp(attribute,"do_ttd")) {
	  run_parameters.do_ttd = atoi(value);
      }

      if (!strcmp(attribute,"do_adjttd")) {
	  run_parameters.do_adjttd = atoi(value);
	  flags[map_variable_to_index("adjttd")]=1;
	  rflags[map_variable_to_index("adjttd_restart")]=1;
      }
      if (!strcmp(attribute,"adjoint_initfile"))
	strcpy(run_parameters.adjoint_initfile,value);

      if (!strcmp(attribute,"ttd"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"ttd_restart"))
	rflags[map_variable_to_index("ttd")] = atoi(value);
      if (!strcmp(attribute,"num_ttd_intervals"))
	run_parameters.num_ttd_intervals = atoi(value);

      if (!strcmp(attribute,"do_n2"))
	run_parameters.do_n2 = atoi(value);
      if (!strcmp(attribute,"n2"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"n2sol"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"n2_restart"))
	rflags[map_variable_to_index("n2")] = atoi(value);

      if (!strcmp(attribute,"do_ar"))
	run_parameters.do_ar = atoi(value);
      if (!strcmp(attribute,"ar"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"arsol"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"ar_restart"))
	rflags[map_variable_to_index("ar")] = atoi(value);

      // Oxygen
      if (!strcmp(attribute,"do_oxygen"))
	run_parameters.do_oxygen = atoi(value);
      if (!strcmp(attribute,"woa_path"))
	strcpy(run_parameters.woa_path,value);
      if (!strcmp(attribute,"oxygen"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"o2sat"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"jo2"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"po4"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"dop"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"jpo4"))
	flags[map_variable_to_index(attribute)] = atoi(value);
      if (!strcmp(attribute,"oxygen_restart"))
	rflags[map_variable_to_index("oxygen")] = atoi(value);
      if (!strcmp(attribute,"po4_restart"))
	rflags[map_variable_to_index("po4")] = atoi(value);
      if (!strcmp(attribute,"dop_restart"))
	rflags[map_variable_to_index("dop")] = atoi(value);

      if(!strcmp(attribute,"do_age_pair"))
	run_parameters.do_age_pair = atoi(value);
      if (!strcmp(attribute,"age_pair_restart")) {
	rflags[map_variable_to_index("age_pair1")] = atoi(value);
	rflags[map_variable_to_index("age_pair2")] = atoi(value);
      }
      if (!strcmp(attribute,"age_pair1"))
	flags[map_variable_to_index("age_pair1")] = atoi(value);
      if (!strcmp(attribute,"age_pair2"))
	flags[map_variable_to_index("age_pair2")] = atoi(value);


      if(!strcmp(attribute,"age_pair_decay")) age_pair_decay = atof(value);
      if(!strcmp(attribute,"age_pair_conversion")) age_pair_conversion = atof(value);

  }
  free(line_read);
  fclose(ptr_file);

  if (run_parameters.conservation_check)
    mTEST = run_parameters.tracer_counter++;

  if (run_parameters.do_age)
    mAGE = run_parameters.tracer_counter++;

  if (run_parameters.do_ttd)
    mTTD = run_parameters.tracer_counter++;

  if (run_parameters.do_adjttd)
    mADJOINT = run_parameters.tracer_counter++;

  if (run_parameters.do_cfcs) {
      mCFC11 = run_parameters.tracer_counter++;
      mCFC12 = run_parameters.tracer_counter++;
      mSF6 = run_parameters.tracer_counter++;
  }
  if (run_parameters.do_n2)
    mN2 = run_parameters.tracer_counter++;

  if (run_parameters.do_ar)
    mAR = run_parameters.tracer_counter++;

  if (run_parameters.do_oxygen) {
      mOXYGEN = run_parameters.tracer_counter++;
      mPHOSPHATE = run_parameters.tracer_counter++;
      mDOP = run_parameters.tracer_counter++;
  }

  if (run_parameters.do_age_pair) {
      mAGE1 = run_parameters.tracer_counter++;
      mAGE2 = run_parameters.tracer_counter++;
  }


}
