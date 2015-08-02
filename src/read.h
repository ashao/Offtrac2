void read_grid();
void read_var3d(char *readpath, char *varname, int imon, double ***readarray);
void read_var2d_time( char inpath[200], int imon, char varname[200], float **data);
void read_D();
void read_uvw(int imon, char *fieldtype, char *readpath);
void read_h(int imon, char *fieldtype, char* readpath, double ***hread);
void read_woa_file(int imon, double ***harray, double ***outarray, char *filename, char *varname, double conv_factor);
