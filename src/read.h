void read_grid();
void read_var3d(char *readpath, char *varname, int imon, double ***readarray);
void read_var2d( char inpath[200], char varname[200], double **data);
void read_D();
void read_uvw(int imon, char *fieldtype, char *readpath);
void read_h(int imon, char *fieldtype, char* readpath, double ***hread);
