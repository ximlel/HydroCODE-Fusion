//extern double * U0;


/* This function reads the initial data file
 * to generate the initial data.
 * It returns 0 if successfully read the file,
 * while returns the index of the wrong entry.
 */
int file_read(FILE * fp, double * U, int num);


void example_io(const char *name, char *add_mkdir, const int i_or_o);

/* This function reads the initial data file. The function 
 * initialize return a pointer pointing to the position of
 * a block of memory consisting (m+1) variables of type
 * double. The value of first of these variables is m.
 * The following m variables are the initial value.
 */
void _1D_initialize(char * name, char * add);


/* This function read the configuration data file,
 * and store the configuration data in the array
 * "config".
 * config[0] is the viscosity efficient
 * config[1] is the length of the time step
 * config[2] is the spatial grid size
 * config[3] is the largest value can be seen as zero
 * config[4] is the number of time steps
 */
void _1D_configurate(double * config, char * name, char * add);


/* This function write the solution into an output file.
 * It is quite simple so there will be no more comments.
 */
void _1D_file_write(int m, int N, double * RHO[], double * U[], double * P[], double * Ene[], double * X[], double * cpu_time, double * config, char * name, char * add_out);
