//extern double * U0;

double str2num(const char * number);


/* This function examine whether a string
 * represents a real number.
 * Transform the string represents a
 * negtive number into a string represents
 * a positive one and return its' sign.
 * It returns 0 if the string do not
 * represents a real number.
 */
int format_string(char * str);


/* This function reads the initial data file
 * to generate the initial data.
 * It returns 0 if successfully read the file,
 * while returns the index of the wrong entry.
 */
int file_read(FILE * fp, double * U, int num);




/* this function counts how many numbers are there
 * the initial data file.
 */
int _1D_file_pre_read(FILE * fp);


/* This function reads the initial data file. The function 
 * initialize return a pointer pointing to the position of
 * a block of memory consisting (m+1) variables of type
 * double. The value of first of these variables is m.
 * The following m variables are the initial value.
 */
void _1D_initialize(char * name, char * addU, char * addP, char * addrho);


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
void _1D_file_write(int m, int N, double * RHO[], double * U[], double * P[], double * Ene[], double * X[], double * cpu_time, double * config, char * scheme);
