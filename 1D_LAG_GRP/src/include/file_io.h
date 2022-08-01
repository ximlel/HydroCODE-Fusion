int file_read(FILE * fp, double * U, const int num);

void example_io(const char * name, char * add_mkdir, const int i_or_o);

void _1D_initialize(const char * name, const char * add_in);

void _1D_configurate(double * config, const char * name, const char * add);

void _1D_file_write(const int m, const int N, 
                    double * RHO[], double * U[], double * P[], double * Ene[], double * X[], 
                    const double * cpu_time, const double * config, const char * name, const char * add_out);
