/**
 * @file file_io.h
 * @brief This file is the header file that controls data input and output.
 * @details This header file declares functions in files '_1D_f_io.c' and 'common.c'.
 */

extern double * U0;   //!< Initial velocity data array pointer.
extern double * P0;   //!< Initial pressure data array pointer.
extern double * RHO0; //!< Initial density  data array pointer.

//! Define the number of configuration parameters.
#ifndef N_CONF
#define N_CONF 7
#endif

int file_read(FILE * fp, double * U, const int num);

void example_io(const char * name, char * add_mkdir, const int i_or_o);

void _1D_initialize(const char * name, const char * add_in);

void _1D_configurate(double * config, const char * name, const char * add);

void _1D_file_write(const int m, const int N, 
                    double * RHO[], double * U[], double * P[], double * Ene[], double * X[], 
                    const double * cpu_time, const double * config, const char * name, const char * add_out);
