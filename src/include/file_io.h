/**
 * @file file_io.h
 * @brief This file is the header file that controls data input and output.
 * @details This header file declares functions in files '_1D_f_io.c' and 'common.c'.
 */

#ifndef FILEIO_H
#define FILEIO_H

// io_control.c
void example_io(const char * name, char * add_mkdir, const int i_or_o);

int flu_var_count(FILE * fp, const char * add);


// _1D_file_in.c
void _1D_initialize(const char * name);


// _1D_file_out.c
void _1D_file_write(const int m, const int N, struct cell_var CV, double * X[], 
                    const double * cpu_time, const char * name);


// config_in.c
void _1D_configurate(const char * name);

#endif
