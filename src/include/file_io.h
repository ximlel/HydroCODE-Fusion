/**
 * @file file_io.h
 * @brief This file is the header file that controls data input and output.
 * @details This header file declares functions in the folder 'file_io'.
 */

#ifndef FILEIO_H
#define FILEIO_H

// io_control.c
void example_io(const char * example, char * add_mkdir, const int i_or_o);

int flu_var_count(FILE * fp, const char * add);
int flu_var_count_line(FILE * fp, const char * add, int * n_x);

int flu_var_read(FILE * fp, double * U, const int num);

// _1D_file_in.c
struct flu_var _1D_initialize(const char * name);
struct flu_var _2D_initialize(const char * name);

// _1D_file_out.c
void _1D_file_write(const int m, const int N, const struct cell_var_stru CV, 
                    double * X[], const double * cpu_time, const char * name, const double * time_plot);
void _2D_file_write(const int n_x, const int n_y, const int N, const struct cell_var_stru CV[],
		    double * X[], double * Y[], const double * cpu_time, const char * name, const double * time_plot);
void _2D_TEC_file_write(const int n_x, const int n_y, const int N, const struct cell_var_stru CV[],
			double * X[], double * Y[], const double * cpu_time, const char * problem, const double * time_plot);

// config_handle.c
void configurate(const char * name);

void config_write(const char * add_out, const double * cpu_time, const char * name);


#endif
