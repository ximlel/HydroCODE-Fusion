/**
 * @file file_io.h
 * @brief This file is the header file that controls data input and output.
 * @details This header file declares functions in the folder 'file_io'.
 */

#ifndef FILEIO_H
#define FILEIO_H

// io_control.c
void example_io(const char * example, char * add_mkdir, const int i_or_o);

int flu_var_count     (FILE * fp, const char * add);
int flu_var_count_line(FILE * fp, const char * add, int * n_x);

int flu_var_read(FILE * fp, double * U, const int num);

void time_plot_read(const char * add_in, int * N_plot, double * time_plot[]);

// file_in.c
struct flu_var initialize_1D(const char * name, int * N_plot, double * time_plot[]);
struct flu_var initialize_2D(const char * name, int * N_plot, double * time_plot[]);

// file_out.c
void file_1D_write          (const int m,                  const int N, const struct cell_var_stru CV, 
                    double * X[], const double * cpu_time, const char * problem, const double time_plot[]);
void file_2D_write          (const int n_x, const int n_y, const int N, const struct cell_var_stru CV[],
		    double ** X, double ** Y, const double * cpu_time, const char * problem, const double time_plot[]);
void file_2D_write_POINT_TEC(const int n_x, const int n_y, const int N, const struct cell_var_stru CV[],
			double ** X, double ** Y, const double * cpu_time, const char * problem, const double time_plot[]);

// file_2D_unstruct_out.c
void file_write_2D_BLOCK_TEC(const struct flu_var FV, const struct mesh_var mv, const char * problem, const double time);
void file_write_3D_VTK      (const struct flu_var FV, const struct mesh_var mv, const char * problem, const double time);


// config_handle.c
void configurate(const char * name);

void config_write(const char * add_out, const double * cpu_time, const char * name);


#endif
