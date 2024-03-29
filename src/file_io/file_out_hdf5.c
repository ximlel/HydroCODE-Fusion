/**
 * @file file_out_hdf5.c
 * @brief This is a set of functions which control the readout of 1-D and 2-D data in HDF5 file format.
 * @attention  Library Dependency: HDF5®
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"
#ifdef HDF5PLOT
#include "hdf5.h"


/* Create the data itself in the dataset.
 * dataset_id = H5Dcreate(loc_id (location id), const char *name (dataset name),
 *                        hid_t dtype_id (data type) hid_t space_id (dataspace id),
 *                        link created_property, dataset created_property, dataset accessed_property);
 */
/* Writes data to the dataset.
 * herr_t write_status = H5Dwrite(dataset_id, memory data type,
 *                                memory_dataspace_id (define the memory dataspace and its selection),
 *                                 - H5S_ALL: the file dataspace is used as the memory dataspace, 
 *                                            the selection in file_dataspace_id is used as the memory dataspace selection.
 *                                file_dataspace_id (define the file dataspace selection),
 *                                 - H5S_ALL: all of the datasapce in the file,
 *                                   defined as all of the dimensional data defined by datasapce in the dataset.
 *                                conversion property of this I/O operation,
 *                                const void * buf (the location of data in memory) );
 */

/**
 * @brief Print out fluid variable 'v' with data array 'v_array'.
 */
#define PRINT_NC(v, v_array)						\
    do {								\
	dataset_id = H5Dcreate(group_id, #v, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); \
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v_array); \
	status = H5Dclose(dataset_id);					\
    } while (0)

/**
 * @brief This function write the 1-D solution into HDF5 output '.h5' files.
 * @param[in] m:   The number of spatial points in the output data.
 * @param[in] N:   The number of time steps in the output data.
 * @param[in] CV:  Structure of grid variable data in computational grid cells.
 * @param[in] X[]: Array of the coordinate data.
 * @param[in] cpu_time:  Array of the CPU time recording.
 * @param[in] problem:   Name of the numerical results for the test problem.
 * @param[in] time_plot: Array of the plotting time recording.
 */
void file_1D_write_HDF5(const int m, const int N, const struct cell_var_stru CV, 
			double * X[], const double * cpu_time, const char * problem, double time_plot[])
{
    char add_out[FILENAME_MAX+40];
    // Get the address of the output data folder of the test example.
    example_io(problem, add_out, 0);
    
    char file_data[FILENAME_MAX+40];
    strcpy(file_data, add_out);
    strcat(file_data, "FLU_VAR.h5");
    
    double *XX = (double*)malloc(m * sizeof(double));
    if(XX == NULL)
	{
	    printf("NOT enough memory! plot X\n");
	    exit(5);
	}

    hid_t file_id, group_id, attr_id, dataspace_id, dataspaceA_id, dataset_id;
    herr_t status;
    const unsigned rank = 1;
    const hsize_t dims[1] = {(hsize_t)m}, dimsA[1] = {1};
    /*
     * file_id = H5Fcreate(const char *filename,
     *                     unsigned overlay_flag,
     *                      - H5F_ACC_TRUNC->can overlay
     *                      - H5F_ACC_EXCL ->cannot overlay, error
     *                     hid_t created_property, hid_t accessed_property);
     */
    file_id = H5Fcreate(file_data, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    /* Create the dataspace information items in the metadata of the dataset.
     * dataspace_id = H5Screate_simple(int rank (spatial dimension),
     *                                 const hsize_t* current_dims (number of elements per dimension),
     *                                 const hsize_t* max_dims, (upper limit on the number of elements per dimension)
     *                                  - NULL: same as current_dim.
     *                                  - H5S_UNLIMITED: no upper limit, but the dataset must be chunked.);
     */
    dataspace_id  = H5Screate_simple(rank, dims,  NULL);
    dataspaceA_id = H5Screate_simple(rank, dimsA, NULL);

    char group_name[14];
    for(int k = 0; k < N; k++)
	{
	    sprintf(group_name, "/T%d", k);
	    /*
	     * new_group_id = H5Gcreate2(group_id, absolute or relative group link name,
	     *                           link created_propertygroup created_property, group accessed_property);
	     */
	    group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    // group2_id = H5Gcreate(group_id, "/MyGroup1/MyGroup2",  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    // group2_id = H5Gcreate(group_id, "./MyGroup2",  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	    /* Create an attribute. */
	    attr_id = H5Acreate(group_id, "time_plot", H5T_NATIVE_FLOAT, dataspaceA_id, H5P_DEFAULT, H5P_DEFAULT);
	    status  = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, time_plot+k);
	    status  = H5Aclose(attr_id);

	    PRINT_NC(RHO, CV.RHO[k]);
	    PRINT_NC(U,   CV.U[k]);
	    PRINT_NC(P,   CV.P[k]);
	    PRINT_NC(E,   CV.E[k]);
#ifdef RADIAL_BASICS
	    PRINT_NC(R, X[k]);
#else
	    for(int j = 0; j < m; ++j)
		XX[j] = 0.5 * (X[k][j] + X[k][j+1]);
	    PRINT_NC(X, XX);
#endif

	    status = H5Gclose(group_id);
	}

    status = H5Sclose(dataspaceA_id);
    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);
    
    free(XX);
    XX = NULL;
    if(status)
        return;
    /*
     * file_id = H5Fopen(const char *filename, 
     *                   unsigned read-write_flag,
     *                    - H5F_ACC_RDWR   read-write
     *                    - H5F_ACC_RDONLY read only
     *                   hid_t accessed_property);
     */
    // file_id = H5Fopen(file_data, H5F_ACC_RDWR, H5P_DEFAULT);
    /*
     * group_id = H5Gopen2(father group_id, absolute or relative group link name,
     *                     group accessed_property);
     */
    // group_id = H5Gopen(file_id, "MyGroup1/MyGroup2", H5P_DEFAULT);
}


/**
 * @brief This function write the 2-D solution into HDF5 output '.h5' files.
 * @param[in] n_x: The number of x-spatial points in the output data.
 * @param[in] n_y: The number of y-spatial points in the output data.
 * @param[in] N:   The number of time steps in the output data.
 * @param[in] CV:  Structure of variable data in computational grid cells.
 * @param[in] X:   Array of the x-coordinate data.
 * @param[in] Y:   Array of the y-coordinate data.
 * @param[in] cpu_time:  Array of the CPU time recording.
 * @param[in] problem:   Name of the numerical results for the test problem.
 * @param[in] time_plot: Array of the plotting time recording.
 */
void file_2D_write_HDF5(const int n_x, const int n_y, const int N, const struct cell_var_stru CV[],
			double ** X, double ** Y, const double * cpu_time, const char * problem, double time_plot[])
{
    char add_out[FILENAME_MAX+40];
    // Get the address of the output data folder of the test example.
    example_io(problem, add_out, 0);
    
    char file_data[FILENAME_MAX+40];
    strcpy(file_data, add_out);
    strcat(file_data, "FLU_VAR.h5");

    double ** XX, ** YY;
    XX = (double **)malloc(n_x * sizeof(double *));
    YY = (double **)malloc(n_x * sizeof(double *));
    if(XX == NULL || YY == NULL)
	{
	    printf("NOT enough memory! plot X or Y\n");
	    exit(5);
	}
    for(int j = 0; j < n_x; ++j)
	{
	    XX[j] = (double *)malloc(n_y * sizeof(double));
	    YY[j] = (double *)malloc(n_y * sizeof(double));
	    if(XX[j] == NULL || YY[j] == NULL)
		{
		    printf("NOT enough memory! plot X[%d] or Y[%d]\n", j, j);
		    exit(5);
		}
	}

    hid_t file_id, group_id, attr_id, dataspace_id, dataspaceA_id, dataset_id;
    herr_t status;
    const unsigned rank = 2;
    const hsize_t dims[2] = {(hsize_t)n_y, (hsize_t)n_x}, dimsA[1] = {1};

    file_id = H5Fcreate(file_data, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dataspace_id  = H5Screate_simple(rank, dims,  NULL);
    dataspaceA_id = H5Screate_simple(rank, dimsA, NULL);

    char group_name[14];
    for(int k = 0; k < N; k++)
	{
	    sprintf(group_name, "/T%d", k);
	    group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	    attr_id = H5Acreate(group_id, "time_plot", H5T_NATIVE_FLOAT, dataspaceA_id, H5P_DEFAULT, H5P_DEFAULT);
	    status  = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, time_plot+k);
	    status  = H5Aclose(attr_id);

	    PRINT_NC(RHO, CV[k].RHO);
	    PRINT_NC(U,   CV[k].U);
	    PRINT_NC(V,   CV[k].V);
	    PRINT_NC(P,   CV[k].P);
	    PRINT_NC(E,   CV[k].E);
	    for(int i = 0; i < n_y; ++i)
		for(int j = 0; j < n_x; ++j)
		    {			    
			XX[j][i] = 0.25*(X[j][i] + X[j][i+1] + X[j+1][i] + X[j+1][i+1]);
			YY[j][i] = 0.25*(Y[j][i] + Y[j][i+1] + Y[j+1][i] + Y[j+1][i+1]);
		    }
	    PRINT_NC(X, XX);
	    PRINT_NC(Y, YY);

	    status = H5Gclose(group_id);
	}

    status = H5Sclose(dataspaceA_id);
    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);

    for(int j = 0; j < n_x; ++j)
	{
	    free(XX[j]);
	    free(YY[j]);
	    XX[j] = NULL;
	    YY[j] = NULL;
	}
    free(XX);
    free(YY);
    XX = NULL;
    YY = NULL;
    if(status)
	return;
}
#endif
