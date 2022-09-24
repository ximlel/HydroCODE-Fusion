#include <stdio.h>
#include <string.h>
#include "hdf5.h"

#include "../include/var_struc.h"
#include "../include/file_io.h"


void file_1D_write_HDF5(const int m, const int N, const struct cell_var_stru CV, 
                    double * X[], const double * cpu_time, const char * problem, const double time_plot[])
{
    char add_out[FILENAME_MAX+40];
    // Get the address of the output data folder of the test example.
    example_io(problem, add_out, 0);
    
    char file_data[FILENAME_MAX+40] = "";
    strcpy(file_data, add_out);
    strcat(file_data, "/FLU_VAR.h5");

    hid_t file_id;
    herr_t status;
    /*
     * file_id = H5Fcreate(const char *filename,
     *                     unsigned overlay_flag,
     *                      - H5F_ACC_TRUNC->can overlay
     *                      - H5F_ACC_EXCL ->cannot overlay, error
     *                     hid_t created_property, hid_t accessed_property);
     */
    file_id = H5Fcreate(file_data, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Fclose(file_id);

    hid_t group_id, group2_id;
    /*
     * new_group_id = H5Gcreate(group_id, absolute or relative group link name,
     *                          link created_property, group created_property, group accessed_property);
     */
    group_id = H5Gcreate(file_id, "/MyGroup1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    group2_id = H5Gcreate(group_id, "/MyGroup1/MyGroup2",  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // group2_id = H5Gcreate(group_id, "MyGroup2",  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // group2_id = H5Gcreate(group_id, "./MyGroup2",  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    unsigned rank = 2;
    hsize_t dims[2];
    dims[0] = rows;
    dims[1] = columns;
    hid_t dataspace_id;
    /* Create the dataspace information items in the metadata of the dataset.
     * dataspace_id = H5Screate_simple(int rank (spatial dimension),
     *                                 const hsize_t* current_dims (number of elements per dimension),
     *                                 const hsize_t* max_dims, (upper limit on the number of elements per dimension)
     *                                  - NULL: same as current_dim.
     *                                  - H5S_UNLIMITED: no upper limit, but the dataset must be chunked.)
     */
    dataspace_id = H5Screate_simple(rank, dims, NULL);

    hid_t dataset_id;
    /* Create the data itself in the dataset.
     * dataset_id = H5Dcreate(loc_id (location id), const char *name (dataset name),
     *                        hid_t dtype_id (data type) hid_t space_id (dataspace id),
     *                        link created_property, dataset created_property, dataset accessed_property)
     */
    dataset_id = H5Dcreate(file_id, "/dset", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Writes data to the dataset.
     * herr_t write_status = H5Dwrite(dataset_id, memory data type,
     *                                memory_dataspace_id (define the memory dataspace and its selection),
     *                                 - H5S_ALL: the file dataspace is used as the memory dataspace, 
     *                                            the selection in file_dataspace_id is used as the memory dataspace selection.
     *                                file_dataspace_id (define the file dataspace selection),
     *                                 - H5S_ALL: all of the datasapce in the file, defined as all of the dimensional data defined by datasapce in the dataset.
     *                                conversion property of this I/O operation,
     *                                const void * buf (the location of data in memory) )
     */
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix[0]);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    status = H5Gclose(group2_id);
    status = H5Gclose(group_id);
    status = H5Fclose(file_id);

    /*
     * file_id = H5Fopen(const char *filename, 
     *                   unsigned read-write_flag,
     *                    - H5F_ACC_RDWR   read-write
     *                    - H5F_ACC_RDONLY read only
     *                   hid_t accessed_property)
     */
    file_id = H5Fopen(file_data, H5F_ACC_RDWR, H5P_DEFAULT);
    /*
     * group_id = H5Gopen2(father group_id, absolute or relative group link name,
     *                     group accessed_property);
     */
    group_id = H5Gopen(file_id, "MyGroup1/MyGroup2", H5P_DEFAULT);
    status = H5Gclose(group_id);
    status = H5Fclose(file_id);
}


void file_2D_write_HDF5(const int n_x, const int n_y, const int N, const struct cell_var_stru CV[],
		    double ** X, double ** Y, const double * cpu_time, const char * problem, const double time_plot[])
{

}
