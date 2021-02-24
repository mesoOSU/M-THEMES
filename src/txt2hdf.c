#include <hdf5.h>
#include <stdio.h>
#include <string.h>

int main() {
	hid_t file_id, euler_space, grain_space, phase_space, euler_set, grain_set, phase_set;
	herr_t status;
	
	hsize_t euler_dims[4]; // dimensions for 4 dimensional array
	hsize_t phase_dims[3]; 
	hsize_t grain_dims[3]; 

	double euler_data[32][32][32][3];
	int grain_data[32][32][32];
	int phase_data[32][32][32];
	double phi;
	double theta;
	double omega;

	int grain;
	int phase;

	FILE *file = fopen("microstructure.txt", "r");
	char buffer[80];

	// read file into memory
	for (int i = 0; i < 32; i++) {
		for (int j = 0; j < 32; j++) {
			for (int k = 0; k < 32; k++) {
				fgets(buffer, 80, file);
				sscanf(buffer,"%lf %lf %lf %*d %*d %*d %d %d", &euler_data[i][j][k][0], &euler_data[i][j][k][1], &euler_data[i][j][k][2], &grain_data[i][j][k], &phase_data[i][j][k]); // read line from buffer 
			}
		}
	}

	// create hdf file
	file_id = H5Fcreate("microstructure.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// write euler angles
	euler_dims[0] = 32;
	euler_dims[1] = 32;
	euler_dims[2] = 32;
	euler_dims[3] = 3;

	euler_space = H5Screate_simple(4, euler_dims, NULL);

	euler_set = H5Dcreate2(file_id, "/euler_angles", H5T_IEEE_F32BE, euler_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(euler_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, euler_data);


	grain_dims[0] = 32;
	grain_dims[1] = 32;
	grain_dims[2] = 32;

	grain_space = H5Screate_simple(3, grain_dims, NULL); 
	grain_set = H5Dcreate2(file_id, "/grain", H5T_STD_I32BE, grain_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(grain_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, grain_data);

	
	phase_dims[0] = 32;
	phase_dims[1] = 32;
	phase_dims[2] = 32;

	phase_space = H5Screate_simple(3, phase_dims, NULL);

	phase_set = H5Dcreate2(file_id, "/phase", H5T_STD_I32BE, phase_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(phase_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, phase_data);


	status = H5Dclose(grain_set);
	status = H5Sclose(grain_space);
	status = H5Dclose(phase_set);
	status = H5Sclose(phase_space);
	status = H5Dclose(euler_set);
	status = H5Sclose(euler_space);
	status = H5Fclose(file_id);

	return 0;
}
