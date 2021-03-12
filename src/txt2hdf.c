#include <hdf5.h>
#include <stdio.h>
#include <string.h>

int main() {
	hid_t file_id, ph_space, th_space, om_space, grain_space, phase_space, ph_set, th_set, om_set, grain_set, phase_set, group_id;
	herr_t status;
	
	hsize_t ph_dims[3]; 
	hsize_t th_dims[3]; 
	hsize_t om_dims[3]; 
	hsize_t phase_dims[3]; 
	hsize_t grain_dims[3]; 

	double ph_data[32][32][32];
	double th_data[32][32][32];
	double om_data[32][32][32];
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
				sscanf(buffer,"%lf %lf %lf %*d %*d %*d %d %d", &ph_data[i][j][k], &th_data[i][j][k], &om_data[i][j][k], &grain_data[i][j][k], &phase_data[i][j][k]); // read line from buffer 
				printf("%lf, %lf, %lf, %d, %d\n", ph_data[i][j][k], th_data[i][j][k], om_data[i][j][k], grain_data[i][j][k], phase_data[i][j][k]);
			}
		}
	}

	// create hdf file
	file_id = H5Fcreate("microstructure.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// create euler group
	group_id = H5Gcreate(file_id, "/euler", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// phi
	ph_dims[0] = 32;
	ph_dims[1] = 32;
	ph_dims[2] = 32;

	ph_space = H5Screate_simple(3, ph_dims, NULL);

	ph_set = H5Dcreate2(file_id, "/euler/phi", H5T_IEEE_F32BE, ph_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(ph_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ph_data);

	// theta
	th_dims[0] = 32;
	th_dims[1] = 32;
	th_dims[2] = 32;

	th_space = H5Screate_simple(3, th_dims, NULL);

	th_set = H5Dcreate2(file_id, "/euler/theta", H5T_IEEE_F32BE, th_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(th_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, th_data);

	// omega
	om_dims[0] = 32;
	om_dims[1] = 32;
	om_dims[2] = 32;

	om_space = H5Screate_simple(3, om_dims, NULL);

	om_set = H5Dcreate2(file_id, "/euler/omega", H5T_IEEE_F32BE, om_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(om_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, om_data);

	// grain
	grain_dims[0] = 32;
	grain_dims[1] = 32;
	grain_dims[2] = 32;

	grain_space = H5Screate_simple(3, grain_dims, NULL); 
	grain_set = H5Dcreate2(file_id, "/grain", H5T_STD_I32BE, grain_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(grain_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, grain_data);

	// phase
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
	status = H5Dclose(ph_set);
	status = H5Sclose(ph_space);
	status = H5Dclose(th_set);
	status = H5Sclose(th_space);
	status = H5Dclose(om_set);
	status = H5Sclose(om_space);
	status = H5Gclose(group_id);
	status = H5Fclose(file_id);

	return 0;
}
