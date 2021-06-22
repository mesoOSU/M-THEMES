#include <hdf5.h>
#include <mpi.h>
#include <unistd.h>
#include "evp.h"

hid_t CreateParallelHDF(char *name) {
	hid_t fapl_id, file_id;
	const int BUFFER_SIZE = 100;
	char fname[BUFFER_SIZE];

	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	snprintf(fname, BUFFER_SIZE, "%s.h5", name);
	file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);

	return file_id;
}/*end CreateParallelHDF()*/

void WriteHDFDataset(hid_t file_id, int dimensions, char *dset_path, hid_t type_id, const void *buffer) {
	hid_t space_id, set_id, xf_id, mem_space; // hdf ids
	hsize_t dims[] = {(unsigned int)CellDim[0], (unsigned int)CellDim[1], (unsigned int)CellDim[2], 3};

	// define hyperslab
	hsize_t slab_count[] = {(unsigned int)CellDim[0]/NumPE, (unsigned int)CellDim[1], (unsigned int)CellDim[2], 3};
	hsize_t slab_offset[] = {mpirank*slab_count[0], 0, 0, 0};

	// define memory space
	mem_space = H5Screate_simple(dimensions, slab_count, NULL);

	// create dataspace
	space_id = H5Screate_simple(dimensions, dims, NULL);
	H5Sselect_hyperslab(space_id, H5S_SELECT_SET, slab_offset, NULL, slab_count, NULL);

	// create dataset
	set_id = H5Dcreate2(file_id, dset_path, type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	xf_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(xf_id, H5FD_MPIO_COLLECTIVE);

	// write
	H5Dwrite(set_id, type_id, mem_space, space_id, xf_id, buffer);

	// close hdf spaces and set
	H5Sclose(space_id);
	H5Sclose(mem_space);
	H5Dclose(set_id);
}/*end WriteHDFDataset()*/

void WriteGrainHDF(char *s, int step)
{
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	hid_t file_id;

	file_id = CreateParallelHDF(s);

	snprintf(dset_name, BUFFER_SIZE, "%s_S%06d", s, step);
	WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_INT, grain_f);
}/*end WriteGrainHDF()*/

#ifdef PF_DRX
void WriteGrainPFHDF(char *s, int step)
{
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	hid_t file_id;

	file_id = CreateParallelHDF(s);

	snprintf(dset_name, BUFFER_SIZE, "%s_S%06d", s, step);
	WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_INT, gID_rex);
}/*end WriteGrainPFHDF()*/
#endif

void WriteDisgradHDF(char *s, int step)
{
	char dset_name[100];
	hid_t file_id;
	real *tmpVector;
	int i;

	file_id = CreateParallelHDF(s);

	AllocMem(tmpVector, lsize, real);

	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = DisGrad[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = DisGrad[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = DisGrad[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = DisGrad[pIDX][0][1];
			}
		}

		sprintf(dset_name, "%01d_S%06d", i+1, step);
		WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);
	}

	H5Fclose(file_id);
	free(tmpVector);
}/*end WriteDisgradHDF()*/

void WriteElsHDF(char *s, int step)
{
	char dset_name[100];
	hid_t file_id;
	real *tmpVector;
	int i;

	file_id = CreateParallelHDF(s);
	AllocMem(tmpVector, lsize, real);

	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = DisGrad[pIDX][i][i] - Eps[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = (DisGrad[pIDX][1][2]+DisGrad[pIDX][2][1])/2.0 - Eps[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = (DisGrad[pIDX][0][2]+DisGrad[pIDX][2][0])/2.0 - Eps[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = (DisGrad[pIDX][0][1]+DisGrad[pIDX][1][0])/2.0 - Eps[pIDX][0][1];
			}
		}

		sprintf(dset_name, "%01d_S%06d", i+1, step);
		WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);
	}

	H5Fclose(file_id);
	free(tmpVector);
}/* WriteElsHDF()*/

void WriteEpsHDF(char *s, int step)
{
	hid_t file_id;
	char dset_name[100];
	real *tmpVector;

	file_id = CreateParallelHDF(s);
	AllocMem(tmpVector, lsize, real);

	for (int i=0;i<6;i++) {
		local_loop{
			if(i<3){
				tmpVector[pIDX] = Eps[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = Eps[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = Eps[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = Eps[pIDX][0][1];
			}
		}

		sprintf(dset_name, "%01d_S%06d", i+1, step);
		WriteHDFDataset(file_id, 3, dset_name, H5T_IEEE_F32BE, tmpVector);
	}

	H5Fclose(file_id);
	free(tmpVector);
}/*end WriteEpsHDF()*/

void WriteSigHDF(char *s, int step)
{
	char dset_name[100];
	hid_t file_id;
	real *tmpVector;
	int i;

	file_id = CreateParallelHDF(s);
	AllocMem(tmpVector, lsize, real);

	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = Sig[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = Sig[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = Sig[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = Sig[pIDX][0][1];
			}
		}

		sprintf(dset_name, "%01d_S%06d", i+1, step);
		WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);
	}

	H5Fclose(file_id);
	free(tmpVector);
}/*end WriteSigHDF()*/

void WriteEdotHDF(char *s, int step)
{
	char dset_name[100];
	hid_t file_id;
	real *tmpVector;
	int i;

	file_id = CreateParallelHDF(s);
	AllocMem(tmpVector, lsize, real);

	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = Edot[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = Edot[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = Edot[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = Edot[pIDX][0][1];
			}
		}

		sprintf(dset_name, "%01d_S%06d", i+1, step);
		WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);
	}

	H5Fclose(file_id);
	free(tmpVector);
}/*end WriteEdotHDF()*/

void WriteTextureHDF(char *s, int step)
{
	// declare vars/allocate memory
	real t1, ph, t2;
	ten2nd sa2xt;
	hid_t file_id;
	char dset_name[100];

	real *angleVector;
	real *grainVector;
	real *phaseVector;

	AllocMem(angleVector, lsize * 3, real);
	AllocMem(grainVector, lsize, real);
	AllocMem(phaseVector, lsize, real);

	// populate vectors
	local_loop{
		T2_loop{
			sa2xt[mi][mj] = TranMat_xt2sa[pIDX][mj][mi];	// transpose
		}
		EulerToTransMatrix(&t1,&ph,&t2,sa2xt,1);
		angleVector[3*pIDX]= t1;
		angleVector[3*pIDX + 1]= ph;
		angleVector[3*pIDX + 2]= t2;
		grainVector[pIDX] = grain_f[pIDX];
		phaseVector[pIDX] = phase_f[pIDX];
	}

	// make file
	file_id = CreateParallelHDF(s);

	// write
	sprintf(dset_name, "Angles_S%06d", step);
	WriteHDFDataset(file_id, 4, dset_name, H5T_NATIVE_DOUBLE, angleVector);

	sprintf(dset_name, "Grains_S%06d", step);
	WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, grainVector);

	sprintf(dset_name, "Phases_S%06d", step);
	WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, phaseVector);

	// free
	H5Fclose(file_id);
	free(angleVector);
	free(grainVector);
	free(phaseVector);

}/*end WriteTextureHDF()*/

void WriteSLIPHDF(char *s, int step)
{
	// init vars/allocate memory
	char dset_name[100];
	hid_t file_id;
	real *SSParentVector[12]; // parent array to store 12 individual SS vectors
	real *grainVector = (real *)malloc(lsize * sizeof(real));

	for (int i = 0; i < 12; i++) {
		real *SSChildVector = (real *)malloc(lsize * sizeof(real));
		SSParentVector[i] = SSChildVector;
	}

	// populate vectors
	local_loop {
		grainVector[pIDX] = grain_f[pIDX];

		for (int i = 0; i < 12; i++) {
			SSParentVector[i][pIDX] = trial_gam_acum[pIDX][i];
		}
	}

	// make file
	file_id = CreateParallelHDF(s);

	// write
	sprintf(dset_name, "Grains_S%06d", step);
	WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, grainVector);
	for (int i = 0; i < 12; i++) {
		sprintf(dset_name, "SS%d_S%04d", i+1, step);

		WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, SSParentVector[i]);
	}

	// free
	H5Fclose(file_id);
	for (int i = 0; i < 12; i++) {
		free(SSParentVector[i]);
	}
}/*end WriteSLIPHDF()*/

void WriteNewPositionHDF(char *s, int step) {
	// init vars/alloc mem
	char dset_name[100];
	hid_t file_id;
	real *SSParentVector[3]; // parent array to store 3 individual SS vectors

	for (int i = 0; i < 3; i++) {
		real *SSChildVector = (real *)malloc(lsize * sizeof(real));
		SSParentVector[i] = SSChildVector;
	}

	// populate vectors
	local_loop {
		for (int i = 0; i < 3; i++) {
			SSParentVector[i][pIDX] = displacement_fluct[pIDX][i];
		}
	}

	// make file
	file_id = CreateParallelHDF(s);

	// write
	for (int i = 0; i < 3; i++) {
		sprintf(dset_name, "SS%d_S%04d", i+1, step);
		WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, SSParentVector[i]);
	}

	// free
	H5Fclose(file_id);
	for (int i = 0; i < 3; i++) {
		free(SSParentVector[i]);
	}
}/*end WriteNewPositionHDF()*/

void WriteDDHDF(char *s, int step)
{
	// init vars/alloc mem
	char dset_name[100];
	hid_t file_id;
	real *grainVector = (real *) malloc(lsize * sizeof(real));
	real *SSParentVector[12];

	for (int i = 0; i < 12; i++) {
		real *SSChildVector = (real *)malloc(lsize * sizeof(real));
		SSParentVector[i] = SSChildVector;
	}

	// populate
	local_loop {
		grainVector[pIDX] = grain_f[pIDX];
		for (int i = 0; i < 12; i++) {
			SSParentVector[i][pIDX] = rho_s[pIDX][i];
		}
	}

	// make file
	file_id = CreateParallelHDF(s);

	// write
	for (int i = 0; i < 12; i++) {
		sprintf(dset_name, "SS%d_S%04d", i+1, step);
		WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, SSParentVector[i]);
	}

	// free
	H5Fclose(file_id);
	for (int i = 0; i < 12; i++) {
		free(SSParentVector[i]);
	}
}

#ifdef DD_BASED_FLAG
void WriteRhoHDF(char *s, char *type, int step)
{
	char fname[100];
	char dset_name[100];
	hid_t file_id;
	int i;
	int jph;

	local_loop{
		jph = phase_f[pIDX];
		rho_tot[pIDX] = 0.0;
		/* SSD */
		if(type[0]=='s'||type[0]=='S'){
			for(i=0;i<nSYS[jph-1];i++){
				rho_tot[pIDX] += rho_s[pIDX][i];
			}
		}
		else if(type[0]=='m'||type[0]=='M'){
			for(i=0;i<nSYS[jph-1];i++){
				rho_tot[pIDX] += rho_m[pIDX][i];
			}
		}
		else if(type[0]=='p'||type[0]=='P'){
			for(i=0;i<nSYS[jph-1];i++){
				rho_tot[pIDX] += rho_P[pIDX][i];
			}
		}
		else if(type[0]=='f'||type[0]=='F'){
			for(i=0;i<nSYS[jph-1];i++){
				rho_tot[pIDX] += rho_F[pIDX][i];
			}
		}
		else if(type[0]=='g'||type[0]=='G'){
			for(i=0;i<nSYS[jph-1];i++){
				rho_tot[pIDX] += sqrt(rho_g1[pIDX][i]*rho_g1[pIDX][i]+rho_g2[pIDX][i]*rho_g2[pIDX][i]+rho_g3[pIDX][i]*rho_g3[pIDX][i]);
			}
		}
		else{
			PError("Wrong type of dislcoations to write!",720);
		}
	}

	// Rho is designed to be called multiple times, so check if HDF file already exists and if so use it
	snprintf(fname, 100, "%s.h5", s);
	if (access(fname, F_OK) != 0) {
		file_id = CreateParallelHDF(s);
	} else {
		file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
	}

	sprintf(dset_name, "%s_S%06d", type, step);
	WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, rho_tot);

	H5Fclose(file_id);
}/*end WriteRhoHDF()*/

void WriteRhoDotHDF(char *s, int step)
{
	char dset_name[100];
	hid_t file_id;
	real *tmpVector;

	AllocMem(tmpVector, lsize, real);

	file_id = CreateParallelHDF(s);

	for (int i=0;i<NSYSMX;i++) {
		// SSD
		local_loop{
			tmpVector[pIDX] = rho_dot_s[pIDX][i];
		}

		sprintf(dset_name, "SSD_slip%02d_S%06d", i+1, step);
		WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);

		// GND-I
		local_loop{
			tmpVector[pIDX] = rho_dot_g1[pIDX][i];
		}

		sprintf(dset_name, "GND1_slip%02d_S%06d", i+1,step);
		WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);

		// GND-II
		local_loop{
			tmpVector[pIDX] = rho_dot_g2[pIDX][i];
		}

		sprintf(dset_name, "GND2_slip%02d_S%06d", i+1,step);
		WriteHDFDataset(file_id, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);
	}

	free(tmpVector);
	H5Fclose(file_id);

	return;
}/*end WriteRhoDotHDF()*/
#endif
