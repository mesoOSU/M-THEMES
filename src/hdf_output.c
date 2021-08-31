#include <hdf5.h>
#include <mpi.h>
#include <unistd.h>
#include "evp.h"

hid_t CreateParallelHDF(char *name) {
	hid_t fapl_id, file_id;

	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);

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

void WriteGrainHDF(hid_t file_id, char *group_name)
{
	WriteHDFDataset(file_id, 3, group_name, H5T_NATIVE_INT, grain_f);
}/*end WriteGrainHDF()*/

#ifdef PF_DRX
void WriteGrainPFHDF(hid_t file_id, char *group_name)
{
	WriteHDFDataset(file_id, 3, group_name, H5T_NATIVE_INT, gID_rex);
}/*end WriteGrainPFHDF()*/
#endif

void WriteDisgradHDF(hid_t file_id, char *group_name)
{
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	real *tmpVector;
	int i;
	hid_t group;

	group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

		snprintf(dset_name, BUFFER_SIZE, "%01d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);
	}

	H5Gclose(group);
	free(tmpVector);
}/*end WriteDisgradHDF()*/

void WriteElsHDF(hid_t file_id, char *group_name)
{
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	real *tmpVector;
	int i;
	hid_t group;

	group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

		snprintf(dset_name, BUFFER_SIZE, "%01d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);
	}

	H5Gclose(group);
	free(tmpVector);
}/* WriteElsHDF()*/

void WriteEpsHDF(hid_t file_id, char *group_name)
{
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	real *tmpVector;
	hid_t group;

	group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

		snprintf(dset_name, BUFFER_SIZE, "%01d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_IEEE_F32BE, tmpVector);
	}

	H5Gclose(group);
	free(tmpVector);
}/*end WriteEpsHDF()*/

void WriteSigHDF(hid_t file_id, char *group_name)
{
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	real *tmpVector;
	int i;
	hid_t group;

	group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

		snprintf(dset_name, BUFFER_SIZE, "%01d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);
	}

	H5Gclose(group);
	free(tmpVector);
}/*end WriteSigHDF()*/

void WriteEdotHDF(hid_t file_id, char *group_name)
{
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	real *tmpVector;
	int i;
	hid_t group;

	group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

		snprintf(dset_name, BUFFER_SIZE, "%01d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);
	}

	H5Gclose(group);
	free(tmpVector);
}/*end WriteEdotHDF()*/

void WriteTextureHDF(hid_t file_id, char *group_name)
{
	// declare vars/allocate memory
	real t1, ph, t2;
	ten2nd sa2xt;
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	hid_t group;

	group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

	// write
	snprintf(dset_name, BUFFER_SIZE, "Euler Angles");
	WriteHDFDataset(group, 4, dset_name, H5T_NATIVE_DOUBLE, angleVector);

	snprintf(dset_name, BUFFER_SIZE, "Grains");
	WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, grainVector);

	snprintf(dset_name, BUFFER_SIZE, "Phases");
	WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, phaseVector);

	// free
	H5Gclose(group);
	free(angleVector);
	free(grainVector);
	free(phaseVector);
}/*end WriteTextureHDF()*/

void WriteSLIPHDF(hid_t file_id, char *group_name)
{
	// init vars/allocate memory
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	real *SSParentVector[12]; // parent array to store 12 individual SS vectors
	real *grainVector = (real *)malloc(lsize * sizeof(real));
	hid_t group;

	group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

	// write
	snprintf(dset_name, BUFFER_SIZE, "Grains");
	WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, grainVector);
	for (int i = 0; i < 12; i++) {
		snprintf(dset_name, BUFFER_SIZE, "SS%d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, SSParentVector[i]);
	}

	// free
	H5Gclose(group);
	for (int i = 0; i < 12; i++) {
		free(SSParentVector[i]);
	}
}/*end WriteSLIPHDF()*/

void WriteNewPositionHDF(hid_t file_id, char *group_name) {
	// init vars/alloc mem
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	real *SSParentVector[3]; // parent array to store 3 individual SS vectors
	hid_t group;

	group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

	// write
	for (int i = 0; i < 3; i++) {
		snprintf(dset_name, BUFFER_SIZE, "SS%d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, SSParentVector[i]);
	}

	// free
	H5Gclose(group);
	for (int i = 0; i < 3; i++) {
		free(SSParentVector[i]);
	}
}/*end WriteNewPositionHDF()*/

void WriteDDHDF(hid_t file_id, char *group_name)
{
	// init vars/alloc mem
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	real *grainVector = (real *) malloc(lsize * sizeof(real));
	real *SSParentVector[12];
	hid_t group;

	group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

	// write
	for (int i = 0; i < 12; i++) {
		snprintf(dset_name, BUFFER_SIZE, "SS%d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, SSParentVector[i]);
	}

	// free
	H5Gclose(group);
	for (int i = 0; i < 12; i++) {
		free(SSParentVector[i]);
	}
}

#ifdef DD_BASED_FLAG
void WriteRhoHDF(hid_t file_id, char *group_name, char *type)
{
	int i;
	int jph;
	hid_t group;

	if (H5Lexists(file_id, group_name, H5P_DEFAULT) == 1) { // check if group already exists
		group = H5Gopen2(file_id, group_name, H5P_DEFAULT);
	} else {
		group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}

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

	WriteHDFDataset(group, 3, type, H5T_NATIVE_DOUBLE, rho_tot);

	H5Gclose(group);
}/*end WriteRhoHDF()*/

void WriteRhoDotHDF(hid_t file_id, char *group_name)
{
	const int BUFFER_SIZE = 100;
	char dset_name[BUFFER_SIZE];
	real *tmpVector;
	hid_t group;

	group = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	AllocMem(tmpVector, lsize, real);

	for (int i=0;i<NSYSMX;i++) {
		// SSD
		local_loop{
			tmpVector[pIDX] = rho_dot_s[pIDX][i];
		}

		snprintf(dset_name, BUFFER_SIZE, "SSD_slip%02d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);

		// GND-I
		local_loop{
			tmpVector[pIDX] = rho_dot_g1[pIDX][i];
		}

		snprintf(dset_name, BUFFER_SIZE, "GND1_slip%02d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);

		// GND-II
		local_loop{
			tmpVector[pIDX] = rho_dot_g2[pIDX][i];
		}

		snprintf(dset_name, BUFFER_SIZE, "GND2_slip%02d", i+1);
		WriteHDFDataset(group, 3, dset_name, H5T_NATIVE_DOUBLE, tmpVector);
	}

	H5Gclose(group);
	free(tmpVector);
}/*end WriteRhoDotHDF()*/
#endif
