/*===========================================
  A phase-field model for recrystallization,
  realized by the function PF_XRX(...).

  Usage of PF_XRX(...):

  Argument list:
  --------------------------------------------
   input
  --------------------------------------------
	1. L_x: grid length of x-dimension
	2. L_y: grid length of y-dimension
	3. L_z: grid length of z-dimension
	4. NumPE: # of PEs used in external MPI
		environment. Currently only only
		NumPE = 2^n, where n is an integer)
	5. mpirank: external MPI rank
	6. E_gb: GB energy in [J/m^2]
	7. M_gb: GB mobility in [m^4/(MJ*s)]
	8. mu: Shear modulus in [GPa]
	9. bb: Burger vector magnitude in [nm]
	10. alpha: Scaler of stored energy in [1]
	11. gID: grain IDs
	12. gID_rex: indicator of recrystallized
			grain (==1) or not (==0)
	13. rho: total dislocation density
	14. Steps: # of PF relaxation steps
  --------------------------------------------
   output
  --------------------------------------------
	15. gID_new: grain IDs after relaxation
	16. grex_new: grain recrystallized
			indicators after relaxation


						-- PY Zhao, 09/16/2014
  ===========================================*/
#ifndef H_PFXRX
#define H_PFXRX
#include<iostream>
#include<fstream>
#include<cmath>
#include<string.h>
#include<stdio.h>
#include<mpi.h>
#include<assert.h>
#include "PFgrid.h"

void PF_XRX(int L_x, int L_y, int L_z, int NumPE, int mpirank,
		double E_gb, double M_gb, double mu, double bb, double alpha,
		int* gID, double* rho, int Steps, int *gID_new, int count,int Flagwrite,double* rho_avg);

double laplacian(const double& right, const double& left, const double& front, const double& back, const double& up,const double& down,const double& right_up, const double& right_down, const double& right_front, const double& right_back, const double& left_up,const double& left_down,const double& left_front, const double& left_back, const double& center_up_front, const double& center_up_back, const double& center_down_front,const double& center_down_back, const double& c);

void Update_gID(PFgrid3D& gd_rex, int *gID, int L_x, int L_y, int L_z, int lnx, int lny, int lnz, int lsize_d ,int istep, int iter_step ,int mpirank);

void Update_gID_final(PFgrid3D& gd_rex,int *gID, int L_x, int L_y, int L_z, int lnx, int lny, int lnz, int mpirank);

void clamping_eta(PFgrid3D& gd_rex, int lnx, int lny, int lnz);

void boundary_conditions(PFgrid3D& gd, int mpirank, int NumPE, int l_x, int l_y, int l_z, int flag);

void Update_grid(PFgrid3D& gd_rex, PFgrid3D& gd_dd, int *gID, int L_x, int L_y, int L_z, int l_x, int l_y, int l_z,
		int mpirank, int NumPE,
		double M_coeff, double kappa, double m, double *E_disl,int iter_step, double *rho_avg, int lsize_d ,int istep, double *DefEng_buf, int *GID); 
	

void GetStatistics(PFgrid3D& gd_rex, int lnx, int lny, int lnz, int Nxyz,int iter_step,int mpirank,int lsize_d,int Flagwrite, int count_1);

void mesh_extend(double *DefEng_buf, int *GID, int mpirank, int NumPE, int L_x, int L_y , int L_z);


void WriteIDMPI(int *ID,int mpirank,int lsize,const char* fn);

struct dd_Info{
	int ID;
	double avgdd;
};

#endif
