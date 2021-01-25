#include "evp.h"
static real NR_quality(voigt sig, voigt xlambda, voigt eps, voigt strain, voigt66 sg, int idx,int jphi)
{
	voigt edotp;
	voigt66 d_edotp;
	voigt tot_eps;	// total strain
  voigt res;
  real f;
  int i,j;
#ifdef DD_BASED_FLAG
#ifdef DD_POWER_LAW
				StrainRate_Orowan_POWER(sig,edotp,d_edotp,idx,jphi);
#else
				StrainRate_Orowan(sig,edotp,d_edotp,idx,jphi);
#endif
#else
				StrainRate_eval(sig,edotp, d_edotp,idx, jphi);
#endif

				/* tot_eps is the total  strain and eps6 is the 
				   current plastic strain, i.e. Eq. 4 */
				for(i=0;i<6;i++){
					tot_eps[i] = eps[i] + edotp[i]*TimeStep;
					for(j=0;j<6;j++){
						tot_eps[i] += sg[i][j]*sig[j];
					}
				}

				// calculate the residual R, Eq. 16
				for(i=0;i<6;i++){
					res[i] = sig[i] - xlambda[i];
					for(j=0;j<6;j++){
						res[i] += C066[i][j]*(tot_eps[j]-strain[j]);
					}
				}
        for(f=0.0,i=0;i<6;i++){
          f += res[i]*res[i];
        }

        return f/2.0;
}/*end NR_quality()*/

#ifdef PF_DRX
static void Gradient_Rho(void)
{
	int px_p, px_m, py_p, py_m, pz_p, pz_m;
	int pIDX_xp, pIDX_xm, pIDX_yp, pIDX_ym, pIDX_zp, pIDX_zm;

	/* Becuase of the slab decomposition, calculating the derivatives
	   requires the data from neighboring PEs */
	if(NumPE>1){	// multiple PEs, use MPI_Send() and MPI_Recv()
		for(px=0;px<1;px++){
			for(py=0;py<CellDim[1];py++){
				for(pz=0;pz<CellDim[2];pz++){
					local_pxp_send[py*CellDim[2]+pz] = rho_tot[pIDX];
				}
			}
		}
		for(px=lnx-1;px<lnx;px++){
			for(py=0;py<CellDim[1];py++){
				for(pz=0;pz<CellDim[2];pz++){
					local_pxm_send[py*CellDim[2]+pz] = rho_tot[pIDX];
				}
			}
		}
		if(mpirank==0){
			MPI_Sendrecv(local_pxm_send,CellDim[1]*CellDim[2],MPI_real,1,0,
					local_pxp_recv,CellDim[1]*CellDim[2],MPI_real,1,3,
					MPI_COMM_WORLD, &status);
			MPI_Sendrecv(local_pxp_send,CellDim[1]*CellDim[2],MPI_real,NumPE-1,1,
					local_pxm_recv,CellDim[1]*CellDim[2],MPI_real,NumPE-1,2*(NumPE-1),
					MPI_COMM_WORLD, &status);
		}
		else if(mpirank==NumPE-1){
			MPI_Sendrecv(local_pxp_send,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*mpirank+1,
					local_pxm_recv,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*(mpirank-1),
					MPI_COMM_WORLD, &status);
			MPI_Sendrecv(local_pxm_send,CellDim[1]*CellDim[2],MPI_real,0,2*mpirank,
					local_pxp_recv,CellDim[1]*CellDim[2],MPI_real,0,1,
					MPI_COMM_WORLD, &status);
		}
		else{
			MPI_Sendrecv(local_pxp_send,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*mpirank+1,
					local_pxm_recv,CellDim[1]*CellDim[2],MPI_real,mpirank-1,2*(mpirank-1),
					MPI_COMM_WORLD, &status);
			MPI_Sendrecv(local_pxm_send,CellDim[1]*CellDim[2],MPI_real,mpirank+1,2*mpirank,
					local_pxp_recv,CellDim[1]*CellDim[2],MPI_real,mpirank+1,2*(mpirank+1)+1,
					MPI_COMM_WORLD, &status);
		}
	}
	else{	// Only one PE (serial run), no need for communication
		for(px=0;px<1;px++){
			for(py=0;py<CellDim[1];py++){
				for(pz=0;pz<CellDim[2];pz++){
					local_pxp_recv[py*CellDim[2]+pz] = rho_tot[pIDX];
				}
			}
		}
		for(px=lnx-1;px<lnx;px++){
			for(py=0;py<CellDim[1];py++){
				for(pz=0;pz<CellDim[2];pz++){
					local_pxm_recv[py*CellDim[2]+pz] = rho_tot[pIDX];
				}
			}
		}

	}

	for(px=0;px<lnx;px++){
		px_p = px+1; px_m = px-1;
		if(px_p>=lnx) px_p = 0;
		if(px_m<0) px_m = lnx-1;
		for(py=0;py<CellDim[1];py++){
			py_p = py+1; py_m = py-1;
			if(py_p>=CellDim[1]) py_p = 0;
			if(py_m<0) py_m = CellDim[1]-1;
			for(pz=0;pz<CellDim[2];pz++){
				pz_p = pz+1; pz_m = pz-1;
				if(pz_p>=CellDim[2]) pz_p = 0;
				if(pz_m<0) pz_m = CellDim[2]-1;
				pIDX_xp = (px_p*CellDim[1]+py)*CellDim[2]+pz;
				pIDX_xm = (px_m*CellDim[1]+py)*CellDim[2]+pz;
				pIDX_yp = (px*CellDim[1]+py_p)*CellDim[2]+pz;
				pIDX_ym = (px*CellDim[1]+py_m)*CellDim[2]+pz;
				pIDX_zp = (px*CellDim[1]+py)*CellDim[2]+pz_p;
				pIDX_zm = (px*CellDim[1]+py)*CellDim[2]+pz_m;


				if(px==0){
					diff_rho[pIDX] =
						local_pxm_recv[py*CellDim[2]+pz]+rho_tot[pIDX_xp]+
						rho_tot[pIDX_ym]+rho_tot[pIDX_yp]+
						rho_tot[pIDX_zm]+rho_tot[pIDX_zp];

					diff_rho[pIDX] = rho_tot[pIDX] - diff_rho[pIDX]/6.0;
				}
				else if(px==lnx-1){
					diff_rho[pIDX] =
						local_pxp_recv[py*CellDim[2]+pz]+rho_tot[pIDX_xm]+
						rho_tot[pIDX_ym]+rho_tot[pIDX_yp]+
						rho_tot[pIDX_zm]+rho_tot[pIDX_zp];

					diff_rho[pIDX] = rho_tot[pIDX] - diff_rho[pIDX]/6.0;
				}
				else{
					diff_rho[pIDX] =
						rho_tot[pIDX_xm]+rho_tot[pIDX_xp]+
						rho_tot[pIDX_ym]+rho_tot[pIDX_yp]+
						rho_tot[pIDX_zm]+rho_tot[pIDX_zp];

					diff_rho[pIDX] = rho_tot[pIDX] - diff_rho[pIDX]/6.0;
				}
			}
		}
	}


	return;
}/*end Gradient_Rho()*/

int NucleationCheckDRX(int STEP, real* kc, real alpha)
{
	int flag,flag_local;
	int i, jph;
	real t1,ph1,t2;
	ten2nd sa2xt;
	static std::vector<G_Info> newID_list;
		static std::vector<N_Info> nucleiID_list;

#ifdef NUCL_RHO_DIFF	
	/* calculate the total dislocation density */
	local_loop{
		jph = phase_f[pIDX];
		rho_tot[pIDX] = 0.0;
		for(i=0;i<nSYS[jph-1];i++){
			rho_tot[pIDX] += rho_s[pIDX][i]+rho_m[pIDX][i]
				+sqrt(pow(rho_g1[pIDX][i],2.0)+pow(rho_g2[pIDX][i],2.0)+pow(rho_g3[pIDX][i],2.0));
		}
	}

	/* update the dislocation difference, the driving force to DRX nucleation */
	Gradient_Rho();
#else	// use the local disl. density to check nucleation
	/* calculate the total dislocation density */
	local_loop{
		jph = phase_f[pIDX];
		rho_tot[pIDX] = 0.0;
		for(i=0;i<nSYS[jph-1];i++){
			rho_tot[pIDX] += rho_s[pIDX][i]+rho_m[pIDX][i]
				+sqrt(pow(rho_g1[pIDX][i],2.0)+pow(rho_g2[pIDX][i],2.0)+pow(rho_g3[pIDX][i],2.0));
		}
		diff_rho[pIDX] = rho_tot[pIDX]/12;
	}
#endif

	/* using the statistical model to calculate the FFT cell transformation rate
	   and compare with a uniformly generated random number over [0,1] */
	flag_local=0;
	local_loop{
#ifndef GB_BULGING_ONLY
	if(Nucl_Static_Flag!=1){
		double dice = gsl_rng_uniform(RandInstance);
		double prob = 1.0-exp(-1.0*pow(diff_rho[pIDX]/kc[pIDX],alpha));
		double age = (TimeTot-t_last[pIDX])/TimeStep;
		if((prob>dice)&&(age>age_drx)){
			gID_rex[pIDX]=1;
			//int index = gID_list.back().ID+newID_list.size()+mpirank*lsize+1;	// the grain IDs may not be consecutive integers from this point
			int index = newID_list.size()+1;	// the grain IDs will keep consecutive
			grain_f[pIDX]=index;
			double ph = 2.0*PI*gsl_rng_uniform(RandInstance);
			double th = PI*gsl_rng_uniform(RandInstance);
			double om = 2.0*PI*gsl_rng_uniform(RandInstance);
			struct G_Info tmp_ID={index,ph,th,om};
			newID_list.push_back(tmp_ID);
			flag_local++;	// flag_local equals to the size of newID_list
		}else{
			gID_rex[pIDX]=0;
		}
	}
	else{
		double age = (TimeTot-t_last[pIDX])/TimeStep;
		if((diff_rho[pIDX]>kappa_drx[pIDX])&&(age>age_drx)){
			gID_rex[pIDX]=1;
			//int index = gID_list.back().ID+newID_list.size()+mpirank*lsize+1;	// the grain IDs may not be consecutive integers from this point
			int index = newID_list.size()+1;	// the grain IDs will keep consecutive
			grain_f[pIDX]=index;
			double ph = 2.0*PI*gsl_rng_uniform(RandInstance);
			double th = PI*gsl_rng_uniform(RandInstance);
			double om = 2.0*PI*gsl_rng_uniform(RandInstance);
			struct G_Info tmp_ID={index,ph,th,om};
			newID_list.push_back(tmp_ID);
			flag_local++;	// flag_local equals to the size of newID_list
      //kc[pIDX] = zeta_kc*k_c;
		}else{
			gID_rex[pIDX]=0;
		}
	}
#else
		if(GB_indicator[pIDX]>0){
		 //   if((mpirank == 0) && (STEP%1600 == 0)) {
		 //       printf("%le\n",diff_rho[pIDX]);
		 //   }
	if(Nucl_Static_Flag!=1){
			double dice = gsl_rng_uniform(RandInstance);
			double prob = 1.0-exp(-1.0*pow(diff_rho[pIDX]/kc[pIDX],alpha));
			double age = (TimeTot-t_last[pIDX])/TimeStep;
			if((prob>dice)&&(age>age_drx)){
				gID_rex[pIDX]=1;
				//int index = gID_list.back().ID+newID_list.size()+mpirank*lsize+1;
				int index = newID_list.size()+1;	// the grain IDs will keep consecutive
				grain_f[pIDX]=index;
				double ph = 2.0*PI*gsl_rng_uniform(RandInstance);
				double th = PI*gsl_rng_uniform(RandInstance);
				double om = 2.0*PI*gsl_rng_uniform(RandInstance);
				struct G_Info tmp_ID={index,ph,th,om};
				newID_list.push_back(tmp_ID);
				flag_local++;	// flag_local equals to the size of newID_list

        /* record the nucleation type */
//        fprintf(fp_stat,"%d\t%d\t%d\n", STEP, sIDX, GB_type[pIDX]);
//        fflush(fp_stat);
			}else{
				gID_rex[pIDX]=0;
			}
	}
	else{
			double age = (TimeTot-t_last[pIDX])/TimeStep;
		//	double dice = gsl_rng_uniform(RandInstance);
			//double prob = (1.0-exp(-pow(e_vm/0.12,30)))*(1.0-exp(-pow(diff_rho[pIDX]/kc[pIDX],alpha)));
		//		double prob = 1.0-exp(-1.0*pow(diff_rho[pIDX]/kc[pIDX],alpha));
		//	if((prob>dice)&&(age>age_drx)){
	//	if(sIDX == 8*8*8 || sIDX == 12*12*12 || sIDX == 14*11*12 || sIDX == 9*10*13){
		if( (diff_rho[pIDX]>kappa_drx[pIDX])){	
				gID_rex[pIDX]=1;
				nuclei[pIDX] = 1;
				swap_check[pIDX]=1;
					T2_loop{
			sa2xt[mi][mj] = TranMat_xt2sa[pIDX][mj][mi];	// transpose
		}
		EulerToTransMatrix(&t1,&ph1,&t2,sa2xt,1);
				//int index = gID_list.back().ID+newID_list.size()+mpirank*lsize+1;	// the grain IDs may not be consecutive integers from this point
				int index = newID_list.size()+1;	// the grain IDs will keep consecutive
				grain_f[pIDX]=index;
				double ph = t1;//2.0*PI*gsl_rng_uniform(RandInstance);
				double th = ph1;//PI*gsl_rng_uniform(RandInstance);
				double om = t2;//2.0*PI*gsl_rng_uniform(RandInstance);
				struct G_Info tmp_ID={index,ph,th,om};
				newID_list.push_back(tmp_ID);
				flag_local++;	// flag_local equals to the size of newID_list
        kc[pIDX]=zeta_kc*k_c;

        /* record the nucleation type */
//        fprintf(fp_stat,"%d\t%d\t%d\n", STEP, sIDX, GB_type[pIDX]);
//        fflush(fp_stat);
			}
			else{
				nuclei[pIDX]=0;
			}
	}
		}	// GB_indicator
		else{
			nuclei[pIDX]=0;
		}
#endif
	}
	MPI_Allreduce(&flag_local, &flag, 1, MPI_INT,
			MPI_SUM, MPI_COMM_WORLD);
	
	if(flag==0){
		return flag;	// equals to the # of newly nucleated grains
	}
	else{
		/* update the grain IDs in the local  newID_list based on the previous node info. */
		int pre_size = 0;
		int local_size = newID_list.size();
		if(NumPE>1){
			if(mpirank==0){
				MPI_Send(&local_size, 1, MPI_INT,
						1, 0, MPI_COMM_WORLD);
			}
			else if(mpirank!=(NumPE-1)){
				MPI_Recv(&pre_size, 1, MPI_INT,
						mpirank-1, mpirank-1, MPI_COMM_WORLD, &status);
                                int z =0;
                                z = pre_size + local_size;
				MPI_Send(&z, 1, MPI_INT,
						mpirank+1, mpirank, MPI_COMM_WORLD);
			}
			else{
				MPI_Recv(&pre_size, 1, MPI_INT,
						mpirank-1, mpirank-1, MPI_COMM_WORLD, &status);
			}
		}
		else{
			pre_size = 0;
		}

		std::vector<G_Info>::iterator it, end;
		end = newID_list.end();
		for(it=newID_list.begin();it!=end;it++){
			(*it).ID += gID_list.size()+pre_size;	// now the new IDs are consecutive
		}
		/* update the DRX IDs */
		i = 0;
		local_loop{
			if(nuclei[pIDX]==1){
				grain_f[pIDX] = newID_list[i].ID;
				i++;
				int X= px+lxs;
				int Y= py;
				int Z= pz; 
				int index2 = grain_f[pIDX];
			//	double ph = newID_list[i].t1;
			//	double th = newID_list[i].Phi;
			//	double om = newID_list[i].t2;
				struct N_Info tmp_ID2={X,Y,Z,index2};
				nucleiID_list.push_back(tmp_ID2);
			}
		}
		assert(i==newID_list.size());
		//cp_22_07_19
				int global_size2;
		MPI_Allreduce(&local_size, &global_size2, 1,
				MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		std::vector<N_Info> global_nucleiID_list(global_size2);	// only store new IDs
		std::vector<int> nsize1(NumPE);
		MPI_Allgather(&local_size,1,MPI_INT,&nsize1[0],1,MPI_INT,MPI_COMM_WORLD);
		std::vector<int> disps1(NumPE);
		for(int i=0;i<NumPE;i++){
			disps1[i]=(i>0)?(disps1[i-1]+nsize1[i-1]):0;
		}
	
		// define G_ID_Entry
		struct N_Info value1 = {0};
		MPI_Datatype N_ID_Entry;
		int count1=4;
		int block_lens1[4]={1,1,1,1};
		MPI_Aint indices1[4];
		MPI_Datatype old_types1[4]={MPI_INT, MPI_INT, MPI_INT, MPI_INT};
		MPI_Address(&value1.X,&indices1[0]);
		MPI_Address(&value1.Y,&indices1[1]);
		MPI_Address(&value1.Z,&indices1[2]);
		MPI_Address(&value1.ID,&indices1[3]);
		// make relative
		for(int i=count1-1;i>=0;i--)
			indices1[i] -= indices1[0];
		MPI_Type_struct(count1,block_lens1,indices1,old_types1,&N_ID_Entry);
		MPI_Type_commit(&N_ID_Entry);
	
		MPI_Allgatherv(&nucleiID_list[0],local_size,N_ID_Entry,
				&global_nucleiID_list[0],&nsize1[0],&disps1[0],N_ID_Entry,MPI_COMM_WORLD);
			
			
			
		//22_07_19
		std::vector<N_Info>::iterator it2, end2;
			end2 =global_nucleiID_list.end();
		for(it2=global_nucleiID_list.begin();it2!=end2;it2++){
		    local_loop{
		     if (pow(px+lxs-(*it2).X,2)+pow(py-(*it2).Y,2)+pow(pz-(*it2).Z,2)<= nucleus_radius ){
		         grain_f[pIDX]=(*it2).ID;
		         gID_rex[pIDX]=1;
                        //  growth[pIDX] = 1;
                       nuclei[pIDX]=1;
                     //  nucleation_count[pIDX] = 1;
                  // for(i=0;i<nSYS[jph-1];i++) {
                            // rho_s[pIDX][i] = 0.1*rho_s[pIDX][i];
                             //rho_m[pIDX][i] = 0.1*rho_m[pIDX][i];
                            // rho_g1[pIDX][i] = 0.1*rho_g1[pIDX][i];
                          //   rho_g2[pIDX][i] = 0.1*rho_g2[pIDX][i];
                        //     rho_g3[pIDX][i] = 0.1*rho_g3[pIDX][i];
                      //   } 
                    //     rho_tot[pIDX] = 0.0;
                  //       for(i=0;i<nSYS[jph-1];i++){
		//	rho_tot[pIDX] += rho_s[pIDX][i]+rho_m[pIDX][i]
		//		+sqrt(pow(rho_g1[pIDX][i],2.0)+pow(rho_g2[pIDX][i],2.0)+pow(rho_g3[pIDX][i],2.0));
		//}

// T2_loop {
   //   Sig[pIDX][mi][mj] *= 0.75;
    // Eps[pIDX][mi][mj] = 0.0;
     //   Edot[pIDX][mi][mj] = 0.0;
     //   DisGrad[pIDX][mi][mj] = 0.0;
     //  VelGrad[pIDX][mi][mj] = 0.0;
     
  //  }
	   // for(int i=0;i<nSYS[jph-1];i++){
   //    gamdot[pIDX][i] = 0.0;
   // }  
                      //   printf("%d %d" , gID_rex[pIDX],grain_f[pIDX]);
	//	         printf("nuclei %d %d %d %d %d\n", px, py, pz, grain_f[pIDX],gID_rex[pIDX]);
		     }      
		    }
		}
		/* update gID_list */

		/* construct global gID_list */
		int global_size;
		MPI_Allreduce(&local_size, &global_size, 1,
				MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		std::vector<G_Info> global_gID_list(global_size);	// only store new IDs
		std::vector<int> nsize(NumPE);
		MPI_Allgather(&local_size,1,MPI_INT,&nsize[0],1,MPI_INT,MPI_COMM_WORLD);
		std::vector<int> disps(NumPE);
		for(int i=0;i<NumPE;i++){
			disps[i]=(i>0)?(disps[i-1]+nsize[i-1]):0;
		}
	
		// define G_ID_Entry
		struct G_Info value = {0};
		MPI_Datatype G_ID_Entry;
		int count=4;
		int block_lens[4]={1,1,1,1};
		MPI_Aint indices[4];
		MPI_Datatype old_types[4]={MPI_INT, MPI_real, MPI_real, MPI_real};
		MPI_Address(&value,&indices[0]);
		MPI_Address(&value.t1,&indices[1]);
		MPI_Address(&value.Phi,&indices[2]);
		MPI_Address(&value.t2,&indices[3]);
		// make relative
		for(int i=count-1;i>=0;i--)
			indices[i] -= indices[0];
		MPI_Type_struct(count,block_lens,indices,old_types,&G_ID_Entry);
		MPI_Type_commit(&G_ID_Entry);
	
		MPI_Allgatherv(&newID_list[0],local_size,G_ID_Entry,
				&global_gID_list[0],&nsize[0],&disps[0],G_ID_Entry,MPI_COMM_WORLD);
		std::vector<G_Info>::iterator it_g, end_g;
		end_g=global_gID_list.end();
		for(it_g=global_gID_list.begin();it_g!=end_g;it_g++){
			end=gID_list.end();
			int check=0;
			for(it=gID_list.begin();it!=end;it++){
			//    printf("ID=%d\n",(*it_g).ID);
				if((*it).ID==(*it_g).ID){
					check=1;
					break;
				}
                             
			}
			if(check==0){
                                // printf("ID_check=%d\n",(*it_g).ID);
				gID_list.push_back((*it_g));	// gID_list is always the same to all PEs.
			}
		}
		/* sort updated gID_list */
		//sort(gID_list.begin(),gID_list.end(),G_Info::before);

		newID_list.clear();
        nucleiID_list.clear();
         global_nucleiID_list.clear();
		return flag;	// equals to the # of newly nucleated grains
	}
}/* NucleationCheckDRX()*/

static int TrilinearInterpolationINT(double x0, double y0, double z0, int nx, int ny, int nz, double dx, double dy, double dz, int *V,
    double x, double y, double z)
{
  int i;
  int idx_x, idx_y, idx_z;
  int idx_px, idx_py, idx_pz;
  double sx, sy, sz;  // reduced coordinates in the "unit" cell
  /* distance to eight lattice points in the following order:
     [000,100,110,010,001,101,111,011] */
  double dist[8];
  int VV[8];
  int min_idx;
  double min_val;

  // assume slab decomposition along x-axis
  if((x>=x0)&&(x<=(x0+(nx-1)*dx))&&(y>=y0)&&(y<=(y0+(ny-1)*dy))&&(z>=z0)&&(z<=(z0+(nz-1)*dz))){
	
	  idx_x = (int)floor(x/dx);
	  idx_y = (int)floor(y/dy);
	  idx_z = (int)floor(z/dz);
	  sx = x/dx-(double)(idx_x);
	  sy = y/dy-(double)(idx_y);
	  sz = z/dz-(double)(idx_z);
    assert(idx_x>=0);
    assert(idx_y>=0);
    assert(idx_z>=0);

    idx_px = idx_x+1;
    if(idx_px>=nx) idx_px = nx-1;
    idx_py = idx_y+1;
    if(idx_py>=ny) idx_py = ny-1;
    idx_pz = idx_z+1;
    if(idx_pz>=nz) idx_pz = nz-1;

	
  /* distance to eight lattice points in the following order:
     [000,100,110,010,001,101,111,011] */
    dist[0] = sx*sx+sy*sy+sz*sz;
    dist[1] = (1-sx)*(1-sx)+sy*sy+sz*sz;
    dist[2] = (1-sx)*(1-sx)+(1-sy)*(1-sy)+sz*sz;
    dist[3] = sx*sx+(1-sy)*(1-sy)+sz*sz;
    dist[4] = sx*sx+sy*sy+(1-sz)*(1-sz);
    dist[5] = (1-sx)*(1-sx)+sy*sy+(1-sz)*(1-sz);
    dist[6] = (1-sx)*(1-sx)+(1-sy)*(1-sy)+(1-sz)*(1-sz);
    dist[7] = sx*sx+(1-sy)*(1-sy)+(1-sz)*(1-sz);

	  VV[0] = V[(idx_x*ny+idx_y)*nz+idx_z];
	  VV[1] = V[((idx_px)*ny+idx_y)*nz+idx_z];
	  VV[3] = V[(idx_x*ny+idx_py)*nz+idx_z];
	  VV[4]= V[(idx_x*ny+idx_y)*nz+idx_pz];
	  VV[7]= V[(idx_x*ny+idx_py)*nz+idx_pz];
	  VV[5]= V[((idx_px)*ny+idx_y)*nz+idx_pz];
	  VV[2]= V[((idx_px)*ny+idx_py)*nz+idx_z];
	  VV[6]= V[((idx_px)*ny+idx_py)*nz+idx_pz];

    min_idx=0;
    min_val = dist[0];
    for(i=1;i<8;i++){
      if(min_val>dist[i]){
        min_idx = i;
        min_val = dist[i];
      }
    }
    return VV[min_idx];
  }else{
    PError("Error in TrilinearInterpolation()",1203);
    return -1;
  }

}/*end TrilinearInterpolationINT()*/

static double TrilinearInterpolation(double x0, double y0, double z0, int nx, int ny, int nz, double dx, double dy, double dz, double *V,
    double x, double y, double z)
{
  int idx_x, idx_y, idx_z;
  int idx_px, idx_py, idx_pz;
  double sx, sy, sz;  // reduced coordinates in the "unit" cell
  double V_000, V_100, V_010, V_001, V_011, V_101, V_110, V_111;
  double V_xyz;

  // assume slab decomposition along x-axis
  if((x>=x0)&&(x<=(x0+(nx-1)*dx))&&(y>=y0)&&(y<=(y0+(ny-1)*dy))&&(z>=z0)&&(z<=(z0+(nz-1)*dz))){
	
	  idx_x = (int)floor(x/dx);
	  idx_y = (int)floor(y/dy);
	  idx_z = (int)floor(z/dz);
	  sx = x/dx-(double)(idx_x);
	  sy = y/dy-(double)(idx_y);
	  sz = z/dz-(double)(idx_z);
    assert(idx_x>=0);
    assert(idx_y>=0);
    assert(idx_z>=0);

    idx_px = idx_x+1;
    if(idx_px>=nx) idx_px = nx-1;
    idx_py = idx_y+1;
    if(idx_py>=ny) idx_py = ny-1;
    idx_pz = idx_z+1;
    if(idx_pz>=nz) idx_pz = nz-1;
	
	  V_000 = V[(idx_x*ny+idx_y)*nz+idx_z];
	  V_100 = V[((idx_px)*ny+idx_y)*nz+idx_z];
	  V_010 = V[(idx_x*ny+idx_py)*nz+idx_z];
	  V_001 = V[(idx_x*ny+idx_y)*nz+idx_pz];
	  V_011 = V[(idx_x*ny+idx_py)*nz+idx_pz];
	  V_101 = V[((idx_px)*ny+idx_y)*nz+idx_pz];
	  V_110 = V[((idx_px)*ny+idx_py)*nz+idx_z];
	  V_111 = V[((idx_px)*ny+idx_py)*nz+idx_pz];
	
	  V_xyz = V_000*(1.0-sx)*(1.0-sy)*(1.0-sz)+
	    V_100*sx*(1.0-sy)*(1.0-sz)+
	    V_010*sy*(1.0-sx)*(1.0-sz)+
	    V_001*sz*(1.0-sy)*(1.0-sx)+
	    V_101*sx*(1.0-sy)*sz+
	    V_011*sy*(1.0-sx)*sz+
	    V_110*sx*(1.0-sz)*sy+
	    V_111*sx*sy*sz;
	  return V_xyz;
  }else{
    PError("Error in TrilinearInterpolation()",1203);
    return -1.0;
  }

}/*end TrilinearInterpolation()*/

void InterpolateGID(void)
{
  int i,j,k,idx;
  int grid_ratio = CellDim_pf[0]/CellDim[0];
  int lnx_pf = lnx*grid_ratio;
  int lxs_pf = lxs*grid_ratio;
  real x0,y0,z0,dx,dy,dz,dx_pf,dy_pf,dz_pf;

  /* constants */
  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;
  dx = 1.0/(CellDim[0]-1);
  dy = 1.0/(CellDim[1]-1);
  dz = 1.0/(CellDim[2]-1);
  dx_pf = 1.0/(CellDim_pf[0]-1);
  dy_pf = 1.0/(CellDim_pf[1]-1);
  dz_pf = 1.0/(CellDim_pf[2]-1);

  AllocMem(table_grain_drx,Nxyz,int);
   AllocMem(table_pf_drx,Nxyz,int);
 // AllocMem(table_rho_drx,Nxyz,real);

  /* Interpolate grain_f */
	MPI_Allgather(grain_f,lsize,MPI_INT,table_grain_drx,lsize,MPI_INT,MPI_COMM_WORLD);
  for(i=0;i<lnx_pf;i++){
    for(j=0;j<CellDim_pf[1];j++){
      for(k=0;k<CellDim_pf[2];k++){
        idx = (i*CellDim_pf[1]+j)*CellDim_pf[2]+k;
        grain_f_drx[idx] = TrilinearInterpolationINT(x0,y0,z0,CellDim[0],CellDim[1],CellDim[2],dx,dy,dz,
            table_grain_drx, x0+(i+lxs_pf)*dx_pf,y0+j*dy_pf,z0+k*dz_pf);
      }
    }
  }
  
  	MPI_Allgather(first_pf,lsize,MPI_INT,table_pf_drx,lsize,MPI_INT,MPI_COMM_WORLD);
  for(i=0;i<lnx_pf;i++){
    for(j=0;j<CellDim_pf[1];j++){
      for(k=0;k<CellDim_pf[2];k++){
        idx = (i*CellDim_pf[1]+j)*CellDim_pf[2]+k;
        first_pf_drx[idx] = TrilinearInterpolationINT(x0,y0,z0,CellDim[0],CellDim[1],CellDim[2],dx,dy,dz,
            table_pf_drx, x0+(i+lxs_pf)*dx_pf,y0+j*dy_pf,z0+k*dz_pf);
      }
    }
  }


  free(table_grain_drx);
 // free(table_rho_drx);
  free(table_pf_drx);
	return;
}/*end InterpolateGID()*/

void InterpolateGrid(void)
{
  int i,j,k,idx;
  int grid_ratio = CellDim_pf[0]/CellDim[0];
  int lnx_pf = lnx*grid_ratio;
  int lxs_pf = lxs*grid_ratio;
  real x0,y0,z0,dx,dy,dz,dx_pf,dy_pf,dz_pf;

  /* constants */
  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;
  dx = 1.0/(CellDim[0]-1);
  dy = 1.0/(CellDim[1]-1);
  dz = 1.0/(CellDim[2]-1);
  dx_pf = 1.0/(CellDim_pf[0]-1);
  dy_pf = 1.0/(CellDim_pf[1]-1);
  dz_pf = 1.0/(CellDim_pf[2]-1);

  AllocMem(table_grain_drx,Nxyz,int);
  AllocMem(table_rho_drx,Nxyz,real);
  // AllocMem(table_dd_average_drx,Nxyz,real);

  /* Interpolate grain_f */
	MPI_Allgather(grain_f,lsize,MPI_INT,table_grain_drx,lsize,MPI_INT,MPI_COMM_WORLD);
  for(i=0;i<lnx_pf;i++){
    for(j=0;j<CellDim_pf[1];j++){
      for(k=0;k<CellDim_pf[2];k++){
        idx = (i*CellDim_pf[1]+j)*CellDim_pf[2]+k;
        grain_f_drx[idx] = TrilinearInterpolationINT(x0,y0,z0,CellDim[0],CellDim[1],CellDim[2],dx,dy,dz,
            table_grain_drx, x0+(i+lxs_pf)*dx_pf,y0+j*dy_pf,z0+k*dz_pf);
      }
    }
  }
  
  /*	MPI_Allgather(dd_index,lsize,MPI_INT,table_grain_drx,lsize,MPI_INT,MPI_COMM_WORLD);
  for(i=0;i<lnx_pf;i++){
    for(j=0;j<CellDim_pf[1];j++){
      for(k=0;k<CellDim_pf[2];k++){
        idx = (i*CellDim_pf[1]+j)*CellDim_pf[2]+k;
        dd_index_drx[idx] = TrilinearInterpolationINT(x0,y0,z0,CellDim[0],CellDim[1],CellDim[2],dx,dy,dz,
            table_grain_drx, x0+(i+lxs_pf)*dx_pf,y0+j*dy_pf,z0+k*dz_pf);
      }
    }
  }
*/  

  /* Interpolate rho_tot */
	MPI_Allgather(rho_tot,lsize,MPI_real,table_rho_drx,lsize,MPI_real,MPI_COMM_WORLD);
  for(i=0;i<lnx_pf;i++){
    for(j=0;j<CellDim_pf[1];j++){
      for(k=0;k<CellDim_pf[2];k++){
        idx = (i*CellDim_pf[1]+j)*CellDim_pf[2]+k;
        rho_tot_drx[idx] = TrilinearInterpolation(x0,y0,z0,CellDim[0],CellDim[1],CellDim[2],dx,dy,dz,
            table_rho_drx, x0+(i+lxs_pf)*dx_pf,y0+j*dy_pf,z0+k*dz_pf);
      }
    }
  }
  	//MPI_Allgather(dd_average,lsize,MPI_real,table_dd_average_drx,lsize,MPI_real,MPI_COMM_WORLD);
/*  for(i=0;i<lnx_pf;i++){
    for(j=0;j<CellDim_pf[1];j++){
      for(k=0;k<CellDim_pf[2];k++){
        idx = (i*CellDim_pf[1]+j)*CellDim_pf[2]+k;
        dd_average_drx[idx] = TrilinearInterpolation(x0,y0,z0,CellDim[0],CellDim[1],CellDim[2],dx,dy,dz,
            table_dd_average_drx, x0+(i+lxs_pf)*dx_pf,y0+j*dy_pf,z0+k*dz_pf);
      }
    }
  } */

  /* Interpolate gID_rex */
/*	MPI_Allgather(gID_rex,lsize,MPI_INT,table_grain_drx,lsize,MPI_INT,MPI_COMM_WORLD);
  for(i=0;i<lnx_pf;i++){
    for(j=0;j<CellDim_pf[1];j++){
      for(k=0;k<CellDim_pf[2];k++){
        idx = (i*CellDim_pf[1]+j)*CellDim_pf[2]+k;
        gID_rex_drx[idx] = TrilinearInterpolationINT(x0,y0,z0,CellDim[0],CellDim[1],CellDim[2],dx,dy,dz,
            table_grain_drx, x0+(i+lxs_pf)*dx_pf,y0+j*dy_pf,z0+k*dz_pf);
      }
    }
  } */

  free(table_grain_drx);
  free(table_rho_drx);
  //free(table_dd_average_drx);

	return;
}/*end InterpolateGrid()*/

real DrivingForceBoundaryMigration(void)
{
  real dd_local, dd_global;
  real dforce;

  Gradient_Rho();
#ifndef GB_BULGING_ONLY
  dd_local=0.0;
  dd_global=0.0;
  local_loop{
    if(dd_local<fabs(diff_rho[pIDX])){
      dd_local = fabs(diff_rho[pIDX]);
    }
  }
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&dd_local, &dd_global, 1, MPI_real,
			MPI_MAX, MPI_COMM_WORLD);
#else
  dd_local=0.0;
  dd_global=0.0;
  local_loop{
		if(GB_indicator[pIDX]>0){
      if(dd_local<fabs(diff_rho[pIDX])){
       dd_local = fabs(diff_rho[pIDX]);
      }
    }
  }
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&dd_local, &dd_global, 1, MPI_real,
			MPI_MAX, MPI_COMM_WORLD);
#endif
  dforce = Scale_Fdeform*dd_global*Shear_G[0]*bb_B[0][0]*bb_B[0][0]*1E-6; // MPa
  if(mpirank==0){
    printf("The driving force for the following PF is %f\n",dforce);
  }

  return dforce;
}/*end DrivingForceBoundaryMigration()*/

void Update_microstructure()
{
  
	std::vector<G_Info>::iterator it, end;
		std::vector<G_Info>::iterator it2, end2;
	static std::vector<G_Info> swapID_list;
	static std::vector<G_Info> startID_list;
	int tmp_flag,stop_check,start_check,size;
	int g_tot_temp,g_check_temp, g_tot, g_check;
	double g_ave;
	real ph,th,om;
	real ph1,th1,om1;
	int id,id1;
    voigt66 c066_local = {0.0};
    //int nph1, nph1_all;
  	voigt aux6;
	ten2nd aux33;
	//voigt66 caux66;
	//ten2nd sa2xt;	// transform matrix (sample -> xtal)
	//ten4th caux3333;
	if(DRX_Interpolation_Flag==1){
	  real x0,y0,z0,dx,dy,dz,dx_pf,dy_pf,dz_pf;
    if(mpirank==0){
      printf("Coarse-graining from PF grid to FFT grid...\n");
    }
	
	  /* constants */
	  x0 = 0.0;
	  y0 = 0.0;
	  z0 = 0.0;
	  dx = 1.0/(CellDim[0]-1);
	  dy = 1.0/(CellDim[1]-1);
	  dz = 1.0/(CellDim[2]-1);
	  dx_pf = 1.0/(CellDim_pf[0]-1);
	  dy_pf = 1.0/(CellDim_pf[1]-1);
	  dz_pf = 1.0/(CellDim_pf[2]-1);
    /* coarse-graining from PF grid to FFT grid */
    AllocMem(table_grain_drx,CellDim_pf[0]*CellDim_pf[1]*CellDim_pf[2],int);
    if(table_grain_drx==0){
      printf("rank=%d: Error in allocating table_grain_drx\n",mpirank);
    }

    /* Interpolate gID_new */
		MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Allgather(gID_new_drx,lsize_pf,MPI_INT,table_grain_drx,lsize_pf,MPI_INT,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
    local_loop{
      gID_new[pIDX] = TrilinearInterpolationINT(x0,y0,z0,CellDim_pf[0],CellDim_pf[1],CellDim_pf[2],dx_pf,dy_pf,dz_pf,
          table_grain_drx, x0+(px+lxs)*dx,y0+py*dy,z0+pz*dz);
	  }

    /* Interpolate grex_new */
/*		MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Allgather(grex_new_drx,lsize_pf,MPI_INT,table_grain_drx,lsize_pf,MPI_INT,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
    local_loop{
      grex_new[pIDX] = TrilinearInterpolationINT(x0,y0,z0,CellDim_pf[0],CellDim_pf[1],CellDim_pf[2],dx_pf,dy_pf,dz_pf,
          table_grain_drx, x0+(px+lxs)*dx,y0+py*dy,z0+pz*dz);
	  }
	    /* Interpolate gID_rex */
	/*	MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Allgather(gID_rex_drx,lsize_pf,MPI_INT,table_grain_drx,lsize_pf,MPI_INT,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
    local_loop{
      gID_rex[pIDX] = TrilinearInterpolationINT(x0,y0,z0,CellDim_pf[0],CellDim_pf[1],CellDim_pf[2],dx_pf,dy_pf,dz_pf,
          table_grain_drx, x0+(px+lxs)*dx,y0+py*dy,z0+pz*dz);
	  } */

    free(table_grain_drx);
    if(mpirank==0){
      printf("Coarse-graining is finished\n");
    }
  }

 //

	/* grain_f is the one used in FFT-CP
	   gID_new is the one updated from PF_XRX() */
	local_loop{
	/* Orientation of recrystallize grains

		   the main thing to update is the TranMat_xt2sa[][][]
		   check update_orient() in constituive.c for detail

		   if nothing is done, the orientation of the new grain
		   remains the same as that of the old grain.
		 
		 */


	grain_f_new[pIDX]=gID_new[pIDX];
	//printf("ID_list=%d\n",grain_f_new[pIDX]);

		/* update dislocation density */
			/* update dislocation density */
    if(gID_rex[pIDX]==1 && grain_f_new[pIDX]>=grain_f[pIDX]){
    //    growth[pIDX]=growth[pIDX]+1;
  //  }
//    if(gID_rex[pIDX]==1){
    	end=gID_list.end();
		//	tmp_flag=0;
			for(it=gID_list.begin();it!=end;it++){
                                 
				if(grain_f_new[pIDX]==(*it).ID){
				//	tmp_flag=1;
					ph1=(*it).t1;
			     	th1=(*it).Phi;
					om1=(*it).t2;
					id1=(*it).ID;
					// printf("ID_list=%d\n",id1);
				//	break;
		//		}
					//for(it_g=startD_list.begin();it_g!=end_g;it_g++){
		//	end2=startID_list.begin();
		//	int check2=0;
		//	for(it2=startID_list.end();it2!=end2;it2--){
			  // printf("ID=%d\n",(*it2).ID);
		//		if((*it).ID==(*it2).ID){
					//check2=1;
		//			break;
		//		}
		//	}
		//	if(it2==end2){
			    struct G_Info tmp_ID={id1,ph1,th1,om1};
                           
			startID_list.push_back(tmp_ID);
//                        printf("ID=%d\n",id1);
			//	startID_list.push_back((*it));	// gID_list is always the same to all PEs.
			}
			}
			}
			
			  if(gID_rex[pIDX]==1 && grain_f_new[pIDX]<grain_f[pIDX]){
                  
                  gID_rex[pIDX]=0;
                  growth[pIDX]=0;
              }
              
             
			//grain_f[pIDX]=grain_f_new[pIDX];
    }



//		if(mpirank==0){
   //     printf("gid=%d\n",id1);
  //  }	
    	
  //  swap_check[pIDX]=0;
  // if(mpirank==0){
     //   printf("swap=%d\n",id1);
  // }

//    }
    //if (flag_stop){
//        growth[pIDX]=0;
  //  }
   // if(grex_new[pIDX]==1){
		//	grain_f[pIDX]=gID_new[pIDX];
			//int jph=phase_f[pIDX];
		//	gID_rex[pIDX] = 1;
  //  }
    
  //  if(grex_new[pIDX] == 0) {
  //      gID_rex[pIDX] = 0;
  //  }
		//	if(swap_check[pIDX] >= 1) {


/*if(swap_check[pIDX]==1) {
    	end=gID_list.end();
		//	tmp_flag=0;
			for(it=gID_list.begin();it!=end;it++){
				if(grain_f[pIDX]==(*it).ID){
					tmp_flag=1;
					ph1=(*it).t1;
					th1=(*it).Phi;
					om1=(*it).t2;
					id1=(*it).ID;
					break;
				}
			}
//		if(mpirank==0){
   //     printf("gid=%d\n",id1);
  //  }	
    	struct G_Info tmp_ID={id1,ph1,th1,om1};
				startID_list.push_back(tmp_ID);
  //  swap_check[pIDX]=0;
  // if(mpirank==0){
     //   printf("swap=%d\n",id1);
  // }
}
if(swap_check[pIDX]>=1){
   swap_check[pIDX] = swap_check[pIDX]+1;
  
    // if(mpirank==0){
      //  printf("swap=%d\n",swap_check[pIDX]);
    //}
} */
			/* update the dislocation content in recrystallized grains */
		/*	if(swap_check[pIDX]==2){
			     if(mpirank==0){
        printf("swap2=%d\n",swap_check[pIDX]);
    }
			    start_check=grain_f[pIDX];
			 //   local_loop{
			        if(grain_f[pIDX]==start_check){
			            
			        
			    
			
			//if(growth[pIDX]<3){
			    int jph = phase_f[pIDX];
			for(int i=0;i<nSYS[jph-1];i++){
				rho_s[pIDX][i] = MAX(rho_s[pIDX][i]*ssdScaler_DRX, rho_SSD_initial);
				trial_rho_s[pIDX][i] = MAX(trial_rho_s[pIDX][i]*ssdScaler_DRX, rho_SSD_initial);
				rho_g1[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
				rho_g2[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
				rho_g3[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
				trial_gam_acum[pIDX][i] *= gndScaler_DRX;	// accumulated shear on each slip system, trial version
				trial_rho_g1[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
				trial_rho_g2[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
				trial_rho_g3[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
			}
			

      /* update the strain states. DRX grain should be undeformed. */

     //T2_loop {
     //Eps[pIDX][mi][mj] = 0.0;
     //   Edot[pIDX][mi][mj] = 0.0;
     //   DisGrad[pIDX][mi][mj] = 0.0;
    //    VelGrad[pIDX][mi][mj] = 0.0;
     
   // }
//	    for(int i=0;i<nSYS[jph-1];i++){
    //    gamdot[pIDX][i] = 0.0;
  //  }
      

//      /* moduli softening */
//      C6_loop{
//        C_gr[pIDX][mi][mj] *= 0.9;
//      }
     

			/* update the grain orientation in recrystallized grains */
		/*	end=gID_list.end();
			tmp_flag=0;
			for(it=gID_list.begin();it!=end;it++){
				if(grain_f[pIDX]==(*it).ID){
					tmp_flag=1;
					ph=(*it).t1;
					th=(*it).Phi;
					om=(*it).t2;
					break;
				}
			}*/
			//if(tmp_flag==0){
			//	PError("Wrong in Update_microstructure during DRX",102);
			//}

		//	ten2nd sa2xt;	// transform matrix (sample -> xtal)
		//	ten4th caux3333;
			//voigt aux6;
			//ten2nd aux33;
			//voigt66 caux66;
			//EulerToTransMatrix(&ph, &th, &om, sa2xt, 2);
		//	T2_loop{
		//		TranMat_xt2sa[pIDX][mi][mj] = sa2xt[mj][mi];
		//	}
		//	Ten4thTransform(caux3333,sa2xt,Cijkl[jph-1],2);	// Cijkl[jph] is in xtal ref. Do inverse transform
	//		chg_basis(aux6, aux33, caux66,caux3333,4);
	//		C6_loop{
	//		    if(C_gr[pIDX][mi][mj] < caux66[mi][mj]){
	//			C_gr[pIDX][mi][mj] = caux66[mi][mj];
	//		    }
//				c066_local[mi][mj] = caux66[mi][mj];
//}


//if(mpirank == 0){
  //  printf("orientation and stiffness tensor not changed after DRX");
//}

			/* update the t_last */
		//	t_last[pIDX] = TimeTot;

//#ifdef LOCAL_RELAX_FIRST

			/* update the local stress to make sure it satifies the
			   constitutive laws. This updated stress will be used to 
			   as an initial guess for stress-equilibrium solver*/
			/* voigt66 sg66;
			ten2nd xlambda_aux, sig_aux, eps_aux, strain_aux;
			voigt xlambda6, sig6, eps6, strain6;
			voigt edotp6;
			ten2nd edotp_aux;
			voigt66 d_edotp66;
			voigt tot_eps;	// total strain
			voigt sig6_old;
			ten4th aux3333;
			voigt66 aux66;
			real signorm, epsnorm;
			real conv_NR;
			real erroral, erral;
		
			voigt res;	// residual R to be nullified
			voigt66	jacob_inv;	// Jacobian of R
			int itmaxal, iterl;
			int i,j, k;

			C6_loop{
				sg66[mi][mj] = C_gr[pIDX][mi][mj];
			}
			LU_inv_66(sg66);
			
			T2_loop{  */
				/* the iteration starts with the current stress field */
				/* xlambda_aux[mi][mj] = Sig[pIDX][mi][mj]*sigScaler_DRX;
				sig_aux[mi][mj] = xlambda_aux[mi][mj];
				sig_aux[mi][mj] = 0.0; */
				/* plastic strain at time t. This is used in Eq. 4 to update total strain
				 from calculated strain rate*/
			//	eps_aux[mi][mj] = Eps[pIDX][mi][mj];
				/* DisGrad stores the updated displacement gradient obtained from Eq. 15.
                            So here strain_aux/strain6 stores the updated total strain*/
			/*	strain_aux[mi][mj] = (DisGrad[pIDX][mi][mj]+DisGrad[pIDX][mj][mi])/2.0;
			}
			chg_basis(xlambda6,xlambda_aux,aux66,aux3333,2);
			chg_basis(sig6,sig_aux,aux66,aux3333,2);
			chg_basis(eps6,eps_aux,aux66,aux3333,2);
			chg_basis(strain6,strain_aux,aux66,aux3333,2);

			signorm = 0.0;
			epsnorm = 0.0;
			T2_loop{
				signorm += xlambda_aux[mi][mj]*xlambda_aux[mi][mj];
				epsnorm += strain_aux[mi][mj]*strain_aux[mi][mj];
			}
			signorm = sqrt(signorm);
			epsnorm = sqrt(epsnorm); */
			/* Here based on the current accumulated plastic shear and
			   total strain which must be compatible to the current
			   strain field and hence unchanged upon DRX nucleation, we
			   use N-R method to solve Eq. 4.
			   The total strain is now fixed to the current value before DRX nucleation.*/

		/*	erroral = 1E-7;
			itmaxal = 100;
			iterl = 0;
			erral = 10*erroral; */

			/* Newton-Raphson method to solve augmented Lagrangians */
		/*	while((iterl<itmaxal)&&(fabs(erral)>fabs(erroral))){
				iterl++;

				for(i=0;i<6;i++) sig6_old[i] = sig6[i]; */

				/* Update the plastic strain rate based on current stress */
//#ifdef DD_POWER_LAW
				/* StrainRate_Orowan_POWER(sig6,edotp6,d_edotp66,pIDX,jph);
#else
				StrainRate_Orowan(sig6,edotp6,d_edotp66,pIDX,jph);
#endif
*/
				/* tot_eps is the total  strain and eps6 is the 
				   current plastic strain, i.e. Eq. 4 */
				/*for(i=0;i<6;i++){
					tot_eps[i] = eps6[i] + edotp6[i]*TimeStep;
					for(j=0;j<6;j++){
						tot_eps[i] += sg66[i][j]*sig6[j];
					}
				}

				// calculate the residual R, which is now simply Eq. 4
				for(i=0;i<6;i++){
					res[i] = sig6[i] - xlambda6[i];
					for(j=0;j<6;j++){
						res[i] += C066[i][j]*(tot_eps[j]-strain6[j]);
					}
				}

				// calculate the Jacobian of R
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
						jacob_inv[i][j] = (real)(i==j);
						for(k=0;k<6;k++){
							// Eq. 18
							jacob_inv[i][j] += C066[i][k]*(sg66[k][j]+d_edotp66[k][j]*TimeStep);
						}
					}
				}
				LU_inv_66(jacob_inv);

				// Newton-Raphson update
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
						sig6[i] -= jacob_inv[i][j]*res[j];
					}
				}

				// convergence check
				conv_NR= 0.0;
				for(i=0;i<6;i++){
					conv_NR += (sig6[i]-sig6_old[i])*(sig6[i]-sig6_old[i]);
				}
				erral = conv_NR/signorm;

			}	// end of while() loop

			chg_basis(sig6,sig_aux,aux66,aux3333,1);
			chg_basis(edotp6,edotp_aux,aux66,aux3333,1);
			// update stress and strain rate fields
			T2_loop{
				Sig[pIDX][mi][mj] = sig_aux[mi][mj];
				Edot[pIDX][mi][mj] = edotp_aux[mi][mj];
			}


			if(mpirank==0){
				printf("pyz: Relax local stress after DRX nucleation: iter = %d\n", iterl);
//				printf("pyz: relaxed stress tensor:\n");
//				PrintTensor(Sig[pIDX]);
//				printf("pyz: relaxed shear rate tensor:\n");
//				PrintTensor(Edot[pIDX]);
			}
#elif defined(LOCAL_RELAX_SECOND)
			voigt66 sg66;
			ten2nd sig_aux, eps_aux, strain_aux;
			voigt sig6, eps6, strain6;
			voigt sig6_old;
			voigt edotp6;
			ten2nd edotp_aux;
			voigt66 d_edotp66;
			ten4th aux3333;
			voigt66 aux66;
			real signorm, epsnorm;
			real conv_NR;
			real erroral, erral;
		
			voigt res;	// residual R to be nullified
			voigt66	jacob_inv;	// Jacobian of R
			int itmaxal, iterl;
			int i,j, k;

			C6_loop{
				sg66[mi][mj] = C_gr[pIDX][mi][mj];
			}
			LU_inv_66(sg66);
			
			T2_loop{
				/* the iteration starts with the current stress field */
			/*	sig_aux[mi][mj] = Sig[pIDX][mi][mj]*sigScaler_DRX;
				/* plastic strain */
			/*	eps_aux[mi][mj] = Eps[pIDX][mi][mj];
				/* keep the current total strain*/
			/*	strain_aux[mi][mj] = (DisGrad[pIDX][mi][mj]+DisGrad[pIDX][mj][mi])/2.0;
			}
			chg_basis(sig6,sig_aux,aux66,aux3333,2);
			chg_basis(strain6,strain_aux,aux66,aux3333,2);
			chg_basis(eps6,eps_aux,aux66,aux3333,2);

			erroral = 1E-7;
			itmaxal = 100;
			iterl = 0;
			erral = 10*erroral;

			/* Newton-Raphson method to solve augmented Lagrangians */
		/*	while((iterl<itmaxal)&&(fabs(erral)>fabs(erroral))){
				iterl++;

				for(i=0;i<6;i++) sig6_old[i] = sig6[i];
				signorm = 0.0;
				for(i=0;i<6;i++){
					signorm += sig6_old[i]*sig6_old[i];
				}
				signorm = sqrt(signorm);

				/* Update the plastic strain rate based on current stress */
/*#ifdef DD_POWER_LAW
				StrainRate_Orowan_POWER(sig6,edotp6,d_edotp66,pIDX,jph);
#else
				StrainRate_Orowan(sig6,edotp6,d_edotp66,pIDX,jph);
#endif

				// calculate the residual R, which is now simply Eq. 4
				for(i=0;i<6;i++){
					res[i] = strain6[i] - eps6[i] - edotp6[i]*TimeStep;
					for(j=0;j<6;j++){
						res[i] -= sg66[i][j]*sig6[j];
					}
				}

				// calculate the Jacobian of R
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
						jacob_inv[i][j] = -1.0*sg66[i][j]-1.0*d_edotp66[i][j]*TimeStep;
					}
				}
				LU_inv_66(jacob_inv);

				// Newton-Raphson update
				for(i=0;i<6;i++){
					for(j=0;j<6;j++){
						sig6[i] -= jacob_inv[i][j]*res[j];
					}
				}

				// convergence check
				conv_NR= 0.0;
				for(i=0;i<6;i++){
					conv_NR += (sig6[i]-sig6_old[i])*(sig6[i]-sig6_old[i]);
				}
				erral = conv_NR/signorm;

			}	// end of while() loop

			chg_basis(sig6,sig_aux,aux66,aux3333,1);
			chg_basis(edotp6,edotp_aux,aux66,aux3333,1);
			// update stress and strain rate fields
			T2_loop{
				Sig[pIDX][mi][mj] = sig_aux[mi][mj];
				Edot[pIDX][mi][mj] = edotp_aux[mi][mj];
			}


			if(mpirank==0){
				printf("pyz: Relax local stress after DRX nucleation: iter = %d\n", iterl);
				printf("pyz: relaxed stress tensor:\n");
				PrintTensor(Sig[pIDX]);
				printf("pyz: relaxed shear rate tensor:\n");
				PrintTensor(Edot[pIDX]);
			}
#else
			T2_loop{
				/* we treat the DRX as a phase-transformation that
				   leads to a plastic-strain-free state */
			/*	   Sig[pIDX][mi][mj] *= sigScaler_DRX;
			}
#endif

		}
			   // }
			} 
//Ten4thTransform(caux3333,sa2xt,Cijkl[jph-1],2);	// Cijkl[jph] is in xtal ref. Do inverse transform
//chg_basis(aux6, aux33, caux66,caux3333,4);
//if(grex_new[pIDX] != 1){
//C6_loop{
//   c066_local[mi][mj] += caux66[mi][mj];
//}
//}

/*if(swap_check[pIDX]>5) {
    	end=gID_list.end();
		//	tmp_flag=0;
			for(it=gID_list.begin();it!=end;it++){
				if(grain_f[pIDX]==(*it).ID){
					tmp_flag=1;
					ph=(*it).t1;
					th=(*it).Phi;
					om=(*it).t2;
					id=grain_f[pIDX];
					break;
				}
			}
    	struct G_Info tmp_ID={id,ph,th,om};
				swapID_list.push_back(tmp_ID);
    swap_check[pIDX]=0;
 //   if(mpirank==0){
   //     printf("swap3=%d\n",swap_check[pIDX]);
//    }
    //stop_check = grain_f[pIDX];
    //local_loop{
  //      if(grain_f[pIDX] == stop_check) {
    //growth[pIDX] = 0;
//    gID_rex[pIDX] = 0;
    //    }
  //  }
//}
	} */

//C6_loop{
//		MPI_Barrier(MPI_COMM_WORLD);
//		MPI_Allreduce(&c066_local[mi][mj], &C066[mi][mj], 1, MPI_real,
//				MPI_SUM, MPI_COMM_WORLD);
//		C066[mi][mj] *= WGT;
//	}
//	printf("Value of C066[1][1] is %le\n",C066[0][0]);

//	C6_loop{
//		S066[mi][mj] = C066[mi][mj];
//	}

//	LU_inv_66(S066);

//	chg_basis(aux6,aux33,C066,C0,3);
//	chg_basis(aux6,aux33,S066,S0,3);	
//}
//exit(0);
 
//	int local_size = swapID_list.size();
	/* construct global swapID_list */
/*		int global_size;
		MPI_Allreduce(&local_size, &global_size, 1,
				MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		std::vector<G_Info> global_swapID_list(global_size);	// only store new IDs
		std::vector<int> nsize(NumPE);
		MPI_Allgather(&local_size,1,MPI_INT,&nsize[0],1,MPI_INT,MPI_COMM_WORLD);
		std::vector<int> disps(NumPE);
		for(int i=0;i<NumPE;i++){
			disps[i]=(i>0)?(disps[i-1]+nsize[i-1]):0;
		}
	
		// define G_ID_Entry
		struct G_Info value = {0};
		MPI_Datatype G_ID_Entry;
		int count=4;
		int block_lens[4]={1,1,1,1};
		MPI_Aint indices[4];
		MPI_Datatype old_types[4]={MPI_INT, MPI_real, MPI_real, MPI_real};
		MPI_Address(&value,&indices[0]);
		MPI_Address(&value.t1,&indices[1]);
		MPI_Address(&value.Phi,&indices[2]);
		MPI_Address(&value.t2,&indices[3]);
		// make relative
		for(int i=count-1;i>=0;i--)
			indices[i] -= indices[0];
		MPI_Type_struct(count,block_lens,indices,old_types,&G_ID_Entry);
		MPI_Type_commit(&G_ID_Entry);
	
		MPI_Allgatherv(&swapID_list[0],local_size,G_ID_Entry,
				&global_swapID_list[0],&nsize[0],&disps[0],G_ID_Entry,MPI_COMM_WORLD);
		std::vector<G_Info>::iterator it_g, end_g;
		end_g=global_swapID_list.end();
		//for(it_g=global_swapID_list.begin();it_g!=end_g;it_g++){
		//	end=swapID_list.begin();
		//	int check=0;
		//	for(it=swapID_list.end();it!=end;it--){
		//		if((*it).ID==(*it_g).ID){
		//			check=1;
		//			break;
		//		}
		//	}
		//	if(check==0){
		//		swapID_list.push_back((*it_g));	// gID_list is always the same to all PEs.
		//	}
		//}
		//	for(it=swapID_list.end();it!=end;it--){
//		for(it_g=global_swapID_list.begin();it_g!=end_g;it_g++){
//			    printf("SwapID_list=%d\n",(*it_g).ID);
//			}
		//	int local_size2 = startID_list.size(); 
	/* construct global swapID_list */
	int local_size2 = startID_list.size();
	int global_size2;
		MPI_Allreduce(&local_size2, &global_size2, 1,
				MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		std::vector<G_Info> global_startID_list(global_size2);	// only store new IDs
		std::vector<int> nsize1(NumPE);
		MPI_Allgather(&local_size2,1,MPI_INT,&nsize1[0],1,MPI_INT,MPI_COMM_WORLD);
		std::vector<int> disps1(NumPE);
		for(int i=0;i<NumPE;i++){
			disps1[i]=(i>0)?(disps1[i-1]+nsize1[i-1]):0;
		}
	
		// define G_ID_Entry
		struct G_Info value1 = {0};
		MPI_Datatype G_ID_Entry1;
		int count1=4;
		int block_lens1[4]={1,1,1,1};
		MPI_Aint indices1[4];
		MPI_Datatype old_types1[4]={MPI_INT, MPI_real, MPI_real, MPI_real};
		MPI_Address(&value1,&indices1[0]);
		MPI_Address(&value1.t1,&indices1[1]);
		MPI_Address(&value1.Phi,&indices1[2]);
		MPI_Address(&value1.t2,&indices1[3]);
		// make relative
		for(int i=count1-1;i>=0;i--)
			indices1[i] -= indices1[0];
		MPI_Type_struct(count1,block_lens1,indices1,old_types1,&G_ID_Entry1);
		MPI_Type_commit(&G_ID_Entry1);
	
		MPI_Allgatherv(&startID_list[0],local_size2,G_ID_Entry1,
				&global_startID_list[0],&nsize1[0],&disps1[0],G_ID_Entry1,MPI_COMM_WORLD);
      // end=global_startID_list.end();
//for(it=global_startID_list.begin();it!=end;it++){
	//	    printf("StartID_list=%d\n",(*it).ID);
	//		}
          
	//std::vector<G_Info>::iterator it_g1, end_g1;
	//	end_g1=global_startID_list.end();
	//for(it_g1=global_startID_list.begin();it_g1!=end_g1;it_g1++){
	//		end=startID_list.begin();
	//		int check=0;
	//		for(it=startID_list.end();it!=end;it--){
	//			if((*it).ID==(*it_g).ID){
	//				check=1;
	//				break;
	//			}
	//		}
	//		if(check==0){
	//			startID_list.push_back((*it_g));	// gID_list is always the same to all PEs.
	//		}
		//	for(it=startID_list.end();it!=end;it--){
		//	    printf("StartID_list=%d\n",(*it_g1).ID);
	//	}
//exit(0);
	//	}

    	end=global_startID_list.end();
			//tmp_flag=0;
			for(it=global_startID_list.begin();it!=end;it++){
                         // printf("StartID_list=%d\n",(*it).ID);
                local_loop{
                    

				if(grain_f_new[pIDX]==(*it).ID){
					gID_rex[pIDX]=1;
                  
					/*if(mpirank==0){
					  printf("growth= %d %d\n",growth[pIDX], (*it).ID);
					}*/
				}
			}
}
 
 

local_loop{
    if(gID_rex[pIDX]==1){
        growth[pIDX]=growth[pIDX]+1;
// printf("growth= %d %d %d\n",growth[pIDX], grain_f_new[pIDX],grain_f[pIDX]);
}
}
//exit(0);
	 	end=gID_list.end();
	for(it=gID_list.begin();it!=end;it++){
	  g_tot_temp = 0;
	  g_check_temp=0;
		g_ave=0.0;
                g_tot=0;
        g_check=0;
	  	local_loop{
		    if(grain_f_new[pIDX]==(*it).ID){
		        g_tot_temp ++;
		        if(growth[pIDX]==1){
		        g_check_temp ++ ;
               // printf("checkid2=%d\n",(*it).ID);
		        }
		    }
		   
		}
	    
            
	    	MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&g_tot_temp, &g_tot, 1, MPI_INT,
				MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&g_check_temp, &g_check, 1, MPI_INT,
				MPI_SUM, MPI_COMM_WORLD);
                                      //  printf("%d %d %d\n",g_tot,g_check,(*it).ID);
                                        if(g_tot!=0){
					g_ave= double ((g_check*100)/g_tot);
                                        }
                                       // exit(0);
      		 if (g_ave < 0.001){
                  //    printf("check = %d\n",(*it).ID);

					local_loop{
				
		   if(grain_f_new[pIDX]==(*it).ID){
		        growth[pIDX]=0;
		        gID_rex[pIDX]=0;
		        //  printf("checkid2=%d\n",(*it).ID);  
			//		dd_index++;
	}
					}
					
									
	}
	
	//exit(0);
        }

          // size =  ((lnx*CellDim[1] + CellDim[1])*CellDim[2] + CellDim[2]);
	  //  printf("file_write_start_fft");
                 /* local_loop{

               printf("gid_rex=%d\n",gID_rex[pIDX]);

                                   }*/
				         /*	MPI_File fp;
	        MPI_Status status;
	        char fname[100]={0};
            sprintf(fname,"nuclei_%06d.iout",istep);
	        MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	        MPI_File_set_view(fp, mpirank*lsize*sizeof(int), MPI_INT,
			MPI_INT, "native", MPI_INFO_NULL);
		MPI_File_write(fp,gID_rex,lsize, MPI_INT, &status);
				MPI_File_close(&fp);*/
	    				   //  }
	//}
        //exit(0);

 local_loop{
//     
      grain_f[pIDX]=grain_f_new[pIDX];
}
//     	//end=global_startID_list.end();
//     	//	end=startID_list.end();
//     	//	exit(0);
// 			//tmp_flag=0;
// 		//	for(it=startID_list.begin();it!=end;it++){
// 		//		if(grain_f[pIDX]==(*it).ID){
// 		if(growth[pIDX]==1 && rho_tot[pIDX] > 50.0 ){
// 					 int jph = phase_f[pIDX];
// 			for(int i=0;i<nSYS[jph-1];i++){
// 				rho_s[pIDX][i] = rho_s[pIDX][i]*ssdScaler_DRX;
// 				trial_rho_s[pIDX][i] = rho_s[pIDX][i]*ssdScaler_DRX;
// 				rho_g1[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 				rho_g2[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 				rho_g3[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 				trial_gam_acum[pIDX][i] *=  gndScaler_DRX;	// accumulated shear on each slip system, trial version
// 				trial_rho_g1[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 				trial_rho_g2[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 				trial_rho_g3[pIDX][i] *= gndScaler_DRX;	// GND is initially zero
// 			}
// 				T2_loop{
// 				/* we treat the DRX as a phase-transformation that
// 				   leads to a plastic-strain-free state */
// 				   Sig[pIDX][mi][mj] *= sigScaler_DRX;
// // Eps[pIDX][mi][mj] *= sigScaler_DRX;
//                                   
//      
// 			}
// 			
// 		
// 					/* if(mpirank==0){
// 					    printf("checkid1=%d",(*it).ID);
// 					} */
// #ifdef LOCAL_RELAX_FIRST
// 
// 			/* update the local stress to make sure it satifies the
// 			   constitutive laws. This updated stress will be used to 
// 			   as an initial guess for stress-equilibrium solver*/
// 			voigt66 sg66;
// 			ten2nd xlambda_aux, sig_aux, eps_aux, strain_aux;
// 			voigt xlambda6, sig6, eps6, strain6;
// 			voigt edotp6;
// 			ten2nd edotp_aux;
// 			voigt66 d_edotp66;
// 			voigt tot_eps;	// total strain
// 			voigt sig6_old;
// 			ten4th aux3333;
// 			voigt66 aux66;
// 			real signorm, epsnorm;
// 			real conv_NR,  conv_istep_e, conv_istep_s;
// 			real erroral, erral;
// 		
// 			voigt res;	// residual R to be nullified
// 			voigt66	jacob_inv;	// Jacobian of R
// 			int itmaxal, iterl;
// 			int i,j, k;
// 
// 			C6_loop{
// 				sg66[mi][mj] = C_gr[pIDX][mi][mj];
// 			}
// 			LU_inv_66(sg66);
// 			
// 			T2_loop{
// 				/* the iteration starts with the current stress field */
// 				xlambda_aux[mi][mj] = Sig[pIDX][mi][mj];
// 				sig_aux[mi][mj] = xlambda_aux[mi][mj];
// 				//sig_aux[mi][mj] = 0.0;
// 				/* plastic strain at time t. This is used in Eq. 4 to update total strain
// 				 from calculated strain rate*/
// 				eps_aux[mi][mj] = Eps[pIDX][mi][mj];
// 				/* DisGrad stores the updated displacement gradient obtained from Eq. 15.
// 				   So here strain_aux/strain6 stores the updated total strain*/
// 				strain_aux[mi][mj] = (DisGrad[pIDX][mi][mj]+DisGrad[pIDX][mj][mi])/2.0;
// 			}
// 			chg_basis(xlambda6,xlambda_aux,aux66,aux3333,2);
// 			chg_basis(sig6,sig_aux,aux66,aux3333,2);
// 			chg_basis(eps6,eps_aux,aux66,aux3333,2);
// 			chg_basis(strain6,strain_aux,aux66,aux3333,2);
// 
// 			signorm = 0.0;
// 			epsnorm = 0.0;
// 			T2_loop{
// 				signorm += xlambda_aux[mi][mj]*xlambda_aux[mi][mj];
// 			//	epsnorm += strain_aux[mi][mj]*strain_aux[mi][mj];
// 			}
// 			signorm = sqrt(signorm);
// 			//epsnorm = sqrt(epsnorm);
// 			/* Here based on the current accumulated plastic shear and
// 			   total strain which must be compatible to the current
// 			   strain field and hence unchanged upon DRX nucleation, we
// 			   use N-R method to solve Eq. 4.
// 			   The total strain is now fixed to the current value before DRX nucleation.*/
// 
// 			erroral = 1E-7;
// 			itmaxal = 100;
// 			iterl = 0;
// 			erral = 10*erroral;
// 
// 			/* Newton-Raphson method to solve augmented Lagrangians */
// 			while((iterl<itmaxal)&&(fabs(erral)>fabs(erroral))){
// 				iterl++;
// 
// 				for(i=0;i<6;i++) sig6_old[i] = sig6[i];
// 
// 				/* Update the plastic strain rate based on current stress */
// #ifdef DD_POWER_LAW
// 				StrainRate_Orowan_POWER(sig6,edotp6,d_edotp66,pIDX,jph);
// #else
// 				StrainRate_Orowan(sig6,edotp6,d_edotp66,pIDX,jph);
// #endif
// 
// 				/* tot_eps is the total  strain and eps6 is the 
// 				   current plastic strain, i.e. Eq. 4 */
// 				for(i=0;i<6;i++){
// 					tot_eps[i] = eps6[i] + edotp6[i]*TimeStep;
// 					for(j=0;j<6;j++){
// 						tot_eps[i] += sg66[i][j]*sig6[j];
// 					}
// 				}
// 
// 				// calculate the residual R, which is now simply Eq. 4
// 				for(i=0;i<6;i++){
// 					res[i] = sig6[i] - xlambda6[i];
// 					for(j=0;j<6;j++){
// 						res[i] += C066[i][j]*(tot_eps[j]-strain6[j]);
// 					}
// 				}
// 
// 				// calculate the Jacobian of R
// 				for(i=0;i<6;i++){
// 					for(j=0;j<6;j++){
// 						jacob_inv[i][j] = (real)(i==j);
// 						for(k=0;k<6;k++){
// 							// Eq. 18
// 							jacob_inv[i][j] += C066[i][k]*(sg66[k][j]+d_edotp66[k][j]*TimeStep);
// 						}
// 					}
// 				}
// 				LU_inv_66(jacob_inv);
// 
// 
// 				// Newton-Raphson update
// 				for(i=0;i<6;i++){
// 					for(j=0;j<6;j++){
// 						sig6[i] -= jacob_inv[i][j]*res[j];
// 					}
// 				}
// 
// 
// 				// convergence check
// 				conv_NR= 0.0;
//                                  conv_istep_e= 0.0;
// 				conv_istep_s = 0.0;
// 				for(i=0;i<6;i++){
// 					conv_NR += (sig6[i]-sig6_old[i])*(sig6[i]-sig6_old[i]);
//                                         conv_istep_e += (tot_eps[i]-strain6[i])*(tot_eps[i]-strain6[i]);
// 					conv_istep_s += (sig6[i]-xlambda6[i])*(sig6[i]-xlambda6[i]);
//                                 
// 				}
// 				erral = conv_NR/signorm;
// 
// 			}	// end of while() loop
// 
// 			chg_basis(sig6,sig_aux,aux66,aux3333,1);
// 			chg_basis(edotp6,edotp_aux,aux66,aux3333,1);
// 			// update stress and strain rate fields
// 			T2_loop{
// 				Sig[pIDX][mi][mj] = sig_aux[mi][mj];
// 				Edot[pIDX][mi][mj] = edotp_aux[mi][mj];
// 			}
//                       //  *Err_s_local += conv_istep_s;
// 			//*Err_e_local += conv_istep_e;
//                  
// #endif
// 
// //			if(mpirank==0){
// //				printf("pyz: Relax local stress after DRX nucleation: iter = %d\n", iterl);
// //				printf("pyz: relaxed stress tensor:\n");
// //				PrintTensor(Sig[pIDX]);
// //				printf("pyz: relaxed shear rate tensor:\n");
// //				PrintTensor(Edot[pIDX]);
// //			}
// 				}
// 			}
			
//}
//exit(0);

	/* update the neighboring gID info */
	Gradient_ExchangeGrainID();


 
 //	swapID_list.clear();
 		startID_list.clear();
 //		 	global_swapID_list.clear();
 	global_startID_list.clear();
	return;
}/*end Update_microstructure()*/

#endif
