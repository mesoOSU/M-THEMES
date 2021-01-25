
#include "evp.h"

using namespace std;

double laplacian(const double& right, const double& left, const double& front, const double& back, const double& up,const double& down,const double& right_up, const double& right_down, const double& right_front, const double& right_back, const double& left_up,const double& left_down,const double& left_front, const double& left_back, const double& center_up_front, const double& center_up_back, const double& center_down_front,const double& center_down_back, const double& c)
{
	return (0.5*(right+left+front+back+up+down)+0.125*(left_up+left_front+left_down+left_back+right_down+right_up+right_front+right_back+center_up_front+center_down_back+center_up_back+center_down_front)-4.5*c); 
}/*end laplacian()*/

void Update_gID(PFgrid3D& gd_rex, int *gID, int L_x, int L_y, int L_z, int lnx, int lny, int lnz)
{
	int i, N, idx_rex;
	double val, max_val_rex;
	for(int px=3;px<lnx;px++){
		for(int py=3;py<lny;py++){
			for(int pz=3;pz<lnz;pz++){
				int pIDX_d = (((px-3)*L_y+py-3)*L_z+pz-3);	// index in original data, excluding extra 4 layers
				PFdata& data_rex = gd_rex[px][py][pz];
			//	PFdata& data_def = gd_def[px][py][pz];
		
				N = data_rex.nonzero();
				idx_rex = 0;
				max_val_rex = 0.0;
				for(i=0;i<N;i++){
					val = data_rex.value(i);
				
					if(max_val_rex < val){
						max_val_rex = val;
						idx_rex = data_rex.index(i);
					}
				}
		
			/*	N = data_def.nonzero();
				idx_def = 0;
				max_val_def = 0.0;
				for(i=0;i<N;i++){
					val = data_def.value(i);
					if(max_val_def < val){
						max_val_def = val;
						idx_def = data_def.index(i);
					}
				}  */
		       
				//if(max_val_rex >= 0.25) {
					gID[pIDX_d] = idx_rex;
			//		if(flag_rex==1){
					   // flag_rex==0;
					//    gID_rex[pIDX_d]=1;
					//}
				
                             //   }
				    
					
			//		if(flag_rex==1){
					    //flag_rex==0;
					  //  gID_rex[pIDX_d]=0;
					//}
				
				}
			}
		}
	

	return;
}/*end Update_gID()*/

void Update_gID_final(PFgrid3D& gd_rex, int *gID, int L_x, int L_y, int L_z, int lnx, int lny, int lnz,int mpirank)
{
	int i, N, idx_rex,c;
	double val, max_val_rex;
	for(int px=3;px<lnx;px++){
		for(int py=3;py<lny;py++){
			for(int pz=3;pz<lnz;pz++){
				int pIDX_d = (((px-3)*L_y+py-3)*L_z+pz-3);	// index in original data, excluding extra 4 layers
				PFdata& data_rex = gd_rex[px][py][pz];
			//	PFdata& data_def = gd_def[px][py][pz];
		  
		  //   c=0;
		   
				N = data_rex.nonzero();
				idx_rex = 0;
				max_val_rex = 0.0;
				for(i=0;i<N;i++){
					val = data_rex.value(i);
			 /* if (mpirank==0){
					if(data_rex.index(i) == 25) {
					   printf("%d %d %d %d %le\n", px-2, py-2, pz-2, pIDX_d, val);
					   c=1;
					}
	   					 } */		
					if(max_val_rex < val){
						max_val_rex = val;
						idx_rex = data_rex.index(i);
					}
				}
		/* if (c==0){
					     printf("%d %d %d %d %le\n", px-2, py-2, pz-2, pIDX_d, 0.0);
					}*/
			/*	N = data_def.nonzero();
				idx_def = 0;
				max_val_def = 0.0;
				for(i=0;i<N;i++){
					val = data_def.value(i);
					if(max_val_def < val){
						max_val_def = val;
						idx_def = data_def.index(i);
					}
				}  */
		
			//	if(grex_new[pIDX_d] == 1){
				  //  if(max_val_rex >= 0.25) {
					gID[pIDX_d] = idx_rex;
                                        if(mpirank==0) {
                                      //   printf("%d\n",idx_rex);   
                                        }
				//	if(flag_rex==1){
					   // flag_rex==0;
				//	    gID_rex[pIDX_d]=1;
				//	}
				//}
				
				  /*  if(max_val_def) {
					gID[pIDX_d] = idx_def;
				//	if(flag_rex==1){
					    //flag_rex==0;
					    gID_rex[pIDX_d]=0;
				//	}
				} */
				}
			}
		}
	

	return;
}/*end Update_gID_final()*/

void clamping_eta(PFgrid3D& gd_rex,int lnx, int lny, int lnz)
{
	int N,N1,j, i,count;
	double sum_val,s;

	for(int px=3;px<lnx;px++){
		for(int py=3;py<lny;py++){
			for(int pz=3;pz<lnz;pz++){
				sum_val = 0.0;
                              //  count =0;
				N = gd_rex[px][py][pz].nonzero();
                                
                                for(i=0;i<N;i++){
				    if(gd_rex[px][py][pz].value(i)>1.0) {
                                        
                        
                                       // count++;
				      gd_rex[px][py][pz].value(i) = 1.0;
				      for(j=0;j<N;j++) {
				          if(j!=i){
				              gd_rex[px][py][pz].value(j) = 0.0;
				          }
				      }
				      
				      
                                    }
                                }
                                
                            //  printf("count = %d\n",count);  
                                
                                
                                
                           // sum_val = 0.0;
				/*count = 0;
				N = gd_rex[px][py][pz].nonzero();
					//N1 = gd_def[px][py][pz].nonzero();
				 for(i=0;i<N;i++){
				    if(gd_rex[px][py][pz].value(i)<0.0) {
				       s = fabs(gd_rex[px][py][pz].value(i));
				        gd_rex[px][py][pz].value(i) = 0.0;
				        for(j=0;j<N;j++) {
				            if(j!=i && gd_rex[px][py][pz].value(j)>0.0 && gd_rex[px][py][pz].value(j)<1.0) {
				                count++;
                                               // if(mpirank==0) {
                                                
                                                    
                 //   printf("count = %d\n", count);
                                                
                                               // }
				            }
				        }
				      
				     for(j=0;j<N;j++) {
				       if(j!=i && gd_rex[px][py][pz].value(j)>0.0 && gd_rex[px][py][pz].value(j)<1.0) {
				           gd_rex[px][py][pz].value(j)  -=  s/count; 
 				    }
				     }
 				    
				    }
				   
					//sum_val += gd_rex[px][py][pz].value(i);
				}  */
				
                               // N = gd_rex[px][py][pz].nonzero();
// 				for(i=0;i<N;i++) {
//                                 if(gd_rex[px][py][pz].value(i) < 0.0) {
//                                     for(j =0;j<N;j++) {
//                                         if(j != i) {
                              
				
				
				N = gd_rex[px][py][pz].nonzero();
				for(i=0;i<N;i++){
                    
					sum_val += gd_rex[px][py][pz].value(i);
				}
				/*if(sum_val < 0.99) {
				printf("%le\n",sum_val);
                                } */
			//	sum_val = sqrt(sum_val);
		
			/*	if(sum_val> 0.0){
					N = gd_rex[px][py][pz].nonzero();
					for(i=0;i<N;i++){
						gd_rex[px][py][pz].value(i) /= sum_val;
					}
				/*	N = gd_def[px][py][pz].nonzero();
					for(i=0;i<N;i++){
						gd_def[px][py][pz].value(i) /= sum_val;
					} */
				} 
                        }
                }
                
        return;
}
			
			 /*   if(sum_val>1.0){
					N = gd_rex[px][py][pz].nonzero();
					//N1 = gd_def[px][py][pz].nonzero();
					for(i=0;i<N;i++){
						gd_rex[px][py][pz].value(i) += (1.0-sum_val)/(N+N1);
					}
				//	N = gd_def[px][py][pz].nonzero();
					for(i=0;i<N1;i++){
						gd_def[px][py][pz].value(i) += (1.0-sum_val)/(N+N1);
					}
				} 
			}
		}
	}

	return;
}/*end clamping_eta()*/

void boundary_conditions(PFgrid3D& gd, int mpirank, int NumPE, int l_x, int l_y, int l_z, int flag)
{
  int x, y, z, N, index, j;
  int next,front;
  double value;
  int count;
  int node_count;
  int sendcutcount;
  int *sendcut;
  int *recvcut;
  int *sendgid;
  int *recvgid;
  double *sendrho;
  double *recvrho;
  int recvcount;
  int recv_node_count;
  MPI_Status status;

  // segment size on each processor
  MPI_Datatype vectortype, oldtypes[2];
  int blockcounts[2];
  MPI_Aint offsets[2], extent;

  typedef struct { 
  double value; 
  int  index; 
  } buff_data;

  buff_data *senddata, *recvdata;

  next=mpirank+1;if(next==NumPE)next=0; 
  front=mpirank-1;if(front==-1)front=NumPE-1;
  
 /* Setup description vector*/ 
  offsets[0] = 0; 
  oldtypes[0] = MPI_DOUBLE; 
  blockcounts[0] = 1; 
  MPI_Type_extent(MPI_DOUBLE, &extent); 
  offsets[1] = 1 * extent; 
  oldtypes[1] = MPI_INT; 
  blockcounts[1] = 1; 

  MPI_Type_struct(2, blockcounts, offsets, oldtypes, &vectortype); 
  MPI_Type_commit(&vectortype); 
  /*End  description vector*/ 

//if(flag<3){
  // The following applies periodic boundary condition along y-axis
  for (x=3;x<l_x-3; x++) {
    for (z=3; z<l_z-3; z++) {
        
         y = l_y-6;
      PFdata& data1 = gd[x][y][z];
      N = data1.nonzero();
      for (j=0; j<N; j++) {
        value = data1.value(j);
		index = data1.index(j);
		gd[x][0][z][index] = value;
      }
         y = l_y-5;
      PFdata& data2 = gd[x][y][z];
      N = data2.nonzero();
      for (j=0; j<N; j++) {
        value = data2.value(j);
		index = data2.index(j);
		gd[x][1][z][index] = value;
      }
      y = l_y-4;
      PFdata& data3 = gd[x][y][z];
      N = data3.nonzero();
      for (j=0; j<N; j++) {
        value = data3.value(j);
		index = data3.index(j);
		gd[x][2][z][index] = value;
      }
      y = l_y-3;
      PFdata& data4 = gd[x][3][z];
      N = data4.nonzero();
      for (j=0; j<N; j++) {
		value = data4.value(j);
		index = data4.index(j);
		gd[x][y][z][index] = value;
      }
      y = l_y-2;
      PFdata& data5 = gd[x][4][z];
      N = data5.nonzero();
      for (j=0; j<N; j++) {
		value = data5.value(j);
		index = data5.index(j);
		gd[x][y][z][index] = value;
      }
      y = l_y-1;
      PFdata& data6 = gd[x][5][z];
      N = data6.nonzero();
      for (j=0; j<N; j++) {
		value = data6.value(j);
		index = data6.index(j);
		gd[x][y][z][index] = value;
      }
    }
  }


  // The following applies periodic boundary condition along z-axis
  for (x=3; x<l_x-3; x++) {
    for (y=0; y<l_y; y++) {
      z = l_z-1;
      PFdata& data1 = gd[x][y][5];
      N = data1.nonzero();
      for (j=0; j<N; j++) {
		value = data1.value(j);
		index = data1.index(j);
		gd[x][y][z][index] = value;
      }
      z = l_z-2;
      PFdata& data2 = gd[x][y][4];
      N = data2.nonzero();
      for (j=0; j<N; j++) {
		value = data2.value(j);
		index = data2.index(j);
		gd[x][y][z][index] = value;
      }
      z = l_z-3;
      PFdata& data3 = gd[x][y][3];
      N = data3.nonzero();
      for (j=0; j<N; j++) {
		value = data3.value(j);
		index = data3.index(j);
		gd[x][y][z][index] = value;
      }
      z = l_z-6;
      PFdata& data4 = gd[x][y][z];
      N = data4.nonzero();
      for (j=0; j<N; j++) {
		value = data4.value(j);
		index = data4.index(j);
		gd[x][y][0][index] = value;
      }
      z = l_z-5;
      PFdata& data5 = gd[x][y][z];
      N = data5.nonzero();
      for (j=0; j<N; j++) {
		value = data5.value(j);
		index = data5.index(j);
		gd[x][y][1][index] = value;
      }
        z = l_z-4;
      PFdata& data6 = gd[x][y][z];
      N = data6.nonzero();
      for (j=0; j<N; j++) {
		value = data6.value(j);
		index = data6.index(j);
		gd[x][y][2][index] = value;
      }
    }
  }
  

  /* prepare sending data*/
  count=0;
  node_count=0;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      count+=gd[3][y][z].nonzero();
      count+=gd[4][y][z].nonzero();
      count+=gd[5][y][z].nonzero();
	  node_count++;
	  node_count++;
      node_count++;
    }
  }
  
  senddata=(buff_data *) malloc(count*(sizeof(buff_data))); 

 // store the size of non-zero eta vector for each node
  sendcutcount=l_y*l_z*3;
  sendcut= (int*) malloc(sendcutcount*sizeof(int));	
  recvcut= (int*) malloc(sendcutcount*sizeof(int));
  

  node_count=0;
  count=0;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      PFdata& data1 = gd[3][y][z];
      sendcut[node_count] = data1.nonzero();
      for (j=0; j<sendcut[node_count]; j++) {
		senddata[count].value = data1.value(j);
		senddata[count].index = data1.index(j);
		count++;
      }
      node_count++;
    }
  }

  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      PFdata& data1 = gd[4][y][z];
      sendcut[node_count] = data1.nonzero();
      for (j=0; j<sendcut[node_count]; j++) {
		senddata[count].value = data1.value(j);
		senddata[count].index = data1.index(j);
		count++;
      }
      node_count++;
    }
  }
  
    for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      PFdata& data1 = gd[5][y][z];
      sendcut[node_count] = data1.nonzero();
      for (j=0; j<sendcut[node_count]; j++) {
		senddata[count].value = data1.value(j);
		senddata[count].index = data1.index(j);
		count++;
      }
      node_count++;
    }
  }
  
  /* end of send data preparation*/

  
  MPI_Sendrecv(sendcut,3*l_y*l_z,MPI_INT,front,13*(front+1),recvcut,3*l_y*l_z,MPI_INT,next,13*(mpirank+1),MPI_COMM_WORLD,&status);

  // count record the total number of nonzero etas in the buff-zone
  MPI_Sendrecv(&count,1,MPI_INT,front,3*(front+1),&recvcount, 1, MPI_INT,next, 3*(mpirank+1),MPI_COMM_WORLD,&status);
  MPI_Sendrecv(&node_count,1,MPI_INT,front,4*(front+1),&recv_node_count, 1, MPI_INT,next, 4*(mpirank+1),MPI_COMM_WORLD,&status);


  recvdata=(buff_data *) malloc(recvcount*(sizeof(buff_data))); 
  MPI_Sendrecv(senddata,count,vectortype,front,11*(front+1),recvdata,recvcount,vectortype,next,11*(mpirank+1),MPI_COMM_WORLD,&status); 
   
  count=0;
  node_count=0;
   
  x = l_x-3;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      N=recvcut[y*l_z+z];
      for (j=0; j<N; j++) {
		gd[x][y][z][recvdata[count].index]=recvdata[count].value;
		count++;
      }
	  node_count++;
    }
  }
  
  x = l_x-2;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      N=recvcut[(y+l_y)*l_z+z];
      for (j=0; j<N; j++) {
		gd[x][y][z][recvdata[count].index]=recvdata[count].value;
		count++;
      }
	  node_count++;
    }
  }

  x = l_x-1;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      N=recvcut[(y+l_y+l_y)*l_z+z];
      for (j=0; j<N; j++) {
		gd[x][y][z][recvdata[count].index]=recvdata[count].value;
		count++;
      }
	  node_count++;
    }
  } 
  

  free (senddata);
  free (recvdata);
  free (sendcut);
  free (recvcut);
   
  // another side of system
  
  count=0;
  node_count=0;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      count+=gd[l_x-6][y][z].nonzero();
      count+=gd[l_x-5][y][z].nonzero();
      count+=gd[l_x-4][y][z].nonzero();
      node_count++;
      node_count++;
    }
  }
  
  senddata=(buff_data *) malloc(count*(sizeof(buff_data)));
  sendcutcount=l_y*l_z*3;
  sendcut= (int*) malloc(sendcutcount*sizeof(int));
  recvcut= (int*) malloc(sendcutcount*sizeof(int));
  
  node_count=0;
  count=0;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      PFdata& data1 = gd[l_x-6][y][z];
      sendcut[node_count] = data1.nonzero();
      for (j=0; j<sendcut[node_count]; j++) {
		senddata[count].value = data1.value(j);
		senddata[count].index = data1.index(j);
		count++;
      }
      node_count++;
    }
  }

  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      PFdata& data1 = gd[l_x-5][y][z];
      sendcut[node_count] = data1.nonzero();
      for (j=0; j<sendcut[node_count]; j++) {
		senddata[count].value = data1.value(j);
		senddata[count].index = data1.index(j);
		count++;
      }
      node_count++;
    }
  }
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      PFdata& data1 = gd[l_x-4][y][z];
      sendcut[node_count] = data1.nonzero();
      for (j=0; j<sendcut[node_count]; j++) {
		senddata[count].value = data1.value(j);
		senddata[count].index = data1.index(j);
		count++;
      }
      node_count++;
    }
  }
  /* end of send data preparation*/

  MPI_Sendrecv(sendcut,3*l_y*l_z,MPI_INT,next,13*(next+1),recvcut,3*l_y*l_z,MPI_INT,front,13*(mpirank+1),MPI_COMM_WORLD,&status);

  MPI_Sendrecv(&count,1,MPI_INT,next,3*(next+1),&recvcount, 1, MPI_INT,front, 3*(mpirank+1),MPI_COMM_WORLD,&status);
  MPI_Sendrecv(&node_count,1,MPI_INT,front,4*(front+1),&recv_node_count, 1, MPI_INT,next, 4*(mpirank+1),MPI_COMM_WORLD,&status);

  recvdata=(buff_data *) malloc(recvcount*(sizeof(buff_data))); 
  MPI_Sendrecv(senddata,count,vectortype,next,11*(next+1),recvdata,recvcount,vectortype,front,11*(mpirank+1),MPI_COMM_WORLD,&status); 
  
  count=0;
  node_count=0;
  x = 0;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      N=recvcut[y*l_z+z];
      for ( j=0; j<N; j++) {
		gd[x][y][z][recvdata[count].index]=recvdata[count].value;
		count++;
      }
	  node_count++;
    }
  }
  
  x = 1;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      N=recvcut[(y+l_y)*l_z+z];
      for ( j=0; j<N; j++) {
		gd[x][y][z][recvdata[count].index]=recvdata[count].value;
		count++;
      }
	  node_count++;
    }
  }
    x = 2;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
      N=recvcut[(y+l_y+l_y)*l_z+z];
      for ( j=0; j<N; j++) {
		gd[x][y][z][recvdata[count].index]=recvdata[count].value;
		count++;
      }
	  node_count++;
    }
  }
  
  free (senddata);
  free (recvdata);
  free (sendcut);
  free (recvcut);
  
 
  
  if(flag==1){
      return;
  }

  // clamp the eta to [0,1]
  double sum_val;
 
  for(int pfx=0;pfx<l_x;pfx++){
	  for(int pfy=0;pfy<l_y;pfy++){
		  for(int pfz=0;pfz<l_z;pfz++){
			  PFdata& data = gd[pfx][pfy][pfz];
			  N = data.nonzero();
			  sum_val = 0.0;
			  for(j=0;j<N;j++){
				  value = data.value(j);
				  sum_val += value; 
			  }
			 // sum_val = sqrt(sum_val);
		
			 /*( if(sum_val>0.0){
				  for(j=0;j<N;j++){
					  data.value(j) /= sum_val;
				  }
			  } */
		  }
	  }
  } 
	
	//

	
	
	
	
	
}/*end boundary_conditions()*/



void mesh_extend(double *DefEng_buf, int *GID, int mpirank, int NumPE, int L_x, int L_y, int L_z)
{
 int x, y, z, N, index, j;
 int l_x,l_y,l_z;
 int pIDX_yd1,pIDX_yd2,pIDX_zd1,pIDX_zd2;
  int next,front;
 // double value;
  int count,scount;
 // int node_count;
  int sendcutcount;
  //int *sendcut;
  //int *recvcut;
  int *sendgid;
  int *recvgid;
  double *sendrho;
  double *recvrho;
  //int recvcount;
//  int recv_node_count;
  MPI_Status status;


  next=mpirank+1;if(next==NumPE)next=0; 
  front=mpirank-1;if(front==-1)front=NumPE-1;
  
 l_x=L_x+6;
 l_y=L_y+6;
 l_z=L_z+6;
    //
    	//cp
//if (flag==3){
    //periodic conditions along y axis
    for (x=3;x<l_x-3; x++) {
    for (z=3; z<l_z-3; z++) {
    //y=l_y-6;
    pIDX_yd1 = ((x)*l_y+(l_y-6))*l_z +z;
     pIDX_yd2 = ((x)*l_y+(0))*l_z +z;
    DefEng_buf[pIDX_yd2]= DefEng_buf[pIDX_yd1];
    GID[pIDX_yd2] = GID [pIDX_yd1];
       pIDX_yd1 = ((x)*l_y+(l_y-5))*l_z +z;
     pIDX_yd2 = ((x)*l_y+(1))*l_z +z;
    DefEng_buf[pIDX_yd2]= DefEng_buf[pIDX_yd1];
      GID[pIDX_yd2] = GID [pIDX_yd1];
       pIDX_yd1 = ((x)*l_y+(l_y-4))*l_z +z;
     pIDX_yd2 = ((x)*l_y+(2))*l_z +z;
    DefEng_buf[pIDX_yd2]= DefEng_buf[pIDX_yd1];
      GID[pIDX_yd2] = GID [pIDX_yd1];
       pIDX_yd1 = ((x)*l_y+(l_y-3))*l_z +z;
     pIDX_yd2 = ((x)*l_y+(3))*l_z +z;
    DefEng_buf[pIDX_yd1]= DefEng_buf[pIDX_yd2];
      GID[pIDX_yd1] = GID [pIDX_yd2];
       pIDX_yd1 = ((x)*l_y+(l_y-2))*l_z +z;
     pIDX_yd2 = ((x)*l_y+(4))*l_z +z;
    DefEng_buf[pIDX_yd1]= DefEng_buf[pIDX_yd2];
      GID[pIDX_yd1] = GID [pIDX_yd2];
       pIDX_yd1 = ((x)*l_y+(l_y-1))*l_z +z;
     pIDX_yd2 = ((x)*l_y+(5))*l_z +z;
    DefEng_buf[pIDX_yd1]= DefEng_buf[pIDX_yd2];
      GID[pIDX_yd1] = GID [pIDX_yd2];
            }
    }
    
    for (x=3; x<l_x-3; x++) {
    for (y=0; y<l_y; y++) {
        pIDX_zd1 = ((x)*l_y+y)*l_z+(l_z-6);
         pIDX_zd2 = ((x)*l_y+y)*l_z+(0);
         DefEng_buf[pIDX_zd2]= DefEng_buf[pIDX_zd1];
         GID[pIDX_zd2] = GID[pIDX_zd1];
         pIDX_zd1 = ((x)*l_y+y)*l_z+(l_z-5);
         pIDX_zd2 = ((x)*l_y+y)*l_z+(1);
         DefEng_buf[pIDX_zd2]= DefEng_buf[pIDX_zd1];
         GID[pIDX_zd2] = GID[pIDX_zd1];
         pIDX_zd1 = ((x)*l_y+y)*l_z+(l_z-4);
         pIDX_zd2 = ((x)*l_y+y)*l_z+(2);
         DefEng_buf[pIDX_zd2]= DefEng_buf[pIDX_zd1];
         GID[pIDX_zd2] = GID[pIDX_zd1];
         pIDX_zd1 = ((x)*l_y+y)*l_z+(l_z-3);
         pIDX_zd2 = ((x)*l_y+y)*l_z+(3);
         DefEng_buf[pIDX_zd1]= DefEng_buf[pIDX_zd2];
         GID[pIDX_zd1] = GID[pIDX_zd2];
         pIDX_zd1 = ((x)*l_y+y)*l_z+(l_z-2);
         pIDX_zd2 = ((x)*l_y+y)*l_z+(4);
         DefEng_buf[pIDX_zd1]= DefEng_buf[pIDX_zd2];
         GID[pIDX_zd1] = GID [pIDX_zd2];
         pIDX_zd1 = ((x)*l_y+y)*l_z+(l_z-1);
         pIDX_zd2 = ((x)*l_y+y)*l_z+(5);
         DefEng_buf[pIDX_zd1]= DefEng_buf[pIDX_zd2];
         GID[pIDX_zd1] = GID [pIDX_zd2];
    }}
    
    //send data along x
    scount=3*l_y*l_z;
     sendrho=(double *) malloc(scount*(sizeof(double)));
     sendgid= (int *) malloc(scount*(sizeof(int)));
      recvrho=(double *) malloc(scount*(sizeof(double)));
     recvgid= (int *) malloc(scount*(sizeof(int)));
     
     count=0;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
    sendrho[count]=DefEng_buf[((3)*l_y+y)*l_z+z];
     sendgid[count]=GID[((3)*l_y+y)*l_z+z];
        count++;
            }}
              for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
    sendrho[count]=DefEng_buf[((4)*l_y+y)*l_z+z];
      sendgid[count]=GID[((4)*l_y+y)*l_z+z];
        count++;
            }}
              for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
    sendrho[count]=DefEng_buf[((5)*l_y+y)*l_z+z];
      sendgid[count]=GID[((5)*l_y+y)*l_z+z];
        count++;
            }}
            
              
  MPI_Sendrecv(sendrho,3*l_y*l_z,MPI_DOUBLE,front,13*(front+1),recvrho,3*l_y*l_z,MPI_DOUBLE,next,13*(mpirank+1),MPI_COMM_WORLD,&status);
  MPI_Sendrecv(sendgid,3*l_y*l_z,MPI_INT,front,13*(front+1),recvgid,3*l_y*l_z,MPI_INT,next,13*(mpirank+1),MPI_COMM_WORLD,&status);
 
             // MPI_Sendrecv(&count,1,MPI_INT,front,3*(front+1),&recvcount, 1, MPI_INT,next, 3*(mpirank+1),MPI_COMM_WORLD,&status);
            //  recvdata=(double *) malloc(count*(sizeof(double)));
    // recvdata2= (int *) malloc(count*(sizeof(int)));
            //   MPI_Sendrecv(senddata,count,MPI_DOUBLE,front,11*(front+1),recvdata,recvcount,MPI_DOUBLE,next,11*(mpirank+1),MPI_COMM_WORLD,&status); 
             //  MPI_Sendrecv(senddata2,count,MPI_INT,front,11*(front+1),recvdata2,recvcount,MPI_INT,next,11*(mpirank+1),MPI_COMM_WORLD,&status); 
count = 0;
 for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
  DefEng_buf[(((l_x-3)*l_y+y)*l_z+z)]=recvrho[count];
  GID[(((l_x-3)*l_y+y)*l_z+z)]=recvgid[count];
   count++;     
    }}
      
     for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
  DefEng_buf[(((l_x-2)*l_y+y)*l_z+z)]=recvrho[count];
  GID[(((l_x-2)*l_y+y)*l_z+z)]=recvgid[count];
   count++;     
    }}
     for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
  DefEng_buf[(((l_x-1)*l_y+y)*l_z+z)]=recvrho[count];
  GID[(((l_x-1)*l_y+y)*l_z+z)]=recvgid[count];
   count++;     
    }}
    
  free (sendrho);
  free(sendgid);
  free (recvrho);
    free(recvgid);
    //anotherside
   
        scount=3*l_y*l_z;
     sendrho=(double *) malloc(scount*(sizeof(double)));
     sendgid= (int *) malloc(scount*(sizeof(int)));
          recvrho=(double *) malloc(scount*(sizeof(double)));
     recvgid= (int *) malloc(scount*(sizeof(int)));
     count=0;
  for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
    sendrho[count]=DefEng_buf[((l_x-6)*l_y+y)*l_z+z];
     sendgid[count]=GID[((l_x-6)*l_y+y)*l_z+z];
        count++;
            }}
              for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
    sendrho[count]=DefEng_buf[((l_x-5)*l_y+y)*l_z+z];
      sendgid[count]=GID[((l_x-5)*l_y+y)*l_z+z];
        count++;
            }}
              for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
    sendrho[count]=DefEng_buf[((l_x-4)*l_y+y)*l_z+z];
      sendgid[count]=GID[((l_x-4)*l_y+y)*l_z+z];
        count++;
            }}
            MPI_Sendrecv(sendrho,3*l_y*l_z,MPI_DOUBLE,next,13*(next+1),recvrho,3*l_y*l_z,MPI_DOUBLE,front,13*(mpirank+1),MPI_COMM_WORLD,&status);
              MPI_Sendrecv(sendgid,3*l_y*l_z,MPI_INT,next,13*(next+1),recvgid,3*l_y*l_z,MPI_INT,front,13*(mpirank+1),MPI_COMM_WORLD,&status);
               
//                MPI_Sendrecv(&count,1,MPI_INT,next,3*(next+1),&recvcount, 1, MPI_INT,front, 3*(mpirank+1),MPI_COMM_WORLD,&status);
//               recvdata=(double *) malloc(count*(sizeof(double)));
//      recvdata2= (int *) malloc(count*(sizeof(int)));
//                 MPI_Sendrecv(senddata,count,vectortype,next,11*(next+1),recvdata,recvcount,vectortype,front,11*(mpirank+1),MPI_COMM_WORLD,&status); 
//    MPI_Sendrecv(senddata2,count,vectortype,next,11*(next+1),recvdata,recvcount,vectortype,front,11*(mpirank+1),MPI_COMM_WORLD,&status); 
  
count = 0;
 for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
  DefEng_buf[(((0)*l_y+y)*l_z+z)]=recvrho[count];
  GID[(((0)*l_y+y)*l_z+z)]=recvgid[count];
   count++;     
    }}
     for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
  DefEng_buf[(((1)*l_y+y)*l_z+z)]=recvrho[count];
  GID[(((1)*l_y+y)*l_z+z)]=recvgid[count];
   count++;     
    }}
     for (y=0; y<l_y; y++) {
    for (z=0; z<l_z; z++) {
  DefEng_buf[(((2)*l_y+y)*l_z+z)]=recvrho[count];
  GID[(((2)*l_y+y)*l_z+z)]=recvgid[count];
   count++;     
    }}
  free(sendrho);
  free(sendgid);
  free(recvrho);
    free(recvgid);
  
//}
}

void Update_grid(PFgrid3D& gd_rex, int *gID, int L_x, int L_y, int L_z, int l_x, int l_y, int l_z,
		int mpirank, int NumPE,
		double M_coeff, double kappa, double m, double *DefEng, double E_disl, int iter_step, double *rho_avg, double *defEng_buf, int *GID)
{
  static PFgrid3D temp_rex(l_x,l_y,l_z);
 // static PFgrid3D temp_def(l_x,l_y,l_z);
  static std::vector<int> nonzero_rex;
//  static std::vector<int> nonzero_def;
  std::vector<int>::iterator it, it1, it2, end;
//  std::vector<int>::iterator it2, end2;
  double f_def;
  int idx, k,g,pidx_b;
  int NN_rex,count,N_rex_nonlocal;
  double sum_rex,sum_clamp_rex, sum1, h_phi_def, f_coeff;
double val, max_val_rex;
int id, p_idx, index;
double *Avg=NULL;
int *count_dd=NULL;



  const int xlim = l_x-3;
  const int ylim = l_y-3;
  const int zlim = l_z-3;

  // Iterate through each grid node...
  for (int x=3; x<xlim; ++x) {
    for (int y=3; y<ylim; ++y) {
      for ( int z=3; z< zlim; ++z){

		  idx = ((x-3)*L_y+(y-3))*L_z+(z-3);
		  pidx_b = ((x)*(L_y+6)+(y))*(L_z+6)+(z);
		//  f_def = DefEng[idx];


	// Limits of this node's neigborhood.
	int xmin = x-1; int xmax = x+1;
	int ymin = y-1; int ymax = y+1;
	int zmin = z-1; int zmax = z+1;

	// Make a list of nonzero order parameters in this neighborhood.
	// the list could be different for eta_rex and eta_def
	for (int i=xmin; i<=xmax; ++i){
	  for (int j=ymin; j<=ymax; ++j){
	    for (int m=zmin; m<=zmax; ++m) {

	      PFdata& data_rex = gd_rex[i][j][m];
	    //  PFdata& data_def = gd_def[i][j][m];
	   //   PFdata& data_dd = gd_dd[i][j][m];
	      int N_rex = data_rex.nonzero();
	    //  int N_def = data_def.nonzero();
       //  int pIDX_d = ((px-2)*L_y+py-2)*L_z+pz-2;
		  // construct nonzero_rex
	      for (int k=0; k<N_rex; ++k){
			if (data_rex.value(k)>EPSILON) {
				int index = data_rex.index(k);	// index is one of the grain IDs in consideration

				// iterate over the current nonzero_rex (list of indices (grain IDs) of nonzero etas)
				end = nonzero_rex.end();
				for (it=nonzero_rex.begin(); it!=end; ++it)
				if ((*it)==index) break;
				if (it==end) nonzero_rex.push_back(index);
			}
		  }

		  // construct nonzero_def
	 /*     for (int k=0; k<N_def; ++k){
			if (data_def.value(k)>EPSILON) {
				int index = data_def.index(k);

				end = nonzero_def.end();
				for (it=nonzero_def.begin(); it!=end; ++it)
				if ((*it)==index) break;
				if (it==end) nonzero_def.push_back(index);
			}
		  } */
		}
	  }
	}


	// Use local variables to avoid excessive subscript operator calls.
	PFdata& update_rex = temp_rex[x][y][z];
//	PFdata& update_def = temp_def[x][y][z];

	/***********************************
	  neighboring grids of _rex
	  **********************************/
	PFdata& right_rex = gd_rex[xmax][y][z];
	PFdata& left_rex = gd_rex[xmin][y][z];
	PFdata& front_rex = gd_rex[x][ymax][z];
	PFdata& back_rex = gd_rex[x][ymin][z];
	PFdata& up_rex = gd_rex[x][y][zmax];
	PFdata& down_rex = gd_rex[x][y][zmin];
	PFdata& self_rex = gd_rex[x][y][z];


	PFdata& left_up_rex = gd_rex[xmin][y][zmax];
	PFdata& left_down_rex = gd_rex[xmin][y][zmin];
	PFdata& left_front_rex = gd_rex[xmin][ymin][z];
	PFdata& left_back_rex = gd_rex[xmin][ymax][z];
	PFdata& center_up_front_rex = gd_rex[x][ymin][zmax];
	PFdata& center_up_back_rex = gd_rex[x][ymax][zmax];
	
	PFdata& center_down_front_rex = gd_rex[x][ymin][zmin];
	PFdata& center_down_back_rex = gd_rex[x][ymax][zmin];
	PFdata& right_up_rex = gd_rex[xmax][y][zmax];
	PFdata& right_down_rex = gd_rex[xmax][y][zmin];
	PFdata& right_front_rex = gd_rex[xmax][ymin][z];
	PFdata& right_back_rex = gd_rex[xmax][ymax][z];
	//PFdata& dd_ave = gd_dd[x][y][z];

	/***********************************
	  neighboring grids of _def
	  **********************************/
/*	PFdata& right_def = gd_def[xmax][y][z];
	PFdata& left_def = gd_def[xmin][y][z];
	PFdata& front_def = gd_def[x][ymax][z];
	PFdata& back_def = gd_def[x][ymin][z];
	PFdata& up_def = gd_def[x][y][zmax];
	PFdata& down_def = gd_def[x][y][zmin];
	PFdata& self_def = gd_def[x][y][z];


	PFdata& left_up_def = gd_def[xmin][y][zmax];
	PFdata& left_down_def = gd_def[xmin][y][zmin];
	PFdata& left_front_def = gd_def[xmin][ymin][z];
	PFdata& left_back_def = gd_def[xmin][ymax][z];
	PFdata& center_up_front_def = gd_def[x][ymin][zmax];
	PFdata& center_up_back_def = gd_def[x][ymax][zmax];
	
	PFdata& center_down_front_def = gd_def[x][ymin][zmin];
	PFdata& center_down_back_def = gd_def[x][ymax][zmin];
	PFdata& right_up_def = gd_def[xmax][y][zmax];
	PFdata& right_down_def = gd_def[xmax][y][zmin];
	PFdata& right_front_def = gd_def[xmax][ymin][z];
	PFdata& right_back_def = gd_def[xmax][ymax][z];     */

	// calculate the sum terms in the deformation energy term
	NN_rex = self_rex.nonzero();
		 Avg = (double *)malloc(NN_rex * sizeof(double));
	 count_dd = (int *)malloc(NN_rex*sizeof(int));
	
f_coeff = 0.0;
        count = 0;
       // sum_clamp_rex = 0.0;
	

// finding dd_avg for all nonzero grains in a cube box
// Limits of this node's neigborhood.
	//int xmin = x-1; int xmax = x+1;
	//int ymin = y-1; int ymax = y+1;
	//int zmin = z-1; int zmax = z+1;

for(k=0;k<NN_rex;k++){
Avg[k] = 0.0;
count_dd[k] = 0;
}
if(NN_rex > 1) {
        for (int i=x-3; i<=x+3; ++i){
	  for (int j=y-3; j<=y+3; ++j){
	    for (int m=z-3; m<=z+3; ++m) {
              //PFdata& data_rex_nonlocal = gd_rex[i][j][m];
             // N_rex_nonlocal = data_rex_nonlocal.nonzero();
              p_idx = ((i)*(L_y+6)+(j))*(L_z+6)+(m);
                //id = 0;
                //max_val_rex = 0.0;
                 //for(k=0;k<N_rex_nonlocal;k++){
                  //val = data_rex_nonlocal.value(k);
                //    if(max_val_rex < val) {
                  //     id = data_rex_nonlocal.index(k);
                //        }
            //    }
                 
                 for(k=0;k<NN_rex;k++) {
                 if(GID[p_idx] == self_rex.index(k)) {
                     Avg[k] += defEng_buf[p_idx]; 
                 count_dd[k] += 1;
                   }

            }
   
        }
          }
        }
}
    
for(k=0;k<NN_rex;k++) {
    index = self_rex.index(k);
    rho_avg[index-1] = 0.0;
   // printf("index = %d\n",index);
    if(count_dd[k] > 0) {
    rho_avg[index-1] = Avg[k]/count_dd[k];
}
}
//printf("rho_avg(1) = %le and rho_avg(45) = %le\n",rho_avg[0],rho_avg[44]);
sum_rex = 0.0;
        sum1 = 0.0;
        h_phi_def = 0.0;
//for(k=0;k<NN_rex;k++){
		
if(NN_rex > 1) {
for(k=0;k<NN_rex;k++){
                if(count_dd[k] > 0 && Avg[k] > 0.0) {
                    sum_rex += (self_rex.value(k))*(self_rex.value(k));
                sum1 += self_rex.value(k)*self_rex.value(k)*self_rex.value(k)*(10.0 - 15.0*self_rex.value(k) + 6.0*self_rex.value(k)*self_rex.value(k));
                h_phi_def += (self_rex.value(k)*self_rex.value(k)*self_rex.value(k)*(10.0 - 15.0*self_rex.value(k) + 6.0*self_rex.value(k)*self_rex.value(k)))*(Avg[k]/count_dd[k]);
	}
}
}
if(NN_rex == 1) {
    for(k=0;k<NN_rex;k++){
    sum_rex += (self_rex.value(k))*(self_rex.value(k));
    sum1 += self_rex.value(k)*self_rex.value(k)*self_rex.value(k)*(10.0 - 15.0*self_rex.value(k) + 6.0*self_rex.value(k)*self_rex.value(k));
     rho_avg[index-1]= defEng_buf[pidx_b];
    }
}
	
//	NN_def = self_def.nonzero();
//	sum_def = 0.0;
//	for(k=0;k<NN_def;k++){
//		sum_def += (self_def.value(k))*(self_def.value(k));
//	}
//	sum_tot = sum_rex + sum_def;
//	sum_tot = sum_tot*sum_tot;

	// Calculate update values for each non-zero order parameter.
//	int t =0;
	/*end = nonzero_rex.end();
	for (it=nonzero_rex.begin(); it!=end; ++it) {
	  int index = (*it);
	  {

	  //cp
      //    if (mpirank==0){
    //    printf("dd_average_rex=%lf %d \n", dd_average[index], index);
      //    }	
      //exit(0);
              sum1 = 0.0;
              for (it1=nonzero_rex.begin(); it1!=end; ++it1) {
	  int index1 = (*it1);
          for (it2=nonzero_rex.begin(); it2!=end; ++it2) {
	  int index2 = (*it2);
          if(index1 != index && index2 != index && index1<index2) {
              sum1 += self_rex[index1]*self_rex[index2];
          }
          }
              }
              
              
      
              double center_rex = self_rex[index];

			  double lap_rex =laplacian(right_rex[index], left_rex[index], front_rex[index], back_rex[index],
	                          up_rex[index], down_rex[index], left_up_rex[index], left_down_rex[index], 
                                  left_front_rex[index], left_back_rex[index],center_up_front_rex[index], 
                                  center_up_back_rex[index], center_down_front_rex[index], center_down_back_rex[index],
                                  right_up_rex[index], right_down_rex[index],
	                          right_front_rex[index], right_back_rex[index], center_rex);
   // if(fabs(lap_rex) > 0.0) {
             // double value_rex = center_rex - 
				 // M_coeff*(0.242*(center_rex*center_rex*center_rex - center_rex + 3.0*center_rex*(sum_rex - center_rex*center_rex))
						//  -kappa*lap_rex - 6.0*center_rex*(1.0-center_rex)*(f_def - E_disl*dd_ave[index]));
						 // -2.0*center_rex*sum_def*(f_def-E_disl*dd_ave[index])/(sum_tot));
              
         sum_clamp_rex +=    - M_coeff*(0.1875E06*(center_rex*center_rex*center_rex - center_rex + 3.0*center_rex*(sum_rex - center_rex*center_rex))
						  -kappa*lap_rex - (6.0*center_rex*(1.0-center_rex) + 2.0*sum1)*(f_def - E_disl*dd_ave[index]));                     
                                  
      //   count++;                         
                                  
   //}
   //  if (value_rex>1E-04) update_rex[index] = value_rex;
   // }
/* if(iter_step%20 == 0.0) {
           	 if(mpirank==0 && index == 25){
                
		printf("%lf\n",value_rex);
		}

}



 
	 //   t++;
          }
	} */


end = nonzero_rex.end();
	for (it=nonzero_rex.begin(); it!=end; ++it) {
	  int index = (*it);
          
              if(index > init_gid) {
                  count ++;
              }
          }
          


		end = nonzero_rex.end();
	for (it=nonzero_rex.begin(); it!=end; ++it) {
	  int index = (*it);
          {
              

//exit(0); 
          
              double h_phi = self_rex[index]*self_rex[index]*self_rex[index]*(10.0 - 15.0*self_rex[index] + 6.0*self_rex[index]*self_rex[index]);
              double dh_phi = self_rex[index]*self_rex[index]*(30.0 - 60.0*self_rex[index] + 30.0*self_rex[index]*self_rex[index]);
               double center_rex = self_rex[index];

			  double lap_rex =laplacian(right_rex[index], left_rex[index], front_rex[index], back_rex[index],
	                          up_rex[index], down_rex[index], left_up_rex[index], left_down_rex[index], 
                                  left_front_rex[index], left_back_rex[index],center_up_front_rex[index], 
                                  center_up_back_rex[index], center_down_front_rex[index], center_down_back_rex[index],
                                  right_up_rex[index], right_down_rex[index],
	                          right_front_rex[index], right_back_rex[index], center_rex);
                if(count == 0 && rho_avg[index-1] > 0.0) {
              double value_rex = center_rex - 
				 M_coeff*(m*(center_rex*center_rex*center_rex - center_rex + 3.0*center_rex*(sum_rex - center_rex*center_rex))
						  -kappa*lap_rex + E_disl*dh_phi*((sum1-h_phi)*rho_avg[index-1] - (h_phi_def - h_phi*rho_avg[index-1]))/(sum1*sum1));   //- (sum_clamp_rex/NN_rex);
   
if(value_rex > 1E-04){
 temp_rex[x][y][z][index] = value_rex;
}
}
                 if(count > 0 && rho_avg[index-1] > 0.0) {
              double value_rex = center_rex - 
				  M_coeff*(m*(center_rex*center_rex*center_rex - center_rex + 3.0*center_rex*(sum_rex - center_rex*center_rex))
						  -kappa*lap_rex + E_disl*dh_phi*((sum1-h_phi)*rho_avg[index-1] -  (h_phi_def - h_phi*rho_avg[index-1]))/(sum1*sum1));   //- (sum_clamp_rex/NN_rex);
             

if(index > init_gid && center_rex < 0.5 && value_rex > 0.5) {
          defEng_buf[pidx_b] = 48.0;
}

if(value_rex > 1E-04) {
temp_rex[x][y][z][index] = value_rex;
}
}
 
     
         
              
              

 //  }
          }
	}

/*if(NN_rex == 1) {
   end = nonzero_rex.end();
	for (it=nonzero_rex.begin(); it!=end; ++it) {
	  int index = (*it);
    temp_rex[x][y][z][index] = self_rex[index];
        }
} */
	

	
	
	
	
free(Avg);
free(count_dd);	
	
nonzero_rex.clear();

	// Calculate update values for each non-zero order parameter.
/*	end = nonzero_def.end();
	for (it=nonzero_def.begin(); it!=end; ++it) {
	  int index = (*it);
          {     
           if (mpirank==0){
   printf("dd_average_def=%lf %d \n", dd_ave[index], index);
       }	
      // exit(0);
              double center_def = self_def[index];
              
			  double lap_def =laplacian(right_def[index], left_def[index], front_def[index], back_def[index],
	                          up_def[index], down_def[index], left_up_def[index], left_down_def[index], 
                                  left_front_def[index], left_back_def[index],center_up_front_def[index], 
                                  center_up_back_def[index], center_down_front_def[index], center_down_back_def[index],
                                  right_up_def[index], right_down_def[index],
	                          right_front_def[index], right_back_def[index], center_def);

              double value_def = center_def - 
				  M_coeff*(0.625*(center_def*center_def*center_def - center_def + 3.0*center_def*(sum_rex + sum_def - center_def*center_def))
						  -kappa*lap_def -6.0*center_def*(1.0-center_def)*(f_def - E_disl*gd_dd[x][y][z][grain_ID]);
						 // +((2.0*center_def/(sqrt(sum_tot))-2.0*center_def*sum_def/(sum_tot)))*(f_def));
              if (value_def>1E-05) update_def[index] = value_def;
     //         	if(mpirank==0){
	//	printf("INDEX = %d , value_def = %lf\n",
	//				index,value_def);
	//	} 
          }
	}
	nonzero_def.clear();  */
	}
   }
  }

//exit(0);
  /* clamping eta */
  //clamping_eta(temp_rex,xlim,ylim,zlim);
  boundary_conditions(temp_rex,mpirank,NumPE,l_x,l_y,l_z,0);
  //boundary_conditions(temp_def,mpirank,NumPE,l_x,l_y,l_z,0);
 //Update_gID(gd_rex, gID, L_x, L_y, L_z, xlim,ylim,zlim);


  // Swap the contents of the original grid with the temporary grid.
  gd_rex.swap(temp_rex);
 // gd_def.swap(temp_def);
  
  // Clear out the contents of the temporary grid.
  for (int x=0; x<l_x; x++){
    for (int y=0; y<l_y; y++){
      for (int z=0; z<l_z; z++){
		temp_rex[x][y][z].clear();
	//	temp_def[x][y][z].clear();
	  }
	}
  }


	return;
}/*end Update_grid()*/

void GetStatistics(PFgrid3D& gd_rex, int lnx, int lny, int lnz, int Nxyz,int iter_step,int mpirank,int lsize_d,int istep, int count_1)
{
	int i,j, N, count,pIDX_d,size;
	double tmp_rex, def_frac_local, def_frac,sum;
	
	double l[lsize_d];

	def_frac_local = 0.0;
	count=0;
	size=lsize_d;
	
/*	for(int px=2;px<lnx;px++){
		for(int py=2;py<lny;py++){
			for(int pz=2;pz<lnz;pz++){
			   int pIDX_d = (((px-2)*(lny+2)+(py-2))*(lnz+2)+pz-2);
			   l[pIDX_d]=0.0;
			}
		}
	}*/
	
	
	//if(iter_step%200 == 0.0){
	//MPI_File fp;
	//MPI_Status status;
//	char fname[80];
	 //char fname[100]={0};
      //sprintf(fname,"Inner_PF_STEP_%06d",iter_step);
	
//	sprintf(fname, "%s.iout", fn);
	//MPI_File_open(MPI_COMM_WORLD, fname,
	//		MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	//MPI_File_set_view(fp, mpirank*lsize_d*sizeof(int), MPI_INT,
	//		MPI_INT, "native", MPI_INFO_NULL);
//}
	for(int px=3;px<lnx;px++){
		for(int py=3;py<lny;py++){
			for(int pz=3;pz<lnz;pz++){
				tmp_rex = 0.0;
				pIDX_d = (((px-3)*(lny-3)+(py-3))*(lnz-3)+pz-3);
				
				N = gd_rex[px][py][pz].nonzero();
			/*	for(i=0;i<N;i++){
					tmp_rex += gd_rex[px][py][pz].value(i)*gd_rex[px][py][pz].value(i);
				
				} */
		
			/*	tmp_def = 0.0;
				N = gd_def[px][py][pz].nonzero();
				for(i=0;i<N;i++){
					tmp_def += gd_def[px][py][pz].value(i)*gd_def[px][py][pz].value(i);
						
				} */
				l[pIDX_d] = 0.0;
                               sum = 0.0;
                               int check =0;
			for(i=0;i<N;i++) {
                           sum += gd_rex[px][py][pz].value(i);
                           if (gd_rex[px][py][pz].index(i) ==1){
                                check ==1;
                           }
                            //if(gd_rex[px][py][pz].index(i) == 9) {
				 for(j=0;j<N;j++){
				    if(j>i) {
				         l[pIDX_d] += gd_rex[px][py][pz].value(i)*gd_rex[px][py][pz].value(j);
				         
				         
				    } 
				 }
  
}
if(sum <0.5) {
printf("sum is %le, %ld\n",sum,pIDX_d);
assert(sum>0.5);
}
				
                        

				//}
				//if(mpirank == 0 && iter_step%100 == 0.0) {
					  //  printf("%d %d %d %le\n",px-2,py-2,pz-2,l[pIDX_d]);
					   //	MPI_File_write(fp,l, lsize_d, MPI_INT, &status);
					//}
				  
					//else if(mpirank == 0 && iter_step%200 == 0.0)  {
					  //  MPI_File_write(fp,0, lsize_d, MPI_INT, &status);
					//}
				
		
			/*	double tmp_f = tmp_def/(tmp_rex+tmp_def);
				def_frac_local += tmp_f;
				if(tmp_f<0.5){
					grex[count]=1;
				}
				else{
					grex[count]=0;
				} */
				count++;
			}
		}
	}
	//if(iter_step%200 == 0.0){
	  //	MPI_File_close(&fp);
	//}
//	MPI_Allreduce(&def_frac_local, &def_frac, 1,
		//	MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				if(count_1%400 == 0 && iter_step%17 == 0.0){
	    printf("file_write_start");
	    				 //    if(l> 0.0 && iter_step%200 == 0.0) {
				         	//FILE *fp;
	      MPI_File fp;
	        MPI_Status status;
//	char fname[80];
	        char fname[100]={0};
            sprintf(fname,"PF_STEP_%06d_%d.iout",istep,iter_step);
//fp = fopen(fname,"w");
                // if(mpirank == 0){
//for(int px=2;px<lnx;px++){
		//for(int py=2;py<lny;py++){
			//for(int pz=2;pz<lnz;pz++){
				//tmp_rex = 0.0;
				//pIDX_d = (((px-2)*(lny-2)+(py-2))*(lnz-2)+pz-2);
				
//fprintf(fp,"%le\n",l[pIDX_d]);
//}////
//}
///}
//}
//fclose(fp);	
//	sprintf(fname, "%s.iout", fn);
	        MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	        MPI_File_set_view(fp, mpirank*size*sizeof(double), MPI_DOUBLE,
			MPI_DOUBLE, "native", MPI_INFO_NULL);
		//	for(int px=2;px<lnx;px++){
		//for(int py=2;py<lny;py++){
		//	for(int pz=2;pz<lnz;pz++){
		//	    	int pIDX_d = (((px-2)*lny+py-2)*lnz+pz-2);
					   // printf("%d %d %d %le\n",px,py,pz,l);
		//			   if (l[pIDX_d]>0.0){
					   	MPI_File_write(fp,l, size, MPI_DOUBLE, &status);
		//			}
		//			else {
		//			    MPI_File_write(fp,0, lsize_d, MPI_INT, &status);
		//			}
				
	  
	//}
	//	}
	//		}
				MPI_File_close(&fp);
	    				   //  }
	} 
//	return (def_frac /= Nxyz);

} /*end GetStatistics()*/

void PF_XRX(int L_x, int L_y, int L_z, int NumPE, int mpirank,
		double E_gb, double M_gb, double mu, double bb, double alpha,
		int* gID, double* rho, int Steps, int *gID_new, int count, int Flagwrite,double* rho_avg)
{
	int iter_step;
	int lsize_d, size_d;	// division of original data
	int l_x, l_y, l_z;	// size of local array, including extra 4 buffer layers
	int lnx, lny, lnz;	// position of local data index
	int DX, Nxyz;
	double L_gb;	// GB width
	double dx;		// voxel length
	double dt;		// physical time of one PF step
	double M_coeff;	// combined coefficient for dynamic equation
	double kappa;	// combined coefficient for gradient term
	double E_disl;	// unit stored energy corresponding to a dislocation density of 1E12/m^2
	double def_frac;	// the average fraction of deformed
	double m;         //potenital_well_height
	
	/* Material properties */
	/* Moelans, et. al., PRB, 88 (2013), 054103 */
	/* set 6*E_gb = L_gb, L_gb = sqrt(9.6)*dx */
	 
	dx = pf_length_scale; //L_gb/sqrt(9.0);  // [mm], same as L0 in FFT grid 
	
        L_gb =3.0*dx; //28.8*E_gb;
        dt = pf_time_step;  //0.0375*L_gb/M_gb; // [s]
	M_coeff = dt*4.0/3.0*M_gb/L_gb; 
	kappa = 3.0/4.0*E_gb*L_gb/dx/dx;
	m=6*E_gb/L_gb;
	E_disl = alpha*mu*bb*bb;	// [KJ/m^3] by multiplying a density of 1E12/m^2 and in shear modulus in GPa
	
//printf("init_gid dt dx m=%d %lf %lf %lf\n",init_gid,dt,dx,m);
//exit(0);
	if(mpirank==0){
		printf("Voxel length = %lf \n",dx);
		printf("Time increment = %lf \n",dt);
		printf("Grain boundary length (reduced): L_gb/dx = %d\n",(int)(L_gb/dx));
	}
//	assert(L_x%NumPE==0);	// only allow an equal distribution among PEs
	// slab-decomposition along x-axis
	DX = L_x/NumPE;
	Nxyz = L_x*L_y*L_z;
	// Add extra two layers as buffers to store the neighboring data
//	l_x = DX+4; l_y = L_y+4; l_z = L_z+4;
    l_x = DX+6; l_y = L_y+6; l_z = L_z+6;    //extra 3 layers
	lsize_d = DX*L_y*L_z;
     size_d=(l_x)*(l_y)*(l_z);
//	lnx = l_x-2; lny =l_y-2; lnz =l_z-2;	// position of local data index
    	lnx = l_x-3; lny =l_y-3; lnz =l_z-3;	// position of local data index
	int GID[(l_x)*(l_y)*(l_z)];
	/*****************************************
	  Allocate PF grids and stored energy
	  ****************************************/
	PFgrid3D grid_rex(l_x, l_y, l_z);
//	PFgrid3D grid_def(l_x, l_y, l_z);
	//PFgrid3D grid_dd(l_x, l_y, l_z);
//	double *DefEng = (double*)malloc(lsize_d*sizeof(double));
    	double *DefEng = (double*)malloc(lsize_d*sizeof(double));
    double *DefEng_buf = (double*)malloc(size_d*sizeof(double));

	/*****************************************
	  initialized PF grids and stored energy
	  using the input gID, gID_rex, and rho
	  ****************************************/
	for(int px=3;px<lnx;px++){
		for(int py=3;py<lny;py++){
			for(int pz=3;pz<lnz;pz++){
				int pIDX_d = ((px-3)*L_y+py-3)*L_z+pz-3;	// index in original data, excluding extra 4 layers
				int pIDX_db = (px*l_y+py)*l_z+pz;
		        int grainID=gID[pIDX_d];
		      //  int rexID=gID_rex[pIDX_d];
		        // need to write grain ID based phase field order parameter
			//	if(rexID==1){
					grid_rex[px][py][pz][grainID] = 1.0;
			//		grid_def[px][py][pz][grainID] = 0.0;
			//	}else{
			//		grid_rex[px][py][pz][grainID] = 0.0;
			//		grid_def[px][py][pz][grainID] = 1.0;
			//	} */
				//if(first_pf[pIDX_d]==1){
                                    // rho_avg[grainID -1] *= 0.25;
				   //grid_dd[px][py][pz][grainID] = 1.2E13;
                                    // rho[pIDX_d] *= 0.2 ;
			//}else{
				   // grid_dd[px][py][pz][grainID] = dd_average[pIDX_d];
			//}
		
				DefEng[pIDX_d] = rho[pIDX_d];
//                  if(grainID==1){
//                       DefEng_buf[pIDX_db]=500;
//                       DefEng[pIDX_d]=500;
//                }
                DefEng_buf[pIDX_db]= rho[pIDX_d];
     //           DefEng[pIDX_d]=1.8;
    //            exit(0);
        //    local_loop{			
				//if (mpirank==0){
           // printf("dd_average= %le %d\n", dd_average[grainID-1], grainID);
          // }
       //     }

          
				// copy gID to gID_new
				gID_new[pIDX_d] = grainID;
                GID[pIDX_db] = grainID;
			}
			
		}
		//exit(0);
	}
 //
// exit(0);
	/******************************************
	  Apply the boundary condition to the
	  initial data
	  *****************************************/
	boundary_conditions(grid_rex,mpirank,NumPE,l_x,l_y,l_z,0);
   //  boundary_conditions(DefEng_buf,mpirank,NumPE,l_x,l_y,l_z,3);
   // boundary_conditions(GID,mpirank,NumPE,l_x,l_y,l_z,3);
    mesh_extend(DefEng_buf,GID,mpirank, NumPE, DX,L_y,L_z);
//	boundary_conditions(grid_def,mpirank,NumPE,l_x,l_y,l_z,0);
	//boundary_conditions(grid_dd,mpirank,NumPE,l_x,l_y,l_z,1);
	// Update_gID(grid_rex, gID_new,L_x, L_y, L_z, lnx,lny,lnz);
       
      //  exit(0);

	/*****************************************
	  microstructure evolution
	  ****************************************/
//	def_frac = GetStatistics(grid_rex,grid_def,lnx, lny, lnz, Nxyz, grex_new);
//	if(mpirank==0){
//		printf("STEP = %d (TIME = %lf): Frac_of_Def = %lf\n",
//				0,0.0,def_frac);
//	}
//	for(iter_step=0;iter_step<Steps;iter_step++){
//		Update_grid(grid_rex, grid_def, gID_new, L_x, L_y, L_z, l_x, l_y, l_z, 
//				mpirank, NumPE,
//				M_coeff, kappa, DefEng);
//		def_frac = GetStatistics(grid_rex,grid_def,lnx, lny, lnz, Nxyz, grex_new);
//		if(mpirank==0){
//			printf("STEP = %d (TIME = %lf): Frac_of_Def = %lf\n",
//					iter_step,iter_step*dt,def_frac);
//		}
//	}

	iter_step=0;
	while(iter_step<Steps){
		
	Update_grid(grid_rex, gID_new, L_x, L_y, L_z, l_x, l_y, l_z, 
				mpirank, NumPE,
				M_coeff, kappa,m, DefEng,E_disl,iter_step,rho_avg,DefEng_buf, GID);
			//	printf("%le\n",rho_avg[24]);
			//	if(iter_step==200){
			//	exit(0);
//	}
Update_gID(grid_rex, gID_new,L_x, L_y, L_z, lnx,lny,lnz);
mesh_extend(DefEng_buf,GID,mpirank, NumPE, DX,L_y,L_z);
              //   if(count%400 == 0) {
		GetStatistics(grid_rex,lnx, lny, lnz, Nxyz,iter_step,mpirank, lsize_d,Flagwrite,count);
             // }
	iter_step++;
		 //   if(iter_step%10 == 0.0) {
			//		    printf("%d %d %d %le\n",px,py,pz,grid_rex[px][py][pz][grain_ID]);
				//	}
		//	printf("STEP = %d (TIME = %lf): Frac_of_Def = %lf\n",
				//	iter_step,iter_step*dt,def_frac);
		 
	//	if(def_frac<0.1){
		    
	//		break;
	//	}
		//cp
	//exit(0);

		//if(FlagWrite>0){
     // char fname[100]={0};
     // sprintf(fname,"Inner_PF_STEP%d_gID_%06d",FlagWrite,iter_step);
    //  WriteIDMPI(gID_new,mpirank,lsize_d,fname);
		//}
	}
	
	
	//	Update_gID_final(grid_rex, gID_new, L_x, L_y, L_z, lnx,lny,lnz,mpirank);
                // exit(0);
		
	/*	if(def_frac <0.1) {
		    for(int px=2;px<lnx;px++){
		for(int py=2;py<lny;py++){
			for(int pz=2;pz<lnz;pz++){
				int pIDX_d = ((px-2)*L_y+py-2)*L_z+pz-2;
				gID_rex[pIDX_d] = 0;
		}
		}
		    }
		} */

/*	for(int px=2;px<lnx;px++){
		for(int py=2;py<lny;py++){
			for(int pz=2;pz<lnz;pz++){
				int pIDX_d = ((px-2)*L_y+py-2)*L_z+pz-2;	// index in original data, excluding extra 4 layers
		        int grainID=gID[pIDX_d];
				if(grid_rex[px][py][pz][grainID]>0.0){
					gID_rex[pIDX_d] = 1.0;
				}else{
					gID_rex[pIDX_d] = 0.0;
				}
			}
		}
	}*/


      
	/*****************************************
	  clear the contents of the PF grids and
	  free the stored energy array
	  ****************************************/
	for(int x=0; x<l_x; x++){
		for(int y=0; y<l_y; y++){
			for(int z=0; z<l_z; z++){
				grid_rex[x][y][z].clear();
			//	grid_def[x][y][z].clear();
				//grid_dd[x][y][z].clear();
			}
		}
	}
	free(DefEng);
	free(DefEng_buf);

	return;
}/*end PF_XRX()*/

void WriteIDMPI(int *ID,int mpirank,int lsize,const char* fn)
{
	MPI_File fp;
	MPI_Status status;
	char fname[80];

	sprintf(fname, "%s.iout", fn);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize*sizeof(int), MPI_INT,
			MPI_INT, "native", MPI_INFO_NULL);
	MPI_File_write(fp, ID, lsize, MPI_INT, &status);
	MPI_File_close(&fp);

	return;
}/*end WriteIDMPI()*/

/*void WriterexMPI(int *ID,int mpirank,int lsize,const char* fn)
{
	MPI_File fp;
	MPI_Status status;
	char fname[80];

	sprintf(fname, "%s.iout", fn);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize*sizeof(int), MPI_INT,
			MPI_INT, "native", MPI_INFO_NULL);
	MPI_File_write(fp, ID, lsize, MPI_INT, &status);
	MPI_File_close(&fp);

	return;
}*/
