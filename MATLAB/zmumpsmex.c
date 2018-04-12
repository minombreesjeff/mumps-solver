#include "mex.h"
#include "zmumps_c.h"

#define dmumps_c       zmumps_c				
#define dmumps_par		 zmumps_par			
#define DMUMPS_STRUC_C ZMUMPS_STRUC_C	
#define DMUMPS_alloc   ZMUMPS_alloc     
#define DMUMPS_free    ZMUMPS_free
#define double2        mumps_double_complex
#define mxREAL2        mxCOMPLEX

#define SYM        (prhs[0])
#define JOB        (prhs[1])
#define ICNTL      (prhs[2])
#define CNTL       (prhs[3])
#define PERM_IN    (prhs[4])
#define COLSCA     (prhs[5])
#define ROWSCA     (prhs[6])
#define RHS        (prhs[7])
#define VAR_SCHUR  (prhs[8])
#define INST       (prhs[9])
#define A_IN       (prhs[10])

#define INFO_OUT  (plhs[0])
#define RINFO_OUT (plhs[1])
#define RHS_OUT   (plhs[2])
#define INST_OUT  (plhs[3])
#define SCHUR_OUT (plhs[4])
#define DEFICIENC (plhs[5])
#define NULLSPACE (plhs[6])
#define PERM_OUT  (plhs[7])
#define UNS_PERM  (plhs[8])


#define MYMALLOC(ptr,l,type)                      \
  if(!(ptr = (type *) malloc(l*sizeof(type)))){   \
    mexErrMsgTxt ("Malloc failed in mumpsmex.c");	\
  }                                               



#define MYFREE(ptr) \
  if(ptr){   \
    free(ptr); \
		ptr = 0;	 \
 }



#define EXTRACT_FROM_MATLAB_TOPTR(mxcomponent,mumpspointer,type,length)          \
  ptr_matlab = mxGetPr(mxcomponent);                                             \
  MYFREE(mumpspointer);			 		         	         \
  if(ptr_matlab[0] != -9999){                                                    \
    MYMALLOC(mumpspointer,length,type);                               		 \
    for(i=0;i<length;i++){                                                       \
      mumpspointer[i] = ptr_matlab[i];                                           \
    }                                                                            \
  }


#define EXTRACT_FROM_MATLAB_TOARR(mxcomponent,mumpsarray,type,length)            \
  ptr_matlab = mxGetPr(mxcomponent);                                             \
  if(ptr_matlab[0] != -9999){                                                    \
    for(i=0;i<length;i++){                                                       \
      if(ptr_matlab[i] != -9998){                                                \
        mumpsarray[i] = ptr_matlab[i];                                           \
      }                                                                          \
    }                                                                            \
  }                                                                              \

#define EXTRACT_FROM_MATLAB_TOVAL(mxcomponent,mumpsvalue)                        \
  ptr_matlab = mxGetPr(mxcomponent);                                             \
  if(ptr_matlab[0] != -9999){                                                    \
      mumpsvalue = ptr_matlab[0];                                                \
  }                                                                              \



#define EXTRACT_FROM_C_TO_MATLAB(mxcomponent,mumpspointer,length)                \
  if(mumpspointer == 0){                                                         \
    mxcomponent = mxCreateDoubleMatrix (1, 1, mxREAL);                           \
    ptr_matlab = mxGetPr (mxcomponent);                                          \
    ptr_matlab[0] = -9999;                                                       \
  }else{                                                                         \
    mxcomponent = mxCreateDoubleMatrix (1,length,mxREAL);                        \
    ptr_matlab = mxGetPr (mxcomponent);                                          \
    for(i=0;i<length;i++){                                                       \
      ptr_matlab[i]=(double)(mumpspointer)[i];                                   \
    }                                                                            \
  }                                                                              \


#define EXTRACT_CMPLX_FROM_MATLAB_TOPTR(mxcomponent,mumpspointer,type,length)    \
  ptr_matlab = mxGetPr(mxcomponent);                                             \
  if(ptr_matlab[0] != -9999){                                                    \
    MYFREE(mumpspointer);	 			         		 \
    MYMALLOC(mumpspointer,length,double2);                                       \
    for(i=0;i<length;i++){                                                       \
      (mumpspointer[i]).r = ptr_matlab[i];				 	 \
    }                                                                            \
    ptr_matlab = mxGetPi(mxcomponent);			 			 \
    if(ptr_matlab){ 								 \
      for(i=0;i<length;i++){							 \
        (mumpspointer[i]).i = ptr_matlab[i];					 \
      }  									 \
    }else{									 \
      for(i=0;i<length;i++){							 \
        (mumpspointer[i]).i = 0.0;						 \
      }										 \
    }	 									 \
  }										 \


#define EXTRACT_CMPLX_FROM_C_TO_MATLAB(mxcomponent,mumpspointer,length)         \
  if(mumpspointer == 0){                                                        \
    mxcomponent = mxCreateDoubleMatrix (1, 1, mxCOMPLEX);                       \
    ptr_matlab = mxGetPr (mxcomponent);                                         \
    ptr_matlab[0] = -9999;							\
    ptr_matlab = mxGetPi (mxcomponent);                                         \
    ptr_matlab[0] = -9999;							\
  }else{                                                                        \
    mxcomponent = mxCreateDoubleMatrix (1,length,mxCOMPLEX);                    \
    ptr_matlab = mxGetPr (mxcomponent);                                         \
    for(i=0;i<length;i++){                                                      \
      ptr_matlab[i] = (mumpspointer[i]).r;					\
    }                                                                           \
    ptr_matlab = mxGetPi (mxcomponent);                                         \
    for(i=0;i<length;i++){							\
      ptr_matlab[i] = (mumpspointer[i]).i;					\
    }										\
  }                                                            

void DMUMPS_free(DMUMPS_STRUC_C **dmumps_par){
  if(*dmumps_par){
  MYFREE( (*dmumps_par)->irn );
  MYFREE( (*dmumps_par)->jcn  );
  MYFREE( (*dmumps_par)->a );
  MYFREE( (*dmumps_par)->irn_loc );
  MYFREE( (*dmumps_par)->jcn_loc );
  MYFREE( (*dmumps_par)->a_loc );
  MYFREE( (*dmumps_par)->eltptr );
  MYFREE( (*dmumps_par)->eltvar );
  MYFREE( (*dmumps_par)->a_elt );
  MYFREE( (*dmumps_par)->perm_in );
  MYFREE( (*dmumps_par)->colsca );
  MYFREE( (*dmumps_par)->rowsca  );
  MYFREE( (*dmumps_par)->nullspace );
  MYFREE( (*dmumps_par)->listvar_schur );
  MYFREE( (*dmumps_par)->sym_perm );
  MYFREE( (*dmumps_par)->uns_perm );
  MYFREE( (*dmumps_par)->irhs_ptr);
  MYFREE( (*dmumps_par)->irhs_sparse);
  MYFREE( (*dmumps_par)->rhs_sparse);
  MYFREE( (*dmumps_par)->rhs);
  MYFREE(*dmumps_par);
  }
}

void DMUMPS_alloc(DMUMPS_STRUC_C **dmumps_par){

  MYMALLOC((*dmumps_par),1,DMUMPS_STRUC_C);
  (*dmumps_par)->irn  = NULL;
  (*dmumps_par)->jcn  = NULL;
  (*dmumps_par)->a  = NULL;
  (*dmumps_par)->irn_loc  = NULL;
  (*dmumps_par)->jcn_loc  = NULL;
  (*dmumps_par)->a_loc  = NULL;
  (*dmumps_par)->eltptr  = NULL;
  (*dmumps_par)->eltvar  = NULL;
  (*dmumps_par)->a_elt  = NULL;
  (*dmumps_par)->perm_in  = NULL;
  (*dmumps_par)->colsca  = NULL;
  (*dmumps_par)->rowsca  = NULL;
  (*dmumps_par)->rhs  = NULL;
  (*dmumps_par)->rhs_sparse = NULL;
  (*dmumps_par)->irhs_sparse = NULL;
  (*dmumps_par)->irhs_ptr = NULL;
  (*dmumps_par)->nullspace  = NULL;
  (*dmumps_par)->listvar_schur  = NULL;
  (*dmumps_par)->schur  = NULL;
  (*dmumps_par)->sym_perm  = NULL;
  (*dmumps_par)->uns_perm  = NULL;
}



void mexFunction(int nlhs, mxArray *plhs[ ],
		 int nrhs, const mxArray *prhs[ ]) { 
  
  int i,j,pos;
  int *ptr_int;
  double *ptr_double;
  double *ptr_matlab;
  double *ptri_matlab;

  /* C pointer for input parameters */
  size_t inst_address;
  int n,m,ne,netrue,inst,job;
  int *irn_in,*jcn_in;
  
  /* variable for multiple and sparse rhs */
  int nbrhs,ldrhs,nz_rhs,posrhs;
  int *irhs_ptr, *irhs_sparse;
  double *rhs_sparse;
  double *im_rhs_sparse;

  DMUMPS_STRUC_C *dmumps_par;
  int dosolve = 0;
  int doanal = 0;
  
  
  EXTRACT_FROM_MATLAB_TOVAL(JOB,job);


  dosolve = (job == 3 || job == 5 || job == 6);
  doanal = (job == 1 || job == 4 || job == 6);

  if(job == -1){
    DMUMPS_alloc(&dmumps_par);
    EXTRACT_FROM_MATLAB_TOVAL(SYM,dmumps_par->sym);
    dmumps_par->job = -1;
    dmumps_par->par = 1;
    dmumps_c(dmumps_par);
    dmumps_par->nz = -1;
    dmumps_par->nz_alloc = -1;
    EXTRACT_CMPLX_FROM_C_TO_MATLAB(SCHUR_OUT,dmumps_par->schur,1);
  }else{
    EXTRACT_FROM_MATLAB_TOVAL(INST,inst_address);
    ptr_int = (int *) inst_address;
    
    dmumps_par = (DMUMPS_STRUC_C *) ptr_int;
    
    if(job == -2){
      dmumps_par->job = -2;
      dmumps_c(dmumps_par);
      DMUMPS_free(&dmumps_par);
    }else{
      
      /* check of input arguments */
      n = mxGetN(A_IN);
      m = mxGetM(A_IN);
			
      if (!mxIsSparse(A_IN) || n != m ){
	mexErrMsgTxt ("Input matrix must be a sparse square matrix");
      }
      jcn_in = mxGetJc(A_IN);
      ne = jcn_in[n];
      irn_in = mxGetIr(A_IN);
      dmumps_par->n = n;
      
      if(dmumps_par->sym != 0){
	netrue = (n+ne)/2;
      }else{
	netrue = ne;
      }
      if(dmumps_par->nz_alloc < netrue || dmumps_par->nz_alloc >= 2*netrue){  
	MYFREE(dmumps_par->jcn);
	MYFREE(dmumps_par->irn);
	MYFREE(dmumps_par->a);
	MYMALLOC((dmumps_par->jcn),netrue,int);
	MYMALLOC((dmumps_par->irn),netrue,int);
	MYMALLOC((dmumps_par->a),netrue,double2);
	dmumps_par->nz_alloc = netrue;
      }


      if(dmumps_par->sym == 0){
	/* if analysis already performed then we only need to read
	   numerical values
	   Note that we suppose that matlab did not change the internal
	   format of the matrix between the 2 calls */
	if(doanal){ 
	  /* || dmumps_par->info[22] == 0 */
	  for(i=0;i<n;i++){
	    for(j=jcn_in[i];j<jcn_in[i+1];j++){
	      (dmumps_par->jcn)[j] = i+1;
	      (dmumps_par->irn)[j] = irn_in[j]+1;
	    }
	  }
	}
	ptr_matlab = mxGetPr(A_IN);
	for(i=0;i<ne;i++){                                                   
	  ((dmumps_par->a)[i]).r = ptr_matlab[i];
	}
	ptr_matlab = mxGetPi(A_IN);
	if(ptr_matlab){
	  for(i=0;i<ne;i++){                                                   
	    ((dmumps_par->a)[i]).i = ptr_matlab[i];
	  }
	}else{
	  for(i=0;i<ne;i++){                                                   
	    ((dmumps_par->a)[i]).i = 0.0;
	  }
	}
	dmumps_par->nz = ne;
      }else{
	/* in the symmetric case we do not need to check doanal */
	pos = 0;
	ptr_matlab = mxGetPr(A_IN);
	ptri_matlab = mxGetPi(A_IN);
	for(i=0;i<n;i++){
	  for(j=jcn_in[i];j<jcn_in[i+1];j++){
	    if(irn_in[j] >= i){
	      if(pos >= netrue){
		mexErrMsgTxt ("Input matrix must be symmetric");}
	      (dmumps_par->jcn)[pos] = i+1;
	      (dmumps_par->irn)[pos] = irn_in[j]+1;
	      ((dmumps_par->a)[pos]).r = ptr_matlab[j];
	      if(ptri_matlab){
		((dmumps_par->a)[pos]).i = ptri_matlab[j];
	      }else{
		((dmumps_par->a)[pos]).i = 0.0;
	      }
	      pos++;
	    }
	  }
	}
	dmumps_par->nz = pos;
      }
    

      EXTRACT_FROM_MATLAB_TOVAL(JOB,dmumps_par->job);
      EXTRACT_FROM_MATLAB_TOARR(ICNTL,dmumps_par->icntl,int,40);
      EXTRACT_FROM_MATLAB_TOARR(CNTL,dmumps_par->cntl,double,5);
      EXTRACT_FROM_MATLAB_TOPTR(PERM_IN,(dmumps_par->perm_in),int,n);

      EXTRACT_CMPLX_FROM_MATLAB_TOPTR(COLSCA,(dmumps_par->colsca),double,n);
      EXTRACT_CMPLX_FROM_MATLAB_TOPTR(ROWSCA,(dmumps_par->rowsca),double,n);

      dmumps_par->size_schur = mxGetN(VAR_SCHUR);
      EXTRACT_FROM_MATLAB_TOPTR(VAR_SCHUR,(dmumps_par->listvar_schur),int,dmumps_par->size_schur);
      if(!dmumps_par->listvar_schur) dmumps_par->size_schur = 0;

      ptr_matlab = mxGetPr (RHS);
      if((!dosolve) || ptr_matlab[0] == -9999){ /* rhs not already provided */
        EXTRACT_CMPLX_FROM_MATLAB_TOPTR(RHS,(dmumps_par->rhs),double,1);
      }else{
	nbrhs = mxGetN(RHS);
	ldrhs = mxGetM(RHS);
	dmumps_par->nrhs = nbrhs;
	dmumps_par->lrhs = ldrhs;
	if(ldrhs != n){
	  mexErrMsgTxt ("Incompatible number of row in RHS");
	}
	if (!mxIsSparse(RHS)){ /* full rhs */
	  dmumps_par->icntl[19] = 0;
	  EXTRACT_CMPLX_FROM_MATLAB_TOPTR(RHS,(dmumps_par->rhs),double,(nbrhs*ldrhs));
	}else{ /* sparse rhs */
	  printf("sparse RHS ldrhs = %d nrhs = %d\n",ldrhs,nbrhs);
	  dmumps_par->icntl[19] = 1;
          irhs_ptr = mxGetJc(RHS);
	  irhs_sparse = mxGetIr(RHS);
	  rhs_sparse = mxGetPr(RHS);
	  im_rhs_sparse = mxGetPi(RHS);

	  nz_rhs = irhs_ptr[nbrhs];

	  MYMALLOC((dmumps_par->irhs_ptr),(nbrhs+1),int);
	  MYMALLOC((dmumps_par->irhs_sparse),nz_rhs,int);
	  MYMALLOC((dmumps_par->rhs_sparse),nz_rhs,double2);
	  MYMALLOC((dmumps_par->rhs),(nbrhs*ldrhs),double2);

	  for(i=0;i<nbrhs;i++){
	    for(j=irhs_ptr[i];j<irhs_ptr[i+1];j++){
	      (dmumps_par->irhs_sparse)[j] = irhs_sparse[j]+1;
	    }
	    (dmumps_par->irhs_ptr)[i] = irhs_ptr[i]+1;
	  }
	  (dmumps_par->irhs_ptr)[nbrhs] = nz_rhs+1;
	  if(im_rhs_sparse){
	    for(i=0;i<nz_rhs;i++){                                                   
	      ((dmumps_par->rhs_sparse)[i]).r = rhs_sparse[i];
	      ((dmumps_par->rhs_sparse)[i]).i = im_rhs_sparse[i];
	    }
	  }else{
	    for(i=0;i<nz_rhs;i++){                                                   
	      ((dmumps_par->rhs_sparse)[i]).r = rhs_sparse[i];
	      ((dmumps_par->rhs_sparse)[i]).i = 0.0;
	    }
	  }
	}
      }
      


      if(dmumps_par->size_schur > 0){
	MYMALLOC((dmumps_par->schur),((dmumps_par->size_schur)*(dmumps_par->size_schur)),double2);
	dmumps_par->icntl[18] = 1;
      }else{
	dmumps_par->icntl[18] = 0;
      }
      dmumps_c(dmumps_par);
    }
  }
  if(nlhs > 0){
    EXTRACT_FROM_C_TO_MATLAB( INFO_OUT  ,(dmumps_par->infog),40);
    EXTRACT_FROM_C_TO_MATLAB( RINFO_OUT ,(dmumps_par->rinfog),20);
    if(dmumps_par->rhs && dosolve){
      RHS_OUT = mxCreateDoubleMatrix (n,nbrhs,mxREAL2);
      ptr_matlab = mxGetPr (RHS_OUT);
      ptri_matlab = mxGetPi (RHS_OUT);
      for(j=0;j<nbrhs;j++){
	posrhs = j*n;
	for(i=0;i<n;i++){
	  ptr_matlab[posrhs+i]= (dmumps_par->rhs)[posrhs+i].r;	
	  ptri_matlab[posrhs+i]= (dmumps_par->rhs)[posrhs+i].i;		
	}              
      }
    }else{
      EXTRACT_CMPLX_FROM_C_TO_MATLAB( RHS_OUT,(dmumps_par->rhs),1);
    }

    ptr_int = (int *)dmumps_par;
    inst_address = (size_t) ptr_int;
    EXTRACT_FROM_C_TO_MATLAB( INST_OUT   ,&inst_address,1); 
    EXTRACT_FROM_C_TO_MATLAB( DEFICIENC  ,&(dmumps_par->deficiency),1);
    EXTRACT_CMPLX_FROM_C_TO_MATLAB( NULLSPACE  ,dmumps_par->nullspace,dmumps_par->deficiency);
    EXTRACT_FROM_C_TO_MATLAB( PERM_OUT   ,dmumps_par->sym_perm,n);
    EXTRACT_FROM_C_TO_MATLAB( UNS_PERM   ,dmumps_par->uns_perm,n);

    if(dmumps_par->size_schur > 0){
      SCHUR_OUT = mxCreateDoubleMatrix(dmumps_par->size_schur,dmumps_par->size_schur,mxREAL2);
      ptr_matlab = mxGetPr (SCHUR_OUT);
      ptri_matlab = mxGetPi (SCHUR_OUT);
      for(i=0;i<dmumps_par->size_schur;i++){
        pos = i*(dmumps_par->size_schur);
	for(j=0;j<dmumps_par->size_schur;j++){
	  ptr_matlab[j+pos] = ((dmumps_par->schur)[j+pos]).r;
	  ptri_matlab[j+pos] = ((dmumps_par->schur)[j+pos]).i;
	}
      }
    }else{
      SCHUR_OUT = mxCreateDoubleMatrix(1,1,mxREAL2);
      ptr_matlab = mxGetPr (SCHUR_OUT);
      ptr_matlab[0] = -9999; 
      ptr_matlab = mxGetPi (SCHUR_OUT);
      ptr_matlab[0] = -9999;
    }
    MYFREE(dmumps_par->schur);
    MYFREE(dmumps_par->irhs_ptr);
    MYFREE(dmumps_par->irhs_sparse);
    MYFREE(dmumps_par->rhs_sparse);
    MYFREE(dmumps_par->rhs);
  }
}
