#include "mex.h"
#include "stack-c.h"
#include <stdlib.h>
#include <stdio.h>
#include "zmumps_c.h"

#define dmumps_c      zmumps_c
#define dmumps_par    zmumps_par
#define DMUMPS_STRUC_C  ZMUMPS_STRUC_C
#define DMUMPS_alloc    ZMUMPS_alloc     
#define DMUMPS_free     ZMUMPS_free
#define double2         mumps_double_complex
#define nb_RHS 11

#define MYFREE(ptr)\
if(ptr){	\
  free(ptr);  \
  ptr=0;} 	\
 

#define EXTRACT_FROM_SCILAB_TOPTR(it,ptr_scilab1,ptr_scilab2,mumpspointer,type,length)	\
if(ptr_scilab1[0] != -9999){                                            		\
  free(mumpspointer);                                                   		\
  mumpspointer = (type *) malloc(length*sizeof(type));  				\
  for(i=0;i<length;i++){                                        			\
  mumpspointer[i] = ptr_scilab1[i];                             			\
  }                                                             			\
}						\

#define EXTRACT_FROM_SCILAB_TOARR(ptr_scilab,mumpsarray,type,length)    \
if(ptr_scilab[0] != -9999){                                            	\
  for(i=0;i<length;i++){                                                \
    if(ptr_scilab[i] != -9998){                                         \
      mumpsarray[i] = ptr_scilab[i];                                    \
      }                                                                 \
   }                                                                    \
   }\

#define EXTRACT_INT_FROM_C_TO_SCILAB(num,ptr_scilab,mumpspointer,length1,length2,one) 	\
if(mumpspointer == 0){  								\
CreateVar(nb_RHS+num,"i",&one,&one,&ptr_scilab);        				\
*istk(ptr_scilab)=-9999;                                				\
}else{                                                  				\
CreateVar(nb_RHS+num,"i",&length1,&length2,&ptr_scilab);				\
int l=length1*length2;  								\
for (i=0;i<l;i++){ 									\
   *istk(ptr_scilab+i)=0;                                  				\
   *istk(ptr_scilab+i)=(mumpspointer)[i];}                      			\
}                                                   				\
LhsVar(num)=nb_RHS+num;                                 			\
    
#define EXTRACT_DOUBLE_FROM_C_TO_SCILAB(num,it,ptr_scilab1,ptr_scilab2,mumpspointer,length1,length2,one)\
if(mumpspointer == 0){                                          					\
CreateVar(nb_RHS+num,"d",&one,&one,&ptr_scilab1);                					\
*stk(ptr_scilab1)=-9999;                                         					\
}else{                                                          					\
CreateVar(nb_RHS+num,"d",&length1,&length2,&ptr_scilab1); 						\
int l=length1*length2;  										\
for (i=0;i<l;i++){                                      						\
   *stk(ptr_scilab1+i)=(mumpspointer)[i];}                       					\
}                                                           					\
LhsVar(num)=nb_RHS+num; 								\
                                                                                

#define EXTRACT_CMPLX_FROM_SCILAB_TOPTR(it,ptr_scilab1,ptr_scilab2,mumpspointer,type,length)    \
if(ptr_scilab1[0] != -9999){                                                    		\
  free(mumpspointer);	 			         					\
  mumpspointer=(double2 *) malloc(length*sizeof(double2));      					\
  for(i=0;i<length;i++){                                                			\
    (mumpspointer[i]).r = ptr_scilab1[i];} 							\
  if(it==1){ 											\
    for(i=0;i<length;i++){									\
      (mumpspointer[i]).i = ptr_scilab2[i];}							\
  }else{											\
    for(i=0;i<length;i++){									\
      (mumpspointer[i]).i = 0.0;}								\
    }	 											\
  }											\

#define EXTRACT_CMPLX_FROM_C_TO_SCILAB(num,it,ptr_scilab1,ptr_scilab2,mumpspointer,length1,length2,one) \
  if(it!=1){												\
    Scierror(999,"Internal error, the variable built must be complex.");}  				\
  if(mumpspointer == 0){                                                        			\
    CreateCVar(nb_RHS+num,"d",&it, &one,&one,&ptr_scilab1,&ptr_scilab2);         			\
    *stk(ptr_scilab1) = -9999;  									\
    *stk(ptr_scilab2) = -9999;  									\
  }else{												\
    int l=length1*length2;										\
    CreateCVar(nb_RHS+num,"d",&it,&length1,&length2,&ptr_scilab1,&ptr_scilab2);          		\
    for(i=0;i<l;i++){                                                      				\
      *stk(ptr_scilab1+i) = (mumpspointer[i]).r; 							\
    }                                                                           			\
    for(i=0;i<l;i++){											\
      *stk(ptr_scilab2+i) = (mumpspointer[i]).i;							\
    }													\
  }													\
  LhsVar(num)=nb_RHS+num;                                                       	\

     


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

  *dmumps_par = (DMUMPS_STRUC_C *) malloc(sizeof(DMUMPS_STRUC_C));
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
  (*dmumps_par)->irhs_ptr  = NULL;
  (*dmumps_par)->irhs_sparse  = NULL;
  (*dmumps_par)->rhs_sparse  = NULL;
  (*dmumps_par)->nullspace  = NULL;
  (*dmumps_par)->listvar_schur  = NULL;
  (*dmumps_par)->schur  = NULL;
  (*dmumps_par)->sym_perm  = NULL;
  (*dmumps_par)->uns_perm  = NULL;
}



 static int dmumpsc(char *fname){

 
  /* RhsVar parameters */
  int njob, mjob, ljob, mint, nint, lint, nsym, msym, lsym, nA, mA, nRHS, mRHS,lRHS, liRHS;
  int nicntl, micntl, licntl, ncntl, mcntl, lcntl, nperm, mperm, lperm;
  int ncols, mcols, lcols, licols, nrows, mrows, lrows, lirows, ns_schu , ms_schu, ls_schu;
  int nv_schu, mv_schu, lv_schu, nschu, mschu, lschu;
  int type_rhs, mtype_rhs, ntype_rhs, ltype_rhs;

  /* LhsVar parameters */
  int linfog, lrinfog, lrhsout,lrhsouti, linstout, lschurout, lschurouti, ldef;
  int lnullsp, lnullspi, lmapp, lsymperm, lunsperm;
  int one=1, temp1=40, temp2=20, temp3, temp4, temp5;
  int it, itRHS; /* parameter for real/complex types */
  
  int i,j,k1,k2, nb_in_row,netrue;
  int *ptr_int;
  double *ptr_double;
  double *ptr_scilab;
  double * ptri_scilab;

  /* C pointer for input parameters */
  int inst_adress;
  int ne,inst;
  int *irn_in,*jcn_in; 
 
  /* Variable for multiple and sparse RHS*/
  int posrhs, posschur, nz_RHS,col_ind,k;
  int *irhs_ptr;
  int *irhs_sparse;
  double *rhs_sparse;
  double *im_rhs_sparse;

  SciSparse A;
  SciSparse RHS_SPARSE;
  DMUMPS_STRUC_C *dmumps_par;
  
  int dosolve=0;
  int doanal = 0;
  /* Check number of input parameters */
  CheckRhs(10,11);

  /* Get job value. njob/mjob are the dimensions of variable job. */
  GetRhsVar(2,"i",&mjob,&njob,&ljob);
  dosolve = (*istk(ljob) == 3 || *istk(ljob) == 5 ||*istk(ljob) == 6);
  doanal = (*istk(ljob) == 1 || *istk(ljob) == 4 || *istk(ljob) == 6);
  if(*istk(ljob) == -1){

    DMUMPS_alloc(&dmumps_par);
    GetRhsVar(1,"i",&msym,&nsym,&lsym);
    dmumps_par->sym=*istk(lsym);
    dmumps_par->job = -1;
    dmumps_par->par = 1;
    dmumps_c(dmumps_par);
    dmumps_par->nz = -1;
    dmumps_par->nz_alloc=-1;
    it=1;
    EXTRACT_CMPLX_FROM_C_TO_SCILAB(5,it,lschurout,lschurouti,(dmumps_par->schur),one,one,one);    

  }else{
    /* Obtain pointer on instance */ 
    GetRhsVar(10,"i",&mint,&nint,&lint);
    inst_adress=*istk(lint); /* EXTRACT_FROM_SCILAB_TOVAL(INST,inst_address); */
    ptr_int = (int *) inst_adress;
    
    dmumps_par = (DMUMPS_STRUC_C *) ptr_int;
    if(*istk(ljob) == -2){
      dmumps_par->job = -2;
      dmumps_c(dmumps_par);
      DMUMPS_free(&dmumps_par);
    }else{
      /* Get the sparse matrix A */
      GetRhsVar(11,"s",&mA,&nA,&A);

      if (nA != mA || mA<1 ){
	Scierror(999,"%s : Bad dimensions for the matrix",fname, 12);
	return 0;
      }
      
      ne=A.nel;
      dmumps_par->n = nA;
	
      if(dmumps_par->sym != 0){
	netrue = (nA+ne)/2;
      }else{
	netrue = ne;
      }

      if(dmumps_par->nz_alloc < netrue ||dmumps_par->nz_alloc >= 2*netrue){  
	MYFREE(dmumps_par->jcn);
	MYFREE(dmumps_par->irn);
	MYFREE(dmumps_par->a);                                         
	
	dmumps_par->jcn = (int*)malloc(netrue*sizeof(int));
	dmumps_par->irn = (int*)malloc(netrue*sizeof(int));
	dmumps_par->a = (double2 *) malloc(netrue*sizeof(double2));
	dmumps_par->nz_alloc = netrue;
      }
      /* Check for symmetry in order to initialize only
       * lower triangle on entry to symmetric MUMPS code */
      if ((dmumps_par->sym)==0){
        /*
	 * Unsymmetric case:
	 * build irn from mnel for MUMPS format
	 * mA : number of rows
	 */

        if(doanal){ 
	  for(i=0;i<ne;i++){
	    (dmumps_par->jcn)[i]=(A.icol)[i];}
	  k1=0;
	  for (k2=1;k2<mA+1;k2++){
	    nb_in_row=0; 
	    while(nb_in_row<(A.mnel)[k2-1]){
	      dmumps_par->irn[k1]=k2; /* matrix indices start at 1 */
	      k1=k1+1;
	      nb_in_row=nb_in_row+1;
	    } 
	  }
        }

	for(i=0;i<ne;i++){
          ((dmumps_par->a)[i]).r = (A.R)[i];}
        if(A.it == 1){
          for(i=0;i<ne;i++){
            ((dmumps_par->a)[i]).i = (A.I)[i];}
        }else{
          for(i=0;i<ne;i++){
            ((dmumps_par->a)[i]).i = 0.0;}
        }
	dmumps_par->nz = ne; 
	}
      else{
	/* symmetric case */
	k1=0;
        i=0;
	for (k2=1;k2<mA+1;k2++){
	  nb_in_row=0; 
	  while(nb_in_row<(A.mnel)[k2-1]){
            if( k2 >= (A.icol)[i]){
	      if(k1>=netrue){	
	 	Scierror(999," The matrix must be symmetric ",fname, 13);
	  	return 0;
	      }
              (dmumps_par->jcn)[k1]=(A.icol)[i];
	      (dmumps_par->irn)[k1]=k2; 
	      (dmumps_par->a)[k1].r=(A.R)[i];
              if(A.it == 1){
                ((dmumps_par->a)[k1]).i = (A.I)[i];}
              else{
                ((dmumps_par->a)[k1]).i = 0.0;}
	      k1=k1+1;}
		
	      nb_in_row=nb_in_row+1;
	      i=i+1;
	     } 
	  }     
	dmumps_par->nz = k1;
      	}  
	
        GetRhsVar(2,"i",&mjob,&njob,&ljob);
	dmumps_par->job=*istk(ljob);
	
	GetRhsVar(3,"i",&micntl,&nicntl,&licntl);
	EXTRACT_FROM_SCILAB_TOARR(istk(licntl),dmumps_par->icntl,int,40);

	GetRhsVar(4,"d",&mcntl,&ncntl,&lcntl);
	EXTRACT_FROM_SCILAB_TOARR(stk(lcntl),dmumps_par->cntl,double,5);
	
        GetRhsVar(5,"i",&mperm, &nperm, &lperm);
	EXTRACT_FROM_SCILAB_TOPTR(A.it,istk(lperm),istk(lperm),(dmumps_par->perm_in),int,nA);

	GetRhsCVar(6,"d",&it,&mcols,&ncols,&lcols,&licols);
        EXTRACT_CMPLX_FROM_SCILAB_TOPTR(it,stk(lcols),stk(licols),(dmumps_par->colsca),double2,nA);
	
        GetRhsCVar(7,"d",&it,&mrows,&nrows,&lrows,&lirows);
        EXTRACT_CMPLX_FROM_SCILAB_TOPTR(it,stk(lrows),stk(lirows),(dmumps_par->rowsca),double2,nA);
       
	
        if(GetType(8)!=5){
          GetRhsCVar(8,"d",&itRHS,&mRHS,&nRHS,&lRHS,&liRHS);
      
          if((!dosolve) || (stk(lRHS)[0]) == -9999){
            EXTRACT_CMPLX_FROM_SCILAB_TOPTR(itRHS,stk(lRHS),stk(liRHS),(dmumps_par->rhs),double2,one);
          }else{
  
            dmumps_par->nrhs = nRHS;
            dmumps_par->lrhs = mRHS;
            if(mRHS!=nA){
              Scierror(999,"Incompatible number of rows in RHS",fname,8);}
            dmumps_par->icntl[19]=0;            
            EXTRACT_CMPLX_FROM_SCILAB_TOPTR(itRHS,stk(lRHS),stk(liRHS),(dmumps_par->rhs),double2,(nRHS*mRHS));
          }
        }else{
          GetRhsVar(8,"s",&mRHS,&nRHS,&RHS_SPARSE);
          if((stk(lRHS)[0]) == -9999){
            EXTRACT_CMPLX_FROM_SCILAB_TOPTR(it,stk(lRHS),stk(liRHS),(dmumps_par->rhs),double2,one);}
          dmumps_par->icntl[19]=1;
          dmumps_par->nrhs = nRHS;
          dmumps_par->lrhs = mRHS;
          nz_RHS=RHS_SPARSE.nel;
          
          irhs_ptr=(int*)malloc((nRHS+1)*sizeof(int));
   
          dmumps_par->irhs_ptr=(int*)malloc((nRHS+1)*sizeof(int));
          dmumps_par->irhs_sparse=(int*)malloc(nz_RHS*sizeof(int));
          dmumps_par->rhs_sparse=(double2*)malloc(nz_RHS*sizeof(double2));
          dmumps_par->rhs=(double2*)malloc((nRHS*mRHS)*sizeof(double2));
          /* transform row-oriented sparse multiple rhs (scilab)
	   * into column-oriented sparse multiple rhs (mumps) */
          k=0;
          for(i=0;i<nRHS+1;i++){
            irhs_ptr[i]=0;
            dmumps_par->irhs_ptr[i]=0;}
          for(i=1;i<mRHS+1;i++){
            for(j=0;j<(RHS_SPARSE.mnel)[i-1];j++){
              col_ind=(RHS_SPARSE.icol)[k];
              k++;
              ((dmumps_par->irhs_ptr)[col_ind])++;
            }
          }  
          (dmumps_par->irhs_ptr)[0]=1;
          irhs_ptr[0]=(dmumps_par->irhs_ptr)[0];
          for(i=1;i<nRHS+1;i++){
            (dmumps_par->irhs_ptr)[i]=(dmumps_par->irhs_ptr)[i]+(dmumps_par->irhs_ptr)[i-1];
            irhs_ptr[i]= (dmumps_par->irhs_ptr)[i];
          }  
          k=RHS_SPARSE.nel-1;
          for(i=mRHS;i>=1;i--){
       
            for(j=0;j<(RHS_SPARSE.mnel)[i-1];j++){
              col_ind=(RHS_SPARSE.icol)[k];
             (dmumps_par->irhs_sparse)[irhs_ptr[col_ind]-2]=i;
              if(RHS_SPARSE.it==1){
                ((dmumps_par->rhs_sparse)[irhs_ptr[col_ind]-2]).r=RHS_SPARSE.R[k];
                ((dmumps_par->rhs_sparse)[irhs_ptr[col_ind]-2]).i=RHS_SPARSE.I[k];
              }else{
                ((dmumps_par->rhs_sparse)[irhs_ptr[col_ind]-2]).r=RHS_SPARSE.R[k];
                ((dmumps_par->rhs_sparse)[irhs_ptr[col_ind]-2]).i=0.0;
              }
              k--;
              irhs_ptr[col_ind]=irhs_ptr[col_ind]-1;
            }
          }
 	MYFREE(irhs_ptr);
	}
	
	GetRhsVar(9,"i",&nv_schu,&mv_schu,&lv_schu);
	dmumps_par-> size_schur=mv_schu;
	EXTRACT_FROM_SCILAB_TOPTR(A.it,istk(lv_schu),istk(lv_schu),(dmumps_par->listvar_schur),int,dmumps_par->size_schur);
	if(!dmumps_par->listvar_schur) dmumps_par->size_schur=0;

        if((dmumps_par->size_schur)>0){
	  MYFREE(dmumps_par->schur);	
          if(!(dmumps_par->schur=(double2 *)malloc((dmumps_par->size_schur*dmumps_par->size_schur)*sizeof(double2)))){
            sciprint("Malloc failed in intmumpsc.c");
          }
          dmumps_par->icntl[18]=1;
        }else{
          dmumps_par->icntl[18]=0;
        }
  
	dmumps_c(dmumps_par);
      }
    }
    
    if(*istk(ljob)==-2){
      return 0;
    }else{    
      
      CheckLhs(9,9);    
      
      EXTRACT_INT_FROM_C_TO_SCILAB(1,linfog,(dmumps_par->infog),one,temp1,one);
      
      EXTRACT_DOUBLE_FROM_C_TO_SCILAB(2,it,lrinfog,lrinfog,(dmumps_par->rinfog),one,temp2,one);
      
      if(dmumps_par->rhs && dosolve){ /* Just to know if solution step was called */
	it=1;
        CreateCVar(nb_RHS+3,"d",&it,&nA,&nRHS,&lrhsout,&lrhsouti);
        for(i=0;i<nRHS;i++){
          posrhs=i*nA;
          for(j=0;j<nA;j++){
            *stk(lrhsout+posrhs+j)=((dmumps_par->rhs)[posrhs+j]).r;
            *stk(lrhsouti+posrhs+j)=((dmumps_par->rhs)[posrhs+j]).i;
          }
        }
        LhsVar(3)=nb_RHS+3;  
      }else{
        it=1;
        temp5=(dmumps_par->n)*(dmumps_par->nrhs);
        EXTRACT_CMPLX_FROM_C_TO_SCILAB(3,it,lrhsout,lrhsouti,(dmumps_par->rhs),one,one,one);
      }

      ptr_int = (int *)dmumps_par;
      inst_adress = (int) ptr_int;
      EXTRACT_INT_FROM_C_TO_SCILAB(4,linstout,&inst_adress,one,one,one);
      
      temp4=dmumps_par->size_schur;
      if(temp4>0){
        CreateCVar(nb_RHS+5,"d",&it,&temp4,&temp4,&lschurout,&lschurouti);
        for(i=0;i<temp4;i++){
          posschur=i*temp4;
          for(j=0;j<temp4;j++){
            *stk(lschurout+posschur+j)=((dmumps_par->schur)[posschur+j]).r;
            *stk(lschurouti+posschur+j)=((dmumps_par->schur)[posschur+j]).i;
          }
        }
LhsVar(5)=nb_RHS+5;
      }else{
        it=1;
          EXTRACT_CMPLX_FROM_C_TO_SCILAB(5,it,lschurout,lschurouti,(dmumps_par->schur),one,one,one);
      }
      
      MYFREE(dmumps_par->schur);
      MYFREE(dmumps_par->irhs_ptr);
      MYFREE(dmumps_par->irhs_sparse);
      MYFREE(dmumps_par->rhs_sparse);
      MYFREE(dmumps_par->rhs);
      
      EXTRACT_INT_FROM_C_TO_SCILAB(6,ldef, &(dmumps_par->deficiency),one,one,one);       
      
      temp3=dmumps_par->deficiency;
      EXTRACT_CMPLX_FROM_C_TO_SCILAB(7,it,lnullsp,lnullspi,(dmumps_par->nullspace),one,temp3,one);       
      
      EXTRACT_INT_FROM_C_TO_SCILAB(8,lsymperm,(dmumps_par->sym_perm),one,nA,one);   
      
      EXTRACT_INT_FROM_C_TO_SCILAB(9,lunsperm,(dmumps_par->uns_perm),one,nA,one); 
          
      return 0;
         
    }
}


static GenericTable Tab[]={
{(Myinterfun) sci_gateway, dmumpsc,"zmumpsc"}
};

int C2F(scizmumps)()
{Rhs = Max(0, Rhs);
(*(Tab[Fin-1].f))(Tab[Fin-1].name,Tab[Fin-1].F);
return 0;
}
