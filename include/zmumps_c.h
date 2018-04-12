/*

   THIS FILE IS PART OF MUMPS VERSION 4.6.3
   This Version was built on Thu Jun 22 13:22:44 2006


  This version of MUMPS is provided to you free of charge. It is public
  domain, based on public domain software developed during the Esprit IV
  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. 
  Since this first public domain version in 1999, the developments are
  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
  INRIA.

  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.

  Up-to-date copies of the MUMPS package can be obtained
  from the Web pages http://www.enseeiht.fr/apo/MUMPS/
  or http://graal.ens-lyon.fr/MUMPS


   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.


  User documentation of any code that uses this software can
  include this complete notice. You can acknowledge (using
  references [1], [2], and [3] the contribution of this package
  in any scientific publication dependent upon the use of the
  package. You shall use reasonable endeavours to notify
  the authors of the package of this publication.

   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent (1998),
   Multifrontal parallel distributed symmetric and unsymmetric solvers,
   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).

   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
   A fully asynchronous multifrontal solver using distributed dynamic
   scheduling, SIAM Journal of Matrix Analysis and Applications,
   Vol 23, No 1, pp 15-41 (2001).

   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
   S. Pralet, Hybrid scheduling for the parallel solution of linear
   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).

*/
/* $Id: zmumps_c.h,v 1.28 2006/06/13 11:06:00 jylexcel Exp $ */
/* Mostly written in march 2002 (JYL) */

#if ( ! defined ZMUMPS_C_H )
#define ZMUMPS_C_H

/* Complex datatypes */
typedef struct {float r,i;} zmumps_complex;
typedef struct {double r,i;} mumps_double_complex;

/* Next line defines F_INT, F_DOUBLE and F_DOUBLE2 */
#include "zmumps_prec.h"
/*
 * Definition of the (simplified)
 * MUMPS C structure
 */
typedef struct
  {
    F_INT sym, par, job;
    F_INT comm_fortran;    /* Fortran communicator */
    F_INT icntl[40];
    F_DOUBLE2 cntl[5];
    F_INT n;
   
    F_INT nz_alloc; /* used in matlab interface to decide if
                       we free + malloc when we have large variation */

    /* Assembled entry */
    F_INT nz; F_INT *irn; F_INT *jcn; F_DOUBLE *a;
    /* Distributed entry */
    F_INT nz_loc; F_INT *irn_loc; F_INT *jcn_loc; F_DOUBLE *a_loc;
    /* Element entry */
    F_INT nelt; F_INT *eltptr; F_INT *eltvar; F_DOUBLE *a_elt;

    /* Ordering, if given by user */
    F_INT *perm_in;

    /* Orderings returned to user */
    /* symmetric permutation */
    F_INT *sym_perm;
    /* column permutation */
    F_INT *uns_perm;

    /* Scaling (input only in this version) */
    F_DOUBLE *colsca; F_DOUBLE *rowsca;
    /* RHS, solution, ouptput data and statistics */
    F_DOUBLE *rhs, *rhs_sparse, *sol_loc;
    F_INT *irhs_sparse, *irhs_ptr, *isol_loc;
    F_INT nrhs, lrhs, nz_rhs, lsol_loc;
  F_INT schur_mloc, schur_nloc, schur_lld;
  F_INT mblock, nblock, nprow, npcol;
    F_INT info[40],infog[40];
    F_DOUBLE2 rinfo[20], rinfog[20];
    /* Null space */
    F_INT deficiency; F_DOUBLE * nullspace; F_INT * mapping;
    /* Schur */
    F_INT size_schur; F_INT *listvar_schur; F_DOUBLE *schur;
    /* Internal parameters */
    F_INT instance_number;
    /* For out-of-core */
    char ooc_tmpdir[151];
    char ooc_prefix[151];
  } ZMUMPS_STRUC_C;


#if defined(UPPER) || defined(_WIN32)
#define zmumps_f77_ ZMUMPS_F77
#define zmumps_affect_mapping_ ZMUMPS_AFFECT_MAPPING
#define zmumps_affect_nullspace_ ZMUMPS_AFFECT_NULLSPACE
#define zmumps_affect_colsca_ ZMUMPS_AFFECT_COLSCA
#define zmumps_affect_rowsca_ ZMUMPS_AFFECT_ROWSCA 
#define zmumps_affect_uns_perm_     ZMUMPS_AFFECT_UNS_PERM
#define zmumps_affect_sym_perm_     ZMUMPS_AFFECT_SYM_PERM
#define zmumps_nullify_c_mapping_   ZMUMPS_NULLIFY_C_MAPPING
#define zmumps_nullify_c_nullspace_ ZMUMPS_NULLIFY_C_NULLSPACE
#define zmumps_nullify_c_sym_perm_  ZMUMPS_NULLIFY_C_SYM_PERM
#define zmumps_nullify_c_uns_perm_  ZMUMPS_NULLIFY_C_UNS_PERM
#define zmumps_nullify_c_colsca_    ZMUMPS_NULLIFY_C_COLSCA
#define zmumps_nullify_c_rowsca_    ZMUMPS_NULLIFY_C_ROWSCA
#elif defined(Add__)
#define zmumps_f77_ zmumps_f77__
#define zmumps_affect_mapping_ zmumps_affect_mapping__
#define zmumps_affect_nullspace_ zmumps_affect_nullspace__
#define zmumps_affect_colsca_ zmumps_affect_colsca__
#define zmumps_affect_rowsca_ zmumps_affect_rowsca__
#define zmumps_affect_uns_perm_     zmumps_affect_uns_perm__     
#define zmumps_affect_sym_perm_     zmumps_affect_sym_perm__     
#define zmumps_nullify_c_mapping_   zmumps_nullify_c_mapping__    
#define zmumps_nullify_c_nullspace_ zmumps_nullify_c_nullspace__  
#define zmumps_nullify_c_sym_perm_  zmumps_nullify_c_sym_perm__   
#define zmumps_nullify_c_uns_perm_  zmumps_nullify_c_uns_perm__   
#define zmumps_nullify_c_colsca_    zmumps_nullify_c_colsca__     
#define zmumps_nullify_c_rowsca_    zmumps_nullify_c_rowsca__    
#elif defined(Add_)
/* Default. Nothing to do */
#else
/* Name without underscore is used */
#define zmumps_f77_ zmumps_f77
#define zmumps_affect_mapping_ zmumps_affect_mapping
#define zmumps_affect_nullspace_ zmumps_affect_nullspace
#define zmumps_affect_colsca_ zmumps_affect_colsca
#define zmumps_affect_rowsca_ zmumps_affect_rowsca
#define zmumps_affect_uns_perm_     zmumps_affect_uns_perm     
#define zmumps_affect_sym_perm_     zmumps_affect_sym_perm     
#define zmumps_nullify_c_mapping_   zmumps_nullify_c_mapping    
#define zmumps_nullify_c_nullspace_ zmumps_nullify_c_nullspace  
#define zmumps_nullify_c_sym_perm_  zmumps_nullify_c_sym_perm   
#define zmumps_nullify_c_uns_perm_  zmumps_nullify_c_uns_perm   
#define zmumps_nullify_c_colsca_    zmumps_nullify_c_colsca     
#define zmumps_nullify_c_rowsca_    zmumps_nullify_c_rowsca
#endif

#if defined(_WIN32)
/* 
 * Next line May be needed depending on your Windows environment:
 * #define MUMPS_CALL __stdcall
 */
#else
#define MUMPS_CALL
#endif

void MUMPS_CALL zmumps_c(ZMUMPS_STRUC_C * zmumps_par);
void MUMPS_CALL zmumps_affect_mapping_(F_INT * f77mapping);
void MUMPS_CALL zmumps_affect_nullspace_(F_DOUBLE * f77nullspace);
void MUMPS_CALL zmumps_affect_uns_perm_(F_INT * f77sym_perm);
void MUMPS_CALL zmumps_affect_sym_perm_(F_INT * f77uns_perm);
void MUMPS_CALL zmumps_nullify_c_mapping_();
void MUMPS_CALL zmumps_nullify_c_nullspace_();
void MUMPS_CALL zmumps_nullify_c_sym_perm_();
void MUMPS_CALL zmumps_nullify_c_uns_perm_();
#ifdef return_scaling
void MUMPS_CALL zmumps_affect_colsca_(F_DOUBLE * f77colsca);
void MUMPS_CALL zmumps_affect_rowsca_(F_DOUBLE * f77rowsca);
void MUMPS_CALL zmumps_nullify_c_colsca_();
void MUMPS_CALL zmumps_nullify_c_rowsca_();
#endif

void MUMPS_CALL zmumps_f77_(F_INT *job, F_INT *sym, F_INT *par, F_INT * comm_fortran,
F_INT *n, F_INT *icntl,
F_DOUBLE2 *cntl, F_INT *nz, F_INT *irn, F_INT *irn_avail,
F_INT *jcn, F_INT *jcn_avail, F_DOUBLE *a, F_INT *a_avail,
  F_INT *nz_loc, F_INT * irn_loc, F_INT * irn_loc_avail,
  F_INT * jcn_loc, F_INT * jcn_loc_avail, F_DOUBLE *a_loc,
  F_INT * a_loc_avail, F_INT * nelt, F_INT * eltptr,
  F_INT * eltptr_avail, F_INT * eltvar, F_INT * eltvar_avail,
  F_DOUBLE * a_elt, F_INT * a_elt_avail, F_INT * perm_in,
  F_INT * perm_in_avail, F_DOUBLE * rhs, F_INT * rhs_avail,
  F_INT * info, F_DOUBLE2 * rinfo, F_INT * infog, F_DOUBLE2 * rinfog,
  F_INT * deficiency, F_INT * size_schur, F_INT * listvar_schur,
  F_INT * listvar_schur_avail, F_DOUBLE * schur,
  F_INT * schur_avail, F_DOUBLE * colsca, F_INT * colsca_avail,
  F_DOUBLE * rowsca, F_INT * rowsca_avail, F_INT * instance_number,
  F_INT * nrhs, F_INT * lrhs, F_DOUBLE * rhs_sparse, F_INT * rhs_sparse_avail,
  F_DOUBLE * sol_loc, F_INT * sol_loc_avail, F_INT * irhs_sparse,
  F_INT * irhs_sparse_avail, F_INT * irhs_ptr, F_INT * irhs_ptr_avail,
  F_INT * isol_loc, F_INT * isol_loc_avail, F_INT * nz_rhs, F_INT * lsol_loc
  ,F_INT * schur_mloc, F_INT *schur_nloc, F_INT * schur_lld,
  F_INT * schur_mblock, F_INT * schur_nblock, F_INT * schur_nprow,
  F_INT * schur_npcol,
  F_INT * ooc_tmpdir,
  F_INT * ooc_prefix,
  F_INT * ooc_tmpdirlen,
  F_INT * ooc_prefixlen
  );
#endif

