/*
 *
 *  This file is part of MUMPS 4.9.1, built on Fri Oct 23 14:06:12 UTC 2009
 *
 *
 *  This version of MUMPS is provided to you free of charge. It is public
 *  domain, based on public domain software developed during the Esprit IV
 *  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL.
 *  Since this first public domain version in 1999, the developments are
 *  supported by the following institutions: CERFACS, CNRS, INPT(ENSEEIHT)-
 *  IRIT, and INRIA.
 *
 *  Current development team includes Patrick Amestoy, Alfredo Buttari,
 *  Abdou Guermouche, Jean-Yves L'Excellent, Bora Ucar.
 *
 *  Up-to-date copies of the MUMPS package can be obtained
 *  from the Web pages:
 *  http://mumps.enseeiht.fr/  or  http://graal.ens-lyon.fr/MUMPS
 *
 *
 *   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 *   EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
 *
 *
 *  User documentation of any code that uses this software can
 *  include this complete notice. You can acknowledge (using
 *  references [1] and [2]) the contribution of this package
 *  in any scientific publication dependent upon the use of the
 *  package. You shall use reasonable endeavours to notify
 *  the authors of the package of this publication.
 *
 *   [1] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
 *   A fully asynchronous multifrontal solver using distributed dynamic
 *   scheduling, SIAM Journal of Matrix Analysis and Applications,
 *   Vol 23, No 1, pp 15-41 (2001).
 *
 *   [2] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
 *   S. Pralet, Hybrid scheduling for the parallel solution of linear
 *   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
 *
 */
#include "mumps_common.h"
/* Special case of mapping and pivnul_list -- allocated from MUMPS */
static MUMPS_INT * MUMPS_MAPPING;
static MUMPS_INT * MUMPS_PIVNUL_LIST;
/* as uns_perm and sym_perm */
static MUMPS_INT * MUMPS_SYM_PERM;
static MUMPS_INT * MUMPS_UNS_PERM;
MUMPS_INT*
mumps_get_mapping()
{
    return MUMPS_MAPPING;
}
void MUMPS_CALL
MUMPS_AFFECT_MAPPING(MUMPS_INT * f77mapping)
{
    MUMPS_MAPPING = f77mapping;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_MAPPING()
{
    MUMPS_MAPPING = 0;
}
MUMPS_INT*
mumps_get_pivnul_list()
{
    return MUMPS_PIVNUL_LIST;
}
void MUMPS_CALL
MUMPS_AFFECT_PIVNUL_LIST(MUMPS_INT * f77pivnul_list)
{
    MUMPS_PIVNUL_LIST = f77pivnul_list;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_PIVNUL_LIST()
{
    MUMPS_PIVNUL_LIST = 0;
}
MUMPS_INT*
mumps_get_sym_perm()
{
    return MUMPS_SYM_PERM;
}
void MUMPS_CALL
MUMPS_AFFECT_SYM_PERM(MUMPS_INT * f77sym_perm)
{
    MUMPS_SYM_PERM = f77sym_perm;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_SYM_PERM()
{
    MUMPS_SYM_PERM = 0;
}
MUMPS_INT*
mumps_get_uns_perm()
{
    return MUMPS_UNS_PERM;
}
void MUMPS_CALL
MUMPS_AFFECT_UNS_PERM(MUMPS_INT * f77uns_perm)
{
    MUMPS_UNS_PERM = f77uns_perm;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_UNS_PERM()
{
    MUMPS_UNS_PERM = 0;
}
