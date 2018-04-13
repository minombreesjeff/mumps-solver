/*
 *
 *  This file is part of MUMPS 4.8.0, built on Fri Jul 25 14:46:02 2008
 *
 *
 *  This version of MUMPS is provided to you free of charge. It is public
 *  domain, based on public domain software developed during the Esprit IV
 *  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL.
 *  Since this first public domain version in 1999, the developments are
 *  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
 *  INRIA.
 *
 *  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
 *  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.
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
 *  references [1], [2], and [3]) the contribution of this package
 *  in any scientific publication dependent upon the use of the
 *  package. You shall use reasonable endeavours to notify
 *  the authors of the package of this publication.
 *
 *   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent,
 *   Multifrontal parallel distributed symmetric and unsymmetric solvers,
 *   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).
 *
 *   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
 *   A fully asynchronous multifrontal solver using distributed dynamic
 *   scheduling, SIAM Journal of Matrix Analysis and Applications,
 *   Vol 23, No 1, pp 15-41 (2001).
 *
 *   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
 *   S. Pralet, Hybrid scheduling for the parallel solution of linear
 *   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
 *
 */
/*    $Id: mumps_io_err.h 5043 2008-07-18 08:56:02Z pcombes $  */
#include <errno.h>
#include "mumps_common.h"
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
# include <pthread.h>
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
/* Exported global variables */
extern char error_str[200];
extern char* mumps_err;
extern int* dim_mumps_err;
extern int err_flag;
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
extern pthread_mutex_t err_mutex;
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
/* Exported functions */
#define MUMPS_LOW_LEVEL_INIT_ERR_STR \
    F_SYMBOL(low_level_init_err_str,LOW_LEVEL_INIT_ERR_STR)
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_ERR_STR( int* dim, char* err_str, mumps_ftnlen l1);
int mumps_io_prop_err_info(int ierr);
int mumps_io_build_err_str(int errnum, int mumps_err,const char* desc,char* buf,int size);
#if ! ( defined(MUMPS_WIN32) || defined(WITHOUT_PTHREAD) )
int mumps_io_init_err_lock();
int mumps_io_destroy_err_lock();
int mumps_check_error_th();
#endif /* ! ( MUMPS_WIN32 || WITHOUT_PTHREAD ) */
