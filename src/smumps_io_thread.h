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

/*    $Id: smumps_io_thread.h,v 1.21 2006/03/27 16:46:59 jylexcel Exp $  */
#ifndef _WIN32  

#include <unistd.h>
#include <pthread.h>
#include <sys/types.h>


#include <sys/time.h>
#include <time.h>

#define MAX_IO 20
#define MAX_FINISH_REQ 40

#define IO_FLAG_STOP 1
#define IO_FLAG_RUN 0

#define IO_READ 1
#define IO_WRITE 0

struct request_io{
  int inode;
  int req_num; /*request number*/
  void* addr;  /*memory address (either source or dest)*/
  int size;    /* size of the requested io (unit=size of elementary mumps data)*/
  int* pos_in_file; /*offset in the swap file*/
  int* num_file; /* file number*/
  int io_type; /*read or write*/
  pthread_cond_t local_cond;
  int int_local_cond;
};

void*  _smumps_async_thread_function_with_sem (void* arg);
void*  _smumps_async_thread_function (void* arg);

int   smumps_wait_req_sem_th(int *request_id);
int   smumps_is_there_finished_request_th(int* flag);
int   smumps_clean_request_th(int* request_id);
int   smumps_test_request_th(int* request_id,int *flag);
int   smumps_wait_request_th(int *request_id);
int   smumps_wait_all_requests_th();
int   smumps_low_level_init_ooc_c_th(int* async, int* ierr);
int   smumps_async_write_th(const int * strat_IO,void * address_block,int * block_size,int * pos_in_file,
                           int * file_number,int * inode,int * request_arg,int * ierr);
int   smumps_async_read_th(const int * strat_IO,void * address_block,int * block_size,int * from_where,
                          int * file_number,int * inode,int * request_arg,int * ierr);
int smumps_clean_io_data_c_th(int *myid);

int _smumps_get_sem(void *arg,int *value);

int _smumps_wait_sem(void *arg,pthread_cond_t *cond);

int _smumps_post_sem(void *arg,pthread_cond_t *cond);

int _smumps_clean_finished_queue_th();

#endif /*_WIN32*/

