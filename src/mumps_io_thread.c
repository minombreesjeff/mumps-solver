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
/*    $Id: mumps_io_thread.c 5043 2008-07-18 08:56:02Z pcombes $  */
#include "mumps_io_basic.h"
#include "mumps_io_err.h"
#include "mumps_io_thread.h"
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
/*******************************************************************************/
/* Algorithm:                                                                  */
/*                                                                             */
/* The communication between io-thread and mumps thread is done with 2 queues: */
/* 1- A queue for posting io-requests (by mumps)                               */
/* 2- A queue for storing the request-id's of finished requests that have not  */
/* been yet checked by mumps.                                                  */
/*                                                                             */
/* mumps-thread:                                                               */
/* =============                                                               */
/*                                                                             */
/* posts requests for io's in the queue (io_queue) when reading and writing    */
/* data is needed by mumps.                                                    */
/* When the queue is full, it waits until there is an available position in    */
/* the queue. The wait, is a succession of sleep periods together with tests   */
/* on the number of free positions in the queue. Note that in this             */
/* case, we try to treat finished requests to avoid deadlocks with the         */
/* io-thread (we want to avoid to have the two queues simultaneously full.     */
/* We provide at the mumps-thread level functions to deal with requests        */
/* and to clean the finished requests queue (clean_finished_queue_th,          */
/* mumps_test_request_th, mumps_wait_request_th, mumps_wait_all_requests_th */
/*                                                                             */
/* io-thread:                                                                  */
/* ==========                                                                  */
/*                                                                             */
/* This thread is created by the mumps-thread when mumps_low_level_init_ooc_c  */
/* is called. It's execution scheme is to loop until the termination of the    */
/* execution (io_flag_stop=1). The loop is described below :                   */
/*                                                                             */
/* While (not termination){                                                    */
/*   if(there is no pending request){                                          */
/*      sleep;                                                                 */
/*	continue with next iteration;                                          */
/*   }else{                                                                    */
/*      treat the io request (do the read or the write);                       */
/*      wait for the moment where there is a free position in the finished     */
/*      request queue (just in case it is full);                               */
/*      put the request in the queue of finished request                       */
/*   }                                                                         */
/* }                                                                           */
/*                                                                             */
/*                                                                             */
/* Two mechanisms for synchronization have been designed: one based on mutexes */
/* and sleeps, and the other one on mutexes and semaphores.                    */
/*                                                                             */
/* Mechanism 1:                                                                */
/* ============                                                                */
/*                                                                             */
/* One mutex (io_mutex) is used to ensure that things are accessed in a good   */
/* way.                                                                        */
/* In the case where one of the queues if full, the mechanism used is to sleep */
/* for a period and to check again until a free position is available.         */
/*                                                                             */
/* Mechanism 2:                                                                */
/* ============                                                                */
/*                                                                             */
/* One mutex (io_mutex) is used to ensure that things are accessed in a good   */
/* way. In addition semaphores are used to synchronise threads when accessing  */
/* the queues. The idea is to initialize a semaphore to the length of the      */
/* queue. Then each time a request is inserted, the value of the semaphore     */
/* is decreased by the thread that does the operation.                         */
/* Before accessing to a queue, the thread executes sem_wait on the semaphore. */
/* Thus, if the value of the semaphore is 0 (there no free positions in the    */
/* associated queue), the thread is suspend until the moment when  the value   */
/* of the semaphore is greater than 0. The advantage of the mechanism is to    */
/* avoid the sleep periods. The semaphores are modified as follow :            */
/*                                                                             */
/* mumps_thread :                                                              */
/* ==============                                                              */
/* 1- For each io_operation call a sem_wait on the sem_nb_free_active_requests to check         */
/* that there is a free position in the queue. In addition do a sem_post       */
/* on sem_io (to tell the other thread that there is a new request).           */
/* 2- When cleaning a finished_io_request (mumps_clean_request) call a         */
/* sem_post on the sem_nb_finish to increase the value of the semaphore        */
/* corresponding to the number of free position in the finished requests       */
/* queue.                                                                      */
/* 3- At the end of the execution call a sem_post on sem_top.                  */
/*                                                                             */
/* io_thread :                                                                 */
/* ===========                                                                 */
/* 1- Check the value of sem_stop using sem_getvalue at each iteration of the  */
/* main loop.                                                                  */
/* 2- At the begining of the loop call a sem_wait on sem_io to know if there   */
/* is something to do.                                                         */
/* 3- for each treated request (from the io_request queue) do a sem_post on    */
/* the sem_nb_free_active_requests semaphore.                                                   */
/* 4- before inserting the treated request in the finished_request             */
/* queue do a sem_wait on sem_nb_finish to check if there available            */
/* positions.                                                                  */
/*                                                                             */
/* Important remark:                                                           */
/* =================                                                           */
/* Before the doing a sem_wait, the mumps_thread always calls                  */
/* mumps_clean_finished_queue_th (to clean the finished_queue). This is done to*/
/* avoid deadlocks (the case where the two queues are full).                   */
/*                                                                             */
/* Iinitial values:                                                            */
/* ================                                                            */
/* sem_io -> 0 (comptabilise le nombre de taches restantes a commencer) ;      */
/* sem_nb_free_active_requests -> MAX_IO (empeche le debordement du nombre de requetes          */
/*                       d'IO empilees) ;                                      */
/* sem_nb_free_finished_requests -> MAX_FINISH_REQ (                                           */
/* sem_stop -> 1 (gere l'arret du thread d'I/O) lorsque MUMPS termine.         */
/*                                                                             */
/*                                                                             */
/*******************************************************************************/
/* Exported global variables */
int io_flag_stop,current_req_num;
pthread_t io_thread,main_thread;
pthread_mutex_t io_mutex;
pthread_cond_t cond_io,cond_nb_free_finished_requests,cond_nb_free_active_requests,cond_stop;
pthread_mutex_t io_mutex_cond;
int int_sem_io,int_sem_nb_free_finished_requests,int_sem_nb_free_active_requests,int_sem_stop;
int with_sem;
struct request_io *io_queue;
int first_active,last_active,nb_active;
int *finished_requests_inode,*finished_requests_id,first_finished_requests,
  last_finished_requests,nb_finished_requests,smallest_request_id;
int mumps_owns_mutex;
int test_request_called_from_mumps;
/* Other global variables */
double inactive_time_io_thread;
int time_flag_io_thread;
struct timeval origin_time_io_thread;
/**
 * Main loop of the io thread when semaphores are not used
 */
void* mumps_async_thread_function (void* arg){   
  struct request_io *current_io_request;
  struct timeval start_time,end_time;
  int ierr,ret_code;
   /* Locks the IO mutex */
   pthread_mutex_lock(&io_mutex);
   /* Loops */
   while (!io_flag_stop){
      if(nb_active==0){
	if(time_flag_io_thread){
	  gettimeofday(&start_time,NULL);
	}
	pthread_mutex_unlock(&io_mutex);
	usleep(50);
	gettimeofday(&end_time,NULL);
	if(time_flag_io_thread){
	  inactive_time_io_thread=inactive_time_io_thread+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
	}else{
	  inactive_time_io_thread=((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)origin_time_io_thread.tv_sec+((double)origin_time_io_thread.tv_usec/1000000));
	}
      }
      else{
	if(!time_flag_io_thread) 
	  time_flag_io_thread=1;
	pthread_mutex_unlock(&io_mutex);
	current_io_request=&io_queue[first_active];
	if(current_io_request->io_type==0){
	  /*write*/	
	  ret_code=mumps_io_do_write_block(current_io_request->addr,
				    &(current_io_request->size),
                                    &(current_io_request->file_type),
				    current_io_request->vaddr,
				    &ierr);
	  if(ret_code<0){
	    break;
	  }
	}
	else{
	  if(current_io_request->io_type==1){
	    /*read*/
	    ret_code=mumps_io_do_read_block(current_io_request->addr,
				     &(current_io_request->size),
                                     &(current_io_request->file_type),
    				     current_io_request->vaddr,
				     &ierr);
	    if(ret_code<0){
	      break;
	    }
	  }
	  else{
	    printf("Internal error in mumps_async_thread_function\n");
	    exit (-3);
	  }
	}
	pthread_mutex_lock(&io_mutex);
	while(nb_finished_requests==MAX_FINISH_REQ-1){
	  pthread_mutex_unlock(&io_mutex);
	  usleep(50);
	  pthread_mutex_lock(&io_mutex);
	}	
	nb_active--;
	if(first_active<MAX_IO-1){
	  first_active++;
	}
	else{
	  first_active=0;
	}	 
	finished_requests_id[last_finished_requests]=current_io_request->req_num;
	finished_requests_inode[last_finished_requests]=current_io_request->inode;
	last_finished_requests=(last_finished_requests+1)%(MAX_FINISH_REQ);
	nb_finished_requests++;
	pthread_mutex_unlock(&io_mutex);  
      }      
      pthread_mutex_lock(&io_mutex);  
   }
   pthread_exit(NULL);
/* Not reached */
   return NULL;
}
/**
 * Main function of the io thread when semaphores are used.
 */
void*  mumps_async_thread_function_with_sem (void* arg){
   struct request_io *current_io_request;
   int ierr,_sem_stop;
   struct timeval start_time,end_time;
   int ret_code;
   for (;;){      
     gettimeofday(&start_time,NULL);
       if(with_sem==2){
	 mumps_wait_sem(&int_sem_io,&cond_io);
       }
     /*     sem_wait(&sem_io);  */
     gettimeofday(&end_time,NULL);
     if(time_flag_io_thread){
       inactive_time_io_thread=inactive_time_io_thread+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
     }else{
       inactive_time_io_thread=((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)origin_time_io_thread.tv_sec+((double)origin_time_io_thread.tv_usec/1000000));
     }
     if(!time_flag_io_thread){
       time_flag_io_thread=1;
     }
     /* Check if the main thread ordered to stop this slave thread */
     /*     sem_getvalue(&sem_stop,&_sem_stop); */
       if(with_sem==2){
	 mumps_get_sem(&int_sem_stop,&_sem_stop);
       }
     if(_sem_stop==IO_FLAG_STOP){
       /* The thread must stop */
       break; /* Breaks the while loop. */
     }
      current_io_request=&io_queue[first_active];
      switch(current_io_request->io_type)
	{
	 case IO_WRITE:
	   ret_code=mumps_io_do_write_block(current_io_request->addr,
				     &(current_io_request->size),
                                     &(current_io_request->file_type),
				     current_io_request->vaddr,
				     &ierr);
	   if(ret_code<0){
	     goto end;
	   }
	   break;
	case IO_READ:
	  ret_code=mumps_io_do_read_block(current_io_request->addr,
				   &(current_io_request->size),
                                   &(current_io_request->file_type),
              			   current_io_request->vaddr,
				   &ierr);
	   if(ret_code<0){
	     goto end;
	   }	  
	   break;
	 default:
	   printf("Error : Mumps_IO : Operation %d is neither READ nor WRITE\n",current_io_request->io_type);
	   exit (-3);
	}
      /* Notify that the IO was performed */
      /* Wait that finished_requests queue could register 
       the notification */
       if(with_sem==2){
	 mumps_wait_sem(&int_sem_nb_free_finished_requests,&cond_nb_free_finished_requests);
       }
      pthread_mutex_lock(&io_mutex);
      /* Updates active queue bounds */
      /* Register the notification in finished_requests queue
       and updates its bounds. */
      finished_requests_id[last_finished_requests]=current_io_request->req_num;
      finished_requests_inode[last_finished_requests]=current_io_request->inode;
      last_finished_requests=(last_finished_requests+1)%(MAX_FINISH_REQ); /* ??? */
      nb_finished_requests++;      
      /* Realeases the lock : ***UNLOCK*** */
      nb_active--;      
      if(first_active<MAX_IO-1){
	 first_active++; 
      }
      else{
	 first_active=0;
      }
      if(with_sem==2){
	mumps_post_sem(&(current_io_request->int_local_cond),&(current_io_request->local_cond));
      }
      pthread_mutex_unlock(&io_mutex);  
      /* Finally increases the number of free active requests.*/
      /*      sem_post(&sem_nb_free_active_requests); */
	  mumps_post_sem(&int_sem_nb_free_active_requests,&cond_nb_free_active_requests);
   }
 end:
   /* The main thread ordered the end of the IO thread (it changed sem_stop).
    We exit. */
   pthread_exit(NULL);
/* Not reached */
   return NULL;
}
int mumps_is_there_finished_request_th(int* flag){
  if(!mumps_owns_mutex) pthread_mutex_lock(&io_mutex);
  if(nb_finished_requests==0){
    *flag=0;
  }else{
    *flag=1;
    /*    printf("finished : %d\n",nb_finished_requests);     */
  }
  if(!mumps_owns_mutex) pthread_mutex_unlock(&io_mutex);
  return 0;
}
int mumps_clean_request_th(int* request_id){
  int ierr;
  ierr=mumps_check_error_th();
  if(ierr!=0){
    return ierr;
  }
  if(!mumps_owns_mutex)pthread_mutex_lock(&io_mutex);
  *request_id=finished_requests_id[first_finished_requests];
  if(smallest_request_id!=finished_requests_id[first_finished_requests]){
    sprintf(error_str,"Internal error in OOC Management layer (mumps_clean_request_th)\n");
    return -91;
  }
  finished_requests_id[first_finished_requests]=-9999;
  first_finished_requests=(first_finished_requests+1)%(MAX_FINISH_REQ);
  nb_finished_requests--;
  /*we treat the io requests in their arrival order => we just have to
    increase smallest_request_id*/
  smallest_request_id++;
  if(!mumps_owns_mutex) pthread_mutex_unlock(&io_mutex);
  if(with_sem) {
      if(with_sem==2){
	mumps_post_sem(&int_sem_nb_free_finished_requests,&cond_nb_free_finished_requests);
      }
  }
  return 0;
}
int mumps_clean_finished_queue_th(){
   /* Cleans the finished request queue. On exit, the queue is empty.*/
   int local_flag;
   int cur_req;
   int loc_owned_mutex=0,ierr;
   if(!mumps_owns_mutex){
      pthread_mutex_lock(&io_mutex);
      mumps_owns_mutex=1;
      loc_owned_mutex=1;
  }
  /* this block of code is designed for avoiding deadlocks between
     the two threads*/
   mumps_is_there_finished_request_th(&local_flag);
   while(local_flag){
     ierr=mumps_clean_request_th(&cur_req);
     if(ierr!=0){
       return ierr;
     }
     mumps_is_there_finished_request_th(&local_flag);
   }
   if((!mumps_owns_mutex)||(loc_owned_mutex)){
      pthread_mutex_unlock(&io_mutex);
      mumps_owns_mutex=0;
   }
   return 0;
}
int mumps_test_request_th(int* request_id,int *flag){
  /* Tests if the request "request_id" has finished. It sets the flag  */
  /* argument to 1 if the request has finished (0 otherwise)           */
  int request_pos;
  int i;
  i=mumps_check_error_th();
  if(i!=0){
    return i;
  }
  pthread_mutex_lock(&io_mutex);
  /*  printf("entering test !!! \n"); */
  if(*request_id < smallest_request_id){
    *flag=1;
    /*    exit (-2); */
  }else{
    if(nb_finished_requests==0){
      *flag=0;
    }else{
      request_pos=(first_finished_requests+nb_finished_requests-1)%(MAX_IO*2);    
      if(*request_id > finished_requests_id[request_pos]){
	/*the request has not been treated yet since it is not in 
	the list of treated requests*/
	i=0;
	/*this loop is only for checking (no special treatment is done*/
	while(i<nb_active){
	  request_pos=(first_active+i)%(MAX_IO);
	  if(io_queue[request_pos].req_num==*request_id){
	    break;
	  }
	  i++;
	}
	if(i==nb_active){
	  sprintf(error_str,"Internal error in OOC Management layer (mumps_test_request_th (1))\n");
	  return -91;
	}
	*flag=0;
      }else{
	/*the request has been treated yet since it is the list of
	  treated requests*/
	i=0;
	while(i<nb_finished_requests){
	  request_pos=(first_finished_requests+i)%(MAX_IO*2);
	  if(finished_requests_id[request_pos]==*request_id){
	    break;
	  }
	  i++;
	}
	if(i==nb_finished_requests){
	  sprintf(error_str,"Internal error in OOC Management layer (mumps_test_request_th (2))\n");
	  return -91;
	}      
	*flag=1;
      }
    }
  }
  /*  printf("je vais essayer de nettoyer %d\n",nb_finished_requests);
    pthread_mutex_unlock(&io_mutex);
    sleep (10); */
  mumps_owns_mutex=1;
  /* printf("ici nb_finished_requests=%d\n",nb_finished_requests);*/
  mumps_clean_finished_queue_th();
  mumps_owns_mutex=0;
  pthread_mutex_unlock(&io_mutex);
  return 0;
}
int mumps_wait_req_sem_th(int *request_id){
  int i,j;
  j=first_active;
  for(i=0;i<nb_active;i++){
    if(io_queue[j].req_num==*request_id) break;
    j=(j+1)%MAX_IO;
  }
  if(i<nb_active){
    mumps_wait_sem(&(io_queue[j].int_local_cond),&(io_queue[j].local_cond));
  }
  return 0;
}
int mumps_wait_request_th(int *request_id){
  /* waits for the termination of the request "request_id" */
  int flag=0,ierr;
  if(with_sem!=2){
    while(!flag){
      ierr=mumps_test_request_th(request_id,&flag);
      if(ierr!=0)return ierr;
    }
  }else{
    ierr=mumps_test_request_th(request_id,&flag);
    if(ierr!=0)return ierr;
    if(!flag){
      mumps_wait_req_sem_th(request_id);
      ierr=mumps_test_request_th(request_id,&flag);
      if(ierr!=0)return ierr;
    }
  }
  return 0;
}
int mumps_wait_all_requests_th(){
  /* waits for the termination of all the oingoing requests */
  int flag,ierr;
  pthread_mutex_lock(&io_mutex);
  while((nb_active!=0)||(nb_finished_requests!=0)){
    pthread_mutex_unlock(&io_mutex);
    mumps_is_there_finished_request_th(&flag);
    if(flag){
      ierr=mumps_clean_request_th(&flag);
      if(ierr!=0){
	return ierr;
      }
    }else{
      usleep(10);
    }
    pthread_mutex_lock(&io_mutex);
  }
  pthread_mutex_unlock(&io_mutex); 
  return 0;
}
int mumps_low_level_init_ooc_c_th(int* async, int* ierr){
  int ret_code;    
  int i;
  /* Computes the number of files needed. Uses ceil value. */
  current_req_num=0;
  with_sem=2;
  first_active=0;
  last_active=0;
  nb_active=0;
  first_finished_requests=0;
  last_finished_requests=0;
  nb_finished_requests=0;
  smallest_request_id=0;
  mumps_owns_mutex=0;
  inactive_time_io_thread=0;
  time_flag_io_thread=0;
  gettimeofday(&origin_time_io_thread,NULL);
  /*  mumps_io_flag_async=*async; */
  if(*async!=IO_ASYNC_TH){
    sprintf(error_str,"Internal error: mumps_low_level_init_ooc_c_th should not to be called with strat_IO=%d\n",*async);
    return -91;
  }
  if(*async){
    pthread_mutex_init(&io_mutex,NULL);
    mumps_io_init_err_lock();
#ifdef WITH_PFUNC
    mumps_io_init_pointers_lock();
#endif
    io_queue=(struct request_io *)malloc(MAX_IO*sizeof(struct request_io));
    if(with_sem==2){
      for(i=0;i<MAX_IO;i++){
	pthread_cond_init(&(io_queue[i].local_cond),NULL);
	io_queue[i].int_local_cond=0;
      }
    }
    finished_requests_id=(int *)malloc(MAX_IO*2*sizeof(int));
    finished_requests_inode=(int *)malloc(MAX_IO*2*sizeof(int));
    for(i=0;i<MAX_IO*2;i++){
      finished_requests_id[i]=-9999;
      finished_requests_inode[i]=-9999;
    }
    if(with_sem){
      switch(with_sem){
      case 2:
	int_sem_io=0;
	int_sem_stop=0;
	int_sem_nb_free_finished_requests=MAX_FINISH_REQ;
	int_sem_nb_free_active_requests=MAX_IO;
	pthread_cond_init(&cond_stop,NULL);
	pthread_cond_init(&cond_io,NULL);
	pthread_cond_init(&cond_nb_free_active_requests,NULL);
	pthread_cond_init(&cond_nb_free_finished_requests,NULL);
	pthread_mutex_init(&io_mutex_cond,NULL);
	break;
      default:
	sprintf(error_str,"Internal error: mumps_low_level_init_ooc_c_th should not to be called with strat_IO=%d\n",*async);
	return -92;
      }
      ret_code=pthread_create(&io_thread,NULL,mumps_async_thread_function_with_sem,NULL);
    }else{
      ret_code=pthread_create(&io_thread,NULL,mumps_async_thread_function,NULL);
    }
    if(ret_code!=0){
      mumps_io_build_err_str(ret_code,-92,"Unable to create I/O thread",error_str,200);
      return -92;      
    }
    main_thread=pthread_self();
  }
  return 0;
}
int mumps_async_write_th(const int * strat_IO, 
			void * address_block,
			int * block_size,
			int * inode,
			int * request_arg,
		        int * type,
			long long vaddr,
			int * ierr){
  int cur_req;
  *ierr=mumps_check_error_th();
  if(*ierr!=0){
    return *ierr;
  }
  if(with_sem){ 
    mumps_clean_finished_queue_th();
      if(with_sem==2){
	mumps_wait_sem(&int_sem_nb_free_active_requests,&cond_nb_free_active_requests);
      }
    /*    sem_wait(&sem_nb_free_active_requests); */
    pthread_mutex_lock(&io_mutex);
  }else{
    pthread_mutex_lock(&io_mutex);
    mumps_owns_mutex=1;
    mumps_clean_finished_queue_th();
    mumps_owns_mutex=0;
    while(nb_active==MAX_IO){
      pthread_mutex_unlock(&io_mutex);
      usleep(10);
      pthread_mutex_lock(&io_mutex);
    }
  }
  if(nb_active<=MAX_IO){
    if(nb_active==0){
      first_active=last_active;
    }
    else{
      last_active=(last_active+1)%MAX_IO;
    }
    cur_req=last_active;
    nb_active++;
    io_queue[cur_req].inode=*inode; 
    io_queue[cur_req].req_num=current_req_num; 
    io_queue[cur_req].addr=address_block; 
    io_queue[cur_req].size=*block_size;  
    io_queue[cur_req].vaddr=vaddr;
    io_queue[cur_req].io_type=0;
    io_queue[cur_req].file_type=*type;
    if(with_sem==2){
      io_queue[cur_req].int_local_cond=0;      
    }
    *request_arg=current_req_num; 
    current_req_num++;
  }else{
    sprintf(error_str,"Internal error in OOC Management layer (mumps_async_write_th)\n");
    return -91;
    /*    exit(-3);*/
  }
  pthread_mutex_unlock(&io_mutex);
  if(with_sem){
    /*    sem_post(&sem_io); */
      if(with_sem==2){
	mumps_post_sem(&int_sem_io,&cond_io);
      }
  }
  return 0;
}
int mumps_async_read_th(const int * strat_IO, 
		       void * address_block,
		       int * block_size,
		       int * inode,
		       int * request_arg,
         	       int * type,
		       long long vaddr,
		       int * ierr){
  int cur_req;  
  *ierr=mumps_check_error_th();
  if(*ierr!=0){
    return *ierr;
  }
  if(with_sem){
    mumps_clean_finished_queue_th();
    /* end of the block*/
      if(with_sem==2){
	mumps_wait_sem(&int_sem_nb_free_active_requests,&cond_nb_free_active_requests);
      }
    /*    sem_wait(&sem_nb_free_active_requests); */
    pthread_mutex_lock(&io_mutex);
  }else{
    pthread_mutex_lock(&io_mutex);
    /* this block of code is designed for avoiding deadlocks between
       the two threads*/
    mumps_owns_mutex=1;
    mumps_clean_finished_queue_th();
    mumps_owns_mutex=0;
    while(nb_active==MAX_IO){
      pthread_mutex_unlock(&io_mutex);
      usleep(10);
      pthread_mutex_lock(&io_mutex);
    }
  }
  if(nb_active<MAX_IO){
    if(nb_active==0){
      first_active=last_active;
    }else{
      last_active=(last_active+1)%MAX_IO;
      /*	if(last_active<MAX_IO-1){
		cur_req=last_active+1;
		last_active++;
		}else{
		cur_req=0;
		last_active=0;
		}*/
    }
    cur_req=last_active;
    nb_active++;
    io_queue[cur_req].inode=*inode;
    io_queue[cur_req].req_num=current_req_num; 
    io_queue[cur_req].addr=address_block; 
    io_queue[cur_req].size=*block_size;  
    io_queue[cur_req].vaddr=vaddr;
    io_queue[cur_req].io_type=1;
    io_queue[cur_req].file_type=*type;
    if(with_sem==2){
      io_queue[cur_req].int_local_cond=0;      
    }
    *request_arg=current_req_num; 
    current_req_num++;
  }else{
    sprintf(error_str,"Internal error in OOC Management layer (mumps_async_read_th)\n");
    return -91;
  }
  if(with_sem){
    /*    sem_post(&sem_io); */
      if(with_sem==2){
	mumps_post_sem(&int_sem_io,&cond_io);
      }
  }
  pthread_mutex_unlock(&io_mutex);
  return 0;
}
int mumps_clean_io_data_c_th(int *myid){
  int i;
  /* cleans the thread/io management data*/
  if(mumps_io_flag_async){
    /*we can use either signals or mutexes for this step */
    if(with_sem){
	if(with_sem==2){
	  mumps_post_sem(&int_sem_stop,&cond_stop);
	  mumps_post_sem(&int_sem_io,&cond_io);
	}
    }else{
      pthread_mutex_lock(&io_mutex);
      io_flag_stop=1;
      pthread_mutex_unlock(&io_mutex);
    }
    pthread_join(io_thread,NULL);
    pthread_mutex_destroy(&io_mutex);
    mumps_io_destroy_err_lock();
#ifdef WITH_PFUNC
    mumps_io_destroy_pointers_lock();
#endif
    if(with_sem){
	if(with_sem==2){
	  pthread_cond_destroy(&cond_stop);
	  pthread_cond_destroy(&cond_io);
	  pthread_cond_destroy(&cond_nb_free_active_requests);
	  pthread_cond_destroy(&cond_nb_free_finished_requests);
	  pthread_mutex_destroy(&io_mutex_cond);
	}
    }
  }
  if(with_sem==2){
    for(i=0;i<MAX_IO;i++){
      pthread_cond_destroy(&(io_queue[i].local_cond));
    }
  }
  free(io_queue);
  free(finished_requests_id);
  free(finished_requests_inode);
  return 0;
}
int mumps_get_sem(void *arg,int *value){
  switch(with_sem){
  case 2:
    pthread_mutex_lock(&io_mutex_cond);
    *value=*((int *)arg);
    pthread_mutex_unlock(&io_mutex_cond);
    break;
  default:
    sprintf(error_str,"Internal error in OOC Management layer (mumps__get_sem)\n");
    return -91;
  }
  return 0;
}
int mumps_wait_sem(void *arg,pthread_cond_t *cond){
  int *tmp_pointer;
  switch(with_sem){
  case 2:
    pthread_mutex_lock(&io_mutex_cond);
    tmp_pointer=(int *)arg;
    while(*tmp_pointer==0){
      pthread_cond_wait(cond, &io_mutex_cond);
    }
    (*tmp_pointer)--;
    pthread_mutex_unlock(&io_mutex_cond);
    break;
  default:
    sprintf(error_str,"Internal error in OOC Management layer (mumps_wait_sem)\n");
    return -91;
  }
  return 0;  
}
int mumps_post_sem(void *arg,pthread_cond_t *cond){
  int *tmp_pointer;
  switch(with_sem){
  case 2:
    pthread_mutex_lock(&io_mutex_cond);
    tmp_pointer=(int *)arg;
    (*tmp_pointer)++;
    if(*tmp_pointer==1){
      pthread_cond_broadcast(cond);
    }
    pthread_mutex_unlock(&io_mutex_cond);
    break;
  default:
    sprintf(error_str,"Internal error in OOC Management layer (mumps_post_sem)\n");
    return -91;
  }
  return 0;  
}
#endif /* _WIN32 && WITHOUT_PTHREAD */
