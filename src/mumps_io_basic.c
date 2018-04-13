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
/*    $Id: mumps_io_basic.c 5083 2008-07-24 14:48:27Z aguermou $  */
#include "mumps_io_basic.h"
#include "mumps_io_err.h"
/* Exported global variables */
#if ! defined (MUMPS_WIN32)
# if defined(WITH_PFUNC) && ! defined (WITHOUT_PTHREAD)
#  include <pthread.h>
pthread_mutex_t mumps_io_pwrite_mutex;
# endif
/* int* mumps_io_pfile_pointer_array; */
/* int* mumps_io_current_file; */
/* #else /\*_WIN32*\/ */
/* FILE** mumps_io_current_file; */
/* FILE** mumps_io_pfile_pointer_array; */
#endif /* MUMPS_WIN32 */
/* mumps_file_struct* mumps_io_pfile_pointer_array;
   mumps_file_struct* mumps_io_current_file; */
mumps_file_type* mumps_files;
/* int mumps_io_current_file_number; */
char* mumps_ooc_file_prefix;
/* char** mumps_io_pfile_name; */
/* int mumps_io_current_file_position; */
/* int mumps_io_write_pos; */
/* int mumps_io_last_file_opened; */
int mumps_elementary_data_size;
int mumps_io_is_init_called;
int mumps_io_myid;
int mumps_io_max_file_size;
/* int mumps_io_nb_file; */
int mumps_io_flag_async;
int mumps_io_k211;
/* int mumps_flag_open;*/
int mumps_directio_flag;
int mumps_io_nb_file_type;
/* Functions */
int mumps_set_file(int type,int file_number_arg){
  /* Defines the pattern for the file name. The last 6 'X' will be replaced
     so as to name were unique */
  char name[150];
  int ret_code, fd;  
  mumps_file_struct  *mumps_io_pfile_pointer_array;
  /*  if ((mumps_files+type)->mumps_io_current_file_number >= ((mumps_files+type)->mumps_io_nb_file)-1){*/
  if (file_number_arg > ((mumps_files+type)->mumps_io_nb_file)-1){
    /* Exception : probably thrown because of a bad estimation
       of number of files. */
    /* We increase the number of file needed and then realloc. */
    ((mumps_files+type)->mumps_io_nb_file)++;
/* #if ! defined(_WIN32) || defined(__MINGW32)  */
/*     mumps_io_pfile_pointer_array=realloc((void *)mumps_io_pfile_pointer_array,mumps_io_nb_file*sizeof(int)); */
/* #else */
/*     mumps_io_pfile_pointer_array=realloc((void *)mumps_io_pfile_pointer_array,mumps_io_nb_file*sizeof(FILE*)); */
/* #endif */
    (mumps_files+type)->mumps_io_pfile_pointer_array=realloc((void *)(mumps_files+type)->mumps_io_pfile_pointer_array,((mumps_files+type)->mumps_io_nb_file)*sizeof(mumps_file_struct));
    if((mumps_files+type)->mumps_io_pfile_pointer_array==NULL){
      sprintf(error_str,"Allocation problem in low-level OOC layer\n");
      return -13;
    }
/*     mumps_io_pfile_name=realloc((void*)mumps_io_pfile_name,(mumps_io_nb_file)*sizeof(char *)); */
/*     if(mumps_io_pfile_name==NULL){ */
/*       sprintf(error_str,"Allocation problem in low-level OOC layer\n"); */
/*       return -13; */
/*     } */
    (((mumps_files+type)->mumps_io_pfile_pointer_array)+((mumps_files+type)->mumps_io_nb_file)-1)->is_opened=0;
  }
  mumps_io_pfile_pointer_array=(mumps_files+type)->mumps_io_pfile_pointer_array;  
  /* Do change the current file */
  ((mumps_files+type)->mumps_io_current_file_number)=file_number_arg;
  if((mumps_io_pfile_pointer_array+file_number_arg)->is_opened!=0){
    ((mumps_files+type)->mumps_io_current_file_number)=file_number_arg;
    /*    printf("je sors %d | %d\n",file_number_arg,((mumps_files+type)->mumps_io_nb_file)-1);*/
    return 0;
  }
/* #if ! defined(_WIN32) || defined(__MINGW32__)*/
/* MinGW does not have a mkstemp function and MinGW defines _WIN32,
 * so we also go in the else branch below with MinGW */
#if ! defined(_WIN32)
  strcpy(name,mumps_ooc_file_prefix);
  fd=mkstemp(name);  
  /* Note that a file name is built by mkstemp and that the file is 
     opened. fd hold the file descriptor to access it.
     We want to close the file that will be opened later
     and might be removed before the end of the processus.
    */
    if(fd<0){
      ret_code=-99;
      return ret_code;
    }
    else{
      close(fd); 
    }
#else
  sprintf(name,"%s_%d",mumps_ooc_file_prefix,((mumps_files+type)->mumps_io_current_file_number)+1);
#endif
/*   *(mumps_io_pfile_pointer_array+mumps_io_current_file_number)=fopen(name,"w+"); */
/*   *(mumps_io_pfile_name+mumps_io_current_file_number)=(char *)malloc((strlen(name)+1)*sizeof(char)); */
/*   if(*(mumps_io_pfile_name+mumps_io_current_file_number)==NULL){ */
/*     sprintf(error_str,"Allocation problem in low-level OOC layer\n"); */
/*     return -13; */
/*   } */
  strcpy((mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->name,name);
  /* See mumps_io_basic.h for comments on the I/O flags passed to open */
#if ! defined(_WIN32) || defined(__MINGW32__)
   (mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->file=open(name,(mumps_files+type)->mumps_flag_open); 
  /* 
CPA: for LU factor file: 
(mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->file= open(name, O_WRONLY | O_CREAT | O_TRUNC, 0666); */
  if((mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->file==-1){
    mumps_io_build_err_str(errno,-90,"Unable to open OOC file",error_str,200);
    return -90;
  }
#else
  (mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->file=fopen(name,(mumps_files+type)->mumps_flag_open);
  if((mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number)->file==NULL){
    /*    mumps_io_build_err_str(errno,-90,"Unable to open OOC file",error_str,200); */
    sprintf(error_str,"Problem while opening OOC file");
    return -90;
   }
#endif
  (mumps_files+type)->mumps_io_current_file=(mumps_io_pfile_pointer_array+(mumps_files+type)->mumps_io_current_file_number);
  ((mumps_files+type)->mumps_io_nb_file_opened)++;
  if((mumps_files+type)->mumps_io_current_file_number>(mumps_files+type)->mumps_io_last_file_opened){
    (mumps_files+type)->mumps_io_last_file_opened=(mumps_files+type)->mumps_io_current_file_number;
  }
  /*  if(*(mumps_io_pfile_pointer_array+mumps_io_current_file_number)==NULL) */
  ((mumps_files+type)->mumps_io_current_file)->write_pos=0;
  ((mumps_files+type)->mumps_io_current_file)->is_opened=1;
  /*  printf("new file created -> num = %d \n", ((mumps_files+type)->mumps_io_last_file_opened));*/
  /*  printf("new file created %d\n",mumps_io_current_file_number);*/
  return 0;
}
void mumps_update_current_file_position(mumps_file_struct* file_arg){
  file_arg->current_pos=file_arg->write_pos;
/*   mumps_io_current_file_position=mumps_io_write_pos; */
}
int mumps_compute_where_to_write(const double to_be_written,const int type,long long vaddr,size_t already_written){
  /* Check if the current file has enough memory to receive the whole block*/
  int ret_code;
  int file;
  mumps_file_struct *current_file;
  long long vaddr_loc;
  int pos;
  /* Virtual address based file management scheme */
  vaddr_loc=vaddr*(long long)mumps_elementary_data_size+(long long)already_written;
  mumps_gen_file_info(vaddr_loc,&pos,&file);
  ret_code=mumps_set_file(type,file);
  if(ret_code<0){
    return ret_code;
  }
  current_file=(mumps_files+type)->mumps_io_current_file;
  current_file->write_pos=pos;
  mumps_update_current_file_position(current_file);
  return 0;
}
int mumps_prepare_pointers_for_write(double to_be_written,int * pos_in_file, int * file_number,const int type,long long vaddr,size_t already_written){
  int ret_code;
  ret_code=mumps_compute_where_to_write(to_be_written,type,vaddr,already_written);
  if(ret_code<0){
    return ret_code;
  }
  *pos_in_file=((mumps_files+type)->mumps_io_current_file)->current_pos;
  /* should be modified to take into account the file arg */
  *file_number=(mumps_files+type)->mumps_io_current_file_number;
  return 0;
}
MUMPS_INLINE int mumps_gen_file_info(long long vaddr, int * pos, int * file){
  *file=(int)(vaddr/(long long)mumps_io_max_file_size);
  *pos=(int)(vaddr%(long long)mumps_io_max_file_size);
  return 0;
}
int mumps_compute_nb_concerned_files(int * block_size, int * nb_concerned_files,long long vaddr){
  int file,pos,available_size;
  long long vaddr_loc;
  vaddr_loc=vaddr*(long long)mumps_elementary_data_size;
  mumps_gen_file_info(vaddr_loc,&pos,&file);
  available_size=mumps_io_max_file_size-pos+1;
  *nb_concerned_files=(int)my_ceil((double)(my_max(0,((*block_size)*(double)(mumps_elementary_data_size))-available_size))/(double)mumps_io_max_file_size)+1;
  return 0;
}
int mumps_io_do_write_block(void * address_block,
              		     int * block_size,
         	             int * type_arg,
		             long long vaddr,
		             int * ierr){   
  /* Type of fwrite : size_t fwrite(const void *ptr, size_t size, 
                                    *size_t nmemb, FILE *stream); */
  size_t write_size;
  int i;
  int nb_concerned_files=0;
  int ret_code,file_number_loc,pos_in_file_loc;
  double to_be_written;
#if ! defined(_WIN32) || defined(__MINGW32__)
  int* file;
#else
  FILE** file;
#endif
  int where;
  void* loc_addr;
  int type;
  size_t already_written=0;
  type=*type_arg;
  loc_addr=address_block;
  mumps_compute_nb_concerned_files(block_size,&nb_concerned_files,vaddr);
  to_be_written=((double)mumps_elementary_data_size)*((double)(*block_size));
  /*  printf("nb_concerned -> %d | %lf \n",nb_concerned_files,to_be_written); */
  for(i=0;i<nb_concerned_files;i++){
#if ( ! defined (_WIN32) || defined(__MINGW32__)) && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
    if(mumps_io_flag_async==IO_ASYNC_TH){
      mumps_io_protect_pointers();
    }
#endif
#endif
    ret_code=mumps_prepare_pointers_for_write(to_be_written,&pos_in_file_loc,&file_number_loc,type,vaddr,already_written);
    if(ret_code<0){
#if ( ! defined (_WIN32) || defined(__MINGW32__)) && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
    if(mumps_io_flag_async==IO_ASYNC_TH){
      mumps_io_unprotect_pointers();
    }
#endif
#endif
      return ret_code;
    }
    if((double)(mumps_io_max_file_size-((mumps_files+type)->mumps_io_current_file)->write_pos)>to_be_written){
      write_size=(size_t)to_be_written;
      already_written=(size_t)to_be_written;
    }else{
      write_size=(size_t)((double)(mumps_io_max_file_size-((mumps_files+type)->mumps_io_current_file)->write_pos));
      already_written=already_written+(size_t)write_size;
    }
#if defined (_WIN32) && ! defined(__MINGW32__)
    write_size=(size_t)(int)((write_size)/mumps_elementary_data_size);
#endif
    file=&(((mumps_files+type)->mumps_io_current_file)->file);
    where=((mumps_files+type)->mumps_io_current_file)->write_pos;
#if ( ! defined (_WIN32) || defined(__MINGW32__) ) && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
    if(mumps_io_flag_async==IO_ASYNC_TH){
      mumps_io_unprotect_pointers();
    }
#endif
#endif
    /* printf("1 write -> size = %d | off = %d | file = %d (%d) \n",(int)write_size,where,*file,((mumps_files+type)->mumps_io_current_file)->write_pos); */
    ret_code=mumps_io_write__(file,loc_addr,write_size,where,type);
    if(ret_code<0){
      return ret_code;
    }
#if ( ! defined (_WIN32) || defined(__MINGW32__) ) && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
    if(mumps_io_flag_async==IO_ASYNC_TH){
      mumps_io_protect_pointers();
    }
#endif
#endif
#if ! defined(_WIN32) || defined(__MINGW32__)
    ((mumps_files+type)->mumps_io_current_file)->write_pos=((mumps_files+type)->mumps_io_current_file)->write_pos+((int)write_size);
    to_be_written=to_be_written-((int)write_size);
    loc_addr=(void*)((size_t)loc_addr+write_size);
/*     mumps_io_write_pos=mumps_io_write_pos+((int)write_size); */
/*     to_be_written=to_be_written-((int)write_size); */
/*     loc_addr=(void*)((size_t)loc_addr+write_size); */
#else
    /* fread and write */
    ((mumps_files+type)->mumps_io_current_file)->write_pos=((mumps_files+type)->mumps_io_current_file)->write_pos+((int)write_size*mumps_elementary_data_size);
    to_be_written=to_be_written-((int)write_size*mumps_elementary_data_size);
    loc_addr=(void*)((size_t)loc_addr+(size_t)((int)write_size*mumps_elementary_data_size));
/*     mumps_io_write_pos=mumps_io_write_pos+((int)write_size*mumps_elementary_data_size); */
/*     to_be_written=to_be_written-((int)write_size*mumps_elementary_data_size); */
/*     loc_addr=(void*)((size_t)loc_addr+(size_t)((int)write_size*mumps_elementary_data_size)); */
#endif
#if ( ! defined (_WIN32) || defined(__MINGW32__) )  && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
    if(mumps_io_flag_async==IO_ASYNC_TH){
      mumps_io_unprotect_pointers();
    }
#endif
#endif
  }
  if(to_be_written!=0){
    sprintf(error_str,"Internal (1) error in low-level I/O operation %lf",to_be_written);
    return -90;
  }
  /* printf("write ok -> %d \n");*/
  return 0;
}
int mumps_io_do_read_block(void * address_block,
	            int * block_size,
                    int * type_arg,
         	    long long vaddr,
                    int * ierr){
  int ret_code;
  size_t size;
#if ! defined(_WIN32) || defined(__MINGW32__)
  int* file;
#else
  FILE** file;
#endif
  double read_size;
  int local_fnum,local_offset;
  void *loc_addr;
  long long vaddr_loc;
  int type;
  type=*type_arg;
  /*  if(((double)(*block_size))*((double)(mumps_elementary_data_size))>(double)mumps_io_max_file_size){
    sprintf(error_str,"Internal error in low-level I/O operation (requested size too big for file system) \n");
    return -90;
    }*/
  if(*block_size==0){
    return 0;
  }
  read_size=(double)mumps_elementary_data_size*(double)(*block_size);
  /*  if((*file_number<0)&&(read_size<(double)mumps_io_max_file_size)){
    sprintf(error_str,"Internal error (1) in low level read op\n");
    return -90;
    }*/
  loc_addr=address_block;
  vaddr_loc=vaddr*(long long)mumps_elementary_data_size;
  while(read_size>0){
    /* Virtual addressing based management stuff */
    local_fnum=(int)(vaddr_loc/(long long)mumps_io_max_file_size);
    local_offset=(int)(vaddr_loc%(long long)mumps_io_max_file_size);
    file=&((((mumps_files+type)->mumps_io_pfile_pointer_array)+local_fnum)->file);
    /* printf("1 read | file -> %d | fnum -> %d | vaddr -> %d \n",*file,local_fnum,(int)vaddr_loc); */
#if ! defined(_WIN32) || defined(__MINGW32__)
    if(read_size+(double)local_offset>(double)mumps_io_max_file_size){
      size=(size_t)mumps_io_max_file_size-(size_t)local_offset;
    }else{
      size=(size_t)read_size;
    }
#else
    if(read_size+(double)local_offset>(double)mumps_io_max_file_size){
      size=((size_t)mumps_io_max_file_size-(size_t)local_offset)/(size_t)mumps_elementary_data_size;
    }else{
      size=(size_t)(read_size/mumps_elementary_data_size);
    }
#endif
    ret_code=mumps_io_read__(file,loc_addr,size,local_offset,type);
    if(ret_code<0){
      return ret_code;
    }
#if defined(_WIN32) && ! defined(__MINGW32__)
    size=size*mumps_elementary_data_size;
#endif
    vaddr_loc=vaddr_loc+(long long)size;
    read_size=read_size-(double)size;
    loc_addr=(void*)((size_t)loc_addr+size);
    local_fnum++;
    local_offset=0;
    if(local_fnum>(mumps_files+type)->mumps_io_nb_file){
      sprintf(error_str,"Internal error (2) in low level read op\n");
      return -90;
    }
  }
  return 0;
}
int mumps_free_file_pointers(int *step){
  int i,j,bound,ierr;
  bound=mumps_io_nb_file_type;
/*   if(*step==0){ */
/*     /\* factorization *\/ */
/*     bound=NB_FILE_TYPE_FACTO; */
/*   }else{ */
/*     /\* solve *\/ */
/*     bound=NB_FILE_TYPE_SOLVE; */
/*   } */
  for(j=0;j<bound;j++){
    for(i=0;i<(mumps_files+j)->mumps_io_nb_file_opened;i++){
#if ! defined(_WIN32) || defined(__MINGW32__)
      ierr=close((((mumps_files+j)->mumps_io_pfile_pointer_array)+i)->file);
      if(ierr==-1){
	mumps_io_build_err_str(errno,-90,"Problem while closing OOC file",error_str,200);
	return -90;
      }
#else
      ierr=fclose((((mumps_files+j)->mumps_io_pfile_pointer_array)+i)->file);
      if(ierr==-1){
	sprintf(error_str,"Problem while opening OOC file\n");
	return -90;
      }    
#endif
      /*     free(*(mumps_io_pfile_name+i)); */
    } 
    free((mumps_files+j)->mumps_io_pfile_pointer_array);
  }
/*   free(mumps_io_pfile_name); */
  free(mumps_files);
/*   Free prefix only for facto  */
  if (*step == 0) free(mumps_ooc_file_prefix);
  return 0;
}
int mumps_init_file_structure(int* _myid, int* total_size_io,int* size_element,int nb_file_type,int * flag_tab){
  /* Computes the number of files needed. Uses ceil value. */
  int ierr;
#if ! defined(_WIN32) || defined(__MINGW32__)
  int k211_loc;
  int mumps_flag_open;
#endif
  int i,nb;
  int mumps_io_nb_file;
  mumps_io_max_file_size=MAX_FILE_SIZE;
  mumps_io_nb_file_type=nb_file_type;
  mumps_io_nb_file=(int)((((double)(*total_size_io)*1000000)*((double)(*size_element)))/(double)mumps_io_max_file_size)+1;
  mumps_directio_flag=0;
#if ! defined(_WIN32) || defined(__MINGW32__)
  mumps_flag_open=0;
#endif
  mumps_io_myid=*_myid;
  mumps_elementary_data_size=*size_element;
  /* Allocates the memory necessary to handle the file pointer array.*/
  mumps_files=(mumps_file_type *)malloc(mumps_io_nb_file_type*sizeof(mumps_file_type));
  if(mumps_files==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  ierr=0;
  for(i=0;i<mumps_io_nb_file_type;i++){
    /* FIXME : We have to compute mumps_io_nb_file in a good way when
       we have a file type for L and another one for U */
    if(flag_tab[i]==0||flag_tab[i]==1){
      nb=mumps_io_nb_file;
    }else{
      nb=1;
    }
    (mumps_files+i)->mumps_io_nb_file=nb;
    switch(flag_tab[i]){
    case 0:
#if ! defined(_WIN32) || defined(__MINGW32__)
      (mumps_files+i)->mumps_flag_open=mumps_flag_open|O_WRONLY;
#else
      strcpy((mumps_files+i)->mumps_flag_open,"wb");
#endif
      break;
    case 1:
#if ! defined(_WIN32) || defined(__MINGW32__)
      (mumps_files+i)->mumps_flag_open=mumps_flag_open|O_RDONLY;
#else
      strcpy((mumps_files+i)->mumps_flag_open,"rb");
#endif
      break;
    case 2:
#if ! defined(_WIN32) || defined(__MINGW32__)
      (mumps_files+i)->mumps_flag_open=mumps_flag_open|O_RDWR;
#else
      strcpy((mumps_files+i)->mumps_flag_open,"rwb");
#endif
      break;
    default:
      printf("unknown value of flag_open");
      ierr=-90;
    }
    if(ierr<0){
      return ierr;
    }
    ierr=mumps_io_alloc_file_struct(&nb,i);
    if(ierr<0){
      return ierr;
    }
    ierr=mumps_set_file(i,0);
    if(ierr<0){
      return ierr;
    }
  }
  /* Init the current file.*/
  return 0;
}
int mumps_init_file_name(char* mumps_dir,char* mumps_file,
			 int* mumps_dim_dir,int* mumps_dim_file,int* _myid){
  int i;
  char *tmp_dir,*tmp_fname;
  char base_name[20];
  int dir_flag=0,file_flag=0;
  char mumps_base[10]="mumps_";
  tmp_dir=(char *)malloc(((*mumps_dim_dir)+1)*sizeof(char));
  if(tmp_dir==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  tmp_fname=(char *)malloc(((*mumps_dim_file)+1)*sizeof(char));
  if(tmp_fname==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  for(i=0;i<*mumps_dim_dir;i++){
    tmp_dir[i]=mumps_dir[i];
  }
  tmp_dir[i]=0;
  for(i=0;i<*mumps_dim_file;i++){
    tmp_fname[i]=mumps_file[i];
  }
  tmp_fname[i]=0;  
  if(strcmp(tmp_dir,UNITIALIZED)==0){
    dir_flag=1;
    free(tmp_dir);
    tmp_dir=getenv("MUMPS_OOC_TMPDIR");
    if(tmp_dir==NULL){
#ifdef SP_
#ifndef CINES_
      tmp_dir=getenv("TMPDIR");
      if(tmp_dir==NULL){
	tmp_dir=MUMPS_OOC_DEFAULT_DIR;
      }
#else
      tmp_dir=MUMPS_OOC_DEFAULT_DIR;       
#endif       
#else
      tmp_dir=MUMPS_OOC_DEFAULT_DIR;
#endif      
    }
  }
  if(strcmp(tmp_fname,UNITIALIZED)==0){
    free(tmp_fname);
    tmp_fname=getenv("MUMPS_OOC_PREFIX");
    file_flag=1;
  }
  if(tmp_fname!=NULL){
#if ! defined(_WIN32) || defined(__MINGW32__)
      sprintf(base_name,"_%s%d_XXXXXX",mumps_base,*_myid);
#else
      sprintf(base_name,"_%s%d",mumps_base,*_myid);
#endif
      mumps_ooc_file_prefix=(char *)malloc((strlen(SEPARATOR)+strlen(tmp_dir)+strlen(tmp_fname)+strlen(base_name)+1+1)*sizeof(char));
      if(mumps_ooc_file_prefix==NULL){
	sprintf(error_str,"Allocation problem in low-level OOC layer\n");
	return -13;
      }
      sprintf(mumps_ooc_file_prefix,"%s%s%s%s",tmp_dir,SEPARATOR,tmp_fname,base_name);
  }else{
#if ! defined(_WIN32) || defined(__MINGW32__)
    sprintf(base_name,"%s%s%d_XXXXXX",SEPARATOR,mumps_base,*_myid);
#else
    sprintf(base_name,"%s%s%d",SEPARATOR,mumps_base,*_myid);
#endif
      mumps_ooc_file_prefix=(char *)malloc((strlen(SEPARATOR)+strlen(tmp_dir)+strlen(base_name)+1)*sizeof(char));
      if(mumps_ooc_file_prefix==NULL){
	sprintf(error_str,"Allocation problem in low-level OOC layer\n");
	return -13;
      }
      sprintf(mumps_ooc_file_prefix,"%s%s%s",tmp_dir,SEPARATOR,base_name);
  }  
  if(!dir_flag){
    free(tmp_dir);
  }
  if(!file_flag){
    free(tmp_fname);
  }
  return 0;
}
int mumps_io_alloc_file_struct(int* nb,int which){
  /* allocate structure for files corresponding to factors */
  int i;
  (mumps_files+which)->mumps_io_pfile_pointer_array=(mumps_file_struct *)malloc((*nb)*sizeof(mumps_file_struct));
  if((mumps_files+which)->mumps_io_pfile_pointer_array==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  for(i=0;i<*nb;i++){
    (((mumps_files+which)->mumps_io_pfile_pointer_array)+i)->is_opened=0;
  }
/*   mumps_io_pfile_name=(char **)malloc((*nb)*sizeof(char *)); */
/*   if(mumps_io_pfile_name==NULL){ */
/* #ifndef _WIN32       */
/*     sprintf(error_str,"Allocation problem in low-level OOC layer\n"); */
/* #else */
/*     sprintf(error_str,"Allocation problem in low-level OOC layer\n"); */
/* #endif     */
/*     return -13; */
/*     /\*    return -1;*\/ */
/*   } */
  (mumps_files+which)->mumps_io_current_file_number = -1;
  (mumps_files+which)->mumps_io_last_file_opened=-1;
  (mumps_files+which)->mumps_io_nb_file_opened=0;
  return 0;
}
int mumps_io_get_nb_files(int* nb_files,int* type){
  *nb_files=((mumps_files+*type)->mumps_io_last_file_opened)+1;
  return 0;
}
int mumps_io_get_file_name(int* indice,char* name,int* length,int* type){
  int i;
  i=(*indice)-1;
  strcpy(name,(((mumps_files+*type)->mumps_io_pfile_pointer_array)+i)->name);
  *length=strlen(name)+1;
  return 0;  
}
int mumps_io_alloc_pointers(int* nb_file_type,int * dim){
  int ierr;
  int i;
  /* This is called by solve step, we have only one type of files */
  mumps_io_nb_file_type=*nb_file_type;
  mumps_files=(mumps_file_type *)malloc(mumps_io_nb_file_type*sizeof(mumps_file_type));
  if(mumps_files==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  /* FIXME : This loop must be modified when NB_FILE_TYPE_SOLVE is
     greater than 1 (the dim parameter is not good) */
  for(i=0;i<mumps_io_nb_file_type;i++){
    (mumps_files+i)->mumps_io_nb_file=*(dim+i);
    ierr=mumps_io_alloc_file_struct(dim+i,i);
    if(ierr<0){
      return ierr;
    }
  }
  return 0;
}
int mumps_io_init_vars(int* myid_arg,int* size_element,int* async_arg){
#if ! defined(_WIN32) || defined(__MINGW32__)
  int k211_loc;
  int mumps_flag_open;
#endif
  int i;
  mumps_io_max_file_size=MAX_FILE_SIZE;
  mumps_directio_flag=0;
#if ! defined(_WIN32) || defined(__MINGW32__)
  mumps_flag_open=0;
#endif
  /* must be changed when we will have more than one file type during solve step */
  for(i=0;i<mumps_io_nb_file_type;i++){
#if ! defined(_WIN32) || defined(__MINGW32__)
    (mumps_files+i)->mumps_flag_open=mumps_flag_open|O_RDONLY;
#else
    strcpy((mumps_files+i)->mumps_flag_open,"rb");
#endif
  }
  mumps_io_myid=*myid_arg;
  mumps_elementary_data_size=*size_element;
  mumps_io_flag_async=*async_arg;
  return 0;
}
int mumps_io_set_file_name(int* indice,char* name,int* length,int* type){
  int i;
  i=(*indice)-1;
/*   *(mumps_io_pfile_name+i)=(char *) malloc((*length)*strlen(name)); */
/*   if(*(mumps_io_pfile_name+i)==NULL){ */
/*     sprintf(error_str,"Allocation problem in low-level OOC layer"); */
/*     return -13; */
/*   } */
  strcpy((((mumps_files+*type)->mumps_io_pfile_pointer_array)+i)->name,name);
  return 0;  
}
int mumps_io_open_files_for_read(){
  int i,j;
  mumps_file_struct  *mumps_io_pfile_pointer_array;
#ifdef IRIX64_
  struct dioattr dio;
#endif
  for(j=0;j<mumps_io_nb_file_type;j++){
    mumps_io_pfile_pointer_array=(mumps_files+j)->mumps_io_pfile_pointer_array;
    for(i=0;i<(mumps_files+j)->mumps_io_nb_file;i++){
#if ! defined(_WIN32) || defined(__MINGW32__)
      (mumps_io_pfile_pointer_array+i)->file=open((mumps_io_pfile_pointer_array+i)->name,(mumps_files+j)->mumps_flag_open);
      if((mumps_io_pfile_pointer_array+i)->file==-1){
	mumps_io_build_err_str(errno,-90,"Problem while opening OOC file",error_str,200);
	return -90;
      }
#else
      (mumps_io_pfile_pointer_array+i)->file=fopen((mumps_io_pfile_pointer_array+i)->name,(mumps_files+j)->mumps_flag_open);      
      if((mumps_io_pfile_pointer_array+i)->file==NULL){
	sprintf(error_str,"Problem while opening OOC file");
	return -90;
      }
      (mumps_io_pfile_pointer_array+i)->is_opened=1;
#endif
    }
  }
  return 0;
}
int mumps_io_set_last_file(int* dim,int* type){
  (mumps_files+*type)->mumps_io_last_file_opened=*dim-1;
  (mumps_files+*type)->mumps_io_nb_file_opened=*dim;
  return 0;
}
#if (! defined(_WIN32) || defined(__MINGW32__)) &&  ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
int mumps_io_protect_pointers(){
  pthread_mutex_lock(&mumps_io_pwrite_mutex);
  return 0;
}
int mumps_io_unprotect_pointers(){
  pthread_mutex_unlock(&mumps_io_pwrite_mutex);
  return 0;
}
int mumps_io_init_pointers_lock(){
  pthread_mutex_init(&mumps_io_pwrite_mutex,NULL);
  return 0;
}
int mumps_io_destroy_pointers_lock(){
  pthread_mutex_destroy(&mumps_io_pwrite_mutex);
  return 0;
}
#endif /*WITH_PFUNC*/
#endif /* _WIN32 && WITHOUT_PTHREAD */
 int mumps_io_read__(void * file,void * loc_addr,size_t size,int local_offset,int type){
  int ret_code;
#if ! defined(_WIN32) || defined(__MINGW32__)
  if(!mumps_directio_flag){
    ret_code=mumps_io_read_os_buff__(file,loc_addr, size,local_offset);
    if(ret_code<0){
      return ret_code;
    }
  }
#else
  ret_code=mumps_io_read_win32__(file,loc_addr, size,local_offset);
  if(ret_code<0){
    return ret_code;
  }
#endif  
  return 0;
}
#if ! defined(_WIN32) || defined(__MINGW32__)
int mumps_io_read_os_buff__(void * file,void * loc_addr,size_t size,int local_offset){
  size_t ret_code;
  /* printf("Read with buff %d %d %d\n",(int) size, local_offset,*((int *)file)); */
#ifdef WITH_PFUNC
  ret_code=pread(*(int *)file,loc_addr,size,local_offset);
#else
  lseek(*(int *)file,(long) local_offset,SEEK_SET);
  ret_code=read(*(int *)file,loc_addr,size);
#endif
  if((int) ret_code==-1){
    mumps_io_build_err_str(errno,-90,"Problem with low level read",error_str,200);
    return -90;
  }
  return 0;
}
#endif
#if defined(_WIN32) && ! defined(__MINGW32__)
int mumps_io_read_win32__(void * file,void * loc_addr,size_t size,int local_offset){
  size_t ret_code;
  fseek(*(FILE **)file,(long) local_offset,SEEK_SET);
  ret_code=fread(loc_addr,mumps_elementary_data_size,size,*(FILE **)file);
  if(ret_code!=size){
    sprintf(error_str,"Problem with I/O operation\n");
    return -90;
  }
  return 0;
}
#endif
int mumps_io_write__(void *file, void *loc_addr, size_t write_size, int where,int type){
  int ret_code;
#if ! defined(_WIN32) || defined(__MINGW32__)
  if(!mumps_directio_flag){
    ret_code=mumps_io_write_os_buff__(file,loc_addr, write_size,where);
    if(ret_code<0){
      return ret_code;
    }
  }
#else
  ret_code=mumps_io_write_win32__(file,loc_addr, write_size,where);
  if(ret_code<0){
    return ret_code;
  }
#endif
  return 0;
}
#if ! defined(_WIN32) || defined(__MINGW32__)
int mumps_io_write_os_buff__(void *file, void *loc_addr, size_t write_size, int where){
  size_t ret_code;
  /* printf("write with buff %d %d %d\n",(int) write_size, where,*((int *)file)); */
#ifdef WITH_PFUNC
  ret_code=pwrite(*(int *)file,loc_addr,write_size,where);
#else
  /*in this case all the I/O's are made by the I/O thread => we don't
    need to protect the file pointer.*/
  lseek(*(int *)file,(long)where,SEEK_SET); 
  ret_code=write(*(int *)file,loc_addr,write_size);
#endif
  if((int)ret_code==-1){
    mumps_io_build_err_str(errno,-90,"Problem with low level write",error_str,200);
    return -90;
  }else{
    if(ret_code!=write_size){
      sprintf(error_str,"Error not enough space on disk \n");
      return -90;
    }
  }
  return 0;
}
#endif
#if defined(_WIN32) && ! defined(__MINGW32__)
int mumps_io_write_win32__(void *file, void *loc_addr, size_t write_size, int where){
  size_t ret_code;
  fseek(*(FILE **)file,(long)where,SEEK_SET);  
  ret_code=fwrite(loc_addr,mumps_elementary_data_size, write_size,*(FILE**)file);
  if(ret_code!=write_size){
    sprintf(error_str,"Problem with I/O operation\n");
    return -90;
  }
  return 0;
}
#endif
