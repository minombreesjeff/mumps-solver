C
C   THIS FILE IS PART OF MUMPS VERSION 4.6.3
C   This Version was built on Thu Jun 22 13:22:44 2006
C
C
C  This version of MUMPS is provided to you free of charge. It is public
C  domain, based on public domain software developed during the Esprit IV
C  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. 
C  Since this first public domain version in 1999, the developments are
C  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
C  INRIA.
C
C  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
C  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.
C
C  Up-to-date copies of the MUMPS package can be obtained
C  from the Web pages http://www.enseeiht.fr/apo/MUMPS/
C  or http://graal.ens-lyon.fr/MUMPS
C
C
C   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
C   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
C
C
C  User documentation of any code that uses this software can
C  include this complete notice. You can acknowledge (using
C  references [1], [2], and [3] the contribution of this package
C  in any scientific publication dependent upon the use of the
C  package. You shall use reasonable endeavours to notify
C  the authors of the package of this publication.
C
C   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent (1998),
C   Multifrontal parallel distributed symmetric and unsymmetric solvers,
C   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).
C
C   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
C   A fully asynchronous multifrontal solver using distributed dynamic
C   scheduling, SIAM Journal of Matrix Analysis and Applications,
C   Vol 23, No 1, pp 15-41 (2001).
C
C   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
C   S. Pralet, Hybrid scheduling for the parallel solution of linear
C   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
C
      SUBROUTINE DMUMPS_301( id)
      USE DMUMPS_STRUC_DEF
      USE DMUMPS_BUFFER
      USE DMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'mpif.h'
#if defined(V_T)
      INCLUDE 'VT.inc'
#endif
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      TYPE (DMUMPS_STRUC), TARGET :: id
      INTEGER MP,LP, MPG
      LOGICAL PROK, PROKG
      INTEGER MTYPE, ICNTL20, ICNTL21
      LOGICAL LSCAL, ERANAL, GIVSOL
      INTEGER ICNTL10, ICNTL11
      INTEGER I,K,JPERM, J, II
      INTEGER IZ, NZ_THIS_BLOCK, IRHS_PTR_BEG, SHIFT_PTR
      INTEGER LIW,LIWW
      INTEGER LA, LWC, ISRW2, LIW_PASSED, LA_PASSED
      INTEGER DMUMPS_LBUF, DMUMPS_LBUF_INT
      INTEGER NITREF, NOITER, SOLVET, KASE, JOBIREF
      DOUBLE PRECISION RSOL(1), ONE
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0, ONE = 1.0D0 )
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS_MUMPS
      DOUBLE PRECISION, ALLOCATABLE :: SAVERHS(:), RW1(:),
     *                                 RW2(:), SRW1(:), SRW3(:),
     *                                 Y(:), W(:), D(:), 
     *                                 RHSSOL(:)
      INTEGER, DIMENSION(:), ALLOCATABLE :: UNS_PERM_INV
      INTEGER LIWK_SOLVE
      INTEGER, ALLOCATABLE :: IW1(:), IWK_SOLVE(:)
      INTEGER, POINTER :: N, NZ, MAXIS, MAXS
      INTEGER, POINTER :: NRHS, LRHS
      DOUBLE PRECISION, DIMENSION(:), POINTER :: CNTL
      INTEGER, DIMENSION (:), POINTER :: KEEP,ICNTL,INFO
      INTEGER*8, DIMENSION (:), POINTER :: KEEP8
      INTEGER, DIMENSION (:), POINTER :: IS
      DOUBLE PRECISION, DIMENSION(:),POINTER::S
      DOUBLE PRECISION, DIMENSION(:),POINTER::   RINFOG
      type scaling_data_t
        SEQUENCE
        DOUBLE PRECISION, dimension(:), pointer :: SCALING
        DOUBLE PRECISION, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      type (scaling_data_t) :: scaling_data
      DOUBLE PRECISION ARRET
      DOUBLE PRECISION DUMMY(1)
      INTEGER IDUMMY(1), JDUMMY(1)
      INTEGER allocok, PERLU
      INTEGER NBRHS, NBRHS_EFF, BEG_RHS, END_RHS, 
     &        IBEG, LD_RHS, KDEC, MINLWC, NBRHSMIN, 
     &        MASTER_ROOT
      INTEGER NB_K133
      INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE_COPY
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS_SPARSE_COPY
#if defined(V_T)
      INTEGER soln_drive_class, glob_comm_ini, perm_scal_ini, soln_dist,
     &        soln_assem, perm_scal_post
#endif
      LOGICAL I_AM_SLAVE, MAXSUPDATED
      INTEGER MAXSSAVED
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275
      N    =>id%N
      NZ   =>id%NZ
      NRHS =>id%NRHS
      LRHS =>id%LRHS
      MAXIS=>id%MAXIS
      MAXS =>id%MAXS
      CNTL =>id%CNTL
      KEEP =>id%KEEP
      KEEP8=>id%KEEP8
      IS   =>id%IS
      S    =>id%S
      ICNTL=>id%ICNTL
      INFO =>id%INFO
      RINFOG =>id%RINFOG
      MAXSSAVED   = MAXS
      MAXSUPDATED = .FALSE.
      NB_K133     = 3     ! Max level of recursivity to evaluate MINLWC
#if defined(V_T)
      CALL VTCLASSDEF( 'Soln driver',soln_drive_class,ierr)
      CALL VTFUNCDEF( 'glob_comm_ini',soln_drive_class,
     &     glob_comm_ini,ierr)
      CALL VTFUNCDEF( 'perm_scal_ini',soln_drive_class,
     &     perm_scal_ini,ierr)
      CALL VTFUNCDEF( 'soln_dist',soln_drive_class,soln_dist,ierr)
      CALL VTFUNCDEF( 'soln_assem',soln_drive_class,soln_assem,ierr)
      CALL VTFUNCDEF( 'perm_scal_post',soln_drive_class,
     &     perm_scal_post,ierr)
#endif
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     *             ( id%MYID .eq. MASTER .AND.
     *               KEEP(46) .eq. 1 ) )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LP  = ICNTL( 1 )
      PROK  = (MP.GT.0)
      PROKG = (MPG.GT.0 .and. id%MYID.eq.MASTER)
      IF ( PROK  ) WRITE(MP,100)
      IF ( PROKG ) WRITE(MPG,100)
      IF (id%MYID .eq. MASTER) THEN
        NBRHS = KEEP(84)
        IF (NBRHS .GT. NRHS ) NBRHS = NRHS
      ENDIF
#if defined(V_T)
      CALL VTBEGIN(glob_comm_ini,ierr)
#endif
      CALL MPI_BCAST(NRHS,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(NBRHS,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
        IF (id%KEEP(201).NE.0) THEN
          IF (id%MYID .eq. MASTER) THEN
             KEEP(107) = MAX(0,KEEP(107))
             IF ((KEEP(107).EQ.0).AND.(KEEP(209).EQ.0).AND.
     &            (KEEP(204).EQ.0).AND.(KEEP(211).NE.1) ) THEN
              KEEP(209) = 1
             ENDIF
          ENDIF
          CALL MPI_BCAST( KEEP(107), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(204), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(208), 2, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
        ENDIF
      IF ( I_AM_SLAVE ) THEN
        MINLWC = NB_K133*KEEP(133)
        IF (id%KEEP(201).NE.0) THEN
           IF ( (KEEP(209).GT.0)  .AND. (KEEP(209).LT.MAXS) ) THEN 
              MAXSUPDATED = .TRUE.
              MAXS = MAX(KEEP(209), MINLWC*NBRHS + KEEP(203)+1)
              IF ( ASSOCIATED (id%S) ) THEN
                 DEALLOCATE (id%S)
                 NULLIFY (id%S)
              ENDIF
           ELSE 
              KEEP(209) = 0
           ENDIF
           IF (.NOT. ASSOCIATED(id%S)) THEN
               ALLOCATE (id%S(MAXS), stat = allocok)
               IF ( allocok .GT. 0 ) THEN
                WRITE(*,*) ' Problem reallocation of MAXS '
                INFO(1) = -13
                INFO(2) = MAXS
               ENDIF
           ENDIF
           S => id%S
        ENDIF
        IF ( KEEP( 38 ) .NE. 0 .OR. KEEP( 20 ) .NE. 0 ) THEN 
          IF ( KEEP( 38 ) .ne. 0 ) THEN
            MASTER_ROOT = DMUMPS_275(id%STEP( KEEP(38)),
     *                    id%PROCNODE_STEPS, id%NSLAVES )
          ELSE
            MASTER_ROOT = DMUMPS_275(id%STEP(KEEP(20)),
     *                    id%PROCNODE_STEPS, id%NSLAVES )
          END IF
          IF ( MASTER_ROOT .eq. id%MYID_NODES ) THEN
              NB_K133 = NB_K133 + 1
          END IF
        ENDIF
        IF(KEEP(201).EQ.0)THEN
           LA  = KEEP(31)
           LWC = MAXS - LA
           IF (LWC.LT.MINLWC) THEN
             INFO(1) = -11
             INFO(2) = MINLWC + KEEP(31) - MAXS
           ENDIF
        ELSE
           IF ( id%MAXS .LT. MINLWC + KEEP(203) ) THEN
             INFO(1) = -11
             INFO(2) = MINLWC + KEEP(203) - id%MAXS
           ENDIF
           LWC = MIN(id%MAXS - KEEP(203),
     *               MAX( INT(DBLE(id%MAXS*0.05)), MINLWC ))
           LA  = id%MAXS-LWC
           IF(LA.GT.(KEEP(31)+KEEP(203)*(KEEP(107)+1)))THEN
              LA=KEEP(31)+KEEP(203)*(KEEP(107)+1)
              LWC=id%MAXS-LA
           ENDIF
        ENDIF
        IF ( LWC .LT. NBRHS*MINLWC ) THEN
           NBRHS =  LWC / MINLWC
        ENDIF
      ENDIF
      CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF (INFO(1) < 0) RETURN
      CALL MPI_ALLREDUCE(NBRHS, NBRHSMIN, 1, MPI_INTEGER, MPI_MIN,
     &                   id%COMM, IERR)
      NBRHS = NBRHSMIN
      IF ( NBRHS .LT. 1) THEN
         INFO(1) = -11
         INFO(2) = 3*KEEP(133)
         RETURN
      ENDIF
      IF ( I_AM_SLAVE ) THEN
         IF (KEEP(201).NE.0) THEN
          CALL DMUMPS_590(LA)
          CALL DMUMPS_586(id)
          IF (INFO(1).LT.0) THEN
             GOTO 111
          ENDIF
        ENDIF
      ENDIF
      ALLOCATE(RHSSOL(N*NBRHS),stat = allocok)
            IF ( allocok .GT. 0 ) THEN
              WRITE(*,*) ' Problem in solve: error allocating RHSSOL'
              INFO(1) = -13
              INFO(2) = N*NBRHS
              GOTO 111
            END IF
      IF (id%MYID .eq. MASTER) THEN
        MTYPE = ICNTL(  9 )
        IF ( PROKG )  WRITE( MPG, 150 )
     *                NRHS, NBRHS, ICNTL(9), ICNTL(10), ICNTL(11),
     *                ICNTL(20), ICNTL(21)
        ICNTL20 = ICNTL(20)
        ICNTL21 = ICNTL(21)
        IF (ICNTL20.ne.0.and.ICNTL20.ne.1) ICNTL20=0
        IF (ICNTL21.ne.0.and.ICNTL21.ne.1) ICNTL21=0
        LSCAL = (((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 6)) .OR. (
     *    KEEP(52) .EQ. -1) .OR. KEEP(52) .EQ. -2)
        ERANAL = ((ICNTL(11) .GT. 0) .OR. (ICNTL(10) .GT. 0))
        IF ( (KEEP(55).eq.0) .AND. KEEP(54).eq.0 .AND. 
     &      .NOT.ASSOCIATED(id%A) ) THEN
          ICNTL10 = 0
          ICNTL11 = 0
          ERANAL = .FALSE.
        ELSE
          ICNTL10 = ICNTL(10)
          ICNTL11 = ICNTL(11)
        ENDIF
        IF ((ERANAL .AND. NBRHS > 1) .OR. ICNTL(21) > 0) THEN
          IF (ICNTL(11) > 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     *     ' WARNING: ICNTL(11) treated as if set to zero'
            ICNTL11=0
          ENDIF
          IF (ICNTL(10) > 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     *     ' WARNING: ICNTL(10) treated as if set to zero'
            ICNTL10=0
          ENDIF
          ERANAL = .FALSE.     
        ENDIF
        IF (ERANAL) THEN
          ALLOCATE(SAVERHS(N*NBRHS),stat = allocok)
            IF ( allocok .GT. 0 ) THEN
              WRITE(*,*) ' Problem in solve: error allocating SAVERHS'
              INFO(1) = -13
              INFO(2) = N*NBRHS
              GOTO 111
            END IF
        ENDIF
      END IF
      CALL MPI_BCAST(ICNTL10,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL11,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL20,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL21,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(ERANAL,1,MPI_LOGICAL,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(LSCAL,1,MPI_LOGICAL,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(MTYPE,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
        id%LBUFR_BYTES = ( 4 + KEEP(133) ) * KEEP(34) +
     *                      KEEP(133) * NBRHS * KEEP(35)
        id%LBUFR = ( id%LBUFR_BYTES * KEEP(34) - 1 ) / KEEP(34)
        IF ( ASSOCIATED (id%BUFR) ) DEALLOCATE(id%BUFR)
        ALLOCATE (id%BUFR(id%LBUFR),stat=allocok)
        IF ( allocok .GT. 0 ) THEN
              WRITE(*,*) ' Problem in solve: error allocating BUFR'
              INFO(1) = -13
              INFO(2) = id%LBUFR
              GOTO 111
        ENDIF
      IF ( I_AM_SLAVE ) THEN
        DMUMPS_LBUF = id%LBUFR_BYTES * id%NSLAVES
        DMUMPS_LBUF_INT = ( 20 + id%NSLAVES * id%NSLAVES  * 4 )
     *                 * KEEP(34)
        CALL DMUMPS_55( DMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          INFO(1) = -13
          INFO(2) = DMUMPS_LBUF_INT
          IF ( LP .GT. 0 ) THEN
            WRITE(LP,*) id%MYID,
     *      ':Error allocating small Send buffer:IERR=',IERR
          END IF
          GOTO 111
        END IF
        CALL DMUMPS_53( DMUMPS_LBUF, IERR )
        IF ( IERR .NE. 0 ) THEN
          INFO(1) = -13
          INFO(2) = DMUMPS_LBUF
          IF ( LP .GT. 0 ) THEN
            WRITE(LP,*) id%MYID,
     *      ':Error allocating Send buffer:IERR=', IERR
          END IF
          GOTO 111
        END IF
      ENDIF
      IF ( id%MYID .NE. MASTER .or.
     *     (id%MYID .EQ. MASTER .AND.
     *      icntl20.ne.0 .and. icntl21.ne.0 )) THEN
        ALLOCATE(RHS_MUMPS(N*NBRHS),stat=ierr)
        IF ( ierr .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=N*NBRHS
          IF (LP > 0)
     *      WRITE(LP,*) 'ERREUR while allocating RHS on a slave'
          GOTO 111
        END IF
      ELSE
        RHS_MUMPS=>id%RHS
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        LIWK_SOLVE = 5 * KEEP(28) + 2
        ALLOCATE ( IWK_SOLVE( LIWK_SOLVE), stat = allocok )
        IF (allocok .GT. 0 ) THEN
         INFO(1)=-13
         INFO(2)=LIWK_SOLVE
         GOTO 111
        END IF
        LIW = KEEP(32)
        LIWW  = MAXIS - LIW
        IF ( LIWW .LT. 0 ) THEN
          IF (LP>0) THEN
            WRITE(LP,*) 'ERROR in DMUMPS_301: LIWW<0'
          ENDIF
          CALL DMUMPS_ABORT()
        END IF
        ISRW2   = LA + 1
        ALLOCATE(SRW3(KEEP(133)), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=KEEP(133)
          GOTO 111
        END IF
      END IF
 111  CONTINUE
#if defined(V_T)
      CALL VTEND(glob_comm_ini,ierr)
#endif
      CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF (INFO(1) .LT.0 ) GOTO 90
      IF ( ICNTL21==1 ) THEN
        IF (LSCAL) THEN
          IF (id%MYID.NE.MASTER) THEN
            IF (MTYPE == 1) THEN
              ALLOCATE(id%COLSCA(id%N),stat=allocok)
            ELSE
              ALLOCATE(id%ROWSCA(id%N),stat=allocok)
            ENDIF
            IF (allocok > 0) THEN
              IF (LP > 0) THEN
                WRITE(LP,*) 'Error allocating temporary scaling array'
              ENDIF
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 40
            ENDIF
          ENDIF
          IF (MTYPE == 1) THEN
              CALL MPI_BCAST(id%COLSCA,id%N,
     *                       MPI_DOUBLE_PRECISION,MASTER,
     *                       id%COMM,IERR)
              scaling_data%SCALING=>id%COLSCA
          ELSE
              CALL MPI_BCAST(id%ROWSCA,id%N,
     *                       MPI_DOUBLE_PRECISION,MASTER,
     *                       id%COMM,IERR)
              scaling_data%SCALING=>id%ROWSCA
          ENDIF
          IF (I_AM_SLAVE) THEN
            ALLOCATE(scaling_data%SCALING_LOC(id%KEEP(89)),
     *               stat=allocok)
            IF (allocok > 0) THEN
              IF (LP > 0) THEN
                WRITE(LP,*) 'Error allocating local scaling array'
              ENDIF
              INFO(1)=-13
              INFO(2)=id%KEEP(89)
              GOTO 40
            ENDIF
          ENDIF
        ENDIF
        IF ( I_AM_SLAVE ) THEN
          LIW_PASSED=MAX(1,LIW)
          IF (KEEP(89) .GT. 0) THEN
            CALL DMUMPS_535( MTYPE, id%ISOL_LOC(1),
     *               id%LSOL_LOC, id%PTLUST_S(1),
     *               id%KEEP(1),id%KEEP8(1),
     *               id%IS(1), LIW_PASSED,id%MYID_NODES,
     *               id%N, id%STEP(1), id%PROCNODE_STEPS(1),
     *               id%NSLAVES, scaling_data, LSCAL )
          ENDIF
          IF (id%MYID.NE.MASTER .AND. LSCAL) THEN
            IF (MTYPE == 1) THEN
              DEALLOCATE(id%COLSCA)
              NULLIFY(id%COLSCA)
            ELSE
              DEALLOCATE(id%ROWSCA)
              NULLIFY(id%ROWSCA)
            ENDIF
          ENDIF
        ENDIF
        IF (KEEP(23) .NE. 0 .AND. MTYPE==1) THEN
          IF (id%MYID.NE.MASTER) THEN
            ALLOCATE(id%UNS_PERM(id%N),stat=allocok)
            IF (allocok > 0) THEN
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 40
            ENDIF
          ENDIF
        ENDIF
 40     CONTINUE
        CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
        IF (INFO(1) .LT.0 ) GOTO 90
        IF (KEEP(23) .NE. 0 .AND. MTYPE==1) THEN
          CALL MPI_BCAST(id%UNS_PERM,id%N,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
          IF (I_AM_SLAVE) THEN
            DO I=1, KEEP(89)
              id%ISOL_LOC(I) = id%UNS_PERM(id%ISOL_LOC(I))
            ENDDO
          ENDIF
          IF (id%MYID.NE.MASTER) THEN
            DEALLOCATE(id%UNS_PERM)
            NULLIFY(id%UNS_PERM)
          ENDIF
        ENDIF
      ENDIF
      DO BEG_RHS=1, NRHS, NBRHS
        NBRHS_EFF    = MIN(NRHS-BEG_RHS+1, NBRHS)
        END_RHS      = BEG_RHS + NBRHS_EFF - 1
        IF (ICNTL20.eq.0 .or. ICNTL21.eq.0) THEN
          IF (id%MYID .eq. MASTER) THEN
            IF (NRHS.GT.1) THEN 
              LD_RHS    = LRHS
            ELSE
              LD_RHS    = N
            ENDIF
            IBEG      = (BEG_RHS-1) * LD_RHS + 1
          ELSE 
            LD_RHS    = N
            IBEG      = 1
          END IF
        ELSE
          LD_RHS = N
          IBEG   = 1
        ENDIF
#if defined(V_T)
      CALL VTBEGIN(perm_scal_ini,ierr)
#endif
      IF (id%MYID .eq. MASTER) THEN
        IF (ICNTL20==1) THEN
          NZ_THIS_BLOCK=id%IRHS_PTR(BEG_RHS+NBRHS_EFF)-
     *                    id%IRHS_PTR(BEG_RHS)
          IF (KEEP(23) .NE. 0 .and. MTYPE .NE. 1) THEN
            ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK),stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=NZ_THIS_BLOCK
              GOTO 30
            endif
          ELSE
            IRHS_SPARSE_COPY
     *      =>
     *            id%IRHS_SPARSE(id%IRHS_PTR(BEG_RHS):
     *                        id%IRHS_PTR(BEG_RHS)+NZ_THIS_BLOCK-1)
          ENDIF
          IF (LSCAL) THEN
            ALLOCATE(RHS_SPARSE_COPY(NZ_THIS_BLOCK), stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=NZ_THIS_BLOCK
              GOTO 30
            endif
          ELSE
            RHS_SPARSE_COPY
     *         => id%RHS_SPARSE(id%IRHS_PTR(BEG_RHS):
     *                       id%IRHS_PTR(BEG_RHS)+NZ_THIS_BLOCK-1)
          ENDIF
        ENDIF
        IF (KEEP(23) .NE. 0) THEN
          IF (MTYPE .NE. 1) THEN
            IF (ICNTL20==0) THEN
              DO K = 1, NBRHS_EFF
               KDEC = IBEG+(K-1)*LD_RHS
               DO I = 1, N
                RHSSOL(I+(K-1)*N)=RHS_MUMPS(I-1+KDEC)
               END DO
               DO I = 1, N
                JPERM = id%UNS_PERM(I)
                RHS_MUMPS(KDEC+I-1) = RHSSOL(JPERM+(K-1)*N)
               END DO
              END DO
            ELSE
              ALLOCATE(UNS_PERM_INV(N),stat=allocok) !FIXME: build it outside...
              if (allocok .GT.0 ) THEN
                INFO(1)=-13
                INFO(2)=N
                GOTO 30
              endif
              DO I = 1, N
                UNS_PERM_INV(id%UNS_PERM(I))=I
              ENDDO
              DO I = id%IRHS_PTR(BEG_RHS),
     *               id%IRHS_PTR(BEG_RHS+NBRHS_EFF)-1
                JPERM = UNS_PERM_INV(id%IRHS_SPARSE(I))
                IRHS_SPARSE_COPY(I-id%IRHS_PTR(BEG_RHS)+1)=JPERM
              ENDDO
              DEALLOCATE(UNS_PERM_INV) !FIXME: free it oustide
            ENDIF
          ENDIF
        ENDIF
        IF (ERANAL) THEN
         IF ( ICNTL20 == 0 ) THEN
          DO K = 1, NBRHS_EFF
            KDEC = IBEG+(K-1)*LD_RHS
            DO I = 1, N
              SAVERHS(I+(K-1)*N) = RHS_MUMPS(KDEC+I-1)
            END DO
          ENDDO
         ENDIF
        ENDIF
        IF (LSCAL) THEN
         IF (ICNTL20==0) THEN
          IF (MTYPE .EQ. 1) THEN
            DO K =1, NBRHS_EFF 
             KDEC = (K-1) * LD_RHS + IBEG - 1
             DO I = 1, N
              RHS_MUMPS(KDEC+I) = RHS_MUMPS(KDEC+I) * id%ROWSCA(I)
             END DO
            ENDDO
          ELSE
            DO K =1, NBRHS_EFF 
             KDEC = (K-1) * LD_RHS + IBEG - 1
             DO I = 1, N
              RHS_MUMPS(KDEC+I) = RHS_MUMPS(KDEC+I) * id%COLSCA(I)
             END DO
            ENDDO
          ENDIF
         ELSE
          KDEC=id%IRHS_PTR(BEG_RHS)
          IF (MTYPE .eq. 1) THEN
            DO IZ=1,NZ_THIS_BLOCK
              I=IRHS_SPARSE_COPY(IZ)
              RHS_SPARSE_COPY(IZ)=id%RHS_SPARSE(KDEC+IZ-1)*id%ROWSCA(I)
            ENDDO
          ELSE
            DO IZ=1,NZ_THIS_BLOCK
              I=IRHS_SPARSE_COPY(IZ)
              RHS_SPARSE_COPY(IZ)=id%RHS_SPARSE(KDEC+IZ-1)*id%COLSCA(I)
            ENDDO
          ENDIF
         ENDIF
        END IF
      ENDIF
#if defined(V_T)
      CALL VTEND(perm_scal_ini,ierr)
#endif
 30   CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF (INFO(1) .LT.0 ) GOTO 90
#if defined(V_T)
      CALL VTBEGIN(soln_dist,ierr)
#endif
      IF (ICNTL20 == 0) THEN
       DO K= 1, NBRHS_EFF
        KDEC = (K-1) * LD_RHS + IBEG 
        CALL MPI_BCAST(RHS_MUMPS(KDEC), N, MPI_DOUBLE_PRECISION,
     *                MASTER,id%COMM,IERR)
       END DO
      ELSE
       CALL MPI_BCAST( NZ_THIS_BLOCK,1, MPI_INTEGER,
     *                 MASTER, id%COMM,IERR)
       IF (id%MYID==MASTER) THEN
         CALL MPI_BCAST(IRHS_SPARSE_COPY,
     *                NZ_THIS_BLOCK,
     *                MPI_INTEGER,
     *                MASTER, id%COMM,IERR)
         CALL MPI_BCAST(RHS_SPARSE_COPY,
     *                NZ_THIS_BLOCK,
     *                MPI_DOUBLE_PRECISION,
     *                MASTER, id%COMM,IERR)
         IRHS_PTR_BEG=BEG_RHS
         CALL MPI_BCAST(id%IRHS_PTR(BEG_RHS),
     *                NBRHS_EFF+1,
     *                MPI_INTEGER,
     *                MASTER, id%COMM,IERR)
       ELSE
         ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK))
         ALLOCATE(RHS_SPARSE_COPY(NZ_THIS_BLOCK))
         ALLOCATE(id%IRHS_PTR(NBRHS_EFF+1))
         CALL MPI_BCAST(IRHS_SPARSE_COPY,
     *                NZ_THIS_BLOCK,
     *                MPI_INTEGER,
     *                MASTER, id%COMM,IERR)
         CALL MPI_BCAST(RHS_SPARSE_COPY,
     *                NZ_THIS_BLOCK,
     *                MPI_DOUBLE_PRECISION,
     *                MASTER, id%COMM,IERR)
         IRHS_PTR_BEG=1
         CALL MPI_BCAST(id%IRHS_PTR(IRHS_PTR_BEG),
     *                NBRHS_EFF+1,
     *                MPI_INTEGER,
     *                MASTER, id%COMM,IERR)
       ENDIF
       SHIFT_PTR=id%IRHS_PTR(IRHS_PTR_BEG)-1
       IF ( I_AM_SLAVE ) THEN
         DO K = 1, NBRHS_EFF
          KDEC = (K-1) * LD_RHS + IBEG - 1
          RHS_MUMPS(KDEC+1:KDEC+id%N)=DBLE(ZERO)
          DO IZ=id%IRHS_PTR(IRHS_PTR_BEG+K-1)-SHIFT_PTR,
     *          id%IRHS_PTR(IRHS_PTR_BEG+K)-1-SHIFT_PTR
            I=IRHS_SPARSE_COPY(IZ)
            RHS_MUMPS(KDEC+I)= RHS_SPARSE_COPY(IZ)
          ENDDO
         ENDDO
       END IF
       IF (id%MYID .ne. MASTER) THEN
         DEALLOCATE(IRHS_SPARSE_COPY)
         DEALLOCATE(RHS_SPARSE_COPY)
         DEALLOCATE(id%IRHS_PTR)
         NULLIFY(IRHS_SPARSE_COPY)
         NULLIFY(RHS_SPARSE_COPY)
         NULLIFY(id%IRHS_PTR)
       ELSE
         IF (KEEP(23).ne.0 .and. MTYPE.ne.1) THEN
           DEALLOCATE(IRHS_SPARSE_COPY)
         ENDIF
         NULLIFY(IRHS_SPARSE_COPY)
         IF (LSCAL) THEN
           DEALLOCATE(RHS_SPARSE_COPY)
         ENDIF
         NULLIFY(RHS_SPARSE_COPY)
       ENDIF
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        LIW_PASSED = MAX( LIW, 1 )
        LA_PASSED  = MAX( LA, 1 )
        CALL DMUMPS_245(id%root, N, S(1), LA_PASSED,
     *    IS(1), LIW_PASSED,
     *    S(ISRW2), LWC, IS( LIW + 1 ),
     *    LIWW, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF, RHSSOL,
     *    id%NA,id%LNA,id%NE_STEPS(1), SRW3, MTYPE,
     *    ICNTL, id%STEP(1), id%FRERE_STEPS(1), 
     *    id%DAD_STEPS(1), id%FILS(1),
     *    id%PTLUST_S(1), id%PTRFAC(1),
     *    IWK_SOLVE, LIWK_SOLVE,
     *    id%PROCNODE_STEPS,
     *    id%NSLAVES, INFO, KEEP,KEEP8,
     *    id%COMM, id%COMM_NODES, id%MYID,
     *    id%MYID_NODES,
     *    id%BUFR,
     *    id%LBUFR, id%LBUFR_BYTES, 
     *
     *    id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1))
      END IF
      CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF (INFO(1).eq.-2) then
        INFO(1)=-11
        write(*,*) ' WARNING : -11 error code obtained in solve'
      END IF
      IF (INFO(1).eq.-3) then
        INFO(1)=-14
        write(*,*) ' WARNING : -14 error code obtained in solve'
      END IF
      IF (INFO(1).LT.0) GO TO 90
#if defined(V_T)
      CALL VTBEGIN(soln_assem,ierr)
#endif
      IF (ICNTL21 == 0) THEN
        LIW_PASSED = MAX( LIW, 1 )
        IF ( .NOT.I_AM_SLAVE ) THEN
        CALL DMUMPS_521(id%NSLAVES,id%N, id%MYID, id%COMM,
     *          MTYPE, RHSSOL, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     *          JDUMMY, id%KEEP(1), id%KEEP8(1), id%PROCNODE_STEPS(1),
     *          IDUMMY, 1,
     *          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES)
        ELSE
        CALL DMUMPS_521(id%NSLAVES,id%N, id%MYID, id%COMM,
     *          MTYPE, RHSSOL, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     *          id%PTLUST_S(1), id%KEEP(1), id%KEEP8(1),
     *          id%PROCNODE_STEPS(1),
     *          IS(1), LIW_PASSED,
     *          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES)
        ENDIF
      ELSE
        IF ( I_AM_SLAVE ) THEN
         LIW_PASSED = MAX( LIW, 1 )
         IF ( KEEP(89) .GT. 0 ) THEN
           CALL DMUMPS_532(id%NSLAVES,
     *          id%N, id%MYID,
     *          MTYPE, RHSSOL, NBRHS_EFF,
     *          id%ISOL_LOC(1),
     *          id%SOL_LOC(1), BEG_RHS, id%LSOL_LOC,
     *          id%PTLUST_S(1), id%PROCNODE_STEPS(1),
     *          id%KEEP(1),id%KEEP8(1),
     *          IS(1), LIW_PASSED,
     *          id%STEP(1), scaling_data, LSCAL )
         ENDIF
        ENDIF
      ENDIF
#if defined(V_T)
      CALL VTEND(soln_assem,ierr)
      CALL VTBEGIN(perm_scal_post,ierr)
#endif
      IF ( id%MYID.eq.MASTER) THEN
        IF (ICNTL21==0.AND.LSCAL) THEN
          IF (MTYPE .EQ. 1) THEN
             DO K= 1, NBRHS_EFF
              KDEC = (K-1) * LD_RHS + IBEG - 1
              DO I = 1, N
                RHS_MUMPS(KDEC+ I) = RHS_MUMPS(KDEC+ I) * id%COLSCA(I)
              END DO
             END DO
          ELSE
             DO K= 1, NBRHS_EFF
              KDEC = (K-1) * LD_RHS + IBEG - 1
              DO I = 1, N
                RHS_MUMPS(KDEC+I) = RHS_MUMPS(KDEC+I) * id%ROWSCA(I)
              END DO
             END DO
          ENDIF
        END IF
      END IF
#if defined(V_T)
      CALL VTEND(perm_scal_post,ierr)
#endif
      IF ( ICNTL10 > 0 .AND. NBRHS_EFF > 1 ) THEN
        DO I = 1, ICNTL10
          write(*,*) 'FIXME: to be implemented'
        END DO
      END IF
      IF (ERANAL) THEN
        IF ((ICNTL10 .GT. 0) .AND. (ICNTL11 .GT. 0)) THEN
          IF (id%MYID .EQ. MASTER) THEN
            GIVSOL = .FALSE.
            IF (MP .GT. 0) WRITE( MP, 170 )
            ALLOCATE(RW1(N),stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=N
              GOTO 776
            endif
            ALLOCATE(RW2(N),stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=N
              GOTO 776
            endif
          END IF
          IF ( KEEP(54) .ne. 0 ) THEN
            ALLOCATE( SRW1( N ), stat =allocok )
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=N
            endif
          END IF
 776      CONTINUE
          CALL DMUMPS_276( ICNTL, INFO,
     *                  id%COMM,id%MYID)
          IF ( INFO(1) .LT. 0 ) GOTO 90
          IF ( KEEP(54) .eq. 0 ) THEN
            IF (id%MYID .EQ. MASTER) THEN
              IF (KEEP(55).EQ.0) THEN
                CALL DMUMPS_278( ICNTL(9), N, NZ, id%A(1),
     *             id%IRN(1), id%JCN(1),
     *             RHS_MUMPS(IBEG), SAVERHS, RW1, RW2, KEEP,KEEP8 )
              ELSE
                CALL DMUMPS_121( ICNTL(9), N, 
     *          id%NELT, id%ELTPTR, 
     *          id%LELTVAR, id%ELTVAR,
     *          id%NA_ELT, id%A_ELT,
     *          RHS_MUMPS(IBEG), SAVERHS, RW1, RW2, KEEP,KEEP8 )
              ENDIF
            END IF
          ELSE
            CALL MPI_BCAST( RHS_MUMPS(IBEG), id%N,
     *              MPI_DOUBLE_PRECISION, MASTER,
     *              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_192( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        RHS_MUMPS(IBEG), SRW1, KEEP(50), MTYPE )
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, RW2,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
              RW2 = SAVERHS - RW2
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc,
     *          SRW1, id%KEEP,id%KEEP8 )
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, RW1,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            DEALLOCATE( SRW1 )
          END IF
          IF ( id%MYID .EQ. MASTER )  THEN
            CALL DMUMPS_205(ICNTL(9),INFO(1),N,NZ,
     *        RHS_MUMPS(IBEG), SAVERHS,RW1,RW2,GIVSOL,
     *        RSOL,RINFOG(4),RINFOG(5),RINFOG(6),MP,ICNTL,
     *        KEEP,KEEP8)
            DEALLOCATE(RW1)
            DEALLOCATE(RW2)
          END IF
        END IF
      IF ( PROK  .AND. ICNTL10 .GT. 0 ) WRITE( MP, 270 )
      IF ( PROKG .AND. ICNTL10 .GT. 0 ) WRITE( MPG, 270 )
      ALLOCATE(Y(N), stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        GOTO 777
      ENDIF
      IF ( id%MYID .EQ. MASTER ) THEN
        ALLOCATE( IW1( 2 * N ),stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=2 * N
          GOTO 777
        ENDIF
        ALLOCATE( D(N),stat =allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=N
          GOTO 777
        ENDIF
        ALLOCATE( W(3*N), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=N
          GOTO 777
        ENDIF
        NITREF = ICNTL10
        JOBIREF= ICNTL11
        IF ( PROKG .AND. ICNTL10 .GT. 0 )
     *    WRITE( MPG, 240) 'MAXIMUM NUMBER OF STEPS =', NITREF
        DO I = 1, N
          D( I ) = ONE
        END DO
      END IF
      ALLOCATE(SRW1(N),stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        GOTO 777
      ENDIF
      KASE = 0
 777  CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF ( INFO(1) .LT. 0 ) GOTO 90
  22    CONTINUE
        IF ( KEEP(54) .eq. 0 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            IF ( KASE .eq. 0 ) THEN
              IF (KEEP(55).NE.0) THEN 
               CALL DMUMPS_119(MTYPE, N, 
     *           id%NELT, id%ELTPTR(1), 
     *           id%LELTVAR, id%ELTVAR(1),
     *           id%NA_ELT, id%A_ELT,
     *           W(N+1), KEEP,KEEP8 )
              ELSE
               IF ( MTYPE .eq. 1 ) THEN
                 CALL DMUMPS_207
     *       ( id%A(1), NZ, N, id%IRN(1), id%JCN(1), W(N+1), KEEP,KEEP8)
               ELSE
                 CALL DMUMPS_207
     *       ( id%A(1), NZ, N, id%JCN(1), id%IRN(1), W(N+1), KEEP,KEEP8)
               END IF
              ENDIF
            ENDIF
          END IF
        ELSE
          IF ( KASE .eq. 0 ) THEN
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              IF ( MTYPE .eq. 1 ) THEN
              CALL DMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc,
     *          SRW1, id%KEEP,id%KEEP8 )
              ELSE
              CALL DMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%JCN_loc, id%IRN_loc,
     *          SRW1, id%KEEP,id%KEEP8 )
              END IF
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, W( N + 1 ),
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY, 
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
          END IF
        END IF
        IF ( id%MYID .eq. MASTER ) THEN
            ARRET = CNTL(2)
            IF (ARRET .LT. 0.0D0) THEN
              ARRET = SQRT(EPSILON(0.0D0))
            END IF
            CALL DMUMPS_206(NZ,N,SAVERHS,RHS_MUMPS(IBEG),
     *      Y, D, W,
     *      IW1, KASE,RINFOG(7),
     *      RINFOG(9), JOBIREF, RINFOG(10), NITREF, NOITER, MP,
     *      KEEP,KEEP8, ARRET )
        END IF
        IF ( KEEP(54) .ne. 0 ) THEN
          CALL MPI_BCAST( KASE, 1, MPI_INTEGER, MASTER,
     *    id%COMM, IERR )
        END IF
        IF ( KEEP(54) .eq. 0 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            IF ( KASE .eq. 14 ) THEN
              IF (KEEP(55).NE.0) THEN
               CALL DMUMPS_122( MTYPE, N, 
     *            id%NELT, id%ELTPTR, id%LELTVAR,
     *            id%ELTVAR, id%NA_ELT, id%A_ELT,
     *            SAVERHS, RHS_MUMPS(IBEG),
     *            Y, W, KEEP(50))
              ELSE
                 IF ( MTYPE .eq. 1 ) THEN
                   CALL DMUMPS_208
     *    (id%A(1), NZ, N, id%IRN(1), id%JCN(1), SAVERHS,
     *    RHS_MUMPS(IBEG), Y, W, KEEP,KEEP8)
                 ELSE
                   CALL DMUMPS_208
     *    (id%A(1), NZ, N, id%JCN(1), id%IRN(1), SAVERHS,
     *    RHS_MUMPS(IBEG), Y, W, KEEP,KEEP8)
                 END IF
              ENDIF
              GOTO 22
            END IF
          END IF
        ELSE
          IF ( KASE.eq.14 ) THEN
            CALL MPI_BCAST( RHS_MUMPS(IBEG), id%N,
     *              MPI_DOUBLE_PRECISION, MASTER,
     *              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_192( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        RHS_MUMPS(IBEG), SRW1, KEEP(50), MTYPE )
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, Y,
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
              Y = SAVERHS - Y
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY, 
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and. id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_193( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        RHS_MUMPS(IBEG), SRW1, KEEP(50), MTYPE )
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, W,
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY, 
     *          id%N, MPI_DOUBLE_PRECISION,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            GOTO 22
          END IF
        END IF
      IF ( id%MYID .eq. MASTER ) THEN
        IF ( KASE .GT. 0 ) THEN
          IF ( MTYPE .EQ. 1 ) THEN
            SOLVET = KASE - 1
          ELSE
            SOLVET = KASE
          END IF
          IF ( LSCAL ) THEN
            IF ( SOLVET .EQ. 1 ) THEN
              DO K = 1, N
                Y( K ) = Y( K ) * id%ROWSCA( K )
              END DO
            ELSE
              DO K = 1, N
                Y( K ) = Y( K ) * id%COLSCA( K )
              END DO
            END IF
          END IF
        END IF
      END IF
      CALL MPI_BCAST( KASE , 1, MPI_INTEGER, MASTER,
     *                id%COMM, IERR)
      CALL MPI_BCAST( SOLVET, 1, MPI_INTEGER, MASTER,
     *                id%COMM, IERR)
      IF ( KASE .GT. 0 ) THEN
        CALL MPI_BCAST( Y, N, MPI_DOUBLE_PRECISION, MASTER,
     *                id%COMM, IERR )
        IF ( I_AM_SLAVE ) THEN
          LIW_PASSED = MAX( LIW, 1 )
          LA_PASSED = MAX( LA, 1 )
          CALL DMUMPS_245( id%root, N,
     *    S(1), LA_PASSED, id%IS( 1 ),
     *    LIW_PASSED, S( ISRW2 ), LWC, id%IS( LIW + 1 ),
     *    LIWW,
     *    Y, N, NBRHS_EFF, RHSSOL, id%NA, id%LNA, id%NE_STEPS,
     *    SRW1, SOLVET, ICNTL,
     *    id%STEP(1), id%FRERE_STEPS(1), id%DAD_STEPS(1), 
     *    id%FILS(1),
     *    id%PTLUST_S(1), id%PTRFAC(1),
     *    IWK_SOLVE, LIWK_SOLVE,
     *    id%PROCNODE_STEPS, id%NSLAVES, INFO, KEEP,KEEP8,
     *    id%COMM,
     *    id%COMM_NODES,
     *    id%MYID, id%MYID_NODES,
     *    id%BUFR, id%LBUFR, id%LBUFR_BYTES , 
     *
     *    id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1))
        END IF
        CALL DMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
        IF (INFO(1).eq.-2) INFO(1)=-12
        IF (INFO(1).eq.-3) INFO(1)=-15
        IF (INFO(1).LT.0) GO TO 90
        LIW_PASSED = MAX( LIW, 1 )
        IF ( .NOT. I_AM_SLAVE ) THEN
        CALL DMUMPS_521(id%NSLAVES,id%N, id%MYID, id%COMM,
     *          MTYPE, RHSSOL, Y, LD_RHS, NBRHS_EFF,
     *          JDUMMY, id%KEEP(1),id%KEEP8(1), id%PROCNODE_STEPS(1),
     *          IDUMMY, 1,
     *          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES)
        ELSE
        CALL DMUMPS_521(id%NSLAVES,id%N, id%MYID, id%COMM,
     *          MTYPE, RHSSOL, Y, LD_RHS, NBRHS_EFF,
     *          id%PTLUST_S(1), id%KEEP(1),id%KEEP8(1),
     *          id%PROCNODE_STEPS(1),
     *          IS(1), LIW_PASSED,
     *          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES)
        ENDIF
        IF ( id%MYID.eq.MASTER) THEN
          IF (LSCAL) THEN
            IF (SOLVET .EQ. 1) THEN
               DO K = 1, N
                 Y(K) = Y(K) * id%COLSCA(K)
               END DO
            ELSE
               DO K = 1, N
                 Y(K) = Y(K) * id%ROWSCA(K)
               END DO
            ENDIF
          END IF
        END IF
        GO TO 22
      ELSEIF ( KASE .LT. 0 ) THEN
           INFO( 1 ) = INFO( 1 ) + 8
      END IF
      IF ( id%MYID .eq. MASTER ) DEALLOCATE(W)
      IF ( id%MYID .eq. MASTER ) DEALLOCATE(IW1)
      IF ( PROKG .AND. NITREF .GT. 0 .AND.
     *id%MYID .EQ. MASTER ) THEN
        WRITE( MPG, 81 ) 
        WRITE( MPG, 141 ) 'NUMBER OF STEPS OF ITERATIVE REFINEMENTS  
     *=', NOITER
      ENDIF
      IF ( id%MYID .EQ. MASTER .AND. NITREF .GT. 0 ) THEN
        id%INFOG(15) = NOITER
      END IF
      IF ( PROK .AND. NITREF .GT.0 ) WRITE( MP, 131 ) 
      IF (ICNTL11 .GT. 0) THEN
        IF ( KEEP(54) .eq. 0 ) THEN
          IF (id%MYID .EQ. MASTER) THEN
            IF (KEEP(55).EQ.0) THEN
              CALL DMUMPS_278( MTYPE, N, NZ, id%A(1),
     *          id%IRN(1), id%JCN(1),
     *          RHS_MUMPS(IBEG), SAVERHS, Y, D, KEEP,KEEP8 )
            ELSE
              CALL DMUMPS_121( MTYPE, N, 
     *          id%NELT, id%ELTPTR, 
     *          id%LELTVAR, id%ELTVAR,
     *          id%NA_ELT, id%A_ELT,
     *          RHS_MUMPS(IBEG), SAVERHS, Y, D, KEEP,KEEP8 )
            ENDIF
          END IF
        ELSE
            CALL MPI_BCAST( RHS_MUMPS(IBEG), id%N,
     *              MPI_DOUBLE_PRECISION, MASTER, 
     *              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_192( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        RHS_MUMPS(IBEG), SRW1, KEEP(50), MTYPE ) 
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, D,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
              D = SAVERHS - D
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL DMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc,
     *          SRW1, id%KEEP,id%KEEP8 )
            ELSE
              SRW1 = DBLE(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, Y,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_PRECISION,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
        END IF
        IF (id%MYID .EQ. MASTER) THEN
         IF ((MPG .GT. 0) .AND. (NITREF .GT. 0)) WRITE( MPG, 65 ) 
         IF ((MPG .GT. 0) .AND. (NITREF .LE. 0)) WRITE( MPG, 170 ) 
         GIVSOL = .FALSE.
         CALL DMUMPS_205(MTYPE,INFO(1),N,NZ,RHS_MUMPS(IBEG),
     *        SAVERHS,Y,D,GIVSOL,
     *        RSOL,RINFOG(4),RINFOG(5),RINFOG(6),MPG,ICNTL,
     *        KEEP,KEEP8)
         IF ( MPG .GT. 0 ) THEN
          WRITE( MPG, 115 )
     *'RINFOG(7):COMPONENTWISE SCALED RESIDUAL(W1)=', RINFOG(7)
          WRITE( MPG, 115 )
     *'------(8):---------------------------- (W2)=', RINFOG(8)
          WRITE( MPG, 115 )
     *'------(9):Upper bound ERROR ...............=', RINFOG(9)
          WRITE( MPG, 115 )
     *'-----(10):CONDITION NUMBER (1) ............=', RINFOG(10)
          WRITE( MPG, 115 )
     *'-----(11):CONDITION NUMBER (2) ............=', RINFOG(11)
         END IF
        END IF ! MASTER
      END IF ! ICNTL11 > 0
      IF (id%myid == MASTER) DEALLOCATE(D)
      DEALLOCATE(Y)
      DEALLOCATE(SRW1)
      END IF
      IF ( id%MYID .EQ. MASTER .and.ICNTL(21)==0) THEN
        IF ((ICNTL(9) .EQ. 1) .AND. (KEEP(23) .NE. 0)) THEN
          ALLOCATE( RW1( N ),stat =allocok )
          IF ( allocok .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=N
            WRITE(*,*) 'could not allocate ', N, 'integers.'
            CALL DMUMPS_ABORT()
          END IF
          DO K = 1, NBRHS_EFF
           KDEC = (K-1)*LD_RHS+IBEG-1
           DO 70 I = 1, N
            RW1(I) = RHS_MUMPS(KDEC+I)
 70        CONTINUE
           DO 80 I = 1, N
            JPERM = id%UNS_PERM(I)
            RHS_MUMPS( KDEC+JPERM ) = RW1( I )
 80        CONTINUE
          END DO
          DEALLOCATE( RW1 )
        END IF
        IF (PROK) WRITE( MP, 120 ) INFO(1), INFO(2)
      END IF
      IF (id%MYID.EQ.MASTER .and.ICNTL(21)==0) THEN
        IF ( INFO(1) .GE. 0 .AND. ICNTL(4).GE.3 .AND. ICNTL(3).GT.0)
     *    THEN
          K = min0(10, N)
          IF (ICNTL(4) .eq. 4 ) K = N
          J = min0(10,NBRHS_EFF)
          IF (ICNTL(4) .eq. 4 ) J = NBRHS_EFF
          DO II=1, J
            WRITE(ICNTL(3),110) BEG_RHS+II-1
            WRITE(ICNTL(3),160)
     *    (RHS_MUMPS(IBEG+(II-1)*LRHS+I-1),I=1,K)
          ENDDO
        END IF
      END IF
      ENDDO
 90   CONTINUE
      IF (KEEP(201).NE.0)THEN
        IF (I_AM_SLAVE) THEN
          CALL DMUMPS_582(IERR)
          IF (IERR.LT.0 .AND. INFO(1) .GE. 0) INFO(1) = IERR
        ENDIF
        CALL DMUMPS_276( ICNTL, INFO,
     *         id%COMM,id%MYID)
      ENDIF
      IF (associated(id%BUFR)) THEN
          DEALLOCATE(id%BUFR)
          NULLIFY(id%BUFR)
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        IF (allocated(IWK_SOLVE)) DEALLOCATE( IWK_SOLVE )
        CALL DMUMPS_57( IERR )
        CALL DMUMPS_59( IERR )
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
        IF (allocated(SAVERHS)) DEALLOCATE( SAVERHS)
        IF (ICNTL20 .ne. 0 .and. ICNTL21.ne.0) THEN
          DEALLOCATE(RHS_MUMPS)
          NULLIFY(RHS_MUMPS)
        ENDIF
      ELSE
        IF (ASSOCIATED(RHS_MUMPS)) THEN
          DEALLOCATE(RHS_MUMPS)
          NULLIFY(RHS_MUMPS)
        END IF
      END IF
      IF (I_AM_SLAVE) THEN
        IF (allocated(SRW3)) DEALLOCATE(SRW3)
        IF (LSCAL .AND. ICNTL21==1) THEN
          DEALLOCATE(scaling_data%SCALING_LOC)
        ENDIF
      ENDIF
      IF (MAXSUPDATED)
     &  id%MAXS= MAXSSAVED ! restore original value of MAXS
      DEALLOCATE (RHSSOL)
      RETURN
 65   FORMAT (//' ERROR ANALYSIS AFTER ITERATIVE REFINEMENT')
 100  FORMAT(//' ****** SOLVE & CHECK STEP ********'/)
 110  FORMAT (//' VECTOR SOLUTION FOR COLUMN ',I12)
 115  FORMAT(1X, A44,1P,D9.2)
 120  FORMAT(//' LEAVING SOLVER WITH:  INFO(1) ............ =',I12/
     *         '                       INFO(2) ............ =',I12)
 150  FORMAT (/' STATISTICS PRIOR SOLVE PHASE     ...........'/
     *        ' NUMBER OF RIGHT-HAND-SIDES                    =',I12/
     *        ' BLOCKING FACTOR FOR MULTIPLE RHS              =',I12/
     *        ' ICNTL (9)                                     =',I12/
     *        '  --- (10)                                     =',I12/
     *        '  --- (11)                                     =',I12/
     *        '  --- (20)                                     =',I12/
     *        '  --- (21)                                     =',I12)
 160  FORMAT (' RHS'/(1X,1P,5D14.6))
 170  FORMAT (//' ERROR ANALYSIS' )
 240  FORMAT (1X, A42,I4)
 270  FORMAT (//' BEGIN ITERATIVE REFINEMENT' )
  81  FORMAT (/' STATISTICS AFTER ITERATIVE REFINEMENT ')
 131  FORMAT (/' END   ITERATIVE REFINEMENT ')
 141  FORMAT(1X, A42,I4)
      END SUBROUTINE DMUMPS_301
      SUBROUTINE DMUMPS_245(root, N, A, LA, IW, LIW, W, LWC, 
     * IWCB,LIWW,RHS,LRHS,NRHS,RHSSOL,NA,LNA,NE_STEPS, W2,
     * MTYPE, ICNTL,
     * STEP, FRERE, DAD, FILS, PTRIST, PTRFAC, IW1,LIW1,
     * PROCNODE_STEPS, SLAVEF,
     * INFO, KEEP,KEEP8, COMM, COMM_NODES, MYID,
     * MYID_NODES,
     * BUFR, LBUFR, LBUFR_BYTES,
     * 
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE
     * )
      USE DMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
#if defined(V_T)
      INCLUDE 'VT.inc'
#endif
      TYPE ( DMUMPS_ROOT_STRUC ) :: root
      INTEGER LA,LWC,N,LIW,MTYPE,LIW1,LIWW,LNA
      INTEGER ICNTL(40),INFO(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IW(LIW),IW1(LIW1),NA(LNA),NE_STEPS(KEEP(28)),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)), FILS(N), PTRIST(KEEP(28)),
     * PTRFAC(KEEP(28)), DAD(KEEP(28))
      INTEGER LRHS, NRHS
      DOUBLE PRECISION    A(LA), W(LWC), RHS(LRHS,NRHS),
     *        RHSSOL(N,NRHS), W2(KEEP(133))
      INTEGER SLAVEF, COMM, COMM_NODES, MYID, MYID_NODES
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER MP, LP, LDIAG, LWC_LOC
      INTEGER K,I
      INTEGER LPOOL,MYLEAF
      INTEGER NSTK_S,IPOOL,PTRICB,PTRACB
      INTEGER POOLSS, MTYPE_LOC
      INTEGER SIZE_ROOT, IPT_RHS_ROOT
      INTEGER IERR
      INTEGER       IOLDPS, IAPOS,
     *              LOCAL_M     ,
     *              LOCAL_N     ,
     *              MASTER_ROOT
#if defined(V_T)
      INTEGER soln_c_class, forw_soln, back_soln, root_soln
#endif
           INTEGER INODE, IPOS, LIELL, NPIV,J1,JJ
      INTEGER IZERO
      DOUBLE PRECISION ZERO
      PARAMETER (IZERO = 0 )
      PARAMETER (ZERO = 0.0D0)
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      EXTERNAL DMUMPS_248, DMUMPS_249
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275
      LP      = ICNTL(1)
      MP      = ICNTL(2)
      LDIAG   = ICNTL(4)
#if defined(V_T)
      CALL VTCLASSDEF( 'Soln_c',soln_c_class,ierr)
      CALL VTFUNCDEF( 'forw_soln',soln_c_class,forw_soln,ierr)
      CALL VTFUNCDEF( 'back_soln',soln_c_class,back_soln,ierr)
      CALL VTFUNCDEF( 'root_soln',soln_c_class,root_soln,ierr)
#endif
      NSTK_S   = 1
      POOLSS = NSTK_S + KEEP(28)
      PTRICB = POOLSS + KEEP(28) + 1
      PTRACB = PTRICB + KEEP(28)
      IPOOL  = PTRACB + KEEP(28)
      LPOOL  = LIW1 - IPOOL + 1
      IF ( LPOOL .ne. KEEP(28) + 1 ) THEN
        WRITE(MP,*) ' warning solve: lpool.NE.n+1 :',lpool,n+1
      ENDIF
      IF (KEEP(201).NE.0) THEN
      CALL DMUMPS_584(PTRFAC,KEEP(28))
      CALL DMUMPS_585(A,LA,PTRFAC,KEEP(28),IERR)
      IF(IERR.LT.0)THEN
         INFO(1)=IERR
         INFO(2)=0
         RETURN
      ENDIF
      ENDIF
      LWC_LOC = LWC
      SIZE_ROOT  = 0
      IPT_RHS_ROOT  = LWC
      IF ( KEEP( 38 ) .NE. 0 .OR. KEEP( 20 ) .NE. 0 ) THEN
          IF ( KEEP( 38 ) .ne. 0 ) THEN
            MASTER_ROOT = DMUMPS_275(STEP( KEEP(38)),
     *                    PROCNODE_STEPS, SLAVEF )
          ELSE
            MASTER_ROOT = DMUMPS_275(STEP(KEEP(20)),
     *                    PROCNODE_STEPS, SLAVEF )
          END IF
          IF ( MASTER_ROOT .eq. MYID_NODES ) THEN
            IF ( KEEP( 38 ) .NE. 0 ) THEN
              SIZE_ROOT    = root%TOT_ROOT_SIZE
            ELSE
              SIZE_ROOT    = IW(PTRIST(STEP(KEEP(20)))+XSIZE)
            END IF
            IPT_RHS_ROOT = LWC -NRHS * SIZE_ROOT + 1
            LWC_LOC = LWC_LOC - NRHS * SIZE_ROOT
          END IF
      END IF
      IF ( KEEP( 50 ) .eq. 0 ) THEN
        MTYPE_LOC = MTYPE
      ELSE
        MTYPE_LOC = 1
      ENDIF
#if defined(V_T)
      CALL VTBEGIN(forw_soln,ierr)
#endif
      CALL DMUMPS_248(N, A(1), LA, IW(1), LIW, W(1),
     *           LWC_LOC, RHS, LRHS, NRHS,
     *           IW1(PTRICB), IWCB, LIWW,
     *           RHSSOL,NE_STEPS, NA, LNA, STEP, FRERE,DAD,FILS,
     *           IW1(NSTK_S),IW1(IPOOL),LPOOL,PTRIST,PTRFAC,
     *           MYLEAF,INFO,
     *           KEEP,KEEP8,
     *           PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           W( IPT_RHS_ROOT ), MTYPE_LOC, 
     * 
     *           ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *           )
#if defined(V_T)
      CALL VTEND(forw_soln,ierr)
#endif
      IF ( INFO(1) .LT. 0 ) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP,*) MYID,
     *    ': ** ERROR RETURN FROM DMUMPS_248,INFO(1:2)=',
     *    INFO(1:2)
        END IF
        RETURN
      END IF
      IF ( KEEP( 38 ) .NE. 0 ) THEN
        IF ( KEEP(60) == 0 ) THEN
        IOLDPS = PTRIST(STEP(KEEP(38)))
        LOCAL_M = IW( IOLDPS + 2 + XSIZE)
        LOCAL_N = IW( IOLDPS + 1 + XSIZE)
       IF(KEEP(201).NE.0) THEN
        IF(.NOT.DMUMPS_SOLVE_IS_INODE_IN_MEM(KEEP(38),
     $       PTRFAC,KEEP(28),A,LA,IERR))THEN
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              RETURN
           ENDIF
           CALL DMUMPS_578(KEEP(38),PTRFAC,
     $          KEEP,KEEP8,A,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              RETURN
           ENDIF
           CALL DMUMPS_577(
     &           A(PTRFAC(IW(IOLDPS+4+ XSIZE))),KEEP(38),IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              RETURN
           ENDIF
        ELSE
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              RETURN
           ENDIF
        ENDIF
       ENDIF
        IAPOS   = PTRFAC(IW( IOLDPS + 4 + XSIZE))
          IF ( .NOT. ASSOCIATED( root%QR_TAU ) )
     *          ALLOCATE( root%QR_TAU(1))
          IF ( .NOT. ASSOCIATED( root%GROW ) )
     *          ALLOCATE( root%GROW( 1 ) )
          IF ( .NOT. ASSOCIATED( root%GSIN ) )
     *          ALLOCATE( root%GSIN( 1 ) )
          IF ( .NOT. ASSOCIATED( root%GCOS ) )
     *          ALLOCATE( root%GCOS( 1 ) )
          CALL MPI_BARRIER(COMM_NODES, IERR)
#if defined(V_T)
      CALL VTBEGIN(root_soln,ierr)
#endif
        IF ( root%yes ) THEN
          CALL DESCINIT( root%DESCB, root%TOT_ROOT_SIZE,
     *      NRHS, root%MBLOCK, root%NBLOCK, 0, 0,
     *      root%CNTXT_BLACS, LOCAL_M, IERR )
          IF (IERR.NE.0) THEN
            WRITE(*,*) 'After DESCINIT, IERR = ', IERR
            CALL DMUMPS_ABORT()
          END IF
#if defined(try_null_space_not_tested)
          CALL DMUMPS_352( NRHS, root%DESCRIPTOR,
     *       root%DESCB,
     *       root%CNTXT_BLACS, LOCAL_M, LOCAL_N,
     *       root%MBLOCK, root%NBLOCK,
     *       root%IPIV, root%LPIV, MASTER_ROOT, MYID_NODES,
     *       COMM_NODES,
     *       W( IPT_RHS_ROOT ),
     *       root%TOT_ROOT_SIZE, A( IAPOS ),
     *       INFO(1), MTYPE, KEEP(50), KEEP(19),
     *       root%QR_TAU, W(1), LWC_LOC, KEEP(17),
     *       root%MAXG, root%GIND, root%GROW, root%GCOS, root%GSIN )
#else
          CALL DMUMPS_286( NRHS, root%DESCRIPTOR, root%DESCB,
     *       root%CNTXT_BLACS, LOCAL_M, LOCAL_N,
     *       root%MBLOCK, root%NBLOCK,
     *       root%IPIV, root%LPIV, MASTER_ROOT, MYID_NODES,
     *       COMM_NODES,
     *       W( IPT_RHS_ROOT ),
     *       root%TOT_ROOT_SIZE, A( IAPOS ),
     *       INFO(1), MTYPE, KEEP(50))
#endif
         ENDIF
      ELSE
        IF ( MYID_NODES .eq.  DMUMPS_275( STEP(KEEP(38)),
     *        PROCNODE_STEPS, SLAVEF ) ) THEN
           W( IPT_RHS_ROOT:IPT_RHS_ROOT+NRHS*SIZE_ROOT- 1) = 0.0D0
        ENDIF
        CALL MPI_BARRIER( COMM_NODES, IERR )
      ENDIF
      ELSE IF ( KEEP(20) .NE. 0 ) THEN
        IF ( MYID_NODES .eq.  DMUMPS_275( STEP(KEEP(20)),
     *        PROCNODE_STEPS, SLAVEF ) ) THEN
#if defined(try_null_space)
         IF ( KEEP(60) .eq. 0 ) THEN
           CALL DMUMPS_296(
     *     SIZE_ROOT,A( PTRFAC(
              IW( PTRIST(STEP(KEEP(20)))+4+XSIZE))),
     *     ROOT,
     *     W( IPT_RHS_ROOT ), W(1),LWC_LOC,KEEP,KEEP8,MTYPE,INFO)
         ELSE
#endif
           W( IPT_RHS_ROOT:IPT_RHS_ROOT+NRHS*SIZE_ROOT- 1) = 0.0D0
#if defined(try_null_space)
         END IF
#endif
        END IF
        CALL MPI_BARRIER( COMM_NODES, IERR )
      ELSE
        CALL MPI_BARRIER( COMM_NODES, IERR )
      END IF
#if defined(V_T)
      CALL VTEND(root_soln,ierr)
#endif
      IF ( INFO(1) .LT. 0 ) RETURN
      IF(KEEP(201).NE.0)THEN
        IF(KEEP(38).NE.0)THEN
           IF(KEEP(60).EQ.0)THEN
              CALL DMUMPS_598(KEEP(38),
     $             PTRFAC,KEEP(28),A,LA,IERR)
              IF(IERR.LT.0)THEN
                 INFO(1)=IERR
                 INFO(2)=0
                 RETURN
              ENDIF
           ENDIF
        ENDIF
        CALL DMUMPS_583(PTRFAC,KEEP(28))
        IF(KEEP(38).NE.0)THEN
           CALL DMUMPS_598(KEEP(38),
     $          PTRFAC,KEEP(28),A,LA,IERR)
        ENDIF 
      ENDIF
      IF ( KEEP( 50 ) .eq. 0 ) THEN
        MTYPE_LOC = MTYPE
      ELSE
        MTYPE_LOC = IZERO
      ENDIF
#if defined(V_T)
      CALL VTBEGIN(back_soln,ierr)
#endif
      CALL DMUMPS_249( N, A, LA, IW, LIW, W(1), LWC_LOC,
     *              RHS, LRHS, NRHS,
     *              RHSSOL, IW1(PTRICB),IW1(PTRACB),IWCB,LIWW,
     *              IW1(POOLSS),W2,
     *              NE_STEPS, NA, LNA, STEP, FRERE,FILS,
     *              IW1(IPOOL),LPOOL,PTRIST,PTRFAC,MYLEAF,INFO,
     *              PROCNODE_STEPS, SLAVEF, COMM_NODES,MYID_NODES,
     *              BUFR, LBUFR, LBUFR_BYTES, KEEP,KEEP8,
     *              W( IPT_RHS_ROOT ),
     *              MTYPE_LOC, 
     *              ISTEP_TO_INIV2, TAB_POS_IN_PERE )
#if defined(V_T)
      CALL VTEND(back_soln,ierr)
#endif
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
        K = MIN0(10,N)
        IF (LDIAG.EQ.4) K = N
        WRITE (MP,99992)
        IF (N.GT.0) WRITE (MP,99993) (RHSSOL(I,1),I=1,K)
        IF (N.GT.0.and.NRHS>1) 
     *              WRITE (MP,99994) (RHSSOL(I,2),I=1,K)
      ENDIF
      RETURN
99993 FORMAT (' RHSSOL (first column)'/(1X,1P,5D14.6))
99994 FORMAT (' RHSSOL (2 nd  column)'/(1X,1P,5D14.6))
99992 FORMAT (//' LEAVING SOLVE (MPI41C) WITH')
      END SUBROUTINE DMUMPS_245
      SUBROUTINE DMUMPS_521(NSLAVES, N, MYID, COMM,
     *           MTYPE, RHSSOL, RHS, LRHS, NRHS, PTRIST,
     *           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, BUFFER,
     *           SIZE_BUF, SIZE_BUF_BYTES )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM, LIW, MTYPE
      INTEGER NRHS, LRHS
      DOUBLE PRECISION RHSSOL(   N, NRHS)
      DOUBLE PRECISION RHS   (LRHS, NRHS)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N)
      INTEGER SIZE_BUF, SIZE_BUF_BYTES
      INTEGER BUFFER(SIZE_BUF_BYTES)
      INTEGER I, J, JJ, J1, J2, INODE, MASTER,
     *        MYID_NODES, TYPE_PARAL, N2RECV
      INTEGER LIELL, IPOS, NPIV
      INTEGER MSGSOU, STATUS(MPI_STATUS_SIZE), IERR
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
      INTEGER RECORD_SIZE_P_1, SIZE1, SIZE2
      INTEGER POS_BUF, N2SEND, IROW
      DOUBLE PRECISION ONE_ROW(NRHS)
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275
      N2SEND=0
      N2RECV=N
      POS_BUF=0
      TYPE_PARAL = KEEP(46)
      IF (NSLAVES > 1 .OR. TYPE_PARAL == 0) THEN
        CALL MPI_PACK_SIZE(2,MPI_INTEGER, COMM, SIZE1, IERR)
        CALL MPI_PACK_SIZE(NRHS,MPI_DOUBLE_PRECISION, COMM,
     *                   SIZE2, IERR)
        RECORD_SIZE_P_1= SIZE1+SIZE2
      ELSE
        RECORD_SIZE_P_1 = -9999
      ENDIF
      I_AM_SLAVE = MYID .ne. 0 .OR. TYPE_PARAL .eq. 1
      IF ( TYPE_PARAL == 1 ) THEN
        MYID_NODES = MYID
      ELSE
        MYID_NODES = MYID-1
      ENDIF
      IF (I_AM_SLAVE) THEN
        POS_BUF = 0
        DO INODE = 1, N
          IF ( STEP(INODE) > 0 ) THEN
            IF (MYID_NODES == DMUMPS_275(STEP(INODE),
     *          PROCNODE_STEPS,NSLAVES)) THEN
              IF ( KEEP(38) .eq. INODE .OR.
     *             KEEP(20) .eq. INODE) THEN
                    IPOS = PTRIST(STEP(INODE)) 
                    LIELL = IW(IPOS+3+XSIZE)
                    NPIV = LIELL
                    IPOS= PTRIST(STEP(INODE))+5+XSIZE
              ELSE
                  IPOS = PTRIST(STEP(INODE)) + 2+ XSIZE
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(STEP(INODE)) + 5 +XSIZE)
              END IF
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              IF (MYID .EQ. MASTER) THEN
               DO JJ=J1,J1+NPIV-1
                RHS(IW(JJ),1:NRHS)=RHSSOL(IW(JJ),1:NRHS)
               ENDDO
               N2RECV=N2RECV-NPIV
              ELSE
               DO JJ=J1,J1+NPIV-1
                CALL DMUMPS_522( IW(JJ), RHSSOL(IW(JJ),1:NRHS) )
               ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        CALL DMUMPS_523()   ! Send remaining rows
      ENDIF
      IF ( MYID .EQ. MASTER ) THEN
       DO WHILE (N2RECV .NE. 0)
        CALL MPI_RECV( BUFFER, SIZE_BUF_BYTES, MPI_PACKED,
     *                 MPI_ANY_SOURCE,
     *                 GatherSol, COMM, STATUS, IERR )
        POS_BUF = 0
        CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     *                   IROW, 1, MPI_INTEGER, COMM, IERR)
        DO WHILE (IROW.NE.0)
          CALL MPI_UNPACK(BUFFER, SIZE_BUF_BYTES, POS_BUF,
     *                   ONE_ROW, NRHS, MPI_DOUBLE_PRECISION,
     *                   COMM, IERR)
          RHS(IROW,1:NRHS)=ONE_ROW
          CALL MPI_UNPACK( BUFFER, SIZE_BUF_BYTES, POS_BUF,
     *                   IROW, 1, MPI_INTEGER, COMM, IERR)
          N2RECV=N2RECV-1
        ENDDO
       ENDDO
      ENDIF
      RETURN
      CONTAINS
        SUBROUTINE DMUMPS_522(IROW, RHS_VAL)
        INTEGER IROW
        DOUBLE PRECISION RHS_VAL(NRHS)
        CALL MPI_PACK(IROW, 1, MPI_INTEGER, BUFFER,
     *                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(RHS_VAL, NRHS, MPI_DOUBLE_PRECISION,
     *                BUFFER, SIZE_BUF_BYTES, POS_BUF, COMM,
     *                IERR)
        N2SEND=N2SEND+1
        IF ( POS_BUF + RECORD_SIZE_P_1 > SIZE_BUF_BYTES ) THEN
          CALL DMUMPS_523()
        END IF
        RETURN
        END SUBROUTINE DMUMPS_522
        SUBROUTINE DMUMPS_523()
        IF (N2SEND .NE. 0) THEN
         CALL MPI_PACK(0, 1, MPI_INTEGER, BUFFER,
     *                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
         CALL MPI_SEND(BUFFER, POS_BUF, MPI_PACKED, MASTER, 
     *                 GatherSol, COMM, IERR)
        ENDIF
        POS_BUF=0
        N2SEND=0
        RETURN
        END SUBROUTINE DMUMPS_523
      END SUBROUTINE DMUMPS_521
      SUBROUTINE DMUMPS_535(MTYPE, ISOL_LOC,
     *             LSOL_LOC, PTRIST, KEEP,KEEP8,
     *             IW, LIW_PASSED, MYID_NODES, N, STEP,
     *             PROCNODE, NSLAVES, scaling_data, LSCAL)
      IMPLICIT NONE
      INTEGER MTYPE, LSOL_LOC, MYID_NODES, N, NSLAVES
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE(KEEP(28))
      INTEGER ISOL_LOC(LSOL_LOC)
      INTEGER LIW_PASSED
      INTEGER IW(LIW_PASSED)
      INTEGER STEP(N)
      LOGICAL LSCAL
      type scaling_data_t
        SEQUENCE
        DOUBLE PRECISION, dimension(:), pointer :: SCALING
        DOUBLE PRECISION, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      type (scaling_data_t) :: scaling_data
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275
      INTEGER ISTEP, K
      INTEGER J1, IPOS, LIELL, NPIV, JJ
      LOGICAL ROOT
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      K=0
      DO ISTEP=1, KEEP(28)
          IF ( MYID_NODES == DMUMPS_275( ISTEP,
     *         PROCNODE, NSLAVES)) THEN
              ROOT=.false.
              IF (KEEP(38).ne.0) ROOT = STEP(KEEP(38))==ISTEP
              IF (KEEP(20).ne.0) ROOT = STEP(KEEP(20))==ISTEP
              IF ( ROOT ) THEN
                    IPOS = PTRIST(ISTEP)+XSIZE
                    LIELL = IW(IPOS+3)
                    NPIV = LIELL
                    IPOS= PTRIST(ISTEP)+5+XSIZE
              ELSE
                  IPOS = PTRIST(ISTEP) + 2 + XSIZE
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 + XSIZE)
              END IF
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              DO JJ=J1,J1+NPIV-1
                  K=K+1
                  ISOL_LOC(K)=IW(JJ)
                  IF (LSCAL) THEN
                    scaling_data%SCALING_LOC(K)=
     *              scaling_data%SCALING(IW(JJ))
                  ENDIF
              ENDDO
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_535
      SUBROUTINE DMUMPS_532(
     *           SLAVEF, N, MYID_NODES,
     *           MTYPE, RHSSOL, NRHS,
     *           ISOL_LOC, SOL_LOC, BEG_RHS, LSOL_LOC,
     *           PTRIST,
     *           PROCNODE_STEPS, KEEP,KEEP8, IW, LIW, STEP,
     *           scaling_data, LSCAL)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      type scaling_data_t
        SEQUENCE
        DOUBLE PRECISION, dimension(:), pointer :: SCALING
        DOUBLE PRECISION, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      TYPE (scaling_data_t) :: scaling_data
      LOGICAL LSCAL
      INTEGER SLAVEF, N, MYID_NODES, LIW, MTYPE, NRHS
      INTEGER LSOL_LOC, BEG_RHS
      INTEGER ISOL_LOC(LSOL_LOC)
      DOUBLE PRECISION SOL_LOC( LSOL_LOC, *)
      DOUBLE PRECISION RHSSOL(   N, NRHS)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N)
      INTEGER JJ, J1, J2, ISTEP, K
      INTEGER IPOS, LIELL, NPIV
      LOGICAL ROOT
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275
      K=0
        DO ISTEP = 1, KEEP(28)
            IF (MYID_NODES == DMUMPS_275(ISTEP,
     *          PROCNODE_STEPS,SLAVEF)) THEN
              ROOT=.false.
              IF (KEEP(38).ne.0) ROOT = STEP(KEEP(38))==ISTEP
              IF (KEEP(20).ne.0) ROOT = STEP(KEEP(20))==ISTEP
              IF ( ROOT ) THEN
                    IPOS = PTRIST(ISTEP) + XSIZE
                    LIELL = IW(IPOS+3)
                    NPIV = LIELL
                    IPOS= PTRIST(ISTEP)+5+XSIZE
              ELSE
                  IPOS = PTRIST(ISTEP) + 2 +XSIZE
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +XSIZE)
              END IF
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              DO JJ=J1,J1+NPIV-1
                K=K+1
                IF (LSCAL) THEN
                  SOL_LOC(K,BEG_RHS:BEG_RHS+NRHS-1) =
     *            scaling_data%SCALING_LOC(K)*RHSSOL(IW(JJ),1:NRHS)
                ELSE
                  SOL_LOC(K,BEG_RHS:BEG_RHS+NRHS-1) =
     *            RHSSOL(IW(JJ),1:NRHS)
                ENDIF
              ENDDO
            ENDIF
        ENDDO
      RETURN
      END SUBROUTINE DMUMPS_532
      SUBROUTINE DMUMPS_248(N, A, LA, IW, LIW, WCB, LWCB,
     *    RHS, LRHS, NRHS,
     *    PTRICB, IWCB, LIWCB, W1, NE_STEPS, NA, LNA, STEP,
     *    FRERE, DAD, FILS,
     *    NSTK_S, IPOOL, LPOOL, PTRIST, PTRFAC, MYLEAF, INFO,
     *    KEEP,KEEP8,
     *    PROCNODE_STEPS,
     *    SLAVEF, COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,
     *    RHS_ROOT, MTYPE, 
     *
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      USE DMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER N, LA, LIW, LWCB, LPOOL, LIWCB, LNA
      INTEGER SLAVEF, MYLEAF, COMM, MYID
      INTEGER INFO( 40 ), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER LRHS, NRHS
      DOUBLE PRECISION A( LA ), RHS( LRHS, NRHS ), WCB( LWCB )
      DOUBLE PRECISION W1( N, NRHS ), RHS_ROOT( * )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER NA( LNA ), NE_STEPS( KEEP(28) )
      INTEGER STEP( N ), FRERE( KEEP(28) ), FILS( N ),
     *        DAD( KEEP(28) )
      INTEGER NSTK_S(KEEP(28)), IPOOL( LPOOL )
      INTEGER PTRIST(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER PTRICB( KEEP(28) )
      INTEGER IW( LIW ), IWCB( LIWCB )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MSGTAG, MSGSOU, DUMMY(1)
      LOGICAL FLAG
      INTEGER NBFIN, MYROOT
      INTEGER POSIWCB,POSWCB,PLEFTWCB
      INTEGER INODE
      INTEGER I, K
      INTEGER III, NBROOT,NBLEAF,LEAF
      LOGICAL BLOQ
      EXTERNAL DMUMPS_275
      INTEGER DMUMPS_275
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      POSIWCB = LIWCB
      POSWCB  = LWCB
      PLEFTWCB= 1
      DO I = 1, KEEP(28)
        NSTK_S(I)   = NE_STEPS(I)
      ENDDO
      PTRICB = 0
      DO K=1, NRHS
      DO 10 I = 1, N
        W1(I,K)     = DBLE(ZERO)
   10 CONTINUE
      ENDDO
      CALL DMUMPS_362(N, LEAF, NBROOT, MYROOT, MYID,
     *     SLAVEF, NA, LNA, KEEP,KEEP8, STEP,
     *     PROCNODE_STEPS, IPOOL, LPOOL)
      NBFIN = SLAVEF
      IF ( MYROOT .EQ. 0 ) THEN
        NBFIN = NBFIN - 1
        DUMMY(1) = 1
        CALL DMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *       RACINE_SOLVE, SLAVEF)
      END IF
      MYLEAF = LEAF - 1
      III    = 1
   50 CONTINUE
      IF (SLAVEF .EQ. 1) THEN
         CALL DMUMPS_574
     &        ( IPOOL(1), LPOOL, III, LEAF, INODE,
     &          KEEP(208), MYID )
        GOTO 60
      ENDIF
      BLOQ = ( III .EQ. LEAF )
      CALL DMUMPS_303( BLOQ, FLAG,
     *     BUFR, LBUFR, LBUFR_BYTES,
     *     MYID, SLAVEF, COMM,
     *     N, NRHS, IPOOL, LPOOL, III, LEAF,
     *     NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     *     IWCB, LIWCB,
     *     WCB, LWCB, POSWCB,
     *     PLEFTWCB, POSIWCB,
     *     PTRICB, INFO, KEEP,KEEP8, STEP,
     *     PROCNODE_STEPS,W1)
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      IF (.not. FLAG) THEN
         IF (III .NE. LEAF) THEN
            CALL DMUMPS_574
     &           (IPOOL(1), LPOOL, III, LEAF, INODE,
     &            KEEP(208), MYID )
            GOTO 60
         ENDIF
      ENDIF
      GOTO 50
 60   CONTINUE
      CALL DMUMPS_302( INODE, BUFR, LBUFR, LBUFR_BYTES,
     *        MSGTAG, MSGSOU, MYID, SLAVEF, COMM,  N,
     *        IPOOL, LPOOL, III, LEAF, NBFIN, NSTK_S,
     *        IWCB, LIWCB, WCB, LWCB, A, LA,
     *        IW, LIW, RHS, LRHS, NRHS, W1,
     *        POSWCB, PLEFTWCB, POSIWCB,
     *        PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     *        FILS, STEP, FRERE, DAD,
     *        MYROOT, INFO, KEEP,KEEP8, RHS_ROOT, MTYPE, 
     *
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *      )
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      GOTO 50
  260 CONTINUE
      CALL DMUMPS_150( MYID,COMM,BUFR,
     *                            LBUFR,LBUFR_BYTES )
      RETURN
      END SUBROUTINE DMUMPS_248
      RECURSIVE SUBROUTINE DMUMPS_323
     *         ( BUFR, LBUFR, LBUFR_BYTES,
     *           MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     *           N, NRHS, IPOOL, LPOOL, III, LEAF,
     *           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     *           PTRFAC, IWCB, LIWCB,
     *           WCB, LWCB, POSWCB,
     *           PLEFTWCB, POSIWCB,
     *           PTRICB,
     *           INFO, KEEP,KEEP8, STEP, PROCNODE_STEPS, W1 )
      USE DMUMPS_OOC 
      USE DMUMPS_BUFFER 
      IMPLICIT NONE
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIW, LA
      INTEGER N, NRHS, LPOOL, III, LEAF, NBFIN
      INTEGER LIWCB, LWCB, POSWCB, PLEFTWCB, POSIWCB
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER*8 KEEP8(150)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ),  NSTK_S( N )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      INTEGER PTRICB(KEEP(28)),PTRIST(KEEP(28)),PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      DOUBLE PRECISION WCB( LWCB ), A( LA ), W1( N, NRHS )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, K, JJ
      INTEGER FINODE, FPERE, LONG, NCB, POSITION, NCV, NPIV
      INTEGER PTRX, PTRY, APOS, PDEST, I
      LOGICAL FLAG
      EXTERNAL DMUMPS_275
      INTEGER  DMUMPS_275
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER( ONE = 1.0D0, ALPHA = -1.0D0 )
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      IF ( MSGTAG .EQ. RACINE_SOLVE ) THEN
          NBFIN = NBFIN - 1
          IF ( NBFIN .eq. 0 ) GOTO 270
      ELSE  IF (MSGTAG .EQ. ContVec ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FINODE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FPERE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         NCB, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         LONG, 1, MPI_INTEGER, COMM, IERR )
          IF ( NCB .eq. 0 ) THEN
             PTRICB(STEP(FINODE)) = -1
             NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
               IPOOL( LEAF ) = FPERE
               LEAF = LEAF + 1
               IF ( LEAF > LPOOL ) THEN
                 WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
                 CALL DMUMPS_ABORT()
               END IF
             END IF
          ELSE
          IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
             PTRICB(STEP(FINODE)) = NCB + 1
          END IF
          IF ( ( POSIWCB - LONG ) .LT. 0 ) THEN
            INFO( 1 ) = -14
            INFO( 2 ) = LONG
            GOTO 260
          END IF
          IF ( POSWCB - PLEFTWCB + 1 .LT. LONG * NRHS) THEN
            INFO( 1 ) = -11
            INFO( 2 ) = PLEFTWCB - POSWCB - 1 + LONG * NRHS
            GOTO 260
          END IF
          IF (LONG .GT. 0) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *          IWCB( 1 ),
     *          LONG, MPI_INTEGER, COMM, IERR )
           DO K = 1, NRHS
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *               WCB( PLEFTWCB ),
     *               LONG, MPI_DOUBLE_PRECISION, COMM, IERR )
            DO I = 1, LONG
              W1(IWCB(I),K) = W1(IWCB(I),K) + WCB(PLEFTWCB+I-1)
            ENDDO
           END DO
           PTRICB(STEP(FINODE)) = PTRICB(STEP(FINODE)) - LONG
          ENDIF
          IF ( PTRICB(STEP(FINODE)) == 1 ) THEN
               NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
          END IF
          IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
            IPOOL( LEAF ) = FPERE
            LEAF = LEAF + 1
            IF ( LEAF > LPOOL ) THEN
              WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
              CALL DMUMPS_ABORT()
            END IF
          ENDIF
        END IF
        ELSEIF ( MSGTAG .EQ. MASTER2SLAVE ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FINODE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FPERE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         NCV, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         NPIV, 1, MPI_INTEGER, COMM, IERR )
          PTRY = PLEFTWCB
          PTRX = PLEFTWCB + NCV * NRHS
          PLEFTWCB = PLEFTWCB + (NPIV + NCV) * NRHS
          IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
            INFO(1) = -11
            INFO(2) = -POSWCB + PLEFTWCB -1
            GO TO 260
          END IF
          DO K=1, NRHS
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         WCB( PTRY + (K-1) * NCV ), NCV,
     *         MPI_DOUBLE_PRECISION, COMM, IERR )
          ENDDO
          IF ( NPIV .GT. 0 ) THEN
            DO K=1, NRHS
              CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *          WCB( PTRX + (K-1)*NPIV ), NPIV,
     *          MPI_DOUBLE_PRECISION, COMM, IERR )
            END DO
          END IF
          IF (KEEP(201).NE.0) THEN
           IF(.NOT.DMUMPS_SOLVE_IS_INODE_IN_MEM(FINODE,PTRFAC,
     $         KEEP(28),A,LA,IERR))THEN
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
             CALL DMUMPS_578(FINODE,PTRFAC,
     $            KEEP,KEEP8,A,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
             CALL DMUMPS_577(
     $            A(PTRFAC(IW( PTRIST( STEP(FINODE)) + 3 +XSIZE))),
     $            FINODE,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
           ELSE
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
           ENDIF
          ENDIF
          APOS = PTRFAC(IW( PTRIST( STEP(FINODE)) + 3 + XSIZE ))
          IF ( NRHS == 1 ) THEN
            CALL DGEMV( 'T', NPIV, NCV, ALPHA, A(APOS), NPIV,
     *                WCB( PTRX ), 1, ONE,
     *                WCB( PTRY ), 1 )
          ELSE
            CALL DGEMM( 'T', 'N', NCV, NRHS, NPIV, ALPHA, A(APOS), NPIV,
     *                WCB( PTRX), NPIV, ONE,
     *                WCB( PTRY), NCV )
          ENDIF
          IF (KEEP(201).NE.0) THEN
             CALL DMUMPS_598(FINODE,PTRFAC,
     $            KEEP(28),A,LA,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
          ENDIF
          PLEFTWCB = PLEFTWCB - NPIV * NRHS
          PDEST = DMUMPS_275( STEP(FPERE),
     *            PROCNODE_STEPS, SLAVEF )
          IF ( PDEST .EQ. MYID ) THEN
            IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
              NCB = IW( PTRIST(STEP(FINODE)) + 2 + XSIZE )
              PTRICB(STEP(FINODE)) = NCB + 1
            END IF
            DO I = 1, NCV
              JJ=IW(PTRIST(STEP(FINODE))+3+I+ XSIZE )
              DO K=1, NRHS
                W1(JJ,K)= W1(JJ,K) + WCB(PTRY+I-1+(K-1)*NCV)
              ENDDO
            END DO
            PTRICB(STEP(FINODE)) =
     *      PTRICB(STEP(FINODE)) - NCV
            IF ( PTRICB( STEP( FINODE ) ) == 1 ) THEN
                 NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
            END IF
            IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
              IPOOL( LEAF ) = FPERE
              LEAF = LEAF + 1
              IF ( LEAF > LPOOL ) THEN
                WRITE(*,*) 'INTERNAL Error 41r: Pool is too small.'
                CALL DMUMPS_ABORT()
              END IF
            ENDIF
          ELSE
 210       CONTINUE
           CALL DMUMPS_78( NRHS, FINODE, FPERE,
     *       IW(PTRIST(STEP( FINODE )) + 2 + XSIZE ), NCV, NCV,
     *       IW(PTRIST(STEP(FINODE))+4+ XSIZE ),
     *       WCB( PTRY ), PDEST, ContVec, COMM, IERR )
            IF ( IERR .EQ. -1 ) THEN
              CALL DMUMPS_303( .FALSE., FLAG,
     *               BUFR, LBUFR, LBUFR_BYTES,
     *               MYID, SLAVEF, COMM,
     *               N, NRHS, IPOOL, LPOOL, III, LEAF,
     *               NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     *               IWCB, LIWCB,
     *               WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     *               PTRICB, INFO, KEEP,KEEP8, STEP,
     *               PROCNODE_STEPS, W1 )
              IF ( INFO( 1 )  .LT. 0 )  GOTO 270
              GOTO 210
            ELSE IF ( IERR .EQ. -2 ) THEN
              INFO( 1 ) = -17
              INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     *                      NCV * KEEP( 35 )
              GOTO 260
            END IF
          END IF
          PLEFTWCB = PLEFTWCB - NCV * NRHS
        ELSEIF ( MSGTAG .EQ. TERREUR ) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GOTO 270
        ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     *      (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
        ELSE
          INFO(1)=-100
          INFO(2)=MSGTAG
          GO TO 260
        ENDIF
        GO TO 270
 260    CONTINUE
        CALL DMUMPS_44( MYID, SLAVEF, COMM )
 270    CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_323
      SUBROUTINE DMUMPS_302( INODE,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     *           N, IPOOL, LPOOL, III, LEAF,
     *           NBFIN, NSTK_S,
     *           IWCB, LIWCB,
     *           WCB, LWCB, A, LA, IW, LIW,
     *           RHS, LRHS, NRHS, W1, POSWCB,
     *           PLEFTWCB, POSIWCB,
     *           PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     *           FILS, STEP, FRERE, DAD,
     *           MYROOT,
     *           INFO, KEEP,KEEP8, RHS_ROOT, MTYPE,
     *
     *           ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *           )
      USE DMUMPS_OOC
      USE DMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER INODE, LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIWCB, LWCB, LIW, LA, POSWCB, PLEFTWCB, POSIWCB
      INTEGER N, LPOOL, III, LEAF, NBFIN
      INTEGER MYROOT
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER*8 KEEP8(150)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ), NSTK_S(KEEP(28))
      INTEGER IWCB( LIWCB ), IW( LIW )
      INTEGER LRHS, NRHS
      DOUBLE PRECISION WCB( LWCB ), A( LA )
      DOUBLE PRECISION RHS(LRHS, NRHS ), W1( N, NRHS ),
     *                 RHS_ROOT( * )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER FILS( N ), STEP( N ), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      EXTERNAL DGEMV, DTRSV, DGEMM, DTRSM, DMUMPS_275
      INTEGER DMUMPS_275
      DOUBLE PRECISION ALPHA,ONE,ZERO
      PARAMETER(ONE=1.0D0, ALPHA=-1.0D0, ZERO=0.0D0)
      INTEGER I, J, K, IPOS, NSLAVES, J1, J2, J3, FPERE, NPIV, NCB,
     *        IERR,
     *        APOS, APOS1, IF, IFR, IPOSCB, APOSCB, LIELL, IN, JJ,
     *        NELIM, PLEFT, PCB_COURANT, PPIV_COURANT
      INTEGER Effective_CB_SIZE, NUPDATE, ISLAVE, PDEST, FirstIndex
      LOGICAL FLAG
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INTEGER APOS2,POSWCB1,POSWCB2,J22
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY( 1 )
      IF ( INODE .eq. KEEP(38 ) .OR. INODE .eq.KEEP( 20 ) ) THEN
        LIELL = IW( PTRIST( STEP(INODE)) + 3 + XSIZE)
        NPIV  = LIELL
        NELIM = 0
        NSLAVES = 0
        IPOS = PTRIST( STEP(INODE)) + 5 + XSIZE
      ELSE
        IPOS = PTRIST(STEP(INODE)) + 2 + XSIZE
        LIELL = IW(IPOS-2)+IW(IPOS+1)
        NELIM = IW(IPOS-1)
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + XSIZE )
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IPOS = IPOS + 1
        IF (KEEP(201).NE.0) THEN
         IF(.NOT.DMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,PTRFAC,
     $       KEEP(28),A,LA,IERR))THEN
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           CALL DMUMPS_578(INODE,PTRFAC,
     $          KEEP,KEEP8,A,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           CALL DMUMPS_577(A(PTRFAC(IW(IPOS))),INODE,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
         ELSE
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
         ENDIF
        ENDIF
        APOS = PTRFAC(IW(IPOS))
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + XSIZE)
        IPOS = IPOS + 1 + NSLAVES
      END IF
      IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
        J1 = IPOS + 1
        J2 = IPOS + LIELL
        J3 = IPOS + NPIV
      ELSE
        J1 = IPOS + LIELL + 1
        J2 = IPOS + 2 * LIELL
        J3 = IPOS + LIELL + NPIV
      END IF
      NCB = LIELL-NPIV
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq. KEEP(20) ) THEN
        IFR = 0
        DO JJ = J1, J3
          J = IW( JJ )
          IFR = IFR + 1
          DO K=1,NRHS
            RHS_ROOT(IFR+NPIV*(K-1)) = RHS(J,K) + W1(J,K)
          END DO
        END DO
        IF ( NPIV .LT. LIELL ) THEN
          WRITE(*,*) ' Internal error in SOLVE_NODE for Root node'
          CALL DMUMPS_ABORT()
        END IF
        MYROOT = MYROOT - 1
        IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
              DUMMY (1) = 1
              CALL DMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID,
     *             COMM, RACINE_SOLVE, SLAVEF)
            ENDIF
        END IF
        GO TO 270
      END IF
      PLEFT    = PLEFTWCB
      PPIV_COURANT = PLEFTWCB
      PCB_COURANT = PLEFT + NPIV*NRHS
      PLEFTWCB = PLEFTWCB + LIELL * NRHS
      IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
        INFO(1) = -11
        INFO(2) = PLEFTWCB - POSWCB - 1
        GO TO 260
      END IF
      IFR = PPIV_COURANT - 1
      DO 130 JJ = J1, J3
        J = IW(JJ)
        IFR = IFR + 1
        DO K=1, NRHS
          WCB(IFR+(K-1)*NPIV) = RHS(J,K) + W1(J,K)
        END DO
  130 CONTINUE
      IFR = PCB_COURANT - 1
      IF (NPIV .LT. LIELL) THEN
        DO 140 JJ = J3 + 1, J2
          J = IW(JJ)
          IFR = IFR + 1
          DO K=1, NRHS
            WCB(IFR+(K-1)*NCB) = W1(J,K)
            W1(J,K)=ZERO
          ENDDO
  140   CONTINUE
      ENDIF
      IF ( NPIV .NE. 0 ) THEN
        IF ( KEEP(50).NE.0) THEN
          IF ( NRHS == 1 ) THEN
            CALL DTRSV( 'U', 'T', 'U', NPIV, A(APOS), NPIV,
     *                   WCB(PPIV_COURANT), 1 )
          ELSE
            CALL DTRSM( 'L','U','T','U', NPIV, NRHS, ONE,
     *                   A(APOS), NPIV, WCB(PPIV_COURANT),
     *                   NPIV )
          ENDIF
        ELSE
          IF ( MTYPE .eq. 1 ) THEN
            IF ( NRHS == 1)  THEN
              CALL DTRSV( 'U', 'T', 'U', NPIV, A(APOS), LIELL, 
     *        WCB(PPIV_COURANT), 1 )
            ELSE
              CALL DTRSM( 'L','U','T','U', NPIV, NRHS, ONE,
     *                   A(APOS), LIELL, WCB(PPIV_COURANT),
     *                   NPIV )
            ENDIF
          ELSE
            IF (NRHS == 1) THEN
               CALL DTRSV( 'L', 'N', 'N', NPIV, A(APOS), LIELL,
     *         WCB(PPIV_COURANT), 1 )
            ELSE
              CALL DTRSM('L','L','N','N',NPIV, NRHS, ONE,
     *                   A(APOS), LIELL, WCB(PPIV_COURANT),
     *                   NPIV)
            ENDIF
          END IF
        END IF
      END IF
      NCB   = LIELL - NPIV
      IF ( MTYPE .EQ. 1 ) THEN
        IF ( KEEP(50) .eq. 0 ) THEN
          APOS1 = APOS  + NPIV * LIELL
        ELSE
          APOS1 = APOS + NPIV * NPIV
        END IF
        IF ( NSLAVES .EQ. 0 .OR. NPIV .eq. 0 ) THEN
          NUPDATE = NCB
        ELSE
          NUPDATE = NELIM
        END IF
      ELSE
        APOS1 = APOS + NPIV
        NUPDATE = NCB
      END IF
      IF ( NPIV .NE. 0 .AND. NUPDATE.NE.0 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          IF ( NRHS == 1 ) THEN
                CALL DGEMV('T', NPIV, NUPDATE, ALPHA, A(APOS1),
     *            NPIV,  WCB(PPIV_COURANT), 1, ONE,
     *            WCB(PCB_COURANT), 1)
          ELSE
             CALL DGEMM('T', 'N', NUPDATE, NRHS, NPIV, ALPHA,
     *            A(APOS1), NPIV, WCB(PPIV_COURANT), NPIV, ONE,
     *            WCB(PCB_COURANT), NCB)
          END IF
        ELSE
          IF ( NRHS == 1 ) THEN
                CALL DGEMV('N',NUPDATE, NPIV, ALPHA, A(APOS1),
     *            LIELL, WCB(PPIV_COURANT), 1,
     *            ONE, WCB(PCB_COURANT), 1 )
          ELSE
             CALL DGEMM('N', 'N', NUPDATE, NRHS, NPIV, ALPHA,
     *            A(APOS1), LIELL, WCB(PPIV_COURANT), NPIV, ONE,
     *            WCB(PCB_COURANT), NCB)
          END IF
        END IF
      END IF
      IFR = PPIV_COURANT - 1
      IF ( KEEP(50) .eq. 0 ) THEN
        DO 150 JJ = J1, J3
          J = IW(JJ)
          IFR = IFR + 1
          DO K=1, NRHS
            RHS(J,K) = WCB(IFR+(K-1)*NPIV)
          ENDDO
  150   CONTINUE
      ELSE
        APOS1 = APOS
        JJ = J1
        DO 
           IF(JJ .GT. J3) EXIT
           J = IW(JJ)
           IFR = IFR + 1
           IF(IW(JJ+LIELL) .GT. 0) THEN
              DO K=1, NRHS
                 RHS( J, K ) = WCB( IFR+(K-1)*NPIV ) * A( APOS1 )
              END DO
              APOS1 = APOS1 + NPIV + 1
              JJ = JJ+1
           ELSE
              J22 = IW(JJ+1)
              APOS2 = APOS1+NPIV+1
              DO K=1, NRHS
                 POSWCB1 = IFR+(K-1)*NPIV
                 POSWCB2 = POSWCB1+1
                 RHS(J,K) = WCB(POSWCB1)*A(APOS1)
     *                + WCB(POSWCB2)*A(APOS1+1)
                 RHS(J22,K) = WCB(POSWCB1)*A(APOS1+1)
     *                + WCB(POSWCB2)*A(APOS2)
              END DO
              APOS1 = APOS2 + NPIV + 1
              JJ = JJ+2
              IFR = IFR+1
           ENDIF
        ENDDO
      END IF
      IF (KEEP(201).NE.0) THEN
        CALL DMUMPS_598(INODE,PTRFAC,KEEP(28),
     $        A,LA,IERR)
        IF(IERR.LT.0)THEN
           INFO(1)=IERR
           INFO(2)=0
           GOTO 260
        ENDIF
      END IF
      FPERE = DAD(STEP(INODE))
      IF ( FPERE .EQ. 0 ) THEN
          MYROOT = MYROOT - 1
          PLEFTWCB = PLEFTWCB - LIELL *NRHS
          IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
              DUMMY (1) = 1
              CALL DMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID,
     *             COMM, RACINE_SOLVE, SLAVEF)
            ENDIF
          END IF
          GO TO 270
      ENDIF
      IF ( NUPDATE .NE. 0 .OR. NCB.eq.0 ) THEN
        IF (DMUMPS_275(STEP(FPERE),PROCNODE_STEPS,
     *      SLAVEF) .EQ. MYID) THEN
         IF ( NCB .ne. 0 ) THEN
          PTRICB(STEP(INODE)) = NCB + 1
          DO 190 I = 1, NUPDATE
            DO K=1, NRHS
             W1( IW(J3 + I), K ) = W1( IW(J3 + I), K )
     *       + WCB(PCB_COURANT + I-1 +(K-1)*NCB)
            ENDDO
  190     CONTINUE
          PTRICB(STEP( INODE )) = PTRICB(STEP( INODE )) - NUPDATE
          IF ( PTRICB(STEP(INODE)) == 1 ) THEN
            NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
            IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN
              IPOOL( LEAF ) = FPERE
              LEAF = LEAF + 1
            ENDIF
          END IF
         ELSE
          PTRICB(STEP( INODE )) = -1
          NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
          IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN
            IPOOL( LEAF ) = FPERE
            LEAF = LEAF + 1
          ENDIF
         ENDIF
        ELSE
 210      CONTINUE
          CALL DMUMPS_78( NRHS, INODE, FPERE, NCB, NCB,
     *                          NUPDATE,
     *                          IW( J3 + 1 ), WCB( PCB_COURANT ),
     *                          DMUMPS_275(STEP(FPERE),
     *                          PROCNODE_STEPS,SLAVEF),
     *                          ContVec,
     *                          COMM, IERR )
          IF ( IERR .EQ. -1 ) THEN
            CALL DMUMPS_303( .FALSE., FLAG,
     *             BUFR, LBUFR, LBUFR_BYTES,
     *             MYID, SLAVEF, COMM,
     *             N, NRHS, IPOOL, LPOOL, III, LEAF,
     *             NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     *             IWCB, LIWCB,
     *             WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     *             PTRICB, INFO, KEEP,KEEP8, STEP,
     *             PROCNODE_STEPS, W1 )
            IF ( INFO( 1 )  .LT. 0 )  GOTO 270
            GOTO 210
          ELSE IF ( IERR .EQ. -2 ) THEN
            INFO( 1 ) = -17
            INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     *                ( NUPDATE + 3 ) * KEEP( 34 )
            GOTO 260
          END IF
        ENDIF
      END IF
      IF ( NSLAVES .NE. 0 .AND. MTYPE .eq. 1
     *     .and. NPIV .NE. 0 ) THEN
        DO ISLAVE = 1, NSLAVES
          PDEST = IW( PTRIST(STEP(INODE)) + 5 + ISLAVE +XSIZE)
          CALL DMUMPS_49( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB - NELIM, 
     &                NSLAVES, 
     &                Effective_CB_Size, FirstIndex )
 222      CALL DMUMPS_72( NRHS,
     *    INODE, FPERE,
     *    Effective_CB_Size, NCB, NPIV,
     *    WCB( PCB_COURANT + NELIM + FirstIndex - 1 ),
     *    WCB( PPIV_COURANT ),
     *    PDEST, COMM, IERR )
          IF ( IERR .EQ. -1 ) THEN
            CALL DMUMPS_303( .FALSE., FLAG,
     *               BUFR, LBUFR, LBUFR_BYTES,
     *               MYID, SLAVEF, COMM,
     *               N, NRHS, IPOOL, LPOOL, III, LEAF,
     *               NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     *               IWCB, LIWCB,
     *               WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     *               PTRICB, INFO, KEEP,KEEP8, STEP,
     *               PROCNODE_STEPS, W1 )
            IF ( INFO( 1 )  .LT. 0 )  GOTO 270
            GOTO 222
          ELSE IF ( IERR .EQ. -2 ) THEN
            INFO( 1 ) = -17
            INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS*KEEP(35) +
     *      ( Effective_CB_Size + 4 ) * KEEP( 34 )
            GOTO 260
          END IF
        END DO
      END IF
      PLEFTWCB = PLEFTWCB - LIELL*NRHS
  270 CONTINUE
      RETURN
  260 CONTINUE
      CALL DMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE DMUMPS_302
      RECURSIVE SUBROUTINE DMUMPS_303( BLOQ, FLAG,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           MYID, SLAVEF, COMM,
     *           N, NRHS, IPOOL, LPOOL, III, LEAF,
     *           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     *           IWCB, LIWCB,
     *           WCB, LWCB, POSWCB,
     *           PLEFTWCB, POSIWCB,
     *           PTRICB, INFO, KEEP,KEEP8, STEP, PROCNODE_STEPS,
     *           W1 )
      LOGICAL BLOQ
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, NRHS, LPOOL, III, LEAF, NBFIN
      INTEGER LIWCB, LWCB, POSWCB, PLEFTWCB, POSIWCB
      INTEGER LIW, LA
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER*8 KEEP8(150)
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL)
      INTEGER NSTK_S( KEEP(28) )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      DOUBLE PRECISION WCB( LWCB ), A( LA ), W1( N, NRHS )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      LOGICAL FLAG
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      INTEGER MSGSOU, MSGTAG, MSGLEN
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
        CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     *                   COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     *                   FLAG, STATUS, IERR )
      END IF
      IF ( FLAG ) THEN
         MSGSOU = STATUS( MPI_SOURCE )
         MSGTAG = STATUS( MPI_TAG )
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL DMUMPS_44( MYID, SLAVEF, COMM )
         ELSE
           CALL MPI_RECV( BUFR, LBUFR_BYTES, MPI_PACKED,
     *                  MSGSOU, MSGTAG, COMM, STATUS, IERR )
           CALL DMUMPS_323( BUFR, LBUFR, LBUFR_BYTES,
     *          MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     *          N, NRHS, IPOOL, LPOOL, III, LEAF,
     *          NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     *          IWCB, LIWCB,
     *          WCB, LWCB, POSWCB,
     *          PLEFTWCB, POSIWCB,
     *          PTRICB, INFO, KEEP,KEEP8, STEP,
     *          PROCNODE_STEPS, W1 )
         END IF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_303
      SUBROUTINE DMUMPS_574
     &           ( IPOOL, LPOOL, III, LEAF, 
     &             INODE, STRATEGIE, MYID )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: STRATEGIE, MYID, LPOOL
      INTEGER IPOOL (LPOOL)
      INTEGER III,LEAF
      INTEGER, INTENT(OUT) :: INODE
         LEAF  = LEAF - 1
         INODE = IPOOL( LEAF )
      RETURN
      END SUBROUTINE DMUMPS_574
      SUBROUTINE DMUMPS_249(N, A, LA, IW, LIW, W, LWC,
     *    RHS, LRHS, NRHS, RHSSOL, 
     *    PTRICB, PTRACB, IWCB, LIWW, POOLSS, W2, 
     *    NE_STEPS, NA, LNA, STEP,
     *    FRERE, FILS, IPOOL, LPOOL, PTRIST, PTRFAC, MYLEAF, INFO, 
     *    PROCNODE_STEPS,
     *    SLAVEF, COMM,MYID, BUFR, LBUFR, LBUFR_BYTES,
     *    KEEP,KEEP8, RHS_ROOT, MTYPE, 
     *
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      USE DMUMPS_OOC
      USE DMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER N,LA,LIW,LIWW,LWC,LPOOL,LNA
      INTEGER SLAVEF,MYLEAF,COMM,MYID
      INTEGER KEEP( 500 )
      INTEGER*8 KEEP8(150)
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER NA(LNA),NE_STEPS(KEEP(28))
      INTEGER IPOOL(LPOOL),POOLSS(N)
      INTEGER INFO(40)
      INTEGER PTRIST(KEEP(28)),PTRFAC(KEEP(28)),
     *        PTRICB(KEEP(28)),PTRACB(KEEP(28))
      INTEGER LRHS, NRHS
      DOUBLE PRECISION A(LA), RHS(LRHS,NRHS), W(LWC),
     *                 RHSSOL(N, NRHS) 
      DOUBLE PRECISION W2(KEEP(133))
      INTEGER IW(LIW),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)),FILS(N)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      DOUBLE PRECISION RHS_ROOT( * )
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      LOGICAL FLAG
      INTEGER POSIWCB,POSWCB,K
      INTEGER APOS,APOSCB,NPIV
      INTEGER IPOS,IPOSCB,LIELL,NELIM,IFR,JJ,I
      INTEGER J1,J2,J,IST,NCB,NBFINF
      INTEGER NBLEAF,INODE,NBROOT,NROOT,NBFILS
      INTEGER IN,IF,LONG,IFRR,POOL_FIRST_POS,TMP
      INTEGER III,IIPOOL,MYLEAFE
      INTEGER IISSAR
      INTEGER NSLAVES
      DOUBLE PRECISION ALPHA,ONE,ZERO
      PARAMETER(ONE=1.0D0, ALPHA=-1.0D0, ZERO=0.0D0)
      LOGICAL BLOQ,DEBUT
      INTEGER PROCDEST, DEST
      INTEGER SSII,POSII, POSINDICES
      INTEGER DUMMY(1)
      INTEGER PLEFTW, LDA
      INTEGER Offset, EffectiveSize, ISLAVE, FirstIndex
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      LOGICAL DEJA_SEND( 0:SLAVEF-1 )
      LOGICAL DMUMPS_283
      INTEGER DMUMPS_330
      EXTERNAL DGEMV, DTRSV, DTRSM, DGEMM,
     *         DMUMPS_283, DMUMPS_330
      PLEFTW = 1
      POSIWCB = LIWW
      POSWCB = LWC
      NBFINF = SLAVEF
      MYLEAFE = MYLEAF
      IF (MYLEAFE .EQ. 0) THEN
        CALL DMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM, FEUILLE,
     *                  SLAVEF)
        NBFINF = NBFINF - 1
        IF (NBFINF .EQ. 0) GOTO 340
      ENDIF
      NROOT = 0
      NBLEAF = NA(1)
      NBROOT = NA(2)
      DO I = 1, NBROOT
        INODE = NA(NBLEAF+I+2)
        IF (DMUMPS_275(STEP(INODE),PROCNODE_STEPS,
     *      SLAVEF) .EQ. MYID) THEN
          NROOT = NROOT + 1
          IPOOL(NROOT) = INODE
        ENDIF
      END DO
      III = 1
      IIPOOL = NROOT + 1
      IISSAR = 1
   50 CONTINUE
      BLOQ = (III .EQ. IIPOOL .AND. IISSAR .LE. 1 )
      CALL DMUMPS_41( BLOQ, FLAG, BUFR, LBUFR,
     *     LBUFR_BYTES, MYID, SLAVEF, COMM,
     *     N, IWCB, LIWW, POSIWCB,
     *     W, LWC, POSWCB,
     *     IIPOOL, III, IISSAR,
     *     NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *     IPOOL, LPOOL, STEP,  FRERE, FILS, PROCNODE_STEPS,
     *     PLEFTW, KEEP,KEEP8,
     *     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *     RHSSOL, RHS, LRHS, NRHS )
      IF ( INFO(1) .LT. 0 ) GOTO 340
      IF ( .NOT. FLAG ) THEN
        IF (III .NE. IIPOOL) THEN
          INODE = IPOOL(III)
          III = III + 1
          GO TO 60
        ELSE IF ( IISSAR .GT. 1 ) THEN
          GO TO 200
        ENDIF
      END IF
      IF ( NBFINF .eq. 0 ) GOTO 340
      GOTO 50
   60 CONTINUE
      IF ( INODE .EQ. KEEP( 38 ) .OR. INODE .EQ. KEEP( 20 ) ) THEN
        IPOS = PTRIST(STEP(INODE))+XSIZE
        LIELL = IW(IPOS+3)
        NPIV  = LIELL
        IPOS =  PTRIST(STEP(INODE)) + 5 + XSIZE
        IF ( MTYPE .EQ. 1 .AND. KEEP(50) .EQ. 0) THEN
          J1   = IPOS + LIELL + 1
          J2   = IPOS + LIELL + NPIV
        ELSE
          J1   = IPOS + 1
          J2   = IPOS + NPIV
        END IF
        IFR  = 0
        DO JJ = J1, J2
          J  = IW( JJ )
          IFR = IFR + 1
          DO K=1,NRHS
            RHSSOL(J,K) = RHS_ROOT(IFR+NPIV*(K-1))
          ENDDO
        END DO 
        IN = INODE
  270   IN = FILS(IN)
        IF (IN .GT. 0) GOTO 270
        IF (IN .EQ. 0) THEN
          MYLEAFE = MYLEAFE - 1
          IF (MYLEAFE .EQ. 0) THEN
            CALL DMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *                       FEUILLE, SLAVEF )
            NBFINF = NBFINF - 1
            IF (NBFINF .EQ. 0) GOTO 340
          ENDIF
          GOTO 50
        ENDIF
        IF   = -IN
        LONG = NPIV
        NBFILS = NE_STEPS(STEP(INODE))
        DEBUT = .TRUE.
        DO I = 0, SLAVEF - 1
          DEJA_SEND( I ) = .FALSE.
        END DO
        POOL_FIRST_POS=IIPOOL
        DO I = 1, NBFILS
          IF (DMUMPS_275(STEP(IF),PROCNODE_STEPS,SLAVEF)
     *       .EQ. MYID) THEN
             IPOOL(IIPOOL) = IF
             IIPOOL = IIPOOL + 1
          ELSE
            PROCDEST = DMUMPS_275(STEP(IF),PROCNODE_STEPS,
     *                 SLAVEF)
            IF (.NOT. DEJA_SEND( PROCDEST ))  THEN
 600          CALL DMUMPS_78( NRHS, IF, 0, 0,
     *           LONG, LONG, IW( J1 ),
     *           RHS_ROOT( 1 ), PROCDEST,
     *           NOEUD, COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                CALL DMUMPS_41(
     *          .FALSE., FLAG,
     *          BUFR, LBUFR, LBUFR_BYTES,
     *          MYID, SLAVEF, COMM,
     *          N, IWCB, LIWW, POSIWCB,
     *          W, LWC, POSWCB,
     *          IIPOOL, III, IISSAR,
     *          NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *          IPOOL, LPOOL, STEP, FRERE, FILS, PROCNODE_STEPS,
     *          PLEFTW, KEEP,KEEP8,
     *          PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *          RHSSOL, RHS, LRHS, NRHS )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 600
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = LONG * KEEP(35) +
     *                      ( LONG + 2 ) * KEEP(34)
                GOTO 330
              END IF
              DEJA_SEND( PROCDEST ) = .TRUE.
            END IF
            IF ( IERR .NE. 0 ) CALL DMUMPS_ABORT()
          ENDIF
          IF = FRERE(STEP(IF))
        ENDDO
        DO I=1,(IIPOOL-POOL_FIRST_POS)/2
           TMP=IPOOL(POOL_FIRST_POS+I-1)
           IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
           IPOOL(IIPOOL-I)=TMP
        ENDDO
        GOTO 50
      END IF
      IF ( DMUMPS_283( STEP(INODE), PROCNODE_STEPS,
     *     SLAVEF ) ) THEN
        POOLSS(IISSAR) = INODE
        IISSAR = IISSAR + 1
        GO TO 50
      END IF
      IPOS = PTRIST(STEP(INODE)) + 2 + XSIZE
      LIELL = IW(IPOS-2)+IW(IPOS+1)
      NELIM = IW(IPOS-1)
      IPOS = IPOS + 1
      NPIV = IW(IPOS)
      IPOS = IPOS + 1
      IF (KEEP(201).NE.0) THEN
       IF(.NOT.DMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,
     $     PTRFAC,KEEP(28),A,LA,IERR))THEN
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
         CALL DMUMPS_578(INODE,PTRFAC,
     $        KEEP,KEEP8,A,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
         CALL DMUMPS_577(A(PTRFAC(IW(IPOS))),INODE,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
       ELSE
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
       ENDIF
      ENDIF
      APOS = PTRFAC(IW(IPOS))
      NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + XSIZE )
      IPOS = IPOS + 1 + NSLAVES
      IFRR = 0
      LONG = 0
      IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
        J1 = IPOS + 1
        J2 = IPOS + NPIV
      ELSE
        J1 = IPOS + LIELL + 1
        J2 = IPOS + NPIV + LIELL
      END IF
      IF ( POSIWCB - LIELL - 2 .LT. 0 .or.
     *     POSWCB - LIELL*NRHS .LT. PLEFTW - 1 ) THEN
        CALL DMUMPS_95( NRHS, N, KEEP(28), IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
        IF ( POSWCB - LIELL*NRHS .LT. PLEFTW - 1 ) THEN
          INFO( 1 ) = -11
          INFO( 2 ) = LIELL * NRHS - POSWCB - PLEFTW + 1
          GOTO 330
        END IF
        IF ( POSIWCB - LIELL - 2 .LT. 0 ) THEN
          INFO( 1 ) = -14
          INFO( 2 ) = LIELL + 2 - POSIWCB
          GO TO 330
        END IF
      END IF
      POSIWCB = POSIWCB - LIELL - 2
      POSWCB  = POSWCB - LIELL*NRHS
      PTRICB(STEP( INODE )) = POSIWCB + 1
      PTRACB(STEP( INODE )) = POSWCB  + 1
      IWCB( PTRICB(STEP( INODE ))     ) = LIELL
      IWCB( PTRICB(STEP( INODE )) + 1 ) = 1
      IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
         POSINDICES = IPOS + LIELL + 1
      ELSE
         POSINDICES = IPOS + 1
      END IF
      IWCB( PTRICB(STEP( INODE )) + 2 :
     * PTRICB(STEP( INODE )) + 1 + LIELL )
     *   = IW( POSINDICES : POSINDICES + LIELL - 1 )
      IFR = PTRACB(STEP( INODE )) - 1
      DO JJ = J1, J2
        J = IW(JJ)
        IFR = IFR + 1
        DO K=1, NRHS
          W(IFR+(K-1)*LIELL) = RHS(J,K)
        ENDDO
      END DO
      IF ( LIELL .GT. NPIV ) THEN
        IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
          J1 = IPOS + LIELL + NPIV + 1
          J2 = IPOS + 2 * LIELL
        ELSE
          J1 = IPOS + NPIV + 1
          J2 = IPOS + LIELL
        END IF
        DO 130 JJ = J1, J2
          J = IW(JJ)
          IFR = IFR + 1
          DO K=1, NRHS
            W(IFR+(K-1)*LIELL) = RHSSOL(J,K)
          ENDDO
  130   CONTINUE
        NCB = LIELL - NPIV
        IF (NPIV .EQ. 0) GOTO 160
        IF ( MTYPE .eq. 1 ) THEN
          IST = APOS + NPIV
          IF (NRHS == 1) THEN
            CALL DGEMV( 'T', NCB, NPIV, ALPHA, A(IST), LIELL,
     *              W(NPIV + PTRACB(STEP(INODE))), 1,
     *              ONE,
     *              W(PTRACB(STEP(INODE))), 1 )
          ELSE
            CALL DGEMM('T','N', NPIV, NRHS, NCB, ALPHA, A(IST), LIELL,
     *              W(NPIV+PTRACB(STEP(INODE))), LIELL, ONE,
     *              W(PTRACB(STEP(INODE))), LIELL)
          ENDIF
        ELSE
          IF ( KEEP(50) .eq. 0 ) THEN
            IST = APOS + NPIV * LIELL
          ELSE
            IST = APOS + NPIV * NPIV
          END IF
          IF ( DMUMPS_330(STEP(INODE),PROCNODE_STEPS,
     *         SLAVEF).eq.2 ) THEN
            IF ( NELIM .ne. 0 ) THEN
              IF ( NRHS == 1 ) THEN
                CALL DGEMV( 'N', NPIV, NELIM, ALPHA,
     *                A( IST ), NPIV,
     *                W( NPIV + PTRACB(STEP(INODE)) ),
     *                1, ONE,
     *                W(PTRACB(STEP(INODE))), 1 )
              ELSE
                CALL DGEMM( 'N', 'N', NPIV, NRHS, NELIM, ALPHA,
     *                A(IST), NPIV, W(NPIV+PTRACB(STEP(INODE))),LIELL,
     *                ONE, W(PTRACB(STEP(INODE))),LIELL)
              END IF
            ENDIF
            Offset        = NPIV + NELIM
            IWCB( PTRICB(STEP( INODE )) + 1 ) =
     *      IWCB( PTRICB(STEP( INODE )) + 1 ) + NSLAVES
            DO ISLAVE = 1, NSLAVES
              CALL DMUMPS_49( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB - NELIM, 
     &                NSLAVES, 
     &                EffectiveSize,
     &                FirstIndex )
 500         DEST = IW( PTRIST(STEP(INODE))+5+ISLAVE+XSIZE)
             CALL DMUMPS_63(NRHS, INODE,
     *             W(Offset+PTRACB(STEP(INODE))), EffectiveSize, LIELL,
     *             DEST,
     *             BACKSLV_MASTER2SLAVE,
     *             COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                CALL DMUMPS_41(
     *          .FALSE., FLAG,
     *          BUFR, LBUFR, LBUFR_BYTES,
     *          MYID, SLAVEF, COMM,
     *          N, IWCB, LIWW, POSIWCB,
     *          W, LWC, POSWCB,
     *          IIPOOL, III, IISSAR,
     *          NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *          IPOOL, LPOOL, STEP, FRERE, FILS,
     *          PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     *          PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *          RHSSOL, RHS, LRHS, NRHS )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 500
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = EffectiveSize * KEEP(35) +
     *                            2 * KEEP(34)
                GOTO 330
              END IF
              Offset = Offset + EffectiveSize
            END DO
            GOTO 50
          ELSE
            IF ( NRHS == 1 ) THEN
              CALL DGEMV( 'N', NPIV, NCB, ALPHA, A( IST ), NPIV,
     *                W( NPIV + PTRACB(STEP(INODE))),
     *                1, ONE,
     *                W(PTRACB(STEP(INODE))), 1 )
            ELSE
                CALL DGEMM( 'N', 'N', NPIV, NRHS, NCB, ALPHA,
     *                A(IST), NPIV, W(NPIV+PTRACB(STEP(INODE))),LIELL,
     *                ONE, W(PTRACB(STEP(INODE))),LIELL)
            END IF
          END IF ! node level =  2
        END IF ! MTYPE == 1
      ENDIF  ! LIELL > NPIV
      IF ( MTYPE .eq. 1 ) THEN
        IF ( NRHS == 1 ) THEN
          CALL DTRSV('L', 'T', 'N', NPIV, A(APOS), LIELL,
     *              W(PTRACB(STEP(INODE))), 1)
        ELSE
          CALL DTRSM('L','L','T','N', NPIV, NRHS, ONE, A(APOS),
     *              LIELL, W(PTRACB(STEP(INODE))), LIELL)
        ENDIF
      ELSE
        IF ( KEEP(50) .EQ. 0 ) THEN
          IF ( NRHS == 1 ) THEN
            CALL DTRSV('U','N','U', NPIV, A(APOS), LIELL,
     *              W(PTRACB(STEP(INODE))), 1)
          ELSE
            CALL DTRSM('L','U','N','U', NPIV, NRHS, ONE, A(APOS),
     *                 LIELL,W(PTRACB(STEP(INODE))),LIELL)
          END IF
        ELSE
          IF ( NRHS == 1 ) THEN
            CALL DTRSV('U','N','U', NPIV, A(APOS), NPIV,
     *              W(PTRACB(STEP(INODE))), 1)
          ELSE
            CALL DTRSM('L','U','N','U',NPIV, NRHS, ONE, A(APOS),
     *           NPIV, W(PTRACB(STEP(INODE))), LIELL)
          END IF
        END IF
      END IF
      IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0) THEN
        J1 = IPOS + LIELL + 1
      ELSE
        J1 = IPOS + 1
      END IF
      DO 150 I = 1, NPIV
        JJ = IW(J1 + I - 1)
        DO K=1, NRHS
          RHSSOL(JJ, K) = W(PTRACB(STEP(INODE))+I-1+(K-1)*LIELL)
        ENDDO
  150 CONTINUE
  160 CONTINUE
      IF (KEEP(201).NE.0) THEN
         CALL DMUMPS_598(INODE,PTRFAC,KEEP(28),
     $        A,LA,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
      ENDIF
      IN = INODE
  170 IN = FILS(IN)
      IF (IN .GT. 0) GOTO 170
      IF (IN .EQ. 0) THEN
        MYLEAFE = MYLEAFE - 1
        IF (MYLEAFE .EQ. 0) THEN
          CALL DMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *                     FEUILLE, SLAVEF )
          NBFINF = NBFINF - 1
          IF (NBFINF .EQ. 0) GOTO 340
        ENDIF
        GOTO 50
      ENDIF
      IF = -IN
      LONG = LONG + NPIV
      NBFILS = NE_STEPS(STEP(INODE))
      DEBUT = .TRUE.
      DO i = 0, SLAVEF - 1
        DEJA_SEND( i ) = .FALSE.
      END DO
      POOL_FIRST_POS=IIPOOL
      DO 190 I = 1, NBFILS
        IF (DMUMPS_275(STEP(IF),PROCNODE_STEPS,
     *      SLAVEF) .EQ. MYID) THEN
           IPOOL(IIPOOL) = IF
           IIPOOL = IIPOOL + 1
           IF = FRERE(STEP(IF))
        ELSE
          PROCDEST = DMUMPS_275(STEP(IF),PROCNODE_STEPS,SLAVEF)
          IF (.not. DEJA_SEND( PROCDEST ))  THEN
 400        CONTINUE
            CALL DMUMPS_78( NRHS, IF, 0, 0, LIELL,
     *         LIELL,
     *         IWCB( PTRICB(STEP( INODE )) + 2),
     *         W   ( PTRACB(STEP( INODE ))), PROCDEST,
     *         NOEUD, COMM, IERR )
            IF ( IERR .EQ. -1 ) THEN
              CALL DMUMPS_41(
     *        .FALSE., FLAG,
     *        BUFR, LBUFR, LBUFR_BYTES,
     *        MYID, SLAVEF, COMM,
     *        N, IWCB, LIWW, POSIWCB,
     *        W, LWC, POSWCB,
     *        IIPOOL, III, IISSAR,
     *        NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *        IPOOL, LPOOL, STEP, FRERE, FILS, PROCNODE_STEPS,
     *        PLEFTW, KEEP,KEEP8,
     *        PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *        RHSSOL, RHS, LRHS, NRHS )
              IF ( INFO( 1 ) .LT. 0 ) GOTO 340
              GOTO 400
            ELSE IF ( IERR .EQ. -2 ) THEN
              INFO( 1 ) = -17
              INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
              GOTO 330
            END IF
            DEJA_SEND( PROCDEST ) = .TRUE.
          END IF
          IF = FRERE(STEP(IF))
        ENDIF
  190 CONTINUE
      DO I=1,(IIPOOL-POOL_FIRST_POS)/2
         TMP=IPOOL(POOL_FIRST_POS+I-1)
         IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
         IPOOL(IIPOOL-I)=TMP
      ENDDO
      IWCB(PTRICB(STEP(INODE))+1) = IWCB(PTRICB(STEP(INODE))+1)-1
      CALL DMUMPS_151(NRHS,N, KEEP(28), IWCB, LIWW, W, LWC,
     *     POSWCB,POSIWCB,PTRICB,PTRACB)
      GOTO 50
  200 CONTINUE
      IISSAR = IISSAR - 1
      INODE  = POOLSS( IISSAR )
      SSII  = IISSAR
      POSII = IISSAR + 1
  230 CONTINUE
      POSII = POSII - 1
      INODE = POOLSS( POSII )
      IPOS = PTRIST(STEP(INODE)) + 2 + XSIZE
      LIELL = IW(IPOS-2)+IW(IPOS+1)
      IPOS = IPOS + 1
      NPIV = IW(IPOS)
      IPOS = IPOS + 1
      IF (KEEP(201).NE.0) THEN
       IF(.NOT.DMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,
     $     PTRFAC,KEEP(28),A,LA,IERR))THEN
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF         
         CALL DMUMPS_578(INODE,PTRFAC,
     $        KEEP,KEEP8,A,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
         CALL DMUMPS_577(A(PTRFAC(IW(IPOS))),INODE,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
       ELSE
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
       ENDIF
      ENDIF
      IF (NPIV .EQ. 0) GOTO 280
      APOS = PTRFAC(IW(IPOS))
      NSLAVES = IW( PTRIST(STEP(INODE)) +5 + XSIZE)
      IPOS = IPOS + 1 + NSLAVES
      IF ( MTYPE .EQ. 1 .OR. KEEP(50).NE.0 ) THEN
        J1 = IPOS + 1
        J2 = IPOS + NPIV
      ELSE
        J1 = IPOS + 1 + LIELL
        J2 = IPOS + NPIV + LIELL
      END IF
      IFR = 0
      IF ( POSWCB .LT. LIELL*NRHS ) THEN
          CALL DMUMPS_95( NRHS, N, KEEP(28), IWCB, LIWW, W, LWC,
     *                 POSWCB, POSIWCB, PTRICB, PTRACB)
          IF ( POSWCB .LT. LIELL*NRHS ) THEN
            INFO(1) = -11
            INFO(2) = LIELL*NRHS - POSWCB
            GOTO 330
          END IF
      END IF
      DO 240 JJ = J1, J2
        J = IW(JJ)
        IFR = IFR + 1
        DO K=1,NRHS
          W(IFR+(K-1)*LIELL) = RHS(J,K)
        ENDDO
  240 CONTINUE
      IF (LIELL .GT. NPIV) THEN
        IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
          J1 = IPOS + LIELL + NPIV + 1
          J2 = IPOS + 2 * LIELL
        ELSE
          J1 = IPOS + NPIV + 1
          J2 = IPOS + LIELL
        END IF
        DO 250 JJ = J1, J2
          J = IW(JJ)
          IFR = IFR + 1
          DO K=1,NRHS
            W(IFR+(K-1)*LIELL) = RHSSOL(J,K)
          ENDDO
  250   CONTINUE
        NCB = LIELL - NPIV
        IF ( MTYPE .EQ. 1) THEN
          IST = APOS + NPIV
          IF ( NRHS == 1 ) THEN
            CALL DGEMV( 'T', NCB, NPIV, ALPHA, A(IST), LIELL,
     *                W(NPIV+1), 1, ONE,
     *                W(1), 1 )
          ELSE
            CALL DGEMM( 'T', 'N', NPIV, NRHS, NCB, ALPHA, A(IST),
     *                LIELL, W(NPIV+1), LIELL, ONE, W(1), LIELL)
          END IF
        ELSE
          IF ( KEEP(50) .eq.0 ) THEN
            IST = APOS + NPIV * LIELL
          ELSE
            IST = APOS + NPIV * NPIV
          END IF
          IF ( NRHS == 1 ) THEN
            CALL DGEMV( 'N', NPIV, NCB, ALPHA, A(IST), NPIV,
     *                W(NPIV+1), 1, ONE,
     *                W(1), 1 )
          ELSE
            CALL DGEMM( 'N','N', NPIV, NRHS, NCB, ALPHA, A(IST),
     *                NPIV, W(NPIV+1), LIELL, ONE, W(1),LIELL)
          END IF
        END IF
      ENDIF
      IF ( MTYPE .EQ. 1) THEN
        IF (NRHS == 1) THEN
            CALL DTRSV('L', 'T', 'N', NPIV, A(APOS), LIELL,
     *               W(1), 1)
        ELSE
          CALL DTRSM('L','L','T','N',NPIV, NRHS, ONE, A(APOS),
     *                LIELL, W(1), LIELL)
        ENDIF
      ELSE
        IF ( KEEP(50) .eq. 0 ) THEN
          LDA = LIELL
        ELSE
          LDA = NPIV
        ENDIF
        IF ( NRHS == 1 ) THEN
            CALL DTRSV('U', 'N', 'U', NPIV, A(APOS), LDA,
     *               W(1), 1)
        ELSE
          CALL DTRSM('L','U','N','U', NPIV, NRHS, ONE, A(APOS),
     *               LDA, W(1), LIELL)
        END IF
      END IF
      IFR = 0
      IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
        J1 = IPOS + LIELL + 1
      ELSE
        J1 = IPOS + 1
      END IF
      J2 = J1 + NPIV - 1
      DO I = 1, NPIV
        JJ = IW(J1 + I - 1)
        DO K=1,NRHS
          RHSSOL(JJ,K) = W(I+(K-1)*LIELL)
        ENDDO
      ENDDO
  280 CONTINUE
      IF (KEEP(201).NE.0) THEN
       CALL DMUMPS_598(INODE,PTRFAC,KEEP(28),A,
     $        LA,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
      ENDIF
      IN = INODE
  290 IN = FILS(IN)
      IF (IN .GT. 0) GOTO 290
      IF (IN .EQ. 0) THEN
        MYLEAFE = MYLEAFE - 1
        IF (MYLEAFE .EQ. 0) THEN
          CALL DMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID,
     *        COMM, FEUILLE, SLAVEF)
          NBFINF = NBFINF - 1
          IF (NBFINF .EQ. 0) GOTO 340
        ENDIF
        IF (SSII .EQ. POSII) GOTO 50
        GOTO 230
      ENDIF
      IF = -IN
      NBFILS = NE_STEPS(STEP(INODE))
      DO 300 I = 1, NBFILS
        POOLSS((POSII-I+1)+NBFILS-I) = IF
        POSII = POSII + 1
        IF = FRERE(STEP(IF))
  300 CONTINUE
      GOTO 230
  330 CONTINUE
      CALL DMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM, TERREUR,
     * SLAVEF)
  340 CONTINUE
      CALL DMUMPS_150( MYID,COMM,BUFR,
     *                            LBUFR,LBUFR_BYTES )
      RETURN
      END SUBROUTINE DMUMPS_249
      RECURSIVE SUBROUTINE DMUMPS_41(
     *     BLOQ, FLAG,
     *     BUFR, LBUFR, LBUFR_BYTES,
     *     MYID, SLAVEF, COMM,
     *     N, IWCB, LIWW, POSIWCB,
     *     W, LWC, POSWCB,
     *     IIPOOL, III, IISSAR,
     *     NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *     IPOOL, LPOOL, STEP, FRERE, FILS, PROCNODE_STEPS,
     *     PLEFTW, KEEP,KEEP8,
     *     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, RHSSOL,
     *     RHS, LRHS, NRHS )
      IMPLICIT NONE
      LOGICAL BLOQ, FLAG
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LIWW
      INTEGER IWCB( LIWW )
      INTEGER LWC
      DOUBLE PRECISION W( LWC )
      INTEGER POSIWCB, POSWCB
      INTEGER IIPOOL, III, IISSAR, LPOOL
      INTEGER IPOOL( LPOOL ), POOLSS( N )
      INTEGER NBFINF, INFO(40)
      INTEGER PLEFTW, KEEP( 500)
      INTEGER*8 KEEP8(150)
      INTEGER PROCNODE_STEPS( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER PTRICB(KEEP(28)), PTRACB(KEEP(28)), STEP( N ), FILS( N )
      INTEGER LIW, LA
      INTEGER PTRIST(KEEP(28)), PTRFAC(KEEP(28)), IW( LIW )
      DOUBLE PRECISION A( LA ), W2( KEEP(133) )
      INTEGER LRHS, NRHS
      DOUBLE PRECISION RHSSOL( N, NRHS ), RHS(LRHS, NRHS)
      INTEGER MYLEAFE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MSGSOU, MSGTAG, MSGLEN
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
        CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     *                   COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     *                   FLAG, STATUS, IERR )
      END IF
      IF (FLAG) THEN
         MSGSOU=STATUS(MPI_SOURCE)
         MSGTAG=STATUS(MPI_TAG)
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL DMUMPS_44( MYID, SLAVEF, COMM )
         ELSE
           CALL MPI_RECV(BUFR, LBUFR_BYTES, MPI_PACKED, MSGSOU,
     *                   MSGTAG, COMM, STATUS, IERR)
           CALL DMUMPS_42( MSGTAG, MSGSOU,
     *                BUFR, LBUFR, LBUFR_BYTES,
     *                MYID, SLAVEF, COMM,
     *                N, IWCB, LIWW, POSIWCB,
     *                W, LWC, POSWCB,
     *                IIPOOL, III, IISSAR,
     *                NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *                IPOOL, LPOOL, STEP,
     *                FRERE, FILS, PROCNODE_STEPS, PLEFTW,
     *                KEEP,KEEP8,
     *                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *                RHSSOL, RHS, LRHS, NRHS )
         END IF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_41
      RECURSIVE SUBROUTINE DMUMPS_42(
     *                MSGTAG, MSGSOU,
     *                BUFR, LBUFR, LBUFR_BYTES,
     *                MYID, SLAVEF, COMM,
     *                N, IWCB, LIWW, POSIWCB,
     *                W, LWC, POSWCB,
     *                IIPOOL, III, IISSAR,
     *                NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *                IPOOL, LPOOL, STEP,
     *                FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     *                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *                RHSSOL, RHS, LRHS, NRHS )
      USE DMUMPS_OOC
      USE DMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MSGTAG, MSGSOU
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LIWW
      INTEGER IWCB( LIWW )
      INTEGER LWC
      DOUBLE PRECISION W( LWC )
      INTEGER POSIWCB, POSWCB
      INTEGER IIPOOL, III, IISSAR, LPOOL
      INTEGER IPOOL( LPOOL ), POOLSS( N )
      INTEGER NBFINF, INFO(40)
      INTEGER PLEFTW, KEEP( 500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRICB(KEEP(28)), PTRACB(KEEP(28)), STEP( N ), FILS( N )
      INTEGER FRERE(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER LIW, LA
      INTEGER IW( LIW ), PTRIST( KEEP(28) ), PTRFAC(KEEP(28))
      DOUBLE PRECISION A( LA ), W2( KEEP(133) )
      INTEGER LRHS, NRHS
      DOUBLE PRECISION RHSSOL(N, NRHS), RHS(LRHS, NRHS)
      INTEGER MYLEAFE, K
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER POSITION, IF, INODE, IERR, LONG, DUMMY(1)
      INTEGER P_UPDATE, P_SOL_MAS, LIELL
      INTEGER NPIV, NPIV_CHECK, NROW_L, APOS, IPOS, NROW_RECU
      INTEGER I, JJ, IN, PROCDEST, LDA
      LOGICAL FLAG
      DOUBLE PRECISION ZERO, ALPHA, ONE
      PARAMETER( ZERO = 0.0D0, ALPHA = -1.0D0, ONE = 1.0D0)
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      LOGICAL DEJA_SEND( 0:SLAVEF-1 )
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275, DTRSV, DTRSM, DGEMV, DGEMM
      IF (MSGTAG .EQ. FEUILLE) THEN
          NBFINF = NBFINF - 1
      ELSE IF (MSGTAG .EQ. NOEUD) THEN
          POSITION = 0
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        INODE, 1, MPI_INTEGER,
     *        COMM, IERR)
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        LONG, 1, MPI_INTEGER,
     *        COMM, IERR)
          IF (   POSIWCB - LONG - 2 .LT. 0
     *      .OR. POSWCB - PLEFTW + 1 .LT. LONG ) THEN
            CALL DMUMPS_95(NRHS, N, KEEP(28), IWCB,
     *      LIWW, W, LWC,
     *      POSWCB, POSIWCB, PTRICB, PTRACB)
            IF ((POSIWCB - LONG - 2 ) .LT. 0) THEN
              INFO(1)=-14
              INFO(2)=-POSIWCB + LONG + 2
              GOTO 260
            END IF
            IF ( POSWCB - PLEFTW + 1 .LT. LONG ) THEN
              INFO(1) = -11
              INFO(2) = LONG + PLEFTW - POSWCB - 1
              GOTO 260
            END IF
          ENDIF
          POSIWCB = POSIWCB - LONG
          POSWCB = POSWCB - LONG
          IF (LONG .GT. 0) THEN
            CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          IWCB(POSIWCB + 1), 
     *          LONG, MPI_INTEGER, COMM, IERR)
            DO K=1,NRHS
             CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          W(POSWCB + 1), LONG, 
     *          MPI_DOUBLE_PRECISION, COMM, IERR)
             DO JJ=0, LONG-1
               RHSSOL(IWCB(POSIWCB+1+JJ),K) = W(POSWCB+1+JJ)
             ENDDO
            ENDDO
            POSIWCB = POSIWCB + LONG
            POSWCB = POSWCB + LONG
          ENDIF
          IF = FRERE( STEP(INODE) )
          DO WHILE ( IF .GT. 0 )
            IF ( DMUMPS_275(STEP(IF),PROCNODE_STEPS,
     *           SLAVEF) .eq. MYID ) THEN
              IPOOL( IIPOOL ) = IF
              IIPOOL = IIPOOL + 1
            END IF
            IF = FRERE( STEP( IF ) )
          END DO
          IPOOL( IIPOOL ) = INODE
          IIPOOL = IIPOOL + 1
      ELSE IF ( MSGTAG .EQ. BACKSLV_MASTER2SLAVE ) THEN
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   INODE, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   NROW_RECU, 1, MPI_INTEGER, COMM, IERR )
        IPOS   = PTRIST( STEP(INODE) ) + XSIZE
        NPIV   = - IW( IPOS     )
        NROW_L =   IW( IPOS + 1 )
        IF (KEEP(201).NE.0) THEN
         IF(.NOT.DMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,
     $       PTRFAC,KEEP(28),A,LA,IERR))THEN
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           CALL DMUMPS_578(INODE,PTRFAC,
     $          KEEP,KEEP8,A,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           CALL DMUMPS_577(A(PTRFAC(IW( IPOS + 3 ))),INODE,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
         ELSE
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
         ENDIF
        ENDIF
        APOS   =   PTRFAC(IW( IPOS + 3 ))
        IF ( NROW_L .NE. NROW_RECU ) THEN
          WRITE(*,*) 'Error1 in 41S : NROW L/RECU=',NROW_L, NROW_RECU
          CALL DMUMPS_ABORT()
        END IF
        LONG = NROW_L + NPIV
        IF ( POSWCB - LONG*NRHS .LT. PLEFTW - 1 ) THEN
           CALL DMUMPS_95(NRHS, N, KEEP(28), IWCB,
     *          LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
           IF ( POSWCB - LONG*NRHS .LT. PLEFTW - 1 ) THEN
             INFO(1) = -11
             INFO(2) = LONG * NRHS- POSWCB
             GOTO 260
           END IF
        END IF
        P_UPDATE  = PLEFTW
        P_SOL_MAS = PLEFTW + NPIV * NRHS
        PLEFTW    = P_SOL_MAS + NROW_L * NRHS
        DO K=1, NRHS
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   W( P_SOL_MAS+(K-1)*NROW_L),NROW_L,
     *                   MPI_DOUBLE_PRECISION,
     *                   COMM, IERR )
        ENDDO
        IF ( NRHS == 1 ) THEN
          CALL DGEMV( 'N', NPIV, NROW_L, ALPHA, A( APOS ), NPIV,
     *              W( P_SOL_MAS ), 1, ZERO,
     *              W( P_UPDATE ), 1 )
        ELSE
          CALL DGEMM( 'N', 'N', NPIV, NRHS, NROW_L, ALPHA, A(APOS),
     *                NPIV, W( P_SOL_MAS ), NROW_L, ZERO, W( P_UPDATE ),
     *                NPIV )
        END IF
        IF (KEEP(201).NE.0) THEN
         CALL DMUMPS_598(INODE,PTRFAC,KEEP(28),
     $          A,LA,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 260
         ENDIF
        ENDIF
        PLEFTW = PLEFTW - NROW_L * NRHS
 100    CONTINUE
        CALL DMUMPS_63( NRHS, INODE, W(P_UPDATE),
     *                               NPIV, NPIV,
     *                                MSGSOU, 
     *                                BACKSLV_UPDATERHS,
     *                                COMM, IERR )
        IF ( IERR .EQ. -1 ) THEN
          CALL DMUMPS_41(
     *     .FALSE., FLAG,
     *     BUFR, LBUFR, LBUFR_BYTES,
     *     MYID, SLAVEF, COMM,
     *     N, IWCB, LIWW, POSIWCB,
     *     W, LWC, POSWCB,
     *     IIPOOL, III, IISSAR,
     *     NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *     IPOOL, LPOOL, STEP,
     *     FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     *     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE,
     *     RHSSOL, RHS, LRHS, NRHS )
          IF ( INFO( 1 ) .LT. 0 ) GOTO 270
          GOTO 100
        ELSE IF ( IERR .EQ. -2 ) THEN
          INFO( 1 ) = -17
          INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
          GOTO 260
        END IF
        PLEFTW = PLEFTW - NPIV * NRHS
      ELSE IF ( MSGTAG .EQ. BACKSLV_UPDATERHS ) THEN
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   INODE, 1, MPI_INTEGER, COMM, IERR )
          IPOS  = PTRIST(STEP(INODE)) + 2 + XSIZE
          LIELL = IW(IPOS-2)+IW(IPOS+1)
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   NPIV, 1, MPI_INTEGER, COMM, IERR )
        DO K=1, NRHS
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   W2, NPIV, MPI_DOUBLE_PRECISION,
     *                   COMM, IERR )
        W( (K-1)*LIELL+PTRACB(STEP(INODE)):
     *     (K-1)*LIELL+PTRACB(STEP(INODE)) + NPIV - 1 )
     *    =   W( (K-1)*LIELL+PTRACB(STEP(INODE)):
     *     (K-1)*LIELL+PTRACB(STEP(INODE)) + NPIV - 1 )
     *    +   W2( 1: NPIV )
        ENDDO
        IWCB( PTRICB(STEP(INODE)) + 1 )
     *    = IWCB( PTRICB(STEP(INODE)) + 1 ) - 1
        IF ( IWCB( PTRICB(STEP(INODE)) + 1 ) .EQ. 1 ) THEN
          IPOS  = PTRIST(STEP(INODE)) + 2 + XSIZE
          LIELL = IW(IPOS-2)+IW(IPOS+1)
          IPOS  = IPOS + 1
          NPIV_CHECK = IW( IPOS )
          IF ( NPIV .NE. NPIV_CHECK ) THEN
            WRITE(*,*) 'Error 2 in usedbympima41s.F'
            CALL DMUMPS_ABORT()
          END IF
          IPOS  = IPOS + 1
          IF (KEEP(201).NE.0) THEN
           IF(.NOT.DMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,
     $         PTRFAC,KEEP(28),A,LA,IERR))THEN
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
             CALL DMUMPS_578(INODE,PTRFAC,
     $            KEEP,KEEP8,A,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
             CALL DMUMPS_577(A(PTRFAC(IW( IPOS ))),INODE,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
           ELSE
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
           ENDIF
          ENDIF
          APOS  = PTRFAC(IW( IPOS ))
          IF ( KEEP(50) .eq. 0 ) THEN
           LDA = LIELL
          ELSE
           LDA = NPIV
          ENDIF
          IF ( NRHS == 1 ) THEN
              CALL DTRSV( 'U', 'N', 'U', NPIV, A(APOS), LDA,
     *                  W(PTRACB(STEP(INODE))),1)
          ELSE
             CALL DTRSM( 'L','U', 'N', 'U', NPIV, NRHS, ONE,
     *                   A(APOS), LDA,
     *                   W(PTRACB(STEP(INODE))),LIELL)
          END IF
        IF (KEEP(201).NE.0)THEN
           CALL DMUMPS_598(INODE,PTRFAC,KEEP(28),
     $          A,LA,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
        ENDIF
          IPOS = PTRIST(STEP(INODE)) + 5 + XSIZE
          IPOS = IW( IPOS ) + IPOS + 1
          DO I = 1, NPIV
            JJ = IW( IPOS + I - 1 )
            DO K=1,NRHS
              RHSSOL( JJ, K ) = W( PTRACB(STEP(INODE))+I-1
     *         + (K-1)*LIELL )
            ENDDO
          END DO
          IN = INODE
  200     IN = FILS(IN)
          IF (IN .GT. 0) GOTO 200
          IF (IN .EQ. 0) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
              CALL DMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *                       FEUILLE, SLAVEF )
              NBFINF = NBFINF - 1
            ENDIF
            IWCB( PTRICB(STEP(INODE)) + 1 )
     *        = IWCB( PTRICB(STEP(INODE)) + 1 ) - 1
            CALL DMUMPS_151(NRHS, N, KEEP(28),
     *          IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
            GOTO 270
          ENDIF
          DO I = 0, SLAVEF - 1
            DEJA_SEND( I ) = .FALSE.
          END DO
          IN = -IN
 300      CONTINUE
            IF (DMUMPS_275(STEP(IN),PROCNODE_STEPS,
     *          SLAVEF) .EQ. MYID) THEN
              IPOOL(IIPOOL ) = IN
              IIPOOL = IIPOOL + 1
            ELSE
              PROCDEST = DMUMPS_275( STEP(IN), PROCNODE_STEPS,
     *                   SLAVEF )
              IF ( .NOT. DEJA_SEND( PROCDEST ) ) THEN
 110            CALL DMUMPS_78( NRHS, IN, 0, 0,
     *          LIELL, LIELL,
     *          IWCB( PTRICB(STEP(INODE))+2 ),
     *          W( PTRACB(STEP(INODE))),
     *          PROCDEST, NOEUD, COMM, IERR )
                IF ( IERR .EQ. -1 ) THEN
                  CALL DMUMPS_41(
     *            .FALSE., FLAG,
     *            BUFR, LBUFR, LBUFR_BYTES,
     *            MYID, SLAVEF, COMM,
     *            N, IWCB, LIWW, POSIWCB,
     *            W, LWC, POSWCB,
     *            IIPOOL, III, IISSAR,
     *            NBFINF, PTRICB, PTRACB, POOLSS, INFO,
     *            IPOOL, LPOOL, STEP,
     *            FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     *            PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, RHSSOL, 
     *            RHS, LRHS, NRHS)
                  IF ( INFO( 1 ) .LT. 0 ) GOTO 270
                  GOTO 110
                ELSE IF ( IERR .eq. -2 ) THEN
                  INFO(1) = -17
                  INFO(2) = LIELL * NRHS * KEEP(35) +
     *                    ( LIELL + 2 ) * KEEP(34)
                  GOTO 260
                END IF
                DEJA_SEND( PROCDEST ) = .TRUE.
              END IF
            END IF
          IN = FRERE( STEP( IN ) )
          IF ( IN .GT. 0 ) GOTO 300
          IWCB( PTRICB(STEP( INODE )) + 1 )
     *      = IWCB( PTRICB(STEP( INODE )) + 1 ) - 1
          CALL DMUMPS_151(NRHS, N, KEEP(28),
     *          IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
        END IF
      ELSE IF (MSGTAG.EQ.TERREUR) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GO TO 270
       ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     *      (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
      ELSE
          INFO(1) = -100
          INFO(2) = MSGTAG
          GOTO 260
      ENDIF
      GO TO 270
 260  CONTINUE
      CALL DMUMPS_44( MYID, SLAVEF, COMM )
 270  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_42
#ifdef try_null_space_not_tested
      SUBROUTINE DMUMPS_352( NRHS,DESCA_PAR, DESCB_PAR,
     *  CNTXT_PAR,LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *  IPIV,LPIV,MASTER_ROOT,MYID,COMM,
     *  RHS_SEQ,SIZE_ROOT,A,INFO,MTYPE,LDLT,QR,
     *  QR_TAU,WK,LWK,DEFICIENCY,MAXG,GIND,GROW,GCOS,GSIN )
      IMPLICIT NONE
      INTEGER NRHS, MTYPE
      INTEGER DESCA_PAR( 9 ), DESCB_PAR( 9 )
      INTEGER LOCAL_M, LOCAL_N, MBLOCK, NBLOCK
      INTEGER CNTXT_PAR, MASTER_ROOT, SIZE_ROOT
      INTEGER MYID, COMM, DEFICIENCY
      INTEGER LPIV, IPIV( LPIV )
      INTEGER INFO(40), LDLT, QR, LWK
      DOUBLE PRECISION WK( LWK )
      DOUBLE PRECISION RHS_SEQ( * )
      DOUBLE PRECISION QR_TAU( * )
      DOUBLE PRECISION A( LOCAL_M, LOCAL_N )
      INTEGER MAXG,GIND
      DOUBLE PRECISION GROW(*),GCOS(*),GSIN(*)
      INCLUDE 'mpif.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER IERR, NPROW, NPCOL, MYROW, MYCOL
      INTEGER allocok
      INTEGER J, K, PRO, LOC, PRO2, LOC2, ISOUR,II
      INTEGER SEQ,NN,LDA
      DOUBLE PRECISION C, S, TEMP(2)
      DOUBLE PRECISION RHS_PAR( LOCAL_M, 1 )
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV_SEQ
      NN=SIZE_ROOT
      LDA=DESCA_PAR(9)
      IF ( NRHS > 1 ) THEN
        WRITE(*,*) "Internal error: experimental nullspace"
     *  "functionalities not compatible with multiple RHS"
        STOP
      ENDIF
      CALL BLACS_GRIDINFO( CNTXT_PAR, NPROW, NPCOL, MYROW, MYCOL )
      IF(MBLOCK.GE.NN) THEN
        SEQ=1
      ELSE
        SEQ=0
      END IF
      IF ( MYID .eq. MASTER_ROOT ) THEN
        ALLOCATE( IPIV_SEQ( SIZE_ROOT), stat=allocok )
        IF ( allocok.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = SIZE_ROOT
          CALL DMUMPS_ABORT()
        ENDIF
      ELSE
        ALLOCATE( IPIV_SEQ( 1 ), stat=allocok )
        IF ( allocok.GT.0 ) THEN
          INFO(1)  = -13
          INFO(2)  = 1
          CALL DMUMPS_ABORT()
        ENDIF
      END IF
      IF(SEQ.EQ.1) THEN
        IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
          CALL DCOPY(NN,RHS_SEQ,1,RHS_PAR,1)
        END IF
      ELSE
        IF(.not.(QR.EQ.1.AND.MTYPE.EQ.0)) THEN
          CALL DMUMPS_290( MYID, SIZE_ROOT, 1, RHS_SEQ,
     *      LOCAL_M, 1,
     *      MBLOCK, NBLOCK, RHS_PAR, MASTER_ROOT,
     *      NPROW, NPCOL, COMM )
         END IF
      END IF
      IF(SEQ.EQ.1) THEN
        IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
        IF(QR.EQ.1.OR.QR.EQ.2.OR.QR.EQ.3.OR.QR.EQ.4) THEN
          IF(MTYPE.EQ.1 ) then
            CALL DORMQR( 'L','T',NN,1,NN,
     *        A,LDA,QR_TAU,
     *        RHS_PAR,LDA,WK(1),LWK,IERR)
            IF ( IERR .LT. 0 ) THEN
              WRITE(*,*) 'Error return from DORMQR: IERR=', IERR
              CALL DMUMPS_ABORT()
            END IF
            DO J=1,GIND
              K=GROW(J)
              C=GCOS(J)
              S=GSIN(J)
              CALL DMUMPS_31(1,C,S,RHS_PAR(K,1),LDA,
     *             RHS_PAR(K+1,1),LDA)
            END DO
          END IF
          IF(MTYPE.eq.1) THEN
            CALL DTRTRS('U','N','N',NN-DEFICIENCY,1,
     *        A,LDA,RHS_PAR,LDA,IERR)
          ELSEIF(mtype.eq.0) then
            IF(myid.eq.master_root) then
              CALL DMUMPS_334(NN,IPIV,RHS_SEQ,WK)
            END IF
            CALL DCOPY(NN,RHS_SEQ,1,RHS_PAR,1)
            CALL DTRTRS('U','T','N',NN-DEFICIENCY,1,
     *        A,LDA,RHS_PAR,LDA,IERR)
            DO J=1,GIND
              K=GROW(J)
                C=GCOS(J)
                S=GSIN(J)
               ii=1
         CALL DMUMPS_31(1,C,S,RHS_PAR(K,1),ii,RHS_PAR(K+1,1),ii)
            END DO
            CALL DORMQR( 'L','N',NN,1,NN,A,LDA,QR_TAU,
     *        RHS_PAR,LDA,WK,LWK,IERR)
            IF ( IERR .LT. 0 ) THEN
              WRITE(*,*) 'Error return from DORMQR: IERR=', IERR
              CALL DMUMPS_ABORT()
            END IF
          ENDIF
        ELSEIF(QR.EQ.5.OR.QR.EQ.6) THEN
         IF ( MTYPE .eq. 1 ) THEN
         CALL DGETRS('N',NN-DEFICIENCY,1,A,LDA,IPIV,
     *     RHS_PAR,LDA,IERR)
         ELSE
          CALL DGETRS('T',NN-DEFICIENCY,1,A,LDA,IPIV,
     *      RHS_PAR,LDA,IERR)
         END IF
        ELSEIF(QR.EQ.7.OR.QR.EQ.8) THEN
         IF ( MTYPE .eq. 1 ) THEN
           CALL DGETRS('N',NN-DEFICIENCY,1,A,LDA,IPIV,
     *       RHS_PAR,LDA,IERR)
         ELSE
           CALL DGETRS('T',NN-DEFICIENCY,1,A,LDA,IPIV,
     *       RHS_PAR,LDA,IERR)
         END IF
        ELSE IF ( LDLT .eq. 0 .OR. LDLT .eq. 2 ) THEN
          IF ( MTYPE .eq. 1 ) THEN
           CALL DGETRS('N',NN-DEFICIENCY,1,A,LDA,IPIV,
     *       RHS_PAR,LDA,IERR)
          ELSE
           CALL DGETRS('T',NN-DEFICIENCY,1,A,LDA,IPIV,
     *       RHS_PAR,LDA,IERR)
          END IF
        ELSE
        END IF
        IF ( IERR .LT. 0 ) THEN
          WRITE(*,*) ' Problem during solve of the root'
          CALL DMUMPS_ABORT()
        END IF
        END IF
        ELSE
        IF(QR.EQ.1.OR.QR.EQ.2.OR.QR.EQ.3.OR.QR.EQ.4) THEN
          IF(MTYPE.EQ.1 ) then
            CALL DMUMPS_262( 'L','T', SIZE_ROOT, 1, SIZE_ROOT,
     *        A, 1, 1, DESCA_PAR, QR_TAU,
     *        RHS_PAR, 1, 1, DESCB_PAR,
     *        WK(1), LWK, INFO(1) )
            IF ( INFO(1) .LT. 0 ) THEN
              WRITE(*,*) 'Error return from PDORMQR: INFO(1)=', INFO(1)
              CALL DMUMPS_ABORT()
            END IF
            DO J=1,GIND
              K=GROW(J)
              IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
                C=GCOS(J)
                S=GSIN(J)
              END IF
              CALL DMUMPS_161(MYROW,MBLOCK,NPCOL,K,PRO)
              CALL DMUMPS_160(MYROW,MBLOCK,NPCOL,K,LOC)
              IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
                CALL MPI_SEND(RHS_PAR(LOC,1),1,
     *            MPI_DOUBLE_PRECISION,
     *            0,163,COMM,IERR)
              END IF
              IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
                ISOUR=PRO*NPCOL+0
                CALL MPI_RECV(TEMP(1),1,MPI_DOUBLE_PRECISION,
     *            ISOUR,163,COMM,STATUS,IERR)
              END IF
              CALL DMUMPS_161(MYROW,MBLOCK,NPCOL,K+1,PRO2)
              CALL DMUMPS_160(MYROW,MBLOCK,NPCOL,K+1,LOC2)
              IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO2) THEN
                CALL MPI_SEND(RHS_PAR(LOC2,1),1,
     *            MPI_DOUBLE_PRECISION,
     *            0,164,COMM,IERR)
              END IF
              IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
                ISOUR=PRO2*NPCOL+0
                CALL MPI_RECV(TEMP(2),1,MPI_DOUBLE_PRECISION,
     *            ISOUR,164,COMM,STATUS,IERR)
              END IF
              IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
                CALL DMUMPS_31(1,C,S,TEMP(1),1,TEMP(2),1)
              END IF
              IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
                ISOUR=PRO*NPCOL+0
                CALL MPI_SEND(TEMP(1),1,
     *            MPI_DOUBLE_PRECISION,
     *            ISOUR,163,COMM,IERR)
              END IF
              IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
                ISOUR=0
                CALL MPI_RECV(RHS_PAR(LOC,1),1,
     *            MPI_DOUBLE_PRECISION,
     *            ISOUR,163,COMM,STATUS,IERR)
              END IF
              IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
                ISOUR=PRO2*NPCOL+0
                CALL MPI_SEND(TEMP(2),1,MPI_DOUBLE_PRECISION,
     *            ISOUR,164,COMM,IERR)
              END IF
              IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO2) THEN
                ISOUR=0
                CALL MPI_RECV(RHS_PAR(LOC2,1),1,
     *            MPI_DOUBLE_PRECISION,
     *            ISOUR,164,COMM,STATUS,IERR)
              END IF
            END DO
          END IF
        IF(MTYPE.eq.1) THEN
          CALL PDTRTRS('U','N','N',SIZE_ROOT-DEFICIENCY,1,
     *      A,1,1,DESCA_PAR,RHS_PAR,1,1,DESCB_PAR,IERR)
          CALL DMUMPS_173(SIZE_ROOT,1,
     *      LOCAL_M,1,MBLOCK,NBLOCK,IPIV,IPIV_SEQ,
     *      MYID,MASTER_ROOT,NPROW,NPCOL,COMM)
        ELSEIF(mtype.eq.0) then
          CALL DMUMPS_173(SIZE_ROOT,1,
     *      LOCAL_M,1,MBLOCK,NBLOCK,IPIV,IPIV_SEQ,
     *      MYID,MASTER_ROOT,NPROW,NPCOL,COMM)
          IF(myid.eq.master_root) then
            CALL DMUMPS_334(SIZE_ROOT,IPIV_SEQ,RHS_SEQ,WK)
          END IF
          CALL DMUMPS_290( MYID, SIZE_ROOT, 1, RHS_SEQ,
     *      LOCAL_M,1,MBLOCK,NBLOCK,RHS_PAR,MASTER_ROOT,
     *      NPROW,NPCOL,COMM)
          CALL PDTRTRS('U','T','N',SIZE_ROOT-DEFICIENCY,1,
     *      A,1,1,DESCA_PAR,RHS_PAR,1,1,DESCB_PAR,IERR)
          DO J=1,GIND
            K=GROW(J)
            IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
              C=GCOS(J)
              S=GSIN(J)
            END IF
            CALL DMUMPS_161(MYROW,MBLOCK,NPCOL,K,PRO)
            CALL DMUMPS_160(MYROW,MBLOCK,NPCOL,K,LOC)
            IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
              CALL MPI_SEND(RHS_PAR(LOC,1),1,
     *          MPI_DOUBLE_PRECISION,0,163,COMM,IERR)
            END IF
            IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
              ISOUR=PRO*NPCOL+0
              CALL MPI_RECV(TEMP(1),1,MPI_DOUBLE_PRECISION,
     *          ISOUR,163,COMM,STATUS,IERR)
            END IF
            CALL DMUMPS_161(MYROW,MBLOCK,NPCOL,K+1,PRO2)
            CALL DMUMPS_160(MYROW,MBLOCK,NPCOL,K+1,LOC2)
            IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO2) THEN
              CALL MPI_SEND(RHS_PAR(LOC2,1),1,
     *          MPI_DOUBLE_PRECISION,0,164,COMM,IERR)
            END IF
            IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
              ISOUR=PRO2*NPCOL+0
              CALL MPI_RECV(TEMP(2),1,MPI_DOUBLE_PRECISION,
     *          ISOUR,164,COMM,STATUS,IERR)
            END IF
            IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
              CALL DMUMPS_31(1,C,S,TEMP(1),1,TEMP(2),1)
            END IF
            IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
              ISOUR=PRO*NPCOL+0
              CALL MPI_SEND(TEMP(1),1,
     *          MPI_DOUBLE_PRECISION,
     *          ISOUR,163,COMM,IERR)
            END IF
            IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
              ISOUR=0
              CALL MPI_RECV(RHS_PAR(LOC,1),1,
     *          MPI_DOUBLE_PRECISION,
     *          ISOUR,163,COMM,STATUS,IERR)
            END IF
            IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
              ISOUR=PRO2*NPCOL+0
              CALL MPI_SEND(TEMP(2),1,
     *          MPI_DOUBLE_PRECISION,ISOUR,164,COMM,IERR)
            END IF
            IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO2) THEN
              ISOUR=0
              CALL MPI_RECV(RHS_PAR(LOC2,1),1,
     *          MPI_DOUBLE_PRECISION,
     *          ISOUR,164,COMM,STATUS,IERR)
            END IF
          END DO
          CALL DMUMPS_262( 'L','N', SIZE_ROOT, 1, SIZE_ROOT,
     *      A, 1, 1, DESCA_PAR, QR_TAU,
     *      RHS_PAR, 1, 1, DESCB_PAR,
     *      WK, LWK, INFO(1))
          IF ( INFO(1) .LT. 0 ) THEN
            WRITE(*,*) 'Error return from PDORMQR: INFO(1)=', INFO(1)
            CALL DMUMPS_ABORT()
          END IF
        ENDIF
      ELSEIF(QR.EQ.5.OR.QR.EQ.6) THEN
       IF ( MTYPE .eq. 1 ) THEN
        CALL PDGETRS('N',SIZE_ROOT-DEFICIENCY,1,A,1,1,DESCA_PAR,IPIV,
     *                RHS_PAR, 1, 1, DESCB_PAR, IERR )
       ELSE
        CALL PDGETRS('T',SIZE_ROOT-DEFICIENCY,1,A,1,1,DESCA_PAR,IPIV,
     *                RHS_PAR, 1, 1, DESCB_PAR, IERR )
       END IF
      ELSEIF(QR.EQ.7.OR.QR.EQ.8) THEN
       IF ( MTYPE .eq. 1 ) THEN
        CALL PDGETRS('N',SIZE_ROOT-DEFICIENCY,1,A,1,1,DESCA_PAR,IPIV,
     *                RHS_PAR, 1, 1, DESCB_PAR, IERR )
       ELSE
        CALL PDGETRS('T',SIZE_ROOT-DEFICIENCY,1,A,1,1,DESCA_PAR,IPIV,
     *                RHS_PAR, 1, 1, DESCB_PAR, IERR )
       END IF
      ELSE IF ( LDLT .eq. 0 .OR. LDLT .eq. 2 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          CALL PDGETRS('N',SIZE_ROOT,1,A,1,1,DESCA_PAR,IPIV,
     *      RHS_PAR,1,1,DESCB_PAR,IERR)
        ELSE
          CALL PDGETRS('T',SIZE_ROOT,1,A,1,1,DESCA_PAR,IPIV,
     *      RHS_PAR, 1, 1, DESCB_PAR,IERR)
        END IF
      ELSE
        CALL PDPOTRS( 'L', SIZE_ROOT, 1, A, 1, 1, DESCA_PAR,
     *    RHS_PAR, 1, 1, DESCB_PAR, IERR )
      END IF
      END IF
      IF ( IERR .LT. 0 ) THEN
        WRITE(*,*) ' Problem during solve of the root'
        CALL DMUMPS_ABORT()
      END IF
      IF(SEQ.EQ.1) THEN
        IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
          CALL DCOPY(NN,RHS_PAR,1,RHS_SEQ,1)
        ELSE
        END IF
      ELSE
        CALL DMUMPS_156( MYID, SIZE_ROOT, 1,
     *    RHS_SEQ, LOCAL_M, 1,
     *    MBLOCK, NBLOCK, RHS_PAR, MASTER_ROOT,
     *    NPROW, NPCOL, COMM )
      END IF
      IF(SEQ.EQ.1) THEN
        IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
          IF (MYID.eq.MASTER_ROOT.AND.
     *      (QR.EQ.1.OR.QR.EQ.2.OR.QR.EQ.3.OR.QR.EQ.4).and.
     *      mtype.eq.1) THEN
            CALL DMUMPS_333(NN,IPIV,RHS_SEQ,WK)
          END IF
        ELSE
        END IF
      ELSE
        IF (MYID.eq.MASTER_ROOT.AND.
     *    (QR.EQ.1.OR.QR.EQ.2.OR.QR.EQ.3.OR.QR.EQ.4).and.
     *    mtype.eq.1) THEN
          CALL DMUMPS_333( SIZE_ROOT, IPIV_SEQ, RHS_SEQ, WK )
        END IF
      END IF
      DEALLOCATE(IPIV_SEQ)
      RETURN
      END SUBROUTINE DMUMPS_352
#endif
      SUBROUTINE DMUMPS_286( NRHS, DESCA_PAR, DESCB_PAR,
     *  CNTXT_PAR,LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *  IPIV,LPIV,MASTER_ROOT,MYID,COMM,
     *  RHS_SEQ,SIZE_ROOT,A,INFO,MTYPE,LDLT )
      IMPLICIT NONE
      INTEGER NRHS, MTYPE
      INTEGER DESCA_PAR( 9 ), DESCB_PAR( 9 )
      INTEGER LOCAL_M, LOCAL_N, MBLOCK, NBLOCK
      INTEGER CNTXT_PAR, MASTER_ROOT, SIZE_ROOT
      INTEGER MYID, COMM
      INTEGER LPIV, IPIV( LPIV )
      INTEGER INFO(40), LDLT
      DOUBLE PRECISION RHS_SEQ( SIZE_ROOT *NRHS)
      DOUBLE PRECISION A( LOCAL_M, LOCAL_N )
      INCLUDE 'mpif.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER IERR, NPROW, NPCOL, MYROW, MYCOL
      INTEGER LOCAL_N_RHS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :,: ) ::RHS_PAR
      EXTERNAL NUMROC
      INTEGER  NUMROC
      INTEGER allocok
      CALL BLACS_GRIDINFO( CNTXT_PAR, NPROW, NPCOL, MYROW, MYCOL )
      LOCAL_N_RHS = NUMROC(NRHS, NBLOCK, MYCOL, 0, NPCOL)
      LOCAL_N_RHS = MAX(1,LOCAL_N_RHS)
      ALLOCATE(RHS_PAR(LOCAL_M, LOCAL_N_RHS),stat=allocok)
      IF (allocok > 0 ) THEN
        WRITE(*,*) ' Problem during solve of the root.'
        WRITE(*,*) ' Reduce number of right hand sides.'
        CALL DMUMPS_ABORT()
      ENDIF
      CALL DMUMPS_290( MYID, SIZE_ROOT, NRHS, RHS_SEQ,
     *      LOCAL_M, LOCAL_N_RHS,
     *      MBLOCK, NBLOCK, RHS_PAR, MASTER_ROOT,
     *      NPROW, NPCOL, COMM )
      IF ( LDLT .eq. 0 .OR. LDLT .eq. 2 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          CALL PDGETRS('N',SIZE_ROOT,NRHS,A,1,1,DESCA_PAR,IPIV,
     *      RHS_PAR,1,1,DESCB_PAR,IERR)
        ELSE
          CALL PDGETRS('T',SIZE_ROOT,NRHS,A,1,1,DESCA_PAR,IPIV,
     *      RHS_PAR, 1, 1, DESCB_PAR,IERR)
        END IF
      ELSE
        CALL PDPOTRS( 'L', SIZE_ROOT, NRHS, A, 1, 1, DESCA_PAR,
     *    RHS_PAR, 1, 1, DESCB_PAR, IERR )
      END IF
      IF ( IERR .LT. 0 ) THEN
        WRITE(*,*) ' Problem during solve of the root'
        CALL DMUMPS_ABORT()
      END IF
      CALL DMUMPS_156( MYID, SIZE_ROOT, NRHS,
     *    RHS_SEQ, LOCAL_M, LOCAL_N_RHS,
     *    MBLOCK, NBLOCK, RHS_PAR, MASTER_ROOT,
     *    NPROW, NPCOL, COMM )
      DEALLOCATE(RHS_PAR)
      RETURN
      END SUBROUTINE DMUMPS_286
