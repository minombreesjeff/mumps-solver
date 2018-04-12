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
      SUBROUTINE SMUMPS_324(A, NFRONT, NPIV, NBROW, K59)
      IMPLICIT NONE
      INTEGER NFRONT, NPIV, NBROW, K59
      REAL A(NFRONT*(NBROW+NPIV))
      INTEGER IOLD, INEW, I , J
      IF ((NPIV.EQ.0).OR.(NFRONT.EQ.NPIV)) GOTO 500
      IOLD = NFRONT + 1
      INEW = NPIV + 1
      IF ( NPIV < K59 ) THEN
        DO I = 1, NBROW+NPIV-1
          DO J = 0, NPIV-1
            A(INEW+J) = A(IOLD+J)
          ENDDO
          INEW = INEW + NPIV
          IOLD = IOLD + NFRONT
        ENDDO
      ELSE
        DO I = 1, NPIV - 2
          DO J = 0, I+1
            A( INEW + J ) = A( IOLD + J )
          END DO
          INEW = INEW + NPIV
          IOLD = IOLD + NFRONT
        END DO
        I = NPIV - 1
        DO J = 0, I
            A( INEW + J ) = A( IOLD + J )
        END DO
        INEW = INEW + NPIV
        IOLD = IOLD + NFRONT
        DO I = NPIV, NBROW+NPIV-1
          DO J = 0, NPIV - 1
            A( INEW + J ) = A( IOLD + J )
          END DO
          INEW = INEW + NPIV
          IOLD = IOLD + NFRONT
        ENDDO
      END IF
 500  RETURN
      END SUBROUTINE SMUMPS_324
      SUBROUTINE SMUMPS_325(A, NPIV, NBROW)
      IMPLICIT NONE
      INTEGER NPIV, NBROW
      REAL A((NBROW+NPIV)*(NBROW+NPIV))
      INTEGER IOLD, INEW, I , J, LDA
      IF ((NPIV.EQ.0).OR.(NBROW.EQ.0)) GOTO 500
      LDA  = NBROW+NPIV
      IOLD = LDA + 1
      INEW = NPIV + 1
      DO I = 1, NBROW+NPIV-1
        DO J = 0, NPIV-1
         A(INEW+J) = A(IOLD+J)
        ENDDO
        INEW = INEW + NPIV
        IOLD = IOLD + LDA
      ENDDO
 500  RETURN
      END SUBROUTINE SMUMPS_325
      SUBROUTINE SMUMPS_308( N, PERM, IW, LIW )
      IMPLICIT NONE
      INTEGER N, LIW
      INTEGER PERM( N ), IW( LIW )
      INTEGER I, SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, LIW - 1
          IF ( PERM( IW( I ) ) .GT. PERM( IW( I + 1 ) ) ) THEN
            DONE = .FALSE.
            SWAP  = IW( I + 1 )
            IW( I + 1 ) = IW( I )
            IW( I ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_308
      SUBROUTINE SMUMPS_309( N, NASS1, PERM, ITLOC,
     *                             SMALL, LSMALL,
     *                             LARGE, LLARGE,
     *                             MERGE, LMERGE )
      IMPLICIT NONE
      INTEGER N, NASS1, LSMALL, LLARGE, LMERGE
      INTEGER PERM( N ), ITLOC( N )
      INTEGER SMALL(LSMALL), LARGE(LLARGE), MERGE(LMERGE)
      INTEGER PSMALL, PLARGE, PMERGE, VSMALL, VLARGE, VMERGE
      PSMALL = 1
      PLARGE = 1
      PMERGE = 1
      DO WHILE ( PSMALL .LE. LSMALL .or. PLARGE.LE. LLARGE )
        IF ( PSMALL .GT. LSMALL ) THEN
          VMERGE = LARGE( PLARGE )
          PLARGE = PLARGE + 1
        ELSE IF ( PLARGE .GT. LLARGE ) THEN
          VMERGE = SMALL( PSMALL )
          PSMALL = PSMALL + 1
        ELSE
          VSMALL = SMALL( PSMALL )
          VLARGE = LARGE( PLARGE )
          IF ( PERM( VSMALL ) .LT. PERM( VLARGE ) ) THEN
            VMERGE = VSMALL
            PSMALL   = PSMALL + 1
          ELSE
            VMERGE = VLARGE
            PLARGE   = PLARGE + 1
          END IF
        END IF
        MERGE( PMERGE ) = VMERGE
        ITLOC( VMERGE ) = PMERGE + NASS1
        PMERGE = PMERGE + 1
      END DO
      PMERGE = PMERGE - 1
      RETURN
      END SUBROUTINE SMUMPS_309
      SUBROUTINE SMUMPS_88( MYID, N, PERM,
     *                             ITLOC,
     *                             IW, LIW, NASS1 )
      IMPLICIT NONE
      INTEGER MYID, N, LIW, NASS1
      INTEGER PERM( N ), ITLOC( N ), IW( LIW )
      INTEGER I, J
      IF ( LIW .eq. 0 ) RETURN
      IF ( NASS1 .eq. 0 ) THEN
      END IF
      DO I = 1, LIW - 1
        IF ( PERM( IW( I ) ) .GE. PERM( IW( I + 1 ) ) ) THEN
          WRITE(*,*) MYID,':Error in CHECK_ORDER_LIST: not sorted'
          WRITE(*,*) MYID,':List = ',IW(1:LIW)
          WRITE(*,*) MYID,':Permuted list =', (PERM(IW(J)),J=1,LIW)
          CALL SMUMPS_ABORT()
        END IF
        IF (ITLOC( IW( I + 1 ) ) .NE. ITLOC( IW( I ) ) + 1 ) THEN
          WRITE(*,*) MYID,':Error in MUMP_CHECK_ORDER_LIST: bad ITLOC'
          WRITE(*,*) MYID,':ITLOC( IW ) =',(ITLOC(IW(J)),J=1,LIW)
          CALL SMUMPS_ABORT() 
        END IF
      END DO
      IF ( ITLOC( IW( 1 ) )  .NE. NASS1 + 1 ) THEN
        WRITE(*,*) MYID,': Bad First ITLOC ITLOC(IW(1)),NASS1+1=',
     *              ITLOC( IW( 1 ) ), NASS1 + 1
        CALL SMUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE SMUMPS_88
      SUBROUTINE SMUMPS_181( N, NBROOT, NBLEAF, POOL, 
     &          LPOOL, PERM, FILS, FRERE, NSTK )
      IMPLICIT NONE
      INTEGER N, NBROOT, NBLEAF, LPOOL
      INTEGER POOL( LPOOL ), PERM( N ), FILS( N ),
     *        FRERE( N ), NSTK( N )
      INTEGER IPERM, INODE, IN
      IPERM = 1
      DO WHILE ( NBLEAF .NE. 0 )
        INODE = POOL( NBLEAF )
        NBLEAF = NBLEAF - 1
        IN = INODE
        DO WHILE ( IN .GT. 0 )
          PERM ( IN ) = IPERM
          IPERM = IPERM + 1
          IN = FILS( IN )
        END DO
        IN = FRERE( INODE )
        DO WHILE( IN .GT. 0 )
          IN = FRERE( IN )
        END DO
        IF ( IN .eq. 0 ) THEN
          NBROOT = NBROOT - 1
        ELSE
          IN = - IN
          NSTK( IN ) = NSTK( IN ) - 1
          IF ( NSTK( IN ) .eq. 0 ) THEN
            NBLEAF = NBLEAF + 1
            POOL( NBLEAF ) = IN
          END IF
        END IF
      END DO
      END SUBROUTINE SMUMPS_181
      SUBROUTINE SMUMPS_81(MYID, INODE, N, IOLDPS,
     *           HF, NFRONT, NFRONT_EFF,
     *           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *           IFSON, STEP, PIMASTER, PTRAIW, IW, LIW, 
     *           INTARR, ITLOC, FILS, FRERE, 
     *           SON_LEVEL2, NIV1, NBPROCFILS, KEEP,KEEP8, IFLAG)
      IMPLICIT NONE
      INTEGER INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     *        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG
      INTEGER KEEP(500), NFRONT_EFF
      INTEGER*8 KEEP8(150)
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PTRAIW(N), IW(LIW), 
     *        ITLOC(N), FILS(N), FRERE(KEEP(28)),
     *        NBPROCFILS(KEEP(28))
      INTEGER INTARR(MAX(1,KEEP(14)))
      LOGICAL SON_LEVEL2, NIV1
      INTEGER NEWEL, INEW, IOLDP2, INEW1,
     *        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     *        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     *        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     *        I, K
      LOGICAL LEVEL1
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      SON_LEVEL2 = .FALSE.
      IOLDP2 = IOLDPS + HF - 1   ! HF includes XSIZE shift
      NEWEL = IOLDP2 + NASS1
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF 
      INEW1 = 1
   50 J1 = PTRAIW(IN) + 2
      JT1 = INTARR(J1)
      INTARR(J1) = INEW1
      ITLOC(JT1) = INEW1
      IW(INEW) = JT1
      INEW = INEW + 1
      INEW1 = INEW1 + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 50
      IF (NUMSTK .NE. 0) THEN
        NTOTFS = NUMORG
        ISON = IFSON
        ICT11 = IOLDP2 + NFRONT
        DO 100 IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK = IW(J2+XSIZE)
          NELIM = IW(J2 + 1+XSIZE)
          NPIVS = IW(J2 + 3+XSIZE)
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON = IW(J2 + 5+XSIZE)
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1    = NSLSON.EQ.0
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          ITRANS = NROWS
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
          ELSE
           IF (LEVEL1) THEN
            NBPROCFILS(STEP(ISON)) = 1
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     *                               NBPROCFILS(STEP(ISON))
          ENDIF
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+XSIZE)
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6 + XSIZE
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          IF (NELIM .EQ. 0) GOTO 70
          DO 60 JJ = J1, J3
            NTOTFS = NTOTFS + 1
            JT1 = IW(JJ)
            IW(ICT11 + NTOTFS) = JT1
            ITLOC(JT1) = NTOTFS
            IW(JJ) = NTOTFS
            IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
   60     CONTINUE
   70     J1 = J3 + 1
          IF (NASS1 .NE. NFRONT) THEN
            DO 80 JJ = J1, J2
              J = IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN
                NEWEL = NEWEL + 1
                NFRONT_EFF = NFRONT_EFF + 1
                IW(NEWEL) = J
                IW(JJ) = NFRONT_EFF
                ITLOC(J) = NFRONT_EFF
              ELSE
                IW(JJ) = ITLOC(J)
              ENDIF
   80       CONTINUE
          ELSE
            DO 90 JJ = J1, J2
              IW(JJ) = ITLOC(IW(JJ))
   90       CONTINUE
          ENDIF
          ISON = FRERE(STEP(ISON))
  100   CONTINUE
      ENDIF
      IBROT = INODE
      DO 120 IORG = 1, NUMORG
        J1 = PTRAIW(IBROT) + 2
        IBROT = FILS(IBROT)
        J2 = J1 + INTARR(J1 - 2) - INTARR(J1 - 1)
        J1 = J1 + 1
        IF (J1 .LE. J2) THEN
          DO 110 JJ = J1, J2
            J = INTARR(JJ)
            IF (ITLOC(J) .EQ. 0) THEN
              NEWEL = NEWEL + 1
              NFRONT_EFF = NFRONT_EFF + 1
              IW(NEWEL) = J
              INTARR(JJ) = NFRONT_EFF
              ITLOC(J) = NFRONT_EFF
            ELSE
              INTARR(JJ) = ITLOC(J)
            ENDIF
  110     CONTINUE
        ENDIF
  120 CONTINUE
      IF (NFRONT.NE.NFRONT_EFF) THEN
       IF (NUMORG.EQ.NASS1) THEN
        IP1 = IOLDPS + HF 
        IP2 = IOLDPS + HF + NFRONT_EFF - 1 
        DO I = IP1, IP2
         IW(I + NFRONT_EFF) = IW(I)
        ENDDO
       ELSE
        IP1 = IOLDPS + NFRONT + HF + NUMORG 
        IP2 = IOLDPS + HF + NFRONT_EFF + NUMORG 
        IW(IP2:IP2+NASS-1) = IW(IP1:IP1+NASS-1)
        IP1 = IOLDPS + NASS1 + HF 
        IP2 = IOLDPS + HF + NFRONT - 1 
        DO I = IP1, IP2
         IW(I + NFRONT_EFF) = IW(I)
        ENDDO
        IP1 = IOLDPS + HF 
        IP2 = IOLDPS + HF + NUMORG - 1 
        DO I = IP1, IP2
          IW(I + NFRONT_EFF) = IW(I)
        ENDDO
       ENDIF
      ELSE
       IP1 = IOLDPS + NASS1 + HF 
       IP2 = IOLDPS + HF + NFRONT - 1 
       DO I = IP1, IP2
        IW(I + NFRONT) = IW(I)
       ENDDO
       IP1 = IOLDPS + HF 
       IP2 = IOLDPS + HF + NUMORG - 1 
       DO I = IP1, IP2
         IW(I + NFRONT) = IW(I)
       ENDDO
      ENDIF
      K1 = IOLDPS + HF + NUMORG 
      K2 = K1 + NFRONT_EFF - 1 + NASS
      DO 150 K = K1, K2
        I = IW(K)
        ITLOC(I) = 0
  150 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_81
      SUBROUTINE SMUMPS_86(MYID, INODE, N, IOLDPS,
     *           HF, NFRONT, NFRONT_EFF, PERM, 
     *           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *           IFSON, STEP, PIMASTER, PTRAIW, IW, LIW, 
     *           INTARR, ITLOC, FILS, FRERE_STEPS, 
     *           SON_LEVEL2, NIV1, NBPROCFILS, KEEP,KEEP8, IFLAG)
      IMPLICIT NONE
      INTEGER INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     *        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG, 
     *        NFRONT_EFF
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PTRAIW(N), IW(LIW), 
     *        ITLOC(N), FILS(N), FRERE_STEPS(KEEP(28)),
     *        NBPROCFILS(KEEP(28)), PERM(N)
      INTEGER INTARR(MAX(1,KEEP(14)))
      LOGICAL SON_LEVEL2, NIV1
      INTEGER NEWEL, IOLDP2, INEW, INEW1,
     *        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     *        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     *        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     *        I, K, JDEBROW, ILOC, NEWEL_SAVE, NEWEL1_SAVE,
     *        LAST_J_ASS, JMIN, MIN_PERM
      LOGICAL LEVEL1
       INTEGER  XSIZE
       PARAMETER (XSIZE=3)
      INTEGER allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTTRI, PTLAST
      ALLOCATE(PTTRI(NUMSTK+1), stat=allocok)
      IF (allocok .GT. 0) THEN
       IFLAG = -13
       GOTO 800
      ENDIF
      ALLOCATE(PTLAST(NUMSTK+1), stat=allocok)
      IF (allocok .GT. 0) THEN
       IFLAG = -13
       GOTO 800
      ENDIF
      SON_LEVEL2 = .FALSE.
      IOLDP2     = IOLDPS + HF - 1
      ICT11      = IOLDP2 + NFRONT
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF
      INEW1 = 1
      JDEBROW = PTRAIW(INODE)+3
      PTTRI(NUMSTK+1)  = JDEBROW
      PTLAST(NUMSTK+1) = JDEBROW + INTARR(JDEBROW-3) - 1
   50 J1 = PTRAIW(IN) + 2
      JT1 = INTARR(J1)
      INTARR(J1) = INEW1
      ITLOC(JT1) = INEW1
      IW(INEW)         = JT1
      IW (INEW+NFRONT) = JT1
      INEW = INEW + 1
      INEW1 = INEW1 + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 50
      NTOTFS = NUMORG
      IF (NUMSTK .NE. 0) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK   = IW(J2    +XSIZE)
          NELIM  = IW(J2 + 1+XSIZE)
          NPIVS  = IW(J2 + 3+XSIZE)
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5+XSIZE)
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1    = NSLSON.EQ.0
          NCOLS  = NPIVS + LSTK
          NROWS  = NCOLS
          ITRANS = NROWS
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
          ELSE
           IF (LEVEL1) THEN
            NBPROCFILS(STEP(ISON)) = 1
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     *                               NBPROCFILS(STEP(ISON))
          ENDIF
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+XSIZE)
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6 + XSIZE
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          IF (NELIM .NE. 0) THEN
            DO JJ = J1, J3
              NTOTFS = NTOTFS + 1
              JT1 = IW(JJ)
              IW(ICT11 + NTOTFS) = JT1
              ITLOC(JT1) = NTOTFS
              IW(JJ) = NTOTFS
              IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
            ENDDO
          ENDIF
          PTTRI(IELL)  = J2+1
          PTLAST(IELL) = J2
          J1 = J3 + 1
          IF (NASS1 .NE. NFRONT) THEN
            DO JJ = J1, J2
              J = IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN 
                PTTRI(IELL) = JJ
                EXIT
              ENDIF
            ENDDO
          ELSE
            DO JJ = J1, J2
              IW(JJ) = ITLOC(IW(JJ))
            ENDDO
          ENDIF
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IF (NFRONT.EQ.NASS1) GOTO 600
 199  CONTINUE
      IF ( PTTRI( NUMSTK + 1 ) .LE. PTLAST( NUMSTK + 1 ) ) THEN
      IF ( ITLOC( INTARR( PTTRI( NUMSTK + 1 ) ) ) .NE. 0 ) THEN
       PTTRI( NUMSTK + 1 ) = PTTRI( NUMSTK + 1 ) + 1
       GOTO 199
      END IF
      END IF
      MIN_PERM = N + 1
      DO IELL = 1, NUMSTK 
        ILOC = PTTRI( IELL )
        IF ( ILOC .LE. PTLAST( IELL ) ) THEN 
         IF ( PERM( IW( ILOC ) ) .LT. MIN_PERM ) THEN
           JMIN     = IW( ILOC )
           MIN_PERM = PERM( JMIN )
         END IF
        END IF
      END DO
      IELL = NUMSTK + 1
      ILOC =  PTTRI( IELL )
      IF ( ILOC .LE. PTLAST( IELL ) ) THEN
        IF ( PERM( INTARR( ILOC ) ) .LT. MIN_PERM ) THEN
         JMIN        = INTARR( ILOC )
         MIN_PERM = PERM( JMIN )
        END IF
      END IF
      NEWEL = IOLDP2 + NASS1 + NFRONT
      DO WHILE ( MIN_PERM .NE. N + 1 )
          NEWEL  = NEWEL + 1
          NFRONT_EFF = NFRONT_EFF + 1
          IW( NEWEL ) = JMIN
          ITLOC( JMIN ) = NFRONT_EFF
          LAST_J_ASS = JMIN
          MIN_PERM = N + 1
          DO IELL = 1,  NUMSTK
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
              IF ( IW( PTTRI( IELL ) ) .eq. LAST_J_ASS )
     *        PTTRI( IELL ) = PTTRI( IELL ) + 1
            ENDIF
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
             IF ( PERM(IW( PTTRI( IELL )) ) .LT. MIN_PERM ) THEN
                JMIN        = IW( PTTRI( IELL ) )
                MIN_PERM = PERM( JMIN )
             END IF
            END IF
          END DO
          IELL = NUMSTK + 1
 145      CONTINUE
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
            IF ( INTARR( PTTRI( IELL ) ) .eq. LAST_J_ASS ) THEN
              PTTRI( IELL ) = PTTRI( IELL ) + 1 
              GOTO 145
            END IF
          END IF
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
            IF (PERM(INTARR( PTTRI(IELL) )) .LT. MIN_PERM) THEN
              JMIN        = INTARR( PTTRI(IELL) )
              MIN_PERM = PERM( JMIN )
            END IF
          END IF
      END DO
      NEWEL_SAVE  = NEWEL
      NEWEL1_SAVE = NFRONT_EFF
      IF (NEWEL1_SAVE.LT.NFRONT) THEN 
       IBROT = FILS( INODE )
       DO IORG = 2, NUMORG
         J1    = PTRAIW(IBROT) + 2
         IBROT = FILS( IBROT )
         J2    = J1 + INTARR(J1 - 2)
         J1    = J1 + 1
         DO JJ = J1, J2
           J     = INTARR( JJ )
           IF ( ITLOC( J ) .eq. 0 ) THEN
            NEWEL  = NEWEL + 1
            NFRONT_EFF = NFRONT_EFF + 1
            IW( NEWEL ) = J
            ITLOC( J ) = NFRONT_EFF
           END IF
         ENDDO
       ENDDO
      ENDIF
      IF ( NEWEL1_SAVE .eq. NFRONT_EFF ) THEN
       IW( IOLDP2+NASS1+1 : IOLDP2+NFRONT_EFF ) = 
     *  IW( ICT11+NASS1+1 : ICT11+NFRONT_EFF ) 
      ELSE
        CALL SMUMPS_308( N, PERM, 
     *           IW( NEWEL_SAVE + 1 ), NFRONT_EFF - NEWEL1_SAVE )
          CALL SMUMPS_309( N, NASS1, PERM, ITLOC,
     *               IW( NEWEL_SAVE + 1), NFRONT_EFF - NEWEL1_SAVE,
     *               IW( ICT11  + NASS1 + 1 ), NEWEL1_SAVE - NASS1,
     *               IW( IOLDP2 + NASS1 + 1 ), NFRONT_EFF - NASS1 )
      END IF
      IF ( NFRONT_EFF .NE. NEWEL1_SAVE ) THEN
        IP1 = IOLDPS + NASS1 + HF    ! XSIZE is included in HF
        IP2 = IOLDPS + HF + NFRONT - 1
        DO I = IP1, IP2
          IW(I + NFRONT) = IW(I)
        ENDDO
      END IF
      IF (NFRONT.NE.NFRONT_EFF) THEN
        IP1 = IOLDPS + NFRONT + HF 
        IP2 = IOLDPS + NFRONT_EFF + HF
        IW (IP2:IP2+NFRONT_EFF) = IW (IP1:IP1+NFRONT_EFF)
      ENDIF
  600 CONTINUE
      IF ((NUMSTK .NE. 0).AND.(NFRONT.GT.NASS1)) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK = IW(J2+XSIZE)
          NELIM = IW(J2 + 1+XSIZE)
          NPIVS = IW(J2 + 3+XSIZE)
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5+XSIZE)
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+XSIZE)
          ENDIF
          HS = NSLSON + 6 +XSIZE
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          J1 = J3 + 1
          DO JJ = J1, J2
              J = IW(JJ)
                IW(JJ) = ITLOC(J)
          ENDDO
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IBROT = INODE
      DO IORG = 1, NUMORG
        J1 = PTRAIW(IBROT) + 2
        IBROT = FILS(IBROT)
        J2 = J1 + INTARR(J1 - 2) - INTARR(J1 - 1)
        J1 = J1 + 1
        IF (J1 .LE. J2) THEN
          DO JJ = J1, J2
            J = INTARR(JJ)
            INTARR(JJ) = ITLOC(J)
          ENDDO
        ENDIF
      ENDDO
        K1 = IOLDPS + HF + NUMORG   ! XSIZE included in HF
        K2 = K1 + NFRONT_EFF - 1 + NASS
        DO K = K1, K2
          I = IW(K)
          ITLOC(I) = 0
        ENDDO
  800 CONTINUE
      IF (ALLOCATED(PTTRI)) DEALLOCATE(PTTRI)
      IF (ALLOCATED(PTLAST)) DEALLOCATE(PTLAST)
      RETURN
      END SUBROUTINE SMUMPS_86
      SUBROUTINE SMUMPS_140( N, INODE, IW, LIW, A, LA,
     *                           IOLDPS, POSELT, IFLAG,
     *                           UU, NNEG, NPVW,
     *                           KEEP,KEEP8,
     *                           MYID, SEUIL, AVOID_DELAYED, ETATASS,
     *     DKEEP,PIVNUL_LIST,LPN_LIST )
      IMPLICIT NONE
      INTEGER N, INODE, LIW, LA, IFLAG, NNEG, NPVW
      INTEGER MYID, SLAVEF, IOLDPS, POSELT
      INTEGER KEEP( 500 )
      INTEGER*8 KEEP8(150)
      REAL UU, SEUIL
      INTEGER POSFAC, NIRBDU, COMP
      REAL A( LA )
      INTEGER IW( LIW )
      LOGICAL AVOID_DELAYED
      INTEGER ETATASS
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(30)
      INTEGER INOPV, IFINB, NFRONT, NPIV, LAELL, NBTLKJ,IBEG_BLOCK
      INTEGER NASS, NEL1, LDA
      REAL UUTEMP
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      EXTERNAL SMUMPS_330, SMUMPS_222, SMUMPS_234, 
     *         SMUMPS_230, SMUMPS_226, 
     *         SMUMPS_237
      INTEGER  SMUMPS_330
      LOGICAL STATICMODE
      REAL SEUIL_LOC
      INTEGER PIVSIZ,IWPOS
      INOPV = 0
      SEUIL_LOC = SEUIL
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
        STATICMODE = .TRUE.
        UUTEMP=UU
        SEUIL_LOC = MAX(SEUIL,EPSILON(SEUIL))
      ELSE
        UUTEMP=UU
      ENDIF
      IBEG_BLOCK = 1
      NFRONT = IW(IOLDPS+XSIZE)
      LDA    = NFRONT
      NASS   = IABS(IW(IOLDPS+2+XSIZE))
      IF (NASS .GT. KEEP(3)) THEN
        NBTLKJ = MIN( KEEP(6), NASS )
      ELSE
        NBTLKJ = MIN( KEEP(5), NASS )
      ENDIF
      IW(IOLDPS+3+XSIZE) = MIN0(NASS,NBTLKJ)
      UUTEMP = UU
 50   CONTINUE
      CALL SMUMPS_222(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *                INOPV, NNEG, IFLAG,IOLDPS,POSELT,UUTEMP,
     *                SEUIL_LOC,KEEP,KEEP8,PIVSIZ,
     *     DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
      IF(KEEP(109).GT. 0) THEN
         IF(PIVNUL_LIST(KEEP(109)).EQ.-1) THEN
            IWPOS = IOLDPS+IW(IOLDPS+1+XSIZE)+6
            PIVNUL_LIST(KEEP(109)) = IW(IWPOS+XSIZE)
         ENDIF
      ENDIF
      IF (IFLAG.LT.0) GOTO 500
      IF (INOPV.EQ.1) THEN
         IF(STATICMODE) THEN
            INOPV = -1
            GOTO 50
         ENDIF
         CALL SMUMPS_237(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *        LDA, IOLDPS,POSELT, KEEP,KEEP8, UUTEMP, ETATASS )
         GOTO 500
      END IF
      IF (INOPV.EQ.2) THEN
         CALL SMUMPS_234(IBEG_BLOCK,
     *            NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            LDA, IOLDPS,POSELT,NBTLKJ,KEEP(4),UUTEMP,
     *            KEEP,KEEP8)
         GOTO 50
      ENDIF
      NPVW = NPVW + PIVSIZ
      IF (NASS.LE.1) THEN
       CALL SMUMPS_230(NFRONT,N,INODE,IW,LIW,A,LA,
     *                 IOLDPS,POSELT)
       IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
       GO TO 500
      ENDIF
       CALL SMUMPS_226(IBEG_BLOCK,
     *             NFRONT, NASS, N,INODE,IW,LIW,A,LA,
     *             LDA, UUTEMP, IOLDPS,POSELT,IFINB,
     *             NBTLKJ,PIVSIZ)
       IF(PIVSIZ .EQ. 2) THEN
          IWPOS = IOLDPS+IW(IOLDPS+1+XSIZE)+6
          IW(IWPOS+NFRONT+XSIZE) = -IW(IWPOS+NFRONT+XSIZE)
       ENDIF
       IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + PIVSIZ
       IF (IFINB.EQ.0) GOTO 50
       NPIV   = IW(IOLDPS+1+XSIZE)
       NEL1   = NASS - NPIV
       CALL SMUMPS_234(IBEG_BLOCK,
     *            NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            LDA, IOLDPS,POSELT,NBTLKJ,KEEP(4),UUTEMP,
     *            KEEP,KEEP8)
       IF (IFINB.EQ.-1) THEN 
         CALL SMUMPS_237(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *         LDA, IOLDPS,POSELT, KEEP,KEEP8, UUTEMP, ETATASS)
         GOTO 500
       ENDIF
      GO TO 50
 500  CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_140
      SUBROUTINE SMUMPS_222 
     *   (NFRONT,NASS,N,INODE,IW,LIW,
     *    A,LA, INOPV,
     *    NNEG,
     *    IFLAG,IOLDPS,POSELT,UU, SEUIL,KEEP,KEEP8,PIVSIZ,
     *     DKEEP,PIVNUL_LIST,LPN_LIST)
#if defined (PROFILE_BLAS_ASS_G)
      USE SMUMPS_LOAD
#endif
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFLAG,INOPV,
     *        IOLDPS, POSELT, NNEG
      INTEGER PIVSIZ,LPIV
      REAL A(LA) 
      REAL UU, UULOC, SEUIL
      INTEGER IW(LIW)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(30)
      include 'mpif.h'
      INTEGER POSPV1,POSPV2,OFFDAG,APOSJ
      INTEGER JMAX
      REAL RMAX,AMAX,TMAX,SWOP,TOL
      REAL DELTA,MAXPIV
      REAL PIVNUL,FIXA
      REAL PIVOT,DETPIV
      PARAMETER(TOL = 1.0E-20)
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INTEGER APOS, LDA
      REAL ZERO,ONE
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,J1,JJ,J2,K
      INTRINSIC MAX
      DATA ZERO /0.0E0/
      DATA ONE /1.0E0/
      PIVNUL = DKEEP(1)
      FIXA = DKEEP(2)
      LDA   = NFRONT
      UULOC = UU
      PIVSIZ = 1
      NPIV    = IW(IOLDPS+1+XSIZE)
      NPIVP1  = NPIV + 1
      NASSW   = IABS(IW(IOLDPS+3+XSIZE))
      IF(INOPV .EQ. -1) THEN
         APOS = POSELT + LDA*NPIV + NPIV
         POSPV1 = APOS
         IF(ABS(A(APOS)).LT.SEUIL) THEN
            IF(REAL(A(APOS)) .GE. ZERO) THEN
               A(APOS) = SEUIL
            ELSE
               A(APOS) = -SEUIL
           NNEG = NNEG+1
            ENDIF
            KEEP(98) = KEEP(98)+1
         ENDIF
         GO TO 420
      ENDIF
      INOPV   = 0
      DO 460 IPIV=NPIVP1,NASSW
         APOS = POSELT + LDA*(IPIV-1) + NPIV
         POSPV1 = APOS + IPIV - NPIVP1
         PIVOT = A(POSPV1)
         IF (UULOC.EQ.ZERO) THEN 
            IF (ABS(A(APOS)).EQ.ZERO) GO TO 630
              IF (A(APOS).LT.ZERO) NNEG = NNEG+1
            GO TO 420
         ENDIF
         AMAX = ZERO
         JMAX = 0
         J1 = APOS
         J2 = POSPV1 - 1
         DO JJ=J1,J2
            IF(ABS(A(JJ)) .GT. AMAX) THEN
               AMAX = ABS(A(JJ))
               JMAX = IPIV - (POSPV1-JJ)
            ENDIF
         ENDDO
         J1 = POSPV1 + LDA
         J2 = NASSW - IPIV
         DO JJ=1,J2
            IF(ABS(A(J1)) .GT. AMAX) THEN
               AMAX = ABS(A(J1))
               JMAX = IPIV + JJ
            ENDIF
            J1 = J1 + LDA
         ENDDO
         RMAX = ZERO
         J2 = NFRONT - NASSW 
         DO JJ=1,J2
            RMAX = MAX(ABS(A(J1)),RMAX)
            J1 = J1 + LDA
         ENDDO
         IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.PIVNUL) THEN
            KEEP(109) = KEEP(109)+1
            PIVNUL_LIST(KEEP(109)) = -1
            IF(FIXA.GT.ZERO) THEN
               IF(REAL(PIVOT) .GE. ZERO) THEN
                  A(POSPV1) = FIXA
               ELSE
                  A(POSPV1) = -FIXA
               ENDIF
            ELSE
               J1 = APOS
               J2 = POSPV1 - 1
               DO JJ=J1,J2
                  A(JJ) = ZERO
               ENDDO
               J1 = POSPV1 + LDA
               J2 = NASSW - IPIV
               DO JJ=1,J2
                  A(J1) = ZERO
                  J1 = J1 + LDA
               ENDDO
               J2 = NFRONT - NASSW 
               DO JJ=1,J2
                  A(J1) = ZERO
                  J1 = J1 + LDA
               ENDDO
               A(POSPV1) = ONE
            ENDIF
            PIVOT = A(POSPV1)
            GO TO 415
         ENDIF
         IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.TOL) THEN
            IF(SEUIL .GT. EPSILON(SEUIL)) THEN
               IF(REAL(PIVOT) .GE. ZERO) THEN
                  A(POSPV1) = SEUIL
               ELSE
                  A(POSPV1) = -SEUIL
           NNEG = NNEG+1
               ENDIF
               PIVOT = A(POSPV1)
               WRITE(*,*) 'WARNING matrix may be singular'
               KEEP(98) = KEEP(98)+1
               GO TO 415
            ENDIF
         ENDIF
         IF (MAX(AMAX,ABS(PIVOT)).LE.TOL) GO TO 460
            IF (ABS(PIVOT).GT.UULOC*MAX(RMAX,AMAX)) THEN
               A(POSPV1) = PIVOT
            IF (A(POSPV1).LT.ZERO) NNEG = NNEG+1
               GO TO 415
            END IF
            IF (AMAX.LE.TOL) GO TO 460
            IF (RMAX.LT.AMAX) THEN
               J1 = APOS
               J2 = POSPV1 - 1
               DO JJ=J1,J2
                  IF(POSPV1-JJ .NE. IPIV-JMAX) THEN
                     RMAX = MAX(RMAX,ABS(A(JJ)))
                  ENDIF
               ENDDO
               J1 = POSPV1 + LDA
               J2 = NASS - IPIV
               DO JJ=1,J2
                  IF(IPIV+JJ .NE. JMAX) THEN
                     RMAX = MAX(ABS(A(J1)),RMAX)
                  ENDIF
                  J1 = J1 + LDA
               ENDDO
           ENDIF
           APOSJ = POSELT + (JMAX-1)*LDA + NPIV
           POSPV2 = APOSJ + JMAX - NPIVP1
           IF (IPIV.LT.JMAX) THEN
              OFFDAG = APOSJ + IPIV - NPIVP1
           ELSE
              OFFDAG = APOS + JMAX - NPIVP1
           END IF
           TMAX = ZERO
           IF(JMAX .LT. IPIV) THEN
              JJ = POSPV2
              DO K = 1, NFRONT-JMAX
                 JJ = JJ+NFRONT
                 IF (JMAX+K.NE.IPIV) THEN
                    TMAX=MAX(TMAX,ABS(A(JJ)))
                 ENDIF
              ENDDO
              DO K =  APOSJ, POSPV2-1
                 TMAX = MAX(TMAX,ABS(A(K)))
              ENDDO
           ELSE
              JJ = POSPV2
              DO K = 1, NFRONT-JMAX
                 JJ = JJ+NFRONT
                 TMAX=MAX(TMAX,ABS(A(JJ)))
              ENDDO
              DO K =  APOSJ, POSPV2-1
                 IF (K.NE.OFFDAG) THEN
                    TMAX = MAX(TMAX,ABS(A(K)))
                 ENDIF
              ENDDO
           ENDIF
           DETPIV = A(POSPV1)*A(POSPV2) - A(OFFDAG)**2
           MAXPIV = MAX(ABS(A(POSPV1)),ABS(A(POSPV2)))
           IF (MAXPIV.EQ.ZERO) MAXPIV = ONE
           IF (ABS(DETPIV)/MAXPIV.LE.TOL) GO TO 460
           IF ((ABS(A(POSPV2))*RMAX+AMAX*TMAX)*UULOC.GT.
     +          ABS(DETPIV)) GO TO 460
           IF ((ABS(A(POSPV1))*TMAX+AMAX*RMAX)*UULOC.GT.
     +          ABS(DETPIV)) GO TO 460
           PIVSIZ = 2
           KEEP(103) = KEEP(103)+1
           IF(DETPIV .LT. ZERO) THEN
             NNEG = NNEG+1
           ELSE IF(A(POSPV2) .LT. ZERO) THEN
             NNEG = NNEG+2
           ENDIF
 415       CONTINUE
           IF (PIVSIZ.EQ.2) THEN
              LPIV = MIN(IPIV,JMAX)
           ELSE
              LPIV = IPIV              
           ENDIF
           DO K=1,PIVSIZ
              IF (LPIV.EQ.NPIVP1) THEN
                 GOTO 416
              ENDIF
              CALL SMUMPS_319( A, LA, IW, LIW,
     *             IOLDPS, NPIVP1, LPIV, POSELT, NASS,
     *             LDA, NFRONT, 1 )
 416          LPIV = MAX(IPIV,JMAX)
              NPIVP1 = NPIVP1+1
           ENDDO
           IF(PIVSIZ .EQ. 2) THEN
              A(POSELT+(LDA+1)*NPIV+1) = DETPIV
           ENDIF
           GOTO 420
  460   CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      PIVSIZ = 0
      IFLAG = -10
  420 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_222
      SUBROUTINE SMUMPS_226(IBEG_BLOCK,
     *     NFRONT,NASS,N,INODE,IW,LIW,
     *     A,LA,LDA, UU,
     *     IOLDPS,POSELT,IFINB,LKJIB,PIVSIZ)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFINB,LKJIB,LDA,
     &        NPBEG, IBEG_BLOCK
      REAL    A(LA)
      REAL   UU
      INTEGER IW(LIW)
      REAL    VALPIV
      INTEGER APOS, IOLDPS, POSELT, K1POS, NCB1
      INTEGER NPIV,JROW2
      INTEGER NEL2,LPOS,NEL
      REAL ONE, ALPHA
      REAL  ZERO
      INTEGER PIVSIZ,NPIV_NEW,OFFDAG_OLD,LPOS1,LPOS2,JJ
      INTEGER POSPV1,POSPV2,OFFDAG,IBEG,IEND,J2,K1,K2,IROW
      REAL SWOP,DETPIV,MULT1,MULT2
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0, ZERO=0.0E0)
      NPIV   = IW(IOLDPS+1+XSIZE)
      NPIV_NEW = NPIV + PIVSIZ
      NEL    = NFRONT - NPIV_NEW
      IFINB  = 0
      JROW2 = IW(IOLDPS+3+XSIZE)
      NPBEG = IBEG_BLOCK
      NEL2   = JROW2 - NPIV_NEW
      IF (NEL2.EQ.0) THEN
        IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      IF(PIVSIZ .EQ. 1) THEN
         APOS   = POSELT + NPIV*(NFRONT + 1)
         VALPIV = ONE/A(APOS)
         A(APOS) = VALPIV
         LPOS   = APOS + LDA
         IF ( UU.NE.ZERO) THEN
            CALL SCOPY(NFRONT-NPIV_NEW, A(LPOS), LDA, A(APOS+1), 1)
         ELSE
            CALL SCOPY(NASS-NPIV_NEW, A(LPOS), LDA, A(APOS+1), 1)
         END IF
         CALL SSYR('U', NEL2, -VALPIV, A(LPOS), LDA, A(LPOS+1), 
     &        LDA)
         IF (UU.NE.ZERO) THEN
            CALL SSCAL(NFRONT-NPIV_NEW, VALPIV, A(LPOS), LDA)
         ELSE
            CALL SSCAL(NASS-NPIV_NEW, VALPIV, A(LPOS), LDA)
         ENDIF
         IF (NEL2.GT.0) THEN
            K1POS = LPOS + NEL2*LDA
            IF ( UU.NE.ZERO ) THEN
               NCB1  = NFRONT - JROW2
            ELSE
               NCB1  = NASS   - JROW2
            END IF
            CALL SGER(NEL2, NCB1 , ALPHA, A(APOS+1), 1, 
     &           A(K1POS), LDA, A(K1POS+1), LDA)
         ENDIF
      ELSE
         POSPV1 = POSELT + NPIV*(NFRONT + 1)
         POSPV2 = POSPV1+NFRONT+1
         OFFDAG_OLD = POSPV2 - 1
         OFFDAG = POSPV1+1
         SWOP = A(POSPV2)
         DETPIV = A(OFFDAG)
         A(POSPV2) = A(POSPV1)/DETPIV
         A(POSPV1) = SWOP/DETPIV
         A(OFFDAG) = -A(OFFDAG_OLD)/DETPIV
         A(OFFDAG_OLD) = ZERO
         LPOS1   = POSPV2 + LDA - 1
         LPOS2   = LPOS1+1
         CALL SCOPY(NFRONT-NPIV_NEW, A(LPOS1), LDA, A(POSPV1+2), 1)
         CALL SCOPY(NFRONT-NPIV_NEW, A(LPOS2), LDA, A(POSPV2+1), 1)
         JJ = POSPV2 + NFRONT-1
         IBEG = JJ + 2
         IEND = IBEG
         DO J2 = 1,NEL2
            K1 = JJ
            K2 = JJ+1
            MULT1 = - (A(POSPV1)*A(K1)+A(OFFDAG)*A(K2))
            MULT2 = - (A(OFFDAG)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1+2
            K2 = POSPV2+1
            DO IROW = IBEG,IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1
               K2 = K2 + 1
            ENDDO
            A(JJ) = -MULT1
            A(JJ+1) = -MULT2
            IBEG = IBEG + NFRONT
            IEND = IEND + NFRONT + 1
            JJ = JJ+NFRONT
         ENDDO
         IEND = IEND-1
         DO J2 = JROW2+1,NFRONT
            K1 = JJ
            K2 = JJ+1
            MULT1 = - (A(POSPV1)*A(K1)+A(OFFDAG)*A(K2))
            MULT2 = - (A(OFFDAG)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1+2
            K2 = POSPV2+1
            DO IROW = IBEG,IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1
               K2 = K2 + 1
            ENDDO
            A(JJ) = -MULT1
            A(JJ+1) = -MULT2
            IBEG = IBEG + NFRONT
            IEND = IEND + NFRONT
            JJ = JJ+NFRONT
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_226
      SUBROUTINE SMUMPS_230(NFRONT,N,INODE,IW,LIW,A,LA,
     &           IOLDPS,POSELT)
      IMPLICIT NONE
      INTEGER NFRONT,N,INODE,LA,LIW
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    VALPIV
      INTEGER APOS, POSELT
      INTEGER IOLDPS,NEL
      INTEGER LPOS,JROW
      INTEGER LDA
      REAL ONE
      PARAMETER(ONE=1.0E0)
        APOS   = POSELT 
        VALPIV = ONE/A(APOS)
        A(APOS) = VALPIV
        NEL    = NFRONT - 1
        IF (NEL.EQ.0) GO TO 500
        LDA    = NFRONT
        LPOS   = APOS + LDA
        CALL SSYR('U',NEL, -VALPIV, 
     &             A(LPOS), NFRONT, A(LPOS+1), NFRONT)
          DO JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT
          END DO
  500   CONTINUE
        RETURN
        END SUBROUTINE SMUMPS_230
      SUBROUTINE SMUMPS_234(IBEG_BLOCK,
     *    NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *    LDA,
     *    IOLDPS,POSELT,LKJIB,LKJIT,UUTEMP,KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,LA,LIW, IBEG_BLOCK
      REAL    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB, INODE, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER POSELT, LDA
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1, NEL11
      INTEGER LBP, HF
      INTEGER LPOS,UPOS,APOS
      INTEGER LKJIT
      INTEGER LKJIBOLD, IROW
      INTEGER I, Block
      INTEGER BLSIZE
      REAL UUTEMP
      REAL ONE, ALPHA
      REAL ZERO
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0, ZERO = 0.0E0)
      LKJIBOLD = LKJIB
      NPIV   = IW(IOLDPS+1+XSIZE)
      JROW2  = IABS(IW(IOLDPS+3+XSIZE))
      NPBEG  = IBEG_BLOCK
      HF     = 6 + IW(IOLDPS+5+XSIZE) + XSIZE
      NEL1   = NASS - JROW2
      LKJIW  = NPIV - NPBEG + 1
      NEL11  = NFRONT - NPIV
      IF ( LKJIW .NE. LKJIB ) THEN
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
          LKJIB       = LKJIB + NONEL
          IW(IOLDPS+3+XSIZE)= MIN0(NPIV+LKJIB,NASS)
          LKJIB       = MIN0(LKJIB, NASS - NPIV)
        ELSE
          LKJIB = NASS - NPIV
          IW(IOLDPS+3+XSIZE) = NASS
        ENDIF
        IBEG_BLOCK = NPIV + 1
      ELSEIF (JROW2.LT.NASS) THEN
          IBEG_BLOCK   = NPIV + 1
          IW(IOLDPS+3+XSIZE) = MIN0(JROW2+LKJIB,NASS)
          LKJIB  = MIN0(LKJIB,NASS-NPIV)
      ENDIF
      IF (LKJIW.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF ( NASS - JROW2 > KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NASS - JROW2
        END IF
        IF ( NASS - JROW2 .GT. 0 ) THEN
#if defined(SAK_BYROW)
         DO IROW = JROW2+1, NASS, BLSIZE
           Block = MIN( BLSIZE, NASS - IROW + 1 )
           LPOS = POSELT + (IROW  - 1) * LDA + NPBEG - 1
           UPOS = POSELT + (NPBEG - 1) * LDA + IROW - 1
           APOS = POSELT + (IROW  - 1) * LDA + JROW2
           CALL SGEMM( 'N','N', IROW + Block - JROW2 - 1, Block, LKJIW,
     *                ALPHA, A( UPOS ), LDA,
     *                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         ENDDO
#else
         DO IROW = JROW2+1, NASS, BLSIZE
          Block = MIN( BLSIZE, NASS - IROW + 1 )
           LPOS = POSELT + ( IROW - 1 ) * LDA + NPBEG - 1
           UPOS = POSELT + ( NPBEG - 1 ) * LDA + IROW - 1
           APOS = POSELT + ( IROW - 1 ) * LDA + IROW - 1
           CALL SGEMM( 'N','N', Block, NASS - IROW + 1, LKJIW,
     *                ALPHA, A( UPOS ), LDA,
     *                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         END DO
#endif
        END IF
       LPOS = POSELT + NASS*LDA + NPBEG - 1
       UPOS = POSELT + (NPBEG-1) * LDA + JROW2
       APOS = POSELT + NASS*LDA + JROW2
       IF ( UUTEMP. NE. ZERO ) THEN
         CALL SGEMM('N', 'N', NEL1, NFRONT-NASS, LKJIW, ALPHA, 
     *              A(UPOS), LDA, A(LPOS), LDA, ONE, 
     *              A(APOS), LDA)
       END IF
      ENDIF
      LBP = JROW2 - NPIV
      IF ( LBP  .GT. 0     .AND.
     *     NEL1 .GE. LKJIB .AND.
     *     NEL1 .GE. LBP + LBP ) THEN
      DO I = 1, LBP
        CALL SMUMPS_319( A, LA, IW, LIW,
     *                 IOLDPS, NPIV+I, NASS-LBP+I, POSELT, NASS,
     *                 LDA, NFRONT, 1 )
      END DO
      END IF
  500 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_234
      SUBROUTINE SMUMPS_321(NFRONT,N,INODE,IW,LIW,A,LA,
     *                          IOLDPS,POSELT)
      IMPLICIT NONE
      INTEGER NFRONT,N,INODE,LA,LIW
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    VALPIV
      INTEGER APOS, POSELT
      INTEGER IOLDPS,NPIV,NEL
      INTEGER LPOS,JROW
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      REAL ONE
      DATA ONE /1.0E0/
        NPIV   = IW(IOLDPS+1+XSIZE)
        NEL    = NFRONT - NPIV - 1
        APOS   = POSELT + (NPIV)*NFRONT + NPIV
        IF (NEL.EQ.0) GO TO 500
        VALPIV = A(APOS)
        LPOS   = APOS + NFRONT
        CALL SSYR('U', NEL, -VALPIV, A(LPOS), NFRONT, A(LPOS+1), 
     &             NFRONT)
        DO JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT
        ENDDO
  500   RETURN
        END SUBROUTINE SMUMPS_321
        SUBROUTINE SMUMPS_319( A, LA, IW, LIW,
     *                       IOLDPS, NPIVP1, IPIV, POSELT, NASS,
     *                       LDA, NFRONT, LEVEL )
        IMPLICIT NONE
      INTEGER LA, LIW, IOLDPS, NPIVP1, POSELT, IPIV
      INTEGER LDA, NFRONT, NASS, LEVEL
      REAL A( LA )
      INTEGER IW( LIW )
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INTEGER ISW, ISWPS1, ISWPS2, APOS, IDIAG, HF
      REAL SWOP
            APOS = POSELT + LDA*(IPIV-1) + NPIVP1-1
            IDIAG = APOS + IPIV - NPIVP1
            HF = 6 + IW( IOLDPS + 5 + XSIZE) + XSIZE
            ISWPS1 = IOLDPS + HF + NPIVP1 - 1
            ISWPS2 = IOLDPS + HF + IPIV - 1
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            ISW = IW(ISWPS1+NFRONT)
            IW(ISWPS1+NFRONT) = IW(ISWPS2+NFRONT)
            IW(ISWPS2+NFRONT) = ISW
            IF ( LEVEL .eq. 2 ) THEN
              CALL SSWAP( NPIVP1 - 1,
     *            A( POSELT + NPIVP1 - 1 ), LDA,
     *            A( POSELT + IPIV - 1 ), LDA )
            END IF
            CALL SSWAP( NPIVP1-1, A( POSELT+(NPIVP1-1) * LDA ), 1,
     *                          A( POSELT + (IPIV-1) * LDA ), 1 )
             CALL SSWAP( IPIV - NPIVP1 - 1,
     *                  A( POSELT + NPIVP1 * LDA + NPIVP1 - 1),
     *                  LDA, A( APOS + 1 ), 1 )
            SWOP = A( IDIAG )
            A( IDIAG  ) = A( POSELT+(NPIVP1-1)*LDA+NPIVP1-1 )
            A( POSELT+(NPIVP1-1)*LDA+NPIVP1-1) = SWOP
            CALL SSWAP( NASS - IPIV, A( APOS + LDA ), LDA,
     *                  A( IDIAG + LDA ), LDA )
            IF ( LEVEL .eq. 1 ) THEN
              CALL SSWAP( NFRONT - NASS,
     *        A( APOS  + ( NASS - IPIV + 1 ) * LDA ), LDA,
     *        A( IDIAG + ( NASS - IPIV + 1 ) * LDA ), LDA )
            END IF
#if defined(PARPIV)
            IF ( LEVEL .eq. 2) THEN
              APOS=POSELT+LDA*LDA-1
              SWOP=A(APOS+NPIVP1)
              A(APOS+NPIVP1)=A(APOS+IPIV)
              A(APOS+IPIV)=SWOP
            ENDIF
#endif
        RETURN
        END SUBROUTINE SMUMPS_319
      SUBROUTINE SMUMPS_237(NFRONT,NASS,N,INODE,
     *    IW,LIW,A,LA,
     *    LDA,
     *    IOLDPS,POSELT,KEEP,KEEP8,UU, ETATASS)
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,INODE,LA,LIW
      REAL    A(LA)
      REAL   UU
      INTEGER IW(LIW) 
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER POSELT, LDA
      INTEGER IOLDPS, ETATASS
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INTEGER BLSIZE, BLSIZE2, Block, IROW, NPIV, I, IROWEND
      INTEGER UPOS, APOS, LPOS
      INTEGER I2, I2END, Block2
      REAL  ONE, ALPHA, BETA
      REAL ZERO
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0,ZERO=0.0E0)
      IF (ETATASS.EQ.1) THEN
        BETA = ZERO
      ELSE
        BETA = ONE
      ENDIF
      IF ( NFRONT - NASS > KEEP(57) ) THEN
        BLSIZE = KEEP(58)
      ELSE
        BLSIZE = NFRONT - NASS
      END IF
      BLSIZE2 = 50
      NPIV = IW( IOLDPS + 1 + XSIZE)
      IF ( NFRONT - NASS .GT. 0 ) THEN
       IF ( UU.eq.ZERO ) THEN
         CALL STRSM( 'L', 'U', 'T', 'U',
     *               NPIV, NFRONT-NPIV, ONE,
     *               A( POSELT ), LDA,
     *               A( POSELT + LDA * NPIV ), LDA )
       ENDIF
       DO IROWEND = NFRONT - NASS, 1, -BLSIZE
        Block = MIN( BLSIZE, IROWEND )
        IROW  = IROWEND - Block + 1
        LPOS = POSELT + NASS*LDA + (IROW-1) * LDA
        APOS = POSELT + NASS*LDA + (IROW-1) * LDA + NASS + IROW - 1
        UPOS = POSELT + NASS
        IF (UU.NE.ZERO) THEN
          UPOS = POSELT + NASS + IROW - 1
        ENDIF
        IF (UU.EQ.ZERO) THEN
         DO I = 1, NPIV
          CALL SCOPY( Block, A(LPOS+I-1), LDA, A(UPOS+(I-1)*LDA),1)
          CALL SSCAL( Block, A(POSELT+(LDA+1)*(I-1)),
     *                A( LPOS + I - 1 ), LDA )
         ENDDO
        ENDIF
        DO I2END = Block, 1, -BLSIZE2
          Block2 = MIN(BLSIZE2, I2END)
          I2 = I2END - Block2+1
          CALL SGEMM('N', 'N', Block2, Block-I2+1, NPIV, ALPHA,
     *               A(UPOS+I2-1), LDA, A(LPOS+(I2-1)*LDA),LDA,
     *               BETA,
     *               A(APOS + I2-1 + (I2-1)*LDA), LDA)
        ENDDO
        IF ( NFRONT - NASS - IROW + 1 - Block > 0 ) THEN
        CALL SGEMM( 'N', 'N', Block, NFRONT-NASS-Block-IROW+1, NPIV,
     *              ALPHA,  A( UPOS ), LDA,
     *              A( LPOS + LDA * Block ), LDA,
     *              BETA,
     *              A( APOS + LDA * Block ), LDA )
        ENDIF
       END DO
      END IF
      RETURN
      END SUBROUTINE SMUMPS_237
      SUBROUTINE SMUMPS_320( BUF, BLOCK_SIZE,
     *                           MYROW, MYCOL, NPROW, NPCOL,
     *                           A, LOCAL_M, LOCAL_N, N, MYID, COMM )
      IMPLICIT NONE
      INTEGER BLOCK_SIZE, NPROW, NPCOL, LOCAL_M, LOCAL_N, N, COMM
      INTEGER MYROW, MYCOL, MYID
      REAL BUF( BLOCK_SIZE * BLOCK_SIZE )
      REAL A( LOCAL_M, LOCAL_N )
      INTEGER NBLOCK, IBLOCK, JBLOCK, IBLOCK_SIZE, JBLOCK_SIZE
      INTEGER ROW_SOURCE, ROW_DEST, COL_SOURCE, COL_DEST
      INTEGER IGLOB, JGLOB
      INTEGER IROW_LOC_SOURCE, JCOL_LOC_SOURCE
      INTEGER IROW_LOC_DEST, JCOL_LOC_DEST
      INTEGER PROC_SOURCE, PROC_DEST
      NBLOCK = ( N - 1 ) / BLOCK_SIZE + 1
      DO IBLOCK = 1, NBLOCK
        IF ( IBLOCK .NE. NBLOCK
     *    ) THEN
          IBLOCK_SIZE = BLOCK_SIZE
        ELSE
          IBLOCK_SIZE = N - ( NBLOCK - 1 ) * BLOCK_SIZE
        END IF
        ROW_SOURCE = MOD( IBLOCK - 1, NPROW ) 
        COL_DEST   = MOD( IBLOCK - 1, NPCOL )
        IGLOB = ( IBLOCK - 1 ) * BLOCK_SIZE + 1
        IROW_LOC_SOURCE = BLOCK_SIZE *
     *                    ( ( IGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     *                  + MOD( IGLOB - 1, BLOCK_SIZE ) + 1
        JCOL_LOC_DEST   = BLOCK_SIZE *
     *                    ( ( IGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     *                  + MOD( IGLOB - 1, BLOCK_SIZE ) + 1
        DO JBLOCK = 1, IBLOCK
          IF ( JBLOCK .NE. NBLOCK
     *      ) THEN
            JBLOCK_SIZE = BLOCK_SIZE
          ELSE
            JBLOCK_SIZE = N - ( NBLOCK - 1 ) * BLOCK_SIZE
          END IF
          COL_SOURCE = MOD( JBLOCK - 1, NPCOL )
          ROW_DEST   = MOD( JBLOCK - 1, NPROW )
          PROC_SOURCE = ROW_SOURCE * NPCOL + COL_SOURCE
          PROC_DEST   = ROW_DEST   * NPCOL + COL_DEST
          IF ( PROC_SOURCE .eq. PROC_DEST ) THEN
           IF ( MYID .eq. PROC_DEST ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            JCOL_LOC_SOURCE = BLOCK_SIZE *
     *                  ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            IROW_LOC_DEST   = BLOCK_SIZE *
     *                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            IF ( IBLOCK .eq. JBLOCK ) THEN
              IF ( IBLOCK_SIZE .ne. JBLOCK_SIZE ) THEN
                WRITE(*,*) MYID,': Error in calling transdiag:unsym'
                CALL SMUMPS_ABORT()
              END IF
              CALL SMUMPS_327( A( IROW_LOC_SOURCE,
     *                 JCOL_LOC_SOURCE),
     *                 IBLOCK_SIZE, LOCAL_M )
            ELSE
              CALL SMUMPS_326(
     *           A( IROW_LOC_SOURCE, JCOL_LOC_SOURCE ),
     *           A( IROW_LOC_DEST, JCOL_LOC_DEST ),
     *           IBLOCK_SIZE, JBLOCK_SIZE, LOCAL_M )
            END IF
           END IF
          ELSE IF (  MYROW .eq. ROW_SOURCE 
     *    .AND. MYCOL .eq. COL_SOURCE ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            JCOL_LOC_SOURCE = BLOCK_SIZE *
     *                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            CALL SMUMPS_293( BUF,
     *           A( IROW_LOC_SOURCE, JCOL_LOC_SOURCE ), LOCAL_M,
     *           IBLOCK_SIZE, JBLOCK_SIZE, COMM, PROC_DEST )
          ELSE IF ( MYROW .eq. ROW_DEST 
     *    .AND.     MYCOL .eq. COL_DEST ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            IROW_LOC_DEST   = BLOCK_SIZE *
     *                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            CALL SMUMPS_281( BUF,
     *           A( IROW_LOC_DEST, JCOL_LOC_DEST ), LOCAL_M,
     *           JBLOCK_SIZE, IBLOCK_SIZE, COMM, PROC_SOURCE )
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_320
      SUBROUTINE SMUMPS_293( BUF, A, LDA, M, N, COMM, DEST )
      IMPLICIT NONE
      INTEGER M, N, LDA, DEST, COMM
      REAL BUF(*), A(LDA,*)
      INTEGER I, IBUF, IERR
      INTEGER J
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      IBUF = 1
      DO J = 1, N
        BUF( IBUF: IBUF + M - 1 ) = A( 1 : M, J )
        DO I = 1, M
        END DO
        IBUF = IBUF + M
      END DO
      CALL MPI_SEND( BUF, M * N, MPI_REAL,
     *     DEST, SYMMETRIZE, COMM, IERR )
      RETURN
      END SUBROUTINE SMUMPS_293
      SUBROUTINE SMUMPS_281( BUF, A, LDA, M, N, COMM, SOURCE )
      IMPLICIT NONE
      INTEGER LDA, M, N, COMM, SOURCE
      REAL BUF(*), A( LDA, *)
      INTEGER I, IBUF, IERR
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      CALL MPI_RECV( BUF(1), M * N, MPI_REAL, SOURCE,
     *               SYMMETRIZE, COMM, STATUS, IERR )
      IBUF = 1
      DO I = 1, M
        CALL SCOPY( N, BUF(IBUF), 1, A(I,1), LDA )
        IBUF = IBUF + N
      END DO
      RETURN
      END
      SUBROUTINE SMUMPS_327( A, N, LDA )
      IMPLICIT NONE
      INTEGER N,LDA
      REAL A( LDA, * )
      INTEGER I, J
      DO I = 2, N
        DO J = 1, I - 1
          A( J, I ) = A( I, J )
        END DO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_327
      SUBROUTINE SMUMPS_326( A1, A2, M, N, LD )
      IMPLICIT NONE
      INTEGER M,N,LD
      REAL A1( LD,* ), A2( LD, * )
      INTEGER I, J
      DO J = 1, N
        DO I = 1, M
          A2( J, I ) = A1( I, J )
        END DO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_326
      RECURSIVE SUBROUTINE SMUMPS_274( 
     *   COMM_LOAD, ASS_IRECV,
     *   BUFR, LBUFR,
     *   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     *   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     *   A, LA, NIRBDU, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     *   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     *   MYID, COMM, IFLAG, IERROR, NBFIN,
     *
     *    PTLUST_S, PTRFAC, root, OPASSW, OPELIW, ITLOC, FILS,  
     *    PTRARW, PTRAIW, INTARR, DBLARR,
     *    ICNTL, KEEP,KEEP8, IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      USE SMUMPS_BUFFER
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER*8 KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, LIW, LA
      INTEGER COMP
      INTEGER NIRBDU, IFLAG, IERROR, POSFAC, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     *        PTRAST(KEEP(28)), NSTK_S(KEEP(28))
      INTEGER NBPROCFILS( KEEP(28) ), STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *        ITLOC(N), FILS(N), ND(KEEP(28))
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER FRERE_STEPS(KEEP(28))
      INTEGER INTARR( MAX(1,KEEP(14)) )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      REAL DBLARR( MAX(1,KEEP(13)) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      LOGICAL FLAG
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER PIVI
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INTEGER POSPV1,POSPV2,OFFDAG,J2,K1,K2,JJ,LPOS1
      REAL MULT1,MULT2
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER INODE, POSITION, NPIV, IERR, IERR_MPI
      INTEGER LAELL, NCOL, POSBLOCFACTO, NROW
      INTEGER IOLDPS, POSELT, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW, DEST
      INTEGER ICT11, LPOS, LPOS2, DPOS, UPOS
      INTEGER I, IPOS, KPOS, IPIV, FPERE, NSLAVES_TOT,
     *        NSLAVES_FOLLOW,  NB_BLOC_FAC
      INTEGER LCONT,NELIM,NASS, LDA, NCOL_TO_SEND,
     *        SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON, SHIFT_VAL_SON
      INTEGER IPOSK, JPOSK, NPIVSENT, Block, IROW, BLSIZE
      INTEGER allocok
      REAL, DIMENSION(:),ALLOCATABLE :: UIP21K
      INTEGER, DIMENSION(:), ALLOCATABLE :: LIST_SLAVES_FOLLOW
      INTEGER ITYPE2
      PARAMETER (ITYPE2=2)
      LOGICAL LASTBL
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      REAL ONE,ALPHA
      REAL VALPIV
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0)
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      FPERE    = -1
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NSLAVES_TOT, 1,
     *                 MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NB_BLOC_FAC, 1,
     *                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LAELL = NPIV * NCOL
      IF ( NPIV.GT.0 ) THEN
       IF ( LRLU .LT. LAELL .OR. IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          IERROR = LAELL - LRLU
          GOTO 700
        END IF
        CALL SMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *       NIRBDU, LRLU, IPTRLU,
     *       IWPOS, IWPOSCB, PTRIST, PTRAST,
     *       STEP, PIMASTER, PAMASTER, ITLOC)
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress ass..blocfacto: LRLU,LRLUS='
     *       ,LRLU,LRLUS
             IFLAG = -9
             IERROR = LAELL -LRLU
             GOTO 700
        END IF
        IF ( IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
          IFLAG = -8
          IERROR = IWPOS + NPIV - 1 - IWPOSCB
          GOTO 700
        END IF
       END IF
       LRLU  = LRLU - LAELL
       LRLUS = LRLUS - LAELL
      ENDIF
      KEEP(67) = MIN(LRLUS, KEEP(67))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL SMUMPS_471(.FALSE.,.FALSE.,
     *                           LA-LRLUS,0,LAELL,KEEP,KEEP8,LRLU)
      IF ( NPIV.GT.0 ) THEN
        IPIV = IWPOS
        IWPOS = IWPOS + NPIV
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 IW( IPIV ), NPIV,
     *                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *              A(POSBLOCFACTO), NPIV*NCOL, MPI_REAL,
     *              COMM, IERR )
      ENDIF
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
        WRITE(*,*) MYID,
     *   ': Internal ERROR 1  in SMUMPS_274',
     *   ' INODE =', INODE,
     *   ' MAITRE_DESC_BANDE not yet received from ', MSGSOU
        CALL SMUMPS_ABORT()
      ENDIF
      DO WHILE ( NBPROCFILS(STEP(INODE)) .NE. 0 ) 
        BLOCKING = .TRUE.
        SET_IRECV=.FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     *    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MPI_ANY_SOURCE, CONTSIPERENIV2,
     *    STATUS,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END  DO
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     *    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     *    STATUS,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC, 
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      IOLDPS = PTRIST(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      LCONT1 = IW( IOLDPS + XSIZE)
      NASS1  = IW( IOLDPS + 1 + XSIZE)
      NROW1  = IW( IOLDPS + 2 + XSIZE)
      NPIV1  = IW( IOLDPS + 3 + XSIZE)
      NSLAV1 = IW( IOLDPS + 5 + XSIZE)
      NSLAVES_FOLLOW = NSLAV1 - 2
      HS     = 6 + NSLAV1 + XSIZE
      NCOL1  = LCONT1 + NPIV1
      IF ( LASTBL ) THEN
        IW( IOLDPS + 6 + XSIZE ) = IW( IOLDPS + 6 + XSIZE ) -
     *  ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) * 
     *  NB_BLOC_FAC
      END IF
      IF (NPIV.EQ.0) GOTO 200
      IF ( NPIV1 + NCOL .NE. NASS1 ) THEN
        WRITE(*,*) 'SymBLFC Error: NPIV1 + NCOL .NE. NASS1 :',
     *             NPIV1,NCOL,NASS1
        CALL SMUMPS_ABORT()
      END IF
      ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
      DO 100 I = 1, NPIV
         PIVI = ABS(IW(IPIV+I-1))
         IF (PIVI.EQ.I) GOTO 100
        ISW = IW(ICT11+I)
        IW(ICT11+I) = IW(ICT11+PIVI)
        IW(ICT11+PIVI) = ISW
        IPOS = POSELT + NPIV1 + I - 1
        KPOS = POSELT + NPIV1 + PIVI - 1
        CALL SSWAP(NROW1, A(IPOS), NCOL1, A(KPOS), NCOL1)
 100  CONTINUE
      ALLOCATE( UIP21K( NPIV * NROW1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IFLAG = -13
        IERROR = NPIV * NROW1
        GOTO 700
      END IF
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
        ALLOCATE( LIST_SLAVES_FOLLOW ( NSLAVES_FOLLOW ),
     *            stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          IFLAG = -13
          IERROR = NSLAVES_FOLLOW
          GOTO 700
        END IF
        LIST_SLAVES_FOLLOW(1:NSLAVES_FOLLOW)=
     *  IW(IOLDPS+8+XSIZE:IOLDPS+7+XSIZE+NSLAVES_FOLLOW)
      END IF
      CALL STRSM( 'L', 'U', 'T', 'U', NPIV, NROW1, ONE,
     *             A( POSBLOCFACTO ), NCOL, A(POSELT+NPIV1), NCOL1 )
      LPOS = POSELT + NPIV1
      UPOS = 1
      DO I = 1, NROW1
        UIP21K( UPOS: UPOS + NPIV-1 ) = A( LPOS: LPOS+NPIV-1)
        LPOS = LPOS + NCOL1
        UPOS = UPOS + NPIV
      END DO
      LPOS = POSELT + NPIV1
      DPOS = POSBLOCFACTO
      I = 1
      DO
         IF(I .GT. NPIV) EXIT
         IF(IW(IPIV+I-1) .GT. 0) THEN
            CALL SSCAL( NROW1, A(DPOS), A(LPOS), NCOL1 )
            LPOS = LPOS + 1
            DPOS = DPOS + NCOL + 1
            I = I+1
         ELSE
            POSPV1 = DPOS
            POSPV2 = DPOS+ NCOL + 1
            OFFDAG = POSPV1+1
            LPOS1 = LPOS
            DO J2 = 1,NROW1
               K1 = JJ
               K2 = JJ+1
               MULT1 = A(POSPV1)*A(LPOS1)+A(OFFDAG)*A(LPOS1+1)
               MULT2 = A(OFFDAG)*A(LPOS1)+A(POSPV2)*A(LPOS1+1)
               A(LPOS1) = MULT1
               A(LPOS1+1) = MULT2
               LPOS1 = LPOS1 + NCOL1
            ENDDO
            LPOS = LPOS + 2
            DPOS = POSPV2 + NCOL + 1
            I = I+2
         ENDIF
      ENDDO
      LPOS2 = POSELT + NPIV1
      UPOS = POSBLOCFACTO+NPIV
      LPOS  = LPOS2 + NPIV
      CALL SGEMM('N','N', NCOL-NPIV,NROW1,NPIV,ALPHA,A(UPOS),NCOL,
     *           A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
      DPOS = POSELT + NCOL1 - NROW1
      IF ( NROW1 .GT. KEEP(7) ) THEN
        BLSIZE = KEEP(8)
      ELSE
        BLSIZE = NROW1
      ENDIF
      IF ( NROW1 .GT. 0 ) THEN
      DO IROW = 1, NROW1, BLSIZE
        Block = MIN( BLSIZE, NROW1 - IROW + 1 )
        DPOS  = POSELT + NCOL1 - NROW1 + ( IROW - 1 ) * ( NCOL1 + 1 )
        LPOS2 = POSELT + NPIV1         + ( IROW - 1 ) * NCOL1
        UPOS  = ( IROW - 1 ) * NPIV + 1
        DO I = 1, Block
          CALL SGEMV( 'T', NPIV, Block-I+1, ALPHA,
     *                A( LPOS2 + (I - 1) * NCOL1 ), NCOL1,
     *                UIP21K( UPOS + NPIV * ( I - 1 ) ), 1, ONE,
     *                A(DPOS+(NCOL1+1)*(I-1)),NCOL1 )
        END DO
       IF ( NROW1-IROW+1-Block .ne. 0 )
     * CALL SGEMM( 'T', 'N', Block, NROW1-IROW+1-Block, NPIV, ALPHA,
     *             UIP21K( UPOS ), NPIV,
     *             A( LPOS2 + Block * NCOL1 ), NCOL1, ONE,
     *             A( DPOS + Block * NCOL1 ), NCOL1 )
      ENDDO
      ENDIF
      FLOP1 = DBLE(NROW1) * DBLE(NPIV) *
     *           ( 2 * NCOL  - NPIV + NROW1 +1 )
      FLOP1 = -FLOP1
      CALL SMUMPS_190( 1, .FALSE., FLOP1, KEEP,KEEP8 )
 200  CONTINUE
      IW(IOLDPS+XSIZE) = IW(IOLDPS+XSIZE) - NPIV
      IW(IOLDPS + 3+XSIZE) = IW(IOLDPS+3+XSIZE) + NPIV
      IF (LASTBL) IW(IOLDPS+1+XSIZE) = IW(IOLDPS + 3+XSIZE)
      LRLU  = LRLU + LAELL
      LRLUS = LRLUS + LAELL
      POSFAC = POSFAC - LAELL
      IWPOS = IWPOS - NPIV
      CALL SMUMPS_471(.FALSE.,.FALSE.,
     *                           LA-LRLUS,0,-LAELL,KEEP,KEEP8,LRLU)
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
         IPOSK = NPIV1 + 1
         JPOSK = NCOL1 - NROW1 + 1
           NPIVSENT = NPIV
          IERR = -1
           DO WHILE ( IERR .eq. -1 )
            CALL SMUMPS_64(
     *                    INODE, NPIVSENT, FPERE,
     *                    IPOSK, JPOSK,
     *                    UIP21K, NROW1,
     *                    NSLAVES_FOLLOW,
     *                    LIST_SLAVES_FOLLOW(1),
     *                    COMM, IERR )
            IF (IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV= .FALSE.
              MESSAGE_RECEIVED = .FALSE.
              CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     *         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *         MPI_ANY_SOURCE, MPI_ANY_TAG,
     *         STATUS, 
     *         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *         IWPOS, IWPOSCB, IPTRLU,
     *         LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *         PTLUST_S, PTRFAC,
     *         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *         IFLAG, IERROR, COMM,
     *         NBPROCFILS,
     *         IPOOL, LPOOL, LEAF,
     *         NBFIN, MYID, SLAVEF,
     *         root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *         INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *         LPTRAR, NELT, FRTPTR, FRTELT, 
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE )
             IF ( IFLAG .LT. 0 ) GOTO 600
            END IF
            IF ( IERR .eq. -2 ) THEN
             IFLAG = -17
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
            END IF
           END DO
           DEALLOCATE( LIST_SLAVES_FOLLOW )
      END IF
      IF ( NPIV .NE. 0 ) DEALLOCATE( UIP21K )
      IF ( LASTBL .AND. IW(PTRIST(STEP(INODE))+6+XSIZE ) .eq. 0
     *   .and. KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 1 )
     *   THEN
         DEST = SMUMPS_275( STEP(INODE), PROCNODE_STEPS, SLAVEF )
         CALL SMUMPS_62( INODE, DEST, END_NIV2_LDLT,
     *                             COMM, IERR )
         IF ( IERR .LT. 0 ) THEN
           write(*,*) ' Internal error in PROCESS_SYM_BLOCFACTO.'
           IFLAG = -99
           GOTO 700
         END IF
      END IF
      IF ( LASTBL  .AND. IW(PTRIST(STEP(INODE))+6+XSIZE ) .eq. 0
     *     .AND. (KEEP(38).EQ.FPERE) ) THEN
       IOLDPS = PTRIST(STEP( INODE ))
       LCONT  = IW(IOLDPS+XSIZE)
       NROW   = IW(IOLDPS+2+XSIZE)
       NPIV   = IW(IOLDPS+3+XSIZE)
       NASS   = IW(IOLDPS+4+XSIZE)
       NELIM  = NASS-NPIV
       NCOL_TO_SEND =  LCONT-NELIM
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5+XSIZE) + XSIZE
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW + NASS
       SHIFT_VAL_SON      = NASS
       LDA                = LCONT + NPIV
       CALL SMUMPS_80( COMM_LOAD, ASS_IRECV, 
     *    N, INODE, FPERE, 
     *    PTRIST, PTRAST, 
     *    root, NROW, NCOL_TO_SEND, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDA, 
     *    ROOT_CONT_STATIC, MYID, COMM,
     *    
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER,
     *    PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, .FALSE., ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
       IF ( IFLAG < 0 ) GOTO 600
       IF (NELIM.EQ.0) THEN
        CALL SMUMPS_314( N, INODE,  
     *    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, ITYPE2
     $       )
       ENDIF
      ENDIF
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_274
      SUBROUTINE SMUMPS_141( COMM_LOAD, ASS_IRECV, 
     *           N, INODE, FPERE, IW, LIW, A, LA,
     *           UU, NOFFW,
     *           NPVW,
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     *             PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED )
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, LA, NOFFW, NPVW
      REAL A( LA )
      REAL UU, SEUIL
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER LPTRAR, NELT
      INTEGER NBFIN, IFLAG, IERROR, LEAF, LPOOL, SLAVEF, 
     *        POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU,
     *        COMP 
      INTEGER NB_BLOC_FAC
      INTEGER ICNTL(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IW( LIW )
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL), ITLOC(N)
      INTEGER PTRARW(LPTRAR), PTRAIW(LPTRAR), ND( KEEP(28) )
      INTEGER FRERE(KEEP(28)), FILS(N) 
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)),
     * PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *        
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)),
     *        NSTK_S(KEEP(28)), NBPROCFILS(KEEP(28)),
     *        PROCNODE_STEPS(KEEP(28)), STEP(N)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      REAL  DBLARR(MAX(1,KEEP(13)))
      LOGICAL AVOID_DELAYED
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEGKJI, NBTLKJ
      INTEGER NASS, IEND, IOLDPS, POSELT, LDAFS,allocok, IBEG_BLOCK
      LOGICAL LASTBL 
      REAL UUTEMP
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      EXTERNAL SMUMPS_223, SMUMPS_235,
     *         SMUMPS_227, SMUMPS_294,
     *         SMUMPS_44
      LOGICAL STATICMODE
      REAL SEUIL_LOC
      INTEGER PIVSIZ,IWPOSPIV
      REAL ONE
      PARAMETER(ONE = 1.0E0)
      INOPV = 0
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
        STATICMODE = .TRUE.
        UUTEMP=UU
        SEUIL_LOC = MAX(SEUIL,EPSILON(SEUIL))
      ELSE
        UUTEMP=UU
      ENDIF
      IBEG_BLOCK=1
      NB_BLOC_FAC = 0
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST( STEP( INODE ))
      NFRONT = IW(IOLDPS+XSIZE)
      NASS   = IABS(IW(IOLDPS+2+XSIZE))
      LDAFS  = NASS
      IF (NASS .GT. KEEP(3)) THEN
        NBTLKJ = MIN( KEEP(6), NASS )
      ELSE
        NBTLKJ = MIN( KEEP(5), NASS )
      ENDIF
      IW(IOLDPS+3+XSIZE) = MIN0(NASS,NBTLKJ)
      ALLOCATE( IPIV( NASS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) MYID,' : FACTO_NIV2 :failed to allocate ',NASS,
     * ' integers'
        IFLAG=-13
        IERROR=NASS
        GO TO 490
      END IF
 50   CONTINUE
      IBEGKJI = IBEG_BLOCK
      CALL SMUMPS_223(
     *                NFRONT,NASS,IBEGKJI, NASS, IPIV,
     *                N,INODE,IW,LIW,A,LA,NOFFW,INOPV,
     *                IFLAG,IOLDPS,POSELT,UU, SEUIL_LOC,
     *                KEEP,KEEP8,PIVSIZ)
      IF (IFLAG.LT.0) GOTO 490
         IF(INOPV.EQ. 1 .AND. STATICMODE) THEN
            INOPV = -1
            GOTO 50
         ENDIF
      IF (INOPV.GE.1) THEN
          LASTBL = (INOPV.EQ.1)
          IEND = IW(IOLDPS+1+XSIZE)
          CALL SMUMPS_294( COMM_LOAD, ASS_IRECV,
     *             N, INODE, FPERE, IW, LIW, 
     *             IOLDPS, POSELT, A, LA, LDAFS,
     *             IBEGKJI, IEND, IPIV, NASS,LASTBL, NB_BLOC_FAC,
     *
     *             COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     *             PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
      ENDIF
      IF (INOPV.EQ.1) GO TO 500
      IF (INOPV.EQ.2) THEN
         CALL SMUMPS_235(IBEG_BLOCK,
     *            NASS,N,INODE,IW,LIW,A,LA,
     *            LDAFS, 
     *            IOLDPS,POSELT,NBTLKJ,KEEP(4),KEEP,KEEP8)
         GOTO 50
      ENDIF
      NPVW = NPVW + PIVSIZ
      IF (NASS.LE.1) THEN
        IFINB = -1
        IF (NASS == 1) A(POSELT)=ONE/A(POSELT)
      ELSE
         CALL SMUMPS_227(IBEG_BLOCK,
     *             NASS, N,INODE,IW,LIW,A,LA,
     *             LDAFS, IOLDPS,POSELT,IFINB,
     *             NBTLKJ,KEEP(4),PIVSIZ)
         IF(PIVSIZ .EQ. 2) THEN
            IWPOSPIV = IOLDPS+XSIZE+IW(IOLDPS+1+XSIZE)+6+
     &                 IW(IOLDPS+5+XSIZE)
            IW(IWPOSPIV+NFRONT) = -IW(IWPOSPIV+NFRONT)
         ENDIF
      ENDIF
      IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + PIVSIZ
       IF (IFINB.EQ.0) GOTO 50
       IF ((IFINB.EQ.1).OR.(IFINB.EQ.-1)) THEN
          LASTBL = (IFINB.EQ.-1) 
          IEND = IW(IOLDPS+1+XSIZE)
          CALL SMUMPS_294(COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE, IW, LIW, 
     *             IOLDPS, POSELT, A, LA, LDAFS, 
     *             IBEGKJI, IEND, IPIV, NASS, LASTBL,NB_BLOC_FAC,
     *
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     *             PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
       ENDIF
       IF (IFINB.EQ.(-1)) GOTO 500
       NPIV   = IW(IOLDPS+1+XSIZE)
       CALL SMUMPS_235(IBEG_BLOCK,
     *            NASS,N,INODE,IW,LIW,A,LA,
     *            LDAFS, 
     *            IOLDPS,POSELT,NBTLKJ,KEEP(4),KEEP,KEEP8)
      GO TO 50
 490  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
 500  CONTINUE
      DEALLOCATE( IPIV )
      RETURN
      END SUBROUTINE SMUMPS_141
      SUBROUTINE SMUMPS_223( NFRONT, NASS,
     *                   IBEGKJI, NASS2, TIPIV,
     *                   N, INODE, IW, LIW,
     *                   A, LA, NNEG, 
     *                   INOPV, IFLAG,
     *                   IOLDPS, POSELT, UU, 
     *                   SEUIL,KEEP,KEEP8,PIVSIZ)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFLAG,INOPV
      INTEGER NASS2, IBEGKJI, NNEG
      INTEGER TIPIV( NASS2 )
      INTEGER PIVSIZ,LPIV
      REAL A(LA) 
      REAL UU, UULOC, SEUIL
      INTEGER IW(LIW) 
      INTEGER   IOLDPS, POSELT
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      include 'mpif.h'
      INTEGER POSPV1,POSPV2,OFFDAG,APOSJ
      INTEGER JMAX
      REAL RMAX,AMAX,TMAX,SWOP,TOL
      REAL DELTA,MAXPIV
      REAL PIVOT,DETPIV
      PARAMETER(TOL = 1.0E-20)
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INTEGER APOSMAX
      INTEGER APOS, LDAFS
      REAL ZERO,ONE
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,J1,JJ,J2,ILOC,K
      INTRINSIC MAX
      DATA ZERO /0.0E0/
      DATA ONE /1.0E0/
        LDAFS = NASS
        UULOC = UU
        PIVSIZ = 1
        NPIV    = IW(IOLDPS+1+XSIZE)
        NPIVP1  = NPIV + 1
        ILOC = NPIVP1 - IBEGKJI + 1
        TIPIV( ILOC ) = ILOC
        NASSW   = IABS(IW(IOLDPS+3+XSIZE))
        APOSMAX = POSELT+LDAFS**2-1
        IF(INOPV .EQ. -1) THEN
           APOS = POSELT + LDAFS*(NPIVP1-1) + NPIV
           POSPV1 = APOS
           IF(ABS(A(APOS)).LT.SEUIL) THEN
              IF(REAL(A(APOS)) .GE. ZERO) THEN
                 A(APOS) = SEUIL
              ELSE
                 A(APOS) = -SEUIL
              ENDIF
              KEEP(98) = KEEP(98)+1
           ENDIF
           GO TO 420
        ENDIF
        INOPV   = 0
        DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + LDAFS*(IPIV-1) + NPIV
            POSPV1 = APOS + IPIV - NPIVP1
            PIVOT = A(POSPV1)
            IF (UULOC.EQ.ZERO) THEN 
              IF (ABS(A(APOS)).EQ.ZERO) GO TO 630
        IF (A(APOS).LT.ZERO) NNEG = NNEG+1
              GO TO 420
            ENDIF
            AMAX = ZERO
            JMAX = 0
            J1 = APOS
            J2 = POSPV1 - 1
            DO JJ=J1,J2
               IF(ABS(A(JJ)) .GT. AMAX) THEN
                  AMAX = ABS(A(JJ))
                  JMAX = IPIV - (POSPV1-JJ)
               ENDIF
            ENDDO
            J1 = POSPV1 + LDAFS
            J2 = NASSW - IPIV
            DO JJ=1,J2
               IF(ABS(A(J1)) .GT. AMAX) THEN
                  AMAX = MAX(ABS(A(J1)),AMAX)
                  JMAX = IPIV + JJ
               ENDIF
               J1 = J1 + LDAFS
            ENDDO
#if defined(PARPIV)
            RMAX = MAX(SEUIL/UULOC,A(APOSMAX+IPIV))
#else
            RMAX = SEUIL/UULOC
#endif
            J2 = NASS - NASSW 
            DO JJ=1,J2
               RMAX = MAX(ABS(A(J1)),RMAX)
               J1 = J1 + LDAFS
            ENDDO
            IF (MAX(AMAX,ABS(PIVOT)).LE.TOL) GO TO 460
            IF (ABS(PIVOT).GT. UULOC*MAX(RMAX,AMAX)
     *           .AND. ABS(PIVOT).GT. SEUIL) THEN
            IF (A(POSPV1).LT.ZERO) NNEG = NNEG+1
               GO TO 415
            END IF
            IF (AMAX.LE.TOL) GO TO 460
            IF (RMAX.LT.AMAX) THEN
               J1 = APOS
               J2 = POSPV1 - 1
               DO JJ=J1,J2
                  IF(POSPV1-JJ .NE. IPIV-JMAX) THEN
                     RMAX = MAX(RMAX,ABS(A(JJ)))
                  ENDIF
               ENDDO
               J1 = POSPV1 + LDAFS
               J2 = NASS - IPIV
               DO JJ=1,J2
                  IF(IPIV+JJ .NE. JMAX) THEN
                     RMAX = MAX(ABS(A(J1)),RMAX)
                  ENDIF
                  J1 = J1 + LDAFS
               ENDDO
            ENDIF            
            APOSJ = POSELT + (JMAX-1)*LDAFS + NPIV
            POSPV2 = APOSJ + JMAX - NPIVP1
            IF (IPIV.LT.JMAX) THEN
               OFFDAG = APOSJ + IPIV - NPIVP1
            ELSE
               OFFDAG = APOS + JMAX - NPIVP1
            END IF
#if defined(PARPIV)
            TMAX = MAX(SEUIL/UULOC,A(APOSMAX+JMAX))
#else
            TMAX = SEUIL/UULOC
#endif
            IF(JMAX .LT. IPIV) THEN
               JJ = POSPV2
               DO K = 1, NASS-JMAX
                  JJ = JJ+NASS
                  IF (JMAX+K.NE.IPIV) THEN
                     TMAX=MAX(TMAX,ABS(A(JJ)))
                  ENDIF
               ENDDO
               DO K =  APOSJ, POSPV2-1
                  TMAX = MAX(TMAX,ABS(A(K)))
               ENDDO
            ELSE
               JJ = POSPV2
               DO K = 1, NASS-JMAX
                  JJ = JJ+NASS
                  TMAX=MAX(TMAX,ABS(A(JJ)))
               ENDDO
               DO K =  APOSJ, POSPV2-1
                  IF (K.NE.OFFDAG) THEN
                     TMAX = MAX(TMAX,ABS(A(K)))
                  ENDIF
               ENDDO
            ENDIF
            DETPIV = A(POSPV1)*A(POSPV2) - A(OFFDAG)**2
            MAXPIV = MAX(ABS(A(POSPV1)),ABS(A(POSPV2)))
            IF (MAXPIV.EQ.ZERO) MAXPIV = ONE
            IF (ABS(DETPIV)/MAXPIV.LE.TOL) GO TO 460
            IF ((ABS(A(POSPV2))*RMAX+AMAX*TMAX)*UULOC.GT.
     +           ABS(DETPIV)) GO TO 460
            IF ((ABS(A(POSPV1))*TMAX+AMAX*RMAX)*UULOC.GT.
     +           ABS(DETPIV)) GO TO 460
           PIVSIZ = 2
           KEEP(105) = KEEP(105)+1
           IF(DETPIV .LT. ZERO) THEN
             NNEG = NNEG+1
           ELSE IF(A(POSPV2) .LT. ZERO) THEN
             NNEG = NNEG+2
           ENDIF
 415       CONTINUE
           LPIV = IPIV
           IF (PIVSIZ.EQ.2) THEN
              LPIV = MIN(IPIV,JMAX)
              TIPIV(ILOC) = -(LPIV - IBEGKJI + 1)
              TIPIV(ILOC+1) = -(MAX(IPIV,JMAX) - IBEGKJI + 1)
           ELSE
              TIPIV(ILOC) = IPIV - IBEGKJI + 1
           ENDIF
           DO K=1,PIVSIZ
              IF (LPIV.EQ.NPIVP1) THEN
                 GOTO 416
              ENDIF
              CALL SMUMPS_319( A, LA, IW, LIW,
     *             IOLDPS, NPIVP1, LPIV, POSELT, NASS,
     *             LDAFS, NFRONT, 2 )
 416          LPIV = MAX(IPIV,JMAX)
              NPIVP1 = NPIVP1+1
           ENDDO
           IF(PIVSIZ .EQ. 2) THEN
              A(POSELT+LDAFS*NPIV+NPIV+1) = DETPIV
           ENDIF
           GOTO 420
  460   CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      IFLAG = -10
  420 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_223
      SUBROUTINE SMUMPS_235(
     *                 IBEG_BLOCK,
     *                 NASS, N, INODE,
     *                 IW, LIW, A, LA,
     *                 LDAFS,
     *                 IOLDPS, POSELT,
     *                 LKJIB, LKJIT, KEEP,KEEP8 )
      IMPLICIT NONE
      INTEGER NASS,N,LA,LIW
      REAL    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB, INODE, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER POSELT, LDAFS, IBEG_BLOCK
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1
      INTEGER HF
      INTEGER LPOS,UPOS,APOS
      INTEGER LKJIT
      INTEGER LKJIBOLD, IROW
      INTEGER J, Block
      INTEGER BLSIZE
      REAL ONE, ALPHA
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0)
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      LKJIBOLD = LKJIB
      NPIV   = IW(IOLDPS+1+XSIZE)
      JROW2  = IABS(IW(IOLDPS+3+XSIZE))
      NPBEG  = IBEG_BLOCK
      HF     = 6 + IW(IOLDPS+5+XSIZE) + XSIZE
      NEL1   = NASS - JROW2
      LKJIW  = NPIV - NPBEG + 1
      IF ( LKJIW .NE. LKJIB ) THEN
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
          LKJIB       = LKJIB + 2 * NONEL
          IW(IOLDPS+3+XSIZE)= MIN0(NPIV+LKJIB,NASS)
          LKJIB       = MIN0(LKJIB, NASS - NPIV)
        ELSE
          LKJIB = NASS - NPIV
          IW(IOLDPS+3+XSIZE) = NASS
        ENDIF
      ELSEIF (JROW2.LT.NASS) THEN
          IW(IOLDPS+3+XSIZE) = MIN0(JROW2+LKJIB,NASS)
      ENDIF
      IBEG_BLOCK = NPIV + 1
      IF (LKJIW.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF ( NASS - JROW2 > KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NASS - JROW2
        END IF
        IF ( NASS - JROW2 .GT. 0 ) THEN
         DO IROW = JROW2+1, NASS, BLSIZE
          Block = MIN( BLSIZE, NASS - IROW + 1 )
          LPOS = POSELT + ( IROW - 1 ) * LDAFS + NPBEG - 1
          UPOS = POSELT + ( NPBEG - 1 ) * LDAFS + IROW - 1
          APOS =  POSELT + (IROW-1) * LDAFS + IROW - 1
          DO J=1, Block
            CALL SGEMV( 'T', LKJIW, Block - J + 1, ALPHA,
     &                  A( LPOS ), LDAFS, A( UPOS ), LDAFS,
     &                  ONE, A( APOS ), LDAFS )
            LPOS = LPOS + LDAFS
            APOS = APOS + LDAFS + 1
            UPOS = UPOS + 1
          END DO
          LPOS = POSELT + ( IROW - 1 + Block ) * LDAFS + NPBEG - 1
          UPOS = POSELT + ( NPBEG - 1 ) * LDAFS + IROW - 1
          APOS = POSELT + ( IROW - 1 + Block ) * LDAFS + IROW - 1
          CALL SGEMM( 'N','N', Block, NASS - IROW + 1 - Block, LKJIW,
     *                ALPHA, A( UPOS ), LDAFS,
     *                A( LPOS ), LDAFS, ONE, A( APOS ), LDAFS )
         END DO
        END IF
      END IF
  500 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_235
      SUBROUTINE SMUMPS_227
     *     ( IBEG_BLOCK, NASS, N, INODE, IW, LIW,
     *     A, LA, LDAFS, 
     *     IOLDPS,POSELT,IFINB,LKJIB,LKJIT,PIVSIZ)
      IMPLICIT NONE
      INTEGER NASS,N,LA,LIW,INODE,IFINB,LKJIB,LDAFS,
     &        NPBEG
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    VALPIV
      INTEGER APOS, IOLDPS, POSELT, K1POS, NCB1
      INTEGER LKJIT, IBEG_BLOCK
      INTEGER NPIV,JROW2
      INTEGER NEL2,LPOS
      REAL ONE, ALPHA
      REAL  ZERO
      INTEGER PIVSIZ,NPIV_NEW,OFFDAG_OLD,LPOS1,LPOS2,JJ
      INTEGER POSPV1,POSPV2,OFFDAG,IBEG,IEND,J2,K1,K2,IROW
      REAL SWOP,DETPIV,MULT1,MULT2
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0, ZERO=0.0E0)
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      NPIV   = IW(IOLDPS+1+XSIZE)
      NPIV_NEW = NPIV + PIVSIZ
      IFINB  = 0
      IF (IW(IOLDPS+3+XSIZE).LE.0) THEN
         IW(IOLDPS+3+XSIZE) = MIN0(NASS,LKJIB)
      ENDIF
      JROW2 = IW(IOLDPS+3+XSIZE)
      NPBEG = IBEG_BLOCK
      NEL2   = JROW2 - NPIV_NEW
      IF (NEL2.EQ.0) THEN
        IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      IF(PIVSIZ .EQ. 1) THEN
         APOS   = POSELT + NPIV*(LDAFS + 1)
         VALPIV = ONE/A(APOS)
         A(APOS) = VALPIV
         LPOS   = APOS + LDAFS
         CALL SCOPY(NASS-NPIV_NEW, A(LPOS), LDAFS, A(APOS+1), 1)
         CALL SSYR('U', NEL2, -VALPIV, A(LPOS), LDAFS, A(LPOS+1), 
     &        LDAFS)
         CALL SSCAL(NASS-NPIV_NEW, VALPIV, A(LPOS), LDAFS)
         IF (NEL2.GT.0) THEN
            K1POS = LPOS + NEL2*LDAFS
            NCB1  = NASS - JROW2
            CALL SGER(NEL2, NCB1 , ALPHA, A(APOS+1), 1, 
     &           A(K1POS), LDAFS, A(K1POS+1), LDAFS)
         ENDIF
      ELSE
         POSPV1 = POSELT + NPIV*(LDAFS + 1)
         POSPV2 = POSPV1+LDAFS+1
         OFFDAG_OLD = POSPV2 - 1
         OFFDAG = POSPV1+1
         SWOP = A(POSPV2)
         DETPIV = A(OFFDAG)
         A(POSPV2) = A(POSPV1)/DETPIV
         A(POSPV1) = SWOP/DETPIV
         A(OFFDAG) = -A(OFFDAG_OLD)/DETPIV
         A(OFFDAG_OLD) = ZERO
         LPOS1   = POSPV2 + LDAFS - 1
         LPOS2   = LPOS1+1
         CALL SCOPY(NASS-NPIV_NEW, A(LPOS1), LDAFS, A(POSPV1+2), 1)
         CALL SCOPY(NASS-NPIV_NEW, A(LPOS2), LDAFS, A(POSPV2+1), 1)
         JJ = POSPV2 + NASS-1
         IBEG = JJ + 2
         IEND = IBEG
         DO J2 = 1,NEL2
            K1 = JJ
            K2 = JJ+1
            MULT1 = - (A(POSPV1)*A(K1)+A(OFFDAG)*A(K2))
            MULT2 = - (A(OFFDAG)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1+2
            K2 = POSPV2+1
            DO IROW = IBEG,IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1
               K2 = K2 + 1
            ENDDO
            A(JJ) = -MULT1
            A(JJ+1) = -MULT2
            IBEG = IBEG + NASS
            IEND = IEND + NASS + 1
            JJ = JJ+NASS
         ENDDO
         IEND = IEND-1
         DO J2 = JROW2+1,NASS
            K1 = JJ
            K2 = JJ+1
            MULT1 = - (A(POSPV1)*A(K1)+A(POSPV1+1)*A(K2))
            MULT2 = - (A(POSPV1+1)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1+2
            K2 = POSPV2+1
            DO IROW = IBEG,IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1
               K2 = K2 + 1
            ENDDO
            A(JJ) = -MULT1
            A(JJ+1) = -MULT2
            IBEG = IBEG + NASS
            IEND = IEND + NASS
            JJ = JJ+NASS
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_227
      RECURSIVE SUBROUTINE SMUMPS_263(
     *   COMM_LOAD, ASS_IRECV,
     *   BUFR, LBUFR,
     *   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     *   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     *   A, LA, NIRBDU, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     *   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     *   MYID, COMM, IFLAG, IERROR, NBFIN,
     *
     *    PTLUST_S, PTRFAC, root, OPASSW, OPELIW, ITLOC, FILS,  
     *    PTRARW, PTRAIW, INTARR, DBLARR,
     *    ICNTL, KEEP,KEEP8, IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     *     )
      USE SMUMPS_BUFFER
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER*8 KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, IPTRLU, 
     &        LRLU, LRLUS, LIW, LA
      INTEGER COMP
      INTEGER NIRBDU, IFLAG, IERROR, POSFAC, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     *        PTRAST(KEEP(28)),
     *        NSTK_S(KEEP(28))
      INTEGER NBPROCFILS(KEEP(28)), STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER ITLOC( N ), FILS( N )
      INTEGER ND( KEEP(28) ), FRERE_STEPS( KEEP(28) )
      INTEGER INTARR( MAX(1,KEEP(14)) )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      REAL  DBLARR( MAX(1,KEEP(13)) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER ITYPE2
      PARAMETER(ITYPE2=2)
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      INTEGER INODE, IPOSK, JPOSK, NCOLU, NPIV, POSITION, IERR
      INTEGER LAELL, POSBLOCFACTO
      INTEGER IOLDPS, POSELT, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NELIM1, NCOL_TO_SEND
      INTEGER LDA
      INTEGER NSLAV1, HS, DEST, NSLAVES_FOLLOW
      INTEGER FPERE
      INTEGER CPOS, LPOS
      INTEGER SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON, SHIFT_VAL_SON
      LOGICAL DYNAMIC
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER allocok
      REAL, ALLOCATABLE, DIMENSION(:) :: UDYNAMIC
      REAL ONE,ALPHA
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0)
      DYNAMIC = .FALSE.
      POSITION  = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, IPOSK, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, JPOSK, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     *                 MPI_INTEGER, COMM, IERR )
      IF ( NPIV .LE. 0 ) THEN
      NPIV = - NPIV
        WRITE(*,*) MYID,':error, received negative NPIV in BLFAC'
        CALL SMUMPS_ABORT()
      END IF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOLU, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LAELL = NPIV * NCOLU
      IF ( LRLU .LT. LAELL ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          IERROR = LAELL - LRLU
          GOTO 700
        END IF
        CALL SMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *        NIRBDU, LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC)
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress ass..blocfacto: LRLU,LRLUS='
     *       ,LRLU,LRLUS
             IFLAG = -9
             IERROR = LAELL -LRLU
             GOTO 700
        END IF
      END IF
      LRLU  = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      KEEP(67) = MIN(LRLUS, KEEP(67))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL SMUMPS_471(.FALSE.,.FALSE.,
     *                           LA-LRLUS,0, LAELL,KEEP,KEEP8,LRLU)
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 A(POSBLOCFACTO), NPIV*NCOLU,
     *                 MPI_REAL,
     *                 COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) DYNAMIC = .TRUE.
      IF ( (PTRIST(STEP( INODE )).NE.0) .AND.
     &  (IPOSK + NPIV -1 .GT. IW(PTRIST(STEP(INODE))+3+XSIZE)) )THEN
        DYNAMIC = .TRUE.
      ENDIF
      IF (DYNAMIC)  THEN
        ALLOCATE(UDYNAMIC(LAELL), stat=allocok)
        if (allocok .GT. 0) THEN
          write(*,*) myid, ' : PB allocation U in blfac_slave '
     *     , LAELL
          IFLAG = -13 
          IERROR = LAELL
          GOTO 700
        endif
        UDYNAMIC(1:LAELL) = A(POSBLOCFACTO:POSBLOCFACTO+LAELL-1)
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        POSFAC = POSFAC - LAELL
      CALL SMUMPS_471(.FALSE.,.FALSE.,
     *          LA-LRLUS,0,-LAELL,KEEP,KEEP8,LRLU)
      ENDIF
      DO WHILE ( PTRIST(STEP(INODE)) .EQ. 0 )
        MSGSOU = SMUMPS_275( STEP(INODE),
     *           PROCNODE_STEPS, SLAVEF )
        SET_IRECV = .FALSE.
        BLOCKING  = .TRUE.
        MESSAGE_RECEIVED = .FALSE.
        CALL SMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MSGSOU, MAITRE_DESC_BANDE,
     *    STATUS, 
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF ( IFLAG .LT. 0 ) GOTO 600
      ENDDO
      DO WHILE ( IPOSK + NPIV -1 .GT.
     *            IW( PTRIST(STEP( INODE )) + 3 +XSIZE) )
        MSGSOU = SMUMPS_275( STEP(INODE), PROCNODE_STEPS, SLAVEF )
        SET_IRECV = .FALSE.
        BLOCKING  = .TRUE.
        MESSAGE_RECEIVED = .FALSE.
        CALL SMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MSGSOU, BLOC_FACTO_SYM,
     *    STATUS, 
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END DO
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL SMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     *    STATUS, 
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      IOLDPS  = PTRIST(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      LCONT1 = IW( IOLDPS + XSIZE )
      NROW1  = IW( IOLDPS + 2  + XSIZE)
      NPIV1  = IW( IOLDPS + 3  + XSIZE)
      NSLAV1 = IW( IOLDPS + 5  + XSIZE)
      HS     = 6 + NSLAV1 + XSIZE
      NCOL1  = LCONT1 + NPIV1
      CPOS = POSELT + JPOSK - 1
      LPOS = POSELT + IPOSK - 1
      IF ( NPIV .GT. 0 ) THEN
       IF (DYNAMIC) THEN
        CALL SGEMM( 'T', 'N', NCOLU, NROW1, NPIV, ALPHA,
     *            UDYNAMIC(1), NPIV,
     *            A( LPOS ), NCOL1, ONE,
     *            A( CPOS ), NCOL1 )
       ELSE
        CALL SGEMM( 'T', 'N', NCOLU, NROW1, NPIV, ALPHA,
     *            A( POSBLOCFACTO ), NPIV,
     *            A( LPOS ), NCOL1, ONE,
     *            A( CPOS ), NCOL1 )
       ENDIF
       FLOP1 = DBLE(NCOLU*NPIV)*DBLE(2*NROW1)
       FLOP1 = -FLOP1
       CALL SMUMPS_190(1, .FALSE., FLOP1, KEEP,KEEP8 )
      ENDIF
      IW( IOLDPS + 6 + XSIZE ) = IW( IOLDPS + 6+ XSIZE ) + 1
      IF (DYNAMIC) THEN
       DEALLOCATE(UDYNAMIC)
      ELSE
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        POSFAC = POSFAC - LAELL
      CALL SMUMPS_471(.FALSE.,.FALSE.,
     *                      LA-LRLUS,0,-LAELL,KEEP,KEEP8,LRLU)
      ENDIF
      NSLAVES_FOLLOW = IW( IOLDPS + 5 +XSIZE ) - 2
      IF ( IW( IOLDPS + 6  +XSIZE) .eq. 0 .and.
     *  KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 1 ) THEN
         DEST = SMUMPS_275( STEP(INODE), PROCNODE_STEPS, SLAVEF )
         CALL SMUMPS_62( INODE, DEST, END_NIV2_LDLT,
     *                             COMM, IERR )
         IF ( IERR .LT. 0 ) THEN
           write(*,*) ' Internal error in PROCESS_BLFAC_SLAVE.'
           IFLAG = -99
           GOTO 700
         END IF
      END IF
      IF ( IW( IOLDPS + 6 + XSIZE ) .eq. 0 
     &         .AND. (KEEP(38).EQ.FPERE) ) THEN
       NPIV1  = IW( IOLDPS + 3 +XSIZE)
       NASS1  = IW(IOLDPS+4+XSIZE)
       NELIM1  = NASS1-NPIV1
       NCOL_TO_SEND =  LCONT1-NELIM1
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5+XSIZE)  + XSIZE
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW1 + NASS1
       SHIFT_VAL_SON      = NASS1
       LDA                = LCONT1 + NPIV1
       CALL SMUMPS_80( COMM_LOAD, ASS_IRECV, 
     *    N, INODE, FPERE, 
     *    PTRIST, PTRAST, 
     *    root, NROW1, NCOL_TO_SEND, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDA, 
     *    ROOT_CONT_STATIC, MYID, COMM,
     *    
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8,
     *    .FALSE., ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
       IF ( IFLAG < 0 ) GOTO 600
       IF (NELIM1.EQ.0) THEN
        CALL SMUMPS_314( N, INODE,  
     *    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, ITYPE2
     $       )
       ENDIF
      ENDIF
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_263
      SUBROUTINE SMUMPS_179(MYID, N, STEP, FRERE, FILS,
     *     NA, LNA, NE, DAD, ND, PROCNODE, SLAVEF,
     *     NRLADU, NIRADU, NIRNEC, NRLNEC,
     *     NRLNEC_ACTIVE,
     *     MAXFR, OPSA,
     *     KEEP,KEEP8, LOCAL_M, LOCAL_N,
     *     SBUF_SEND, SBUF_REC, OPS_SUBTREE, NSTEPS,
     *     NMB_PAR2, ISTEP_TO_INIV2, CANDIDATES, 
     *     IFLAG, IERROR
     $     ,MAX_FRONT_SURFACE_LOCAL
     *     ,MAX_SIZE_FACTOR
     $     )
      IMPLICIT NONE
      INTEGER  MYID, N, LNA, IFLAG, IERROR
      INTEGER  NRLADU, NIRADU, NRLNEC, NIRNEC
      INTEGER NRLADU_CURRENT, NRLADU_ROOT_3, NRLNEC_ACTIVE
      INTEGER  MAXFR, NSTEPS
      INTEGER STEP(N)
      INTEGER FRERE(NSTEPS), FILS(N), NA(LNA), NE(NSTEPS),
     *        ND(NSTEPS), PROCNODE(NSTEPS), DAD(NSTEPS)
      INTEGER  SLAVEF, KEEP(500), LOCAL_M, LOCAL_N
      INTEGER*8 KEEP8(150)
      INTEGER  SBUF_SEND, SBUF_REC
      INTEGER  NMB_PAR2
      INTEGER  ISTEP_TO_INIV2( KEEP(71) )
      INTEGER  CANDIDATES( SLAVEF+1, MAX(NMB_PAR2,1) )
      INTEGER  IDUMMY
      INTEGER MAX_SIZE_FACTOR
      REAL OPSA
      REAL OPS_SUBTREE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LSTKR, TNSTK, IPOOL, 
     *                                      LSTKI 
      INTEGER I,ISTKR,ITOP,NELIM,NFR
      INTEGER ISTKI, ISTKIM, STKI, INDX
      INTEGER K,LSTK,NSTK, IFATH
      INTEGER INODE, LEAF, NBLEAF, NBROOT, III, IN
      INTEGER LEVEL, MAXTEMPCB, MAXITEMPCB
      LOGICAL UPDATE, UPDATEF, MASTER, MASTERF 
      INTEGER LEVELF, NCB, SIZECB, SIZECBI, SIZECBINFR
      INTEGER NBROWMAX, NSLAVES, TEMP, NSLAVES_LOC,
     *         NELIMF, NFRF, NCBF, CBMAXR, CBMAXS,
     *         NSLAVESF, NBROWMAXF, NCOLMAX, LKJIB,
     *         LKJIBT, KMAX, NBR
      INTEGER LWK_QR, LIWK_QR, SIZE_ROOT
      INTEGER IROOT
      INTEGER allocok
      LOGICAL ROOT_OWNER, COMPRESSCB
      DOUBLE PRECISION OPS_NODE
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      INTEGER MAX_FRONT_SURFACE_LOCAL
       INTEGER WHAT
      INTRINSIC MIN, INT, REAL
      INTEGER SMUMPS_275, SMUMPS_330, SMUMPS_497
      INTEGER SMUMPS_52
      EXTERNAL SMUMPS_503, SMUMPS_52
      LOGICAL SMUMPS_170, SMUMPS_358
      EXTERNAL SMUMPS_275, SMUMPS_330, 
     *         SMUMPS_170, SMUMPS_358
      logical :: force_cand , concerned, UPDATES, STACKCB, MASTERSON
      integer :: nmb_of_cand,istat, IFSON, LEVELSON
      integer, allocatable :: cand_of_node(:)
      COMPRESSCB=(XSIZE.GE.3)
      MAX_FRONT_SURFACE_LOCAL=0
      MAX_SIZE_FACTOR=0
      ALLOCATE( LSTKR(NSTEPS), TNSTK(NSTEPS), IPOOL(NSTEPS),
     *          LSTKI(NSTEPS) , stat=allocok)
      if (allocok .GT. 0) THEN
        IFLAG  =-7
        IERROR = 4*NSTEPS
        RETURN
      endif
      IDUMMY = 0
      LKJIB = MAX(KEEP(5),KEEP(6))
      TNSTK = NE
      LEAF = NA(1)+1
      IPOOL(1:LEAF-1) = NA(3:3+LEAF-2)
      NBROOT = NA(2)
      ISTKR  = 0
      ISTKI  = 0
      ISTKIM = 0
      OPSA   = 0.0E0
      OPS_SUBTREE = 0.0E0
      NRLADU = 0
      NIRADU = 0
      NRLADU_CURRENT = 0
      NRLADU_ROOT_3 = 0
      NRLNEC_ACTIVE = 0
      NRLNEC = 0
      NIRNEC = 0
      MAXFR  = 0
      ITOP = 0
      MAXTEMPCB  = 0
      MAXITEMPCB = 0
      SBUF_SEND  = 1
      SBUF_REC   = 1
      IF (KEEP(38) .NE. 0 .AND. KEEP(60).EQ.0) THEN
        INODE  = KEEP(38)
        NRLADU = LOCAL_M*LOCAL_N
        NRLADU_ROOT_3 = LOCAL_M*LOCAL_N
        NRLNEC_ACTIVE = NRLADU_CURRENT
         MAX_SIZE_FACTOR=MAX(MAX_SIZE_FACTOR,LOCAL_M*LOCAL_N)
        NRLNEC = NRLADU
        IF (SMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
     *                                       .EQ. MYID) THEN
          NIRADU = XSIZE+6+2*ND(STEP(INODE))
        ELSE
          NIRADU = XSIZE+6
        ENDIF
        NIRNEC = NIRADU
      ENDIF
      IF((KEEP(24).eq.0).OR.(KEEP(24).eq.1)) THEN
         force_cand=.FALSE.           
      ELSE
         force_cand=(mod(KEEP(24),2).eq.0)
      END IF
      III = 1
 90   CONTINUE
        IF (LEAF .NE. 1) THEN
           LEAF = LEAF - 1
           INODE = IPOOL(LEAF)
        ELSE 
           WRITE(MYID+6,*) ' ERROR 1 in file mpima41md '
           CALL SMUMPS_ABORT()
        ENDIF
 95     CONTINUE 
        MASTER = SMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
     *           .EQ. MYID
        LEVEL  = SMUMPS_330(STEP(INODE),PROCNODE,SLAVEF)
        IF (KEEP(60).GT.1) THEN
         IF (MASTER .AND. INODE.EQ.KEEP(38)) THEN
          NIRADU = NIRADU+XSIZE+6+2*ND(STEP(INODE))
         ENDIF
        ENDIF
      UPDATE=.FALSE.
      if(.NOT.force_cand) then
         UPDATE = ( (MASTER.AND.(LEVEL.NE.3) ).OR. LEVEL.EQ.2 )
      else
         if(master.and.(level.ne.3)) then
            UPDATE = .TRUE.
         else if(level.eq.2) then
            if ( SMUMPS_358( MYID, SLAVEF, inode,
     *           NMB_PAR2, ISTEP_TO_INIV2, KEEP(71), STEP, N,
     *           CANDIDATES, KEEP(24) ) ) then
              UPDATE = .TRUE.
            end if
         end if
      end if
        NFR    = ND(STEP(INODE))
        NSTK   = NE(STEP(INODE))
        NELIM = 0 
        IN = INODE
 100    NELIM = NELIM + 1 
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
        IFSON = -IN
        IF (LEVEL.EQ.3) THEN
         IF ( 
     *     KEEP(60).LE.1
     *      ) THEN
           NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR+
     *                 LOCAL_M*LOCAL_N)
           NRLADU_CURRENT = LOCAL_M*LOCAL_N
           NRLNEC_ACTIVE = MAX0(NRLNEC_ACTIVE,NRLADU_ROOT_3 + 
     *                        NRLADU_CURRENT+ISTKR)
         ENDIF
         IF (MASTER) THEN 
          IF (NFR.GT.MAXFR) MAXFR = NFR
         ENDIF
        ENDIF
        NCB      = NFR-NELIM
        SIZECBINFR = NCB*NCB
        IF ((LEVELF.NE.3).AND.(COMPRESSCB)) THEN
          SIZECB    = (NCB*(NCB+1))/2
        ELSE
          SIZECB    = SIZECBINFR
        ENDIF
        IF(KEEP(86).EQ.1)THEN
           IF(MASTER.AND.(.NOT.SMUMPS_170(STEP(INODE),
     $          PROCNODE,SLAVEF)))THEN
              IF(LEVEL.EQ.1)THEN
                 MAX_FRONT_SURFACE_LOCAL=MAX(MAX_FRONT_SURFACE_LOCAL,
     $                NFR*NFR)
              ELSEIF(LEVEL.EQ.2)THEN
                 IF(KEEP(50).EQ.0)THEN
                    MAX_FRONT_SURFACE_LOCAL=MAX(MAX_FRONT_SURFACE_LOCAL,
     $                   NFR*NELIM)
                 ELSE
                    MAX_FRONT_SURFACE_LOCAL=MAX(MAX_FRONT_SURFACE_LOCAL,
     $                   NELIM*NELIM)
#if defined(PARPIV)
                    IF (KEEP(50).EQ.2) THEN
                    MAX_FRONT_SURFACE_LOCAL=MAX(MAX_FRONT_SURFACE_LOCAL,
     $                   NELIM*(NELIM+1))
                    ENDIF
#endif
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
        SIZECBI  = 2* NCB  + 6 + 3
        IF (LEVEL.EQ.2) THEN
          IF (MASTER) THEN
             SBUF_SEND = MAX(SBUF_SEND, NELIM*LKJIB+NELIM+6)
          ELSE if (UPDATE) then
             SBUF_REC = MAX( SBUF_REC, NELIM*LKJIB+NELIM+6 )
             IF (KEEP(50).EQ.1) THEN
               LKJIBT  = LKJIB
             ELSE
               LKJIBT  = NELIM
             ENDIF
             IF (KEEP(48) .EQ. 5) THEN
                  WHAT = 4 ! Compute NBROWMAX only
                  IF (FORCE_CAND) THEN
                    NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     $                  ISTEP_TO_INIV2(STEP(INODE)))
                  ELSE
                    NSLAVES_LOC=SLAVEF-1
                  ENDIF
             ELSE
                  WHAT = 1 ! Compute NBROWMAX only
                  NSLAVES_LOC=SLAVEF
             ENDIF
             CALL SMUMPS_503(WHAT, KEEP,KEEP8,
     &       NCB, NFR, NSLAVES_LOC, NBROWMAX, IDUMMY)
             SBUF_SEND = MAX(SBUF_SEND,
     *                       LKJIBT*NBROWMAX+6)
             SBUF_REC = MAX( SBUF_REC, NBROWMAX*LKJIBT+6 )
          ENDIF
        ENDIF
        IFATH = DAD(STEP(INODE))
        IF ( UPDATE ) THEN
          IF ( (MASTER) .AND. (LEVEL.EQ.1) ) THEN
             NIRADU = NIRADU + 6 + 2*NFR + XSIZE
             NRLADU = NRLADU + NELIM*NFR
             NRLADU_CURRENT = NELIM*NFR
             MAX_SIZE_FACTOR=MAX(MAX_SIZE_FACTOR,NELIM*NFR)
            SIZECBI    = 2* NCB  + 6 + 3
          ELSEIF (LEVEL.EQ.2) THEN
            IF (MASTER) THEN
             NIRADU = NIRADU + 6 + SLAVEF - 1 + 2*NFR + XSIZE
             NRLADU = NRLADU + NELIM*NELIM
             NRLADU_CURRENT = NELIM*NELIM
             MAX_SIZE_FACTOR=MAX(MAX_SIZE_FACTOR,NELIM*NELIM)
             SIZECB     = 0
             SIZECBINFR = 0
             SIZECBI    = NCB + 5 +  SLAVEF - 1
            ELSE
             IF (KEEP(48) .EQ. 5) THEN
               WHAT = 5
               IF (FORCE_CAND) THEN
                 NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     $                 ISTEP_TO_INIV2(STEP(INODE)))
               ELSE
                 NSLAVES_LOC=SLAVEF-1
               ENDIF
             ELSE
               WHAT = 2 ! Both NBROWMAX and SIZECB
               NSLAVES_LOC=SLAVEF
             ENDIF
             CALL SMUMPS_503(WHAT, KEEP,KEEP8,
     &         NCB, NFR, NSLAVES_LOC, NBROWMAX, SIZECB )
             SIZECBINFR = SIZECB
             NIRADU   = NIRADU+4+NELIM+NBROWMAX 
             NRLADU   = NRLADU + NELIM*NBROWMAX
             NRLADU_CURRENT = NELIM*NBROWMAX
             MAX_SIZE_FACTOR=MAX(MAX_SIZE_FACTOR,NELIM*NBROWMAX)
             NSLAVES = SMUMPS_52( SLAVEF,
     *       KEEP(48), KEEP(10), KEEP(50), NFR, NCB)
             SIZECBI      = 4 + NBROWMAX + NCB + NSLAVES + 1
            ENDIF
          ENDIF
          NIRNEC = MAX0(NIRNEC,NIRADU+ISTKI+SIZECBI+MAXITEMPCB)
#if defined(PARPIV)
          IF (MASTER .AND. KEEP(50).EQ.2) THEN
             NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR+SIZECBINFR+
     *            MAX0(SIZECB,MAXTEMPCB)+NELIM )
          ELSE
             NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR+SIZECBINFR+
     *            MAX0(SIZECB,MAXTEMPCB) )
          ENDIF
#else
          NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR+SIZECBINFR+
     *         MAX0(SIZECB,MAXTEMPCB) )
#endif
         NRLNEC_ACTIVE = MAX0(NRLNEC_ACTIVE,NRLADU_CURRENT+
     *                        NRLADU_ROOT_3+ISTKR+SIZECBINFR+
     *                        MAX0(SIZECB,MAXTEMPCB) )
          IF (NFR.GT.MAXFR) MAXFR = NFR
          IF (NSTK.GT.0) THEN
           DO 70 K=1,NSTK
            LSTK = LSTKR(ITOP)
            ISTKR = ISTKR - LSTK
            STKI = LSTKI( ITOP )
            ISTKI = ISTKI - STKI
            ITOP = ITOP - 1
            IF (ITOP.LT.0) THEN
             write(MYID+6,*) ' ERROR 2 in file mpima41md'
             CALL SMUMPS_ABORT()
            ENDIF
   70      CONTINUE
          ENDIF  
      ELSE IF (LEVEL.NE.3) THEN
         DO WHILE (IFSON.GT.0) 
            UPDATES=.FALSE.
            MASTERSON = SMUMPS_275(STEP(IFSON),PROCNODE,SLAVEF)
     *                  .EQ.MYID
            LEVELSON  = SMUMPS_330(STEP(IFSON),PROCNODE,SLAVEF)
            if(.NOT.force_cand) then
               UPDATES =((MASTERSON.AND.(LEVELSON.NE.3)).OR. 
     *                   LEVELSON.EQ.2)
            else
               if(MASTERSON.and.(LEVELSON.ne.3)) then
                  UPDATES = .TRUE.
               else if(LEVELSON.eq.2) then
                  if ( SMUMPS_358( MYID, SLAVEF, IFSON,
     *                 NMB_PAR2, ISTEP_TO_INIV2, KEEP(71), STEP, N,
     *                 CANDIDATES, KEEP(24)) ) then
                    UPDATES = .TRUE.
                  end if
               end if
            end if
            IF (UPDATES) THEN
              LSTK = LSTKR(ITOP)
              ISTKR = ISTKR - LSTK
              STKI = LSTKI( ITOP )
              ISTKI = ISTKI - STKI
              ITOP = ITOP - 1
              IF (ITOP.LT.0) THEN
                write(*,*) MYID,': ERROR 2 in mpima41md. ITOP = ',ITOP
               CALL SMUMPS_ABORT()
              ENDIF
            ENDIF
            IFSON = FRERE(STEP(IFSON)) ! process next son
         END DO
        ENDIF
        IF (MASTER.OR.LEVEL.EQ.3) THEN
         IF (
     *        ( (INODE.NE.KEEP(20)).OR.(KEEP(60).EQ.0) ) 
     *       .AND.
     *        ( (INODE.NE.KEEP(38)).OR.(KEEP(60).LE.1) ) 
     *      )
     *   THEN
            CALL SMUMPS_511(NFR, NELIM, NELIM,KEEP(50),
     *           1,OPS_NODE)
         ELSE
           OPS_NODE = 0.0E0
         ENDIF
         IF ( LEVEL .eq. 3 ) THEN
            OPSA = OPSA + REAL(OPS_NODE) / REAL( SLAVEF )
         ELSE IF (MASTER) THEN
            OPSA = OPSA + REAL(OPS_NODE)
         ENDIF
         IF (SMUMPS_170(STEP(INODE),
     *   PROCNODE, SLAVEF) .OR. NE(STEP(INODE))==0) THEN
           IF (LEVEL == 1) THEN
             OPS_SUBTREE = OPS_SUBTREE + REAL(OPS_NODE)
           ELSE
             CALL SMUMPS_511(NFR, NELIM, NELIM,KEEP(50),
     *           1,OPS_NODE)
             OPS_SUBTREE = OPS_SUBTREE + REAL(OPS_NODE)
           ENDIF
         ENDIF
        ENDIF
        IF (IFATH .EQ. 0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 115
         GOTO 90
        ELSE
         LEVELF = SMUMPS_330(STEP(IFATH),PROCNODE,SLAVEF)
         MASTERF= SMUMPS_275(STEP(IFATH),PROCNODE,SLAVEF)
     *            .EQ. MYID
         UPDATEF=.FALSE.
         if(.NOT.force_cand) then
            UPDATEF= ((MASTERF.AND.(LEVELF.NE.3)).OR.LEVELF.EQ.2)
         else
            if(masterf.and.(levelf.ne.3)) then
               UPDATEF = .TRUE.
            else if(levelf.eq.2) then
               if ( SMUMPS_358( MYID, SLAVEF, ifath,
     *              NMB_PAR2, ISTEP_TO_INIV2, KEEP(71), STEP, N,
     *              CANDIDATES, KEEP(24) )) THEN
                 UPDATEF = .TRUE.
               end if
            end if
         end if
         concerned  = UPDATEF .OR. UPDATE
         NCB    = NFR-NELIM
         IF (LEVELF.EQ.3) THEN
           TEMP       = MIN(SIZECB, LOCAL_M*LOCAL_N) 
           MAXTEMPCB  = MAX(MAXTEMPCB, TEMP)
           MAXITEMPCB = MAX(MAXITEMPCB,SIZECBI)
           TEMP       = MIN(SIZECB, LOCAL_M*LOCAL_N) + SIZECBI
           SBUF_REC   = MAX(SBUF_REC, TEMP)
         ENDIF
        IF (concerned) THEN
         IN = IFATH
         NELIMF = 0
         DO WHILE (IN.GT.0)
            IN = FILS(IN)
            NELIMF = NELIMF+1
         ENDDO
         NFRF = ND(STEP(IFATH))
         NCBF = NFRF - NELIMF
         IF (LEVELF.EQ.2) THEN
           IF (KEEP(48) .EQ. 5) THEN
               WHAT = 4
               IF (FORCE_CAND) THEN
                 NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     $                 ISTEP_TO_INIV2(STEP(IFATH)))
               ELSE
                 NSLAVES_LOC=SLAVEF-1
               ENDIF
             ELSE
               WHAT = 1 ! NBROWMAX only
               NSLAVES_LOC=SLAVEF
           ENDIF
           CALL SMUMPS_503(WHAT,
     *     KEEP,KEEP8, NCBF, NFRF, NSLAVES_LOC, NBROWMAXF,IDUMMY)
           NCOLMAX  = NCB
           IF (UPDATE) THEN
              IF(MASTERF)THEN
                 NBR = MIN(NBROWMAXF,NCB)
              ELSE
                 NBR = MAX(MIN(NELIMF,NCB),MIN(NBROWMAXF,NCB))
              ENDIF
              CBMAXS = NBR*NCOLMAX - (NBR*(NBR-1))/2
           ELSE
              CBMAXS = 0
           END IF
           IF (MASTERF) THEN
                CBMAXR = MIN(NELIMF,NCB)*MIN(NCOLMAX,NELIMF)
                IF (LEVEL.EQ.1) THEN
                  IF (.NOT.COMPRESSCB) SIZECB = (NCB*(NCB+1))/2
                  CBMAXR = MIN(CBMAXR, SIZECB)
                ENDIF
           ELSE IF (UPDATEF) THEN
                 NBR = MIN(NBROWMAXF,NCB)
                 IF (LEVEL.EQ.1) THEN
                  CBMAXS = NBR*NCOLMAX - (NBR*(NBR-1))/2
                 ELSE
                  CBMAXR = NBR*NCOLMAX
                 ENDIF
           else
              cbmaxr = 0
           ENDIF
         ELSEIF (LEVELF.EQ.3) THEN
           CBMAXR = LOCAL_M*LOCAL_N
           CBMAXS = LOCAL_M*LOCAL_N
         ELSE
           IF (LEVEL.EQ.1) THEN
             CBMAXR = (NFRF*(NFRF+1))/2
             CBMAXS = (NFRF*(NFRF+1))/2
            ELSE
             CBMAXR = NFRF*NFRF
             CBMAXS = NFRF*NFRF
            ENDIF
         ENDIF
         IF (UPDATE) THEN
           SIZECB = MIN(CBMAXS, SIZECB)
           IF ( .not. ( LEVELF .eq. 1 .AND. UPDATEF ) )THEN
              SBUF_SEND = MAX(SBUF_SEND, SIZECB+SIZECBI)
           ENDIF
         ENDIF
         STACKCB = .FALSE.
         IF (UPDATEF) THEN 
          STACKCB = .TRUE.
          SIZECBI = 2 * NFR + 6 + XSIZE
          IF (LEVEL.EQ.1) THEN
             IF ((LEVELF.NE.3).AND.COMPRESSCB) THEN
               SIZECB  = (NCB*(NCB+1))/2
             ELSE
               SIZECB  = NCB*NCB
             ENDIF
             IF (MASTER) THEN
               SIZECBI = 2+ XSIZE
             ELSE IF (LEVELF.EQ.1) THEN
              SIZECBI = 2 * NCB + 6 + 3+ XSIZE
              SIZECB  = MIN(CBMAXR,SIZECB)
              SBUF_REC = MAX(SBUF_REC, SIZECBI+SIZECB)
             ELSE 
              SIZECB = MIN(SIZECB,CBMAXR)
              MAXTEMPCB  = MAX(MAXTEMPCB, SIZECB)
              MAXITEMPCB = MAX(MAXITEMPCB, 2 * NCB + 6 + 3)
              SBUF_REC = MAX(SBUF_REC, SIZECB +2 * NCB + 6 + 3)
              SIZECBI= 0
              SIZECB = 0
             ENDIF
          ELSE 
             IF (KEEP(48) .EQ. 5) THEN
               WHAT = 5! Both NBROWMAX and SIZECB max are computed
               IF (FORCE_CAND) THEN
                 NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     $                 ISTEP_TO_INIV2(STEP(INODE)))
               ELSE
                 NSLAVES_LOC=SLAVEF-1
               ENDIF
             ELSE
               WHAT = 2 ! Both NBROWMAX and SIZECB max are computed
               NSLAVES_LOC=SLAVEF
             ENDIF
             CALL SMUMPS_503( WHAT, KEEP,KEEP8,
     &       NCB, NFR, NSLAVES_LOC, NBROWMAX, SIZECB )
             SIZECB  = MIN(CBMAXR,SIZECB)
             MAXTEMPCB  = MAX(MAXTEMPCB, SIZECB)
             MAXITEMPCB = MAX(MAXITEMPCB, NBROWMAX + NCB + 6 )
             SBUF_REC = MAX(SBUF_REC,SIZECB+NBROWMAX + NCB + 6)
             IF (MASTER) THEN
              SIZECBI =  NCB + 5 +  SLAVEF - 1+ XSIZE
              SIZECB  = 0
             ELSE IF (UPDATE) THEN
              SIZECBI =  NFR + 6 + SLAVEF - 1 + XSIZE
              NSLAVES = SMUMPS_52(
     &            SLAVEF, KEEP(48), KEEP(10), KEEP(50),
     &            NFR, NCB)
              SIZECBI = SIZECBI+NBROWMAX+NFR+6+NSLAVES+1+XSIZE
             ELSE
              SIZECB  = 0 
              SIZECBI = 0
             ENDIF
          ENDIF
         ELSE
           IF (LEVELF.NE.3) THEN
               STACKCB = .TRUE.
               SIZECB  = 0 
               SIZECBI = 0
               IF ( (LEVEL.EQ.1) .AND. (LEVELF.NE.1) ) THEN
                  IF (COMPRESSCB) THEN 
                      SIZECB  = (NCB*(NCB+1))/2
                  ELSE
                      SIZECB  = NCB*NCB
                  ENDIF
                  SIZECBI = 2 * NCB + 6 + 3 + XSIZE
               ELSE IF (LEVEL.EQ.2) THEN
                IF (MASTER) THEN
                 SIZECBI =  NCB + 5 +  SLAVEF - 1 + XSIZE
                ELSE 
                 IF (KEEP(48) .EQ. 5) THEN
                  WHAT = 5! Both NBROWMAX and SIZECB max are computed
                  IF (FORCE_CAND) THEN
                    NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     $                  ISTEP_TO_INIV2(STEP(INODE)))
                  ELSE
                    NSLAVES_LOC=SLAVEF-1
                  ENDIF
                 ELSE
                  WHAT = 2 ! Both NBROWMAX and SIZECB max are computed
                  NSLAVES_LOC=SLAVEF
                 ENDIF
                 CALL SMUMPS_503(WHAT, KEEP,KEEP8,
     &           NCB, NFR, NSLAVES_LOC, NBROWMAX, SIZECB )
                 SIZECBI =  NFR + 6 + SLAVEF - 1 + XSIZE
                 SIZECBI = SIZECBI + NBROWMAX + NFR + 6 + XSIZE
                ENDIF 
               ENDIF
           ENDIF
         ENDIF
         IF (STACKCB) THEN
          IF (FRERE(STEP(INODE)).EQ.0) THEN
            write(*,*) ' ERROR 3 in LDLT_mpima41md'
            CALL SMUMPS_ABORT()
          ENDIF
           ITOP = ITOP + 1
           LSTKI(ITOP) = SIZECBI
           ISTKI=ISTKI+LSTKI(ITOP)
           ISTKIM=MAX( ISTKI, ISTKIM )
           LSTKR(ITOP) = SIZECB
           ISTKR = ISTKR + LSTKR(ITOP)
           NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR+MAXTEMPCB)
           NIRNEC = MAX0(NIRNEC,NIRADU+ISTKI+MAXITEMPCB)
         ENDIF ! end of (stackcb)
        ENDIF ! end of (concerned) 
         TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
         IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
         ELSE
            GOTO 90
         ENDIF
        ENDIF
  115 CONTINUE
      ISTKIM = ISTKIM / 2 + MAXFR + N
      IF (KEEP(38) .ne. 0 ) THEN
        IF (SMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
     *                                       .EQ. MYID) THEN
          ISTKIM = ISTKIM + ND(STEP( KEEP(38) ))
        END IF
      END IF
      IF ( NRLNEC - NRLADU .LT. ISTKIM ) THEN
        NRLNEC = NRLADU + ISTKIM
      END IF
      IF ( NIRNEC - NIRADU .LT. ISTKIM ) THEN
        NIRNEC = NIRADU + ISTKIM
      END IF
      IF ( KEEP(53) .NE. 0 ) THEN
        IF ( KEEP(38) .ne. 0 ) THEN
          IROOT = KEEP( 38 )
        ELSE
          IROOT = KEEP( 20 )
        END IF
        ROOT_OWNER = ( MYID .eq.
     *  SMUMPS_275( STEP(IROOT), PROCNODE, SLAVEF ) )
        SIZE_ROOT = ND(STEP(IROOT))
        CALL SMUMPS_279( KEEP(51), KEEP(51), SIZE_ROOT,
     *       LOCAL_M, LOCAL_N, ROOT_OWNER, KEEP,KEEP8,
     *       LIWK_QR, LWK_QR )
        IF ( NRLNEC-NRLADU .LT. LWK_QR ) THEN
          NRLNEC = NRLADU + LWK_QR
        END IF
        IF ( NIRNEC-NIRADU .LT. LIWK_QR ) THEN
          NIRNEC = NIRNEC + LIWK_QR
        END IF
      END IF
      NRLNEC = MAX(NRLNEC, NRLADU+4*KEEP(127)*KEEP(84))
      NRLNEC_ACTIVE = MAX(NRLNEC_ACTIVE, MAX_SIZE_FACTOR+
     *                    4*KEEP(127)*KEEP(84))
      SBUF_REC = MAX(SBUF_REC, MAXTEMPCB+MAXITEMPCB )
      SBUF_REC = SBUF_REC   + 10
      SBUF_SEND = SBUF_SEND + 10
#if defined(PARPIV)
      IF(KEEP(50) .EQ. 2) THEN
         SBUF_REC = SBUF_REC+KEEP(108)+1
         SBUF_SEND = SBUF_SEND+KEEP(108)+1
      ENDIF
#endif
      IF (SLAVEF.EQ.1) THEN 
         SBUF_REC = 1
         SBUF_SEND= 1
      ENDIF
      DEALLOCATE( LSTKR, TNSTK, IPOOL,
     *          LSTKI )
      RETURN
      END SUBROUTINE SMUMPS_179
      SUBROUTINE SMUMPS_137( INODE, N, PROCNODE_STEPS,
     *           SLAVEF,
     &           ND, FILS, FRERE_STEPS, STEP, PIMASTER,
     *           KEEP28, KEEP50,
     *           FLOP1,
     &           IW, LIW )
      IMPLICIT NONE
      INTEGER INODE, N, KEEP50, LIW, SLAVEF, KEEP28
      INTEGER PROCNODE_STEPS(KEEP28), ND(KEEP28),
     &        FILS(N), FRERE_STEPS(KEEP28),
     *        STEP(N), 
     * PIMASTER(KEEP28),
     *  IW( LIW )
      DOUBLE PRECISION FLOP1
      INTEGER NUMORG, IN, NASS, IFSON, NUMSTK, NFRONT, NPIV, NCB, 
     &        LEVEL, ISON
      LOGICAL SMUMPS_170
      INTEGER SMUMPS_330
      EXTERNAL SMUMPS_170, SMUMPS_330
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0E0)
      INTEGER XSIZE
      PARAMETER (XSIZE=3)
      FLOP1 = 0.0E0
      IF (SMUMPS_170(STEP(INODE),
     *    PROCNODE_STEPS, SLAVEF) ) RETURN
      IN     = INODE
      NUMORG = 0
   10 NUMORG = NUMORG + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 10
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      IF (ISON .EQ. 0) GOTO 30
   20 NUMSTK = NUMSTK + 1
      NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 +XSIZE)
      ISON = FRERE_STEPS(STEP(ISON))
      IF (ISON .GT. 0) GOTO 20
   30 NFRONT = ND(STEP(INODE)) + NASS
      NPIV  = NASS + NUMORG
      NCB   = NFRONT - NPIV
      LEVEL = SMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      CALL SMUMPS_511(NFRONT,NPIV,NPIV,KEEP50,LEVEL,FLOP1)
      RETURN
      END SUBROUTINE SMUMPS_137
      SUBROUTINE SMUMPS_511(NFRONT,NPIV,NASS,
     *                                 KEEP50,LEVEL,COST)
      IMPLICIT NONE
      INTEGER NFRONT,NPIV,KEEP50,LEVEL, NASS
      DOUBLE PRECISION COST
      IF (KEEP50.EQ.0) THEN
        IF (LEVEL.EQ.1) THEN
          COST = DBLE(2) * DBLE(NFRONT) * DBLE(NPIV) *
     *      DBLE(NFRONT - NPIV - 1) +
     *      DBLE(NPIV) * DBLE(NPIV + 1) * DBLE(2 * NPIV + 1)
     *          / DBLE(3)
          COST = COST + DBLE(2 * NFRONT - NPIV - 1)
     *      * DBLE(NPIV) /DBLE(2)
        ELSEIF (LEVEL.EQ.2) THEN
          COST = DBLE(2*NASS)*DBLE(NFRONT) -
     *          DBLE(NASS+NFRONT)*DBLE(NPIV+1)
          COST = DBLE(NPIV)*COST +
     *     DBLE(2 * NASS - NPIV - 1) * DBLE(NPIV) / DBLE(2) +
     *     DBLE(NPIV) * DBLE(NPIV + 1) *
     *     DBLE(2 * NPIV + 1) /DBLE(3)
        ENDIF
      ELSE
        IF (LEVEL.EQ.1) THEN
          COST = DBLE(NPIV) * (
     *          DBLE( NFRONT ) * DBLE( NFRONT ) +
     *          DBLE( NFRONT ) - (
     *          DBLE( NFRONT)*DBLE(NPIV) + DBLE(NPIV+1)
     *          )) +( DBLE(NPIV)*DBLE(NPIV+1)
     *          *DBLE(2*NPIV+1))/ DBLE(6)
        ELSE IF (LEVEL.EQ.3.AND.KEEP50.EQ.2) THEN
          COST = DBLE(2) * DBLE(NFRONT) * DBLE(NPIV) *
     *      DBLE(NFRONT - NPIV - 1) +
     *      DBLE(NPIV) * DBLE(NPIV + 1) *
     *      DBLE(2 * NPIV + 1) / DBLE(3)
          COST = COST + DBLE(2 * NFRONT - NPIV - 1)
     *         * DBLE(NPIV) / DBLE(2)
        ELSE
          COST = DBLE(NPIV) * (
     *          DBLE( NASS ) * DBLE( NASS ) + DBLE( NASS )
     *        - DBLE( NASS*NPIV +  NPIV + 1 )    )
     *        + ( DBLE(NPIV)*DBLE(NPIV+1)*DBLE(2*NPIV+1) )
     *        / DBLE( 6 )
        ENDIF
      ENDIF
      END SUBROUTINE SMUMPS_511
