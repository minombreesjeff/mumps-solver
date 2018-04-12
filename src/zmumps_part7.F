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
#if defined(try_null_space)
      SUBROUTINE ZMUMPS_295
     *(NN,A,ROOT,WR03,LWR03,KEEP,KEEP8,RIGHT_SV,INFO,LP)
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      INCLUDE 'mpif.h'
      TYPE ( ZMUMPS_ROOT_STRUC ) :: root
      INTEGER NN,LWR03,LP
      COMPLEX*16 A(NN*NN)
      INTEGER INFO(2),KEEP(500)
      INTEGER*8 KEEP8(150)
      COMPLEX*16 WR03(LWR03),RIGHT_SV(NN)
      INTEGER LDLT,QR,DEFICIENCY
      COMPLEX*16, DIMENSION(:), POINTER :: GROW, GCOS, GSIN
      INTEGER, DIMENSION(:), POINTER :: PERM_ROW
      INTEGER, DIMENSION(:), POINTER :: PERM_COL
      COMPLEX*16, DIMENSION(:), POINTER :: ELELG
      COMPLEX*16, DIMENSION(:), POINTER :: EUELG
      INTEGER, DIMENSION(:), POINTER :: ELROW
      INTEGER, DIMENSION(:), POINTER :: EUROW
      INTEGER, DIMENSION(:), POINTER :: PTREL
      INTEGER, DIMENSION(:), POINTER :: PTREU
      COMPLEX*16, DIMENSION(:), POINTER :: DL
      COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: AL
      COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: AU
      INTEGER EUIND,ELIND,DIML,DIMU,ROT_STEPS,KK
      INTEGER LSHIFT,USHIFT,LUSHIFT,IND_ROT,L,IND1,IND2
      INTEGER JJ,KND,NULL_IND
      DOUBLE PRECISION XMAX,YMAX
      DOUBLE PRECISION TEMP,SV_MU,ABSYL,ABSXK
      COMPLEX*16 TEMP1,TEMP2,ELG,SMIN
      DOUBLE PRECISION ABSVK
      INTEGER ELIND_SAVE,DIMLU,NLUPDATE,NUUPDATE,NULL_MAX
      INTEGER ELG_MAX,NUM_RANK,EST_DEFICIENCY
      INTEGER I,J,K,LDA,GIND,IERR,IND,II,JND,JND1
      INTEGER MAXG,IND_FROM,IND_TO,NSHORT
      INTEGER SMP_DEFICIENCY,ALLOCOK,MAXDEF
      LOGICAL GAP_FIND_IS_ON
      COMPLEX*16 ALPHA,SUM,VMAX, ONE
      DOUBLE PRECISION TOL,EPS,GAP_SIZE,AVG
      COMPLEX*16 C,S
      DOUBLE PRECISION TOL16,TOL13,TOL11,ZERO,TOL9
      DOUBLE PRECISION TOL8,TOL12,TOL10
      PARAMETER(TOL9=1.0D-9,TOL8=1.0D-8,TOL12=1.0D-12,TOL10=1.0D-10)
      PARAMETER(TOL16=1.0D-16,TOL13=1.0D-13,TOL11=1.1D-11)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      PARAMETER(SV_MU=0.1D0,EPS=1.1D-16,GAP_SIZE=1.0D2)
      LOGICAL SPECDEF
      COMPLEX*16 ZMUMPS_104,ZMUMPS_106
      IF (  LWR03 .LT. 3 * NN ) THEN
        INFO(1) = -9
        INFO(2) = 3 * NN - LWR03
        IF ( LP .GT. 0 )
     *  WRITE(LP,*) 'ZMUMPS_295: LWR03 too small'
        RETURN
      END IF
      LDLT=KEEP(50)
      QR=KEEP(19)
      IF ( QR .eq. 9 .OR. QR .eq. 10 ) THEN
        QR = QR-8
        SPECDEF = .TRUE.
      ELSE
        SPECDEF = .FALSE.
      END IF
      MAXDEF=KEEP(21)
      IF(MAXDEF.LE.0 .AND. .NOT. SPECDEF) MAXDEF=NN
      IF(MAXDEF.GT.NN)THEN
          MAXDEF=NN
          WRITE(*,*) 'Warning: deficiency no larger than ', NN
      END IF
      IF (LDLT.NE.0) THEN
        CALL ZMUMPS_297(NN,A)
      END IF
      DEFICIENCY=0
      KEEP(17)=DEFICIENCY
      LDA=NN
      IF(QR.EQ.1.OR.QR.EQ.2.OR.QR.EQ.3.OR.QR.EQ.4) THEN
        GIND=0
        root%GIND=GIND
        IF(ASSOCIATED(root%QR_TAU)) DEALLOCATE(root%QR_TAU)
        ALLOCATE(root%QR_TAU(NN),stat=IERR )
        IF(IERR.GT.0) THEN
          INFO(1)=-13
          INFO(2)=NN
          IF ( LP .GT. 0 )
     *    WRITE(LP,*) 'ZMUMPS_295: Problem allocating QR/QP'
          RETURN
        END IF
        IF (ASSOCIATED(root%IPIV)) DEALLOCATE(root%IPIV)
        ALLOCATE(root%IPIV(NN), stat=IERR)
        IF (IERR.GT.0) THEN
          INFO(1)=-13
          INFO(2)=NN
          IF ( LP .GT. 0 )
     *    WRITE(LP,*) 'ZMUMPS_295:Problem allocating IPIV'
          RETURN
        END IF
        IF(QR.EQ.1.OR.QR.EQ.2) then
          root%IPIV = 0
          CALL ZGEQPF(NN,NN,A(1),LDA,root%IPIV,
     *      root%QR_TAU(1),WR03,IERR)
          IF(IERR.NE.0) THEN
            WRITE(*,*) ' Problem in ZGEQPF : IERR = ', IERR
            CALL ZMUMPS_ABORT()
          END IF
          IF(MAXDEF.EQ.0) GO TO 5488
          TOL=1.1*SQRT(DBLE(NN))*EPS*ZMUMPS_104(NN,A,RIGHT_SV)
          TEMP=A(LDA*NN)
          TEMP=ABS(A(LDA*NN))
          IF(TEMP.GE.10*TOL) THEN
            GAP_FIND_IS_ON=.FALSE.
          ELSE
            GAP_FIND_IS_ON=.TRUE.
            AVG=0.02*TEMP
            SUM= DCMPLX(ZERO)
            GO TO 132
          END IF 
          TEMP1=ABS(A(LDA*NN))
          DO I=1,MIN(NN-1,23)
            TEMP=A(LDA*NN-I*(NN+1))
            TEMP=ABS(A(LDA*NN-I*(NN+1)))
            IF(ABS(TEMP1).LE.10E-4*TEMP) THEN
              GAP_FIND_IS_ON=.TRUE.
              AVG=0.02*ABS(A(LDA*NN))
              SUM= DCMPLX(ZERO)
              GO TO 132
            END IF
            TEMP1=TEMP
          END DO
          GAP_FIND_IS_ON=.FALSE.
 132      CONTINUE
          IF(GAP_FIND_IS_ON) THEN
            DO IND=0,NN-1
              TEMP=A(LDA*NN-IND*(NN+1))
              IF(ABS(TEMP).GT.AVG*GAP_SIZE+EPS) GO TO 1487
              SUM=SUM+ABS(TEMP)
              AVG=SUM/(IND+1)
              IF(IND.GE.MAXDEF) GOTO 1487
            END DO
          END IF
 1487     CONTINUE
          DEFICIENCY=IND
          IF(GAP_FIND_IS_ON) GO TO 5488
          DEFICIENCY=0
 1488     CONTINUE
          IND=(NN-DEFICIENCY)*LDA-DEFICIENCY
          ALPHA=A(IND)
          TOL=1.8*SQRT(DBLE(NN))*EPS*ZMUMPS_104(NN,A,RIGHT_SV)
          IF(
     *          (.NOT. SPECDEF .AND. ABS(ALPHA).LE.TOL )
     *     .OR. ( SPECDEF .AND. DEFICIENCY .NE. MAXDEF )
     *     ) THEN
            DEFICIENCY=DEFICIENCY+1
            IF(DEFICIENCY.LT.MAXDEF) GO TO 1488
          END IF
 5488     CONTINUE
          IF(DEFICIENCY.GT.0) THEN
            CALL ZTRTRS('U','N','N',NN-DEFICIENCY,DEFICIENCY,
     *      A,LDA,A(LDA*(NN-DEFICIENCY)+1),LDA,IERR)
            IF ( IERR .NE. 0 ) THEN
              WRITE(*,*) ' Problem in ZTRTRS : IERR = ', IERR
              CALL ZMUMPS_ABORT()
            END IF
          END IF
          root%GIND = GIND
#ifdef try_null_space_not_tested
        ELSEIF(QR.EQ.3.OR.QR.EQ.4) THEN
          CALL ZMUMPS_175(NN,0,ROOT%IPIV,1)
          CALL ZGEQPF(NN,NN,A,LDA,root%IPIV(1),
     *      root%QR_TAU,WR03,IERR)
          IF(IERR.NE.0) THEN
            WRITE(*,*) ' Problem in ZGEQPF : IERR = ', IERR
            CALL ZMUMPS_ABORT()
          END IF
          CALL ZMUMPS_117(3*NN,ZERO,WR03,1)
          GIND=0
          TOL=SQRT(DBLE(NN))*EPS*ZMUMPS_104(NN,A,RIGHT_SV)
          DEFICIENCY=0
          CALL ZMUMPS_98(NN,A,NN,SMP_DEFICIENCY,TOL13)
          NSHORT=NN-SMP_DEFICIENCY
          IND_FROM=LDA*(NN-SMP_DEFICIENCY)+1
          CALL ZTRSM('Left','Upper','No transpose',
     *      'Non-unit',NSHORT,
     *      SMP_DEFICIENCY,ONE,A,LDA,A(IND_FROM),NN)
          I=NN
          NULL_IND=1+(I-1)*LDA
          IND=NN*LDA
          JND=(SMP_DEFICIENCY-1)*NN+1
          DO II=1,SMP_DEFICIENCY
            IND=IND-LDA-1
            I=I-1
          END DO
          NULL_IND=1+(I-1)*LDA
          NULL_IND=NULL_IND-SMP_DEFICIENCY*LDA
          DEFICIENCY=SMP_DEFICIENCY
          I=NN-SMP_DEFICIENCY
          CALL ZMUMPS_317(NN-SMP_DEFICIENCY,LDA,A,WR03(1),
     *     WR03(NN+1),WR03(2*NN+1),SMIN,RIGHT_SV)
          IF(ABS(SMIN).GT.TOL) THEN
            GO TO 211
          END IF
          root%MAXG = 5000
          IF (ASSOCIATED(root%GROW)) DEALLOCATE(root%GROW)
          ALLOCATE(root%GROW(root%MAXG),stat=allocok)
          IF(allocok.GT.0) THEN
            INFO(1) = -13
            INFO(2) = NN
            IF (LP.GT.0) WRITE(LP,*)
     *      'ZMUMPS_295: Problem allocating GROW'
            RETURN
          ENDIF
          IF (ASSOCIATED(root%GCOS)) DEALLOCATE(root%GCOS)
          ALLOCATE(root%GCOS(root%MAXG),stat=allocok)
          IF(allocok.GT.0) THEN
            INFO(1) = -13
            INFO(2) = NN
            IF (LP.GT.0) WRITE(LP,*)
     *      'ZMUMPS_295: Problem allocating GCOS'
            RETURN
          ENDIF
          IF (ASSOCIATED(root%GSIN)) DEALLOCATE(root%GSIN)
          ALLOCATE(root%GSIN(root%MAXG),stat=allocok)
          IF(allocok.GT.0) THEN
            INFO(1) = -13
            INFO(2) = NN
            IF (LP.GT.0) WRITE(LP,*)
     *      'ZMUMPS_295: Problem allocating GSIN'
            RETURN
          ENDIF
          MAXG = root%MAXG
          GROW => root%GROW
          GCOS => root%GCOS
          GSIN => root%GSIN
 111      CONTINUE
          DEFICIENCY=DEFICIENCY+1
          I=I-1
          CALL ZCOPY(I,RIGHT_SV,1,A(NULL_IND),1)
          NULL_IND=NULL_IND-LDA
          CALL ZTRTRS('U','N','N',
     *      NN-DEFICIENCY,1,
     *      A,LDA,A(1+LDA*(NN-DEFICIENCY)),
     *      LDA,IERR)
          VMAX=ZMUMPS_106(I+1,RIGHT_SV)
          DO K=I+1,1,-1
            ABSVK=ABS(RIGHT_SV(K))
            IF(ABS(VMAX)*SV_MU.LT.ABSVK) GO TO 50
          END DO
 50       CONTINUE
          IF(K.LT.I+1) THEN
            KND=(K-1)*LDA+1
            CALL ZCOPY(K-1,A((K-1)*LDA+1),1,WR03,1)
            DO II=K+1,I+1
              CALL ZCOPY(K-1,A(KND+LDA),1,A(KND),1)
              KND=KND+LDA
            END DO
            CALL ZCOPY(K-1,WR03,1,A(I*LDA+1),1)
            CALL ZMUMPS_171(K,I+1,root%IPIV)
            IND=1+LDA*(I+1)
            DO II=I+2,NN
              IND1=K
              IND2=I
              CALL ZMUMPS_105(IND1,IND2,A(IND))
              IND=IND+LDA
            END DO
            IF(MAXG.LT.GIND+I+1-K) CALL
     *         ZMUMPS_298(0,2,'SPARSLAB','SLABQR1',
     *        ' Missing space to store Givens rotations',1601,2)
            DO J=K,I
             JND1=(J-1)*NN+J
              JND=J*NN+J
            CALL ZMUMPS_159(A(JND),A(JND+1),C,S)
            TEMP1=A(JND1)
            TEMP2= DCMPLX(ZERO)
            CALL ZMUMPS_31(I+1-J,C,S,A(JND),NN,A(JND+1),NN)
            CALL ZMUMPS_31(1,C,S,TEMP1,1,TEMP2,1)
            KND=(J-1)*NN+J
            A(KND)=A(KND+LDA)
            KND=KND+LDA
            DO II=J+2,I+1
              CALL ZCOPY(2,A(KND+LDA),1,A(KND),1)
              KND=KND+LDA
            END DO
            KND=(I-1)*LDA+J
            A(KND)=TEMP1
            A(KND+1)=TEMP2
            GIND=GIND+1
            GROW(GIND)=J
            GCOS(GIND)=C
            GSIN(GIND)=S
            END DO
          END IF
          CALL ZMUMPS_317(I,LDA,A,WR03(1),
     *     WR03(NN+1),WR03(2*NN+1),SMIN,RIGHT_SV)
          IF(ABS(SMIN).LT.TOL.AND.NN-I.LE.MAXDEF) THEN
            GO TO 111
          END IF
          NUM_RANK=I
 211      CONTINUE
          root%GIND = GIND
#endif
        END IF
#ifdef try_null_space_not_tested
      ELSEIF(QR.EQ.5.OR.QR.EQ.6) THEN
        IF (ASSOCIATED(root%IPIV)) DEALLOCATE(root%IPIV)
        ALLOCATE(root%IPIV(NN), stat=IERR)
        IF (IERR.GT.0) THEN
          INFO(1)=-13
          INFO(2)=NN
          WRITE(*,*) 'Problem allocating IPIV in root'
          RETURN
        END IF
        CALL ZMUMPS_175(NN,0,root%IPIV(1),1)
        CALL ZGETRF(NN,NN,A,LDA,root%IPIV(1),IERR)
        DEFICIENCY=0
        TOL=10*SQRT(DBLE(NN))*EPS*ZMUMPS_104(NN,A,RIGHT_SV)
        CALL ZMUMPS_98(NN,A,NN,SMP_DEFICIENCY,TOL)
        DEFICIENCY=0
        CALL ZMUMPS_98(NN,A,NN,SMP_DEFICIENCY,TOL)
        SMP_DEFICIENCY=MIN(SMP_DEFICIENCY,MAXDEF)
        NSHORT=NN-SMP_DEFICIENCY
        IND_FROM=LDA*(NN-SMP_DEFICIENCY)+1
        CALL ZTRSM('Left','Upper','No transpose',
     *    'Non-unit',NSHORT,
     *    SMP_DEFICIENCY,ONE,A,LDA,A(IND_FROM),NN)
        I=NN
        NULL_IND=1+(I-1)*LDA
        IND=NN*LDA
        JND=(SMP_DEFICIENCY-1)*NN+1
        DO II=1,SMP_DEFICIENCY
          IND=IND-LDA-1
          I=I-1
        END DO
        NULL_IND=1+(I-1)*LDA
        NULL_IND=NULL_IND-SMP_DEFICIENCY*LDA
        DEFICIENCY=SMP_DEFICIENCY
        I=NN-SMP_DEFICIENCY
      ELSEIF(QR.EQ.7.OR.QR.EQ.8) THEN
        IF (ASSOCIATED(root%IPIV)) DEALLOCATE(root%IPIV)
        ALLOCATE(root%IPIV(NN), stat=IERR)
        IF (IERR.GT.0) THEN
          INFO(1)=-13
          INFO(2)=NN
          WRITE(*,*) 'Problem allocating IPIV in root'
          RETURN
        END IF
        CALL ZMUMPS_175(NN,0,root%IPIV(1),1)
        CALL ZGETRF(NN,NN,A,LDA,root%IPIV(1),IERR)
        DEFICIENCY=0
        CALL ZMUMPS_98(NN,A,NN,SMP_DEFICIENCY,TOL13)
        TOL=SQRT(DBLE(NN))*EPS*ZMUMPS_104(NN,A,RIGHT_SV)
        TOL=TOL8
        DEFICIENCY=0
        CALL ZMUMPS_98(NN,A,NN,SMP_DEFICIENCY,TOL10)
        SMP_DEFICIENCY=MIN(SMP_DEFICIENCY,MAXDEF)
        NSHORT=NN-SMP_DEFICIENCY
        IND_FROM=LDA*(NN-SMP_DEFICIENCY)+1
        CALL ZTRSM('Left','Upper','No transpose',
     *    'Non-unit',NSHORT,
     *    SMP_DEFICIENCY,ONE,A,LDA,A(IND_FROM),NN)
        I=NN
        NULL_IND=1+(I-1)*LDA
        IND=NN*LDA
        JND=(SMP_DEFICIENCY-1)*NN+1
        DO II=1,SMP_DEFICIENCY
          IND=IND-LDA-1
          I=I-1
        END DO
        NULL_IND=1+(I-1)*LDA
        NULL_IND=NULL_IND-SMP_DEFICIENCY*LDA
        DEFICIENCY=SMP_DEFICIENCY
        I=NN-SMP_DEFICIENCY
        EST_DEFICIENCY=1
        root%ELG_MAX = 5000
        ELG_MAX=root%ELG_MAX
        root%NULL_MAX = 20
        NULL_MAX=root%NULL_MAX
        IF (ASSOCIATED(root%PERM_COL)) DEALLOCATE(root%PERM_COL)
        ALLOCATE(root%PERM_COL(NN),stat=allocok)
        IF(allocok.GT.0) THEN
          INFO(1) = -13
          INFO(2) = NN
          IF ( LP .GT. 0 )
     *      WRITE(LP,*)
     *      ' ZMUMPS_295: Problem allocating PERM_COL'
          RETURN
        ENDIF
        CALL ZMUMPS_174(NN,0,root%PERM_COL)
        PERM_COL => root%PERM_COL
        IF (ASSOCIATED(root%PERM_ROW)) DEALLOCATE (root%PERM_ROW)
        ALLOCATE(root%PERM_ROW(NN),stat=allocok)
        IF(allocok.GT.0) THEN
          INFO(1) = -13
          INFO(2) = NN
          IF ( LP .GT. 0 )
     *      WRITE(LP,*)
     *      ' ZMUMPS_295: Problem allocating PERM_ROW'
          RETURN
        ENDIF
        CALL ZMUMPS_174(NN,0,root%PERM_ROW)
        PERM_ROW => root%PERM_ROW
        IF (ASSOCIATED(root%ELROW)) DEALLOCATE(root%ELROW)
        ALLOCATE(root%ELROW(root%ELG_MAX),stat=allocok)
        IF(allocok.GT.0) THEN
          INFO(1) = -13
          INFO(2) = root%ELG_MAX
          IF ( LP .GT. 0 )
     *      WRITE(LP,*)
     *      ' ZMUMPS_295: Problem allocating ELROW'
          RETURN
        ENDIF
        ELROW => root%ELROW
        IF (ASSOCIATED(root%ELELG)) DEALLOCATE(root%ELELG)
        ALLOCATE(root%ELELG(root%ELG_MAX),stat=allocok)
        IF(allocok.GT.0) THEN
          INFO(1) = -13
          INFO(2) = root%ELG_MAX
          IF ( LP .GT. 0 )
     *      WRITE(LP,*)
     *      ' ZMUMPS_295: Problem allocating ELELG'
          RETURN
        ENDIF
        ELELG => root%ELELG
        IF (ASSOCIATED(root%EUROW)) DEALLOCATE(root%EUROW)
        ALLOCATE(root%EUROW(root%ELG_MAX),stat=allocok)
        IF(allocok.GT.0) THEN
          INFO(1) = -13
          INFO(2) = root%ELG_MAX
          IF ( LP .GT. 0 )
     *      WRITE(LP,*)
     *      ' ZMUMPS_295: Problem allocating EUROW'
          RETURN
        ENDIF
        EUROW => root%EUROW
        IF (ASSOCIATED(root%EUELG)) DEALLOCATE(root%EUELG)
        ALLOCATE(root%EUELG(root%ELG_MAX),stat=allocok)
        IF(allocok.GT.0) THEN
          INFO(1) = -13
          INFO(2) = root%ELG_MAX
          IF ( LP .GT. 0 )
     *      WRITE(LP,*)
     *      ' ZMUMPS_295: Problem allocating EUELG'
          RETURN
        ENDIF
        EUELG => root%EUELG
        IF (ASSOCIATED(root%DL)) DEALLOCATE(root%DL)
        ALLOCATE(root%DL(root%ELG_MAX),stat=allocok)
        IF(allocok.GT.0) THEN
          INFO(1) = -13
          INFO(2) = root%ELG_MAX
          IF ( LP .GT. 0 )
     *      WRITE(LP,*)
     *      ' ZMUMPS_295: Problem allocating DL'
          RETURN
        ENDIF
        DL => root%DL
        IF (ASSOCIATED(root%PTREL)) DEALLOCATE(root%PTREL)
        ALLOCATE(root%PTREL(root%NULL_MAX),stat=allocok)
        IF(allocok.GT.0) THEN
          INFO(1) = -13
          INFO(2) = root%NULL_MAX
          IF ( LP .GT. 0 )
     *      WRITE(LP,*)
     *      ' ZMUMPS_295: Problem allocating PTREL'
          RETURN
        ENDIF
        PTREL => root%PTREL
        IF (ASSOCIATED(root%PTREU)) DEALLOCATE(root%PTREU)
        ALLOCATE(root%PTREU(root%ELG_MAX),stat=allocok)
        IF(allocok.GT.0) THEN
          INFO(1) = -13
          INFO(2) = root%ELG_MAX
          IF ( LP .GT. 0 )
     *      WRITE(LP,*)
     *      ' ZMUMPS_295: Problem allocating PTREU'
          RETURN
        ENDIF
        PTREU => root%PTREU
        ELIND=0
        EUIND=0
        PTREL(1)=1
        PTREU(1)=1
        NLUPDATE=0
        NUUPDATE=0
        CALL ZMUMPS_168(NN,A,LDA,DEFICIENCY,5,WR03(NN+1),
     *    PERM_COL,ELROW,ELELG,DL,
     *    EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WR03,SMIN)
        CALL ZMUMPS_182(NN,A,LDA,DEFICIENCY,WR03(NN+1),RIGHT_SV,
     *    PERM_COL,ELROW,ELELG,DL,
     *    EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WR03,SMIN)
        XMAX=ZMUMPS_106(I,WR03(NN+1))
        YMAX=ZMUMPS_106(I,RIGHT_SV)
        DO L=I,1,-1
          ABSYL=ABS(RIGHT_SV(L))
          IF(YMAX*SV_MU.LT.ABSYL) GO TO 60
        END DO
 60     CONTINUE
        L=I
        ROT_STEPS=SMP_DEFICIENCY
        IND_ROT=0
        DO II=1,ROT_STEPS
          DO K=I,1,-1
            ABSXK=ABS(WR03(NN+K))
            IF(XMAX*SV_MU.LT.ABSXK) GO TO 40
          END DO
 40       CONTINUE
          IND_ROT=IND_ROT+1
        END DO
        ROT_STEPS=0
        KK=L
        DO II=1,ROT_STEPS
        END DO
        LUSHIFT=KK-1
        DIMLU=I-LUSHIFT
        allocate(AL(DIMLU,DIMLU),stat=IERR)
        CALL ZMUMPS_117(DIMLU*DIMLU,ZERO,AL,1)
        IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *    ' Allocation error',8,2)
        allocate(AU(DIMLU,DIMLU),stat=IERR)
        CALL ZMUMPS_117(DIMLU*DIMLU,ZERO,AU,1)
        IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *    ' Allocation error',8,2)
        IF(KK.LT.I) THEN
          CALL ZMUMPS_117(DIMLU*DIMLU,ZERO,AU,1)
          DO J=KK,I
            DO II=KK,J
              TEMP=A(LDA*(J-1)+II)
              AU(II-LUSHIFT,J-LUSHIFT)=TEMP
            END DO
          END DO
          NUUPDATE=1
          IF(.NOT.ROT_STEPS.EQ.0) THEN
            DO JJ=1,ROT_STEPS
              K=I
              LSHIFT=K-1
              DIML=I-LSHIFT
              NLUPDATE=NLUPDATE+1
              CALL ZCOPY(I-DIML,A(K),LDA,WR03,1)
              DO II=K+1,I
                CALL ZCOPY(I-DIML,A(II),LDA,
     *            A(II-1),LDA)
              END DO
              CALL ZCOPY(I-DIML,WR03,1,A(NN),LDA)
              DO J=K,I
                AL(J-LSHIFT,J-LSHIFT)= ONE
                DO II=J+1,I
                  TEMP=A(LDA*(J-1)+II)
                  AL(II-LSHIFT,J-LSHIFT)=TEMP
                END DO
              END DO
              CALL ZCOPY(DIML,AL(K-LSHIFT,1),DIMLU,WR03,1)
              DO II=K-LSHIFT+1,I-LSHIFT
                CALL ZCOPY(DIML,AL(II,1),DIMLU,AL(II-1,1),DIMLU)
              END DO
              CALL ZCOPY(DIML,WR03,1,AL(I-LSHIFT,1),DIMLU)
              CALL ZMUMPS_171(K,I,PERM_ROW)
              IF(ELG_MAX.LT.ELIND+I-K) CALL ZMUMPS_298(0,2,
     *          'SPARSLAB','SLABRR3',
     *          ' Missing space to store elementary matrices',1604,2)
              ELIND_SAVE=ELIND+1
              DO J=K-LSHIFT,I-LSHIFT-1
                CALL ZMUMPS_177(J,AL,DIMLU,LSHIFT,ELIND,ELROW,ELELG)
              END DO
              ELIND=ELIND+1
              ELROW(ELIND)=elrow(ptrel(nlupdate))
              PTREL(NLUPDATE+1)=ELIND+1
              DO J=PTREL(NLUPDATE),ELIND-1
                II=ELROW(J)
                ELROW(J)=0
                ELG=ELELG(J)
                IF(II.GT.0) THEN
                  CALL ZMUMPS_29(DIMLU,-ELG,AU(II-LUSHIFT,1),
     *              DIMLU,AU(II-LUSHIFT+1,1),DIMLU)
                ELSE
                  CALL ZMUMPS_33(DIMLU,ELG,AU(-II-LUSHIFT,1),
     *              DIMLU,AU(-II-LUSHIFT+1,1),DIMLU)
                END IF
              END DO
              DO II=1,DIML
                DL(ELIND_SAVE)= ONE/AL(II,II)
                ELIND_SAVE=ELIND_SAVE+1
                CALL ZSCAL(DIML-II+1,ONE/AL(II,II),AL(II,II),1)
              END DO
              IND=DIMLU
              DO II=DIML,1,-1
                ELIND_SAVE=ELIND_SAVE-1
                CALL ZSCAL(DIMLU,ONE/DL(ELIND_SAVE),AU(IND,1),DIMLU)
                DL(ELIND_SAVE)= ONE
                IND=IND-1
              END DO
            DO J=K,I
              DO II=J+1,I
                TEMP=AL(II-LSHIFT,J-LSHIFT)
                A(LDA*(J-1)+II)=TEMP
              END DO
            END DO
            NLUPDATE=0
            ELIND=0
            END DO
            DO J=DIMLU-DIMLU+1,DIMLU
              DO II=DIMLU-1,J,-1 ! DANO NAZPET Z DIMLU-2
                TEMP1=AU(II,J)
                TEMP2=AU(II+1,J)
                IF(ABS(TEMP1).EQ.ZERO.AND.ABS(TEMP2).EQ.ZERO) THEN
                ELSEIF(ABS(TEMP1).GT.ABS(TEMP2)) THEN
                  CALL ZMUMPS_157(TEMP1,TEMP2,ELG)
                  CALL ZMUMPS_29(DIMLU-J+1,ELG,AU(II,J),
     *              DIMLU,AU(II+1,J),DIMLU)
                  EUIND=EUIND+1
                  EUROW(EUIND)=II+LUSHIFT
                  EUELG(EUIND)=ELG
                ELSE
                  CALL ZMUMPS_157(AU(II+1,J),AU(II,J),ELG)
                  CALL ZMUMPS_30(DIMLU-J+1,ELG,AU(II,J),
     *              DIMLU,AU(II+1,J),DIMLU)
                  EUIND=EUIND+1
                  EUROW(EUIND)=-(II+LUSHIFT)
                  EUELG(EUIND)=ELG
                END IF
              END DO
            END DO
            deallocate(AL)
          END IF
        END IF ! End if QR options:
        IF(L.LT.I) THEN
          USHIFT=L-1
          DIMU=I-USHIFT
          KND=(L-1)*NN+1
          CALL ZCOPY(I-DIMLU,A((L-1)*LDA+1),
     *      1,WR03,1)
          DO II=L+1,I
            CALL ZCOPY(I-DIMLU,A(KND+NN),
     *        1,A(KND),1)
            KND=KND+NN
          END DO
          CALL ZCOPY(I-DIMLU,WR03,1,A((I-1)*LDA+1),1)
          CALL ZCOPY(DIMLU,AU(1,DIMLU-DIMU+1),1,WR03,1)
          DO II=L-USHIFT+DIMLU-DIMU+1,I-USHIFT+DIMLU-DIMU
            CALL ZCOPY(DIMLU,AU(1,II),1,AU(1,II-1),1)
          END DO
          CALL ZCOPY(DIMLU,WR03,1,AU(1,I-USHIFT+DIMLU-DIMU),1)
          CALL ZMUMPS_171(L,I,PERM_COL)
          IND=(NULL_MAX-NN+I-1)*NN+1
          DO II=NULL_MAX-NN+I,NULL_MAX
            IND1=L
            IND2=I
            CALL ZMUMPS_105(IND1,IND2,A(IND))
            IND=IND+NN
          END DO
          IF(ELG_MAX.LT.EUIND+I-L)
     *      CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *      ' Missing space to store elementary matrices',1604,2)
          DO J=L-USHIFT+DIMLU-DIMU,I-USHIFT-1+DIMLU-DIMU ! DANO ZPET
            CALL ZMUMPS_331(J,AU,DIMLU,LUSHIFT,EUIND,EUROW,EUELG)
          END DO
        END IF
        PTREU(NUUPDATE+1)=EUIND+1
        IF(KK.LT.I) THEN
          DO J=KK,I
            DO II=KK,J
              TEMP=AU(II-LUSHIFT,J-LUSHIFT)
              A(LDA*(J-1)+II)=TEMP
            END DO
          END DO
        END IF
        deallocate(AU)
      CALL ZMUMPS_117(DEFICIENCY,ZERO,WR03(2*NN-DEFICIENCY+1),1)
      CALL ZMUMPS_117(DEFICIENCY,ZERO,RIGHT_SV(NN-DEFICIENCY+1),1)
      CALL ZMUMPS_168(NN,A,LDA,DEFICIENCY,5,WR03(NN+1),
     *  PERM_COL,ELROW,ELELG,DL,
     *  EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WR03,SMIN)
       CALL ZMUMPS_182(NN,A,LDA,DEFICIENCY,WR03(NN+1),RIGHT_SV,
     *  PERM_COL,ELROW,ELELG,DL,
     *  EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WR03,SMIN)
      IF(ABS(SMIN).LE.TOL.AND.DEFICIENCY.LE.MAXDEF) THEN
        IF(NN-I.GT.NULL_MAX)
     *    CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *    ' Missing space to store null space basis',1602,2)
      ELSE
        GO TO 210
      END IF
      CALL ZMUMPS_334(NN,PERM_COL,RIGHT_SV,WR03(NN+1))
      ALPHA=RIGHT_SV(I)
      ALPHA=-ONE/ALPHA
      CALL ZSCAL(I,ALPHA,RIGHT_SV,1)
 110  CONTINUE
        NULL_IND=1+(I-1)*LDA
        CALL ZMUMPS_117(I-1,ZERO,A(NULL_IND),1)
        CALL ZCOPY(I-1,RIGHT_SV,1,A(NULL_IND),1)
        NULL_IND=NULL_IND-LDA
        DEFICIENCY=DEFICIENCY+1
        I=I-1
        YMAX=ZMUMPS_106(I,RIGHT_SV)
        DO L=I,1,-1
          ABSYL=ABS(RIGHT_SV(L))
          IF(YMAX*SV_MU.LT.ABSYL) GO TO 65
        END DO
 65     CONTINUE
        L=I
        KK=L
        LUSHIFT=KK-1
        DIMLU=I-LUSHIFT
        allocate(AU(DIMLU,DIMLU),stat=IERR)
        CALL ZMUMPS_117(DIMLU*DIMLU,ZERO,AU,1)
        IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *    ' Allocation error',8,2)
        IF(L.LT.I) THEN
          CALL ZMUMPS_117(DIMLU*DIMLU,ZERO,AU,1)
          DO J=KK,I
            DO II=KK,J
              TEMP=A(LDA*(J-1)+II)
              AU(II-LUSHIFT,J-LUSHIFT)=TEMP
            END DO
          END DO
          USHIFT=L-1
          DIMU=I-USHIFT
          KND=(L-1)*NN+1
          CALL ZCOPY(I-DIMLU,A((L-1)*LDA+1),1,WR03,1)
          DO II=L+1,I
            CALL ZCOPY(I-DIMLU,
     *        A(KND+LDA),1,A(KND),1)
            KND=KND+LDA
          END DO
          CALL ZCOPY(I-DIMLU,WR03,1,A((I-1)*LDA+1),1)
          CALL ZCOPY(DIMLU,AU(1,DIMLU-DIMU+1),1,WR03,1)
          DO II=L-USHIFT+DIMLU-DIMU+1,I-USHIFT+DIMLU-DIMU
            CALL ZCOPY(DIMLU,AU(1,II),1,AU(1,II-1),1)
          END DO
          CALL ZCOPY(DIMLU,WR03,1,AU(1,I-USHIFT+DIMLU-DIMU),1)
          CALL ZMUMPS_171(L,I,PERM_COL)
          IND=(NULL_MAX-NN+I-1)*NN+1
          DO II=NULL_MAX-NN+I,NULL_MAX
            IND1=L
            IND2=I
            CALL ZMUMPS_105(IND1,IND2,A(IND))
            IND=IND+NN
          END DO
          IF(ELG_MAX.LT.EUIND+I-L) CALL ZMUMPS_298
     *      (0,2,'SPARSLAB','SLABRR3',
     *      ' Missing space to store elementary matrices',1604,2)
          DO J=L-USHIFT+DIMLU-DIMU,I-USHIFT-1+DIMLU-DIMU ! DANO ZPET
            CALL ZMUMPS_331(J,AU,DIMLU,LUSHIFT,EUIND,EUROW,EUELG)
          END DO
          PTREU(NUUPDATE+1)=EUIND+1
          DO J=KK,I
            DO II=KK,J
              TEMP=AU(II-LUSHIFT,J-LUSHIFT)
              A(LDA*(J-1)+II)=TEMP
            END DO
          END DO
        END IF
        deallocate(AU)
        CALL ZMUMPS_117(DEFICIENCY,ZERO,WR03(2*NN-DEFICIENCY+1),1)
        CALL ZMUMPS_117(DEFICIENCY,ZERO,RIGHT_SV(NN-DEFICIENCY+1),1)
        CALL ZMUMPS_317(I,NN,A,WR03(NN+1),WR03,
     *              WR03(2*NN+1),SMIN,RIGHT_SV)
        CALL ZCOPY(NN,RIGHT_SV,1,WR03(NN+1),1)
        CALL ZMUMPS_333(NN,PERM_COL,WR03(NN+1),WR03)
        IF(ABS(SMIN).LT.TOL.AND.NN-I.LE.MAXDEF) THEN
          IF(NN-I.GT.NULL_MAX)
     *      CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *      ' Missing space to store null space basis',1602,2)
          GO TO 110
        END IF
 210    CONTINUE
        NUM_RANK=I
        root%ELIND = ELIND
        root%EUIND = EUIND
        root%NLUPDATE = NLUPDATE
        root%NUUPDATE = NUUPDATE
#endif
      END IF
      KEEP(17)=DEFICIENCY
      RETURN
      END SUBROUTINE ZMUMPS_295
      SUBROUTINE ZMUMPS_296
     *  (NN,A,ROOT,RHS,WK,LWK,KEEP,KEEP8,MTYPE,INFO)
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      INTEGER NN,LWK
      TYPE ( ZMUMPS_ROOT_STRUC ) :: root
      COMPLEX*16 A(*)
      INTEGER INFO(2),KEEP(*)
      COMPLEX*16 WK(*),RHS(*)
      INTEGER I,J,K,LDA,GIND,ALLOCOK,IERR,MTYPE_SAVE
      INTEGER LDLT,QR,MTYPE,DEFICIENCY
      INTEGER ELIND,EUIND,NLUPDATE,NUUPDATE
      COMPLEX*16 C,S
      DOUBLE PRECISION TOL16,TOL13,TOL11
      PARAMETER(TOL16=1.0D-16,TOL13=1.0D-13,TOL11=1.1D-11)
      COMPLEX*16, DIMENSION(:), POINTER :: GROW, GCOS, GSIN
      INTEGER, DIMENSION(:), POINTER :: PERM_ROW
      INTEGER, DIMENSION(:), POINTER :: PERM_COL
      COMPLEX*16, DIMENSION(:), POINTER :: ELELG
      COMPLEX*16, DIMENSION(:), POINTER :: EUELG
      INTEGER, DIMENSION(:), POINTER :: ELROW
      INTEGER, DIMENSION(:), POINTER :: EUROW
      INTEGER, DIMENSION(:), POINTER :: PTREL
      INTEGER, DIMENSION(:), POINTER :: PTREU
      COMPLEX*16, DIMENSION(:), POINTER :: DL
      LDLT=KEEP(50)
      QR=KEEP(19)
      IF ( QR.eq.9 .OR. QR.eq.10) THEN
        QR = QR - 8
      END IF
      DEFICIENCY=KEEP(17)
      GIND=root%GIND
      MTYPE_SAVE=MTYPE
      IF(MTYPE.EQ.1) THEN
        MTYPE=0
      ELSEIF(MTYPE.EQ.0) THEN
        MTYPE=1
      END IF
      LDA=NN
      IERR = 0
      IF(QR.EQ.1.OR.QR.EQ.2.OR.QR.EQ.3.OR.QR.EQ.4) THEN
        IF(MTYPE.EQ.1) THEN
          CALL ZUNMQR('L','T',NN,1,NN,A,LDA,root%QR_TAU,
     *      RHS,LDA,WK(1),LWK,IERR)
          IF(IERR.LT.0) THEN
            WRITE(*,*) 'Error return from ZUNMQR: IERR=', IERR
            RETURN
          END IF
#if defined(try_null_space_not_tested)
          DO J=1,GIND
            K=root%GROW(J)
            C=root%GCOS(J)
            S=root%GSIN(J)
            CALL ZMUMPS_31(1,C,S,RHS(K),LDA,RHS(K+1),LDA)
          END DO
#endif
          CALL ZTRTRS('U','N','N',NN-DEFICIENCY,1,
     *      A,LDA,RHS,LDA,IERR)
          RHS( NN - DEFICIENCY + 1: NN ) = 0.0D0
        ELSEIF(MTYPE.EQ.0) THEN
          CALL ZMUMPS_334(NN,root%IPIV,RHS,WK)
          CALL ZTRTRS('U','T','N',NN-DEFICIENCY,1,
     *      A,LDA,RHS,LDA,IERR)
          RHS( NN - DEFICIENCY + 1: NN ) = 0.0D0
          IF(IERR.LT.0) THEN
              WRITE(*,*) 'Error return from ZTRTRS: IERR=', IERR
              RETURN
          END IF
#if defined(try_null_space_not_tested)
          DO J=1,GIND
            K=root%GROW(J)
            C=root%GCOS(J)
            S=root%GSIN(J)
            CALL ZMUMPS_31(1,C,S,RHS(K),1,RHS(K+1),1)
          END DO
#endif
          CALL ZUNMQR( 'L','N',NN,1,NN,A,LDA,root%QR_TAU,
     *      RHS,LDA,WK,LWK,IERR)
          IF(IERR.LT.0) THEN
              WRITE(*,*) 'Error return from ZUNMQR: IERR=', IERR
              RETURN
          END IF
        ENDIF
#if defined(try_null_space_not_tested)
      ELSEIF(QR.EQ.5.OR.QR.EQ.6) THEN
        IF(MTYPE.EQ.1) THEN
          CALL ZGETRS('N',NN-DEFICIENCY,1,A,LDA,root%IPIV,
     *      RHS,LDA,IERR)
        ELSEIF(MTYPE.EQ.0) THEN
          CALL ZGETRS('T',NN-DEFICIENCY,1,A,LDA,root%IPIV,
     *      RHS,LDA,IERR)
        END IF
      ELSEIF(QR.EQ.7.OR.QR.EQ.8) THEN
        PERM_COL => root%PERM_COL
        PERM_ROW => root%PERM_ROW
        ELROW => root%ELROW
        ELELG => root%ELELG
        EUROW => root%EUROW
        EUELG => root%EUELG
        DL => root%DL
        PTREL => root%PTREL
        PTREU => root%PTREU
        ELIND = root%ELIND
        EUIND = root%EUIND
        NLUPDATE = root%NLUPDATE
        NUUPDATE =root%NUUPDATE
        IF ( MTYPE .eq. 1 ) THEN
          CALL ZMUMPS_304(NN-DEFICIENCY,A,LDA,DL,NLUPDATE,
     *      PTREL,ELROW,ELELG,
     *      NUUPDATE,PTREU,EUROW,EUELG,
     *      root%IPIV,PERM_ROW,PERM_COL,
     *      RHS,WK)
        ELSE
          CALL ZMUMPS_306(NN-DEFICIENCY,A,LDA,DL,NLUPDATE,
     *      PTREL,ELROW,ELELG,
     *      NUUPDATE,PTREU,EUROW,EUELG,
     *      root%IPIV,PERM_ROW,PERM_COL,
     *      RHS,WK)
        END IF
#endif
      END IF
      IF ( IERR .LT. 0 ) THEN
        WRITE(*,*) ' Problem during solve of the root: IERR =',IERR
        RETURN
      END IF
      IF((QR.EQ.1.OR.QR.EQ.2.OR.QR.EQ.3.OR.QR.EQ.4).and.
     *  mtype.eq.1) THEN
        CALL ZMUMPS_333(NN,root%IPIV,RHS,WK)
      END IF
      MTYPE=MTYPE_SAVE
      RETURN
      END SUBROUTINE ZMUMPS_296
      SUBROUTINE ZMUMPS_297(N,A)
      INTEGER N
      COMPLEX*16 A( N, N )
      INTEGER I,J
      DO I = 2, N
        DO J = 1, I - 1
          A( I, J ) = A( J, I )
        END DO
      END DO
      RETURN
      END SUBROUTINE ZMUMPS_297
#endif
      SUBROUTINE ZMUMPS_625()
      RETURN
      END SUBROUTINE ZMUMPS_625
      SUBROUTINE ZMUMPS_85(NN,A,LDA,NULL,
     *  COLNO,DEFICIENCY,QR_TAU,IPIV,WR01,LWK)
      IMPLICIT NONE
      INTEGER NN,COLNO,DEFICIENCY,LDA,LWK
      INTEGER IPIV(NN)
      COMPLEX*16 A(NN*NN),NULL(NN),WR01(LWK)
      COMPLEX*16 QR_TAU(NN)
      INTEGER I,J,K,IND,JND,IERR
      COMPLEX*16 TEMP
      DOUBLE PRECISION ZERO,TEMP1
      PARAMETER(ZERO=0.0D0)
      IND=(NN-COLNO)*LDA+1
      CALL ZCOPY(NN-DEFICIENCY,A(IND),1,NULL,1)
      DO I=NN,NN-DEFICIENCY+1,-1
        NULL(I)= DCMPLX(ZERO)
      END DO
      NULL(NN-COLNO+1)=-1.0
      TEMP1=ZERO
      JND=1
      DO I=1,NN
        TEMP=0.0D0
        IND=(I-1)*LDA+I
        DO J=JND,NN
          TEMP=TEMP+A(IND)*NULL(J)
          IND=IND+LDA
        END DO
        TEMP1=MAX(TEMP1,ABS(TEMP))
        JND=JND+1
      END DO
       CALL ZMUMPS_333(NN,IPIV,null,wr01)
      RETURN
      END
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_82(NN,A,LDA,NULL,
     *  COLNO,DEFICIENCY,IPIV)
      IMPLICIT NONE
      INTEGER NN,COLNO,DEFICIENCY,LDA,LWK
      INTEGER IPIV(*)
      COMPLEX*16 A(*),NULL(*),WR01
      INTEGER I,J,K,IND,JND,IERR
      COMPLEX*16 TEMP
      DOUBLE PRECISION ZERO,ONE,TEMP1
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      IND=(NN-COLNO)*LDA+1
      CALL ZCOPY(NN-DEFICIENCY,A(IND),1,NULL,1)
      DO I=NN,NN-DEFICIENCY+1,-1
        NULL(I)= DCMPLX(ZERO)
      END DO
      NULL(NN-COLNO+1)=-1.0D0
      TEMP1=ZERO
      JND=1
      DO I=1,NN
        TEMP=0.0D0
        IND=(I-1)*LDA+I
        DO J=JND,NN
          TEMP=TEMP+A(IND)*NULL(J)
          IND=IND+LDA
        END DO
        TEMP1=MAX(TEMP1,ABS(TEMP))
        JND=JND+1
      END DO
      CALL ZLASWP(1,NULL,LDA,1,LDA,IPIV,-1)
      RETURN
      END
#endif
      SUBROUTINE ZMUMPS_84( id )
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (ZMUMPS_STRUC), TARGET :: id
      INCLUDE 'mpif.h'
      INTEGER LIW, LA
      INTEGER INULL_SOL, INULL, IW2
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      COMPLEX*16 ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER MTYPE
      INTEGER i,j,K
      COMPLEX*16 temp
      INTEGER LIWK_SOLVE, LWK_SOLVE
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWK_SOLVE
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: WK_SOLVE
      DOUBLE PRECISION DZNRM2
      EXTERNAL DZNRM2
      INTEGER allocok
      LOGICAL I_AM_SLAVE
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     *             ( id%MYID .eq. MASTER .AND.
     *               id%KEEP(46) .eq. 1 ) )
      IF ( I_AM_SLAVE ) THEN
        LA = id%KEEP(31)
        LIW = id%KEEP(32)
        IF ( id%KEEP(50) .eq. 0 ) THEN
          MTYPE = 1
        ELSE
          MTYPE = 0
        END IF
        LIWK_SOLVE = id%KEEP(28) * 4 + 2
        ALLOCATE( IWK_SOLVE( LIWK_SOLVE ),
     *            stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = LIWK_SOLVE
        END IF
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
        IF ( associated( id%NULL_SPACE ) )
     *       DEALLOCATE( id%NULL_SPACE )
        ALLOCATE( id%NULL_SPACE( id%N,
     *            id%KEEP(17) ),
     *            stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          WRITE(*,*) ' Error : unable to allocate null space basis'
          id%INFO(1) = -13
          id%INFO(2) = id%N * id%KEEP(17)
        END IF
      END IF
      CALL ZMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      IF ( I_AM_SLAVE ) THEN
        LWK_SOLVE = id%N * 3
      ELSE
        LWK_SOLVE = id%N * 3
      END IF
      ALLOCATE( WK_SOLVE( LWK_SOLVE ), stat=allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) 'Error : unable to allocate WK_SOLVE'
        id%INFO(1) = -13
        id%INFO(2) = LWK_SOLVE
      END IF
      CALL ZMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      INULL_SOL = 1
      INULL = INULL_SOL + id%N
      IW2 = INULL + id%N
      DO i = 1, id%KEEP(17)
        IF ( I_AM_SLAVE ) THEN
        CALL ZMUMPS_247( i, id%root, id%N,
     *  id%S(1), LA, id%IS(1), LIW,
     *  id%S(LA + 1), id%MAXS - LA,
     *  id%IS(LIW + 1), id%MAXIS - LIW,
     *  WK_SOLVE(INULL), WK_SOLVE(INULL_SOL),
     *  id%NA, id%LNA, id%NE_STEPS(1),
     *  WK_SOLVE(IW2), MTYPE, id%ICNTL,
     *  id%STEP(1), id%FRERE_STEPS(1),id%FILS(1),
     *  id%PTLUST_S(1),id%PTRFAC(1), IWK_SOLVE, LIWK_SOLVE,
     *  id%PROCNODE_STEPS, id%NSLAVES,
     *  id%INFO, id%KEEP,id%KEEP8,
     *  id%COMM,
     *  id%COMM_NODES, id%MYID,
     *  id%MYID_NODES,
     *  id%BUFR, id%LBUFR, id%LBUFR_BYTES, 
     *  id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1)
     *  )
      ELSE
        WK_SOLVE(INULL_SOL:INULL_SOL+id%N-1) = ZERO
      END IF
      CALL ZMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      IF ( id%MYID .eq. MASTER ) THEN
        CALL MPI_REDUCE( WK_SOLVE( INULL_SOL ),
     *       id%NULL_SPACE( 1, i ),
     *       id%N, MPI_DOUBLE_COMPLEX,
     *       MPI_SUM, MASTER, id%COMM, IERR )
        CALL ZMUMPS_259(id%NULL_SPACE(1,i), id%N)
      ELSE
        CALL MPI_REDUCE( WK_SOLVE( INULL_SOL ),
     *       WK_SOLVE(IW2),
     *       id%N, MPI_DOUBLE_COMPLEX,
     *       MPI_SUM, MASTER, id%COMM, IERR )
      END IF
      IF (id%MYID .eq. MASTER) THEN
        IF  ( ((id%KEEP(52) .GT. 0)
     *   .AND. (id%KEEP(52) .LE. 6)) .OR.
     *        ( id%KEEP(52) .EQ. -1) .OR.
     *        ( id%KEEP(52) .EQ. -2) ) THEN
          DO K = 1, id%N
          id%NULL_SPACE(K,i) = id%NULL_SPACE(K,i) *
     *                                id%COLSCA( K )
          END DO
        END IF
        IF ( id%KEEP(23) .NE. 0 ) THEN
          WK_SOLVE( IW2: IW2 + id%N - 1 ) =
     *    id%NULL_SPACE( 1: id%N, i )
          DO K = 1, id%N
            id%NULL_SPACE( id%UNS_PERM( K ), i ) = 
     *      WK_SOLVE( IW2 + K - 1 )
          END DO
        END IF
      END IF
      END DO
      IF ( I_AM_SLAVE )
     *     DEALLOCATE( IWK_SOLVE )
      DEALLOCATE( WK_SOLVE )
      RETURN
      END SUBROUTINE ZMUMPS_84
      SUBROUTINE ZMUMPS_279( MBLOCK, NBLOCK, SIZE_ROOT,
     *           LOCAL_M, LOCAL_N, ROOT_OWNER, KEEP,KEEP8,
     *           LIWK_QR, LWK_QR )
      INTEGER MBLOCK, NBLOCK, SIZE_ROOT, LOCAL_M, LOCAL_N
      INTEGER LIWK_QR, LWK_QR,QR,PAR_ROOT
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      LOGICAL ROOT_OWNER
      QR=KEEP(19)
      PAR_ROOT=KEEP(38)
      LIWK_QR = 0
      IF(PAR_ROOT.NE.0) THEN
        IF(QR.GE.1.AND.QR.LE.4) THEN
          LWK_QR = MAX(5*MBLOCK*NBLOCK+
     *     5*MAX(LOCAL_M,LOCAL_N)+6+SIZE_ROOT,4*SIZE_ROOT+5)
        ELSE
          LWK_QR = MAX(5*MBLOCK*NBLOCK+
     *     5*MAX(LOCAL_M,LOCAL_N)+6+SIZE_ROOT,4*SIZE_ROOT+5)
        END IF
      ELSE
        IF(ROOT_OWNER) THEN
          IF(QR.GE.1.AND.QR.LE.4) THEN
            LWK_QR=10000+5*SIZE_ROOT
          ELSE
            LWK_QR=10000+5*SIZE_ROOT
          END IF
        ELSE
          LWK_QR=SIZE_ROOT+1
        END IF
      END IF
      RETURN
      END
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_103(MYID, M, N, ASEQ,
     *                    LOCAL_M, LOCAL_N,
     *                    MBLOCK, NBLOCK,
     *                    APAR,
     *                    MASTER_ROOT,
     *                    NPROW, NPCOL,
     *                    COMM)
      IMPLICIT NONE
      INTEGER MYID, MASTER_ROOT, COMM
      INTEGER M, N
      INTEGER NPROW, NPCOL
      INTEGER LOCAL_M, LOCAL_N
      INTEGER MBLOCK, NBLOCK
      COMPLEX*16 APAR( LOCAL_M, LOCAL_N )
      COMPLEX*16 ASEQ( M, N )
      INCLUDE 'mpif.h'
      INTEGER I, J, SIZE_IBLOCK, SIZE_JBLOCK, ISOUR, IROW, ICOL
      INTEGER IBLOCK, JBLOCK, II, JJ, KK
      INTEGER IAPAR, JAPAR, IERR,idist
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMPLEX*16 WK( MBLOCK * NBLOCK )
      LOGICAL JUPDATE
        IAPAR = 1
        JAPAR = 1
        DO J=1,N,NBLOCK
          SIZE_JBLOCK=NBLOCK
          IF (J+NBLOCK>N) THEN
            SIZE_JBLOCK=N-J+1
          END IF
          JUPDATE=.FALSE.
          DO I=1,M,MBLOCK
            SIZE_IBLOCK=MBLOCK
            IF (I+MBLOCK>M) THEN
              SIZE_IBLOCK=M-I+1
            END IF
            IBLOCK=I/MBLOCK
            JBLOCK=J/NBLOCK
            IROW=MOD(IBLOCK,NPROW)
            ICOL=MOD(JBLOCK,NPCOL)
            ISOUR=IROW*NPCOL+ICOL
            IF ( ISOUR .NE. MASTER_ROOT ) THEN
              IF ( MYID .EQ. MASTER_ROOT ) THEN
                CALL MPI_RECV(WK,SIZE_IBLOCK*SIZE_JBLOCK,
     *                         MPI_DOUBLE_COMPLEX,
     *                         ISOUR,128,COMM,STATUS,IERR)
                KK=1
                DO JJ=J,J+SIZE_JBLOCK-1
                DO II=I,I+SIZE_IBLOCK-1
                  ASEQ(II,JJ)=WK(KK)
                 KK=KK+1
                END DO
                END DO
              ELSE IF ( MYID .EQ. ISOUR ) THEN
                KK=1
                DO JJ=JAPAR,JAPAR+SIZE_JBLOCK-1
                DO II=IAPAR,IAPAR+SIZE_IBLOCK-1
                  WK(KK)=APAR(II,JJ)
                  KK=KK+1
                END DO
                END DO
                CALL MPI_SEND( WK( 1 ),
     *                         SIZE_IBLOCK*SIZE_JBLOCK,
     *                         MPI_DOUBLE_COMPLEX,
     *                         MASTER_ROOT,128,COMM,IERR)
                JUPDATE = .TRUE.
                IAPAR = IAPAR + SIZE_IBLOCK
              END IF
            ELSE IF ( MYID.EQ. MASTER_ROOT ) THEN
              ASEQ(I:I+SIZE_IBLOCK-1,J:J+SIZE_JBLOCK-1)
     *        = APAR( IAPAR:IAPAR+SIZE_IBLOCK-1,
     *                JAPAR:JAPAR+SIZE_JBLOCK-1)
              JUPDATE = .TRUE.
              IAPAR = IAPAR + SIZE_IBLOCK
            END IF
          END DO
          IF ( JUPDATE ) THEN
            IAPAR = 1
            JAPAR = JAPAR + SIZE_JBLOCK
          END IF
        END DO
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_171(START,STOP,IX)
      INTEGER START,STOP
      INTEGER IX(*)
      INTEGER I,ITEMP
      IF(STOP.LE.START) RETURN
      ITEMP=IX(START)
      DO I=START+1,STOP
        IX(I-1)=IX(I)
      END DO
      IX(STOP)=ITEMP
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_159(X,Y,CS,SN)
      COMPLEX*16 X,Y
      DOUBLE PRECISION CS,SN
      COMPLEX*16 TEMP
      DOUBLE PRECISION ZERO,ONE
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      IF(Y.EQ.ZERO) THEN
        CS= ONE
        SN= ZERO
      ELSE
        IF(ABS(Y).GT.ABS(X)) THEN
          TEMP=X/Y
          SN= ONE/SQRT(ONE+TEMP*TEMP)
          CS=SN*TEMP
        ELSE
          TEMP=Y/X
          CS= ONE/SQRT(ONE+TEMP*TEMP)
          SN=CS*TEMP
        END IF
      END IF
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_32(LEN,CS,SN,X,INCX,Y,INCY)
      INTEGER LEN,INCX,INCY
      DOUBLE PRECISION CS,SN
      COMPLEX*16 X(*),Y(*)
      INTEGER I,J,K
      DOUBLE PRECISION TEMP
      I=1
      J=1
      DO K=1,LEN
        TEMP=CS*X(I)-SN*Y(J)
        Y(J)=SN*X(I)+CS*Y(J)
        X(I)=TEMP
        I=I+INCX
        J=J+INCY
      END DO
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_160(MYX,XBLOCK,NPY,GLO,LOC)
      IMPLICIT NONE
      INTEGER MYX,XBLOCK,NPY,LOC,GLO
      INTEGER IBLOCK
      IBLOCK=(GLO-1)/XBLOCK
      LOC=GLO-XBLOCK*IBLOCK+XBLOCK*(IBLOCK/NPY)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_161(MYX,XBLOCK,NPY,GLO,PRO)
      IMPLICIT NONE
      INTEGER MYX,XBLOCK,NPY,GLO,PRO
      INTEGER IBLOCK
      IBLOCK=(GLO-1)/XBLOCK
      PRO=MOD(IBLOCK,NPY)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_191(MYX,XBLOCK,NPY,LOC,GLO)
      IMPLICIT NONE
      INTEGER MYX,XBLOCK,NPY,LOC,GLO
      INTEGER IBLOCK,JBLOCK
      IBLOCK=LOC/XBLOCK
      JBLOCK=MOD(LOC,XBLOCK)
      IF(JBLOCK.NE.0) THEN
        GLO=NPY*IBLOCK*XBLOCK+MYX*XBLOCK+JBLOCK
      ELSE
        GLO=NPY*(IBLOCK-1)*XBLOCK+(MYX+1)*XBLOCK
      END IF
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_113(MYID,MYROW,MYCOL,NPROW,NPCOL,MBLOCK,
     *  NBLOCK,COL,
     *  LOCAL_M,DV,DMAX,GMAX,IMAX,RMAX,SV_MU,COMM)
      IMPLICIT NONE
      INTEGER LOCAL_M,MBLOCK,NBLOCK,GMAX,RMAX
      INTEGER MYID,MYROW,MYCOL,NPROW,NPCOL,IMAX,COL
      INTEGER COMM
      COMPLEX*16 DMAX,SV_MU
      COMPLEX*16 DV(*)
      INTEGER I,J,IERR,JBLOCK,COLPROC,GROUP,COLGROUP
      INTEGER COLCOMM,SIZE,LEN,MYRANK
      COMPLEX*16 TEMP
      DOUBLE PRECISION ZERO,ABSDVI
      COMPLEX*16 VALUE
      COMPLEX*16 TUPLEIN(2,1),TUPLEOUT(2,1)
      INTEGER RANGES(3,NPROW)
      PARAMETER(ZERO=0.0D0)
      INCLUDE 'mpif.h'
      MYRANK=MYROW*NPROW+MYCOL
      VALUE=SV_MU*DMAX
      JBLOCK=(COL-1)/NBLOCK
      COLPROC=MOD(JBLOCK,NPROW)
      CALL MPI_COMM_GROUP(COMM,GROUP,IERR)
      DO I=1,NPROW
        RANGES(1,I)=(I-1)*NPCOL+1
        RANGES(2,I)=(I-1)*NPCOL+NPCOL-1
        RANGES(3,I)=1
      END DO
      CALL MPI_GROUP_RANGE_EXCL(GROUP,NPROW,RANGES,COLGROUP,IERR)
      CALL MPI_COMM_CREATE(COMM,COLGROUP,COLCOMM,IERR)
      CALL MPI_GROUP_SIZE(COLGROUP,SIZE,IERR)
      CALL MPI_REDUCE(LOCAL_M,LEN,1,MPI_INTEGER,
     *  MPI_MAX,0,COLCOMM,IERR)
      LEN=LEN+1
      CALL MPI_BCAST(LEN,1,MPI_INTEGER,0,COLCOMM,IERR)
      IF(MYCOL.EQ.COLPROC) THEN
        IMAX=0
        DO I=LOCAL_M,1,-1
          ABSDVI=ABS(DV(I))
          IF(ABSDVI.GE.ABS(VALUE)) THEN
            IMAX=I
            GO TO 50
          END IF
        END DO
 50     CONTINUE
        CALL ZMUMPS_191(MYROW,MBLOCK,NPCOL,IMAX,GMAX)
        TUPLEIN(1,1)=DCMPLX(GMAX)
        IMAX=MYRANK*LEN+IMAX
        TUPLEIN(2,1)=DCMPLX(IMAX)
        CALL MPI_REDUCE(TUPLEIN,TUPLEOUT,1,MPI_2DOUBLE_PRECISION,
     *    MPI_MAXLOC,0,COLCOMM,IERR)
        IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
          GMAX=TUPLEOUT(1,1)
          J=TUPLEOUT(2,1)
          IMAX=MOD(J,LEN)
          RMAX=J/LEN
        END IF
      END IF
      CALL MPI_COMM_FREE(COLCOMM,IERR)
      CALL MPI_GROUP_FREE(GROUP,IERR)
      CALL MPI_GROUP_FREE(COLGROUP,IERR)
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_107(MYID,MYROW,MYCOL,
     *  NPROW,NPCOL,NBLOCK,COL,
     *  LOCAL_M,DV,DMAX,IMAX,RMAX,COMM)
      IMPLICIT NONE
      INTEGER LOCAL_M,NBLOCK,RMAX
      INTEGER MYID,MYROW,MYCOL,NPROW,NPCOL,IMAX,COL
      INTEGER COMM
      COMPLEX*16 DMAX
      COMPLEX*16 DV(*)
      INTEGER I,J,IERR,JBLOCK,COLPROC,GROUP,COLGROUP
      INTEGER COLCOMM,SIZE,LEN,MYRANK
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION ZERO,ABSDVI
      COMPLEX*16 TUPLEIN(2,1),TUPLEOUT(2,1)
      INTEGER RANGES(3,NPROW)
      PARAMETER(ZERO=0.0D0)
      INCLUDE 'mpif.h'
      MYRANK=MYROW*NPROW+MYCOL
      JBLOCK=(COL-1)/NBLOCK
      COLPROC=MOD(JBLOCK,NPROW)
      CALL MPI_COMM_GROUP(COMM,GROUP,IERR)
      DO I=1,NPROW
        RANGES(1,I)=(I-1)*NPCOL+1
        RANGES(2,I)=(I-1)*NPCOL+NPCOL-1
        RANGES(3,I)=1
      END DO
      CALL MPI_GROUP_RANGE_EXCL(GROUP,NPROW,RANGES,COLGROUP,IERR)
      CALL MPI_COMM_CREATE(COMM,COLGROUP,COLCOMM,IERR)
      CALL MPI_GROUP_SIZE(COLGROUP,SIZE,IERR)
      CALL MPI_REDUCE(LOCAL_M,LEN,1,MPI_INTEGER,
     *  MPI_MAX,0,COLCOMM,IERR)
      LEN=LEN+1
      CALL MPI_BCAST(LEN,1,MPI_INTEGER,0,COLCOMM,IERR)
      IF(MYCOL.EQ.COLPROC) THEN
        TEMP=ZERO
        IMAX=1
        DO I=1,LOCAL_M
          ABSDVI=ABS(DV(I))
          IF(ABSDVI.GT.TEMP) THEN
            IMAX=I
            TEMP=ABSDVI
          END IF
        END DO
        IMAX=MYRANK*LEN+IMAX
        TUPLEIN(1,1)=TEMP
        TUPLEIN(2,1)=DCMPLX(IMAX)
        CALL MPI_REDUCE(TUPLEIN,TUPLEOUT,1,MPI_2DOUBLE_PRECISION,
     *    MPI_MAXLOC,0,COLCOMM,IERR)
        IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
          DMAX=TUPLEOUT(1,1)
          J=TUPLEOUT(2,1)
          IMAX=MOD(J,LEN)
          RMAX=J/LEN
        END IF
      END IF
      CALL MPI_COMM_FREE(COLCOMM,IERR)
      CALL MPI_GROUP_FREE(GROUP,IERR)
      CALL MPI_GROUP_FREE(COLGROUP,IERR)
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_318(MYID,NR,DEF,LOCAL_M,LOCAL_N,A,DESCA,
     *  NPROW,NPCOL,CNTXT_BLACS,COMM,SMIN,VAUX,ZERO_DIAG,
     *  INFO1,INFO2,TOL)
      IMPLICIT NONE
      INTEGER NR,DEF,LOCAL_M,LOCAL_N,INFO1,INFO2
      INTEGER NPROW,NPCOL,MYROW,MYCOL,MYID,CNTXT_BLACS
      INTEGER COMM
      INTEGER DESCA(*)
      DOUBLE PRECISION SMIN
      COMPLEX*16 A(*),VAUX(NR,*)
      INTEGER I,ROWPROC
      INTEGER LPIV, IERR, allocok
      INTEGER JJ,ROW,IDIST
      integer rsrc,csrc,lrindx,lcindx
      COMPLEX*16 qr_diag_tol
      integer DESCP(9),iarow,iacol,mblock,nblock,csrca,rsrca
      integer indxg2p,DIM,J,K,ZERO_DIAG
      INTEGER POS,ISOUR,RSRCP,CSRCP,IVROW,IVCOL
      INTEGER RSRCU,CSRCU,DESCU(9)
      INTEGER RSRCV,CSRCV,DESCV(9),IPCOL,IPROW,COLPROC
      COMPLEX*16 RROOT,MU,ETA,BETA,XJ
      DOUBLE PRECISION C,S
      COMPLEX*16 LAMBDMAX,VNEW,GAMMA
      DOUBLE PRECISION PNOR, VNOR, ZMUMPS_277, UTP, PHI
      DOUBLE PRECISION ZERO,NEARLY_ZERO,TOL,EPS,ONE
      COMPLEX*16 ALPHA
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,EPS=1.1D-15)
      COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: V_PAR
      COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: U_PAR
      COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: P_PAR
      INCLUDE 'mpif.h'
      SMIN=0.0D0
      NEARLY_ZERO=TOL
      ALLOCATE( V_PAR(1,LOCAL_N+5), stat=allocok )
      IF ( allocok.GT.0 ) THEN
        INFO1 = -13
        INFO2 = NR
        STOP
      ENDIF
      ALLOCATE( U_PAR(1,LOCAL_N+5), stat=allocok )
      IF ( allocok.GT.0 ) THEN
        INFO1 = -13
        INFO2 = NR
        STOP
      ENDIF
      ALLOCATE( P_PAR(1,LOCAL_N+5), stat=allocok )
      IF ( allocok.GT.0 ) THEN
        INFO1 = -13
        INFO2 = NR
        STOP
      ENDIF
           mblock=DESCA(5)
           nblock=DESCA(6)
           rsrca=DESCA(7)
           csrca=DESCA(8)
        CALL DESCINIT( DESCP,1,NR+5,
     *       MBLOCK,NBLOCK,0,0,
     *      CNTXT_BLACS,1,IERR)
        IF ( IERR .ne. 0 ) THEN
          WRITE(*,*) 'After DESCINIT, IERR = ', IERR
          STOP
        END IF
        RSRCP=0
        CSRCP=0
        CALL DESCINIT( DESCV,NR,1,
     *       MBLOCK,NBLOCK,0,0,
     *      CNTXT_BLACS,NR,IERR)
        IF ( IERR .ne. 0 ) THEN
          WRITE(*,*) 'After DESCINIT, IERR = ', IERR
          STOP
        END IF
        RSRCV=0
        CSRCV=0
        CALL DESCINIT( DESCU,1,NR+5,
     *       MBLOCK,NBLOCK,0,0,
     *      CNTXT_BLACS,1,IERR)
        IF ( IERR .ne. 0 ) THEN
          WRITE(*,*) 'After DESCINIT, IERR = ', IERR
          STOP
        END IF
        RSRCU=0
        CSRCU=0
      CALL BLACS_GRIDINFO(CNTXT_BLACS,NPROW,NPCOL,MYROW,MYCOL)
           call infog2l(NR-DEF,
     *         NR-DEF,
     *         DESCA,nprow,npcol,
     *        myrow,mycol,lrindx,lcindx,rsrc,csrc)
           DIM=NR-DEF+1
           CALL ZMUMPS_117(LOCAL_N+5,ZERO,V_PAR,1)
           CALL ZMUMPS_117(LOCAL_N+5,ZERO,U_PAR,1)
           CALL ZMUMPS_117(LOCAL_N+5,ZERO,P_PAR,1)
           ZERO_DIAG=0
           IAROW=INDXG2P(DIM,MBLOCK,MYROW,0,NPROW)
           IACOL=INDXG2P(DIM,NBLOCK,MYCOL,0,NPCOL)
           IF(MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
             CALL INFOG2L(DIM,DIM,DESCA,NPROW,NPCOL,MYROW,MYCOL,
     *          LRINDX,LCINDX,RSRC,CSRC)
             GAMMA=A((LCINDX-1)*LOCAL_M+LRINDX)
           END IF
           ISOUR=IAROW*NPCOL+IACOL
           CALL MPI_BCAST(GAMMA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *         COMM,IERR)
           IF(ABS(GAMMA).LE.NEARLY_ZERO) THEN
             ZERO_DIAG=1
             I=DIM
             DO J=DIM,1,-1
               IAROW=INDXG2P(J,MBLOCK,MYROW,0,NPROW)
               IACOL=INDXG2P(J,NBLOCK,MYCOL,0,NPCOL)
               IF(MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
                 CALL INFOG2L(J,J,DESCA,NPROW,NPCOL,MYROW,MYCOL,
     *              LRINDX,LCINDX,RSRC,CSRC)
                 GAMMA=A((LCINDX-1)*LOCAL_M+LRINDX)
               END IF
               ISOUR=IAROW*NPCOL+IACOL
               CALL MPI_BCAST(GAMMA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *           COMM,IERR)
               IF(ABS(GAMMA).LE.NEARLY_ZERO) I=J
             END DO
          CALL ZMUMPS_110(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *      A,0,VAUX,DIM,DIM,COLPROC,MYROW,
     *      MYCOL,NPROW,NPCOL,
     *      MYID,COMM)
           IF(MYCOL.EQ.0) THEN
           END IF
           CALL PZTRTRS('U','N','N',
     *       I-1,1,
     *       A,1,1,DESCA,VAUX,1,
     *       1,DESCV,IERR)
           DIM=I
           GO TO 400
           END IF
           CALL ZMUMPS_111(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *       A,0,P_PAR,1,DIM,ROWPROC,MYROW,MYCOL,NPROW,NPCOL,
     *       MYID,COMM)
           iarow=indxg2p(1,mblock,myrow,rsrca,nprow)
           iacol=indxg2p(1,nblock,mycol,csrca,npcol)
           IF(MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
             CALL INFOG2L(1,1,DESCA,NPROW,NPCOL,MYROW,MYCOL,
     *          LRINDX,LCINDX,RSRC,CSRC)
             ALPHA= ONE/A((LCINDX-1)*LOCAL_M+LRINDX)
             CALL ZGEBS2D(CNTXT_BLACS,'R',' ',1,1,
     *         A((LCINDX-1)*LOCAL_M+LRINDX),1)
           ELSE
             CALL ZGEBR2D(CNTXT_BLACS,'R',' ',1,1,
     *         ALPHA,1,MYROW,IACOL)
             ALPHA= ONE/ALPHA
           END IF
           IVROW=INDXG2P(1,MBLOCK,MYROW,0,NPROW)
           IVCOL=INDXG2P(1,NBLOCK,MYCOL,0,NPCOL)
           IF(MYROW.EQ.IVROW.AND.MYCOL.EQ.IVCOL) THEN
             CALL INFOG2L(1,1,DESCP,NPROW,NPCOL,MYROW,MYCOL,
     *          LRINDX,LCINDX,RSRCP,CSRCP)
             V_PAR(LRINDX,LCINDX)=ALPHA
             VNOR=ABS(ALPHA)
           END IF
           ISOUR=IVROW*NPCOL+IVCOL
           CALL MPI_BCAST(VNOR,1,MPI_DOUBLE_PRECISION,ISOUR,
     *         COMM,IERR)
           CALL PZSCAL(DIM,ALPHA,P_PAR,1,1,DESCP,1)
           DO J=2,DIM-1
           call mpi_barrier(comm,ierr)
           CALL ZMUMPS_111(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *       A,0,U_PAR,J,DIM,ROWPROC,MYROW,MYCOL,NPROW,NPCOL,
     *       MYID,COMM)
             IF(MYROW.EQ.0) THEN
             CALL PDDOT(DIM-J+1,UTP,P_PAR,1,J+1,DESCP,1,
     *         U_PAR,1,J+1,DESCU,1)
             END IF
           ISOUR=0
           CALL MPI_BCAST(UTP,1,MPI_DOUBLE_PRECISION,ISOUR,
     *         COMM,IERR)
             IAROW=INDXG2P(J,MBLOCK,MYROW,0,NPROW)
             IACOL=INDXG2P(J,NBLOCK,MYCOL,0,NPCOL)
           IF(MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
             CALL INFOG2L(J,J,DESCA,NPROW,NPCOL,MYROW,MYCOL,
     *          LRINDX,LCINDX,RSRC,CSRC)
             GAMMA=A((LCINDX-1)*LOCAL_M+LRINDX)
           END IF
           ISOUR=IAROW*NPCOL+IACOL
           CALL MPI_BCAST(GAMMA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *         COMM,IERR)
           IPROW=INDXG2P(1,MBLOCK,MYROW,0,NPROW)
           IPCOL=INDXG2P(J,NBLOCK,MYCOL,0,NPCOL)
           IF(MYROW.EQ.IPROW.AND.MYCOL.EQ.IPCOL) THEN
             CALL INFOG2L(1,J,DESCP,NPROW,NPCOL,MYROW,MYCOL,
     *          LRINDX,LCINDX,RSRCP,CSRCP)
             XJ=P_PAR(LRINDX,LCINDX)
           END IF
           ISOUR=IPROW*NPCOL+IPCOL
           CALL MPI_BCAST(XJ,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *         COMM,IERR)
             IF(MYROW.EQ.0) THEN
             CALL PDZNRM2(DIM-J+1,PHI,U_PAR,1,J+1,DESCU,1)
             END IF
           ISOUR=0
           CALL MPI_BCAST(PHI,1,MPI_DOUBLE_PRECISION,ISOUR,
     *         COMM,IERR)
             PHI= ONE+PHI*PHI
             IF(MYROW.EQ.0) THEN
             CALL PDZNRM2(DIM-J+1,PNOR,P_PAR,1,J+1,DESCP,1)
             END IF
           ISOUR=0
           CALL MPI_BCAST(PNOR,1,MPI_DOUBLE_PRECISION,ISOUR,
     *         COMM,IERR)
             ALPHA=XJ*PHI-GAMMA*UTP
             IF(ALPHA.EQ.ZERO) THEN
               BETA=GAMMA*GAMMA*(VNOR*VNOR+PNOR*PNOR)-XJ*XJ*PHI
               IF(ABS(BETA).GT.ABS(PHI)) THEN
                 S= ONE
                 C= ZERO
               ELSE
                 S= ZERO
                 C= ONE
               END IF
             ELSE
               BETA=GAMMA*GAMMA*(VNOR*VNOR+PNOR*PNOR)+(XJ*XJ-1)*PHI
               BETA=BETA-2*XJ*GAMMA*UTP
               ETA=BETA/(2.0D0*ALPHA)
               IF(ABS(ALPHA).GT.0) THEN
                 MU=ETA+ZMUMPS_277(DCMPLX(ETA),DCMPLX(ONE))
               ELSE
                 MU=ETA-ZMUMPS_277(DCMPLX(ETA),DCMPLX(ONE))
               END IF
               RROOT=ZMUMPS_277(DCMPLX(MU),DCMPLX(ONE))
               S=MU/RROOT
               C=-ONE/RROOT
             END IF
             IF(MYROW.EQ.0) THEN
             CALL PZSCAL(J-1,S,V_PAR,1,1,DESCP,1)
             END IF
             VNEW=(C-S*XJ)/GAMMA
             IVROW=INDXG2P(1,MBLOCK,MYROW,0,NPROW)
             IVCOL=INDXG2P(J,NBLOCK,MYCOL,0,NPCOL)
             IF(MYROW.EQ.IVROW.AND.MYCOL.EQ.IVCOL) THEN
               CALL INFOG2L(1,J,DESCP,NPROW,NPCOL,MYROW,MYCOL,
     *           LRINDX,LCINDX,RSRCP,CSRCP)
                 V_PAR(LRINDX,LCINDX)=VNEW
             END IF
             VNOR=ZMUMPS_277(DCMPLX(S*VNOR),VNEW)
             IF(MYROW.EQ.0) THEN
             CALL PZSCAL(DIM-J,S,P_PAR,1,J+1,DESCP,1)
             CALL PZAXPY(DIM-J,VNEW,U_PAR,1,J+1,DESCU,1,P_PAR,1,J+1,
     *         DESCP,1)
             END IF
           END DO
             IAROW=INDXG2P(DIM,MBLOCK,MYROW,0,NPROW)
             IACOL=INDXG2P(DIM,NBLOCK,MYCOL,0,NPCOL)
           IF(MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
             CALL INFOG2L(DIM,DIM,DESCA,NPROW,NPCOL,MYROW,MYCOL,
     *          LRINDX,LCINDX,RSRC,CSRC)
             GAMMA=A((LCINDX-1)*LOCAL_M+LRINDX)
             IF(ABS(GAMMA).LE.TOL) THEN
               GAMMA=EPS
               A((LCINDX-1)*LOCAL_M+LRINDX)=GAMMA
             END IF
           END IF
           ISOUR=IAROW*NPCOL+IACOL
           CALL MPI_BCAST(GAMMA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *         COMM,IERR)
           IPROW=INDXG2P(1,MBLOCK,MYROW,0,NPROW)
           IPCOL=INDXG2P(DIM,NBLOCK,MYCOL,0,NPCOL)
           IF(MYROW.EQ.IPROW.AND.MYCOL.EQ.IPCOL) THEN
             CALL INFOG2L(1,DIM,DESCP,NPROW,NPCOL,MYROW,MYCOL,
     *          LRINDX,LCINDX,RSRCP,CSRCP)
             ALPHA=P_PAR(LRINDX,LCINDX)
           END IF
           ISOUR=IPROW*NPCOL+IPCOL
           CALL MPI_BCAST(ALPHA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *         COMM,IERR)
           IF(ABS(ALPHA).EQ.0) THEN
             BETA=GAMMA*GAMMA*VNOR*VNOR
             IF(ABS(BETA).GT.ONE) THEN
               S= ONE
               C= ZERO
               LAMBDMAX=BETA
             ELSE
               S= ZERO
               C= ONE
               LAMBDMAX= DCMPLX(ONE)
             END IF
           ELSE
             BETA=GAMMA*GAMMA*VNOR*VNOR+ALPHA*ALPHA-ONE
             ETA=BETA/(2*ALPHA)
             IF(ABS(ALPHA).GT.0) THEN
               MU=ETA+ZMUMPS_277(DCMPLX(ETA),DCMPLX(ONE))
             ELSE
               MU=ETA-ZMUMPS_277(DCMPLX(ETA),DCMPLX(ONE))
             END IF
              LAMBDMAX=ALPHA*MU+ONE
             RROOT=ZMUMPS_277(DCMPLX(MU),DCMPLX(ONE))
             S=MU/RROOT
             C=-ONE/RROOT
           END IF
           IF(MYROW.EQ.0) THEN
           CALL PZSCAL(DIM-1,S,V_PAR,1,1,DESCP,1)
           END IF
           VNEW=(C-S*ALPHA)/GAMMA
           IVROW=INDXG2P(1,MBLOCK,MYROW,0,NPROW)
           IVCOL=INDXG2P(DIM,NBLOCK,MYCOL,0,NPCOL)
           IF(MYROW.EQ.IVROW.AND.MYCOL.EQ.IVCOL) THEN
             CALL INFOG2L(1,DIM,DESCP,NPROW,NPCOL,MYROW,MYCOL,
     *         LRINDX,LCINDX,RSRCP,CSRCP)
               V_PAR(LRINDX,LCINDX)=VNEW
           END IF
           VNOR=SQRT(LAMBDMAX)/ABS(GAMMA)
           VNOR= ONE/VNOR
           IF(MYROW.EQ.0) THEN
           CALL PZSCAL(DIM,VNOR,V_PAR,1,1,DESCP,1)
           END IF
            CALL ZMUMPS_112(DIM,LOCAL_M,V_PAR,VAUX,
     *        MYID,NPCOL,NPROW,
     *        MBLOCK,NBLOCK,1,1,COMM)
           CALL PZTRTRS('U','N','N',
     *       DIM,1,
     *       A,1,1,DESCA,VAUX,1,
     *       1,DESCV,IERR)
           call mpi_barrier(comm,ierr)
           IF(MYCOL.EQ.0) THEN
             CALL PDZNRM2(DIM,SMIN,VAUX,1,1,DESCV,1)
           END IF
           ISOUR=0
           CALL MPI_BCAST(SMIN,1,MPI_DOUBLE_PRECISION,ISOUR,
     *         COMM,IERR)
           call mpi_barrier(comm,ierr)
           SMIN= ONE/SMIN
           IF(MYCOL.EQ.0) then
           CALL PZSCAL(DIM,SMIN,VAUX,1,1,DESCV,1)
           END IF
           IVROW=INDXG2P(DIM,MBLOCK,MYROW,0,NPROW)
           IVCOL=INDXG2P(1,NBLOCK,MYCOL,0,NPCOL)
           IF(MYROW.EQ.IVROW.AND.MYCOL.EQ.IVCOL) THEN
             CALL INFOG2L(DIM,1,DESCV,NPROW,NPCOL,MYROW,MYCOL,
     *         LRINDX,LCINDX,RSRCP,CSRCP)
               ALPHA=-ONE/VAUX(LRINDX,LCINDX)
           END IF
      ISOUR=IVROW*NPCOL+IVCOL
      CALL MPI_BCAST(ALPHA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *  COMM,IERR)
      CALL PZSCAL(DIM,ALPHA,VAUX,1,1,DESCV,1)
      IF(MYCOL.EQ.0) THEN
      END IF
 400  CONTINUE
      deallocate(P_PAR)
      deallocate(U_PAR)
      deallocate(V_PAR)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_112(GLOBAL_M,LOCAL_M,DVFROM,DVTO,
     *  MYID,NPCOL,NPROW,
     *  MBLOCK,NBLOCK,ROW,COLUMN,COMM)
      IMPLICIT NONE
      INTEGER GLOBAL_M
      INTEGER LOCAL_M,MYID,NPCOL,NPROW,MBLOCK,NBLOCK,ROW,COLUMN
      INTEGER COMM
      COMPLEX*16 DVFROM(*),DVTO(*)
      INCLUDE 'mpif.h'
      INTEGER SOUR,DEST,ISOUR,IDEST,IBLOCK,JBLOCK
      INTEGER PTR_FROM_ROW,PTR_TO_COL,I,SIZE
      INTEGER COLPROC,ROWPROC,TAG,IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      PARAMETER (TAG=141)
      IBLOCK=(ROW-1)/MBLOCK
      ROWPROC=MOD(IBLOCK,NPCOL)
      JBLOCK=(COLUMN-1)/NBLOCK
      COLPROC=MOD(JBLOCK,NPROW)
      SOUR=NPCOL*ROWPROC
      DEST=COLPROC
      ISOUR=SOUR
      IDEST=DEST
      PTR_FROM_ROW=1
      PTR_TO_COL=1
      DO I=1,GLOBAL_M,MBLOCK
        SIZE=MIN(MBLOCK,IABS(GLOBAL_M-I+1))
        IF(ISOUR.EQ.IDEST.AND.ISOUR.EQ.MYID) THEN
          CALL ZCOPY(SIZE,DVFROM(PTR_FROM_ROW),1,
     *      DVTO(PTR_TO_COL),1)
          PTR_FROM_ROW=PTR_FROM_ROW+NBLOCK
          PTR_TO_COL=PTR_TO_COL+MBLOCK
        ELSE
        IF(MYID.EQ.ISOUR) THEN
          CALL MPI_SEND(DVFROM(PTR_FROM_ROW),SIZE,
     *      MPI_DOUBLE_COMPLEX,
     *      IDEST,TAG,COMM,IERR)
          PTR_FROM_ROW=PTR_FROM_ROW+NBLOCK
        END IF
        IF(MYID.EQ.IDEST) THEN
          CALL MPI_RECV(DVTO(PTR_TO_COL),SIZE,
     *      MPI_DOUBLE_COMPLEX,ISOUR,TAG,COMM,STATUS,IERR)
          PTR_TO_COL=PTR_TO_COL+MBLOCK
        END IF
        END IF
        ISOUR=SOUR+MOD(ISOUR-SOUR+1,NPCOL)
        IDEST=DEST+MOD(IDEST-DEST+NPCOL,NPROW*NPCOL)
      END DO
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_115(GLOBAL_M,GLOBAL_N,
     *  LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,DVFROM,DVTO,ROW,
     *  MYID,MASTER_ROOT,NPROW,NPCOL,COMM)
      IMPLICIT NONE
      INTEGER GLOBAL_M,GLOBAL_N,LOCAL_M,LOCAL_N,MBLOCK,NBLOCK
      INTEGER MYID,MASTER_ROOT,NPROW,NPCOL,COMM,ROW
      COMPLEX*16 DVFROM(LOCAL_M,LOCAL_N)
      COMPLEX*16 DVTO(GLOBAL_M,GLOBAL_N)
      INCLUDE 'mpif.h'
      INTEGER I,J,SIZE_IBLOCK,SIZE_JBLOCK,ISOUR,IROW,ICOL
      INTEGER IBLOCK,JBLOCK,II,JJ,KK,LL,IDIST
      INTEGER PTR_FROM_ROW,PTR_FROM_COL,IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMPLEX*16 DWORK(MBLOCK*NBLOCK)
      LOGICAL ROW_READY
      PTR_FROM_ROW=1
      PTR_FROM_COL=1
      DO I=1,GLOBAL_M,MBLOCK
        IDIST=ROW-I
        SIZE_IBLOCK=MIN(MBLOCK,IABS(GLOBAL_M-I+1))
        ROW_READY= .FALSE.
        DO J=1,GLOBAL_N,NBLOCK
          SIZE_JBLOCK=MIN(NBLOCK,IABS(GLOBAL_N-J+1))
          IBLOCK=I/MBLOCK
          JBLOCK=J/NBLOCK
          IROW=MOD(IBLOCK,NPROW)
          ICOL=MOD(JBLOCK,NPCOL)
          ISOUR=IROW*NPCOL+ICOL
          IF(ISOUR.NE.MASTER_ROOT) THEN
            IF(MYID.EQ.MASTER_ROOT) THEN
              IF(I+SIZE_IBLOCK.GT.ROW) THEN
                CALL MPI_RECV(DWORK, SIZE_IBLOCK*SIZE_JBLOCK,
     *            MPI_INTEGER,ISOUR,139,COMM,STATUS,IERR)
                KK=1
                DO II=I,I+SIZE_JBLOCK-1
                  DVTO(II,1)=DWORK(KK)
                  KK=KK+1
                END DO
              END IF
            ELSE IF(MYID.EQ.ISOUR) THEN
              IF(I+SIZE_IBLOCK.GT.ROW) THEN
                KK=1
                DO II=PTR_FROM_COL,PTR_FROM_COL+SIZE_JBLOCK-1
                  DWORK(KK)=DVFROM(PTR_FROM_ROW+IDIST,II)
                  KK=KK+1
                END DO
                CALL MPI_SEND(DWORK(1),SIZE_IBLOCK*SIZE_JBLOCK,
     *            MPI_INTEGER,MASTER_ROOT,139,COMM,IERR)
              END IF
              ROW_READY= .TRUE.
              PTR_FROM_COL=PTR_FROM_COL+SIZE_JBLOCK
            END IF
          ELSE IF(MYID.EQ.MASTER_ROOT) THEN
            IF(I+SIZE_IBLOCK.GT.ROW) THEN
              KK=PTR_FROM_COL
              DO II=I,I+SIZE_JBLOCK-1
                DVTO(II,1)=DVFROM(PTR_FROM_ROW+IDIST,KK)
                KK=KK+1
              END DO
            END IF
            ROW_READY= .TRUE.
            PTR_FROM_COL=PTR_FROM_COL+SIZE_JBLOCK
          END IF
        END DO
        IF(ROW_READY) THEN
           PTR_FROM_COL=1
           PTR_FROM_ROW=PTR_FROM_ROW+SIZE_IBLOCK
        END IF
      ENDDO
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_108(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *  DVFROM,DVTO,COLUMN,COLPROC,MYCOL,NPROW,NPCOL)
      IMPLICIT NONE
      INTEGER LOCAL_M,LOCAL_N,MBLOCK,NBLOCK
      INTEGER MYCOL,MASTER_ROOT,NPROW,NPCOL,COLUMN,COLPROC
      COMPLEX*16 DVFROM(LOCAL_M,LOCAL_N)
      COMPLEX*16 DVTO(LOCAL_M,1)
      INCLUDE 'mpif.h'
      INTEGER I,J,SIZE_IBLOCK
      INTEGER II,IDIST,JBLOCK,COLLOC
      JBLOCK=COLUMN/NBLOCK
      COLPROC=MOD(JBLOCK,NPROW)
      IF(MYCOL.EQ.COLPROC) THEN
        COLLOC=COLUMN-NBLOCK*JBLOCK
        DO I=1,LOCAL_M,MBLOCK
          DO J=1,LOCAL_N,NBLOCK
            IDIST=COLLOC
           IF(J+NBLOCK.GT.COLLOC) THEN
              SIZE_IBLOCK=MIN(MBLOCK,IABS(LOCAL_M-I+1))
              DO II=I,I+SIZE_IBLOCK-1
                DVTO(II,1)=DVFROM(II,IDIST)
              END DO
              GO TO 50
           END IF
          END DO
        END DO
      END IF
 50   CONTINUE
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_110(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *  DVFROM,TOCOL,DVTO,COLUMN,LENCOL,COLPROC,MYROW,
     *  MYCOL,NPROW,NPCOL,
     *  MYID,COMM)
      IMPLICIT NONE
      INTEGER LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,LENCOL,MYROW
      INTEGER MYCOL,MASTER_ROOT,NPROW,NPCOL,COLUMN,COLPROC
      INTEGER TOCOL,COMM,MYID
      COMPLEX*16 DVTO(LOCAL_M,1)
      COMPLEX*16 DVFROM(LOCAL_M,LOCAL_N)
      INCLUDE 'mpif.h'
      INTEGER I,J,SIZE_IBLOCK,KK
      INTEGER II,IDIST,JBLOCK,COLLOC,JJ,LLENCOL
      INTEGER ISOUR,IDEST,IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMPLEX*16 WK(LOCAL_M)
      JBLOCK=(COLUMN-1)/NBLOCK
      COLPROC=MOD(JBLOCK,NPCOL)
      II=LENCOL/MBLOCK
      JJ=II/NPROW
      KK=MOD(II,NPROW)
      IF(MYROW.LT.KK) THEN
       LLENCOL=(JJ+1)*MBLOCK
      ELSEIF(MYROW.GT.KK) THEN
       LLENCOL=JJ*MBLOCK
      ELSEIF(MYROW.EQ.KK) THEN
       LLENCOL=JJ*MBLOCK+MOD(LENCOL,MBLOCK)
      END IF
      IF(MYCOL.EQ.COLPROC.OR.MYCOL.EQ.TOCOL) THEN
        COLLOC=COLUMN-NBLOCK*JBLOCK
        DO I=1,MIN(LLENCOL,LOCAL_M),MBLOCK
          DO J=1,LOCAL_N,NBLOCK
            IDIST=COLLOC+NBLOCK*(JBLOCK/NPCOL)
           IF(J+NBLOCK.GT.COLLOC) THEN
              SIZE_IBLOCK=MIN(MBLOCK,IABS(LOCAL_M-I+1))
              SIZE_IBLOCK=MIN(SIZE_IBLOCK,IABS(LLENCOL-I+1))
              IF(COLPROC.EQ.TOCOL) THEN
                DO II=I,I+SIZE_IBLOCK-1
                  DVTO(II,1)=DVFROM(II,IDIST)
                END DO
              ELSEIF(MYCOL.EQ.COLPROC) THEN
                IDEST=MYROW*NPROW+TOCOL
                KK=1
                DO II=I,I+SIZE_IBLOCK-1
                  WK(KK)=DVFROM(II,IDIST)
                  KK=KK+1
                END DO
                CALL MPI_SEND(WK,SIZE_IBLOCK,
     *            MPI_DOUBLE_COMPLEX,IDEST,160,COMM,IERR)
              ELSEIF(MYCOL.EQ.TOCOL) THEN
                ISOUR=MYROW*NPROW+COLPROC
                CALL MPI_RECV(WK,SIZE_IBLOCK,
     *            MPI_DOUBLE_COMPLEX,ISOUR,160,COMM,STATUS,IERR)
                KK=1
                DO II=I,I+SIZE_IBLOCK-1
                  DVTO(II,1)=WK(KK)
                  KK=KK+1
                END DO
              END IF
           END IF
            GO TO 50
          END DO
 50       CONTINUE
        END DO
      END IF
      call mpi_barrier(comm,ierr)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_109(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *  DVFROM,FROMCOL,DVTO,COLUMN,LENCOL,COLPROC,MYROW,
     *  MYCOL,NPROW,NPCOL,
     *  MYID,COMM)
      IMPLICIT NONE
      INTEGER LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,LENCOL,MYROW
      INTEGER MYCOL,MASTER_ROOT,NPROW,NPCOL,COLUMN,COLPROC
      INTEGER FROMCOL,COMM,MYID
      COMPLEX*16 DVFROM(LOCAL_M,1)
      COMPLEX*16 DVTO(LOCAL_M,LOCAL_N)
      INCLUDE 'mpif.h'
      INTEGER I,J,SIZE_IBLOCK,KK
      INTEGER II,IDIST,JBLOCK,COLLOC,JJ,LLENCOL
      INTEGER ISOUR,IDEST,IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMPLEX*16 WK(LOCAL_M)
      JBLOCK=(COLUMN-1)/NBLOCK
      COLPROC=MOD(JBLOCK,NPCOL)
      II=LENCOL/MBLOCK
      JJ=II/NPROW
      KK=MOD(II,NPROW)
      IF(MYROW.LT.KK) THEN
       LLENCOL=(JJ+1)*MBLOCK
      ELSEIF(MYROW.GT.KK) THEN
       LLENCOL=JJ*MBLOCK
      ELSEIF(MYROW.EQ.KK) THEN
       LLENCOL=JJ*MBLOCK+MOD(LENCOL,MBLOCK)
      END IF
      IF(MYCOL.EQ.COLPROC.OR.MYCOL.EQ.FROMCOL) THEN
        COLLOC=COLUMN-NBLOCK*JBLOCK
        DO I=1,MIN(LLENCOL,LOCAL_M),MBLOCK
          DO J=1,LOCAL_N,NBLOCK
            IDIST=COLLOC+NBLOCK*(JBLOCK/NPCOL)
           IF(J+NBLOCK.GT.COLLOC) THEN
              SIZE_IBLOCK=MIN(MBLOCK,IABS(LOCAL_M-I+1))
              SIZE_IBLOCK=MIN(SIZE_IBLOCK,IABS(LLENCOL-I+1))
              IF(COLPROC.EQ.FROMCOL) THEN
                DO II=I,I+SIZE_IBLOCK-1
                  DVTO(II,IDIST)=DVFROM(II,1)
                END DO
              ELSEIF(MYCOL.EQ.FROMCOL) THEN
                IDEST=MYROW*NPROW+COLPROC
                KK=1
                DO II=I,I+SIZE_IBLOCK-1
                  WK(KK)=DVFROM(II,1)
                  KK=KK+1
                END DO
                CALL MPI_SEND(WK,SIZE_IBLOCK,
     *            MPI_DOUBLE_COMPLEX,IDEST,160,COMM,IERR)
              ELSEIF(MYCOL.EQ.COLPROC) THEN
                ISOUR=MYROW*NPROW+FROMCOL
                CALL MPI_RECV(WK,SIZE_IBLOCK,
     *            MPI_DOUBLE_COMPLEX,ISOUR,160,COMM,STATUS,IERR)
                KK=1
                DO II=I,I+SIZE_IBLOCK-1
                  DVTO(II,IDIST)=WK(KK)
                  KK=KK+1
                END DO
              END IF
           END IF
            GO TO 50
          END DO
 50       CONTINUE
        END DO
      END IF
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_111(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *  DVFROM,ROWTO,DVTO,ROW,LENROW,ROWPROC,MYROW,MYCOL,
     *  NPROW,NPCOL,MYID,COMM)
      IMPLICIT NONE
      INTEGER LOCAL_M,LOCAL_N,MBLOCK,NBLOCK
      INTEGER MYROW,MASTER_ROOT,NPROW,NPCOL,ROW,ROWPROC
      INTEGER ROWTO,MYCOL,MYID,COMM,LENROW
      COMPLEX*16 DVFROM(LOCAL_M,LOCAL_N)
      COMPLEX*16 DVTO(1,LOCAL_N)
      INCLUDE 'mpif.h'
      INTEGER I,J,SIZE_JBLOCK,IERR,II
      INTEGER JJ,IDIST,IBLOCK,ROWLOC,LLENROW
      INTEGER ISOUR,IDEST,KK,COLTO,JBLOCK
      LOGICAL COLUMN_READY
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMPLEX*16 WK(LOCAL_N)
      II=LENROW/NBLOCK
      JJ=II/NPCOL
      KK=MOD(II,NPCOL)
      IF(MYCOL.LT.KK) THEN
        LLENROW=(JJ+1)*NBLOCK
      ELSEIF(MYCOL.GT.KK) THEN
        LLENROW=JJ*NBLOCK
      ELSEIF(MYCOL.EQ.KK) THEN
        LLENROW=JJ*NBLOCK+MOD(LENROW,NBLOCK)
      END IF
      IBLOCK=(ROW-1)/MBLOCK
      ROWPROC=MOD(IBLOCK,NPROW)
      IF(MYROW.EQ.ROWPROC.OR.MYROW.EQ.ROWTO) THEN
        ROWLOC=ROW-MBLOCK*IBLOCK
        DO J=1,MIN(LOCAL_N,LLENROW),NBLOCK
          JBLOCK=(J-1)/NBLOCK
          COLTO=MOD(JBLOCK,NPCOL)
          COLUMN_READY=.FALSE.
          DO I=1,LOCAL_M,MBLOCK
            IDIST=ROWLOC+MBLOCK*(IBLOCK/NPROW)
            IF(I+MBLOCK.GT.IDIST) THEN
              SIZE_JBLOCK=MIN(NBLOCK,IABS(LOCAL_N-J+1))
              SIZE_JBLOCK=MIN(SIZE_JBLOCK,IABS(LLENROW-J+1))
              IF(ROWPROC.EQ.ROWTO) THEN
                DO JJ=J,J+SIZE_JBLOCK-1
                  DVTO(1,JJ)=DVFROM(IDIST,JJ)
                END DO
              ELSEIF(MYROW.EQ.ROWPROC) THEN
                IDEST=ROWTO*NPROW+MYCOL
                KK=1
                DO JJ=J,J+SIZE_JBLOCK-1
                  WK(KK)=DVFROM(IDIST,JJ)
                  KK=KK+1
                END DO
                CALL MPI_SEND(WK,SIZE_JBLOCK,
     *            MPI_DOUBLE_COMPLEX,IDEST,170,COMM,IERR)
                COLUMN_READY=.TRUE.
              ELSEIF(MYROW.EQ.ROWTO) THEN
                ISOUR=ROWPROC*NPROW+MYCOL
                CALL MPI_RECV(WK,SIZE_JBLOCK,
     *            MPI_DOUBLE_COMPLEX,ISOUR,170,COMM,STATUS,IERR)
                KK=1
                DO JJ=J,J+SIZE_JBLOCK-1
                  DVTO(1,JJ)=WK(KK)
                  KK=KK+1
                END DO
              END IF
              GO TO 50
            END IF
          END DO
 50       CONTINUE
        END DO
      END IF
      RETURN
      END
#endif
#ifdef try_null_space_not_tested
      SUBROUTINE ZMUMPS_114(GLOBAL_M,GLOBAL_N,
     *  LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,DVFROM,DVTO,COLUMN,
     *  MYID,MASTER_ROOT,NPROW,NPCOL,COMM)
      IMPLICIT NONE
      INTEGER GLOBAL_M,GLOBAL_N,LOCAL_M,LOCAL_N,MBLOCK,NBLOCK
      INTEGER MYID,MASTER_ROOT,NPROW,NPCOL,COMM,COLUMN
      COMPLEX*16 DVFROM(LOCAL_M,LOCAL_N)
      COMPLEX*16 DVTO(GLOBAL_M,GLOBAL_N)
      INCLUDE 'mpif.h'
      INTEGER I,J,SIZE_IBLOCK,SIZE_JBLOCK,ISOUR,IROW,ICOL
      INTEGER IBLOCK,JBLOCK,II,JJ,KK,LL,IDIST
      INTEGER PTR_FROM_ROW,PTR_FROM_COL,IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMPLEX*16 DWORK(MBLOCK*NBLOCK)
      LOGICAL COLUMN_READY
      PTR_FROM_ROW=1
      PTR_FROM_COL=1
      DO J=1,GLOBAL_N,NBLOCK
        IDIST=COLUMN-J
        SIZE_JBLOCK=MIN(NBLOCK,IABS(GLOBAL_N-J+1))
        COLUMN_READY= .FALSE.
        DO I=1,GLOBAL_M,MBLOCK
          SIZE_IBLOCK=MIN(MBLOCK,IABS(GLOBAL_M-I+1))
          IBLOCK=I/MBLOCK
          JBLOCK=J/NBLOCK
          IROW=MOD(IBLOCK,NPROW)
          ICOL=MOD(JBLOCK,NPCOL)
          ISOUR=IROW*NPCOL+ICOL
          IF(ISOUR.NE.MASTER_ROOT) THEN
            IF(MYID.EQ.MASTER_ROOT) THEN
              IF(J+SIZE_JBLOCK.GT.COLUMN) THEN
                CALL MPI_RECV(DWORK(1),SIZE_IBLOCK*SIZE_JBLOCK,
     *            MPI_DOUBLE_COMPLEX,ISOUR,129,COMM,STATUS,IERR)
                KK=1
                DO II=I,I+SIZE_IBLOCK-1
                  DVTO(II,1)=DWORK(KK)
                  KK=KK+1
                END DO
              END IF
            ELSE IF(MYID.EQ.ISOUR) THEN
              IF(J+SIZE_JBLOCK.GT.COLUMN) THEN
                KK=1
                DO II=PTR_FROM_ROW,PTR_FROM_ROW+SIZE_IBLOCK-1
                  DWORK(KK)=DVFROM(II,PTR_FROM_COL+IDIST)
                  KK=KK+1
                END DO
                CALL MPI_SEND(DWORK(1),SIZE_IBLOCK*SIZE_JBLOCK,
     *            MPI_DOUBLE_COMPLEX,MASTER_ROOT,129,COMM,IERR)
              END IF
              COLUMN_READY= .TRUE.
              PTR_FROM_ROW=PTR_FROM_ROW + SIZE_IBLOCK
            END IF
          ELSE IF(MYID.EQ.MASTER_ROOT) THEN
            IF(J+SIZE_JBLOCK.GT.COLUMN) THEN
              KK=PTR_FROM_ROW
              DO II=I,I+SIZE_IBLOCK-1
                DVTO(II,1)=DVFROM(KK,PTR_FROM_COL+IDIST)
                KK=KK+1
              END DO
            END IF
            COLUMN_READY= .TRUE.
            PTR_FROM_ROW=PTR_FROM_ROW+SIZE_IBLOCK
          END IF
        END DO
        IF(COLUMN_READY) THEN
           PTR_FROM_ROW=1
           PTR_FROM_COL=PTR_FROM_COL + SIZE_JBLOCK
        END IF
      ENDDO
      RETURN
      END
#endif
#ifdef try_null_space_not_tested
      SUBROUTINE ZMUMPS_154(MYID, M, N, ASEQ,
     *                    LOCAL_M, LOCAL_N,
     *                    MBLOCK, NBLOCK,
     *                    APAR,
     *                    MASTER_ROOT,
     *                    NPROW, NPCOL,
     *                    COMM,column)
      IMPLICIT NONE
      INTEGER MYID, MASTER_ROOT, COMM
      INTEGER M, N,column
      INTEGER NPROW, NPCOL
      INTEGER LOCAL_M, LOCAL_N
      INTEGER MBLOCK, NBLOCK
      COMPLEX*16 APAR( LOCAL_M, LOCAL_N )
      COMPLEX*16 ASEQ( M, N )
      INCLUDE 'mpif.h'
      INTEGER I, J, SIZE_IBLOCK, SIZE_JBLOCK, ISOUR, IROW, ICOL
      INTEGER IBLOCK, JBLOCK, II, JJ, KK
      INTEGER IAPAR, JAPAR, IERR,idist
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMPLEX*16 WK( MBLOCK * NBLOCK )
      LOGICAL JUPDATE
        IAPAR = 1
        JAPAR = 1
        DO J = 1, N, NBLOCK
          idist=column-j
          SIZE_JBLOCK = NBLOCK
          IF ( J + NBLOCK > N ) THEN
            SIZE_JBLOCK = N - J + 1
          END IF
          JUPDATE = .FALSE.
          DO I = 1, M, MBLOCK
            SIZE_IBLOCK = MBLOCK
            IF ( I + MBLOCK > M ) THEN
              SIZE_IBLOCK = M - I + 1
            END IF
            IBLOCK = I / MBLOCK
            JBLOCK = J / NBLOCK
            IROW = MOD ( IBLOCK, NPROW )
            ICOL = MOD ( JBLOCK, NPCOL )
            ISOUR = IROW * NPCOL + ICOL
            IF ( ISOUR .NE. MASTER_ROOT ) THEN
              IF ( MYID .EQ. MASTER_ROOT ) THEN
                if(j+size_jblock.gt.column) then
                CALL MPI_RECV( WK, SIZE_IBLOCK,
     *                         MPI_DOUBLE_COMPLEX,
     *                         ISOUR, 128, COMM, STATUS, IERR )
                KK=1
                DO JJ=1,1
                DO II=I,I+SIZE_IBLOCK-1
                  ASEQ(II,JJ)=WK(KK)
                  KK=KK+1
                END DO
                END DO
                end if
              ELSE IF ( MYID .EQ. ISOUR ) THEN
                if(j+size_jblock.gt.column) then
                KK=1
                do JJ=japar+idist,japar+idist
                DO II=IAPAR,IAPAR+SIZE_IBLOCK-1
                  WK(KK)=APAR(II,JJ)
                  KK=KK+1
                END DO
                END DO
                CALL MPI_SEND( WK( 1 ),
     *                         SIZE_IBLOCK,
     *                         MPI_DOUBLE_COMPLEX,
     *                         MASTER_ROOT,128,COMM,IERR)
                end if
                JUPDATE = .TRUE.
                IAPAR = IAPAR + SIZE_IBLOCK
              END IF
            ELSE IF ( MYID.EQ. MASTER_ROOT ) THEN
              if(j+size_jblock.gt.column) then
              ASEQ(I:I+SIZE_IBLOCK-1,1:1)
     *        = APAR( IAPAR:IAPAR+SIZE_IBLOCK-1,
     *                JAPAR+idist:JAPAR+idist )
              end if
              JUPDATE = .TRUE.
              IAPAR = IAPAR + SIZE_IBLOCK
            END IF
          END DO
          IF ( JUPDATE ) THEN
            IAPAR = 1
            JAPAR = JAPAR + SIZE_JBLOCK
          END IF
        END DO
      RETURN
      END SUBROUTINE ZMUMPS_154
#endif
#ifdef try_null_space_not_tested
      SUBROUTINE ZMUMPS_172(GLOBAL_M,GLOBAL_N,
     *  LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,IVFROM,IVTO,
     *  MYID,MASTER_ROOT,NPROW,NPCOL,COMM)
      IMPLICIT NONE
      INTEGER GLOBAL_M,GLOBAL_N,LOCAL_M,LOCAL_N,MBLOCK,NBLOCK
      INTEGER MYID,MASTER_ROOT,NPROW,NPCOL,COMM
      INTEGER IVFROM(LOCAL_M,LOCAL_N)
      INTEGER IVTO(GLOBAL_M,GLOBAL_N)
      INCLUDE 'mpif.h'
      INTEGER I,J,SIZE_IBLOCK,SIZE_JBLOCK,IDEST,IROW,ICOL
      INTEGER IBLOCK,II,JJ,KK,LL
      INTEGER PTR_FROM_ROW,PTR_FROM_COL,IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      INTEGER IWORK(MBLOCK*NBLOCK)
      LOGICAL COLUMN_READY
      PTR_FROM_ROW=1
      PTR_FROM_COL=1
      call mpi_barrier(comm,ierr)
      DO J=1,GLOBAL_N,NBLOCK
        SIZE_JBLOCK=MIN(NBLOCK,IABS(GLOBAL_N-J+1))
        COLUMN_READY= .FALSE.
        DO I=1,GLOBAL_M,MBLOCK
          SIZE_IBLOCK=MIN(MBLOCK,IABS(GLOBAL_M-I+1))
          IBLOCK=I/MBLOCK
          IDEST=MOD(IBLOCK,NPCOL)
          IF(IDEST.NE.MASTER_ROOT) THEN
            IF(MYID.EQ.MASTER_ROOT) THEN
              KK=1
              DO JJ=J,J+SIZE_JBLOCK-1
                DO II=I,I+SIZE_IBLOCK-1
                  IWORK(KK)=IVFROM(II,JJ)
                  KK=KK+1
                END DO
              END DO
              CALL MPI_SEND(IWORK, SIZE_IBLOCK*SIZE_JBLOCK,
     *          MPI_INTEGER,IDEST,129,COMM,IERR)
            ELSE IF(MYID.EQ.IDEST) THEN
              CALL MPI_RECV(IWORK(1),SIZE_IBLOCK*SIZE_JBLOCK,
     *          MPI_INTEGER,MASTER_ROOT,129,COMM,STATUS,IERR)
              KK=1
              DO JJ=PTR_FROM_COL,PTR_FROM_COL+SIZE_JBLOCK-1
                DO II=PTR_FROM_ROW,PTR_FROM_ROW+SIZE_IBLOCK-1
                  IVTO(II,JJ)=IWORK(KK)
                  KK=KK+1
                END DO
              END DO
              COLUMN_READY= .TRUE.
              PTR_FROM_ROW=PTR_FROM_ROW+SIZE_IBLOCK
            END IF
          ELSE IF(MYID.EQ.MASTER_ROOT) THEN
            KK=PTR_FROM_ROW
            LL=PTR_FROM_COL
            DO II=I,I+SIZE_IBLOCK-1
              DO JJ=j,J+SIZE_JBLOCK-1
                IVTO(KK,LL)=IVFROM(II,JJ)
                LL=LL+1
              END DO
              KK=KK+1
              LL=PTR_FROM_COL
            END DO
            COLUMN_READY= .TRUE.
            PTR_FROM_ROW=PTR_FROM_ROW+SIZE_IBLOCK
          END IF
        END DO
        IF(COLUMN_READY) THEN
           PTR_FROM_ROW=1
           PTR_FROM_COL=PTR_FROM_COL + SIZE_JBLOCK
        END IF
      ENDDO
      call mpi_barrier(comm,ierr)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_173(GLOBAL_M,GLOBAL_N,
     *  LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,IVFROM,IVTO,
     *  MYID,MASTER_ROOT,NPROW,NPCOL,COMM)
      IMPLICIT NONE
      INTEGER GLOBAL_M,GLOBAL_N,LOCAL_M,LOCAL_N,MBLOCK,NBLOCK
      INTEGER MYID,MASTER_ROOT,NPROW,NPCOL,COMM
      INTEGER IVFROM(LOCAL_M,LOCAL_N)
      INTEGER IVTO(GLOBAL_M,GLOBAL_N)
      INCLUDE 'mpif.h'
      INTEGER I,J,SIZE_IBLOCK,SIZE_JBLOCK,ISOUR,IROW,ICOL
      INTEGER IBLOCK,II,JJ,KK,LL
      INTEGER PTR_FROM_ROW,PTR_FROM_COL,IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      INTEGER IWORK(MBLOCK*NBLOCK)
      LOGICAL COLUMN_READY
      PTR_FROM_ROW=1
      PTR_FROM_COL=1
      DO J=1,GLOBAL_N,NBLOCK
        SIZE_JBLOCK=MIN(NBLOCK,IABS(GLOBAL_N-J+1))
        COLUMN_READY= .FALSE.
        DO I=1,GLOBAL_M,MBLOCK
          SIZE_IBLOCK=MIN(MBLOCK,IABS(GLOBAL_M-I+1))
          IBLOCK=I/MBLOCK
          ISOUR=MOD(IBLOCK,NPCOL)
          IF(ISOUR.NE.MASTER_ROOT) THEN
            IF(MYID.EQ.MASTER_ROOT) THEN
              CALL MPI_RECV(IWORK, SIZE_IBLOCK*SIZE_JBLOCK,
     *          MPI_INTEGER,ISOUR,129,COMM,STATUS,IERR)
              KK=1
              DO JJ=J,J+SIZE_JBLOCK-1
                DO II=I,I+SIZE_IBLOCK-1
                  IVTO(II,JJ)=IWORK(KK)
                  KK=KK+1
                END DO
              END DO
            ELSE IF(MYID.EQ.ISOUR) THEN
              KK=1
              DO JJ=PTR_FROM_COL,PTR_FROM_COL+SIZE_JBLOCK-1
                DO II=PTR_FROM_ROW,PTR_FROM_ROW+SIZE_IBLOCK-1
                  IWORK(KK)=IVFROM(II,JJ)
                  KK=KK+1
                END DO
              END DO
              CALL MPI_SEND(IWORK(1),SIZE_IBLOCK*SIZE_JBLOCK,
     *          MPI_INTEGER,MASTER_ROOT,129,COMM,IERR)
              COLUMN_READY= .TRUE.
              PTR_FROM_ROW=PTR_FROM_ROW + SIZE_IBLOCK
            END IF
          ELSE IF(MYID.EQ.MASTER_ROOT) THEN
            KK=PTR_FROM_ROW
            LL=PTR_FROM_COL
            DO II=I,I+SIZE_IBLOCK-1
              DO JJ=j,J+SIZE_JBLOCK-1
                IVTO(II,JJ)=IVFROM(KK,LL)
                LL=LL+1
              END DO
              KK=KK+1
              LL=PTR_FROM_COL
            END DO
            COLUMN_READY= .TRUE.
            PTR_FROM_ROW=PTR_FROM_ROW+SIZE_IBLOCK
          END IF
        END DO
        IF(COLUMN_READY) THEN
           PTR_FROM_ROW=1
           PTR_FROM_COL=PTR_FROM_COL + SIZE_JBLOCK
        END IF
      ENDDO
      RETURN
      END
#endif
#ifdef try_null_space_not_tested
      SUBROUTINE ZMUMPS_155(MYID, M, N, AISEQ,
     *                    LOCAL_M, LOCAL_N,
     *                    MBLOCK, NBLOCK,
     *                    AIPAR,
     *                    MASTER_ROOT,
     *                    NPROW, NPCOL,
     *                    COMM)
      IMPLICIT NONE
      INTEGER MYID, MASTER_ROOT, COMM
      INTEGER M, N
      INTEGER NPROW, NPCOL
      INTEGER LOCAL_M, LOCAL_N
      INTEGER MBLOCK, NBLOCK
      INTEGER AIPAR( LOCAL_M, LOCAL_N )
      INTEGER AISEQ( M, N )
      INCLUDE 'mpif.h'
      INTEGER I, J, SIZE_IBLOCK, SIZE_JBLOCK, ISOUR, IROW, ICOL
      INTEGER IBLOCK, JBLOCK, II, JJ, KK
      INTEGER IAPAR, JAPAR, IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      INTEGER WK( MBLOCK * NBLOCK )
      LOGICAL JUPDATE
        IAPAR = 1
        JAPAR = 1
         j=1
          SIZE_JBLOCK = NBLOCK
          IF ( J + NBLOCK > N ) THEN
            SIZE_JBLOCK = N - J + 1
          END IF
          size_jblock=min(nblock,abs(n-j+1))
          JUPDATE = .FALSE.
          DO I = 1, M, MBLOCK
            SIZE_IBLOCK = MBLOCK
            IF ( I + MBLOCK > M ) THEN
              SIZE_IBLOCK = M - I + 1
            END IF
            size_iblock=min(mblock,abs(m-i+1))
            IBLOCK = I / MBLOCK
            JBLOCK = J / NBLOCK
            IROW = MOD ( IBLOCK, NPROW )
            ICOL = MOD ( JBLOCK, NPCOL )
            ISOUR = IROW * NPCOL + ICOL
             isour=mod(iblock,npcol)
            IF ( ISOUR .NE. MASTER_ROOT ) THEN
              IF ( MYID .EQ. MASTER_ROOT ) THEN
                CALL MPI_RECV( WK, SIZE_IBLOCK*SIZE_JBLOCK,
     *                         MPI_INTEGER,
     *                         ISOUR, 129, COMM, STATUS, IERR )
                KK=1
                DO JJ=J,J+SIZE_JBLOCK-1
                DO II=I,I+SIZE_IBLOCK-1
                  AISEQ(II,JJ)=WK(KK)
                  KK=KK+1
                END DO
                END DO
              ELSE IF ( MYID .EQ. ISOUR ) THEN
                KK=1
                DO JJ=JAPAR,JAPAR+SIZE_JBLOCK-1
                DO II=IAPAR,IAPAR+SIZE_IBLOCK-1
                  WK(KK)=AIPAR(II,JJ)
                  KK=KK+1
                END DO
                END DO
                CALL MPI_SEND( WK( 1 ),
     *                         SIZE_IBLOCK*SIZE_JBLOCK,
     *                         MPI_INTEGER,
     *                         MASTER_ROOT,129,COMM,IERR)
                JUPDATE = .TRUE.
                IAPAR = IAPAR + SIZE_IBLOCK
              END IF
            ELSE IF ( MYID.EQ. MASTER_ROOT ) THEN
              AISEQ(I:I+SIZE_IBLOCK-1,J:J+SIZE_JBLOCK-1)
     *        = AIPAR( IAPAR:IAPAR+SIZE_IBLOCK-1,
     *                JAPAR:JAPAR+SIZE_JBLOCK-1 )
              JUPDATE = .TRUE.
              IAPAR = IAPAR + SIZE_IBLOCK
            END IF
          END DO
      RETURN
      END SUBROUTINE ZMUMPS_155
#endif
      SUBROUTINE ZMUMPS_333(N,PERM,X,RN01)
      INTEGER N,PERM(N),I
      COMPLEX*16 RN01(N),X(N)
      DO 100 I=1,N
      RN01(PERM(I))=X(I)
100   CONTINUE
      DO 200 I=1,N
      X(I)=RN01(I)
200   CONTINUE
      RETURN
      END
#if defined(try_null_space_not_tested)
      DOUBLE PRECISION FUNCTION ZMUMPS_277(X,Y)
      COMPLEX*16 X,Y
      COMPLEX*16 TEMP,TEMP1
      DOUBLE PRECISION ZERO
      PARAMETER(ZERO=0.0D0)
      TEMP=MAX(ABS(X),ABS(Y))
      IF (TEMP.EQ.ZERO) THEN
        ZMUMPS_277= ZERO
      ELSE
        TEMP1=ABS(X)/TEMP
        ZMUMPS_277=TEMP1*TEMP1
        TEMP1=ABS(Y)/TEMP
        ZMUMPS_277=ZMUMPS_277+TEMP1*TEMP1
        ZMUMPS_277=TEMP*SQRT(ZMUMPS_277)
      END IF
      RETURN
      END
#endif
#ifdef try_null_space_not_tested
      SUBROUTINE ZMUMPS_261
     *  ( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     $                    IB, JB, DESCB, INFO )
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX*16   A( * ), B( * )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX*16   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      LOGICAL            NOTRAN
      INTEGER            IAROW, IBROW, ICOFFA, ICTXT, IROFFA, IROFFB,
     $                   MYCOL, MYROW, NPCOL, NPROW
      INTEGER            DESCIP( DLEN_ ), IDUM1( 1 ), IDUM2( 1 )
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, PCHK2MAT,
     $                   PZLAPIV, PZTRSM, PXERBLA
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, LSAME, NUMROC
      INTRINSIC          ICHAR, MOD
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(700+CTXT_)
      ELSE
         NOTRAN = LSAME( TRANS, 'N' )
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 2, NRHS, 3, IB, JB, DESCB, 12, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                       NPROW )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
               INFO = -1
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(700+NB_)
            ELSE IF( IROFFB.NE.0 .OR. IBROW.NE.IAROW ) THEN
               INFO = -10
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(1200+NB_)
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -(1200+CTXT_)
            END IF
         END IF
         IF( NOTRAN ) THEN
            IDUM1( 1 ) = ICHAR( 'N' )
         ELSE IF( LSAME( TRANS, 'T' ) ) THEN
            IDUM1( 1 ) = ICHAR( 'T' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'C' )
         END IF
         IDUM2( 1 ) = 1
         CALL PCHK2MAT( N, 2, N, 2, IA, JA, DESCA, 7, N, 2, NRHS, 3,
     $                  IB, JB, DESCB, 12, 1, IDUM1, IDUM2, INFO )
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZGETRS', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
      CALL DESCSET( DESCIP, DESCA( M_ ) + DESCA( MB_ )*NPROW, 1,
     $              DESCA( MB_ ), 1, DESCA( RSRC_ ), MYCOL, ICTXT,
     $              DESCA( MB_ ) + NUMROC( DESCA( M_ ), DESCA( MB_ ),
     $              MYROW, DESCA( RSRC_ ), NPROW ) )
      IF( NOTRAN ) THEN
         CALL PZTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $                NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
      ELSE
         write(*,*) ' not to be used'
         stop
      END IF
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_169(MYID,MYROW,MYCOL,LOCAL_M,LOCAL_N,
     *  MBLOCK,NBLOCK,NPROW,NPCOL,
     *  N,A,
     *  DESCA,DEFICIENCY,MAXSTEPS,U,DESCU,
     *  PERM_COL,ELROW,ELELG,DL,
     *  EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WORK,SMIN)
      IMPLICIT NONE
      INTEGER MYID,MYROW,MYCOL,LOCAL_M,LOCAL_N
      INTEGER MBLOCK,NBLOCK,NPROW,NPCOL
      INTEGER N,LDA,DEFICIENCY,MAXSTEPS
      COMPLEX*16 SMIN
      COMPLEX*16 A(*),U(LOCAL_M+1,1)
      INTEGER NLUPDATE,NUUPDATE,PTREL(*),PTREU(*)
      INTEGER PERM_COL(*)
      INTEGER DESCA(*),DESCU(*)
      INTEGER ELROW(*),EUROW(*)
      COMPLEX*16 ELELG(*),EUELG(*),DL(*),WORK(*)
      INTEGER I
      INTEGER NSTEPS,LRINDX,LCINDX
      INTEGER RSRCP,CSRCP,IVROW,IVCOL,DIM
      INTEGER INDXG2P
      COMPLEX*16 ALPHA
      DOUBLE PRECISION ONE,UNORM
      PARAMETER(ONE=1.0D0)
      DOUBLE PRECISION DZNRM2
      DIM=N-DEFICIENCY
      LDA=DESCA(9)
      NSTEPS=0
      CALL ZMUMPS_117(N-DEFICIENCY,ONE,U,1)
      DO I=N-DEFICIENCY+1,N
        IVROW=INDXG2P(I,MBLOCK,MYROW,0,NPROW)
        IVCOL=INDXG2P(1,NBLOCK,MYCOL,0,NPCOL)
         IF(MYROW.EQ.IVROW.AND.MYCOL.EQ.IVCOL) THEN
           CALL INFOG2L(I,1,DESCU,NPROW,NPCOL,MYROW,MYCOL,
     *       LRINDX,LCINDX,RSRCP,CSRCP)
             ALPHA=U(LRINDX,LCINDX)
             U(LRINDX,LCINDX)=0.0D0
         END IF
      END DO
 100  CONTINUE
      IF(MYCOL.EQ.0) THEN
        CALL PDZNRM2(DIM,UNORM,U,1,1,DESCU,1)
      END IF
      IF(MYCOL.EQ.0) THEN
        ALPHA= ONE/UNORM
        CALL PZSCAL(DIM,ALPHA,U,1,1,DESCU,1)
      END IF
      NSTEPS=NSTEPS+1
      IF(NSTEPS.LT.MAXSTEPS) GO TO 100
      UNORM=DZNRM2(N,U,1)
      CALL ZSCAL(N,ONE/UNORM,U,1)
      SMIN=SQRT(ONE/UNORM)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_174(N,IALPHA,IX)
      INTEGER N,IALPHA
      INTEGER IX(*)
      INTEGER I
      DO I=1,N
        IX(I)=IALPHA+I
      END DO
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_98(N,A,LDA,DEFICIENCY,TOL)
      INTEGER N,LDA,DEFICIENCY
      DOUBLE PRECISION TOL
      COMPLEX*16 A(*)
      INTEGER I,J,IND,JSTRT,JSTOP
      COMPLEX*16 TEMP
      IND=(N-1)*LDA+N
      DO I=N,1,-1
        JSTRT=IND
        JSTOP=IND+(N-I)*N
        DO J=JSTRT,JSTOP,N
          TEMP=A(J)
          IF(ABS(TEMP).GT.TOL) THEN
            DEFICIENCY=N-I
            GO TO 100
          END IF
        END DO
        IND=IND-LDA-1
      END DO
      DEFICIENCY=N
 100  CONTINUE
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_177(J,AL,DIML,SHIFT,ELCOUNT,ELROW,ELELG)
      INTEGER J,DIML,ELCOUNT,SHIFT
      INTEGER ELROW(*)
      COMPLEX*16 ELELG(*),AL(DIML,*)
      COMPLEX*16 TEMP1,TEMP2,ELG
      DOUBLE PRECISION ZERO
      PARAMETER(ZERO=0.0D0)
      TEMP1=AL(J,J)
      TEMP2=AL(J,J+1)
      IF(TEMP1.EQ.ZERO.AND.TEMP2.EQ.ZERO) THEN
      ELSEIF(ABS(TEMP1).GT.ABS(TEMP2)) THEN
        CALL ZMUMPS_157(AL(J,J),AL(J,J+1),ELG)
        CALL ZMUMPS_29(DIML-J+1,ELG,AL(J,J),1,AL(J,J+1),1)
        ELCOUNT=ELCOUNT+1
        ELROW(ELCOUNT)=J+SHIFT
        ELELG(ELCOUNT)=ELG
      ELSE
        CALL ZMUMPS_157(AL(J,J+1),AL(J,J),ELG)
        CALL ZMUMPS_30(DIML-J+1,ELG,AL(J,J),1,AL(J,J+1),1)
        ELCOUNT=ELCOUNT+1
        ELROW(ELCOUNT)=-(J+SHIFT)
        ELELG(ELCOUNT)=ELG
      END IF
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_331(J,AU,DIMU,SHIFT,EUCOUNT,EUROW,EUELG)
      INTEGER J,DIMU,EUCOUNT,SHIFT
      INTEGER EUROW(*)
      COMPLEX*16 EUELG(*),AU(DIMU,*)
      COMPLEX*16 TEMP1,TEMP2,ELG
      DOUBLE PRECISION ZERO
      PARAMETER(ZERO=0.0D0)
      TEMP1=AU(J,J)
      TEMP2=AU(J+1,J)
      IF(TEMP1.EQ.ZERO.AND.TEMP2.EQ.ZERO) THEN
      ELSEIF(ABS(TEMP1).GT.ABS(TEMP2)) THEN
        CALL ZMUMPS_157(AU(J,J),AU(J+1,J),ELG)
        CALL ZMUMPS_29(DIMU-J+1,ELG,AU(J,J),DIMU,AU(J+1,J),DIMU)
        EUCOUNT=EUCOUNT+1
        EUROW(EUCOUNT)=J+SHIFT
        EUELG(EUCOUNT)=ELG
      ELSE
        CALL ZMUMPS_157(AU(J+1,J),AU(J,J),ELG)
        CALL ZMUMPS_30(DIMU-J+1,ELG,AU(J,J),DIMU,AU(J+1,J),DIMU)
        EUCOUNT=EUCOUNT+1
        EUROW(EUCOUNT)=-(J+SHIFT)
        EUELG(EUCOUNT)=ELG
      END IF
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_304(N,AA,LDA,
     *  DL,NLUPDATE,PTREL,ELROW,ELELG,
     *  NUUPDATE,PTREU,EUROW,EUELG,PERM,PERM_ROW,PERM_COL,
     *  RHS,WORK)
      INTEGER N,LDA,NLUPDATE,NUUPDATE
      INTEGER ELROW(*),EUROW(*),PTREL(*),PTREU(*),PERM(*)
      INTEGER PERM_COL(*),PERM_ROW(*)
      COMPLEX*16 ELELG(*),EUELG(*),RHS(*),WORK(*),DL(*),AA(*)
      INTEGER II,IND,J,K,JSTRT,JSTOP
      COMPLEX*16 ELG
      COMPLEX*16 ONE
      PARAMETER(ONE=1.0D0)
      CALL ZLASWP(1,RHS,LDA,1,N,PERM,1)
      CALL ZMUMPS_334(N,PERM_ROW,RHS,WORK)
      CALL ZTRSM('Left','Lower','No transpose','Unit',N,1,
     *  ONE,AA,LDA,RHS,LDA)
      DO II=NLUPDATE,1,-1
        JSTRT=PTREL(II)
        JSTOP=PTREL(II+1)-1
        K=ELROW(JSTRT)
        IF(K.EQ.0) THEN
          K=ELROW(JSTOP)
        END IF
        IND=JSTRT
        DO J=K,K+JSTOP-JSTRT
          RHS(J)=RHS(J)*DL(IND)
          IND=IND+1
        END DO
        DO J=JSTOP,JSTRT,-1
          K=ELROW(J)
          IF(J.EQ.JSTOP) K=0
          IF(K.GT.0) THEN
            ELG=ELELG(J)
            CALL ZMUMPS_28(1,ELG,RHS(K),N,RHS(K+1),N)
          ELSEIF(K.LT.0) THEN
            ELG=ELELG(J)
            CALL ZMUMPS_30(1,ELG,RHS(-K),N,RHS(-K+1),N)
          ELSE
          END IF
        END DO
      END DO
      DO II=1,NUUPDATE
        JSTRT=PTREU(II)
        JSTOP=PTREU(II+1)-1
        DO J=JSTRT,JSTOP
          K=EUROW(J)
          ELG=EUELG(J)
          IF(K.GT.0) THEN
            CALL ZMUMPS_29(1,ELG,RHS(K),N,RHS(K+1),N)
          ELSE
            CALL ZMUMPS_30(1,ELG,RHS(-K),N,RHS(-K+1),N)
          END IF
        END DO
      END DO
      CALL ZTRSM('Left','Upper','No transpose','Non-unit',N,
     *  1,ONE,AA,LDA,RHS,LDA)
      CALL ZMUMPS_333(N,PERM_COL,RHS,WORK)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_306(N,AA,LDA,
     *  DL,NLUPDATE,PTREL,ELROW,ELELG,
     *  NUUPDATE,PTREU,EUROW,EUELG,PERM,PERM_ROW,PERM_COL,
     *  RHS,WORK)
      INTEGER N,LDA,NLUPDATE,NUUPDATE
      INTEGER ELROW(*),EUROW(*),PTREL(*),PTREU(*),PERM(*)
      INTEGER PERM_COL(*),PERM_ROW(*)
      COMPLEX*16 ELELG(*),EUELG(*),RHS(*),WORK(*),DL(*),AA(*)
      INTEGER II,IND,J,K,JSTRT,JSTOP
      COMPLEX*16 ELG
      COMPLEX*16 ONE
      PARAMETER(ONE=1.0D0)
      CALL ZMUMPS_334(N,PERM_COL,RHS,WORK)
      CALL ZTRSM('Left','Upper','Transpose','Non-unit',N,
     *  1,ONE,AA,LDA,RHS,LDA)
      DO II=NUUPDATE,1,-1
        JSTRT=PTREU(II)
        JSTOP=PTREU(II+1)-1
        DO J=JSTOP,JSTRT,-1
          K=EUROW(J)
          ELG=EUELG(J)
          IF(K.GT.0) THEN
            CALL ZMUMPS_28(1,ELG,RHS(K),N,RHS(K+1),N)
          ELSEIF(K.LT.0) THEN
            CALL ZMUMPS_30(1,ELG,RHS(-K),N,RHS(-K+1),N)
          END IF
        END DO
      END DO
      DO II=1,NLUPDATE
        JSTRT=PTREL(II)
        JSTOP=PTREL(II+1)-1
        DO J=JSTRT,JSTOP
          K=ELROW(J)
          IF(J.EQ.JSTOP) K=0
          IF(K.GT.0) THEN
            ELG=ELELG(J)
            CALL ZMUMPS_29(1,ELG,RHS(K),N,RHS(K+1),N)
          ELSEIF(K.LT.0) THEN
            ELG=ELELG(J)
            CALL ZMUMPS_30(1,ELG,RHS(-K),N,RHS(-K+1),N)
          END IF
        END DO
        K=ELROW(JSTRT)
        IF(K.EQ.0) THEN
          K=ELROW(JSTOP)
        END IF
        IND=JSTRT
        DO J=K,K+JSTOP-JSTRT
          RHS(J)=RHS(J)*DL(IND)
          IND=IND+1
        END DO
      END DO
      CALL ZTRSM('Left','Lower','Transpose','Unit',N,1,
     *  ONE,AA,LDA,RHS,LDA)
      CALL ZMUMPS_333(N,PERM_ROW,RHS,WORK)
      CALL ZLASWP(1,RHS,LDA,1,lda,PERM,-1)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_33(LEN,ELG,X,INCX,Y,INCY)
      INTEGER LEN,INCX,INCY
      COMPLEX*16 ELG
      COMPLEX*16 X(*),Y(*)
      INTEGER I,J,K
      COMPLEX*16 TEMP
      I=1
      J=1
      DO K=1,LEN
        TEMP=-ELG*X(I)+Y(J)
        Y(J)=X(I)
        X(I)=TEMP
        I=I+INCX
        J=J+INCY
      END DO
      RETURN
      END SUBROUTINE ZMUMPS_33
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_30(LEN,ELG,X,INCX,Y,INCY)
      INTEGER LEN,INCX,INCY
      COMPLEX*16 ELG
      COMPLEX*16 X(*),Y(*)
      INTEGER I,J,K
      COMPLEX*16 TEMP
      I=1
      J=1
      DO K=1,LEN
        TEMP=ELG*Y(J)+X(I)
        X(I)=Y(J)
        Y(J)=TEMP
        I=I+INCX
        J=J+INCY
      END DO
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_29(LEN,ELG,X,INCX,Y,INCY)
      INTEGER LEN,INCX,INCY
      COMPLEX*16 ELG
      COMPLEX*16 X(*),Y(*)
      INTEGER I,J,K
      COMPLEX*16 TEMP
      I=1
      J=1
      DO K=1,LEN
        TEMP=ELG*X(I)+Y(J)
        Y(I)=TEMP
        I=I+INCX
        J=J+INCY
      END DO
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_28(LEN,ELG,X,INCX,Y,INCY)
      INTEGER LEN,INCX,INCY
      COMPLEX*16 ELG
      COMPLEX*16 X(*),Y(*)
      INTEGER I,J,K
      COMPLEX*16 TEMP
      I=1
      J=1
      DO K=1,LEN
        TEMP=ELG*Y(J)+X(I)
        X(I)=TEMP
        I=I+INCX
        J=J+INCY
      END DO
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_157(X,Y,ELG)
      COMPLEX*16 X,Y,ELG
      DOUBLE PRECISION ZERO,ONE
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      IF(X.EQ.ZERO) THEN
        WRITE(*,*) 'Zero multiplier'
        STOP
      ELSE
        ELG=-Y/X
      END IF
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_299(EUNIT,LUNIT,LIBNAM,
     *  FILNAM,STRING,ERR_NO,ERR_LEV,
     *  STRING1,STRING2,ERR_INF1)
      INTEGER EUNIT,LUNIT,ERR_NO,ERR_LEV,ERR_INF1
      CHARACTER*(*) LIBNAM,FILNAM,STRING,STRING1,STRING2
      IF(EUNIT.EQ.0) THEN
        WRITE(*,1) '       Error No. ',ERR_NO
        WRITE(*,3) ' -- ',STRING,STRING1
        WRITE(*,4) ' -- -- ',STRING2,ERR_INF1
        WRITE(*,2) '         found in routine  ',FILNAM
        WRITE(*,2) '         called within package  ',LIBNAM
      ELSE
        WRITE(EUNIT,1) '       Error No. ',ERR_NO
        WRITE(EUNIT,3) ' -- ',STRING,STRING1
        WRITE(EUNIT,4) ' -- -- ',STRING2,ERR_INF1
        WRITE(EUNIT,2) '         found in routine  ',FILNAM
        WRITE(EUNIT,2) '         called within package  ',LIBNAM
      END IF
      IF(LUNIT.NE.EUNIT) THEN
        WRITE(LUNIT,1) '       Error No. ',ERR_NO
        WRITE(LUNIT,3) ' -- ',STRING,STRING1
        WRITE(LUNIT,4) ' -- -- ',STRING2,ERR_INF1
        WRITE(LUNIT,2) '         found in routine  ',FILNAM
        WRITE(LUNIT,2) '         called within package  ',LIBNAM
      END IF
      IF(ERR_LEV.EQ.2) THEN
        STOP
      END IF
 1    FORMAT((A),I7)
 2    FORMAT((A),2X,(A))
 3    FORMAT(9X,(A),(A),1X,(A))
 4    FORMAT(9X,(A),(A),I7)
      RETURN
      END
#endif
      SUBROUTINE ZMUMPS_334(N,PERM,X,RN01)
      INTEGER N,PERM(N),I
      COMPLEX*16 RN01(N),X(N)
      DO 100 I=1,N
      RN01(I)=X(PERM(I))
100   CONTINUE
      DO 200 I=1,N
      X(I)=RN01(I)
200   CONTINUE
      RETURN
      END
#if defined(try_null_space_not_tested)
      FUNCTION ZMUMPS_332(N,X)
      INTEGER N
      COMPLEX*16 X(N),ZMUMPS_332
      COMPLEX*16 POM
      INTEGER I
      DOUBLE PRECISION ZERO,DEN
      PARAMETER (ZERO=0.0D0)
      DEN=ZERO
      DO 1 I=1,N
      DEN=MAX(DEN,ABS(X(I)))
    1 CONTINUE
      POM= DCMPLX(ZERO)
      IF (DEN.GT.ZERO) THEN
      DO 2 I=1,N
      POM=POM+(X(I)/DEN)**2
    2 CONTINUE
      ENDIF
      ZMUMPS_332=DEN*SQRT(POM)
      RETURN
      END
#endif
      DOUBLE PRECISION FUNCTION ZMUMPS_104(N,AA,WR01)
      INTEGER N
      COMPLEX*16 AA(N*N),WR01(N)
      INTEGER I,J,JSTRT
      COMPLEX*16 ZERO
      COMPLEX*16 VALUE
      PARAMETER(ZERO=0.0D0)
      VALUE=ZERO
      CALL ZMUMPS_117(N,ZERO,WR01,1)
      JSTRT=1
      DO I=1,N
        DO J=JSTRT,JSTRT+N-1
          WR01(I)=WR01(I)+ABS(AA(J))
        END DO
        JSTRT=JSTRT+N
      END DO
      DO I=1,N
        VALUE=MAX(ABS(VALUE),ABS(WR01(I)))
      END DO
      ZMUMPS_104=VALUE
      RETURN
      END
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_298(EUNIT,LUNIT,LIBNAM,FILNAM,
     *           STRING,ERR_NO,ERR_LEV)
      INTEGER EUNIT,LUNIT,ERR_NO,ERR_LEV
      CHARACTER*(*) LIBNAM,FILNAM,STRING
      IF(EUNIT.EQ.0) THEN
        WRITE(*,1) '       Error No. ',ERR_NO
        WRITE(*,3) ' -- ',STRING
        WRITE(*,2) '         found in routine  ',FILNAM
        WRITE(*,2) '         called within package  ',LIBNAM
      ELSE
        WRITE(EUNIT,1) '       Error No. ',ERR_NO
        WRITE(EUNIT,3) ' -- ',STRING
        WRITE(EUNIT,2) '         found in routine  ',FILNAM
        WRITE(EUNIT,2) '         called within package  ',LIBNAM
      END IF
      IF(LUNIT.NE.EUNIT) THEN
        WRITE(LUNIT,1) '       Error No. ',ERR_NO
        WRITE(LUNIT,3) ' -- ',STRING
        WRITE(LUNIT,2) '         found in routine  ',FILNAM
        WRITE(LUNIT,2) '         called within package  ',LIBNAM
      END IF
      IF(ERR_LEV.EQ.2) THEN
        STOP
      END IF
 1    FORMAT((A),I7)
 2    FORMAT((A),2X,(A))
 3    FORMAT(9X,(A),(A))
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_168(N,A,LDA,DEFICIENCY,MAXSTEPS,U,
     *  PERM_COL,ELROW,ELELG,DL,
     *  EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WORK,SMIN)
      INTEGER N,LDA,DEFICIENCY,MAXSTEPS
      COMPLEX*16 SMIN
      COMPLEX*16 A(LDA,*),U(*)
      INTEGER NLUPDATE,NUUPDATE,PTREL(*),PTREU(*)
      INTEGER PERM_COL(*)
      INTEGER ELROW(*),EUROW(*)
      COMPLEX*16 ELELG(*),EUELG(*),DL(*),WORK(*)
      INTEGER NSTEPS
      COMPLEX*16 UNORM
      DOUBLE PRECISION ONE
      PARAMETER(ONE=1.0D0)
      DOUBLE PRECISION DZNRM2
      NSTEPS=0
      CALL ZMUMPS_117(N-DEFICIENCY,ONE,U,1)
      CALL ZMUMPS_117(DEFICIENCY,0.0D0,U(N-DEFICIENCY+1),1)
 100  CONTINUE
      UNORM=DZNRM2(N-DEFICIENCY,U,1)
      CALL ZSCAL(N-DEFICIENCY,ONE/UNORM,U,1)
      CALL ZMUMPS_305(N,DEFICIENCY,A,N,
     *  DL,NLUPDATE,PTREL,ELROW,ELELG,
     *  NUUPDATE,PTREU,EUROW,EUELG,PERM_COL,
     *  U,WORK)
      CALL ZMUMPS_307(N,DEFICIENCY,A,N,
     *  DL,NLUPDATE,PTREL,ELROW,ELELG,
     *  NUUPDATE,PTREU,EUROW,EUELG,PERM_COL,
     *  U,WORK)
      NSTEPS=NSTEPS+1
      IF(NSTEPS.LT.MAXSTEPS) GO TO 100
      UNORM=DZNRM2(N,U,1)
      CALL ZSCAL(N,ONE/UNORM,U,1)
      SMIN=SQRT(ONE/UNORM)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      FUNCTION ZMUMPS_106(N,X)
      INTEGER N
      COMPLEX*16 X(N)
      DOUBLE PRECISION ZMUMPS_106
      INTEGER I
      DOUBLE PRECISION TEMP
      COMPLEX*16 XI
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      TEMP=ZERO
      DO I=1,N
        XI=X(I)
        TEMP=MAX(TEMP,ABS(XI))
      END DO
      ZMUMPS_106=TEMP
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_105(START,STOP,X)
      INTEGER START,STOP
      COMPLEX*16 X(*)
      INTEGER I
      COMPLEX*16 TEMP
      IF(STOP.LE.START) RETURN
      TEMP=X(START)
      DO I=START+1,STOP
        X(I-1)=X(I)
      END DO
      X(STOP)=TEMP
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_317(N,LDR,R,V,P,U,SMIN,VMIN)
      INTEGER N,LDR
      COMPLEX*16 R(LDR,*),V(*),P(*),U(*)
      COMPLEX*16 SMIN,VMIN(*)
      INTEGER I
      COMPLEX*16 V1,XI,ALPHA,BETA,GAMMA
      DOUBLE PRECISION C,S,ETA,MU
      DOUBLE PRECISION PNOR,DNOR,VNOR,PHI,UTP
      COMPLEX*16 LAMBDMAX,ROOT
      DOUBLE PRECISION ZERO,ONE
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION DDOT,DZNRM2,ZMUMPS_277
      V(1)= ONE/R(1,1)
      VNOR=ABS(V(1))
      V1=V(1)
      CALL ZCOPY(N,R(1,1),LDR,P,1)
      CALL ZSCAL(N,V1,P,1)
      DO I=2,N-1
        CALL ZCOPY(N-I,R(I,I+1),LDR,U(I+1),1)
        UTP=DDOT(N-I,U(I+1),1,P(I+1),1)
        GAMMA=R(I,I)
        XI=P(I)
        DNOR=DZNRM2(N-I,U(I+1),1)
        PHI=1+DNOR*DNOR
        PNOR=DZNRM2(N-I,P(I+1),1)
        ALPHA=XI*PHI-GAMMA*UTP
        IF(ALPHA.EQ.ZERO) THEN
          BETA=GAMMA*GAMMA*(VNOR*VNOR+PNOR*PNOR)-XI*XI*PHI
          IF(ABS(BETA).GT.ABS(PHI)) THEN
            S= ONE
            C= ZERO
          ELSE
            S= ZERO
            C= ONE
          END IF
        ELSE
          BETA=GAMMA*GAMMA*(VNOR*VNOR+PNOR*PNOR)+(XI*XI-1)*PHI
          BETA=BETA-2*XI*GAMMA*UTP
          ETA=BETA/(2.0D0*ALPHA)
          IF(ABS(ALPHA).GT.0) THEN
            MU=ETA+ZMUMPS_277(DCMPLX(ETA),DCMPLX(ONE))
          ELSE
            MU=ETA-ZMUMPS_277(DCMPLX(ETA),DCMPLX(ONE))
          END IF
          ROOT=ZMUMPS_277(DCMPLX(MU),DCMPLX(ONE))
          S=MU/ROOT
          C=-ONE/ROOT
        END IF
        CALL ZSCAL(I-1,S,V,1)
        V(I)=(C-S*XI)/GAMMA
        VNOR=ZMUMPS_277(DCMPLX(S*VNOR),(C-S*XI)/GAMMA)
        CALL ZMUMPS_116(N-I,S,P(I+1),1,V(I),U(I+1),1,P(I+1),1)
      END DO
      ALPHA=P(N)
      GAMMA=R(N,N)
      IF(ALPHA.EQ.0) THEN
        BETA=GAMMA*GAMMA*VNOR*VNOR
        IF(ABS(BETA).GT.ONE) THEN
          S= ONE
          C= ZERO
          LAMBDMAX=BETA
        ELSE
          S= ZERO
          C= ONE
          LAMBDMAX= DCMPLX(ONE)
        END IF
      ELSE
        BETA=GAMMA*GAMMA*VNOR*VNOR+ALPHA*ALPHA-ONE
        ETA=BETA/(2*ALPHA)
        IF(ABS(ALPHA).GT.0) THEN
          MU=ETA+ZMUMPS_277(DCMPLX(ETA),DCMPLX(ONE))
        ELSE
          MU=ETA-ZMUMPS_277(DCMPLX(ETA),DCMPLX(ONE))
        END IF
        LAMBDMAX=ALPHA*MU+ONE
        ROOT=ZMUMPS_277(DCMPLX(MU),DCMPLX(ONE))
        S=MU/ROOT
        C=-ONE/ROOT
      END IF
      CALL ZSCAL(N-1,S,V,1)
      V(N)=(C-S*ALPHA)/GAMMA
      VNOR=SQRT(LAMBDMAX)/ABS(GAMMA)
      CALL ZCOPY(N,V,1,VMIN,1)
      CALL ZSCAL(N,ONE/VNOR,VMIN,1)
      CALL ZTRSM('Left','Upper','No transpose','Non-unit',N,
     *  1,ONE,R,LDR,VMIN,N)
      SMIN= ONE/DZNRM2(N,VMIN,1)
      CALL ZSCAL(N,SMIN,VMIN,1)
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_305(N,DEFICIENCY,AA,LDA,
     *  DL,NLUPDATE,PTREL,ELROW,ELELG,
     *  NUUPDATE,PTREU,EUROW,EUELG,PERM_COL,
     *  RHS,WORK)
      INTEGER N,LDA,NLUPDATE,NUUPDATE,DEFICIENCY
      INTEGER ELROW(*),EUROW(*),PTREL(*),PTREU(*)
      INTEGER PERM_COL(*)
      COMPLEX*16 ELELG(*),EUELG(*),RHS(*),WORK(*),DL(*),AA(*)
      INTEGER II,IND,J,K,JSTRT,JSTOP
      COMPLEX*16 ELG
      COMPLEX*16 ONE
      PARAMETER(ONE=1.0D0)
      CALL ZTRSM('Left','Lower','No transpose','Unit',N-DEFICIENCY,1,
     *  ONE,AA,LDA,RHS,LDA)
      DO II=NLUPDATE,1,-1
        JSTRT=PTREL(II)
        JSTOP=PTREL(II+1)-1
        K=ELROW(JSTRT)
        IF(K.EQ.0) THEN
          K=ELROW(JSTOP)
        END IF
        IND=JSTRT
        DO J=K,K+JSTOP-JSTRT
          RHS(J)=RHS(J)*DL(IND)
          IND=IND+1
        END DO
        DO J=JSTOP,JSTRT,-1
          K=ELROW(J)
          IF(J.EQ.JSTOP) K=0
          IF(K.GT.0) THEN
            ELG=ELELG(J)
            CALL ZMUMPS_28(1,ELG,RHS(K),N,RHS(K+1),N)
          ELSEIF(K.LT.0) THEN
            ELG=ELELG(J)
            CALL ZMUMPS_30(1,ELG,RHS(-K),N,RHS(-K+1),N)
          ELSE
          END IF
        END DO
      END DO
      DO II=1,NUUPDATE
        JSTRT=PTREU(II)
        JSTOP=PTREU(II+1)-1
        DO J=JSTRT,JSTOP
          K=EUROW(J)
          ELG=EUELG(J)
          IF(K.GT.0) THEN
            CALL ZMUMPS_29(1,ELG,RHS(K),N,RHS(K+1),N)
          ELSE
            CALL ZMUMPS_30(1,ELG,RHS(-K),N,RHS(-K+1),N)
          END IF
        END DO
      END DO
      CALL ZTRSM('Left','Upper','No transpose','Non-unit',
     *  N-DEFICIENCY, 1,ONE,AA,LDA,RHS,LDA)
      CALL ZMUMPS_333(N,PERM_COL,RHS,WORK)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_307(N,DEFICIENCY,AA,LDA,
     *  DL,NLUPDATE,PTREL,ELROW,ELELG,
     *  NUUPDATE,PTREU,EUROW,EUELG,PERM_COL,
     *  RHS,WORK)
      INTEGER N,LDA,NLUPDATE,NUUPDATE,DEFICIENCY
      INTEGER ELROW(*),EUROW(*),PTREL(*),PTREU(*)
      INTEGER PERM_COL(*)
      COMPLEX*16 ELELG(*),EUELG(*),RHS(*),WORK(*),DL(*),AA(*)
      INTEGER II,IND,J,K,JSTRT,JSTOP
      COMPLEX*16 ELG
      COMPLEX*16 ONE
      PARAMETER(ONE=1.0D0)
      CALL ZMUMPS_334(N,PERM_COL,RHS,WORK)
      CALL ZTRSM('Left','Upper','Transpose','Non-unit',N-DEFICIENCY,
     *  1,ONE,AA,LDA,RHS,LDA)
      DO II=NUUPDATE,1,-1
        JSTRT=PTREU(II)
        JSTOP=PTREU(II+1)-1
        DO J=JSTOP,JSTRT,-1
          K=EUROW(J)
          ELG=EUELG(J)
          IF(K.GT.0) THEN
            CALL ZMUMPS_28(1,ELG,RHS(K),N,RHS(K+1),N)
          ELSEIF(K.LT.0) THEN
            CALL ZMUMPS_30(1,ELG,RHS(-K),N,RHS(-K+1),N)
          END IF
        END DO
      END DO
      DO II=1,NLUPDATE
        JSTRT=PTREL(II)
        JSTOP=PTREL(II+1)-1
        DO J=JSTRT,JSTOP
          K=ELROW(J)
          IF(J.EQ.JSTOP) K=0
          IF(K.GT.0) THEN
            ELG=ELELG(J)
            CALL ZMUMPS_29(1,ELG,RHS(K),N,RHS(K+1),N)
          ELSEIF(K.LT.0) THEN
            ELG=ELELG(J)
            CALL ZMUMPS_30(1,ELG,RHS(-K),N,RHS(-K+1),N)
          END IF
        END DO
        K=ELROW(JSTRT)
        IF(K.EQ.0) THEN
          K=ELROW(JSTOP)
        END IF
        IND=JSTRT
        DO J=K,K+JSTOP-JSTRT
          RHS(J)=RHS(J)*DL(IND)
          IND=IND+1
        END DO
      END DO
      CALL ZTRSM('Left','Lower','Transpose','Unit',N-DEFICIENCY,1,
     *  ONE,AA,LDA,RHS,LDA)
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_116(N,DA,DX,INCX,DB,DY,INCY,DZ,INCZ)
      INTEGER N,INCX,INCY,INCZ
      COMPLEX*16 DA,DB,DX(*),DY(*),DZ(*)
      INTEGER IX,IY,IZ,I,M,MP1
      DOUBLE PRECISION ZERO
      PARAMETER(ZERO=0.0D0)
      IF(N.LE.0) RETURN
      IF(DA.EQ.ZERO) THEN
        CALL ZCOPY(N,DY,1,DZ,1)
        CALL ZSCAL(N,DB,DZ,1)
      END IF
      IF(.NOT.(INCX.EQ.1.AND.INCY.EQ.1.AND.INCZ.EQ.1)) THEN
        IX=1
        IY=1
        IZ=1
        IF(INCX.LT.0) IX=(-N+1)*INCX+1
        IF(INCY.LT.0) IY=(-N+1)*INCY+1
        IF(INCZ.LT.0) IZ=(-N+1)*INCZ+1
        DO I=1,N
          DZ(IZ)=DB*DY(IY)+DA*DX(IX)
          IX=IX+INCX
          IY=IY+INCY
          IZ=IZ+INCZ
        END DO
      ELSE
        M=MOD(N,4)
        IF(M.NE.0) THEN
          DO I=1,M
            DZ(I)=DB*DY(I)+DA*DX(I)
          END DO
          IF(N.LT.4) RETURN
        END IF
        MP1=M+1
        DO I=MP1,N,4
          DZ(I)=DB*DY(I)+DA*DX(I)
          DZ(I+1)=DB*DY(I+1)+DA*DX(I+1)
          DZ(I+2)=DB*DY(I+2)+DA*DX(I+2)
          DZ(I+3)=DB*DY(I+3)+DA*DX(I+3)
        END DO
      END IF
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_182(N,A,LDA,DEFICIENCY,U,V,
     *  PERM_COL,ELROW,ELELG,DL,
     *  EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WORK,SMIN)
      INTEGER N,LDA,DEFICIENCY
      COMPLEX*16 SMIN
      COMPLEX*16 A(LDA,*),U(*),V(*)
      INTEGER NLUPDATE,NUUPDATE,PTREL(*),PTREU(*)
      INTEGER PERM_COL(*)
      INTEGER ELROW(*),EUROW(*)
      COMPLEX*16 ELELG(*),EUELG(*),DL(*),WORK(*)
      COMPLEX*16 VNORM
      DOUBLE PRECISION ZERO,ONE
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION DZNRM2
      CALL ZCOPY(N,U,1,V,1)
      CALL ZMUMPS_305(N,DEFICIENCY,A,N,
     *  DL,NLUPDATE,PTREL,ELROW,ELELG,
     *  NUUPDATE,PTREU,EUROW,EUELG,PERM_COL,
     *  V,WORK)
      VNORM=DZNRM2(N,V,1)
      CALL ZSCAL(N,ONE/VNORM,V,1)
      SMIN= ONE/VNORM
      RETURN
      END
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_262( SIDE, TRANS, M, N, K, A, IA, JA,
     *                    DESCA, TAU,
     *                    C, IC, JC, DESCC, WORK, LWORK, INFO )
      CHARACTER          SIDE, TRANS
      INTEGER            IA, IC, INFO(2), JA, JC, K, LWORK, M, N
      INTEGER            DESCA( * ), DESCC( * )
      COMPLEX*16   A( * ), C( * ), TAU( * ), WORK( * )
      INTEGER MINUS_ONE
      PARAMETER( MINUS_ONE = -1 )
      INTEGER LREQ, allocok
      COMPLEX*16, ALLOCATABLE, DIMENSION (:) :: WORKSPACE
      CALL PZUNMQR( SIDE, TRANS, M, N, K, A, IA, JA,
     *                    DESCA, TAU,
     *                    C, IC, JC, DESCC, WORK, MINUS_ONE, INFO(1) )
      LREQ = INT( WORK(1) )
      IF ( LWORK .LT. LREQ ) THEN
        ALLOCATE( WORKSPACE( LREQ ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = LREQ
        ELSE
          CALL PZUNMQR( SIDE, TRANS, M, N, K, A, IA, JA,
     *                    DESCA, TAU,
     *                    C, IC, JC, DESCC, WORKSPACE, LREQ, INFO(1) )
          DEALLOCATE( WORKSPACE )
        END IF
      ELSE
        CALL PZUNMQR( SIDE, TRANS, M, N, K, A, IA, JA,
     *                    DESCA, TAU,
     *                    C, IC, JC, DESCC, WORK, LWORK, INFO(1))
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_262
#endif
      SUBROUTINE ZMUMPS_259(V,N)
      INTEGER N
      COMPLEX*16 V( N )
      INTEGER I
      DOUBLE PRECISION TEMP
      TEMP = 0.0D0
      DO I = 1, N
        IF ( ABS( V( I ) ) .GT. TEMP ) TEMP = ABS( V( I ) )
      END DO
      TEMP = 1.0D0 / TEMP
      DO I = 1, N
        V( I ) = TEMP * V( I )
      END DO
      RETURN
      END SUBROUTINE ZMUMPS_259
      SUBROUTINE ZMUMPS_247(null_basis_col,root,
     * N, A, LA, IW, LIW, W, LWC,
     * IWCB,LIWW,NULL,NULLSOL,NA, LNA, NE_STEPS, W2, MTYPE, ICNTL,
     * STEP, FRERE,FILS, PTRIST, PTRFAC, IW1,LIW1,
     * PROCNODE_STEPS, SLAVEF,
     * INFO, KEEP,KEEP8, COMM, COMM_NODES, MYID,
     * MYID_NODES, BUFR, LBUFR, LBUFR_BYTES, 
     *
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      TYPE ( ZMUMPS_ROOT_STRUC ) :: root
      INTEGER LA,LWC,N,LIW,MTYPE,LIW1,LIWW,null_basis_col,LNA
      INTEGER ICNTL(40),INFO(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IW(LIW),IW1(LIW1),NA(LNA),NE_STEPS(KEEP(28)),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)), FILS(N), PTRIST(KEEP(28))
      INTEGER PTRFAC(KEEP(28))
      COMPLEX*16    A(LA), W(LWC), NULL(N),
     *        NULLSOL(*), W2(N)
      INTEGER SLAVEF, COMM, COMM_NODES, MYID, MYID_NODES
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER XSIZE
      PARAMETER(XSIZE=3)
      INTEGER MP, LDIAG, LWC_LOC
      INTEGER K,I,NN
      INTEGER LPOOL,MYLEAF
      INTEGER IPOOL,PTRICB,PTRACB
      INTEGER POOLSS,QR
      INTEGER SIZE_ROOT, IPT_NULL_ROOT
      INTEGER IERR,allocok
      INTEGER       IOLDPS, IAPOS,
     *              LOCAL_M     ,
     *              LOCAL_N     ,
     *              MASTER_ROOT
      INTEGER IROOTPOS
      COMPLEX*16 ZERO
      PARAMETER(ZERO=0.0D0)
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: NULL_PAR
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV_SEQ
      INCLUDE 'mpif.h'
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
      CALL ZMUMPS_158
     *              ( N, MYID_NODES, STEP, PROCNODE_STEPS,
     *                KEEP(28), SLAVEF,
     *                NA, LNA, MYLEAF )
      MP      = ICNTL(2)
      LDIAG   = ICNTL(4)
      POOLSS = 1
      PTRICB = POOLSS + KEEP(28)+1
      PTRACB = PTRICB + KEEP(28)
      IPOOL  = PTRACB + KEEP(28)
      LPOOL  = LIW1 - IPOOL + 1
      IF ( LPOOL .ne. KEEP(28) + 1 ) THEN
        WRITE(MP,*) ' warning solve: lpool.NE.n+1 :',lpool,n+1
      ENDIF
      LWC_LOC = LWC
      SIZE_ROOT  = 0
      IPT_NULL_ROOT  = LWC
      IF ( KEEP( 38 ) .NE. 0 .OR. KEEP(20) .NE. 0 ) THEN
          IF ( KEEP(38) .NE. 0 ) THEN
            MASTER_ROOT = ZMUMPS_275( STEP(KEEP(38)),
     *                    PROCNODE_STEPS, SLAVEF )
          ELSE
            MASTER_ROOT = ZMUMPS_275( STEP(KEEP(20)),
     *                    PROCNODE_STEPS, SLAVEF )
          END IF
          IF ( MASTER_ROOT .eq. MYID_NODES ) THEN
            IF ( KEEP(38) .NE. 0 ) THEN
              SIZE_ROOT    = root%TOT_ROOT_SIZE
            ELSE
              SIZE_ROOT    = IW( PTRIST(STEP( KEEP(20) ) +XSIZE))
            END IF
            IPT_NULL_ROOT = LWC - SIZE_ROOT + 1
            LWC_LOC = LWC_LOC - SIZE_ROOT
          END IF
      END IF
      DO i=1, N
         NULLSOL(i)= ZERO
         NULL   (i)= ZERO
      END DO
#ifdef try_null_space_not_tested
      IF ( KEEP( 38 ) .NE. 0 .and. root%yes ) THEN
        IOLDPS = PTRIST(STEP( KEEP(38) ))
        LOCAL_M = IW( IOLDPS + 2 +XSIZE )
        LOCAL_N = IW( IOLDPS + 1 +XSIZE )
        IAPOS   = PTRFAC(IW( IOLDPS + 4 + XSIZE ))
        CALL ZCOPY(root%TOT_ROOT_SIZE,NULL,1,W(IPT_NULL_ROOT),1)
      END IF
#endif
#ifdef try_null_space_not_tested
      IF ( KEEP(38) .NE. 0 ) THEN
        IF ( MYID_NODES .eq. MASTER_ROOT ) THEN
          ALLOCATE( IPIV_SEQ(root%TOT_ROOT_SIZE), stat=allocok )
          IF ( allocok.GT.0 ) THEN
            INFO(1) = -13
            INFO(2) = SIZE_ROOT
            STOP
          ENDIF
        ELSE
          ALLOCATE( IPIV_SEQ( 1 ), stat=allocok )
          IF ( allocok.GT.0 ) THEN
            INFO(1)  = -13
            INFO(2) = 1
            STOP
          ENDIF
        END IF
        ALLOCATE( NULL_PAR(LOCAL_M), stat=allocok )
        IF ( allocok.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = SIZE_ROOT
          STOP
        ENDIF
        do i=1,local_m
           null_par(i)=0.0d0
        end do
       QR=KEEP(19)
       IF ( QR .eq. 9 .or. QR.eq.10) THEN
         QR = QR - 8
       END IF
       IF(QR.EQ.2.OR.QR.EQ.4) THEN
        CALL ZMUMPS_114(root%TOT_ROOT_SIZE,root%TOT_ROOT_SIZE,
     *     LOCAL_M,LOCAL_N,root%MBLOCK,root%NBLOCK,A(IAPOS),
     *     W(IPT_NULL_ROOT),root%TOT_ROOT_SIZE-NULL_BASIS_COL+1,
     *     MYID_NODES,MASTER_ROOT,root%NPROW,root%NPCOL,COMM_NODES)
       CALL ZMUMPS_173(root%TOT_ROOT_SIZE,1,
     *   LOCAL_M,1,root%MBLOCK,root%NBLOCK,root%IPIV,IPIV_SEQ,
     *   MYID_NODES,MASTER_ROOT,root%NPROW,root%NPCOL,COMM_NODES)
          if(myid_nodes.eq.master_root) then
          do i=root%TOT_ROOT_SIZE,root%TOT_ROOT_SIZE-keep(17)+1,-1
            w(ipt_null_root+i-1)=0.0d0
          end do
          w(ipt_null_root+root%tot_root_size-null_basis_col)=-1.0d0
           CALL ZMUMPS_333(root%TOT_ROOT_SIZE,IPIV_SEQ,
     *                       W(IPT_NULL_ROOT),W2 )
          end if
         ELSEIF(QR.EQ.6.OR.QR.EQ.8) THEN
           CALL ZMUMPS_114(root%TOT_ROOT_SIZE,root%TOT_ROOT_SIZE,
     *       LOCAL_M,LOCAL_N,root%MBLOCK,root%NBLOCK,A(IAPOS),
     *       W(IPT_NULL_ROOT),root%TOT_ROOT_SIZE-NULL_BASIS_COL+1,
     *       MYID_NODES,MASTER_ROOT,root%NPROW,root%NPCOL,COMM_NODES)
          if(myid_nodes.eq.master_root) then
          do i=root%TOT_ROOT_SIZE,root%TOT_ROOT_SIZE-KEEP(17)+1,-1
            w(ipt_null_root+i-1)=0.0d0
          end do
          w(ipt_null_root+root%tot_root_size-null_basis_col)=-1.0d0
          end if
         END IF
        CALL ZMUMPS_249(N, A, LA, IW, LIW, W(1), LWC_LOC, NULL,
     *              N, 1, NULLSOL,
     *              IW1(PTRICB),IW1(PTRACB),IWCB,LIWW,
     *              IW1(POOLSS),W2,
     *              NE_STEPS, NA, LNA, STEP, FRERE,FILS,
     *              IW1(IPOOL),LPOOL,PTRIST,PTRFAC,MYLEAF,INFO,
     *              PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     *              BUFR, LBUFR, LBUFR_BYTES, KEEP,KEEP8,
     *              W( IPT_NULL_ROOT ),
     *              MTYPE, 
     *              ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        IF (INFO(1).eq.-2) INFO(1)=-11
        IF (INFO(1).eq.-3) INFO(1)=-14
        deallocate(IPIV_SEQ)
        deallocate(null_par)
      ELSE
#endif
        IF ( MYID_NODES .eq. MASTER_ROOT ) THEN
          IROOTPOS = IW( PTRIST(STEP( KEEP(20) )) + 4 + XSIZE )
        ELSE
          IROOTPOS = 1
        END IF
        IF ( MYID_NODES .eq. MASTER_ROOT ) THEN
          QR=KEEP(19)
          IF (QR.eq.9 .or. QR.eq.10) QR = QR - 8
#if defined(try_null_space_not_tested)
          IF(QR.EQ.2.OR.QR.EQ.4) THEN
#endif
            CALL ZMUMPS_85(SIZE_ROOT,A(IROOTPOS),SIZE_ROOT,
     *      W(IPT_NULL_ROOT),null_basis_col,KEEP(17),root%QR_TAU,
     *      root%IPIV,W(1), LWC_LOC )
#if defined(try_null_space_not_tested)
          ELSE
            CALL ZMUMPS_82(SIZE_ROOT,A(IROOTPOS),SIZE_ROOT,
     *      W(IPT_NULL_ROOT),null_basis_col,KEEP(17),
     *      root%IPIV)
          END IF
#endif
        END IF
        CALL ZMUMPS_249(N, A, LA, IW, LIW, W(1), LWC_LOC, NULL,
     *              N, 1, NULLSOL,
     *              IW1(PTRICB),IW1(PTRACB),IWCB,LIWW,
     *              IW1(POOLSS),W2,
     *              NE_STEPS, NA, LNA, STEP, FRERE,FILS,
     *              IW1(IPOOL),LPOOL,PTRIST,PTRFAC,MYLEAF,INFO,
     *              PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     *              BUFR, LBUFR, LBUFR_BYTES, KEEP,KEEP8,
     *              W( IPT_NULL_ROOT ),
     *              MTYPE, 
     *              ISTEP_TO_INIV2, TAB_POS_IN_PERE)
#ifdef try_null_space_not_tested
      END IF
#endif
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
        K = MIN0(10,N)
        IF (LDIAG.EQ.4) K = N
        WRITE (MP,99992)
        IF (N.GT.0) WRITE (MP,99993) (NULLSOL(I),I=1,K)
      ENDIF
      RETURN
99993 FORMAT (' NULSOL'/(1X,1P,5D14.6))
99992 FORMAT (//' LEAVING SOLVE (MPI41C) WITH')
      END SUBROUTINE ZMUMPS_247
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_31(LEN,CS,SN,X,INCX,Y,INCY)
      INTEGER LEN,INCX,INCY
      DOUBLE PRECISION CS,SN
      DOUBLE PRECISION X(*),Y(*)
      INTEGER I,J,K
      DOUBLE PRECISION TEMP
      I=1
      J=1
      DO K=1,LEN
        TEMP=CS*X(I)+SN*Y(J)
        Y(J)=-SN*X(I)+CS*Y(J)
        X(I)=TEMP
        I=I+INCX
        J=J+INCY
      END DO
      RETURN
      END
#endif
      SUBROUTINE ZMUMPS_117(N,ALPHA,DX,INCX)
      INTEGER N,INCX
      COMPLEX*16 ALPHA,DX(*)
      INTEGER I,IX,M,MP1
      INTRINSIC MOD
      IF(N.LE.0) RETURN
      IF(INCX.NE.1) THEN
        IX=1
        IF(INCX.LT.0) IX=(-N+1)*INCX+1
        DO I=1,N
          DX(IX)=ALPHA
          IX=IX+INCX
        END DO
        RETURN
      ELSE
        M=MOD(N,7)
        IF(M.NE.0) THEN
          DO I=1,M
            DX(I)=ALPHA
          END DO
          IF(N.LT.7) RETURN
        END IF
        MP1=M+1
        DO I=MP1,N,7
          DX(I)=ALPHA
          DX(I+1)=ALPHA
          DX(I+2)=ALPHA
          DX(I+3)=ALPHA
          DX(I+4)=ALPHA
          DX(I+5)=ALPHA
          DX(I+6)=ALPHA
        END DO
        RETURN
      END IF
      END
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_175(N,IALPHA,IX,INCX)
      INTEGER N,IALPHA,INCX
      INTEGER IX(*)
      INTEGER I,IND,M,MP1
      INTRINSIC MOD
      IF(N.LE.0) RETURN
      IF(INCX.NE.1) THEN
        IND=1
        IF(INCX.LT.0) IND=(-N+1)*INCX+1
        DO I=1,N
          IX(IND)=IALPHA
          IND=IND+INCX
        END DO
        RETURN
      ELSE
        M=MOD(N,7)
        IF(M.NE.0) THEN
          DO I=1,M
            IX(I)=IALPHA
          END DO
          IF(N.LT.7) RETURN
        END IF
        MP1=M+1
        DO I=MP1,N,7
          IX(I)=IALPHA
          IX(I+1)=IALPHA
          IX(I+2)=IALPHA
          IX(I+3)=IALPHA
          IX(I+4)=IALPHA
          IX(I+5)=IALPHA
          IX(I+6)=IALPHA
        END DO
        RETURN
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_175
#endif
#if defined(try_null_space_not_tested)
      SUBROUTINE ZMUMPS_147( MYID, root, N, IROOT,
     *           COMM, IW, LIW, IFREE,
     *           A, LA, PTRAST, PTLUST_S, PTRFAC,
     *           STEP, INFO, LDLT, QR,
     *           WK, LWK, KEEP,KEEP8, WK1,LWK1,RV,VAUX_SEQ )
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      INCLUDE 'mpif.h'
      TYPE ( ZMUMPS_ROOT_STRUC ) :: root
      INTEGER N, IROOT, COMM, LIW, LA, MYID, LWK, LIWK, IFREE,LWK1
      COMPLEX*16 WK( LWK ),WK1(LWK1,1)
      COMPLEX*16 RV(LWK1,1),VAUX_SEQ(LWK1)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRAST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *STEP(N), IW( LIW )
      INTEGER INFO( 2 ), LDLT, QR
      COMPLEX*16 A( LA )
      INTEGER IOLDPS, IAPOS
      INTEGER XSIZE
      PARAMETER(XSIZE=3)
      INTEGER LOCAL_M, LOCAL_N, LPIV, IERR, allocok,i
      INTEGER JJ,ROW,NR,IDIST
      integer deficiency,rsrc,csrc,lrindx,lcindx
      integer nprow,npcol,myrow,mycol,rowproc
      integer iarow,iacol,MBLOCK,NBLOCK,csrca,rsrca
      integer indxg2p,master_root,DIM,J,K
      INTEGER POS,ISOUR,RSRCP,CSRCP,DEF,IVROW,IVCOL
      INTEGER COLPROC,PRO2,LOC2,IDEST
      INTEGER RSRCV,CSRCV,DESCV(9),DESCU(9)
      INTEGER MAXG,GIND,ZERO_DIAG
      INTEGER STATUS(MPI_STATUS_SIZE)
      INTEGER IMAX,GMAX,RMAX,PRO,LOC,ENVELOPE,ADD_COLUMN
      INTEGER CURRENT_DEFICIENCY
      COMPLEX*16 VMAX
      COMPLEX*16 SMIN,RROOT,ETA,C,S,BETA,PHI,XJ
      COMPLEX*16 LAMBDMAX,VNEW,PYTHAG,PNOR,VNOR,GAMMA
      COMPLEX*16 TEMP, ONE
      DOUBLE PRECISION ZERO, SV_MU, EPS
      DOUBLE PRECISION RANK_REVEALING_TOL_1
      DOUBLE PRECISION RANK_REVEALING_TOL_2
      DOUBLE PRECISION TOL,MU,TOL1
      COMPLEX*16 ALPHA,UTP,ANORM,RCOND,PDLANGE
      COMPLEX*16 ZMUMPS_104
      INTEGER SMP_DEFICIENCY,NN
      INTEGER ELIND,EUIND,NLUPDATE,NUUPDATE
      INTEGER NULL_MAX,ELG_MAX,EST_DEFICIENCY
      INTEGER RSRCU,CSRCU,LDA,MAXDEF
      LOGICAL NOT_YET
      PARAMETER(NOT_YET=.FALSE.)
      INTEGER IND,NULL_IND,JND,NUM_RANK
      INTEGER IND1,IND2,KND,USHIFT
      INTEGER KK,LUSHIFT,LSHIFT,IND_ROT
      INTEGER DIML,DIMU,DIMLU,ELIND_SAVE,ROT_STEPS
      INTEGER L,II,IND_TO,IND_FROM,NSHORT
      COMPLEX*16 TEMP1,TEMP2,ELG
      DOUBLE PRECISION ABSYL,ABSXK,YMAX,XMAX
      DOUBLE PRECISION ABSVK
      COMPLEX*16 ZMUMPS_106
      DOUBLE PRECISION EPS10,EPS12,EPS14,EPS9
      PARAMETER(EPS10=1.0D-10,EPS12=1.0D-12,EPS14=1.0D-14)
      PARAMETER(EPS9=1.0D-9)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      PARAMETER(NULL_MAX=20,ELG_MAX=10000)
      PARAMETER(SV_MU=0.1D0,EPS=1.1D-16)
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VMIN
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV_SEQ
      INTEGER, DIMENSION(:), ALLOCATABLE :: PERM_ROW
      INTEGER, DIMENSION(:), ALLOCATABLE :: PERM_COL
      INTEGER, DIMENSION(:), ALLOCATABLE :: PERM_INIT
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: ELELG
      INTEGER, DIMENSION(:), ALLOCATABLE :: ELROW
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: EUELG
      INTEGER, DIMENSION(:), ALLOCATABLE :: EUROW
      INTEGER, DIMENSION(:), ALLOCATABLE :: PTREL
      INTEGER, DIMENSION(:), ALLOCATABLE :: PTREU
      INTEGER, DIMENSION(:), ALLOCATABLE :: WN01
      INTEGER, DIMENSION(:), ALLOCATABLE :: LROT
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: AR
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: RIGHT_SV
      COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: AL
      COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: AU
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: DL
      COMPLEX*16, DIMENSION(:), POINTER :: GROW, GCOS, GSIN
      IF ( root%yes ) THEN
        IOLDPS  = PTLUST_S(STEP(IROOT)) + XSIZE
        IAPOS   = PTRAST(STEP(IROOT))
        LOCAL_M = IW( IOLDPS + 2 )
        LOCAL_N = IW( IOLDPS + 1 )
        IAPOS = PTRFAC(IW ( IOLDPS + 4 ))
        IF ( LDLT.EQ.0 .OR. LDLT.EQ.2 .OR. QR.ne.0 ) THEN
         LPIV = LOCAL_M + root%MBLOCK
        ELSE
         LPIV = 1
        END IF
        IF (ASSOCIATED( root%IPIV )) DEALLOCATE(root%IPIV)
        root%LPIV = LPIV
        ALLOCATE( root%IPIV( LPIV ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = LPIV
          WRITE(*,*) MYID,': problem allocating IPIV(',LPIV,') in root'
          CALL ZMUMPS_ABORT()
        END IF
        CALL DESCINIT( root%DESCRIPTOR, root%TOT_ROOT_SIZE,
     *      root%TOT_ROOT_SIZE, root%MBLOCK, root%NBLOCK,
     *      0, 0, root%CNTXT_BLACS, LOCAL_M, IERR )
        CALL DESCINIT( root%DESCB, root%TOT_ROOT_SIZE,
     *      1, root%MBLOCK, root%NBLOCK, 0, 0,
     *      root%CNTXT_BLACS, LOCAL_M, IERR )
        IF (IERR.NE.0) THEN
          WRITE(*,*) 'After DESCINIT, IERR = ', IERR
          CALL ZMUMPS_ABORT()
        END IF
        IF(LDLT.EQ.2.OR.(LDLT.NE.0.AND.QR.NE.0)) THEN
          IF(root%MBLOCK.NE.root%NBLOCK) THEN
            WRITE(*,*) ' Error: symmetrization only works for'
            WRITE(*,*) ' square block sizes, MBLOCK/NBLOCK=',
     *      root%MBLOCK, root%NBLOCK
            CALL ZMUMPS_ABORT()
          END IF
          IF(root%yes) THEN
            IF ( LWK .LT.
     *        MIN(root%MBLOCK * root%NBLOCK,root%TOT_ROOT_SIZE*
     *        root%TOT_ROOT_SIZE ) ) THEN
              WRITE(*,*) 'Not enough workspace for symmetrization.'
              CALL ZMUMPS_ABORT()
            END IF
            CALL ZMUMPS_320( WK, root%MBLOCK,
     *      root%MYROW, root%MYCOL, root%NPROW, root%NPCOL,
     *      A( IAPOS ), LOCAL_M, LOCAL_N,
     *      root%TOT_ROOT_SIZE, MYID, COMM )
          END IF
        END IF
        DEFICIENCY=0
        KEEP(17)=DEFICIENCY
        TOL=EPS12
        RANK_REVEALING_TOL_1=EPS12
        NN=root%TOT_ROOT_SIZE
        MBLOCK=root%DESCRIPTOR(5)
        NBLOCK=root%DESCRIPTOR(6)
        rsrca=root%DESCRIPTOR(7)
        csrca=root%DESCRIPTOR(8)
        LDA=root%DESCRIPTOR(9)
        MAXDEF=KEEP(21)
        IF(MAXDEF.LE.0) MAXDEF=NN
        GIND=0
        root%GIND=GIND
        CALL BLACS_GRIDINFO(root%DESCRIPTOR(2),NPROW,NPCOL,
     *    MYROW,MYCOL)
        IF(QR.EQ.7) QR=5
        IF(QR.EQ.8) QR=6
        IF(QR.EQ.1.OR.QR.EQ.2.OR.QR.EQ.3.OR.QR.EQ.4) THEN
          GIND=0
          root%GIND=GIND
          IF(ASSOCIATED( root%QR_TAU ) ) DEALLOCATE( root%QR_TAU )
          ALLOCATE( root%QR_TAU( root%TOT_ROOT_SIZE ),
     *      stat=allocok )
          IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=LPIV
            WRITE(*,*) MYID,': problem allocating Householder coeffs'
            CALL ZMUMPS_ABORT()
          END IF
          IF(QR.EQ.1.OR.QR.EQ.2) then
            IF(root%MBLOCK.GE.NN.AND.NOT_YET) THEN
              IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
                 CALL ZGEQPF(NN,NN,A(IAPOS),LDA,root%IPIV(1),
     *             root%QR_TAU,WK,INFO)
                IF ( IERR .NE. 0 ) THEN
                  WRITE(*,*) ' Problem in ZGEQPF : IERR = ', IERR
                  CALL ZMUMPS_ABORT()
                END IF
 1488           CONTINUE
                IND=IAPOS-1+(NN-DEFICIENCY)*LDA-DEFICIENCY
                ALPHA=A(IND)
                IF(abs(alpha).le.RANK_REVEALING_TOL_1) THEN
                  deficiency=deficiency+1
                  go to 1488
                END IF
                IF(DEFICIENCY.GT.0) THEN
                  CALL ZTRTRS('U','N','N',
     *              NN-DEFICIENCY,DEFICIENCY,
     *              A(iapos),LDA,a(iapos+LDA*(NN-DEFICIENCY)),
     *              LDA,IERR)
                  IF ( IERR .NE. 0 ) THEN
                    WRITE(*,*) ' Problem in ZTRTRS : IERR = ', IERR
                    CALL ZMUMPS_ABORT()
                  END IF
                END IF
              END IF
            ELSE
            CALL PZGEQPF(NN,NN,A(IAPOS),
     *      1, 1, root%DESCRIPTOR, root%IPIV(1),
     *      root%QR_TAU, WK, LWK, IERR )
            IF ( IERR .NE. 0 ) THEN
              WRITE(*,*) ' Problem in PZGEQPF : IERR = ', IERR
              CALL ZMUMPS_ABORT()
            END IF
 148        CONTINUE
            call infog2l(root%TOT_ROOT_SIZE-deficiency,
     *        root%TOT_ROOT_SIZE-deficiency,
     *        root%DESCRIPTOR,nprow,npcol,
     *        myrow,mycol,lrindx,lcindx,rsrc,csrc)
            if(myrow.EQ.rsrc.and.mycol.EQ.csrc) then
              alpha=a(local_m*(lcindx-1)+lrindx+iapos-1)
            END IF
            ISOUR=RSRC*NPCOL+CSRC
            CALL MPI_BCAST(ALPHA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *        COMM,IERR)
            IF(abs(alpha).le.RANK_REVEALING_TOL_1) THEN
              deficiency=deficiency+1
              IF(DEFICIENCY.LE.MAXDEF) GO TO 148
            END IF
            call MPI_BCAST(deficiency,1,MPI_INTEGER,rsrc*nprow+csrc,
     *        COMM,IERR)
            IF(DEFICIENCY.GT.0) THEN
              CALL PZTRTRS('U','N','N',
     *          root%TOT_ROOT_SIZE-DEFICIENCY,DEFICIENCY,
     *          A(iapos),1,1,root%DESCRIPTOR,a(iapos),1,
     *          root%TOT_ROOT_SIZE-deficiency+1,root%DESCRIPTOR,IERR)
            END IF
            END IF
          ELSEIF(QR.EQ.3.OR.QR.EQ.4) THEN
            IF(root%MBLOCK.GE.NN.AND.NOT_YET) THEN
              IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
                CALL ZGEQPF(NN,NN,A(IAPOS),LDA,root%IPIV(1),
     *            root%QR_TAU,WK,INFO)
                IF(IERR.NE.0) THEN
                  WRITE(*,*) ' Problem in ZGEQPF : IERR = ', IERR
                  CALL ZMUMPS_ABORT()
                END IF
              root%MAXG = 5000
              ALLOCATE(root%GROW(root%MAXG),stat=allocok)
              IF(allocok.GT.0) THEN
                INFO(1) = -13
                INFO(2) = N
                CALL ZMUMPS_ABORT()
              ENDIF
              ALLOCATE(root%GCOS(root%MAXG),stat=allocok)
              IF(allocok.GT.0) THEN
                INFO(1) = -13
                INFO(2) = N
                CALL ZMUMPS_ABORT()
              ENDIF
              ALLOCATE(root%GSIN(root%MAXG),stat=allocok)
              IF(allocok.GT.0) THEN
                INFO(1) = -13
                INFO(2) = N
                CALL ZMUMPS_ABORT()
              ENDIF
              MAXG = root%MAXG
              GROW => root%GROW
              GCOS => root%GCOS
              GSIN => root%GSIN
              ALLOCATE( IPIV_SEQ(root%TOT_ROOT_SIZE), stat=allocok )
              IF(allocok.GT.0) THEN
                INFO(1) = -13
                INFO(2) = N
                CALL ZMUMPS_ABORT()
              ENDIF
            CALL ZMUMPS_117(LOCAL_M+1,ZERO,WK1,1)
            CALL ZMUMPS_117(LOCAL_M+1,ZERO,WK1(LOCAL_M+2,1),1)
            ALLOCATE( AR(NN*NN), stat=allocok )
            IF ( allocok.GT.0 ) THEN
              INFO(1) = -13
              INFO(2) = N
              CALL ZMUMPS_ABORT()
            ENDIF
            CALL ZCOPY(NN*NN,A(IAPOS),1,AR,1)
            GIND=0
            CALL ZCOPY(NN,root%IPIV,1,IPIV_SEQ,1)
          TOL=SQRT(DBLE(NN))*EPS*ZMUMPS_104(NN,A(IAPOS),VAUX_SEQ)
              TOL=EPS9
              TOL1=EPS9
              DEFICIENCY=0
              CALL ZMUMPS_98(NN,A(IAPOS),NN,SMP_DEFICIENCY,TOL1)
            DEF=SMP_DEFICIENCY
            DIM=root%TOT_ROOT_SIZE-DEF+1
            NR=root%TOT_ROOT_SIZE
      NSHORT=NN-DEF
      IND_FROM=LDA*(NN-DEF)+1
      IND_TO=1
      IND=0
      DO I=1,SMP_DEFICIENCY
        CALL ZCOPY(NSHORT,A(IAPOS+IND_FROM-1),
     *    1,VAUX_SEQ(IND_TO),1)
        CALL ZSCAL(NSHORT,ONE,VAUX_SEQ(IND_TO),1)
        IND_TO=IND_TO+NSHORT
        CALL ZMUMPS_117(SMP_DEFICIENCY,ZERO,VAUX_SEQ(IND_TO),1)
        VAUX_SEQ(IND_TO+IND)= ONE
        IND=IND+1
        IND_FROM=IND_FROM+LDA
        IND_TO=IND_TO+SMP_DEFICIENCY
      END DO
      CALL ZTRSM('Left','Upper','No transpose','Non-unit',NSHORT,
     *  SMP_DEFICIENCY,ONE,A(IAPOS),LDA,VAUX_SEQ,NN)
      I=NN
      NULL_IND=IAPOS+(I-1)*LDA
      IND=NN*LDA
      JND=(SMP_DEFICIENCY-1)*NN+1
      DO II=1,SMP_DEFICIENCY
        IND=IND-LDA-1
        CALL ZCOPY(NSHORT,VAUX_SEQ(JND),1,A(NULL_IND),1)
        JND=JND-NN
        NULL_IND=NULL_IND-LDA
        I=I-1
      END DO
      DEFICIENCY=SMP_DEFICIENCY
            CALL ZMUMPS_317(NN-DEF,LDA,A(IAPOS),WK1,
     *        WK1(LOCAL_M+2,1),VAUX_SEQ,
     *        SMIN,RV)
            IF(ABS(SMIN).GT.TOL) THEN
              DEF=DEF-1
              GO TO 211
            ELSE
              DIM=root%TOT_ROOT_SIZE-DEF+1
              DEFICIENCY=DEFICIENCY+1
              GO TO 111
            END IF
 111        CONTINUE
        I=I-1
        CALL ZCOPY(I,RV,1,A(NULL_IND),1)
        NULL_IND=NULL_IND-LDA
                  CALL ZTRTRS('U','N','N',
     *              NN-DEFICIENCY,DEFICIENCY,
     *              A(iapos),LDA,a(iapos+LDA*(NN-DEFICIENCY)),
     *              LDA,IERR)
        VMAX=ZMUMPS_106(I+1,RV)
        DO K=I+1,1,-1
          ABSVK=ABS(RV(K,1))
           IF(K.EQ.I-5) GO TO 50
        END DO
 50     CONTINUE
        IF(K.LT.I+1) THEN
          KND=(K-1)*LDA+1
          CALL ZCOPY(NN,AR((K-1)*LDA+1),1,WK1,1)
          DO II=K+1,I+1
            CALL ZCOPY(NN,AR(KND+LDA),1,AR(KND),1)
            KND=KND+LDA
          END DO
          CALL ZCOPY(NN,WK1,1,AR(I*LDA+1),1)
          CALL ZMUMPS_171(K,I+1,root%IPIV)
          IND=IAPOS+LDA*(I+1)
          DO II=I+2,NN
            IND1=K
            IND2=I
            CALL ZMUMPS_105(IND1,IND2,A(IND))
            IND=IND+LDA
          END DO
          IF(MAXG.LT.GIND+I+1-K)
     *       CALL ZMUMPS_298(0,2,'SPARSLAB','SLABQR1',
     *      ' Missing space to store Givens rotations',1601,2)
          DO J=K,I
            JND=(J-1)*NN+J
            CALL ZMUMPS_159(AR(JND),AR(JND+1),C,S)
            CALL ZMUMPS_31(NN-K+1,C,S,AR(JND),LDA,AR(JND+1),LDA)
            GIND=GIND+1
            GROW(GIND)=J
            GCOS(GIND)=C
            GSIN(GIND)=S
          END DO
        END IF
        CALL ZMUMPS_317(I-1,LDA,AR,VAUX_SEQ,WK1,
     *   WK1(LOCAL_M+2,1),SMIN,RV)
        IF(ABS(SMIN).LT.TOL) THEN
          GO TO 111
        END IF
        NUM_RANK=I
        IND=IAPOS
        K=1
        J=0
        DO I=1,NN-DEFICIENCY
          J=J+1
          CALL ZCOPY(J,AR(K),1,A(IND),1)
          IND=IND+LDA
          K=K+LDA
        END DO
 211        CONTINUE
            root%GIND = GIND
            ELSE
              END IF
            ISOUR=0
            CALL MPI_BCAST(DEFICIENCY,1,MPI_INTEGER,ISOUR,
     *        COMM,IERR)
            keep(17)=deficiency
            ELSE
            CALL PZGEQPF(NN,NN,A(IAPOS),
     *      1, 1, root%DESCRIPTOR, root%IPIV(1),
     *      root%QR_TAU, WK, LWK, IERR )
            IF ( IERR .NE. 0 ) THEN
              WRITE(*,*) ' Problem in PZGEQPF : IERR = ', IERR
              CALL ZMUMPS_ABORT()
            END IF
            IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
              root%MAXG = 5000
              ALLOCATE(root%GROW(root%MAXG),stat=allocok)
              IF(allocok.GT.0) THEN
                INFO(1) = -13
                INFO(2) = N
                CALL ZMUMPS_ABORT()
              ENDIF
              ALLOCATE(root%GCOS(root%MAXG),stat=allocok)
              IF(allocok.GT.0) THEN
                INFO(1) = -13
                INFO(2) = N
                CALL ZMUMPS_ABORT()
              ENDIF
              ALLOCATE(root%GSIN(root%MAXG),stat=allocok)
              IF(allocok.GT.0) THEN
                INFO(1) = -13
                INFO(2) = N
                CALL ZMUMPS_ABORT()
              ENDIF
              MAXG = root%MAXG
              GROW => root%GROW
              GCOS => root%GCOS
              GSIN => root%GSIN
            END IF
            IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
              ALLOCATE( IPIV_SEQ(root%TOT_ROOT_SIZE), stat=allocok )
              IF(allocok.GT.0) THEN
                INFO(1) = -13
                INFO(2) = N
                CALL ZMUMPS_ABORT()
              ENDIF
            ELSE
              ALLOCATE(IPIV_SEQ(1),stat=allocok )
              IF(allocok.GT.0) THEN
                INFO(1)=-13
                INFO(2)=N
                CALL ZMUMPS_ABORT()
              ENDIF
            END IF
            CALL ZMUMPS_117(LOCAL_M+1,ZERO,WK1,1)
            CALL ZMUMPS_117(LOCAL_M+1,ZERO,WK1(LOCAL_M+2,1),1)
            RSRCV=0
            CSRCV=0
            CALL DESCINIT( DESCV,root%TOT_ROOT_SIZE+1,1,
     *        root%MBLOCK,root%NBLOCK,0,0,
     *        root%CNTXT_BLACS,root%TOT_ROOT_SIZE+1,IERR)
            IF ( IERR .ne. 0 ) THEN
              WRITE(*,*) 'After DESCINIT, IERR = ', IERR
              CALL ZMUMPS_ABORT()
            END IF
            GIND=0
            CALL ZMUMPS_173(root%TOT_ROOT_SIZE,1,
     *        LOCAL_M,1,MBLOCK,NBLOCK,root%IPIV,IPIV_SEQ,
     *        MYID,0,NPROW,NPCOL,COMM)
            DEF=1
            DIM=root%TOT_ROOT_SIZE-DEF+1
            NR=root%TOT_ROOT_SIZE
            CALL ZMUMPS_318(MYID,root%TOT_ROOT_SIZE,DEF,
     *        LOCAL_M,LOCAL_N, A(IAPOS),
     *        root%DESCRIPTOR,NPROW,NPCOL,
     *        root%CNTXT_BLACS,COMM,SMIN,RV,ZERO_DIAG,
     *        INFO(1),INFO(2),TOL)
            IF(ZERO_DIAG.EQ.1) GO TO 550
            CALL MPI_BARRIER(COMM,IERR)
            IF(ABS(SMIN).GT.TOL) THEN
              DEF=DEF-1
              GO TO 200
            ELSE
              DIM=root%TOT_ROOT_SIZE-DEF+1
              GO TO 100
            END IF
 100        CONTINUE
        call mpi_barrier(comm,ierr)
            CALL ZMUMPS_107(MYID,MYROW,MYCOL,NPROW,NPCOL,NBLOCK,1,
     *        LOCAL_M,RV,VMAX,IMAX,RMAX,COMM)
            CALL ZMUMPS_113(MYID,MYROW,MYCOL,NPROW,NPCOL,MBLOCK,
     *        NBLOCK,1,
     *        LOCAL_M,RV,VMAX,GMAX,IMAX,RMAX,SV_MU,COMM)
            CALL MPI_BCAST(GMAX,1,MPI_INTEGER,0,COMM,IERR)
            IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
              CALL ZMUMPS_171(GMAX,DIM,IPIV_SEQ)
            END IF
            CALL MPI_BARRIER(COMM,IERR)
            CALL ZMUMPS_161(MYROW,MBLOCK,NPCOL,GMAX,PRO)
            CALL ZMUMPS_160(MYROW,MBLOCK,NPCOL,GMAX,LOC)
            IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
              CALL MPI_SEND(RV(LOC,1),1,MPI_DOUBLE_COMPLEX,
     *          0,165,COMM,IERR)
            END IF
            IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
              ISOUR=PRO*NPCOL+0
              CALL MPI_RECV(TEMP,1,MPI_DOUBLE_COMPLEX,
     *          ISOUR,165,COMM,STATUS,IERR)
            END IF
            DO J=GMAX+1,DIM
              CALL ZMUMPS_161(MYROW,MBLOCK,NPCOL,J,PRO)
              CALL ZMUMPS_160(MYROW,MBLOCK,NPCOL,J,LOC)
              CALL ZMUMPS_161(MYROW,MBLOCK,NPCOL,J-1,PRO2)
              CALL ZMUMPS_160(MYROW,MBLOCK,NPCOL,J-1,LOC2)
              IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
                IDEST=PRO2*NPCOL+0
                CALL MPI_SEND(RV(LOC,1),1,MPI_DOUBLE_COMPLEX,
     *            IDEST,166,COMM,IERR)
              END IF
              IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO2) THEN
                ISOUR=PRO*NPCOL+0
                CALL MPI_RECV(RV(LOC2,1),1,MPI_DOUBLE_COMPLEX,
     *            ISOUR,166,COMM,STATUS,IERR)
              END IF
            END DO
            IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
              IDEST=PRO*NPCOL+0
              CALL MPI_SEND(TEMP,1,MPI_DOUBLE_COMPLEX,
     *          IDEST,167,COMM,IERR)
            END IF
            IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
              CALL MPI_RECV(RV(LOC,1),1,MPI_DOUBLE_COMPLEX,
     *          0,167,COMM,STATUS,IERR)
            END IF
            IVROW=INDXG2P(DIM,MBLOCK,MYROW,0,NPROW)
            IVCOL=INDXG2P(1,NBLOCK,MYCOL,0,NPCOL)
            IF(MYROW.EQ.IVROW.AND.MYCOL.EQ.IVCOL) THEN
              CALL INFOG2L(DIM,1,DESCV,NPROW,NPCOL,MYROW,MYCOL,
     *          LRINDX,LCINDX,RSRCP,CSRCP)
                ALPHA=-ONE/RV(LRINDX,LCINDX)
            END IF
            ISOUR=IVROW*NPCOL+IVCOL
            CALL MPI_BCAST(ALPHA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *        COMM,IERR)
            CALL PZSCAL(DIM,ALPHA,RV,1,1,DESCV,1)
            CALL ZMUMPS_110(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *        A(IAPOS),0,WK1,GMAX,GMAX,COLPROC,MYROW,
     *        MYCOL,NPROW,NPCOL,
     *        MYID,COMM)
            DO I=GMAX+1,DIM
              CALL ZMUMPS_110(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *          A(IAPOS),0,WK1(LOCAL_M+2,1),I,I,COLPROC,MYROW,
     *          MYCOL,NPROW,NPCOL,
     *          MYID,COMM)
              DO J=GMAX,I
                CALL ZMUMPS_161(MYROW,MBLOCK,NPCOL,J,PRO)
                CALL ZMUMPS_160(MYROW,MBLOCK,NPCOL,J,LOC)
                IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
                  CALL MPI_SEND(WK1(LOCAL_M+1+LOC,1),
     *              1,MPI_DOUBLE_COMPLEX,
     *              0,163,COMM,IERR)
                END IF
                IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
                  ISOUR=PRO*NPCOL+0
                  CALL MPI_RECV(VAUX_SEQ(J),1,MPI_DOUBLE_COMPLEX,
     *              ISOUR,163,COMM,STATUS,IERR)
                END IF
              END DO
              DO J=1,GIND
                K=GROW(J)
                IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
                  C=GCOS(J)
                  S=GSIN(J)
                  CALL ZMUMPS_31(1,C,S,VAUX_SEQ(K),1,
     *              VAUX_SEQ(K+1),1)
                END IF
              END DO
              IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
                CALL ZMUMPS_159(VAUX_SEQ(I-1),VAUX_SEQ(I),C,S)
                CALL ZMUMPS_31(1,C,S,VAUX_SEQ(I-1),1,
     *            VAUX_SEQ(I),1)
              END IF
              GIND=GIND+1
              GROW(GIND)=I-1
              IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
                GCOS(GIND)=C
                GSIN(GIND)=S
              END IF
              DO J=GMAX,I
                CALL ZMUMPS_161(MYROW,MBLOCK,NPCOL,J,PRO)
                CALL ZMUMPS_160(MYROW,MBLOCK,NPCOL,J,LOC)
                IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
                  ISOUR=PRO*NPCOL+0
                  CALL MPI_SEND(VAUX_SEQ(J),1,MPI_DOUBLE_COMPLEX,
     *              ISOUR,164,COMM,IERR)
                END IF
                IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
                  ISOUR=0
                  CALL MPI_RECV(WK1(LOCAL_M+1+LOC,1),
     *              1,MPI_DOUBLE_COMPLEX,
     *              ISOUR,164,COMM,STATUS,IERR)
                END IF
              END DO
              CALL ZMUMPS_109(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *          WK1(LOCAL_M+2,1),0,A(IAPOS),I-1,I-1,COLPROC,MYROW,
     *          MYCOL,NPROW,NPCOL,
     *          MYID,COMM)
            END DO
            DO J=GMAX,DIM
              CALL ZMUMPS_161(MYROW,MBLOCK,NPCOL,J,PRO)
              CALL ZMUMPS_160(MYROW,MBLOCK,NPCOL,J,LOC)
              IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
                CALL MPI_SEND(WK1(LOC,1),1,MPI_DOUBLE_COMPLEX,
     *            0,163,COMM,IERR)
              END IF
              IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
                ISOUR=PRO*NPCOL+0
                CALL MPI_RECV(VAUX_SEQ(J),1,MPI_DOUBLE_COMPLEX,
     *            ISOUR,163,COMM,STATUS,IERR)
              END IF
            END DO
            IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
              DO J=1,GIND
                K=GROW(J)
                C=GCOS(J)
                S=GSIN(J)
                CALL ZMUMPS_31(1,C,S,VAUX_SEQ(K),1,VAUX_SEQ(K+1),1)
              END DO
            END IF
            DO J=GMAX,DIM
              CALL ZMUMPS_161(MYROW,MBLOCK,NPCOL,J,PRO)
              CALL ZMUMPS_160(MYROW,MBLOCK,NPCOL,J,LOC)
              IF(MYCOL.EQ.0.AND.MYROW.EQ.0) THEN
                ISOUR=PRO*NPCOL+0
                CALL MPI_SEND(VAUX_SEQ(J),1,MPI_DOUBLE_COMPLEX,
     *            ISOUR,164,COMM,IERR)
              END IF
              IF(MYCOL.EQ.0.AND.MYROW.EQ.PRO) THEN
                ISOUR=0
                CALL MPI_RECV(WK1(LOC,1),1,MPI_DOUBLE_COMPLEX,
     *            ISOUR,164,COMM,STATUS,IERR)
              END IF
            END DO
            CALL ZMUMPS_109(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *        WK1,0,A(IAPOS),DIM,DIM,COLPROC,MYROW,
     *        MYCOL,NPROW,NPCOL,
     *        MYID,COMM)
 550        CONTINUE
            CALL ZMUMPS_109(LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *        RV,0,A(IAPOS),DIM,DIM-1,COLPROC,
     *        MYROW,MYCOL,NPROW,NPCOL,MYID,COMM)
            CALL MPI_BARRIER(COMM,IERR)
            DEF=DEF+1
            CALL ZMUMPS_318(MYID,NR,DEF,LOCAL_M,LOCAL_N,
     *        A(IAPOS),root%DESCRIPTOR,NPROW,NPCOL,
     *        root%CNTXT_BLACS,COMM,SMIN,RV,ZERO_DIAG,
     *        INFO(1),INFO(2),TOL)
            IF(DEFICIENCY.LE.MAXDEF) THEN
              IF(ABS(SMIN).GT.TOL) THEN
                DEF=DEF-1
                GO TO 200
              ELSEIF(ZERO_DIAG.EQ.1) THEN
                DIM=root%TOT_ROOT_SIZE-DEF+1
                GO TO 550
              ELSE
                DIM=root%TOT_ROOT_SIZE-DEF+1
                GO TO 100
              END IF
            ELSE
              DEF=DEF-1
              GO TO 200
            END IF
            CALL MPI_BARRIER(COMM,IERR)
 200        CONTINUE
            CALL ZMUMPS_172(root%TOT_ROOT_SIZE,1,
     *        LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,IPIV_SEQ,root%IPIV,
     *        MYID,0,NPROW,NPCOL,COMM)
            DEFICIENCY=DEF
            CALL MPI_BCAST(DEFICIENCY,1,MPI_INTEGER,RSRCA*NPROW+CSRCA,
     *        COMM,IERR)
            IF(MYROW.EQ.0.and.MYCOL.EQ.0) root%GIND = GIND
          END IF
          END IF
          keep(17)=deficiency
          IF ( IERR .EQ. -10 ) THEN
            WRITE(*,*) 'Not enough workspace LWK=',LWK
            CALL ZMUMPS_ABORT()
          ELSE IF ( IERR .LT. 0 ) THEN
            WRITE(*,*) 'Error code on return from PZGEQPF:',IERR
            CALL ZMUMPS_ABORT()
          END IF
        ELSEIF (LDLT.EQ.0.OR.LDLT.EQ.2) THEN
          IF(QR.EQ.0) THEN
          CALL PZGETRF( root%TOT_ROOT_SIZE, root%TOT_ROOT_SIZE,
     *      A( IAPOS ),
     *      1, 1, root%DESCRIPTOR, root%IPIV(1), IERR )
            IF ( IERR .GT. 0 ) THEN
              INFO(1)=-10
              INFO(2)=IERR-1
            END IF
          ELSEIF(QR.EQ.5.OR.QR.EQ.6) THEN
            IF(root%MBLOCK.GE.NN.AND.NOT_YET) THEN
              IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
            LDA=root%DESCRIPTOR(9)
      CALL ZGETRF(NN,NN,A(IAPOS),LDA,root%IPIV(1),IERR)
              TOL=SQRT(DBLE(NN))*EPS*ZMUMPS_104(NN,A(IAPOS),WK)
              TOL1=1.0D-8
              DEFICIENCY=0
      CALL ZMUMPS_98(NN,A(IAPOS),NN,SMP_DEFICIENCY,TOL1)
      allocate(RIGHT_SV(MAX(NN,NN*SMP_DEFICIENCY)),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      NSHORT=NN-SMP_DEFICIENCY
      IND_FROM=LDA*NSHORT+1
      IND_TO=1
      IND=0
      DO I=1,SMP_DEFICIENCY
        CALL ZCOPY(NSHORT,A(IAPOS+IND_FROM-1),
     *    1,RIGHT_SV(IND_TO),1)
        CALL ZSCAL(NSHORT,-ONE,RIGHT_SV(IND_TO),1)
        IND_TO=IND_TO+NSHORT
        CALL ZMUMPS_117(SMP_DEFICIENCY,ZERO,RIGHT_SV(IND_TO),1)
        RIGHT_SV(IND_TO+IND)= ONE
        IND=IND+1
        IND_FROM=IND_FROM+LDA
        IND_TO=IND_TO+SMP_DEFICIENCY
      END DO
      CALL ZTRSM('Left','Upper','No transpose','Non-unit',NSHORT,
     *  SMP_DEFICIENCY,ONE,A(IAPOS),LDA,RIGHT_SV,NN)
      I=NN
      NULL_IND=IAPOS+(I-1)*LDA
      IND=NN*LDA
      JND=(SMP_DEFICIENCY-1)*NN+1
      DO II=1,SMP_DEFICIENCY
        IND=IND-LDA-1
        CALL ZCOPY(NN,RIGHT_SV(JND),1,A(NULL_IND),1)
        JND=JND-NN
        NULL_IND=NULL_IND-LDA
        I=I-1
      END DO
      DEFICIENCY=SMP_DEFICIENCY
              ELSE
              END IF
            ELSE
          CALL PZGETRF( root%TOT_ROOT_SIZE, root%TOT_ROOT_SIZE,
     *      A( IAPOS ),
     *      1, 1, root%DESCRIPTOR, root%IPIV(1), IERR )
 149        CONTINUE
            call infog2l(root%TOT_ROOT_SIZE-deficiency,
     *        root%TOT_ROOT_SIZE-deficiency,
     *        root%DESCRIPTOR,nprow,npcol,
     *        myrow,mycol,lrindx,lcindx,rsrc,csrc)
            if(myrow.EQ.rsrc.and.mycol.EQ.csrc) then
              alpha=a(local_m*(lcindx-1)+lrindx+iapos-1)
            END IF
            call mpi_barrier(comm,ierr)
            ISOUR=RSRC*NPCOL+CSRC
            CALL MPI_BCAST(ALPHA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *        COMM,IERR)
            call mpi_barrier(comm,ierr)
              if(abs(alpha).le.RANK_REVEALING_TOL_1) then
                deficiency=deficiency+1
                IF(DEFICIENCY.LE.MAXDEF) GO TO 149 
              end if
            call MPI_BCAST(deficiency,1,MPI_INTEGER,rsrc*nprow+csrc,
     *        COMM,IERR)
            IF(DEFICIENCY.GT.0) THEN
              CALL ZMUMPS_261('N',
     *          root%TOT_ROOT_SIZE-DEFICIENCY,DEFICIENCY,
     *          A(iapos),1,1,root%DESCRIPTOR,
     *          root%IPIV,a(iapos),1,
     *          root%TOT_ROOT_SIZE-deficiency+1,root%DESCRIPTOR,IERR)
            END IF
            END IF
            keep(17)=deficiency
          ELSEIF(QR.EQ.7.OR.QR.EQ.8) THEN
            NN=root%TOT_ROOT_SIZE
            CALL PZGETRF(NN,NN,A(IAPOS),
     *        1,1,root%DESCRIPTOR,root%IPIV(1),IERR)
            LDA=root%DESCRIPTOR(9)
            IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
            END IF
            IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
              allocate(PERM_INIT(NN),stat=IERR)
              IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *          ' Allocation error',8,2)
              CALL ZMUMPS_174(NN,0,PERM_INIT)
              TOL=SQRT(DBLE(NN))*EPS*ZMUMPS_104(NN,A(IAPOS),WK)
              TOL1=1.0D-8
              DEFICIENCY=0
      CALL ZMUMPS_98(NN,A(IAPOS),NN,SMP_DEFICIENCY,TOL1)
      allocate(RIGHT_SV(MAX(NN,NN*SMP_DEFICIENCY)),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      NSHORT=NN-SMP_DEFICIENCY
      IND_FROM=LDA*NSHORT+1
      IND_TO=1
      IND=0
      DO I=1,SMP_DEFICIENCY
        CALL ZCOPY(NSHORT,A(IAPOS+IND_FROM-1),
     *    1,RIGHT_SV(IND_TO),1)
        CALL ZSCAL(NSHORT,-ONE,RIGHT_SV(IND_TO),1)
        IND_TO=IND_TO+NSHORT
        CALL ZMUMPS_117(SMP_DEFICIENCY,ZERO,RIGHT_SV(IND_TO),1)
        RIGHT_SV(IND_TO+IND)= ONE
        IND=IND+1
        IND_FROM=IND_FROM+LDA
        IND_TO=IND_TO+SMP_DEFICIENCY
      END DO
      CALL ZTRSM('Left','Upper','No transpose','Non-unit',NSHORT,
     *  SMP_DEFICIENCY,ONE,A(IAPOS),LDA,RIGHT_SV,NN)
      SMP_DEFICIENCY=0 !!
      I=NN
      NULL_IND=IAPOS+(I-1)*LDA
      IND=NN*LDA
      JND=(SMP_DEFICIENCY-1)*NN+1
      DO II=1,SMP_DEFICIENCY
        IND=IND-LDA-1
        CALL ZCOPY(NN,RIGHT_SV(JND),1,A(NULL_IND),1)
        JND=JND-NN
        NULL_IND=NULL_IND-LDA
        I=I-1
      END DO
      DEFICIENCY=SMP_DEFICIENCY
      EST_DEFICIENCY=1
      allocate(WN01(NN),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      allocate(PERM_COL(NN),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      CALL ZMUMPS_174(NN,0,PERM_COL)
      allocate(PERM_ROW(NN),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      CALL ZMUMPS_174(NN,0,PERM_ROW)
      allocate(ELROW(ELG_MAX),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      allocate(ELELG(ELG_MAX),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      allocate(EUROW(ELG_MAX),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      allocate(EUELG(ELG_MAX),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      allocate(PTREL(NULL_MAX),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      allocate(PTREU(NULL_MAX),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      allocate(DL(ELG_MAX),stat=IERR)
      IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *  ' Allocation error',8,2)
      ELIND=0
      EUIND=0
      PTREL(1)=1
      PTREU(1)=1
      NLUPDATE=0
      NUUPDATE=0
      CALL ZMUMPS_168(NN,A(IAPOS),LDA,DEFICIENCY,5,WK(LOCAL_M+1),
     *  PERM_COL,ELROW,ELELG,DL,
     *  EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WK,SMIN)
      CALL ZMUMPS_182(NN,A(IAPOS),LDA,DEFICIENCY,
     *  WK(LOCAL_M+1),RIGHT_SV,
     *  PERM_COL,ELROW,ELELG,DL,
     *  EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WK,SMIN)
        MU=0.1D0
        XMAX=ZMUMPS_106(I,WK(LOCAL_M+1))
        YMAX=ZMUMPS_106(I,RIGHT_SV)
        DO L=I,1,-1
          ABSYL=ABS(RIGHT_SV(L))
          IF(YMAX*MU.LT.ABSYL) GO TO 60
        END DO
 60     CONTINUE
        L=I
        ROT_STEPS=SMP_DEFICIENCY
        CALL ZMUMPS_175(NN,0,WN01,1)
        allocate(LROT(NN),stat=IERR)
        IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *    ' Allocation error',8,2)
        IND_ROT=0
        DO II=1,ROT_STEPS
          DO K=I,1,-1
            ABSXK=ABS(WK(LOCAL_M+K))
            IF(XMAX*MU.LT.ABSXK.AND.WN01(K).EQ.0) GO TO 40
          END DO
 40       CONTINUE
          IND_ROT=IND_ROT+1
          LROT(IND_ROT)=K
          WN01(K)=IND_ROT
        END DO
        ROT_STEPS=3
        LROT(1)=I-6
        LROT(2)=I-2
        LROT(3)=I-7
        rot_steps=0
        KK=L
        DO II=1,ROT_STEPS
          KK=MIN(KK,LROT(II))
        END DO
        LUSHIFT=KK-1
        DIMLU=I-LUSHIFT
        allocate(AL(DIMLU,DIMLU),stat=IERR)
        CALL ZMUMPS_117(DIMLU*DIMLU,ZERO,AL,1)
        IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *    ' Allocation error',8,2)
        allocate(AU(DIMLU,DIMLU),stat=IERR)
        CALL ZMUMPS_117(DIMLU*DIMLU,ZERO,AU,1)
        IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *    ' Allocation error',8,2)
        IF(KK.LT.I) THEN
          CALL ZMUMPS_117(DIMLU*DIMLU,ZERO,AU,1)
          DO J=KK,I
            DO II=KK,J
              TEMP=A(IAPOS-1+LDA*(J-1)+II)
              AU(II-LUSHIFT,J-LUSHIFT)=TEMP
            END DO
          END DO
          NUUPDATE=1
          IF(.NOT.(ROT_STEPS.EQ.0.OR.
     *      (ROT_STEPS.EQ.1.AND.LROT(1).EQ.I))) THEN
            DO JJ=1,ROT_STEPS
              K=LROT(JJ)
              LSHIFT=K-1
              DIML=I-LSHIFT
              NLUPDATE=NLUPDATE+1
              CALL ZCOPY(I-DIML,A(IAPOS-1+K),LDA,WK,1)
              DO II=K+1,I
                CALL ZCOPY(I-DIML,A(IAPOS-1+II),LDA,
     *            A(IAPOS-1+II-1),LDA)
              END DO
              CALL ZCOPY(I-DIML,WK,1,A(IAPOS-1+N),LDA)
              DO J=K,I
                AL(J-LSHIFT,J-LSHIFT)= ONE
                DO II=J+1,I
                  TEMP=A(IAPOS-1+LDA*(J-1)+II)
                  AL(II-LSHIFT,J-LSHIFT)=TEMP
                END DO
              END DO
              CALL ZCOPY(DIML,AL(K-LSHIFT,1),DIMLU,WK,1)
              DO II=K-LSHIFT+1,I-LSHIFT
                CALL ZCOPY(DIML,AL(II,1),DIMLU,AL(II-1,1),DIMLU)
              END DO
              CALL ZCOPY(DIML,WK,1,AL(I-LSHIFT,1),DIMLU)
              CALL ZMUMPS_171(K,I,PERM_ROW)
              IF(ELG_MAX.LT.ELIND+I-K) CALL ZMUMPS_298(0,2,
     *          'SPARSLAB','SLABRR3',
     *          ' Missing space to store elementary matrices',1604,2)
              ELIND_SAVE=ELIND+1
              DO J=K-LSHIFT,I-LSHIFT-1
                CALL ZMUMPS_177(J,AL,DIMLU,LSHIFT,ELIND,ELROW,ELELG)
              END DO
              ELIND=ELIND+1
              ELROW(ELIND)=elrow(ptrel(nlupdate))
              PTREL(NLUPDATE+1)=ELIND+1
              DO J=PTREL(NLUPDATE),ELIND-1
                II=ELROW(J)
                ELROW(J)=0
                ELG=ELELG(J)
                IF(II.GT.0) THEN
                  CALL ZMUMPS_29(DIMLU,-ELG,AU(II-LUSHIFT,1),
     *              DIMLU,AU(II-LUSHIFT+1,1),DIMLU)
                ELSE
                  CALL ZMUMPS_33(DIMLU,ELG,AU(-II-LUSHIFT,1),
     *              DIMLU,AU(-II-LUSHIFT+1,1),DIMLU)
                END IF
              END DO
              DO II=1,DIML
                DL(ELIND_SAVE)= ONE/AL(II,II)
                ELIND_SAVE=ELIND_SAVE+1
                CALL ZSCAL(DIML-II+1,ONE/AL(II,II),AL(II,II),1)
              END DO
              IND=DIMLU
              DO II=DIML,1,-1
                ELIND_SAVE=ELIND_SAVE-1
                CALL ZSCAL(DIMLU,ONE/DL(ELIND_SAVE),AU(IND,1),DIMLU)
                DL(ELIND_SAVE)= ONE
                IND=IND-1
              END DO
            DO J=K,I
              DO II=J+1,I
                TEMP=AL(II-LSHIFT,J-LSHIFT)
                A(IAPOS-1+LDA*(J-1)+II)=TEMP
              END DO
            END DO
            NLUPDATE=0
            ELIND=0
            END DO
            DO J=DIMLU-DIMLU+1,DIMLU
              DO II=DIMLU-1,J,-1 ! DANO NAZPET Z DIMLU-2
                TEMP1=AU(II,J)
                TEMP2=AU(II+1,J)
                IF(TEMP1.EQ.ZERO.AND.TEMP2.EQ.ZERO) THEN
                ELSEIF(ABS(TEMP1).GT.ABS(TEMP2)) THEN
                  CALL ZMUMPS_157(TEMP1,TEMP2,ELG)
                  CALL ZMUMPS_29(DIMLU-J+1,ELG,AU(II,J),
     *              DIMLU,AU(II+1,J),DIMLU)
                  EUIND=EUIND+1
                  EUROW(EUIND)=II+LUSHIFT
                  EUELG(EUIND)=ELG
                ELSE
                  CALL ZMUMPS_157(AU(II+1,J),AU(II,J),ELG)
                  CALL ZMUMPS_30(DIMLU-J+1,ELG,AU(II,J),
     *              DIMLU,AU(II+1,J),DIMLU)
                  EUIND=EUIND+1
                  EUROW(EUIND)=-(II+LUSHIFT)
                  EUELG(EUIND)=ELG
                END IF
              END DO
            END DO
            deallocate(AL)
          END IF
        END IF
        IF(L.LT.I) THEN
          USHIFT=L-1
          DIMU=I-USHIFT
          KND=(L-1)*N+1
          CALL ZCOPY(I-DIMLU,A(IAPOS-1+(L-1)*LDA+1),
     *      1,WK,1)
          DO II=L+1,I
            CALL ZCOPY(I-DIMLU,A(IAPOS-1+KND+N),
     *        1,A(IAPOS-1+KND),1)
            KND=KND+N
          END DO
          CALL ZCOPY(I-DIMLU,WK,1,A(IAPOS-1+(I-1)*LDA+1),1)
          CALL ZCOPY(DIMLU,AU(1,DIMLU-DIMU+1),1,WK,1)
          DO II=L-USHIFT+DIMLU-DIMU+1,I-USHIFT+DIMLU-DIMU
            CALL ZCOPY(DIMLU,AU(1,II),1,AU(1,II-1),1)
          END DO
          CALL ZCOPY(DIMLU,WK,1,AU(1,I-USHIFT+DIMLU-DIMU),1)
          CALL ZMUMPS_171(L,I,PERM_COL)
          IND=(NULL_MAX-NN+I-1)*NN+1
          DO II=NULL_MAX-NN+I,NULL_MAX
            IND1=L
            IND2=I
            CALL ZMUMPS_105(IND1,IND2,A(IND))
            IND=IND+N
          END DO
          IF(ELG_MAX.LT.EUIND+I-L) 
     *      CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *      ' Missing space to store elementary matrices',1604,2)
          DO J=L-USHIFT+DIMLU-DIMU,I-USHIFT-1+DIMLU-DIMU ! DANO ZPET
            CALL ZMUMPS_331(J,AU,DIMLU,LUSHIFT,EUIND,EUROW,EUELG)
          END DO
        END IF
        PTREU(NUUPDATE+1)=EUIND+1
        IF(KK.LT.I) THEN
          DO J=KK,I
            DO II=KK,J
              TEMP=AU(II-LUSHIFT,J-LUSHIFT)
              A(IAPOS-1+LDA*(J-1)+II)=TEMP
            END DO
          END DO
        END IF
        deallocate(AU)
      CALL ZMUMPS_117(DEFICIENCY,ZERO,WK(LOCAL_M+NN-DEFICIENCY+1),1)
      CALL ZMUMPS_117(DEFICIENCY,ZERO,RIGHT_SV(NN-DEFICIENCY+1),1)
      CALL ZMUMPS_168(NN,A(IAPOS),LDA,DEFICIENCY,5,WK(LOCAL_M+1),
     *  PERM_COL,ELROW,ELELG,DL,
     *  EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WK,SMIN)
       CALL ZMUMPS_182(NN,A(IAPOS),LDA,DEFICIENCY,
     *  WK(LOCAL_M+1),RIGHT_SV,
     *  PERM_COL,ELROW,ELELG,DL,
     *  EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WK,SMIN)
      IF(ABS(SMIN).LE.TOL1) THEN
        IF(NN-I.GT.NULL_MAX)
     *    CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *    ' Missing space to store null space basis',1602,2)
      ELSE
        GO TO 210
      END IF
      CALL ZMUMPS_334(NN,PERM_COL,RIGHT_SV,WK(LOCAL_M+1))
      ALPHA=RIGHT_SV(I)
      ALPHA=-ONE/ALPHA
      CALL ZSCAL(I,ALPHA,RIGHT_SV,1)
 110  CONTINUE
        NULL_IND=IAPOS+(I-1)*LDA
        CALL ZMUMPS_117(I-1,ZERO,A(NULL_IND),1)
        CALL ZCOPY(I-1,RIGHT_SV,1,A(NULL_IND),1)
        NULL_IND=NULL_IND-LDA
        DEFICIENCY=DEFICIENCY+1
        I=I-1
        MU=0.1D0
        YMAX=ZMUMPS_106(I,RIGHT_SV)
        DO L=I,1,-1
          ABSYL=ABS(RIGHT_SV(L))
          IF(YMAX*MU.LT.ABSYL) GO TO 65
        END DO
 65     CONTINUE
        L=I
        KK=L
        LUSHIFT=KK-1
        DIMLU=I-LUSHIFT
        allocate(AU(DIMLU,DIMLU),stat=IERR)
        CALL ZMUMPS_117(DIMLU*DIMLU,ZERO,AU,1)
        IF(IERR.NE.0) CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *    ' Allocation error',8,2)
        IF(L.LT.I) THEN
          CALL ZMUMPS_117(DIMLU*DIMLU,ZERO,AU,1)
          DO J=KK,I
            DO II=KK,J
              TEMP=A(IAPOS-1+LDA*(J-1)+II)
              AU(II-LUSHIFT,J-LUSHIFT)=TEMP
            END DO
          END DO
          USHIFT=L-1
          DIMU=I-USHIFT
          KND=(L-1)*N+1
          CALL ZCOPY(I-DIMLU,A(IAPOS-1+(L-1)*LDA+1),1,WK,1)
          DO II=L+1,I
            CALL ZCOPY(I-DIMLU,
     *        A(IAPOS-1+KND+LDA),1,A(IAPOS-1+KND),1)
            KND=KND+LDA
          END DO
          CALL ZCOPY(I-DIMLU,WK,1,A(IAPOS-1+(I-1)*LDA+1),1)
          CALL ZCOPY(DIMLU,AU(1,DIMLU-DIMU+1),1,WK,1)
          DO II=L-USHIFT+DIMLU-DIMU+1,I-USHIFT+DIMLU-DIMU
            CALL ZCOPY(DIMLU,AU(1,II),1,AU(1,II-1),1)
          END DO
          CALL ZCOPY(DIMLU,WK,1,AU(1,I-USHIFT+DIMLU-DIMU),1)
          CALL ZMUMPS_171(L,I,PERM_COL)
          IND=(NULL_MAX-NN+I-1)*N+1
          DO II=NULL_MAX-NN+I,NULL_MAX
            IND1=L
            IND2=I
            CALL ZMUMPS_105(IND1,IND2,A(IND))
            IND=IND+N
          END DO
          IF(ELG_MAX.LT.EUIND+I-L)
     *      CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *      ' Missing space to store elementary matrices',1604,2)
          DO J=L-USHIFT+DIMLU-DIMU,I-USHIFT-1+DIMLU-DIMU ! DANO ZPET
            CALL ZMUMPS_331(J,AU,DIMLU,LUSHIFT,EUIND,EUROW,EUELG)
          END DO
          PTREU(NUUPDATE+1)=EUIND+1
          DO J=KK,I
            DO II=KK,J
              TEMP=AU(II-LUSHIFT,J-LUSHIFT)
              A(IAPOS-1+LDA*(J-1)+II)=TEMP
            END DO
          END DO
        END IF
        deallocate(AU)
      CALL ZMUMPS_117(DEFICIENCY,ZERO,WK(LOCAL_M+NN-DEFICIENCY+1),1)
      CALL ZMUMPS_117(DEFICIENCY,ZERO,RIGHT_SV(NN-DEFICIENCY+1),1)
      CALL ZMUMPS_317(I,NN,A(IAPOS),WK(LOCAL_M+1),WK,
     *  WK(LOCAL_M+1),SMIN,RIGHT_SV)
      CALL ZCOPY(NN,RIGHT_SV,1,WK(LOCAL_M+1),1)
      CALL ZMUMPS_333(NN,PERM_COL,WK(LOCAL_M+1),WK)
      IF(ABS(SMIN).LT.TOL1) THEN
        IF(NN-I.GT.NULL_MAX)
     *    CALL ZMUMPS_298(0,2,'SPARSLAB','SLABRR3',
     *    ' Missing space to store null space basis',1602,2)
        GO TO 110
      END IF
 210  CONTINUE
      NUM_RANK=I
            ELSE
          END IF
            CALL MPI_BARRIER(COMM,IERR)
            IF(root%MBLOCK.LT.NN.AND.NOT_YET) THEN
             call BLACS_GRIDINFO(root%DESCRIPTOR(2),nprow,npcol,
     *         myrow,mycol)
 150       CONTINUE
           call infog2l(root%TOT_ROOT_SIZE-deficiency,
     *         root%TOT_ROOT_SIZE-deficiency,
     *         root%DESCRIPTOR,nprow,npcol,
     *        myrow,mycol,lrindx,lcindx,rsrc,csrc)
           if(myrow.EQ.rsrc.and.mycol.EQ.csrc) then
             alpha=a(local_m*(lcindx-1)+lrindx+iapos-1)
           END IF
           call mpi_barrier(comm,ierr)
           ISOUR=RSRC*NPCOL+CSRC
           CALL MPI_BCAST(ALPHA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *       COMM,IERR)
           call mpi_barrier(comm,ierr)
             if(abs(alpha).le.RANK_REVEALING_TOL_1) then
               deficiency=deficiency+1
                IF(DEFICIENCY.LE.MAXDEF) GO TO 150
             end if
           call MPI_BCAST(deficiency,1,MPI_INTEGER,rsrc*nprow+csrc,
     *       COMM,IERR)
           IF(DEFICIENCY.GT.0) THEN
             CALL ZMUMPS_261('N',
     *         root%TOT_ROOT_SIZE-DEFICIENCY,DEFICIENCY,
     *         A(iapos),1,1,root%DESCRIPTOR,
     *         root%IPIV,a(iapos),1,
     *         root%TOT_ROOT_SIZE-deficiency+1,root%DESCRIPTOR,IERR)
           END IF
           KEEP(17)=DEFICIENCY
          ANORM=PDLANGE('1',root%TOT_ROOT_SIZE-DEFICIENCY,
     *      root%TOT_ROOT_SIZE-DEFICIENCY,
     *      A(IAPOS),1,1,root%DESCRIPTOR,WK)
          LIWK=LIW-IFREE+1
          CALL PZGECON('1',root%TOT_ROOT_SIZE-DEFICIENCY,
     *      A(IAPOS),1,1,
     *      root%DESCRIPTOR,ANORM,RCOND,WK,LWK,IW(IFREE),
     *      LIWK,IERR)
          RANK_REVEALING_TOL_2=ANORM*(ONE/RCOND)
          RANK_REVEALING_TOL_2=
     *      RANK_REVEALING_TOL_1*(root%TOT_ROOT_SIZE)**0.5D0*
     *      RANK_REVEALING_TOL_2
           RANK_REVEALING_TOL_2=MIN(RANK_REVEALING_TOL_2,1.0D-8)
      EST_DEFICIENCY=1
        RSRCU=0
        CSRCU=0
        CALL DESCINIT(DESCU,root%TOT_ROOT_SIZE+1,1,
     *    root%MBLOCK,root%NBLOCK,0,0,root%CNTXT_BLACS,
     *    root%TOT_ROOT_SIZE+1,IERR)
        IF(IERR.NE.0) THEN
          WRITE(*,*) ' After DESCINIT, IERR = ',IERR
          CALL ZMUMPS_ABORT()
        END IF
      IF(MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
      allocate(PERM_COL(root%TOT_ROOT_SIZE),stat=IERR)
        IF (IERR.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = N
          CALL ZMUMPS_ABORT()
        ENDIF
      CALL ZMUMPS_174(root%TOT_ROOT_SIZE,0,PERM_COL)
      allocate(PERM_ROW(root%TOT_ROOT_SIZE),stat=IERR)
        IF (IERR.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = N
          CALL ZMUMPS_ABORT()
        ENDIF
      CALL ZMUMPS_174(root%TOT_ROOT_SIZE,0,PERM_ROW)
      allocate(ELROW(ELG_MAX),stat=IERR)
        IF (IERR.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = N
          CALL ZMUMPS_ABORT()
        ENDIF
      allocate(ELELG(ELG_MAX),stat=IERR)
        IF (IERR.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = N
          CALL ZMUMPS_ABORT()
        ENDIF
      allocate(EUROW(ELG_MAX),stat=IERR)
        IF (IERR.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = N
          CALL ZMUMPS_ABORT()
        ENDIF
      allocate(EUELG(ELG_MAX),stat=IERR)
        IF (IERR.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = N
          CALL ZMUMPS_ABORT()
        ENDIF
      allocate(PTREL(NULL_MAX),stat=IERR)
        IF (IERR.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = N
          CALL ZMUMPS_ABORT()
        ENDIF
      allocate(PTREU(NULL_MAX),stat=IERR)
        IF (IERR.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = N
          CALL ZMUMPS_ABORT()
        ENDIF
      allocate(DL(ELG_MAX),stat=IERR)
        IF (IERR.GT.0 ) THEN
          INFO(1) = -13
          INFO(2) = N
          CALL ZMUMPS_ABORT()
        ENDIF
      ELIND=0
      EUIND=0
      PTREL(1)=1
      PTREU(1)=1
      NLUPDATE=0
      NUUPDATE=0
      END IF
      CALL ZMUMPS_169(MYID,MYROW,MYCOL,LOCAL_M,LOCAL_N,
     *  MBLOCK,NBLOCK,NPROW,NPCOL,
     *  root%TOT_ROOT_SIZE,A(IAPOS),
     *  root%DESCRIPTOR,DEFICIENCY,5,WK(LOCAL_M+1),DESCU,
     *  PERM_COL,ELROW,ELELG,DL,
     *  EUROW,EUELG,NLUPDATE,NUUPDATE,PTREL,PTREU,WK,SMIN)
           ENVELOPE=root%TOT_ROOT_SIZE-DEFICIENCY
           CURRENT_DEFICIENCY=0
 151       CONTINUE
           ADD_COLUMN=0
           DO I=ENVELOPE,root%TOT_ROOT_SIZE-DEFICIENCY
             CALL INFOG2L(root%TOT_ROOT_SIZE-DEFICIENCY,
     *         root%TOT_ROOT_SIZE-DEFICIENCY,
     *         root%DESCRIPTOR,NPROW,NPCOL,
     *         MYROW,MYCOL,LRINDX,LCINDX,RSRC,CSRC)
             if(myrow.EQ.rsrc.and.mycol.EQ.csrc) then
               alpha=a(local_m*(lcindx-1)+lrindx+iapos-1)
             END IF
             ISOUR=RSRC*NPCOL+CSRC
             CALL MPI_BCAST(ALPHA,1,MPI_DOUBLE_COMPLEX,ISOUR,
     *         COMM,IERR)
             IF(ABS(ALPHA).LE.RANK_REVEALING_TOL_2) then
               ADD_COLUMN=1
             END IF
           END DO
           IF(ADD_COLUMN.EQ.1) THEN
             CURRENT_DEFICIENCY=CURRENT_DEFICIENCY+1
             ENVELOPE=ENVELOPE-1
             GO TO 151
           END IF
           END IF
      ISOUR=0
      CALL MPI_BCAST(DEFICIENCY,1,MPI_INTEGER,ISOUR,
     *  COMM,IERR)
           KEEP(17)=DEFICIENCY
          END IF
        ELSE
          CALL PZPOTRF('L',root%TOT_ROOT_SIZE,A(IAPOS),
     *      1,1,root%DESCRIPTOR,IERR)
            IF ( IERR .GT. 0 ) THEN
              INFO(1)=-10
              INFO(2)=IERR-1
            END IF
          END IF
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_147
#endif
      SUBROUTINE ZMUMPS_146( MYID, root, N, IROOT,
     *           COMM, IW, LIW, IFREE,
     *           A, LA, PTRAST, PTLUST_S, PTRFAC,
     *           STEP, INFO, LDLT, QR,
     *           WK, LWK, KEEP,KEEP8)
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      INCLUDE 'mpif.h'
      TYPE ( ZMUMPS_ROOT_STRUC ) :: root
      INTEGER N, IROOT, COMM, LIW, LA, MYID, LWK, LIWK, IFREE
      COMPLEX*16 WK( LWK )
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRAST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *STEP(N), IW( LIW )
      INTEGER INFO( 2 ), LDLT, QR
      COMPLEX*16 A( LA )
      INTEGER IOLDPS, IAPOS
      INTEGER LOCAL_M, LOCAL_N, LPIV, IERR, allocok,i
      INTEGER XSIZE
      PARAMETER(XSIZE=3)
        IF ( .NOT. root%yes ) RETURN
        IF ( KEEP(60) .NE. 0 ) THEN
          IF ((LDLT == 1 .OR. LDLT == 2) .AND. KEEP(60) == 3 ) THEN
            CALL ZMUMPS_320( WK, root%MBLOCK,
     *      root%MYROW, root%MYCOL, root%NPROW, root%NPCOL,
     *      root%SCHUR_POINTER(1),
     *      root%SCHUR_LLD, root%SCHUR_NLOC,
     *      root%TOT_ROOT_SIZE, MYID, COMM )
          ENDIF
        RETURN
        ENDIF
        IOLDPS  = PTLUST_S(STEP(IROOT))+XSIZE
        IAPOS   = PTRAST(STEP(IROOT))
        LOCAL_M = IW( IOLDPS + 2 )
        LOCAL_N = IW( IOLDPS + 1 )
        IAPOS = PTRFAC(IW ( IOLDPS + 4 ))
        IF ( LDLT.EQ.0 .OR. LDLT.EQ.2 .OR. QR.ne.0 ) THEN
         LPIV = LOCAL_M + root%MBLOCK
        ELSE
         LPIV = 1
        END IF
        IF (ASSOCIATED( root%IPIV )) DEALLOCATE(root%IPIV)
        root%LPIV = LPIV
        ALLOCATE( root%IPIV( LPIV ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = LPIV
          WRITE(*,*) MYID,': problem allocating IPIV(',LPIV,') in root'
          CALL ZMUMPS_ABORT()
        END IF
        CALL DESCINIT( root%DESCRIPTOR, root%TOT_ROOT_SIZE,
     *      root%TOT_ROOT_SIZE, root%MBLOCK, root%NBLOCK,
     *      0, 0, root%CNTXT_BLACS, LOCAL_M, IERR )
        IF ( LDLT.EQ.2 ) THEN
            IF(root%MBLOCK.NE.root%NBLOCK) THEN
              WRITE(*,*) ' Error: symmetrization only works for'
              WRITE(*,*) ' square block sizes, MBLOCK/NBLOCK=',
     *        root%MBLOCK, root%NBLOCK
              CALL ZMUMPS_ABORT()
            END IF
            IF ( LWK .LT.
     *        MIN(root%MBLOCK * root%NBLOCK,root%TOT_ROOT_SIZE*
     *        root%TOT_ROOT_SIZE ) ) THEN
              WRITE(*,*) 'Not enough workspace for symmetrization.'
              CALL ZMUMPS_ABORT()
            END IF
            CALL ZMUMPS_320( WK, root%MBLOCK,
     *      root%MYROW, root%MYCOL, root%NPROW, root%NPCOL,
     *      A( IAPOS ), LOCAL_M, LOCAL_N,
     *      root%TOT_ROOT_SIZE, MYID, COMM )
        END IF
        IF (LDLT.EQ.0.OR.LDLT.EQ.2) THEN
          CALL PZGETRF( root%TOT_ROOT_SIZE, root%TOT_ROOT_SIZE,
     *      A( IAPOS ),
     *      1, 1, root%DESCRIPTOR, root%IPIV(1), IERR )
          IF ( IERR .GT. 0 ) THEN
              INFO(1)=-10
              INFO(2)=IERR-1
          END IF
        ELSE
          CALL PZPOTRF('L',root%TOT_ROOT_SIZE,A(IAPOS),
     *      1,1,root%DESCRIPTOR,IERR)
            IF ( IERR .GT. 0 ) THEN
              INFO(1)=-10
              INFO(2)=IERR-1
            END IF
        END IF
        RETURN
      END SUBROUTINE ZMUMPS_146
      SUBROUTINE ZMUMPS_556(
     *     N,PIV,FRERE,FILS,NFSIZ,IKEEP,
     *     NCST,KEEP,KEEP8,id)
      USE ZMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (ZMUMPS_STRUC) :: id
      INTEGER N,NCST
      INTEGER PIV(N),FRERE(N),FILS(N),NFSIZ(N),IKEEP(N,3)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER I,P11,P1,P2,K1,K2,NLOCKED
      LOGICAL V1,V2
      NCST = 0
      NLOCKED = 0
      P11 = KEEP(93)
      DO I=KEEP(93)-1,1,-2
         P1 = PIV(I)
         P2 = PIV(I+1)
         K1 = IKEEP(P1,1)
         IF(K1 .GT. 0) THEN
            V1 = (ABS(id%A(K1)*(id%ROWSCA(P1)**2)).GE.1.0D-1)
         ELSE
            V1 = .FALSE.
         ENDIF
         K2 = IKEEP(P2,1)
         IF(K2 .GT. 0) THEN
            V2 = (ABS(id%A(K2)*(id%ROWSCA(P2)**2)).GE.1.0D-1)
         ELSE
            V2 = .FALSE.
         ENDIF
         IF(V1 .AND. V2) THEN
            PIV(P11) = P1
            P11 = P11 - 1
            PIV(P11) = P2
            P11 = P11 - 1
         ELSE IF(V1) THEN
            NCST = NCST+1
            FRERE(NCST) = P1
            NCST = NCST+1
            FRERE(NCST) = P2
         ELSE IF(V2) THEN
            NCST = NCST+1
            FRERE(NCST) = P2                
            NCST = NCST+1
            FRERE(NCST) = P1
         ELSE
            NLOCKED = NLOCKED + 1
            FILS(NLOCKED) = P1
            NLOCKED = NLOCKED + 1
            FILS(NLOCKED) = P2                   
         ENDIF
      ENDDO
      DO I=1,NLOCKED
         PIV(I) = FILS(I)
      ENDDO
      KEEP(94) = KEEP(94) + KEEP(93) - NLOCKED
      KEEP(93) = NLOCKED
      DO I=1,NCST
         NLOCKED = NLOCKED + 1
         PIV(NLOCKED) = FRERE(I)
      ENDDO
      DO I=1,KEEP(93)/2
         NFSIZ(I) = 0
      ENDDO
      DO I=(KEEP(93)/2)+1,(KEEP(93)/2)+NCST,2
         NFSIZ(I) = I+1
         NFSIZ(I+1) = -1
      ENDDO
      DO I=(KEEP(93)/2)+NCST+1,(KEEP(93)/2)+KEEP(94)
         NFSIZ(I) = 0
      ENDDO
      END SUBROUTINE ZMUMPS_556
      SUBROUTINE ZMUMPS_548(N,PE,NV,WORK)
      IMPLICIT NONE
      INTEGER N
      INTEGER PE(N),NV(N),WORK(N)
      INTEGER I,FATHER,LEN,K,NEWSON,NEWFATHER
      DO I=1,N
         IF(NV(I) .GT. 0) CYCLE
         LEN = 1
         WORK(LEN) = I
         FATHER = -PE(I)
         DO
            IF(NV(FATHER) .GT. 0) THEN
               NEWSON = FATHER
               EXIT
            ENDIF
            LEN = LEN + 1
            WORK(LEN) = FATHER
            NV(FATHER) = 1
            FATHER = -PE(FATHER)
         ENDDO
         NEWFATHER = -PE(FATHER)
         PE(WORK(LEN)) = -NEWFATHER
         PE(NEWSON) = -WORK(1)
      ENDDO      
      END SUBROUTINE ZMUMPS_548
      SUBROUTINE ZMUMPS_549(N,PE,PERM,NFILS,WORK)
      IMPLICIT NONE
      INTEGER N
      INTEGER PE(N),PERM(N),NFILS(N),WORK(N)
      INTEGER I,FATHER,STKLEN,STKPOS,PERMPOS,CURVAR
      NFILS = 0
      DO I=1,N
         FATHER = -PE(I)
         IF(FATHER .NE. 0) NFILS(FATHER) = NFILS(FATHER) + 1
      ENDDO
      STKLEN = 0
      PERMPOS = 1
      DO I=1,N
         IF(NFILS(I) .EQ. 0) THEN
            STKLEN = STKLEN + 1
            WORK(STKLEN) = I
            PERM(I) = PERMPOS
            PERMPOS = PERMPOS + 1
         ENDIF
      ENDDO
      DO STKPOS = 1,STKLEN
         CURVAR = WORK(STKPOS)
         FATHER = -PE(CURVAR)
         DO
            IF(FATHER .EQ. 0) EXIT
            IF(NFILS(FATHER) .EQ. 1) THEN
               PERM(FATHER) = PERMPOS
               FATHER = -PE(FATHER)
               PERMPOS = PERMPOS + 1
            ELSE
               NFILS(FATHER) = NFILS(FATHER) - 1
               EXIT
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE ZMUMPS_549
      SUBROUTINE ZMUMPS_567(N,NCMP,N11,N22,PIV,
     *     PE,PE_EXP)
      IMPLICIT NONE
      INTEGER N11,N22,N,NCMP
      INTEGER PIV(N),PE(NCMP),PE_EXP(N)
      INTEGER I,J,N2,K,THEROOT,FATHER,K1,K2
      N2 = N22/2
#if (SP_DEBUG >= 5)
      WRITE(*,*) 'PIV',PIV
      WRITE(*,*) 'PE',PE
#endif
      IF(N22+N11 .EQ. N) THEN
         THEROOT = 0
      ELSE
         THEROOT = PIV(N)
         PE_EXP(THEROOT) = 0
      ENDIF
      K = 1
      DO I=1,N2
         K1 = PIV(K)
         K = K+1
         K2 = PIV(K)
         K = K+1
         FATHER = -PE(I)
         IF(FATHER .EQ. 0) THEN
            FATHER = THEROOT
         ELSE IF(FATHER .LE. N2) THEN
            FATHER = PIV(2*FATHER)
         ELSE
            FATHER = PIV(N2+FATHER)
         ENDIF
         PE_EXP(K1) = -K2
         PE_EXP(K2) = -FATHER
      ENDDO
#if(SP_DEBUG >=1)
      IF(K .NE. N22+1) THEN
         WRITE(*,*) 'K NE N22+1',K,N22+1
         CALL ZMUMPS_ABORT()
      ENDIF
#endif
      DO I=N2+1,NCMP
         K1 = PIV(K)
         K = K+1
         FATHER = -PE(I)
         IF(FATHER .EQ. 0) THEN
            FATHER = THEROOT
         ELSE IF(FATHER .LE. N2) THEN
            FATHER = PIV(2*FATHER)
         ELSE
            FATHER = PIV(N2+FATHER)
         ENDIF
         PE_EXP(K1) = -FATHER
      ENDDO
      DO I=N11+N22+1,N-1
         K1 = PIV(I)
         PE(K1) = -THEROOT
      ENDDO
#if (SP_DEBUG >= 5)
      WRITE(*,*) 'PE_EXP',PE_EXP
#endif
      END SUBROUTINE ZMUMPS_567
      SUBROUTINE ZMUMPS_566(N,LIW,COLPTR,IW,N11,N22,PIV,
     *     PE,NV,WORK)
      IMPLICIT NONE
      INTEGER N11,N22,N,LIW
      INTEGER PIV(N),NV(N),PE(N),COLPTR(N+1),IW(LIW),WORK(N)
      INTEGER I,J,K,COL,FLAG,THEROOT
      FLAG = 1
      WORK = 0
      THEROOT = PIV(N)
      DO I = N11+N22+1,N
         J = PIV(I)
         WORK(THEROOT) = FLAG
         DO K=COLPTR(J),COLPTR(J+1)-1
            COL = IW(K)
            DO
               IF(NV(COL) .GT. 0) EXIT
               COL = -PE(COL)
            ENDDO
            DO
               IF(WORK(COL) .EQ. FLAG) EXIT
               NV(COL) = NV(COL) + 1
               WORK(COL) = FLAG
               COL = -PE(COL)
            ENDDO
         ENDDO
         IF(FLAG+1 .GT. FLAG) THEN
            FLAG = FLAG+1
         ELSE
            WORK = 0
            FLAG = 1
         ENDIF
      ENDDO
#if (SP_DEBUG >= 2)
      WRITE(*,*) 'PIV',PIV
      WRITE(*,*) 'PE',PE
      WRITE(*,*) 'NV',NV
#endif
      END SUBROUTINE ZMUMPS_566
      SUBROUTINE ZMUMPS_565(N,NCMP,N11,N22,PIV,
     *     PE,NV,PE_EXP,NV_EXP)
      IMPLICIT NONE
      INTEGER N11,N22,N,NCMP
      INTEGER PIV(N),NV(NCMP),PE(NCMP),NV_EXP(N),PE_EXP(N)
      INTEGER I,J,N2,K,THEROOT,FATHER,K1,K2
      N2 = N22/2
#if (SP_DEBUG >= 5)
      WRITE(*,*) 'PIV',PIV
      WRITE(*,*) 'PE',PE
      WRITE(*,*) 'NV',NV
#endif
      IF(N22+N11 .EQ. N) THEN
         THEROOT = 0
      ELSE
         THEROOT = PIV(N)
         PE_EXP(THEROOT) = 0
         NV_EXP(THEROOT) = N-N22-N11
      ENDIF
      K = 1
      DO I=1,N2
         K1 = PIV(K)
         K = K+1
         K2 = PIV(K)
         K = K+1
         NV_EXP(K1) = 0
         NV_EXP(K2) = NV(I)
         FATHER = -PE(I)
         IF(FATHER .EQ. 0) THEN
            FATHER = THEROOT
         ELSE IF(FATHER .LE. N2) THEN
            FATHER = PIV(2*FATHER)
         ELSE
            FATHER = PIV(N2+FATHER)
         ENDIF
         PE_EXP(K1) = -K2
         PE_EXP(K2) = -FATHER
      ENDDO
#if(SP_DEBUG >=1)
      IF(K .NE. N22+1) THEN
         WRITE(*,*) 'K NE N22+1',K,N22+1
         CALL ZMUMPS_ABORT()
      ENDIF
#endif
      DO I=N2+1,NCMP
         K1 = PIV(K)
         K = K+1
         NV_EXP(K1) = NV(I)
         FATHER = -PE(I)
         IF(FATHER .EQ. 0) THEN
            FATHER = THEROOT
         ELSE IF(FATHER .LE. N2) THEN
            FATHER = PIV(2*FATHER)
         ELSE
            FATHER = PIV(N2+FATHER)
         ENDIF
         PE_EXP(K1) = -FATHER
      ENDDO
      DO I=N11+N22+1,N-1
         K1 = PIV(I)
         NV(K1) = 0
         PE(K1) = -THEROOT
      ENDDO
#if (SP_DEBUG >= 5)
      WRITE(*,*) 'PE_EXP',PE_EXP
      WRITE(*,*) 'NV_EXP',NV_EXP
#endif
      END SUBROUTINE ZMUMPS_565
      SUBROUTINE ZMUMPS_550(N,NCMP,N11,N22,PIV,
     *     PERM,IPERM)
      IMPLICIT NONE
      INTEGER N11,N22,N,NCMP
      INTEGER PIV(N),PERM(N),IPERM(N)
      INTEGER CMP_POS,EXP_POS,I,J,N2,K
      N2 = N22/2
#if(SP_DEBUG >= 5)
      WRITE(*,*) 'PIV',PIV
      WRITE(*,*) 'PERM',PERM(1:NCMP)
      WRITE(*,*) 'IPERM',IPERM(1:NCMP)
#endif
      EXP_POS = 1
      DO CMP_POS=1,NCMP
         J = IPERM(CMP_POS)
         IF(J .LE. N2) THEN
            K = 2*J-1
            I = PIV(K)
            PERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
            K = K+1
            I = PIV(K)
            PERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
         ELSE
            K = N2 + J
            I = PIV(K)
            PERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
         ENDIF
      ENDDO
      DO K=N22+N11+1,N
         I = PIV(K)
         PERM(I) = EXP_POS
         EXP_POS = EXP_POS+1
      ENDDO
#if(SP_DEBUG >= 5)
      WRITE(*,*) 'PERM',PERM
#endif
      END SUBROUTINE ZMUMPS_550
      SUBROUTINE ZMUMPS_547(
     *     N,NZ, IRN, ICN, PIV,
     *     NCMP, IW, LW, IPE, LEN, IQ, 
     *     FLAG, ICMP, IWFR,
     *     IERROR, KEEP,KEEP8, ICNTL)
      IMPLICIT NONE
      INTEGER N,NZ,NCMP,LW,IWFR,IERROR
      INTEGER ICNTL(40),KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IRN(NZ),ICN(NZ),IW(LW),PIV(N),IPE(N+1)
      INTEGER LEN(N),IQ(N),FLAG(N),ICMP(N) 
      INTEGER MP,N11,N22,NDUP
      INTEGER I,K,J,N1,LAST,K1,K2,L
      MP = ICNTL(2)
      IERROR = 0
      N22 = KEEP(93)
      N11 = KEEP(94)
      NCMP = N22/2 + N11
      DO I=1,NCMP
         IPE(I) = 0
      ENDDO
      K = 1
      DO I=1,N22/2
         J = PIV(K)
         ICMP(J) = I
         K = K + 1
         J = PIV(K)
         ICMP(J) = I
         K = K + 1
      ENDDO
      K = N22/2 + 1
      DO I=N22+1,N22+N11
         J = PIV(I)
         ICMP(J) = K
         K = K + 1
      ENDDO
      DO I=N11+N22+1,N
         J = PIV(I)
         ICMP(J) = 0
      ENDDO
      DO K=1,NZ
         I = IRN(K)
         J = ICN(K)
         I = ICMP(I)
         J = ICMP(J)
         IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *        .OR.(J.LT.1)) THEN
            IERROR = IERROR + 1
         ELSE
            IF (I.NE.J) THEN
               IPE(I) = IPE(I) + 1
               IPE(J) = IPE(J) + 1
            ENDIF
         ENDIF
      ENDDO
      IQ(1) = 1
      N1 = NCMP - 1
      IF (N1.GT.0) THEN
         DO I=1,N1
            IQ(I+1) = IPE(I) + IQ(I)
         ENDDO 
      ENDIF
      LAST = MAX(IPE(NCMP)+IQ(NCMP)-1,IQ(NCMP))
      DO I = 1,NCMP
         FLAG(I) = 0
         IPE(I)  = IQ(I)
      ENDDO
      DO K=1,LAST
        IW(K) = 0
      ENDDO
      IWFR = LAST + 1
      DO K=1,NZ
         I = IRN(K)
         J = ICN(K)
         I = ICMP(I)
         J = ICMP(J)
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
            IF ((I.GE.1).AND.(J.LE.N)) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1 
            ENDIF
          ELSE
            IF ((J.GE.1).AND.(I.LE.N)) THEN
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1
            ENDIF
          ENDIF
         ENDIF
      ENDDO
      NDUP = 0
      DO I=1,NCMP
         K1 = IPE(I) 
         K2 = IQ(I) -1
         IF (K1.GT.K2) THEN
            LEN(I) = 0
            IQ(I)  = 0
         ELSE
            DO K=K1,K2
               J     = -IW(K)
               IF (J.LE.0) GO TO 250
               L     = IQ(J) 
               IQ(J) = L + 1
               IF (FLAG(J).EQ.I) THEN
                  NDUP = NDUP + 1
                  IW(L) = 0
                  IW(K) = 0
               ELSE
                  IW(L)   = I
                  IW(K)   = J
                  FLAG(J) = I
               ENDIF
            ENDDO
 250        IQ(I) = IQ(I) - IPE(I)
            IF (NDUP.EQ.0) LEN(I) = IQ(I)
         ENDIF
      ENDDO
#if (SP_DEBUG >= 2)
      WRITE(*,*) 'NDUP =',NDUP
#endif
      IF (NDUP.NE.0) THEN
         IWFR = 1
         DO I=1,NCMP
            K1 = IPE(I) 
            IF (IQ(I).EQ.0) THEN
               LEN(I) = 0
               IPE(I) = IWFR
               CYCLE
            ENDIF
            K2 = K1 + IQ(I) - 1
            L = IWFR
            IPE(I) = IWFR
            DO K=K1,K2
               IF (IW(K).NE.0) THEN
                  IW(IWFR) = IW(K)
                  IWFR     = IWFR + 1
               ENDIF
            ENDDO
            LEN(I) = IWFR - L 
         ENDDO
      ENDIF
      IPE(NCMP+1) = IPE(NCMP) + LEN(NCMP)
      IWFR = IPE(NCMP+1)
#if(SP_DEBUG >= 5)
      WRITE(*,*) 'PIV',PIV
      WRITE(*,*) 'ICMP',ICMP
      WRITE(*,*) 'IPE',IPE(1:NCMP+1)
      WRITE(*,*) 'LEN',LEN(1:NCMP)
      WRITE(*,*) 'IW',IW(1:IWFR-1)
#endif
      RETURN
      END SUBROUTINE ZMUMPS_547
      SUBROUTINE ZMUMPS_551(
     *     N, NE, IP, IRN, SCALING,LSC,CPERM, DIAG,
     *     ICNTL, WEIGHT,MARKED,FLAG,
     *     PIV_OUT, INFO)
      IMPLICIT NONE
      INTEGER N, NE, ICNTL(10), INFO(10),LSC,LWEIGHT
      INTEGER CPERM(N),PIV_OUT(N),IP(N+1),IRN(NE), DIAG(N)
      DOUBLE PRECISION SCALING(LSC),WEIGHT(N+2)
      INTEGER MARKED(N),FLAG(N)
      INTEGER NUM1,NUM2,NUMTOT,PATH_LENGTH,NLAST
      INTEGER I,J,K,L,CUR_EL,CUR_EL_PATH,CUR_EL_PATH_NEXT,BEST_BEG
      INTEGER L1,L2,PTR_SET1,PTR_SET2,TUP,T22,INTER,MERGE
      DOUBLE PRECISION BEST_SCORE,CUR_VAL,TMP,VAL
      DOUBLE PRECISION INITSCORE, ZMUMPS_UPDATESCORE, 
     *     ZMUMPS_UPDATE_INVERSE, ZMUMPS_METRIC2x2
#if(SP_DEBUG >= 1)
      INTERFACE
         FUNCTION ZMUMPS_UPDATESCORE(A,B,T)
         DOUBLE PRECISION A,B
         INTEGER T
         END FUNCTION
         FUNCTION ZMUMPS_UPDATE_INVERSE(A,B,T)
         DOUBLE PRECISION A,B
         INTEGER T
         END FUNCTION ZMUMPS_UPDATE_INVERSE
         FUNCTION ZMUMPS_METRIC2x2(CUR_EL,CUR_EL_PATH,
     *        SET1,SET2,L1,L2,VAL,DIAG,N,FLAG,FLAGON,T)
         INTEGER CUR_EL,CUR_EL_PATH,L1,L2,N
         INTEGER SET1(L1),SET2(L2),DIAG(N),FLAG(N)
         DOUBLE PRECISION VAL
         LOGICAL FLAGON
         INTEGER T
         END FUNCTION ZMUMPS_METRIC2x2
      END INTERFACE
#endif
      LOGICAL VRAI,FAUX,MAX_CARD_DIAG,USE_SCALING
      INTEGER JOB_DEF,DONE
      INTEGER SUM,PROD,STRUCT,MA47,MAGNITUDE
      INTEGER ZERO,ONE
      PARAMETER (JOB_DEF = 5, DONE = -1,
     *     SUM = 1, PROD = 2, STRUCT=1, MA47=2, MAGNITUDE=3,
     *     VRAI = .TRUE., FAUX = .FALSE.)
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      MAX_CARD_DIAG = .TRUE.
      NUM1 = 0
      NUM2 = 0
      NUMTOT = 0
      NLAST = N
      INFO = 0      
      MARKED = 1
      FLAG = 0
      VAL = ONE
      IF(LSC .GT. 1) THEN
         USE_SCALING = .TRUE.
      ELSE
         USE_SCALING = .FALSE.
      ENDIF
      TUP = ICNTL(2)
      IF(TUP .EQ. SUM) THEN
        INITSCORE = ZERO
      ELSE
        INITSCORE = ONE
      ENDIF
      IF(ICNTL(2) .GT. 2 .OR. ICNTL(2) .LE. 0) THEN
         WRITE(*,*)
     *        'ERROR: WRONG VALUE FOR ICNTL(2) = ',ICNTL(2)
         INFO(1) = -1
         RETURN
      ENDIF
      T22 = ICNTL(1)
      IF(ICNTL(1) .LT. 0 .OR. ICNTL(1) .GT. 2) THEN
         WRITE(*,*)
     *        'ERROR: WRONG VALUE FOR ICNTL(1) = ',ICNTL(1)
         INFO(1) = -1
         RETURN
      ENDIF
      DO CUR_EL=1,N
         IF(MARKED(CUR_EL) .LE. 0) THEN
            CYCLE
         ENDIF
         IF(CPERM(CUR_EL) .LT. 0) THEN
            MARKED(CUR_EL) = -1
            CYCLE
         ENDIF
         PATH_LENGTH = 2
         CUR_EL_PATH = CPERM(CUR_EL)
         IF(CUR_EL_PATH .EQ. CUR_EL) THEN
            MARKED(CUR_EL) = -1
            CYCLE
         ENDIF
         MARKED(CUR_EL) = 0
         WEIGHT(1) = INITSCORE
         WEIGHT(2) = INITSCORE
         L1 = IP(CUR_EL+1)-IP(CUR_EL)
         L2 = IP(CUR_EL_PATH+1)-IP(CUR_EL_PATH)
         PTR_SET1 = IP(CUR_EL)
         PTR_SET2 = IP(CUR_EL_PATH)
         IF(USE_SCALING) THEN
            VAL = -SCALING(CUR_EL_PATH) - SCALING(CUR_EL+N)
         ENDIF
         CUR_VAL = ZMUMPS_METRIC2x2(
     *        CUR_EL,CUR_EL_PATH,
     *        IRN(PTR_SET1),IRN(PTR_SET2),
     *        L1,L2,
     *        VAL,DIAG,N,FLAG,FAUX,T22)
         WEIGHT(PATH_LENGTH+1) = 
     *        ZMUMPS_UPDATESCORE(WEIGHT(1),CUR_VAL,TUP)
         DO
            IF(CUR_EL_PATH .EQ. CUR_EL) EXIT
            PATH_LENGTH = PATH_LENGTH+1
            MARKED(CUR_EL_PATH) = 0
            CUR_EL_PATH_NEXT = CPERM(CUR_EL_PATH)
            L1 = IP(CUR_EL_PATH+1)-IP(CUR_EL_PATH)
            L2 = IP(CUR_EL_PATH_NEXT+1)-IP(CUR_EL_PATH_NEXT)
            PTR_SET1 = IP(CUR_EL_PATH)
            PTR_SET2 = IP(CUR_EL_PATH_NEXT)
            IF(USE_SCALING) THEN
               VAL = -SCALING(CUR_EL_PATH_NEXT) 
     *              - SCALING(CUR_EL_PATH+N)
            ENDIF
            CUR_VAL = ZMUMPS_METRIC2x2(
     *           CUR_EL_PATH,CUR_EL_PATH_NEXT,
     *           IRN(PTR_SET1),IRN(PTR_SET2),
     *           L1,L2,
     *           VAL,DIAG,N,FLAG,VRAI,T22)
            WEIGHT(PATH_LENGTH+1) = 
     *           ZMUMPS_UPDATESCORE(WEIGHT(PATH_LENGTH-1),CUR_VAL,TUP)
            CUR_EL_PATH = CUR_EL_PATH_NEXT
         ENDDO
         IF(MOD(PATH_LENGTH,2) .EQ. 1) THEN
            IF(WEIGHT(PATH_LENGTH+1) .GE. WEIGHT(PATH_LENGTH)) THEN
               CUR_EL_PATH = CPERM(CUR_EL)
            ELSE
               CUR_EL_PATH = CUR_EL
            ENDIF
            DO I=1,(PATH_LENGTH-1)/2
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
            NUMTOT = NUMTOT + PATH_LENGTH - 1
         ELSE
            IF(MAX_CARD_DIAG) THEN
               CUR_EL_PATH = CPERM(CUR_EL)
               IF(DIAG(CUR_EL) .NE. 0) THEN 
                  BEST_BEG = CUR_EL_PATH
                  GOTO 1000
               ENDIF
               DO I=1,(PATH_LENGTH/2)
                  CUR_EL_PATH_NEXT = CPERM(CUR_EL_PATH)
                  IF(DIAG(CUR_EL_PATH) .NE. 0) THEN 
                     BEST_BEG = CUR_EL_PATH_NEXT
                     GOTO 1000
                  ENDIF
               ENDDO
            ENDIF
            BEST_BEG = CUR_EL
            BEST_SCORE = WEIGHT(PATH_LENGTH-1)
            CUR_EL_PATH = CPERM(CUR_EL)
            DO I=1,(PATH_LENGTH/2)-1
               TMP = ZMUMPS_UPDATESCORE(WEIGHT(PATH_LENGTH),
     *              WEIGHT(2*I-1),TUP)
               TMP = ZMUMPS_UPDATE_INVERSE(TMP,WEIGHT(2*I),TUP)
               IF(TMP .GT. BEST_SCORE) THEN
                  BEST_SCORE = TMP
                  BEST_BEG = CUR_EL_PATH
               ENDIF
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               TMP = ZMUMPS_UPDATESCORE(WEIGHT(PATH_LENGTH+1),
     *              WEIGHT(2*I),TUP)
               TMP = ZMUMPS_UPDATE_INVERSE(TMP,WEIGHT(2*I+1),TUP)
               IF(TMP .GT. BEST_SCORE) THEN
                  BEST_SCORE = TMP
                  BEST_BEG = CUR_EL_PATH
               ENDIF
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
 1000       CUR_EL_PATH = BEST_BEG
            DO I=1,(PATH_LENGTH/2)-1
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
            NUMTOT = NUMTOT + PATH_LENGTH - 2
            MARKED(CUR_EL_PATH) = -1
         ENDIF
      ENDDO
      DO I=1,N
         IF(MARKED(I) .LT. 0) THEN
            IF(DIAG(I) .EQ. 0) THEN
               PIV_OUT(NLAST) = I
               NLAST = NLAST - 1
            ELSE
               NUM1 = NUM1 + 1
               PIV_OUT(NUM2+NUM1) = I
               NUMTOT = NUMTOT + 1
            ENDIF
         ENDIF
      ENDDO
      INFO(2) = NUMTOT
      INFO(3) = NUM1
      INFO(4) = NUM2
      RETURN
      END SUBROUTINE ZMUMPS_551
      FUNCTION ZMUMPS_UPDATESCORE(A,B,T)
      IMPLICIT NONE
      DOUBLE PRECISION ZMUMPS_UPDATESCORE
      DOUBLE PRECISION A,B
      INTEGER T
      INTEGER SUM,PROD
      PARAMETER(SUM = 1,PROD = 2)
      IF(T .EQ. SUM) THEN
         ZMUMPS_UPDATESCORE = A+B
      ELSE
         ZMUMPS_UPDATESCORE = A*B
      ENDIF
      END FUNCTION ZMUMPS_UPDATESCORE
      FUNCTION ZMUMPS_UPDATE_INVERSE(A,B,T)
      IMPLICIT NONE
      DOUBLE PRECISION ZMUMPS_UPDATE_INVERSE
      DOUBLE PRECISION A,B
      INTEGER T
      INTEGER SUM,PROD
      PARAMETER(SUM = 1,PROD = 2)
      IF(T .EQ. SUM) THEN
         ZMUMPS_UPDATE_INVERSE = A-B
      ELSE
         ZMUMPS_UPDATE_INVERSE = A/B
      ENDIF
      END FUNCTION ZMUMPS_UPDATE_INVERSE
      FUNCTION ZMUMPS_METRIC2x2(CUR_EL,CUR_EL_PATH,
     *     SET1,SET2,L1,L2,VAL,DIAG,N,FLAG,FLAGON,T)
      IMPLICIT NONE
      DOUBLE PRECISION ZMUMPS_METRIC2x2
      INTEGER CUR_EL,CUR_EL_PATH,L1,L2,N
      INTEGER SET1(L1),SET2(L2),DIAG(N),FLAG(N)
      DOUBLE PRECISION VAL
      LOGICAL FLAGON
      INTEGER T
      INTEGER I,INTER,MERGE
      INTEGER STRUCT,MA47,MAGNITUDE
      PARAMETER(STRUCT=0,MA47=1,MAGNITUDE=2)
      IF(T .EQ. STRUCT) THEN
         IF(.NOT. FLAGON) THEN
            DO I=1,L1
               FLAG(SET1(I)) = CUR_EL
            ENDDO            
         ENDIF
         INTER = 0
         DO I=1,L2
            IF(FLAG(SET2(I)) .EQ. CUR_EL) THEN
               INTER = INTER + 1
               FLAG(SET2(I)) = CUR_EL_PATH
            ENDIF
         ENDDO
         MERGE = L1 + L2 - INTER
         ZMUMPS_METRIC2x2 = DBLE(INTER) / DBLE(MERGE)
      ELSE IF (T .EQ. MA47) THEN
         MERGE = 3
         IF(DIAG(CUR_EL) .NE. 0) MERGE = 2
         IF(DIAG(CUR_EL_PATH) .NE. 0) MERGE = MERGE - 2
         IF(MERGE .EQ. 0) THEN
            ZMUMPS_METRIC2x2 = L1+L2-2
            ZMUMPS_METRIC2x2 = -(ZMUMPS_METRIC2x2**2)/2
         ELSE IF(MERGE .EQ. 1) THEN
            ZMUMPS_METRIC2x2 = - DBLE(L1+L2-4) * DBLE(L1-2)
         ELSE IF(MERGE .EQ. 2) THEN
            ZMUMPS_METRIC2x2 = - DBLE(L1+L2-4) * DBLE(L2-2)
         ELSE
            ZMUMPS_METRIC2x2 = - DBLE(L1-2) * DBLE(L2-2)
         ENDIF
      ELSE
         ZMUMPS_METRIC2x2 = VAL
      ENDIF
      END FUNCTION 
      SUBROUTINE ZMUMPS_622(NA, NCMP,
     *      PERM,IPERM, 
     *      LISTVAR_SCHUR, SIZE_SCHUR, AOTOA)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: SIZE_SCHUR, LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, INTENT(IN):: NA, NCMP
      INTEGER, INTENT(IN):: AOTOA(NCMP), IPERM(NCMP)
      INTEGER, INTENT(OUT):: PERM(NA) ! not used a input
      INTEGER CMP_POS, IO, I, K, IPOS
      DO CMP_POS=1, NCMP
        IO              = IPERM(CMP_POS)
        PERM(AOTOA(IO)) = CMP_POS
      ENDDO
      IPOS = NCMP
      DO K =1,  SIZE_SCHUR
        I       = LISTVAR_SCHUR(K)
        IPOS    = IPOS+1
        PERM(I) = IPOS
      ENDDO
      RETURN
      END SUBROUTINE ZMUMPS_622
      SUBROUTINE ZMUMPS_623
     * (NA,N,NZ, IRN, ICN, IW, LW, IPE, LEN,
     * IQ, FLAG, IWFR,
     * NRORM, NIORM, IFLAG,IERROR, ICNTL, 
     * symmetry, SYM, MedDens, NBQD, AvgDens, 
     * LISTVAR_SCHUR, SIZE_SCHUR, ATOAO, AOTOA)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: NA,N,NZ,LW
      INTEGER, INTENT(IN)  :: SIZE_SCHUR, LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, INTENT(IN)  :: IRN(NZ), ICN(NZ) 
      INTEGER, INTENT(IN)  :: ICNTL(40), SYM
      INTEGER, INTENT(OUT) :: IFLAG,IERROR,NRORM,NIORM,IWFR
      INTEGER, INTENT(OUT) :: AOTOA(N)
      INTEGER, INTENT(OUT) :: ATOAO(NA)  ! only used internallly
      INTEGER, INTENT(OUT) :: LEN(N), IPE(N+1)
      INTEGER, INTENT(OUT) :: symmetry, 
     &                        MedDens, NBQD, AvgDens
      INTEGER, INTENT(OUT)  :: FLAG(N), IW(LW), IQ(N)
      INTEGER MP, MPG, NAO
      INTEGER I,K,J,N1,LAST,NDUP,K1,K2,L
      INTEGER NBERR, THRESH, IAO
      INTEGER NZOFFA, NDIAGA
      DOUBLE PRECISION RSYM
      INTRINSIC NINT
      ATOAO(1:NA) = 0
      DO I = 1, SIZE_SCHUR
        ATOAO(LISTVAR_SCHUR(I)) = -1
      ENDDO
      IAO = 0  
      DO I= 1, NA
        IF (ATOAO(I).LT.0) CYCLE
        IAO = IAO +1   
        ATOAO(I)   = IAO
        AOTOA(IAO) = I
      ENDDO
      MP = ICNTL(2)
      MPG= ICNTL(3)
      NIORM  = 3*N
      NDIAGA = 0
      IERROR = 0
      IPE(1:N+1) = 0
      DO K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     *                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
        ELSE
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  ! in schur
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1
           IPE(J) = IPE(J) + 1
           NIORM  = NIORM + 1
          ELSE
           NDIAGA = NDIAGA + 1
          ENDIF
        ENDIF
      ENDDO
      NZOFFA  = NIORM - 3*N
      IF (IERROR.GE.1) THEN
         NBERR  = 0
         IFLAG  = IFLAG+1
         IF ((MP.GT.0).AND.(ICNTL(4).GE.2))  THEN 
          WRITE (MP,99999) 
          DO 70 K=1,NZ
           I = IRN(K)
           J = ICN(K)
           IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     *                            .OR.(J.LT.1)) THEN
            NBERR = NBERR + 1
            IF (NBERR.LE.10)  THEN
               IF (MOD(K,10).GT.3 .OR. MOD(K,10).EQ.0 .OR.
     *             (10.LE.K .AND. K.LE.20)) THEN
                 WRITE (MP,'(I8,A,I8,A,I8,A)')
     *             K,'th entry (in row',I,' and column',J,') ignored'
               ELSE
                 IF (MOD(K,10).EQ.1) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'st entry (in row',I,' and column',J,') ignored'
                 IF (MOD(K,10).EQ.2) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'nd entry (in row',I,' and column',J,') ignored'
                 IF (MOD(K,10).EQ.3) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'rd entry (in row',I,' and column',J,') ignored'
               ENDIF
            ELSE
               GO TO 100
            ENDIF
           ENDIF
   70     CONTINUE
         ENDIF
      ENDIF
  100 NRORM = NIORM - 2*N
      IQ(1) = 1
      N1 = N - 1
      IF (N1.GT.0) THEN
        DO 110 I=1,N1
            IQ(I+1) = IPE(I) + IQ(I) 
  110   CONTINUE
      ENDIF
      LAST = MAX(IPE(N)+IQ(N)-1,IQ(N))
      FLAG(1:N) = 0
      IPE(1:N)  = IQ(1:N)
      IW(1:LAST) = 0
      IWFR = LAST + 1
      DO 200 K=1,NZ
         I = IRN(K)
         J = ICN(K)
         IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     *                          .OR.(J.LT.1)) CYCLE
         I = ATOAO(I)
         J = ATOAO(J)
         IF ((I.LT.0).OR.(J.LT.0)) CYCLE  ! in schur
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1 
          ELSE
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1
          ENDIF
         ENDIF
  200 CONTINUE
      NDUP = 0
      DO 260 I=1,N
        K1 = IPE(I) 
        K2 = IQ(I) -1
        IF (K1.GT.K2) THEN
         LEN(I) = 0
         IQ(I)  = 0
        ELSE
         DO 240 K=K1,K2
           J     = -IW(K)
           IF (J.LE.0) GO TO 250
           L     = IQ(J) 
           IQ(J) = L + 1
           IF (FLAG(J).EQ.I) THEN
            NDUP = NDUP + 1
            IW(L) = 0
            IW(K) = 0
           ELSE
            IW(L)   = I
            IW(K)   = J
            FLAG(J) = I
           ENDIF
  240    CONTINUE
  250    IQ(I) = IQ(I) - IPE(I)
         IF (NDUP.EQ.0) LEN(I) = IQ(I)
        ENDIF
  260 CONTINUE
      IF (NDUP.NE.0) THEN
       IWFR = 1
       DO 280 I=1,N
         IF (IQ(I).EQ.0) THEN
             LEN(I) = 0
            IPE(I) = IWFR
            GOTO 280
         ENDIF
         K1 = IPE(I) 
         K2 = K1 + IQ(I) - 1
         L = IWFR
         IPE(I) = IWFR
         DO 270 K=K1,K2
           IF (IW(K).NE.0) THEN
            IW(IWFR) = IW(K)
            IWFR     = IWFR + 1
           ENDIF
  270    CONTINUE
         LEN(I) = IWFR - L 
  280  CONTINUE
      ENDIF
      IPE(N+1) = IPE(N) + LEN(N)
      IWFR = IPE(N+1)
      IF (SYM.EQ.0) THEN
      RSYM =  DBLE(NDIAGA+2*NZOFFA - (IWFR-1))/
     *            DBLE(NZOFFA+NDIAGA) 
      symmetry = NINT (100.0*RSYM)
         IF (MPG .GT. 0)
     &  write(MPG,'(A,I5)') 
     &  ' ... Structural symmetry (in percent)=', symmetry
        IF (MP.GT.0 .AND. MPG.NE.MP)
     &  write(MP,'(A,I5)') 
     &  ' ... Structural symmetry (in percent)=', symmetry
      ELSE
       symmetry = 100
      ENDIF
      AvgDens = NINT(DBLE(IWFR-1)/DBLE(N))
      THRESH  = AvgDens*50 - AvgDens/10 + 1
      NBQD    = 0
      IF (N.GT.2) THEN
        IQ(1:N) = 0
        DO I= 1, N
          K = MAX(LEN(I),1)
          IQ(K) = IQ(K) + 1
          IF (K.GT.THRESH) NBQD = NBQD+1
        ENDDO
        K = 0
        MedDens = 0
        DO WHILE (K .LT. (N/2))
         MedDens = MedDens + 1
         K       = K+IQ(MedDens)
        ENDDO
      ELSE
        MedDens = AvgDens
      ENDIF
         IF (MPG .GT. 0)
     &  write(MPG,'(A,3I5)') 
     &  ' Density: NBdense, Average, Median   =',
     &  NBQD, AvgDens, MedDens
        IF (MP.GT.0 .AND. MPG.NE.MP)
     &  write(MP,'(A,3I5)') 
     &  ' Density: NBdense, Average, Median   =',
     &  NBQD, AvgDens, MedDens
      RETURN
99999 FORMAT (/'*** Warning message from analysis routine ***')
      END SUBROUTINE ZMUMPS_623
