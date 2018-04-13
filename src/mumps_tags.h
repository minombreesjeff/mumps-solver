C
C  This file is part of MUMPS 4.8.3, built on Wed Sep 24 17:11:30 UTC 2008
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
C  from the Web pages:
C  http://mumps.enseeiht.fr/  or  http://graal.ens-lyon.fr/MUMPS
C
C
C   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
C   EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
C
C
C  User documentation of any code that uses this software can
C  include this complete notice. You can acknowledge (using
C  references [1], [2], and [3]) the contribution of this package
C  in any scientific publication dependent upon the use of the
C  package. You shall use reasonable endeavours to notify
C  the authors of the package of this publication.
C
C   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent,
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
      INTEGER ARROWHEAD, ARR_INT, ARR_REAL, ELT_INT, ELT_REAL
      PARAMETER ( ARROWHEAD = 20,
     *            ARR_INT = 29,
     *            ARR_REAL = 30,
     *            ELT_INT = 31,
     *            ELT_REAL = 32 )
      INTEGER COLLECT_NZ, COLLECT_IRN, COLLECT_JCN
      PARAMETER( COLLECT_NZ  = 35,
     *           COLLECT_IRN = 36,
     *           COLLECT_JCN = 37 )
      INTEGER RACINE,
     *        NOEUD,
     *        TERREUR,
     *        MAITRE_DESC_BANDE,
     *        MAITRE2,
     *        BLOC_FACTO,
     *        CONTRIB_TYPE2,
     *        MAPLIG,
     *        FACTOR
      PARAMETER ( RACINE            = 2,
     *            NOEUD             = 3,
     *            MAITRE_DESC_BANDE = 4,
     *            MAITRE2           = 5,
     *            BLOC_FACTO        = 6,
     *            CONTRIB_TYPE2     = 7,
     *            MAPLIG            = 8,
     *            FACTOR            = 9,
     *            TERREUR           = 99 )
      INTEGER ROOT_NELIM_INDICES,
     *        ROOT_CONT_STATIC,
     *        ROOT_NON_ELIM_CB,
     *        ROOT_2SLAVE,
     *        ROOT_2SON
       PARAMETER( ROOT_NELIM_INDICES = 15,
     *        ROOT_CONT_STATIC       = 16,
     *        ROOT_NON_ELIM_CB       = 17,
     *        ROOT_2SLAVE            = 18,
     *        ROOT_2SON              = 19 )
      INTEGER RACINE_SOLVE,
     *        ContVec,
     *        Master2Slave,
     *        GatherSol,
     *        ScatterRhsI,
     *        ScatterRhsR
      PARAMETER( RACINE_SOLVE = 10,
     *           ContVec      = 11,
     *           Master2Slave = 12,
     *           GatherSol    = 13,
     *           ScatterRhsI  = 54,
     *           ScatterRhsR  = 55)
      INTEGER FEUILLE,
     *        BACKSLV_UPDATERHS,
     *        BACKSLV_MASTER2SLAVE
      PARAMETER( FEUILLE = 21,
     *           BACKSLV_UPDATERHS = 22,
     *           BACKSLV_MASTER2SLAVE = 23 )
      INTEGER SYMMETRIZE
      PARAMETER ( SYMMETRIZE = 24 )
      INTEGER BLOC_FACTO_SYM,
     *        BLOC_FACTO_SYM_SLAVE, END_NIV2_LDLT
      PARAMETER ( BLOC_FACTO_SYM = 25,
     *            BLOC_FACTO_SYM_SLAVE = 26 , 
     *            END_NIV2_LDLT = 33 )
      INTEGER UPDATE_LOAD
      PARAMETER ( UPDATE_LOAD = 27 )
      INTEGER DEFIC_TAG
      PARAMETER(  DEFIC_TAG = 28 )
      INTEGER TAG_SCHUR
      PARAMETER( TAG_SCHUR = 38 )
      INTEGER TAG_DUMMY
      PARAMETER( TAG_DUMMY = 39 )
      INTEGER ZERO_PIV
      PARAMETER( ZERO_PIV = 40 )
