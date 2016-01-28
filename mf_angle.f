!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate the angle between 3 atoms 

      SUBROUTINE MF_INTL_ANGLE(STN)
!
      USE specify
      USE structure
      USE build 
      USE elements  
!
      Implicit none
!
      INTEGER :: STN,ACNT,FID,SZ,EL,AN,I,NAr,B
     & ,ios
      REAL*8 :: ANGSCL
      CHARACTER(IDSZ) :: Ai,Bi,Ci,Ao,Bo,Co
     & ,AT,GT,RT,SM
      CHARACTER(IDSZ) ,ALLOCATABLE :: IDb(:) 
      CHARACTER(CHSZ) :: ESTAT
      LOGICAL :: ADD 
!     
!     Set the local variables
!
      ACNT = ATOM_ANGL_CNT(STN)               ! number of angles to calculate 
      NAr = NA(STN) 
      FID  = 80                       ! Initial file number 
      ANGSCL = 1.0d0
      SZ = INT(180.0d0*ANGSCL) + 2   ! 0 - 181
!
!     Verbose output
!
      IF (VERB ) THEN
        WRITE(*,*) "Calc atomic angles "
        WRITE(*,*) "  Structure ",STN," with ", ACNT," angles"
      ENDIF
!
!     Allocate arrays 
!
      ALLOCATE( ATOM_THSUM(ACNT),ATOM_THDIS(0:SZ,ACNT) )
      ATOM_THSUM(ACNT) = 0.0d0
      ATOM_THDIS(:,:) = 0
      ALLOCATE( IDb(ACNT),BCNT(ACNT), ATOMB_NB(NAr,ACNT) )

      DO AN=1,ACNT
        IDb(AN) =  ATOM_ANAX(2,AN,STN)
      ENDDO  
!
!     Intialize
!
      BCNT(:) = 0 
      ATOM_THSUM(:) = 0
      ATOM_THDIS(:,:) = 0
!
!     Find central atom types and record atom #
!
      DO I = 1,NA(STN)
          EL =  ELN(I,STN) 
          AT =  TYP(I,STN)
          GT =  GRID(I,STN)
          RT =  RID(I,STN)
          SM = ATSYM(EL) 
          DO AN=1,ACNT 
           Bo = IDb(AN) 
           ADD = .FALSE.
           IF( AT .EQ.Bo ) ADD=.TRUE.
           IF( GT .EQ.Bo ) ADD=.TRUE.
           IF( RT .EQ.Bo ) ADD=.TRUE.
           IF( SM .EQ.Bo ) ADD=.TRUE.
           IF( ADD ) THEN
             B = BCNT(AN)
             B = B + 1
             ATOMB_NB(B,AN) = I
             BCNT(AN) = B
           ENDIF
          ENDDO
      ENDDO 
      DEALLOCATE( IDb )!
!
      RETURN 
      END SUBROUTINE MF_INTL_ANGLE 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate the angle between 3 atoms 

      SUBROUTINE MF_CALC_ANGLE(STN)
!
      USE specify
      USE structure
      USE build
      USE elements  
!
      Implicit none
!
      INTEGER :: STN,I,J,ACNT,Ai,Bi,Ci,C,D,AN,EL,B
     & ,ios,NBi,NBf,NB,NAi
      REAL*8, DIMENSION(NDIM) :: DF_AB,DF_BC
      REAL*8 :: SQ_AB,SQ_BC,SQ_ABBC,MAG_ABBC,MAG_AB,MAG_BC
     & ,COSTH,THETA,THETAd,LVr(NDIM,NDIM)
      CHARACTER(IDSZ) ::  Ao,Co,AT,GT,RT,SM
      LOGICAL :: ADD 
       CHARACTER(CHSZ) :: ESTAT
!     
!     Set the local variables
!
      ACNT = ATOM_ANGL_CNT(STN)               ! number of angles to calculate 
 !
!     Triclinc PBC's
      LVr(:,:) = LV(:,:,STN)
      NAi = NA(STN)
      ALLOCATE( Fr(NDIM,NAi) ) 
!
!     Loop over all atomic angles 
!
      DO AN=1,ACNT 
!
!       Asign local id's for side atoms 
!
        Ao = ATOM_ANAX(1,AN,STN)
        Co = ATOM_ANAX(3,AN,STN)
!
!       Loop over central atoms
!
        DO B=1,BCNT(AN)
!
!        For centeral atom find side atoms 
! 
         Bi = ATOMB_NB(B,AN)
!
!        Loop over nieghbors 
!
         AI = 0
         Ci = 0
         NBi = NINDX(Bi,STN)
         NBf = NINDX(Bi+1,STN) -1


         IF(VERB ) WRITE(*,*) 'Looping over ',Bi,'s'
     &  ,NBf-NBi+1,' nieghbors '


         DO NB = NBi,NBf
          J = NLIST(NB,STN)
          EL = ELN(J,STN) 
          AT = TYP(J,STN)
          GT = GRID(J,STN)
          RT = RID(J,STN)
          SM = ATSYM(EL)

! debug
          WRITE(*,*) Bi,J,EL,AT,GT
 
          ADD = .FALSE.
          IF( AT .EQ.Ao ) ADD=.TRUE.
          IF( GT .EQ.Ao ) ADD=.TRUE.
          IF( RT .EQ.Ao ) ADD=.TRUE.
          IF( SM .EQ.Ao ) ADD=.TRUE.
          IF( ADD ) Ai  = J
          ADD = .FALSE.
          IF( AT .EQ.Co ) ADD=.TRUE.
          IF( GT .EQ.Co ) ADD=.TRUE.
          IF( RT .EQ.Co ) ADD=.TRUE.
          IF( SM .EQ.Co ) ADD=.TRUE.
          IF( ADD ) Ci  = J

         ENDDO 
         IF ( Ai.EQ.0 .OR. Ci.EQ.0 ) THEN
           WRITE(ESTAT,*) 'For ',Bi,' can not find nb ',Ao,Co
           CALL PRNERROR(-12,ESTAT) 
         ENDIF 
!
!        Calculate angle 
!
         SQ_AB = 0
         SQ_BC = 0
         SQ_ABBC = 0
         DO D=1,NDIM
          DF_AB(D) = Fr(D,Ai) - Fr(D,Bi)
          DF_BC(D) = Fr(D,Ci) - Fr(D,Bi)
!
!         Apply PBC's
!
          IF( PBCS(D,STN) ) THEN
            DF_AB(D) = DF_AB(D) - ANINT(DF_AB(D))
            DF_BC(D) = DF_BC(D) - ANINT(DF_BC(D))
          ENDIF
          SQ_AB = SQ_AB + DF_AB(D)*DF_AB(D)
          SQ_BC = SQ_BC + DF_BC(D)*DF_BC(D)
          SQ_ABBC = SQ_ABBC + DF_AB(D)*DF_BC(D)
         ENDDO
         MAG_AB = SQRT( SQ_AB )
         MAG_BC = SQRT( SQ_BC )
         COSTH = SQ_ABBC/MAG_AB/MAG_BC

         THETA = ACOS(COSTH)
         THETAd = THETA*360.0d0/(2.0d0*PI)      

         IF(VERB) THEN
           WRITE(*,*) Ai,Bi,Ci,THETAd  
         ENDIF

         ATOM_THSUM(AN) = ATOM_THSUM(AN)  + THETAd
         C =  ANINT( THETAd )
         ATOM_THDIS(C,AN) = ATOM_THDIS(C,AN) + 1
        ENDDO ! BCNT
!
      ENDDO 
!
      DEALLOCATE( Fr ) 
      RETURN
      END SUBROUTINE MF_CALC_ANGLE 

!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate the angle between 3 atoms 

      SUBROUTINE MF_FIN_ANGLE(STN)
!
      USE specify
      USE structure
      USE build 
!
      Implicit none
!
      INTEGER :: STN,I,ACNT,FID,AN,F,SZ
     & ,NPNTS,C,CNT
      CHARACTER(IDSZ) :: Ao,Bo,Co
      REAL*8 :: ANGSCL,AVEC,THETAD,TCNTS,CNTf,PROB 
 
!     
!     Set the local variables
!
      FID  = 80                       ! Initial file number 
      ANGSCL = 1.0d0
      ACNT = ATOM_ANGL_CNT(STN)
      SZ = INT(180.0d0*ANGSCL) + 2   ! 0 - 181

      DO AN=1,ACNT 
!
!       Open ouput files
!

        IF(VERB) THEN 
          WRITE(*,*) 'Writing angle dist ',AN,' to file '
     &  ,ATOM_ANFILE(AN,STN)
        ENDIF 

       Ao = ATOM_ANAX(1,AN,STN)
       Co = ATOM_ANAX(2,AN,STN)
       Co = ATOM_ANAX(3,AN,STN)

        F = FID + AN 
        OPEN(UNIT=F,FILE=ATOM_ANFILE(AN,STN),STATUS='UNKNOWN')
        WRITE(F,*) "# Angular distribution"
        WRITE(F,*) "# Scale ", 1/ANGSCL
        WRITE(F,'("# Atom names:",3A4)' ) Ao,Bo,Co
!
!       Calc avwerage 
!
        NPNTS = 0
        DO C=1,SZ 
          NPNTS  = NPNTS + ATOM_THDIS(C,AN)
        ENDDO
!       Calculate average
        TCNTS = FLOAT(NPNTS)
        AVEC = ATOM_THSUM(AN)/TCNTS 
        WRITE(F,*) "# total vectors found", NPNTS      
        WRITE(F,*) "# average cos \theta", AVEC 
!     
        DO C=0,SZ
          CNT = ATOM_THDIS(C,AN)
          CNTf = FLOAT(CNT) 
          THETAd = FLOAT( C )
          PROB = CNTf/TCNTS 
          WRITE(F,3601) THETAd,PROB,CNT
        ENDDO
!
!     Clsoe ouput files
!
        CLOSE(F)
      ENDDO     ! ANCNT         

!
!     deallocate arrays
!
      DEALLOCATE( ATOM_THSUM,ATOM_THDIS,BCNT,ATOMB_NB )

!
      RETURN 
 3601 FORMAT(2F12.6,I6)
      END SUBROUTINE MF_FIN_ANGLE 
