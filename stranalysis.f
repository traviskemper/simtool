!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Return atomic number based on atomic symbol
!
      SUBROUTINE getanumb(SYMr,ELNr)
      USE elements
!
      IMPLICIT none
!
      INTEGER :: EL,STN,ELNr
      CHARACTER(ATSZ) :: SYMr
!
      ELNr = 0
      DO EL=0,NELM
        IF( ADJUSTL(SYMr).EQ.ADJUSTL(ATSYM(EL)) ) THEN
           ELNr = EL
          EXIT
        ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE  getanumb
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate number of elements
!
      SUBROUTINE CALCNELM(STN)
      USE specify
      USE structure 
      USE elements 
!
      IMPLICIT none
      LOGICAL :: EFIND
      INTEGER :: STN,I ,EL,ETOT,E,ELr,CNT,N,Ef
      CHARACTER(CHSZ) :: ESTAT
!
      ETOT = 0
!     Initialize element list
      ETOT = 1
      I = 1
      EL = ELN(I,STN)
      NEl(ETOT,STN) = EL
      CNT = 1
      ELCNT(ETOT,STN) = CNT
      ELIST(CNT,ETOT,STN) = I
      ELREF(I,STN) = ETOT
!
      DO I=2,NA(STN)
        EL = ELN(I,STN)
        EFIND=.FALSE.
        DO E=1,ETOT
          ELr = NEL(E,STN)
          IF( ELr.EQ.EL) THEN
            EFIND = .TRUE.
            Ef = E
            EXIT
          ENDIF
        ENDDO
        IF( EFIND ) THEN
            CNT =  ELCNT(Ef,STN)+ 1
            ELCNT(Ef,STN) = CNT
            ELIST(CNT,Ef,STN) = I
            ELREF(I,STN) = Ef
        ELSE
            ETOT = ETOT + 1
            IF( ETOT.GT.NELM) THEN
              WRITE(ESTAT,'(2I4,"Found :",5I4)') EL,Elr,NEL(1:5,STN)
              CALL prnerror(-1,ESTAT)
            ENDIF
            NEl(ETOT,STN) = EL
            CNT = 1
            ELCNT(ETOT,STN) = CNT
            ELIST(CNT,ETOT,STN) = I
            ELREF(I,STN) = ETOT 
        ENDIF    
      ENDDO
      ELT(STN) = ETOT
      IF( VERB) THEN
        WRITE(6,*) "Final structure",STN," has",ELT(STN)," elements"
        DO E=1,ETOT           
          WRITE(*,*) E,NEL(E,STN)," has ",ELCNT(E,STN)," atoms"
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE  CALCNELM
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 10/24/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Calculate lattice constants and angles based on lattice vectors
!
      SUBROUTINE calclc(STN)
      USE structure
      USE specify 
!
      IMPLICIT none
!
      INTEGER :: N,EL,D,STN
      REAL*8 :: LVr(NDIM,NDIM),LVSQ,DT(3),LCr(NDIM)
     &  ,COSTH,THETA,THETAd
!      
      IF(VERB)WRITE(6,*) " Calculating lattice constants ",
     & " based on lattice vectors"
!     Calc lattice constants and angles
!     cos Oij = i*j/( |i| |j| )

      LVr(:,:) = LV(:,:,STN)
!
      DT(:) = 0.0d0 
      DO D = 1,NDIM
        LVSQ=LVr(D,1)*LVr(D,1)+LVr(D,2)*LVr(D,2)+LVr(D,3)*LVr(D,3) 
        LCr(D) = SQRT( LVSQ )
!       a*b
        DT(1) = DT(1) + LVr(1,D)*LVr(2,D)
!       a*c
        DT(2) = DT(2) + LVr(1,D)*LVr(3,D)
!       b*c      
        DT(3) = DT(3) + LVr(2,D)*LVr(3,D)
      ENDDO
!     Gamma  a-b
      COSTH = DT(1)/LCr(1)/LCr(2)
      THETA = ACOS(COSTH)
      THETAd = THETA*360/(2.0d0*PI)
      LA(3,STN) = THETAd
!     Beta a-c
      COSTH = DT(2)/LCr(1)/LCr(3)
      THETA = ACOS(COSTH)
      THETAd = THETA*360/(2.0d0*PI)
      LA(2,STN) = THETAd
!     alpha b-c
      COSTH = DT(3)/LCr(2)/LCr(3)
      THETA = ACOS(COSTH)
      THETAd = THETA*360/(2.0d0*PI)
      LA(1,STN) = THETAd
!     Pass to global
      LC(:,STN) = LCr(:)
!
      END SUBROUTINE calclc
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 4/11/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Calculate lattice constants and angles based on lattice vectors

      SUBROUTINE calclv(STN)
      USE structure
      USE const
      USE specify 
!
      IMPLICIT none
 !
      INTEGER :: N,EL,STN
      REAL*8 :: RCONV,alpha,beta,gamma,A,B,C
     & ,ax,ay,az,bx,by,bz,cx,cy,cz

      IF(VERB)WRITE(6,*) " Calculating lattice vectors ",
     & " based on lattice constants"
      RCONV = 2.0d0*PI/360.0d0
      alpha = LA(1,STN)*RCONV
      beta  = LA(2,STN)*RCONV
      gamma = LA(3,STN)*RCONV
      A = LC(1,STN)
      B = LC(2,STN)
      C = LC(3,STN)
 !
      ax = A
      bx = cos( gamma )* B      
      by = sin( gamma )* B      
      cx = C*cos(beta)
      cy = C*(cos(alpha)- cos(beta)*cos(gamma))/sin(gamma)
      cz = sqrt( C*C - cx*cx - cy*cy )
      IF( ax**2 .LT. ZCUT ) ax = 0.0d0 
      IF( bx**2 .LT. ZCUT ) bx = 0.0d0 
      IF( by**2 .LT. ZCUT ) by = 0.0d0 
      IF( cx**2 .LT. ZCUT ) cx = 0.0d0 
      IF( cy**2 .LT. ZCUT ) cy = 0.0d0 
      IF( cz**2 .LT. ZCUT ) cz = 0.0d0 
!     Set lattice vectors
      LV(:,:,STN) = 0.0d0
      LV(1,1,STN) = ax
      LV(2,1,STN) = bx
      LV(2,2,STN) = by
      LV(3,1,STN) = cx
      LV(3,2,STN) = cy
      LV(3,3,STN) = cz
!     
      IF(VERB) THEN
        WRITE(6,*) "Lattice vectors ",LV(:,:,STN)
     & ,"calculated based on lattice constants",LC(:,STN),LA(:,STN)
      ENDIF
      RETURN
      END SUBROUTINE calclv
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/02/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!    Calculate basic properties 

      SUBROUTINE CALCVD(STN)
      USE structure
      USE const
      USE elements 
      USE potential 
      USE build  
      USE specify 
!
      IMPLICIT none
!
      INTEGER :: I,STN,EL,D
      REAL*8  :: VSQ,VTOT,VAVE,LVr(NDIM,NDIM),A,B,C
!     V = ( b x c) * a
      LVr(:,:) = LV(:,:,STN)
      A = LVr(2,2)*LVr(3,3) - LVr(2,3)*LVr(3,2)
      B = LVr(2,1)*LVr(3,3) - LVr(2,3)*LVr(3,1)
      C = LVr(2,1)*LVr(3,2) - LVr(2,2)*LVr(3,1)
      VOL(STN) =  LVr(1,1)*A - LVr(1,2)*B + LVr(1,3)*C
!     V = ( a x b ) * c
      A = LVr(1,2)*LVr(2,3) - LVr(1,3)*LVr(2,2)
      B = LVr(1,1)*LVr(2,3) - LVr(1,3)*LVr(2,1)
      C = LVr(1,1)*LVr(2,2) - LVr(1,2)*LVr(2,1)
      VOL(STN) =  LVr(3,1)*A - LVr(3,2)*B + LVr(3,3)*C

! Total mass
      TMAS(STN) = 0.0d0
      TCHG(STN) = 0.0d0 
      VTOT = 0.0d0
      DO I=1,NA(STN)
        TMAS(STN) = TMAS(STN) + AMAS(I,STN)
        TCHG(STN) = TCHG(STN) + ACHG(I,STN) 
        VSQ = 0.0d0
        DO D=1,NDIM
          VSQ = VSQ + R1(D,I,STN)
        ENDDO
        VTOT = VTOT + SQRT(VSQ)
!        KE = 
      ENDDO  
      VAVE = VTOT/NA(STN)
!
      DEN(STN) = TMAS(STN)/VOL(STN)/AVO*10.0d0
!     Maxmin density
!      IF( 
      IF( MKSLAB(STN).OR.MKVAC(1)) THEN 
        CALL MAXMIN(STN)
        SHT(STN) = MMN(2,SDIM,STN) - MMN(1,SDIM,STN)
        LVr(SDIM,:) =0.0d0
        LVr(SDIM,SDIM) = SHT(STN)
        A = LVr(1,2)*LVr(2,3) - LVr(1,3)*LVr(2,2)
        B = LVr(1,1)*LVr(2,3) - LVr(1,3)*LVr(2,1)
        C = LVr(1,1)*LVr(2,2) - LVr(1,2)*LVr(2,1)
        SVOL(STN) = LVr(SDIM,1)*A-LVr(SDIM,2)*B+LVr(SDIM,3)*C
        SDEN(STN) = TMAS(STN)/SVOL(STN)/AVO*10.0d0
      ENDIF
!     P = V/(nkT)
!      NpernmK = GPRES*MPatoJ/BOLZ
!      NGAS  = INT( NpernmK*VVOL*Atonm*Atonm*Atonm/GASTEMP )
 
      RETURN
      END SUBROUTINE CALCVD 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 4.0 04/27/2013 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!    Calculate volume
!
      SUBROUTINE CALCMAS(STN)
      USE structure
      USE const
      USE elements 
      USE potential 
!
      IMPLICIT none
!
      INTEGER :: I,STN,EL
      REAL*8  :: MS
!
      DO I=1,NA(STN) 
        EL = ELN( I,STN )
        MS = AMASS( EL )
        AMAS(I,STN) = MS        
      ENDDO 

      RETURN
      END SUBROUTINE CALCMAS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Return atomic number based on atomic symbol
!
      SUBROUTINE FINDELN(STN)
      USE elements
      USE structure 
      USE potential 
      USE specify  
!
      IMPLICIT none
!
      INTEGER :: STN,I,ELNr
      REAL*8 :: MASr
!
      DO I =1, NA(STN)
        MASr = AMAS(I,STN)
        ELNr = -1
        CALL GETMASAN(MASr,ELNr)
        IF( VERB.AND.ELNr.EQ.-1 ) THEN
            WRITE(*,*) " Atom ",I," of str ",STN
            WRITE(*,*) " with mass of ",MASr
            WRITE(*,*) " could not be found"
            ELNr = 0 
        ENDIF
        ELN(I,STN) = ELNr
      ENDDO
!
      RETURN
      END SUBROUTINE  FINDELN 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Return atomic number based on atomic symbol
!
      SUBROUTINE GETMASAN(MASr,ELNr)
      USE elements
!
      IMPLICIT none
!
      INTEGER :: ELNr,MRi,MRr,EL
      REAL*8 :: MASr,MASi,CPSC
!
      CPSC = 1.0d0 
      MRr = INT( MASr*CPSC )
      DO EL=0,NELM
        MASi = AMASS(EL) 
        MRi = INT( MASi*CPSC ) 
        IF( MRr.EQ.MRi ) THEN
           ELNr = EL
           EXIT
        ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE  GETMASAN
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate the angle between a vector and an axis

      SUBROUTINE ANGAXM(STN)
      USE specify
      USE structure
      USE build 
!
      Implicit none
!
      INTEGER :: STN,I,ELi,ELj,D,C,SZ,CNT,NPNTS,J,Jo,Jf,Nj,Ni,NAr
     & ,MTOT,MNi,Mi,Mf,ACNT,FID,F,AN ,R,RRANG,RSHFT ,TACNT 
      REAL*8 :: X,Y,Z,AA,AB,AC,BA,BB,BC,CA,CB,CC,Fi(NDIM),Fj(NDIM)
     & ,MLV(NDIM),DF(NDIM),DR(NDIM),RSQ
     & ,RAi(NDIM),RCENT(NDIM),LVr(NDIM,NDIM),HDR(NDIM)
     & ,MGR,ADB,COSTH,ANGSCL,AVEC,THETA,THETAd
     & ,RZ,TASUM,TAAVE 
      INTEGER, ALLOCATABLE :: THDIS(:,:),DPTHCNT(:,:)
      REAL*8, ALLOCATABLE :: RAVE(:,:,:),THSUM(:),REFV(:,:),MREF(:)
     & ,DPTHAVE(:,:),PROF(:) 
      CHARACTER(CHSZ) :: ESTAT,FL
      CHARACTER(IDSZ) :: Ai,Bi,Ao,Bo
      CHARACTER(16*11) :: TAOUT 
!     
!     Set the local variables
      ACNT = ANCNT(STN)               ! number of angles to calculate 
      FID  = 50                       ! Initial file number 
      MTOT = MOLCNT(STN)
!     Verbose output
      IF (VERB ) THEN
        WRITE(*,*) "Starting angle calc within molecules"
        WRITE(*,*) "  Structure ",STN," with ", MTOT," molecules"
      ENDIF
!     Set scale for bins
        CALL MAXMIN(STN)
        RSHFT = ANINT( MMN(1,SDIM,STN)*RSCALE(STN) ) - RSCALE(STN) 
        RRANG = ANINT( MMN(2,SDIM,STN)*RSCALE(STN) ) -  RSHFT 

       ANGSCL = 1.0d0
!      SZ = INT(2.0d0*ANGSCL) + 1   !-1 to 1 for cos
      SZ = INT(180.0d0*ANGSCL) + 2   ! 0 - 181
!
!     Allocate count arrays
! 
      ALLOCATE( THDIS(0:SZ,ACNT),RAVE(0:SZ,NDIM,ACNT)
     & ,REFV(NDIM,ACNT),THSUM(ACNT),MREF(ACNT)
     & ,DPTHCNT(RRANG,ACNT),DPTHAVE(RRANG,ACNT) 
     & ,PROF(ACNT) )

      THDIS(:,:) = 0.0d0
      THSUM(:) = 0.0d0
      RAVE(:,:,:) = 0.0d0

      DPTHCNT(:,:) =  0
      DPTHAVE(:,:) =  0.0d0 
         
!     Open ouput files
      DO AN=1,ACNT 
        F = FID + AN 
        Ao = ANAX(1,AN,STN)
        Bo = ANAX(2,AN,STN)
        REFV(:,AN) = ANGX(:,AN,STN)
        MREF(AN) = 0.0d0
        DO D =1,NDIM
          MREF(AN) = MREF(AN) + REFV(D,AN)*REFV(D,AN)
        ENDDO
        OPEN(UNIT=F,FILE=ANGAXFL(AN,STN),STATUS='UNKNOWN')
        WRITE(F,*) "# Angular distribution"
        WRITE(F,*) "# Scale ", 1/ANGSCL
        WRITE(F,'("# Atom names:",2A4)' ) Ao,Bo
        WRITE(F,'("#   vector :",3F16.3)' ) REFV(:,AN)
      ENDDO 
      F = FID + ACNT + 1
!     WRITE(FL,'Prof-',A) = ANGAXFL(AN,STN)
      OPEN(UNIT=F,FILE=DPANGAXFL(STN),STATUS='UNKNOWN')
      WRITE(F,*) "# Depth distribution of angles in",SDIM
      WRITE(F,*) "# Scale ", RSCALE(STN),RSHFT
      WRITE(F,'("# Atom names:",2A4)' ) Ao,Bo
      WRITE(F,'("#   vector :",3F16.3)' ) REFV(:,AN)
!     Triclinc PBC's
      NAr = NA(STN)
      LVr(:,:) = LV(:,:,STN)
!     Resize lattice vector in vacuum direction 
      IF( MKVAC(1) ) THEN
        LVr(SDIM,SDIM)=MMN(2,SDIM,STN)-MMN(1,SDIM,STN)+MINBF(STN)*2
      ENDIF 
!     Calculate fractional coordinates 
      ALLOCATE( Fr(NDIM,NAr) ) 
      CALL fracR(STN,LVr,1)
!
!     Loop over all molecules
      DO MNi=1,MTOT 
       Mi = MPNT(MNi,STN)
       Mf = MPNT(MNi+1,STN)-1           
!
!      Get center of mass in surface direction 
!
       RZ = MCOMS(SDIM,MNi,STN)*RSCALE(STN) 
       R = ANINT( RZ ) - RSHFT 

       DO Ni=Mi,Mf 
        I = MOLST(Ni,STN)
        Ai =  GRID(I,STN) 
        RAi(:) = R0(:,I,STN)
        Fi(:) = Fr(:,I)
!       Loop over other atoms in molecule
        DO Nj=Mi,Mf
         J = MOLST(Nj,STN)
         Bi = GRID(J,STN)
         DO AN=1,ACNT 
          Ao = ANAX(1,AN,STN)
          Bo = ANAX(2,AN,STN)
          F = FID + AN
          IF( Ai.EQ.Ao.AND.Bi.EQ.Bo)  THEN
             Fj(:) = Fr(:,J)
             DO D=1,NDIM
               DF(D) = Fj(D) - Fi(D)
               IF( PBCS(D,STN) )  DF(D) = DF(D) - ANINT(DF(D))
             ENDDO
!
!             Calculate cartion coordinates 
!
             CALL  REALF(LVr,DF,DR) 
!             DR(1)=DF(1)*LVr(1,1) + DF(2)*LVr(2,1)+ DF(3)*LVr(3,1)
!             DR(2)=DF(1)*LVr(1,2) + DF(2)*LVr(2,2) + DF(3)*LVr(3,2)
!             DR(3)=DF(1)*LVr(1,3) + DF(2)*LVr(2,3) + DF(3)*LVr(3,3)
             HDR(1) = DR(1)/2.0d0
             HDR(2) = DR(2)/2.0d0
             HDR(3) = DR(3)/2.0d0
             RCENT(:) = RAi(:) + HDR(:)
             RSQ = 0.0d0
             ADB = 0.0d0
             DO D=1,NDIM
               RSQ = RSQ + DR(D)*DR(D)
               ADB = ADB + DR(D)*REFV(D,AN)
             ENDDO
             MGR = SQRT(RSQ)
             COSTH = ADB/MGR/MREF(AN) 
!             C = NINT( (COSTH + 1.d0)*ANGSCL ) + 1
!             IF( C.GT.SZ .OR. C.LT.1 ) THEN
!               WRITE(ESTAT,'(I12)') C
!               CALL PRNERROR(-2,ESTAT)
!             ENDIF
!             THSUM = THSUM + COSTH
             THETA = ACOS(COSTH)
             THETAd = THETA*360.0d0/(2.0d0*PI)
             THSUM(AN) = THSUM(AN) + THETAd
             C =  ANINT( THETAd )
             THDIS(C,AN) = THDIS(C,AN) + 1
             RAVE(C,:,AN) = RAVE(C,:,AN) + RCENT(:)
             DPTHCNT(R,AN) =  DPTHCNT(R,AN) + 1 
             DPTHAVE(R,AN) =  DPTHAVE(R,AN) + THETAd 
             IF( VERB ) THEN
                WRITE(F,3602) I,J,THETAd,RCENT,MGR,Ai,Bi,MNi,Mf-Mi+1
             ENDIF
          ENDIF
         ENDDO  ! ANCNT 
        ENDDO
       ENDDO
      ENDDO
!     Print distribution 
      DO AN=1,ACNT 
        F = FID + AN
        NPNTS = 0
        DO C=1,SZ 
          NPNTS  = NPNTS + THDIS(C,AN)
        ENDDO
!       Calculate average
        AVEC = THSUM(AN)/FLOAT(NPNTS)
        WRITE(F,*) "# total vectors found", NPNTS      
        WRITE(F,*) "# average cos \theta", AVEC 
!     
        DO C=0,SZ
          CNT = THDIS(C,AN)
          THETAd = FLOAT( C )
          IF( CNT.GT.0) THEN
            DO D=1,NDIM
              RCENT(D) = RAVE(C,D,AN)/FLOAT(CNT)
            ENDDO
          ELSE
            RCENT(:) = (/0.0d0,0.0d0,0.0d0/)
          ENDIF
          WRITE(F,3601) THETAd,CNT,RCENT(:)
        ENDDO
        CLOSE(F)
      ENDDO     ! ANCNT 
!
!     Write out depth profile of angles
!
      F = FID + ACNT + 1
      WRITE(F,3604) ,
      DO R=1,RRANG
        RZ = FLOAT(R + RSHFT)*RSCALE(STN) 
        TACNT = 0 
        TASUM = 0.0d0 
        DO AN=1,ACNT
          CNT = DPTHCNT(R,AN) 
          IF( CNT .GT. 0) THEN
            THETAd = DPTHAVE(R,AN)
            THETA = THETAd / (FLOAT(CNT))
          ELSE
            THETA = 0.0d0
          ENDIF
          PROF(AN) = THETA 
!
!         Track average of all considered angles
!
          TACNT = TACNT + CNT 
          TASUM = TASUM + THETAd 
        ENDDO ! ACNT 
        TAAVE = 0.0d0 
        IF( TACNT.GT.0)TAAVE = TASUM/FLOAT(TACNT)
        WRITE(TAOUT,*) PROF(:)
        WRITE(F,3603) RZ,TAAVE,TACNT,TAOUT
      ENDDO !R=1,RRANG 
      CLOSE(F)

      DEALLOCATE(THDIS,RAVE,Fr,REFV,THSUM,MREF
     & ,DPTHCNT,DPTHAVE,PROF)
!
      RETURN
 3601 FORMAT(F12.6,I6,3F18.6)
 3602 FORMAT(" # ",2I8,5F16.6," ",2A6," ",2I6)
 3603 FORMAT(2F10.3,I6,A)
 3604 FORMAT(" Hieght in ",I6,'tot ave ang | average angles of 1:',I6)
      END SUBROUTINE ANGAXM
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/02/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate dipole 
!
      SUBROUTINE TOTDIP(STN)
      USE structure
      USE const
      USE elements 
      USE potential 
      USE specify 
!
      IMPLICIT none
!
      INTEGER :: I,STN,D
      REAL*8  :: TDIPr,DIPr(NDIM),TDSQr
!
      CALL CENTERMASS(STN)
      DIPr(:) = 0.0d0 
      DO I=1,NA(STN) 
        DO D=1,NDIM
          DIPr(D) = DIPr(D)+ R0(D,I,STN)*ACHG(I,STN)
       ENDDO 
      ENDDO 

      TDSQr = 0.0d0  
      DO D=1,NDIM
          TDSQr = TDSQr+DIPr(D)*DIPr(D)
      ENDDO 
      TDIPr = SQRT(TDSQr) 
      DIP(STN) = TDIPr*DCON   ! D/(eA) 
      IF( VERB) THEN
        WRITE(6,*) " Dipole x y z "
        WRITE(6,*) DIPr(:)
      ENDIF
!
      RETURN
      END SUBROUTINE TOTDIP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/02/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate electric field at each molecule's cmos
!      E = 1/(4pi e)  sum i-N qi/r^2 dot hat(r)
!
      SUBROUTINE MOLEFLD(STN)
      USE structure
      USE const
      USE elements 
      USE potential 
      USE specify 
!
      IMPLICIT none 
!
      INTEGER :: STN,NAr,D,I,MNi,Is,Js,Ks
     & ,EFRANG,EPRANG,EFi,EPi,EPMINi,EPMAXi
     & ,EFx,EFy,EFz,EFMINi,EFMAXi 
      REAL*8 :: R(NDIM),LVr(NDIM,NDIM),RSQ,CH
     & ,Fm(NDIM),Fi(NDIM),FD(NDIM),DR(NDIM),DF(NDIM)
     & ,Eo(NDIM),Vo,Ei(NDIM),Vi,MAGR,UR(NDIM),FPE,Ef,Etot
     & ,Ir,Jr,Kr,SDF(NDIM),EFr
     & ,EFMAX,EPMAX,EFBIN,EPBIN,EPMIN,EFMIN
     & ,EFSUM,EFXSUM,EFYSUM,EFZSUM,EPSUM
     & ,EFAVE,EFXAVE,EFYAVE,EFZAVE,EPAVE
      INTEGER, ALLOCATABLE :: EFDIS(:),EPDIS(:)
     & ,EFDISX(:),EFDISZ(:),EFDISY(:)
      CHARACTER(CHSZ) :: FINA
!
!
!     Calculate 1/(4 pi e_o) (Nm/C) ^-9
!
      FPE = ELCH/(4.0d0*PI*EPSN )*100
!
      NAr = NA(STN)
!     Triclinc PBC's
      LVr(:,:) = LV(:,:,STN)
      ALLOCATE( Fr(NDIM,NAr) ) 
      CALL fracR(STN,LVr,1)
!
!
!
!
      EFMAX = 0.0d0
      EFMIN = 1.0d16
      EPMAX = 0.0d0
      EPMIN = 1.0d16
! 
!     Open molecule electric field / potential output file
!  
      OPEN(UNIT=60,FILE=MEFLFL(STN),STATUS='UNKNOWN')
      WRITE(60,6000) 
!
!     Loop over all molecules 
      DO MNi=1,MOLCNT(STN) 
        R(:) = MCOMS(:,MNi,STN)
        Eo(:) = 0.0d0 
        Vo  = 0.0d0 
        CALL SPFRAC(LVr,Fm,R)
        DO I=1,NAr
          IF( MOLN(I,STN).NE.MNi) THEN
            Fi(:) = Fr(:,I)
            CH = ACHG(I,STN)
            DO D=1,NDIM
               DF(D) = Fm(D) - Fi(D)
               IF( PBCS(D,STN) )  DF(D) = DF(D) - ANINT(DF(D))
            ENDDO
            DO Is=-1,1
              DO Js=-1,1
                DO Ks=-1,1
                  Ir = FLOAT( Is )
                  Jr = FLOAT( Js )
                  Kr = FLOAT( Ks )
                  SDF(1) = DF(1) + Ir
                  SDF(2) = DF(2) + Jr
                  SDF(3) = DF(3) + Kr
                  CALL REALF(LVr,SDF,DR)
                 
                  RSQ = 0.0d0
                  DO D=1,NDIM
                    RSQ = RSQ + DR(D)*DR(D)
                  ENDDO  
                  MAGR = SQRT( RSQ) 
                  DO D=1,NDIM
                    UR(D) = DR(D)/MAGR
                    Eo(D) = Eo(D) + CH/RSQ*UR(D)
                  ENDDO
                  Vo = Vo + CH/MAGR 
                ENDDO
              ENDDO
            ENDDO    
          ENDIF 
        ENDDO ! NAr
!
!       Electric field  (N/C) 
!
        Ei(1) = FPE*Eo(1)*1.0d11 
        Ei(2) = FPE*Eo(2)*1.0d11
        Ei(3) = FPE*Eo(3)*1.0d11
!
!       Electronic potential  (V)
!
        Vi = FPE*Vo*10  

        Etot = Ei(1)*Ei(1)+Ei(2)*Ei(2)+Ei(3)*Ei(3)
        Ef = SQRT( Etot )

!
!       Write potential and efield for each molecule
!
        WRITE(60,6001) MNi,R(:),Ef,Ei(:),Vi

        MOLEF(:,MNi,STN) = Ei(:)
        MOLEP(MNi,STN) = Vi
!
!       Track maximum and minimum 
!

        IF( EF.GT.EFMAX ) EFMAX = Ef
        IF( Ei(1).GT.EFMAX ) EFMAX = Ei(1)
        IF( Ei(2).GT.EFMAX ) EFMAX = Ei(2)
        IF( Ei(3).GT.EFMAX ) EFMAX = Ei(3)
        IF( EF.LT.EFMIN ) EFMIN = Ef
        IF( Ei(1).LT.EFMIN ) EFMIN = Ei(1)
        IF( Ei(2).LT.EFMIN ) EFMIN = Ei(2)
        IF( Ei(3).LT.EFMIN ) EFMIN = Ei(3)
        IF( Vi.GT.EPMAX ) EPMAX = Vi
        IF( Vi.LT.EPMIN ) EPMIN = Vi

!
      ENDDO ! MNi


!debug 
      WRITE(60,6002) EFMAX,EFMIN
      WRITE(60,6003) EPMAX,EPMIN

      CLOSE(60)


      EFBIN = 1.d7
      EPBIN = 0.1d0 

!
!     Calculate distribution 
!
      EFMINi = INT(EFMIN/EFBIN) - 1
      EFMAXi = INT(EFMAX/EFBIN) + 1
      EFMAX = FLOAT( EFMAXi)*EFBIN
      EFMIN = FLOAT( EFMINi)*EFBIN
      EFRANG = INT( (EFMAX-EFMIN)/ EFBIN ) + 1

!     Set potentials to integer values 
      EPMINi = INT( EPMIN ) - 1
      EPMAXi = INT( EPMAX ) + 1
      EPMIN = FLOAT( EPMINi)
      EPMAX = FLOAT( EPMAXi)
      EPRANG = INT( (EPMAX - EPMIN )/ EPBIN ) + 1
      ALLOCATE( EFDIS(0:EFRANG)
     & ,EFDISX(0:EFRANG),EFDISY(0:EFRANG),EFDISZ(0:EFRANG)
     & ,EPDIS(0:EPRANG))
      EFDIS(:) = 0
      EFDISX(:) = 0
      EFDISY(:) = 0
      EFDISZ(:) = 0
      EPDIS(:) = 0
!
!
!     Open distribution files
!
      WRITE(FINA,6100) STN
      OPEN(UNIT=61,FILE=FINA,STATUS='UNKNOWN')
      WRITE(61,*) 
      WRITE(61,6002) EFMAX,EFMIN
      WRITE(61,6102)   EFBIN,EFRANG
!
      WRITE(FINA,6200) STN
      OPEN(UNIT=62,FILE=FINA,STATUS='UNKNOWN')
      WRITE(62,6003) EPMAX,EPMIN
      WRITE(62,6202)   EPBIN,EPRANG

!     Loop over all molecules
!
      DO MNi=1,MOLCNT(STN) 
        Eo(:) = MOLEF(:,MNi,STN)
        Etot = Eo(1)*Eo(1)+Eo(2)*Eo(2)+Eo(3)*Eo(3)
        Ef = SQRT( Etot ) - EFMIN
!        
        Ei(1) = Eo(1) - EFMIN
        Ei(2) = Eo(2) - EFMIN
        Ei(3) = Eo(3) - EFMIN
!
        EFi = ANINT( Ef/EFBIN  ) 
        EFx = ANINT( Ei(1)/EFBIN  ) 
        EFy = ANINT( Ei(2)/EFBIN  ) 
        EFz = ANINT( Ei(3)/EFBIN  )

        Vo = MOLEP(MNi,STN) - EPMIN
        EPi = ANINT( Vo/ EPBIN  )

! debug
        WRITE(506,'(5E12.2)') Ef,Ei(:),Vo
        WRITE(506,'(5I12)') EFi,EFx,EFy,EFz,EPi
 

        EFDIS(EFi) = EFDIS(EFi) + 1
        EFDISX(EFx) = EFDISX(EFx) + 1
        EFDISY(EFY) = EFDISY(EFY) + 1
        EFDISZ(EFZ) = EFDISZ(EFZ) + 1
!
        EPDIS(EPi) = EPDIS(EPi) + 1
!
      ENDDO
!
!     Print distributions
! 
      WRITE(61,6110) 
      DO EFi=0,EFRANG
       EFr = FLOAT( EFi )*EFBIN + EFMIN
!      EFr = FPE*EF*1.0d11 
       WRITE(61,6101) EFr,EFDIS(EFi),EFDISX(EFi),EFDISY(EFi),EFDISZ(EFi)
      ENDDO 
      CLOSE(61) 
      WRITE(62,6210) 
      DO EPi=0,EPRANG
       Vi = FLOAT( EPi )*EPBIN + EPMIN
!       Electronic potential  (V)
!       Vi = FPE*Vo*10  
       WRITE(62,6201) Vi,EPDIS(EPi)
      ENDDO 
      CLOSE(62) 

      DEALLOCATE( Fr,EFDIS,EPDIS
     & ,EFDISX,EFDISY,EFDISZ ) 
!
      RETURN   
 6000 FORMAT('# Mol#'
     & ,'| centermass_x | centermass_y | centermass_z  '
     & ,'! E_tot            '
     & ,'| E_x              | E_y               | E_x               '
     & ,'| V ')
 6001 FORMAT(I6,3F12.4,5E16.4)
 6002 FORMAT('# Efield Max: ',E16.4,' Min: ',E16.4)
 6003 FORMAT('# V Max: ',E16.4,' Min: ',E16.4)
 6111 FORMAT(I6,5E16.4,3F12.4)
 6100 FORMAT('mE-',I1,'.dat')
 6200 FORMAT('mV-',I1,'.dat')
 6101 FORMAT(E16.4,4I6)
 6201 FORMAT(F12.4,I6)
 6102 FORMAT('# Bin size : ',E16.4,' Range: ',I12)
 6202 FORMAT('# Bin size : ',F12.4,' Range: ',I12)
 6110 FORMAT('# Efield   ! N_tot            '
     & ,'| N_x              | N_y               | N_x               ')
 6210 FORMAT('# Potential   ! N              ')

      END SUBROUTINE MOLEFLD 

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/02/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate electric field at each atom
!      E = 1/(4pi e)  sum i-N qi/r^2 dot hat(r) 
!
      SUBROUTINE ATOMEFLD(STN)
!
      USE structure
      USE const
      USE elements 
      USE potential 
      USE specify 
!
      IMPLICIT none 
!
      INTEGER :: STN,NAr,D,I,J,Is,Js,Ks
     & ,EFRANG,EPRANG,EFi,EPi,EPMINi,EPMAXi
     & ,EFx,EFy,EFz,EFMINi,EFMAXi
      REAL*8 :: R(NDIM),LVr(NDIM,NDIM),RSQ,CHi,CHj
     & ,Fj(NDIM),Fi(NDIM),FD(NDIM),DR(NDIM),DF(NDIM)
     & ,Eo(NDIM),Vo,Ei(NDIM),Vi,MAGR,UR(NDIM),FPE,Ef,Etot
     & ,Ir,Jr,Kr,SDF(NDIM),EFr
     & ,EFMAX,EPMAX,EFBIN,EPBIN,EPMIN,EFMIN
      REAL*8, ALLOCATABLE :: EFDIS(:),EPDIS(:)
     & ,EFDISX(:),EFDISZ(:),EFDISY(:)
      CHARACTER(CHSZ) :: FINA
!
!
!     Calculate 1/(4 pi e_o) (Nm/C) ^-9
!
      FPE = ELCH/(4.0d0*PI*EPSN )*100

!
      NAr = NA(STN)
!     Triclinc PBC's
      LVr(:,:) = LV(:,:,STN)
      ALLOCATE( Fr(NDIM,NAr) ) 
      CALL fracR(STN,LVr,1)
!
!
!
      EFBIN = 1.d11
      EPBIN = 1.d0 
!
      EFMAX = 0.0d0
      EFMIN = 1.0d16
      EPMAX = 0.0d0
      EPMIN = 1.0d16
!
!     Loop over all atomic pairs 
      DO I=1,NAr
        R(:) = R0(:,I,STN)
        CHi = ACHG(I,STN)
        Eo(:) = 0.0d0 
        Vo  = 0.0d0 
        CALL SPFRAC(LVr,Fi,R)
        DO J=1,NAr
          IF( I.NE.J) THEN 
            CHj = ACHG(J,STN)
            Fj(:) = Fr(:,J)
            DO D=1,NDIM
               DF(D) = Fj(D) - Fi(D)
               IF( PBCS(D,STN) )  DF(D) = DF(D) - ANINT(DF(D))
            ENDDO
            DO Is=-1,1
              DO Js=-1,1
                DO Ks=-1,1
                  Ir = FLOAT( Is )
                  Jr = FLOAT( Js )
                  Kr = FLOAT( Ks )
                  SDF(1) = DF(1) + Ir
                  SDF(2) = DF(2) + Jr
                  SDF(3) = DF(3) + Kr
                  CALL REALF(LVr,SDF,DR)
                 
                  RSQ = 0.0d0
                  DO D=1,NDIM
                    RSQ = RSQ + DR(D)*DR(D)
                  ENDDO  
                  MAGR = SQRT( RSQ) 
                  DO D=1,NDIM
                    UR(D) = DR(D)/MAGR
                    Eo(D) = Eo(D) + CHj/RSQ*UR(D)
                  ENDDO
                  Vo = Vo + CHj/MAGR 

                ENDDO
              ENDDO
            ENDDO    
          ENDIF 
        ENDDO ! NAr
!
!       Electric field  (N/C) 
!
        Ei(1) = FPE*Eo(1)*1.0d11 
        Ei(2) = FPE*Eo(2)*1.0d11
        Ei(3) = FPE*Eo(3)*1.0d11
!
!       Electronic potential  (V)
!
        Vi = FPE*Vo*10  

        Etot = Ei(1)*Ei(1)+Ei(2)*Ei(2)+Ei(3)*Ei(3)
        Ef = SQRT( Etot )

        ATOMEF(:,I,STN) = Ei(:)
        ATOMEP(I,STN) = Vi
        IF( EF.GT.EFMAX ) EFMAX = Ef
        IF( Ei(1).GT.EFMAX ) EFMAX = Ei(1)
        IF( Ei(2).GT.EFMAX ) EFMAX = Ei(2)
        IF( Ei(3).GT.EFMAX ) EFMAX = Ei(3)
        IF( EF.LT.EFMIN ) EFMIN = Ef
        IF( Ei(1).LT.EFMIN ) EFMIN = Ei(1)
        IF( Ei(2).LT.EFMIN ) EFMIN = Ei(2)
        IF( Ei(3).LT.EFMIN ) EFMIN = Ei(3)
        IF( Vi.GT.EPMAX ) EPMAX = Vi
        IF( Vi.LT.EPMIN ) EPMIN = Vi

!
      ENDDO ! MNi

!
!
!     Calculate distribution 
!
      EFMINi = INT(EFMIN/EFBIN) - 1
      EFMAXi = INT(EFMAX/EFBIN) + 1
      EFMAX = FLOAT( EFMAXi)*EFBIN
      EFMIN = FLOAT( EFMINi)*EFBIN
      EFRANG = INT( (EFMAX-EFMIN)/ EFBIN ) + 1

!     Set potentials to integer values 
      EPMINi = INT( EPMIN ) - 1
      EPMAXi = INT( EPMAX ) + 1
      EPMIN = FLOAT( EPMINi)
      EPMAX = FLOAT( EPMAXi)
      EPRANG = INT( (EPMAX - EPMIN )/ EPBIN ) + 1
      ALLOCATE( EFDIS(0:EFRANG)
     & ,EFDISX(0:EFRANG),EFDISY(0:EFRANG),EFDISZ(0:EFRANG)
     & ,EPDIS(0:EPRANG))
      EFDIS(:) = 0
      EPDIS(:) = 0
!


!debug 
      WRITE(709,*) EFMAX,EFMIN
      WRITE(709,*) EFBIN,EFRANG
      WRITE(709,*) EPMAXi,EPMINi
      WRITE(709,*) EPMAX,EPMIN
      WRITE(709,*) EPBIN,EPRANG


!     Loop over all atoms
!
      DO I=1,NAr
        Eo(:) = ATOMEF(:,I,STN)
        Etot = Eo(1)*Eo(1)+Eo(2)*Eo(2)+Eo(3)*Eo(3)
        Ef = SQRT( Etot ) - EFMIN
!
        Ei(1) = Eo(1) - EFMIN
        Ei(2) = Eo(2) - EFMIN
        Ei(3) = Eo(3) - EFMIN
!
        EFi = ANINT( Ef/EFBIN  ) 
        EFx = ANINT( Ei(1)/EFBIN  ) 
        EFy = ANINT( Ei(2)/EFBIN  ) 
        EFz = ANINT( Ei(3)/EFBIN  ) 

        Vo = ATOMEP(I,STN) - EPMIN
        EPi = ANINT( Vo/ EPBIN  )
!
        EFDIS(EFi)  = EFDIS(EFi) + 1
        EFDISX(EFx) = EFDISX(EFx) + 1
        EFDISY(EFY) = EFDISY(EFY) + 1
        EFDISZ(EFZ) = EFDISZ(EFZ) + 1
        EPDIS(EPi)  = EPDIS(EPi) + 1
!
      ENDDO
!
!     Print distributions
! 
      WRITE(FINA,6100) STN
      OPEN(UNIT=61,FILE=FINA,STATUS='UNKNOWN')
      WRITE(61,*) 
      DO EFi=0,EFRANG
       EFr = FLOAT( EFi )*EFBIN + EFMIN
!      EFr = FPE*EF*1.0d11 
       WRITE(61,6101) EFr,EFDISX(EFi),EFDISY(EFi),EFDISZ(EFi),EFDIS(EFi)
      ENDDO 
      CLOSE(61) 
!
      WRITE(FINA,6200) STN
      OPEN(UNIT=62,FILE=FINA,STATUS='UNKNOWN')
      DO EPi=0,EPRANG
       Vi = FLOAT( EPi )*EPBIN + EPMIN
!       Electronic potential  (V)
!       Vi = FPE*Vo*10  
       WRITE(62,6201) Vi,EPDIS(EPi)
      ENDDO 
      CLOSE(62) 

      DEALLOCATE( Fr,EFDIS,EPDIS
     & ,EFDISX,EFDISY,EFDISZ ) 
!
      RETURN
 6111 FORMAT(I6,5E16.4,3F12.4)
 6100 FORMAT('aE-',I1,'.dat')
 6200 FORMAT('aV-',I1,'.dat')
 6101 FORMAT(E16.4,4F12.4)
 6201 FORMAT(2F12.4)
      END SUBROUTINE ATOMEFLD 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/02/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate dipole for each molecule 
      SUBROUTINE CMOLDIP(STN)
      USE structure
      USE const
      USE elements 
      USE potential 
      USE specify 
!
      IMPLICIT none 
!
      INTEGER :: STN,MNi,Mi,Mf,Ni,I,D
     & ,RSHFT,DSHFT,DRANG,RRANG,R,DCNT,DIM 
     & ,Dti,Dxi,Dyi,Dzi 
      INTEGER, ALLOCATABLE :: MDDIS(:,:)
     & ,MDDISt(:),MDDISx(:),MDDISy(:),MDDISz(:)  
      REAL*8 :: DIPr(NDIM),DIPCON(NDIM) 
     & ,RZ,DZ,DSUM,DAVE,DTOT,Dt,Dx,Dy
      CHARACTER(12*51) :: MDSTC
       CHARACTER(CHSZ) :: ESTAT
!
!     Open molecular dipole file 
!
      OPEN(UNIT=51,FILE=CMDIPFL(STN),STATUS='UNKNOWN')
      WRITE(51,500) 
!
!      If distribution needs to be print open distribution file
!
      IF( CMLDPD(STN) ) THEN
        OPEN(UNIT=50,FILE=MOLDIPAXFL(STN),STATUS='UNKNOWN')
        OPEN(UNIT=52,FILE=MOLDIPFL(STN),STATUS='UNKNOWN')
!
!       Set spatial range for dipole distribution
!
        DIM = MDDIM(STN) 
        CALL MAXMIN(STN)
        RSHFT = ANINT( MMN(1,DIM,STN)*RSCALE(STN) ) - RSCALE(STN) 
        RRANG = ANINT( MMN(2,DIM,STN)*RSCALE(STN) ) -  RSHFT 
!
!
        DRANG = ANINT( 2.0d0*DMAX(STN)*DSCALE(STN) + 1.d0 )    !Max dybe of molecule 
        DSHFT = ANINT(DMAX(STN)*DSCALE(STN) + 1.0d0)
        ALLOCATE( MDDIS(RRANG,DRANG),MDDISt(DRANG)
     &  ,MDDISx(DRANG),MDDISy(DRANG),MDDISz(DRANG)  ) 
        MDDIS(:,:) = 0  
        MDDISt(:) =  0
        MDDISx(:) = 0
        MDDISy(:) = 0
        MDDISz(:) = 0
      ENDIF
!
!     Loop over all molecules
!
      DO MNi=1,MOLCNT(STN) 
        Mi = MPNT(MNi,STN)
        Mf = MPNT(MNi+1,STN)-1           
        DIPr(:) = 0.0d0 
        DO Ni=Mi,Mf 
          I = MOLST(Ni,STN)
          DO D=1,NDIM
             DIPr(D) = DIPr(D)+ R0(D,I,STN)*ACHG(I,STN)
          ENDDO 
        ENDDO ! Ni=Mi,Mf 
        DSUM  = 0.0d0
        DO D=1,NDIM
          DIPCON(D) = DIPr(D)*DCON  
          MOLDIP(D,MNi,STN) =  DIPCON(D)
          DSUM = DSUM +  DIPCON(D)*DIPCON(D)
        ENDDO 
        DTOT = SQRT(DSUM)
        IF(  CMLDPD(STN) ) THEN 
!
!         Save general dipole distribution 
!
          
          Dt = DTOT*DSCALE(STN)  
          Dx = DIPCON(1)*DSCALE(STN)  
          Dy = DIPCON(2)*DSCALE(STN)  
          Dz = DIPCON(3)*DSCALE(STN)  
          Dti = ANINT( Dt ) + DSHFT 
          Dxi = ANINT( Dx ) + DSHFT 
          Dyi = ANINT( Dy ) + DSHFT 
          DZi = ANINT( DZ ) + DSHFT 
          MDDISt(Dti) = MDDISt(Dti) + 1
          MDDISx(Dxi) = MDDISx(Dxi) + 1
          MDDISy(Dyi) = MDDISy(Dyi) + 1
          MDDISz(Dzi) = MDDISz(Dzi) + 1
!
!         Print dipole of each molecule 
!
          WRITE(51,501) MNi,DTOT,DIPCON(:),MCOMS(:,MNi,STN) 
!
!         Save dipole distribution along specified axis
!
          RZ = MCOMS(DIM,MNi,STN)*RSCALE(STN) 
          R = ANINT( RZ ) - RSHFT 
          DZ = DIPCON(DIM)*DSCALE(STN)  
          D = ANINT( DZ ) + DSHFT 
          IF( D.GT.DRANG) THEN
            WRITE(ESTAT,*)  D,' greater than ',DMAX(STN),' for mole',MNi
            CALL  prnerror(81,ESTAT)
          ENDIF 
          MDDIS(R,D) = MDDIS(R,D) + 1
        ENDIF 
      ENDDO   ! MNi=1,MTOT 
!
      IF( CMLDPD(STN) ) THEN
!        WRITE(50,503) -1*DRANG/2,DRANG/ANINT(DSCALE(STN)),DRANG/2
        WRITE(52,520)  
        WRITE(52,504)  DSCALE(STN),DRANG,DSHFT
!
!       Print molecular dipole histogram
!
        DCNT = 0
        DO D=1,DRANG 
          Dti = D - DSHFT
          Dt  = FLOAT(Dti)/DSCALE(STN)
          WRITE(52,5202) Dt,MDDISt(D),MDDISx(D),MDDISy(D),MDDISz(D)  
        ENDDO  ! DRANG
!
!        Print dipole distribution along axis 
!
        WRITE(50,520)  
        WRITE(50,504)  DSCALE(STN),DRANG,DSHFT
        DO R=1,RRANG
          DSUM = 0.0d0
          DCNT = 0
          DO D=1,DRANG 
             DCNT = DCNT + MDDIS(R,D) 
             DZ = FLOAT(D - DSHFT)/DSCALE(STN) 
             DSUM =  DSUM + DZ*MDDIS(R,D) 
          ENDDO  ! DRANG
          DAVE = 0.0d0  
          IF( DCNT.GT.0) DAVE = DSUM / FLOAT( DCNT ) 
          RZ = FLOAT(R + RSHFT)*RSCALE(STN) 
          WRITE(MDSTC,*) MDDIS(R,:) 
          WRITE(50,502) RZ,DAVE,MDSTC        
        ENDDO  !RRANG 
        CLOSE(50)
        DEALLOCATE( MDDIS,MDDISt,MDDISx,MDDISy,MDDISz ) 
      ENDIF 
!
      CLOSE(51)
      RETURN
 500  FORMAT('#  Molecule # | Total (D) |   mu_x      |       mu_y     ',
     & ' |       mu_z  | centermass_x | centermass_y | centermass_z  ')
 520  FORMAT('#  Dipole | N_total | N_x | N_y | N_z ')
 501  FORMAT(I6,7F18.3)
 502  FORMAT(' ',2F18.3,A)
 503  FORMAT('# Centermass_z |  <mu_z>  |  histogram '
     & ,I6,':',I6,':',I6)
 504  FORMAT('#  Dipole scale :',F10.3,' range :',I6,' shift ;',I6)
 5202 FORMAT(F12.4,4I6)
      END SUBROUTINE CMOLDIP 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/02/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate dipole momment distribution for an index
!
      SUBROUTINE CINDPROP(STN)
      USE structure
      USE const
      USE specify 
      USE potential 
!
      IMPLICIT none
!
      INTEGER :: STN,N,CMX,I,C,D,IND
      INTEGER, ALLOCATABLE :: INDX(:)
      REAL*8 :: TOTMr,RMAS(NDIM),MS,DIPr(NDIM),TCHr,CHr
!
      DO N=1,IDPN(STN)
        ALLOCATE( INDX(NA(STN)) ) 
        IF( IDPID(N,STN).EQ.1 ) INDX(:) = MOLN(:,STN)
        IF( IDPID(N,STN).EQ.2 ) INDX(:) = RN(:,STN)  
        IF( IDPID(N,STN).EQ.3 ) INDX(:) = CHN(:,STN)      
        IF( IDPID(N,STN).EQ.4 ) INDX(:) = THRM(:,STN)
!
!       Loop over all atoms find max index value
        CMX = 0
        DO I=1,NA(STN)
          IF( INDX(I).GT.CMX )CMX= INDX(I) 
        ENDDO 
!       Loop over all index values and cal dipole moment of each group 
        DO C=1,CMX
          TOTMr = 0.0d0
          RMAS(:) = 0.0d0
!         Find center of mass of group 
          DO I=1,NA(STN)
            IF( INDX(I).EQ.C) THEN
              MS = AMAS(I,STN ) 
              TOTMr = TOTMr + MS 
              DO D=1,NDIM
               RMAS(D) = RMAS(D) + MS*R0(D,I,STN)
              ENDDO 
            ENDIF 
          ENDDO 
          COMAS(1) = RMAS(1)/TOTMr
          COMAS(2) = RMAS(2)/TOTMr
          COMAS(3) = RMAS(3)/TOTMr
!         Find center dipole of center of mass structure 
          DIPr(:) = 0.0d0 
          TCHr = 0.0d0 
          DO I=1,NA(STN)
            IF( INDX(I).EQ.C) THEN
              CHR = ACHG(I,STN)
              TCHr = TCHr + CHR 
              DO D=1,NDIM
               DIPr(D) = DIPr(D)+ (R0(D,I,STN)-COMAS(D) )*CHR 
              ENDDO 
            ENDIF 
          ENDDO 
          GPROP(1:3,C,STN) = COMAS(:) 
          GPROP(4,C,STN) =  TOTMr 
          GPROP(5,C,STN) =  TCHr  
          DO D=1,NDIM 
            IND = D + 5
            GPROP(IND,C,STN) = DIPr(D)*DCON   
          ENDDO
        ENDDO !
        DEALLOCATE( INDX ) 
      ENDDO  !
!
      RETURN
      END SUBROUTINE  CINDPROP
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate the angle between a vector and an axis

      SUBROUTINE ANGAXM_TRI(STN)
      USE specify
      USE structure
      USE build 
!
      Implicit none
!
      INTEGER :: STN,I,ELi,ELj,D,C,SZ,CNT,NPNTS,J,Jo,Jf,Nj,Ni,NAr
     & ,MTOT,MNi,Mi,Mf,ACNT,FID,F,AN ,R,RRANG,RSHFT ,TACNT,K,Nk
      REAL*8 :: X,Y,Z,AA,AB,AC,BA,BB,BC,CA,CB,CC,Fi(NDIM),Fj(NDIM)
     & ,MLV(NDIM),DF(NDIM),DR(NDIM),RSQ,Fk(NDIM),Dv(NDIM)
     & ,RAi(NDIM),RCENT(NDIM),LVr(NDIM,NDIM),HDR(NDIM)
     & ,MGR,ADB,COSTH,ANGSCL,AVEC,THETA,THETAd
     & ,RZ,TASUM,TAAVE,COSTH2
      INTEGER, ALLOCATABLE :: THDIS(:,:),DPTHCNT(:,:)
      REAL*8, ALLOCATABLE :: RAVE(:,:,:),THSUM(:),REFV(:,:),MREF(:)
     & ,DPTHAVE(:,:),PROF(:) 
      CHARACTER(16*11) :: TAOUT 
      CHARACTER(CHSZ) :: ESTAT,FL
      CHARACTER(IDSZ) :: Ai,Bi,Ao,Bo,Co,Ci
!     
!     Set the local variables
      ACNT = TANCNT(STN)               ! number of angles to calculate 
      FID  = 50                       ! Initial file number 
      MTOT = MOLCNT(STN)
!     Verbose output
      IF (VERB ) THEN
        WRITE(*,*) "Starting 3 point angle calc within molecules"
        WRITE(*,*) "  Structure ",STN," with ", MTOT," molecules"
        WRITE(*,*) " Will calculate ",ACNT," angles "
      ENDIF
!     Set scale for bins
        CALL MAXMIN(STN)
        RSHFT = ANINT( MMN(1,SDIM,STN)*RSCALE(STN) ) - RSCALE(STN) 
        RRANG = ANINT( MMN(2,SDIM,STN)*RSCALE(STN) ) -  RSHFT 

       ANGSCL = 1.0d0
!      SZ = INT(2.0d0*ANGSCL) + 1   !-1 to 1 for cos
      SZ = INT(180.0d0*ANGSCL) + 2   ! 0 - 181
!
!     Allocate count arrays
! 
      ALLOCATE( THDIS(0:SZ,ACNT),RAVE(0:SZ,NDIM,ACNT)
     & ,REFV(NDIM,ACNT),THSUM(ACNT),MREF(ACNT)
     & ,DPTHCNT(RRANG,ACNT),DPTHAVE(RRANG,ACNT) 
     & ,PROF(ACNT) )

      THDIS(:,:) = 0.0d0
      THSUM(:) = 0.0d0
      RAVE(:,:,:) = 0.0d0

      DPTHCNT(:,:) =  0
      DPTHAVE(:,:) =  0.0d0 

!     Open ouput files
      DO AN=1,ACNT 
        F = FID + AN 
        Ao = TANAX(1,AN,STN)
        Bo = TANAX(2,AN,STN)
        Co = TANAX(3,AN,STN)
        REFV(:,AN) = TANGX(:,AN,STN)
        MREF(AN) = 0.0d0
        DO D =1,NDIM
          MREF(AN) = MREF(AN) + REFV(D,AN)*REFV(D,AN)
        ENDDO
        OPEN(UNIT=F,FILE=TANGAXFL(AN,STN),STATUS='UNKNOWN')
        IF(VERB) THEN
        WRITE(*,*) "opening TANGAXFL(AN,STN)",AN,TANGAXFL(AN,STN)
        ENDIF
        WRITE(F,*) "# Angular distribution"
        WRITE(F,*) "# Scale ", 1/ANGSCL
        WRITE(F,'("# Atom names:",3A4)' ) Ao,Bo,Co
        WRITE(F,'("#   vector :",3F16.3)' ) REFV(:,AN)
      ENDDO 
      F = FID + ACNT + 1
!     WRITE(FL,'Prof-',A) = TANGAXFL(AN,STN)
      OPEN(UNIT=F,FILE=DPTANGAXFL(STN),STATUS='UNKNOWN')
        IF(VERB) THEN
        WRITE(*,*) "opening DPTANGAXFL(STN)",DPTANGAXFL(STN)
        ENDIF
      WRITE(F,*) "# Depth distribution of angles in",SDIM
      WRITE(F,*) "# Scale ", RSCALE(STN),RSHFT
      WRITE(F,'("# Atom names:",3A4)' ) Ao,Bo,Co
      WRITE(F,'("#   vector :",3F16.3)' ) REFV(:,AN)
!     Triclinc PBC's
      NAr = NA(STN)
      LVr(:,:) = LV(:,:,STN)
!     Resize lattice vector in vacuum direction 
      IF( MKVAC(1) ) THEN
        LVr(SDIM,SDIM)=MMN(2,SDIM,STN)-MMN(1,SDIM,STN)+MINBF(STN)*2
      ENDIF 
!     Calculate fractional coordinates 
      ALLOCATE( Fr(NDIM,NAr) ) 
      CALL fracR(STN,LVr,1)
!
!     Loop over all molecules
      DO MNi=1,MTOT 
       Mi = MPNT(MNi,STN)
       Mf = MPNT(MNi+1,STN)-1           
!
!      Get center of mass in surface direction 
!
       RZ = MCOMS(SDIM,MNi,STN)*RSCALE(STN) 
       R = ANINT( RZ ) - RSHFT 
!
       DO Ni=Mi,Mf 
        I = MOLST(Ni,STN)
        Ai =  GRID(I,STN) 
        RAi(:) = R0(:,I,STN)
        Fi(:) = Fr(:,I)
!       Loop over other atoms in molecule
        DO Nj=Mi,Mf
         J  = MOLST(Nj,STN)
         Bi = GRID(J,STN)
!        Loop over other atoms in molecule
         DO Nk=Mi,Mf
          K  = MOLST(Nk,STN)
          Ci = GRID(K,STN)
          DO AN=1,ACNT 
          F = FID + AN 
          Ao = TANAX(1,AN,STN)
          Bo = TANAX(2,AN,STN)
          Co = TANAX(3,AN,STN)

          F = FID + AN

          IF( Ai.EQ.Ao.AND.Bi.EQ.Bo.AND.Ci.EQ.Co)  THEN
             Fj(:) = Fr(:,J) 
             Fk(:) = Fr(:,K)
!
!              
!            V = Vk - 1/2(Vi -Vj)
!
             DO D=1,NDIM
               DF(D) =  Fk(D) - 0.5d0*Fi(D) - 0.5d0*Fj(D)
               IF( PBCS(D,STN) )  DF(D) = DF(D) - ANINT(DF(D))
             ENDDO
!
!             Calculate cartion coordinates 
!
             CALL  REALF(LVr,DF,DR) 
             HDR(1) = DR(1)/2.0d0
             HDR(2) = DR(2)/2.0d0
             HDR(3) = DR(3)/2.0d0
             RCENT(:) = RAi(:) + HDR(:)
!
!            cos(O) = a*b/( ||a|| ||b||)
!               ADB  = DR*REFV : (a*bv)
!               ||b|| = sqrt( DR*DR )
!

             RSQ = 0.0d0
             ADB = 0.0d0
             DO D=1,NDIM
               RSQ = RSQ + DR(D)*DR(D)
               ADB = ADB + DR(D)*REFV(D,AN)
             ENDDO
             MGR = SQRT(RSQ)
             COSTH = ADB/MGR/MREF(AN) 


!                Using the law of consines 
c$$$             AA = 0.0d0
c$$$             BB = 0.0d0 
c$$$             CC = 0.0d0 
c$$$             DO D=1,NDIM
c$$$               X = DR(D) - REFV(D,AN)   
c$$$               CC = CC + X*X
c$$$               BB = BB + DR(D)*DR(D)
c$$$               AA = AA + REFV(D,AN)*REFV(D,AN)
c$$$             ENDDO
c$$$             COSTH2 = AA + BB - CC
c$$$             AA=SQRT(AA)
c$$$             BB=SQRT(BB)
c$$$             COSTH2 = COSTH2/( 2.0d0*AA*BB)

!             C = NINT( (COSTH + 1.d0)*ANGSCL ) + 1
!             IF( C.GT.SZ .OR. C.LT.1 ) THEN
!               WRITE(ESTAT,'(I12)') C
!               CALL PRNERROR(-2,ESTAT)
!             ENDIF
!             THSUM = THSUM + COSTH
             THETA = ACOS(COSTH)
             THETAd = THETA*180.0d0/PI
             THSUM(AN) = THSUM(AN) + THETAd
             C =  ANINT( THETAd )
             THDIS(C,AN) = THDIS(C,AN) + 1
             RAVE(C,:,AN) = RAVE(C,:,AN) + RCENT(:)
             DPTHCNT(R,AN) =  DPTHCNT(R,AN) + 1 
             DPTHAVE(R,AN) =  DPTHAVE(R,AN) + THETAd 
             IF( VERB ) THEN
                WRITE(F,3602) I,J,THETAd,RCENT,MGR,Ai,Bi,MNi,Mf-Mi+1
             ENDIF


! debug  
            WRITE(155,*) COSTH2
            WRITE(155,*) COSTH,THETAd


          ENDIF
         ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!     Print distribution 
      DO AN=1,ACNT 
        F = FID + AN
        NPNTS = 0
        DO C=1,SZ 
          NPNTS  = NPNTS + THDIS(C,AN)
        ENDDO
!       Calculate average
        AVEC = THSUM(AN)/FLOAT(NPNTS)
        WRITE(F,*) "# total vectors found", NPNTS      
        WRITE(F,*) "# average cos \theta", AVEC 
!     
        DO C=0,SZ
          CNT = THDIS(C,AN)
          THETAd = FLOAT( C )
          IF( CNT.GT.0) THEN
            DO D=1,NDIM
              RCENT(D) = RAVE(C,D,AN)/FLOAT(CNT)
            ENDDO
          ELSE
            RCENT(:) = (/0.0d0,0.0d0,0.0d0/)
          ENDIF
          WRITE(F,3601) THETAd,CNT,RCENT(:)
        ENDDO
        CLOSE(F)
      ENDDO     ! TANCNT 
!
!     Write out depth profile of angles
!
      F = FID + ACNT + 1
      WRITE(F,3604) ,
      DO R=1,RRANG
        RZ = FLOAT(R + RSHFT)*RSCALE(STN) 
        TACNT = 0 
        TASUM = 0.0d0 
        DO AN=1,ACNT
          CNT = DPTHCNT(R,AN) 
          IF( CNT .GT. 0) THEN
            THETAd = DPTHAVE(R,AN)
            THETA = THETAd / (FLOAT(CNT))
          ELSE
            THETA = 0.0d0
          ENDIF
          PROF(AN) = THETA 
!
!         Track average of all considered angles
!
          TACNT = TACNT + CNT 
          TASUM = TASUM + THETAd 
        ENDDO ! ACNT 
        TAAVE = 0.0d0 
        IF( TACNT.GT.0)TAAVE = TASUM/FLOAT(TACNT)
        WRITE(TAOUT,*) PROF(:)
        WRITE(F,3603) RZ,TAAVE,TACNT,TAOUT
      ENDDO !R=1,RRANG 
      CLOSE(F)

      DEALLOCATE(THDIS,RAVE,Fr,REFV,THSUM,MREF
     & ,DPTHCNT,DPTHAVE,PROF)
!
      RETURN
 3601 FORMAT(F12.6,I6,3F18.6)
 3602 FORMAT(" # ",2I8,5F16.6," ",2A6," ",2I6)
 3603 FORMAT(2F10.3,I6,A)
 3604 FORMAT(" Hieght in ",I6,'tot ave ang | average angles of 1:',I6)
      END SUBROUTINE ANGAXM_TRI
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++! 



      SUBROUTINE get_bonds(STN)
      USE structure
      USE specify
      USE potential 
!
      INTEGER :: I,J,K,L,NAr,Ni,Nf,NJi,NJf,NJ,STN
     & ,NNAB,SN
     & ,NBl,NANGl,NDIHl
!
!     Set max molecules to # of atoms for now
      NAr= NA(STN)
      NBl   = 0 
      NANGl = 0
      NDIHl = 0
      
!
!
!     Loop over all atoms
      DO I=1,NAr
         Ni = NINDX(I,STN) 
         Nf = NINDX(I+1,STN)-1
         
         IF( DEBUG ) WRITE(200,*) '  I  ',I

!
!        Loop over all pairs  
!
         DO N = Ni,Nf
            J = NLIST(N,STN)
            IF( DEBUG ) WRITE(200,*) '   J  ',J
!
!           Find bonds 
!
            IF( J .GT. I ) THEN
               NBl = NBl + 1
               BNDI(NBl,STN) = I
               BNDJ(NBl,STN) = J
               !BTYP(NBl,STN) = TP
               !BVAL(NBl,STN) = VAL
               !BCNST(NBl,STN) = CONSTi 
            ENDIF 
!
!           Find angles 
!        
            DO SN = Ni,Nf
               K = NLIST(SN,STN)
               IF( K .GT. J ) THEN
                  NANGl = NANGl + 1
                  ANGI(NANGl,STN) = I
                  ANGJ(NANGl,STN) = J
                  ANGK(NANGl,STN) = K
                  !   ATYP(NANGl,STN) = TP
                  !   AVAL(NANGl,STN) = VAL
                  !   ACNST(NANGl,STN) = CONSTi 
               ENDIF 

!
!              Find dihedrals
!
               IF( J .GT. I .AND. K.NE.J ) THEN
                  NJi = NINDX(J,STN) 
                  NJf = NINDX(J+1,STN)-1
                  DO NJ = NJi,NJf
                     L =  NLIST(NJ,STN)
                     IF( L.NE.I .AND. L.NE.K  ) THEN
                        NDIHl = NDIHl + 1 
                        DIHI(NDIHl,STN) = I
                        DIHJ(NDIHl,STN) = J 
                        DIHK(NDIHl,STN) = K 
                        DIHL(NDIHl,STN) = L 
                        !DTYP(NDIHl,STN) = TP
                        IF(DEBUG) WRITE(601,*) K,I,J,L 
                     ENDIF
                  ENDDO 
               ENDIF
            ENDDO 
         ENDDO
         NNAB = Nf - Ni + 1
         IF( DEBUG) WRITE(401,*) I,NNAB           
      ENDDO

!     Update global #'s
      NB(STN) = NBl
      NANG(STN) = NANGl
      NDIH(STN) = NDIHl

!     Verbose output
      IF(VERB) WRITE(*,*) ' Bonds found for structure ',
     & STN,' are ',NBl

      END SUBROUTINE get_bonds 
