!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Add velocities to specified structure

      SUBROUTINE ADDVEL(STN)
      USE structure
      USE specify 
!
      IMPLICIT none
!
      INTEGER :: STN,I
!      
      IF( VERB) THEN
         WRITE(6,*)"Adding velocity",VEL(:,STN),"to str",STN
      ENDIF
      DO I=1,NA(STN)
        R1(:,I,STN) =  R1(:,I,STN) + VEL(:,STN)
      ENDDO
      IF( VERB) THEN
         WRITE(6,*)"Final velocity",R1(:,NA(STN),STN),"to str",STN
      ENDIF
!
      RETURN
      END SUBROUTINE ADDVEL
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 10/24/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Center structure

      SUBROUTINE centstr(STN)
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: N,EL,D,STN
      REAL*8  :: SCENT(NDIM),CENT(NDIM)
     &  ,aA,AB,AC,BA,BB,BC,CA,CB,CC 
!     Find coordinates of specified lattice vector
!     Set lattice vectors
      AA = LV(1,1,STN)
      AB = LV(1,2,STN)
      AC = LV(1,3,STN)
      BA = LV(2,1,STN)
      BB = LV(2,2,STN)
      BC = LV(2,3,STN)
      CA = LV(3,1,STN)
      CB = LV(3,2,STN)
      CC = LV(3,3,STN)

      CENT(1)=AA*FCENT(1)+AB*FCENT(2)+AC*FCENT(3)
      CENT(2)=BA*FCENT(1)+BB*FCENT(2)+BC*FCENT(3)
      CENT(3)=CA*FCENT(1)+CB*FCENT(2)+CC*FCENT(3)
!
      CALL MAXMIN(STN)
!     Find center of strcuture
      SCENT(1) = MMN(1,1,STN) + (MMN(2,1,STN)-MMN(1,1,STN) )/2.0d0
      SCENT(2) = MMN(1,2,STN) + (MMN(2,2,STN)-MMN(1,2,STN) )/2.0d0
      SCENT(3) = MMN(1,3,STN) + (MMN(2,3,STN)-MMN(1,3,STN) )/2.0d0
!
c$$$      IF(VERB) THEN
c$$$        WRITE(6,'("Centering structure ",I6)') STN
c$$$        WRITE(6,'("  Max :",3F8.2)') MMN(2,:,STN)
c$$$        WRITE(6,'("  Min :",3F8.2)') MMN(1,:,STN)
c$$$        WRITE(6,'("  Structure center ",3F8.2)') SCENT(:)
c$$$        WRITE(6,'("  will be moved to ",3F8.2)') CENT(:)
c$$$      ENDIF
      DO N=1,NA(STN)
        DO D=1,NDIM
          R0(D,N,STN) =  R0(D,N,STN) - SCENT(D) + CENT(D)
        ENDDO        
      ENDDO 
!
!     Apply PBC's 
      !CALL PBC(STN) 
      CALL MAXMIN(STN)

      END SUBROUTINE centstr

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 5/02/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Center of mass
!
      SUBROUTINE centmolmass(MN,STN)
!
      USE structure
      USE elements 
      USE potential 
      USE specify 
!
      IMPLICIT none 
!
      INTEGER :: I,EL,D,STN,MN,Mi,Mf,N
      REAL*8  :: MS,RMAS(NDIM),TOTMr
!
      TOTMr = 0.0d0
      RMAS(:) = 0.0d0
!
      Mi = MPNT(MN,STN)
      Mf = MPNT(MN+1,STN)-1          
      DO N=Mi,Mf
        I = MOLST(N,STN) 
        MS = AMAS(I,STN ) 
        DO D =1,3        
          RMAS(D) = RMAS(D) + MS*R0(D,I,STN)
        ENDDO
        TOTMr = TOTMr + MS 
      ENDDO
!
      COMAS(1) = RMAS(1)/TOTMr
      COMAS(2) = RMAS(2)/TOTMr
      COMAS(3) = RMAS(3)/TOTMr
!
      IF( VERB ) THEN 
        WRITE(*,*) " Center of mass of ",STN," mol",MN
        WRITE(*,*) " with  ",Mf-Mi+1," atoms "
        WRITE(*,*) COMAS(:) 
      ENDIF 
      END SUBROUTINE centmolmass
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 10/24/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Center structure

      SUBROUTINE centmol(MN,STN)
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: N,EL,D,STN,MN
      REAL*8  :: SCENT(NDIM),CENT(NDIM),FCr(NDIM) 
     &  ,LVr(NDIM,NDIM)
!     Find coordinates of specified lattice vector
!     Set lattice vectors
      LVr(:,:) = LV(:,:,STN)
! 
!     Calcuate shift to specified fractional center
!
      FCr(:) = MCCENT(:,STN)
      CENT(1)=LVr(1,1)*FCr(1)+LVr(2,1)*FCr(2)+LVr(3,1)*FCr(3)
      CENT(2)=LVr(1,2)*FCr(1)+LVr(2,2)*FCr(2)+LVr(3,2)*FCr(3)
      CENT(3)=LVr(1,3)*FCr(1)+LVr(2,3)*FCr(2)+LVr(3,3)*FCr(3)
!
      IF( MLMNCENT( STN)  ) THEN 
!        Use maxmin to center 
        CALL MAXMINMOL(MN,STN)
!       Find center of strcuture
        SCENT(1) = MMN(1,1,STN) + (MMN(2,1,STN)-MMN(1,1,STN) )/2.0d0
        SCENT(2) = MMN(1,2,STN) + (MMN(2,2,STN)-MMN(1,2,STN) )/2.0d0
        SCENT(3) = MMN(1,3,STN) + (MMN(2,3,STN)-MMN(1,3,STN) )/2.0d0
        IF(VERB ) WRITE(*,*) " Using max min of molecule ",MN 
      ELSE
!       Find center of molecule 
        CALL centmolmass(MN,STN) 
        SCENT(:) = COMAS(:)
        IF(VERB ) WRITE(*,*) " Using center of mass of molecule ",MN
      ENDIF
!
      IF(VERB) THEN
        WRITE(6,'("Centering structure ",I6)') STN
        WRITE(6,'("  Structure center ",3F8.2)') SCENT(:)
        WRITE(6,'("  will be moved to ",3F8.2)') CENT(:)
      ENDIF
      DO N=1,NA(STN)
        DO D=1,NDIM
          R0(D,N,STN) =  R0(D,N,STN) - SCENT(D) + CENT(D)
        ENDDO        

      ENDDO 
!     Apply PBC's 
      !CALL PBC(STN) 


!

      END SUBROUTINE centmol
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Remove specified atom types 
      SUBROUTINE RMATOMTYP(STN)
      USE structure
      USE potential 
      USE specify
      USE CONST
      USE elements 

      IMPLICIT none
      INTEGER :: STN,ACNT,Io,A,NATR,EL
      CHARACTER(IDSZ) :: AT,GT,AN,RT
      LOGICAL :: ADD

        NATR = NATPR(STN)
        ACNT = 0
        DO Io=1,NA(STN)
          EL =  ELN(Io,STN) 
          AT =  TYP(Io,STN)
          GT =  GRID(Io,STN)
          RT =  RID(Io,STN)
          AN = ATSYM(EL) 
          ADD = .TRUE.
          DO A=1,NATR
            IF( AT .EQ.RATOMTP(A,STN) ) ADD=.FALSE.
            IF( GT .EQ.RATOMTP(A,STN) ) ADD=.FALSE.
            IF( RT .EQ.RATOMTP(A,STN) ) ADD=.FALSE.
            IF( AN .EQ.RATOMTP(A,STN) ) ADD=.FALSE.
          ENDDO
          IF( ADD ) THEN
            ACNT = ACNT + 1    
            CALL ATPASS(Io,STN,ACNT,STN)
          ENDIF
        ENDDO
!
!
!
        IF( VERB ) THEN
        WRITE(6,*) 'For stn',STN,' of ', NA(STN),', '
     &  ,NA(STN)-ACNT
     &  ,' atoms hav e been removed to create new structure with '
     &  ,ACNT,' aotms'
        ENDIF

        NA(STN) = ACNT

      RETURN
      END SUBROUTINE RMATOMTYP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Cut out molecules

      SUBROUTINE CUTMOL(STN)
      USE specify
      USE structure
      USE potential
      USE const
!
      Implicit none
!
      INTEGER :: STN,ACNT,MN,MNo,MNf,N,Mi,Mf,Io
      REAL*8 :: LVr(NDIM,NDIM),F(NDIM),DV,CFCr(2,NDIM)
      LOGICAL :: ADD(NDIM)
!     Apply PBC's 
      !CALL PBC(STN) 
!

      ACNT = 0
!
!     Loop over all molecules
!   
      MNo =  CMOLN(1,STN) 
      MNf =  CMOLN(1,STN) 
!
      DO MN = MNo,MNf
         Mi = MPNT(MN,STN)
         Mf = MPNT(MN+1,STN)-1          
         DO N=Mi,Mf
            Io = MOLST(N,STN) 
            ACNT = ACNT + 1              
            CALL ATPASS(Io,STN,ACNT,STN)
         ENDDO
      ENDDO
!     Update number of atoms and latice vectors
      NA(STN) = ACNT

!    
!     Verbose
!
      IF(VERB) THEN
        WRITE(*,*) " Molecules " , MNo," - ", MNf," have been removed"
        WRITE(*,*) "   creating a new structure ",STN," of ", NA(STN)
      ENDIF
!
      IF( CLVD(STN) )  CALL CALCVD(STN)
      CALL CALCNELM(STN)
      CALL MAXMIN(STN)
!
      RETURN
      END SUBROUTINE CUTMOL
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Cut out region based on fractional coordinates

      SUBROUTINE CUTFRAC(STN)
      USE specify
      USE structure
      USE potential
      USE const
!
      Implicit none
!
      INTEGER :: STN,Io,ACNT,NAr,I,J
      REAL*8 :: LVr(NDIM,NDIM),F(NDIM),DV,CFCr(2,NDIM)
      LOGICAL :: ADD(NDIM)
!     Apply PBC's 
      !CALL PBC(STN) 
!
      NAr = NA(STN)
      CFCr(:,:) = CUTFC(:,:,STN) 
!     Triclinc PBC's
      LVr(:,:) = LV(:,:,STN)
      ALLOCATE( Fr(NDIM,NAr),MLVr(NDIM) ) 
      CALL fracR(STN,LVr,1)
!
      MLVr(:) = 0.0d0
      DO I=1,NDIM
        DO J=1,NDIM
          MLVr(I) = MLVr(I) + LVr(I,J)*LVr(I,J)
        ENDDO
        MLVr(I) = sqrt( MLVr(I) )
      ENDDO
      IF ( VERB ) THEN
        WRITE(*,*) " Cutting structure",STN
        WRITE(*,*) " at :"
        WRITE(*,*)  CFCr(:,1)
        WRITE(*,*)  CFCr(:,2)
        WRITE(*,*)  CFCr(:,3)
        WRITE(*,*) "   Latice vectors:"
        WRITE(*,'(3F12.4)') LVr(1,:)
        WRITE(*,'(3F12.4)') LVr(2,:)
        WRITE(*,'(3F12.4)') LVr(3,:)
        WRITE(*,*) "   Magnitudes of latice vectors:"
        WRITE(*,'(3F12.4)') MLVr(:)
      ENDIF

      ACNT = 0
      DO Io = 1,NAr
        F(:)  = Fr(:,Io)
        ADD(:) = .FALSE.
      IF(F(1).GT.CFCr(1,1).AND.F(1).LT.CFCr(2,1))ADD(1)=.TRUE.
      IF(F(2).GT.CFCr(1,2).AND.F(2).LT.CFCr(2,2))ADD(2)=.TRUE.
      IF(F(3).GT.CFCr(1,3).AND.F(3).LT.CFCr(2,3))ADD(3)=.TRUE.
        IF( ADD(1).AND.ADD(2).AND.ADD(3) ) THEN
            ACNT = ACNT + 1              
            CALL ATPASS(Io,STN,ACNT,STN)
        ENDIF
      ENDDO
!     Update number of atoms and latice vectors
      NA(STN) = ACNT
      DO I=1,NDIM
        DV =  (CFCr(2,I)-CFCr(1,I))
        DO J=1,NDIM
          LVr(I,J) = DV*LVr(I,J)
        ENDDO
      ENDDO
      LV(:,:,STN)  =LVr(:,:) 
!
      IF ( VERB ) THEN
        WRITE(*,*) "   Updated latice vectors:"
        WRITE(*,'(3F12.4)') LVr(1,:)
        WRITE(*,'(3F12.4)') LVr(2,:)
        WRITE(*,'(3F12.4)') LVr(3,:)
      ENDIF
!
      DEALLOCATE( Fr,MLVr )
!     Apply PBC's 
      CALL PBC(STN) 
!     Update properties 
        IF( CLLC(STN) )  CALL CALCLC(STN)
        IF( CLLV(STN) )  CALL CALCLV(STN)
        IF( CLMAS(STN) )  CALL CALCMAS(STN) 
        IF( CLVD(STN) )  CALL CALCVD(STN)
        CALL CALCNELM(STN)
        CALL MAXMIN(STN)
!
      RETURN
      END SUBROUTINE CUTFRAC      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 12/7/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Center structure
!
      SUBROUTINE pbc(STN)
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: STN,I,D,NAr,J 
      REAL*8  :: F(3),FT,LVr(NDIM,NDIM)
!     Triclinc PBC's
      NAr = NA(STN)
      LVr(:,:) = LV(:,:,STN)
      ALLOCATE( Fr(NDIM,NAr) ) 
      CALL fracR(STN,LVr,1)

      DO I=1,NA(STN)
        F(:) = Fr(:,I)
!       Reduce fractional coordinates to less than 1


        DO D=1,NDIM 
          FT = F(D)

! debug 
!           WRITE(990,*) D,FT
!  end debug 


          IF( PBCS(D,STN) )  F(D) = FT - ANINT( FT )
! debug 
!           WRITE(990,*) D,F(D),ANINT(FT)
!  end debug 

        ENDDO
!       Get new coordinate based on fractional positions
        R0(1,I,STN) = LVr(1,1)*F(1) + LVr(1,2)*F(2)+LVr(1,3)*F(3)
        R0(2,I,STN) = LVr(2,1)*F(1) + LVr(2,2)*F(2)+LVr(2,3)*F(3)
        R0(3,I,STN) = LVr(3,1)*F(1) + LVr(3,2)*F(2)+LVr(3,3)*F(3)
      ENDDO
      DEALLOCATE( Fr )
!
      RETURN
      END SUBROUTINE pbc

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Cut out region based on fractional coordinates

      SUBROUTINE FIXREGION(STN)
      USE specify
      USE structure
      USE potential
      USE const
!
      Implicit none
!
      INTEGER :: STN,Io,CNT,NAr,I,J
      REAL*8 :: LVr(NDIM,NDIM),F(NDIM),DV,CFCr(2,NDIM)
      LOGICAL :: ADD(NDIM)
!
 
      NAr = NA(STN)
      CFCr(:,:) = RFIXC(:,:,STN) 
!     Triclinc PBC's
      LVr(:,:) = LV(:,:,STN)
      ALLOCATE( Fr(NDIM,NAr) ) 
      CALL fracR(STN,LVr,1)
!
      IF ( VERB ) THEN
        WRITE(*,*) " Fixing region ",STN
        WRITE(*,*) " at :"
        WRITE(*,*)  CFCr(:,1)
        WRITE(*,*)  CFCr(:,2)
        WRITE(*,*)  CFCr(:,3)
        WRITE(*,*) "   Latice vectors:"
        WRITE(*,'(3F12.4)') LVr(1,:)
        WRITE(*,'(3F12.4)') LVr(2,:)
        WRITE(*,'(3F12.4)') LVr(3,:)
      ENDIF

      CNT = 0 
      DO Io = 1,NAr
        F(:)  = Fr(:,Io)
        ADD(:) = .FALSE.
      IF(F(1).GT.CFCr(1,1).AND.F(1).LT.CFCr(2,1))ADD(1)=.TRUE.
      IF(F(2).GT.CFCr(1,2).AND.F(2).LT.CFCr(2,2))ADD(2)=.TRUE.
      IF(F(3).GT.CFCr(1,3).AND.F(3).LT.CFCr(2,3))ADD(3)=.TRUE.
        IF( ADD(1).AND.ADD(2).AND.ADD(3) ) THEN
             CNT = CNT + 1
             FTAG(:,Io,STN) = FRTAG(:,STN)
        ENDIF
      ENDDO

      DEALLOCATE( Fr )
      IF(VERB) THEN
         WRITE(*,*) " Number of fix tags updated ",CNT
      ENDIF
!
      RETURN
      END SUBROUTINE FIXREGION      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 10/24/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Calculate maximum and minimum of molecule 

      SUBROUTINE MAXMINMOL(MN,STN)
      USE structure
      USE specify 
!
      IMPLICIT none
!
      INTEGER :: I,D,STN,MN,Mi,Mf,N
!      
      MMN(1,:,STN) = (/1.d16,1.0d16,1.0d16/)
      MMN(2,:,STN) = (/-1.0d16,-1.0d16,-1.0d16/)
!     
      Mi = MPNT(MN,STN)
      Mf = MPNT(MN+1,STN)-1          
      DO N=Mi,Mf
        I = MOLST(N,STN) 
       DO D =1,3        
          IF ( R0(D,I,STN).LT.MMN(1,D,STN)) MMN(1,D,STN)=R0(D,I,STN)
          IF ( R0(D,I,STN).GT.MMN(2,D,STN)) MMN(2,D,STN)=R0(D,I,STN)
        ENDDO
      ENDDO
!
      END SUBROUTINE MAXMINMOL 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!    Copy region based on fractional coordinates
!
      SUBROUTINE CPFAC(STN)
      USE specify
      USE structure
      USE potential
      USE build 
      USE const
!
      Implicit none
!
      INTEGER :: STN,Io,ACNT,NAr,I,J,STNi 
      REAL*8 :: LVr(NDIM,NDIM),F(NDIM),DV,CFCr(2,NDIM),R(NDIM)
      LOGICAL :: ADD(NDIM)
!     Copy region of one structure to new structure 
      STNi = STCPN(STN) 
!
!     Apply PBC's 
      !CALL PBC(STN) 
!
      NAr = NA(STN)
      CFCr(1:2,1:3) = RCPRG(1:2,1:3,STN) 
!     Triclinc PBC's
      LVr(:,:) = LV(:,:,STN)
      ALLOCATE( Fr(NDIM,NAr),MLVr(NDIM) ) 
      CALL fracR(STN,LVr,1)
!
      MLVr(:) = 0.0d0
      DO I=1,NDIM
        DO J=1,NDIM
          MLVr(I) = MLVr(I) + LVr(I,J)*LVr(I,J)
        ENDDO
        MLVr(I) = sqrt( MLVr(I) )
      ENDDO
      IF ( VERB ) THEN
        WRITE(*,*) " Copying structure",STN," to ",STNi
        WRITE(*,*) " at :"
        WRITE(*,*)  CFCr(:,1)
        WRITE(*,*)  CFCr(:,2)
        WRITE(*,*)  CFCr(:,3)
        WRITE(*,*) "and displacing ",CPRDS(:,STN) 
        WRITE(*,*) "   Latice vectors:"
        WRITE(*,'(3F12.4)') LVr(1,:)
        WRITE(*,'(3F12.4)') LVr(2,:)
        WRITE(*,'(3F12.4)') LVr(3,:)
        WRITE(*,*) "   Magnitudes of latice vectors:"
        WRITE(*,'(3F12.4)') MLVr(:)
      ENDIF

      ACNT = 0 
      DO Io = 1,NAr
        F(:)  = Fr(:,Io)
        ADD(:) = .FALSE.
      IF(F(1).GT.CFCr(1,1).AND.F(1).LT.CFCr(2,1))ADD(1)=.TRUE.
      IF(F(2).GT.CFCr(1,2).AND.F(2).LT.CFCr(2,2))ADD(2)=.TRUE.
      IF(F(3).GT.CFCr(1,3).AND.F(3).LT.CFCr(2,3))ADD(3)=.TRUE.
        IF( ADD(1).AND.ADD(2).AND.ADD(3) ) THEN
            ACNT = ACNT + 1              
            CALL ATPASS(Io,STN,ACNT,STNi)
             F(1) = F(1) + CPRDS(1,STN)
             F(2) = F(2) + CPRDS(2,STN)
             F(3) = F(3) + CPRDS(3,STN)
             CALL REALF(LVr,F,R)
             R0(:,ACNT,STNi) = R(:)
!
        ENDIF
      ENDDO
!     Update number of atoms and latice vectors
      NA(STNi) = ACNT
      LV(:,:,STNi)  =LVr(:,:) 
!
      DEALLOCATE( Fr,MLVr )
!     Apply PBC's 
      CALL PBC(STN) 
!     Update properties 
      IF( MTOAN(STNi) )    CALL FINDELN(STNi)
      IF( CLLC(STNi) )  CALL CALCLC(STNi)
      IF( CLLV(STNi) )  CALL CALCLV(STNi)
      IF( CLMAS(STNi) )  CALL CALCMAS(STNi) 
      IF( CLVD(STNi) )  CALL CALCVD(STNi)
      CALL CALCNELM(STNi)
      CALL MAXMIN(STNi)
      CALL WRITE_STRP(STN)
!
      RETURN
      END SUBROUTINE     CPFAC  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 05/19/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      Asign charge groups
!
      SUBROUTINE CHRGNB(STN)
      USE specify
      USE structure
      USE build
!
      IMPLICIT none
!
      INTEGER :: STN,I,J,CHG,NCNT,No,Nf,N,NAr 
      REAL*8 :: RSQ(NDIM),DR(NDIM),Ri(NDIM)
      LOGICAL :: SQNB(NDIM)
      LOGICAL, ALLOCATABLE :: CHFIND(:) 
      
!     Locals
      NAr = NA(STN) 
!     Initialize all atoms as unfound
      ALLOCATE( CHFIND(NAr) )
      CHFIND(:) = .TRUE.
!
      NCNT = 0
      CHG   = 1
!
      DO I=1,NAr
        IF( CHFIND(I) ) THEN
          NCNT = NCNT + 1
          IF( NCNT.GT.CHGMX) THEN
            NCNT = 0
            CHG = CHG + 1
          ENDIF
          CHN(I,STN) = CHG
          CHFIND(I) = .FALSE.
          No = NINDX(I,STN) 
          Nf = NINDX(I+1,STN) -1 
          DO N=No,Nf
            J=NLIST(N,STN)
            NCNT = NCNT + 1
            IF( NCNT.GT.CHGMX) THEN
               NCNT = 0
               CHG = CHG + 1
            ENDIF
            CHN(J,STN) = CHG
            CHFIND(J) = .FALSE.
          ENDDO ! N=No,Nf
        ENDIF ! CHFIND(I)
      ENDDO ! I   
!
      IF( VERB) THEN
      WRITE(6,*) ' Structure ',STN,' has ',CHG,' charge groups'
      ENDIF 
      DEALLOCATE( CHFIND) 
      RETURN
      END SUBROUTINE CHRGNB
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Rename 
!
      SUBROUTINE rename_TYP(STN)
      USE specify
      USE structure 
!
      IMPLICIT none
      INTEGER :: STN,I,N
!
      DO I=1,NA(STN)
        DO N=1,RENM(STN)
          IF( GRID(I,STN).EQ.RNIDI(N,STN)) GRID(I,STN) = RNIDF(N,STN)
          IF( TYP(I,STN).EQ.RNIDI(N,STN)) TYP(I,STN) = RNIDF(N,STN)
          IF( RID(I,STN).EQ.RNIDI(N,STN)) RID(I,STN) = RNIDF(N,STN)
        ENDDO 
      ENDDO
!
      END SUBROUTINE rename_TYP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Hydrogen terminate carbon bonds
!
      SUBROUTINE hterm(STN)
      USE specify
      USE structure 
      USE elements
!
      IMPLICIT none
      INTEGER :: STN,I,NB,Ni,Nf,N,ACNT,J,EL,NAr
      REAL*8 :: RH(NDIM),Ri(NDIM),RNSQ,RMAG,BL(NELM)
     &  ,RNBL(NDIM)
!
      BL(:) = 1.5d0
      BL(1) = 1.1d0
      SDIM = 3      ! Set surface dimension in Z
!
      NAr = NA(STN)
      DO I=1,NAr
        IF( ELN(I,STN) .EQ. 6) THEN
          Ri(:) = R0(:,I,STN)
          Ni = NINDX(I,STN)
          Nf = NINDX(I+1,STN) -1
          NB = Nf - Ni + 1
          RH(:) =  0.0d0
          WRITE(*,*) I,NB
          DO N =Ni,Nf
            J = NLIST(N,STN)
            RH(:) = RH(:) + Ri(:) - R0(:,J,STN) 
          ENDDO
          RNSQ = RH(1)*RH(1) + RH(2)*RH(2)+RH(3)*RH(3)
!         Does not work due pbc's
!          RMAG = SQRT(RNSQ)
!          RNBL(:) = RH(:)/RMAG*BL(1)
          RNBL(:) = 0.0d0
          RMAG = SQRT( RH(SDIM)*RH(SDIM))
          RNBL(SDIM) = RH(SDIM)/RMAG*BL(1)
!
          IF( NB.LT.4) THEN
            ACNT = ACNT + 1
            R0(:,ACNT,STN) = Ri(:) + RNBL(:)
            CALL ATPASS(I,STN,ACNT,STN)             
            ELN(ACNT,STN) = 1 
          ENDIF
        ENDIF
      ENDDO
      NA(STN) = ACNT
      END SUBROUTINE hterm
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! shift 
!
      SUBROUTINE SHIFT(STN)
      USE specify
      USE structure 
!
      IMPLICIT none
      INTEGER :: STN,I,N,D
      REAL*8 :: VEC(NDIM)
!
      VEC(:) = SHVEC(:,STN)
      IF(VERB) WRITE(6,*)"Structure",STN,'will be shifted',VEC(:)
!
      DO I=1,NA(STN)
        DO D=1,NDIM
          R0(D,I,STN) = R0(D,I,STN) -  VEC(D)
        ENDDO
      ENDDO
!
      END SUBROUTINE SHIFT 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! shift max to zero  
!
      SUBROUTINE ZSHIFT(STN)
      USE specify
      USE structure 
!
      IMPLICIT none
      INTEGER :: STN,I,N,D
      REAL*8 :: SHIFT

      D = SZDIR(STN)
      CALL MAXMIN(STN)
      SHIFT = MMN(2,D,STN)
      IF(VERB) WRITE(6,*)"Structure",STN,'will be shifted',SHIFT
!
      DO I=1,NA(STN)
        R0(D,I,STN) = R0(D,I,STN) - SHIFT 
      ENDDO
!
      END SUBROUTINE ZSHIFT 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Find molecule neighbors 
      SUBROUTINE WHOLEMOL(STN)
      USE const
      USE structure
      USE specify 
!
      IMPLICIT none
!
      INTEGER :: STN,MN,Mi,Mf,N,MNB,MNBo,Ni,Nf,Nb,I,J
     & ,MNn,MICK,MCK 
      LOGICAL :: ADDM 
      CHARACTER(ERSZ) :: ESTAT
!
!     Set initial 
!
      MNB = 0 



      END SUBROUTINE  WHOLEMOL 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Center structure
!
      SUBROUTINE molpbc(STN)
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: STN,I,D,NAr,J,MN,Mi,Mf,N
      REAL*8  :: F(NDIM),FT(NDIM),LVr(NDIM,NDIM),R(NDIM)
     & ,Fi
      LOGICAL  :: TRANS 
!     Triclinc PBC's
      NAr = NA(STN)
      LVr(:,:) = LV(:,:,STN)
      ALLOCATE( Fr(NDIM,NAr) ) 
      CALL fracR(STN,LVr,1)
!
!     Get center of mass
!
      CALL centmass(STN)
!
!     Loop over all molecules
!
      DO MN=1,MOLCNT(STN)
        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1          
        R(:) = MCOMS(:,MN,STN) 
        CALL SPFRAC(LVr,F,R)
        TRANS = .FALSE.  
!
!       Find out if comas outside box
!
        DO D=1,NDIM 
           IF( PBCS(D,STN).AND.ANINT( F(D)).NE.0 ) TRANS = .TRUE. 
           FT(D) = ANINT( F(D) )
        ENDDO

        CALL REALF(LVr,FT,R)

!
!       If comas ouside box shift back into box 
!
        IF( TRANS ) THEN
!
!     Loop over all atoms in molecule 
!
          DO N=Mi,Mf
            I = MOLST(N,STN) 
            DO D=1,NDIM
              Fi = Fr(D,I)
              F(D) = Fi -  FT(D)
            ENDDO
!
!           Get new coordinate based on fractional positions
!
            CALL REALF(LVr,F,R)
            R0(:,I,STN) = R(:)

          ENDDO 
        ENDIF !       If comas ouside box shift back into box 
      ENDDO !     Loop over all molecules

      DEALLOCATE( Fr )
!
      RETURN
      END SUBROUTINE molpbc
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Cut out neighbors of molecule 
!
      SUBROUTINE molnbcut(STN)
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: STN,MN,NBS,Ni,Nf,I,N,Nb
     & ,ACNT,Mi,Mf,NBcnt,J,JMi,JMf,JN,JI
!
!     Initialize count 
!
      ACNT = 0
!
      MN  = CMOLN(1,STN) 
      NBS = CMOLN(2,STN) 

!
!     Verbose
!
      IF(VERB) WRITE(6,*) 'Cuttig ', NBS ,' around ',MN
!
!     Add atoms of molecule i to updated list 
!
      Mi = MPNT(MN,STN)
      Mf = MPNT(MN+1,STN)-1          
      DO N=Mi,Mf
         I = MOLST(N,STN) 
         ACNT = ACNT + 1              
         CALL ATPASS(I,STN,ACNT,STN)
      ENDDO
!
      CALL BUILDMNB(STN)
!      CALL CENTMOL(MN,STN)
!      CALL molPBC(STN)          !Apply PBC's 
!
!     Loop over molecule neigbor list of atom 
!
      NBcnt = 0
!
      !DO WHILE ( NBcnt .LT. NBS )
         
         Ni = NINDXML(MN,STN)           
         Nf = NINDXML(MN+1,STN)-1

         IF(VERB)THEN
            WRITE(6,*)' Mol ',MN,' has ',Nf-Ni-1,' neighbros',Ni,Nf
         ENDIF

         DO Nb=Ni,Nf
            NBcnt = NBcnt + 1

            IF( DEBUG) WRITE(162,*) Nb,NBcnt 

            IF( NBcnt .LE. NBS ) THEN 
               J = NLISTML(Nb,STN)
               IF(VERB) WRITE(6,*) ' adding mol ',J,' count ',NBcnt 
               JMi = MPNT(J,STN)
               JMf = MPNT(J+1,STN)-1          
               DO JN=JMi,JMf
                 JI = MOLST(JN,STN) 
                 ACNT = ACNT + 1              
                 CALL ATPASS(JI,STN,ACNT,STN)
              ENDDO
            ENDIF 
         ENDDO
      !ENDDO
!
!     Update number of atoms and latice vectors
!
      NA(STN) = ACNT
!
      IF(VERB) WRITE(6,*) ' Mol nb cut finished with ', ACNT,' atoms'
!
      IF( CLVD(STN) )  CALL CALCVD(STN)
      CALL CALCNELM(STN)
      CALL MAXMIN(STN)

      RETURN
      END SUBROUTINE molnbcut 

