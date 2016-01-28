!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Programer:
!     Travis Kemper
!     Department of Materials Sxience and Engineering
!     University of FLorida
!     traviskemper@ufl.edu
!
!     Version 1.0  4/11/11 T. W. Kemper
!     Version 2.0  5/02/11 T. W. Kemper
!     GaTech
!     Version 3.0 12/06/11 T. W. Kemper
!     Version 4.0 12/07/11 T. W. Kemper
!       :add vacuum option
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Randomly distrubute molecules
!
      SUBROUTINE addmol 
!
      USE structure
      USE specify 
      USE potential 
      USE build 
      USE elements 
! 
      IMPLICIT none
!
      INTEGER :: N,M,IM,NG,II,J,K,STNo,STNi,So,Si,D,RNREF,MREF,I
     & ,LDIM(NDIM),NMAX,MN,ROTCNT,MINEL(2),RSTN
      INTEGER :: ATEMPS,NATOM,SFIX,IE,NEWCNT,NGTOT,ELii,ELim 
      INTEGER, ALLOCATABLE :: NREF(:)
      REAL*8 :: RNi,VAC(NDIM,NDIM),RTIM,DRSQ,DELR(NDIM),MINOV,SZSQ(NDIM)
     &   ,DRSQI, X,Y,Z,AA,AB,AC,BA,BB,BC,CA,CB,CC,RTL(3),RIL(3),F(3),FT
     &   ,Vj,Vk,MD,MDSQ,MFRAC,SZF,SZSTR(NDIM),SHIFT,MINBFSQ
     &   ,MINDR(NDIM),DELSQ(NDIM),CBUF,VDBUF,FRZL,ANa,ANb
     &   ,MOL_HIEGHT,MAXH
      REAL*8, ALLOCATABLE :: RT(:,:),RPNT(:,:)
      INTEGER, ALLOCATABLE :: MCNT(:)
      LOGICAL :: OVERLAP,MKNEW,ADDMOLi,NWROT
!
      NMAX = 0 
      MAXH = 0.0d0 
!
!     Recenter strucutre at origin 
!
      FCENT(:) = 0.0d0 
      ALLOCATE( NREF(STOTo) )
!
!     Initalize box as last random structure lattice vectors
!     
      NGTOT = 0 
      DO STNo=1,STOTo
        STNi = STNo + STOTo 
        N = NA(STNo)*NDEP(STNo)
!
!       Use another structure as box size 
!
        IF( UDBOX(STNo) ) THEN 
          RSTN = RSTBOX(STNo) 
          LV(:,:,STNi) = LV(:,:,RSTN)
          VAC(:,:) = LV(:,:,RSTN)
!         Set  intial hieght to max of reference structure 
          CALL MAXMIN(RSTN)
          MAXH  = MMN(2,SDIM,RSTN) 
        ENDIF
        IF( VERB ) THEN
            WRITE(*,*) ' Randbox lattice vectors',STNi
            WRITE(*,*) VAC(1,:)
            WRITE(*,*) VAC(2,:)
            WRITE(*,*) VAC(3,:)
            WRITE(*,*) ''
        ENDIF
!
!       If molecules need to be added increase NMAX
!
        IF( ADDDMOL(STNo) ) THEN
          CALL MOLCHECK(STNo)
          IF( N.GT.NMAX ) NMAX=N
          NGTOT = NGTOT + NDEP(STNo) 
        ENDIF 
!
      ENDDO
      IF(VERB) THEN
        WRITE(*,'(" Total number for molecules to add ",I6)')NGTOT
      ENDIF

!
!     Allocate and intialize
!

      ALLOCATE( RT(NMAX,NDIM),RPNT(NGTOT,NDIM),MCNT(STOT) )
      MCNT(:) = 0
      NEWCNT = 0
      RNREF = 0
      MREF = 0
      NREF(:) = 0
      MN   = 0 
!
!     Set final values of # of atoms for random atoms
      DO STNo=1,STOTo
         IF( ADDDMOL(STNo) ) THEN
           STNi = STNo + STOTo
           NA(STNi)   = 0
           NB(STNi)   = 0
!           NPR(STNi)  = 0
           NANG(STNi) = 0
           NDIH(STNi) = 0 
!          Local mol count 
           MCNT(STNi)  = 0
         ENDIF
      ENDDO
!
      ADDMOLi = .TRUE.
      DO WHILE ( ADDMOLi )
        ADDMOLi = .FALSE.
!       Loop over all strucutres
        DO STNo=1,STOTo
!        Only use the ones with randimize turned on
         IF( ADDDMOL(STNo) ) THEN
           STNi = STNo + STOTo
           CALL CENTERMASS(STNo)
           IF ( MCNT(STNi).LT.NDEP(STNo) ) THEN
            ADDMOLi = .TRUE.
            OVERLAP = .TRUE.
            ATEMPS = 0
            MN = MN + 1 
!
            IF ( MN.GT. NGTOT ) THEN
                WRITE(*,*) " Mol count ",MN
     &          ," GT  total # to be added ",NGTOT
                WRITE(*,*)  ' STNo ',  STNo              
                WRITE(*,*)  ' MCNT(STNi) ',    MCNT(STNi)            
!                WRITE(*,*)  '  ',                 
                   WRITE(*,*) " For attempt ",ATEMPS," for mol ",MN 
                   WRITE(*,*) " Rotation atemps ",ROTCNT 
                   WRITE(*,*) "  MKNEW ", MKNEW
                   WRITE(*,*) "   ADDMOLi ", ADDMOLi
                   WRITE(*,*) "    NEWCNT  ", NEWCNT
                   WRITE(*,*) "     OVERLAP ", OVERLAP
                   WRITE(*,*) "       NWROT ", NWROT
                STOP 
            ENDIF
!
!           Rotate structure 
!
            IF( RANROT(STNo) ) THEN
!                  Get random #'s 
                   SEED = SEED*SEED
                   SEED = SEED - SEED/2 
                   IF( VERB) WRITE(*,*) ' SEED = ',SEED 
                   ANa =  ran(seed)*( 2*PI )
                   SEED = SEED - SEED/2 
                   IF( VERB) WRITE(*,*) ' SEED = ',SEED 
                   ANb = ran(seed)*( 2*PI )
                   IF(VERB) WRITE(*,*) " angle theta,phi",ANa,ANb
                   CALL CENTSTR(STNo)    
                   CALL ROTATE(STNo,ANa,ANb)
            ENDIF 
!
!           Generate random point
!
            RPNT(MN,:) = 0.0d0
            SEED = SEED*SEED
            DO M =1,NDIM
                SEED = SEED - SEED/2 
                RNi =  ran(seed) 
                F(M) = RNi
                DO D=1,NDIM
                   IF(M.NE.SDIM) RPNT(MN,D) = RPNT(MN,D) +  RNi*VAC(M,D) 
                ENDDO
            ENDDO
!
!           Add hieght and buffer distance
!
            IF( UDBOX(STNo) ) THEN 
               RSTN = RSTBOX(STNo)
               MAXH  = MMN(2,SDIM,RSTN)
               CALL MAXMIN(STNo)
               RPNT(MN,SDIM)=MAXH-MMN(1,SDIM,STNo)+DEPBF(STNo)
            ENDIF
!
            IF( VERB) THEN
              WRITE(*,*)' MAXH  MMN(1,SDIM,RSTN) DEPBF(STNo) MOL_HIEGHT'
              WRITE(*,*)MAXH,MMN(1,SDIM,STNo),DEPBF(STNo),RPNT(MN,SDIM)
              WRITE(*,'("Random point",3F10.3)') RPNT(MN,:)  
            ENDIF
!
!
!           Shift structure to random point 
!   
            DO IM =1,NA(STNo)
                   RT(IM,:) = R0(:,IM,STNo) + RPNT(MN,:) 
                   R0(:,IM,STNo) = RT(IM,:) 
            ENDDO 
!
!           Add to global  atoms        
!
            DO I = 1,NA(STNo)
                N = NA(STNi) + I
                CALL ATPASS(I,STNo,N,STNi)       
            ENDDO
!
!           Update total number of atoms NP
!
            NA(STNi)   = NA(STNi)   + NA(STNo)
            NB(STNi)   = NB(STNi)   + NB(STNo)
!             NPR(STNi)  = NPR(STNi)  + NPR(STNo)
            NANG(STNi) = NANG(STNi) + NANG(STNo)
            NDIH(STNi) = NDIH(STNi) + NDIH(STNo)
!
!           Save last numbers of some counts as a reference
!
            NREF(STNo)  = NA(STNi)
            RNREF = RN(NREF(STNo),STNi)    
            MREF  = MOLN(NREF(STNo),STNi)
            MCNT(STNi)  = MCNT(STNi)  + 1
            CALL MAXMIN(STNi)
            IF( VERB) THEN
                 WRITE(*,*) "Added ",MCNT(STNi)," of str#",STNi
                WRITE(*,*) " min of added mole ",MMN(1,:,STNi)
                WRITE(*,*) " max of added mole ",MMN(2,:,STNi)
                 WRITE(6,*) ' which now has ',NA(STNi)
     &           ,' atoms from',NA(STNo)
                WRITE(6,'(A,2I6)') ' Ref :',RNREF,MREF
                WRITE(6,*) 
            ENDIF
           ENDIF    ! NCNT(STN).LT.NGAS(STN) 
         ENDIF      ! ADDDM 
        ENDDO       ! STN
      ENDDO        ! ADDMOLi 

      IF( VERB ) THEN
        WRITE(6,*) 'Final vacuum space:'
        WRITE(6,'(3F8.2)') VAC(1,:)
        WRITE(6,'(3F8.2)') VAC(2,:)
        WRITE(6,'(3F8.2)') VAC(3,:)
        CALL MAXMIN(STNi)
        WRITE(*,*) " min of added mole ",MMN(1,:,STNi)
        WRITE(*,*) " max of added mole ",MMN(2,:,STNi)
      ENDIF
!
      !
      DEALLOCATE(RT,NREF,RPNT,MCNT )
!
      RETURN

      END SUBROUTINE addmol

!
!
! 
      SUBROUTINE ADDTEMP(STN)
      USE specify
      USE structure
      
!
      IMPLICIT none 
!
      INTEGER I,STN,MN
      REAL*8 :: GASTEMP,KE,MVELSQ,MVEL
!
!     Update molecule list 
!
      IF( MOLCH(STN)  ) CALL MOLCHECK(STN) 
!
!     Get molmasses 
!
      CALL centmass(STN)
!
      GASTEMP = DEPTEMP(STN)
      DO MN = 1,MOLCNT(STN)
     
        KE = 3.0d0*BOLZeV*GASTEMP/2.d0
        MVELSQ = 2.0d0*KE/ECONV/MOLMAS(MN,STN)
        MVEL = -1.0d0*SQRT( MVELSQ )

        MVEL = -1.0d0

        DO I=1,NA(STN)
          R1(SDIM,I,STN) =  MVEL
        ENDDO
      ENDDO
        
      IF( VERB) THEN
         WRITE(6,*)"Final velocity",R1(:,NA(STN),STN),"to str",STN
      ENDIF

      END SUBROUTINE ADDTEMP
 
