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
      SUBROUTINE mkrandom 
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
     &   ,MAXH,DELTS
      REAL*8, ALLOCATABLE :: RT(:,:),RPNT(:,:)
      INTEGER, ALLOCATABLE :: MCNT(:)
      LOGICAL :: OVERLAP,MKNEW,ADDMOL,NWROT
!     
      SZF = 1/8.0d0
!
      NMAX = 0 
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
        N = NA(STNo)*NGAS(STNo)
        IF( MKRAND(STNo) ) THEN
          CALL MOLCHECK(STNo)
          IF( N.GT.NMAX ) NMAX=N
          VAC(:,:) =  LV(:,:,STNi)
          SFIX  = STNo
          NGTOT = NGTOT + NGAS(STNo) 
        ENDIF 
      ENDDO
      IF(VERB) THEN
        WRITE(*,'(" Total number for molecules to add ",I6)')NGTOT
      ENDIF
      DO STNo=1,STOTo 
!       Check if fixed unit cell is used
        STNi = STNo + STOTo 
        IF( FIXLV(STNo) ) THEN
          SFIX = STNo
          IF( VERB ) THEN
            WRITE(*,*) "Using lattice vectors ",STNo ,STNi
            WRITE(*,*) " except ",SDIM
          ENDIF
          DO M =1,NDIM
            DO D=1,NDIM
              IF(M.NE.SDIM) VAC(M,D) =  LV(M,D,STNi)
            ENDDO
          ENDDO
          IF( VERB ) THEN
            WRITE(*,*) ' Fix lattice vectors'
            WRITE(*,*) VAC(1,:)
            WRITE(*,*) VAC(2,:)
            WRITE(*,*) VAC(3,:)
            WRITE(*,*) ''
          ENDIF
        ENDIF
        IF( RANBOX(STNo) ) THEN
!         Set vacuum box to maximum and minimum of specified structure
          CALL MAXMIN(STNi)
          IF( FIXLV(STNo) ) THEN
            VAC(SDIM,SDIM) = MMN(2,SDIM,STNi) - MMN(1,SDIM,STNi)
            IF(VERB) THEN
              WRITE(6,*) " Setting vacuum dimension",SDIM
              WRITE(6,*) " to max/min of structure",STNi
              WRITE(6,*)MMN(2,SDIM,STNi),MMN(1,SDIM,STNi),VAC(SDIM,SDIM)
            ENDIF
          ELSE
            VAC(:,:) = 0
            VAC(1,1) = MMN(2,1,STNi) - MMN(1,1,STNi)
            VAC(2,2) = MMN(2,2,STNi) - MMN(1,2,STNi)
            VAC(3,3) = MMN(2,3,STNi) - MMN(1,3,STNi)
          ENDIF
!          Not sure 
!          DO I=1,NA(STNi)
!            R0i(STN,I,SDIM) = R0i(STN,I,SDIM) - MMNi(STN,1,SDIM)
!          ENDDO
          IF( VERB ) THEN
            WRITE(*,*) ' Randbox lattice vectors',STNi
            WRITE(*,*) VAC(1,:)
            WRITE(*,*) VAC(2,:)
            WRITE(*,*) VAC(3,:)
            WRITE(*,*) ''
          ENDIF
        ENDIF 
!
!       Use another structure as box size 
!
        IF( UDBOX(STNo) ) THEN 
             RSTN = RSTBOX(STNo) + STOTo
!        
!            Add hieght and buffer distance
!
             CALL MAXMIN(RSTN)
             MAXH = MMN(2,SDIM,RSTN) 
             VAC(:,:) = LV(:,:,RSTN)
        ENDIF
!
!       Reset box size 
        DO D =1,NDIM
            IF( RANFRX(D,STNo) ) THEN
              FRZL = RANFRSZ(D,STNo) 
              IF( VERB ) WRITE(*,*) "Resetting box in dim",D
     &        , "to",FRZL,STNo
              VAC(D,D) = RANFRSZ(D,STNo)  

            ENDIF
        ENDDO
!
!       Set molecule position for depostion 
!          in that molecule is placed beyond cut off 
!
        
      ENDDO


      IF(VERB)THEN
        DO STNo=1,STOTo
          IF( MKRAND(STNo).OR. RANBOX(STNo) ) THEN
            WRITE(*,'(A,I6)') 'Structure ',STNo
     &      ," will be used in neighbor check"
            IF(MKRAND(STNo)) THEN
              WRITE(6,*) NGAS(STNo)," of str",STNo,"will be added"
            ENDIF
          ENDIF
        ENDDO
        WRITE(*,'(3F12.4)') VAC(1,:)
        WRITE(*,'(3F12.4)') VAC(2,:)
        WRITE(*,'(3F12.4)') VAC(3,:)
        WRITE(*,'(A,F8.2)') 'Using min buffer of',SQRT(MBUF)
      ENDIF

      ALLOCATE( RT(NMAX,NDIM),RPNT(NGTOT,NDIM),MCNT(STOT) )
      MCNT(:) = 0
!     Triclinc PBC's
      AA = VAC(1,1)
      AB = VAC(1,2)
      AC = VAC(1,3)
      BA = VAC(2,1)
      BB = VAC(2,2)
      BC = VAC(2,3)
      CA = VAC(3,1)
      CB = VAC(3,2)
      CC = VAC(3,3)
!
!     Record original vacuum volume  
!
      Vk = VAC(1,1)*VAC(2,2)*VAC(3,3)
      Vj = Vk
      MKNEW = .TRUE.
      NEWCNT = 0
      DO WHILE ( MKNEW ) 
       MKNEW = .FALSE.
       RNREF = 0
       MREF = 0
       NREF(:) = 0
       MN   = 0 
!
!      Set final values of # of atoms for random atoms
       DO STNo=1,STOTo
         IF( MKRAND(STNo) ) THEN
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
       ADDMOL = .TRUE.
       DO WHILE ( ADDMOL )
        ADDMOL = .FALSE.
!       Loop over all strucutres
        DO STNo=1,STOTo
!        Only use the ones with randimize turned on
         IF( MKRAND(STNo) ) THEN
           STNi = STNo + STOTo


           CALL CENTERMASS(STNo)
           IF ( MCNT(STNi).LT.NGAS(STNo) ) THEN
            ADDMOL = .TRUE.
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
                   WRITE(*,*) "   ADDMOL ", ADDMOL
                   WRITE(*,*) "    NEWCNT  ", NEWCNT
                   WRITE(*,*) "     OVERLAP ", OVERLAP
                   WRITE(*,*) "       NWROT ", NWROT
                STOP 
            ENDIF
            DO WHILE ( OVERLAP )
               OVERLAP = .FALSE.
               ATEMPS = ATEMPS + 1
               MINOV  = AA*BB*CC
               MINDR(:) = MINOV
!              Generate random point
               RPNT(MN,:) = 0.0d0
               SEED = SEED*SEED
               DO M =1,NDIM
                 SEED = SEED - SEED/2 
                IF( VERB) WRITE(*,*) ' SEED = ',SEED 
                 RNi =  ran(seed) 
                 F(M) = RNi
                 DO D=1,NDIM
                   RPNT(MN,D) = RPNT(MN,D) +  RNi*VAC(M,D) 
                 ENDDO
               ENDDO
               IF( VERB ) WRITE(6,*) "Random point ",RPNT(MN,:)
     &          ,' used for structure ',STNo,' molecule ',MN
!
!              Remove half the size of the added molecule from the top of the box
!                to keep the molecule out of the vacuum space 
!              Place molecule
!
               NWROT = .TRUE.
               ROTCNT = 0
               DO WHILE ( NWROT ) 
                 NWROT  = .FALSE.
                 ROTCNT = ROTCNT + 1 
                 ATEMPS = ATEMPS + 1

                 IF( RANROT(STNo) ) THEN
!                  Get random #'s 
                   SEED = SEED*SEED
                   SEED = SEED - SEED/2 
                   IF( VERB) WRITE(*,*) ' SEED = ',SEED 
                   ANa =  ran(seed)*( 2*PI )
                   SEED = SEED - SEED/2 
                   IF( VERB) WRITE(*,*) ' SEED = ',SEED 
                   ANb = ran(seed)*( 2*PI )
                    CALL ROTATE(STNo,ANa,ANb)
                 ENDIF 
 
!
!!!                RPNT(MN,SDIM) = RPNT(MN,SDIM)
!!!     &             -ANINT( RPNT(MN,SDIM)/VAC(SDIM,SDIM))*SZSTR(STN)
!                Check proximity to prvioulsy generated points
                 LDIM(:) = 0
c$$$                 DO M = 1,MN-1
c$$$                   DO D=1,NDIM
c$$$                     MD = RPNT(M,D) - RPNT(MN,D)
c$$$                     MDSQ = MD*MD
c$$$                     IF( MDSQ.LT.SZSQ(D) ) THEN
c$$$                       LDIM(D) = 1
c$$$                     ENDIF
c$$$                   ENDDO
c$$$                   IF( LDIM(1)*LDIM(2)*LDIM(3).GT.0 ) THEN
c$$$                     IF(VERB) THEN
c$$$                   WRITE(*,'( "Random pnt",I6," too close to ",I6)')MN,M
c$$$                   WRITE(*,*) LDIM(:)
c$$$                     ENDIF
c$$$                     OVERLAP = .TRUE.
c$$$                   ENDIF 
c$$$                 ENDDO                  
!                Update global  posisions 
                 CALL CENTSTR(STNo)    
                 CALL MAXMIN(STNo)
                 SZSTR(:) = ( MMN(2,:,STNo) - MMN(1,:,STNo) )
                 DO D = 1,NDIM
                   SZSQ(D) = SZSTR(D)*SZSTR(D)*SZF
                 ENDDO

!
!               Use another structure as box size 
!
                IF( UDBOX(STNo) ) THEN 
!                  DELTS = MAXH-MMN(1,SDIM,STNo)+DEPBF(STNo)
                  RPNT(MN,SDIM) = MAXH +0.50d0*SZSTR(SDIM)+DEPBF(STNo)
!
                  IF( VERB) THEN
                   WRITE(*,*)
     &            ' MAXH  MMN(1,SDIM,RSTN) DEPBF(STNo) MOL_HIEGHT'
                   WRITE(*,*) ' Max hieght of sub ',MAXH
                   WRITE(*,*) 'delta sdim',DELTS
                   WRITE(*,*) MMN(1,SDIM,STNo),DEPBF(STNo),RPNT(MN,SDIM)
                   WRITE(*,*) 'size sdim',SZSTR(SDIM)
                  ENDIF
                 ENDIF
!

                 DO IM =1,NA(STNo)
!                  Skip this check if rpnt check already failed 
!                   IF( OVERLAP ) EXIT
                   RT(IM,:) = R0(:,IM,STNo) + RPNT(MN,:) 
                   R0(:,IM,STNo) = RT(IM,:) 
                 ENDDO 

                 DO IM =1,NA(STNo)
!                  Skip this check if rpnt check already failed 
!                   IF( OVERLAP ) EXIT
                   ELim  = ELN(IM,STNo)
!                  Apply PBC
                   X = RT(IM,1)
                   Y = RT(IM,2)
                   Z = RT(IM,3)
               F(1)=-(-BC*CB*X+BB*CC*X+BC*CA*Y-BA*CC*Y-BB*CA*Z+BA*CB*Z)/
     &         (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC)
               F(2)=-(AC*CB*X-AB*CC*X-AC*CA*Y+AA*CC*Y+AB*CA*Z-AA*CB*Z)/
     &         (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC)
               F(3)=-(-AC*BB*X+AB*BC*X+AC*BA*Y-AA*BC*Y-AB*BA*Z+AA*BB*Z)/
     &         (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC)
!                  Reduce fractional coordinates to less than 1
                   DO D=1,2 !Warning only in xy plane
                     FT = F(D)
!                    Acount for postive nature of lattice 
                     IF( F(D).LT.0.0d0 ) FT = F(D) - 1.0d0
                     F(D) = F(D) -INT( FT )
                   ENDDO
!                  Get new coordinate based on fractional positions
                   RTL(1) = AA*F(1) + BA*F(2)+CA*F(3)
                   RTL(2) = AB*F(1) + BB*F(2)+CB*F(3)
                   RTL(3) = AC*F(1) + BC*F(2)+CC*F(3)
!                  Check proximity
                   DO So=1,STOTo
                    Si = So + STOTo
!                   Check if structure is included in random generation
                    IF( MKRAND(So) .OR. RANBOX(So) ) THEN

                     DO II = 1,NA(Si)
                       ELii  = ELN(II,Si)
!
!                      Apply PBC
!
                       X = R0(1,II,Si)
                       Y = R0(2,II,Si)
                       Z = R0(3,II,Si)
               F(1)=-(-BC*CB*X+BB*CC*X+BC*CA*Y-BA*CC*Y-BB*CA*Z+BA*CB*Z)/
     &         (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC)
               F(2)=-(AC*CB*X-AB*CC*X-AC*CA*Y+AA*CC*Y+AB*CA*Z-AA*CB*Z)/
     &         (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC)
               F(3)=-(-AC*BB*X+AB*BC*X+AC*BA*Y-AA*BC*Y-AB*BA*Z+AA*BB*Z)/
     &         (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC)
!
!                      Reduce fractional coordinates to less than 1
!
                       DO D=1,2    !Warning only in xy plane
                         FT = F(D)
                         IF( F(D).LT.0.0d0 ) FT = F(D) - 1.0d0
                         F(D) = F(D) -INT( FT )
                       ENDDO
!                      Get new coordinate based on fractional positions
                       RIL(1) = AA*F(1) + BA*F(2)+CA*F(3)
                       RIL(2) = AB*F(1) + BB*F(2)+CB*F(3)
                       RIL(3) = AC*F(1) + BC*F(2)+CC*F(3)
!                      Set buffer cylinder to covalent radi 
!                       CYBUF = ( CRDI(Elii) + CRDI(ELim) )**2
!                      Set buffer cylinder to van der Waals radi 
                       VDBUF = ( VDWRDI(Elii)+VDWRDI(ELim) )**2

                       DRSQ = 0.0d0
                       DO M = 1,NDIM
                         DELR(M) =  RTL(M) - RIL(M)
                         DELSQ(M) = DELR(M)*DELR(M)  
                         DRSQ = DRSQ +   DELSQ(M) 
                       ENDDO
!                      Track min speeration in SDIM for minization
                       IF( DELSQ(1).LT.CYBUF.AND.DELSQ(2).LT.CYBUF)
     &                    THEN
                           IF(DELR(3).LT.MINDR(SDIM) ) THEN
                            MINEL(1) = ELim
                            MINEL(2) = ELii
                            MINDR(SDIM) = DELR(SDIM)


                           ENDIF
                       ENDIF
                       IF ( DRSQ.LT.VDBUF ) THEN
                         NWROT =.TRUE.
                         EXIT
                       ELSE
                         IF( DRSQ.LT.MINOV ) MINOV = DRSQ 
                       ENDIF
                     ENDDO  !NAi
                    ENDIF   !MKRAND(Si)
                    IF( NWROT ) EXIT
                   ENDDO  !STN
                   IF( NWROT ) EXIT
                 ENDDO    ! NA(STNo) 
!
!                Minimize speration between two structures
!
                 IF( MINRD(STNo).AND.MCNT(STNi).LE.MINML(STNo)  ) THEN
                  ELim = MINEL(1)
                  ELii = MINEL(2)
                  VDBUF = ( VDWRDI(Elii)+VDWRDI(ELim) )*CRBUF                
                  IF( .NOT.OVERLAP.AND.MINDR(SDIM).GT.VDBUF ) THEN 
                   Si = MINRS(STNo) 
                   CALL MAXMIN(STNo)
                   CALL MAXMIN(Si)
                   SHIFT=VDBUF - MINDR(SDIM) 
!                  Only allow negative shifts
                   IF( SHIFT.LT.0.0d0) THEN
                     DO I=1,NA(STNo)
                       R0(SDIM,I,STNo)  = R0(SDIM,I,STNo) + SHIFT
                     ENDDO
!                    Check to make sure no overlap after shift
                   ENDIF
                   IF(VERB) THEN
                     WRITE(*,*) "Minimizing seperations between",Si
     &               ," and ",STNo
                     WRITE(*,*) " min seperation squared",MINOV
                     WRITE(*,*) " mole ",NA(STNo),' atoms'
                     WRITE(*,*) " min of mole ",MMN(1,:,STNo)
                     WRITE(*,*) " max of mole ",MMN(2,:,STNo)
                     WRITE(*,*) " min seperation in ",SDIM
     &                          ," ",MINDR(SDIM) 
                     WRITE(*,*) " shift is ",SHIFT
     &               ," to have a distance of",VDBUF
     &               ," based on the VDW radi",VDWRDI(Elii),VDWRDI(ELim)
                   ENDIF  
                  ELSE
                     IF(.NOT.OVERLAP.AND.VERB) THEN
                       WRITE(*,*) "No shift need"
     &                  ,MINDR(SDIM),"LT",VDBUF
                     ENDIF
                  ENDIF
                 ENDIF
!
                 IF ( ROTCNT.GT.ATMAX ) THEN
                   IF( VERB) WRITE(6,603)  MCNT(STNi),ROTCNT
                   OVERLAP = .TRUE.
                   EXIT
                 ENDIF
               ENDDO       ! NWROT
!!! Warning Error  ATEMPS has been hacked seems to work need to scale all max attempts 
!
               IF ( ATEMPS.GT.ATMAX*1000 ) THEN 
                 IF( VERB) WRITE(6,604)  MCNT(STNi),ROTCNT,NEWCNT 
                 MKNEW = .TRUE. 
                 ADDMOL = .FALSE. 
                 OVERLAP = .TRUE.
                 NWROT = .FALSE.
                 NEWCNT = NEWCNT + 1
!                Adjust size factor to determine good point if number of rotations too low
                 IF ( ROTCNT.LT.ATMAX/2 ) THEN
                   SZF = SZF - FLOAT(1/ATMAX)
                 ENDIF 
!                Calculate percent of molecule created
                 MFRAC = MCNT(STNi)/NGAS(STNo)
                 IF ( NEWCNT.GT.ATMAX.OR.MFRAC.LT..75 ) THEN
!                   If less than 3/4 are created or more than max attempts
                   IF( EXPRAN(STNo) ) THEN
                     Vk = Vj
                     VAC(1,:) = VAC(1,:) + VACEXP(1,:,STNo) 
                     VAC(2,:) = VAC(2,:) + VACEXP(2,:,STNo) 
                     VAC(3,:) = VAC(3,:) + VACEXP(3,:,STNo) 
                     AA = VAC(1,1)
                     AB = VAC(1,2)
                     AC = VAC(1,3)
                     BA = VAC(2,1)
                     BB = VAC(2,2)
                     BC = VAC(2,3)
                     CA = VAC(3,1)
                     CB = VAC(3,2)
                     CC = VAC(3,3)
                     Vj = AA*BB*CC
                     IF(VERB) THEN
                       WRITE(6,*)'Expanding vacuum space from'
                       WRITE(6,'(F16.2," to ",F16.2)') Vk,Vj
                       WRITE(6,*) VAC(1,1:3)
                       WRITE(6,*) VAC(2,1:3)
                       WRITE(6,*) VAC(3,1:3)
                       WRITE(6,*) " while adding in molecule:",MN 
 
                     ENDIF
                   ELSE
                     WRITE(6,*) 'Try inceasing box size'
                     WRITE(6,*) 'decreasing the number of molecules'
                     WRITE(6,*) 'or use the the "expandrand" function' 
                     WRITE(6,*) 'created structure in ',OUTXYZ
                     CALL WRITE_XYZ(STNi) 
                     STOP
                   ENDIF
                 ENDIF
               ENDIF
               IF ( MKNEW ) EXIT
             ENDDO  ! OVERLAP
             IF ( MKNEW ) EXIT
!
c$$$!            Add potential paramters for new molecule
c$$$             DO D=1,NC(STNo) 
c$$$                I =  NC(STNi) + D
c$$$                CONIi(STN,I) = CONIo(STN,D)+ NREF(STNo)
c$$$                CONJi(STN,I) = CONJo(STN,D)+ NREF(STNo)
c$$$                CTYPi(STN,I) = CTYPo(STN,D)
c$$$                CDISIi(STN,I) = CDISIo(STN,D)
c$$$                CDISJi(STN,I) = CDISJo(STN,D)
c$$$             ENDDO
c$$$             DO D=1,NB(STNo) 
c$$$                I =  NB(STNi) + D
c$$$                BNDIi(STN,I) = BNDIo(STN,D) + NREF(STN)
c$$$                BNDJi(STN,I) = BNDJo(STN,D) + NREF(STN)
c$$$                BTYPi(STN,I) = BTYPo(STN,D) 
c$$$             ENDDO
c$$$             DO D=1,NPR(STNo)
c$$$                 I = NPRi(STN ) + D
c$$$                 PRSIi(STN,I)  = PRSIo(STN,D)  + NREF(STN)
c$$$                 PRSJi(STN,I)  = PRSJo(STN,D)  + NREF(STN)
c$$$                 PRTYPi(STN,I) = PRTYPo(STN,D)
c$$$             ENDDO
c$$$             DO D=1,NANG(STNo) 
c$$$                I = NANG(STNi) + D
c$$$                ANGIi(STN,I)  = ANGIo(STN,D) + NREF(STN)
c$$$                ANGJi(STN,I)  = ANGJo(STN,D) + NREF(STN)
c$$$                ANGKi(STN,I)  = ANGKo(STN,D) + NREF(STN)
c$$$                ANTYPi(STN,I) = ANTYPo(STN,D)
c$$$             ENDDO
c$$$             DO D=1,NDIH(STNo)
c$$$                I = NDIH(STNi) + D
c$$$                DIHIi(STN,I)  = DIHIo(STN,D) + NREF(STN)
c$$$                DIHJi(STN,I)  = DIHJo(STN,D) + NREF(STN)
c$$$                DIHKi(STN,I)  = DIHKo(STN,D) + NREF(STN)
c$$$                DIHLi(STN,I)  = DIHLo(STN,D) + NREF(STN)
c$$$                DHTYPi(STN,I) = DHTYPo(STN,D)
c$$$             ENDDO
c$$$ 
!            Add to global  atoms        


             DO I = 1,NA(STNo)
                N = NA(STNi) + I
                CALL ATPASS(I,STNo,N,STNi)       
             ENDDO
!
!            Update total number of atoms NP
!
             NA(STNi)   = NA(STNi)   + NA(STNo)
             NB(STNi)   = NB(STNi)   + NB(STNo)
!             NPR(STNi)  = NPR(STNi)  + NPR(STNo)
             NANG(STNi) = NANG(STNi) + NANG(STNo)
             NDIH(STNi) = NDIH(STNi) + NDIH(STNo)
!
!            Save last numbers of some counts as a reference
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
                 WRITE(6,'(A,I6,A,F10.2)')
     &           '  Atempt :',ATEMPS,' min ',SQRT(MINOV)
                WRITE(6,'(A,2I6)') ' Ref :',RNREF,MREF
                WRITE(6,*) 
             ENDIF
           ENDIF    ! NCNT(STN).LT.NGAS(STN) 
           IF ( MKNEW ) EXIT
         ENDIF      ! MKRAND
         IF ( MKNEW ) EXIT
 ! debug
         IF( VERB ) WRITE(*,*) "Starting to check next str ",STNo+1
        ENDDO       ! STN
        IF ( MKNEW ) EXIT
! debug
         IF( VERB ) WRITE(*,*) "Starting new molecule ",MN+1
       ENDDO        ! ADDMOL 
! debug
       IF( VERB ) WRITE(*,*) "Starting new random structure ",NEWCNT
      ENDDO         ! MKNEQ
!     Update global lattice vector
      LV(:,:,STNi) = VAC(:,:)
!
      IF( VERB ) THEN
        WRITE(6,*) ' Final structure has: ',NA(STNi)
     &           ,' atoms from',NA(STNo)
        WRITE(6,*) 'Final vacuum space:'
        WRITE(6,'(3F8.2)') VAC(1,:)
        WRITE(6,'(3F8.2)') VAC(2,:)
        WRITE(6,'(3F8.2)') VAC(3,:)
        CALL MAXMIN(STNi)
        WRITE(*,605) " min of added mole ",MMN(1,:,STNi)
        WRITE(*,605) " max of added mole ",MMN(2,:,STNi)
      ENDIF
!
      !
      DEALLOCATE(RT,NREF,RPNT,MCNT )
!
      RETURN
 601  FORMAT("Molecule ",I6," not added, attempts ",I6," gt ",I6)
 602  FORMAT("Molecule ",I6," not added, new guesses ",I6," gt ",I6)
 603  FORMAT(" Rotated  ",2I6," times trying new point")
 604  FORMAT("Molecule ",I6," not added,rotations:"
     & ,I6," system restarts",I6)
 605  FORMAT(A,3F12.3)
      END SUBROUTINE mkrandom
