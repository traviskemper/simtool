!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  NanoMD MDbuild makenano ...                                 !
!     Minipulate molecules and polymers                        !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Main file
!
      INCLUDE 'module_const.f'
      INCLUDE 'module_elements.f'
      INCLUDE 'module_structure.f'
      INCLUDE 'module_specify.f'
      INCLUDE 'module_potential.f'
      INCLUDE 'module_build.f'

      PROGRAM mdtool
      USE specify
      USE structure 
      USE build 
CC      USE potential
!
      IMPLICIT none
      INTEGER :: STNo,STNi,MN
      REAL*8 ::  ANa,ANb 
      LOGICAL :: CLDIPf
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Read in settings and allocate arrays 
      STNo = 1
      CALL READIN
C      IF( RANSD ) CALL SETSEED
      CALL ALLOCATE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Read in structural information 
      DO STNo = 1, STOTo- STADD 
!       Read in structure
        IF( RXYZ(STNo) )    CALL READ_XYZ(STNo)
        IF( RCRD(STNo) )    CALL READ_CRD(STNo)
        IF( RGRO(STNo) )    CALL READ_GRO(STNo)
        IF( RPDB(STNo) )    CALL READ_PDB(STNo)
        IF( RCAR(STNo) )    CALL READ_CAR(STNo)
        IF( RTOP(STNo) )    CALL READ_TOP(STNo)
        IF( RTNK(STNo) )    CALL READ_TNK(STNo)
        IF( RTPA(STNo) )    CALL READ_TPA(STNo)
        IF( RXMOL(STNo) )    CALL READ_XMOL(STNo)
        IF( RTBML(STNo) )    CALL READ_TBML(STNo)
        IF( MTOAN(STNo) )    CALL FINDELN(STNo)
        IF( RCHR(STNo)  )    CALL READ_CHR(STNo)
!
C        IF( RSURF(STNo))    CALL READ_SURF(STNo)
C        IF( RRTP(STNo) )    CALL READ_RTP(STNo)
C        IF( MKCHAIN(STNo) ) CALL BUILDCHAIN(STNo)
!        IF( ACON(STNo) )  CALL ADDCONST(STNo) 
!       Modify initial inputs
        IF( CPREG(STNo) ) CALL CPFAC(STNo) 
        IF( CFRAC(STNo) ) CALL CUTFRAC(STNo) 
        IF( CMOL(STNo) )  CALL CUTMOL(STNo)

        IF( RMATTP(STNo) )  CALL RMATOMTYP(STNo)
        IF( MKCENT(STNo) )  CALL CENTSTR(STNo)
        IF( SZMAX(STNo) ) CALL ZSHIFT(STNo)
        IF( SHFT(STNo) ) CALL SHIFT(STNo)

!       Calculate properties
        IF( CLLC(STNo) )  CALL CALCLC(STNo)
        IF( CLLV(STNo) )  CALL CALCLV(STNo)
        IF( CLMAS(STNo) )  CALL CALCMAS(STNo) 
        IF( CLVD(STNo) )  CALL CALCVD(STNo)

        CALL CALCNELM(STNo)
        CALL MAXMIN(STNo)
        IF( NBLST       )  CALL BUILDNBL(STNo)
        IF( MOLCH(STNo)  ) CALL MOLCHECK(STNo)
        IF ( GBONDS(STNo) ) CALL GET_BONDS(STNo)
        CALL CENTMASS(STNo) 
        CALL MASSPEC(STNo)

        CALL BUILDMNB(STNo)
        IF( CLDIP(STNo))  CALL TOTDIP(STNo)
        IF( CHR_UPDATE(STNo) ) CALL READ_MCHR(STNo)
!      
!        Asign tinker ID numbers
!
            IF( TNKOPLS(STNo)  ) CALL ATYPNOPLSAA(STNo)
            IF( TNKAMOEBA(STNo) ) CALL ATYPNAMOEBA(STNo) 
            IF( TNKADD(STNo) )   CALL ATYPNADD(STNo)  

        IF( CUTMOLNB(STNo) )  CALL molnbcut(STNo)

! Should operate on final structures !!!!!!!!!!
        IF( ANGAXS(STNo) ) CALL ANGAXM(STNo) 
        IF( TANGAXS(STNo) ) CALL ANGAXM_TRI(STNo) 
        IF( PRDIS(STNo)  ) CALL PRDIST(STNo)
        IF ( MKMRDF(STNo) ) CALL MOLRDF(STNo)
        IF( CMDIP(STNo))  CALL CMOLDIP(STNo)
        IF( IPROP(STNo) ) CALL CINDPROP(STNo)
        IF( MEFL(STNo) ) CALL MOLEFLD(STNo)
        IF( ATEFL(STNo) ) CALL ATOMEFLD(STNo)
        IF( GMMIN(STNo) ) CALL MOL_MIN(STNo)
!
!       Multi frame analysis
!
!       Initialize analysis
        IF( CALC_ATOM_ANG(STNo) ) CALL MF_INTL_ANGLE(STNo)
!
!       Read through frames
!
        IF( RXMOL(STNo) ) CALL MF_READ_XMOL(STNo) 
!
!       Print analysis results 
!
        IF( CALC_ATOM_ANG(STNo) ) CALL MF_FIN_ANGLE(STNo)
!
!       Center based on molecule center of mass
        MN= MICENT(STNo)    
        IF( MLCENT(STNo) )  CALL CENTMOL(MN,STNo)        
c$$$!        IF( MOLCH(STNo)  )  CALL MOLCHECK(STNo,0) !need to fic nblisto
c$$$!       Find inion cut out 
c$$$        IF( CUTONION(STNo) ) CALL FINDCUT(NA,STNo)
!       Add velocity 
        IF( MKVEL(STNo) )    CALL ADDVEL(STNo)
        IF( REFIX(STNo) )    CALL FIXREGION(STNo)
        IF( MKDTMP(STNo) )   CALL ADDTEMP(STNo)
c$$$!       Calculate the number of final atoms
c$$$        CALL CALCMXN(STNo,NA)
      ENDDO 
!     Write properties of read in structures

       CALL WRITE_PROPo
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Opperate on orignal structures to produce final structures
!
      IF( WRTOUT ) THEN 
          CALL WRITE_PROC 
          DO STNo = 1,STOTo
            STNi = STOTo + STNo 
!           Reset atomic properties
            IF( STTOG ) CALL TTOG(STNo)
!           Find inion cut out
C           IF ( CUTONION(STNo) ) CALL CUTREGION(STNo)
C            IF(.NOT.MKSUP(STNo) ) 
            STNi = STNo + STOTo
            CALL STRPASS(STNo,STNi)
!           Make supercell
           IF( MKSUP(STNo) ) CALL SUPERCELL(STNo)
           IF( MKSLAB(STNo) ) CALL SLAB(STNi)
           IF( MKVAC(1) ) THEN
                LV(SDIM,SDIM,STNi) = VACL
                LC(SDIM,STNi) = VACL 
            ENDIF 
          ENDDO
!         Build randomized structures simultaneously
          IF( CRAND )  CALL  MKRANDOM
!         Add deposited molecule
          IF( ADDDM ) CALL ADDMOL 
          CLDIPf = .FALSE. 
          DO STNi=STOTo+1,STOT-1 
            STNo = STNi - STOTo 
!           Hydrogen terminate structures
            IF( MKHT(STNo)  ) CALL HTERM(STNi)
!           Rename atom id's and residue id's
            IF( REN(STNo)  ) CALL RENAME_TYP(STNi)
            IF( NBCHRG(STNo)) CALL CHRGNB(STNi)
C                  IF( CHFT(STNi)  ) CALL FXTRN(STN,1)
C                  IF( CFRAC(STNi) ) CALL CUTFRAC(STNi)
            IF(ROTA(STNo) ) THEN 
               ANa = RANA(STNo)*(PI/180.0d0)
               ANb = RANB(STNo)*(PI/180.0d0)
               CALL ROTATE(STNi,ANa,ANb)
            ENDIF
!
!           Reset tinker ID numbers
!
            IF( MTNKADD(STNo) ) CALL  MOL_ATYPNADD(STNi)
            IF ( TNKNMD(STNo) ) CALL ATYPTNKMOD(STNi) 
            IF ( MTNKNMD(STNo) ) CALL MTYPTNKMOD(STNi) 
!
!           Calculate final properties
            IF( UPBC(STNo) ) CALL molPBC(STNi)          !Apply PBC's 
            ! IF( UPBC(STNo) ) CALL PBC(STNi)          !Apply PBC's 
C                  IF( NBLST ) CALL FINDNBi(STNi)
C                  IF( MOLCH(STNi)  ) CALL MOLCHECK(STN,1)
C                  IF( ANGAXS(STNi) ) CALL ANGAXM(STNi)
C                  IF( MXNDEN(STNi) ) CALL MXMNDENS(STNi)
C                  CALL FINDIMP(STNi)
C                  CALL CMDIST(STNi) 
C                  CALL CMPRDIST(STNi)
            IF( CLVD(STNo) )  CALL CALCVD(STNi)
            CALL MAXMIN(STNi)
            CALL CALCMAS(STNi)         
            CALL CALCVD(STNi)
            IF( CLDIP(STNo) ) THEN
             CALL TOTDIP(STNi)
             CLDIPf = .TRUE. 
            ENDIF 
            IF( NBLST       )  CALL BUILDNBL(STNi)
            CALL BUILDMNB(STNi)
            IF( WAMONIOM(STNo) ) CALL WRTMOLCOM(STNi) 
            IF( WAMONIOMM(STNo) ) CALL MOLCOM(STNi) 
            IF( WONIOMSP(STNo) ) THEN
              CALL SPCOM(STNi)  
              CALL WRITE_SPCOM(STNi) 
            ENDIF 
          ENDDO
!         Write properties of read in structures
          CALL WRITE_PROPi
!    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!           Set final structure 
!    
!           Combine  all intermdiate structures
            IF( MKSTACK ) CALL STACK
            CALL STRREND
            STNi = STOT
            CALL CALCMAS(STNi)         
            CALL CALCVD(STNi)
            CALL CALCNELM(STNi)
            CALL MAXMIN(STNi)
            IF(  NBLST       )  CALL BUILDNBL(STNi)
            IF(  MOLCH(STNi)  ) CALL MOLCHECK(STNi)
            IF ( GBONDS(STNi) ) CALL GET_BONDS(STNi)

            CALL CENTMASS(STNi) 
            CALL MASSPEC(STNi)
            CALL BUILDMNB(STNi)
 
            IF( CLDIPf ) CALL TOTDIP(STNi)
            WRITE(6,602) 
            CALL WRITE_STRP(STNi)

!           Write final structure 
            IF( WEMB ) THEN
               CALL DEFINE_EMBED(STOT) 
            ELSE
               IF( WCOM  ) CALL WRITE_COM(STOT)
               IF( WTBML ) CALL WRITE_TBML(STOT)
               IF( WXYZ  ) CALL WRITE_XYZ(STOT) 
            ENDIF 
            IF( WCRD ) CALL WRITE_CRD(STOT) 
            IF( WGRO ) CALL WRITE_GRO(STOT)
            IF( WPDB ) CALL WRITE_PDB(STOT)
            IF( WCAR ) CALL WRITE_CAR(STOT)
            IF( WTOP ) CALL WRITE_TOP(STOT)
            IF( WTNK ) CALL WRITE_TNK(STOT)
            IF( WLMP ) CALL WRITE_LMP(STOT)
            IF( WMCOMS ) CALL WRITE_MCOMS(STOT)
!
      ENDIF ! WRTOUT 
!
      CALL DEALL
 602  FORMAT('--  Combined structure  --')
      END PROGRAM mdtool
!
      INCLUDE 'subroutines.f'
