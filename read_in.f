!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Read data from standard input
!
      SUBROUTINE readin
      USE structure 
      USE specify
      USE elements  
      USE potential 
      USE build 
!
      IMPLICIT none
!
      INTEGER:: istat,rstat,S,I,STN,J,N,FMAX,Io,Jo,To
     & ,CNMX,STMX,CNT(12),PNUMB
      REAL*8 :: V(NDIM),R,VEC(NDIM),AA,AB,MS,CR,VR
      LOGICAL :: LIN,LG1(NDIM)
      CHARACTER(CHSZ) :: kywrd,var,VARi,ID,Ao,Bo,Co,Cc,Dc
     &    ,FMT1,FMT2,FMT3,FMT4,FMT5,ESTAT,TP
!
      WRITE(6,*) '!++++++++++++++++++++++++++++++++++++++++++++++!'
      WRITE(6,*) '!    Programer:                                !'
      WRITE(6,*) '!      Travis Kemper                           !'
      WRITE(6,*) '!      tkemper3@mail.gatech.edu                !'
      WRITE(6,*) '!++++++++++++++++++++++++++++++++++++++++++++++!'
      WRITE(6,*) '!    Version:                                  !'
      WRITE(6,*) '!     5.334  12.21.2012                        !'
      WRITE(6,*) '!++++++++++++++++++++++++++++++++++++++++++++++!'
!
!     Determine the number of structure in input file
      STOTo = 0
      STADD = 0 
      CNT(:) = 0 
!     Run modifications and output geometry  
      WRTOUT = .FALSE. 
      READ(5,*,iostat=istat) kywrd
      BACKSPACE(5)

      DO WHILE (ISTAT.EQ.0)
        IF( kywrd .EQ. 'verbose' ) THEN
          READ(5,*) 
          VERB =.TRUE.
         ELSEIF( kywrd .EQ. 'help' ) THEN
          READ(5,*) 
          HELP =  .TRUE.
        ELSEIF( kywrd .EQ. 'debug' ) THEN
          READ(5,*) 
          DEBUG =  .TRUE.
        ELSEIF( kywrd .EQ. 'inxyz' ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(36,ESTAT)
          IF( S.GT.STOTo ) STOTo = S
        ELSEIF( kywrd .EQ. 'inxmol' ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(61,ESTAT)
          IF( S.GT.STOTo ) STOTo = S
        ELSEIF( kywrd .EQ. 'inpdb'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(37,ESTAT)
          IF( S.GT.STOTo ) STOTo = S
        ELSEIF( kywrd .EQ. 'ingro'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(39,ESTAT)
          IF( S.GT.STOTo ) STOTo = S
        ELSEIF( kywrd .EQ. 'incrd'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(38,ESTAT)
          IF( S.GT.STOTo ) STOTo = S
        ELSEIF( kywrd .EQ. 'incar'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(40,ESTAT)
          IF( S.GT.STOTo ) STOTo = S
        ELSEIF( kywrd .EQ. 'intnk'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(42,ESTAT)
          IF( S.GT.STOTo ) STOTo = S
        ELSEIF( kywrd .EQ. 'inlammps'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(57,ESTAT)
          IF( S.GT.STOTo ) STOTo = S
        ELSEIF( kywrd .EQ. 'inturbomol'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(58,ESTAT)
          IF( S.GT.STOTo ) STOTo = S
        ELSEIF( kywrd .EQ. 'inchr'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(1,ESTAT)
          IF( S.GT.STOTo ) STOTo = S
         ELSEIF( kywrd .EQ. 'outxyz' ) THEN
           READ(5,*) 
           WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'outpdb' ) THEN
           READ(5,*) 
           WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'outcrd' ) THEN
            READ(5,*) 
            WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'outgro' ) THEN
           READ(5,*) 
            WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'outtop' ) THEN
           READ(5,*) 
           WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'outcar' ) THEN
            READ(5,*) 
           WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'outcom' ) THEN
           READ(5,*) 
           WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'outtnk' ) THEN
           READ(5,*) 
           WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'outmcoms' ) THEN
           READ(5,*) 
           WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'outlammps' ) THEN
           READ(5,*) 
           WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'outturbomol' ) THEN
           READ(5,*) 
           WRTOUT = .TRUE. 
        ELSEIF( kywrd .EQ. 'copyregion'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(91,ESTAT)
          STADD = STADD + 1
        ELSEIF( kywrd .EQ. 'rmatomtype'  ) THEN
          READ(5,*,iostat=rstat) var,S,Ao
           IF( rstat.NE.0) CALL PRNERROR(25,ESTAT)
          CNT(1) = CNT(1) + 1
        ELSEIF( kywrd .EQ. 'vecangle'  ) THEN
          READ(5,*) 
          CNT(2) = CNT(2) + 1
       ELSEIF( kywrd .EQ. 'vectriangle'  ) THEN
          READ(5,*) 
          CNT(11) = CNT(11) + 1
       ELSEIF( kywrd .EQ. 'atomrdftypes' ) THEN
          READ(5,*) 
          CNT(3) = CNT(3) + 1
        ELSEIF( kywrd .EQ. 'molrdfanumbs' ) THEN
          READ(5,*) 
          CNT(4) = CNT(4) + 1
        ELSEIF( kywrd .EQ. 'addsurfatom' ) THEN
          READ(5,*) 
          CNT(5) = CNT(5) + 1
        ELSEIF( kywrd(1:11) .EQ. 'rename' ) THEN
          READ(5,*) 
          CNT(6) = CNT(6) + 1
        ELSEIF( kywrd .EQ. 'stack' ) THEN
           READ(5,*) 
          CNT(7) = CNT(7) + 1
        ELSEIF( kywrd .EQ. 'indexprop' ) THEN
          READ(5,*) 
          CNT(9) = CNT(9) + 1
        ELSEIF( kywrd .EQ. 'moleculecom' ) THEN
          READ(5,*) 
          CNT(10) = CNT(10) + 1
        ELSEIF( kywrd .EQ. 'atom_angle_dist' ) THEN
          READ(5,*) 
          CNT(12) = CNT(12) + 1
         ELSE
          READ(5,*,iostat=istat)
        ENDIF
        READ(5,*,iostat=istat) kywrd
        BACKSPACE(5)
      ENDDO
      IF(VERB) THEN
        WRITE(6,*) "Number of  structures found in input file",STOTo
        WRITE(6,*) "Counts ",CNT(:) 
        IF( DEBUG) WRITE(6,*) "Debug mode is on, ",
     &        "results are in debug.out"
      ENDIF
      REWIND(5)
!     Add one for the final structure
      STOTo = STOTo + STADD 
      STOT = STOTo*2 + 1 
!     Allocate # of structure dependent arrays 
      ALLOCATE( NA(STOT),LV(NDIM,NDIM,STOT),LC(NDIM,STOT),LA(NDIM,STOT)
     & ,VOL(STOT),DEN(STOT),TMAS(STOT),MMN(2,NDIM,STOT)
     & ,TITLE(STOT),INXYZ(STOT),INPDB(STOT),INTRF(STOT),INTOP(STOT)
     & ,INGRO(STOT),INCRD(STOT),INCAR(STOT),INCOM(STOT)
     & ,INTBML(STOT),RTBML(STOT)
     & ,RXYZ(STOT),RPDB(STOT),RTOP(STOT),RGRO(STOT),RCAR(STOT)
     & ,RCRD(STOT),RTRF(STOT),RCOM(STOT)
     & ,INXMOL(STOT),RXMOL(STOT) ! xmol in
     & ,ELT(STOT),ELCNT(NELM,STOT),NEL(NELM,STOT)
     & ,NB(STOT),NANG(STOT),NDIH(STOT),NIMP(STOT),NMOL(STOT)
     & ,NBTYP(STOT),NATYP(STOT),NDTYP(STOT),NITYP(STOT)
     & ,MKRAND(STOT),FIXLV(STOT),RANBOX(STOT)
     & ,RANROT(STOT),EXPRAN(STOT),MINRD(STOT),MINRS(STOT)
     & ,MINML(STOT),MINBF(STOT),MTOAN(STOT)
     & ,NGAS(STOT),VACEXP(NDIM,NDIM,STOT)
     & ,CLLC(STOT),CLLV(STOT),CLVD(STOT),CLMAS(STOT)
     & ,MKCENT(STOT),MKVEL(STOT),VEL(NDIM,STOT)
     & ,MKSUP(STOT),SC(NDIM,STOT),LCRD(STOT),LVRD(STOT)
     & ,RANFRX(NDIM,STOT),RANFRSZ(NDIM,STOT),MOLCH(STOT)
     & ,NATPR(STOT), RATOMTP(CNT(1),STOT),RMATTP(STOT)                 ! rmatomtype
     & ,CFRAC(STOT),CUTFC(2,NDIM,STOT)                                 ! cutregion
     & ,MICENT(STOT),MCCENT(NDIM,STOT),MLCENT(STOT),MLMNCENT(STOT)     ! molcenter
     & ,REFIX(STOT),RFIXC(2,NDIM,STOT),FRTAG(NDIM,STOT)                ! fixregion 
     & ,ANAX(2,CNT(2),STOT),ANGX(NDIM,CNT(2),STOT),ANGAXS(STOT)        ! angle distribution 
     & ,ANGAXFL(CNT(2),STOT),ANCNT(STOT),DPANGAXFL(STOT)               ! angle distribution 
     & ,CPREG(STOT), CPRDS(NDIM,STOT),RCPRG(2,NDIM,STOT),STCPN(STOT)   ! copyregion 
     & ,PRDIS(STOT),UDIS(NDIM,CNT(3),STOT),RDIS(2,CNT(3),STOT)         ! atomic rdf 
     & ,PRCNT(STOT),PRDISFL(CNT(3),STOT),RMIAM(STOT)                   ! atomic rdf 
     & ,MKMRDF(STOT),NMOLRi(CNT(4),STOT),NMOLRj(CNT(4),STOT)           ! molecular rdf 
     & ,NMCNT(STOT),MPRDISFL(CNT(4),STOT),UMDIS(NDIM,CNT(4),STOT)      ! molecular rdf 
     & ,UPBC(STOT),PBCS(NDIM,STOT)                                     ! PBC's 
     & ,MKSLAB(STOT),SCUT(STOT)                                        ! Surface
     & ,SATOMN(STOT),SATOMID(CNT(5),STOT),SATOMV(NDIM,CNT(5),STOT)     ! Add surface atoms 
     & ,REN(STOT),RENM(STOT),RNIDI(CNT(6),STOT),RNIDF(CNT(6),STOT)               ! Relabel 
     & ,NBCHRG(STOT)                                                   ! Use naborlist for charge groups
     & ,STKI(CNT(7)),STKJ(CNT(7)),STKBUF(CNT(7)),SHFIX(CNT(7))                 ! Stack geometeries  
     & ,MKHT(STOT)                                                     ! Hterm 
     & ,MKVAC(1),SVOL(STOT),SDEN(STOT),SHT(STOT)                      ! surface
     & ,ROTA(STOT),RANA(STOT),RANB(STOT)                    ! rotate all 
     & ,TCHG(STOT),DIP(STOT),CLDIP(STOT),CMDIP(STOT)                ! total molecular dipole
     & ,IPROP(STOT),IDPN(STOT),IDPID(CNT(9),STOT)         !  index properties 
     & ,SZMAX(STOT),SZDIR(STOT),SHFT(STOT),SHVEC(NDIM,STOT)    ! shift 
     & ,MDDIM(STOT),MOLDIPFL(STOT),CMLDPD(STOT),CMDIPFL(STOT) ! Molecular dipole distribution 
     & ,MOLDIPAXFL(STOT)
     & ,RSCALE(STOT),DSCALE(STOT),DMAX(STOT)     ! Molecular dipole distribution 
     & ,WAMONIOM(STOT) ,WAMONIOMA(STOT),WAMONIOMM(STOT),WONIOMSP(STOT) !     Print oniom inputs   
     & ,ONRAD(2,STOT),ONCENT(NDIM,STOT),USONIOMFX(STOT)  !  QM/MM lists 
     & ,ONIONQM(STOT),ONIONMM(STOT),ONIONFX(STOT ) ! QM/MM lists 
     & ,NONMOL(STOT),WONIONMOL(CNT(10),STOT)       ! QM/MM molecules 
     & ,MEFL(STOT),MEFLFL(STOT),ATEFL(STOT)        ! Elctric field at molecule     
     & ,RTNK(STOT),INTNK(STOT),RTPA(STOT),INTPA(STOT)   ! atom type number for tinker
     & ,RLMP(STOT),INLMP(STOT)
     & ,QMNB(STOT)                                ! QM molecule nieghbors to be added to QM region
     & ,TANAX(3,CNT(11),STOT),TANGX(NDIM,CNT(11),STOT),TANGAXS(STOT)        ! 3 vec angle di/updastribution 
     & ,TANGAXFL(CNT(11),STOT),TANCNT(STOT),DPTANGAXFL(STOT)               ! 3 vec angle distribution 
     & ,GMMIN(STOT)   ! minimize each molecule with gromacs 
     & ,CHR_UPDATE(STOT) !     Update charges 
     & ,MFRAM_XYZ(STOT)   !     Multi fram analysis
     & ,ATOM_NB(3,STOT) !     Atomic analysis arrays
     & ,CALC_ATOM_ANG(STOT),ATOM_ANAX(3,CNT(12),STOT)
     & ,ATOM_ANGL_CNT(STOT),ATOM_ANFILE(CNT(12),STOT)    !     Atomic angle distribution 
     & ,RDBCNST(STOT),RDACNST(STOT),RDDCNST(STOT),RDICNST(STOT)      ! Read in constants from top file 
     & ,UDBOX(STOT),NDEP(STOT),RSTBOX(STOT),DEPBF(STOT)  ! Depostion 
     & ,MKDTMP(STOT),DEPTEMP(STOT),ADDDMOL(STOT)               ! Depostion temperature  
     & ,TNKAMOEBA(STOT),TNKOPLS(STOT),TNKADD(STOT),TNKADDI(STOT)   ! Tinker types
     & ,TNKNMD(STOT),TNKNUMBM(3,STOT),MTNKNMD(STOT),MTNKNUMBM(3,STOT)                    ! modify tinker numbers 
     & ,MTNKADD(STOT),MTNKADDI(STOT)   ! Tinker types
     & ,CMOL(STOT),CMOLN(2,STOT)   ! Cut molecule 
     & ,INCHR(STOT),RCHR(STOT)   ! charges 
     & ,GBONDS(STOT),GANGLES(STOT),GDIHS(STOT)   ! force field info 
     & ,CUTMOLNB(STOT)   ! cut out 
     & )
!     Set for default names 
      NATPR(:) = CNT(1)
      ANCNT(:) = CNT(2)
      PRCNT(:) = CNT(3) 
      NMCNT(:) = CNT(4)
      TANCNT(:) = CNT(11)
      ATOM_ANGL_CNT(:) =CNT(12)
!
!     Set defaults 
      CALL DEFAULT
!
!     Counts for key words which can repeat 
      ANCNT(:) = 0 
      NATPR(:) = 0 
      PRCNT(:) = 0 
      NMCNT(:) = 0 
      RENM(:) = 0 
      SATOMN(:) = 0 
      TANCNT(:) = 0 
      ATOM_ANGL_CNT(:) = 0
!
      READ(5,*,iostat=istat) kywrd
      BACKSPACE(5)
      DO WHILE (ISTAT.EQ.0)
        IF( kywrd .EQ. 'verbose' ) THEN
          READ(5,*) 
          VERB =.TRUE.
        ELSEIF( kywrd .EQ. 'help' ) THEN
          READ(5,*) 
          HELP =  .TRUE.
        ELSEIF( kywrd .EQ. 'debug' ) THEN
          READ(5,*) 
          DEBUG =  .TRUE.
        ELSEIF( kywrd .EQ. 'inxyz' ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(36,ESTAT)
          INXYZ(S) =  vari
          RXYZ(S)  = .TRUE.   
        ELSEIF( kywrd .EQ. 'inpdb'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
          IF( HELP ) WRITE(*,*) ' '
           IF( rstat.NE.0) CALL PRNERROR(37,ESTAT)
          INPDB(S)  =  vari
          RPDB(S)   = .TRUE.
        ELSEIF( kywrd .EQ. 'incrd'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(38,ESTAT)
          INCRD(S)  =  vari
          RCRD(S)   = .TRUE.
        ELSEIF( kywrd .EQ. 'ingro' ) THEN 
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(39,ESTAT)
          INGRO(S) =  vari
          RGRO(S)  = .TRUE.    
        ELSEIF( kywrd .EQ. 'incar'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(40,ESTAT)
          INCAR(S) = vari
          RCAR(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'intop' ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(41,ESTAT)
          INTOP(S) = vari  
          RTOP(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'intnk' ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(42,ESTAT)
          INTNK(S) = vari  
          RTNK(S) = .TRUE.
        ELSEIF( kywrd .EQ. 'inchr' ) THEN
          READ(5,*,iostat=rstat) var,S,vari
          IF( rstat.NE.0) CALL PRNERROR(1,ESTAT)
          INCHR(S) = vari
          RCHR(S) = .TRUE.
        ELSEIF( kywrd .EQ. 'tnkoplsaa' ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(43,ESTAT)
          TNKOPLS(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'tnkamoeba' ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(44,ESTAT)
          TNKAMOEBA(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'tnkadd' ) THEN
          READ(5,*,iostat=rstat) var,S,AA
           IF( rstat.NE.0) CALL PRNERROR(45,ESTAT)
          TNKADDI(S) = AA
          TNKADD(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'tnkmoladd' ) THEN
          READ(5,*,iostat=rstat) var,S,AA
           IF( rstat.NE.0) CALL PRNERROR(46,ESTAT)
          MTNKADDI(S) = AA
          MTNKADD(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'intpa' ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(47,ESTAT)
          INTPA(S) = vari  
          RTPA(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'outxyz' ) THEN
          READ(5,*,iostat=rstat) var,vari
           IF( rstat.NE.0) CALL PRNERROR(48,ESTAT)
          OUTXYZ = vari
          WXYZ = .TRUE.
        ELSEIF( kywrd .EQ. 'outtnk' ) THEN
          READ(5,*,iostat=rstat) var,vari
           IF( rstat.NE.0) CALL PRNERROR(49,ESTAT)
          OUTTNK = vari
          WTNK = .TRUE.
        ELSEIF( kywrd .EQ. 'outpdb' ) THEN
          READ(5,*,iostat=rstat) var,OUTPDB
           IF( rstat.NE.0) CALL PRNERROR(50,ESTAT)
          WPDB = .TRUE.
        ELSEIF( kywrd .EQ. 'outcrd' ) THEN
          READ(5,*,iostat=rstat) var,OUTCRD
           IF( rstat.NE.0) CALL PRNERROR(51,ESTAT)
          WCRD = .TRUE.
        ELSEIF( kywrd .EQ. 'outgro' ) THEN
          READ(5,*,iostat=rstat) var,OUTGRO
           IF( rstat.NE.0) CALL PRNERROR(52,ESTAT)
           WGRO = .TRUE.
         ELSEIF( kywrd .EQ. 'outtop' ) THEN
          READ(5,*,iostat=rstat) var,OUTTOP
           IF( rstat.NE.0) CALL PRNERROR(53,ESTAT)
          WTOP = .TRUE.
        ELSEIF( kywrd .EQ. 'outcar' ) THEN
          READ(5,*,iostat=rstat) var,OUTCAR
           IF( rstat.NE.0) CALL PRNERROR(54,ESTAT)
          WCAR = .TRUE.
         ELSEIF( kywrd .EQ. 'outcom' ) THEN
          READ(5,*,iostat=rstat) var,OUTCOM
           IF( rstat.NE.0) CALL PRNERROR(55,ESTAT)
          WCOM = .TRUE.
        ELSEIF( kywrd .EQ. 'outmcoms' ) THEN
          READ(5,*,iostat=rstat) var,OUTMCOMS
           IF( rstat.NE.0) CALL PRNERROR(56,ESTAT)
          WMCOMS = .TRUE.
        ELSEIF( kywrd .EQ. 'inlammps'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(57,ESTAT)
          INLMP(S) = vari  
          RLMP(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'inturbomol'  ) THEN
          READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(58,ESTAT)
          INTBML(S) = vari  
          RTBML(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'outlammps' ) THEN
          READ(5,*,iostat=rstat) var,OUTLMP
           IF( rstat.NE.0) CALL PRNERROR(59,ESTAT)
           WLMP  = .TRUE. 
        ELSEIF( kywrd .EQ. 'embedded' ) THEN
          READ(5,*,iostat=rstat) var
          IF( rstat.NE.0) CALL PRNERROR(60,ESTAT)
          WEMB = .TRUE.
        ELSEIF( kywrd .EQ. 'outturbomol' ) THEN
          READ(5,*,iostat=rstat) var,OUTTBMOL
          IF( rstat.NE.0) CALL PRNERROR(60,ESTAT)
          WTBML = .TRUE.
        ELSEIF( kywrd .EQ. 'inxmol' ) THEN
           READ(5,*,iostat=rstat) var,S,vari
           IF( rstat.NE.0) CALL PRNERROR(61,ESTAT)
           INXMOL(S) = vari 
           RXMOL(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'set_ttog' ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(62,ESTAT)
          STTOG = .TRUE.
        ELSEIF( kywrd .EQ. 'seed' ) THEN
           READ(5,*,iostat=rstat) var,seed
           IF( rstat.NE.0) CALL PRNERROR(63,ESTAT)
        ELSEIF( kywrd .EQ. 'randomguess' ) THEN
          READ(5,*,iostat=rstat) var,ATMAX
           IF( rstat.NE.0) CALL PRNERROR(65,ESTAT)
        ELSEIF( kywrd .EQ. 'random' ) THEN
           READ(5,*,iostat=rstat) var,S,N
           IF( rstat.NE.0) CALL PRNERROR(65,ESTAT)
           NGAS(S) = N 
           MKRAND(S) = .TRUE. 
           CRAND = .TRUE.
       ELSEIF( kywrd .EQ. 'fixlatvec' ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(70,ESTAT)
          FIXLV(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'rotate' ) THEN
           READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(67,ESTAT)
           RANROT(S) = .TRUE.
        ELSEIF( kywrd .EQ. 'expandrand' ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(68,ESTAT)
          READ(5,*) VACEXP(1,:,S)
          IF( rstat.NE.0) CALL PRNERROR(57,ESTAT)
          READ(5,*) VACEXP(2,:,S)
          IF( rstat.NE.0) CALL PRNERROR(57,ESTAT)
          READ(5,*) VACEXP(3,:,S)
          IF( rstat.NE.0) CALL PRNERROR(57,ESTAT)
          EXPRAN(S) =.TRUE.
        ELSEIF( kywrd.EQ. 'maxminrand' ) THEN 
          READ(5,*,iostat=rstat) var,S
          IF( rstat.NE.0) CALL PRNERROR(58,ESTAT)
          RANBOX(S) =.TRUE.
        ELSEIF( kywrd .EQ. 'latconst' ) THEN
          READ(5,*,iostat=rstat) var,S,V(:),VEC(:)
           IF( rstat.NE.0) CALL PRNERROR(69,ESTAT)
          LC(1:3,S) = V(1:3)
          LA(1:3,S) = VEC(:)
          LCRD(S) =.TRUE.
        ELSEIF( kywrd.EQ. 'minbuffrand' ) THEN 
          READ(5,*,iostat=rstat) var,S,J,N
          IF( rstat.NE.0) CALL PRNERROR(60,ESTAT)
          MINRD(S) =.TRUE.
          MINML(S) =  N
          MINRS(S) = J
        ELSEIF( kywrd .EQ. 'latvec' ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(70,ESTAT)
          READ(5,*,iostat=rstat) LV(1,1:3,S)
          IF( rstat.NE.0) CALL PRNERROR(62,ESTAT)
          READ(5,*,iostat=rstat) LV(2,1:3,S)
          IF( rstat.NE.0) CALL PRNERROR(62,ESTAT)
          READ(5,*,iostat=rstat) LV(3,1:3,S)
          IF( rstat.NE.0) CALL PRNERROR(62,ESTAT)
          LVRD(S) =.TRUE.
        ELSEIF( kywrd.EQ. 'resetranbox' ) THEN 
!         Use set box size to certain value
          READ(5,*,iostat=rstat) var,S,Io,R
          IF( rstat.NE.0) CALL PRNERROR(63,ESTAT)
          RANFRX(Io,S) =.TRUE.
          RANFRSZ(Io,S) = R
        ELSEIF( kywrd .EQ. 'center' ) THEN
          READ(5,*,iostat=rstat) var,S,V(1:3)
           IF( rstat.NE.0) CALL PRNERROR(71,ESTAT)
          FCENT(1:3) = V(1:3)
          MKCENT(S) =.TRUE.
        ELSEIF( kywrd .EQ. 'vacuum' ) THEN
           READ(5,*,iostat=rstat) var,SDIM,VACL
           IF( rstat.NE.0) CALL PRNERROR(72,ESTAT)
           MKVAC(1) = .TRUE.
        ELSEIF( kywrd .EQ. 'velocity' ) THEN
          READ(5,*,iostat=rstat) var,S,VEC(:)
           IF( rstat.NE.0) CALL PRNERROR(73,ESTAT)
          VEL(:,S) = VEC(:)
          MKVEL(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'supercell' ) THEN
          READ(5,*,iostat=rstat) var,S,V(:)
           IF( rstat.NE.0) CALL PRNERROR(74,ESTAT)
          SC(1:3,S) = V(:)
          MKSUP(S) =.TRUE.
        ELSEIF( kywrd .EQ. 'mass2atnumb' ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(75,ESTAT)
          MTOAN(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'molcheck' ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(76,ESTAT)
          MOLCH(S) = .TRUE.
        ELSEIF( kywrd .EQ. 'molcmcenter' ) THEN
          READ(5,*,iostat=rstat) var,S,I,V(1:3)
           IF( rstat.NE.0) CALL PRNERROR(77,ESTAT)
          MICENT(S) = I
          MCCENT(1:3,S) = V(1:3)
          MLCENT(S) =.TRUE.
        ELSEIF( kywrd .EQ. 'molmxmncenter' ) THEN
          READ(5,*,iostat=rstat) var,S,I,V(1:3)
           IF( rstat.NE.0) CALL PRNERROR(78,ESTAT)
          MICENT(S) = I
          MCCENT(1:3,S) = V(1:3)
          MLCENT(S) =.TRUE.
          MLMNCENT(S) =.TRUE.
        ELSEIF( kywrd.EQ. 'rmatomtype' ) THEN
           READ(5,*,iostat=rstat) var,S,Ao
           IF( rstat.NE.0) CALL PRNERROR(25,ESTAT)
            N =  NATPR(S) 
            N = N + 1
            NATPR(S) = N
            RATOMTP(N,S) = Ao
            RMATTP(S) = .TRUE.
        ELSEIF( kywrd .EQ. 'cutregion' ) THEN 
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(79,ESTAT)
          DO I=1,NDIM
            READ(5,*,iostat=rstat) V(1),V(2)
            IF( rstat.NE.0) CALL PRNERROR(74,ESTAT)
            CUTFC(:,I,S) = V(1:2)
          ENDDO
          CFRAC(S) = .TRUE.
        ELSEIF( kywrd .EQ. 'cut_mol' ) THEN 
          READ(5,*,iostat=rstat) var,S,Io,Jo
           IF( rstat.NE.0) CALL PRNERROR(80,ESTAT)
          CMOLN(1,S) = Io
          CMOLN(2,S) = Jo
          CMOL(S) = .TRUE.
       ELSEIF( kywrd .EQ. 'cut_molneighbors' ) THEN 
          READ(5,*,iostat=rstat) var,S,Io,Jo
           IF( rstat.NE.0) CALL PRNERROR(80,ESTAT)
          CMOLN(1,S) = Io
          CMOLN(2,S) = Jo
          CUTMOLNB(S) = .TRUE.
       ELSEIF( kywrd .EQ. 'pbc' ) THEN
          READ(5,*,iostat=rstat) var,S,LG1(:)
           IF( rstat.NE.0) CALL PRNERROR(81,ESTAT)
          PBCS(:,S) = LG1(:)
          UPBC(:) = .TRUE. 
       ELSEIF( kywrd .EQ. 'fixregion' ) THEN
          READ(5,*,iostat=rstat) var,S,Ao,Bo,Cc
           IF( rstat.NE.0) CALL PRNERROR(82,ESTAT)
          FRTAG(:,S) = (/ Ao,Bo,Cc/)
          DO I=1,NDIM
            READ(5,*,iostat=rstat) V(1),V(2)
            IF( rstat.NE.0) CALL PRNERROR(78,ESTAT)
            RFIXC(:,I,S) = V(1:2)
          ENDDO
          REFIX(S) = .TRUE.
        ELSEIF( kywrd .EQ. 'molrdfanumbs' ) THEN
          READ(5,*,iostat=rstat) var,S,I,J,LG1(:) 
           IF( rstat.NE.0) CALL PRNERROR(83,ESTAT)
          N = NMCNT(S) 
          N = N + 1
          NMOLRi(N,S) = I   
          NMOLRj(N,S) = J   
          UMDIS(:,N,S) = LG1(:)
          MKMRDF(S) = .TRUE. 
          NMCNT(S)  = N 
        ELSEIF( kywrd .EQ. 'molrdffile'  ) THEN
          READ(5,*,iostat=rstat) var,S,ID
           IF( rstat.NE.0) CALL PRNERROR(84,ESTAT)
          N = NMCNT(S) 
          MPRDISFL(N,S) = ID
        ELSEIF( kywrd .EQ. 'vecangle' ) THEN
           READ(5,*,iostat=rstat) var,S,Ao,Bo,VEC(:)
           IF( rstat.NE.0) CALL PRNERROR(85,ESTAT)
           N = ANCNT(S) 
           N = N + 1
           ANCNT(S) = N
           ANAX(1,N,S) = Ao
           ANAX(2,N,S) = Bo
           ANGX(:,N,S) = VEC(:)
           ANGAXS(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'vecanglefile'  ) THEN
          READ(5,*,iostat=rstat) var,S,ID
           IF( rstat.NE.0) CALL PRNERROR(86,ESTAT)
           N = ANCNT(S)  
          ANGAXFL(N,S) = ID
        ELSEIF( kywrd .EQ. 'depthvecanglefile'  ) THEN
          READ(5,*,iostat=rstat) var,S,ID
           IF( rstat.NE.0) CALL PRNERROR(87,ESTAT)
          DPANGAXFL(S) = ID
        ELSEIF( kywrd .EQ. 'vectriangle' ) THEN
           READ(5,*,iostat=rstat) var,S,Ao,Bo,Co,VEC(:)
           IF( rstat.NE.0) CALL PRNERROR(88,ESTAT)
           N = TANCNT(S) 
           N = N + 1
           TANCNT(S) = N
           TANAX(1,N,S) = Ao
           TANAX(2,N,S) = Bo
           TANAX(3,N,S) = Co
           TANGX(:,N,S) = VEC(:)
           TANGAXS(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'vectrianglefile'  ) THEN
          READ(5,*,iostat=rstat) var,S,ID
           IF( rstat.NE.0) CALL PRNERROR(89,ESTAT)
           N = TANCNT(S)  
           TANGAXFL(N,S) = ID
        ELSEIF( kywrd .EQ. 'depthvectrianglefile'  ) THEN
          READ(5,*,iostat=rstat) var,S,ID
           IF( rstat.NE.0) CALL PRNERROR(90,ESTAT)
          DPTANGAXFL(S) = ID 
        ELSEIF( kywrd .EQ. 'copyregion'  ) THEN
          READ(5,*,iostat=rstat) var,S,VEC(:)
           IF( rstat.NE.0) CALL PRNERROR(91,ESTAT)
          CPRDS(:,S) = VEC(:) 
          DO I=1,NDIM
            READ(5,*,iostat=rstat) V(1),V(2)
            IF( rstat.NE.0) CALL PRNERROR(88,ESTAT)
            RCPRG(:,I,S) = V(1:2)
          ENDDO
          CPREG(S) = .TRUE.
        ELSEIF( kywrd .EQ. 'atomrdftypes' ) THEN
          READ(5,*,iostat=rstat) var,S,Ao,Bo,LG1(:)
           IF( rstat.NE.0) CALL PRNERROR(92,ESTAT)
          N = PRCNT(S) 
          N = N + 1
          RDIS(1,N,S) = Ao
          RDIS(2,N,S) = Bo
          UDIS(:,N,S) = LG1(:)
          PRDIS(S) =.TRUE.
          PRCNT(S) = N
        ELSEIF( kywrd .EQ. 'rdfdistfile'  ) THEN
          READ(5,*,iostat=rstat) var,S,ID
           IF( rstat.NE.0) CALL PRNERROR(93,ESTAT)
          N = PRCNT(S) 
          PRDISFL(N,S) = ID
       ELSEIF  ( kywrd .eq. 'nointramol' ) THEN
          READ(5,*,iostat=rstat) var,S
          IF( rstat.NE.0) CALL PRNERROR(91,ESTAT)
          RMIAM(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'slab' ) THEN
           READ(5,*,iostat=rstat) var,S,SDIM
           IF( rstat.NE.0) CALL PRNERROR(94,ESTAT)
           MKSLAB(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'rename' ) THEN
           READ(5,*,iostat=rstat) var,S,Ao,Bo
           IF( rstat.NE.0) CALL PRNERROR(95,ESTAT)
           N = RENM(S) + 1
           RNIDI(N,S) = Ao
           RNIDF(N,S) = Bo
           REN(S) = .TRUE. 
           RENM(S) = N
        ELSEIF( kywrd .EQ. 'nbcharges'  ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(96,ESTAT)
           NBCHRG(S) = .TRUE.
        ELSEIF( kywrd .EQ. 'addsurfatom' ) THEN
           READ(5,*,iostat=rstat) var,S,ID,VEC(:)
           IF( rstat.NE.0) CALL PRNERROR(97,ESTAT)
           N = SATOMN(S)  + 1
           SATOMID(N,S)  = ID
           SATOMV(:,N,S) = VEC(1:3)
           SATOMN(S) = N
        ELSEIF( kywrd .EQ. 'stack' ) THEN
           READ(5,*,iostat=rstat) var,N
           IF( rstat.NE.0) CALL PRNERROR(98,ESTAT)
           MKSTACK = .TRUE.
           NSTK = N
           DO S=1,N
             READ(5,*,iostat=rstat) I,J,R
             IF( rstat.NE.0) CALL PRNERROR(97,ESTAT)
             STKI(S)    =  I
             STKJ(S)    =  J 
             STKBUF(S)  =  R
           ENDDO
       ELSEIF( kywrd .EQ. 'hterm' ) THEN 
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(99,ESTAT)
          MKHT(S)  = .TRUE.
        ELSEIF( kywrd .EQ. 'fixstack' ) THEN
           READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(100,ESTAT)
           SHFIX(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'totaldipole' ) THEN
           READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(101,ESTAT)
           CLDIP(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'rotateall' ) THEN
           READ(5,*,iostat=rstat) var,S,AA,AB
           IF( rstat.NE.0) CALL PRNERROR(102,ESTAT)
           RANA(S) = AA
           RANB(S) = AB
           ROTA(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'moldipole' ) THEN
           READ(5,*,iostat=rstat) var,S,Ao
           IF( rstat.NE.0) CALL PRNERROR(103,ESTAT)
           CMDIP(S) = .TRUE. 
           CMDIPFL(S) = Ao 
        ELSEIF( kywrd .EQ. 'moldipdist' ) THEN
           READ(5,*,iostat=rstat) var,S,Ao
           IF( rstat.NE.0) CALL PRNERROR(104,ESTAT)
           MOLDIPFL(S)= Ao 
           CMLDPD(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'moldipaxdist' ) THEN
           READ(5,*,iostat=rstat) var,S,Io,Ao
           IF( rstat.NE.0) CALL PRNERROR(105,ESTAT)
           MDDIM(S) = Io 
           MOLDIPAXFL(S)= Ao 
           CMLDPD(S) = .TRUE.
         ELSEIF( kywrd .EQ. 'moldipbin' ) THEN
           READ(5,*,iostat=rstat) var,S,AA,AB
           IF( rstat.NE.0) CALL PRNERROR(106,ESTAT)
           DSCALE(S) =  AA
           DMAX(S) = AB
        ELSEIF( kywrd .EQ. 'indexprop' ) THEN
           READ(5,*,iostat=rstat) var,S,Io
           IF( rstat.NE.0) CALL PRNERROR(107,ESTAT)
           N = IDPN(S)
           N = N + 1
           IDPID(N,STOT) = Io 
           IDPN(S) = N
           IPROP(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'nonblist' ) THEN
           READ(5,*,iostat=rstat) var
           IF( rstat.NE.0) CALL PRNERROR(108,ESTAT)
           NBLST = .FALSE. 
        ELSEIF( kywrd .EQ. 'zeromax' ) THEN
           READ(5,*,iostat=rstat) var,S,Io
           IF( rstat.NE.0) CALL PRNERROR(109,ESTAT)
           SZDIR(S) = Io
           SZMAX(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'shift' ) THEN
           READ(5,*,iostat=rstat) var,S,VEC(:)
           IF( rstat.NE.0) CALL PRNERROR(110,ESTAT)
           SHVEC(:,S) = VEC(:)
           SHFT(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'allmolcom' ) THEN
           READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(111,ESTAT)
           WAMONIOM(S) = .TRUE. 
           WAMONIOMA(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'moleculecom' ) THEN
           READ(5,*,iostat=rstat) var,S,I
           IF( rstat.NE.0) CALL PRNERROR(112,ESTAT)
!           WAMONIOM(S) = .TRUE. 
           WAMONIOMM(S) = .TRUE. 
           N = NONMOL(S)
           N = N + 1
           WONIONMOL(N,S) = I 
           NONMOL(S) = N
        ELSEIF( kywrd .EQ. 'qmsphere' ) THEN
           READ(5,*,iostat=rstat) var,S,R,AA,VEC
           IF( rstat.NE.0) CALL PRNERROR(113,ESTAT)
           WONIOMSP(S) = .TRUE. 
           ONRAD(1,S) = R
           ONRAD(2,S) = AA
           ONCENT(:,S) = VEC(:) 
           USONIOMFX(S) = .TRUE.

! debug hack
           WONIOMSP(STOT) = .TRUE. 
           ONRAD(1,STOT) = R
           ONRAD(2,STOT) = AA
           ONCENT(:,STOT) = VEC(:) 
           USONIOMFX(STOT) = .TRUE.


        ELSEIF( kywrd .EQ. 'qmnieghbors' ) THEN
           READ(5,*,iostat=rstat) var,S,J
           IF( rstat.NE.0) CALL PRNERROR(114,ESTAT)
           QMNB(S) = J

! debug hack
           QMNB(STOT) = J

        ELSEIF( kywrd .EQ. 'atomefield' ) THEN
           READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(115,ESTAT)
           ATEFL(S) = .TRUE. 
        ELSEIF( kywrd .EQ. 'molefield' ) THEN
           READ(5,*,iostat=rstat) var,S,Ao
           IF( rstat.NE.0) CALL PRNERROR(116,ESTAT)
           MEFL(S) = .TRUE. 
           MEFLFL(S) = Ao 
        ELSEIF( kywrd .EQ. 'delR' ) THEN
           READ(5,*,iostat=rstat) var,S,AA
           IF( rstat.NE.0) CALL PRNERROR(117,ESTAT)
           RSCALE(S) = AA
        ELSEIF( kywrd .EQ. 'deldip' ) THEN
           READ(5,*,iostat=rstat) var,S,AA,AB
           IF( rstat.NE.0) CALL PRNERROR(118,ESTAT)
           DSCALE(S) = AA
           DMAX(S) = AB 
      ELSEIF( kywrd .EQ. 'noqmmmfixed' ) THEN
           READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(119,ESTAT)
           USONIOMFX(S) = .FALSE.
      ELSEIF( kywrd .EQ. 'g_molmin' ) THEN
           READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(120,ESTAT)
           GMMIN(S) = .TRUE.
      ELSEIF( kywrd .EQ. 'update_charges' ) THEN
           READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(121,ESTAT)
           CHR_UPDATE(S) = .TRUE.
      ELSEIF( kywrd .EQ. 'xyz_trj' ) THEN
           READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(122,ESTAT)
           MFRAM_XYZ(S) = .TRUE.
      ELSEIF( kywrd .EQ. 'atom_angle_dist' ) THEN
           READ(5,*,iostat=rstat) var,S,Ao,Bo,Co
           IF( rstat.NE.0) CALL PRNERROR(123,ESTAT)
           N = ATOM_ANGL_CNT(S) 
           N = N + 1
           ATOM_ANGL_CNT(S) = N
           ATOM_ANAX(1,N,S) = Ao 
           ATOM_ANAX(2,N,S) = Bo
           ATOM_ANAX(3,N,S) = Co
           CALC_ATOM_ANG(S) = .TRUE. 
      ELSEIF( kywrd .EQ. 'atom_angle_file'  ) THEN
          READ(5,*,iostat=rstat) var,S,ID
           IF( rstat.NE.0) CALL PRNERROR(124,ESTAT)
           N = ATOM_ANGL_CNT(S)  
           ATOM_ANFILE(N,S) = ID
      ELSEIF( kywrd .EQ. 'cov_buffer'  ) THEN
          READ(5,*,iostat=rstat) var,AA
           IF( rstat.NE.0) CALL PRNERROR(125,ESTAT)
           CRBUF = 1.0d0 + AA 
      ELSEIF( kywrd .EQ. 'read_bconst'  ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(127,ESTAT)
          RDBCNST(S) = .TRUE.
      ELSEIF( kywrd .EQ. 'read_bconst'  ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(127,ESTAT)
          RDBCNST(S) = .TRUE.
      ELSEIF( kywrd .EQ. 'read_aconst'  ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(128,ESTAT)
          RDACNST(S) = .TRUE.
      ELSEIF( kywrd .EQ. 'read_dconst'  ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(129,ESTAT)
          RDDCNST(S) = .TRUE.
      ELSEIF( kywrd .EQ. 'read_iconst'  ) THEN
          READ(5,*,iostat=rstat) var,S
           IF( rstat.NE.0) CALL PRNERROR(130,ESTAT)
          RDICNST(S) = .TRUE. 
      ELSEIF( kywrd .EQ. 'deposit' ) THEN
           READ(5,*,iostat=rstat) var,S,N,SDIM,AB
           IF( rstat.NE.0) CALL PRNERROR(131,ESTAT)
           NDEP(S)  =   N 
           ADDDM = .TRUE. 
           ADDDMOL(S) = .TRUE. 
           DEPBF(S) =   AB
      ELSEIF( kywrd.EQ. 'depbox' ) THEN 
!         Use maximum and minmum of other structures as box size
          READ(5,*,iostat=rstat) var,S,N
          IF( rstat.NE.0) CALL PRNERROR(50,ESTAT)
          UDBOX(S) = .TRUE.
          RSTBOX(S) = N
      ELSEIF( kywrd.EQ. 'deptemperature' ) THEN 
          READ(5,*,iostat=rstat) var,S,AA
          IF( rstat.NE.0) CALL PRNERROR(50,ESTAT)
          DEPTEMP(S) = AA
          MKDTMP(S) = .TRUE. 
      ELSEIF( kywrd .EQ. 'tnknumbmod'  ) THEN
          READ(5,*,iostat=rstat) var,S,Io,Jo,To 
           IF( rstat.NE.0) CALL PRNERROR(132,ESTAT)
          TNKNMD(S) = .TRUE. 
          TNKNUMBM(1,S) = Io
          TNKNUMBM(2,S) = Jo
          TNKNUMBM(3,S) = To
      ELSEIF( kywrd .EQ. 'tnkmolmod'  ) THEN
          READ(5,*,iostat=rstat) var,S,Io,To 
           IF( rstat.NE.0) CALL PRNERROR(133,ESTAT)
          MTNKNMD(S) = .TRUE. 
          MTNKNUMBM(1,S) = Io
          MTNKNUMBM(2,S) = To
      ELSEIF( kywrd .EQ. 'set_oniomcm2'  ) THEN
          READ(5,*,iostat=rstat) var,S,LAYER_CHR(1),LAYER_MLT(1)
     &  ,LAYER_CHR(2),LAYER_MLT(2)
           IF( rstat.NE.0) CALL PRNERROR(130,ESTAT)
          ONIOMLAYER = 2
      ELSEIF( kywrd .EQ. 'set_oniomcm3'  ) THEN
          READ(5,*,iostat=rstat) var,S,LAYER_CHR(1),LAYER_MLT(1)
     &  ,LAYER_CHR(2),LAYER_MLT(2),LAYER_CHR(3),LAYER_MLT(3)
          IF( rstat.NE.0) CALL PRNERROR(130,ESTAT)
          ONIOMLAYER = 3
      ELSEIF( kywrd .EQ. 'set_oniom'  ) THEN
          READ(5,*,iostat=rstat) var,S
          IF( rstat.NE.0) CALL PRNERROR(130,ESTAT)
          READ(5,'(A100)') ONIOM_METH
      ELSEIF( kywrd .EQ. 'set_ptype'  ) THEN
          READ(5,*,iostat=rstat) var,PNUMB,TP,MS,CR,VR
          IF( rstat.NE.0) CALL PRNERROR(130,ESTAT)
          ATSYM( PNUMB ) = TP
          AMASS( PNUMB ) = MS
          CRDI( PNUMB )  = CR
          VDWRDI( PNUMB ) = VR
      ELSEIF( kywrd .EQ. 'get_bonds'  ) THEN
          READ(5,*,iostat=rstat) var,S
          IF( rstat.NE.0) CALL PRNERROR(130,ESTAT)
          GBONDS(S) = .TRUE.
      ELSEIF( kywrd .EQ. 'get_angles'  ) THEN
          READ(5,*,iostat=rstat) var,S
          IF( rstat.NE.0) CALL PRNERROR(130,ESTAT)
          GANGLES(S) = .TRUE.
      ELSEIF( kywrd .EQ. 'get_dihedrals'  ) THEN
          READ(5,*,iostat=rstat) var,S
          IF( rstat.NE.0) CALL PRNERROR(130,ESTAT)
          GDIHS(S) = .TRUE.
      ELSEIF  ( kywrd .eq. '!' ) THEN
          READ(5,*)
      ELSEIF  ( kywrd .eq. 'C ' ) THEN
          READ(5,*)
        ELSE 
          IF (VERB) WRITE(6,*) kywrd,' unknown' 
          READ(5,*,iostat=istat)
        ENDIF
        READ(5,*,iostat=istat) kywrd
        BACKSPACE(5) 
      ENDDO
!     Check input 
      DO STN = 1,STOTo
          IF( RXYZ(STN) ) THEN
            IF( (.NOT.LVRD(STN)).AND.(.NOT.LCRD(STN))) THEN
            WRITE(6,*) 'Warning: Reading xyz with no lattice info'
            WRITE(6,*)'Specify either latvec or latconst and latangle'
            ENDIF
          ENDIF
          IF( RPDB(STN) ) CLLV(STN) = .TRUE.
          IF( RCAR(STN) ) CLLC(STN) = .TRUE.
          IF( RGRO(STN) ) CLLC(STN) = .TRUE.
          IF( RTOP(STN) ) CLMAS(STN) = .FALSE.
          IF( LCRD(STN) ) THEN
             CLLC(STN) = .FALSE.
             CLLV(STN) = .TRUE.
          ENDIF
          IF( LVRD(STN) )  THEN
             CLLV(STN) = .FALSE.
             CLLC(STN) = .TRUE.
          ENDIF
!         Check redundancies 
          IF( CLLV(STN).AND.CLLC(STN) ) THEN
             CALL PRNERROR(73,ESTAT)
          ENDIF 
!         Turn on mass check if top and gro are read in
          IF( RTOP(STN).AND.RGRO(STN) ) MTOAN(STN) = .TRUE. 
!         Recalc molecules and neighbors if atom is removed
          IF( RMATTP(S) ) THEN
             NBLST = .TRUE.
             MOLCH(S) = .TRUE. 
          ENDIF 
!         
      ENDDO
!
      IF ( HELP ) THEN
        DO STN = 1,STOTo
            IF( VERB ) THEN
              WRITE(*,*) ' Verbose output will '
            ELSEIF( DEBUG ) THEN
              WRITE(*,*) ' Debug will output fort.* ',
     &         'files for debugging ',
     &         ' grep -in "debug" to find where data is generated '
            ELSEIF( RXYZ(STN) ) THEN
              WRITE(*,*) '  '
            ELSEIF( RXYZ(STN) ) THEN
              WRITE(*,*) '  '
            ELSEIF( RXYZ(STN) ) THEN
              WRITE(*,*) '  '
            ENDIF
        ENDDO  
      ENDIF 

!
  101 FORMAT(I6,A6,' atoms will be added to',3I6)
      END SUBROUTINE readin
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Write out inital structure properties 
!
      SUBROUTINE  WRITE_PROPo
      USE structure 
      USE specify
      USE potential
      USE build 
!
      IMPLICIT none
!
      INTEGER :: STN,N
      WRITE(6,*) "-------------------------------------------"
!     Summereize input information
      WRITE(6,622) SEED
      WRITE(6,601) STOTo
      DO STN=1,STOTo
        WRITE(6,602) STN
       IF (RXYZ(STN) ) THEN
          WRITE(6,603) INXYZ(STN)
        ENDIF
        IF (RPDB(STN)) THEN
          WRITE(6,603) INPDB(STN)
        ENDIF
        CALL WRITE_STRP(STN)
        IF( PRCNT(STN).GT.0 ) THEN          
          WRITE(6,630) PRCNT(STN)
          DO N=1,PRCNT(STN)
            WRITE(6,631) N,RDIS(1,N,STN),RDIS(2,N,STN) 
            WRITE(6,632) UDIS(:,N,STN)
            WRITE(6,633)  PRDISFL(N,STN)
          ENDDO 
        ENDIF
        IF( RMIAM(STN) ) THEN          
          WRITE(6,634)  
        ENDIF 
        IF( REN(STN) ) THEN
          DO N =1,RENM(STN)
             WRITE(6,*) " Changing ",RNIDI(N,STN)," to ", RNIDF(N,STN)
          ENDDO 
        ENDIF 
        IF( SZMAX(STN)) THEN
          WRITE(6,*) 'Shifting the maximum in the '
     &   ,SZDIR(STN),' to zero'
        ENDIF
!       Find inion cut out
c$$$        IF ( CUTONION(STN) ) THEN
c$$$          WRITE(6,620) NAo(STN)/52,ONi(STN)
c$$$          WRITE(6,621) ONCUT(STN)
c$$$        ENDIF 
c$$$        IF( ADTNC(STN).GT.0 ) THEN
c$$$          WRITE(6,623) ADTNC(STN)
c$$$          DO I=1,ADTNC(STN)
c$$$            WRITE(6,624) ADTCONI(STN,I),ADTCONJ(STN,I)
c$$$          ENDDO
c$$$        ENDIF
c$$$        IF( CFRAC(STN) ) THEN
c$$$          WRITE(6,628) 
c$$$          WRITE(6,629) CUTFC(STN,1,:)
c$$$          WRITE(6,629) CUTFC(STN,2,:)
c$$$          WRITE(6,629) CUTFC(STN,3,:)
c$$$        ENDIF
      ENDDO
!
 623  FORMAT("    ",I12," constraints will be added")
 624  FORMAT("         Bassed on atom ",A6," having added atom ",A6)
 601  FORMAT(I12,' structures will be read in')
 602  FORMAT('Structure ',I12,' ')
 603  FORMAT('    Sturcture read from ',A)
 604  FORMAT('    Concetivity read from ',A)
 620  FORMAT('      Found',I12,' molecules around atom',I12)
 621  FORMAT('      Using a cut off of ',F12.4)
 622  FORMAT('A random seed number of ',I12,' will be used')
 628  FORMAT("    Removing atoms not in fractional regions:")
 629  FORMAT("      ",2F16.6)
 630  FORMAT(' ',I12,' atomic RDFs will be printed')
 631  FORMAT('    REF ',I12,' is between ',2A6)
 632  FORMAT('      Dimensions ',3L)
 633  FORMAT('      and will be writen to file :',A40)
 634  FORMAT('    Intramolecular atoms will be neglected for RDFs ')
!
      RETURN 
      END SUBROUTINE WRITE_PROPo
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Write out inital structure properties 
!
      SUBROUTINE write_propi
      USE structure 
      USE specify
      USE potential
!
      IMPLICIT none
!
      INTEGER :: STN,I
!     Summereize final str information
      WRITE(6,*) "-------------------------------------------"
      DO STN=STOTo+1,STOT-1
        WRITE(6,602) STN-STOTo
        CALL WRITE_STRP(STN)
c$$$!       Find inion cut out
c$$$        IF ( CUTONION(STN) ) THEN
c$$$          WRITE(6,620) NAi(STN)/52,ONi(STN)
c$$$          WRITE(6,621) ONCUT(STN)
c$$$        ENDIF 
c$$$        IF( ADTNC(STN).GT.0 ) THEN
c$$$          WRITE(6,623) ADTNC(STN)
c$$$          DO I=1,ADTNC(STN)
c$$$            WRITE(6,624) ADTCONI(STN,I),ADTCONJ(STN,I)
c$$$          ENDDO
c$$$        ENDIF

      ENDDO
!
 623  FORMAT("    ",I12," constraints will be added")
 624  FORMAT("         Bassed on atom ",A6," having added atom ",A6)
 601  FORMAT(I12,' structures will be read in')
 602  FORMAT('Final structure ',I12,' ')
 620  FORMAT('      Found',I12,' molecules around atom',I12)
 621  FORMAT('      Using a cut off of ',F12.4)
 622  FORMAT('A random seed number of ',I12,' will be used')
!
      RETURN 
      END SUBROUTINE write_propi
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Write out structure properties 

      SUBROUTINE WRITE_STRP(STN)
      USE structure 
      USE potential 
      USE specify 
      USE build 
!
      IMPLICIT none
      INTEGER :: STN,E
!
        WRITE(6,605) NA(STN)
!        WRITE(6,625) NCo(STN)
        WRITE(6,606) NB(STN)
!        WRITE(6,607) NPRo(STN)
        WRITE(6,608) NANG(STN)
        WRITE(6,609) NDIH(STN)
        WRITE(6,610)
        WRITE(6,611) LV(1,1:3,STN)
        WRITE(6,611) LV(2,1:3,STN)
        WRITE(6,611) LV(3,1:3,STN)
        WRITE(6,612) 
        WRITE(6,613) LC(1:3,STN)
        WRITE(6,614) LA(1:3,STN)
        WRITE(6,615) VOL(STN)
        WRITE(6,616) TMAS(STN)
        WRITE(6,617) DEN(STN)
        WRITE(6,618) MMN(2,:,STN)
        WRITE(6,619) MMN(1,:,STN) 
        WRITE(6,627) ELT(STN)
        DO E=1,ELT(STN)
          WRITE(6,626) NEL(E,STN),ELCNT(E,STN)
        ENDDO
        WRITE(6,633) TCHG(STN) 
        IF(  CLDIP(STN)) WRITE(6,634) DIP(STN) 
        IF( MKSLAB(STN).OR.MKVAC(1) ) THEN 
          WRITE(6,630) SHT(STN)
          WRITE(6,631) SVOL(STN)
          WRITE(6,632) SDEN(STN)
        ENDIF
        IF( WTBML ) THEN
          WRITE(6,635)
          IF( QMNB(STN).GT.0 )WRITE(6,636)
          IF( WONIOMSP(STN) ) WRITE(6,637) ONRAD(2,STN),ONRAD(1,STN)
        ENDIF
      RETURN
 605  FORMAT('    ',I12,' atoms')
 607  FORMAT('    ',I12,' pairs')
 625  FORMAT('    ',I12,' constraints')
 606  FORMAT('    ',I12,' bonds')
 608  FORMAT('    ',I12,' angles')
 609  FORMAT('    ',I12,' dihedrals')
 610  FORMAT('    Lattice vectors:')
 611  FORMAT('    ',3F12.4)
 612  FORMAT('    Lattice constants:')
 613  FORMAT('      A,B,C           :',3F12.4)
 614  FORMAT('      alpha,beta,gamma:',3F12.4)
 615  FORMAT('      Volume :',F22.4,' A^3')
 616  FORMAT('      Mass   :',F22.4,' AMU')
 617  FORMAT('      Density:',F22.4,' g/cm^3')
 618  FORMAT('      Max(A):',3F12.4)
 619  FORMAT('      Min(A):',3F12.4)
 626  FORMAT("      Atomic #",I12," has  ",I12," atoms")
 627  FORMAT("      Number of elements",I12)
 630  FORMAT('      Surface Height :',F22.4,' A')
 631  FORMAT('      Surface Volume :',F22.4,' A^3')
 632  FORMAT('      Surface Density:',F22.4,' g/cm^3')
 633  FORMAT("      Charge :",F22.4)
 634  FORMAT('      Total dipole moment :',F22.4,' D')
 635  FORMAT('      Writing QM/MM files for each molecule')
 636  FORMAT('        with ',I6
     & ,' neighbor shells included in the QM region ')
 637  FORMAT('        a sphere of ',F12.3
     & ,' will cut and included  in the MM region'
     & ,' with molecule greater than ',F12.3,' being fixed')
      RETURN 
      END SUBROUTINE WRITE_STRP
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Write out inital structure properties 
!
      SUBROUTINE  WRITE_PROC 
      USE structure 
      USE specify
      USE potential
      USE build 
!
      IMPLICIT none
!
      INTEGER :: STN,I,D
!
      DO STN = 1, STOTo
        IF (RXYZ(STN) ) THEN
          WRITE(6,603) INXYZ(STN)
        ENDIF
        IF (RPDB(STN)) THEN
          WRITE(6,603) INPDB(STN)
        ENDIF
        IF (RCRD(STN)) THEN
          WRITE(6,603) INCRD(STN)
        ENDIF
        IF (RGRO(STN)) THEN
          WRITE(6,603) INGRO(STN)
        ENDIF
       IF (RCAR(STN)) THEN
          WRITE(6,603) INCAR(STN)
        ENDIF
        IF ( RTOP(STN) ) THEN
          WRITE(6,604) INTOP(STN)
        ENDIF
        IF( MKRAND(STN )) THEN
          WRITE(6,607) STN
          IF( MINRD(STN) ) THEN
             WRITE(6,608) MINBF(STN),STN,MINRS(STN),SDIM,MINML(STN)
          ENDIF
          DO D=1,NDIM
          IF( RANFRX(D,STN)) WRITE(6,*)"Resetting box in",D
     &    ,"to",RANFRSZ(D,STN)
          ENDDO 
        ENDIF 

      ENDDO 
        IF (WXYZ ) THEN
          WRITE(6,605) OUTXYZ
        ENDIF
        IF (WPDB) THEN
          WRITE(6,605) OUTPDB
        ENDIF
        IF (WCRD) THEN
          WRITE(6,605) OUTCRD
        ENDIF
        IF (WGRO) THEN
          WRITE(6,605) OUTGRO
        ENDIF
       IF (WCAR) THEN
          WRITE(6,605) OUTCAR
        ENDIF
        IF ( WTOP ) THEN
          WRITE(6,606) OUTTOP
        ENDIF
 
!
      RETURN 
 603  FORMAT('    Sturcture read from ',A)
 604  FORMAT('    Concetivity read from ',A)
 605  FORMAT('    Sturcture output ',A)
 606  FORMAT('    Concetivity output ',A)
 607  FORMAT("      A random gas will be created based "
     & ," on structure",I6)
 608  FORMAT("        A minimum seperation of ",F8.2
     & , "between ",I6," and ",I6," will be inforced in the "
     & ,I6," dimension. for the first",I6," molecules ")
      END SUBROUTINE WRITE_PROC
