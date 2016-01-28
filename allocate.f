!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/02/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Allocate arrays
!
      SUBROUTINE ALLOCATE
      USE specify
      USE structure
      USE elements
      USE potential
      USE build 
!
      IMPLICIT none
      INTEGER :: NAr,istat,NBMX,STN,E,io,tstat,rstat
     & ,I,J,K,L,Ia,In,Id,Ib,Ii,NNMX,NDMX,NBLMX,ML,STNi
     & ,NBMULT ,CPCNT,NBMULTLJ,NBLMXLJ,NBMULTML,NBLMXML
     & ,NBLJSTOT,NBSTOT,Imol 
      INTEGER, ALLOCATABLE :: ELISTL(:)
      REAL*8 :: XX,YY,ZZ
      CHARACTER(CHSZ) :: RLINE,LLINE,LI,var,kywd,FIL2
      CHARACTER(IDSZ) :: AT
!
      NAMX = 0
      Ia = 0
      Ib = 0 
      In = 0
      Id = 0
      Imol = 0 
      CPCNT = 0 
      DO STN=1,STOTo - STADD 
         IF( RXYZ(STN) ) THEN
           OPEN(UNIT=11,FILE=INXYZ(STN),STATUS='unknown')
           READ(11,*) NAr
           CLOSE(11) 
           NA(STN) = NA(STN) + NAr
         ENDIF
         IF( RXMOL(STN) ) THEN
           OPEN(UNIT=19,FILE=INXMOL(STN),STATUS='unknown')
           READ(19,*) NAr
           CLOSE(19) 
           NA(STN) = NA(STN) + NAr
         ENDIF
         IF( RPDB(STN) ) THEN
           OPEN(unit=12,FILE=INPDB(STN),STATUS='unknown')
    !      Find number of atoms 
           NAr =  NA(STN)
           READ(12,*,iostat=istat) RLINE
           BACKSPACE(12)
           DO WHILE (ISTAT.EQ.0)
             IF( RLINE(1:6) .EQ. 'ATOM' ) THEN
               NAr = NAr + 1
             ENDIF
             IF( RLINE(1:6) .EQ. 'KATM' ) THEN
               NAr = NAr + 1
             ENDIF
              IF( RLINE(1:6) .EQ. 'HETATM' ) THEN
               NAr = NAr + 1
             ENDIF
             READ(12,*,iostat=istat) RLINE
    !        BACKSPACE(12)
           ENDDO
           CLOSE(12)
           NA(STN) = NAr
         ENDIF
         IF( RGRO(STN) ) THEN
           OPEN(UNIT=11,FILE=INGRO(STN),STATUS='unknown')
           READ(11,*) 
           READ(11,*) NAr
           CLOSE(11) 
           NA(STN) = NA(STN) + NAr 
         ENDIF
         IF( RTBML(STN) ) THEN
          NAr = 0
          OPEN (UNIT=11,FILE=INTBML(STN),STATUS='unknown')     
          READ(11,*,iostat=istat) RLINE
          DO WHILE (ISTAT.EQ.0)
            LLINE = ADJUSTL(RLINE)
            IF( LLINE.EQ.'$coord' ) THEN 
               READ(11,*,iostat=rstat) XX,YY,ZZ,AT
               DO WHILE ( RSTAT.EQ.0 )
                  NAr = NAr + 1
                  READ(11,*,iostat=rstat) XX,YY,ZZ,AT
               ENDDO           
               BACKSPACE(11)
               READ(11,*,iostat=istat) RLINE      
               LLINE = ADJUSTL(RLINE)
            ELSE
              READ(11,*,iostat=istat) RLINE      
              LLINE = ADJUSTL(RLINE)
            ENDIF
          ENDDO
          CLOSE(11) 
          NA(STN) = NA(STN) + NAr 
         ENDIF

         IF( RCAR(STN) ) THEN
           NAr = NA(STN)
           OPEN(UNIT=16,FILE=INCAR(STN),STATUS='unknown')
           DO E=1,6
               READ(16,*)
           ENDDO
           DO E=NELM,1,-1
             ALLOCATE( ELISTL(E) )
             READ(16,*,iostat=io) ELISTL(:)
             IF( io.EQ.0) THEN
               NAr = NAr + SUM(ELISTL)
               EXIT
             ENDIF
             BACKSPACE(16)
             BACKSPACE(16)
             DEALLOCATE( ELISTL  )
           ENDDO
           CLOSE(16)
           NA(STN) = NAr 
         ENDIF
         IF( RCRD(STN) ) THEN
           OPEN(UNIT=11,FILE=INCRD(STN),STATUS='unknown')
           READ(11,*) 
           READ(11,*) NAr
           READ(11,*) 
           CLOSE(11) 
           NA(STN) = NA(STN) + NAr 
         ENDIF
         IF( RCOM(STN) ) THEN
           WRITE(*,*) ' Gaussian input read in not supported'
           STOP 
         ENDIF
         IF( RTOP(STN)  ) THEN 
            Ia = 0
            Ib = NB(STN)
            In = NANG(STN)
            Id = NDIH(STN)
            Imol = NMOL(STN)
!
            OPEN (UNIT=16,FILE=INTOP(STN),STATUS='unknown')     
            READ(16,*,iostat=tstat) RLINE
            BACKSPACE(16)
            DO WHILE (TSTAT.EQ.0)
               LLINE = ADJUSTL(RLINE)
               LI = LLINE(1:1)      
               IF( LLINE.EQ.'#include' ) THEN 
                   READ(16,*,iostat=rstat) var,FIL2
                   IF(VERB)WRITE(*,*) 'Opeing sub top file  :',FIL2
                   OPEN (UNIT=17,FILE=FIL2,STATUS='unknown')    
                   READ(17,*,iostat=istat) RLINE
                   BACKSPACE(17)
                   DO WHILE (ISTAT.EQ.0)
                     LLINE = ADJUSTL(RLINE)
                     LI = LLINE(1:1) 
                     IF ( LI.EQ.'[' ) THEN
                        Io=SCAN( RLINE,'[' )      
                        Ii=SCAN( RLINE,']')      
                        IF ( Ii.LT.Io ) THEN
                           READ(17,*) var,KYWD 
                        ELSE 
                           KYWD = LLINE(Io+1:Ii-1)
                           READ(17,*)
                        ENDIF
                        IF( KYWD.EQ.'moleculetype' ) THEN
                          Imol  = Imol  + 1 
                        ENDIF
                     ELSE
                        READ(17,*)
                     ENDIF
                     READ(17,*,iostat=istat) RLINE
                     BACKSPACE(17)
                  ENDDO
               ELSEIF ( LI.EQ.'[' ) THEN
                  Io=SCAN( RLINE,'[' )      
                  Ii=SCAN( RLINE,']')      
                  KYWD = LLINE(Io+1:Ii-1)
                  IF ( Ii.LE.Io ) THEN
                     READ(16,*) var,KYWD 
                  ELSE 
                     READ(16,*)
                  ENDIF 
                 IF( KYWD.EQ.'moleculetype' ) THEN
                     READ(16,*) 
                     Imol  = Imol  + 1 
                 ELSEIF ( KYWD.EQ.'atoms' ) THEN          
                     RSTAT = 0 
                     READ(16,*,iostat=rstat) I,AT
                     DO WHILE ( RSTAT.EQ.0 )
                        Ia = Ia + 1
                        READ(16,*,iostat=rstat) I,AT
                     ENDDO           
                     BACKSPACE(16)   
                  ELSEIF( KYWD.EQ.'bonds' ) THEN
                     READ(16,*,iostat=rstat) I,J
                     DO WHILE ( RSTAT.EQ.0 )
                        Ib = Ib + 1
                       READ(16,*,iostat=rstat) I,J
                     ENDDO           
                     BACKSPACE(16)   
                  ELSEIF( KYWD.EQ.'angles' ) THEN
                     READ(16,*,iostat=rstat) I,J,K
                     DO WHILE ( RSTAT.EQ.0 )
                        In = In + 1
                        READ(16,*,iostat=rstat) I,J,K
                     ENDDO           
                     BACKSPACE(16)   
                  ELSEIF( KYWD.EQ.'dihedrals' ) THEN
                     READ(16,*,iostat=rstat) I,J,K,L
                     DO WHILE ( RSTAT.EQ.0 )
                         Id = Id + 1
                        READ(16,*,iostat=rstat) I,J,K,L
                     ENDDO           
                     BACKSPACE(16)   
                  ELSEIF( KYWD.EQ.'constraints' ) THEN
                  ELSEIF( KYWD.EQ.'exclusions' ) THEN
                  ELSEIF( KYWD.EQ.'virtual' ) THEN
                  ENDIF 
               ELSE
                  READ(16,*)
               ENDIF
               READ(16,*,iostat=tstat) RLINE
               BACKSPACE(16)
            ENDDO
            CLOSE(16)
            NB(STN)   =  Ib
            NANG(STN) =  In
            NDIH(STN) =  Id
            NMOL(STN) =  Imol 
        ENDIF 
!     Added structure props
        IF( CPREG(STN) ) THEN
           CPCNT = CPCNT + 1
           NA(STOTo - CPCNT +1 ) = NA(STN) 
           STCPN(STN) = STOTo - CPCNT+1 
         ENDIF
      ENDDO 
!     Find max number of atoms
      NAMX = 0
      NBMX = 0
      NNMX = 0
      NDMX = 0
      DO STN=1,STOTo
        NAr = NA(STN) 
        IF( MKRAND(STN) ) NAr = NAr*NGAS(STN)
        IF( MKSUP(STN) ) THEN
         ML =  CEILING( SC(1,STN)*SC(2,STN)*SC(3,STN) )
         NAr = NAr*ML
        ENDIF
        IF( MKHT(STN) ) NAr = NAr*4              ! If all C then max C*3 H
        IF( SATOMN(STN).GT.0 ) NAr = NAr*SATOMN(STN)  !if 1 layer max x2 atoms 
        IF( MKSLAB(STN) ) NAr = NAr*2            !if 1 layer max x2 atoms         
        NAMX = NAMX + NAr
!       IF( NAr .GT. NAMX) NAMX = NAr
        IF( NB(STN).GT.NBMX) NBMX = NB(STN)
        IF( NANG(STN).GT.NNMX) NNMX = NANG(STN)
        IF( NDIH(STN).GT.NDMX) NDMX = NDIH(STN)
      ENDDO
!     
      IF(VERB) WRITE(6,*)" Max atoms",NAMX 

!  debug
        NBMX = NAMX*12
        NNMX = NBMX*2
        NDMX = NNMX*4
!     
      ALLOCATE( 
     &  ELN(NAMX,STOT)
     & ,TYP(NAMX,STOT),TYP_IND(NAMX,STOT)
     & ,TYP_REF(NAMX,STOT),TYP_CNT(NAMX,STOT)
     & ,REND_IND(NAMX,STOT)
     & ,RID(NAMX,STOT),RN(NAMX,STOT)
     & ,MOLN(NAMX,STOT),ACHG(NAMX,STOT),CHN(NAMX,STOT),GRID(NAMX,STOT)
     & ,THRM(NAMX,STOT)
     & ,R0(NDIM,NAMX,STOT),R1(NDIM,NAMX,STOT),R2(NDIM,NAMX,STOT)
     & ,FTAG(NDIM,NAMX,STOT),ELIST(NAMX,NELM,STOT)
     & ,BNDI(NBMX,STOT),BNDJ(NBMX,STOT)
     & ,BTYP(NBMX,STOT),BVAL(NBMX,STOT),BCNST(NBMX,STOT)
     & ,ANGI(NNMX,STOT),ANGJ(NNMX,STOT),ANGK(NNMX,STOT)
     & ,ATYP(NNMX,STOT),AVAL(NNMX,STOT),ACNST(NNMX,STOT)
     & ,DIHI(NDMX,STOT),DIHJ(NDMX,STOT)
     & ,DIHK(NDMX,STOT),DIHL(NDMX,STOT)
     & ,DTYP(NDMX,STOT),DVAL(NDMX,STOT),DCNST(NDMX,STOT)
     & ,IMPI(NDMX,STOT),IMPJ(NDMX,STOT)
     & ,IMPK(NDMX,STOT),IMPL(NDMX,STOT)
     & ,ITYP(NDMX,STOT),IVAL(NDMX,STOT),ICNST(NDMX,STOT)
     & ,AMAS(NAMX,STOT),ELREF(NAMX,STOT)
     & ,MPNT(NAMX+1,STOT),MOLST(NAMX+1,STOT),MOLCNT(STOT) !molecule mass
     & ,MCOMS(NDIM,NAMX+1,STOT) , MOLMAS(NAMX+1,STOT)       !molecule mass
     & ,MOLCHR(NAMX+1,STOT)     !molecule charge
     & ,MOLDIP(NDIM,NAMX+1,STOT)                   ! molecular dipole 
     & ,MOLEF(NDIM,NAMX+1,STOT),MOLEP(NAMX+1,STOT) ! molecular efield and E pot
     & ,ATOMEF(NDIM,NAMX+1,STOT),ATOMEP(NAMX+1,STOT) ! molecular efield and E pot
     & ,GPROP(8,NAMX,STOT)        ! index properties 
     & ,ONMTAG(NAMX,STOT)  ! oniom tag 
     & ,TYPN(NAMX,STOT)    ! atom type number for tinker
     &  )
!
!     Assume max 12 neighbors 
      NBMULT = 24
      NBLMX =NAMX*NBMULT
      
      ALLOCATE ( NBTOT(STOT),NINDX(NAMX+1,STOT),NLIST(NBLMX,STOT)
     & ,NBMAX(STOT) )
      NBTOT(:) = 0
      NINDX(:,:) = 0
      NLIST(:,:) = 0
      NBMAX(:) = 1 
!     LJ nieghbor list
      NBMULTLJ = 64
      NBLMXLJ =NAMX*NBMULTLJ
      ALLOCATE (NBTOTLJ(STOT),NINDXLJ(NAMX+1,STOT),NLISTLJ(NBLMXLJ,STOT)
     & ,NBMAXLJ(STOT) )
      NBTOTLJ(:) = 0
      NINDXLJ(:,:) = 0
      NLISTLJ(:,:) = 0
      NBMAXLJ(:) = 1 
!     ML nieghbor list
      NBLMXML =NAMX*NBMULTLJ
      ALLOCATE (NBTOTML(STOT),NINDXML(NAMX+1,STOT),NLISTML(NBLMXML,STOT)
     & ,NBMAXML(STOT) )
      NBTOTML(:) = 0
      NINDXML(:,:) = 0
      NLISTML(:,:) = 0
      NBMAXML(:) = NBLMXML 
      NBLJSTOT = 0
      NBSTOT = 0
!     Set  NBmax 
      DO STN=1,STOTo
        STNi = STN + STOTo 
        NAr = NA(STN) 
        IF( MKRAND(STN) ) NAr = NAr*NGAS(STN)
        IF( MKSUP(STN) ) THEN
         ML =  CEILING( SC(1,STN)*SC(2,STN)*SC(3,STN) )
         NAr = NAr*ML
        ENDIF
        IF(VERB) WRITE(*,*) STN," max neighbors",NAr*NBMULT
        NBMAX(STNi) = NAr*NBMULT
        NBMAX(STN) = NAr*NBMULT 
        NBMAXLJ(STNi) = NAr*NBMULTLJ
        NBMAXLJ(STN) = NAr*NBMULTLJ 
        NBSTOT  = NBSTOT + NAr*NBMULT 
        NBLJSTOT = NBLJSTOT + NAr*NBMULTLJ 
!       Molecules
!          redundant 1
        MOLCNT(STN) = 1
        MPNT(1,STN) = Nar
        DO I =1,NAr
         MOLST(I,STN) = I
        ENDDO 
        MOLST(NAr+1,STN) = NAr+1
      ENDDO
!     Set for final structure
      NBMAX(STOT) = NBSTOT
      NBMAXLJ(STOT) = NBLJSTOT
!
!     Reset number of atoms in each structure 
      NA(:)     = 0          ! Number of atoms
      NB(:)     = 0          ! Number of bonds
      NANG(:)   = 0          ! Number of angles
      NDIH(:)   = 0          ! Number of dihedrals 
      NIMP(:)   = 0          ! Number of dihedrals 
      ELN(:,:)  = 0          ! Atomic #
      TYP(:,:)  = ' C'       ! Force field type
      RID(:,:)  = 'MOL'      ! Residue ID
      RN(:,:)   = 1          ! Residue #
      MOLN(:,:) = 1          ! Molecule #
      ACHG(:,:) = 0.0d0      ! Charge #
      CHN(:,:)  = 1          ! Charge group #
      GRID(:,:) = ' C'       ! GROMACS ID
      THRM(:,:) = 0          ! Thermostat flag (REBO) 
      R0(:,:,:)   =  0.d0    ! coodinates
      R1(:,:,:)   =  0.d0    ! velocity
      R2(:,:,:)   =  0.d0    ! acceleration 
      FTAG(:,:,:) =  'T'     ! fix tag
      TYPN(:,:) = 0          ! atom type #
!     Potential 
        BNDI(:,:) = 1
        BNDJ(:,:) = 1
        ANGI(:,:) = 1
        ANGJ(:,:) = 1
        ANGK(:,:) = 1
        DIHI(:,:) = 1
        DIHJ(:,:) = 1
        DIHK(:,:) = 1
        DIHL(:,:) = 1
        ! types 
        BTYP(:,:) = 1
        ATYP(:,:) = 1 
        DTYP(:,:) = 1
        ITYP(:,:) = 1
        ! ro
        BVAL(:,:) = -1.0d0 
        AVAL(:,:) = -1.0d0 
        DVAL(:,:) = -1.0d0 
        IVAL(:,:) = -1.0d0 
        ! constant 
        BCNST(:,:) = '0'
        ACNST(:,:) = '0'
        DCNST(:,:) = '0'
        ICNST(:,:) = '0'
   
!     Molecules
!          redundant 1
      MPNT(:,:) = 1
      MOLST(:,:) = 1
      MOLN(:,:) = 1
      MOLCNT(:) = 1
      MCOMS(:,:,:) = 0.0d0 
      MOLMAS(:,:)  = 0.0d0
      MOLDIP(:,:,:) = 0.0d0 
      MOLEF(:,:,:) = 0.0d0 
      MOLEP(:,:) = 0.0d0 
!     Index properties 
      GPROP(:,:,:) = 0.0d0 
!
!      QM/MM
!
      ALLOCATE( QMLIST(NAMX,STOT),MMLIST(NAMX,STOT),FXLIST(NAMX,STOT) )  !  QM/MM lists 
      QMLIST(:,:)  = 0
      MMLIST(:,:)  = 0
      FXLIST(:,:)  = 0
!      
! 
      END SUBROUTINE ALLOCATE
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Deallocate allocateable arrays
!
      SUBROUTINE deall
      USE structure
      USE specify
      USE potential
      USE build 
!
      IMPLICIT none
!
       DEALLOCATE( NA,LV,LC,LA
     & ,VOL,DEN,TMAS
     & ,TYP_IND,TYP_REF,TYP_CNT,REND_IND
     & ,TITLE,INXYZ,INPDB,INTRF,INTOP
     & ,INGRO,INCRD,INCAR,INCOM,THRM,INXMOL,INTBML,RTBML
     & ,RXYZ,RPDB,RTOP,RGRO,RCAR,RCRD,RTRF,RCOM,RANROT
     & ,ELT,ELCNT,NEL,NB,NANG,NDIH,NIMP,ELREF,NMOL
     & ,NGAS,MKRAND,FIXLV,RANBOX,MMN
     & ,EXPRAN,VACEXP,CLLC,CLLV,CLVD,CLMAS,MKCENT,MKVEL
     & ,MKSUP,SC,LCRD,LVRD,MTOAN,DPANGAXFL 
     & ,NBTYP,NATYP,NDTYP,NITYP
     & )

      DEALLOCATE( 
     &  ELN,TYP,RID,RN,MOLN,ACHG,CHN,GRID,R0,R1,R2,FTAG,ELIST
     & ,BTYP,BVAL,BCNST                             ! bond constants
     & ,ATYP,AVAL,ACNST                             ! angle constants
     & ,DTYP,DVAL,DCNST                             ! dihedral constants 
     & ,ITYP,IVAL,ICNST                             ! improper constants 
     & ,ANGI,ANGJ,ANGK,DIHI,DIHJ,DIHK,DIHL,AMAS,NBTOT,NINDX,NLIST
     & ,IMPI,IMPJ,IMPK,IMPL
     & ,MPNT,MOLST,MOLCNT,MOLCH,NBMAX               
     & ,NATPR,RATOMTP,RMATTP                        ! rmatomtype
     & ,CFRAC,CUTFC                                 ! cutregion
     & ,MICENT,MCCENT,MLCENT ,MLMNCENT              ! molcenter
     & ,REFIX,RFIXC,FRTAG                           ! fixregion 
     & ,MCOMS,MOLMAS ,MOLDIP                               ! molecule mass
     & ,CPREG,CPRDS,RCPRG                           ! copyregion 
     & ,ANAX,ANGX,ANGAXS,ANGAXFL,ANCNT              ! angle distribution 
     & ,PRDIS,UDIS,RDIS,PRCNT,PRDISFL,RMIAM         ! atomic rdf 
     & ,MKMRDF,NMOLRi,NMOLRj,NMCNT,MPRDISFL,UMDIS   ! molecular rdf 
     & ,PBCS,UPBC                                   ! PBC's 
     & ,MKSLAB,SCUT                                 ! Surface
     & ,SATOMN,SATOMID,SATOMV                       ! Add surface atoms 
     & ,RENM,RNIDI,RNIDF                            ! Relabel 
     & ,NBCHRG                                      ! Use naborlist for charge groups
     & ,STKI,STKJ,STKBUF,SHFIX                      ! Stack geometeries  
     & ,MKHT                                        ! Hterm 
     & ,MKVAC,SVOL,SDEN,SHT                         ! surface 
     & ,TCHG,DIP,CLDIP,CMDIP      ! dipoles  
     & ,ROTA,RANA,RANB     ! rotate all 
     & ,IPROP,IDPN,IDPID,GPROP   ! index properties 
     & ,SZMAX,SZDIR,SHFT,SHVEC  ! shift 
     & ,MDDIM,MOLDIPFL,CMLDPD,CMDIPFL,MOLDIPAXFL ! molecular dipole distribution 
     & ,RSCALE,DSCALE,DMAX             ! Molecular dipole distribution  
     & ,NBTOTLJ,NINDXLJ,NLISTLJ,NBMAXLJ   ! LJ neighbor list 
     & ,NBTOTML,NINDXML,NLISTML,NBMAXML   ! ML neighbor list 
     & ,WAMONIOM,ONMTAG,WAMONIOMA,WAMONIOMM,WONIOMSP     !     Print oniom inputs   
     & ,ONRAD,ONCENT,USONIOMFX               !  QM/MM lists 
     & ,QMLIST,MMLIST,FXLIST       !  QM/MM lists 
     & ,ONIONQM,ONIONMM,ONIONFX    !  QM/MM lists          
     & ,NONMOL,WONIONMOL           ! QM/MM molecules 
     & ,MEFL,MEFLFL,ATEFL                ! elctric field at molecule     
     & ,RTNK,INTNK,TYPN,RTPA,INTPA
     & ,MOLEF,MOLEP,ATOMEF,ATOMEP ! molecular efield and E pot
     & ,QMNB
     & ,TANAX,TANGX,TANGAXS,TANGAXFL,TANCNT,DPTANGAXFL        ! 3 vec angle distribution 
     & ,GMMIN       !     minimize each molecule with gromacs 
     & ,CHR_UPDATE  !     Update charges 
     & ,MFRAM_XYZ   !     Multi fram analysis
     & ,ATOM_NB     !     Atomic analysis arrays
     & ,CALC_ATOM_ANG,ATOM_ANAX,ATOM_ANGL_CNT,ATOM_ANFILE    ! Atomic angle distribution 
     & ,RDBCNST,RDACNST,RDDCNST,RDICNST                      ! Read in constants from top file 
     & ,UDBOX,NDEP,RSTBOX,DEPBF,MKDTMP,DEPTEMP,ADDDMOL         ! Depostion 
     & ,TNKAMOEBA,TNKOPLS,TNKADD,TNKADDI                      ! TInker types 
     & ,TNKNMD,TNKNUMBM,MTNKNMD,MTNKNUMBM  ! modify tinker numbers 
     & ,MTNKADD,MTNKADDI   ! Tinker types
     & ,CMOL,CMOLN    ! Cut molecule 
     & ,INCHR,RCHR   ! charges 
     & ,MOLCHR
     & ,GBONDS,GANGLES,GDIHS   ! force field info 
     & ,CUTMOLNB   ! cut out 
     &  )
!
      END SUBROUTINE deall
