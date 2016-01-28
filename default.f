!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Default values 

      SUBROUTINE default 
      USE specify
      USE structure
      USE elements  
      USE potential
      USE build 
!
      IMPLICIT none
!
      INTEGER :: ELj,ELi,STN,AN,D
      CHARACTER(CHSZ) :: F,ESTAT

! 
       HELP = .FALSE.
       DEBUG = .FALSE.
!
!     Property calculations 
      CLLC(:) = .FALSE.
      CLLV(:) = .FALSE.
      CLVD(:) = .TRUE.
      CLMAS(:) = .TRUE.
!     Set atomic radi
      CRBUF = 1.25d0         ! covalent bond buffer distance 
      LJBUF = 1.25d0         ! LJ bond buffer distance 
!     Initialize particle properties 
      ATSYM(:) = 'VS'
      CRDI(:) = 1.0d0
      VDWRDI(:) = 1.0d0
      AMASS(:) = 1.0d0
      CALL  CRADI
      CALL  VDWRADI
      CALL  ATOMSYM
      CALL  ATOMMAS
!     Set zero for virtual  sites 
      ATSYM(0)   = 'VS'
      AMASS(0)  = 0.0d0
      VDWRDI(0) = 2.0d0
      CRDI(0) = 1.0d0 
!
      TITLE(:) = "Unknown"
! Input output files
      INXYZ(:)  = 'in.xyz'
      INPDB(:)  = 'in.pdb'
      INCRD(:)  = 'in.crd'
      INTRF(:)  = 'in.trf'
      INTOP(:)  = 'in.top'
      INGRO(:)  = 'in.gro'
      INCAR(:)  = 'in.car'
      INCOM(:)  = 'in.com'
      INTNK(:)  = 'in.tnk'
      INTPA(:)  = 'in.tpa'
      RTNK(:)   = .FALSE.
      RTPA(:)  = .FALSE. 
      OUTXYZ = 'out.xyz' 
      OUTPDB = 'out.pdb'
      OUTCOM = 'out.com'
      OUTGRO = 'out.gro'
      OUTTRF = 'out.trf'
      OUTCAR = 'out.car'
      OUTTOP = 'out.top'
      OUTCRD = 'out.crd'
! Control statements
      VERB = .FALSE. 
      RXYZ(:) = .FALSE.
      RPDB(:) = .FALSE.
      RTOP(:) = .FALSE.
      RGRO(:) = .FALSE.
      RCRD(:) = .FALSE.
      RCAR(:) = .FALSE.
      RTRF(:) = .FALSE.
      RCOM(:) = .FALSE.
      RTNK(:) = .FALSE. 
      STTOG =  .FALSE.
      LCRD(:) = .FALSE.
      LVRD(:) = .FALSE.
      MTOAN(:) = .FALSE.
      RMATTP(:)= .FALSE.       ! rmatomtype
      CFRAC(:) = .FALSE.       ! cutregion
      REFIX(:) = .FALSE.       ! fixregion 
      ANAX(:,:,:) = ""                 ! angle distribution 
      ANGX(:,:,:)  =  0.0d0            ! angle distribution 
      ANGAXS(:)  = .FALSE.           ! angle distribution 
      DO STN=1,STOTo
       DO AN=1,ANCNT(STN)
         IF( STN.GE.10 .OR. AN.GE.10 ) CALL prnerror(-8,ESTAT)
         WRITE(F,'(A6,I1,A1,I1,A4)') 'angdis',STN,'-',AN,'.dat'
         ANGAXFL(AN,STN) =  F     ! angle distribution 
     
       ENDDO
       WRITE(F,'(A6,I1,A4)') 'angdpt',STN,'.dat'
       DPANGAXFL(STN) = F
      ENDDO 
!
      WCRD = .FALSE.
      WXYZ = .FALSE.
      WPDB = .FALSE.
      WTOP = .FALSE.
      WGRO = .FALSE.
      WTRF = .FALSE.
      WCOM = .FALSE.
      WCAR = .FALSE.   
      WTNK = .FALSE. 
      WEMB = .FALSE.

!     Initialize
      NA(:) = 0
      NB(:) = 0
      NANG(:) = 0
      NDIH(:) = 0
      NIMP(:) = 0
      NMOL(:) = 0 
      NBTYP(:) = 0
      NATYP(:) = 0
      NDTYP(:) = 0
      NITYP(:) = 0
!
      VOL(:) = 0.0d0
      DEN(:) = 0.0d0
      TMAS(:) = 0.0d0
      DO STN=1,STOT
        MMN (1,:,STN) = (/1.d16,1.0d16,1.0d16/)
        MMN(2,:,STN) = (/-1.0d16,-1.0d16,-1.0d16/)
      ENDDO
!     Set default box size 
      LV(:,:,:) = 0.0d0
      LV(1,1,:) = 100.0d0
      LV(2,2,:) = 100.0d0
      LV(3,3,:) = 100.0d0
      LC(:,:) = 100.0d0
      LA(:,:) = 90.0d0

      VEL(:,:) = 0.0d0
      MCCENT(:,:) = 0.0d0
      CUTFC(:,:,:) = 0.0d0   ! cutregion
      MICENT(:) = 0   
      NATPR(:)  = 0       ! rmatomtype
      RATOMTP(:,:) = ""       ! rmatomtype 
      RFIXC(:,:,:)  = 0       ! fixregion
      FRTAG(:,:) = "T"          ! fixregion
!         
!     Neighbor list
      NBLST = .TRUE. 
!     Random molecule builder 
      CRAND = .FALSE. 
      SEED = 759847358
      MBUF = 3.0d0
      ATMAX = 50 
      NGAS(:) = 0
      MKRAND(:) = .FALSE.
      RANBOX(:) = .FALSE.
      FIXLV(:) = .FALSE. 
      RANROT(:) = .FALSE. 
      EXPRAN(:) = .FALSE. 
      VACEXP(:,:,:) = 0.0d0
      MINRD(:) = .FALSE.
      MINRS(:) = 0 
      MINBF(:)  = 3.0d0
      CYBUF  = 1.0d0
      RANFRX(:,:) = .FALSE. 
      RANFRSZ(:,:) = 0.0d0
!    Surface
      SDIM = 3  ! z-axis 
      SVOL(:) = 0.0d0 
      SDEN(:) = 0.0d0 
      SHT(:) = 0.0d0 
!     Structure modifications
      MKCENT(:) = .FALSE.
      MKVEL(:) = .FALSE.
      VEL(:,:)  = 0.0d0
!     Z-shift
      SZMAX(:) = .FALSE.
      SZDIR(:) = 3
      SHFT(:)  = .FALSE.
      DO STN=1,STOTo
       DO D = 1,NDIM
        SHVEC(D,STN) = 0.0d0
       ENDDO
      ENDDO
!
!     Suppercell
      SC(:,:) = 1.0d0
      MKSUP(:) = .FALSE.
!     Molecules
      MOLCH(:)  = .FALSE. 
! copyregion 
      CPREG(:) = .FALSE.       
      CPRDS(:,:) = 0.d0          ! copyregion        
      RCPRG(:,:,:) = 0.0d0         ! copyregion 
      STCPN(:) = 1 
!     Molecular rdf 
      NMCNT(:) = 0
      MPRDISFL(:,:) = 'molrdf.dat'
      UMDIS(:,:,:) = .TRUE.  
      MKMRDF(:) = .FALSE.      
      NMOLRi(:,:) = 0       
      NMOLRj(:,:) = 0    
!     Atomic rdf
      PRDISFL(:,:) = 'rdf.dat'
      PRCNT(:) = 0
      PRDIS(:) = .FALSE.       
      UDIS(:,:,:) = .TRUE.  
      RDIS(:,:,:) = ''
      RMIAM(:) = .FALSE.
!     Periodic boundary conditions 
      PBCS(:,:) = .FALSE.
      UPBC(:) = .FALSE. 
!     Surface
      MKSLAB(:) = .FALSE.
      SCUT(:) = 0.1d0      ! Fractional coordinate below which will be translated up to make slab
!     Add surface atoms 
      SATOMID(:,:)  = 'SURF'
      SATOMV(:,:,:)  = 0.0d0
!     Relabel 
      REN(:) = .FALSE. 
      RNIDI(:,:) = 'LAB'
      RNIDF(:,:) = 'LAB'
      RENM(:) = 0
!     Use naborlist for charge groups
      NBCHRG(:)= .FALSE. 
      CHGMX = 30
!     Stack 
      MKSTACK = .FALSE.
      NSTK = 0 
      STKI(:) = 0
      STKJ(:) = 0
      STKBUF(:) = 3.5d0 
      SHFIX(:) = .FALSE. 
!     Hterm
      MKHT(:) = .FALSE. 
!     Molecule center 
      MLCENT(:) = .FALSE. 
      
!     Electronic
      CLDIP(:) = .FALSE.
      CMDIP(:) = .FALSE. 
      TCHG(:)  = 0.0d0 
      DIP(:)   = 0.0d0 
!     Rotate 
      ROTA(:)  = .FALSE.
      RANA(:) = 0.0d0 
      RANB(:) = 0.0d0 
!     Index properties 
      IPROP(:)  = .FALSE.
      IDPN(:) = 0
      IDPID(:,:) = 0 
!     molecular dipole distribution 
      MDDIM(:) = 3
      MOLDIPFL(:) = 'mDIS.dat'
      MOLDIPAXFL(:) = 'mDIS-AX.dat'
      CMLDPD(:) =.FALSE.
      CMDIPFL(:) = 'Mdip.dat'
      RSCALE(:) = 1.d0
      DSCALE(:) =  1.0d0 
      DMAX(:) = 15.0d0 
!
!     Field electric field and potential at each molecule
      MEFL(:) = .FALSE. 
      MEFLFL(:) = 'Mefield.dat'
!     Print oniom inputs   
      WAMONIOM(:)  = .FALSE. 
      WAMONIOMA(:) = .FALSE. 
      WAMONIOMM(:) = .FALSE. 
      WONIOMSP(:)  = .FALSE. 
      ONRAD(:,:)   = 10.0d0 
      ONCENT(:,:)  = 0.0d0 
      ONIONQM(:) = 0  
      ONIONMM(:) = 0   
      ONIONFX(:) = 0   
      USONIOMFX(:) = .FALSE.
      NONMOL(:) = 0 
      WONIONMOL(:,:)  = 0        ! QM/MM molecules 
      QMNB(:) = 0
!     Tri vector angle dis 
      TANAX(:,:,:) = ""                 ! angle distribution 
      TANGX(:,:,:)  =  0.0d0            ! angle distribution 
      TANGAXS(:)  = .FALSE.           ! angle distribution 
      DO STN=1,STOTo
       DO AN=1,TANCNT(STN)
         IF( STN.GE.10 .OR. AN.GE.10 ) CALL prnerror(-8,ESTAT)
         WRITE(F,'(A6,I1,A1,I1,A4)') 'angdis',STN,'-',AN,'.dat'
         TANGAXFL(AN,STN) =  F     ! angle distribution 
     
       ENDDO
       WRITE(F,'(A6,I1,A4)') 'angdpt',STN,'.dat'
       DPTANGAXFL(STN) = F
      ENDDO
!     Molecule min
      GMMIN(:) = .FALSE.  
!     Charge update
      CHR_UPDATE(:) = .FALSE.
!     Atomic angle distribution
      CALC_ATOM_ANG(:)  = .FALSE.      
      DO STN=1,STOTo
       DO AN=1,ATOM_ANGL_CNT(STN)
         IF( STN.GE.10 .OR. AN.GE.10 ) CALL prnerror(-8,ESTAT)
         WRITE(F,'(A6,I1,A1,I1,A4)') 'atom_angdis',STN,'-',AN,'.dat'
         ATOM_ANFILE(AN,STN) =  F     ! angle distribution 
       ENDDO
      ENDDO
      ATOM_ANAX(:,:,:) = ''
!     Xmol read in
      RXMOL(:) = .FALSE.
      MFRAM_XYZ(:) = .FALSE.
!     ff constants read in 
      RDBCNST(:) = .FALSE.
      RDACNST(:) = .FALSE.
      RDDCNST(:) = .FALSE.
      RDICNST(:) = .FALSE. 
!     Read turbomol  
      RTBML(:) = .FALSE. 
!     Depostion
      ADDDM = .FALSE. 
      ADDDMOL(:) = .FALSE. 
      UDBOX(:) = .FALSE. 
      MKDTMP(:) = .FALSE. 
      NDEP(:)   = 0
      RSTBOX(:)  = 1
      DEPBF(:)   = 3.0d0 
      DEPTEMP(:) = 300.0d0
!
!     Tinker types
!
      TNKAMOEBA(:) = .FALSE.  
      TNKOPLS(:) = .FALSE.
      TNKADD(:) = .FALSE. 
      TNKADDI (:) =  200  
      TNKNMD(:) = .FALSE. 
      MTNKNMD(:) = .FALSE. 
      TNKNUMBM(:,:) = 0 
      MTNKNUMBM(:,:) = 0 
!
! 
!     Cut mol
      CMOL = .FALSE.
      CMOLN(:,:) = 0
!
!
      RCHR(:)  = .FALSE.
      INCHR(:) = 'in.chr'
      OUTCHRS  = 'out.chr'
!
!
      ONIOMLAYER = 2
      LAYER_MLT(:) = 1
      LAYER_CHR(:) = 0
      ONIOM_METH ='# oniom(HF/6-31G*:uff)=embedcharge nosymm  '
!     Find information
      GBONDS(:)  = .FALSE.
      GANGLES(:) = .FALSE.
      GDIHS(:)   = .FALSE.
!
!
      RETURN
      END SUBROUTINE default 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Covalent radi from: http://www.ccdc.cam.ac.uk/products/csd/radii/table.php4
!
      SUBROUTINE ATOMSYM
      USE elements 
      ATSYM(0 ) =   'VS' ! 
      ATSYM(1 ) =   'H' ! 
      ATSYM(2 ) =   'He' ! 
      ATSYM(3 ) =   'Li' ! 
      ATSYM(4 ) =   'Be' ! 
      ATSYM(5 ) =   'B' ! 
      ATSYM(6 ) =   'C' ! 
      ATSYM(7 ) =   'N' ! 
      ATSYM(8 ) =   'O' ! 
      ATSYM(9 ) =   'F' ! 
      ATSYM(10 ) =   'Ne' ! 
      ATSYM(11 ) =   'Na' ! 
      ATSYM(12 ) =   'Mg' ! 
      ATSYM(13 ) =   'Al' ! 
      ATSYM(14 ) =   'Si' ! 
      ATSYM(15 ) =   'P' ! 
      ATSYM(16 ) =   'S' ! 
      ATSYM(17 ) =   'Cl' ! 
      ATSYM(18 ) =   'Ar' ! 
      ATSYM(19 ) =   'K' ! 
      ATSYM(20 ) =   'Ca' ! 
      ATSYM(21 ) =   'Sc' ! 
      ATSYM(22 ) =   'Ti' ! 
      ATSYM(23 ) =   'V' ! 
      ATSYM(24 ) =   'Cr' ! 
      ATSYM(25 ) =   'Mn' ! 
      ATSYM(26 ) =   'Fe' ! 
      ATSYM(27 ) =   'Co' ! 
      ATSYM(28 ) =   'Ni' ! 
      ATSYM(29 ) =   'Cu' ! 
      ATSYM(30 ) =   'Zn' ! 
      ATSYM(31 ) =   'Ga' ! 
      ATSYM(32 ) =   'Ge' ! 
      ATSYM(33 ) =   'As' ! 
      ATSYM(34 ) =   'Se' ! 
      ATSYM(35 ) =   'Br' ! 
      ATSYM(36 ) =   'Kr' ! 
      ATSYM(37 ) =   'Rb' ! 
      ATSYM(38 ) =   'Sr' ! 
      ATSYM(39 ) =   'Y' ! 
      ATSYM(40 ) =   'Zr' ! 
      ATSYM(41 ) =   'Nb' ! 
      ATSYM(42 ) =   'Mo' ! 
      ATSYM(43 ) =   'Tc' ! 
      ATSYM(44 ) =   'Ru' ! 
      ATSYM(45 ) =   'Rh' ! 
      ATSYM(46 ) =   'Pd' ! 
      ATSYM(47 ) =   'Ag' ! 
      ATSYM(48 ) =   'Cd' ! 
      ATSYM(49 ) =   'In' ! 
      ATSYM(50 ) =   'Sn' ! 
      ATSYM(51 ) =   'Sb' ! 
      ATSYM(52 ) =   'Te' ! 
      ATSYM(53 ) =   'I' ! 
      ATSYM(54 ) =   'Xe' ! 
      ATSYM(55 ) =   'Cs' ! 
      ATSYM(56 ) =   'Ba' ! 
      ATSYM(71 ) =   'Lu' ! 
      ATSYM(72 ) =   'Hf' ! 
      ATSYM(73 ) =   'Ta' ! 
      ATSYM(74 ) =   'W' ! 
      ATSYM(75 ) =   'Re' ! 
      ATSYM(76 ) =   'Os' ! 
      ATSYM(77 ) =   'Ir' ! 
      ATSYM(78 ) =   'Pt' ! 
      ATSYM(79 ) =   'Au' ! 
      ATSYM(80 ) =   'Hg' ! 
      ATSYM(81 ) =   'Tl' ! 
      ATSYM(82 ) =   'Pb' ! 
      ATSYM(83 ) =   'Bi' ! 
      ATSYM(84 ) =   'Po' ! 
      ATSYM(85 ) =   'At' ! 
      ATSYM(86 ) =   'Rn' ! 
      ATSYM(58 ) =   'Ce' ! 
      ATSYM(66 ) =   'Dy' ! 
      ATSYM(68 ) =   'Er' ! 
      ATSYM(63 ) =   'Eu' ! 
      ATSYM(64 ) =   'Gd' ! 
      ATSYM(67 ) =   'Ho' ! 
      ATSYM(57 ) =   'La' ! 
      ATSYM(60 ) =   'Nd' ! 
      ATSYM(61 ) =   'Pm' ! 
      ATSYM(59 ) =   'Pr' ! 
      ATSYM(62 ) =   'Sm' ! 
      ATSYM(65 ) =   'Tb' ! 
      ATSYM(69 ) =   'Tm' ! 
      ATSYM(70 ) =   'Yb' ! 
      ATSYM(87 ) =   'Fr' ! 
      ATSYM(88 ) =   'Ra' ! 
      ATSYM(104 ) =   'Rf' ! 
      ATSYM(105 ) =   'Db' ! 
      ATSYM(106 ) =   'Sg' ! 
      ATSYM(107 ) =   'Bh' ! 
      ATSYM(108 ) =   'Hs' ! 
      ATSYM(109 ) =   'Mt' ! 
      ATSYM(110 ) =   'Ds' ! 
      ATSYM(89 ) =   'Ac' ! 
      ATSYM(95 ) =   'Am' ! 
      ATSYM(97 ) =   'Bk' ! 
      ATSYM(98 ) =   'Cf' ! 
      ATSYM(96 ) =   'Cm' ! 
      ATSYM(99 ) =   'Es' ! 
      ATSYM(100 ) =   'Fm' ! 
      ATSYM(101 ) =   'Md' ! 
      ATSYM(102 ) =   'No' ! 
      ATSYM(93 ) =   'Np' ! 
      ATSYM(91 ) =   'Pa' ! 
      ATSYM(94 ) =   'Pu' ! 
      ATSYM(90 ) =   'Th' ! 
      ATSYM(92 ) =   ' U' ! 

      END SUBROUTINE ATOMSYM 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      http://www.ccdc.cam.ac.uk/products/csd/radii/table.php4

      SUBROUTINE ATOMMAS 
      USE elements 
       AMASS(1 ) =  1.008d0 !  H
       AMASS(2 ) =  4.003d0 !  He
       AMASS(3 ) =  6.941d0 !  Li
       AMASS(4 ) =  9.012d0 !  Be
       AMASS(5 ) =  10.811d0 !  B
       AMASS(6 ) =  12.011d0 !  C
       AMASS(7 ) =  14.007d0 !  N
       AMASS(8 ) =  15.999d0 !  O
       AMASS(9 ) =  18.998d0 !  F
       AMASS(10 ) =  20.180d0 !  Ne
       AMASS(11 ) =  22.991d0 !  Na
       AMASS(12 ) =  24.305d0 !  Mg
       AMASS(13 ) =  26.982d0 !  Al
       AMASS(14 ) =  28.086d0 !  Si
       AMASS(15 ) =  30.974d0 !  P
       AMASS(16 ) =  32.066d0 !  S
       AMASS(17 ) =  35.453d0 !  Cl
       AMASS(18 ) =  39.948d0 !  Ar
       AMASS(19 ) =  39.098d0 !  K
       AMASS(20 ) =  40.078d0 !  Ca
       AMASS(21 ) =  44.956d0 !  Sc
       AMASS(22 ) =  47.867d0 !  Ti
       AMASS(23 ) =  50.942d0 !  V
       AMASS(24 ) =  51.996d0 !  Cr
       AMASS(25 ) =  54.938d0 !  Mn
       AMASS(26 ) =  55.845d0 !  Fe
       AMASS(27 ) =  58.933d0 !  Co
       AMASS(28 ) =  58.693d0 !  Ni
       AMASS(29 ) =  63.546d0 !  Cu
       AMASS(30 ) =  65.390d0 !  Zn
       AMASS(31 ) =  69.723d0 !  Ga
       AMASS(32 ) =  72.610d0 !  Ge
       AMASS(33 ) =  74.922d0 !  As
       AMASS(34 ) =  78.960d0 !  Se
       AMASS(35 ) =  79.904d0 !  Br
       AMASS(36 ) =  83.800d0 !  Kr
       AMASS(37 ) =  85.468d0 !  Rb
       AMASS(38 ) =  87.620d0 !  Sr
       AMASS(39 ) =  88.906d0 !  Y
       AMASS(40 ) =  91.224d0 !  Zr
       AMASS(41 ) =  92.906d0 !  Nb
       AMASS(42 ) =  95.940d0 !  Mo
       AMASS(43 ) =  98.0d0 !  Tc
       AMASS(44 ) =  101.070d0 !  Ru
       AMASS(45 ) =  102.906d0 !  Rh
       AMASS(46 ) =  106.420d0 !  Pd
       AMASS(47 ) =  107.868d0 !  Ag
       AMASS(48 ) =  112.411d0 !  Cd
       AMASS(49 ) =  114.818d0 !  In
       AMASS(50 ) =  118.71d0 !  Sn
       AMASS(51 ) =  121.760d0 !  Sb
       AMASS(52 ) =  127.600d0 !  Te
       AMASS(53 ) =  126.904d0 !  I
       AMASS(54 ) =  131.290d0 !  Xe
       AMASS(55 ) =  132.905d0 !  Cs
       AMASS(56 ) =  137.327d0 !  Ba
       AMASS(71 ) =  174.967d0 !  Lu
       AMASS(72 ) =  178.490d0 !  Hf
       AMASS(73 ) =  180.948d0 !  Ta
       AMASS(74 ) =  183.840d0 !  W
       AMASS(75 ) =  186.207d0 !  Re
       AMASS(76 ) =  190.230d0 !  Os
       AMASS(77 ) =  192.217d0 !  Ir
       AMASS(78 ) =  195.078d0 !  Pt
       AMASS(79 ) =  196.967d0 !  Au
       AMASS(80 ) =  200.590d0 !  Hg
       AMASS(81 ) =  204.383d0 !  Tl
       AMASS(82 ) =  207.200d0 !  Pb
       AMASS(83 ) =  208.980d0 !  Bi
       AMASS(84 ) =  210.0d0 !  Po
       AMASS(85 ) =  210.0d0 !  At
       AMASS(86 ) =  222.0d0 !  Rn
       AMASS(58 ) =  140.116d0 !  Ce
       AMASS(66 ) =  162.500d0 !  Dy
       AMASS(68 ) =  167.260d0 !  Er
       AMASS(63 ) =  151.964d0 !  Eu
       AMASS(64 ) =  157.250d0 !  Gd
       AMASS(67 ) =  164.930d0 !  Ho
       AMASS(57 ) =  138.906d0 !  La
       AMASS(60 ) =  144.240d0 !  Nd
       AMASS(61 ) =  145.0d0 !  Pm
       AMASS(59 ) =  140.908d0 !  Pr
       AMASS(62 ) =  150.360d0 !  Sm
       AMASS(65 ) =  158.925d0 !  Tb
       AMASS(69 ) =  168.934d0 !  Tm
       AMASS(70 ) =  173.040d0 !  Yb
       AMASS(87 ) =  223.0d0 !  Fr
       AMASS(88 ) =  226.0d0 !  Ra
       AMASS(104 ) =  261.0d0 !  Rf
       AMASS(105 ) =  262.0d0 !  Db
       AMASS(106 ) =  266.0d0 !  Sg
       AMASS(107 ) =  264.0d0 !  Bh
       AMASS(108 ) =  269.0d0 !  Hs
       AMASS(109 ) =  268.0d0 !  Mt
       AMASS(110 ) =  271.0d0 !  Ds
       AMASS(89 ) =  227.0d0 !  Ac
       AMASS(95 ) =  243.0d0 !  Am
       AMASS(97 ) =  247.0d0 !  Bk
       AMASS(98 ) =  251.0d0 !  Cf
       AMASS(96 ) =  247.0d0 !  Cm
       AMASS(99 ) =  252.0d0 !  Es
       AMASS(100 ) =  257.0d0 !  Fm
       AMASS(101 ) =  258.0d0 !  Md
       AMASS(102 ) =  259.0d0 !  No
       AMASS(93 ) =  237.0d0 !  Np
       AMASS(91 ) =  231.036d0 !  Pa
       AMASS(94 ) =  244.0d0 !  Pu
       AMASS(90 ) =  232.038d0 !  Th
       AMASS(92 ) =  238.029d0 !  U
       RETURN
       END SUBROUTINE ATOMMAS 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Covalent radi from: http://www.ccdc.cam.ac.uk/products/csd/radii/table.php4

      SUBROUTINE CRADI
      USE elements 
!
      CRDI(1 ) =  0.23d0 !  H
      CRDI(2 ) =  1.50d0 !  He
      CRDI(3 ) =  1.28d0 !  Li
      CRDI(4 ) =  0.96d0 !  Be
      CRDI(5 ) =  0.83d0 !  B
      CRDI(6 ) =  0.68d0 !  C
      CRDI(7 ) =  0.68d0 !  N
      CRDI(8 ) =  0.68d0 !  O
      CRDI(9 ) =  0.64d0 !  F
      CRDI(10 ) =  1.50d0 !  Ne
      CRDI(11 ) =  1.66d0 !  Na
      CRDI(12 ) =  1.41d0 !  Mg
      CRDI(13 ) =  1.21d0 !  Al
      CRDI(14 ) =  1.20d0 !  Si
      CRDI(15 ) =  1.05d0 !  P
      CRDI(16 ) =  1.02d0 !  S
      CRDI(17 ) =  0.99d0 !  Cl
      CRDI(18 ) =  1.51d0 !  Ar
      CRDI(19 ) =  2.03d0 !  K
      CRDI(20 ) =  1.76d0 !  Ca
      CRDI(21 ) =  1.70d0 !  Sc
      CRDI(22 ) =  1.60d0 !  Ti
      CRDI(23 ) =  1.53d0 !  V
      CRDI(24 ) =  1.39d0 !  Cr
      CRDI(25 ) =  1.61d0 !  Mn
      CRDI(26 ) =  1.52d0 !  Fe
      CRDI(27 ) =  1.26d0 !  Co
      CRDI(28 ) =  1.24d0 !  Ni
      CRDI(29 ) =  1.32d0 !  Cu
      CRDI(30 ) =  1.22d0 !  Zn
      CRDI(31 ) =  1.22d0 !  Ga
      CRDI(32 ) =  1.17d0 !  Ge
      CRDI(33 ) =  1.21d0 !  As
      CRDI(34 ) =  1.22d0 !  Se
      CRDI(35 ) =  1.21d0 !  Br
      CRDI(36 ) =  1.50d0 !  Kr
      CRDI(37 ) =  2.20d0 !  Rb
      CRDI(38 ) =  1.95d0 !  Sr
      CRDI(39 ) =  1.90d0 !  Y
      CRDI(40 ) =  1.75d0 !  Zr
      CRDI(41 ) =  1.64d0 !  Nb
      CRDI(42 ) =  1.54d0 !  Mo
      CRDI(43 ) =  1.47d0 !  Tc
      CRDI(44 ) =  1.46d0 !  Ru
      CRDI(45 ) =  1.45d0 !  Rh
      CRDI(46 ) =  1.39d0 !  Pd
      CRDI(47 ) =  1.45d0 !  Ag
      CRDI(48 ) =  1.44d0 !  Cd
      CRDI(49 ) =  1.42d0 !  In
      CRDI(50 ) =  1.39d0 !  Sn
      CRDI(51 ) =  1.39d0 !  Sb
      CRDI(52 ) =  1.47d0 !  Te
      CRDI(53 ) =  1.40d0 !  I
      CRDI(54 ) =  1.50d0 !  Xe
      CRDI(55 ) =  2.44d0 !  Cs
      CRDI(56 ) =  2.15d0 !  Ba
      CRDI(71 ) =  1.87d0 !  Lu
      CRDI(72 ) =  1.75d0 !  Hf
      CRDI(73 ) =  1.70d0 !  Ta
      CRDI(74 ) =  1.62d0 !  W
      CRDI(75 ) =  1.51d0 !  Re
      CRDI(76 ) =  1.44d0 !  Os
      CRDI(77 ) =  1.41d0 !  Ir
      CRDI(78 ) =  1.36d0 !  Pt
      CRDI(79 ) =  1.50d0 !  Au
      CRDI(80 ) =  1.32d0 !  Hg
      CRDI(81 ) =  1.45d0 !  Tl
      CRDI(82 ) =  1.46d0 !  Pb
      CRDI(83 ) =  1.48d0 !  Bi
      CRDI(84 ) =  1.40d0 !  Po
      CRDI(85 ) =  1.21d0 !  At
      CRDI(86 ) =  1.50d0 !  Rn
      CRDI(58 ) =  2.04d0 !  Ce
      CRDI(66 ) =  1.92d0 !  Dy
      CRDI(68 ) =  1.89d0 !  Er
      CRDI(63 ) =  1.98d0 !  Eu
      CRDI(64 ) =  1.96d0 !  Gd
      CRDI(67 ) =  1.92d0 !  Ho
      CRDI(57 ) =  2.07d0 !  La
      CRDI(60 ) =  2.01d0 !  Nd
      CRDI(61 ) =  1.99d0 !  Pm
      CRDI(59 ) =  2.03d0 !  Pr
      CRDI(62 ) =  1.98d0 !  Sm
      CRDI(65 ) =  1.94d0 !  Tb
      CRDI(69 ) =  1.90d0 !  Tm
      CRDI(70 ) =  1.87d0 !  Yb
      CRDI(87 ) =  2.60d0 !  Fr
      CRDI(88 ) =  2.21d0 !  Ra
      CRDI(89 ) =  2.15d0 !  Ac
      CRDI(95 ) =  1.80d0 !  Am
      CRDI(97 ) =  1.54d0 !  Bk
      CRDI(98 ) =  1.83d0 !  Cf
      CRDI(96 ) =  1.69d0 !  Cm
      CRDI(99 ) =  1.50d0 !  Es
      CRDI(93 ) =  1.90d0 !  Np
      CRDI(91 ) =  2.00d0 !  Pa
      CRDI(94 ) =  1.87d0 !  Pu
      CRDI(90 ) =  2.06d0 !  Th
      CRDI(92 ) =  1.96d0 !  U
      CRDI( 100:NELM  ) =  1.50d0
!
      RETURN
      END SUBROUTINE  CRADI 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Covalent radi from: http://www.ccdc.cam.ac.uk/products/csd/radii/table.php4
!
      SUBROUTINE VDWRADI
      USE elements 
!
      VDWRDI(1 ) =  1.09d0 !  H
      VDWRDI(2 ) =  1.40d0 !  He
      VDWRDI(3 ) =  1.82d0 !  Li
      VDWRDI(4 ) =  2.00d0 !  Be
      VDWRDI(5 ) =  2.00d0 !  B
      VDWRDI(6 ) =  1.70d0 !  C
      VDWRDI(7 ) =  1.55d0 !  N
      VDWRDI(8 ) =  1.52d0 !  O
      VDWRDI(9 ) =  1.47d0 !  F
      VDWRDI(10 ) =  1.54d0 !  Ne
      VDWRDI(11 ) =  2.27d0 !  Na
      VDWRDI(12 ) =  1.73d0 !  Mg
      VDWRDI(13 ) =  2.00d0 !  Al
      VDWRDI(14 ) =  2.10d0 !  Si
      VDWRDI(15 ) =  1.80d0 !  P
      VDWRDI(16 ) =  1.80d0 !  S
      VDWRDI(17 ) =  1.75d0 !  Cl
      VDWRDI(18 ) =  1.88d0 !  Ar
      VDWRDI(19 ) =  2.75d0 !  K
      VDWRDI(20 ) =  2.00d0 !  Ca
      VDWRDI(21 ) =  2.00d0 !  Sc
      VDWRDI(22 ) =  2.00d0 !  Ti
      VDWRDI(23 ) =  2.00d0 !  V
      VDWRDI(24 ) =  2.00d0 !  Cr
      VDWRDI(25 ) =  2.00d0 !  Mn
      VDWRDI(26 ) =  2.00d0 !  Fe
      VDWRDI(27 ) =  2.00d0 !  Co
      VDWRDI(28 ) =  1.63d0 !  Ni
      VDWRDI(29 ) =  1.40d0 !  Cu
      VDWRDI(30 ) =  1.39d0 !  Zn
      VDWRDI(31 ) =  1.87d0 !  Ga
      VDWRDI(32 ) =  2.00d0 !  Ge
      VDWRDI(33 ) =  1.85d0 !  As
      VDWRDI(34 ) =  1.90d0 !  Se
      VDWRDI(35 ) =  1.85d0 !  Br
      VDWRDI(36 ) =  2.02d0 !  Kr
      VDWRDI(37 ) =  2.00d0 !  Rb
      VDWRDI(38 ) =  2.00d0 !  Sr
      VDWRDI(39 ) =  2.00d0 !  Y
      VDWRDI(40 ) =  2.00d0 !  Zr
      VDWRDI(41 ) =  2.00d0 !  Nb
      VDWRDI(42 ) =  2.00d0 !  Mo
      VDWRDI(43 ) =  2.00d0 !  Tc
      VDWRDI(44 ) =  2.00d0 !  Ru
      VDWRDI(45 ) =  2.00d0 !  Rh
      VDWRDI(46 ) =  1.63d0 !  Pd
      VDWRDI(47 ) =  1.72d0 !  Ag
      VDWRDI(48 ) =  1.58d0 !  Cd
      VDWRDI(49 ) =  1.93d0 !  In
      VDWRDI(50 ) =  2.17d0 !  Sn
      VDWRDI(51 ) =  2.00d0 !  Sb
      VDWRDI(52 ) =  2.06d0 !  Te
      VDWRDI(53 ) =  1.98d0 !  I
      VDWRDI(54 ) =  2.16d0 !  Xe
      VDWRDI(55 ) =  2.00d0 !  Cs
      VDWRDI(56 ) =  2.00d0 !  Ba
      VDWRDI(71 ) =  2.00d0 !  Lu
      VDWRDI(72 ) =  2.00d0 !  Hf
      VDWRDI(73 ) =  2.00d0 !  Ta
      VDWRDI(74 ) =  2.00d0 !  W
      VDWRDI(75 ) =  2.00d0 !  Re
      VDWRDI(76 ) =  2.00d0 !  Os
      VDWRDI(77 ) =  2.00d0 !  Ir
      VDWRDI(78 ) =  1.72d0 !  Pt
      VDWRDI(79 ) =  1.66d0 !  Au
      VDWRDI(80 ) =  1.55d0 !  Hg
      VDWRDI(81 ) =  1.96d0 !  Tl
      VDWRDI(82 ) =  2.02d0 !  Pb
      VDWRDI(83 ) =  2.00d0 !  Bi
      VDWRDI(84 ) =  2.00d0 !  Po
      VDWRDI(85 ) =  2.00d0 !  At
      VDWRDI(86 ) =  2.00d0 !  Rn
      VDWRDI(58 ) =  2.00d0 !  Ce
      VDWRDI(66 ) =  2.00d0 !  Dy
      VDWRDI(68 ) =  2.00d0 !  Er
      VDWRDI(63 ) =  2.00d0 !  Eu
      VDWRDI(64 ) =  2.00d0 !  Gd
      VDWRDI(67 ) =  2.00d0 !  Ho
      VDWRDI(57 ) =  2.00d0 !  La
      VDWRDI(60 ) =  2.00d0 !  Nd
      VDWRDI(61 ) =  2.00d0 !  Pm
      VDWRDI(59 ) =  2.00d0 !  Pr
      VDWRDI(62 ) =  2.00d0 !  Sm
      VDWRDI(65 ) =  2.00d0 !  Tb
      VDWRDI(69 ) =  2.00d0 !  Tm
      VDWRDI(70 ) =  2.00d0 !  Yb
      VDWRDI(87 ) =  2.00d0 !  Fr
      VDWRDI(88 ) =  2.00d0 !  Ra
      VDWRDI(89 ) =  2.00d0 !  Ac
      VDWRDI(95 ) =  2.00d0 !  Am
      VDWRDI(97 ) =  2.00d0 !  Bk
      VDWRDI(98 ) =  2.00d0 !  Cf
      VDWRDI(96 ) =  2.00d0 !  Cm
      VDWRDI(99 ) =  2.00d0 !  Es
      VDWRDI(93 ) =  2.00d0 !  Np
      VDWRDI(91 ) =  2.00d0 !  Pa
      VDWRDI(94 ) =  2.00d0 !  Pu
      VDWRDI(90 ) =  2.00d0 !  Th
      VDWRDI(92 ) =  1.86d0 !  U
      VDWRDI(100:NELM ) =  2.00d0
!
      RETURN
      END SUBROUTINE  VDWRADI 
