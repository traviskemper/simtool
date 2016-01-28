!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Module of specifications  
!
      MODULE specify
      USE const
!
!     Print help info
      LOGICAL :: HELP,DEBUG
!     Run modifications and output geometry 
      LOGICAL :: WRTOUT 
!     Zero cut off
      REAL*8 :: ZCUT 
      PARAMETER( ZCUT = 0.00000001d0 ) 
!     Input names
      CHARACTER(CHSZ), ALLOCATABLE, DIMENSION(:) :: 
     & TITLE,INXYZ,INPDB,INTRF,INTOP,INGRO,INCRD,INCAR,INCOM,INLMP 
     & ,INXMOL,INTBML,INCHR
!      Output name
      CHARACTER(CHSZ) :: OUTXYZ,OUTPDB,OUTTOP,OUTCOM,OUTGRO,OUTCRD
     &   ,OUTTRF,OUTCAR,OUTMCOMS,OUTLMP,OUTTBMOL,OUTXMOL,OUTCHRS
     & ,OUTEMTBMOL 
!      Read in control variables  
       LOGICAL, ALLOCATABLE, DIMENSION(:) :: 
     &  RXYZ,RPDB,RTOP,RGRO,RCAR,RCRD,RTRF,RCOM,LCRD,LVRD,RLMP,RXMOL
     & ,RTBML,RCHR
!      Output control variables 
       LOGICAL :: VERB,WXYZ,WPDB,WTOP,WCOM,WGRO,WCRD,WCAR,WTRF,CRAND
     & ,WMCOMS,WLMP,WTBML,WEMB
!      Atom property controls
       LOGICAL :: STTOG 
!      Property calculations
       LOGICAL, ALLOCATABLE  :: CLLC(:),CLLV(:),CLVD(:),CLMAS(:)
     & ,MKMRDF(:),MKRDF(:)
!      Modify structure controls
       LOGICAL, ALLOCATABLE  :: MKCENT(:),MKVEL(:) ,MKVAC(:),MTOAN(:)
     & ,MLCENT(:),RMATTP(:),CFRAC(:),REFIX(:),MLMNCENT(:)
      INTEGER, ALLOCATABLE  :: MICENT(:),NATPR(:) 
      CHARACTER(IDSZ), ALLOCATABLE  :: RATOMTP(:,:),FRTAG(:,:)
      REAL*8 :: FCENT(NDIM),VACL
      REAL*8, ALLOCATABLE :: VEL(:,:),MCCENT(:,:),CUTFC(:,:,:)
     &  ,RFIXC(:,:,:)
!     molecular rdf
      INTEGER, ALLOCATABLE  :: NMCNT(:) 
      CHARACTER(CHSZ), ALLOCATABLE  :: MPRDISFL(:,:) 
      LOGICAL, ALLOCATABLE :: UMDIS(:,:,:)
       INTEGER,  ALLOCATABLE  :: NMOLRi(:,:),NMOLRj(:,:)

!     Angle  distribution    
      CHARACTER(CHSZ), ALLOCATABLE :: ANGAXFL(:,:),DPANGAXFL(:)
      INTEGER, ALLOCATABLE :: ANCNT(:) 
      CHARACTER(CHSZ), ALLOCATABLE :: ANAX(:,:,:)
      LOGICAL, ALLOCATABLE :: ANGAXS(:) 
      REAL*8,  ALLOCATABLE :: ANGX(:,:,:)
!     Average vector angle  distribution    
      CHARACTER(CHSZ), ALLOCATABLE :: TANGAXFL(:,:),DPTANGAXFL(:)
      INTEGER, ALLOCATABLE :: TANCNT(:) 
      CHARACTER(CHSZ), ALLOCATABLE :: TANAX(:,:,:)
      LOGICAL, ALLOCATABLE :: TANGAXS(:) 
      REAL*8,  ALLOCATABLE :: TANGX(:,:,:)
!     Radial distrobution function 
      INTEGER, ALLOCATABLE ::  PRCNT(:)  
      CHARACTER(CHSZ), ALLOCATABLE :: PRDISFL(:,:) 
      LOGICAL, ALLOCATABLE :: PRDIS(:),UDIS(:,:,:),RMIAM(:)
      CHARACTER(IDSZ), ALLOCATABLE  :: RDIS(:,:,:)
!     PBC's
      LOGICAL, ALLOCATABLE :: UPBC(:) 
!     Relabel 
      LOGICAL, ALLOCATABLE :: REN(:)
      INTEGER, ALLOCATABLE :: RENM(:) 
      CHARACTER(IDSZ), ALLOCATABLE :: RNIDI(:,:),RNIDF(:,:)
!     Use naborlist for charge groups
      LOGICAL, ALLOCATABLE ::  NBCHRG(:)
      INTEGER :: CHGMX
!     Stack geometeries 
      LOGICAL :: MKSTACK 
      INTEGER :: NSTK
      INTEGER, ALLOCATABLE :: STKI(:),STKJ(:)
      REAL*8,  ALLOCATABLE :: STKBUF(:)
      LOGICAL, ALLOCATABLE :: SHFIX(:) 
!     Rotate 
      LOGICAL, ALLOCATABLE :: ROTA(:) 
      REAL*8,  ALLOCATABLE :: RANA(:),RANB(:) 
!     Property of index 
      LOGICAL, ALLOCATABLE :: IPROP(:) 
      INTEGER, ALLOCATABLE :: IDPN(:) ,IDPID(:,:)
      REAL*8,  ALLOCATABLE :: GPROP(:,:,:)
!     Z Shift 
      LOGICAL, ALLOCATABLE ::  SZMAX(:)
      INTEGER, ALLOCATABLE :: SZDIR(:)
!     Shift 
      LOGICAL, ALLOCATABLE ::  SHFT(:)
      REAL*8 , ALLOCATABLE ::  SHVEC(:,:)
!     Molecular dipole distribution 
      LOGICAL, ALLOCATABLE :: CMLDPD(:) 
      INTEGER, ALLOCATABLE :: MDDIM(:)
      CHARACTER(CHSZ), ALLOCATABLE ::  MOLDIPFL(:),CMDIPFL(:)
     & ,MOLDIPAXFL(:)
      REAL*8, ALLOCATABLE :: RSCALE(:),DSCALE(:),DMAX(:) 
!     Print QM/MM inputs
      LOGICAL, ALLOCATABLE :: WAMONIOM(:)
     & ,WAMONIOMA(:),WAMONIOMM(:),WONIOMSP(:),USONIOMFX(:)
      REAL*8, ALLOCATABLE  :: ONRAD(:,:),ONCENT(:,:)
      INTEGER, ALLOCATABLE  ::  QMLIST(:,:),MMLIST(:,:),FXLIST(:,:) 
     & ,ONIONQM(:),ONIONMM(:),ONIONFX(:) 
     & ,NONMOL(:),WONIONMOL(:,:),QMNB(:) 
!     Field electric field and potential at each molecule
      LOGICAL, ALLOCATABLE :: MEFL(:),ATEFL(:) 
      CHARACTER(CHSZ), ALLOCATABLE :: MEFLFL(:)
!     Minimize each molecule 
      LOGICAL, ALLOCATABLE :: GMMIN(:)
      LOGICAL :: MSLET
      INTEGER :: MANUMB 
!     Update charges 
      LOGICAL, ALLOCATABLE :: CHR_UPDATE(:)
!     Multi fram analysis
      LOGICAL, ALLOCATABLE ::  MFRAM_XYZ(:)
!     Atomic angle distribution 
      CHARACTER(CHSZ), ALLOCATABLE :: ATOM_ANFILE(:,:)
      INTEGER, ALLOCATABLE :: ATOM_ANGL_CNT(:), ATOM_THDIS(:,:)
      REAL*8, ALLOCATABLE :: ATOM_THSUM(:)
      CHARACTER(CHSZ), ALLOCATABLE :: ATOM_ANAX(:,:,:)
      LOGICAL, ALLOCATABLE :: CALC_ATOM_ANG(:) 
!     Atomic analysis arrays
      INTEGER, ALLOCATABLE :: ATOM_NB(:,:)
     & ,BCNT(:), ATOMB_NB(:,:) 
!     Depostion of molecules
      LOGICAL, ALLOCATABLE :: ADDDMOL(:),UDBOX(:),MKDTMP(:)
      LOGICAL :: ADDDM
      INTEGER, ALLOCATABLE :: NDEP(:),RSTBOX(:)
      REAL*8, ALLOCATABLE :: DEPBF(:),DEPTEMP(:)
!     Tinker style
      LOGICAL, ALLOCATABLE :: TNKAMOEBA(:),TNKOPLS(:),TNKADD(:)
     & ,TNKNMD(:),MTNKNMD(:),MTNKADD(:)
      INTEGER, ALLOCATABLE :: TNKADDI(:),TNKNUMBM(:,:),MTNKNUMBM(:,:)
     & ,MTNKADDI(:)
!
!     Cut out moleucle
      LOGICAL, ALLOCATABLE :: CMOL(:)
      INTEGER, ALLOCATABLE :: CMOLN(:,:)
!     Get bonds
      LOGICAL, ALLOCATABLE :: GBONDS(:),GANGLES(:),GDIHS(:)
!     Cut out neighbors
      LOGICAL, ALLOCATABLE :: CUTMOLNB(:)

      END MODULE specify 
