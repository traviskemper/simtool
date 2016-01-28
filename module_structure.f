!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Module of structure information 
!
      MODULE structure
      USE const
!

!     Variables dependent on the number of structures 
      INTEGER :: STOT,STOTo,STADD,NAMX
      INTEGER, ALLOCATABLE :: NA(:)
      REAL*8, ALLOCATABLE ::  LC(:,:),LA(:,:),LV(:,:,:),VOL(:),DEN(:)
     & ,TMAS(:)
!     Variables dependent on the number of atoms
      INTEGER, ALLOCATABLE :: ELN(:,:),RN(:,:),MOLN(:,:),THRM(:,:)
     & ,CHN(:,:),TYP_CNT(:,:),TYP_REF(:,:)
      REAL*8, ALLOCATABLE :: R0(:,:,:),ACHG(:,:)
     & ,R1(:,:,:),R2(:,:,:)
      CHARACTER(IDSZ), ALLOCATABLE ::  RID(:,:)
     &  ,FTAG(:,:,:),GRID(:,:)
      CHARACTER(3), ALLOCATABLE ::  TYP(:,:),TYP_IND(:,:)
      INTEGER, ALLOCATABLE :: REND_IND(:,:)
!
!     Element list
      INTEGER, ALLOCATABLE :: ELIST(:,:,:),NEL(:,:)
     & ,ELT(:),ELCNT(:,:),ELREF(:,:)
!     Periodic boundary conditions 
      LOGICAL, ALLOCATABLE  :: PBCS(:,:) 
!     Fractional coordinates
      REAL*8, ALLOCATABLE :: FR(:,:),MLVr(:)
!     Neighbor list
      LOGICAL ::  NBLST
      REAL*8 :: CRBUF 
      INTEGER, ALLOCATABLE ::  NBTOT(:),NINDX(:,:),NLIST(:,:),NBMAX(:) 
!     Center of mass
      REAL*8 :: COMAS(NDIM) 
!     Max.min
      REAL*8, ALLOCATABLE :: MMN(:,:,:)
!     Surface 
      INTEGER :: SDIM 
      REAL*8, ALLOCATABLE :: SVOL(:),SDEN(:),SHT(:)
!     Molecule list
      INTEGER, ALLOCATABLE :: MPNT(:,:),MOLST(:,:),MOLCNT(:)
      REAL*8, ALLOCATABLE :: MCOMS(:,:,:),MOLMAS(:,:),MOLDIP(:,:,:)
     &  ,MOLEF(:,:,:),MOLEP(:,:),ATOMEF(:,:,:),ATOMEP(:,:),MOLCHR(:,:)
      LOGICAL, ALLOCATABLE ::  MOLCH(:) 
!     Electronic properties      
      LOGICAL, ALLOCATABLE ::  CLDIP(:),CMDIP(:)
      REAL*8, ALLOCATABLE ::  TCHG(:),DIP(:)
!     LJ neighbor list
      REAL*8 :: LJBUF 
      INTEGER, ALLOCATABLE ::  NBTOTLJ(:),NINDXLJ(:,:)
     & ,NLISTLJ(:,:),NBMAXLJ(:) 
!     ML neighbor list
      
      INTEGER, ALLOCATABLE ::  NBTOTML(:),NINDXML(:,:)
     & ,NLISTML(:,:),NBMAXML(:) 
!     Tinker structure files
      LOGICAL  :: WTNK  
      LOGICAL, ALLOCATABLE :: RTNK(:),RTPA(:) 
      CHARACTER(CHSZ) ::  OUTTNK
      CHARACTER(CHSZ), ALLOCATABLE  ::  INTNK(:),INTPA(:) 
      INTEGER, ALLOCATABLE :: TYPN(:,:)
!     QM/MM
      INTEGER :: LAYER_MLT(3),LAYER_CHR(3),ONIOMLAYER
      CHARACTER(IDSZ), ALLOCATABLE ::  ONMTAG(:,:)
      CHARACTER(100) ::  ONIOM_METH
!

!
      END MODULE structure 
