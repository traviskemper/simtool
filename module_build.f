!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/7/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Module for structure building

      MODULE build
      USE const
!
!     Random molecule builder 
      REAL*8 :: MBUF,CYBUF
      INTEGER :: ATMAX,SEED
      INTEGER, ALLOCATABLE  :: NGAS(:) ,MINRS(:),MINML(:) 
      LOGICAL, ALLOCATABLE :: MKRAND(:),FIXLV(:),RANBOX(:),RANROT(:) 
     & ,EXPRAN(:) ,MINRD(:),MKSUP(:),RANFRX(:,:) 
      REAL*8, ALLOCATABLE :: VACEXP(:,:,:),MINBF(:),SC(:,:),RANFRSZ(:,:)
!     Copy region
      LOGICAL , ALLOCATABLE ::  CPREG(:)
      INTEGER, ALLOCATABLE :: STCPN(:)
      REAL*8, ALLOCATABLE :: CPRDS(:,:),RCPRG(:,:,:) 
!     Surface
      LOGICAL , ALLOCATABLE ::  MKSLAB(:)
      REAL*8,ALLOCATABLE :: SCUT(:) 
!     Add surface atoms 
      INTEGER , ALLOCATABLE :: SATOMN(:)
      CHARACTER(IDSZ), ALLOCATABLE  ::  SATOMID(:,:) 
      REAL*8, ALLOCATABLE :: SATOMV(:,:,:)
!     H term
      LOGICAL , ALLOCATABLE :: MKHT(:) 
!
      END MODULE build

