!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/02/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Potential information
!
      MODULE potential
      USE const
      INTEGER,ALLOCATABLE :: NB(:)
     & ,BNDI(:,:),BNDJ(:,:),NBTYP(:)
     & ,BTYP(:,:),ATYP(:,:),DTYP(:,:),ITYP(:,:)
     & ,NANG(:),ANGI(:,:),ANGJ(:,:),ANGK(:,:),NATYP(:)
     & ,NDIH(:),DIHI(:,:),DIHJ(:,:),DIHK(:,:),DIHL(:,:),NDTYP(:)
     & ,NIMP(:),IMPI(:,:),IMPJ(:,:),IMPK(:,:),IMPL(:,:),NITYP(:)
     & ,NMOL(:) 
      REAL*8,ALLOCATABLE  ::  AMAS(:,:) 
     & ,BVAL(:,:),AVAL(:,:),DVAL(:,:),IVAL(:,:) 
      CHARACTER(CHSZ), ALLOCATABLE ::  
     & BCNST(:,:),ACNST(:,:),DCNST(:,:),ICNST(:,:)
      LOGICAL,ALLOCATABLE :: RDBCNST(:),RDACNST(:),RDDCNST(:),RDICNST(:)

!
      END MODULE potential
 
