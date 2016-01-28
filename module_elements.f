!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Module of elements  
!
      MODULE elements
      INTEGER :: ATSZ,NELM
      PARAMETER ( ATSZ = 2, NELM =110 )
      CHARACTER(ATSZ), DIMENSION(0:NELM) :: ATSYM 
      REAL*8 ::  AMASS(0:NELM),CRDI(0:NELM),VDWRDI(0:NELM)
!
      END MODULE 
