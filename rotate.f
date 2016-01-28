!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 10/24/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Center of mass
!
      SUBROUTINE centermass(STN)
!
      USE structure
      USE elements 
      USE potential 
      USE specify 
!
      IMPLICIT none 
!
      INTEGER :: I,EL,D,STN
      REAL*8  :: MS,RMAS(NDIM),TOTMr
!
      TOTMr = 0.0d0
      RMAS(:) = 0.0d0
!
      DO I=1,NA(STN)
        MS = AMAS(I,STN ) 
        DO D =1,3        
          RMAS(D) = RMAS(D) + MS*R0(D,I,STN)
        ENDDO
        TOTMr = TOTMr + MS 
      ENDDO
!
      COMAS(1) = RMAS(1)/TOTMr
      COMAS(2) = RMAS(2)/TOTMr
      COMAS(3) = RMAS(3)/TOTMr
!
!      IF( VERB ) THEN 
!        WRITE(*,*) " Center of mass of ",STN," :"
!        WRITE(*,*) COMAS(:) 
!      ENDIF 
      END SUBROUTINE centermass
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 10/24/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Randomize orientation 
!
      SUBROUTINE rotate(STN,ANa,ANb)
!
      USE structure
      USE specify
      USE const
      USE build 
!
      IMPLICIT none 
!
      INTEGER :: I,STN
      REAL*8  :: DP,DT,CY,SY,CZ,SZ,ANa,ANb
     &           ,xd,yd,zd
!
      CALL centermass(STN)  

!
      dp = ANa 
      dt = Anb 
!
      cy = cos(dp)
      sy = sin(dp)
      cz = cos(dt)
      sz = sin(dt)
      DO I=1,NA(STN)
        xd = R0(1,I,STN) - COMAS(1)
        yd = R0(2,I,STN) - COMAS(2)
        zd = R0(3,I,STN) - COMAS(3)
        R0(1,I,STN) =  cy*cz*xd - sz*cy*yd + sy*zd + COMAS(1)
        R0(2,I,STN) =  sz*xd    + cz*yd            + COMAS(2)
        R0(3,I,STN) = -sy*cz*xd + sy*sz*yd + cy*zd + COMAS(3)
      ENDDO
!      IF( VERB ) THEN
!        WRITE(*,*) " Rotating ",STN,dp,dt
!      ENDIF

!
      END SUBROUTINE rotate
