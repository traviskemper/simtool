!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 10/24/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Calculate maximum and minimum

      SUBROUTINE maxmin(STN)
      USE structure
      USE specify 
!
      IMPLICIT none
!
      INTEGER :: I,EL,D,STN
!      
      MMN(1,:,STN) = (/1.d16,1.0d16,1.0d16/)
      MMN(2,:,STN) = (/-1.0d16,-1.0d16,-1.0d16/)
!     
      DO I=1,NA(STN)
        DO D=1,NDIM
          IF ( R0(D,I,STN).LT.MMN(1,D,STN)) MMN(1,D,STN)=R0(D,I,STN)
          IF ( R0(D,I,STN).GT.MMN(2,D,STN)) MMN(2,D,STN)=R0(D,I,STN)
        ENDDO
      ENDDO 
!
      END SUBROUTINE maxmin
