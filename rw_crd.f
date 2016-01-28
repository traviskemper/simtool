! Read in REBO structure file coord.d
!
      SUBROUTINE read_crd(STN)
      USE specify
      USE structure 
!
      IMPLICIT none
      INTEGER :: STN,I,K,T,ELr,NAr,Kr,MNr
      REAL*8 :: Xr,Yr,Zr
!
      OPEN(UNIT=18,FILE=INCRD(STN),STATUS='unknown')

      IF( VERB ) THEN
        WRITE(*,*) 'Reading in  REBO input file :',INCRD(STN),STN
      ENDIF
      READ(18,'(A100)') TITLE(STN)
      READ(18,*) NAr
      READ(18,*)
      READ(18,*)
      READ(18,*) LC(:,STN)
      LA(:,STN) = (/ 90.0,90.0,90.0/)
      DO I=NA(STN) + 1,NA(STN) + NAr
        READ(18,*) K,ELr,Xr,Yr,Zr,T
        ELN(I,STN) = ELr
        R0(:,K,STN) = (/Xr,Yr,Zr /) 
        READ(18,*) Kr,R1(:,I,STN)
        READ(18,*) Kr,R2(:,I,STN)
        READ(18,*) 
      ENDDO
      CLOSE(18)
      NA(STN) = NA(STN) + NAr
!
      END SUBROUTINE READ_CRD
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Write REBO structure file coord.d
!
      SUBROUTINE write_crd(STN)
      USE specify
      USE structure 
!
      IMPLICIT none
      INTEGER :: STN,I,K,T,NAr
      REAL*8 :: R3(NDIM)
!
!     Set locals
      NAr = NA(STN)
      R3 = (/ 0.d0,0.d0,0.d0 /)
      OPEN(UNIT=59,FILE=OUTCRD,STATUS='unknown')
      IF( VERB ) THEN
        WRITE(*,*) 'Writing coord.d',OUTCRD,STN
      ENDIF
      WRITE(59,'(A100)') TITLE(STN)
      WRITE(59,*) NAr
      WRITE(59,*)
      WRITE(59,*)
      WRITE(59,*) LC(:,STN)
      DO I =1,NAr       
        WRITE(59,350) I,ELN(I,STN),R0(:,I,STN),THRM(I,STN)
        WRITE(59,360) I,R1(1:3,I,STN)
        WRITE(59,360) I,R2(1:3,I,STN)
        WRITE(59,360) I,R3
      ENDDO
      CLOSE(59)
  350 FORMAT(2I7,3E20.11,I3)
  360 FORMAT(I7,3E20.11)
!
      END SUBROUTINE WRITE_CRD
