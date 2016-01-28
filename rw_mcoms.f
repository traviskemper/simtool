C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 10/21/2011 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Read molecule center of mass 

      SUBROUTINE READ_MCOMS(STN)
      USE structure
      USE specify
      USE elements  
!
      IMPLICIT none
!
      INTEGER :: I,NAr,ELr,STN,io
      CHARACTER(ATSZ) :: AT
      REAL*8 :: XX,YY,ZZ
!      
!
      RETURN
 602  FORMAT("Reading xyz file ",A22," structure #:",I8)
      END SUBROUTINE READ_MCOMS

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 10/21/2011 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Write xyz file

      SUBROUTINE WRITE_MCOMS(STN)
      USE specify
      USE structure
      USE elements  
!
      IMPLICIT none
!
      INTEGER :: STN,MNi,MAS
      REAL*8 :: FPE,Ef,Etot,Eo(NDIM),Ei(NDIM),Vo,Vi,Dipl(NDIM)
!
!     Calculate 1/(4 pi e_o) (Nm/C) ^-9
!
      FPE = ELCH/(4.0d0*PI*EPSN )*100
!
!     CHARACTER(ATSZ) :: Atp
!      
!     Verbose output
      IF(VERB) WRITE(6,'("Writing xyz file ",A70)')OUTMCOMS
!     Write file
      OPEN(UNIT=24,FILE=OUTMCOMS,STATUS='unknown')
      WRITE(24,*) MOLCNT(STN) 
      WRITE(24,'(A70)') TITLE(1)
      WRITE(24,2410) 

!     Loop over all molecules
      DO MNi=1,MOLCNT(STN) 
          MAS = INT( MOLMAS(MNi,STN) )
! 
!       Dipole moment
!
        Dipl(:) = MOLDIP(:,MNi,STN)
!
!       Electric field  (N/C) 
!
        Eo(:) = MOLEF(:,MNi,STN)
        Etot = Eo(1)*Eo(1)+Eo(2)*Eo(2)+Eo(3)*Eo(3)
        Ef = SQRT( Etot )
!
!       Electronic potential  (V)
!
         Vi = MOLEP(MNi,STN)
!         Vi = FPE*Vo*10  
!
         WRITE(24,2401) MNi,MAS,MCOMS(:,MNi,STN),Dipl(:),Ei(:),Ef,Vi
      ENDDO
      CLOSE(24)
!
      RETURN
 2410 FORMAT(' # Molmas | COM_x COM_y COM_z | '
     & ,' Mol # | E_x E_y E_z |E| (N/C) | V  (V) ')
 2401 FORMAT(2I6,3F12.3,3F12.3,5E16.4)
      END SUBROUTINE WRITE_MCOMS
