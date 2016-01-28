C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 10/21/2011 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Read xyz file

      SUBROUTINE READ_XYZ(STN)
      USE structure
      USE specify
      USE elements  
!
      IMPLICIT none
!
      INTEGER :: I,NAr,ELr,STN,io,LTH
      CHARACTER(ATSZ) :: AT
      REAL*8 :: XX,YY,ZZ
      LOGICAL :: ATREAD
!      
!     Verbose output
      IF(VERB) WRITE(6,602) INXYZ(STN),STN
!
      OPEN(UNIT=11,FILE=INXYZ(STN),STATUS='unknown')
      READ(11,*) NAr
      READ(11,'(A100)') TITLE(STN)
!
!     Read in first line to determine if atomic # or atomic symbols are being used
!

      ATREAD = .FALSE.
      READ(11,*) AT
      BACKSPACE(11) 
      DO I = 1, NELM
         IF ( TRIM(AT) .EQ. TRIM(ATSYM(I)) ) ATREAD  = .TRUE.
      ENDDO

!     Debug 
      IF( DEBUG ) THEN
        OPEN(unit=101,file='debug.out',status='unknown')
        WRITE(101,*) 'in rw_xyz'
        WRITE(101,*) ' AT ', AT
        WRITE(101,*) ' ATREAD ',ATREAD 
        CLOSE(101) 
      ENDIF 
! 
!     Read in values 
!
      DO I=NA(STN) + 1,NA(STN) + NAr
        IF( ATREAD) THEN
          READ(11,*) AT,XX,YY,ZZ
          CALL getanumb(AT,ELr)
        ELSE
          READ(11,*) ELr,XX,YY,ZZ
        ENDIF 
        ELN(I,STN) = ELr
        R0(:,I,STN) = (/XX,YY,ZZ/)
      ENDDO
      CLOSE(11)
      NA(STN) = NA(STN) + NAr

!
      RETURN
 602  FORMAT("Reading xyz file ",A22," structure #:",I8)
      END SUBROUTINE READ_XYZ

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 10/21/2011 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Write xyz file

      SUBROUTINE WRITE_XYZ(STN)
      USE specify
      USE structure
      USE elements  
!
      IMPLICIT none
!
      INTEGER :: I,STN,NAr,ELr
      CHARACTER(ATSZ) :: Atp
!      
!     Verbose output
      IF(VERB) WRITE(6,'("Writing xyz file ",A70)')OUTXYZ
!     Write file
      NAr = NA(STN)
      OPEN(UNIT=11,FILE=OUTXYZ,STATUS='unknown')
      WRITE(11,*) NAr
      WRITE(11,'(A100)') TITLE(1)
      DO I=1,NAr
        Elr = ELN(I,STN)
        Atp = ATSYM(ELr)
        WRITE(11,1102) Atp,R0(:,I,STN)
      ENDDO 
      CLOSE(11)
!
      RETURN
 1101 FORMAT(I6,3F12.3)
 1102 FORMAT(A6,"  ",3F24.8)
      END SUBROUTINE WRITE_XYZ
