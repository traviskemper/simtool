C     Version 1.0 10/21/2011 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Read xyz file

      SUBROUTINE READ_XMOL(STN)
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
!     Verbose output
      IF(VERB) WRITE(6,602) INXMOL(STN),STN
!
      OPEN(UNIT=19,FILE=INXMOL(STN),STATUS='unknown')
      READ(19,*) NAr
      READ(19,'(A100)') TITLE(STN)
      DO I=NA(STN) + 1,NA(STN) + NAr
         ELr = -1
        READ(19,'(A,3F8.2)',iostat=io) AT,XX,YY,ZZ
        CALL getanumb(AT,ELr)
        IF( io.NE.0 ) THEN 
          BACKSPACE(19) 
          READ(19,*) ELr,XX,YY,ZZ
!          READ(11,*,iostat=io) ELr,XX,YY,ZZ
 !        IF(io.EQ.0) CALL prnerror(79)
        ENDIF      
        ELN(I,STN) = ELr
        R0(:,I,STN) = (/XX,YY,ZZ/)
      ENDDO
      CLOSE(19)
      NA(STN) = NA(STN) + NAr
!
      RETURN
 602  FORMAT("Reading xmol file ",A22," structure #:",I8)
      END SUBROUTINE READ_XMOL
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Calculate the angle between 3 atoms 

      SUBROUTINE MF_READ_XMOL(STN)
!
      USE specify
      USE structure
      USE build 
!
      Implicit none
!
      INTEGER :: STN,I,LSTEP,LCHK,NAi
     & ,ios
      REAL*8 :: LVr(NDIM,NDIM),TTIME,DELTA,Ri(NDIM)
      CHARACTER(IDSZ) :: NATOM 
      CHARACTER(CHSZ) :: ESTAT
!
!     Set/initialize  local variables    
!
      LSTEP=0
      LCHK = 1
      TTIME = 0
      DELTA = 1
!
!     Verbose output 
!
      IF(VERB) WRITE(6,*) 'Starting xmol analysis'
 !
!     Triclinc PBC's
      LVr(:,:) = LV(:,:,STN)
      NAi = NA(STN)
  !!!    ALLOCATE( Fr(NDIM,NAi) ) 
     
!
!     Open xmol file 
!
      OPEN(UNIT=19,FILE=INXMOL(STN),STATUS='unknown')
!
      DO
        READ(19,*,IOSTAT=ios) NAi
        IF ( IOS .GT. 0 ) THEN
          WRITE(ESTAT,*)' Atoms',NAi,LSTEP
          CALL PRNERROR(90,ESTAT)
        ELSE IF (ios .LT. 0 ) THEN
          EXIT
        ELSE
          LSTEP = LSTEP + 1
          TTIME=TTIME+DELTA
!
!         Read positions from xmol file
!
          READ(19,*) 
          DO I=1,NAi
            READ(19,*) NATOM,Ri(:)
            R0(:,I,STN) = Ri(:)
          ENDDO
        ENDIF ! IOS
!
!       Calculate fractional coordinates 
!
  !!!     CALL fracR(STN,LVr,1)        
        CALL BUILDNBL(STN)
!
!       Run analysis 
!       
        IF( CALC_ATOM_ANG(STN) ) CALL MF_CALC_ANGLE(STN)
      ENDDO 
!
    !!!  DEALLOCATE( Fr ) 
      CLOSE(19)
      RETURN 
      END SUBROUTINE MF_READ_XMOL
