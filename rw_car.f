
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 02/06/2011 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Read  POSCAR VASP file

      SUBROUTINE READ_CAR(STN)
      USE structure
      USE specify
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: N,ELNr,STN,NELr,io,E,I
      INTEGER, ALLOCATABLE :: ELISTL(:)
      CHARACTER(ATSZ), ALLOCATABLE :: ATLISTL(:)
      CHARACTER(ATSZ) :: SYMr
      REAL*8 :: XX,YY,ZZ,LMULT,LVL(NDIM,NDIM),LVr(NDIM,NDIM)
!     Verbose output
      IF(VERB) WRITE(6,602) INCAR(STN),STN
!
      OPEN(UNIT=16,FILE=INCAR(STN),STATUS='unknown')
      READ(16,'(A)')  TITLE(STN)
      READ(16,*)  LMULT
      READ(16,*)  LVr(1,1:3)
      READ(16,*)  LVr(2,1:3)
      READ(16,*)  LVr(3,1:3)
      READ(16,*)  
      DO E=NELM,1,-1
        ALLOCATE( ELISTL(E) )
        READ(16,*,iostat=io) ELISTL(:)
        IF( io.EQ.0) THEN
          NELr = E
          EXIT
        ENDIF
        BACKSPACE(16)
        BACKSPACE(16)
        DEALLOCATE( ELISTL  )
      ENDDO
!     Go back and read atom types
      BACKSPACE(16)
      BACKSPACE(16)
      ALLOCATE( ATLISTL(NELr) )
      READ(16,*)ATLISTL(:)
      READ(16,*)
      IF(VERB) THEN
       WRITE(6,*) " # of elements in ",STN," is ",NELr
       WRITE(6,*) ATLISTL(:)
      ENDIF
!
!     Calculate lattice vectors
      LVL(:,:) = LVr(:,:)*LMULT
!     Skip over specifications
      READ(16,*,iostat=io) XX
      DO while (io.NE.0)
       READ(16,*,iostat=io) XX
      ENDDO
      BACKSPACE(16)
!
!     Loop over fractional coordinates 
      I = NA(STN)
      DO E=1,NELr
        SYMr = ADJUSTR( ATLISTL(E) )
        CALL GETANUMB(SYMr,ELNr)   
        IF(VERB) WRITE(6,*) " Element ",E," #",ELNr," is ",SYMr
        DO N=1,ELISTL(E)
          I = I + 1
          ELN(I,STN) = ELNr
          READ(16,*) XX,YY,ZZ,FTAG(1:NDIM,I,STN)
          R0(1,I,STN) = LVL(1,1)*XX+LVL(1,2)*YY+LVL(1,3)*ZZ
          R0(2,I,STN) = LVL(2,1)*XX+LVL(2,2)*YY+LVL(2,3)*ZZ
          R0(3,I,STN) = LVL(3,1)*XX+LVL(3,2)*YY+LVL(3,3)*ZZ
        ENDDO
      ENDDO
      CLOSE(16) 
!     Save lattice vectors
      LV(:,:,STN)  = LVL(:,:)
!     Save total number of atoms   
      NA(STN) = I
!
      RETURN
 602  FORMAT("Reading CAR file ",A22," structure #:",I8)
      END SUBROUTINE READ_CAR
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
!     Version 1.0 02/06/2011 T. W. Kemper                      C
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
!  Write  POSCAR VASP file

      SUBROUTINE WRITE_CAR(STN)
      USE structure
      USE specify
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: N,EL,STN,NELr,io,E,I,NAr,J
      INTEGER, ALLOCATABLE :: ELISTL(:)
      REAL*8 :: LMULT,LVr(NDIM,NDIM)
      CHARACTER(500) :: ELSTRNG
      CHARACTER(ATSZ), ALLOCATABLE, DIMENSION(:)  :: ATSTRNG
!
!     Verbose output
      IF(VERB) WRITE(6,'("Writing car file ",A70)')OUTCAR
!     Set local variables for output
      NAr = NA(STN) 
      LMULT = 1.0d0
      NELr = ELT(STN)
      ALLOCATE( ATSTRNG(NELr) )
      DO N=1,NELr
          E = NEL(N,STN)
          ATSTRNG(N) = ATSYM(E)
      ENDDO
!     Triclinc PBC's
      LVr(:,:) = LV(:,:,STN)
      ALLOCATE( Fr(NDIM,NAr) ) 
      CALL fracR(STN,LVr,1)
!
      OPEN(UNIT=60,FILE=OUTCAR,STATUS='unknown')
      WRITE(60,601)  TITLE(STN)
      WRITE(60,602)  LMULT
      WRITE(60,603)  LVr(1,1:3)
      WRITE(60,603)  LVr(2,1:3)
      WRITE(60,603)  LVr(3,1:3)
      WRITE(ELSTRNG,*)  ATSTRNG(:)
      WRITE(60,604)  ELSTRNG
      WRITE(ELSTRNG,*)  ELCNT(1:NELr,STN)
      WRITE(60,605)  ELSTRNG
      WRITE(60,606) 
      WRITE(60,607) 
      DO E=1,ELT(STN)
        DO N=1,ELCNT(E,STN)
         I = ELIST(N,E,STN)
         WRITE(60,608) Fr(:,I),FTAG(:,I,STN)
        ENDDO
      ENDDO
      WRITE(60,*)""
      DO I =1,NAr
          WRITE(60,609) 0.0,0.0,0.0
      ENDDO
!
      DEALLOCATE( ATSTRNG,Fr )
      RETURN
 601  FORMAT(A)
 602  FORMAT(F18.6)
 603  FORMAT(3F18.6)
 604  FORMAT(A)
 605  FORMAT(A)
 606  FORMAT("Selective dynamics")
 607  FORMAT("Direct")
 608  FORMAT(3F18.6,"  ",3A5)
 609  FORMAT(3F18.6)
      END SUBROUTINE WRITE_CAR
