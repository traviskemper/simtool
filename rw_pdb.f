!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Change a xyz file into a .gro for GROMACS
!
      SUBROUTINE read_pdb(STN)
      USE specify
      USE structure 
      USE elements  
!
      IMPLICIT none
!
      INTEGER:: I,J,N,RESNR,NMOLr,ISTAT,CNT,NCHNr,STN,ELNr
      REAL*8 :: KK,LL,R(3),LCr(NDIM),LAr(NDIM),LCCK
      CHARACTER(IDSZ)  :: CR1,CR2,CR3,RIDr,AIDr
      CHARACTER(ATSZ) :: SYMr
      CHARACTER(10) :: POTID(100),ID
      CHARACTER(6) ::  KEYWORD 
      CHARACTER(CHSZ) :: RLINE
      CHARACTER(1) :: CHID
!
!     Verbose output
      IF(VERB) WRITE(6,602) INPDB(STN),STN
!
      OPEN(unit=53,FILE=INPDB(STN),STATUS='unknown')
!     Read in pdb information      
      RESNR = 0
      NMOLr = 1   !# of molecules
      NCHNr = 1   !# of chains 
      I     = NA(STN)
      READ(53,*,iostat=istat) KEYWORD
      BACKSPACE(53)
      DO WHILE (ISTAT.EQ.0)
        IF( KEYWORD(1:6) .EQ. 'REMARK' ) THEN
          READ(53,'(A4,A70)') CR1,TITLE(STN)
        ELSEIF( KEYWORD(1:6) .EQ. 'CRYST1' ) THEN
          READ(53,*) CR1,LCr(:),LAr(:)
        ELSEIF( KEYWORD(1:6) .EQ. 'MODEL' ) THEN
          READ(53,'(A70)') TITLE(STN) 
        ELSEIF( KEYWORD(1:6) .EQ. 'ENDMDL' ) THEN
          READ(53,*)
          NMOLr = NMOLr + 1 
        ELSEIF( KEYWORD(1:6) .EQ. 'TER' ) THEN   
          READ(53,*)       
          NCHNr = NCHNr + 1
        ELSEIF( KEYWORD(1:6) .EQ. 'ATOM' ) THEN
          READ(53,*,iostat=istat) CR1,N,AIDr,RIDr,CHID,J,R(:),KK,LL,SYMr
          IF(istat.NE.0) THEN
            BACKSPACE(53)
            READ(53,*,iostat=istat)  CR1,N,AIDr,RIDr,J,R(:),KK,LL,SYMr
          ENDIF
          CALL GETANUMB(SYMr,ELNr)   
          I = I + 1 
          ELN(I,STN) = ELNr
          R0(:,I,STN) = R(:)
          TYP(I,STN) = ADJUSTL(AIDr)
          GRID(I,STN) = ADJUSTL(AIDr)
          RID(I,STN) = ADJUSTL(RIDr)
          RN(I,STN)  = J
          MOLN(I,STN) =  NMOLr

       ELSEIF( KEYWORD(1:4) .EQ. 'KATM' ) THEN
          READ(53,*,iostat=istat)  CR1,N,AIDr,RIDr,CHID,J,R(:),KK,LL 
          IF(istat.NE.0) THEN
            BACKSPACE(53)
            READ(53,*,iostat=istat)  CR1,N,AIDr,RIDr,J,R(:),KK,LL 
          ENDIF
          IF( istat.NE.0) THEN
            BACKSPACE(53)
            READ(53,'(A60)') RLINE 
            CALL prnerror(82,RLINE)
         ENDIF
         CALL GETANUMB(SYMr,ELNr)   
          IF( AIDr.EQ."Bau" ) ELNr = 79
          IF( AIDr.EQ."VS" ) ELNr = 0
          IF( AIDr(1:1).EQ."C" ) ELNr = 6
          IF( AIDr(1:1).EQ."H" ) ELNr = 1
          IF( AIDr(1:1).EQ."O" ) ELNr = 8
          IF( AIDr(1:1).EQ."S" ) ELNr = 16
 
          I = I + 1 
          ELN(I,STN) = ELNr
          R0(:,I,STN) = R(:)
          TYP(I,STN) = ADJUSTL(AIDr)
          GRID(I,STN) = ADJUSTL(AIDr)
          RID(I,STN) = 'MOL' !ADJUSTL(RIDr)
          RN(I,STN)  = J
          MOLN(I,STN) =  NMOLr
       ELSEIF( KEYWORD(1:6) .EQ. 'HETATM' ) THEN
          READ(53,*,iostat=istat)  CR1,N,AIDr,J,R(:),SYMr
          IF( istat.NE.0) THEN
            BACKSPACE(53)
            READ(53,'(A60)') RLINE 
            CALL prnerror(82,RLINE)
          ENDIF
          I = I + 1 
          CALL GETANUMB(SYMr,ELNr)   
          ELN(I,STN) = ELNr
          R0(:,I,STN) = R(:)
          TYP(I,STN) = ADJUSTL(AIDr)
          GRID(I,STN) = ADJUSTL(AIDr)
          RID(I,STN) = 'MOL'
          RN(I,STN)  = J
          MOLN(I,STN) =  NMOLr
        ELSE
          WRITE(6,*) KEYWORD,' unknown'
          READ(53,*,iostat=istat) KEYWORD
        ENDIF
        READ(53,*,iostat=istat) KEYWORD
        BACKSPACE(53)
      ENDDO
      CLOSE(53)
!     Check lattice constants against user specified 
      IF ( LCRD(STN)) THEN
!       Check user input make user not 0 
        LCCK = LC(1,STN)*LC(2,STN)*LC(3,STN) 
        LCCK = LCCK*LA(1,STN)*LA(2,STN)*LA(3,STN) 
        IF( LCCK.EQ.0.0d0 ) THEN
          IF(VERB)WRITE(6,*)"User defined LC have zeros, will us pdb LC"
          LC(:,STN) = LCr(:)
          LA(:,STN) = LAr(:)
        ENDIF
      ELSE
          LC(:,STN) = LCr(:)
          LA(:,STN) = LAr(:)
      ENDIF
!
      NA(STN) = NA(STN) + I

      RETURN
 602  FORMAT("Reading pdb file ",A22," structure #:",I8)
      END SUBROUTINE read_pdb
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Write pdb file
!
      SUBROUTINE WRITE_PDB(STN)
      USE specify
      USE structure
      USE elements
!
      IMPLICIT none
!
      INTEGER ::I,EL,STN
      REAL*8 :: KK,LL
      CHARACTER(6) :: KWD,KWDi
      CHARACTER(3) :: TYPr 
!
!     Verbose output
      IF(VERB) WRITE(6,'("Writing pdb file ",A70)')OUTPDB
!     Write file
      OPEN(UNIT=52,FILE=OUTPDB,STATUS='unknown')
      KK = 1.0d0
      LL = 0.0d0        
      KWD = ''
      KWD = 'REMARK'
      KWD = ADJUSTL(KWD)
      WRITE(52, '(A6,A100)') KWD,TITLE(STN)
!
      KWD = 'CRYST1'
      KWD = ADJUSTL(KWD)
      WRITE(52,5202) KWD,LC(:,STN),LA(:,STN)
!
      KWD = 'ATOM'
      KWD = ADJUSTL(KWD)
      DO I=1,NA(STN) - 1
         EL=ELN(I,STN)
         TYPr = ADJUSTR( TYP(I,STN) )
         WRITE(52,5211) KWD,I,TYPr,RID(I,STN)
     &                  ,RN(I,STN),R0(:,I,STN)
     &   ,KK,LL,ATSYM(EL)         
         IF( MOLN(I,STN) .NE. MOLN(I,STN) ) THEN
           KWDi = 'ENDMOL'
           KWDi = ADJUSTL(KWDi)
           WRITE(52,5203) KWDi
         ENDIF         
      ENDDO
! Write last atom
      I=NA(STN)
      EL=ELN(I,STN)
      TYPr = ADJUSTR( TYP(I,STN) )
      WRITE(52,5211) KWD,I,TYPr,RID(I,STN)
     &                  ,RN(I,STN),R0(:,I,STN)
     &   ,KK,LL,ATSYM(EL)                 
      KWD = 'ENDMOL'
      KWD = ADJUSTL(KWD)
      WRITE(52,5203) KWD
!
      RETURN
 5201 FORMAT(A4,' ',I6,1X,A4,A5,I5,'    ',3F8.3,2F6.2,A10)
! 5201 FORMAT(A6,I6,2X,A4,A4,I5,'    ',3F8.3,2F6.2,A10)
 5202 FORMAT(A6,3F9.3,3F7.2,' P-1')
!           1-6 7-11  13-16 18-20 23-26 31-38 39-46 47-54 
 5211 FORMAT(A6,I5,' ',A4,' ',A3,'  ',I4,'    ',3F8.3,2F6.2,A10)
 5203 FORMAT(A6)
      END SUBROUTINE write_pdb
