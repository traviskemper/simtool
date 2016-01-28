C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 10/21/2011 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Read tnk file

      SUBROUTINE READ_TNK(STN)
      USE structure
      USE specify
      USE elements  
!
      IMPLICIT none
!
      INTEGER :: I,NAr,ELr,STN,io,N
      CHARACTER(ATSZ) :: AT
      REAL*8 :: XX,YY,ZZ
      CHARACTER(CHSZ) :: ESTAT
!      
!     Verbose output
      IF(VERB) WRITE(6,602) INTNK(STN),STN
!
      OPEN(UNIT=11,FILE=INTNK(STN),STATUS='unknown')
      READ(11,*) NAr
      READ(11,'(A100)') TITLE(STN)
      DO I=NA(STN) + 1,NA(STN) + NAr
         ELr = -1
        
        READ(11,*,iostat=io) N,AT,XX,YY,ZZ,TYP
        CALL getanumb(AT,ELr)

        IF( ELr .LT.0 ) THEN 
          BACKSPACE(11) 
          READ(11,*) ELr,XX,YY,ZZ
!          READ(11,*,iostat=io) ELr,XX,YY,ZZ
 !        IF(io.EQ.0) CALL prnerror(79,ESTAT)
        ENDIF      
        ELN(I,STN) = ELr
        R0(:,I,STN) = (/XX,YY,ZZ/)
      ENDDO
      CLOSE(11)
      NA(STN) = NA(STN) + NAr
!
      RETURN
 602  FORMAT("Reading xyz file ",A22," structure #:",I8)
      END SUBROUTINE READ_TNK

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 10/21/2011 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Write xyz file

      SUBROUTINE WRITE_TNK(STN)
      USE specify
      USE structure
      USE elements  
!
      IMPLICIT none
!
      INTEGER :: I,STN,NAr,ELr,Ni,Nj,TYPNr,STNo,Io
      CHARACTER(ATSZ) :: Atp,GTYPr
      CHARACTER(40) :: WNB

! debug
! 
! debug
!         DO I=1,NA(STN)
!           WRITE(109,*) I,TYPN(I,STN)
!        ENDDO 

!     Verbose output
      IF(VERB) WRITE(6,'("Writing xyz file ",A70)')OUTTNK
!     Write file
      NAr = NA(STN)
      OPEN(UNIT=21,FILE=OUTTNK,STATUS='unknown')
      WRITE(21,*) NAr
      WRITE(21,'(A100)') TITLE(1)
      DO I=1,NAr
        Elr = ELN(I,STN)
        Atp = ATSYM(ELr)
        TYPNr = TYPN(I,STN)
        GTYPr = GRID(I,STN)
        Ni = NINDX(I,STN)
        Nj = NINDX(I+1,STN) - 1
        WRITE(WNB,*) NLIST(Ni:Nj,STN) 
        WRITE(21,2102) I,Atp,R0(:,I,STN),TYPNr,WNB
! 
! debug
!        WRITE(110,*) I,TYPNr
! 
!        WRITE(21,2103) I,Atp,R0(:,I,STN),GTYPr,WNB
      ENDDO 
      CLOSE(21)
!
      RETURN
 2101 FORMAT(I6,3F12.3)
 2102 FORMAT(I6,A6,"  ",3F12.3,I6,A)
 2103 FORMAT(I6,A6,"  ",3F12.3,A6,A)
      END SUBROUTINE WRITE_TNK



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      Read top file 
!
      SUBROUTINE READ_TPA(STN)
      USE specify
      USE structure
      USE potential
      USE elements  
!
      IMPLICIT none

      INTEGER :: STN,Io,Ii,istat,rstat
     & ,RNr,CN,N,I,J,K,L,MN,Mi,Mf,MNi,AMN
     & ,Ia,Ib,In,Id,A,M,Mo,ATN
     & ,MCNTl,MOLCNTr,CNT,Ni,Nf 
      INTEGER, ALLOCATABLE :: MPNTl(:)
      REAL*8 :: CH,MA 
      CHARACTER(CHSZ) :: RLINE,LLINE,LI,var,kywd,FIL2,ESTAT
      CHARACTER(IDSZ) :: AT,RD,GT
      LOGICAL :: STRN
!
!     Initialize 
      Ia = 0
      Ib = 0
      In = 0 
      Id = 0
      MCNTl = 0 
      CNT = 0 
      MOLCNTr = 0 
      ALLOCATE( MPNTl(NAMX) )
!
      IF(VERB) WRITE(6,*)' Readin in tinker ff numbers from ',INTPA(STN)
!
      OPEN (UNIT=22,FILE=INTPA(STN),STATUS='unknown')     
      READ(22,*,iostat=istat) RLINE
      BACKSPACE(22)
      DO WHILE (ISTAT.EQ.0)
        LLINE = ADJUSTL(RLINE)
        LI = LLINE(1:1) 
        IF( LLINE.EQ.'#include' ) THEN 
            READ(22,*,iostat=rstat) var,FIL2
            IF(VERB)WRITE(*,*) 'Opeing sub top file  :',FIL2
            OPEN (UNIT=23,FILE=FIL2,STATUS='unknown')    
            READ(23,*,iostat=istat) RLINE 
            BACKSPACE(23)
            DO WHILE (ISTAT.EQ.0)
                 LLINE = ADJUSTL(RLINE)
                 LI = LLINE(1:1) 
                 IF ( LI.EQ.'[' ) THEN
                    Io=SCAN( RLINE,'[' )      
                    Ii=SCAN( RLINE,']')      
                    IF ( Ii.LT.Io ) THEN
                       READ(23,*) var,KYWD 
                    ELSE 
                       KYWD = LLINE(Io+1:Ii-1)
                       READ(23,*)
                    ENDIF
                    IF( KYWD.EQ.'moleculetype' ) THEN
                      MCNTl = MCNTl + 1     
                      MPNTl(MCNTl) = Ia + 1
                    ELSEIF ( KYWD.EQ.'atoms' ) THEN          
                       READ(23,*,iostat=rstat) N,ATN
                       DO WHILE ( RSTAT.EQ.0 )
                          Ia = Ia + 1
                          TYPN(Ia,STN) = ATN
                         READ(23,*,iostat=rstat) N,ATN 
                       ENDDO         
                       BACKSPACE(23)   
                   ENDIF 
                 ELSE
                    READ(23,*)
                 ENDIF
                 READ(23,*,iostat=istat) RLINE
                 BACKSPACE(23)
              ENDDO
              CLOSE(23)
              IF( VERB ) THEN
                   WRITE(*,*) " Finished reading ", FIL2
                   WRITE(*,*) "   molecule count ",MCNTl
                    WRITE(*,*) "   atom count ",Ia
              ENDIF
        ELSEIF ( LI.EQ.'[' ) THEN
          Io=SCAN( RLINE,'[' )      
          Ii=SCAN( RLINE,']')      
          IF ( Ii.LT.Io ) THEN
             READ(22,*) var,KYWD 
          ELSE 
             KYWD = LLINE(Io+1:Ii-1)
             READ(22,*)
          ENDIF
          IF( KYWD.EQ.'moleculetype' ) THEN
            MCNTl = MCNTl + 1     
            MPNTl(MCNTl) = Ia + 1
          ELSEIF ( KYWD.EQ.'atoms' ) THEN          
            READ(22,*,iostat=rstat) N,ATN
            DO WHILE ( RSTAT.EQ.0 )
              Ia = Ia + 1
              TYPN(Ia,STN)  = ATN
              READ(22,*,iostat=rstat) N,ATN
            ENDDO        
            BACKSPACE(22)   
              IF( VERB ) THEN
                   WRITE(*,*) " Finished reading atoms section "
                   WRITE(*,*) "   molecule count ",MCNTl
                    WRITE(*,*) "   atom count ",Ia
              ENDIF
!
          ELSEIF( KYWD.EQ.'molecules' ) THEN
            MPNTl(MCNTl+1) = Ia + 1
            CNT = 0 
            MOLCNTr = 0 
            MCNTl = 0 
            READ(22,*,iostat=rstat) var,AMN 
            DO WHILE ( RSTAT.EQ.0 )
              MCNTl = MCNTl + 1
              Ni = MPNTl(MCNTl)
              Nf = MPNTl(MCNTl+1) - 1 
              IF(VERB) THEN
                WRITE(*,*) "Replicating mol",MCNTl,AMN
                WRITE(*,*) "   atom count ",Nf - Ni + 1
              ENDIF
              DO MN=1,AMN 
                  DO Io =Ni,Nf
                    CNT = CNT + 1
                    IF( CNT .GT. NAMX ) CALL prnerror(-4,ESTAT)
                    TYPN(CNT,STN)  = TYPN(Io,STN)         ! Force field type number (tinker) 

! debug
                WRITE(88,*) CNT,TYPN(CNT,STN)

                  ENDDO 
              ENDDO 
              READ(22,*,iostat=rstat) var,AMN 
            ENDDO           
            BACKSPACE(22)   

          ENDIF 
       ELSE
         READ(22,*)
       ENDIF
       READ(22,*,iostat=istat) RLINE
       BACKSPACE(22)
      ENDDO
      CLOSE(22)

      IF(VERB) THEN
        WRITE(*,*) " Atoms read in",Ia
        WRITE(*,*) " Molecules read in",MOLCNTr
      ENDIF
!
      DEALLOCATE( MPNTl )
!
      RETURN
      END SUBROUTINE READ_TPA
