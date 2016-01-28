!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 12/01/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Read top reference file
!
      SUBROUTINE READ_TRF(STN)
      USE specify
      USE structure
      USE potential
!
      IMPLICIT none
      INTEGER :: I,N,STN,NA,ISTAT,RN,EL
      REAL*8 :: MSt,CHR,CHRT
      CHARACTER(CHSZ) :: RD,RLINE,var
      CHARACTER(5) :: AN,AT
!
      OPEN (UNIT=16,FILE=INTOP(STN),STATUS='unknown')
      READ(16,*,iostat=istat) RLINE
      BACKSPACE(16)
      DO WHILE (ISTAT.EQ.0)
       IF ( RLINE(1:5).EQ.'atoms' ) THEN
         READ(16,*) var,NA
         IF( NA.NE.NA(STN) )THEN
           WRITE(6,*) 'Error number of atoms in str and top'
           STOP
         ENDIF
         CHRT = 0.0d0 
         DO I=1,NA(STN) !WHILE (ISTAT2.EQ.0 ) THEN
           READ(16,*) N,AT,RN,RD,AN,RN(I,STN)
     &    ,CHR,MSt,EL
           ELN(I,STN) = EL
           RID(I,STN) = RD
           ATYP(I,STN) = AT
           ANM(I,STN) = AN
           MASS(EL) = MSt
           CHRG(I,STN) = CHR
           CHRT = CHRT + CHR
         ENDDO
         IF( VERB ) THEN
           WRITE(6,*) NAo(STN),'atoms read from ',INTOP(STN) 
           WRITE(6,*) ' Total charge ',CHRT
         ENDIF
       ELSEIF( RLINE(1:11).EQ.'constraints' ) THEN
         READ(16,*) var,NCo(STN)
         DO I=1,NCo(STN) !WHILE (ISTAT2.EQ.0 ) THEN
           READ(16,*) CONI(I,STN) ,CONJ(I,STN),CTYP(I,STN)
     &               ,CDISI(I,STN),CDISJ(I,STN)
        ENDDO
         IF( VERB ) WRITE(6,*) NCo(STN),' constraints read from '
     &             ,INTOP(STN)
       ELSEIF( RLINE(1:5).EQ.'bonds' ) THEN
         READ(16,*) var,NBo(STN)
         DO I=1,NBo(STN) !WHILE (ISTAT2.EQ.0 ) THEN
           READ(16,*) BNDI(I,STN),BNDJ(I,STN),BTYP(I,STN)
         ENDDO
         IF( VERB ) WRITE(6,*) NBo(STN),' bonds read from ',INTOP(STN)
       ELSEIF( RLINE(1:5).EQ.'pairs' ) THEN
         READ(16,*) var,NPRo(STN)
         DO I=1,NPRo(STN) !WHILE (ISTAT2.EQ.0 ) THEN
           READ(16,*) PRSI(I,STN),PRSJ(I,STN),PRTYP(I,STN)
         ENDDO
       ELSEIF( RLINE(1:6).EQ.'angles' ) THEN
         READ(16,*) var,NANGo(STN)
         IF(VERB)WRITE(6,*)' Reading',NANGo(STN),' angles '
         DO I=1,NANGo(STN) !WHILE (ISTAT2.EQ.0 ) THEN
           READ(16,*) ANGI(I,STN),ANGJ(I,STN),ANGK(I,STN)
     &     ,ANTYP(I,STN)
         ENDDO
         IF(VERB)WRITE(6,*) NANGo(STN),' angles read from ',INTOP(STN)
       ELSEIF( RLINE(1:9).EQ.'dihedrals' ) THEN
         READ(16,*) var,NDIHo(STN)
         DO I=1,NDIHo(STN) !WHILE (ISTAT2.EQ.0 ) THEN
           READ(16,*)DIHI(I,STN),DIHJ(I,STN),DIHK(I,STN)
     &    ,DIHL(I,STN),DHTYP(I,STN)
         ENDDO
         IF(VERB)WRITE(6,*) NDIH(STN),' dih read from ',INTOP(STN)
       ELSE
         READ(16,*)
       ENDIF
       READ(16,*,iostat=istat) RLINE
       BACKSPACE(16)
      ENDDO
!
      RETURN
      END SUBROUTINE READ_TOP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/02/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Write .top file for GROMACS input
!
      SUBROUTINE WRITE_TOP(STN)
      USE specify
      USE structure
      USE potential 
!
      IMPLICIT none
!
      INTEGER   :: I,J,N,EL,D,STN,RN
      CHARACTER(CHSZ) :: NBID(12)
      CHARACTER(5) :: RD,AT,AN
!
      WRITE(6,*) 'Writing xyz coordinates',OUTXYZ
      OPEN(UNIT=53,FILE=OUTTOP,STATUS='unknown')
      OPEN(UNIT=59,FILE='ref.top',STATUS='unknown')

      WRITE(53,'(A)') '; Include forcefield parameters'
      WRITE(53,'(A)') '#include "./mAlq3-13.ff/forcefield.itp"'
      WRITE(53,'(A)') '  '
      WRITE(53,'(A)') '[ moleculetype ] '
      WRITE(53,'(A)') '; Name            nrexcl '
      WRITE(53,'(A)') 'Protein             3    '
      WRITE(53,'(A)') ''
      WRITE(53,'(A)') '[ atoms ]'
      WRITE(59,'(A,I)') ' atoms ',NAi(STN)
      DO I=1,NAi(STN)
        RD = ADJUSTL( RIDi(STN,I))
        AT = ADJUSTR( ATYPi(STN,I) )
        AN = ADJUSTR( ANMi(STN,I) )
        EL = ELNi(STN,I)
        RN = RNi(STN,I)
        WRITE(53,5303) I,AT,RN,RD,AN,RNi(STN,I)
     &    ,CHRGi(STN,I),MASS(EL)
        WRITE(59,5903) I,AT,RN,RD,AN,RNi(STN,I)
     &    ,CHRGi(STN,I),MASS(EL),EL
      ENDDO
      WRITE(53,*) ''
      WRITE(59,*) ''
      IF( NCi(STN).GT.0) THEN
       WRITE(53,'(A)') '[ constraints ]'
       WRITE(59,'(A,I)') ' constraints',NCi(STN)
       DO I=1,NCi(STN)
        WRITE(53,5304)  CONIi(STN,I) ,CONJi(STN,I),CTYPi(STN,I)
     &               ,CDISIi(STN,I),CDISJi(STN,I)
        WRITE(59,5304)  CONIi(STN,I) ,CONJi(STN,I),CTYPi(STN,I)
     &               ,CDISIi(STN,I),CDISJi(STN,I)
       ENDDO 
      ENDIF
      IF( NBi(STN).GT.0) THEN
       WRITE(53,'(A)') '[ bonds ]'
       WRITE(59,'(A,I)') 'bonds ',NBi(STN)
       DO I=1,NBi(STN)
        WRITE(53,*) BNDIi(STN,I),BNDJi(STN,I),BTYPi(STN,I)
        WRITE(59,*) BNDIi(STN,I),BNDJi(STN,I),BTYPi(STN,I)
       ENDDO 
      ENDIF
      WRITE(53,*) ''
      WRITE(59,*) ''
      IF( NPRi(STN).GT.0) THEN
       WRITE(53,'(A)') '[ pairs ]'
       WRITE(59,'(A,I)') ' pairs ',NPRi(STN)
       DO I=1,NPRi(STN)
        WRITE(53,*) PRSIi(STN,I),PRSJi(STN,I),PRTYPi(STN,I)
        WRITE(59,*) PRSIi(STN,I),PRSJi(STN,I),PRTYPi(STN,I)
       ENDDO 
      ENDIF
      WRITE(53,*) ''
      WRITE(59,*) ''
      IF( NANGi(STN).GT.0) THEN
       WRITE(53,'(A)') '[ angles ]'
       WRITE(59,'(A,I)') ' angles ',NANGi(STN)
       DO I=1,NANGi(STN)
        WRITE(53,*)ANGIi(STN,I),ANGJi(STN,I),ANGKi(STN,I),ANTYPi(STN,I)
        WRITE(59,*)ANGIi(STN,I),ANGJi(STN,I),ANGKi(STN,I),ANTYPi(STN,I)
       ENDDO 
      ENDIF
      WRITE(53,*) ''
      WRITE(59,*) ''
      IF( NDIHi(STN).GT.0) THEN
       WRITE(53,'(A)') '[ dihedrals ]'
       WRITE(59,'(A,I)') ' dihedrals ',NDIHi(STN)
       DO I=1,NDIHi(STN)
        WRITE(53,*)DIHIi(STN,I),DIHJi(STN,I),DIHKi(STN,I)
     &   ,DIHLi(STN,I),DHTYPi(STN,I)
        WRITE(59,*)DIHIi(STN,I),DIHJi(STN,I),DIHKi(STN,I)
     &   ,DIHLi(STN,I),DHTYPi(STN,I)
       ENDDO 
      ENDIF
      WRITE(59,*) ''
      WRITE(53,'(A)') ' '
      WRITE(53,'(A)') '; Include Position restraint file '
!      WRITE(53,'(A)') '#ifdef POSRES '
!      WRITE(53,'(A)') '#include "posre.itp" '
!      WRITE(53,'(A)') '#endif '
      WRITE(53,'(A)') ' '
      WRITE(53,'(A)') '[ system ] '
      WRITE(53,'(A)') '; Name '
      WRITE(53,'(A)') 'Protein '
      WRITE(53,'(A)') ' ' 
      WRITE(53,'(A)') '[ molecules ] '
      WRITE(53,'(A)') '; Compound        #mols '
      WRITE(53,'(A)') 'Protein             1   '
      WRITE(53,'(A)') ' '
      WRITE(53,'(A)') ' '
      CLOSE(53)
      CLOSE(59) 
!
      RETURN
 5301 FORMAT(4I8,' 2')
 5302 FORMAT(4I8,' 2')
 5303 FORMAT(I6,A8,I6,2A8,I6,F16.8,F10.3)
 5903 FORMAT(I6,A8,I6,2A8,I6,F16.8,F10.3,I6)
 5304 FORMAT(2I8,I4,3F10.4)
      END SUBROUTINE WRITE_TOP
