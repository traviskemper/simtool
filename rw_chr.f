C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 09/21/2012 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Read chr file

      SUBROUTINE READ_MCHR(STN)
      USE structure
      USE specify
      USE elements  
!
      IMPLICIT none
!
      INTEGER :: STN,MN,MTOT,N,I,K,Mf,Mi
      REAL*8 :: DMB,CHR
      CHARACTER(CHSZ) :: ESTAT,QFILE
      CHARACTER(IDSZ) :: AT
!
!     Set local variables 
!
      MTOT = MOLCNT(STN)
!
!       Verbose
!
        IF( VERB ) WRITE(6,*) ' Reading in charge files for all'
     & ,MTOT,' in structure',STN
!

      DO MN=1,MTOT

        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1  
!
!       Open char file 
! 
        IF( MN.LT.10) THEN
          WRITE(QFILE,7101) MN
        ELSEIF (MN.LT.100) THEN
          WRITE(QFILE,7102) MN
        ELSEIF (MN.LT.1000) THEN
          WRITE(QFILE,7103) MN
        ELSEIF (MN.LT.10000) THEN
          WRITE(QFILE,7104) MN
        ELSEIF (MN.LT.100000) THEN
          WRITE(QFILE,7105) MN
        ELSEIF (MN.LT.1000000) THEN
          WRITE(QFILE,7106) MN
        ELSE
          WRITE(ESTAT,*)  1000000
          CALL prnerror(-11,ESTAT)
        ENDIF

!
!       Read in char file 
!
        OPEN(UNIT=71,FILE=QFILE,STATUS='UNKNOWN')

        DO N=Mi,Mf
          I = MOLST(N,STN) 
          READ(71,*) K,AT,DMB,CHR
          ACHG(I,STN) = CHR 
        ENDDO 
        CLOSE(71)
      ENDDO ! MTOT
!
      RETURN
 7101 FORMAT('q_0_',I1,'.chr')
 7102 FORMAT('q_0_',I2,'.chr')
 7103 FORMAT('q_0_',I3,'.chr')
 7104 FORMAT('q_0_',I4,'.chr')
 7105 FORMAT('q_0_',I5,'.chr')
 7106 FORMAT('q_0_',I6,'.chr')
        END SUBROUTINE READ_MCHR

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C     Version 1.0 09/21/2012 T. W. Kemper                      C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C  Read chr file

      SUBROUTINE READ_CHR(STN)
      USE structure
      USE specify
      USE elements  
!
      IMPLICIT none
!
      INTEGER :: STN,I,J,NAr
      REAL*8 :: CHR
      CHARACTER(CHSZ) :: TYPEi

      OPEN(UNIT=72,FILE=INCHR(STN),STATUS='unknown')
      READ(72,*) NAr
      READ(72,*) 
      DO I=1,NAr
        READ(72,*) J,TYPEi,ACHG(I,STN)
      ENDDO
      CLOSE(72)

      END SUBROUTINE READ_CHR
