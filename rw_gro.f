!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Read in gromacs structure file .gro
!
      SUBROUTINE read_gro(STN)
      USE specify
      USE structure 
      USE elements
!
      IMPLICIT none
!
      INTEGER :: STN,I,N,RNr,NAr,rstat
      REAL*8  :: Rr1,Rr2,Rr3,Rr4,Rr5,Rr6,Rr7,Rr8,Rr9
      CHARACTER(70) :: RLINE,vi
      CHARACTER(IDSZ) :: RD,AN
!
!     Verbose output
      IF(VERB) WRITE(6,602) INGRO(STN),STN

!
      OPEN(UNIT=17,FILE=INGRO(STN),STATUS='unknown')
      READ(17,'(A100)') TITLE(STN)
      READ(17,*) NAr 
      DO I=NA(STN)+1,NA(STN) + NAr
        READ(17,5802,iostat=rstat)RNr,RD,AN,N,Rr1,Rr2,Rr3,Rr4,Rr5,Rr6
        IF( rstat.NE.0 ) THEN
          BACKSPACE(17)
          READ(17,5803,iostat=rstat)RNr,RD,AN,N,Rr1,Rr2,Rr3
          IF(  rstat.NE.0 ) THEN
             WRITE(*,*) "error in ",INGRO(STN)
     &        ," will try unformated read in "
            BACKSPACE(17)
            READ(17,*,iostat=rstat)RNr,RD,AN,N,Rr1,Rr2,Rr3
          ENDIF
          Rr4 = 0 
          Rr5 = 0 
          Rr6 = 0
        ENDIF
            R0(1,I,STN) = Rr1*10.0d0
            R0(2,I,STN) = Rr2*10.0d0
            R0(3,I,STN) = Rr3*10.0d0
            R1(1,I,STN) = Rr4*10.0d0
            R1(2,I,STN) = Rr5*10.0d0
            R1(3,I,STN) = Rr6*10.0d0
            RN(I,STN)   = RNr
            RID(I,STN)  = RD
            GRID(I,STN)  = AN
            MOLN(I,STN) = RNr
      ENDDO
      Rr4 = 0.0d0
      Rr5 = 0.0d0
      Rr6 = 0.0d0
      Rr7 = 0.0d0
      Rr8 = 0.0d0
      Rr9 = 0.0d0
      READ(17,*,iostat=rstat) Rr1,Rr2,Rr3,Rr4,Rr5,Rr6,Rr7,Rr8,Rr9
      IF( rstat.NE.0) THEN
         IF(VERB) WRITE(*,*) INGRO(STN)," turncated ",Rr1,Rr2,Rr3
          Rr4 = 0.0d0 
          Rr5 = 0.0d0 
          Rr6 = 0.0d0 
          Rr7 = 0.0d0 
          Rr8 = 0.0d0 
          Rr9 = 0.0d0 
      ENDIF
      LV(1,1,STN) = Rr1*10.d0
      LV(2,2,STN) = Rr2*10.d0
      LV(3,3,STN) = Rr3*10.d0
      LV(1,2,STN) = Rr4*10.d0
      LV(1,3,STN) = Rr5*10.d0
      LV(2,1,STN) = Rr6*10.d0
      LV(2,3,STN) = Rr7*10.d0
      LV(3,1,STN) = Rr8*10.d0
      LV(3,2,STN) = Rr9*10.d0
      CLOSE(17)
      NA(STN) = NA(STN) + NAr
!
      RETURN
 602  FORMAT("Reading gro file ",A22," structure #:",I8)
 5802 FORMAT(I5,A5,A5,I5,6F8.3)
 5803 FORMAT(I5,A5,A5,I5,3F8.3)
      END SUBROUTINE READ_GRO

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      Write gromacs structure file 

      SUBROUTINE WRITE_GRO(STN)
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: I,STN,NAr
      CHARACTER(5) :: RD,AN
      REAL*8 :: Rr1,Rr2,Rr3,Rr4,Rr5,Rr6,Rr7,Rr8,Rr9
!
!     Verbose output
      IF(VERB) WRITE(6,'("Writing gro file ",A70)') OUTGRO
!     Write file
      NAr = NA(STN)
      OPEN(UNIT=58,FILE=OUTGRO,STATUS='unknown')
      WRITE(58,'(A100)') TITLE(1)
      WRITE(58,*) NAr
      DO I=1,NAr
        RD = ADJUSTL(RID(I,STN))
        AN = ADJUSTR( GRID(I,STN) )
        Rr1 = R0(1,I,STN)/10.0d0
        Rr2 = R0(2,I,STN)/10.0d0
        Rr3 = R0(3,I,STN)/10.0d0
        Rr4 = R1(1,I,STN)/10.0d0
        Rr5 = R1(2,I,STN)/10.0d0
        Rr6 = R1(3,I,STN)/10.0d0
!        WRITE(58,5802) RNi(STN,I),RD,AN,I,Rr1,Rr2,Rr3,Rr4,Rr5,Rr6
        WRITE(58,5802) RN(I,STN),RD,AN,I,Rr1,Rr2,Rr3,Rr4,Rr5,Rr6
!        WRITE(58,*) MOLNi(STN,I),RD,AN,I,Rr1,Rr2,Rr3,Rr4,Rr5,Rr6
       ENDDO 
       Rr1 = LV(1,1,STN)/10.d0
       Rr2 = LV(2,2,STN)/10.d0
       Rr3 = LV(3,3,STN)/10.d0
       Rr4 = LV(1,2,STN)/10.d0
       Rr5 = LV(1,3,STN)/10.d0
       Rr6 = LV(2,1,STN)/10.d0
       Rr7 = LV(2,3,STN)/10.d0
       Rr8 = LV(3,1,STN)/10.d0
       Rr9 = LV(3,2,STN)/10.d0
       WRITE(58,5803) Rr1,Rr2,Rr3,Rr4,Rr5,Rr6,Rr7,Rr8,Rr9
       CLOSE(58)

      RETURN
 5801 FORMAT(I5,A4,A6,I5,"  ",F6.3,"  ",F6.3,"  ",F6.3)
 5802 FORMAT(I5,A5,A5,I5,6F8.3)
 5803 FORMAT(9F12.4)
      END SUBROUTINE WRITE_GRO
