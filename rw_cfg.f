!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
!  Read cfg file
!
      SUBROUTINE READ_CFG(STN)
      USE structure
      USE specify
!
      IMPLICIT none
!
      INTEGER :: I 
      REAL*8 :: SCL,A11,A12,A13,A21,A22,A23,A31,A32,A33
      CHARACTER(70) :: RLINE
!

      OPEN(UNIT=18,FILE=INCFG(STN),STATUS='unknown')
      READ(53,'(A70)',iostat=istat) RLINE
      BACKSPACE(53)
      DO WHILE (ISTAT.EQ.0)
        IF( RLINE() .EQ."Number of particles = " )THEN
          READ(18,*) Aa,Ab,Ac,Ad,NAo(STN)
        ELSEIF(RLINE(1:3).EQ."A =") THEN
          READ(18,*) Aa,Aa,SCL
        ELSEIF(RLINE(1:7).EQ."H0(1,1)") READ(18,*) Aa,Ab,A11
        ELSEIF(RLINE(1:7).EQ."H0(1,2)") READ(18,*) Aa,Ab,A12
        ELSEIF(RLINE(1:7).EQ."H0(1,3)") READ(18,*) Aa,Ab,A13
        ELSEIF(RLINE(1:7).EQ."H0(2,1)") READ(18,*) Aa,Ab,A21
        ELSEIF(RLINE(1:7).EQ."H0(2,2)") READ(18,*) Aa,Ab,A22
        ELSEIF(RLINE(1:7).EQ."H0(2,3)") READ(18,*) Aa,Ab,A23
        ELSEIF(RLINE(1:7).EQ."H0(3,1)") READ(18,*) Aa,Ab,A31
        ELSEIF(RLINE(1:7).EQ."H0(3,2)") READ(18,*) Aa,Ab,A32
        ELSEIF(RLINE(1:7).EQ."H0(3,3)") READ(18,*) Aa,Ab,A33
        ELSEIF( RLINE is A ) READ(18,*) Aa
        ELSEIF( RLINE is INT ) READ(18,*) mass
        ELSEIF( RLINE is float ) THEN
          I = I + 1
          READ(18,*) RF(1:3)
          R0o(STN,I,1) = RF(1)*A11 + RF(2)*A21 + RF(3)*A31
          R0o(STN,I,2) = RF(1)*A12 + RF(2)*A22 + RF(3)*A32
          R0o(STN,I,3) = RF(1)*A13 + RF(2)*A23 + RF(3)*A33
        ENDIF
      ENDDO 
      LVo(STN,1,1) = A11
      LVo(STN,1,1) = A11
      LVo(STN,1,1) = A11
      LVo(STN,1,1) = A11
      LVo(STN,1,1) = A11
      LVo(STN,1,1) = A11
      LVo(STN,1,1) = A11
      CLOSE(18)
      END SUBROUTINE READ_CFG

      

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
!  Write cfg file
!
      SUBROUTINE setcfg(STN)
!
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: STN
      REAL*8 :: CBOX(NDIM)
!
      N = 1
!
      IF ( N.GT.999999) THEN
        WRITE(6,*) 'Warning more than 999999 cfg files are needed'
        WRITE(6,*) ' change filename format in cfg.f'
        STOP
      ENDIF
!     Set box size 
      CBOX(:) = LC(STN,:)
!
!     Create Cfgs folder
      CALL SYSTEM('sh -c "if [ ! -d Cfgs ]; then mkdir Cfgs ; fi"')
!      
      RETURN
      END SUBROUTINE setcfg 
! 

! Create cfg files for atomeye 
! 
      SUBROUTINE write_cfg(STN)
!
!
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: STN
      REAL*8 :: CBOX(NDIM,NDIM)
!
      CBOX(:,:) = LVi(STN,:,:)
!
      OPEN(UNIT=LP,FILE=OUTCFG,STATUS='UNKNOWN')
      write(LP,10)'Number of particles = ',NAi(SNT)
      write(LP,11)'A = 1.0'
      write(LP,12)'H0(1,1) =', cbox(1,1)
      write(LP,12)'H0(1,2) =', cbox(1,2) 
      write(LP,12)'H0(1,3) =', cbox(1,3) 
      write(LP,12)'H0(2,1) =', cbox(2,1) 
      write(LP,12)'H0(2,2) =', cbox(2,2)
      write(LP,12)'H0(2,3) =', cbox(2,3) 
      write(LP,12)'H0(3,1) =', cbox(3,1) 
      write(LP,12)'H0(3,2) =', cbox(3,2) 
      write(LP,12)'H0(3,3) =', cbox(3,3)
      write(LP,11)'.NO_VELOCITY.'
      write(LP,11)'entry_count = 13'
      write(LP,11)'auxiliary[0]= Atom #'
      write(LP,11)'auxiliary[1]= Molecule number [#]'
      write(LP,11)'auxiliary[2]= Force [A/fs]'
      write(LP,11)'auxiliary[3]= Groups [#]'
      write(LP,11)'auxiliary[4]= Surface atom [#]'
!     Loop over all atom types 
      DO K = 1,NTYPES
        IF( NOEL(K).GT.0 ) THEN
          WRITE(LP,*) XMASS(K)
          WRITE(LP,*) ATYPE(K)
          DO I=1,NAi(STN)
            EL = ELNi(STN,I)
            IF(  .EQ. 
          ENDDO      
        ENDIF
      ENDDO




      N = 1 

!
      INTEGER :: N,I,K,LP,D,J
      REAL*8, DIMENSION(3) :: R0F,R1V,R2A
      REAL*8 :: RNPSQ,PRNRNP,VSQ,ASQ,MFRAC,MBF,XX,TEMPKi
      CHARACTER *79 FILENAME
!
      N = LSTEP / NSTR + 1
      LP = 19
!     Calculate buffer space between beam molecules under substrate
      MFRAC=RMN(1,DEPD)/CBOX(DEPD) + 0.5d0
      MBF= MFRAC/NC/NPM
!     Write output
      WRITE (FILENAME, '("Cfgs/", I6.6, ".cfg")'),N
      OPEN(UNIT=LP,FILE=FILENAME,STATUS='UNKNOWN')
      write(LP,10)'Number of particles = ',NP
      write(LP,11)'A = 1.0'
      write(LP,12)'H0(1,1) =', cbox(1)
      write(LP,12)'H0(1,2) =', 0.
      write(LP,12)'H0(1,3) =', 0.
      write(LP,12)'H0(2,1) =', 0.
      write(LP,12)'H0(2,2) =', cbox(2)
      write(LP,12)'H0(2,3) =', 0.
      write(LP,12)'H0(3,1) =', 0.
      write(LP,12)'H0(3,2) =', 0.
      write(LP,12)'H0(3,3) =', cbox(3)
      write(LP,11)'.NO_VELOCITY.'
      write(LP,11)'entry_count = 13'
      write(LP,11)'auxiliary[0]= Atom #'
      write(LP,11)'auxiliary[1]= Thermostat'
      write(LP,11)'auxiliary[2]= Temperature [K]'
      write(LP,11)'auxiliary[3]= Molecule number [#]'
      write(LP,11)'auxiliary[4]= Reaction number [#]'
      write(LP,11)'auxiliary[5]= Force [A/fs]'
      write(LP,11)'auxiliary[6]= Ave_velocity [A/fs]'
      write(LP,11)'auxiliary[7]= Ave_acceleration [A/fs]'
      write(LP,11)'auxiliary[8]= Groups [#]'
      write(LP,11)'auxiliary[9]= Surface atom [#]'
!     Loop over all atom types 
      DO K = 1,NTYPES
        IF( NOA(K).GT.0 ) THEN
          WRITE(LP,*) XMASS(K)
          WRITE(LP,*) ATYPE(K)
          DO I=1,NP
            IF ( KTYPE(I).EQ.K) THEN
!             Calculate fractional coordinates  
              VSQ = 0.0d0
              ASQ = 0.0d0         
              DO D =1,3
                R0F(D) = R0(I,D)/CBOX(D) + .50d0
                R1V(D) = R1(I,D)/DELTA
                VSQ = VSQ + R1V(D)**2
                R2A(D) = R2(I,D)/DELTSQ
                ASQ = ASQ + R2A(D)**2
              ENDDO 
!             Adjust beam molecules
              IF(R0F(DEPD).GT.1.d0) THEN
                MFRAC = MFRAC - MBF
                R0F(DEPD) = MFRAC
                R1V(1:3) = (/0.0d0,0.0d0,0.d0/)
                R2A(1:3) = (/0.0d0,0.0d0,0.d0/)
                ASQ = 0.0d0
              ENDIF
!             Caclulate force 
              RNPSQ = 0.0d0
              DO  J=1,3
                RNPSQ = RNPSQ + RNP(I,J)**2
              ENDDO
              PRNRNP = SQRT(RNPSQ)
              XX = VSQ*XMASS(KTYPE(I))
              TEMPKi = ENPR*XX/6.0d0/DELTSQ*ECONV
              WRITE(LP,13)R0F(1:3),I,ITR(I),TEMPKi,IMOLi(I),RCLIST(I)
     &          ,PRNRNP,SQRT(VSQ),SQRT(ASQ),GRPi(I),SURFI(I)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      CLOSE(LP)
!
10    format(a,1x,i6)
11    format(a)
12    format(a,1x,f26.6)
13    format(3(E12.4,1x),I8,I5,E12.4,1x,2I8,3(E12.4,1x),2I8)
14    format(3E20.11)
      RETURN
      END SUBROUTINE cfg


! Create cfg files for atomeye 
! 
      SUBROUTINE cfg_group
!
      USE STRUCTR 
      USE POTS 
      USE MDSTP
      USE SPECIF
      USE BEAM
      USE ANALYSIS
!
      IMPLICIT none 
!
      INTEGER :: N,I,K,LP,D,J,SUBCNT
      REAL*8, DIMENSION(3) :: R0F,R1V,R2A
      REAL*8 :: RNPSQ,PRNRNP,VSQ,ASQ,MFRAC,MBF,XX,TEMPKi
      CHARACTER *79 FILENAME
!
      N = LSTEP / NSTR + 1
      LP = 101
!     Calculate buffer space between beam molecules under substrate
      MFRAC=RMN(1,DEPD)/CBOX(DEPD) + 0.5d0
      MBF= MFRAC/NC/NPM
!     Write output
      SUBCNT = 0
      DO I=1,NP
        IF( GRPi(I).NE.0 ) SUBCNT = SUBCNT +1
      ENDDO 
      OPEN(UNIT=LP,FILE='group.cfg',STATUS='UNKNOWN')
      write(LP,10)'Number of particles = ',SUBCNT
      write(LP,11)'A = 1.0'
      write(LP,12)'H0(1,1) =', cbox(1)
      write(LP,12)'H0(1,2) =', 0.
      write(LP,12)'H0(1,3) =', 0.
      write(LP,12)'H0(2,1) =', 0.
      write(LP,12)'H0(2,2) =', cbox(2)
      write(LP,12)'H0(2,3) =', 0.
      write(LP,12)'H0(3,1) =', 0.
      write(LP,12)'H0(3,2) =', 0.
      write(LP,12)'H0(3,3) =', cbox(3)
      write(LP,11)'.NO_VELOCITY.'
      write(LP,11)'entry_count = 13'
      write(LP,11)'auxiliary[0]= Atom #'
      write(LP,11)'auxiliary[1]= Thermostat'
      write(LP,11)'auxiliary[2]= Temperature [K]'
      write(LP,11)'auxiliary[3]= Molecule number [#]'
      write(LP,11)'auxiliary[4]= Reaction number [#]'
      write(LP,11)'auxiliary[5]= Force [A/fs]'
      write(LP,11)'auxiliary[6]= Ave_velocity [A/fs]'
      write(LP,11)'auxiliary[7]= Ave_acceleration [A/fs]'
      write(LP,11)'auxiliary[8]= Groups [#]'
      write(LP,11)'auxiliary[9]= Surface atom [#]'
!     Loop over all atom types 
      DO K = 1,NTYPES
        IF( NOA(K).GT.0 ) THEN
          WRITE(LP,*) XMASS(K)
          WRITE(LP,*) ATYPE(K)
          DO I=1,NP
            IF ( KTYPE(I).EQ.K.AND.GRPi(I).NE.0) THEN
!             Calculate fractional coordinates  
              VSQ = 0.0d0
              ASQ = 0.0d0         
              DO D =1,3
                R0F(D) = R0(I,D)/CBOX(D) + .50d0
                R1V(D) = R1(I,D)/DELTA
                VSQ = VSQ + R1V(D)**2
                R2A(D) = R2(I,D)/DELTSQ
                ASQ = ASQ + R2A(D)**2
              ENDDO 
!             Adjust beam molecules
              IF(R0F(DEPD).GT.1.d0) THEN
                MFRAC = MFRAC - MBF
                R0F(DEPD) = MFRAC
                R1V(1:3) = (/0.0d0,0.0d0,0.d0/)
                R2A(1:3) = (/0.0d0,0.0d0,0.d0/)
                ASQ = 0.0d0
              ENDIF
!             Caclulate force 
              RNPSQ = 0.0d0
              DO  J=1,3
                RNPSQ = RNPSQ + RNP(I,J)**2
              ENDDO
              PRNRNP = SQRT(RNPSQ)
              XX = VSQ*XMASS(KTYPE(I))
              TEMPKi = ENPR*XX/6.0d0/DELTSQ*ECONV
              WRITE(LP,13)R0F(1:3),I,ITR(I),TEMPKi,IMOLi(I),RCLIST(I)
     &          ,PRNRNP,SQRT(VSQ),SQRT(ASQ),GRPi(I),SURFI(I)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      CLOSE(LP)
!
10    format(a,1x,i6)
11    format(a)
12    format(a,1x,f26.6)
13    format(3(E12.4,1x),I8,I5,E12.4,1x,2I8,3(E12.4,1x),2I8)
14    format(3E20.11)
      RETURN
      END SUBROUTINE cfg_group
