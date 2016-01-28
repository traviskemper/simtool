!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/09/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Write pdb file
!
      SUBROUTINE WRITE_COM(STN)
      USE specify
      USE structure
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: I,EL,STN,GFIX,NAr 
      REAL*8 :: CG
      CHARACTER(IDSZ) :: AT,TP
!
      NAr = NA(STN) 
      OPEN(UNIT=54,FILE=OUTCOM,STATUS='unknown')
      WRITE(54,5403) TITLE(STN)
      WRITE(54,*) '%mem=2GB'
      WRITE(54,5400) 
      WRITE(54,*) ''
      WRITE(54,5402) TITLE(STN)
      WRITE(54,*) ''
      WRITE(54,*) '0 1  0 1 0 1'
      DO I=1,NAr
         EL=ELN(I,STN)
         AT = ATSYM(EL)
         TP = TYP(I,STN)
         CG= ACHG(I,STN)
         GFIX = 0 
         IF( FTAG(3,I,STN) .EQ. "F" ) GFIX = -1
         WRITE(54,5401) AT,TP,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
      ENDDO
      WRITE(54,*) ''
      WRITE(54,*) ''
!      WRITE(54,*) '--Link1--'
!      WRITE(54,*) '%chk=',TITLE
!      WRITE(54,*) '%mem=2GB'
!      WRITE(54,*) '# uff NoSymm OPT'
      WRITE(54,*) ''
      WRITE(54,*) ''
!
      RETURN
 5400 FORMAT('# oniom(b3lyp/6-31g(d):hf/sto-3g:uff)=embedcharge 
     & nosymm  pop=mk sp')
 5401 FORMAT(A2,"-",A2,"-",F6.3,I3,3F9.3," ",A5," ")
 5402 FORMAT(A20)
 5403 FORMAT('%chk=',A20)
      END SUBROUTINE WRITE_COM
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Read gaussian input file
!
      SUBROUTINE READ_COM(STN)
!
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: STN 
!
 
  
      END SUBROUTINE READ_COM
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Read gaussian input file
!
      SUBROUTINE READ_ONIOM(STN)
!
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: STN 
!
 
  
      END SUBROUTINE READ_ONIOM
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Write oniom com 
!
      SUBROUTINE WRTMOLCOM(STN)
!
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: STN,MN,Mi,Mf,N,I,Ni,Nf,J,JMi,JMf,JN,JI,Nb
!
      ONMTAG(:,STN) = 'L'
      FTAG(3,:,STN) = 'F'

!
      MN = 138 
!     DO MN=1,MOLCNT(STN) 
        TITLE(STN) = 'mol_1'
        OUTCOM = 'mol_1.com'
        OUTXYZ = 'mol_1.xyz'
        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1  

        CALL CENTMOL(MN,STN)
        CALL molPBC(STN)          !Apply PBC's 
!        CALL WHOLEMOL(STN)

! debug 
       WRITE(68,*) ' Mol ',MN,' has ',Mf-Mi+1,' atoms'
        
        DO N=Mi,Mf
          I = MOLST(N,STN) 
          ONMTAG(I,STN) = 'H'
          FTAG(3,I,STN) = 'T'
       ENDDO
!
!         Loop over molecule neigbor list of atom 
!
        
          Ni = NLISTML(MN,STN)           
          Nf = NLISTML(MN+1,STN)-1
! debug 
       WRITE(68,*) Nf-Ni+1,' nb molecules'

          DO Nb=Ni,Nf
              J = NLISTML(Nb,STN)
              JMi = MPNT(J,STN)
              JMf = MPNT(J+1,STN)-1          

! debug 
       WRITE(68,*) 'NB mole',J,' has ',JMf-JMi+1,' atoms'
        

              DO JN=JMi,JMf
                JI = MOLST(JN,STN) 
                ONMTAG(JI,STN) = 'M'
                FTAG(3,JI,STN) = 'T'
              ENDDO
          ENDDO 
          CALL WRITE_COM(STN)
          CALL WRITE_XYZ(STN)
!     ENDDO   
      END SUBROUTINE WRTMOLCOM
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Write oniom com 
!
      SUBROUTINE MOLCOM(STN)
!
      USE specify
      USE structure
      USE potential 
!
      IMPLICIT none
!
      INTEGER :: STN,Mi,NQM,Mr,Mo,STNi,N,M,Mf,MN,D,I
      REAL*8  :: MS,RMAS(NDIM),TOTMr
      LOGICAL :: NEWM 
!
      STNi = STN - STOTo 

!     Set some local variables
!
      NQM = 0
      DO M = 1,NONMOL(STNi)
        Mi =  WONIONMOL(M,STNi) 
         NEWM = .TRUE. 
        IF( NQM.GT.0) THEN
            DO Mo=1,NQM
              Mr = QMLIST(Mo,STN)
              IF(Mi.EQ.Mr) NEWM = .FALSE.
            ENDDO
        ENDIF
        IF( NEWM ) THEN
            NQM = NQM + 1
            QMLIST(NQM,STN) = Mi 
        ENDIF 
      ENDDO 
!
!     Calcualte center of mass of specified molecules 
!
      TOTMr = 0.0d0
      RMAS(:) = 0.0d0
      DO M=1,NQM
        MN = QMLIST(M,STN)
        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1          
        DO N=Mi,Mf
          I = MOLST(N,STN) 
          MS = AMAS(I,STN ) 
          DO D =1,NDIM        
            RMAS(D) = RMAS(D) + MS*R0(D,I,STN)
          ENDDO
          TOTMr = TOTMr + MS 
        ENDDO
      ENDDO
!      IF( MOLCOMADDNB(STN) ) THEN
!
!      ENDIF 
!
!     Reset oniom sphere's to center of mass of cluster
!
      ONCENT(1,STNi) = RMAS(1)/TOTMr
      ONCENT(2,STNi) = RMAS(2)/TOTMr
      ONCENT(3,STNi) = RMAS(3)/TOTMr
!
      ONIONQM(STN) = NQM

!     Verbose output 
!
      IF(VERB) THEN
        WRITE(6,*) 'molecules found in region 1',NQM
        WRITE(6,*) 'the center of QM/MM sphere has been shifted to:'
        WRITE(6,*) ONCENT(1,STNi)
        WRITE(6,*) ONCENT(2,STNi)
        WRITE(6,*) ONCENT(3,STNi)
      ENDIF 
!
      END SUBROUTINE MOLCOM 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Write oniom com 
!
      SUBROUTINE SPCOM(STN)
!
      USE specify
      USE structure
!
      IMPLICIT none
!
      INTEGER :: STN,D,Mi,NQM,NMM,NFX,Mr,I,Mo,STNi
      REAL*8 :: LVr(NDIM,NDIM),RSQ,DR,RQM,RMM 
      LOGICAL :: NEWM 
      STNi = STN - STOTo 

!     Set some local variables
!
      RQM = ONRAD(1,STNi)
      RMM = ONRAD(2,STNi)
!
!     Intialize
!
      ONMTAG(:,STN) = 'L'
      FTAG(3,:,STN) = 'F'
      NQM = ONIONQM(STN) 
      NMM = 0
      NFX = 0 

!
!         Center structure around specified point
!
!
      DO I=1,NA(STN)
        RSQ = 0.0d0 
        DO D=1,NDIM
          R0(D,I,STN) = R0(D,I,STN) - ONCENT(D,STNi)
        ENDDO           
      ENDDO 
!
!     Apply PBC's 
!
      CALL molPBC(STN)        
!
!
!     Loop over all atoms and add molecules in radius to list 
!
      DO I=1,NA(STN)
         RSQ = 0.0d0 
         DO D=1,NDIM
          RSQ = RSQ + R0(D,I,STN)*R0(D,I,STN)
        ENDDO 
!
        DR = SQRT(RSQ)
        Mi = MOLN(I,STN) 
        NEWM = .TRUE.
        IF( DR.LT.RQM) THEN
          IF( NQM.GT.0) THEN
            DO Mo=1,NQM
              Mr = QMLIST(Mo,STN)
              IF(Mi.EQ.Mr) NEWM = .FALSE.
            ENDDO
          ENDIF
          IF( NEWM ) THEN
            NQM = NQM + 1
            QMLIST(NQM,STN) = Mi 
          ENDIF 
        ELSEIF (DR.LT.RMM) THEN
          IF( NMM.GT.0) THEN
            DO Mo=1,NMM
              Mr = MMLIST(Mo,STN)
              IF(Mi.EQ.Mr) NEWM = .FALSE.
            ENDDO
          ENDIF
          IF( NEWM ) THEN
            NMM = NMM + 1
            MMLIST(NMM,STN) = Mi 
          ENDIF 
        ELSE
          IF( NFX.GT.0) THEN
            DO Mo=1,NFX
              Mr = FXLIST(Mo,STN)
              IF(Mi.EQ.Mr) NEWM = .FALSE.
            ENDDO
          ENDIF
          IF( NEWM ) THEN
            NFX = NFX + 1
            FXLIST(NFX,STN) = Mi 
          ENDIF 
        ENDIF
      ENDDO      
!
!     Verbose output 
!
      IF(VERB) THEN
        WRITE(6,*) 'molecules found in region 1',NQM
        WRITE(6,*) 'molecules found in region 2',NMM
        WRITE(6,*) 'molecules found in region 3',NFX
      ENDIF 


!
      ONIONQM(STN) = NQM
      ONIONMM(STN) = NMM 
      ONIONFX(STN) = NFX 
!
      RETURN 
      END SUBROUTINE SPCOM 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/09/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Write gaussian oniom file
!
      SUBROUTINE WRITE_SPCOM(STN)
      USE specify
      USE structure
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: I,EL,STN,GFIX,NAr,M,MN,Mi,Mf,N 
      REAL*8 :: CG
      CHARACTER(IDSZ) :: AT,TP
!
      NAr = NA(STN) 
!      OPEN(UNIT=54,FILE=OUTSPCOM(STN),STATUS='unknown')
      OPEN(UNIT=55,FILE='sp.xyz',STATUS='unknown')
      WRITE(55,*) NA(STN)
      WRITE(55,*) 
!
      OPEN(UNIT=54,FILE='sp.com',STATUS='unknown')
      WRITE(54,5403) TITLE(STN)  
      WRITE(54,*) '%mem=10GB'
      WRITE(54,5400) 
      WRITE(54,*) ''
      WRITE(54,5402) TITLE(STN)
      WRITE(54,*) ''
      WRITE(54,*) '0 1 0 1 0 1 0 1 0 1 '
      DO M=1,ONIONQM(STN)
        MN = QMLIST(M,STN)
        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1          
        DO N=Mi,Mf
          I = MOLST(N,STN) 
          EL=ELN(I,STN)
          AT = ATSYM(EL)
          TP = TYP(I,STN)
          CG= ACHG(I,STN)
          GFIX = 0 
!          IF( FTAG(3,I,STN) .EQ. "F" ) 
          GFIX = 0
          ONMTAG(I,STN) = 'H'
          IF( CG.LT.0.0d0) THEN 
            WRITE(54,5401) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ELSE
            WRITE(54,5411) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ENDIF 
          WRITE(55,5501) R0(:,I,STN) 
        ENDDO 
      ENDDO
      DO M=1,ONIONMM(STN)
        MN = MMLIST(M,STN)
        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1          
        DO N=Mi,Mf
          I = MOLST(N,STN) 
          EL=ELN(I,STN)
          AT = ATSYM(EL)
          TP = TYP(I,STN)
          CG= ACHG(I,STN)
          GFIX = 0 
!          IF( FTAG(3,I,STN) .EQ. "F" ) GFIX = -1
          GFIX = 0
          ONMTAG(I,STN) = 'L'
          IF( CG.LT.0.0d0) THEN 
            WRITE(54,5401) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ELSE
            WRITE(54,5411) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ENDIF 
          WRITE(55,5502) R0(:,I,STN) 
        ENDDO 
      ENDDO
      IF( USONIOMFX(STN) ) THEN
       DO M=1,ONIONFX(STN)
        MN = FXLIST(M,STN)
        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1          
        DO N=Mi,Mf
          I = MOLST(N,STN) 
          EL=ELN(I,STN) 
          AT = ATSYM(EL)
          TP = TYP(I,STN)
          CG= ACHG(I,STN)
          GFIX = 0 
!          IF( FTAG(3,I,STN) .EQ. "F" ) GFIX = -1
          GFIX = -1
          ONMTAG(I,STN) = 'L'
          IF( CG.LT.0.0d0) THEN 
            WRITE(54,5401) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ELSE
            WRITE(54,5411) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ENDIF 
          WRITE(55,5503) R0(:,I,STN) 
        ENDDO 
       ENDDO
      ENDIF
      WRITE(54,*) ''
      WRITE(54,*) ''
!      WRITE(54,*) '--Link1--'
!      WRITE(54,*) '%chk=',TITLE
!      WRITE(54,*) '%mem=2GB'
!      WRITE(54,*) '# uff NoSymm OPT'
      WRITE(54,*) ''
      WRITE(54,*) ''
!
      CLOSE(54)
      CLOSE(55)

      RETURN
 5400 FORMAT('# oniom(b3lyp/6-31+g(d,p):hf/sto-3g:uff)=embedcharge 
     & nosymm  polar pop=(full,mk) sp guess=save ')
 5401 FORMAT(A2,"-",A2,"_R-",F6.3,I3,3F9.3," ",A5," ")
 5411 FORMAT(A2,"-",A2,"_R-",F6.3,I3,3F9.3," ",A5," ")
 5402 FORMAT(A20)
 5403 FORMAT('%chk=',A20)
 5501 FORMAT(" 50 ",3F9.3)
 5502 FORMAT(" 6 ",3F9.3)
 5503 FORMAT(" 1 ",3F9.3)
      END SUBROUTINE WRITE_SPCOM 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/09/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Write gaussian oniom file
!
      SUBROUTINE WRITE_ONIOM(STN,MN)
      USE specify
      USE structure
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: I,EL,STN,GFIX,MIp,MNp,MN,Mi,Mf,N 
     & ,NQM,NMM,NFX,NAM,LTH 
      REAL*8 :: CG
      CHARACTER(CHSZ) :: FQM,FMM,CORD,CONT,FRM
      CHARACTER(IDSZ) :: AT,TP
!
!
!     Write variables for file names 
!
      WRITE(FQM,565) MN
      WRITE(FMM,575) MN
      WRITE(CORD,585) MN
      WRITE(CONT,595) MN
!
!     Print QM 
!
      NQM =  ONIONQM(STN)
      NMM =  ONIONMM(STN)
      NFX =  ONIONFX(STN)
!
      Mi = MPNT(MN,STN)
      Mf = MPNT(MN+1,STN)-1          
      NAM = Mf - Mi + 1
!
!     Verbose
!
      IF( VERB ) THEN 
          WRITE(6,*) ' Writing files for molecule:',MN
          WRITE(6,*) '   gaussian file ',CORD
          IF(USONIOMFX(STN))WRITE(6,*)'Fixed atoms will be included' 
      ENDIF


!
!     Open xyz for qm region
!
      OPEN(UNIT=56,FILE=FQM,STATUS='UNKNOWN')
      WRITE(56,*) NAM 
!     add efiel and e potential to comment line
      WRITE(56,5603) MOLEF(:,MN,STN),MOLEP(MN,STN)
!     Open xyz for qm region
!
      OPEN(UNIT=57,FILE=FMM,STATUS='UNKNOWN')
      WRITE(57,*) NAM 
!     add efiel and e potential to comment line
      WRITE(57,5603) MOLEF(:,MN,STN),MOLEP(MN,STN)
!
!     Open gaussian oniom input file
!
      OPEN(UNIT=54,FILE=CORD,STATUS='UNKNOWN')
      WRITE(54,5451) MN
      WRITE(54,50) 
      WRITE(54,5452) 
!      WRITE(54,5400) 
!      WRITE(54,5440) ONIOM_METH
      WRITE(54,5461)
      WRITE(54,5462)
      WRITE(54,5463)


      WRITE(54,*) ''
      WRITE(54,5402) MN
      WRITE(54,*) ''
      IF( ONIOMLAYER .EQ. 2) THEN
        WRITE(54,5433) LAYER_CHR(1),LAYER_MLT(1)
     &  ,LAYER_CHR(2),LAYER_MLT(2),LAYER_CHR(1),LAYER_MLT(1)
      ELSE 
        WRITE(54,5433)  LAYER_CHR(1),LAYER_MLT(1)
     &  ,LAYER_CHR(2),LAYER_MLT(2),LAYER_CHR(3),LAYER_MLT(3)
      ENDIF 
!
!     Open gaussian single molecule input file
!
      OPEN(UNIT=55,FILE=CONT,STATUS='UNKNOWN')
      WRITE(55,5551) MN
      WRITE(55,50)
      WRITE(55,5452) 
      WRITE(55,5500) 
      WRITE(55,*) ''
      WRITE(55,5402) MN
      WRITE(55,*) ''
      WRITE(55,5431) LAYER_CHR(3),LAYER_MLT(3)
      DO MIp=1,NQM
        MNp = QMLIST(MIp,STN)
        Mi = MPNT(MNp,STN)
        Mf = MPNT(MNp+1,STN)-1          
        DO N=Mi,Mf
          I = MOLST(N,STN) 
          EL=ELN(I,STN)
          AT = ATSYM(EL)
          AT = TRIM(AT)
          TP = TYP(I,STN)
          CG= ACHG(I,STN)
          GFIX = 0 
!          IF( FTAG(3,I,STN) .EQ. "F" ) 
          GFIX = 0
          ONMTAG(I,STN) = 'H'
          LTH = LEN_TRIM(AT)
          IF(DEBUG ) THEN
            WRITE(101,*) CG,EL,LTH
            WRITE(101,*) AT
            WRITE(101,'(A1)') AT
          ENDIF
          
          IF( CG.LT.0.0d0  ) THEN
            WRITE(54,403) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ELSEIF(  CG.GE.0.0d0  ) THEN
            WRITE(54,405) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ELSE
            WRITE(*,*) ' ERORROR ERORROR ERORROR '
            STOP 
          ENDIF


          WRITE(55,5511) AT,R0(:,I,STN)
          WRITE(56,5501) R0(:,I,STN) 
        ENDDO 
        IF ( DEBUG ) STOP 
      ENDDO
      DO MIp=1,ONIONMM(STN)
        MNp = MMLIST(MIp,STN)
        Mi = MPNT(MNp,STN)
        Mf = MPNT(MNp+1,STN)-1          
        DO N=Mi,Mf
          I = MOLST(N,STN) 
          EL=ELN(I,STN)
          AT = ATSYM(EL)
          TP = TYP(I,STN)
          CG= ACHG(I,STN)
          GFIX = 0 
!          IF( FTAG(3,I,STN) .EQ. "F" ) GFIX = -1
          GFIX = 0
          ONMTAG(I,STN) = 'L'
          LTH = LEN_TRIM(AT)

          IF( CG.LT.0.0d0  ) THEN
            WRITE(54,403) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ELSEIF(  CG.GE.0.0d0  ) THEN
            WRITE(54,405) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ELSE
            WRITE(*,*) ' ERORROR ERORROR ERORROR '
            STOP 
          ENDIF

          WRITE(57,5502) R0(:,I,STN) 
        ENDDO 
      ENDDO
      IF( USONIOMFX(STN) ) THEN
       DO MIp=1,ONIONFX(STN)
        MNp = FXLIST(MIp,STN)
        Mi = MPNT(MNp,STN)
        Mf = MPNT(MNp+1,STN)-1          
        DO N=Mi,Mf
          I = MOLST(N,STN) 
          EL=ELN(I,STN) 
          AT = ATSYM(EL)
          TP = TYP(I,STN)
          CG= ACHG(I,STN)
          GFIX = 0
!          IF( FTAG(3,I,STN) .EQ. "F" ) GFIX = -1
          GFIX = -1
          ONMTAG(I,STN) = 'L'
          LTH = LEN_TRIM(AT)
          IF( CG.LT.0.0d0  ) THEN
            WRITE(54,403) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ELSEIF(  CG.GE.0.0d0  ) THEN
            WRITE(54,405) AT,AT,CG,GFIX,R0(:,I,STN),ONMTAG(I,STN)
          ELSE
            WRITE(*,*) ' ERORROR ERORROR ERORROR '
            STOP 
          ENDIF

          WRITE(57,5503) R0(:,I,STN) 
        ENDDO 
       ENDDO
      ENDIF
      WRITE(54,*) ''
      WRITE(54,*) ''
      WRITE(55,*) ''
      WRITE(55,*) ''


      CLOSE(54)
      CLOSE(55)
      CLOSE(56)
      CLOSE(57)

      RETURN
 
 565  FORMAT(I8,'_qm.xyz')
 575  FORMAT(I8,'_mm.xyz')
 585  FORMAT(I8,'_oniom.com')
 595  FORMAT(I8,'_sm.com')
 5451 FORMAT('%chk=',I3,'_CHR_oniom.chk')
 5551 FORMAT('%chk=',I3,'_CHR_sm.chk')
 5452 FORMAT('%nproc=4')
! 5400 FORMAT('# oniom(b3lyp/6-31+g(d,p):uff)=embedcharge 
!     & nosymm  polar pop=(full,mk) sp guess=save ')
 5400 FORMAT('# oniom(HF/6-31G*:uff)=embedcharge 
     & nosymm  polar pop=(full,mk) sp guess=save ')
 5500 FORMAT('# HF/6-31G* nosymm  polar pop=(full,mk) sp guess=save ')

 401  FORMAT(A2,"-",A1,"_-",F9.6,I3,3F12.6," ",A5," ")
 4011 FORMAT(A2,"-",A1,"_-",F8.6,I3,3F12.6," ",A5," ")
 402  FORMAT(A2,"-",A1,"_R-",F9.6,I3,3F12.6," ",A5," ")
 403  FORMAT(A2,"-",A2,"_R-",F9.6,I3,3F12.6," ",A5," ")
 404  FORMAT(A2,"-",A1,"_R-",F8.6,I3,3F12.6," ",A5," ")
 405  FORMAT(A2,"-",A2,"_R-",F8.6,I3,3F12.6," ",A5," ")

 5411 FORMAT(A2,"-",A2,"_R-",F9.6,I3,3F12.6," ",A5," ")

 5511 FORMAT(A2," ",3F9.3)
 5501 FORMAT(" 50 ",3F9.3)
 5502 FORMAT(" 6 ",3F9.3)
 5503 FORMAT(" 1 ",3F9.3)
 5603 FORMAT("Efield and V",4E16.4)
 5402 FORMAT("  qm/mm geometry for molecule:",I6)
 5431 FORMAT(2I4)
 5432 FORMAT(4I4)
 5433 FORMAT(6I4)
 5440 FORMAT(A100)
 5461 FORMAT("#p oniom(uwb97xd/6-31+g(d,p):uff)=embedcharge nosymm")
 5462 FORMAT("iop(6/7=3,3/107=0232000000,3/108=0232000000) ")
 5463 FORMAT("scf=(maxcycle=200,tight) sp ")
  50  FORMAT('%mem=2GB')
      END SUBROUTINE WRITE_ONIOM
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
