!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Determine molecules
      SUBROUTINE molcheck(STN)
      USE structure
      USE specify
      USE elements 
!
      IMPLICIT none 
!

      INTEGER :: STN,NAr
     & ,NMMAX,MN,NBr
      INTEGER :: CNT,MOLCNTr,EL,I,N,MI,MF,CREF,MJ,Ni,Nf,A,B,Elb
      LOGICAL :: TERM
      LOGICAL, ALLOCATABLE :: MOL(:)
      CHARACTER(CHSZ) :: ESTAT

!     Set max molecules to # of atoms for now
      NAr= NA(STN)
      NMMAX = NAr
      NBr = NBTOT(STN)-1
!     Set max molecules to # of atoms for now
      NMMAX = NAr
!     Initialize 
      ALLOCATE( MOL(NAr) ) 
      MOL(:) = .FALSE.
      CNT     = 0
      MOLCNTr = 0

!
!     Verbose
!
      IF( VERB ) THEN
         WRITE(6,*)'Starting molecule check of ',STN,' of ',NAr,' atoms'
      ENDIF
!
!
!     Loop over all atoms
      DO I=1,NAr
        EL = ELN(I,STN)
!       If not part of a molecule and not a hydrogen
        IF( .NOT.MOL(I).AND.EL.NE.1) THEN
          MOLCNTr = MOLCNTr + 1
          CNT = CNT  + 1
          IF( CNT .GT. NMMAX ) CALL prnerror(-4,ESTAT)
          MOLST(CNT,STN) = I
          MPNT(MOLCNTr,STN) = CNT
          MOLN(I,STN) = MOLCNTr
          MOL(I) = .TRUE.
!         While termination is not present add non hydrogen neighbors to 
          TERM = .TRUE.
          DO WHILE ( TERM ) 
            MI = MPNT(MOLCNTr,STN)
            MF = CNT
            CREF = CNT 
!           Loop over all atoms in molecule (A)
            DO MJ = MI,MF
              A = MOLST(MJ,STN)
              Ni = NINDX(A,STN) 
              Nf = NINDX(A+1,STN)-1
!             Loop over all neighbors (B) of atom A
              IF(NF.GE.Ni) THEN
               DO  N = Ni,Nf
                B = NLIST(N,STN)
                ELb = ELN(B,STN)
!               If not already part of molecule and not a hydrogen
                IF( .NOT.MOL(B).AND.Elb.NE.1 ) THEN
                  CNT = CNT  + 1
                  MOL(B) = .TRUE.
                  MOLST(CNT,STN) = B
                  MOLN(B,STN) = MOLCNTr       
                ENDIF
               ENDDO
              ENDIF
            ENDDO
!           If no atoms are added to the molecule stop
            IF( CNT.EQ.CREF) TERM = .FALSE.             
          ENDDO
!         Add hydorgens
          MI = MPNT(MOLCNTr,STN)
          MF = CNT
          DO MJ = MI,MF
            A = MOLST(MJ,STN)
            Ni = NINDX(A,STN) 
            Nf = NINDX(A+1,STN)-1
!           Loop over al l neighbors (B) of atom A
            IF(NF.GT.Ni) THEN
             DO  N = Ni,Nf
              B = NLIST(N,STN)
              ELb = ELN(B,STN)
!             If not already part of molecule and a hydrogen
              IF( .NOT.MOL(B).AND.Elb.EQ.1 ) THEN
                CNT = CNT  + 1
                MOL(B) = .TRUE.
                MOLST(CNT,STN) = B
                MOLN(B,STN) = MOLCNTr                  
              ENDIF
             ENDDO
            ENDIF
          ENDDO
        ENDIF     
      ENDDO
!     Find any H2
      DO  I = 1,NAr
!       If not part of molecule (H2), hopefully
        IF( .NOT.MOL(I) ) THEN
          MOLCNTr = MOLCNTr + 1
          CNT = CNT + 1
          IF( CNT.GT.NAr) CALL prnerror(-4,ESTAT)
          MOLST(CNT,STN) = I
          MPNT(MOLCNTr,STN) = CNT
          MOLN(I,STN)  = MOLCNTr
          MOL(I) = .TRUE.
          Ni = NINDX(I,STN) 
          Nf = NINDX(I+1,STN)-1
!         Loop over all neighbors (B) of atom I
          IF(NF.GT.Ni) THEN
           DO  N = Ni,Nf
            B = NLIST(N,STN)
            ELb = ELN(B,STN)
            IF( Elb .NE. 1 ) CALL prnerror(-4,ESTAT)
            CNT = CNT + 1
            MOLST(CNT,STN) = B
            MOLN(B,STN) = MOLCNTr
            MOL(B) = .TRUE.
           ENDDO
          ENDIF
        ENDIF
      ENDDO
!     Set end of molecule list
      MPNT(MOLCNTr + 1,STN ) = CNT + 1
      MOLCNT(STN) = MOLCNTr

!
      IF( VERB ) WRITE(*,*) MOLCNTr," molecules found in ",STN

      DEALLOCATE(MOL)
!     Sort list of atoms
      CALL SORTMOLLIST(STN)
!
      RETURN
      END SUBROUTINE molcheck 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Determine molecules
      SUBROUTINE masspec(STN)
      USE structure
      USE specify

      IMPLICIT none 
!

      INTEGER :: STN,MN,Mi,Mf
!

      OPEN(UNIT=60,FILE='massspec.dat',STATUS='unknown')
      WRITE(60,*) "# Mass spec:",STN
      DO MN=1,MOLCNT(STN)
         Mi = MPNT(MN,STN)
         Mf = MPNT(MN+1,STN)-1
         WRITE(60,601) MN,Mf-Mi+1,MCOMS(:,MN,STN), MOLCHR(MN,STN)
      ENDDO
      CLOSE(60)
!
      RETURN
 601  FORMAT(I6,I8,4F12.2)
      END SUBROUTINE masspec
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Calculate center of mass
      SUBROUTINE centmass(STN)
!
      USE const
      USE structure
      USE specify 
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: STN,N,I,EL,MN,Mi,Mf,D
      REAL*8 :: TMASr,MAS,TMA(NDIM),TCHRr,CHRr
      CHARACTER(ERSZ) :: ESTAT
!
!     Open 
      IF( VERB) WRITE(6,*) "Starting center of mass calc"
!     Loop over all molecules
      DO MN=1,MOLCNT(STN)
        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1          
        TMASr = 0.0d0
        TCHRr  = 0.0d0
        TMA(:) = 0.0d0
        IF(DEBUG) WRITE(203,*) MN,Mi,Mf,Mf-Mi+1
        DO N=Mi,Mf
           I = MOLST(N,STN) 
           EL = ELN(I,STN)
           MAS = AMASS(EL)
           CHRr = ACHG(I,STN)
           TMASr = TMASr + MAS
           TCHRr = TCHRr + CHRr
           DO D =1,NDIM
             TMA(D) = TMA(D) + R0(D,I,STN)*MAS
           ENDDO           
         ENDDO
         DO D =1,NDIM
           MCOMS(D,MN,STN) = TMA(D)/TMASr
         ENDDO
         MOLCHR(MN,STN) = TCHRr 
         MOLMAS(MN,STN) = TMASr
         IF ( VERB ) THEN
          WRITE(6,601) MN,MOLMAS(MN,STN),MCOMS(:,MN,STN),MOLCHR(MN,STN)
         ENDIF
      ENDDO          
!
      RETURN
 601  FORMAT(I6,F20.6,4F20.6) 
      END SUBROUTINE centmass
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Find molecule neighbors 
      SUBROUTINE BUILDMNB(STN)
      USE const
      USE structure
      USE specify 
!
      IMPLICIT none
!
      INTEGER :: STN,MN,Mi,Mf,N,MNB,MNBo,Ni,Nf,Nb,I,J
     & ,MNn,MICK,MCK 
      LOGICAL :: ADDM 
      CHARACTER(ERSZ) :: ESTAT
!
!     Set initial 
!
      MNB = 0 
!
      IF(VERB) THEN
         WRITE(6,*) 'Creating mol nb list for ',STN,' with '
     &  ,MOLCNT(STN),' molecules'
      ENDIF
!
!     Loop over all molecules
!
      DO MN=1,MOLCNT(STN)
        MNBo = MNB + 1
        NINDXML(MN,STN) = MNBo
        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1   
!
!       Loop over all atoms in molecule 
!       
        DO N=Mi,Mf
          I = MOLST(N,STN) 
!
!         Loop over LJ neigbor list of atom 
!
          Ni = NINDXLJ(I,STN)           
          Nf = NINDXLJ(I+1,STN)-1
!
          DO Nb=Ni,Nf
              J = NLISTLJ(Nb,STN)
!
!             If LJ neighbor not molecule neighbor list or same molecule 
!
              MNn = MOLN(J,STN)
              ADDM = .TRUE.
              IF ( MNn.EQ.MN ) THEN
                 ADDM = .FALSE.
              ELSEIF( MNB.GE.MNBo) THEN
                DO MICK = MNBo,MNB 
                  MCK = NLISTML(MICK,STN)
                  IF ( MCK.EQ.MNn ) ADDM = .FALSE.
                ENDDO
              ENDIF
!
!             add to molecule neighbor list 
!
              IF( ADDM ) THEN
                MNB = MNB + 1 
                NLISTML(MNB,STN) = MNn
              ENDIF
          ENDDO  !          Loop over LJ neigbor list of atom 
        ENDDO !       Loop over all atoms in molecule 
        IF( VERB) WRITE(6,*) 'Mol:',MN,' has ',MNB-MNBo+1,' mol nb' 
      ENDDO !     Loop over all molecules
      MNBo = MNB + 1
      NINDXML(MN,STN) = MNBo


      RETURN
 601  FORMAT(I6,F20.6,3F20.6) 
      END SUBROUTINE BUILDMNB 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Find molecule neighbors 
      SUBROUTINE SORTMOLLIST(STN)
!
      USE const
      USE structure
      USE specify 
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: STN,N,I,MN,Mi,Mf,NAr,IND
      INTEGER, ALLOCATABLE :: LISTr(:),ORD(:)

      
!     Loop over all molecules
      DO MN=1,MOLCNT(STN)
        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1          
        ! 
        !  Create a temp array of atom #'s for molecule
        !
        NAr = Mf - Mi + 1
        ALLOCATE ( LISTr(NAr),ORD(NAr) ) 
        IND = 0
        DO N=Mi,Mf
           IND = IND + 1
           LISTr(IND) = MOLST(N,STN) 
        ENDDO ! atom # N
        ! Sort temp array
        CALL QSORTI( ORD,NAr,LISTr )
        ! 
        ! Update global list 
        !  
        I = 0
        DO N=Mi,Mf
           I = I + 1
           IND = ORD(I)
           MOLST(N,STN) = LISTr(IND)

           IF( DEBUG ) THEN
              WRITE(105,*) I,LISTr(I),LISTr(IND)
           ENDIF 

        ENDDO ! atom # N
        DEALLOCATE ( LISTr,ORD ) 

      ENDDO   ! molecule # MN 


      RETURN 
      END SUBROUTINE SORTMOLLIST 

