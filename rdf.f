!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Calculate distance between residues 
!
      SUBROUTINE prdist(STN)
      USE specify
      USE structure 
      USE build 
      USE elements 
!     
      IMPLICIT none
      INTEGER :: STN,I,J,L,NPNTS,D,NAr,STNo,EL
     & ,SZ,DNB,RCNT,F,FID,ELi,ELj
      INTEGER, ALLOCATABLE:: RDDST(:,:),NDEN_SUM(:)
      REAL*8 :: Ra(NDIM),DR,CUTSQ,DELR,Rin,Rout,DVi,DN,NDEN,GRN,GRM
     & ,RSQ,LVr(NDIM,NDIM) ,DM
     & ,RDSCALE,RDCUT
     & ,Fi(NDIM),Fj(NDIM),DF(NDIM),DRr(NDIM)
      REAL*8 , ALLOCATABLE :: SPV(:),EX(:) ,NDENS(:),MDENS(:),VCUT(:)
     & ,MDEN_SUM(:),RDMDST(:,:)
      CHARACTER(CHSZ) :: Ai,Bi,Ao,Bo
      CHARACTER(CHSZ), ALLOCATABLE :: RFID1(:,:),RFID2(:,:),RDISr(:,:)
      LOGICAL :: UDIM(NDIM)
!     
!     Set the local variables
      STNo = STN !- STOTo 
      NAr = NA(STN)
      RCNT = PRCNT(STNo) 
      FID  = 40                      ! Initial file number 
      RDSCALE = 10.0d0
      RDCUT   = 25.0d0
      CUTSQ = RDCUT*RDCUT
      SZ = ANINT(RDCUT*RDSCALE)
      ALLOCATE( RDDST(0:SZ,RCNT), RDMDST(0:SZ,RCNT)
     & ,SPV(RCNT),EX(RCNT)
     & ,RFID1(NAr,RCNT),RFID2(NAr,RCNT),RDISr(2,RCNT) 
     & ,NDENS(RCNT),MDENS(RCNT),VCUT(RCNT) 
     & ,NDEN_SUM(RCNT),MDEN_SUM(RCNT) )
!     
      RDDST(:,:) = 0
      DR = 1.0d0/RDSCALE
!     
!
!     Triclinc PBC's
      LVr(:,:) = LV(:,:,STN)
!     Resize lattice vector in vacuum direction 
      IF( MKVAC(1) ) THEN
        LVr(SDIM,SDIM)=MMN(2,SDIM,STN)-MMN(1,SDIM,STN)+MINBF(STN)*2
      ENDIF
!     Calculate fractional coordinates 
      ALLOCATE( Fr(NDIM,NAr) ) 
      CALL fracR(STN,LVr,1)
!
!     Verbose output
      IF (VERB ) THEN
        WRITE(*,*) "Starting atomic rdf"
        WRITE(*,*) "  Structure ",STN," with ", RCNT," rdf's to print"
        WRITE(*,*) 'Using latice vectors '
        WRITE(*,*) LVr(1,:)
        WRITE(*,*) LVr(2,:)
        WRITE(*,*) LVr(3,:)
      ENDIF
!
!     Open ouput file
!
      DO DNB = 1,RCNT 
        F = FID + DNB 
!       Determine ID specification 
        RDISr(1,DNB) = RDIS(1,DNB,STNo)
        RDISr(2,DNB) = RDIS(2,DNB,STNo)
        DO I=1,NA(STN)
          IF( RDISr(1,DNB).EQ.TYP(I,STN) ) THEN
            RFID1(:,DNB) =  TYP(:,STN) ! Force field type
            IF(VERB)WRITE(6,*)' Using Force field type',TYP(I,STN)
     &       ,' for ID 1'
          ENDIF
          IF( RDISr(1,DNB).EQ.RID(I,STN) ) THEN
            RFID1(:,DNB) =  RID(:,STN) ! Residue ID
            IF(VERB)WRITE(6,*)' Using Residue ID',RID(I,STN)
     &       ,' for ID 1'
          ENDIF
          IF( RDISr(1,DNB).EQ.GRID(I,STN)) THEN
            RFID1(:,DNB) = GRID(:,STN) ! GROMACS ID
            IF(VERB)WRITE(6,*)' Using GROMACS ID ',GRID(I,STN)
     &       ,' for ID 1'
          ENDIF
          EL =  ELN(I,STN) 
          IF( RDISr(1,DNB).EQ.ATSYM(EL)) THEN
            RFID1(I,DNB) = ATSYM(EL)  ! Using atomic sym
            IF(VERB)WRITE(6,*)' Using atomic symbol  ',ATSYM(EL)
     &       ,' for ID 1'
          ENDIF
!
          IF( RDISr(2,DNB).EQ.TYP(I,STN) ) THEN
            RFID2(:,DNB) =  TYP(:,STN) ! Force field type
            IF(VERB)WRITE(6,*)' Using Force field type',TYP(I,STN)
     &       ,' for ID 2'
          ENDIF
          IF( RDISr(2,DNB).EQ.RID(I,STN) ) THEN
            RFID2(:,DNB) =  RID(:,STN) ! Residue ID
            IF(VERB)WRITE(6,*)' Using Residue ID',RID(I,STN)
     &       ,' for ID 2'
          ENDIF
          IF( RDISr(2,DNB).EQ.GRID(I,STN)) THEN
            RFID2(:,DNB) = GRID(:,STN) ! GROMACS ID
            IF(VERB)WRITE(6,*)' Using GROMACS ID ',GRID(I,STN)
     &       ,' for ID 2'
          ENDIF
          EL =  ELN(I,STN) 
          IF( RDISr(2,DNB).EQ.ATSYM(EL)) THEN
            RFID2(I,DNB) = ATSYM(EL)  ! Using atomic sym
            IF(VERB)WRITE(6,*)' Using atomic symbol  ',ATSYM(EL)
     &       ,' for ID 1'
          ENDIF
!
        ENDDO  
        Ao = RDISr(1,DNB)
        Bo = RDISr(2,DNB)
        UDIM(:) = UDIS(:,DNB,STNo)
        OPEN(UNIT=F,FILE=PRDISFL(DNB,STNo),STATUS='UNKNOWN')
        WRITE(F,*) "# radial distribution function"
        WRITE(F,*) "# Scale ", RDSCALE
        WRITE(F,*) "# cut off " , RDCUT
        WRITE(F,3500) Ao,Bo
      
!       Find dimensionality of distribution
        IF( UDIM(1).AND.UDIM(2).AND.UDIM(3)) THEN
!         3D spherical 
          IF(VERB) WRITE(6,*) " 3D spherical distibution "
          WRITE(F,*) "# 3D spherical distibution "
          SPV(DNB)=4.0d0*PI/3.0d0
          EX(DNB) = 3
        ELSEIF(  UDIM(1).AND.UDIM(2).AND. .NOT.UDIM(3) ) THEN
!         2D circular
          IF(VERB) WRITE(6,*) " 2D cicular distibution "
          WRITE(F,*) "# 2D cicular distibution "
          SPV(DNB) = 2.0d0*PI
          EX(DNB) = 2
        ELSEIF(  UDIM(1).AND. .NOT.UDIM(2).AND. UDIM(3) ) THEN
!         2D circular
          IF(VERB) WRITE(6,*) " 2D cicular distibution "
          WRITE(F,*) "# 2D cicular distibution "
          SPV(DNB) = 2.0d0*PI
          EX(DNB) = 2
        ELSEIF(   .NOT.UDIM(1).AND.UDIM(2).AND.UDIM(3) ) THEN
!         2D circular
          IF(VERB) WRITE(6,*) " 2D cicular distibution "
          WRITE(F,*) "# 2D cicular distibution "
          SPV(DNB) = 2.0d0*PI
          EX(DNB) = 2
        ELSE 
  !          1D
          IF(VERB) WRITE(6,*) " 1D planer distibution "
          WRITE(F,*) "# 1D  planer distibution "
          SPV(DNB) = 1.0d0
          DO D=1,NDIM
            IF( .NOT.UDIM(D) ) SPV(DNB)= SPV(DNB)*LVr(D,D)
          ENDDO
          EX(DNB)=1
        ENDIF
        WRITE(F,*) "#  with an area of ",SPV(DNB),"*dr"
!
!      Calculate volume 
!
        VCUT(DNB)  = SPV(DNB)*RDCUT**EX(DNB)
      ENDDO 
!     
!     Zero sums 
!
      NDEN_SUM(:) =  0 
      MDEN_SUM(:) =  0.0d0 
      RDMDST(:,:) = 0.0d0
      RDDST(:,:) = 0

!
!     Loop over all atomic pairs 
!

      DO I=1,NA(STN)-1
        Ra(:) = R0(:,I,STN) 
        Fi(:) = Fr(:,I)
        ELi = ELN(I,STN) 
        DO DNB = 1,RCNT 
           Ai = RFID1(I,DNB)  
           Ao = RDISr(1,DNB)    ! RFID(1,DNB)
           Bo = RDISr(2,DNB)    ! RFID(2,DNB)
           IF( (Ai.EQ.Ao) .OR.( Ai.EQ.Bo) )THEN
              NDEN_SUM(DNB) =  NDEN_SUM(DNB) + 1
              MDEN_SUM(DNB) =  MDEN_SUM(DNB) + AMASS(ELi) 
           ENDIF
        ENDDO 
        DO J=I+1,NA(STN)
          Fj(:) = Fr(:,J)
!         Loop over rdf's
          DO DNB= 1,RCNT 
!           Check for intra molecular interactions
            IF( .NOT.(RMIAM(STNo).AND.MOLN(I,STN).EQ.MOLN(J,STN))) THEN 
             Ai = RFID1(I,DNB)  
             Bi = RFID2(J,DNB)  
             Ao = RDISr(1,DNB) ! RFID(1,DNB)
             Bo = RDISr(2,DNB) ! RFID(2,DNB)
             IF((Ai.EQ.Ao.AND.Bi.EQ.Bo).OR.(Bi.EQ.Ao.AND.Ai.EQ.Bo))THEN
                UDIM(:) = UDIS(:,DNB,STNo)
                RSQ = 0.0d0
                DO D = 1,NDIM
                  IF( UDIM(D)) THEN
                    DF(D) = Fj(D) - Fi(D)
                    IF( PBCS(D,STN) )  DF(D) = DF(D) - ANINT(DF(D))
                  ENDIF
                ENDDO
!
!               Calculate cartion coordinates 
!
!                CALL REALF(LVr,DF,DRr)
                DRr(1) = LVr(1,1)*DF(1) + LVr(2,1)*DF(2)+LVr(3,1)*DF(3)
                DRr(2) = LVr(1,2)*DF(1) + LVr(2,2)*DF(2)+LVr(3,2)*DF(3)
                DRr(3) = LVr(1,3)*DF(1) + LVr(2,3)*DF(2)+LVr(3,3)*DF(3)
                RSQ = 0.0d0
                DO D=1,NDIM
                   RSQ = RSQ + DRr(D)*DRr(D)
                ENDDO
                IF( RSQ.LE.CUTSQ ) THEN
                   ELj = ELN(J,STN) 
                   DELR=SQRT(RSQ)
                   L=INT(DELR*RDSCALE)
                   RDDST(L,DNB) = RDDST(L,DNB) + 2 
                   RDMDST(L,DNB)=RDMDST(L,DNB)+AMASS(ELj)+AMASS(ELi) 
!
                   NDEN_SUM(DNB) = NDEN_SUM(DNB) + 2 
                   MDEN_SUM(DNB)=MDEN_SUM(DNB)+AMASS(ELj)+AMASS(ELi) 
!
                   IF(VERB)THEN
                    WRITE(F,3502)I,J,DELR,MOLN(I,STN),MOLN(J,STN) 
                   ENDIF 
                ENDIF 
              ENDIF ! (Ai.EQ.Ao.AND.Bi.EQ.Bo).OR.(Bi.EQ.Ao.AND.Ai.EQ.Bo)
            ENDIF   !  RMIAM(STNo) .AND. MOLN(I,STN).NE. MOLN(J,STN)
          ENDDO     ! DNB 
        ENDDO       ! J 
      ENDDO         ! I

!     Print distribution 
!     Loop over rdf's
      DO DNB= 1,RCNT 
        F = FID + DNB 
        NPNTS = NDEN_SUM(DNB)
!       Use density of sphere formormalization 
        NDENS(DNB) = FLOAT(NPNTS)/VCUT(DNB)
        MDENS(DNB) = MDEN_SUM(DNB)/VCUT(DNB)
!       Verbose info
        WRITE(F,*) 'For rdf ',DNB
        WRITE(F,*) " # total nieghbors found", NPNTS
        WRITE(F,*) ' per unit of ',VCUT(DNB),' A ',EX(DNB)
        WRITE(F,*) ' number density is ', NDEN_SUM(DNB),NDENS(DNB)
        WRITE(F,*) ' mass density is ',MDENS(DNB),MDENS/AVO*10.0d0
!     
        DO L=0,SZ-1
          Rin=FLOAT(L)/RDSCALE
          Rout = Rin+DR
          DVi = SPV(DNB)*(Rout**EX(DNB) - Rin**EX(DNB))
          DN = FLOAT(RDDST(L,DNB))
          DM = RDMDST(L,DNB)
          GRN = DN/DVi/NDENS(DNB)
          GRM = DM/DVi/MDENS(DNB)
          WRITE(F,3501) Rin,Rout,DVi,DN,GRN,DM,GRM
        ENDDO
        CLOSE(F)
      ENDDO ! DNB 
!
      DEALLOCATE( Fr ) 
      DEALLOCATE( RDDST,RDMDST,SPV,EX,RFID1,RFID2,RDISr
     & ,NDENS,MDENS,VCUT,NDEN_SUM,MDEN_SUM )

!
      RETURN
!
 3500 FORMAT("# Atom names: ",2A6)
 3501 FORMAT(7F16.3)
 3502 FORMAT("#",2I12,F6.2,2I6)
!
      END SUBROUTINE PRDIST
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Calculate distance between points  
!
      SUBROUTINE MOLRDF(STN)
      USE specify
      USE structure 
      USE build 
!     
      IMPLICIT none
!
      INTEGER :: STN,SZ,D,L,NPNTS,STNo,RCNT,I,J
     & ,MNi,MIi,MJi,Mf,MNAi,MNAj,MNj,F,FID,DNB
      REAL*8 :: RDSCALE,RDCUT,CUTSQ,DR,Rin,Rout,DVi,DN,GRN,DM,GRM
     & ,LVr(NDIM,NDIM),Ri(NDIM),Rj(NDIM),AVESM,MASSi
     & ,RSQ,DF(NDIM),DRr(NDIM) ,DELR,NDEN,Fi(NDIM),Fj(NDIM)
      REAL*8 , ALLOCATABLE :: SPV(:),EX(:),RDMDST(:,:)
     & ,NDENS(:),MDENS(:),VCUT(:)
     & ,MDEN_SUM(:)
      INTEGER, ALLOCATABLE:: RDDST(:,:),NDEN_SUM(:)
      LOGICAL :: UDIM(NDIM)
!
      STNo = STN           !- STOTo
      RCNT = NMCNT(STNo) 
      RDSCALE = 1.0d0
      RDCUT   = 25.0d0
      CUTSQ = RDCUT*RDCUT
      SZ = ANINT(RDCUT*RDSCALE)
      ALLOCATE( RDMDST(0:SZ,RCNT),RDDST(0:SZ,RCNT),SPV(RCNT),EX(RCNT)
     & ,NDENS(RCNT),MDENS(RCNT),VCUT(RCNT)
     & ,NDEN_SUM(RCNT),MDEN_SUM(RCNT) )
!     
!     Set the local variables
      FID  = 40                      ! Initial file number 
      DR = 1.0d0/RDSCALE
!
!     Triclinc PBC's
      LVr(:,:) = LV(:,:,STN)
!     Resize lattice vector in vacuum direction 
      IF( MKVAC(1) ) THEN
        LVr(SDIM,SDIM)=MMN(2,SDIM,STN)-MMN(1,SDIM,STN)+MINBF(STN)*2
      ENDIF
!      ALLOCATE( Fr(NDIM,NAr) ) 
!      CALL fracR(STN,LVr,1)
!     Verbose output
      IF (VERB ) THEN
        WRITE(*,*) "Starting molecule rdf"
        WRITE(*,*)"  Structure ",STN," with ", MOLCNT(STN)," molecules"
        WRITE(*,*)"  and ",RCNT," rdf's to print"
      ENDIF
!
!     Open ouput file
      DO DNB = 1,RCNT 
        F = FID + DNB 
        OPEN(UNIT=F,FILE=MPRDISFL(DNB,STNo),STATUS='UNKNOWN')
        WRITE(F,*) "# Pair distribution function of mol cent mass"
        WRITE(F,*) "# Scale ", RDSCALE
        WRITE(F,*) "# cut off " , RDCUT
!       Find dimensionality of distribution
        UDIM(:) = UMDIS(:,DNB,STNo)
        IF( UDIM(1).AND.UDIM(2).AND.UDIM(3)) THEN
!         3D spherical 
          IF(VERB) WRITE(6,*) " 3D spherical distibution "
          WRITE(F,*) "# 3D spherical distibution "
          SPV(DNB)=4.0d0*PI/3.0d0
          EX(DNB) = 3
        ELSEIF(  UDIM(1).AND.UDIM(2).AND. .NOT.UDIM(3) ) THEN
!         2D circular
          IF(VERB) WRITE(6,*) " 2D cicular distibution "
          WRITE(F,*) "# 2D cicular distibution "
          SPV(DNB) = 2.0d0*PI
          EX(DNB) = 2
        ELSEIF(  UDIM(1).AND. .NOT.UDIM(2).AND. UDIM(3) ) THEN
!         2D circular
          IF(VERB) WRITE(6,*) " 2D cicular distibution "
          WRITE(F,*) "# 2D cicular distibution "
          SPV(DNB) = 2.0d0*PI
          EX(DNB) = 2
        ELSEIF(   .NOT.UDIM(1).AND.UDIM(2).AND.UDIM(3) ) THEN
!         2D circular
          IF(VERB) WRITE(6,*) " 2D cicular distibution "
          WRITE(F,*) "# 2D cicular distibution "
          SPV(DNB) = 2.0d0*PI
          EX(DNB) = 2
        ELSE 
  !     1D
          IF(VERB) WRITE(6,*) " 1D planer distibution "
          WRITE(F,*) "# 1D  planer distibution "
          SPV(DNB) = 1.0d0
          DO D=1,NDIM
            IF( .NOT.UDIM(D) ) SPV(DNB)= SPV(DNB)*LVr(D,D)
          ENDDO
          EX(DNB)=1
        ENDIF
        WRITE(F,*) "#  with an area of ",SPV(DNB),"*dr"
!
!      Calculate volume 
!
        VCUT(DNB)  = SPV(DNB)*RDCUT**EX(DNB)
      ENDDO 
!     
!     Zero sums 
!
      NDEN_SUM(:) =  0 
      MDEN_SUM(:) =  0.0d0 
      RDMDST(:,:) = 0.0d0
      RDDST(:,:) = 0
!
!     Loop over all molecules 
!
      DO MNi=1,MOLCNT(STN)- 1
        MIi = MPNT(MNi,STN)
        Mf = MPNT(MNi+1,STN)-1          
        MNAi = Mf-MIi+1
        MASSi = MOLMAS(MNi,STN)
        DO DNB =1,RCNT 
           IF( MNAi.EQ.NMOLRi(DNB,STNo).OR.MNAi.EQ.NMOLRj(DNB,STNo))THEN
              NDEN_SUM(DNB) =  NDEN_SUM(DNB) + 1
              MDEN_SUM(DNB) =  MDEN_SUM(DNB) + MASSi 
           ENDIF
        ENDDO
        ! DO MNj= MNi+1,MOLCNT(STN)
        DO MNj= 1,MOLCNT(STN)
          MJi = MPNT(MNj,STN)
          Mf = MPNT(MNj+1,STN)-1          
          MNAj = Mf-MJi+1
          DO DNB =1,RCNT 
           IF((MNAi.EQ.NMOLRi(DNB,STNo).AND.MNAj.EQ.NMOLRj(DNB,STNo)) 
     & .OR.(MNAj.EQ.NMOLRi(DNB,STNo).AND.MNAi.EQ.NMOLRj(DNB,STNo)) )THEN
            IF (MNi .NE. MNj )THEN
              Ri(1) = MCOMS(1,MNi,STN)
              Ri(2) = MCOMS(2,MNi,STN)
              Ri(3) = MCOMS(3,MNi,STN)
              CALL SPFRAC(LVr,Fi,Ri)
              Rj(1) = MCOMS(1,MNj,STN)
              Rj(2) = MCOMS(2,MNj,STN)
              Rj(3) = MCOMS(3,MNj,STN)
              CALL SPFRAC(LVr,Fj,Rj)
              RSQ = 0.0d0
              UDIM(:) = UMDIS(:,DNB,STNo)
              DO D = 1,NDIM
                IF( UDIM(D)) THEN
                  DF(D) = Fj(D) - Fi(D)
                  IF( PBCS(D,STN) )  DF(D) = DF(D) - ANINT(DF(D))  
                ENDIF
              ENDDO
!
!             Calculate cartion coordinates 
!
              CALL REALF(LVr,DF,DRr)
!              DRr(1) = LVr(1,1)*DF(1) + LVr(2,1)*DF(2)+LVr(3,1)*DF(3)
!              DRr(2) = LVr(1,2)*DF(1) + LVr(2,2)*DF(2)+LVr(3,2)*DF(3)
!              DRr(3) = LVr(1,3)*DF(1) + LVr(2,3)*DF(2)+LVr(3,3)*DF(3)
              RSQ = 0.0d0
              DO D=1,NDIM
                RSQ = RSQ + DRr(D)*DRr(D)
              ENDDO

              IF( RSQ.LE.CUTSQ ) THEN
                DELR=SQRT(RSQ)
                L=INT(DELR*RDSCALE)

                RDDST(L,DNB) = RDDST(L,DNB) + 2 
                RDMDST(L,DNB)=RDMDST(L,DNB)+MASSi+MOLMAS(MNj,STN)
!
                NDEN_SUM(DNB) = NDEN_SUM(DNB) + 2 
                MDEN_SUM(DNB)=MDEN_SUM(DNB)+MASSi+MOLMAS(MNj,STN)
!
                IF( VERB ) THEN 
                  WRITE(F,301) MNi,MNj,Ri(:),DELR
     &            ,MOLST(MIi,STN),MOLST(MJi,STN)
                ENDIF
              ENDIF ! RSQ.LE.CUTSQ
             ENDIF
            ENDIF  ! check # of atoms in molecules 
          ENDDO   ! DNB 
        ENDDO    ! MNj
      ENDDO     ! MNi 
!     Print distribution 
!     Loop over rdf's
      DO DNB= 1,RCNT 
        F = FID + DNB 
!       Use density of sphere formormalization 
        NPNTS = NDEN_SUM(DNB)
        NDENS(DNB) = FLOAT(NPNTS)/VCUT(DNB)
        MDENS(DNB) = MDEN_SUM(DNB)/VCUT(DNB)
!       Verbose info
        WRITE(F,*) 'For rdf ',DNB
        WRITE(F,*) " # total nieghbors found", NPNTS
        WRITE(F,*) ' per unit of ',VCUT(DNB),' A ',EX(DNB)
        WRITE(F,*) ' number density is ', NDEN_SUM(DNB),NDENS(DNB)
        WRITE(F,*) ' mass density is ',MDENS(DNB),MDENS/AVO*10.0d0
!     
        DO L=0,SZ-1
          Rin=FLOAT(L)/RDSCALE
          Rout = Rin+DR
          DVi = SPV(DNB)*(Rout**EX(DNB) - Rin**EX(DNB))
          DN = FLOAT(RDDST(L,DNB))
          DM = RDMDST(L,DNB)
          GRN = DN/DVi/NDENS(DNB)
          GRM = DM/DVi/MDENS(DNB)
          WRITE(F,3801) Rin,Rout,DVi,DN,GRN,DM,GRM

        ENDDO
        CLOSE(F)
      ENDDO 
!
      DEALLOCATE( RDMDST,RDDST,SPV,EX,NDENS,MDENS,VCUT
     & ,NDEN_SUM,MDEN_SUM )

      RETURN
 3801 FORMAT(7F16.3)
 301  FORMAT("#",2I12,4F12.2,2I12)
      END SUBROUTINE molrdf

