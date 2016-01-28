!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/09/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Write turbomol file
!
      SUBROUTINE DEFINE_EMBED(STN)
!
      USE specify
      USE structure
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: STN,MN,Mi,Mf,N,I,Ni,Nf,J,JMi,JMf,JN,JI,Nb,MNc
     & ,MTOT,NQM,NMM,NFX,NSHELLS,D,MIND,MNN,MIP,MNP,EL,WRT,ELr
      CHARACTER(ATSZ) :: Atp
      REAL*8 :: RMIN,RMAX,R(NDIM),RSQ,DR,LVr(NDIM,NDIM),FCr(NDIM)
     & ,CENT(NDIM)
      LOGICAL, ALLOCATABLE :: QMADDED(:)
!
!     Set some local variables
!
      MTOT = MOLCNT(STN)
      ALLOCATE( QMADDED(MTOT))

      IF( WONIOMSP(STN) ) THEN
        RMIN = ONRAD(1,STN)
        RMAX = ONRAD(2,STN)
      ENDIF
!
!     Intialize
!
      FTAG(:,:,STN) = 'F'
!
!       Verbose
!
        IF( VERB ) WRITE(6,*) ' Writing Qm/mm files for all'
     & ,MTOT,' in structure',STN
!
!
!       Print xyz for reference
!
      OPEN(UNIT=11,FILE=OUTXYZ,STATUS='unknown')
      WRITE(11,*) NA(STN)
      WRITE(11,'(A100)') TITLE(1)
      DO MN=1,MTOT
       Mi = MPNT(MN,STN)
       Mf = MPNT(MN+1,STN)-1
       DO N=Mi,Mf
           I = MOLST(N,STN) 
           Elr = ELN(I,STN)
           Atp = ATSYM(ELr)
           WRITE(11,1102) Atp,R0(:,I,STN)
        ENDDO
      ENDDO
      CLOSE(11)
!
!
!
      DO MN=1,MTOT
!
!       Set QM molecule list
!
        NQM = 1
        NMM = 0
        NFX = 0
        ONIONQM(STN) = NQM
        QMLIST(NQM,STN) = MN
        QMADDED(:) = .TRUE.
!
        CALL CENTMOL(MN,STN)
        CALL molPBC(STN)          !Apply PBC's 
!
!       Set specified neighbors to QM
!

        IF( QMNB(STN).GT.0) THEN     
          NSHELLS = 0
!
!         Loop over specified # of neighbor shells
!

          DO WHILE ( NSHELLS.LE.QMNB(STN) ) 
!
!           Loop over molecule specified in list
!
            DO MIND = 1,ONIONQM(STN) 
              MNn = QMLIST(MIND,STN)
!
!             Check if already added 
!
              IF( QMADDED(MNn) ) THEN

               QMADDED(MNn) = .FALSE.

               Ni = NINDXML(MNn,STN)           
               Nf = NINDXML(MNn+1,STN)-1
!
               DO Nb=Ni,Nf
                J = NLISTML(Nb,STN)
                IF( QMADDED(J) )THEN
                 QMADDED(J) = .FALSE.
                 JMi = MPNT(J,STN)
                 JMf = MPNT(J+1,STN)-1        
! 
!                Add molecule to list to loop over
!  
                 NQM =  ONIONQM(STN)
                 NQM = NQM + 1
                 ONIONQM(STN) = NQM
                 QMLIST(NQM,STN) = J

                ENDIF 
               ENDDO
             ENDIF ! QMADDED
            ENDDO !MIND
            NSHELLS = NSHELLS + 1
          ENDDO  !WHILE # shells
!
        ELSE
          QMADDED(MN) = .FALSE.
        ENDIF  ! if add qm neighbors
!
!       Cut sphere out of MD system       
!
        IF( WONIOMSP(STN) ) THEN
!
!         Update center of mass
!
          CALL centmass(STN)
          LVr(:,:) = LV(:,:,STN)
!
          FCr(:) = MCCENT(:,STN)

          CALL REALF(LVr,FCr,CENT)

   !       CENT(1)=LVr(1,1)*FCr(1)+LVr(2,1)*FCr(2)+LVr(3,1)*FCr(3)
   !       CENT(2)=LVr(1,2)*FCr(1)+LVr(2,2)*FCr(2)+LVr(3,2)*FCr(3)
   !       CENT(3)=LVr(1,3)*FCr(1)+LVr(2,3)*FCr(2)+LVr(3,3)*FCr(3)
!
!         Loop over all molecule
!
          DO MNc = 1, MOLCNT(STN) 
!
!           IF not added to another list
!
            IF( QMADDED(MNc) ) THEN
!
!              Calculate distance from center molecule 
!
               RSQ = 0.0d0           
               DO D=1,NDIM
                 R(D) = MCOMS(D,MNc,STN)  - CENT(D)
                 RSQ = RSQ + R(D)*R(D) 
               ENDDO 
               DR = SQRT(RSQ)               
               IF( DR.LT.RMIN ) THEN
                 NMM = NMM + 1
                 MMLIST(NMM,STN) = MNc
               ELSEIF( DR.LT.RMAX) THEN 
                 NFX = NFX + 1
                 FXLIST(NFX,STN) = MNc
               ENDIF
            ENDIF 

          ENDDO ! MNc cut out molecules
!
!       If no cut out 
!
        ELSE 
!         Loop over all molecule
!
          DO MNc = 1, MOLCNT(STN) 
!
!           IF not added to another list
!
            IF( QMADDED(MNc) ) THEN
                 NMM = NMM + 1
                 MMLIST(NMM,STN) = MNc
            ENDIF 
!
          ENDDO ! MNc cut out molecules

        ENDIF
!
!       Verbose output 
!
        IF(VERB) THEN
          WRITE(6,*) ' For molecule ',MN
          WRITE(6,*) 'molecules found in region 1',NQM
          WRITE(6,*) 'molecules found in region 2',NMM
          WRITE(6,*) 'molecules found in region 3',NFX
        ENDIF 
! 
!       Save for  print out 
!
        ONIONMM(STN) = NMM
        ONIONFX(STN) = NFX
!
!        CALL WRITE_PEECM(STN,MN)
        IF( WCOM )  CALL WRITE_ONIOM(STN,MN)
        IF( WTBML ) CALL WRITE_TBEMBED(STN,MN)
     
      ENDDO ! MN

      DEALLOCATE( QMADDED )

      RETURN 
 5601 FORMAT(A6,3F24.6)
 1102 FORMAT(A6,"  ",3F24.8)
      END SUBROUTINE DEFINE_EMBED 
