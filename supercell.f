!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 10/21/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Create super cell


      SUBROUTINE supercell(STNo)
      USE structure
      USE potential 
      USE specify 
      USE build 
!
      IMPLICIT none
!
      INTEGER :: MCNTr,A,STNo,STNi,NAr,NAi,I,J
     & ,Il,Kl,LL,Is,Js,Ks,RNr,MOLr,GCHRr
      REAL*8 :: LVr(NDIM,NDIM),Fs(NDIM),Ir,Jr,Kr,ANa,ANb 
      REAL*8, allocatable :: RI(:,:),RF(:,:) 
!
!     Initialize
      STNi = STNo + STOTo 
      NAr = NA(STNo)
      NAi = NA(STNo)*SC(1,STNo)*SC(2,STNo)*SC(3,STNo) 
      A = 0 !Number of translated atoms
      MCNTr = 0
      RNr = 0 
      GCHRr = 0

!     Translate to fractional coordinates 
      LVr(:,:) = LV(:,:,STNo)
      ALLOCATE(Rf(NDIM,NAi),Ri(NDIM,NAr)  )
      Ri(:,:) = R0(:,:,STNo)
!     Exand along lattice vectors
      Il = INT ( SC(1,STNo)-1 )
      Kl = INT ( SC(2,STNo)-1 )
      Ll = INT ( SC(3,STNo)-1 )
      IF(VERB ) THEN
        WRITE(6,*) "Creating supercell of ",SC(:,STNo)
       ENDIF
       DO  Is=0,Il
        DO Js =0,Kl
          DO Ks =0,Ll
            IF( RANROT(STNo) ) THEN 
!                   Get random #'s 
                    SEED = SEED*SEED
                    SEED = SEED - SEED/2 
                    ANa = ( 2*PI )*ran(seed) 
                    ANb = ( 2*PI )*ran(seed) 
                 CALL ROTATE(STNo,ANa,ANb)
            ENDIF 
            DO I=1,NAr
              A = A + 1
              CALL ATPASS(I,STNo,A,STNi)
              Ir = FLOAT( Is )
              Jr = FLOAT( Js )
              Kr = FLOAT( Ks )
              Rf(1,A)=Ri(1,I)+Ir*LVr(1,1)+Jr*LVr(2,1)+Kr*LVr(3,1)
              Rf(2,A)=Ri(2,I)+Ir*LVr(1,2)+Jr*LVr(2,2)+Kr*LVr(3,2)
              Rf(3,A)=Ri(3,I)+Ir*LVr(1,3)+Jr*LVr(2,3)+Kr*LVr(3,3)
              RN(A,STNi) = RN(I,STNo) + RNr 
              CHN(A,STNi) = CHN(I,STNo) + GCHRr 
              MOLN(A,STNi) = MOLN(I,STNo) +  MCNTr
            ENDDO
            RNr   = RN(A,STNi)
            GCHRr = CHN(A,STNi) 
            MCNTr = MOLN(A,STNi)
          ENDDO
        ENDDO
      ENDDO
!     Update global
      NA(STNi) = A
      R0(:,1:NAi,STNi) = Rf(:,1:A) 
!     
      DO I=1,NDIM
        LC(I,STNi) = LC(I,STNo)* SC(I,STNo) 
        DO J=1,NDIM
          LV(I,J,STNi) = LV(I,J,STNo)* SC(I,STNo) 
        ENDDO
      ENDDO
      DEALLOCATE( RF,Ri )
      RETURN 
 5601 FORMAT(2I6," 1 0.07 0.07 ")
      END SUBROUTINE supercell

     
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 05/29/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  add top layer of cystal to complete slab
!
      SUBROUTINE SLAB(STNi)
      USE structure
      USE potential 
      USE specify 
      USE build 
!
      IMPLICIT none
!
      INTEGER :: NAr,ACNT,Is,Js,Ks,I,J,STNi,MCNTr,STNo,D,Ir
      REAL*8 :: SCUTr,LVr(NDIM,NDIM),Fs(NDIM),Rr(NDIM)
!
!     Set local variables 
      STNo = STNi - STOTo 
      NAr = NA(STNi) 
      SCUTr = SCUT(STNo)
      ACNT = NAr
      MCNTr = MOLCNT(STNi)
!     Translate to fractional coordinates 
      LVr(:,:) = LV(:,:,STNi)
      ALLOCATE( Fr(NDIM,NAr) )
      CALL fracR(STNi,LVr,1)
!
!     Set lattice direction to expand 
      Is = 0 !SC(1) - 1
      Js = 0 !SC(2) - 1
      Ks = SC(SDIM,STNo) 
!
      IF(VERB ) THEN
        WRITE(6,*) "Creating slab of", STNi
     &  ," with positions below ",SCUTr ," and s dim",SDIM
     &  ," with multiple ",Ks 
     &  ," with ",SATOMN(STNo)," surface atoms "
        WRITE(*,*) ' Lattice vectors'
        WRITE(*,*) LVr(1,:)
        WRITE(*,*) LVr(2,:)
        WRITE(*,*) LVr(3,:)
        WRITE(*,*) ''
      ENDIF 


      DO I=1,NAr
        IF( R0(SDIM,I,STNi).LT.SCUTr ) THEN
          ACNT = ACNT + 1
          CALL ATPASS(I,STNi,ACNT,STNi) 
          Fs(1) = Fr(1,I) + FLOAT( Is )
          Fs(2) = Fr(2,I) + FLOAT( Js )
          Fs(3) = Fr(3,I) + FLOAT( Ks )
          CALL REALF(LVr,Fs,Rr)
          R0(:,ACNT,STNi) = Rr(:)      

!          MOLN(Ii,STNii) = MCNTr + 1          ! Molecule #
!          CHN(Ii,STNii)  = CHN(Io,STNio)           ! Charge group #
          Ir = ACNT 
          DO  J=1,SATOMN(STNo)
              ACNT = ACNT + 1
              IF(VERB) THEN
                WRITE(6,*) 'Adding surface atom ',J,SATOMID(J,STNo)
                WRITE(6,*) ' to atom ',Ir,' by '
                WRITE(6,'(3F6.2)') SATOMV(:,J,STNo) 
                WRITE(6,'(3F6.2)') R0(:,Ir,STNi)
              ENDIF
              CALL ATPASS(Ir,STNi,ACNT,STNi) 
!             Write added surface atoms  
              DO D=1,NDIM
                R0(D,ACNT,STNi) = R0(D,Ir,STNi) + SATOMV(D,J,STNo)
              ENDDO
              TYP(ACNT,STNi)  = SATOMID(J,STNo) 
              GRID(ACNT,STNi)  = SATOMID(J,STNo) 
              RID(ACNT,STNi)  = "AuI"  !RIDi(STNi,N)
              RN(ACNT,STNi)   = RN(J,STNo)+ MCNTr
              ELN(ACNT,STNi)  = 0     ! ELNi(STNi,N)
              MOLN(ACNT,STNi) = MCNTr + 1
              ACHG(ACNT,STNi) = 0.0d0  ! CHRGi(STNi,N) 
            ENDDO 
        ENDIF 
      ENDDO 
!     Update global atom list
      IF( VERB) THEN
         WRITE(6,*) ACNT-NAr," atoms added"
      ENDIF 
      NA(STNi) = ACNT
      DEALLOCATE( Fr )
!
      RETURN
      END SUBROUTINE SLAB
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
