!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Pass properties from one structure to another 
!
      SUBROUTINE ATPASS(Io,STNo,Ii,STNi)
      USE structure
      USE specify
      USE potential 
!
      IMPLICIT none
!
      INTEGER :: STNo,STNi,Io,Ii
!
      ELN(Ii,STNi)  = ELN(Io,STNo)           ! Atomic #
      ELREF(Ii,STNi)  = ELREF(Io,STNo)       !  element ref for element count
      TYP(Ii,STNi)  = TYP(Io,STNo)           ! Force field type
      TYPN(Ii,STNi)  = TYPN(Io,STNo)         ! Force field type number (tinker) 
      RID(Ii,STNi)  = RID(Io,STNo)           ! Residue ID
      RN(Ii,STNi)   = RN(Io,STNo)            ! Residue #
      MOLN(Ii,STNi) = MOLN(Io,STNo)          ! Molecule #
      ACHG(Ii,STNi) = ACHG(Io,STNo)          ! Charge #
      CHN(Ii,STNi)  = CHN(Io,STNo)           ! Charge group #
      GRID(Ii,STNi) = GRID(Io,STNo)          ! GROMACS ID
      THRM(Ii,STNi)   = THRM(Io,STNo)        ! Thermostat flag (REBO)
      R0(:,Ii,STNi)   = R0(:,Io,STNo)        ! coodinates 
      R1(:,Ii,STNi)   = R1(:,Io,STNo)        ! velocity
      R2(:,Ii,STNi)   = R2(:,Io,STNo)        ! acceleration 
      FTAG(:,Ii,STNi) = FTAG(:,Io,STNo)      ! fix tag
      ONMTAG(Ii,STNi) = ONMTAG(Io,STNo)      ! oniom  tag
      AMAS(Ii,STNi) = AMAS(Io,STNo)          ! Atomic mass
!
      RETURN
      END SUBROUTINE atpass
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Pass properties from one structure to another 

      SUBROUTINE STRPASS(STNo,STNi)
      USE structure
      USE specify
      USE potential
!
      IMPLICIT none
!
      INTEGER :: STNo,STNi,NAr,ACNT,ELr,I

!     Set local variables 
      NAr=Na(STNo)
!     Initialize
      ACNT = 0

      DO I=1,NAr
        ACNT = ACNT + 1
C        OTOI(I) = ACNT
        CALL ATPASS(I,STNo,I,STNi)
      ENDDO
      TITLE(STNi) = TITLE(STNo)
      NA(STNi)      = NA(STNo)               ! Number of atoms
!     Structure 
      LV(:,:,STNi) = LV(:,:,STNo)
      LC(:,STNi) = LC(:,STNo)
      LA(:,STNi) = LA(:,STNo)
!     Element list info
      ELT(STNi) = ELT(STNo)
      NEL(:,STNi) = NEL(:,STNo)
      ELCNT(:,STNi) = ELCNT(:,STNo)
      ELIST(:,:,STNi) = ELIST(:,:,STNo)
!     Neighbor list
      NINDX(:,STNi) = NINDX(:,STNo)
      NLIST(:,STNi) = NLIST(:,STNo)
!     Molecule list 
      MOLCH(STNi) = MOLCH(STNo)
      MOLCHR(:,STNi) = MOLCHR(:,STNo)
      MOLMAS(:,STNi) = MOLMAS(:,STNo)
      MOLCNT(STNi)    = MOLCNT(STNo)
      MPNT(:,STNi)    = MPNT(:,STNo)
      MCOMS(:,:,STNi) = MCOMS(:,:,STNo)
      MOLST(:,STNi)   = MOLST(:,STNo) 
!
!       Bond, angle and dih get
      GBONDS(STNi) = GBONDS(STNo)
!      
!     PBC's
      PBCS(:,STNi) = PBCS(:,STNo)
      IF(DEBUG)  WRITE(120,*) STNi, PBCS(:,STNi) 
!    
!     Connection count 
!
      NB(STNi) = NB(STNo)
           BNDI(:,STNi) = BNDI(:,STNo)
           BNDJ(:,STNi) = BNDJ(:,STNo)      
           BTYP(:,STNi) = BTYP(:,STNo) 
           BVAL(:,STNi) = BVAL(:,STNo) 
           BCNST(:,STNi) = BCNST(:,STNo) 


      NANG(STNi) = NANG(STNo)
           ANGI(:,STNi) = ANGI(:,STNo)
           ANGJ(:,STNi) = ANGJ(:,STNo)           
           ANGK(:,STNi) = ANGK(:,STNo)      
           ATYP(:,STNi) = ATYP(:,STNo) 
           AVAL(:,STNi) = AVAL(:,STNo) 
           ACNST(:,STNi) = ACNST(:,STNo) 
     
      NDIH(STNi) = NDIH(STNo)
           DIHI(:,STNi) = DIHI(:,STNo)
           DIHJ(:,STNi) = DIHJ(:,STNo)           
           DIHK(:,STNi) = DIHK(:,STNo)           
           DIHL(:,STNi) = DIHL(:,STNo)           
           DTYP(:,STNi) = DTYP(:,STNo) 
           DVAL(:,STNi) = DVAL(:,STNo) 
           DCNST(:,STNi) = DCNST(:,STNo) 

      NIMP(STNi) = NIMP(STNo)
           IMPI(:,STNi) = IMPI(:,STNo)
           IMPJ(:,STNi) = IMPJ(:,STNo)           
           IMPK(:,STNi) = IMPK(:,STNo)           
           IMPL(:,STNi) = IMPL(:,STNo)           
           ITYP(:,STNi) = ITYP(:,STNo) 
           IVAL(:,STNi) = IVAL(:,STNo) 
           ICNST(:,STNi) = ICNST(:,STNo) 
!
      RETURN
      END SUBROUTINE strpass

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Reduce  structures into final structure

      SUBROUTINE STRREND
      USE structure
      USE specify
      USE build 
      USE potential
!
      IMPLICIT none
!
      INTEGER :: STNo,STNi,NAr,ACNT,ELr,I,S,NBl,NAl,NDl,NIl
     & ,Ii,Ji,Ki,Li ,Io,Jo,Ko,Lo,TP,MAXCHN_STN,MAXCHN_L,CHN_L
      REAL*8 :: VAL 
      CHARACTER(CHSZ) :: ESTAT,CONSTi
!
      IF(VERB)WRITE(6,*) "-------------------------------------------"
!
!     Initialize
      ACNT = 0
      NBl = 0
      NAl = 0
      NDl = 0
      NIl = 0
      IF(VERB)WRITE(6,*) "   Rendering final structure ",STOT
!
!     Loop over all modified structures            
      MAXCHN_L = 0
      DO STNi=STOTo+1,STOT-1
!
!       Add  atoms to final 
!
        NAr=NA(STNi)
        IF(VERB) THEN 
          WRITE(6,*) "Adding ",STNi
     & ," to final structure, with",NAr," atoms"
        ENDIF
        MAXCHN_STN = MAXCHN_L
        DO I=1,NAr
          ACNT = ACNT + 1
          CALL ATPASS(I,STNi,ACNT,STOT)
          CHN_L =  MAXCHN_STN + CHN(I,STNi) ! Charge group #
          IF(  CHN_L .GT. MAXCHN_L )  MAXCHN_L = CHN_L 
          CHN(ACNT,STOT)  = CHN_L
          REND_IND(I,STNi) = ACNT
        ENDDO
      ENDDO
!     Loop over all modified structures            
      DO STNi=STOTo+1,STOT-1
        IF(VERB) THEN
           WRITE(*,*) 'Render structure ',STNi,' into ',STOT
        ENDIF
!       Molecule list 
        
        MOLCNT(STOT)    = MOLCNT(STNi)
        MPNT(:,STOT)    = MPNT(:,STNi)
        MCOMS(:,:,STOT) = MCOMS(:,:,STNi)
        MOLST(:,STOT)   = MOLST(:,STNi) 
        MOLCH(STOT) = MOLCH(STNi)
!
!       Bond, angle and dih get
        GBONDS(STOT) = GBONDS(STNi)
!
!       Add  bonds to final 
!    
        IF(VERB) THEN
           WRITE(*,*) 'Render  ',NB(STNi) ,' bonds '
        ENDIF
        DO I=1,NB(STNi)
           NBl = NBl + 1
           Ii = BNDI(I,STNi)
           Ji = BNDJ(I,STNi)   
!
           TP  = BTYP(I,STNi)
           VAL = BVAL(I,STNi)
           CONSTi = BCNST(I,STNi)
!
           IF( Ii.GT.NA(STNi) .OR. Ji.GT.NA(STNi) ) THEN
              WRITE(ESTAT,*) STNI,Ii," ",Ji," ",NBL,NA(STNi) 
              CALL PRNERROR(-9,ESTAT)
           ENDIF
           Io = REND_IND(Ii,STNi)
           Jo = REND_IND(Ji,STNi)
!
           BNDI(NBl,STOT) = Io
           BNDJ(NBl,STOT) = Jo      
!
           BTYP(NBl,STOT) = TP
           BVAL(NBl,STOT) = VAL
           BCNST(NBl,STOT) = CONSTi 
        ENDDO
!
!       Add  angles to final 
!    
        DO I=1,NANG(STNi)
           NAl = NAl + 1
           Ii = ANGI(I,STNi)
           Ji = ANGJ(I,STNi)   
           Ki = ANGK(I,STNi)   
           TP  = ATYP(I,STNi)
           VAL = AVAL(I,STNi)
           CONSTi = ACNST(I,STNi)

           IF( Ii.GT.NA(STNi) .OR. Ji.GT.NA(STNi)
     &        .OR. Ki.GT.NA(STNi) ) THEN
              WRITE(ESTAT,*) Ii," ",Ji," ",Ki,NA(STNi)
              CALL PRNERROR(-9,ESTAT)
           ENDIF

           Io = REND_IND(Ii,STNi)
           Jo = REND_IND(Ji,STNi)
           Ko = REND_IND(Ki,STNi)
           ANGI(NAl,STOT) = Io
           ANGJ(NAl,STOT) = Jo      
           ANGK(NAl,STOT) = Ko      
           ATYP(NAl,STOT) = TP
           AVAL(NAl,STOT) = VAL
           ACNST(NAl,STOT) = CONSTi 
        ENDDO
!
!       Add  dih to final 
!    
        DO I=1,NDIH(STNi)
           NDl = NDl + 1
           Ii = DIHI(I,STNi)
           Ji = DIHJ(I,STNi)   
           Ki = DIHK(I,STNi)   
           Li = DIHL(I,STNi)   

           TP  = DTYP(I,STNi)
           VAL = DVAL(I,STNi)
           CONSTi = DCNST(I,STNi)
           
           Io = REND_IND(Ii,STNi)
           Jo = REND_IND(Ji,STNi)
           Ko = REND_IND(Ki,STNi)
           Lo = REND_IND(Li,STNi)
           DIHI(NDl,STOT) = Io
           DIHJ(NDl,STOT) = Jo      
           DIHK(NDl,STOT) = Ko      
           DIHL(NDl,STOT) = Lo      

           IF(DEBUG) WRITE(602,*)Io,Jo,Ko,Lo 

                  
           DTYP(NDl,STOT) = TP
           DVAL(NDl,STOT) = VAL
           DCNST(NDl,STOT) = CONSTi 
        ENDDO
!
!       Add  imp to final 
!    
        DO I=1,NIMP(STNi)
           NIl = NIl + 1
           Ii = IMPI(I,STNi)
           Ji = IMPJ(I,STNi)   
           Ki = IMPK(I,STNi)   
           Li = IMPL(I,STNi)   

           TP  = ITYP(I,STNi)
           VAL = IVAL(I,STNi)
           CONSTi = ICNST(I,STNi)
           
           Io = REND_IND(Ii,STNi)
           Jo = REND_IND(Ji,STNi)
           Ko = REND_IND(Ki,STNi)
           Lo = REND_IND(Li,STNi)
           IMPI(NIl,STOT) = Io
           IMPJ(NIl,STOT) = Jo      
           IMPK(NIl,STOT) = Ko      
           IMPL(NIl,STOT) = Lo      
                  
           ITYP(NIl,STOT) = TP
           IVAL(NIl,STOT) = VAL
           ICNST(NIl,STOT) = CONSTi 
       
        ENDDO
!
        TITLE(STOT) = TITLE(STNi)
!       Structure 
        LV(:,:,STOT) = LV(:,:,STNi)
        LC(:,STOT) = LC(:,STNi)
        LA(:,STOT) = LA(:,STNi)
!       Element list info
        ELT(STOT) = ELT(STNi)
        NEL(:,STOT) = NEL(:,STNi)
        ELCNT(:,STOT) = ELCNT(:,STNi)
        ELIST(:,:,STOT) = ELIST(:,:,STNi)
!       Neighbor list
        NINDX(:,STOT) = NINDX(:,STNi)
        NLIST(:,STOT) = NLIST(:,STNi)
!       PBC's
        PBCS(:,STOT) = PBCS(:,STNi)
      ENDDO
!
!     Set final counts
!
      NA(STOT) = ACNT 
      NB(STOT) = NBl
      NANG(STOT) = NAl
      NDIH(STOT) = NDl
      NIMP(STOT) = NIl
!      
!     Check for fixed lattice vector 
      DO S = 1,STOTo 
        IF( FIXLV(S) ) STNi = S 
      ENDDO
!     Asign final lattice vector 
      TITLE(STOT) = TITLE(STNi)
      LV(:,:,STOT) = LV(:,:,STNi)
      LC(:,STOT) = LC(:,STNi)
      LA(:,STOT) = LA(:,STNi)
!
      CALL CALCNELM(STNi)
!
      IF(VERB)WRITE(6,*) "-------------------------------------------"
      RETURN
      END SUBROUTINE STRREND


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/14/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  creat an interface between structures

      SUBROUTINE STACK
      USE structure
      USE potential 
      USE specify 
!
      IMPLICIT none
!
      INTEGER :: I,EL,D,Si,STK,Sj
      REAL*8 :: SHIFT 
    
!
      DO STK=1,NSTK
        Si= STKI(STK)
        Sj = STKJ(STK)
        SHIFT = 0.0d0 
        IF(  .NOT.SHFIX(STK) ) THEN
!         If other structure has already been shifted relative to Si keep same shift value
          IF(VERB) WRITE(6,601) Sj,Si
          CALL MAXMIN(Si)
          CALL MAXMIN(Sj)
          SHIFT = STKBUF(STK) - (MMN(1,SDIM,Sj) - MMN(2,SDIM,Si) )
        ENDIF
        IF( VERB) THEN
          WRITE(*,*) 'Stacking ',Sj,' on ',Si, 'in ',SDIM
          WRITE(*,*)' Min',Sj,MMN(Sj,1,SDIM) 
          WRITE(*,*)' MAx',Si,MMN(Si,2,SDIM) 
          WRITE(*,*)" Fix stack",SHFIX(STK)
          WRITE(*,*)' Add buffer',STKBUF(STK)
          WRITE(*,*)' Shift',SHIFT
        ENDIF
        DO I=1,NA(Sj)
            R0(SDIM,I,Sj)  = R0(SDIM,I,Sj) + SHIFT
        ENDDO
      ENDDO
      RETURN
601   FORMAT("Calculating shift for ",2I8," stack")
      END SUBROUTINE STACK
!

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      Set atomic properties 
      

      SUBROUTINE TTOG(STN)
      USE structure
!
      IMPLICIT none
!
      INTEGER :: STN,I 
!
      DO I=1,NA(STN)
         GRID(I,STN) = TYP(I,STN)
      ENDDO 
      RETURN
      END SUBROUTINE TTOG
