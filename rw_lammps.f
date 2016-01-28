!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      Read lammps data file 
!
      SUBROUTINE READ_LMP(STN)
      USE specify
      USE structure
      USE potential
      USE elements  
!
      IMPLICIT none
!
      INTEGER   :: STN

!
      RETURN
      END SUBROUTINE READ_LMP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 06/13/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Write lammps data file 
!
      SUBROUTINE WRITE_LMP(STN)
!
      USE specify
      USE structure
      USE potential 
      USE elements 
!
      IMPLICIT none
!
      INTEGER   :: STN,I,MN,EL,EIND,ATYPE_CNT
     & ,T,T_REF,Ip,Jp,Kp,Lp,Ti,Tj,Tk,Tl
      REAL*8 :: xlo,xhi,ylo,yhi,zlo,zhi,Q
      CHARACTER(3) ATi,ATr
      LOGICAL :: NEW_TYPE
!
!
!
      ATYPE_CNT = 1
      I = 1
      ATi = TYP(I,STN)
      TYP_IND(I,STN) = ATi  
      TYP_REF(I,STN) = ATYPE_CNT
      TYP_CNT(ATYPE_CNT,STN) = 1
      DO I =2,NA(STN)
         ATi = TYP(I,STN)
         NEW_TYPE = .TRUE.
         DO T=1,ATYPE_CNT
           ATr = TYP_IND(T,STN)
           IF( ATi.EQ.ATr) THEN
             NEW_TYPE = .FALSE.
             T_ref = T
           ENDIF
         ENDDO 
         IF( NEW_TYPE ) THEN
           ATYPE_CNT = ATYPE_CNT + 1
           TYP_IND(ATYPE_CNT,STN) = ATi
           TYP_REF(I,STN) = ATYPE_CNT
           TYP_CNT(ATYPE_CNT,STN) = 1
         ELSE
           TYP_CNT(T_ref,STN) = TYP_CNT(T_ref,STN) + 1
            TYP_REF(I,STN) = T_ref
         ENDIF
      ENDDO
!
      IF ( DEBUG ) THEN 
        DO I =1,NA(STN)
          WRITE(102,*)  I,TYP_REF(I,STN)
        ENDDO
      ENDIF 

!
      IF( VERB ) THEN
         WRITE(*,*) 'FF types found',ATYPE_CNT
         DO T=1,ATYPE_CNT
           WRITE(*,*) T,TYP_IND(T,STN),TYP_CNT(T,STN)
         ENDDO
      ENDIF
!
      WRITE(6,*) 'Writing LAMMPS data file ',OUTLMP
      OPEN(UNIT=56,FILE=OUTLMP,STATUS='unknown')
      WRITE(56,5601)TITLE
      WRITE(56,*)
      WRITE(56,5602) NA(STN)
      WRITE(56,5603) NB(STN)
      WRITE(56,5604) NANG(STN)
      WRITE(56,5605) NDIH(STN)
      WRITE(56,5606) NIMP(STN)
      WRITE(56,*)
      WRITE(56,5607) ELT(STN)
      WRITE(56,5608) NBTYP(STN)
      WRITE(56,5609) NATYP(STN)
      WRITE(56,5610) NDTYP(STN)
      WRITE(56,5611) NITYP(STN)
      WRITE(56,*)
!
!     Calculate box size
!
      xlo = -1.0d0*LC(1,STN)/2.0d0
      xhi = LC(1,STN)/2.0d0
      ylo = -1.0d0*LC(2,STN)/2.0d0
      yhi = LC(2,STN)/2.0d0
      zlo = -1.0d0*LC(3,STN)/2.0d0
      zhi = LC(3,STN)/2.0d0
      WRITE(56,5612) xlo,xhi
      WRITE(56,5613) ylo,yhi
      WRITE(56,5614) zlo,zhi
      WRITE(56,*)
!
!     Masses
!
      WRITE(56,*)
      WRITE(56,5615) 
      WRITE(56,*)
      DO T=1,ATYPE_CNT
       ATr = TYP_IND(T,STN)
       DO I=1,NA(STN)
         ATi= TYP(I,STN)
         IF( ATi.EQ.ATr) EL = ELN(I,STN)
       ENDDO
       WRITE(56,5620) T,AMASS(EL)
      ENDDO 
!
!     Pair Coeffs
!
      WRITE(56,*)
      WRITE(56,5616) 
      WRITE(56,*)
!      DO C=1,1 !ELT(STN)
!       WRITE(56,5620) C,
!      ENDDO 
!
!     Bond Coeffs
!
      WRITE(56,*)
      WRITE(56,5617) 
      WRITE(56,*)
!      DO EL=1,ELT(STN)
!    WRITE(56,5620) EL,AMASS(EL)
!      ENDDO 
!
!     Angle Coeffs
!
      WRITE(56,*)
      WRITE(56,5618) 
      WRITE(56,*)
!      DO EL=1,ELT(STN)
!       WRITE(56,5620) EL,AMASS(EL)
!      ENDDO 
!
!     Dihedral Coeffs
!
      WRITE(56,*)
      WRITE(56,5619) 
      WRITE(56,*)
!      DO EL=1,ELT(STN)
!       WRITE(56,5620) EL,AMASS(EL)
!      ENDDO 

!
!     Atomic positions  
!
      WRITE(56,*)
      WRITE(56,5630) 
      WRITE(56,*)
      DO I=1,NA(STN)
        MN = MOLN(I,STN)       ! Molecule #
        EL = ELREF(I,STN)      ! Element/atom type 
        Q = ACHG(I,STN)        ! Charge 
        ! WRITE(56,5631) I,MN,EL,Q,R0(:,I,STN)
        WRITE(56,5641) I,MN,EL,R0(:,I,STN)
      ENDDO
!
!     Bonds 
!  
      WRITE(56,*)
      WRITE(56,5632) 
      WRITE(56,*)
      DO I=1,NB(STN)
        Ip = BNDI(I,STN)
        Jp = BNDJ(I,STN)
        Ti = TYP_REF(Ip,STN)
        Tj = TYP_REF(Jp,STN)
        ! WRITE(56,5633)  I,Ti,Tj,Ip,Jp
        WRITE(56,5643)  I,Ti,Ip,Jp
      ENDDO 
!
!     Angles 
!  
      WRITE(56,*)
      WRITE(56,5634) 
      WRITE(56,*)
      DO I=1,NANG(STN)
        Ip = ANGI(I,STN)
        Jp = ANGJ(I,STN)
        Kp = ANGK(I,STN)
        Ti = TYP_REF(Ip,STN)
        Tj = TYP_REF(Jp,STN)
        Tk = TYP_REF(Kp,STN)
!        WRITE(56,5635)I,Ti,Tj,Tk,Kp,Ip,Jp
        WRITE(56,5645)I,Ti,Kp,Ip,Jp
      ENDDO 
!
!     Dihedrals 
!  
      WRITE(56,*)
      WRITE(56,5636) 
      WRITE(56,*)
      DO I=1,NDIH(STN)
        Ip = DIHI(I,STN)
        Jp = DIHJ(I,STN)
        Kp = DIHK(I,STN)
        Lp = DIHL(I,STN)
        Ti = TYP_REF(Ip,STN)
        Tj = TYP_REF(Jp,STN)
        Tk = TYP_REF(Kp,STN)
        Tl = TYP_REF(Lp,STN)
!        WRITE(56,5637)I,Ti,Tj,Tk,Tl,Kp,Ip,Jp,Lp
        WRITE(56,5647)I,Ti,Kp,Ip,Jp,Lp
      ENDDO 
!
!
!     Impropers
!  
      WRITE(56,*)
      WRITE(56,5638) 
      WRITE(56,*)
      DO I=1,NIMP(STN)
        Ip = IMPI(I,STN)
        Jp = IMPJ(I,STN)
        Kp = IMPK(I,STN)
        Lp = IMPL(I,STN)
        Ti = TYP_REF(Ip,STN)
        Tj = TYP_REF(Jp,STN)
        Tk = TYP_REF(Kp,STN)
        Tl = TYP_REF(Lp,STN)
        WRITE(56,5639)I,Ti,Tj,Tk,Tl,Ip,Jp,Kp,Lp
      ENDDO 
!
      WRITE(56,*) ''
      CLOSE(56)
!
      RETURN
 5601 FORMAT('LAMMPS data file:',A70)
 5602 FORMAT(I16,'  atoms ')
 5603 FORMAT(I16,'  bonds ')
 5604 FORMAT(I16,'  angles')
 5605 FORMAT(I16,'  dihedrals ')
 5606 FORMAT(I16,'  impropers ')
 5607 FORMAT(I16,'  atom types ')
 5608 FORMAT(I16,'  bond types ')
 5609 FORMAT(I16,'  angle types')
 5610 FORMAT(I16,'  dihedral types ')
 5611 FORMAT(I16,'  improper types ')
 5612 FORMAT(2F26.16,' xlo xhi')
 5613 FORMAT(2F26.16,' ylo yhi')
 5614 FORMAT(2F26.16,' zlo zhi')
 5615 FORMAT(' Masses ')
 5616 FORMAT(' Pair Coeffs ')
 5617 FORMAT(' Bond Coeffs ')
 5618 FORMAT(' Angle Coeffs ')
 5619 FORMAT(' Dihedral Coeffs  ')
 5620 FORMAT(I16,F26.16)
 5630 FORMAT('Atoms ')
 5631 FORMAT(3I8,F12.4,3F16.8,' 0 0 0 ')
 5641 FORMAT(3I8,3F16.8)
 5632 FORMAT(' Bonds ')
 5633 FORMAT(I8,4I12)
 5643 FORMAT(I8,3I12)
 5634 FORMAT(' Angles')
 5635 FORMAT(I8,6I6)
 5645 FORMAT(I8,4I6)
 5636 FORMAT(' Dihedrals ')
 5637 FORMAT(I8,8I6)
 5647 FORMAT(I8,5I6)
 5638 FORMAT(' Impropers ')
 5639 FORMAT(I8,8I6)
      END SUBROUTINE WRITE_LMP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
