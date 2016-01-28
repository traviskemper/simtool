!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/29/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     switch from gromacs ID's to tinker #'s 
!
      SUBROUTINE ATYPNOPLSAA(STN)
      USE structure
      USE const
      USE specify 
      USE potential 
!
      IMPLICIT none
!
      INTEGER :: STN,I
!
      DO I=1,NA(STN)
        IF( TYP(I,STN).EQ."CR" ) TYPN(I,STN) =  72  !  "HID & HIE CE1" 6 12.011 3
        IF( TYP(I,STN).EQ."CB" ) TYPN(I,STN) =  50  !  "Adenine C4" 6 12.011 3  
        IF( TYP(I,STN).EQ."CT") TYPN(I,STN) =  13  !  "N-Me-HIS CB" 6 12.011 4  
        IF( TYP(I,STN).EQ. 'CA' ) TYPN(I,STN) =  38  !  "Aromatic C" 6 12.011 3  
        IF( TYP(I,STN).EQ. 'C!' ) TYPN(I,STN) =  76  !  "Biphenyl C1" 6 12.011 3  
        IF( TYP(I,STN).EQ. 'HA' ) TYPN(I,STN) =  39  !  "Aromatic H-C" 1 1.008 1  
        IF( TYP(I,STN).EQ. 'NA' ) TYPN(I,STN) =  47  !  "Adenine & Guanine N9" 7 14.007 2
        IF( TYP(I,STN).EQ. 'NB' ) TYPN(I,STN) =  51  !  "Adenine & Guanine N7" 7 14.007 2
      ENDDO

      RETURN
      END SUBROUTINE  ATYPNOPLSAA 

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/29/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     switch from gromacs ID's to tinker #'s  
!
      SUBROUTINE ATYPNAMOEBA(STN)
      USE structure
      USE const
      USE specify
      USE potential
!
      IMPLICIT none
!
      INTEGER :: STN,I,TNCNT(1000)
!
      DO I=1,NA(STN)
        IF( TYP(I,STN).EQ. 'CA' ) TYPN(I,STN) =  213  !  "Benzene C" 6 12.011 3
        IF( TYP(I,STN).EQ. 'HA' ) TYPN(I,STN) =  214  !  "Benzene HC" 1 1.008 1
        IF( TYP(I,STN).EQ. "CR" ) TYPN(I,STN) =  331  !  "HID & HIE CE1" 6 12.011 3
        IF( TYP(I,STN).EQ. "CB" ) TYPN(I,STN) =  285 !  "Adenine C4" 6 12.011 3  
        IF( TYP(I,STN).EQ. 'C!' ) TYPN(I,STN) =  332 !  "Biphenyl C1" 6 12.011 3  
        IF( TYP(I,STN).EQ. 'NA' ) TYPN(I,STN) =  254 !  "Adenine & Guanine N9" 7 14.007 2
        IF( TYP(I,STN).EQ. 'NB' ) TYPN(I,STN) =  258 !  "Adenine & Guanine N7" 7 14.007 2
       ENDDO
       IF(VERB) THEN
       WRITE(*,*) 'Atom type #s set for amoeba for structure ',STN
       ENDIF
!
      RETURN
      END SUBROUTINE  ATYPNAMOEBA
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/29/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     switch from gromacs ID's to tinker #'s  
!
      SUBROUTINE ATYPNADD(STN)
      USE structure
      USE const
      USE specify
      USE potential 
!
      IMPLICIT none
!
      INTEGER :: STN,I,C
 
!   
!     Set intital atom type number
!
      C = TNKADDI(STN)
!
      DO I=1,NA(STN)
        TYPN(I,STN) = C  ! Just add to make some ##
        C = C + 1
      ENDDO
!
      RETURN
      END SUBROUTINE  ATYPNADD
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 5.0 06/29/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     switch from gromacs ID's to tinker #'s  
!
      SUBROUTINE MOL_ATYPNADD(STN)
      USE structure
      USE const
      USE specify
      USE potential 
!
      IMPLICIT none
!
      INTEGER :: STN,STNo,MNi,Mi,Mf,MA,Ni,I,Co,C,Io
!
      STNo = STN - STOTo 
!   
!     Set intital atom type number
!
      Co =  MTNKADDI(STNo)   
!
!     Verbose output
! 
      IF( VERB ) THEN
       WRITE(*,*) ' Tinker reference numbers ',STNo
     & ,' will be modified starting from ',Co
     & ,' for all ',MOLCNT(STN),'molecule '
      ENDIF
      DO MNi=1,MOLCNT(STN)
        Mi = MPNT(MNi,STN)
        Mf = MPNT(MNi+1,STN)-1           
        MA = Mf-Mi+1
        C = Co 
!
!       Find the lowest atom #
!  
        Io = 1000000000
        DO Ni=Mi,Mf 
          I = MOLST(Ni,STN)
          IF ( I .LT. Io ) Io = I
        ENDDO 

        DO I=Io,Io+MA-1
          TYPN(I,STN) =  C 
          C = C + 1
        ENDDO

      ENDDO
      IF(VERB) THEN
      WRITE(*,*) 'Atom type #s set for amoeba for structure ',STN
      ENDIF
!
      RETURN
      END SUBROUTINE  MOL_ATYPNADD
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      SUBROUTINE ATYPTNKMOD(STN)
!
      USE structure
      USE const
      USE specify
      USE potential 
!
      IMPLICIT none
!
      INTEGER :: STN,I,Io,Jo,To,C
!
      Io = TNKNUMBM(1,STN)
      Jo = TNKNUMBM(2,STN)
      To = TNKNUMBM(3,STN)

      C = 0
      DO I=Io,Jo
        TYPN(I,STN) = To + C 
        C = C+1
      ENDDO
!
      RETURN 
      END SUBROUTINE ATYPTNKMOD
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      SUBROUTINE MTYPTNKMOD(STN)
!
      USE structure
      USE const
      USE specify
      USE potential 
!
      IMPLICIT none
!
      INTEGER :: STN,STNo,I,MNi,C,Ni,MA,Mf,Mi,Io
!
      STNo = STN - STOTo
      MNi = MTNKNUMBM(1,STNo)
      C = MTNKNUMBM(2,STNo)

      CALL BUILDNBL(STN)
      CALL MOLCHECK(STN)

      Mi = MPNT(MNi,STN)
      Mf = MPNT(MNi+1,STN)-1           
      MA = Mf-Mi+1
!     Find the lowest atom #
      Io = 1000000000
      DO Ni=Mi,Mf 
        I = MOLST(Ni,STN)
        IF ( I .LT. Io ) Io = I
      ENDDO 
 
!
!     Verbose output
! 
      IF( VERB ) THEN
       WRITE(*,*) ' Tinker reference numbers '
     & ,' will be modified starting from ',C
     & ,' for molecule ',MNi,' with ',MA,' atoms'
     & ,' starting from atom ',Io
      ENDIF

      DO I=Io,Io+MA-1
        TYPN(I,STN) =  C 
        C = C + 1
      ENDDO
     
!
      RETURN 
      END SUBROUTINE MTYPTNKMOD
