!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Take read in file  o and minimize each structure 
! 
      SUBROUTINE MOL_MIN(STN)
!
      USE specify
      USE structure
      USE elements 

      IMPLICIT none
!
      INTEGER :: STN,MN,Mi,Mf,N,I,Ni,Nf,STNi,ACNT,MTOT
      CHARACTER(CHSZ) :: OUTGROf 
       

!
!     Set some local variables
!
      MANUMB = 81 
      MSLET = .TRUE.
      MTOT = MOLCNT(STN)
      STNi = STN +  STOTo    ! Use final structure place holder 
!
!     Verbose
!
      IF( VERB) THEN
        WRITE(*,*) 'Minimizing each molecule in ',STN
        WRITE(*,*) 'Using ',STNi,' as  temp'
        IF( MSLET ) THEN
         WRITE(*,*) 'Only molecule with ',MANUMB,' will be minimized'
        ENDIF
      ENDIF
 
!
!     Store the global 
!
      OUTGROf = OUTGRO 
      OUTGRO = 'mol.gro'       
      INGRO(STNi) = 'min.gro' 
!
!     Structure 
      LV(:,:,STNi) = LV(:,:,STN)
      LC(:,STNi) = LC(:,STN)
      LA(:,STNi) = LA(:,STN)
!
      DO MN=1,MTOT
!     
        Mi = MPNT(MN,STN)
        Mf = MPNT(MN+1,STN)-1  
        ACNT = Mf - Mi + 1
        IF( MSLET.AND.ACNT.EQ.MANUMB ) THEN
!
!         Verbose
!
          IF( VERB) WRITE(*,*) 'Minimizing ',MN,' with ',ACNT,'atoms'
!
!         Store each molecule in STNi
!
          ACNT = 0
          DO N=Mi,Mf
            ACNT = ACNT + 1
            I = MOLST(N,STN)
            CALL ATPASS(I,STN,ACNT,STNi)
          ENDDO
          NA(STNi) = ACNT 
!
!         Write out gromacs input file 
!
          CALL WRITE_GRO(STNi)
!
!         Run gromacs to minimize structure 
!
          CALL RUN_GRO(STNi)
! 
!         Read in optimized gromacs structure 
! 
          NA(STNi) = 0
          CALL READ_GRO(STNi)
!
!         Pass optimized cordinates to STN
!
          ACNT = 0
          DO N=Mi,Mf
            ACNT = ACNT + 1
            I = MOLST(N,STN)
            CALL ATPASS(ACNT,STNi,I,STN)
          ENDDO
        ENDIF 
      ENDDO
!
!     Reset global output file 
!
      OUTGRO = OUTGROf
!
      END SUBROUTINE MOL_MIN

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      
! 
      SUBROUTINE RUN_GRO(STN)
!
      IMPLICIT none
!
      INTEGER :: STN
!     
      CALL SYSTEM('sh -c "./min.sh "')
!
      END SUBROUTINE RUN_GRO
