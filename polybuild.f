!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  
!
      SUBROUTINE polybuild
!
      IMPLICIT none
!
      INTEGER :: STNo 

      DO STNo = 1,STOTo
        IF( MKPOLY(STNo) ) THEN
          !
          ! Analyize input structure as monomer unit 
          !           
          CALL MAXMIN(STNi)
          CALL MOLCHECK(STNo)
          IF( MOLCNT(STNo).NE.1 ) THEN
           WRITE(*,*) 'Warning more than one molecule in monomer'
          ENDIF 
          MERBOX(:,STNo) = MN(2,:,STNo) - MMN(1,:,STNo)
          !
          ! Verbose output 
          !
          IF( VERB ) THEN
             WRITE(6,*) ' Genorating polymer based on strucutre',STNo
             WRITE(6,*) ' Monomer unit :'
             WRITE(6,*)  MERBOX(1,STNo)
             WRITE(6,*)  MERBOX(2,STNo)
             WRITE(6,*)  MERBOX(3,STNo)
          ENDIF 
          !
          ! Genorate structure 
          !
          IF( POLY_TYPE(STNo) .EQ. 1) CALL POLY_GRID(STNo)
        ENDIF 
      ENDDO ! STNo 
      
      END SUBROUTINE polybuild

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  
!
      SUBROUTINE POLY_GRID(STNo)
!
      IMPLICIT none
!
      INTEGER :: STNo 


      END SUBROUTINE POLY_GRID


      LOGICAL :: POLY_BUILD 
      ALLOCATABLE,  REAL*8 :: MERBOX(:,:)
      ALLOCATABLE, LOGICAL :: MKPOLY(:)





