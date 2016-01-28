!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Calculate fractional coordinates from mathmatica 
!  J -> -(-BC CB X + BB CC X + BC CA Y - BA CC Y - BB CA Z + 
!       BA CB Z)/(AC BB CA - AB BC CA - AC BA CB + AA BC CB + 
!       AB BA CC - AA BB CC), 
!  L -> -(AC CB X - AB CC X - AC CA Y + AA CC Y + AB CA Z - 
!       AA CB Z)/(AC BB CA - AB BC CA - AC BA CB + AA BC CB + 
!       AB BA CC - AA BB CC), 
!  M -> -(-AC BB X + AB BC X + AC BA Y - AA BC Y - AB BA Z + 
!       AA BB Z)/(AC BB CA - AB BC CA - AC BA CB + AA BC CB + 
!       AB BA CC - AA BB CC)}}


      SUBROUTINE fracR(STN,LVr,N)
      USE const
      USE structure
      USE specify
!
      Implicit none
!
      INTEGER :: STN,I,N,NAr
      REAL*8 :: LVr(NDIM,NDIM)
     &  ,X,Y,Z,AA,AB,AC,BA,BB,BC,CA,CB,CC 
      REAL*8, ALLOCATABLE :: Rr(:,:)
      CHARACTER(12) :: ESTAT
!
!
      NAr= NA(STN)
      ALLOCATE( Rr(NDIM,NAr) ) 
      DO I=1,NAr
          Rr(:,I) = R0(:,I,STN)
      ENDDO
!     Set lattice vectors
      AA = LVr(1,1)
      AB = LVr(1,2)
      AC = LVr(1,3)
      BA = LVr(2,1)
      BB = LVr(2,2)
      BC = LVr(2,3)
      CA = LVr(3,1)
      CB = LVr(3,2)
      CC = LVr(3,3)
! 
      DO I=1,NAr
        X = Rr(1,I)
        Y = Rr(2,I)
        Z = Rr(3,I)
      Fr(1,I)=-((-BC*CB*X+BB*CC*X+BC*CA*Y-BA*CC*Y-BB*CA*Z+BA*CB*Z)/
     &     (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
      Fr(2,I)=-((AC*CB*X-AB*CC*X-AC*CA*Y+AA*CC*Y+AB*CA*Z-AA*CB*Z)/
     &     (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
      Fr(3,I)=-((-AC*BB*X+AB*BC*X+AC*BA*Y-AA*BC*Y-AB*BA*Z+AA*BB*Z)/
     &      (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
!
      ENDDO
      DEALLOCATE( Rr )
!
      RETURN
      END SUBROUTINE fracR
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 04/16/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     translate fractional coordinates to real coordinates 
!      
      SUBROUTINE REALF(LVr,F,R)
      USE structure
!
      IMPLICIT none
!
      REAL*8 :: LVr(NDIM,NDIM),F(NDIM)
     &  ,AA,AB,AC,BA,BB,BC,CA,CB,CC,R(NDIM) 
!
!     Set lattice vectors
      AA = LVr(1,1)
      AB = LVr(1,2)
      AC = LVr(1,3)
      BA = LVr(2,1)
      BB = LVr(2,2)
      BC = LVr(2,3)
      CA = LVr(3,1)
      CB = LVr(3,2)
      CC = LVr(3,3)
!
      R(1) = F(1)*AA + BA*F(2)+CA*F(3) 
      R(2) = F(1)*AB + BB*F(2)+CB*F(3)
      R(3) = F(1)*AC + BC*F(2)+CC*F(3)
!
      RETURN
      END SUBROUTINE REALF


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Find fractional coordinate of single point 
!
      SUBROUTINE SPFRAC(LVr,F,R)
      USE structure
!
      IMPLICIT none
!
      REAL*8 :: LVr(NDIM,NDIM),F(NDIM),X,Y,Z
     & ,AA,AB,AC,BA,BB,BC,CA,CB,CC,R(NDIM) 
!
!     Set lattice vectors
      AA = LVr(1,1)
      AB = LVr(1,2)
      AC = LVr(1,3)
      BA = LVr(2,1)
      BB = LVr(2,2)
      BC = LVr(2,3)
      CA = LVr(3,1)
      CB = LVr(3,2)
      CC = LVr(3,3)
!
      X=R(1)
      Y=R(2)
      Z=R(3)
!
      F(1)=-((-BC*CB*X+BB*CC*X+BC*CA*Y-BA*CC*Y-BB*CA*Z+BA*CB*Z)/
     &     (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
      F(2)=-((AC*CB*X-AB*CC*X-AC*CA*Y+AA*CC*Y+AB*CA*Z-AA*CB*Z)/
     &      (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
      F(3)=-((-AC*BB*X+AB*BC*X+AC*BA*Y-AA*BC*Y-AB*BA*Z+AA*BB*Z)/
     &      (AC*BB*CA-AB*BC*CA-AC*BA*CB+AA*BC*CB+AB*BA*CC-AA*BB*CC))
!
!
!
      RETURN
      END SUBROUTINE SPFRAC
