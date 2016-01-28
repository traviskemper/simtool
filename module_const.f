!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 10/21/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Module of constants 
!
      MODULE const
!
      REAL*8 :: PI,BOLZ,AVO,ECONV,EPSI,DCON,EPSN,ELCH,ANGTOBOR
      PARAMETER ( PI=3.141592654d0
     &  ,BOLZ=1.380662d0
     &  ,AVO=6.02201415d0
     &  ,EPSI=11604.50D0  
     &  ,ECONV=(1.0D0/(AVO*BOLZ*EPSI))*1.0D+07
     & ,DCON = 4.803204191d0                                            ! Debye/(eA)  
     & ,EPSN = 8.85d0
     & ,ELCH =  1.60217646d0
     & ,ANGTOBOR = 1.889726133d0
     &  ,BOLZeV = 8.6173430d-5
     &  )
      INTEGER :: CHSZ,NDIM,IDSZ,ERSZ
      PARAMETER ( CHSZ = 70 )
      PARAMETER ( NDIM = 3 )
      PARAMETER ( IDSZ = 5 )             !for pdb and gro formats 
      PARAMETER ( ERSZ = 12 )             !for pdb and gro formats 
!
      END MODULE const 
