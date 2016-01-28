!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Print error statements
       SUBROUTINE prnerror(ER,ESTAT)
       USE const
       USE structure
       USE specify
       USE elements  
!
       Implicit none
!
       INTEGER :: ER,S
       CHARACTER(CHSZ) :: RLINE,ESTAT
! 
       IF( ER.EQ.1 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment verbose "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " verbose  " 
         STOP   
       ELSEIF( ER.EQ.2 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment help "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " help  " 
         STOP   
       ELSEIF( ER.EQ.3 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment debug "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " debug  " 
         STOP   
       ELSEIF( ER.EQ.4 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment inxyz "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " inxyz (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.5 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment inxmol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " inxmol  (structure #) "
         STOP   
       ELSEIF( ER.EQ.6 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment inpdb "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " inpdb (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.7 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment ingro "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " ingro  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.8 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment incrd "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " incrd (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.9 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment incar "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " incar (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.10 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment intnk "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " intnk  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.11 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment inlammps "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "inlammps (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.12 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment inturbomol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "inturbomol (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.13 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outxyz "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outxyz (file name)"
         STOP   
       ELSEIF( ER.EQ.14 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outpdb "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outpdb  (file name)"
         STOP   
       ELSEIF( ER.EQ.15 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outcrd "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outcrd (file name)"
         STOP   
       ELSEIF( ER.EQ.16 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outgro "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outgro (file name)"
         STOP   
       ELSEIF( ER.EQ.17 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outtop "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outtop (file name)"
         STOP   
       ELSEIF( ER.EQ.18 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outcar "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outcar  (file name)"
         STOP   
       ELSEIF( ER.EQ.19 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outcom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "outcom  (file)"
         STOP   
       ELSEIF( ER.EQ.20 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outtnk "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outtnk (file)"
         STOP   
       ELSEIF( ER.EQ.21 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outmcoms "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "outmcoms  (file)"
         STOP   
       ELSEIF( ER.EQ.22 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outlammps "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "outlammps  (file)"
         STOP   
       ELSEIF( ER.EQ.23 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outturbomol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "outturbomol  (file)"
         STOP   
       ELSEIF( ER.EQ.24 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment copyregion "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*)  "copyregion  (structure #) (fraction shift) "
         STOP   
       ELSEIF( ER.EQ.25 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment rmatomtype "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "Error in control statment rmatomtype"
         STOP   
       ELSEIF( ER.EQ.26 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment vecangle "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "vecanglefile  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.27 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment vectriangle "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "depthvectrianglefile  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.28 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment atomrdftypes "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "atomrdftypes  (structure #) (atom name 1)"
         STOP   
       ELSEIF( ER.EQ.29 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment molrdfanumbs "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "molrdfanumbs  (structure #) "  
         STOP   
       ELSEIF( ER.EQ.30 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment addsurfatom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "Error in control statment addsurfatom"
         STOP   
       ELSEIF( ER.EQ.31 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment stack "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " fixstack  (structure #)"
         STOP   
       ELSEIF( ER.EQ.32 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment indexprop "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " indexprop (structure #) (index #) " 
         STOP   
       ELSEIF( ER.EQ.33 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment moleculecom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " moleculecom (structure #) (molecule #) " 
         STOP   
       ELSEIF( ER.EQ.34 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment atom_angle_dist "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " atom_angle_dist (structure #) (Type a/b/c) " 
         STOP   
       ELSEIF( ER.EQ.35 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment verbose "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " verbose  " 
         STOP   
       ELSEIF( ER.EQ.36 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment inxyz "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " inxyz (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.37 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment inpdb "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " inpdb (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.38 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment incrd "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " incrd (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.39 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment ingro "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " ingro  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.40 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment incar "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " incar (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.41 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment intop "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "Error in control statment intop"
         STOP   
       ELSEIF( ER.EQ.42 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment intnk "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " intnk  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.43 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment tnkoplsaa "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " tnkoplsaa  (structure #) "
         STOP   
       ELSEIF( ER.EQ.44 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment tnkamoeba "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " tnkamoeba  (structure #) "
         STOP   
       ELSEIF( ER.EQ.45 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment tnkadd "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " tnkadd  (structure #)  (intial atom #)"
         STOP   
       ELSEIF( ER.EQ.46 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment tnkmoladd "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " tnkmoladd  (structure #)  (intial atom #)"
         STOP   
       ELSEIF( ER.EQ.47 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment intpa "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " intpa (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.48 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outxyz "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outxyz (file name)"
         STOP   
       ELSEIF( ER.EQ.49 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outtnk "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outtnk (file)"
         STOP   
       ELSEIF( ER.EQ.50 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outpdb "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outpdb  (file name)"
         STOP   
       ELSEIF( ER.EQ.51 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outcrd "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outcrd (file name)"
         STOP   
       ELSEIF( ER.EQ.52 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outgro "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outgro (file name)"
         STOP   
       ELSEIF( ER.EQ.53 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outtop "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outtop (file name)"
         STOP   
       ELSEIF( ER.EQ.54 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outcar "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " outcar  (file name)"
         STOP   
       ELSEIF( ER.EQ.55 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outcom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "outcom  (file)"
         STOP   
       ELSEIF( ER.EQ.56 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outmcoms "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "outmcoms  (file)"
         STOP   
       ELSEIF( ER.EQ.57 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment inlammps "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "inlammps (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.58 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment inturbomol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "inturbomol (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.59 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outlammps "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "outlammps  (file)"
         STOP   
       ELSEIF( ER.EQ.60 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment outturbomol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "outturbomol  (file)"
         STOP   
       ELSEIF( ER.EQ.61 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment inxmol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " inxmol  (structure #) "
         STOP   
       ELSEIF( ER.EQ.62 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment set_ttog "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " set_ttog (structure #) " 
         STOP   
       ELSEIF( ER.EQ.63 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment seed "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " seed (integer)"
         STOP   
       ELSEIF( ER.EQ.64 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment randomguess "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " randomguess (integer)"
         STOP   
       ELSEIF( ER.EQ.65 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment random "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*)" random molecule generation in ref to (ref str#)"
         STOP   
       ELSEIF( ER.EQ.66 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment fixlatvec "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " fixlatvec (structure #)"
         STOP   
       ELSEIF( ER.EQ.67 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment rotate "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " rotate (structure #)"
         STOP   
       ELSEIF( ER.EQ.68 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment expandrand "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "Error in control statment expandrand"
         STOP   
       ELSEIF( ER.EQ.69 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment latconst "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
           WRITE(*,*) " latconst S# a b c alpha beta gamma "
         STOP   
       ELSEIF( ER.EQ.70 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment latvec "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " fixlatvec (structure #)"
         STOP   
       ELSEIF( ER.EQ.71 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment center "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "molcenter  (structure #) (mol #) "
         STOP   
       ELSEIF( ER.EQ.72 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment vacuum "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " vacuum  (structure #) (box size)"
         STOP   
       ELSEIF( ER.EQ.73 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment velocity "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " velocity (structure #)  (vector)"
         STOP   
       ELSEIF( ER.EQ.74 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment supercell "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " supercell (structure #) (3 real #)"
         STOP   
       ELSEIF( ER.EQ.75 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment mass2atnumb "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "mass2atnumb  (structure #) "
         STOP   
       ELSEIF( ER.EQ.76 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment molcheck "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "Error in control statment molcheck"
         STOP   
       ELSEIF( ER.EQ.77 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment molcmcenter "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "molcmcenter  (structure #) (molecule #) (x/y/z)"
         STOP   
       ELSEIF( ER.EQ.78 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment molmxmncenter "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*)"molmxmncenter (structure #) (molecule #) (x/y/z)"
         STOP   
       ELSEIF( ER.EQ.79 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment cutregion "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "Error in control statment cutregion "
         STOP   
       ELSEIF( ER.EQ.80 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment cut_mol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "cut_mol  (structure #) (mol i) (mol j)"
         STOP   
       ELSEIF( ER.EQ.81 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment pbc "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
            WRITE(6,*) "Use pbc to control which are used  "
         STOP   
       ELSEIF( ER.EQ.82 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment fixregion "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "fixregion   (structure #)  "
         STOP   
       ELSEIF( ER.EQ.83 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment molrdfanumbs "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "molrdfanumbs  (structure #) "  
         STOP   
       ELSEIF( ER.EQ.84 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment molrdffile "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "molrdffile  (structure #)  (mol #)"
         STOP   
       ELSEIF( ER.EQ.85 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment vecangle "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "vecanglefile  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.86 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment vecanglefile "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "vecanglefile  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.87 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment depthvecanglefile "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "depthvecanglefile  (structure #) (ID)"
         STOP   
       ELSEIF( ER.EQ.88 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment vectriangle "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "depthvectrianglefile  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.89 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment vectrianglefile "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "depthvectrianglefile  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.90 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment depthvectrianglefile "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "depthvectrianglefile  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.91 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment copyregion "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*)  "copyregion  (structure #) (fraction shift) "
         STOP   
       ELSEIF( ER.EQ.92 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment atomrdftypes "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "atomrdftypes  (structure #) (atom name 1)"
         STOP   
       ELSEIF( ER.EQ.93 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment rdfdistfile "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "rdfdistfile  (structure #) (file name)"
         STOP   
       ELSEIF( ER.EQ.94 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment slab "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " slab  (structure #) (dimension 1-x,2-y,3-z)"
         STOP   
       ELSEIF( ER.EQ.95 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment rename "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " rename (structure #) (atom name) (residue)"
         STOP   
       ELSEIF( ER.EQ.96 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment nbcharges "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " nbcharges (structure #) "
         STOP   
       ELSEIF( ER.EQ.97 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment addsurfatom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) "Error in control statment addsurfatom"
         STOP   
       ELSEIF( ER.EQ.98 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment stack "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " fixstack  (structure #)"
         STOP   
       ELSEIF( ER.EQ.99 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment hterm "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " hterm (structure #) "
         STOP   
       ELSEIF( ER.EQ.100 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment fixstack "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " fixstack  (structure #)"
         STOP   
       ELSEIF( ER.EQ.101 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment totaldipole "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*)  "totaldipole  (structure #)  "
         STOP   
       ELSEIF( ER.EQ.102 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment rotateall "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*)"rotateall (structure #) (min angle) (max angle)"
         STOP   
       ELSEIF( ER.EQ.103 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment moldipole "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " moldipole (structure #) (file name) " 
         STOP   
       ELSEIF( ER.EQ.104 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment moldipdist "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " moldipdist (structure #) (file name) " 
         STOP   
       ELSEIF( ER.EQ.105 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment moldipaxdist "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*)"moldipaxdist (structure #) "
     &   ,"(dimensions #) (file name)"
         STOP   
       ELSEIF( ER.EQ.106 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment moldipbin "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " moldipbin (structure #) (scale A) (max ) " 
         STOP   
       ELSEIF( ER.EQ.107 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment indexprop "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " indexprop (structure #) (index #) " 
         STOP   
       ELSEIF( ER.EQ.108 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment nonblist "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " nonblist  " 
         STOP   
       ELSEIF( ER.EQ.109 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment zeromax "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " zeromax (structure #) (dimension 1=x,2=y,3=z)" 
         STOP   
       ELSEIF( ER.EQ.110 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment shift "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*)  "copyregion  (structure #) (fraction shift)"
         STOP   
       ELSEIF( ER.EQ.111 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment allmolcom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " allmolcom (structure #) " 
         STOP   
       ELSEIF( ER.EQ.112 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment moleculecom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " moleculecom (structure #) (molecule #) " 
         STOP   
       ELSEIF( ER.EQ.113 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment qmsphere "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " qmsphere (structure #) (inner radius) "
     & ,"(outer raduis) (center x/y/z) " 
         STOP   
       ELSEIF( ER.EQ.114 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment qmnieghbors "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*)"qmnieghbors (structure #) (# of neighbor shells)"
         STOP   
       ELSEIF( ER.EQ.115 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment atomefield "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " atomefield (structure #) " 
         STOP   
       ELSEIF( ER.EQ.116 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment molefield "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " molefield (structure #) (file name) " 
         STOP   
       ELSEIF( ER.EQ.117 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment delR "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " delR (structure #) (scale ) " 
         STOP   
       ELSEIF( ER.EQ.118 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment deldip "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " deldip (structure #) (dscale) (dmax) " 
         STOP   
       ELSEIF( ER.EQ.119 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment noqmmmfixed "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " noqmmmfixed (structure #) " 
         STOP   
       ELSEIF( ER.EQ.120 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment g_molmin "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " g_molmin (structure #) " 
         STOP   
       ELSEIF( ER.EQ.121 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment update_charges "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " update_charges (structure #) " 
         STOP   
       ELSEIF( ER.EQ.122 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment xyz_trj "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " xyz_trj (structure #) " 
         STOP   
       ELSEIF( ER.EQ.123 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment atom_angle_dist "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " atom_angle_dist (structure #) (Type a/b/c) " 
         STOP   
       ELSEIF( ER.EQ.124 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment atom_angle_file "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " atom_angle_file (structure #) (file) " 
         STOP   
       ELSEIF( ER.EQ.125 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment cov_buffer "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " cov_buffer (buffer) " 
         STOP   
       ELSEIF( ER.EQ.126 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment read_bconst "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " read_bconst (structure #) " 
         STOP   
       ELSEIF( ER.EQ.127 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment read_bconst "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " read_bconst (structure #) " 
         STOP   
       ELSEIF( ER.EQ.128 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment read_aconst "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " read_aconst (structure #) " 
         STOP   
       ELSEIF( ER.EQ.129 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment read_dconst "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " read_dconst (structure #) " 
         STOP   
       ELSEIF( ER.EQ.130 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment read_iconst "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " read_iconst (structure #) " 
         STOP   
       ELSEIF( ER.EQ.131 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment deposit "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " deposit (structure #) (# of molecules) "
     & ,"(dimension) (buffer distance) " 
         STOP   
       ELSEIF( ER.EQ.132 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment tnknumbmod "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " tnknumbmod (structure #) (atom min) "
     & ,"(atom max) (intial #) " 
         STOP   
       ELSEIF( ER.EQ.133 ) THEN  
!        Error in reading in control states
         BACKSPACE(5) 
         READ(5,'(A)') RLINE
         WRITE(6,*) "Error in control statment tnkmolmod "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " tnkmolmod (structure #) (mol # min) (mol # max)" 
         STOP   
      ENDIF
      
      END SUBROUTINE prnerror

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
