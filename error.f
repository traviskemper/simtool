ELSEIF( ER.EQ.1 ) THEN  
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
         WRITE(6,*) " indexprop (structure #) (index) " 
         STOP   
ELSEIF( ER.EQ.1 ) THEN  
!        Error in reading in control states
         BACKSPACE         BACKSREAD(5,         BACKSPACE    WRITE(6,*) "Error in control statment verbose "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " verbose  " 
         STOP   
ELSEIF( ER.EQ.2 ) THEN  
!        Error in reading in control states
         BACKSPACE         BACKSREAD(5,         BACKSPACE    WRITE(6,*) "Error in control statment help "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " help  " 
         STOP   
ELSEIF( ER.EQ.3 ) THEN  
!        Error in reading in control states
         BACKSPACE         BACKSREAD(5,         BACKSPACE    WRITE(6,*) "Error in control statment debug "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " debug  " 
         STOP   
ELSEIF( ER.EQ.4 ) THEN  
!        Error in reading in control states
         BACKSPACE         BACKSREAD(5,         BACKSPACE    WRITE(6,*) "Error in control statment inxyz "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " inxyz (structure #) (files #) " 
         STOP   
ELSEIF( ER.EQ.1 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment verbose "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " verbose  " 
         STOP   
ELSEIF( ER.EQ.2 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment help "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " help  " 
         STOP   
ELSEIF( ER.EQ.3 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment debug "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " debug  " 
         STOP   
ELSEIF( ER.EQ.4 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment inxyz "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " inxyz (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.5 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment inxmol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " inxmol  (structure #) "
         STOP   
ELSEIF( ER.EQ.6 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment inpdb "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " inpdb (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.7 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment ingro "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " ingro  (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.8 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment incrd "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " incrd (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.9 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment incar "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " incar (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.10 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment intnk "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " intnk  (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.11 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment inlammps "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "inlammps (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.12 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment inturbomol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "inturbomol (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.13 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outxyz "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outxyz (file name)"
         STOP   
ELSEIF( ER.EQ.14 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outpdb "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outpdb  (file name)"
         STOP   
ELSEIF( ER.EQ.15 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outcrd "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outcrd (file name)"
         STOP   
ELSEIF( ER.EQ.16 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outgro "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outgro (file name)"
         STOP   
ELSEIF( ER.EQ.17 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outtop "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outtop (file name)"
         STOP   
ELSEIF( ER.EQ.18 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outcar "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outcar  (file name)"
         STOP   
ELSEIF( ER.EQ.19 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outcom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "outcom  (file)"
         STOP   
ELSEIF( ER.EQ.20 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outtnk "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outtnk (file)"
         STOP   
ELSEIF( ER.EQ.21 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outmcoms "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "outmcoms  (file)"
         STOP   
ELSEIF( ER.EQ.22 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outlammps "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "outlammps  (file)"
         STOP   
ELSEIF( ER.EQ.23 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment outturbomol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "outturbomol  (file)"
         STOP   
ELSEIF( ER.EQ.24 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment copyregion "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*)  "copyregion  (structure #) (fraction shift) "
         STOP   
ELSEIF( ER.EQ.25 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment rmatomtype "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "Error in control statment rmatomtype"
         STOP   
ELSEIF( ER.EQ.26 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment vecangle "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "vecanglefile  (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.27 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment vectriangle "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "depthvectrianglefile  (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.28 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment atomrdftypes "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "atomrdftypes  (structure #) (atom name 1)"
         STOP   
ELSEIF( ER.EQ.29 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment molrdfanumbs "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "molrdfanumbs  (structure #) "  
         STOP   
ELSEIF( ER.EQ.30 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment addsurfatom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "Error in control statment addsurfatom"
         STOP   
ELSEIF( ER.EQ.31 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment stack "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " fixstack  (structure #)"
         STOP   
ELSEIF( ER.EQ.32 ) THEN  
!        Error in reading in control states
                  (5) 
         READ(         READ(          WRITE(6,*) "Error in control statment indexprop "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " indexprop (structure #) (index) " 
         STOP   
ELSEIF( ER.EQ.1 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment verbose "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " verbose  " 
         STOP   
ELSEIF( ER.EQ.2 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment help "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " help  " 
         STOP   
ELSEIF( ER.EQ.3 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment debug "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
         WRITE(6,*) " debug  " 
         STOP   
ELSEIF( ER.EQ.4 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment inxyz "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " inxyz (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.5 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment inxmol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " inxmol  (structure #) "
         STOP   
ELSEIF( ER.EQ.6 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment inpdb "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " inpdb (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.7 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment ingro "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " ingro  (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.8 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment incrd "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " incrd (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.9 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment incar "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " incar (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.10 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment intnk "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " intnk  (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.11 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment inlammps "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "inlammps (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.12 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment inturbomol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "inturbomol (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.13 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outxyz "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outxyz (file name)"
         STOP   
ELSEIF( ER.EQ.14 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outpdb "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outpdb  (file name)"
         STOP   
ELSEIF( ER.EQ.15 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outcrd "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outcrd (file name)"
         STOP   
ELSEIF( ER.EQ.16 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outgro "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outgro (file name)"
         STOP   
ELSEIF( ER.EQ.17 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outtop "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outtop (file name)"
         STOP   
ELSEIF( ER.EQ.18 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outcar "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outcar  (file name)"
         STOP   
ELSEIF( ER.EQ.19 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outcom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "outcom  (file)"
         STOP   
ELSEIF( ER.EQ.20 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outtnk "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " outtnk (file)"
         STOP   
ELSEIF( ER.EQ.21 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outmcoms "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "outmcoms  (file)"
         STOP   
ELSEIF( ER.EQ.22 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outlammps "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "outlammps  (file)"
         STOP   
ELSEIF( ER.EQ.23 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment outturbomol "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "outturbomol  (file)"
         STOP   
ELSEIF( ER.EQ.24 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment copyregion "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*)  "copyregion  (structure #) (fraction shift) "
         STOP   
ELSEIF( ER.EQ.25 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment rmatomtype "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "Error in control statment rmatomtype"
         STOP   
ELSEIF( ER.EQ.26 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment vecangle "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "vecanglefile  (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.27 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment vectriangle "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "depthvectrianglefile  (structure #) (file name)"
         STOP   
ELSEIF( ER.EQ.28 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment atomrdftypes "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "atomrdftypes  (structure #) (atom name 1)"
         STOP   
ELSEIF( ER.EQ.29 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment molrdfanumbs "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "molrdfanumbs  (structure #) "  
         STOP   
ELSEIF( ER.EQ.30 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment addsurfatom "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) "Error in control statment addsurfatom"
         STOP   
ELSEIF( ER.EQ.31 ) THEN  
!        Error in reading in control states
                                READ(                        WRITE(6,*) "Error in control statment stack "
         WRITE(6,*) RLINE
         WRITE(6,*) "The proper format should be: "
                  WRITE(6,*) " fixstack  (structure #)"
         STOP   
