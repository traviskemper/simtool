!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      Read top file 
!
      SUBROUTINE READ_TOP(STN)
      USE specify
      USE structure
      USE potential
      USE elements  
!
      IMPLICIT none

      INTEGER :: STN,Io,Ii,istat,rstat
     & ,RNr,CN,N,I,J,K,L,MN,Mi,Mf,MNi,AMN
     & ,Ia,Ib,In,Id,A,M,Mo
     & ,MCNTl,MOLCNTr,CNT,Ni,Nf 
     & ,NBl,NANGl,NDIHl,NIMPl
     & ,Jp,Kp,Lp,Jo,Ko,Lo,Ip
     & ,TP,MOLMX
      INTEGER, ALLOCATABLE :: MPNTl(:)
      REAL*8 :: CH,MA,VAL 
      CHARACTER(CHSZ) :: RLINE,LLINE,LI,var,kywd,FIL2,ESTAT
      CHARACTER(IDSZ) :: AT,RD,GT
      CHARACTER(CHSZ) :: CONSTi
      CHARACTER(3), ALLOCATABLE ::  TYPl(:) 
      CHARACTER(IDSZ), ALLOCATABLE ::  RIDl(:)
      INTEGER, ALLOCATABLE :: RNl(:)
     & ,NMB(:),NMA(:),NMD(:),NMI(:)
     & ,BNDIl(:,:),BNDJl(:,:)
     & ,BTYPl(:,:),ATYPl(:,:),DTYPl(:,:),ITYPl(:,:)
     & ,ANGIl(:,:),ANGJl(:,:),ANGKl(:,:)
     & ,DIHIl(:,:),DIHJl(:,:),DIHKl(:,:),DIHLl(:,:)
     & ,IMPIl(:,:),IMPJl(:,:),IMPKl(:,:),IMPLl(:,:)
 
      REAL*8, ALLOCATABLE :: ACHGl(:),CHNl(:),AMASl(:) 
     & ,BVALl(:,:),AVALl(:,:),DVALl(:,:),IVALl(:,:)
      CHARACTER(IDSZ), ALLOCATABLE ::  GRIDl(:)
      CHARACTER(CHSZ), ALLOCATABLE :: BCNSTl(:,:)
     & ,ACNSTl(:,:),DCNSTl(:,:),ICNSTl(:,:)
      LOGICAL :: STRN,COMNT
!
!     Initialize 
      Ia = 0
      Ib = 0
      In = 0 
      Id = 0
      MCNTl = 0 
      CNT = 0 
      MOLCNTr = 0 
!
!     Get number of molecules to read in 
!
      MOLMX = NMOL(STN) 

! hack 
!     MOLMX = 2

!
!     Allocate local arrays 
!
      ALLOCATE( MPNTl(NAMX) 
     & ,TYPl(NAMX),RIDl(NAMX),RNl(NAMX),ACHGl(NAMX) 
     & ,CHNl(NAMX),GRIDl(NAMX),AMASl(NAMX)         
     & ,NMB(12*NAMX),NMA(12*NAMX),NMD(12*NAMX),NMI(12*NAMX)
     & ,BNDIl(12*NAMX,MOLMX),BNDJl(12*NAMX,MOLMX)
     & ,BTYPl(12*NAMX,MOLMX),BVALl(12*NAMX,MOLMX),BCNSTl(12*NAMX,MOLMX)
     & ,ANGIl(12*NAMX,MOLMX),ANGJl(12*NAMX,MOLMX),ANGKl(12*NAMX,MOLMX)
     & ,ATYPl(12*NAMX,MOLMX),AVALl(12*NAMX,MOLMX),ACNSTl(12*NAMX,MOLMX)
     & ,DIHIl(12*NAMX,MOLMX),DIHJl(12*NAMX,MOLMX)
     & ,DIHKl(12*NAMX,MOLMX),DIHLl(12*NAMX,MOLMX)
     & ,DTYPl(12*NAMX,MOLMX)
     & ,IMPIl(12*NAMX,MOLMX),IMPJl(12*NAMX,MOLMX)
     & ,IMPKl(12*NAMX,MOLMX),IMPLl(12*NAMX,MOLMX)
     & )
!
!      Intial global  count 
!
       NB(STN) = 0
       NANG(STN) = 0
       NDIH(STN) = 0
       NIMP(STN) = 0
!
!      Intial local molecule count 
!
       NMB(:) = 0
       NMA(:) = 0
       NMD(:) = 0
       NMI(:) = 0

!
      OPEN (UNIT=16,FILE=INTOP(STN),STATUS='unknown')     
      READ(16,*,iostat=istat) RLINE
      BACKSPACE(16)
      DO WHILE (ISTAT.EQ.0)
        LLINE = ADJUSTL(RLINE)
        LI = LLINE(1:1) 
        IF( LLINE.EQ.'#include' ) THEN 
            READ(16,*,iostat=rstat) var,FIL2
            IF(VERB)WRITE(*,*) 'Opeing sub top file  :',FIL2
            OPEN (UNIT=17,FILE=FIL2,STATUS='unknown')    
            READ(17,*,iostat=istat) RLINE
            BACKSPACE(17)
            DO WHILE (ISTAT.EQ.0)
                 LLINE = ADJUSTL(RLINE)
                 LI = LLINE(1:1) 
                 IF ( LI.EQ.'[' ) THEN
                    Io=SCAN( RLINE,'[' )      
                    Ii=SCAN( RLINE,']')      
                    IF ( Ii.LT.Io ) THEN
                       READ(17,*) var,KYWD 
                    ELSE 
                       KYWD = LLINE(Io+1:Ii-1)
                       READ(17,*)
                    ENDIF
                    IF( KYWD.EQ.'moleculetype' ) THEN
                      MCNTl = MCNTl + 1     
                      MPNTl(MCNTl) = Ia + 1
                    ELSEIF ( KYWD.EQ.'atoms' ) THEN          
                       READ(17,*,iostat=rstat) N,AT,RNr,RD,GT,CN,CH,MA
                       IF( rstat.NE.0) THEN
                         BACKSPACE(17) 
                         READ(17,*) AT
                         IF(AT.EQ.";") THEN
                         READ(17,*,iostat=rstat) N,AT,RNr,RD,GT,CN,CH,MA
                         ENDIF 
                       ENDIF 
                       DO WHILE ( RSTAT.EQ.0 )
                          Ia = Ia + 1
                          TYPl(Ia)  = AT
                          RNl(Ia)   = RNr
                          RIDl(Ia)  = RD
                          GRIDl(Ia) = GT         
                          CHNl(Ia)  = CN
                          ACHGl(Ia) = CH
                          AMASl(Ia) = MA 
                         READ(17,*,iostat=rstat) N,AT,RNr,RD,GT,CN,CH,MA
                       ENDDO         
                       BACKSPACE(17)   
                    ELSEIF( KYWD.EQ.'bonds' ) THEN
                      Ib = NMB( MCNTl )
!
!                     Read in line 
!
                      IF( RDBCNST(STN) ) THEN 
                          READ(17,*,iostat=rstat) I,J,TP,VAL,CONSTi
                                ELSE 
                          READ(17,*,iostat=rstat) I,J,TP
                            VAL = -1.0d0
                            CONSTi = '1.0' 
                      ENDIF 
!
!                     Check for comment line 
!
                      IF( rstat.NE.0) THEN
                        BACKSPACE(17) 
                        COMNT = .TRUE.
                        DO WHILE ( COMNT ) 
                          READ(17,*) AT
                          IF( AT(1:1) .NE. ";" )  COMNT = .FALSE. 
                        ENDDO 
!
!                       If comment line then read in the next line 
!
                        BACKSPACE(17) 
                        IF( RDBCNST(STN) ) THEN 
                          READ(17,*,iostat=rstat) I,J,TP,VAL,CONSTi
                        ELSE 
                          READ(17,*,iostat=rstat) I,J,TP
                          VAL = -1.0d0
                          CONSTi = '1.0' 
                        ENDIF 
                      ENDIF 

!
!                      
!        
                      DO WHILE ( RSTAT.EQ.0 )
                        Ib = Ib + 1
                        IF( I.GT.NAMX .OR. J.GT.NAMX ) THEN
                             WRITE(ESTAT,*) I," ",J," gt  ",NAMX
                             CALL PRNERROR(-10,ESTAT)
                        ENDIF              
                        IF( I.LT.0 .OR. J.LT.0 ) THEN
                             WRITE(ESTAT,*) I," ",J," lt  0"
                             CALL PRNERROR(-10,ESTAT)
                        ENDIF        
                        BNDIl(Ib,MCNTl) = I 
                        BNDJl(Ib,MCNTl) = J
!
!                       Save type , value and constant 
!
                        BTYPl(Ib,MCNTl) = TP
                        BVALl(Ib,MCNTl) = VAL
                        BCNSTl(Ib,MCNTl) = CONSTi 

!
!                       Read in line 
!
                        IF( RDBCNST(STN) ) THEN 
                          READ(17,*,iostat=rstat) I,J,TP,VAL,CONSTi
                                ELSE 
                          READ(17,*,iostat=rstat) I,J,TP
                            VAL = -1.0d0
                            CONSTi = '1.0' 
                        ENDIF 
!
!                       Check for comment line 
!
                        IF( rstat.NE.0) THEN
                          BACKSPACE(17) 
                          COMNT = .TRUE.
                          DO WHILE ( COMNT ) 
                            READ(17,*) AT
                            IF( AT(1:1) .NE. ";" )  COMNT = .FALSE. 
                          ENDDO 
!
!                         If comment line then read in the next line 
! 
                          BACKSPACE(17) 
                          IF( RDBCNST(STN) ) THEN 
                            READ(17,*,iostat=rstat) I,J,TP,VAL,CONSTi
                          ELSE 
                            READ(17,*,iostat=rstat) I,J,TP
                            VAL = -1.0d0
                            CONSTi = '1.0' 
                          ENDIF 
                        ENDIF 

                      ENDDO           
                      NMB( MCNTl ) = Ib 
                      BACKSPACE(17)   



                    ELSEIF( KYWD.EQ.'angles' ) THEN
                      In = NMA(  MCNTl )
!
!                     Read in line 
!         
                      IF( RDACNST(STN) ) THEN 
                         READ(17,*,iostat=rstat) K,I,J,TP,VAL,CONSTi
                      ELSE 
                         READ(17,*,iostat=rstat) K,I,J,TP
                         VAL = -1.0d0
                         CONSTi = '1.0' 
                      ENDIF 

!
!                     Check for comment line 
!
                      IF( rstat.NE.0) THEN
                           BACKSPACE(17) 
                           COMNT = .TRUE.
                           DO WHILE ( COMNT ) 
                             READ(17,*) AT
                             IF( AT(1:1) .NE. ";" )  COMNT = .FALSE. 
                           ENDDO 
!
!                          If comment line then read in the next line 
!
                           BACKSPACE(17) 
                           IF( RDACNST(STN) ) THEN 
                             READ(17,*,iostat=rstat) K,I,J,TP,VAL,CONSTi
                           ELSE 
                             READ(17,*,iostat=rstat) K,I,J,TP
                             VAL = -1.0d0
                             CONSTi = '1.0' 
                           ENDIF 
                      ENDIF 
!         
!                     Read to end of bond section 
!

                      DO WHILE ( RSTAT.EQ.0 )
                        In = In + 1
                        ANGIl(In,MCNTl) = I 
                        ANGJl(In,MCNTl) = J
                        ANGKl(In,MCNTl) = K 
!
!                       Save type , value and constant 
!
                        ATYPl(In,MCNTl) = TP
                        AVALl(In,MCNTl) = VAL
                        ACNSTl(In,MCNTl) = CONSTi 
!
!
!                       Read in line 
!         
                        IF( RDACNST(STN) ) THEN 
                           READ(17,*,iostat=rstat) K,I,J,TP,VAL,CONSTi
                        ELSE 
                           READ(17,*,iostat=rstat) K,I,J,TP
                           VAL = -1.0d0
                           CONSTi = '1.0' 
                        ENDIF 

!
!                       Check for comment line 
!
                        IF( rstat.NE.0) THEN
                           BACKSPACE(17) 
                           COMNT = .TRUE.
                           DO WHILE ( COMNT ) 
                             READ(17,*) AT
                             IF( AT(1:1) .NE. ";" )  COMNT = .FALSE. 
                           ENDDO 
!
!                          If comment line then read in the next line 
!
                           BACKSPACE(17) 
                           IF( RDACNST(STN) ) THEN 
                             READ(17,*,iostat=rstat) K,I,J,TP,VAL,CONSTi
                           ELSE 
                             READ(17,*,iostat=rstat) K,I,J,TP
                             VAL = -1.0d0
                             CONSTi = '1.0' 
                           ENDIF 
                        ENDIF 
                      ENDDO           
                      NMA( MCNTl ) = In 
                      BACKSPACE(17)   
                    ELSEIF( KYWD.EQ.'dihedrals' ) THEN
                       Id = NMD(  MCNTl )
                       READ(17,*,iostat=rstat) K,I,J,L,TP
                       DO WHILE ( RSTAT.EQ.0 )
                        Id = Id + 1
                        DIHIl(Id,MCNTl) = I 
                        DIHJl(Id,MCNTl) = J
                        DIHKl(Id,MCNTl) = K 
                        DIHLl(Id,MCNTl) = L

!                       Save type , value and constant 
!
                        DTYPl(Id,MCNTl) = TP
   !                     DVALl(Id,MCNTl) = VAL
   !                     DCNSTl(Id,MCNTl) = CONSTi 

                          READ(17,*,iostat=rstat) K,I,J,L,TP
                       ENDDO           
                       NMD( MCNTl ) = Id
                       BACKSPACE(17)   
                    ELSEIF( KYWD.EQ.'constraints' ) THEN
!     
                    ELSEIF( KYWD.EQ.'exclusions' ) THEN
!     
                    ELSEIF( KYWD.EQ.'virtual' ) THEN
!
                   ENDIF 
                 ELSE
                    READ(17,*)
                 ENDIF
                 READ(17,*,iostat=istat) RLINE
                 BACKSPACE(17)
              ENDDO
              CLOSE(17)
              IF( VERB ) THEN
                   WRITE(*,*) " Finished reading ", FIL2
                   WRITE(*,*) "   molecule count ",MCNTl
                    WRITE(*,*) "   atom count ",Ia
                    WRITE(*,*) "   bond count ",Ib
                    WRITE(*,*) "   angle count ",In
                    WRITE(*,*) "   dih count ",Id
                    WRITE(*,*) "   imp count not counted yet"
              ENDIF
        ELSEIF ( LI.EQ.'[' ) THEN
          Io=SCAN( RLINE,'[' )      
          Ii=SCAN( RLINE,']')      
          IF ( Ii.LT.Io ) THEN
             READ(16,*) var,KYWD 
          ELSE 
             KYWD = LLINE(Io+1:Ii-1)
             READ(16,*)
          ENDIF
          IF( KYWD.EQ.'moleculetype' ) THEN
            MCNTl = MCNTl + 1     
            MPNTl(MCNTl) = Ia + 1
          ELSEIF ( KYWD.EQ.'atoms' ) THEN          
            READ(16,*,iostat=rstat) N,AT,RNr,RD,GT,CN,CH,MA
            DO WHILE ( RSTAT.EQ.0 )
              Ia = Ia + 1
              TYPl(Ia)  = AT
              RNl(Ia)   = RNr
              RIDl(Ia)  = RD
              GRIDl(Ia) = GT         
              CHNl(Ia)  = CN
              ACHGl(Ia) = CH
              AMASl(Ia) = MA 
               READ(16,*,iostat=rstat) N,AT,RNr,RD,GT,CN,CH,MA
            ENDDO        
            BACKSPACE(16)   
              IF( VERB ) THEN
                   WRITE(*,*) " Finished reading atoms section "
                   WRITE(*,*) "   molecule count ",MCNTl
                    WRITE(*,*) "   atom count ",Ia
              ENDIF
          ELSEIF( KYWD.EQ.'bonds' ) THEN
              Ib = NMB( MCNTl )
!
!             Read in line 
!
              IF( RDBCNST(STN) ) THEN 
                READ(16,*,iostat=rstat) I,J,TP,VAL,CONSTi
              ELSE 
                READ(16,*,iostat=rstat) I,J,TP
                VAL = -1.0d0
                CONSTi = '1.0' 
              ENDIF 

!
!             Check for comment line 
!
              IF( rstat.NE.0) THEN
                  BACKSPACE(16) 
                  COMNT = .TRUE.
                  DO WHILE ( COMNT ) 
                    READ(16,*) AT
                    IF( AT(1:1) .NE. ";" )  COMNT = .FALSE. 
                  ENDDO 
!
!                 If comment line then read in the next line 
!
                  BACKSPACE(16) 
                  IF( RDBCNST(STN) ) THEN 
                    READ(16,*,iostat=rstat) I,J,TP,VAL,CONSTi
                  ELSE 
                    READ(16,*,iostat=rstat) I,J,TP
                    VAL = -1.0d0
                    CONSTi = '1.0' 
                  ENDIF 
              ENDIF 
!
!             Read to end of bond section 
!
              DO WHILE ( RSTAT.EQ.0 )
                Ib = Ib + 1
                IF( I.GT.NAMX .OR. J.GT.NAMX ) THEN
                   WRITE(ESTAT,*) I," ",J," gt  ",NAMX
                   CALL PRNERROR(-10,ESTAT)
                ENDIF              
                IF( I.LT.0 .OR. J.LT.0 ) THEN
                   WRITE(ESTAT,*) I," ",J," lt  0"
                   CALL PRNERROR(-10,ESTAT)
                ENDIF              
                BNDIl(Ib,MCNTl) = I 
                BNDJl(Ib,MCNTl) = J
!
!               Save type , value and constant 
!
                BTYPl(Ib,MCNTl) = TP
                BVALl(Ib,MCNTl) = VAL
                BCNSTl(Ib,MCNTl) = CONSTi 
!
!
!               Read in line 
!
                IF( RDBCNST(STN) ) THEN 
                  READ(16,*,iostat=rstat) I,J,TP,VAL,CONSTi
                ELSE 
                  READ(16,*,iostat=rstat) I,J,TP
                    VAL = -1.0d0
                    CONSTi = '1.0' 
                ENDIF 
!
!               Check for comment line 
!
                IF( rstat.NE.0) THEN
                  BACKSPACE(16) 
                  COMNT = .TRUE.
                  DO WHILE ( COMNT ) 
                    READ(16,*) AT
                    ! skip comment lines 
                    IF( AT(1:1) .NE. ";" )   COMNT = .FALSE.
                  ENDDO 
!
!                 If comment line then read in the next line 
!

                  BACKSPACE(16) 
                  IF( RDBCNST(STN) ) THEN 
                    READ(16,*,iostat=rstat) I,J,TP,VAL,CONSTi
                  ELSE 
                    READ(16,*,iostat=rstat) I,J,TP
                    VAL = -1.0d0
                    CONSTi = '1.0' 
                  ENDIF 
                ENDIF 
!
              ENDDO           
              NMB( MCNTl ) = Ib 
            BACKSPACE(16)   
          ELSEIF( KYWD.EQ.'angles' ) THEN
             In = NMA(  MCNTl )
!
!             Read in line 
!
              IF( RDACNST(STN) ) THEN 
                READ(16,*,iostat=rstat) K,I,J,TP,VAL,CONSTi
              ELSE 
                READ(16,*,iostat=rstat) K,I,J,TP
                VAL = -1.0d0
                CONSTi = '1.0' 
              ENDIF 

!
!             Check for comment line 
!
              IF( rstat.NE.0) THEN
                  BACKSPACE(16) 
                  COMNT = .TRUE.
                  DO WHILE ( COMNT ) 
                    READ(16,*) AT
                    IF( AT(1:1) .NE. ";" )  COMNT = .FALSE. 
                  ENDDO 
!
!                 If comment line then read in the next line 
!
                  BACKSPACE(16) 
                  IF( RDACNST(STN) ) THEN 
                    READ(16,*,iostat=rstat) K,I,J,TP,VAL,CONSTi
                  ELSE 
                    READ(16,*,iostat=rstat) K,I,J,TP
                    VAL = -1.0d0
                    CONSTi = '1.0' 
                  ENDIF 
              ENDIF 
!
!             Read to end of bond section 
!

             DO WHILE ( RSTAT.EQ.0 )
                In = In + 1
                ANGIl(In,MCNTl) = I 
                ANGJl(In,MCNTl) = J
                ANGKl(In,MCNTl) = K
!
!               Save type , value and constant 
!
                ATYPl(In,MCNTl) = TP
                AVALl(In,MCNTl) = VAL
                ACNSTl(In,MCNTl) = CONSTi 
!
!               Read in line 
!
                IF( RDACNST(STN) ) THEN 
                    READ(16,*,iostat=rstat) K,I,J,TP,VAL,CONSTi
                ELSE 
                    READ(16,*,iostat=rstat) K,I,J,TP
                    VAL = -1.0d0
                    CONSTi = '1.0' 
                ENDIF 
!
!               Check for comment line 
!
                IF( rstat.NE.0) THEN
                  BACKSPACE(16) 
                  COMNT = .TRUE.
                  DO WHILE ( COMNT ) 
                    READ(16,*) AT
                    ! skip comment lines 
                    IF( AT(1:1) .NE. ";" )   COMNT = .FALSE.
                  ENDDO 
!
!                 If comment line then read in the next line 
!

                  BACKSPACE(16) 
                  IF( RDACNST(STN) ) THEN 
                    READ(16,*,iostat=rstat) K,I,J,TP,VAL,CONSTi
                  ELSE 
                    READ(16,*,iostat=rstat) K,I,J,TP
                    VAL = -1.0d0
                    CONSTi = '1.0' 
                  ENDIF 
                ENDIF 
 
             ENDDO           
             NMA( MCNTl ) = In 
             BACKSPACE(16)   
         ELSEIF( KYWD.EQ.'dihedrals' ) THEN
             Id = NMD(  MCNTl )
             READ(16,*,iostat=rstat) K,I,J,L,TP
             DO WHILE ( RSTAT.EQ.0 )
                Id = Id + 1
                DIHIl(Id,MCNTl) = I 
                DIHJl(Id,MCNTl) = J
                DIHKl(Id,MCNTl) = K 
                DIHLl(Id,MCNTl) = L

!               Save type , value and constant 
!
                DTYPl(Id,MCNTl) = TP

                READ(16,*,iostat=rstat) K,I,J,L,TP 
             ENDDO           
             NMD( MCNTl ) = Id
             BACKSPACE(16)   
          ELSEIF( KYWD.EQ.'constraints' ) THEN
!
          ELSEIF( KYWD.EQ.'exclusions' ) THEN
! 
          ELSEIF( KYWD.EQ.'virtual' ) THEN
!
          ELSEIF( KYWD.EQ.'molecules' ) THEN
            MPNTl(MCNTl+1) = Ia + 1
            CNT = 0 
            MOLCNTr = 0 
            MCNTl = 0 
            READ(16,*,iostat=rstat) var,AMN 
            DO WHILE ( RSTAT.EQ.0 )
              MCNTl = MCNTl + 1
              Ni = MPNTl(MCNTl)
              Nf = MPNTl(MCNTl+1) - 1 
              IF(VERB) THEN
                WRITE(*,*) "Replicating mol",MCNTl,AMN
                WRITE(*,*) "   atom count ",Nf - Ni + 1
                WRITE(*,*) "   atom numbers ",Nf,Ni 
                WRITE(*,*) "   Bonds numbers ",NMB( MCNTl )
                WRITE(*,*) "   Angles numbers ",NMA( MCNTl )
                WRITE(*,*) "   Dih numbers ",NMD( MCNTl )
                WRITE(*,*) "   Imp numbers ",NMI( MCNTl )
              ENDIF
!
!             Loop over local molecule list 
!
              DO MN=1,AMN 
                  MOLCNTr = MOLCNTr + 1
                  MPNT(MOLCNTr,STN) = CNT + 1
                  I = 0
                  DO Io =Ni,Nf
                    I = I + 1
                    CNT = CNT + 1
                    IF( CNT .GT. NAMX ) CALL prnerror(-4,ESTAT)
                    MOLST(CNT,STN) = CNT
                    MOLN(CNT,STN) = MOLCNTr
                    TYP(CNT,STN)  = TYPl(Io)           ! Force field type
                    RID(CNT,STN)  = RIDl(Io)           ! Residue ID
                    RN(CNT,STN)   = RNl(Io)            ! Residue #
                    ACHG(CNT,STN) = ACHGl(Io)          ! Charge #
                    CHN(CNT,STN)  = CHNl(Io)           ! Charge group #
                    GRID(CNT,STN) = GRIDl(Io)          ! GROMACS ID
                    AMAS(CNT,STN) = AMASl(Io)          ! Atomic mass
                    REND_IND(I,STN) = CNT
                 ENDDO 
!
!                 Add molecule bonds to global 
!              
                  NBl = NB(STN) 
                  DO Ib =1,NMB( MCNTl )
                     NBl = NBl + 1
                     Io = BNDIl(Ib,MCNTl)
                     Jo = BNDJl(Ib,MCNTl)
                     TP  = BTYPl(Ib,MCNTl)
                     VAL = BVALl(Ib,MCNTl)
                     CONSTi = BCNSTl(Ib,MCNTl)
   
                     IF( Io.GT.NAMX .OR. Jo.GT.NAMX ) THEN
                            WRITE(ESTAT,*) " Io ",Io," ",Jo," gt  ",NAMX
                            CALL PRNERROR(-10,ESTAT)
                     ENDIF
                     Ip = REND_IND(Io,STN)
                     Jp = REND_IND(Jo,STN)

                     IF( Ip.GT.NAMX .OR. Jp.GT.NAMX ) THEN
                            WRITE(ESTAT,*) " P ",Ip," ",Jp," gt  ",NAMX
                            CALL PRNERROR(-10,ESTAT)
                     ENDIF
!
!                    Save to global values 
!
                     BNDI(NBl,STN) = Ip
                     BNDJ(NBl,STN) = Jp
                     BTYP(NBl,STN) = TP
                     BVAL(NBl,STN) = VAL
                     BCNST(NBl,STN) = CONSTi 
                  ENDDO
                  NB(STN) = NBl
!    
!                  Add molecule angles to global 
!              
                  NANGl = NANG(STN) 
                  DO Ib =1,NMA( MCNTl )
                     NANGl = NANGl + 1
                     Io = ANGIl(Ib,MCNTl)
                     Jo = ANGJl(Ib,MCNTl)
                     Ko = ANGKl(Ib,MCNTl)
                     TP  = ATYPl(Ib,MCNTl)
                     VAL = AVALl(Ib,MCNTl)
                     CONSTi = ACNSTl(Ib,MCNTl)

                     Ip = REND_IND(Io,STN)
                     Jp = REND_IND(Jo,STN)
                     Kp = REND_IND(Ko,STN)

                     ANGI(NANGl,STN) = Ip
                     ANGJ(NANGl,STN) = Jp
                     ANGK(NANGl,STN) = Kp
                     ATYP(NANGl,STN) = TP
                     AVAL(NANGl,STN) = VAL
                     ACNST(NANGl,STN) = CONSTi 

                  ENDDO
                  NANG(STN) = NANGl
!
!                  Add molecule dih to global 
!              
                  NDIHl = NDIH(STN) 
                  DO Ib =1,NMD( MCNTl )
                     NDIHl = NDIHl + 1
                     Io = DIHIl(Ib,MCNTl)
                     Jo = DIHJl(Ib,MCNTl)
                     Ko = DIHKl(Ib,MCNTl)
                     Lo = DIHLl(Ib,MCNTl)
                     TP  = DTYPl(Ib,MCNTl)

                     Ip = REND_IND(Io,STN)
                     Jp = REND_IND(Jo,STN)
                     Kp = REND_IND(Ko,STN)
                     Lp = REND_IND(Lo,STN)
                     DIHI(NDIHl,STN) = Ip
                     DIHJ(NDIHl,STN) = Jp 
                     DIHK(NDIHl,STN) = Kp 
                     DIHL(NDIHl,STN) = Lp 

                     DTYP(NDIHl,STN) = TP

                  ENDDO
                  NDIH(STN) = NDIHl
!
!                  Add molecule imp to global 
!              
                  NIMPl = NIMP(STN) 
                  DO Ib =1,NMI( MCNTl )
                     NIMPl = NIMPl + 1
                     Io = IMPIl(Ib,MCNTl)
                     Jo = IMPJl(Ib,MCNTl)
                     Ko = IMPKl(Ib,MCNTl)
                     Lo = IMPLl(Ib,MCNTl)
                     Ip = REND_IND(Io,STN)
                     Jp = REND_IND(Jo,STN)
                     Kp = REND_IND(Ko,STN)
                     Lp = REND_IND(Lo,STN)
                     IMPI(NIMPl,STN) = Ip
                     IMPJ(NIMPl,STN) = Jp 
                     IMPK(NIMPl,STN) = Kp 
                     IMPL(NIMPl,STN) = Lp 
                  ENDDO
                  NIMP(STN) = NIMPl
!
              ENDDO 
              READ(16,*,iostat=rstat) var,AMN 
            ENDDO           
            BACKSPACE(16)   
            MOLCNT(STN) = MOLCNTr
            MPNT(MOLCNTr+1,STN) = CNT + 1

      IF(DEBUG) THEN
         WRITE(204,*) MOLCNTr,CNT
         WRITE(204,*) MPNT(MOLCNTr,STN) , MPNT(MOLCNTr + 1,STN) 
      ENDIF

 !
          ENDIF 
       ELSE
         READ(16,*)
       ENDIF
       READ(16,*,iostat=istat) RLINE
       BACKSPACE(16)
      ENDDO
      CLOSE(16)

      IF(VERB) THEN
        WRITE(*,*) " Molecules read in",MOLCNTr
        WRITE(*,*) " Atoms read in",CNT 
        WRITE(*,*) " Bonds read in",NB(STN)
        WRITE(*,*) " Angles read in",NANG(STN)
        WRITE(*,*) " Dih read in",NDIH(STN)
        WRITE(*,*) " Imp read in (not working!!)",NIMP(STN)
      ENDIF

!
      DEALLOCATE( MPNTl
     & ,TYPl,RIDl,RNl,ACHGl 
     & ,CHNl,GRIDl,AMASl         
     & ,NMB,NMA,NMD,NMI
     & ,BNDIl,BNDJl,ANGIl,ANGJl,ANGKl
     & ,DIHIl,DIHJl,DIHKl,DIHLl
     & ,IMPIl,IMPJl,IMPKl,IMPLl
     & ,BTYPl,BVALl,BCNSTl
     & ,ATYPl,AVALl,ACNSTl
     & ,DTYPl
     & )
!
      RETURN
! 161  FORMAT(2I,I,F,A20)

      END SUBROUTINE READ_TOP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 06/13/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Write .top file for GROMACS input
!
      SUBROUTINE WRITE_TOP(STN)
      USE specify
      USE structure
      USE potential 
      USE elements 
!
      IMPLICIT none
!
      INTEGER   :: I,J,N,EL,D,STN,RNB,CN,Ii,Ji,Ki,Li
     & ,TP
      REAL*8 :: CH,MS,VAL 
      CHARACTER(CHSZ) :: NBID(12),CONSTi
      CHARACTER(5) :: RD,AT,AN
!
      WRITE(6,*) 'Writing xyz coordinates',OUTXYZ
      OPEN(UNIT=53,FILE=OUTTOP,STATUS='unknown')

      WRITE(53,'(A)') ''
      WRITE(53,'(A)') '[ atoms ]'
!
!     Loop over all atoms in structure 
!
      DO I=1,NA(STN)
        RD = ADJUSTL( RID(I,STN))                 ! Residue ID
        AT = ADJUSTR( TYP(I,STN) )         ! Force field type
        AN = ADJUSTR( GRID(I,STN) )        ! GROMACS ID 
        EL = ELN(I,STN)        ! Atomic #
        RNB = RN(I,STN)        ! Residue #
        CN = CHN(I,STN)        ! Charge group #
        CH = ACHG(I,STN)       ! Charge 
        MS = AMASS(EL)         ! Mass 
! hack debug
          RNB = MOLN(I,STN)

!       Atom #,
        WRITE(53,5303) I,AT,RNB,RD,AN,CN,CH,MS
      ENDDO


!
!       bonds 
!    
      WRITE(53,'(A)') ''
      WRITE(53,'(A)') '[ bonds ]'
      IF(VERB) THEN
           WRITE(*,*) 'Render  ',NB(STN) ,' bonds '
      ENDIF
      DO I=1,NB(STN)
        Ii = BNDI(I,STN)
        Ji = BNDJ(I,STN)
        TP  = BTYP(I,STN)
        VAL = BVAL(I,STN)
        CONSTi = BCNST(I,STN)
        IF ( VAL .GT. 0 ) THEN
          WRITE(53,5305) Ii,Ji,TP,VAL,CONSTi 
        ELSE
          WRITE(53,5315) Ii,Ji,TP
        ENDIF
      ENDDO

!
!     angles 
!    
      WRITE(53,'(A)') ''
      WRITE(53,'(A)') '[ angles ]'
      DO I=1,NANG(STN)
        Ii = ANGI(I,STN)
        Ji = ANGJ(I,STN)   
        Ki = ANGK(I,STN)   
        TP  = ATYP(I,STN)
        VAL = AVAL(I,STN)
        CONSTi = ACNST(I,STN)
        IF ( VAL .GT. 0 ) THEN
          WRITE(53,5306) Ki,Ii,Ji,TP,VAL,CONSTi 
        ELSE
          WRITE(53,5316) Ki,Ii,Ji,TP
        ENDIF
      ENDDO
!
!         dih 
!    
      WRITE(53,'(A)') ''
      WRITE(53,'(A)') '[ dihedrals ]'
      DO I=1,NDIH(STN)
           Ii = DIHI(I,STN)
           Ji = DIHJ(I,STN)   
           Ki = DIHK(I,STN)   
           Li = DIHL(I,STN)   
           TP  = DTYP(I,STN)
           WRITE(53,5317) Ki,Ii,Ji,Li,TP
      ENDDO
!
!       imp 
!     
      IF( NIMP(STN) .GT. 0 ) THEN
        WRITE(53,'(A)') ''
        WRITE(53,'(A)') '[ impropers ]'
        DO I=1,NIMP(STN)
          Ii = IMPI(I,STN)
          Ji = IMPJ(I,STN)   
          Ki = IMPK(I,STN)   
          Li = IMPL(I,STN)  
          TP  = ITYP(I,STN)
          WRITE(53,5317) Ii,Ji,Ki,Li,TP
        ENDDO
      ENDIF
!
!
!      
      WRITE(53,*) ''
      CLOSE(53)
!
      RETURN
 5301 FORMAT(4I8,' 2')
 5302 FORMAT(4I8,' 2')
 5303 FORMAT(I6,A8,I6,2A8,I6,F16.8,F10.3)
 5903 FORMAT(I6,A8,I6,2A8,I6,F16.8,F10.3,I6)
 5304 FORMAT(2I8,I4,3F10.4)
 5315 FORMAT(3I8)
 5305 FORMAT(3I8,F12.4,"  ",A20)
 5316 FORMAT(4I8)
 5306 FORMAT(4I8,F12.4,"  ",A20)
 5317 FORMAT(5I8)
 5307 FORMAT(5I8,F12.4,"  ",A20)
      END SUBROUTINE WRITE_TOP
c$$$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
c$$$!     Version 1.0 11/02/2011 T. W. Kemper                      !
c$$$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
c$$$!  Write .top file for GROMACS input
c$$$!
c$$$      SUBROUTINE WRITE_TOP(STN)
c$$$      USE specify
c$$$      USE structure
c$$$      USE potential 
c$$$!
c$$$      IMPLICIT none
c$$$!
c$$$      INTEGER   :: I,J,N,EL,D,STN,RN
c$$$      CHARACTER(CHSZ) :: NBID(12)
c$$$      CHARACTER(5) :: RD,AT,AN
c$$$!
c$$$      WRITE(6,*) 'Writing xyz coordinates',OUTXYZ
c$$$      OPEN(UNIT=53,FILE=OUTTOP,STATUS='unknown')
c$$$      OPEN(UNIT=59,FILE='ref.top',STATUS='unknown')
c$$$
c$$$      WRITE(53,'(A)') '; Include forcefield parameters'
c$$$      WRITE(53,'(A)') '#include "./mAlq3-13.ff/forcefield.itp"'
c$$$      WRITE(53,'(A)') '  '
c$$$      WRITE(53,'(A)') '[ moleculetype ] '
c$$$      WRITE(53,'(A)') '; Name            nrexcl '
c$$$      WRITE(53,'(A)') 'Protein             3    '
c$$$      WRITE(53,'(A)') ''
c$$$      WRITE(53,'(A)') '[ atoms ]'
c$$$      WRITE(59,'(A,I)') ' atoms ',NAi(STN)
c$$$      DO I=1,NAi(STN)
c$$$        RD = ADJUSTL( RIDi(I,STN))
c$$$        AT = ADJUSTR( ATYPi(I,STN) )
c$$$        AN = ADJUSTR( ANMi(I,STN) )
c$$$        EL = ELNi(I,STN)
c$$$        RN = RNi(I,STN)
c$$$        WRITE(53,5303) I,AT,RN,RD,AN,RNi(I,STN)
c$$$     &    ,CHRGi(I,STN),MASS(EL)
c$$$        WRITE(59,5903) I,AT,RN,RD,AN,RNi(I,STN)
c$$$     &    ,CHRGi(I,STN),MASS(EL),EL
c$$$      ENDDO
c$$$      WRITE(53,*) ''
c$$$      WRITE(59,*) ''
c$$$      IF( NCi(STN).GT.0) THEN
c$$$       WRITE(53,'(A)') '[ constraints ]'
c$$$       WRITE(59,'(A,I)') ' constraints',NCi(STN)
c$$$       DO I=1,NCi(STN)
c$$$        WRITE(53,5304)  CONIi(I,STN) ,CONJi(I,STN),CTYPi(I,STN)
c$$$     &               ,CDISIi(I,STN),CDISJi(I,STN)
c$$$        WRITE(59,5304)  CONIi(I,STN) ,CONJi(I,STN),CTYPi(I,STN)
c$$$     &               ,CDISIi(I,STN),CDISJi(I,STN)
c$$$       ENDDO 
c$$$      ENDIF
c$$$      IF( NBi(STN).GT.0) THEN
c$$$       WRITE(53,'(A)') '[ bonds ]'
c$$$       WRITE(59,'(A,I)') 'bonds ',NBi(STN)
c$$$       DO I=1,NBi(STN)
c$$$        WRITE(53,*) BNDIi(I,STN),BNDJi(I,STN),BTYPi(I,STN)
c$$$        WRITE(59,*) BNDIi(I,STN),BNDJi(I,STN),BTYPi(I,STN)
c$$$       ENDDO 
c$$$      ENDIF
c$$$      WRITE(53,*) ''
c$$$      WRITE(59,*) ''
c$$$      IF( NPRi(STN).GT.0) THEN
c$$$       WRITE(53,'(A)') '[ pairs ]'
c$$$       WRITE(59,'(A,I)') ' pairs ',NPRi(STN)
c$$$       DO I=1,NPRi(STN)
c$$$        WRITE(53,*) PRSIi(I,STN),PRSJi(I,STN),PRTYPi(I,STN)
c$$$        WRITE(59,*) PRSIi(I,STN),PRSJi(I,STN),PRTYPi(I,STN)
c$$$       ENDDO 
c$$$      ENDIF
c$$$      WRITE(53,*) ''
c$$$      WRITE(59,*) ''
c$$$      IF( NANGi(STN).GT.0) THEN
c$$$       WRITE(53,'(A)') '[ angles ]'
c$$$       WRITE(59,'(A,I)') ' angles ',NANGi(STN)
c$$$       DO I=1,NANGi(STN)
c$$$        WRITE(53,*)ANGIi(I,STN),ANGJi(I,STN),ANGKi(I,STN),ANTYPi(I,STN)
c$$$        WRITE(59,*)ANGIi(I,STN),ANGJi(I,STN),ANGKi(I,STN),ANTYPi(I,STN)
c$$$       ENDDO 
c$$$      ENDIF
c$$$      WRITE(53,*) ''
c$$$      WRITE(59,*) ''
c$$$      IF( NDIHi(STN).GT.0) THEN
c$$$       WRITE(53,'(A)') '[ dihedrals ]'
c$$$       WRITE(59,'(A,I)') ' dihedrals ',NDIHi(STN)
c$$$       DO I=1,NDIHi(STN)
c$$$        WRITE(53,*)DIHIi(I,STN),DIHJi(I,STN),DIHKi(I,STN)
c$$$     &   ,DIHLi(I,STN),DHTYPi(I,STN)
c$$$        WRITE(59,*)DIHIi(I,STN),DIHJi(I,STN),DIHKi(I,STN)
c$$$     &   ,DIHLi(I,STN),DHTYPi(I,STN)
c$$$       ENDDO 
c$$$      ENDIF
c$$$      WRITE(59,*) ''
c$$$      WRITE(53,'(A)') ' '
c$$$      WRITE(53,'(A)') '; Include Position restraint file '
c$$$!      WRITE(53,'(A)') '#ifdef POSRES '
c$$$!      WRITE(53,'(A)') '#include "posre.itp" '
c$$$!      WRITE(53,'(A)') '#endif '
c$$$      WRITE(53,'(A)') ' '
c$$$      WRITE(53,'(A)') '[ system ] '
c$$$      WRITE(53,'(A)') '; Name '
c$$$      WRITE(53,'(A)') 'Protein '
c$$$      WRITE(53,'(A)') ' ' 
c$$$      WRITE(53,'(A)') '[ molecules ] '
c$$$      WRITE(53,'(A)') '; Compound        #mols '
c$$$      WRITE(53,'(A)') 'Protein             1   '
c$$$      WRITE(53,'(A)') ' '
c$$$      WRITE(53,'(A)') ' '
c$$$      CLOSE(53)
c$$$      CLOSE(59) 
c$$$!
c$$$      RETURN
c$$$ 5301 FORMAT(4I8,' 2')
c$$$ 5302 FORMAT(4I8,' 2')
c$$$ 5303 FORMAT(I6,A8,I6,2A8,I6,F16.8,F10.3)
c$$$ 5903 FORMAT(I6,A8,I6,2A8,I6,F16.8,F10.3,I6)
c$$$ 5304 FORMAT(2I8,I4,3F10.4)
c$$$      END SUBROUTINE WRITE_TOP
