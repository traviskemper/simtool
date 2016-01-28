!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 10/09/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Read turbomol file
!
      SUBROUTINE READ_TBML(STN)
!
      USE specify
      USE structure
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: istat,rstat,I,Elr,STN
      REAL*8 :: XX,YY,ZZ,XXa,YYa,ZZa
      CHARACTER(IDSZ) :: AT
      CHARACTER(CHSZ) :: RLINE,LLINE 
!
      I = 0 
!
      OPEN (UNIT=11,FILE=INTBML(STN),STATUS='unknown')     
      READ(11,*,iostat=istat) RLINE


      DO WHILE (ISTAT.EQ.0)
        LLINE = ADJUSTL(RLINE)
        IF( LLINE.EQ.'$coord' ) THEN 
           READ(11,*,iostat=rstat) XX,YY,ZZ,AT


           DO WHILE ( RSTAT.EQ.0 )
              I = I + 1
              XXa = XX / ANGTOBOR
              YYa = YY / ANGTOBOR
              ZZa = ZZ / ANGTOBOR
              CALL CASCVT(AT,IDSZ)
              CALL getanumb(AT,ELr)
              ELN(I,STN) = ELr
              R0(:,I,STN) = (/XXa,YYa,ZZa/)
              READ(11,*,iostat=rstat) XX,YY,ZZ,AT
           ENDDO            
           BACKSPACE(11)
           READ(11,*,iostat=istat) RLINE      
           LLINE = ADJUSTL(RLINE)

        ELSE
          READ(11,*,iostat=istat) RLINE      


        ENDIF
      ENDDO
      NA(STN) = I
!
      IF( VERB) THEN
        WRITE(*,*) I,' atoms read in from ',INTBML(STN)
      ENDIF
!

      CLOSE(11)
!
      END SUBROUTINE READ_TBML 


      SUBROUTINE cascvt(strng,lngt)
C CASe ConVerTer
C subroutine to convert all upper case to lower case
C strng: input and output string 
C lngt: number of characters in the string to convert

      USE const

      INTEGER lngt,i
      CHARACTER(IDSZ) strng
      INTEGER jgap

      jgap=ICHAR('A')-ICHAR('a')
C this could be replaced by jgap=32, since the ASCII code for "a" is
C 97 and for "A" is 65, using the ICHAR function ensures robustness
C in case other systems are used (like EBCDIC--hey I'm an old timer)

      IF(lngt.GT.0) THEN
C  if string is empty, nothing needs to be changed
        DO i=1,lngt
         IF(strng(i:i).LE.'z') THEN
           IF(strng(i:i).GE.'a') THEN
            strng(i:i)=CHAR(ICHAR(strng(i:i))+jgap)
           ENDIF
          ENDIF
        ENDDO
      ENDIF
!
      RETURN
      END
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 1.0 11/09/2011 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Write turbomol file
!
      SUBROUTINE WRITE_TBML(STN)
!
      USE specify
      USE structure
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: STN

      RETURN 
      END SUBROUTINE WRITE_TBML
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Write turbomol file PEECM 
!
      SUBROUTINE WRITE_PEECM(STN,MN)
!
      USE specify
      USE structure
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: STN,MN,MIP,Mi,Mf,N,EL,I,D,MNp
     & ,NQM,NMM,NFX,NPER
      REAL*8 :: R(NDIM),Q
      CHARACTER(CHSZ) :: FQM,FMM,CONT,CORD
      CHARACTER(IDSZ) :: AT
!
      WRITE(FQM,565) MN
      WRITE(FMM,575) MN
      WRITE(CONT,585) MN
      WRITE(CORD,595) MN
!
!       Verbose
!
        IF( VERB ) WRITE(6,*) ' Writing files for molecule:',MN
!
!
!       Print QM 
!
        NQM =  ONIONQM(STN)
        NMM =  ONIONMM(STN)
        NFX =  ONIONFX(STN)
!
!       Open xyz for qm region
!
        OPEN(UNIT=56,FILE=FQM,STATUS='UNKNOWN')
        WRITE(56,*) NQM*81
!       add efiel and e potential to comment line
        WRITE(56,5603) MOLEF(:,MN,STN),MOLEP(MN,STN)
!
!       Open embed section of control file 
!
        OPEN(UNIT=58,FILE=CONT,STATUS='UNKNOWN')
        WRITE(58,5801) 
!       find the number of pcbs used 
        NPER = 0
          DO D = 1,NDIM
            IF( PBCS(NDIM,STN) ) NPER = NPER + 1
          ENDDO
          WRITE(58,5802)  !NPER
          WRITE(58,5803) 
          WRITE(58,5804) LC(:,STN),LA(:,STN)
          WRITE(58,5805)
          DO MIp = 1,NQM
            MNp = QMLIST(MIp,STN)
            Mi = MPNT(MNp,STN)
            Mf = MPNT(MNp+1,STN)-1          
            DO N=Mi,Mf
              I = MOLST(N,STN) 
              EL = ELN(I,STN)
              AT = ATSYM(EL)
              R(:) = R0(:,I,STN)
              WRITE(56,5601) AT,R(:)
              AT = GRID(I,STN)
              R(:) = R(:)*ANGTOBOR
              WRITE(58,5601) AT,R(:)
             ENDDO
          ENDDO
          CLOSE(56)
          OPEN(UNIT=57,FILE=FMM,STATUS='UNKNOWN')
          WRITE(57,*) NMM*81 + NFX*81
          WRITE(57,*)
          DO MIp = 1,NMM
            MNp = MMLIST(MIp,STN)
            Mi = MPNT(MNp,STN)
            Mf = MPNT(MNp+1,STN)-1          
            DO N=Mi,Mf
              I = MOLST(N,STN) 
              EL = ELN(I,STN)
              AT = ATSYM(EL)
              R(:) = R0(:,I,STN)
              WRITE(57,5601) AT,R(:)
              AT = GRID(I,STN)
              R(:) = R(:)*ANGTOBOR
              WRITE(58,5601) AT,R(:)
            ENDDO
          ENDDO
          DO MIp = 1,NFX
            MNp = FXLIST(MIp,STN)
            Mi = MPNT(MNp,STN)
            Mf = MPNT(MNp+1,STN)-1          
            DO N=Mi,Mf
              I = MOLST(N,STN) 
              EL = ELN(I,STN)
              AT = ATSYM(EL)
              R(:) = R0(:,I,STN)
              WRITE(57,5601) AT,R(:)
              AT = GRID(I,STN)
              R(:) = R(:)*ANGTOBOR
              WRITE(58,5601) AT,R(:)
            ENDDO
          ENDDO
          WRITE(58,*)"end" ! end content
          CLOSE(57)
          WRITE(58,*) "cluster"
          OPEN(UNIT=59,FILE=CORD,STATUS='UNKNOWN')
          WRITE(59,*) "$coord"
          DO MIp = 1,NQM
            MNp = QMLIST(MIp,STN)
            Mi = MPNT(MNp,STN)
            Mf = MPNT(MNp+1,STN)-1          
            DO N=Mi,Mf
              I = MOLST(N,STN) 
              AT = GRID(I,STN)
              R(:) = R0(:,I,STN)
              R(:) = R(:)*ANGTOBOR
              WRITE(58,5601) AT,R(:)
              EL = ELN(I,STN)
              AT = ATSYM(EL)          
              WRITE(59,5602) R(:),AT
             ENDDO
          ENDDO       
          WRITE(59,*) "$end"
          CLOSE(59)
          WRITE(58,*)"end" ! end cluster
 
! debug hack
          WRITE(58,*) "charges"
            Mi = MPNT(MN,STN)
            Mf = MPNT(MN+1,STN)-1          
            DO N=Mi,Mf
              I = MOLST(N,STN) 
              AT = GRID(I,STN)
              Q = ACHG(I,STN)
              WRITE(58,601) AT,Q
             ENDDO
         WRITE(58,*)"end" ! end charges
      WRITE(58,*)"$end" ! end control file
      CLOSE(58)
! 
      RETURN 
 5801 FORMAT("$embed")
 5802 FORMAT(" periodic  2")
 5803 FORMAT(" cell ang")
 5804 FORMAT(6F10.4)
 5805 FORMAT("content")
 601  FORMAT(A6,F24.6)
 5601 FORMAT(A6,3F24.6)
 5602 FORMAT(3F24.6,A6)
 5603 FORMAT("Efield and V",4E16.4)
 561  FORMAT('qm_',I1,'.xyz')
 571  FORMAT('mm_',I1,'.xyz')
 562  FORMAT('qm_',I2,'.xyz')
 572  FORMAT('mm_',I2,'.xyz')
 563  FORMAT('qm_',I3,'.xyz')
 573  FORMAT('mm_',I3,'.xyz')
 564  FORMAT('qm_',I4,'.xyz')
 574  FORMAT('mm_',I4,'.xyz')
 565  FORMAT(I8,'_qm.xyz')
 575  FORMAT(I8,'_mm.xyz')
 581  FORMAT('embed_',I1)
 582  FORMAT('embed_',I2)
 583  FORMAT('embed_',I3)
 584  FORMAT('embed_',I4)
 585  FORMAT(I8,'_embed.xyz')
 591  FORMAT('coord_',I1)
 592  FORMAT('coord_',I2)
 593  FORMAT('coord_',I3)
 594  FORMAT('coord_',I4)
 595  FORMAT(I8,'_coord.xyz')
      END SUBROUTINE WRITE_PEECM
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Write turbomol file PEECM 
!
      SUBROUTINE WRITE_TBEMBED(STN,MN)
!
      USE specify
      USE structure
      USE elements 
!
      IMPLICIT none
!
      INTEGER :: STN,MN,NQM,NMM,NFX,Mi,Mf,N,EL,I,D,MNp,MIp
      REAL*8 :: R(NDIM),Q
      CHARACTER(CHSZ) :: FQM,FMM,CONT,CORD
      CHARACTER(IDSZ) :: AT
!
!       Print QM 
!
        NQM =  ONIONQM(STN)
        NMM =  ONIONMM(STN)
        NFX =  ONIONFX(STN)
!
!
!
      IF( MN.LT.10 ) THEN
        WRITE(FQM,5651) MN   ! qm.xyz 
        WRITE(FMM,5751) MN   ! mm.xyz
        WRITE(CONT,5851) MN  ! embed
        WRITE(CORD,5951) MN  ! coord
      ELSEIF( MN.LT.100 ) THEN
        WRITE(FQM,5652) MN   ! qm.xyz 
        WRITE(FMM,5752) MN   ! mm.xyz
        WRITE(CONT,5852) MN  ! embed
        WRITE(CORD,5952) MN  ! coord
      ELSEIF( MN.LT.1000 ) THEN
        WRITE(FQM,5653) MN   ! qm.xyz 
        WRITE(FMM,5753) MN   ! mm.xyz
        WRITE(CONT,5853) MN  ! embed
        WRITE(CORD,5953) MN  ! coord
      ELSEIF( MN.LT.10000 ) THEN
        WRITE(FQM,5654) MN   ! qm.xyz 
        WRITE(FMM,5754) MN   ! mm.xyz
        WRITE(CONT,5854) MN  ! embed
        WRITE(CORD,5954) MN  ! coord
      ENDIF
!
!       Verbose
!
        IF( VERB ) WRITE(6,*) ' Writing files for molecule:',MN
!
!
!       Open xyz for qm region
!
        OPEN(UNIT=56,FILE=FQM,STATUS='UNKNOWN')
        WRITE(56,*) NQM*81
!       add efiel and e potential to comment line
        WRITE(56,5603) MOLEF(:,MN,STN),MOLEP(MN,STN)
!
!       Print out quantum region 
!

          OPEN(UNIT=59,FILE=CORD,STATUS='UNKNOWN')
          WRITE(59,*) "$coord"

          DO MIp = 1,NQM
            MNp = QMLIST(MIp,STN)
            Mi = MPNT(MNp,STN)
            Mf = MPNT(MNp+1,STN)-1          
            DO N=Mi,Mf
              I = MOLST(N,STN) 
              EL = ELN(I,STN)
              AT = ATSYM(EL)
              R(:) = R0(:,I,STN)
              WRITE(56,5601) AT,R(:)
              R(:) = R(:)*ANGTOBOR
              WRITE(59,5602) R(:),AT
             ENDDO
          ENDDO
          WRITE(59,*) "$end"
          CLOSE(59)
          CLOSE(56)



          OPEN(UNIT=57,FILE=FMM,STATUS='UNKNOWN')
          WRITE(57,*) NMM*81 + NFX*81
          WRITE(57,*)

          OPEN(UNIT=58,FILE=CONT,STATUS='UNKNOWN')
          WRITE(58,5801) 

          DO MIp = 1,NMM
            MNp = MMLIST(MIp,STN)
            Mi = MPNT(MNp,STN)
            Mf = MPNT(MNp+1,STN)-1          
            DO N=Mi,Mf
              I = MOLST(N,STN)
              EL = ELN(I,STN)
              AT = ATSYM(EL)
              R(:) = R0(:,I,STN)
              WRITE(57,5601) AT,R(:)
              AT = GRID(I,STN)
              Q = ACHG(I,STN)
              R(:) = R(:)*ANGTOBOR
              WRITE(58,5802) R(:),Q
            ENDDO
          ENDDO
          DO MIp = 1,NFX
            MNp = FXLIST(MIp,STN)
            Mi = MPNT(MNp,STN)
            Mf = MPNT(MNp+1,STN)-1          
            DO N=Mi,Mf
              I = MOLST(N,STN) 
              EL = ELN(I,STN)
              AT = ATSYM(EL)
              R(:) = R0(:,I,STN)
              WRITE(57,5601) AT,R(:)
              AT = GRID(I,STN)
              Q = ACHG(I,STN)
              R(:) = R(:)*ANGTOBOR
              WRITE(58,5802) R(:),Q
            ENDDO
          ENDDO
          CLOSE(57)
          CLOSE(58)


 5651 FORMAT(I1,'_qm.xyz')
 5751 FORMAT(I1,'_mm.xyz')
 5851 FORMAT(I1,'_embed.xyz')
 5951 FORMAT(I1,'_coord.xyz')
 5652 FORMAT(I2,'_qm.xyz')
 5752 FORMAT(I2,'_mm.xyz')
 5852 FORMAT(I2,'_embed.xyz')
 5952 FORMAT(I2,'_coord.xyz')
 5653 FORMAT(I3,'_qm.xyz')
 5753 FORMAT(I3,'_mm.xyz')
 5853 FORMAT(I3,'_embed.xyz')
 5953 FORMAT(I3,'_coord.xyz')
 5654 FORMAT(I4,'_qm.xyz')
 5754 FORMAT(I4,'_mm.xyz')
 5854 FORMAT(I4,'_embed.xyz')
 5954 FORMAT(I4,'_coord.xyz')
 5601 FORMAT(A6,3F24.6)
 5602 FORMAT(3F24.6,A6)
 5603 FORMAT("Efield and V",4E16.4)
 5801 FORMAT("$point_charges")
 5802 FORMAT(3F24.6,F24.8)
      END SUBROUTINE WRITE_TBEMBED
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!





