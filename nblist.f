!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!     Version 2.0 04/12/2012 T. W. Kemper                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  Construct nieghbor list 
!
      SUBROUTINE BUILDNBL(STN)
      USE const
      USE specify
      USE structure
      USE elements 
      USE potential 
!
      IMPLICIT none
!
      INTEGER ::  STN,I,NAr,NBr,ELi,ELj,NBa,NNB,BMXr,Imx
     & ,J,D,NBrLJ
      REAL*8  :: DF(NDIM),DR(NDIM),RSQ,F1(NDIM),F2(NDIM)
     &   ,LVr(NDIM,NDIM),RC,RCSQ
      CHARACTER(CHSZ) :: ESTAT
!
!     Triclinc PBC's
      NAr = NA(STN)
      LVr(:,:) = LV(:,:,STN)
      ALLOCATE( Fr(NDIM,NAr) ) 
      CALL fracR(STN,LVr,1)
!
      NBr = 0
      BMXr = 0
      Imx = 0 
      NBrLJ = 0
!
      DO I=1,NAr
!
         NBa = NBr + 1
         NINDX(I,STN) = NBr + 1
         NINDXLJ(I,STN) = NBrLJ + 1
         ELi = ELN(I,STN)
         F1(:) = Fr(:,I)
         DO J=1,NAr
           IF( I.NE.J ) THEN
             ELj = ELN(J,STN)
!            Apply PBC
             F2(:) = Fr(:,J)
             DO D=1,NDIM
               DF(D) = F2(D) - F1(D)
               IF( PBCS(D,STN) )  DF(D) = DF(D) - ANINT(DF(D))
             ENDDO
             CALL REALF(LVr,DF,DR)
             RSQ = 0.0d0
             DO D=1,NDIM
               RSQ = RSQ + DR(D)*DR(D)
             ENDDO
             RC =  ( VDWRDI(Eli) + VDWRDI(ELj) )*LJBUF

             RCSQ = RC*RC
             IF(RSQ.LT.RCSQ) THEN 
               NBrLJ = NBrLJ + 1
               IF( NBrLJ.GT.NBMAXLJ(STN)) THEN
                  WRITE(*,*) ' For str',STN,"atom",I
     &           ," NBrLJ ",NBrLJ,"gt",NBMAXLJ(STN)
                  CALL prnerror(-6,ESTAT) 
               ENDIF
               NLISTLJ(NBrLJ,STN) = J 
             ENDIF
             RC =  ( CRDI(Eli) + CRDI(ELj) )*CRBUF
             RCSQ = RC*RC
             IF(RSQ.LT.RCSQ) THEN 
               NBr = NBr + 1
               IF( NBr.GT.NBMAX(STN)) THEN
                  WRITE(*,*) ' For str',STN,"atom",I
     &           ," NBr ",NBr,"gt",NBMAX(STN)
                  NNB = NBr - NBa
                  IF( NNB.GT.BMXr ) THEN
                    BMXr  = NNB
                   Imx = I
                  ENDIF
                  WRITE(*,*) "MAx niehgbors", BMXr, "on",Imx
                  Eli = ELN(Imx,STN)
                  WRITE(*,*) " with cutoff of",CRDI(Eli)*CRBUF
                  WRITE(*,*) " element ",Eli," mass ",AMAS(Imx,STN)
                  CALL prnerror(-6,ESTAT) 
               ENDIF
               NLIST(NBr,STN) = J 
             ENDIF 
           ENDIF
         ENDDO 
         NNB = NBr - NBa
         IF( NNB.GT.BMXr ) THEN
           Imx = I
           BMXr = NNB
         ENDIF
      ENDDO
!
      NINDX(NAr+1,STN) =NBr+1
      NBTOT(STN) = NBr
!
      NINDXLJ(NAr+1,STN) =NBrLJ+1
      NBTOTLJ(STN) = NBrLJ
!
      DEALLOCATE( Fr ) 
      IF( VERB) THEN
        WRITE(*,*) " ",NBr," found in structure ",STN  
      ENDIF 
!
      RETURN
      END SUBROUTINE BUILDNBL
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
