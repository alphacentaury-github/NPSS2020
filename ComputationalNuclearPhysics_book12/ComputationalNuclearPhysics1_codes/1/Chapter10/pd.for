      PROGRAM PD
C
C========================================== fitting  version  ========C
C                                                                     C
C     FDU0 Package                                                    C
C     ---- -------                                                    C
C                                                                     C
C     Operator Matrix Elements for the Hamiltonian                    C
C                                                                     C
C=====================================================================C
C
      IMPLICIT  REAL*8 (A-H,O-Z)
      PARAMETER (MAXN=21,MAXNS=84,MAXCFP=7762,MAXP=11221)
      PARAMETER (MAXA=18,MAXPAIR=5000,MAXMULT=5000)
      INTEGER   B,R,A,AP
      INTEGER   JSIZE(0:MAXN,0:2*MAXN),OMI
      INTEGER   CFPEND,CFPINDEX(MAXCFP),PAIREND,PAIRINDEX(MAXPAIR)
      INTEGER   MULTEND,MULTINDEX(MAXMULT)
      REAL*8    CFPTAB(MAXCFP),PTAB(MAXP),PAPER(MAXA,MAXA,0:3,0:2*MAXN)
      REAL*4    PAIRTAB(MAXPAIR),MULTTAB(MAXMULT)
      INTEGER   PEND,PINDEX(MAXP)
      LOGICAL   LTRIANGLE
      CHARACTER ACTIVE*3
      EQUIVALENCE (PAIRTAB,MULTTAB),(PAIRINDEX,MULTINDEX)
      EQUIVALENCE (PTAB,CFPTAB)
C
      INDP(AP,JP,R,A,J)=20000000*J+200000*JP+10000*R+100*A+AP
      INDCFP(A,J,R,B,L)=20000000*J+200000*A+10000*R+100*L+B
      INDPAIR(AP,J,R,A)=1000000*R+10000*J+100*AP+A
      PHASE(I)=1-2*MOD(ABS(I),2)
      LTRIANGLE(J1,J2,J3)=(ABS(J1-J2).LE.J3).AND.(J1+J2.GE.J3)
C
 203  PRINT*,'SO8 OR SP6?'
      READ (*,'(A)') ACTIVE
      IF(ACTIVE.EQ.'SO8') THEN
       OPEN(1,FILE='SO8CFP.DAT',STATUS='OLD',FORM='UNFORMATTED')
       OPEN(3,FILE='SO8P.DAT',STATUS='OLD',FORM='UNFORMATTED')
       OPEN(15,FILE='SO8JSIZE.DAT',STATUS='OLD',FORM='UNFORMATTED')
       OPEN(2,FILE='SO8PAIR.DAT',STATUS='NEW',FORM='UNFORMATTED')
       OPEN(4,FILE='SO8MULT.DAT',STATUS='NEW',FORM='UNFORMATTED')
       IA=4
      ELSE IF(ACTIVE.EQ.'SP6') THEN
       OPEN(1,FILE='SP6CFP.DAT',STATUS='OLD',FORM='UNFORMATTED')
       OPEN(3,FILE='SP6P.DAT',STATUS='OLD',FORM='UNFORMATTED')
       OPEN(15,FILE='SP6JSIZE.DAT',STATUS='OLD',FORM='UNFORMATTED')
       OPEN(2,FILE='SP6PAIR.DAT',STATUS='NEW',FORM='UNFORMATTED')
       OPEN(4,FILE='SP6MULT.DAT',STATUS='NEW',FORM='UNFORMATTED')
       IA=3
      ELSE
       PRINT*,'RE-ENTER ACTIVE MODE'
       GOTO 203
      END IF
 10   READ(15,END=100) OMI,N
      IF(OMI.GT.MAXN) THEN
       PRINT*,'OMI=',OMI
       STOP 'MAXN TOO SMALL'
      END IF
      BACKSPACE(15)
      DO N=1,OMI
       READ(15)
       READ(15) (JSIZE(N,J),J=0,2*N)
      END DO
      DO N=1,OMI
       DO J=0,2*N
        IF(JSIZE(N,J).GT.MAXA) THEN
         PRINT*,'JSIZE=',JSIZE(N,J)
         STOP 'MAXA TOO SMALL'
        END IF
       END DO
      END DO
      JSIZE(0,0)=1
      N=0
      PAIREND=2
      PAIRINDEX(1)=INDPAIR(1,0,0,1)
      PAIRINDEX(2)=INDPAIR(1,0,2,1)
      PAIRTAB(1)=0.
      PAIRTAB(2)=0.
      WRITE(2) OMI,N,PAIREND
      WRITE(2) (PAIRINDEX(I),PAIRTAB(I),I=1,PAIREND)
      SCALE=1E30
      DO N=1,OMI
       READ(1) OMI,N_READ,CFPEND
       IF(CFPEND.GT.MAXCFP) THEN
        PRINT*,'CFPEND=',CFPEND
        STOP 'MAXCFP TOO SMALL'
       END IF
       READ(1) (CFPINDEX(I),CFPTAB(I),I=1,CFPEND)
       DO 1 J=0,2*N
        IF(JSIZE(N,J).LE.0) GOTO 1
        DO AP=1,JSIZE(N,J)
         DO A=1,JSIZE(N,J)
          DO R=0,2,2
           PAPER(AP,A,R,J)=0.
          END DO
         END DO
        END DO
  1    END DO
       DO I1=1,CFPEND
        CALL DECFP(CFPINDEX(I1),AP,J,R,B,JP)
        CFP1=CFPTAB(I1)
        DO A=1,JSIZE(N,J)
         INDEX=INDCFP(A,J,R,B,JP)
         CFP2=CFPTAB(ISECH(CFPINDEX,CFPEND,INDEX,1))
         PAPER(AP,A,R,J)=PAPER(AP,A,R,J)+CFP1*SCALE*CFP2
        END DO
       END DO
       PAIREND=0
       DO R=0,2,2
        DO 2 J=0,2*N
         IF(JSIZE(N,J).LE.0) GOTO 2
         DO AP=1,JSIZE(N,J)
          DO A=1,JSIZE(N,J)
           PAIREND=PAIREND+1
           IF(PAIREND.GT.MAXPAIR) THEN
            PRINT*,'PAIREND=',PAIREND
            STOP 'MAXPAIR TOO SAMLL'
           END IF
           PAIRTAB(PAIREND)=PAPER(AP,A,R,J)/SCALE
           PAIRINDEX(PAIREND)=INDPAIR(AP,J,R,A)
          END DO
         END DO
  2     END DO
       END DO
       WRITE(2) OMI,N,PAIREND
       WRITE(2) (PAIRINDEX(I),PAIRTAB(I),I=1,PAIREND)
      END DO     !N
C                                            multipoles
      N=0
      MULTEND=IA
      DO R=0,IA-1
       MULTINDEX(R+1)=INDPAIR(1,0,R,1)
       MULTTAB(R+1)=0.
      END DO
      WRITE(4) OMI,N,MULTEND
      WRITE(4) (MULTINDEX(I),MULTTAB(I),I=1,MULTEND)
      DO N=1,OMI
       READ(3) OMI,N_READ,PEND
       IF(PEND.GT.MAXP) THEN
        PRINT*,'PEND=',PEND
        STOP 'MAXP TOO SMALL'
       END IF
       READ(3) (PINDEX(I),PTAB(I),I=1,PEND)
       DO 11 J=0,2*N
        IF(JSIZE(N,J).LE.0) GOTO 11
        DO AP=1,JSIZE(N,J)
         DO A=1,JSIZE(N,J)
          DO R=0,IA-1
           PAPER(AP,A,R,J)=0.
          END DO
         END DO
        END DO
  11   END DO
       DO I1=1,PEND
        CALL DEP(PINDEX(I1),AP,J,R,B,JP)
        P1=PTAB(I1)
        DO A=1,JSIZE(N,J)
         INDEX=INDP(A,J,R,B,JP)
         P2=PTAB(ISECH(PINDEX,PEND,INDEX,1))
         PAPER(AP,A,R,J)=PAPER(AP,A,R,J)+P1*SCALE*P2
        END DO
       END DO
       MULTEND=0
       DO R=0,IA-1
        DO 21 J=0,2*N
         IF(JSIZE(N,J).LE.0) GOTO 21
         DO AP=1,JSIZE(N,J)
          DO A=1,JSIZE(N,J)
           MULTEND=MULTEND+1
           IF(MULTEND.GT.MAXMULT) THEN
            PRINT*,'MULTEND=',MULTEND
            STOP 'MAXMULT TOO SAMLL'
           END IF
           MULTTAB(MULTEND)=PAPER(AP,A,R,J)/SCALE
           MULTINDEX(MULTEND)=INDPAIR(AP,J,R,A)
          END DO
         END DO
  21    END DO
       END DO
       WRITE(4) OMI,N,MULTEND
       WRITE(4) (MULTINDEX(I),MULTTAB(I),I=1,MULTEND)
      END DO !N
      GOTO 10
 100  END

      SUBROUTINE DECFP(INDEX,A,J,R,B,L)
C=====================================================================C
C     Gets CFP                                                        C
C=====================================================================C
      IMPLICIT INTEGER (A-Z)
C
      IND=INDEX
      J=IND/20000000
      IND=IND-20000000*J
      A=IND/200000
      IND=IND-200000*A
      R=IND/10000
      IND=IND-10000*R
      L=IND/100
      B=IND-100*L
      RETURN
      END

      SUBROUTINE DEP(INDEX,AP,JP,R,A,J)
C=====================================================================C
C     Form an index                                                   C
C=====================================================================C
      IMPLICIT INTEGER (A-Z)
C
      IND=INDEX
      J=IND/20000000
      IND=IND-20000000*J
      JP=IND/200000
      IND=IND-200000*JP
      R=IND/10000
      IND=IND-10000*R
      A=IND/100
      AP=IND-100*A
      RETURN
      END
