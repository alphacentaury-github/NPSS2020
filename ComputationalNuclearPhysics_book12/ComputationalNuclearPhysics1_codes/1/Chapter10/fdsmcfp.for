      PROGRAM FDSMCFP
C
C========================================== fitting  version  ========C
C                                                                     C
C     FDU0 Package                                                    C
C     ---- -------                                                    C
C                                                                     C
C     Pair Coefficients of Fractional Parentage                       C
C     ---- ------------    ---------- ---------                       C
C                                                                     C
C=====================================================================C
C
      IMPLICIT  REAL*8 (A-H,O-Z)
      PARAMETER (MAXN=21,MAXNS=84,MAXCFP=7762,MAXP=11221) ! parameters based on
      PARAMETER (MAXA=18)                                 ! maxn=21,active='sp6'
      INTEGER   B,BP,R,RP,JSIZE(0:MAXN,0:2*MAXN),OMI,S,C,T
      INTEGER   A,AP,U
      INTEGER   CFPEND,CFPINDEX(MAXCFP)
      INTEGER   CFPPEND,CFPPINDEX(MAXCFP)
      REAL*8    CFPTAB(MAXCFP),CFPPTAB(MAXCFP)
      REAL*8    WTAB(MAXCFP),PTAB(MAXP),PPTAB(MAXP)
      INTEGER   STATEB(MAXNS),STATEL(MAXNS),STATER(MAXNS)
      INTEGER   PEND,PPEND,PINDEX(MAXP),PPINDEX(MAXP)
      INTEGER   OMITAB(10)
      REAL*8    K(0:3,0:3,0:3)
      REAL*8    OVERLAP(MAXNS,MAXNS),E(MAXNS)
      REAL*8    Z(MAXNS,MAXNS)
      REAL*8    APA(MAXA,MAXA)
      LOGICAL   LTRIANGLE
      CHARACTER ACTIVE*3
C
      INDP(AP,JP,R,A,J)=20000000*J+200000*JP+10000*R+100*A+AP
      PHASE(I)=1-2*MOD(ABS(I),2)
      LTRIANGLE(J1,J2,J3)=(ABS(J1-J2).LE.J3).AND.(J1+J2.GE.J3)
C
 203  PRINT*,'SO8 OR SP6?'
      READ(*,'(A)') ACTIVE
      IF(ACTIVE.EQ.'SO8') THEN
       OPEN(2,FILE='SO8CFP.DAT',FORM='UNFORMATTED',STATUS='NEW')
       OPEN(4,FILE='SO8P.DAT',FORM='UNFORMATTED',STATUS='NEW')
       OPEN(16,FILE='SO8JSIZE.TAB',STATUS='NEW')
       OPEN(8,FILE='SO8JSIZE.DAT',FORM='UNFORMATTED',STATUS='NEW')
       IA=4
      ELSE IF(ACTIVE.EQ.'SP6') THEN
       OPEN(2,FILE='SP6CFP.DAT',FORM='UNFORMATTED',STATUS='NEW')
       OPEN(4,FILE='SP6P.DAT',FORM='UNFORMATTED',STATUS='NEW')
       OPEN(16,FILE='SP6JSIZE.TAB',STATUS='NEW')
       OPEN(8,FILE='SP6JSIZE.DAT',FORM='UNFORMATTED',STATUS='NEW')
       IA=3
      ELSE
       PRINT*,'RE-ENTER ACTIVE MODE'
       GOTO 203
      END IF
      DO I1=0,3
       DO I2=0,3
        DO I3=0,3
         K(I1,I2,I3)=(-1)**(IA-1)*DSQRT(IA*(2*I1+1)*(2*I2+1.0D0))*
     &               SIXJ(2*I1+1,2*I2+1,2*I3+1,IA,IA,IA)
        END DO
       END DO
      END DO
      IF(IA.EQ.4) THEN
       OMITAB(1)=2
       OMITAB(2)=6
       OMITAB(3)=10
      ELSE
       OMITAB(1)=6
       OMITAB(2)=15
       OMITAB(3)=21
      END IF
      DO I_OMI=1,3
       OMI=OMITAB(I_OMI)
       DO N=0,MAXN
        DO J=0,2*MAXN
         JSIZE(N,J)=0
        END DO
       END DO
       JSIZE(0,0)=1
       PEND=1
       L=0
       B=1
       LP=0
       BP=1
       T=0
       INDEX=INDP(BP,LP,T,B,L)
       PINDEX(1)=INDEX
       PTAB(1)=0.
C                                                         ! main loop
       DO N=0,OMI-1
        IF(N.GT.MAXN) THEN
         PRINT*,'N=',N
         STOP 'MAXN TOO SMALL'
        END IF
        CFPPEND=0
        DO 10 J=0,2*(N+1)
         PRINT*,'OMI=',OMI,' N=',N+1,' J=',J
         NS=0
         DO R=0,2,2
          DO 1 L=0,2*N
           IF(.NOT.LTRIANGLE(R,L,J)) GOTO 1
           IF(JSIZE(N,L).LE.0) GOTO 1
           DO B=1,JSIZE(N,L)
            NS=NS+1
            IF(NS.GT.MAXNS) THEN
             PRINT*,'NS=',NS
             STOP 'MAXNS TOO SMALL'
            END IF
            STATEB(NS)=B
            STATEL(NS)=L
            STATER(NS)=R
           END DO
 1        END DO
         END DO
         PRINT*,'NS=',NS
         IF(NS.LE.0) GOTO 10
         DO IP=1,NS
          BP=STATEB(IP)
          LP=STATEL(IP)
          RP=STATER(IP)
          DO I=1,IP
           B=STATEB(I)
           L=STATEL(I)
           R=STATER(I)
           IF((BP.EQ.B).AND.(LP.EQ.L).AND.(RP.EQ.R)) THEN
            OV=OMI
           ELSE
            OV=0.0
           END IF
           IF(N.GT.0) THEN
            DO 2 S=0,2*(N-1)
             IF(JSIZE(N-1,S).LE.0) GOTO 2
             IF(.NOT.LTRIANGLE(LP,R,S)) GOTO 2
             IF(.NOT.LTRIANGLE(L,RP,S)) GOTO 2
             F=PHASE(R+RP-L-LP)*DSQRT((2*L+1.0)*(2*LP+1.0D0))
             F=F*SIXJ(2*RP+1,2*S+1,2*L+1,2*R+1,2*J+1,2*LP+1)
             DO C=1,JSIZE(N-1,S)
              INDEX=20000000*LP+200000*BP+10000*R+100*S+C
              CFP1=CFPTAB(ISECH(CFPINDEX,CFPEND,INDEX,1))
              INDEX=20000000*L+200000*B+10000*RP+100*S+C
              CFP2=CFPTAB(ISECH(CFPINDEX,CFPEND,INDEX,1))
              OV=OV+F*CFP1*CFP2
             END DO
  2         END DO
            DO 3 T=0,IA-1
             IF(.NOT.LTRIANGLE(T,RP,R)) GOTO 3
             IF(.NOT.LTRIANGLE(LP,T,L)) GOTO 3
             F=K(RP,R,T)*PHASE(L+T+J)*DSQRT((2*T+1)*(2*LP+1.0D0))
             F=F*SIXJ(2*L+1,2*T+1,2*LP+1,2*RP+1,2*J+1,2*R+1)
             INDEX=INDP(BP,LP,T,B,L)
             P=PTAB(ISECH(PINDEX,PEND,INDEX,1))
             OV=OV-2*F*P
  3         END DO
           END IF
           OVERLAP(I,IP)=OV
           OVERLAP(IP,I)=OV
          END DO
         END DO
C        CALL EIGRS(OVER,NS,1,E,Z,MAXNS,WK,IER)
         CALL DEVCSF(NS,OVERLAP,MAXNS,E,Z,MAXNS)
         DO II=1,NS
          SUM=0.
          DO I=1,NS
           SUM=SUM+Z(I,II)**2
          END DO
          SUM=SQRT(SUM)
          DO I=1,NS
           Z(I,II)=Z(I,II)/SUM
          END DO
         END DO
         PRINT*,'EIGEN VALUES:'
         PRINT '(1X,3F20.15)',(E(I),I=1,NS)
         A=0
         DO 4 I=1,NS
          IF(ABS(E(I)).LT.1E-4) GOTO 4
          A=A+1
          IF(A.GT.MAXA) THEN
           PRINT*,'A=',A
           STOP 'MAXA TOO SMALL'
          END IF
          F=DSQRT(E(I))
          DO II=1,NS
           Z(II,I)=Z(II,I)/F
          END DO
          DO II=1,NS
           L=STATEL(II)
           B=STATEB(II)
           R=STATER(II)
           CFP=0.
           DO III=1,NS
            CFP=CFP+Z(III,I)*OVERLAP(II,III)
           END DO
           CFPPEND=CFPPEND+1
           IF(CFPPEND.GT.MAXCFP) THEN
            PRINT*,'CFPPEND=',CFPPEND
            STOP 'MAXCFP TOO SMALL'
           END IF
           INDEX=20000000*J+200000*A+10000*R+100*L+B
           WTAB(CFPPEND)=Z(II,I)
           CFPPTAB(CFPPEND)=CFP
           CFPPINDEX(CFPPEND)=INDEX
          END DO
  4      END DO              !good states
         JSIZE(N+1,J)=A
  10    END DO               !J
        WRITE(16,'(1X,A,I2,3X,A,I2)') '      OMI=',OMI,'N=',N+1
        WRITE(16,'(1X,A,10(25I3,:,/,3X))') 'J=',(J,J=0,2*(N+1))
        WRITE(16,'(1X,A,10(25I3,:,/,3X))')
     &                        'W=',(JSIZE(N+1,J),J=0,2*(N+1))
        WRITE(16,*)
        WRITE(8) OMI,N+1
        WRITE(8) ( JSIZE(N+1,J), J=0,2*(N+1) )
        PRINT*,'CFPPEND=',CFPPEND
        WRITE(2) OMI,N+1,CFPPEND
        WRITE(2) (CFPPINDEX(I),CFPPTAB(I),I=1,CFPPEND)
        CFPEND=CFPPEND
        DO I=1,CFPEND
         CFPTAB(I)=CFPPTAB(I)
         CFPINDEX(I)=CFPPINDEX(I)
        END DO
        PRINT*,'MULTIPOLES...'
        PPEND=0
        DO 11 J=0,2*(N+1)
         IF(JSIZE(N+1,J).LE.0) GOTO 11
         DO 12 JP=0,2*(N+1)
          IF(JSIZE(N+1,JP).LE.0) GOTO 12
          DO 13 R=0,IA-1
           IF(.NOT.LTRIANGLE(JP,R,J)) GOTO 13
           DO AP=1,JSIZE(N+1,JP)
            DO A=1,JSIZE(N+1,J)
             APA(AP,A)=0.
            END DO
           END DO
           DO 14 L=0,2*N
            IF(JSIZE(N,L).LE.0) GOTO 14
            DO 15 I=0,2,2
             IF(.NOT.LTRIANGLE(L,I,J)) GOTO 15
             DO 16 LP=0,2*N
              IF(JSIZE(N,LP).LE.0) GOTO 16
              IF(.NOT.LTRIANGLE(R,L,LP)) GOTO 16
              IF(.NOT.LTRIANGLE(JP,I,LP)) GOTO 16
              F=SIXJ(2*R+1,2*L+1,2*LP+1,2*I+1,2*JP+1,2*J+1)
              F=F*PHASE(R+I-LP+J)*DSQRT((2*LP+1.0D0)*(2*J+1))
              DO B=1,JSIZE(N,L)
               DO A=1,JSIZE(N+1,J)
                INDEX=20000000*J+200000*A+10000*I+100*L+B
                W=WTAB(ISECH(CFPPINDEX,CFPPEND,INDEX,1))
                DO BP=1,JSIZE(N,LP)
                 DO AP=1,JSIZE(N+1,JP)
                  INDEX=20000000*JP+200000*AP+10000*I+100*LP+BP
                  CFP=CFPPTAB(ISECH(CFPPINDEX,CFPPEND,INDEX,1))
                  INDEX=INDP(BP,LP,R,B,L)
                  P=PTAB(ISECH(PINDEX,PEND,INDEX,1))
                  APA(AP,A)=APA(AP,A)+W*F*CFP*P
                 END DO
                END DO
               END DO
              END DO
16           END DO
             DO 17 U=0,2,2
              IF(.NOT.LTRIANGLE(JP,U,L)) GOTO 17
              IF(.NOT.LTRIANGLE(R,I,U)) GOTO 17
              F=SIXJ(2*L+1,2*I+1,2*J+1,2*R+1,2*JP+1,2*U+1)
              F=F*PHASE(L+JP)*DSQRT((2*J+1.0D0)*(2*U+1))*K(R,I,U)
              DO A=1,JSIZE(N+1,J)
               DO B=1,JSIZE(N,L)
                INDEX=20000000*J+200000*A+10000*I+100*L+B
                W=WTAB(ISECH(CFPPINDEX,CFPPEND,INDEX,1))
                DO AP=1,JSIZE(N+1,JP)
                 INDEX=20000000*JP+200000*AP+10000*U+100*L+B
                 CFP=CFPPTAB(ISECH(CFPPINDEX,CFPPEND,INDEX,1))
                 APA(AP,A)=APA(AP,A)+W*F*CFP
                END DO
               END DO
              END DO
 17          END DO
 15         END DO   !I
 14        END DO    !L
           DO A=1,JSIZE(N+1,J)
            DO AP=1,JSIZE(N+1,JP)
             PPEND=PPEND+1
             IF(PPEND.GT.MAXP) THEN
              PRINT*,'PPEND=',PPEND
              STOP 'MAXP TOO SMALL'
             END IF
             PPTAB(PPEND)=APA(AP,A)
             INDEX=INDP(AP,JP,R,A,J)
             PPINDEX(PPEND)=INDEX
            END DO
           END DO
13        END DO      !R
12       END DO       !JP
11      END DO        !J
        PRINT*,'PPEND=',PPEND
        WRITE(4) OMI,N+1,PPEND
        WRITE(4) (PPINDEX(I),PPTAB(I),I=1,PPEND)
        PEND=PPEND
        DO I=1,PEND
         PTAB(I)=PPTAB(I)
         PINDEX(I)=PPINDEX(I)
        END DO
       END DO         !N
      END DO          !OMI
      END
