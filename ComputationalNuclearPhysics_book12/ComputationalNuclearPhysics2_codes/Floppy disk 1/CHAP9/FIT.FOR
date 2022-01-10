C     PROGRAM FIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(5,6),UR(50),UTR(50),RK(50)
      REAL AG(99),AK(99),AI(99,5),GRR(99),CH(99)
C
      AN=0.61D0
      R0=1.3D0
      IZ=0
C     >open files<
      OPEN(UNIT=5,FILE='FOLDING.DAT')
      OPEN(UNIT=35,FILE='FIT.DAT')
c     >read input parameters<
   10 READ (5,21,END=70) AP,AD,VP
      WRITE (6,21) AP,AD,VP
   21 FORMAT (3F10.3)
C
      READ (5,20,END=500) AGA,AKA,ZGA,ZKA,ENERLS,LANF,LEND
   20 FORMAT (5F10.5,2I10)
C     >write input parameters on auxiliary file 35 <
      WRITE (35,20) AGA,AKA,ZGA,ZKA,ENERLS,LANF,LEND
      GR=R0*(AGA**(1.D0/3.D0)+AKA**(1.D0/3.D0))
      IZ=IZ+1
      AG(IZ)=AGA
      AK(IZ)=AKA
      GRR(IZ)=GR
      IR=27
      DR=0.5D0
C
C     >read results from folding <
      READ (5,30) (UR(I),I=1,IR),RA
C
C     >read results from folding (exchanged)<
      READ (5,30) (UTR(I),I=1,IR)
   30 FORMAT (5D15.7)
C
      R=RA
      DO 50 I=1,IR
      RK(I)=(R-GR)/AN-DLOG(1.+DEXP((R-GR)/AN))
C
C     >u(r) symmetrized<
      UR(I)=(UR(I)+UTR(I))*0.5D0
   50 R=R+DR
C
C     >transfer of input parameters<
      WRITE (6,40) AP,AD,VP,AGA,AKA,GR,RA
   40 FORMAT ('1AP=',F5.2,3X,'AD=',F5.2,3X,'VP=',F5.1,
     1 6X,'AGA=',F6.2,3X,'AKA=',F6.2,6X,'GROSSR=',F5.2,6X,'R(STA) =',
     2 F5.1/)
C
C     >least sqare fit<
      DO 200 L=1,5
      DO 180 J=1,5
      R=RA
      S=0.D0
      DO 170 I=1,IR
      S=S+RK(I)**2*(R-GR)**(J-1)*(R-GR)**(L-1)/UR(I)**2
  170 R=R+DR
  180 A(L,J)=S
      R=RA
      S=0.D0
      DO 190 I=1,IR
      S=S+RK(I)*(R-GR)**(L-1)/UR(I)
  190 R=R+DR
  200 A(L,6)=S
C
C     >solving for least square with gauss-jordan procedure<
      CALL GAJO (A,5,6,5,6,*150)
C
      R=RA
      UMAX=0.D0
      CHI=0.D0
      DO 220 I=1,IR
      UTR(I)=0.D0
      DO 205 J=1,5
  205 UTR(I)=UTR(I)+A(J,6)*(R-GR)**(J-1)
      UTR(I)=UTR(I)*RK(I)
C     >chi = least sqare<
      CHI=CHI+(UR(I)-UTR(I))**2/UR(I)**2
      UTR(I)=DABS(UTR(I))
      UR(I)=DABS(UR(I))
      IF (UTR(I) .GT. UMAX) UMAX=UTR(I)
      IF (UR(I) .GT. UMAX) UMAX=UR(I)
  220 R=R+DR
      WRITE (6,210) AN,(J,A(J,6),J=1,5),CHI
  210 FORMAT (  ' ANY=',F5.2,5X,
     1             5('A',I1,'=',F10.4,5X),'CHI**2=',D12.4
     2//     ' PLOT U(R) = *   AND   UT(R) = +'/)
C
C     >transfer results to auxiliary file 35<
      WRITE (35,225) AN,(A(J,6),J=1,5),R0
  225 FORMAT (7F10.4)
      CH(IZ)=CHI
      DO 226 J=1,5
  226 AI(IZ,J)=A(J,6)
C
C     >plotting and comparing original-fitted potential<
      CALL PLOT (UR,UTR,UMAX,IR)
C
      UMAX=DLOG(40.D0)
      DO 240 I=1,IR
      UR(I)=DLOG(UR(I))
  240 UTR(I)=DLOG(UTR(I))
      WRITE (6,230)
  230 FORMAT(/' PLOT LN(U(R)) = *   AND   LN(UT(R)) = +',5X,'4 DECADES'
     1/)
      CALL PLOT (UR,UTR,UMAX,IR)
C
      GOTO 10
C
 70   CLOSE(UNIT=5)
      CLOSE(UNIT=35)
      STOP
C
  500 WRITE (6,510) (J,J=1,5),(AG(I),AK(I),(AI(I,J),J=1,5),GRR(I),CH(I),
     1 I=1,IZ)
  510 FORMAT ('1  AGA',6X,'AKA', 9X,5('A',I1,11X),'R',11X,'CHI'//(1X,
     1 F6.2,3X,F6.2,3X,5(F10.4,3X),3X,F5.2,6X,E12.4))
C
C     >error exit<
  150 WRITE (6,160)
  160 FORMAT (//' GAJO FAILED'//)
      GOTO 70
      END
C
C     >gauss-jordan routine<
      SUBROUTINE GAJO (A,N,NP1,IDIM1,IDIM2,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER Z,S
      DIMENSION A(IDIM1,IDIM2),B(100,100),Z(100),S(100),X(100)
C
      DO 1 I=1,N
      Z(I)=0
    1 S(I)=0
      DO 20 J=1,N
      PJ=0.D0
      DO 10 I=1,N
      DO 10 K=1,N
      IF (S(K).NE.0 .OR. Z(I).NE.0) GOTO 10
      IF (DABS(A(I,K)).LE.PJ) GOTO 10
      PJ=DABS(A(I,K))
      L=I
      M=K
   10 CONTINUE
C
      IF (PJ.EQ.0.D0) RETURN 1
      S(M)=L
      Z(L)=M
      NN=N
      IF (M.EQ.N) NN=N-1
C
      A(L,M)=1./A(L,M)
      DO 20 K=1,N
      IF (K.EQ.M) GOTO 20
      A(L,K)=-A(L,K)*A(L,M)
      DO 15 I=1,N
      IF (I.EQ.L) GOTO 15
      A(I,K)=A(I,K)+A(I,M)*A(L,K)
      IF (K.EQ.NN) A(I,M)=A(I,M)*A(L,M)
   15 CONTINUE
   20 CONTINUE
C
      DO 30 K=1,N
      KK=S(K)
      DO 30 I=1,N
   30 B(I,KK)=A(I,K)
      DO 40 I=1,N
      II=Z(I)
      DO 40 K=1,N
   40 A(II,K)=B(I,K)
C
      IF (N.EQ.NP1) RETURN
C
      DO 41 I=1,N
      X(I)=0.D0
      DO 41 J=1,N
   41 X(I)=X(I)+A(I,J)*A(J,NP1)
      DO 42 I=1,N
   42 A(I,NP1)=X(I)
      RETURN
      END
C
      SUBROUTINE PLOT (UR,UTR,UMAX,IR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION UR(IR),UTR(IR)
      CHARACTER*1 IZ(115)
      CHARACTER*1 IIZ(5)
      DATA IIZ /'*','+','0',' ','.'/
C
      UMIN=UR(1)
      DO 10 I=1,IR
      IF (UR(I) .LT. UMIN) UMIN=UR(I)
      IF (UTR(I) .LT. UMIN) UMIN=UTR(I)
   10 CONTINUE
      DO 20 I=1,115
   20 IZ(I)=IIZ(5)
      WRITE (6,30) UMIN,UMAX,IZ
   30 FORMAT (6X,'UTR',9X,'MIN=',D14.7,79X,'MAX=',D14.7/18X,115A1)
      DU=UMAX-UMIN
C
      DO 80 I=1,IR
      DO 40 J=2,114
   40 IZ(J)=IIZ(4)
      IZ(1)=IIZ(5)
      IZ(115)=IIZ(5)
      I1=IFIX(SNGL((UR(I)-UMIN)/DU*114.+0.5D0))+1
      IF (I1 .GT. 115) GOTO 50
      IZ(I1)=IIZ(1)
   50 I2=IFIX(SNGL((UTR(I)-UMIN)/DU*114.+0.5D0))+1
      IF (I2 .GT.115) GOTO 60
      IZ(I2)=IIZ(2)
      IF (I1 .EQ. I2) IZ(I2)=IIZ(3)
   60 WRITE (6,70) UTR(I),IZ
   70 FORMAT (D15.7,3X,115A1)
   80 CONTINUE
C
      DO 90 I=1,115
   90 IZ(I)=IIZ(5)
      WRITE (6,100) IZ
  100 FORMAT (18X,115A1)
      RETURN
      END
