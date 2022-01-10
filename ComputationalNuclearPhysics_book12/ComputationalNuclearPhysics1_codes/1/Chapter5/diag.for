      SUBROUTINE EIGSYM(NM,LP,MV,R,ROOT,EIGV,A,B,W,Q,WW,RM,LIG,NITER)
C O.S.	AUGUST 1988 : precision improved 
C	Note : in eigv space should be reserved upto EIGV(LP+1,NEIG+1)
C	         with NEIG=ABS(MV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION R(LP,2),EIGV(LP,1),ROOT(1),A(2),B(2),W(1),Q(1)
      DIMENSION WW(1),RM(1),LIG(1)
      REAL R,ROOT,EIGV
      IF(NM-2)13,12,1
    1 CALL TRIDI(LP,NM,R,A,B,W,Q,Q)
    2 B(1)=0.
      NEIG=ABS(MV)
      CALL SYMQR(A,B,NM,EIGV(1,NEIG),W,MV,ROOT)
C	get eigenvectors
      IF(NEIG.EQ.0) GOTO 14
      ITER=1
      IF((NITER.GE.0) .AND. (NITER.LT.3)) ITER=NITER
      DO 11 I=1,NEIG
      NUMBER=NM-NEIG+I+1
      CALL VECTOR(A,B,NM,LP,R,EIGV(1,I),Q,W,WW,RM,LIG,
     & NUMBER,EIGV(1,I),EIGV(1,NEIG),ITER)
   11 CONTINUE
      GO TO 14
C	special treatment for 2x2 matrix
   12 A(1)=R(1,1)
      A(2)=R(2,2)
      B(2)=R(1,2)
      GO TO 2
C	special treatment for 1x1 matrix
   13 A(1)=R(1,1)
      B(1)=0.
      ROOT(1)=A(1)
      EIGV(1,1)=1.
   14 RETURN
      END
*DECK SYMQR
      SUBROUTINE SYMQR(A,B,N,W,E,MV,ROOT)
C O.S.	AUGUST 1988 : precision improved 
C	W contains on return the ordered eigenvalues in double precision
C	ROOT contains on return the ordered eigenval. in single precision
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(1),B(1),W(1),E(1),ROOT(1)
      REAL ROOT
      ILIM = 1000
      PRECS=1.D-15
	ONE=1.0D0
      BASE=2.D0
      HOV=BASE**50
      SQRT2=SQRT(2.D0)
      N1=N-1
C 
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     GET EIGENVALUES OF TRIDIAGONAL FORM BY KAHAN-VARAH Q-R METHOD
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C 
      TOL = PRECS/(10.*FLOAT(N))
      BMAX = 0.
      TMAX = 0.
      W(N+1) = 0.
      DO 300  I = 1,N
      E(I)=A(I)
      BMAX = MAX(BMAX,ABS(B(I)))
  300 TMAX = MAX(BMAX,ABS(A(I)),TMAX)
      SCALE = ONE
      IF (BMAX.EQ.0.)  GO TO 460
      DO 310  I = 1,ILIM
      IF (SCALE*TMAX.GT.HOV)  GO TO 320
  310 SCALE = SCALE*BASE
  320 DO 330  I = 1,N
      E(I) = A(I)*SCALE
      W(I) = (B(I)*SCALE)**2
  330 CONTINUE
      DELTA = TMAX*SCALE*TOL
      EPS = DELTA**2
      K = N
C		Start iteration
  350 L = K
      IF (L.LE.0)  GO TO 460
      L1 = L - 1
      DO 360  I = 1,L
      K1 = K
      K = K - 1
      IF (W(K1).LT.EPS)  GO TO 380
  360 CONTINUE
  380 IF (K1.NE.L)  GO TO 400
      W(L) = 0.
      GO TO 350
  400 T = E(L) - E(L1)
      X = W(L)
      S = SQRT(X)
      Y = 0.5D0*T
      IF (ABS(T).GT.DELTA)  S = (X/Y)/(ONE+SQRT(ONE+X/Y**2))
      E1 = E(L) + S
      E2 = E(L1)- S
      IF (K1.NE.L1)  GO TO 430
      E(L)  = E1
      E(L1) = E2
      W(L1) = 0.
      GO TO 350
  430 RAMBDA = E1
      IF (ABS(T).LT.DELTA.AND.ABS(E2).LT.ABS(E1))  RAMBDA = E2
      S = 0.
      C = ONE
      GG = E(K1)-RAMBDA
      GO TO 450
  440 C = F/T
      S = X/T
      X = GG
      GG = C*(E(K1)-RAMBDA) - S*X
      E(K) = (X-GG) + E(K1)
  450 IF (ABS(GG).LT.DELTA) GG = GG + SIGN(C*DELTA,GG)
      F = GG**2/C
      K = K1
      K1 = K + 1
      X = W(K1)
      T = X + F
      W(K) = S*T
      IF (K.LT.L)  GO TO 440
      E(K) = GG + RAMBDA
      GO TO 350
C		rescale energies
  460 DO 470  I = 1,N
  470 E(I) = E(I)/SCALE
C	order eigenvalues, array E
      DO 510 I=2,N
      N1=I-1
      DO 502 K=1,N1
      IF(E(I)-E(K))502,502,504
  502 CONTINUE
      GO TO 510
  504 X=E(I)
      MM=K+N1
      DO 506 J=K,N1
      MN=MM-J
  506 E(MN+1)=E(MN)
      E(K)=X
  510 CONTINUE
      NM=N
      J=NM
      K=1
      DO 5 I=1,NM
      IF(ABS(E(K))-ABS(E(J)))3,3,4
    3 W(I)=E(J)
      J=J-1
      GO TO 5
    4 W(I)=E(K)
      K=K+1
    5 CONTINUE
C	reverse order eigenvalues if mv<0
      IF(ISIGN(1,MV))6,9,9
    6 DO 7 I=1,NM
    7 E(I)=W(I)
      DO 8 I=1,NM
    8 W(I)=E(NM-I+1)
    9 CONTINUE
      DO 90 I=1,NM
      ROOT(I)=W(I)
   90 CONTINUE
C	Shift energies to end of space for eigenvectors, remember
C	 than on call W==EIGV(1,NEIG) !!!!
      NEIG=ABS(MV)
      DO 91 I=NEIG,1,-1
      W(NM-NEIG+I+1)=W(I)
   91 CONTINUE
      RETURN
      END
      SUBROUTINE TRIDI(LP,NM,R,A,B,W,Q,P)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION R(LP,1),A(1),B(1),W(1),Q(1),P(1)
      REAL R
	ONE=1.0D0
      W(NM+1)=0.  
      Q(NM+1)=0.
      NN=NM-2
      DO 1 I=1,NM
    1 A(I)=R(I,I)
      DO 16 IR=1,NN
      W(IR)=0.
      S=0.
      IRR=IR+1
      DO 3 I=IRR,NM
    3 S=S+R(I,IR)*R(I,IR)
      S=SQRT(S)
      W(IR+1)=R(IR+1,IR)+S*SIGN(1.,R(IR+1,IR))
      IRR=IR+2
      DO 4 I=IRR,NM
    4 W(I)=R(I,IR)
      TOKR=S*S+S*ABS(R(IR+1,IR))
      IF(TOKR)5,14,5
    5 DO 10 I=IR,NM
C     *********
      D=0.
      DO 7 J=IR,I
    7 D=D+R(I,J)*W(J)
C      *********
      IF(I-NM)8,10,10
    8 JK=I+1
C     *********
      DO 9 J=JK,NM
    9 D=D+R(J,I)*W(J)
C     *********
   10 P(I)=D/TOKR
      FAK=0.
      DO 11 I=IR,NM
   11 FAK=FAK+W(I)*P(I)
      FAK=FAK*.5D0/TOKR
      DO 12 I=IR,NM
   12 Q(I)=P(I)-FAK*W(I)
      IRR=IR+1
      DO 13 I=IRR,NM
      DO 13 J=I,NM
   13 R(J,I)=R(J,I)-W(I)*Q(J)-W(J)*Q(I)
      GO TO 15
   14 B(IR+1)=0.
   15 X=A(IR)
      A(IR)=R(IR,IR)
      R(IR,IR)=X
      B(IR+1)=-S*SIGN(1.,R(IR+1,IR))
   16 R(IR+1,IR)=W(IR+1)
      X=A(NM-1)
      A(NM-1)=R(NM-1,NM-1)
      R(NM-1,NM-1)=X
      X=A(NM)
      A(NM)=R(NM,NM)
      R(NM,NM)=X
      B(NM)=R(NM,NM-1)
      RETURN
      END
      SUBROUTINE VECTOR(A,B,NM,LP,R,VEC,U,V,W,RM,LIG,
C O.S.	AUGUST 1988 : precision improved 
     & NUMBER,VCC,ENERGS,ITER)
C	NOTE: 	in the neigth-1 and neigth call the array ENERGS will
C	   	  be overwritten.
C	solves : 0=B(I)*VCC(I-1)+(A(I)-E)*VCC(I)+B(I+1)*VCC(I+1)
C		    with  E=ENERGS(NUMBER)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION R(LP,1),A(1),B(2),U(1),V(1),VEC(1),W(1),RM(1),LIG(1)
      DIMENSION VCC(1),ENERGS(1)
      REAL R,VEC
      DATA SMALL /1.0D-15/
      E=ENERGS(NUMBER)
      NN=NM-1
      P=A(1)-E
      Q=B(2)
      B(NM+1)=0.
      DO 1 I=1,NM
    1 LIG(I)=1
      DO 6 I=1,NN
      IF(ABS(B(I+1))-ABS(P))2,3,5
    2 U(I)=P
      V(I)=Q
      W(I)=0.
      RM(I+1)=B(I+1)/P
      P=A(I+1)-E-RM(I+1)*Q
      Q=B(I+2)
      GO TO 6
    3 IF(ABS(P)-SMALL)4,5,5
    4 P=SMALL
C    REVISED(700826) TO PROVIDE THE RIGHT BRANCH
      GO TO 2
    5 U(I)=B(I+1)
      V(I)=A(I+1)-E
      W(I)=B(I+2)
      RM(I+1)=P/U(I)
      P=Q-RM(I+1)*V(I)
      Q=-RM(I+1)*W(I)
      LIG(I)=-1
    6 IF(ABS(P).LT.SMALL) P=SMALL
C
      W(NM-1)=0.
      W(NM)=0.
      V(NM)=0.
      U(NM)=P
      IF(ABS(U(NM))-SMALL)7,8,8
    7 U(NM)=SMALL
    8 NEM=NM-1
      VCC(NM)=1.D0/U(NM)
      VCC(NM+1)=0.
C      VCC(NM+2)=0.
      X=VCC(NM)*VCC(NM)
      DO 9 J=2,NM
      I=NM-J+1
      VCC(I)=(1.D0-VCC(I+1)*V(I)-VCC(I+2)*W(I))/U(I)
    9 X=VCC(I)*VCC(I)+X
C	renormalize
      X=1.D0/SQRT(X)
      DO 10 I=1,NM
   10 VCC(I)=VCC(I)*X
C
      NJ=-ITER
C  nj<0 to obtain better precision
C	The number of iteration steps is 1-nj with nj<=0
   11 NJ=NJ+1
      DO 13 I=2,NM
      IF(LIG(I-1))12,12,13
   12 VE=VCC(I-1)
      VCC(I-1)=VCC(I)
      VCC(I)=VE
   13 VCC(I)=VCC(I)-RM(I)*VCC(I-1)
C	reiterate eigenvector for improved precision 
      VCC(NM)=VCC(NM)/U(NM)
      X=VCC(NM)*VCC(NM)
      DO 14 J=2,NM
      I=NM-J+1
      VCC(I)=(VCC(I)-VCC(I+1)*V(I)-VCC(I+2)*W(I))/U(I)
   14 X=VCC(I)*VCC(I)+X
C	renormalize
      X=1.D0/SQRT(X)
      DO 15 I=1,NM
   15 VCC(I)=VCC(I)*X
      IF(NJ-1)11,16,16
C
   16 NN=NM-2
C	When two eigenvalues are degenerate VCC(I0) and VCC(I0+1) are
C	  equal to zero and there is an arbitrary sign for VCC(I) , I>I0
CTEST      Y=0
CTEST      DO I=2,NM
CTEST      X=B(I)*VCC(I-1)+(A(I)-E)*VCC(I)
CTEST      U(I)=X
CTEST      X=ABS(X+B(I+1)*VCC(I+1))
CTEST      IF(X.GT.Y) Y=X
CTEST      ENDDO
CTEST      WRITE(2,*) 'X',(U(I),I=2,NM)
CTEST      WRITE(2,*) 'Y',Y
      DO 20 IJ=1,NN
      I=NEM-IJ
      NT=I+1
      FAK=0.
C     *********
      DO 17 J=NT,NM
   17 FAK=FAK+VCC(J)*R(J,I)
C     *********
      Z=B(NT)*R(NT,I)
      IF(Z)18,20,18
   18 FAK=-FAK/Z
      DO 19 J=NT,NM
   19 VCC(J)=VCC(J)-FAK*R(J,I)
   20 CONTINUE
      X=0.
C	renormalize
      DO 21 I=1,NM
   21 X=X+VCC(I)*VCC(I)
      X=SQRT(X)
      DO 22 I=1,NM
   22 VEC(I)=VCC(I)/X
      RETURN
      END
