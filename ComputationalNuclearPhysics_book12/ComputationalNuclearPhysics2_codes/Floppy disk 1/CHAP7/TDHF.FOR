      PROGRAM TDHF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  Parameters:
C    NPMAX - maximum number of points in spatial grid
C    NMAX  - maximum number of wave functions
C
      PARAMETER (NPMAX=300,NMAX=15,NPMAX2=NPMAX*2)
C
C  Principal global variables
C
C   PSIR,PSII - real and imaginary parts of wave functions
C   PSI2      - modulus squared of wave functions
C   OCC       - occupation numbers
C   X         - physical coordinates of grid points
C   RHO       - total density
C   W         - mean potential
C   WT        - potential term for the total energy evaluation
C   PHIR,PHII - real and imaginary parts of wave functions at half time
C   PHI2      - modulus squared of half time wave functions
C
C   N         - actual number of wave functions propagated
C   NP        - actual number of spatial grid points
C   NOCC      - number of occupied wave functions (static problem)
C   T         - physical time
C   IT        - number of current time step
C   DAMP      - damping factor for static iterations
C
      DIMENSION PSII(0:NPMAX,NMAX),PSIR(0:NPMAX,NMAX),OCC(NMAX),
     1  RHO(0:NPMAX),W(0:NPMAX),WT(0:NPMAX),PSI2(0:NPMAX,NMAX),
     2  PHII(0:NPMAX,NMAX),PHIR(0:NPMAX,NMAX),PHI2(0:NPMAX,NMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/WF/ATOT,EF,OCC,PSIR,PSII,PSI2
      COMMON/MEAN/RHO,W,WT
      COMMON/FORCE/T0,T3,A,V0
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      COMMON/CONST/H2M,PI,TWOPI,HBARC
      CHARACTER*30 FNAME1,FNAME2
C
C  Define two constants used in various places
C
C  H2M is hbar**2/2m, units MeV*fm**2, pi, twopi are self-evident
C  The nucleon mass is assumed to be 938 MeV/c**2
C
      HBARC=197.32D0
      H2M=20.75436162D0
      PI=4.D0*ATAN(1.D0)
      TWOPI=2.D0*PI
C
C  Input phase
C
      READ(*,*) STATIC,NP,ITMAX,IPLOT,DT,DX
      IF(NP.GT.NPMAX) THEN
        WRITE(*,*) ' NP=',NP,' > allowed value:',NPMAX
        STOP
      ENDIF
      IF(STATIC) THEN
        READ(*,*) FNAME1,ATOT,ERR,N,DAMP
        IF(N.GT.NMAX) THEN
           WRITE(*,*) ' N=',N,' > allowed value:',NMAX
           STOP
        ENDIF
      ELSE
        READ(*,*) FNAME1,FNAME2,ECOLL,DIST
      ENDIF
      READ(*,*) T0,T3,A,V0
C
C  Compute a few constant values
C
      CALL INIT
C
C  Static branch of the code
C
      IF(STATIC) THEN
C
C  Initialization phase
C
        CALL INITWF
        CALL ORTHOG
        CALL REOCC(ERRN)
        CALL DENSIT(PSI2)
C
C  Iteration
C
        DO100 IT=1,ITMAX
C
C  Compute density, statistics and do plot
C
          IF(MOD(IT,IPLOT).EQ.1) THEN
            CALL STAT
            CALL PLOT
          ENDIF
C
C  Iterate, orthogonalize, and redetermine occupation
C
          CALL STATIT
          CALL ORTHOG
          CALL DENSIT(PSI2)
          CALL REOCC(ERRN)
C
C  Check for convergence and write output file if finished
C
          IF(MOD(IT,IPLOT).EQ.1) WRITE(*,*) ' Step ',IT,' Error=',ERRN
          IF(ERRN.LT.ERR) THEN
            WRITE(*,*) ' Convergence was reached for ',IT,
     1        ' iterations.'
            WRITE(*,*) ' Writing file ',FNAME1
            OPEN(UNIT=10,FILE=FNAME1,FORM='UNFORMATTED')
            WRITE(10) NOCC,NP
            WRITE(10) (OCC(I),I=1,NOCC)
            WRITE(10) ((PSIR(I,J),I=0,NP),J=1,NOCC)
            CLOSE(UNIT=10)
            CALL STAT
            STOP
          ENDIF
          T=T+DT
100     CONTINUE
C
C  Dynamic part of the code
C
      ELSE
C
C  Determine initial fragment positions from distance and
C   insert boosted wave functions into grid
C
        IC1=NP/2-DIST/(2.D0*DX)
        IC2=NP-IC1
        N=0
        CALL BOOST(FNAME1,IC1,ECOLL)
        CALL BOOST(FNAME2,IC2,ECOLL)
C
C  Start time stepping
C
        DO 200 IT=1,ITMAX
C
C  Compute density, statistics, and do plot
C
          CALL DENSIT(PSI2)
          IF(MOD(IT,IPLOT).EQ.1) THEN
            CALL STAT
            CALL PLOT
          ENDIF
C
C  First propagate the wave functions to half time in order to estimate
C    the density. Half-time values are stored in the PHI arrays.
C
          DTOLD=DT
          DT=0.5*DT
          CALL ALGAM
          CALL PROPWF(PHIR,PHII,PHI2)
C
C  Now compute half step density and propagate to full step
C
          CALL DENSIT(PHI2)
          DT=DTOLD
          CALL ALGAM
          CALL PROPWF(PSIR,PSII,PSI2)
          T=T+DT
200     CONTINUE
      ENDIF
      STOP
      END
C*************************************************************************
C
      SUBROUTINE INIT
C
C  Initialize various quantities
C
C*************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NPMAX2=2*NPMAX)
      DIMENSION RHO(0:NPMAX),W(0:NPMAX),WT(0:NPMAX),QQ(0:25)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/MEAN/RHO,W,WT
      COMMON/FORCE/T0,T3,A,V0
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      COMMON/CONST/H2M,PI,TWOPI,HBARC
      COMMON/POTCOE/QQ
C
C  Compute coordinate grid values and derived grid quantities
C
      XCM=0.5D0*DX*(NP+1)
      XDT=DX**2*HBARC/(DT*H2M)
      XDT4=4.0*XDT
      DO 10 I=0,NP
10    X(I)=DX*I-XCM
C
C  Compute factors QQ to be used in Yukawa potential.
C
      Y=DX/A
      YH=0.5D0*Y
      EPY=EXP(Y)
      EMY=1.D0/EPY
      Q=(2.D0/Y**2+1.D0/Y+0.25D0)*EXP(-YH)-2.D0/Y**2
      B0=(-4.D0/Y**2+1.5D0)*SINH(YH)+2.D0*COSH(YH)/Y
      BP=(2.D0/Y**2+1.D0/Y+0.25D0)*SINH(YH)-(1.D0/Y+0.5D0)*COSH(YH)
      BM=(2.D0/Y**2-1.D0/Y+0.25D0)*SINH(YH)-(1.D0/Y-0.5D0)*COSH(YH)
      QQ(0)=2.*A*(Q+EMY*BM+1.D0-EXP(-YH))
      QQ(1)=A*EMY*(B0+EPY*BP+EMY*BM-(Q+BP)*EPY)
      DO 20 I=2,25
20      QQ(I)=A*EXP(-I*Y)*(B0+EPY*BP+EMY*BM)
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE DENSIT(PHI2)
C
C  Compute the total nuclear density RHO by summing up over N single
C     particle states occupied with probability OCC.
C  From RHO then compute the mean potential W.
C
C  Note that the coefficients for the Yukawa potential have 
C     been precomputed by INIT and stored in QQ.
C
C*************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15,NPMAX2=NPMAX*2)
      DIMENSION PSII(0:NPMAX,NMAX),PSIR(0:NPMAX,NMAX),OCC(NMAX),
     1  RHO(0:NPMAX),W(0:NPMAX),WT(0:NPMAX),PSI2(0:NPMAX,NMAX),
     2  PHI2(0:NPMAX,NMAX),QQ(0:25)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/WF/ATOT,EF,OCC,PSIR,PSII,PSI2
      COMMON/MEAN/RHO,W,WT
      COMMON/FORCE/T0,T3,A,V0
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      COMMON/CONST/H2M,PI,TWOPI,HBARC
      COMMON/POTCOE/QQ
      DO 10 I=0,NP
        RHO(I)=0.D0
        DO 10 J=1,N
          RHO(I)=RHO(I)+OCC(J)*PHI2(I,J)
10    CONTINUE
      DO 30 I=0,NP
        WNEW=0.75D0*T0*RHO(I)+0.1875D0*T3*RHO(I)**2
        WT(I)=0.375D0*T0*RHO(I)+0.0625D0*T3*RHO(I)**2
        WYUK=0.D0
        DO 20 II=MAX(0,I-25),MIN(NP,I+25)
20        WYUK=WYUK+RHO(II)*QQ(ABS(I-II))
        WNEW=WNEW+TWOPI*A*A*V0*WYUK
        WT(I)=WT(I)+0.5D0*TWOPI*A*A*V0*WYUK
C
C  In the static case, dampen the iteration by averaging with the
C    previous potential
C
        IF(STATIC.AND.IT.GT.1) THEN
          W(I)=DAMP*W(I)+(1.0D0-DAMP)*WNEW
        ELSE
          W(I)=WNEW
        ENDIF
30    CONTINUE
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE ALGAM
C
C  Subroutine to calculate alpha and gamma coefficients from the
C     mean field W.
C  AZR and AZI denote the real and imaginary parts of the A0 
C    coefficients.
C
C*************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15,NPMAX2=2*NPMAX)
      DIMENSION W(0:NPMAX),WT(0:NPMAX),GAMR(0:NPMAX),GAMI(0:NPMAX),
     1  AZR(0:NPMAX),RHO(0:NPMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      COMMON/MEAN/RHO,W,WT
      COMMON/GAM/GAMR,GAMI
      COMMON/CONST/H2M,PI,TWOPI,HBARC
C
C     Calculate real and imaginary parts AZR and AZI
C
      DO 10 I=0,NP
         AZR(I)=-2.D0-W(I)*DX**2/H2M
10    CONTINUE
      AZI=2.D0*XDT
C
C     Evaluate Gamma at NP-1 to begin the backward recursion
C
      DENOM=AZR(NP-1)**2+AZI**2
      GAMR(NP-1)=-AZR(NP-1)/DENOM
      GAMI(NP-1)=AZI/DENOM
C
C     Backward recursion
C
      DO 20 I=NP-1,1,-1
         DENOM=(GAMR(I)+AZR(I-1))**2+(GAMI(I)+AZI)**2
         GAMR(I-1)=-(GAMR(I)+AZR(I-1))/DENOM
20       GAMI(I-1)=(GAMI(I)+AZI)/DENOM
      RETURN
      END
C************************************************************************
      SUBROUTINE PROPWF(PHIR,PHII,PHI2)
C
C  Subroutine to propagate wave function to new time by solving
C     the system of equations.   
C
C  Input: coefficients GAMR, GAMI and old wave functions PSIR, PSII.
C
C  Output: new propagated wave functions in the PHI arrays.
C          The time step is given through XDT4.
C
C************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15)
      DIMENSION BETR(0:NPMAX),BETI(0:NPMAX),GAMR(0:NPMAX),GAMI(0:NPMAX),
     1 PSI2(0:NPMAX,NMAX),PSIR(0:NPMAX,NMAX),PSII(0:NPMAX,NMAX),
     2  PHII(0:NPMAX,NMAX),PHIR(0:NPMAX,NMAX),PHI2(0:NPMAX,NMAX),
     3 OCC(NMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      COMMON/WF/ATOT,EF,OCC,PSIR,PSII,PSI2
      COMMON/GAM/GAMR,GAMI
C
C  Start loop over wave functions
C
      DO 30 J=1,N
C
C     Set initial values for backward recursion, then do recursion
C
        BETR(NP-1)=0 
        BETI(NP-1)=0
        DO 10 IX=NP-1,1,-1
           BETR(IX-1)=GAMR(IX)*(BETR(IX)+XDT4*PSII(IX,J))
           BETR(IX-1)=BETR(IX-1)-GAMI(IX)*(BETI(IX)-XDT4*PSIR(IX,J))
           BETI(IX-1)=GAMR(IX)*(BETI(IX)-XDT4*PSIR(IX,J))
           BETI(IX-1)=BETI(IX-1)+GAMI(IX)*(BETR(IX)+XDT4*PSII(IX,J))
10      CONTINUE
C
C     Now do the forward recursion, with initial chi = 0.
        CHIR=0
        CHII=0
        DO 20 IX=1,NP-1
           TEMP=CHIR
           CHIR=GAMR(IX)*CHIR-GAMI(IX)*CHII+BETR(IX-1)
           CHII=GAMR(IX)*CHII+GAMI(IX)*TEMP+BETI(IX-1)
           PHIR(IX,J)=CHIR-PSIR(IX,J)
           PHII(IX,J)=CHII-PSII(IX,J)
           PHI2(IX,J)=PHIR(IX,J)**2+PHII(IX,J)**2
20      CONTINUE
30    CONTINUE
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE SPSTAT
C
C  Compute various quantities characterizing the s. p. states
C
C*************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15,NPMAX2=2*NPMAX)
      DIMENSION PSI2(0:NPMAX,NMAX),PSIR(0:NPMAX,NMAX),
     1  PSII(0:NPMAX,NMAX),RHO(0:NPMAX),W(0:NPMAX),WT(0:NPMAX),
     2  OCC(NMAX),GAMR(0:NPMAX),GAMI(0:NPMAX)
      DIMENSION ANORM(NMAX),ENERGY(NMAX),CENTER(NMAX),SIGMA(NMAX),
     1 ALEFT(NMAX),EKIN(NMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      COMMON/WF/ATOT,EF,OCC,PSIR,PSII,PSI2
      COMMON/GAM/GAMR,GAMI
      COMMON/SP/ANORM,ENERGY,CENTER,SIGMA,ALEFT,EKIN
      COMMON/MEAN/RHO,W,WT
      COMMON/CONST/H2M,PI,TWOPI,HBARC
C
C  First compute characteristics for the individual wave functions,
C    norm, energy, average x-value (CENTER), mean square fluctuation
C    in space (SIGMA= <X**2>-<X>**2 ), and probability for finding
C    the nucleon to the left of x=0.
C
      DO 20 J=1,N
        ANORM(J)=0.D0
        ENERGY(J)=0.D0
        CENTER(J)=0.D0
        SIGMA(J)=0.D0
        ALEFT(J)=0.D0
        EKIN(J)=0.D0
        DO 10 I=1,NP-1
          ANORM(J)=ANORM(J)+PSI2(I,J)
          EK=
     1    -H2M*PSIR(I,J)*(PSIR(I+1,J)+PSIR(I-1,J)-2.D0*PSIR(I,J))/DX**2
     2    -H2M*PSII(I,J)*(PSII(I+1,J)+PSII(I-1,J)-2.D0*PSII(I,J))/DX**2
          ENERGY(J)=ENERGY(J)+W(I)*PSI2(I,J)+EK
          CENTER(J)=CENTER(J)+X(I)*PSI2(I,J)
          SIGMA(J)=SIGMA(J)+X(I)**2*PSI2(I,J)
          IF(X(I).LT.0.D0) ALEFT(J)=ALEFT(J)+PSI2(I,J)
          EKIN(J)=EKIN(J)+EK
10      CONTINUE
        ANORM(J)=ANORM(J)*DX
        ENERGY(J)=ENERGY(J)*DX
        CENTER(J)=CENTER(J)*DX
        SIGMA(J)=SQRT(SIGMA(J)*DX-CENTER(J)**2)
        ALEFT(J)=ALEFT(J)*DX
        EKIN(J)=EKIN(J)*DX
20    CONTINUE
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE STAT
C
C  Compute various quantities characterizing the total system
C
C*************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15)
      DIMENSION PSI2(0:NPMAX,NMAX),PSIR(0:NPMAX,NMAX),
     1  PSII(0:NPMAX,NMAX),RHO(0:NPMAX),W(0:NPMAX),WT(0:NPMAX),
     2  OCC(NMAX),GAMR(0:NPMAX),GAMI(0:NPMAX)
      DIMENSION ANORM(NMAX),ENERGY(NMAX),CENTER(NMAX),SIGMA(NMAX),
     1 ALEFT(NMAX),EKIN(NMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      COMMON/WF/ATOT,EF,OCC,PSIR,PSII,PSI2
      COMMON/MEAN/RHO,W,WT
      COMMON/GAM/GAMR,GAMI
      COMMON/SP/ANORM,ENERGY,CENTER,SIGMA,ALEFT,EKIN
      COMMON/CONST/H2M,PI,TWOPI,HBARC
C
C  First compute single particle statistics
C
      CALL SPSTAT
      WRITE(*,'(A,I5,A,F10.3)') 
     1  ' Statistics for time step no.',IT,' at time ',T,' fm/c'
      WRITE(*,'(1X,A1,5A12,A12)')
     1 '#','OCC','NORM','ENERGY','CENTER','SIGMA','ALEFT'
      DO 10 J=1,N
        WRITE(*,'(I3,6F12.5)')
     1    J,OCC(J),ANORM(J),ENERGY(J),CENTER(J),SIGMA(J),ALEFT(J)
10    CONTINUE
C
C  Now compute total averages
C
      ATOT=0.D0
      ETOT=0.D0
      XAVG=0.D0
      SIGTOT=0.D0
      TLEFT=0.D0
      DO30 J=1,N
        ATOT=ATOT+OCC(J)*ANORM(J)
        ETOT=ETOT+OCC(J)*EKIN(J)+0.5D0*H2M*PI*OCC(J)**2
        XAVG=XAVG+OCC(J)*CENTER(J)
        SIGTOT=SIGTOT+OCC(J)*SIGMA(J)
        TLEFT=TLEFT+OCC(J)*ALEFT(J)
30    CONTINUE
      DO 40 I=1,NP
 40   ETOT=ETOT+RHO(I)*WT(I)*DX
      TRIGHT=ATOT-TLEFT
      WRITE(*,'(1X,A,F15.6,T40,A,F15.6)')
     1 'Total mass:',ATOT,' Total energyv:',ETOT,
     2 'Center of mass:',XAVG,'Fluctuation:',SIGTOT,
     3 'Left hand mass:',TLEFT,'Right hand mass:',TRIGHT
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE BOOST(FNAME,IC,EC)
C
C  Read in a set of wave functions from the file with name given
C    by FNAME and place it with its center at IC, also boost its
C    translational energy to correspond to EC.
C    N is updated to reflect the number of wave functions read.
C
C*************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15,NPMAX2=NPMAX*2)
      DIMENSION PSII(0:NPMAX,NMAX),PSIR(0:NPMAX,NMAX),OCC(NMAX),
     1  RHO(0:NPMAX),W(0:NPMAX),WT(0:NPMAX),PSI2(0:NPMAX,NMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/WF/ATOT,EF,OCC,PSIR,PSII,PSI2
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      COMMON/CONST/H2M,PI,TWOPI,HBARC
      CHARACTER*30 FNAME
C
      OPEN(UNIT=10,FILE=FNAME,FORM='UNFORMATTED')
C
C  Read dimensions of wave functions stored in the file
C
      READ(10) NF,NPF
      IF(N+NF.GT.NMAX) THEN
        WRITE(*,*) ' ***** Error: too many wave functions to be read',
     1  ' from file ',FNAME
        WRITE(*,*) ' Previously read:',N,', number in current file:',
     1   NF,', maximum allowed:',NMAX
        STOP
      ENDIF
C
C  Clear wave function arrays
C
      DO 10 J=1,NF
      DO 10 I=0,NPMAX
        PSIR(I,J+N)=0.D0
        PSII(I,J+N)=0.D0
10      PSI2(I,J+N)=0.D0
C
C  Read real part using PSI2 as intermediate storage, then put it into
C    PSIR at correct position. Also correct normalization.
C
      READ(10) (OCC(J+N),J=1,NF)
      READ(10) ((PSI2(I,J+N),I=0,NPF),J=1,NF)
      DO 30 J=1,NF
        AN=0.D0
        ID=IC-NPF/2
        DO 20 I=MAX(0,-ID),MIN(NPF,NPMAX-ID)
          PSIR(I+ID,J+N)=PSI2(I,J+N)
20        AN=AN+PSI2(I,J+N)**2
        AN=1.D0/SQRT(AN*DX)
        DO 30 I=MAX(0,-ID),MIN(NPF,NPMAX-ID)
30        PSIR(I+ID,J+N)=AN*PSIR(I+ID,J+N)
C
C  Boost wave functions; first compute wave vector. If the center
C   for the wave function is above NP/2, the vector is set negative
C   to get the correct collision behaviour.
C
      AK=SQRT(EC/H2M)
      IF(IC.GT.NP/2) AK=-AK
C
      DO 40 I=0,NP
        COSK=COS(AK*X(I))
        SINK=SIN(AK*X(I))
        DO 40 J=1,NF
          PSI2(I,J+N)=PSIR(I,J+N)**2
          PSII(I,J+N)=PSIR(I,J+N)*SINK
40        PSIR(I,J+N)=PSIR(I,J+N)*COSK
C
C  Finally, update number of wave functions assigned 
C    and close input file.
C
      N=N+NF
      CLOSE(UNIT=10)
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE ORTHOG
C
C  Orthonormalize the wave functions
C  Note that since this is used only in the static calculation,
C    it is assumed that the wave functions are real throughout.
C
C*************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15)
      DIMENSION PSII(0:NPMAX,NMAX),PSIR(0:NPMAX,NMAX),OCC(NMAX),
     1  PSI2(0:NPMAX,NMAX),ANORM(NMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/WF/ATOT,EF,OCC,PSIR,PSII,PSI2
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      DO50 K=1,N
C
C  First compute norm of one vector
C
      ANORM(K)=0.D0
      DO 10 I=0,NP
10    ANORM(K)=ANORM(K)+PSIR(I,K)**2
C
C  Now get overlap with other vectors and orthogonalize them
C
      DO30 J=K+1,N
        OLAP=0.D0
        BNORM=0.D0
        DO 20 I=0,NP
          OLAP=OLAP+PSIR(I,K)*PSIR(I,J)
20        BNORM=BNORM+PSIR(I,J)**2
        FACTOR=OLAP/SQRT(ANORM(K)*BNORM)
        DO 30 I=0,NP
30        PSIR(I,J)=PSIR(I,J)-FACTOR*PSIR(I,K)
C    
C  Finally normalize the orthogonalized wave function and compute
C    PSI2
C    
        FACTOR=1.D0/SQRT(ANORM(K)*DX)
        DO 40 I=0,NP
          PSIR(I,K)=FACTOR*PSIR(I,K)
40        PSI2(I,K)=PSIR(I,K)**2
50    CONTINUE
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE STATIT
C
C  Perform a static iteration step on the wave function
C
C*************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15,NPMAX2=2*NPMAX)
      DIMENSION PSII(0:NPMAX,NMAX),PSIR(0:NPMAX,NMAX),OCC(NMAX),
     1  PSI2(0:NPMAX,NMAX),ANORM(NMAX),
     2  RHO(0:NPMAX),W(0:NPMAX),WT(0:NPMAX)
      DIMENSION PSIO(0:NPMAX,NMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/WF/ATOT,EF,OCC,PSIR,PSII,PSI2
      COMMON/MEAN/RHO,W,WT
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      COMMON/CONST/H2M,PI,TWOPI,HBARC
C
C  Operate on the psi's with the Hamiltonian, replacing the values
C    immediately.
C
      DTH=H2M*DT/(DX**2*HBARC)
      WMIN=+1.E10
      DO 1 I=1,NP-1
1     IF(W(I).LT.WMIN) WMIN=W(I)
      DO10 J=1,N
      PT=0.D0
      PSIO(0,J)=0.D0
      PSIO(NP,J)=0.D0
      DO 10 I=1,NP-1
        PSIO(I,J)=PSIR(I,J)
        PS=(1.D0-DT*(W(I)-WMIN)/HBARC)*PSIR(I,J)
     1     +DTH*(PT+PSIR(I+1,J)-2.D0*PSIR(I,J))
        PT=PSIR(I,J)
10      PSIR(I,J)=PS
C
C  Damp iteration by averaging normalized wf with previous one
C
      DO 30 J=1,N
                  FNORM=0.D0
                  DO 20 I=0,NP
20                        FNORM=FNORM+PSIR(I,J)**2
                  FNORM=1.D0/SQRT(FNORM)
      DO 30 I=1,NP-1
30         PSIR(I,J)=DAMP*PSIO(I,J)+(1.0D0-DAMP)*PSIR(I,J)/FNORM
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE REOCC(ERRN)
C
C  During static calculation, adjust the Fermi energy and the number
C    of occupied states such as to get the correct total nucleon 
C    density
C
C*************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15)
      DIMENSION PSI2(0:NPMAX,NMAX),PSIR(0:NPMAX,NMAX),
     1  PSII(0:NPMAX,NMAX),RHO(0:NPMAX),W(0:NPMAX),WT(0:NPMAX),
     2  OCC(NMAX),GAMR(0:NPMAX),GAMI(0:NPMAX)
      DIMENSION ANORM(NMAX),ENERGY(NMAX),CENTER(NMAX),SIGMA(NMAX),
     1 ALEFT(NMAX),EOLD(NMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/WF/ATOT,EF,OCC,PSIR,PSII,PSI2
      COMMON/SP/ANORM,ENERGY,CENTER,SIGMA,ALEFT
      COMMON/CONST/H2M,PI,TWOPI,HBARC
C
C  Start by computing the characteristics of the s. p. states
C    First save energies to use for convergence check
C
      DO 77 J=1,N
77    EOLD(J)=ENERGY(J)
      CALL SPSTAT
C
C  Now reorder the states to do the occupation correctly
C    Note that the other quantities computed in SPSTAT are not
C    reordered, since they will be recomputed before use
C
      DO 2 K=1,N-1
        DO 2 J=K+1,N
          IF(ENERGY(J).LT.ENERGY(K)) THEN
            SAV=ENERGY(K)
            ENERGY(K)=ENERGY(J)
            ENERGY(J)=SAV
            DO 1 I=0,NP
              SAV=PSIR(I,K)
              PSIR(I,K)=PSIR(I,J)
1             PSIR(I,J)=SAV
          ENDIF
2     CONTINUE
C
C  The number of occupied states is raised until one is found
C    that satisfies the conditions.
C
      DO20 NOCC=1,N-1
        EF=ATOT*H2M*PI
        DO 10 J=1,NOCC
10      EF=EF+ENERGY(J)
        EF=EF/NOCC
        IF(EF.GT.ENERGY(NOCC).AND.EF.LT.ENERGY(NOCC+1)) GOTO 30
20    CONTINUE
      WRITE(*,*) ' Cannot find correct occupation for single ',
     1 'particle spectrum:'
      WRITE(*,'(10F8.3)') (ENERGY(J),J=1,N)
      STOP
C
C  In case of success, compute occupation numbers
30    DO 40 J=1,NOCC
40      OCC(J)=(EF-ENERGY(J))/(H2M*PI)
      DO 50 J=NOCC+1,N
50      OCC(J)=0.D0
C
C  Check convergence
C
      ERRN=0.D0
      DO 60 J=1,NOCC
60    ERRN=ERRN+OCC(J)*ABS(ENERGY(J)-EOLD(J))
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE INITWF
C
C  Initialize wave functions for static calculation by
C    harmonic oscillator states. They are selected such that
C    the scale of the oscillator corresponds to the nuclear radius.
C
C*************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15,NPMAX2=NPMAX*2)
      DIMENSION PSII(0:NPMAX,NMAX),PSIR(0:NPMAX,NMAX),OCC(NMAX),
     1  PSI2(0:NPMAX,NMAX),H(NPMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/WF/ATOT,EF,OCC,PSIR,PSII,PSI2
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
C
      R0=ATOT/0.30D0
      DO 30 I=0,NP
        XI=1.4*X(I)/R0
        EFAC=EXP(-0.5D0*XI*XI)
C  Generate Hermite polynomials by recursion
        H(1)=1.D0
        H(2)=2.D0*XI
        DO 10 J=3,N
10      H(J)=2.D0*XI*H(J-1)-2.D0*(J-2)*H(J-2)
        DO 20 J=1,N
        PSIR(I,J)=H(J)*EFAC
        IF(I.EQ.0.OR.I.EQ.NP) PSIR(I,J)=0.D0
20      PSII(I,J)=0.D0
30    CONTINUE
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE PLOT
C
C     Subroutine for plotting the probability distribution and potential
C       on the printer.
C     This is done by inserting appropriate characters into an array SCREEN
C       that acts like a video screen in allowing direct access to
C       all positions. The array is then dumped for printing.
C
C     Local parameters control the size of the plot in order to
C       be able to do screen as well as printer plots.
C      IPY - number of lines
C      IPXMAX - maximum number of character positions across
C      IDISTX - resulting spacing between x-values
C
C************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPMAX=300,NMAX=15,NPMAX2=NPMAX*2)
      PARAMETER (IPXMAX=70,IPY=25)
      DIMENSION RHO(0:NPMAX),W(0:NPMAX),WT(0:NPMAX)
      LOGICAL STATIC
      COMMON/DIMEN/DAMP,N,NP,NOCC,IT,T,STATIC
      COMMON/MEAN/RHO,W,WT
      COMMON/GRID/DX,DT,XDT,XDT4,X(0:NPMAX),XCM
      CHARACTER*1 SCREEN(120,60),BUF*5
C
C     Clear screen
C
      DO 10 I=1,120
        DO 10 J=1,60
          SCREEN(I,J)=' '
10    CONTINUE
C
C     Determine x axis length
C
      IPX=IPXMAX
      IDISTX=(NP+IPX)/IPX
      IPX=(NP+1)/IDISTX
C
C     Draw x-axis
C
      DO 20 I=1,IPX
20    SCREEN(I,2)='-'
C
C     Plot tick marks and indicator for midpoint
C
      DO 30 I=0,4
        J=MAX(I*(IPX/4),1)
        SCREEN(J,2)='+'
        WRITE(BUF,'(F5.1)') X(J*IDISTX)
        DO 30 K=1,5
          SCREEN(J+K-1,1)=BUF(K:K)
30    CONTINUE
      J=NP/(2*IDISTX)
      SCREEN(J,2)='I'
C
C     First calculate scaling parameters
C
      CALL SCALE(W,WMIN,WMAX,NP)
      CALL SCALE(RHO,RHOMIN,RHOMAX,NP)
      WSCALE=(IPY-2)/(WMAX-WMIN)
      RSCALE=(IPY-2)/(RHOMAX-RHOMIN)
C
C     Plot both potential and RHO into SCREEN.
C
      DO 50 I=0,NP,IDISTX
        J=(W(I)-WMIN)*WSCALE+2
        SCREEN(I/IDISTX,J)='o'
        J=(RHO(I)-RHOMIN)*RSCALE+2
        SCREEN(I/IDISTX,J)='*'
50    CONTINUE
C
C     Print SCREEN 
C
      WRITE(6,'(1H1)')
      DO 60 I=IPY,1,-1
60    WRITE(6,'(120A1)') (SCREEN(K,I),K=1,IPX)
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE SCALE(Y,YMIN,YMAX,N)
C
C  Find maximum and minimum in an array Y
C
C*************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(0:N)
      YMIN=Y(0)
      YMAX=Y(0)
      DO 40 I=1,N
        IF(Y(I).GT.YMAX) YMAX=Y(I)
        IF(Y(I).LT.YMIN) YMIN=Y(I)
40    CONTINUE
      RETURN
      END

