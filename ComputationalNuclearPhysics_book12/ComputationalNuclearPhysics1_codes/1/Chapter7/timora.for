C ****************************************************************
C
C     RELATIVISTIC MEAN FIELD CALCULATIONS FOR FINITE NUCLEI
C     PROGRAM TIMORA.FOR
C     LAST MODIFIED 8/17/89
C
C      AUTHOR:
C
C      CHARLES J. HOROWITZ
C      DEPT OF PHYSICS
C      INDIANA UNIVERSITY
C      BLOOMINGTON, IN 47405
C
C      PHONE   (812) 855-2959
C      BITNET:   CHARLIE@IUCF
C
C      REFERENCE: NUCL. PHYS. A368 (1981) 503
C
C ****************************************************************
C
C      INPUT SHOULD BE PLACED IN FILE TIMORA.INP, AFTER EXTRACTING
C      FROM TIMORA.DAT.                                                 
C
C      OUTPUT INCLUDES:
C      1)SINGLE PARTICLE SPECTRA
C      2)SELF-CONSISTENT MESON FIELDS (OPTICAL POTENTIALS)
C      3)DIRAC WAVEFUNCTIONS (UPPER AND LOWER COMPONENTS)
C      4)NUCLEAR PROTON AND NEUTRON SCALAR AND BARYON DENSITIES
C
C ***************************************************************
C
C VARIABLES:
C     HH=STEP SIZE IN fm
C     XMM=NUCLEON MASS IN MeV
C     IMAX=MAXIMUM NUMBER OF GRID POINTS (<=600)
C
C FIELDS(J,I)
C    I=1 TO 4
C       1 : SCALAR MESON FIELD IN MeV
C       2 : VECTOR
C       3 : RHO
C       4 : PHOTON (COULOMB)
C    J=1 TO IMAX
C       J TH GRID POINT
C       1 : HH fm
C       IMAX: IMAX*HH fm
C
C XRHO(J,I)
C   I = 1 TO 4  CORRESPONDING SOURCE DENSITIES TIMES R**2
C       1 : SCALAR DENSITY * X**2 (IN fm**-1)
C       2 : BARYON
C       3 : ISOVECTOR BARYON
C       4 : PROTON BARYON
C   X= RADIUS = J*HH fm
C
C G(J) = LARGE COMPONENT OF WAVEFUNCTION
C F(J) = SMALL
C
C UNITS ARE IN MeV AND fm AND HBAR*C = 197.33 MeV-fm
C
C **********************************************************************
C
C
      PROGRAM TIMORA
C
      COMMON /GRID/XMM,HH,IMAX
      COMMON /SHELLS/ STATE(7,40),LEVEL,EIGEN(40)
      COMMON /MESONS/ FIELDS(600,4),G2(4),V(600)
      COMMON /DENSIT/ XRHO(600,4),XRHOSP(600)
      COMMON /GREEN/ GIN(600,4),GOUT(600,4),XMASS(3)
      DIMENSION EIGENC(40)
      CHARACTER*1 XLABP,XLABN,XLABEL
      DATA XLABP,XLABN/'p','n'/
C
C
      OPEN(UNIT=23,FILE='TIMORA.LOG',STATUS='UNKNOWN')
      OPEN(UNIT=9,FILE='TIMORA.INP',STATUS='OLD')
      OPEN(UNIT=11,FILE='TIMORA.DEN',STATUS='UNKNOWN')
C
      WRITE(11,1)
1     FORMAT('  r(fm)     rho_p       rho_n      rhos_p',
     1  '      rhos_n (1/fm^3)')
C
      WRITE(6,'(''       THIS IS PROGRAM TIMORA (1989)   '')')
      WRITE(23,'(''       THIS IS PROGRAM TIMORA (1989)   '')')
      READ(9,*) XMM,HH,IMAX
      READ(9,*)XMASS,G2,S0,V0,XR,A,CONVRG,DCONVR,B0,A0,ISTATE
      DO 7 I=1,ISTATE
7         READ(9,11) STATE(1,I),STATE(2,I),STATE(3,I),STATE(4,I),
     1    STATE(5,I),STATE(6,I),STATE(7,I)
      DO 8 I=ISTATE+1,40
8     STATE(1,I)=-2.
11    FORMAT(1X,4F8.3,2A3,F10.3)
      WRITE(23,12) XMM,HH,IMAX,XMASS,G2,S0,V0,XR,A,CONVRG,DCONVR,B0,A0
      WRITE(6,12) XMM,HH,IMAX,XMASS,G2,S0,V0,XR,A,CONVRG,DCONVR,B0,A0
C
12    FORMAT(1X,/,1X,'Nucleon mass (MeV), dr (fm), Grid points',/,1X,
     1  F6.1,F8.4,I5,/,1X,
     2  'Meson masses (MeV): ms, mv, mrho'
     3  ,/,1X,3F12.4,/' Coupling constants: gs^2, gv^2,',
     4  ' grho^2, alpha',/,
     5  1X,4F12.6,/,' Initial potential params: S(0), V(0), XR, A',
     6  /,1X,4F12.6,/,' Convergence tolerances: Hartree, Dirac;',
     7  '  Field strengths: B(0), A(0)',/,
     8  1X, 4F12.4,//, ' Level      degeneracy     Kappa    eigenvalue
     9 isospin  match radius ')
C
      DO 20 I=1,40
          IF (STATE(1,I).LT.0.)GOTO 20
          WRITE(6,22)STATE(5,I),STATE(6,I),STATE(1,I),
     1    STATE(2,I),STATE(3,I),STATE(4,I),STATE(7,I)
          WRITE(23,22)STATE(5,I),STATE(6,I),STATE(1,I),
     1    STATE(2,I),STATE(3,I),STATE(4,I),STATE(7,I)
20    CONTINUE
22    FORMAT(1X,2A3,5F12.2)
C
C initial guess for fields
C
      DO 90 I=1,IMAX
          X=HH*FLOAT(I)
          FIELDS(I,1)=S0/(1.+EXP((X-XR)/A))
          FIELDS(I,2)=V0/(1.+EXP((X-XR)/A))
          FIELDS(I,3)=B0/(1.+EXP((X-XR)/A))
90    FIELDS(I,4)=A0/(1.+EXP((X-XR)/A))
C
C calculate meson greens functions
C
      CALL MESINT
      ITERAT=1
      IFLAG=0
      GOTO 115
C
C come here for next iteration
C
100   ITERAT=ITERAT+1
C
C equate new to old eigenvalue
C
      DO 110 I=1,40
          EIGENC(I)=0.                                                  
          STATE(3,I)=EIGEN(I)
110   CONTINUE
C
C zero densities
C
115   DO 120 I=1,IMAX
          XRHOSP(I)=0.                                                  
          DO 120 L=1,4
120   XRHO(I,L)=0.
      WRITE(6,1130)ITERAT
1130  FORMAT(/' NOW WORKING ON ITERATION ',I3)
      WRITE(23,130)ITERAT
130   FORMAT(///' ITERATION ',I3,' CALCULATIONS *****************')
      IF (IFLAG.NE.1) WRITE(23,132)
      IF (IFLAG.NE.1) GOTO 1333
      READ(9,*)ISTATE
      DO 131 I=1,ISTATE
131   READ(9,11) STATE(1,I),STATE(2,I),STATE(3,I),STATE(4,I),
     1 STATE(5,I),STATE(6,I),STATE(7,I)
      DO 8888 I=ISTATE+1,40
8888  STATE(1,I)=-2.
C
      WRITE(6,134)
      IF (ISTATE.NE.0.) WRITE(23,1134)
      DO 133 I=1,40
          IF(STATE(1,I).LT.0.) GOTO 133
          WRITE(23,22)STATE(5,I),STATE(6,I),STATE(1,I),
     1    STATE(2,I),STATE(3,I),STATE(4,I),STATE(7,I)
133   CONTINUE
C
134   FORMAT(///' CONVERGED *********************************')
1134  FORMAT(//' Wave functions for the following states: '/
     1  ' Level       degeneracy   Kappa     eigenvalue   isospin',
     2  '  match radius')
132   FORMAT(/' Level  isospin  E(MeV)    Delta E'/1X,34('-'))
C
1333  DO 140 LEVEL=1,40
C
          IF(STATE(1,LEVEL).LT.0.)GOTO 140
          CALL DIRAC(IFLAG,DCONVR)
          EIGENC(LEVEL)=EIGEN(LEVEL)-STATE(3,LEVEL)
          XLABEL=XLABP
          IF (STATE(4,LEVEL).LT.0.)XLABEL=XLABN
          IF(IFLAG.NE.1)WRITE(23,135)STATE(5,LEVEL),STATE(6,LEVEL),
     1    XLABEL,EIGEN(LEVEL),EIGENC(LEVEL)
140   CONTINUE
135   FORMAT(1X,A3,A3,2X,A3,F11.5,F11.6)
C
      IF (IFLAG.EQ.1) STOP
C
C test for convergence
C
      IFLAG=1
      DO 150 I=1,40
150   IF(ABS(EIGENC(I)).GT.CONVRG)IFLAG=0
C
C solve meson field equations
C
      CALL MESON
C
      IF(IFLAG.EQ.1)WRITE(23,134)
159   WRITE(23,160)
160   FORMAT(//,' r         RhoB        RhoS    ',
     1      '    Phi       V_0       B_0       A_0',/,
     2         ' (fm)     (1/fm^3)    (1/fm^3) ',
     3      '   (MeV)     (MeV)     (MeV)     (MeV)',/,1X,70('-'))
      J=25
      IF(IFLAG.EQ.1)J=5
      DO 170 I=J,IMAX,J
          X=FLOAT(I)*HH
          XBS=XRHO(I,1)/X/X
          XB=XRHO(I,2)/X/X
170   WRITE(23,180) X,XB,XBS,FIELDS(I,1),FIELDS(I,2),
     1 FIELDS(I,3),FIELDS(I,4)
180   FORMAT(1X,F5.2,2E12.4,4F10.3)
C
C calculate rms radius and sum up protons, neutrons, and energies
C
      E=0.
      XPRO=0.
      XNU=0.
      DO 200 I=1,40
          IF (STATE(1,I).LT.0.)GOTO 200
          IF(STATE(4,I).GT.0.)XPRO=XPRO+STATE(1,I)
          IF(STATE(4,I).LT.0.)XNU=XNU+STATE(1,I)
          E=E+STATE(1,I)*EIGEN(I)
200   CONTINUE
C
      XPRMS=0.
      EA=0.
      XNRMS=0.
      EPHI=0.
      EV=0.
      EB=0.
      DO 220 I=1,IMAX
          X=FLOAT(I)*HH
C
C  output densities to unit 11 if converged
C
          X2=X*X
          XP=XRHO(I,4)/X2
          XPS=XRHOSP(I)/X2
          XN=XRHO(I,2)/X2-XP
          XNS=XRHO(I,1)/X2-XPS
          IF (IFLAG.EQ.1) WRITE(11,250) x,xp,xn,xps,xns
          EB=EB+FIELDS(I,3)*XRHO(I,3)
          EV=EV+FIELDS(I,2)*XRHO(I,2)
          EPHI=EPHI+FIELDS(I,1)*XRHO(I,1)
          EA=EA+FIELDS(I,4)*XRHO(I,4)
          XPRMS=XPRMS+X*X*XRHO(I,4)
220   XNRMS=XNRMS+X*X*(XRHO(I,2)-XRHO(I,4))
250   FORMAT(1X,F5.2,4E12.5)
      FACTOR=2.*3.1415926*HH
      EB=-FACTOR*EB
      EV=-FACTOR*EV
      EPHI=EPHI*FACTOR
      E=(E+EPHI+EV+EB-EA*FACTOR)/(XPRO+XNU)
      XPRMS=SQRT(XPRMS/XPRO*4.*3.1415926*HH)
      XNRMS=SQRT(XNRMS/XNU*4.*3.1415926*HH)
      WRITE(23,230)XPRMS,XNRMS,E
      write(6,230) xprms,xnrms,e
230   FORMAT(/' Proton rms radius =',F8.4,' fm   ',
     1        ' Neutron rms radius =',F8.4,' fm  '/
     2  ' Energy/nucleon =',F10.3,' MeV')
      GOTO 100
      END
C ******************************************************************
C
      SUBROUTINE DIRAC(IFLAG,DCONVR)
C
C SOLVES DIRAC EQUATION INWARD FROM RMAX TO RMATCH AND THEN
C OUTWARD FROM ZERO TO RMATCH.  
C EIGENVALUE ADJUSTED WITH DELTA E = -G(MP)*(F^+ - F^-)
C (G AT MATCH POINT TIMES DISCONTINUITY IN LOWER COMPONENT).            
C
C
      COMMON /GRID/XMM,HH,IMAX
      COMMON /DENSIT/ XRHO(600,4),XRHOSP(600)
      COMMON /SHELLS/ STATE(7,40),LEVEL,EIGEN(40)
      COMMON /MESONS/ FIELDS(600,4),G2(4),V(600)
      DIMENSION G(600),F(600)
C
      IF (STATE(3,LEVEL).GT.0.) STATE(3,LEVEL)=EIGEN(LEVEL-1)
      EIGEN(LEVEL)=STATE(3,LEVEL)
      ITURN=1
C combine vector and rho fields plus coulomb
      DO 50 I=1,IMAX
50    V(I)=FIELDS(I,2)+STATE(4,LEVEL)*(FIELDS(I,3)+
     1 FIELDS(I,4))+FIELDS(I,4)/2.
C
C small r solutions
C
5     G(1)=10.*(HH)**(-STATE(2,LEVEL))
      F(1)=(V(1)-FIELDS(1,1)-EIGEN(LEVEL))*G(1)*HH/197.33
      F(1)=F(1)/(1.-2.*STATE(2,LEVEL))
      IF (STATE(2,LEVEL).LT.0.)GOTO 10
      G(1)=10.*(HH)**(1.+STATE(2,LEVEL))
      F(1)=197.33/HH*G(1)*(2.*STATE(2,LEVEL)+1.)
      F(1)=F(1)/(EIGEN(LEVEL)-V(1)-FIELDS(1,1)+2.*XMM)
C
10    YG=G(1)
      JMATCH=IFIX(STATE(7,LEVEL)/HH+.5)
      YF=F(1)
      X=HH
C
      DO 100 I=2,JMATCH
70        CALL RUNGKU(YG,YF,X,HH)
          G(I)=YG
100   F(I)=YF
C
C store endpoint values for matching
C
      YFS=YF
      YGS=YG
C
C At large r, we do not use coulomb wavefunctions
C but rather asymptotic expansions in 1/r. These appear
C to be just fine for bound state coulomb wavefunctions.  
C
      ALPHA=SQRT(-EIGEN(LEVEL)*(EIGEN(LEVEL)+2.*XMM))/197.33
      RMAX=IMAX*HH
      G(IMAX)=EXP(-RMAX*ALPHA)
      X0=-SQRT(-EIGEN(LEVEL)/(EIGEN(LEVEL)+2.*XMM))
      XK=STATE(2,LEVEL)
      X=RMAX
      V0=12.*V(IMAX)/197.33
      E=197.33/2./(EIGEN(LEVEL)+2.*XMM)
      X1=E*(2.*XK+V0*(X0+1./X0))
      X2=E*(2.*V0+(2.*XK+1.)*X1/X0)-X1*X1/2./X0
      X3=E*((2.*XK+2.)*X2/X0+V0*(2.*X2+X1*X1/X0))-X1*X2/X0
      X4=E*((2.*XK+3.)*X3/X0+V0*2.*(X3+X1*X2/X0))-
     1 (2.*X1*X3+X2*X2)/X0/2.
      F(IMAX)=(X0+(X1+(X2+(X3+X4/X)/X)/X)/X)*G(IMAX)
      YG=G(IMAX)
      YF=F(IMAX)
      JTOP=IMAX-JMATCH
      DO 200 J=1,JTOP
          CALL RUNGKU(YG,YF,X,-HH)
          I=IMAX-J
          G(I)=YG
200   F(I)=YF
C
C match solutions
C
      SCALE=YG / YGS
C
      JTOP=JMATCH-1
      DO 250 I=1,JTOP
          F(I)=F(I)*SCALE
250   G(I)=G(I)*SCALE
C
C normalization integral
C
      XNORM=0.
      DO 300 I=1,IMAX
300   XNORM=XNORM+F(I)**2+G(I)**2
      XNORM=XNORM*HH
C
C adjust eigenvalue
C
      DELTAE=-G(JMATCH)*(F(JMATCH)-YFS*SCALE)*197.33/XNORM
      IF (ITURN.EQ.1) DELTAE=DELTAE/2.
      EIGEN(LEVEL)=EIGEN(LEVEL)+DELTAE
      IF(EIGEN(LEVEL).GT.0.)EIGEN(LEVEL)=-4./FLOAT(ITURN)
      ITURN=ITURN+1
      IF(ITURN.GT.50)GOTO 350
      IF (ABS(DELTAE).GT.DCONVR)GOTO 5
350   IF(ITURN.GT.50)WRITE(23,360)
360   FORMAT(' NO DIRAC CONVERGENCE AFTER 50 TRIES,
     1   LEVEL MAY BE UNBOUND')
C
C sum up the densities
C RHO = (2J+1)/NORM(F**2+G**2)/(4PI*X**2)
C
      FACTOR=STATE(1,LEVEL)/XNORM/4./3.1415926
      DO 400 I=1,IMAX
          XRHOSP(I)=XRHOSP(I)+FACTOR*(STATE(4,LEVEL)+.5)*
     1       (G(I)**2-F(I)**2)
          XRHO(I,1)=XRHO(I,1)+FACTOR*(G(I)**2-F(I)**2)
          XRHO(I,2)=XRHO(I,2)+FACTOR*(G(I)**2+F(I)**2)
          XRHO(I,3)=XRHO(I,3)+FACTOR*STATE(4,LEVEL)*
     1       (G(I)**2+F(I)**2)
400   XRHO(I,4)=XRHO(I,4)+FACTOR*(STATE(4,LEVEL)+.5)*
     1  (G(I)**2+F(I)**2)
C
      IF(IFLAG.NE.1)RETURN
C
C Now that we have converged, print out wavefunctions.  
C
      WRITE(23,410) STATE(5,LEVEL),STATE(6,LEVEL),EIGEN(LEVEL)
     1       ,STATE(4,LEVEL)
410   FORMAT(//' Dirac wave fn for level ',2A3,
     1 ':  energy =',1F10.4,
     2 ' MeV,  isospin=',F3.1/1X,
     3    3('  r(fm)    G       F    ')/1X,72('-'))
      XNORM=SQRT(XNORM)
      IMAX3=IMAX/3
      DO 450 I=5,IMAX3,5
          X=FLOAT(I)*HH
          X1=X+IMAX3*HH
          X2=X1+IMAX3*HH
          YF=F(I)/XNORM
          YF1=F(I+IMAX3)/XNORM
          YF2=F(I+IMAX3)/XNORM
          YG=G(I)/XNORM
          YG1=G(I+IMAX3)/XNORM
          YG2=G(I+IMAX3)/XNORM
450   WRITE(23,460)X,YG,YF,X1,YG1,YF1,X2,YG2,YF2
460   FORMAT(1X,3(1X,F5.1,2F8.4,2X))
      RETURN
      END
C ********************************************************************
C
      SUBROUTINE RUNGKU(G,F,X,H)
      G1=H*GUNT(X,G,F,FUNT)
      F1=H*FUNT
C
      G2=H*GUNT(X+H/2.,G+G1/2.,F+F1/2.,FUNT)
      F2=H*FUNT
C
      G3=H*GUNT(X+H/2.,G+G2/2.,F+F2/2.,FUNT)
      F3=H*FUNT
C
      G4=H*GUNT(X+H,G+G3,F+F3,FUNT)
      F4=H*FUNT
C
      G=G+(G1+2.*(G2+G3)+G4)/6.
      F=F+(F1+2.*(F2+F3)+F4)/6.
      X=X+H
C
      RETURN
      END
C ***********************************************************************
C
      FUNCTION GUNT(X,G1,F1,FUNT)
      COMMON /GRID/XMM,HH,IMAX
      COMMON /SHELLS/ STATE(7,40),LEVEL,EIGEN(40)
C
      CALL FIELD(X,S,V)
      GUNT=F1*(EIGEN(LEVEL)+XMM*2.-S-V)/197.33-STATE(2,LEVEL)*G1/X
      FUNT=G1*(V-S-EIGEN(LEVEL))/197.33+STATE(2,LEVEL)*F1/X
      RETURN
      END
C ***********************************************************************
C
      SUBROUTINE FIELD(X,S,V1)
      COMMON /GRID/XMM,HH,IMAX
      COMMON /MESONS/FIELDS(600,4),G2(4),V(600)
      I=IFIX(X/HH+.1)
      S=FIELDS(I,1)
      V1=V(I)
      IF ((X/HH-FLOAT(I)).LT..3)RETURN
      S=(S+FIELDS(I+1,1))*.5
      V1=(V1+V(I+1))*.5
      RETURN
      END
C ***********************************************************************
C
      SUBROUTINE MESINT
C
C CALCULATES MESON GREEN'S FUNCTIONS
C
C XMASS(I)= MS, MV, MRHO
C
C
      COMMON /GRID/XMM,HH,IMAX
      COMMON /GREEN/ GIN(600,4),GOUT(600,4),XMASS(3)
C
      DO 100 I=1,3
100   XMASS(I)=XMASS(I)/197.33
      DO 200 J=1,IMAX
      X=HH*FLOAT(J)
          DO 190 I=1,3
             EX=EXP(XMASS(I)*X)
             GIN(J,I)=.5/XMASS(I)/X*EX
190       GOUT(J,I)=1./X/EX
      GIN(J,4)=1.
      GOUT(J,4)=1./X
200   CONTINUE
      RETURN
      END
C *****************************************************************
C
      SUBROUTINE MESON
C
C CALCULATES MESON FIELDS
C
C
      COMMON /GRID/XMM,HH,IMAX
      COMMON /GREEN/ GIN(600,4),GOUT(600,4),XMASS(3)
      COMMON /MESONS/ FIELDS(600,4),G2(4),V(600)
      COMMON /DENSIT/ XRHO(600,4),XRHOSP(600)
      DIMENSION XI1(600),XI2(600),FIN(600),FOUT(600),FH(600)
C
      DO 100 L=1,4
          DO 20 I=1,IMAX
             FIN(I)=GIN(I,L)*XRHO(I,L)
20        FOUT(I)=GOUT(I,L)*XRHO(I,L)
C
          CALL HALF(FIN,FH)
          XI1(1)=(4.*FH(1)+FIN(1))/6.  
          DO 30 I=2,IMAX
30        XI1(I)=XI1(I-1)+(FIN(I-1)+4.*FH(I)+FIN(I))/6.  
C
          CALL HALF(FOUT,FH)
          XI2(IMAX)=(4.*FH(IMAX)+FOUT(IMAX))/6.  
          DO 40 I=IMAX-1,1,-1
40        XI2(I)=XI2(I+1)+(FOUT(I+1)+4.*FH(I+1)+FOUT(I))/6.  
          XI20=XI2(1)+(4.*FH(1)+FOUT(1))/6.  
C
          IF(L.NE.4)XX=-.5/XMASS(L)
          IF(L.EQ.4)XX=0.  
          XI20=XI20*XX
C
          DO 100 I=1,IMAX
             FIELDS(I,L)=GOUT(I,L)*(XI1(I)+XI20)+GIN(I,L)*XI2(I)
100   FIELDS(I,L)=FIELDS(I,L)*HH*G2(L)*197.33
C
      RETURN
      END
C *****************************************************************
C
      SUBROUTINE HALF(Y,YH)
C
C     INTERPOLATES TO FIND INTERMEDIATE Y VALUES
C
      DIMENSION Y(600), YH(600)
      COMMON /GRID/ XMM,HH,IMAX
      DATA FAC1,FAC2/ -0.0625,0.5625/
      DO 10 I=3,IMAX-1
10    YH(I)=FAC1*(Y(I-2)+Y(I+1))+FAC2*(Y(I)+Y(I-1))
      YH(2)=(3.*Y(1)+6.*Y(2)-Y(3))/8.
      YH(IMAX)=(3.*Y(IMAX)+6.*Y(IMAX-1)-Y(IMAX-2))/8.
C     For the last point, assume y(0)=0 and y goes like x**2.          
      YH(1)=Y(1)/4.
      RETURN
      END
C *****************************************************************
C *****************************************************************
