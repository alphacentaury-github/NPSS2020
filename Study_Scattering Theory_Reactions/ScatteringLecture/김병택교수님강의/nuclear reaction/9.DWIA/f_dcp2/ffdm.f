      SUBROUTINE FFCALDM
CCCCCC
      PARAMETER(NHS=9,NPS=9,LXA=150,LXB=150,LPH=10,LMM=4,LTL=9,LM1=5)
      PARAMETER(NXA=420,NXB=140,NIN=300,N1X=140,KGG=125,KAA=40)
      PARAMETER(KFM=10,KCC=80,LMJ=10,KJL=600)
CCCCCC      
      COMMON/CNST  /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE 
      COMMON/DWCC/  LDWMIR(4),LDWMXR(4),WN(4),LDWSTR(4),CE(4),CFUNIR(4),
     1              ECM(4),XMESR(4),NXMXR(4),NXMIR(4),NOEXP,MXSTEP,
     2              TMASA,TMASB,TMASX,TMASY,PMASA,PMASB,PMASX,PMASY,
     3              TZA,TZB,TZX,TZY,
     4              PZA,PZB,PZX,PZY,XBARD(4),XBARID(4),WNXD(20),EXPD(20)
      COMMON/CNTRL/ KTRLD(9),KTLOUT(24),NERUN,KRTYPE
      COMMON/FFCC/  VRANG(6),VSTR(12,6),VV(300,24),XMESH,XMESC,
     1              XMESHD,XMESHE,RHED(300),
     2              NBCMI,NBCMX,NONB,NBSTEP,LASTEP,LBSTEP,MAMAX,MXMAX,
     3              MP1MAX,NOLA,NOLB,INTRAN,LMI,LMX,KCETN(2),
     4              NCMAX,NCMIN,LDMAX,JATRTY,NONA,NONAH,NONAR,NASTEP,
     5              NHDMX,NHEMX,NGAUSR,NBSTPD,NBSTPE,N1STEP,JATW,ISATW
      COMMON/SPSTAT/EET,SPECTA(LPH),ESP(LPH),ESH(LPH),NSP(LPH),LSP(LPH),
     1              JTWP(LPH),NSH(LPH),LSH(LPH),JTWH(LPH),IPAIR,NOSP,
     2              NOSH,NPST(LPH),NHST(LPH),NPAIRD,L1R(8),ISR(8),
     3              ITR(8),NLSMAX,JT,KPARIT,IST,NOLTR,N1,LTR(8),
     4              ITP(LPH),ITH(LPH),ITZP(LPH),ITZH(LPH)
      COMMON/CDENS/ RHOD(NXA,NIN,LMM),USAVEX(4*NXA),SPEC(4),
     1              DENSTY(NXA),LAMMXD(3),LBXD(4),NOSTX,NDMY
      COMMON/CBST/  USAVP(NPS*NXA),USAVH(NHS*NXA),CY(LXA)
      COMPLEX       CY
      COMMON/CCKF/  L1KFD(KFM),ISKFD(KFM),KKFD(KFM),LTKFD(KFM),KFMAX,
     1              ITKFD(KFM),NLSKFD(KFM),NLTKFD(KFM),LTMIM,LTMAX
      COMMON/CXLM/  XLMRI(LTL,LM1,NXA),LLP1MX
      COMPLEX       XLMRI 
      COMMON/CGGR/  GGRI(KCC,N1X,LM1)
      COMPLEX       GGRI 
      COMMON/OVDE/  OVDD(20,20),OVED(20,20),OVDDP(10,10),OVDEP(10,10)   
      COMPLEX       OVDD,OVED,OVDDP,OVDEP                               
      COMMON/SPECFC/ALPHA(30,6)
      COMMON/CTFAC/ TFAC(3,10,10),MJMAX                                 
      COMMON/CGFAC /KAMAX,KGMAX,KCMAX,JLSMX,KBBSD(10),
     1              KBC(KCC),LAM2PC(KCC),LPP1C(KCC),LLP1D(KCC),
     2              NINTC(KCC),KAD(KGG),KCD(KGG),NLSD(KGG),LAMD(KAA),
     3              CLEBQ(KCC,5),CLEBF(KAA,5,5),COEF(KGG),CCDD(KJL) 
      COMMON/CSOUR/ SOURCE(NXA,LPH,LMJ)   
      COMPLEX       SOURCE
      DIMENSION     P(130,30),PP(130,10)
      DIMENSION     RTS(96),WGT(96)
      DIMENSION     L9(10)
      DIMENSION     AF(500),AFI(500),XAF(500),XIAF(500),FSP(500),
     1              FSPI(500),Q1(500),AU(500)
      DIMENSION     QA(KFM,NIN,LM1)
      COMPLEX       QA      
      DIMENSION     RHOT(NIN,LMM),TRHO(NXA,5),PD(LMM*LMM)
      DIMENSION     TRHOPH(NXA,5,LPH)                                   
      COMPLEX       ZERO,TTI,URI1,URI2,URI3
      DOUBLE PRECISION   FFMU,RROOT,FFMU2
CCCCCC
      WRITE(8,830)
  830 FORMAT(//1H ,'FFCALD IS CALLED')
      WRITE(8,831)
  831 FORMAT(1H ,5X,'INPUT INFORMATION FOR DIRECT FORM FACTOR')
      WRITE(8,832) MXMAX
  832 FORMAT(1H ,5X,'MXMAX=',5I5)
      WRITE(8,833) NXMIR(4),NXMXR(4),N1STEP,XMESR(4)
  833 FORMAT(1H ,5X,'N1MIN,N1MAX,N1STEP,XEMSH=',3I5,F10.3)
      IF (N1STEP.NE.1) WRITE(8,834)
  834 FORMAT(1H ,5X,'**INTERPOLATION IN R_1 IS MADE')
      WRITE(8,835) NXMIR(1),NXMXR(1),NASTEP,XMESR(1)
      WRITE(8,836) INTRAN,XMESHD
 836  FORMAT(1H ,5X,'INTRAN,XMESHD=',I5,F10.3)
 835  FORMAT(1H ,5X,'NAMIN,NAMAX,NASTEP,XEMSH=',3I5,F10.3)
CCCCCC 
CCCCCC     SECTION 1. EFFECTIVE INTERACTION
CCCCCC	
      CALL EFFINT(0)
CCCCCC
CCCCCC     SECTION 2. PROJECTILE DENSITY
CCCCCC
      CALL PDENST(0)
CCCCCC
CCCCCC     SECTION 3. TARGET DENSITY
CCCCCC
      ZERO=(0.0,0.0)
      TTI=(0.0,1.0)
      KMAS=PMASA+0.001
      IF (KTRLD(1).EQ.3) KMAS=1
      IF (KMAS.EQ.1) N1MAX=NBCMX
      LAP1MI=LDWMIR(1)+1
      LAP1MX=LDWMXR(1)+1
      LASTEP=LDWSTR(1)
      LBP1MI=LDWMIR(2)+1
      LBP1MX=LDWMXR(2)+1
      LBSTEP=LDWSTR(2)
      NOLA=(LAP1MX-LAP1MI)/LASTEP+1
      NOLB=(LBP1MX-LBP1MI)/LBSTEP+1
      XMESA=XMESR(1)
      NMXA=NXMXR(1)
      NMIA=NXMIR(1)
      XMESX=XMESR(3)
      NMXX=NXMXR(3)
      NMIX=NXMIR(3)
      XMES1=XMESR(4)
      N1MAX=NXMXR(4)
      IF (KMAS.EQ.1) N1MAX=NBCMX
      N1MIN=NXMIR(4)
      XMESB=XMESR(2)
      NMXB=NXMXR(2)
      NMIB=NXMIR(2)
      NONB=(NBCMX-NBCMI)/NBSTEP+1
      NBCMX=(NONB-1)*NBSTEP+NBCMI
      NONA=NMXA-NMIA+1
      NONX=NMXX-NMIX+1
      NONY=N1MAX-N1MIN+1
      XMESX3=XMESX/3.0
      XMESH3=XMESH/3.0
      XMES13=XMES1/3.0
      XMESA3=XMESA/3.0                                                  
      R2MX=XMESX*FLOAT(NMXX)
      IF(KMAS.EQ.1) R2MX=XMESH*FLOAT(INTRAN) !NUCLEON SCATTERING
      R2MXSQ=R2MX*R2MX
CCCCCC
CCCCC        CALCULATIONS OF TRANSITION DENSITIES
CCCCC
      LAM1MX=1
      LAM2MX=3
      KBMAXT=NLSMAX
      DO 50 NLS=1,NLSMAX
      NPAIR=NPAIRD                                                      
      DO 50 N1=1,N1MAX
      TRHO(N1,NLS)=0.0
      DO 45 I=1,NPAIR                                                   
      TRHOPH(N1,NLS,I)=0.0                                              
 45   CONTINUE                                                          
 50   CONTINUE
      JTTW=JT*2
      ITZ1=PZA-PZB+0.01
      DO 80 NLS=1,NLSMAX
      IS=ISR(NLS)
      ISTW=IS*2
      L1=L1R(NLS)
      IT1=ITR(NLS)
      L1TW=L1*2
      NPAIR=NPAIRD
      DO 56 I=1,NPAIR
      NP1=NPST(I)
      NH1=NHST(I)
      A1=SPECTA(I)
      N1BASE=(NP1-1)*N1MAX
      N2BASE=(NH1-1)*N1MAX
      LP=LSP(NP1)
      LH=LSH(NH1)
      K1=LP+LH+L1
      IF(MOD(K1,2).NE.0) GO TO 56
      LPTW=LP*2
      LHTW=LH*2
      JPTW=JTWP(NP1)
      JHTW=JTWH(NH1)
      L9(1)=LPTW
      L9(2)=1
      L9(3)=LHTW
      L9(4)=1
      L9(5)=JPTW
      L9(6)=JHTW
      L9(7)=L1TW
      L9(8)=ISTW
      L9(9)=JTTW
      CALL NINEJ(L9,FACLOG,U9)
      T1=U9
      H1=(JPTW+1)*(JHTW+1)*(L1TW+1)*(ISTW+1)
      U9=T1*SQRT(H1)
      PHASE=1.0
      IF(IT1.EQ.1.AND.ITZ1.EQ.0.AND.ITZP(NP1).EQ.-1)  PHASE=-1.0
      IA=LPTW
      IB=LHTW
      IC=L1TW
      H1=FLOAT((IA+1)*(IB+1))/FLOAT(IC+1)
      HAT=SQRT(H1)
      CALL CLEBZ(IA,IB,IC,FACLOG,RAC)
      C1=RAC
      K1=-KPARIT
      S1=1.0
      FAC1=A1*U9*S1*C1*HAT
      NCORR=LP
      PHACO=1-2*MOD(NCORR,2)
      FAC1=FAC1*PHACO
      FAC1P=U9*S1*C1*HAT*PHASE*PHACO
      DO 40 N1=1,N1MAX
      N1M=N1+N1BASE
      N2M=N1+N2BASE
      TRHO(N1,NLS)=TRHO(N1,NLS)+USAVP(N1M)*USAVH(N2M)*FAC1*PHASE
      TRHOPH(N1,NLS,I)=TRHOPH(N1,NLS,I)+USAVH(N2M)*FAC1*PHASE           
   40 CONTINUE
   56 CONTINUE
CCCCCC
      IF (KTLOUT(3).NE.0) THEN
      WRITE(8,871)
      WRITE(8,872) NLS,(TRHO(N1,NLS),N1=10,N1MAX,10)
      ENDIF
  871 FORMAT(/1H ,'TRHO(N1,NLS)= (AT N1=10,NIMAX,10)')
  872 FORMAT(1H ,'NLS=',I3,10E10.3/(8X,10E10.3))
CCCCCC
   80 CONTINUE
CCCCCC
CCCCCC     SECTION 4. CG COEF. IN TRANSITION AMPLITUDES, Eq.(12)
CCCCCC
      MP1MAX=MIN(LTR(NOLTR),MXMAX)+1
      WRITE(8,873) JT, MJMAX
      LTP1MX=LTR(NOLTR)+1
      MJP1MX=MJMAX+1
      JTTW=JT*2
      IA=2
      IB=JTTW
      DO 85 LTP1=1,LTP1MX
      LT=LTP1-1
      LTTW=LT*2
      IC=LTTW
      DO 84 MJP1=1,MJP1MX
      MJ=MJP1-1
      MJTW=MJ*2
      IE=MJTW
      DO 83 MS=1,3
      ID=(MS-2)*2
      T1=2.0
      IF(MS.EQ.2) T1=1.0
      IFF=IE+ID
      CALL CLEB(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)
      TFAC(MS,MJP1,LTP1)=RAC/SQRT(T1)
   83 CONTINUE
      WRITE(8,874) LT,MJ,(TFAC(MS,MJP1,LTP1),MS=1,3)
   84 CONTINUE
   85 CONTINUE
  873 FORMAT(//1H , 'CG COEF. IN TRANSITION AMPLITUDES, EQ.(12) (MAIN)'/
     1    6X, 'JT, MJMAX=',2I3 )
  874 FORMAT(1H ,'TFAC for LT,MJ=',2I3,3E10.3)
CCCCCC
CCCCCC      SECTION 5. A FACTOR WITH M AND SUMS OVER LT AND L1
CCCCCC
      L1MAX=L1R(1)
      DO 90 NLS=1,NLSMAX
      IF(L1MAX.LT.L1R(NLS)) L1MAX=L1R(NLS)
   90 CONTINUE
      M2MAX=LAM2MX
      LPP1MX=L1MAX+1
      KF=0
      LTMIN=LTR(1)
      LTMAX=LTR(NOLTR)
      NLSKMX=NLSMAX*2
      DO 95 NLT=1,NOLTR
      LT=LTR(NLT)
      LTTW=LT*2
      DO 94 NLSK=1,NLSKMX
      K1=(NLSK-1)/NLSMAX
      K=K1*2
      NLS=NLSK-NLSMAX*K1
      IS=ISR(NLS)
      L1=L1R(NLS)
      KTWP1=K*2+1
      IF(LT.GT.(K+L1).OR.LT.LT.IABS(K-L1)) GO TO 94
      IF(K.EQ.2.AND.IS.EQ.0)               GO TO 94
      KF=KF+1
      ITKFD(KF)=ITR(NLS)
      L1KFD(KF)=L1
      ISKFD(KF)=IS
      KKFD(KF)=K
      LTKFD(KF)=LT
      NLSKFD(KF)=NLS
      NLTKFD(KF)=NLT
      T1=LT*2+1
      T2=KTWP1
      FAC1=PI*2.0*SQRT(T1)/T2
      IF(MOD(LT+L1+K,2).NE.0) FAC1=0.0
      M1X=MIN0(K+1,L1+1)
      F1=1.0
      DO 92 M1=1,M1X
      IA=LT*2
      IB=L1*2
      IC=K *2
      ID=0
      IE=(M1-1)*2
      IFF=IE
      CALL CLEB(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)
      IF(M1.GT.1) F1=2.0
      CLEBQ(KF,M1)=RAC*FAC1*F1*ALPHA(NLSK,NLT)
   92 CONTINUE
   94 CONTINUE
   95 CONTINUE
      KFMAX=KF
CCCCCC
      IF (KTLOUT(3).NE.0) THEN
      WRITE(8,875) KFMAX,LTMIN,LTMAX,LAM2MX,LPP1MX,MP1MAX
  875 FORMAT(/1H ,'KFMAX,LTMIN,LTMAX,LAM2MX,LPP1MX,MP1MAX=',6I4)
      ENDIF
CCCCC
CCCCC     SECTION 6. PREPARE GAUSSIAN INTEGRATION
CCCCC
      NGAUS=NGAUSR
      CALL GAUSF(NGAUS,RTS,WGT)
CCCCC
CCCCC     SECTION 7. DIRECT FORM FACTOR 
CCCCC                [BIG N1-LOOP STARTS(UP TO 400)]
CCCCC
      DO 400 N1=N1STEP,N1MAX,N1STEP
      N1M=N1/N1STEP
      R1=XMES1*N1
      R1SQ=R1*R1
      DO 102 KF=1,KFMAX
      DO 102 NH=1,INTRAN
      DO 102 M1=1,MP1MAX
      QA(KF,NH,M1)=0.0
  102 CONTINUE
CCCCC
CCCCC             [NB-LOOP(UP TO 300)]
CCCCC
      F1=2.0
      XMESB3=XMESB/3.0
      DO 300 NB=NMIA,NMXA
      RB=XMESB*FLOAT(NB)
      RBSQ=RB*RB
      AB1=ABS(R1-RB)
      IF(AB1.GE.R2MX) GO TO 300
      F1=6.0-F1
      FAC1=XMESB3*F1
      R1RB2=R1*RB*2.0
CCCCC
      IF (KMAS.EQ.1) GO TO 200
CCCCC
      FMUI=(R1SQ+RBSQ-R2MXSQ)/R1RB2
      IF(FMUI.GT.1.0) GO TO 300
      IF(FMUI.LT.-1.0) FMUI=-1.0
      RANGE=1.0-FMUI
      RANG2=RANGE*0.5
CCCCCC
CCCCCC      [GAUSSIAN INTEGRATION OF EQ.(26), (UP TO 190)]
CCCCCC
      DO 190 NG=1,NGAUS
      ROOT=1.0-(1.0-RTS(NG))*RANG2
      WGT1=WGT(NG)*RANG2
      FAC2=WGT1*FAC1
      R2SQ=R1SQ+RBSQ-R1RB2*ROOT
      R2=SQRT(R2SQ)
      X2=R2/XMESX
      N2=X2
      IF(N2.LT.1) N2=1
      IF(N2.GT.NMXX-1) GO TO 190
      FN2=N2
      FMU=(R1*ROOT-RB)/R2
      FFMU=FMU  
      CALL YLCAL(FFMU,LAM2MX,LAM2MX,P)
      J=0
      DO 124 K1=1,LAM2MX
      DO 124 M1=1,K1
      J=J+1
      PD(J)=P(K1,M1)
  124 CONTINUE
      RROOT=ROOT 
      CALL YLCAL(RROOT,LPP1MX,LPP1MX,P)
      DO 125 LAM2P1=1,LAM2MX
      DO 125 NH=1,INTRAN
      RHOT1=RHOD(N2  ,NH,LAM2P1)
      RHOT2=RHOD(N2+1,NH,LAM2P1)
      RHOT(NH,LAM2P1)=RHOT1+(RHOT2-RHOT1)*(X2-FN2)
  125 CONTINUE
CCCCCC
      DO 130 KF=1,KFMAX
      NLS=NLSKFD(KF)
      NLT=NLTKFD(KF)
      LT=LTR(NLT)
      LTP1=LT+1
      S1=(-1.0)**((LT-KPARIT)/2)  
      LP1=L1KFD(KF)+1
      KP1=KKFD(KF)+1
      M2X=MIN0(KP1,LP1)
      DO 128 M2=1,M2X
      C1=CLEBQ(KF,M2)
      K2=(KP1*(KP1-1))/2+M2
      FAC3=FAC2*PD(K2)*P(LP1,M2)*C1*S1                               
      DO 127 M1=1,MP1MAX
      DO 126 NH=1,INTRAN
      QA(KF,NH,M1)=QA(KF,NH,M1)+XLMRI(LTP1,M1,NB)*RHOT(NH,KP1)*FAC3
  126 CONTINUE
  127 CONTINUE
  128 CONTINUE
  130 CONTINUE
  190 CONTINUE
      GO TO 300
CCCCC
CCCCC     NUCLEON-NUCLEUS SCATTERING
CCCCC
  200 DO 240 NH=1,INTRAN
      RH=XMESH*FLOAT(NH)
      RHSQ=RH*RH
      FMU=(RBSQ+R1SQ-RHSQ)/R1RB2
      IF(ABS(FMU).GT.1.0) GO TO 240
      WGT1=-1.0/(RB*R1*RH)
      FMU2=(R1*FMU-RB)/RH
      IF(ABS(FMU2).GT.1.0) GO TO 240
      KP1=3
      FFMU2=FMU2  
      CALL YLCAL(FFMU2,KP1,KP1,P)
      DO 212 LAM2P1=1,KP1
      DO 212 M1    =1,LAM2P1
  212 PP(LAM2P1,M1)=P(LAM2P1,M1)
      FFMU=FMU
      CALL YLCAL(FFMU,LPP1MX,LPP1MX,P)
      DO 230 KF=1,KFMAX
      NLS=NLSKFD(KF)
      NLT=NLTKFD(KF)                                                    
      LT=LTR(NLT)                                                       
      LTP1=LT+1                                                         
      S1=(-1.0)**((LT-KPARIT)/2)                                        
      LP1=L1KFD(KF)+1
      KP1=KKFD(KF)+1
      H1=KP1*2-1
      FAC2=FAC1*WGT1*SQRT(H1)
      M2X=MIN0(KP1,LP1)
      DO 228 M2=1,M2X
      C1=CLEBQ(KF,M2)
      FAC3=FAC2*PP(KP1,M2)*P(LP1,M2)*C1*S1                              
      DO 227 M1=1,MP1MAX                                                
      QA(KF,NH,M1)=QA(KF,NH,M1)+XLMRI(LTP1,M1,NB)*FAC3                  
  227 CONTINUE                                                          
  228 CONTINUE
  230 CONTINUE
  240 CONTINUE
CCCCC
CCCCC      RB LOOP ENDS
CCCCC
  300 CONTINUE
CCCCC
      DO 315 M1=1,MP1MAX                                                
      DO 315 KF=1,KFMAX                                                 
  315 GGRI(KF,N1M,M1)=0.0
      DO 320 KF=1,KFMAX
      IT=ITKFD(KF)
      L1=L1KFD(KF)
      LT=LTKFD(KF)
      H1=1.0/FLOAT(LT*2+1)
      S1=(-1.0)**L1
      H2=SQRT(H1)*S1
      K=KKFD(KF)
      NI1=IT*3+(K/2)+IS+1
      NI2=NI1+6
      F1=2.0
      DO 318 M1=1,MP1MAX                                                
      DO 318 NH=1,INTRAN
      RHO1=XMESH*FLOAT(NH)
      RHO1SQ=RHO1*RHO1
      URI2=QA(KF,NH,M1)
      F1=6.0-F1
      FAC1=XMESH3*F1                                                    
      URI1=VV(NH,NI1)+VV(NH,NI2)*TTI
      GGRI(KF,N1M,M1)=GGRI(KF,N1M,M1)+URI2*H2*
     1        RHO1SQ*URI1*FAC1
  318 CONTINUE
  320 CONTINUE
CCCCCC
CCCCC      R1 LOOP ENDS
CCCCC
  400 CONTINUE
CCCCCC
      IF (KTLOUT(3).GE.1) THEN
C      DO 422 M1=1,MP1MAX
      WRITE(8,881) 
      DO 420 N1=1,N1M
      WRITE(8,882)  N1, (GGRI(KF,N1,1),KF=1,KFMAX)
  420 CONTINUE
  422 CONTINUE
      ENDIF
  881 FORMAT(/1H ,'DIRECT FORM FACTORS BEFORE INTERPOLATION, M1=1')
  882 FORMAT(1H ,'N1=',I5,10E10.3)
CCCCC
CCCCC     SECTION 8. OVERLAP INTEGRAL 
CCCCC
      MJP1MX=MJMAX+1                                                    
      MJP1X=MJP1MX*3                                                    
      DO 505 I=1,NPAIR                                                  
      DO 505 M1=1,MJP1X                                                 
      DO 505 N1=1,N1MAX                                                 
      SOURCE(N1,I,M1)=0.0                                               
  505 CONTINUE                                                          
      DO 510 NLT=1,NOLTR
      DO 510 M1=1,MP1MAX
      OVDD(NLT,M1)=ZERO
  510 CONTINUE                       
      DO 520 N1=N1STEP,N1MAX,N1STEP
      N1M=N1/N1STEP
      XAF(N1M)=XMES1*N1
  520 CONTINUE
      NMAX1=N1M
      DO 530 N1=1,N1MAX
      XIAF(N1)=XMES1*N1
  530 CONTINUE
      NGINT=N1MAX
      MJP1MX=MJMAX+1
CCCCCC                                                 
      DO 600 KF=1,KFMAX
      NLT=NLTKFD(KF)
      NLS=NLSKFD(KF)
      DO 555 M1=1,MP1MAX                                                
      DO 540 N1M=1,NMAX1
      AF (N1M)=REAL (GGRI(KF,N1M,M1))
      AFI(N1M)=AIMAG(GGRI(KF,N1M,M1))
  540 CONTINUE
      CALL DSPLS3(XAF,AF ,NMAX1,XIAF,FSP ,NGINT,Q1,AU,1)
      CALL DSPLS3(XAF,AFI,NMAX1,XIAF,FSPI,NGINT,Q1,AU,1)
      DO 542 N1=1,N1MAX                                                 
      QA(M1,N1,1)=FSP(N1)+FSPI(N1)*TTI                                  
  542 CONTINUE                                                          
      F1=2.0
      DO 550 N1=1,N1MAX
      F1=6.0-F1
      R1=XMES1*FLOAT(N1)
      R1SQ=R1*R1
      OVDD(NLT,M1)=OVDD(NLT,M1)+(FSP(N1)+FSPI(N1)*TTI)*
     1                      TRHO(N1,NLS)*R1SQ*F1*XMES13
  550 CONTINUE
  555 CONTINUE
      LT=LTR(NLT)
      LTP1=LT+1                                                         
      DO 580 MJP1=1,MJP1MX                                              
      ML1=MJP1-2                                                        
      ML2=MJP1                                                          
      S1=1.0                                                            
      IF(ML1.LT.0) THEN                                                 
        ML1=-ML1                                                        
        S1=1-MOD(IABS(ML1),2)*2                                         
      ENDIF                                                             
      ML1P1=ML1+1                                                       
      ML2P1=ML2+1                                                       
      MJP11=MJP1+MJP1MX                                                 
      MJP12=MJP1+MJP1MX*2                                               
      DO 570 I=1,NPAIR                                                  
      DO 560 N1=1,N1MAX                                                 
      R1=XMES1*N1                                                       
      URI1=TRHOPH(N1,NLS,I)*QA(ML1P1,N1,1)*S1*R1                        
      URI2=TRHOPH(N1,NLS,I)*QA(ML2P1,N1,1)*R1                           
      URI3=TRHOPH(N1,NLS,I)*QA(MJP1,N1,1)*R1                            
      SOURCE(N1,I,MJP1 )=SOURCE(N1,I,MJP1 )+                            
     1    TFAC(1,MJP1,LTP1)*URI1-TFAC(3,MJP1,LTP1)*URI2                 
      SOURCE(N1,I,MJP11)=SOURCE(N1,I,MJP11)+                            
     1   (TFAC(1,MJP1,LTP1)*URI1+TFAC(3,MJP1,LTP1)*URI2)*(-TTI)         
      SOURCE(N1,I,MJP12)=SOURCE(N1,I,MJP12)+                            
     1    TFAC(2,MJP1,LTP1)*URI3                                        
  560 CONTINUE      
  570 CONTINUE      
  580 CONTINUE      
  600 CONTINUE
CCCCCC
      DO 710 LXYZ=1,3                                                   
      DO 710 MJP1=1,MJP1MX                                              
      OVDDP(LXYZ,MJP1)=0.0                                              
  710 CONTINUE                                                          
      DO 750 LXYZ=1,3                                                   
      DO 740 MJP1=1,MJP1MX                                              
      MJP1P=(LXYZ-1)*MJP1MX+MJP1                                        
      DO 730 I=1,NPAIR                                                  
      NP1=NPST(I)                                                       
      N1BASE=(NP1-1)*N1MAX                                              
      F1=2.0                                                            
      DO 720 N1=1,N1MAX                                                 
      F1=6.0-F1                                                         
      R1=XMES1*N1                                                       
      R1SQ=R1*R1                                                        
      T1=R1  *F1*XMES13                                                 
      N1M=N1+N1BASE                                                     
      OVDDP(LXYZ,MJP1)=OVDDP(LXYZ,MJP1)+SOURCE(N1,I,MJP1P)*USAVP(N1M)*T1
  720 CONTINUE                                                          
  730 CONTINUE                                                          
  740 CONTINUE                                                          
  750 CONTINUE                                                          
CCCCCC
      WRITE(8,890)
      DO 760 NLT=1,NOLTR                                                
      LT=LTR(NLT)                                                       
      WRITE(8,892) LT,(OVDD(NLT,M1),M1=1,MP1MAX)                        
  760 CONTINUE
      IF (KTRLD(6).NE.0) THEN      
      DO 770 LXYZ=1,3                                                   
      WRITE(8,893) LXYZ,(OVDDP(LXYZ,M1),M1=1,MJP1MX)                    
  770 CONTINUE 
      ENDIF                                                         
  890 FORMAT(/1H ,'DIRECT TRANSITION AMPLITUDES')
  892 FORMAT(1H ,'LT=',I2,10E10.3)                         
  893 FORMAT(/1H ,5X,'LXYZ=',I3,6E10.3/(9X,6E10.3))
CCCCCC                         
      RETURN
      END


      SUBROUTINE XLMCAL
CCCCCC
      PARAMETER(NHS=9,NPS=9,LXA=150,LXB=150,LPH=10,LMM=4,LTL=9,LM1=5)
      PARAMETER(NXA=420,NXB=140,NIN=300)
      PARAMETER(KFM=10)
CCCCCC	
      COMMON/CNST  /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
      COMMON/SPSTAT/EET,SPECTA(LPH),ESP(LPH),ESH(LPH),NSP(LPH),LSP(LPH),
     1              JTWP(LPH),NSH(LPH),LSH(LPH),JTWH(LPH),IPAIR,NOSP,
     2              NOSH,NPST(LPH),NHST(LPH),NPAIRD,L1R(8),ISR(8),
     3              ITR(8),NLSMAX,JT,KPARIT,IST,NOLTR,N1,LTR(8),
     4              ITP(LPH),ITH(LPH),ITZP(LPH),ITZH(LPH)
      COMMON/DWCC/  LDWMIR(4),LDWMXR(4),WN(4),LDWSTR(4),CE(4),CFUNIR(4),
     1              ECM(4),XMESR(4),NXMXR(4),NXMIR(4),NOEXP,MXSTEP,
     2              TMASA,TMASB,TMASX,TMASY,PMASA,PMASB,PMASX,PMASY,
     3              TZA,TZB,TZX,TZY,
     4              PZA,PZB,PZX,PZY,XBARD(4),XBARID(4),WNXD(20),EXPD(20)
      COMMON/CNTRL/ KTRLD(9),KTLOUT(24),NERUN,KRTYPE
      COMMON/FFCC/  VRANG(6),VSTR(12,6),VV(300,24),XMESH,XMESC,
     1              XMESHD,XMESHE,RHED(300),
     2              NBCMI,NBCMX,NONB,NBSTEP,LASTEP,LBSTEP,MAMAX,MXMAX,
     3              MP1MAX,NOLA,NOLB,INTRAN,LMI,LMX,KCETN(2),
     4              NCMAX,NCMIN,LDMAX,JATRTY,NONA,NONAH,NONAR,NASTEP,
     5              NHDMX,NHEMX,NGAUSR,NBSTPD,NBSTPE,N1STEP,JATW,ISATW
      COMMON/CDENS/ RHOD(NXA,NIN,LMM),USAVEX(4*NXA),SPEC(4),
     1              DENSTY(NXA),LAMMXD(3),LBXD(4),NOSTX,NDMY
      COMMON/CCKF/  L1KFD(KFM),ISKFD(KFM),KKFD(KFM),LTKFD(KFM),KFMAX,
     1              ITKFD(KFM),NLSKFD(KFM),NLTKFD(KFM),LTMIM,LTMAX
      COMMON/CXLM/  XLMRI(LTL,LM1,NXA),LLP1MX
      COMPLEX       XLMRI
      COMMON/ANGLC/ NTHEB,NDTHE,THEB,DTHEB
      COMMON/DISW / DISWB(LXB*NXB),DISWA(LXA*NXA)
      COMPLEX       DISWA,DISWB
      DIMENSION     P(130,30),PP(130,10)
      COMPLEX       URI1,TTI
      DOUBLE PRECISION TTH,TT1
CCCCCC
      TTI=(0.0,1.0)
      IF(LDMAX.GT.14)   LDMAX=14
      NDWB=15
      LAP1MI=LDWMIR(1)+1
      LAP1MX=LDWMXR(1)+1
      LBP1MI=LDWMIR(2)+1
      LBP1MX=LDWMXR(2)+1
      LASTEP=LDWSTR(1)
      IF(LASTEP.EQ.0)  LASTEP=1
      NOLA=(LAP1MX-LAP1MI)/LASTEP+1
      LAP1MX=LAP1MI+(NOLA-1)*LASTEP
      SQ4PAI=SQRT(PI*4.0)
      KTOUT7=0
      NXMAX=NXMXR(1)
      NXMIN=NXMIR(1)
      NONAW=NXMAX-NXMIN+1
      NONBW=NXMXR(2)-NXMIR(2)+1
      NXMAX2=NXMAX+2
      NXMX=NXMAX2
      NXMN=1
      XMESA=XMESR(1)
      XMESB=XMESR(2)
      CONST=(SQ4PAI**2)/(WN(1)*WN(2))
      LTMAX=LTR(NOLTR)
      LTMIN=LTR(1)
      MP1MAX=MIN0(LTMAX+1,MXMAX+1)
      LBLTMX=NOLTR*NOLB
CCCCC
      TH=THEB*PI/180.0
      TTH=TH
      TT1=DCOS(TTH)
      CALL YLCAL(TT1,LBP1MX,MP1MAX,P)
      L1RMX=L1R(NLSMAX)                                                 
      LLP1MX=L1RMX+LAMMXD(1)+1                                          
      IF(LLP1MX.GT.9) LLP1MX=9                                          
      DO 20 LLP1=1,LLP1MX                                               
      DO 20 M1=1,MP1MAX                                                 
      DO 20 NB=1,NXMAX                                                  
      XLMRI(LLP1,M1,NB)=0.0                                             
   20 CONTINUE
CCCCC
      JLS=0
      LAMI=LAP1MI-1
      LBMI=LBP1MI-1
      ID=0
      LAMOUT=(NOLA+1)/2
      LTOUT=LTKFD(1)
      LAM=0
      DO 78 LAP1=LAP1MI,LAP1MX,LASTEP
         LA=LAP1-1
         LAM=LAM+1
         LATW=LA*2
         IA=LATW
         H1=LATW+1
         HAT=SQRT(H1)
      DO 77 LLP1=1,LLP1MX                                               
         LL=LLP1-1               
         LLTW=LL*2         
         IC=LLTW                                                        
         MP1MXT=MIN0(MP1MAX,LLP1)                                        
      DO 75 LBP1=LBP1MI,LBP1MX
         IF(LBP1.GT.(LBP1MX-2)) GO TO 75
         LBM=LBP1-LBP1MI+1
         LB=LBP1-1
         IF(LB.LT.IABS(LA-LL)) GO TO 75                                    
         IF(LB.GT.LA+LL      )  GO TO 75                                   
         IF(MOD(LA+LB+LL,2).EQ.1) GO TO 75                                 
         JLS=JLS+1
         K1=IABS(LA-LB+LL)/2                                               
         S1=(-1.0)**K1
         IB=LB*2
         H2=IC+1
         H3=(IA+1)*(IB+1)
         H4=SQRT(H3/H2)
         CALL CLEBZ(IA,IB,IC,FACLOG,RAC)
         C1=RAC*H4
         NA2=(LAM-1)*NONAW                                                 
         NB2=(LB-LBMI)*NONBW                                               
         DO 74 M1=1,MP1MXT
            IE=(M1-1)*2
            IFF=IE
            CALL CLEB(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)
            C2=RAC*HAT*S1*C1 *P(LBP1,M1)                                
            DO 70 NB=NXMIN,NXMAX                                       
               NA1=NA2+NB
               NB1=NB2+NB
               RB=XMESB*FLOAT(NB)
               URI1=DISWA(NA1)*DISWB(NB1)/(RB*RB)
               XLMRI(LLP1,M1,NB)=XLMRI(LLP1,M1,NB)+URI1*C2*CONST       
   70       CONTINUE
   74    CONTINUE
   75 CONTINUE
   77 CONTINUE
   78 CONTINUE
      JLSMX=JLS
CCCCCC
      IF (KTLOUT(3).GE.2) THEN
      WRITE(8,879) JLSMX,NOLTR,LLP1MX,MXMAX,MP1MAX
      LT=LTR(1)                                                         
      DO 80 LLP1=1,LLP1MX
      DO 80 NB=10,NXMAX,10      
      WRITE (8,880) NB,LLP1,(XLMRI(LLP1,M1,NB),M1=1,MP1MAX)
   80 CONTINUE
      ENDIF	
  879 FORMAT(1H //1H ,'DISTORTION FACTORS (XLMCAL)'
     1      /1H ,5X,'JLSMX,NOLTR,LLP1MX,MXMAX,MP1MAX=',6I6) 
  880 FORMAT(1H ,'NB,LLP1=',2I3,10E10.3)
CCCCCC
      RETURN
      END