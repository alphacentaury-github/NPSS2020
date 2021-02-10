      SUBROUTINE FFCALEM
CCCCCC
      PARAMETER(NHS=9,NPS=9,LXA=150,LXB=150,LPH=10,LMM=4,LTL=9)
      PARAMETER(NXA=420,NXB=140,NIN=300,NOU=300,N1X=140)
      PARAMETER(KFM=10,KAA=40,KBB=30,KCC=80,KGG=125,KJL=600,KTEN=5000)
      PARAMETER(LMJ=10,LRP=11,L1X=10,L2X=60,LM1=5,LM2=30,LAB=2000)
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
      COMMON/CNORE/ FACNR,LRP1MX
      COMMON/CDENS/ RHOE(NXA,NIN,LMM),USAVEX(4*NXA),SPEC(4),
     1              DENSTY(NXA),LAMMXD(3),LBXD(4),NOSTX,NDMY
      COMMON/SPSTAT/EET,SPECTA(LPH),ESP(LPH),ESH(LPH),NSP(LPH),LSP(LPH),
     1              JTWP(LPH),NSH(LPH),LSH(LPH),JTWH(LPH),IPAIR,NOSP,
     2              NOSH,NPST(LPH),NHST(LPH),NPAIRD,L1R(8),ISR(8),
     3              ITR(8),NLSMAX,JT,KPARIT,IST,NOLTR,N1,LTR(8),
     4              ITP(LPH),ITH(LPH),ITZP(LPH),ITZH(LPH)
      COMMON/CCKF/  L1KFD(KFM),ISKFD(KFM),KKFD(KFM),LTKFD(KFM),KFMAX,
     1              ITKFD(KFM),NLSKFD(KFM),NLTKFD(KFM),LTMIM,LTMAX
      COMMON/CBST/  USAVP(NPS*NXA),USAVH(NHS*NXA),CY(LXA)
      COMPLEX       CY
      COMMON/CXLM/  XLMRI(LTL,LM1,NXA),LLP1MX
      COMPLEX       XLMRI
      COMMON/OVDE/  OVDD(20,20),OVED(20,20),OVDDP(10,10),OVDEP(10,10)   
      COMPLEX       OVDD,OVED,OVDDP,OVDEP                               
      COMMON/SPECFC/ALPHA(30,6)
      COMMON/CGFAC /KAMAX,KGMAX,KCMAX,JLSMX,KBBSD(10),
     1              KBC(KCC),LAM2PC(KCC),LPP1C(KCC),LLP1D(KCC),
     2              NINTC(KCC),KAD(KGG),KCD(KGG),NLSD(KGG),LAMD(KAA),
     3              CLEBQ(KCC,5),CLEBF(KAA,5,5),COEF(KGG),CCDD(KJL) 
      COMMON/CGGR/  GGRIE(L1X,N1X)
      COMPLEX       GGRIE
      COMMON/CTFAC/ TFAC(3,10,10),MJMAX                                 
      COMMON/CTRHO/ TRHOPH(N1X,NIN,KBB,LPH)
      COMMON/CSOUR/ SOURCE(NXA,LPH,LMJ)   
      COMPLEX       SOURCE
      COMMON/DISW / DISWB(LXB*NXB),DISWA(LXA*NXA)
      COMPLEX       DISWA,DISWB
      COMMON/TREDEN/GWT(N1X,NIN,LM2),TRHO(N1X,NIN,KBB)	                       
      DIMENSION     P(130,30),PP(130,10)
      DIMENSION     RTS(96),WGT(96)
      DIMENSION     AF(500),AFI(500),XAF(500),XIAF(500),FSP(500),
     1              FSPI(500),Q1(500),AU(500)
      DIMENSION     AK(500),AI(500),NI1R(500),NI2R(500)
      DIMENSION     QAT(KCC,NOU)
      DIMENSION     QA(KCC,NOU,L1X),RHOT(NOU,LMM)
      DIMENSION     GGRI(KAA,NOU,3,L1X),GGRIPH(KAA,NOU,3,L2X)
      DIMENSION     GGRIEP(L1X,NOU,LPH),KKA(KAA)
      COMPLEX       GGRI,GGRIPH,GGRIEP,QA
      COMPLEX       ZERO,TTI,URI1,URI2,URI3,TRI1
      DIMENSION     FJLD(NIN,LRP),SOURCT(NXA,LPH,LMJ)
      COMPLEX       SOURCT
      DOUBLE PRECISION FFMU,FFMU2,RROOT
CCCCCC
      WRITE(8,800)
 800  FORMAT(//1H ,'FFCALE IS CALLED')
      WRITE(8,801)
 801  FORMAT(/1H ,5X,'INPUT INFORMATION FOR EXCHANGE FORM FACTOR')
      IDUM=1
      WRITE(8,802) (LAMMXD(I),I=1,3),MXMAX,MAMAX
 802  FORMAT(1H ,5X,'LAM1MX,LAM2MX,LAMDMX,MXMAX,MAMAX=',6I5)
      WRITE(8,803) NXMIR(1),NXMXR(1),NASTEP,XMESR(1)
 803  FORMAT(1H ,5X,'NAMIN,NAMAX,NASTEP,XEMSH=',3I5,F10.3)
      IF (NASTEP.NE.1) WRITE(8,804)
 804  FORMAT(1H ,5X,'**INTERPOLATION IN R_a IS MADE')
      WRITE(8,805) NBCMI,NBCMX,NBSTEP,XMESR(2)
 805  FORMAT(1H ,5X,'NBCMI,NBCMX,NBSTPE,XMESH=',3I5,F10.3)
C      IF (NBSTPE.NE.1) WRITE(8,806)
 806  FORMAT(1H ,5X,'**INTERPOLATION IN R_b IS MADE')
      WRITE(8,807) NXMIR(3),NXMXR(3),IDUM,XMESR(3)
 807  FORMAT(1H ,5X,'N2MIN,N2MAX,N2STEP,XEMSH=',3I5,F10.3)
      WRITE(8,808) NXMIR(4),NXMXR(4),N1STEP,XMESR(4)
 808  FORMAT(1H ,5X,'N1MIN,N1MAX,N1STEP,XEMSH=',3I5,F10.3)
      IF (N1STEP.NE.1) WRITE(8,809)
 809  FORMAT(1H ,5X,'**INTERPOLATION IN R_1 IS MADE')
      WRITE(8,810) INTRAN,XMESHE
 810  FORMAT(1H ,5X,'INTRAN,XMESHE=',I5,F10.3)
      WRITE(8,811) (KCETN(N),N=1,2)
 811  FORMAT(/1H ,5X,'KCETN(1),KCETN(2)=',2I3)
      IF (KCETN(1).EQ.1) WRITE(9,812)
 812  FORMAT(1H ,5X,'CENTRAL EXCHANGE FF IS NOT CONSIDERED')
      IF (KCETN(2).EQ.1) WRITE(9,813)
 813  FORMAT(1H ,5X,'TENSOR EXCHANGE FF IS NOT CONSIDERED')
      WRITE(8,814) LRP1MX,FACNR                                      
 814  FORMAT(/1H ,5X,'LRP1MX,FACNR      =',I6,F10.5)
      IF (LRP1MX.EQ.1)WRITE(8,815)
 815  FORMAT(1H ,5X,'NO RECOIL APPROXIMATION IS MADE')
      KMAS=PMASA+0.001
      IF (KTRLD(1).EQ.3) KMAS=1
CCCCCC
CCCCCC      SECTION 1. CHECK THE FIELD LENGTH OF VARIABLES
CCCCCC
      CALL DCHECK2
CCCCCC
CCCCCC      SECTION 2. CALCULATE THE EFFECTIVE INTERACTION FOR EXCHANGE FF
CCCCCC
      CALL EFFINT(1)
CCCCCC
CCCCCC      SECTION 3. CALCULATE THE PROJ. AND TARGET DENSITIES
CCCCCC
      CALL TRECALM
CCCCCC      
      ZERO=(0.0,0.0)
      TTI=(0.0,1.0)
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
      N1MIN=NXMIR(4)
      XMESB=XMESR(2)
      NMXB=NXMXR(2)
      NMIB=NXMIR(2)
      NONB=(NBCMX-NBCMI)/NBSTEP+1
      NBCMX=(NONB-1)*NBSTEP+NBCMI
      KFMAX=NOLTR
      LTMIN=LTR(1)
      LTMAX=LTR(NOLTR)
      LTMXP1=LTMAX+1
      MP1MAX=MIN0(LTMAX+1,MXMAX+1)
CCCCCC
      XMES13=XMES1/3.0
      R2MX=XMESX*FLOAT(NMXX)
      R2MXSQ=R2MX*R2MX
      ALPH=PMASA
      ALPHSQ=ALPH*ALPH
      RHOMX=RHED(INTRAN)
      RHODMX=RHOMX/ALPH
      RHDXSQ=RHODMX**2
      NONBNH=NONB*INTRAN
      NONAM=(NONA-1)/NASTEP+1
      LAM1MX=LAMMXD(1)+1
      LAM2MX=LAMMXD(2)+1
      KBMAXT=(NLSMAX*((LAM1MX+1)**2))/2
      NON1=N1MAX/N1STEP
      NLTM1MX=NOLTR*MP1MAX
      NPAIR=NPAIRD                                                      
      L1MAX=L1R(1)
      DO 92 NLS=1,NLSMAX
        IF(L1MAX.LT.L1R(NLS)) L1MAX=L1R(NLS)
   92 CONTINUE
      M2MAX=LAM2MX
      LPP1MX=L1MAX+LAM1MX
      LAMP1X=LAMMXD(3)+1
      DO 94 NLTM1=1,NLTM1MX                                             
      DO 94 N1M=1,NON1                                                  
        GGRIE(NLTM1,N1M)=(0.0,0.0)
      DO 93 I=1,NPAIR                                                   
        GGRIEP(NLTM1,N1M,I)=(0.0,0.0)                                   
 93   CONTINUE    
 94   CONTINUE
CCCCCC
CCCCCC      SECTION 4. CALCULATE THE BESSEL FUNCTIONS
CCCCCC                 [FOR RECOIL EFFECT IN THE PLANE WAVE APPROX.]
CCCCCC
      IF (FACNR.EQ.0.0) LRP1MX=1
      IF (KTLOUT(3).NE.0) THEN
      WRITE(8,820) LRP1MX
      ENDIF
      DO 96 NH=1,INTRAN
        R1=RHED(NH)
        XM1=R1*FACNR*WN(1)/PMASA                                        
        IF (FACNR.NE.0.0) THEN
        X1=XM1
        CALL BESSEL(X1,LRP1MX,AK,AI)
        DO 95 LRP1=1,LRP1MX
          FJLD(NH,LRP1)=AK(LRP1)
   95   CONTINUE
        ELSE
        FJLD(NH,1)=1.0
        ENDIF
      IF (KTLOUT(3).GE.2) THEN
      IF (MOD(NH,10).EQ.0) THEN
      WRITE(8,821) R1,XM1,(FJLD(NH,LRP1),LRP1=1,LRP1MX)
      ENDIF
      ENDIF
   96 CONTINUE
  820 FORMAT(/1H ,'BESSEL FUNCTIONS FOR PLANE WAVE APPRO. UP TO L=',I3)
  821 FORMAT(1H ,'R1,X1,FJ=',2F10.3,5E10.3/(29X,5E10.3))
CCCCCC
CCCCCC      SECTION 5. CALCULATE THE CENTRAL AND TENSOR EXCHANGE FF
CCCCCC
CCCCCC            KCT=1 IS FOR CENTRAL EXCHANGE FORM FACTOR   
CCCCCC                  A  NARROW RH-RANGE SHOULD BE TAKEN
CCCCCC            KCT=2 IS FOR TENSOR  EXCHANGE FORM FACTOR        
CCCCCC                  A  LONG   RH-RANGE SHOULD BE TAKEN
CCCCCC
      DO 1000 KCT=1,2
        IF(KCETN(KCT).EQ.1) GO TO 1000
        KCT1=KCT
CCCCCC
        IF(KCT.EQ.1.AND.KTRLD(1).EQ.2) LAMMXD(3)=0
        IF(KCT.EQ.2.AND.KTRLD(1).EQ.2) LAMMXD(3)=2
CCCCCC
        NGAUS=NGAUSR
        CALL GAUSF(NGAUS,RTS,WGT)
CCCCCC
	  IF (KCT1.EQ.1) WRITE(8,838) KCT1
	  IF (KCT1.EQ.2) WRITE(8,839) KCT1
        WRITE(8,837) LAM1MX,LAM2MX,LAMMXD(1),LAMMXD(2),NGAUS
  837 FORMAT(1H ,5X,'(FFCALE)LAM1MX,LAM2MX,LAMMXD(1,2),NGAUS=',5I3)
  838 FORMAT(/1H ,5X,'CENTRAL EXCHANGE FORM FACTORS, KCT=', I3)
  839 FORMAT(/1H ,5X,'TENSOR EXCHANGE FORM FACTORS, KCT=', I3)
CCCCCC
        CALL GFACM(KCT1)
CCCCCC
CCCCCC      SECTION 6. CALCULATE G-FACTORS 
CCCCCC                   [BIG N1-LOOP STARTS (UP TO 900)]
CCCCCC
      DO 100 KA=1,KAMAX
        KKA(KA)=1
  100 CONTINUE
      DO 900 N1=N1STEP,N1MAX,N1STEP
        N1M=N1/N1STEP
        R1=XMES1*N1
        R1SQ=R1*R1      
CCCCCC
CCCCCC      CALCULATION OF QA
CCCCCC
        DO 102 KC=1,KCMAX
        DO 102 NH=1,INTRAN
        DO 102 NLTM1=1,NLTM1MX                                           
          QA(KC,NH,NLTM1)=0.0                                            
  102   CONTINUE
CCCCCC
      F1=2.0
      XMESB3=XMESB/3.0
      DO 300 NB=NBCMI,NBCMX
        RB=XMESB*NB
        RBSQ=RB*RB
        AB1=ABS(R1-RB)
        IF(AB1.GE.R2MX) GO TO 300
        F1=6.0-F1
        FAC1=XMESB3*F1
        R1RB2=R1*RB*2.0
        FMUI=(R1SQ+RBSQ-R2MXSQ)/R1RB2
        IF (FMUI.GE.1.0) GO TO 300
        IF (FMUI.LT.-1.0) FMUI=-1.0
        RANGE=1.0-FMUI
        RANG2=RANGE*0.5
        DO 110 KC=1,KCMAX                                               
        DO 110 NH=1,INTRAN                                              
          QAT(KC,NH)=0.0                                                 
  110   CONTINUE                                                         
CCCCCC
        IF (KMAS.NE.1) THEN                                         
CCCCCC
        DO 140 NG=1,NGAUS
          ROOT=1.0-(1.0-RTS(NG))*RANG2
          WGT1=WGT(NG)*RANG2
          FAC2=WGT1*FAC1
          R2SQ=R1SQ+RBSQ-R1RB2*ROOT
          R2=SQRT(R2SQ)
          X2=R2/XMESX
          N2=X2
          IF (N2.LT.1) N2=1
          IF (N2.GT.NMXX-1) GO TO 140
          FN2=N2
          FMU=(R1*ROOT-RB)/R2
	    FFMU=FMU
          CALL YLCAL(FFMU,LAM2MX,LAM2MX,P)
          DO 111 LAM2P1=1,LAM2MX
          DO 111 M1=1,LAM2P1
            PP(LAM2P1,M1)=P(LAM2P1,M1)
  111     CONTINUE
          RROOT=ROOT
          CALL YLCAL(RROOT,LPP1MX,M2MAX,P)
          DO 112 LAM2P1=1,LAM2MX
          DO 112 NH=1,INTRAN
            RHOT1=RHOE(N2  ,NH,LAM2P1)
            RHOT2=RHOE(N2+1,NH,LAM2P1)
            RHOT(NH,LAM2P1)=RHOT1+(RHOT2-RHOT1)*(X2-FN2)
  112     CONTINUE
          DO 130 KC=1,KCMAX
          KB=KBC(KC)
          LAM2P1=LAM2PC(KC)
          LPP1=LPP1C(KC)
            DO 128 M2=1,LAM2P1
              C1=CLEBQ(KC,M2)
              FAC3=FAC2*PP(LAM2P1,M2)*P(LPP1,M2)*C1
              DO 126 NH=1,INTRAN
                T1=RHOT(NH,LAM2P1)*FAC3                                
                QAT(KC,NH)=QAT(KC,NH)+T1                               
  126         CONTINUE
  128       CONTINUE
  130     CONTINUE
  140   CONTINUE
CCCCCC
        ELSE
CCCCC
CCCCC   NUNCLEON-NUCLEUS SCATTERING
CCCCC 
        DO 160 NH=1,INTRAN                                              
          RH=XMESH*FLOAT(NH)                                            
          RHSQ=RH*RH                                                    
          FMU=(RBSQ+R1SQ-RHSQ)/R1RB2                                    
          IF (ABS(FMU).GT.1.0)  GO TO 160                               
          WGT1=1.0/(RB*R1*RH)                                           
          FMU2=(R1*FMU-RB)/RH                                           
          IF (ABS(FMU2).GT.1.0) GO TO 160                               
          KP1=3
          FFMU2=FMU2                                                    
          CALL YLCAL(FFMU2,KP1,KP1,P)                                   
          DO 142 LAM2P1=1,KP1                                           
          DO 142 M2X   =1,LAM2P1                                        
            PP(LAM2P1,M2X)=P(LAM2P1,M2X)                                
  142     CONTINUE
          FFMU=FMU 
          CALL YLCAL(FFMU,LPP1MX,LPP1MX,P)                              
          DO 150 KC=1,KCMAX                                             
            KB=KBC(KC)                                                  
            LAM2P1=LAM2PC(KC)                                           
            H1SQ=LAM2P1*2-1                                             
            H1=SQRT(H1SQ)                                               
            FAC2=FAC1*WGT1*H1                                           
            S1=(-1.0)**(LAM2P1-1)                                       
            LPP1=LPP1C(KC)                                              
            DO 148 M2X=1,LAM2P1                                         
              C1=CLEBQ(KC,M2X)                                          
              FAC3=FAC2*PP(LAM2P1,M2X)*P(LPP1,M2X)*C1*S1                
              QAT(KC,NH)=QAT(KC,NH)+FAC3                       
  148       CONTINUE
  150     CONTINUE
  160   CONTINUE
        ENDIF
        DO 193 KC=1,KCMAX
          LLP1=LLP1D(KC)                                               
          DO 192 NLTM1=1,NLTM1MX                                       
            NLT=(NLTM1-1)/MP1MAX+1                                     
            M1=NLTM1-(NLT-1)*MP1MAX                                    
            DO 191 NH=1,INTRAN
              QA(KC,NH,NLTM1)=QA(KC,NH,NLTM1)
     1			           +QAT(KC,NH)*XLMRI(LLP1,M1,NB)
  191       CONTINUE
  192     CONTINUE
  193   CONTINUE
  300 CONTINUE
CCCCCC
CCCCCC      SECTION 7. CALCULATE THE G-FACTOR
CCCCCC
      DO 311 NLS=1,NLSMAX
      DO 311 KA=1,KAMAX
      DO 311 NH=1,INTRAN
      DO 311 NLTM1=1,NLTM1MX
        GGRI(KA,NH,NLS,NLTM1)=(0.0,0.0)                              
  311 CONTINUE
      DO 312 NLS=1,NLSMAX
      DO 312 KA=1,KAMAX
      DO 312 NH=1,INTRAN
      DO 312 NLTM1=1,NLTM1MX
      DO 312 I=1,NPAIR
        I1=(I-1)*NLTM1MX+NLTM1
        GGRIPH(KA,NH,NLS,I1)=(0.0,0.0)                               
  312 CONTINUE
CCCCCC
      DO 325 NLTM1=1,NLTM1MX                                         
	    DO 320 KG=1,KGMAX
          NLS=NLSD(KG)                                             
          KA =KAD (KG)                                             
          IF (KA.EQ.0) GO TO 320
          COEF1=COEF(KG)
          KKA(KA)=0
          KC=KCD(KG)
          KB=KBC(KC)                                               
          NI1=NINTC(KC)+KCT-1
          NI2=NI1+6
          NI1R(KA)=NI1
          NI2R(KA)=NI2
          DO 318 NH=1,INTRAN
            URI1=VV(NH,NI1)+VV(NH,NI2)*TTI
            TRI1=QA(KC,NH,NLTM1)*COEF1*URI1*TRHO(N1M,NH,KB)          
            GGRI(KA,NH,NLS,NLTM1)=GGRI(KA,NH,NLS,NLTM1)+TRI1         
            DO 316 I=1,NPAIR                                         
              TRI1=QA(KC,NH,NLTM1)*COEF1*URI1*TRHOPH(N1M,NH,KB,I)    
              I1=(I-1)*NLTM1MX+NLTM1                                 
              GGRIPH(KA,NH,NLS,I1)=GGRIPH(KA,NH,NLS,I1)+TRI1         
  316       CONTINUE                                                 
  318     CONTINUE
CCCCCC
          IF (KTLOUT(3).GT.2) THEN
          IF (NB.EQ.40) THEN
          WRITE(8,840) KG,KA,LLP1D(KC),KC,COEF(KG),QA(KC,1,1),
     1                                         GGRI(KA,1,NLS,1)
          ENDIF
	    ENDIF
CCCCCC
  320   CONTINUE
  325 CONTINUE
CCCCCC                                                            
      IF(KTLOUT(3).GT.2) THEN
         IF(MOD(NB,10).EQ.0) THEN
         IPRINT=MIN0(INTRAN,30) 
         DO 350 NLS=1,NLSMAX
         DO 350 KA=1,KAMAX
            LAM=LAMD(KA)
            IF(LAM.EQ.0)
     1      WRITE(8,841) NB,NLS,LAM,(GGRI(KA,NH,NLS,1),NH=1,IPRINT)
  350    CONTINUE
         ENDIF
	ENDIF
  840 FORMAT(1H ,'KG,KA,NLT,KC,COEF,QA,GGRI=',4I4,5E10.3)
  841 FORMAT(1H ,'NB,NLS,LAM=',3I3,10E10.3/(21X,10E10.3))
CCCCCC
      CONST=SQRT(PI*4.0)/3.
      KTRP2=(KCT-1)*2+1
      JLS=0
      DO 430 NLTM1=1,NLTM1MX                                        
        NLT=(NLTM1-1)/MP1MAX+1                                      
        LT=LTR(NLT)                                                 
        DO 428 NLS=1,NLSMAX
          L1=L1R(NLS)                                               
          L1TWP1=L1*2+1                                             
          DO 426 KA=1,KAMAX
            LAM=LAMD(KA)
            LL=KA*2-(LAM+2)*LAM+L1-2                                
            LTTWP1=LTR(NLT)*2+1
            H1SQ=1.0                                    
            H1=SQRT(H1SQ)*CONST
            NLSK=NLSMAX*(KCT-1)+NLS
            H1=H1*ALPHA(NLSK,NLT)
            NMAX1=INTRAN
            KT=(KCT-1)*2
            DO 424 LRP1=1,LRP1MX                                         
              LR=LRP1-1                                                 
              IF (LR.GT.KT+LAM) GO TO 424                               
              IF (LR.LT.IABS(KT-LAM)) GO TO 424                         
              IF (MOD(KT+LAM+LR,2).NE.0) GO TO 424
              IF (FACNR.EQ.0.0.AND.KT.NE.LAM) GO TO 424                      
              JLS=JLS+1                                                 
              C1=CCDD(JLS)                                              
              DO 402 NH=1,INTRAN
                RHO1=RHED(NH)
                XAF(NH)=RHO1
                C2=C1*FJLD(NH,LRP1)                                     
                AF(NH)=  GGRI(KA,NH,NLS,NLTM1)*XAF(NH)*C2              
                AFI(NH)=-GGRI(KA,NH,NLS,NLTM1)*XAF(NH)*C2*TTI          
  402         CONTINUE
CCCCCC
                F1=2.0
                H2=H1*XMESHE
                DO 407 NH=1,NHEMX 
                  F1=6.0-F1
                  XIAF(NH)=NH*XMESHE
                  R1=XIAF(NH)**KTRP2  
                  URI1=(AF(NH) +AFI(NH)* TTI)
                  TRI1=URI1*R1*H2
                  GGRIE(NLTM1,N1M)=GGRIE(NLTM1,N1M)+TRI1*F1         
  407           CONTINUE
CCCCCC
              DO 420 I=1,NPAIR                                      
                I1=(I-1)*NLTM1MX+NLTM1                              
                DO 412 NH=1,INTRAN
                  RHO1=RHED(NH)
                  XAF(NH)=RHO1
                  C2=C1*FJLD(NH,LRP1)                                    
                  AF(NH)=  GGRIPH(KA,NH,NLS,I1)*XAF(NH)*C2              
                  AFI(NH)=-GGRIPH(KA,NH,NLS,I1)*XAF(NH)*C2*TTI          
  412           CONTINUE
                F1=2.0
                H2=H1*XMESHE
                DO 417 NH=1,NHEMX 
                  F1=6.0-F1
                  XIAF(NH)=NH*XMESHE
                  R1=XIAF(NH)**KTRP2
                  URI1=(AF(NH) +AFI(NH)* TTI)
                  TRI1=URI1*R1*H2*F1                                 
                  GGRIEP(NLTM1,N1M,I)=GGRIEP(NLTM1,N1M,I)+TRI1         
  417           CONTINUE
  420         CONTINUE
  424       CONTINUE                                                   
  426     CONTINUE
  428   CONTINUE
  430 CONTINUE
CCCCCC        
      IF (KTLOUT(3).GE.2) THEN 
      WRITE(8,851) N1M,(GGRIE(NLTM1,N1M),NLTM1=1,NLTM1MX)
  851 FORMAT(1H ,'N1M=',I3,6E10.3/(8X,6E10.3))                  
      ENDIF
CCCCCC
  900 CONTINUE
 1000 CONTINUE
CCCCCC
CCCCCC      END OF TWO BIG LOOPS (KCT,N1)
CCCCCC
CCCCCC      SECTION 8. CALCULATE OVERLAP INTEGRALS AND SOURCE FUNTIONS
CCCCCC
      DO 510 NLT=1,NOLTR
      DO 510 M1=1,MP1MAX
        OVED(NLT,M1)=ZERO
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
      DO 600 M1=1,MP1MAX
      DO 600 NLT=1,NOLTR
        NLTM1=(NLT-1)*MP1MAX+M1
        DO 540 N1M=1,NMAX1
          AF (N1M)=REAL (GGRIE(NLTM1,N1M))
          AFI(N1M)=AIMAG(GGRIE(NLTM1,N1M))
  540   CONTINUE
        CALL DSPLS3(XAF,AF ,NMAX1,XIAF,FSP ,NGINT,Q1,AU,1)
        CALL DSPLS3(XAF,AFI,NMAX1,XIAF,FSPI,NGINT,Q1,AU,1)
        F1=2.0
        DO 550 N1=1,N1MAX
          F1=6.0-F1
          R1=XMES1*N1
          R1SQ=R1*R1
          URI1=(FSP(N1)+FSPI(N1)*TTI)
          TRI1=URI1*R1SQ
          OVED(NLT,M1)=OVED(NLT,M1)+TRI1*f1*xmes13
  550   CONTINUE
  600 CONTINUE     
CCCCCC
      MJP1MX=MJMAX+1                                                
      MJP1X=MJP1MX*3                                                
      DO 605 I=1,NPAIR                                             
      DO 605 M1=1,MJP1X                                            
      DO 605 N1=1,N1MAX                                            
        SOURCT(N1,I,M1)=0.0                                         
  605 CONTINUE                                                      
      DO 610 LXYZ=1,3                                              
      DO 610 MJP1=1,MJP1MX                                         
        OVDEP(LXYZ,MJP1)=0.0                                        
  610 CONTINUE                                                      
CCCCCC     
      DO 650 MJP1=1,MJP1MX                                            
        ML1=MJP1-2                                                     
        ML2=MJP1                                                       
        S1=1.0                                                         
        IF (ML1.LT.0) THEN                                             
          ML1=-ML1                                                     
          S1=1-MOD(IABS(ML1),2)*2                                      
        ENDIF                                                          
        ML1P1=ML1+1                                                    
        ML2P1=ML2+1                                                    
        MJP11=MJP1+MJP1MX                                              
        MJP12=MJP1+MJP1MX*2                                            
        DO 648 NLT=1,NOLTR                                            
          LT=LTR(NLT)                                                  
          LTP1=LT+1                                                    
          TF1=TFAC(1,MJP1,LTP1)                                        
          TF2=TFAC(3,MJP1,LTP1)                                        
          TF3=TFAC(2,MJP1,LTP1)                                        
          NLTM1=(NLT-1)*MP1MAX+ML1P1                                   
          NLTM2=(NLT-1)*MP1MAX+ML2P1                                   
          NLTM3=(NLT-1)*MP1MAX+MJP1                                    
          DO 646 I=1,NPAIR                                            
          DO 646 N1=1,NMAX1                                           
            URI1=GGRIEP(NLTM1,N1,I)*S1                                 
            URI2=GGRIEP(NLTM2,N1,I)                                    
            URI3=GGRIEP(NLTM3,N1,I)                                    
            SOURCT(N1,I,MJP1 )=SOURCT(N1,I,MJP1 )+TF1*URI1-TF2*URI2    
            SOURCT(N1,I,MJP11)=SOURCT(N1,I,MJP11)
     1              			+(TF1*URI1+TF2*URI2)*(-TTI)
            SOURCT(N1,I,MJP12)=SOURCT(N1,I,MJP12)+TF3*URI3             
  646     CONTINUE                                                     
  648   CONTINUE                                                       
  650 CONTINUE
CCCCCC                                                                         
      DO 680 LXYZ=1,3                                                 
        DO 678 MJP1=1,MJP1MX                                          
          MJP1P=(LXYZ-1)*MJP1MX+MJP1                                   
          DO 676 I=1,NPAIR                                            
            NP1=NPST(I)                                                
            N1BASE=(NP1-1)*N1MAX                                       
            DO 660 N1M=1,NMAX1                                        
              AF(N1M) =REAL (SOURCT(N1M,I,MJP1P))                      
              AFI(N1M)=AIMAG(SOURCT(N1M,I,MJP1P))                      
  660       CONTINUE                                                   
            CALL DSPLS3(XAF,AF ,NMAX1,XIAF,FSP ,NGINT,Q1,AU,1)
            CALL DSPLS3(XAF,AFI,NMAX1,XIAF,FSPI,NGINT,Q1,AU,1)
            DO 662 N1=1,N1MAX                                         
              R1=XMES1*N1                                              
              SOURCT(N1,I,MJP1P)=(FSP(N1)+FSPI(N1)*TTI)*R1             
  662       CONTINUE                                                   
            F1=2.0                                                     
            DO 664 N1=1,N1MAX                                         
              F1=6.0-F1                                                
              R1=XMES1*N1                                              
              R1SQ=R1*R1                                               
              T1=R1  *F1*XMES13                                        
              N1M=N1+N1BASE                                            
              OVDEP(LXYZ,MJP1)=OVDEP(LXYZ,MJP1)
     1           		  +SOURCT(N1,I,MJP1P)*USAVP(N1M)*T1
  664       CONTINUE                                                   
  676     CONTINUE                                                     
  678   CONTINUE                                                       
  680 CONTINUE                                                         
CCCCCC
      WRITE(8,869) 
      DO 710 NLT=1,NOLTR                                              
        LT=LTR(NLT)                                                    
        WRITE(8,868) LT,(OVED(NLT,M1),M1=1,MP1MAX)               
  710 CONTINUE
      DO 720 LXYZ=1,3                                                 
        WRITE(8,870) LXYZ,(OVDEP(LXYZ,M1),M1=1,MJP1MX)               
  720 CONTINUE
  868 FORMAT(1H ,'LT  =',I3,6E10.3/(9X,6E10.3))
  869 FORMAT(/1H ,'EXCHANGE TRANSITION AMPLITUDES')
  870 FORMAT(1H ,'LXYZ=',I3,6E10.3/(9X,6E10.3))
CCCCCC     
      IF(KTRLD(6).EQ.1) THEN
      DO 730 LXYZ=1,3
      DO 730 MJP1=1,MJP1MX
        MJP1P=(LXYZ-1)*MJP1MX+MJP1
        DO 726 I=1,NPAIR
          DO 722 N1=1,N1MAX
            SOURCE(N1,I,MJP1P)=SOURCE(N1,I,MJP1P)+SOURCT(N1,I,MJP1P)
  722     CONTINUE
  726   CONTINUE
  730 CONTINUE
      ENDIF
CCCCCC
      RETURN 
      END
 

      SUBROUTINE TRECALM
CCCCCC
      PARAMETER(NHS=9,NPS=9,LXA=150,LXB=150,LPH=10,LMM=4,LTL=9)
      PARAMETER(NXA=420,NXB=140,NIN=300,NOU=300,N1X=140)
      PARAMETER(KFM=10,KAA=40,KBB=30,KCC=80,KGG=125,KJL=600,KTEN=5000)
      PARAMETER(LMJ=10,LRP=11,L1X=10,L2X=60,LM1=5,LM2=30,LAB=2000)
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
      COMMON/CGGR/  GGRIE(L1X,N1X)
      COMPLEX       GGRIE
      COMMON/SPECFC/ALPHA(30,6)
      COMMON/CGFAC /KAMAX,KGMAX,KCMAX,JLSMX,KBBSD(10),
     1              KBC(KCC),LAM2PC(KCC),LPP1C(KCC),LLP1D(KCC),
     2              NINTC(KCC),KAD(KGG),KCD(KGG),NLSD(KGG),LAMD(KAA),
     3              CLEBQ(KCC,5),CLEBF(KAA,5,5),COEF(KGG),CCDD(KJL) 
      COMMON/CTRHO/ TRHOPH(N1X,NIN,KBB,LPH)
      COMMON/DISW / DISWB(LXB*NXB),DISWA(LXA*NXA)
      COMPLEX       DISWA,DISWB
      COMMON/TREDEN/GWT(N1X,NIN,LM2),TRHO(N1X,NIN,KBB)	                       
      DIMENSION     L9(10)
CCCCC
CCCCC           CALCULATIONS OF TRANSITION DENSITIES
CCCCC
      IF (KTLOUT(3).NE.0) THEN
      WRITE(8,810)
  810 FORMAT(/1H ,'TRECAL IS CALLED')
      ENDIF
      KMAS=PMASA+0.001
      IF (KTRLD(1).EQ.3) KMAS=1
CCCCCC
      CALL PDENST(1)
CCCCCC
      N1MAX=NXMXR(4)
      N1STEPT=N1STEP
      NON1=(N1MAX-1)/N1STEPT+1
      XMES1=XMESR(4)
      LAM1MX=LAMMXD(1)+1
      LAM2MX=LAMMXD(2)+1
      KBMAXT=(NLSMAX*((LAM1MX+1)**2))/2
CCCCCC
      DO 35 N1P=1,NON1
      DO 35 NH=1,INTRAN
      DO 35 KB=1,KBMAXT
      TRHO(N1P,NH,KB)=0.0
 35   CONTINUE
      NPAIR=NPAIRD                                             
      DO 36 KB=1,KBMAXT                                        
      DO 36 I=1,NPAIR                                          
      DO 36 NH=1,INTRAN                                        
      DO 36 N1P=1,NON1                                         
   36 TRHOPH(N1P,NH,KB,I)=0.0                                  
      NGWT=LAM1MX*(LAM1MX+1)/2
      JTTW=JT*2
      ITZ1=PZA-PZB+0.01
      KB=0
      MEXP=0
      DO 80 NLS=1,NLSMAX
      IS=ISR(NLS)
      ISTW=IS*2
      L1=L1R(NLS)
      IT1=ITR(NLS)
      L1TW=L1*2
      NPAIR=NPAIRD
      KBBSD(NLS)=KB
      DO 78 LAM1P1=1,LAM1MX
      LAM1=LAM1P1-1
      LAM1TW=LAM1*2
      LPP=L1-LAM1-1
      NLPPX=2*LAM1+1
      DO 76 NLPP=1,NLPPX
      LPP=LPP+1
      LPPTW=LPP*2
      IF(MOD(LAM1+LPP,2).NE.KPARIT) GO TO 76                 
      IF(LPP.GT.(L1+LAM1))          GO TO 76                 
      IF(LPP.LT.IABS(L1-LAM1))      GO TO 76                 
      IF(LPP.LT.0) GO TO 76
      KB=KB+1
      MEXP=MEXP+1
      IRU=0
      DO 70 I=1,NPAIR
      NP1=NPST(I)
      NH1=NHST(I)
      A1=SPECTA(I)
      IRU=IRU+1
      N1BASE=(NP1-1)*N1MAX
      N2BASE=(NH1-1)*N1MAX
      LP=LSP(NP1)
      LH=LSH(NH1)
      LPTW=LP*2
      LHTW=LH*2
      JPTW=JTWP(NP1)
      JHTW=JTWH(NH1)
      NNN=NH1
      LLL=LH
      LLP=LP
      LE1=LLL
      LLLTW=LLL*2
      LLPTW=LLP*2
      IF(IRU.EQ.1) GO TO 53
      II=I-1
      NEXP=NHST(II)
      IF(NH1.EQ. NEXP) GO TO 59
   53 CALL DENSTM(NNN,LLL)
   59 CONTINUE
      PHASE=1.0
      IF(IT1.EQ.1.AND.ITZ1.EQ.0.AND.ITZP(NP1).EQ.-1) PHASE=-1.0
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
      K11=LP+LH+KPARIT
      K1=-KPARIT
      IF(MOD(K11,2).NE.0) GO TO 70
      S1=1.0
      LHP=LLL-LAM1-2
      NLH=LAM1P1
      DO 68 NLH1=1,LAM1P1
      LHP=LHP+2
      S2=(-1.0)**LHP
      NGW=(LAM1*LAM1P1)/2+NLH1
      IF(LHP.LT.0) GO TO 68
      IF(LHP.GT.(LLP+LPP).OR.LHP.LT.IABS(LLP-LPP)) GO TO 68
      LHPTW=LHP*2
      IA=LPPTW
      IB=LHPTW
      IC=LLPTW
      CALL CLEBZ(IA,IB,IC,FACLOG,RAC) 
      C1=RAC
      H1=(LHPTW+1)*(LLLTW+1)*(LPPTW+1)
      C2=SQRT(H1)
      IA=LPPTW
      IB=LHPTW
      IC=L1TW
      ID=LLLTW
      IE=LLPTW
      IFF=LAM1TW
      CALL RAC7(IA,IB,IC,ID,IE,IFF,FACLOG,RAC) 
      FAC1=A1*U9*S1*S2*C1*C2*RAC
      NCORR=LP
      PHACO=1-2*MOD(NCORR,2)
      IF(KMAS.EQ.1) PHACO=PHACO*(-1.0) !NUCLEON SCATTERING
      FAC1=FAC1*PHACO
      DO 66 NH=1,INTRAN
      DO 66 N1P=1,NON1
      N1=N1P*N1STEPT
      N1M=N1+N1BASE
      N2M=N1+N2BASE
      U2=USAVP(N1M)
      TRHO(N1P,NH,KB)=TRHO(N1P,NH,KB)+U2*GWT(N1P,NH,NGW)*FAC1*PHASE
      TRHOPH(N1P,NH,KB,I)=TRHOPH(N1P,NH,KB,I)                          
     1                    +GWT(N1P,NH,NGW)*FAC1*PHASE                  
   66 CONTINUE
   68 CONTINUE
   70 CONTINUE
   76 CONTINUE
   78 CONTINUE
   80 CONTINUE
      KBMAX=KB
CCCCCC
      IF (KTLOUT(3).NE.0) THEN
      ISTP=(INTRAN+9)/10
      DO 88 KB=1,KBMAX
      WRITE (8,881) KB,INTRAN,ISTP
      DO 86 N1=1,5,4
      WRITE (8,882) N1,(TRHO(N1,NH,KB),NH=1,INTRAN,ISTP)
   86 CONTINUE
   88 CONTINUE
      ENDIF
  881 FORMAT(/1H ,'TRANSITION DENSITY FOR KB=',I3,2X,'INTRAN,ISTP=',
     1      2I3/1H )
  882 FORMAT(1H ,'TRHO AT N1=',I4,10E10.3)
CCCCCC
      RETURN
      END


      SUBROUTINE DENSTM(NNN,LLL)
CCCCCC
      PARAMETER(NHS=9,NPS=9,LXA=150,LXB=150,LPH=10,LMM=4,LTL=9)
      PARAMETER(NXA=420,NXB=140,NIN=300,NOU=300,N1X=140)
      PARAMETER(KFM=10,KAA=40,KBB=30,KCC=80,KGG=125,KJL=600,KTEN=5000)
      PARAMETER(LMJ=10,LRP=11,L1X=10,L2X=60,LM1=5,LM2=30,LAB=2000)
CCCCCC
      COMMON/CNST  /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
      COMMON/BSX/   URSAVE(500),ENEPT,VNEPT,VSOR,DFNR,DFNSO,RZR,RZSO,
     1              RZC,XMES2,Q,FTR,NODER,LBTR,JBTRTW,NXRAWF,KTRL4,
     2              KOPOT,LBTRY,JBTWY
      COMMON/DWCC/  LDWMIR(4),LDWMXR(4),WN(4),LDWSTR(4),CE(4),CFUNIR(4),
     1              ECM(4),XMESR(4),NXMXR(4),NXMIR(4),NOEXP,MXSTEP,
     2              TMASA,TMASB,TMASX,TMASY,PMASA,PMASB,PMASX,PMASY,
     3              TZA,TZB,TZX,TZY,
     4              PZA,PZB,PZX,PZY,XBARD(4),XBARID(4),WNXD(20),EXPD(20)
      COMMON/OPTL/  THMIN,THMAX,THIND
      COMMON/CKXKB/ SX,SA,TX,TA,EB,ROOTD(16),COSMI,PB,PX,PBTX,PXSX,
     1              WNB,WNX,RMASA,RMASB,RANG2,Q3,
     2              VC,VXA,VBL,EBL,RX,CSTHBL
      COMMON/POTCC/ VA,VB,VX,VY,WA,WB,WX,WY ,WAS,WBS,WXS,WYS ,ARA,ARB,
     1              ARX,ARY,AIA,AIB,AIX,AIY,AISA,AISB,AISX,AISY,RZRA,
     2              RZRB,RZRX,RZRY,RZIA,RZIB,RZIX,RZIY,RZISA,RZISB,
     3              RZISX,RZISY,RZCA,RZCB,RZCX,RZCY,IDA,IDB,IDX,IDY,
     4              VRIT(900),NXMN,NXMX
      COMPLEX       VRIT
      COMMON/ANGLC/ NTHEB,NOTHE,THEB,DTHEB
      COMMON/CNTRL/ KTRLD(9),KTLOUT(24),NERUN,KRTYPE
      COMMON/FFCC/  VRANG(6),VSTR(12,6),VV(300,24),XMESH,XMESC,
     1              XMESHD,XMESHE,RHED(300),
     2              NBCMI,NBCMX,NONB,NBSTEP,LASTEP,LBSTEP,MAMAX,MXMAX,
     3              MP1MAX,NOLA,NOLB,INTRAN,LMI,LMX,KCETN(2),
     4              NCMAX,NCMIN,LDMAX,JATRTY,NONA,NONAH,NONAR,NASTEP,
     5              NHDMX,NHEMX,NGAUSR,NBSTPD,NBSTPE,N1STEP,JATW,ISATW
      COMMON/CDENS/ RHOD(NXA,NIN,LMM),USAVEX(4*NXA),SPEC(4),
     1              DENSTY(NXA),LAMMXD(3),LBXD(4),NOSTX,NDMY
      COMMON/SPSTAT/EET,SPECTA(LPH),ESP(LPH),ESH(LPH),NSP(LPH),LSP(LPH),
     1              JTWP(LPH),NSH(LPH),LSH(LPH),JTWH(LPH),IPAIR,NOSP,
     2              NOSH,NPST(LPH),NHST(LPH),NPAIRD,L1R(8),ISR(8),
     3              ITR(8),NLSMAX,JT,KPARIT,IST,NOLTR,N1,LTR(8),
     4              ITP(LPH),ITH(LPH),ITZP(LPH),ITZH(LPH)
      COMMON/CBST/  USAVP(NPS*NXA),USAVH(NHS*NXA),CY(LXA)
      COMPLEX       CY
      COMMON/DISW / DISWB(LXB*NXB),DISWA(LXA*NXA)
      COMPLEX       DISWA,DISWB
      COMMON/TREDEN/GWT(N1X,NIN,LM2),TRHO(N1X,NIN,KBB)	                       
      DIMENSION     P(130,30),PP(130,10)
      DIMENSION     RTS(96),WGT(96)
      DIMENSION     SSUM(20),HATD(20),C(100),PD(40,10,10)
      DOUBLE PRECISION RROOT,FFMU
CCCCCC
      KMAS=PMASA+0.001
      IF (KTRLD(1).EQ.3) KMAS=1
      LAMMAX=LAMMXD(1)
      SQPAI=SQRT(PI)
      TWPAI=PI*2.0
      SQRT2=1.0/SQRT(2.0)
CCCCCC
      NGAUS=NGAUSR
      CALL GAUSF(NGAUS,RTS,WGT)
CCCCCC
      XMES2=XMESR(4)
      I=1
      LAMP1X=LAMMAX+1
      N2MAX=NXMXR(4)
      XMES2P=XMESR(4)
      N2STEP=N1STEP
      NON2=(N2MAX-1)/N2STEP+1
      JLSMAX=(LAMP1X*(LAMP1X+1))/2
      N2MXM1=N2MAX-1
      N2MXM2=N2MAX-2
      DO 13 N2=1,NON2
      DO 13 NH=1,INTRAN
      DO 13 JLS=1,JLSMAX
        GWT(N2,NH,JLS)=0.0
   13 CONTINUE
CCCCCC
      NN=NNN
      LBTR=LLL
      LP1MX=LBTR+1
      MP1MX=LP1MX
      FLBTW1=LBTR*2+1
      K=0
      JLS=0
      DO 20 LAMP1=1,LAMP1X
        LAMD=LAMP1-1
        NLDMX=LAMP1
        HLAMD=TWPAI/FLBTW1
        LD=LBTR-LAMD-2
        DO 19 NLD=1,NLDMX
          JLS=JLS+1
          LD=LD+2
          T2=FLOAT(LD*2+1)
          H1=SQRT(ABS(T2))
          HAT=H1*SQRT(FLBTW1)
          HATD(JLS)=HAT
          LDP1=LD+1
          IA=LD*2
          IB=LAMD*2
          IC=LBTR*2
          ID=0
          T1=1.0
          DO 18 M1=1,LAMP1
            K=K+1
            C(K)=0.0
            IF (LD.LT.0) GO TO 18
            IE=(M1-1)*2
            IFF=IE
            CALL CLEB(IA,IB,IC,ID,IE,IFF,FACLOG,RAC) 
            S1=1.D0
            IF (KMAS.EQ.1) S1=1-2*MOD(M1-1,2) !NUCLEON SCATTERING
            C(K)=RAC*T1*H1*HLAMD*S1
            T1=2.0
   18     CONTINUE
   19   CONTINUE 
   20 CONTINUE
      S1=1.0
      IF (KMAS.EQ.1) S1=-1 !NUCLEON SCATTERING
      DO 25 NG=1,NGAUS
        ROOT=RTS(NG)
        RROOT=ROOT
        CALL YLCAL(RROOT,LAMP1X,LAMP1X,P)
        DO 24 LAMP1=1,LAMP1X
        DO 24 M1=1,   LAMP1
          PD(NG,LAMP1,M1)=P(LAMP1,M1)
   24   CONTINUE
   25 CONTINUE
      N2BASE=(NN-1)*N2MAX
      RH=0.0
      DO 80 NH=1,INTRAN
        RH=RHED(NH)
        RHSQ=RH*RH
        DO 70 N2P=1,NON2
          N2PP=N2P*N2STEP
          R2P=XMES2P*N2PP
          R2PSQ=R2P*R2P
          R2PRH=R2P*RH*2.0*S1
          DO 30 JLS=1,JLSMAX
            SSUM(JLS)=0.0
   30     CONTINUE
          DO 60 NG=1,NGAUS
            ROOT=RTS(NG)
            R2=SQRT(RHSQ+R2PSQ+R2PRH*ROOT)
            X2=R2/XMES2
            N2=X2
            IF (N2.LT.1) N2=1
            IF (N2.GT.N2MXM2) GO TO 60
            FN2=N2
            P2=X2-FN2
            P21=P2-1.0
            P22=P2-2.0
            A0=P21*P22*0.5
            A1=-P2*P22
            A2=P2*P21*0.5
            N3=N2+N2BASE
            U2=USAVH(N3)*A0+USAVH(N3+1)*A1+USAVH(N3+2)*A2
            U3=U2*WGT(NG)
            FMU=(R2P+RH*ROOT*S1)/R2
            IF (ABS(FMU).GT.1.0) GO TO 60
	      FFMU=FMU
            CALL YLCAL(FFMU,LP1MX,MP1MX,P)
            J=0
            K=0
            DO 40 LAMP1=1,LAMP1X
            DO 40 NLD=1,LAMP1
              J=J+1
              DO 39 M1=1,LAMP1
                K=K+1
              SSUM(J)=SSUM(J)+C(K)*U3*P(LP1MX,M1)*PD(NG,LAMP1,M1)
   39         CONTINUE
   40       CONTINUE
   60     CONTINUE
          JLS=0
          DO 65 LAMP1=1,LAMP1X
          DO 65 NLD=1,LAMP1
            JLS=JLS+1
            GWT(N2P,NH,JLS)=SSUM(JLS)
   65     CONTINUE
   70   CONTINUE
   80 CONTINUE
CCCCCC
      RETURN
      END

      
      SUBROUTINE GFACM(KCT)
CCCCCC
      PARAMETER(NHS=9,NPS=9,LXA=150,LXB=150,LPH=10,LMM=4,LTL=9)
      PARAMETER(NXA=420,NXB=140,NIN=300,NOU=300,N1X=140)
      PARAMETER(KFM=10,KAA=40,KBB=30,KCC=80,KGG=125,KJL=600,KTEN=5000)
      PARAMETER(LMJ=10,LRP=11,L1X=10,L2X=60,LM1=5,LM2=30,LAB=2000)
CCCCCC	
      COMMON/CNST  /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
      COMMON/DWCC/  LDWMIR(4),LDWMXR(4),WN(4),LDWSTR(4),CE(4),CFUNIR(4),
     1              ECM(4),XMESR(4),NXMXR(4),NXMIR(4),NOEXP,MXSTEP,
     2              TMASA,TMASB,TMASX,TMASY,PMASA,PMASB,PMASX,PMASY,
     3              TZA,TZB,TZX,TZY,
     4              PZA,PZB,PZX,PZY,XBARD(4),XBARID(4),WNXD(20),EXPD(20)
      COMMON/POTCC/ VA,VB,VX,VY,WA,WB,WX,WY ,WAS,WBS,WXS,WYS ,ARA,ARB,
     1              ARX,ARY,AIA,AIB,AIX,AIY,AISA,AISB,AISX,AISY,RZRA,
     2              RZRB,RZRX,RZRY,RZIA,RZIB,RZIX,RZIY,RZISA,RZISB,
     3              RZISX,RZISY,RZCA,RZCB,RZCX,RZCY,IDA,IDB,IDX,IDY,
     4              VRIT(900),NXMN,NXMX
      COMPLEX       VRIT
      COMMON/CNTRL/ KTRLD(9),KTLOUT(24),NERUN,KRTYPE
      COMMON/FFCC/  VRANG(6),VSTR(12,6),VV(300,24),XMESH,XMESC,
     1              XMESHD,XMESHE,RHED(300),
     2              NBCMI,NBCMX,NONB,NBSTEP,LASTEP,LBSTEP,MAMAX,MXMAX,
     3              MP1MAX,NOLA,NOLB,INTRAN,LMI,LMX,KCETN(2),
     4              NCMAX,NCMIN,LDMAX,JATRTY,NONA,NONAH,NONAR,NASTEP,
     5              NHDMX,NHEMX,NGAUSR,NBSTPD,NBSTPE,N1STEP,JATW,ISATW
      COMMON/CXLM/  XLMRI(LTL,LM1,NXA),LLP1MX
      COMPLEX       XLMRI
      COMMON/CNORE/ FACNR,LRP1MX
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
      COMMON/CGFAC /KAMAX,KGMAX,KCMAX,JLSMX,KBBSD(10),
     1              KBC(KCC),LAM2PC(KCC),LPP1C(KCC),LLP1D(KCC),
     2              NINTC(KCC),KAD(KGG),KCD(KGG),NLSD(KGG),LAMD(KAA),
     3              CLEBQ(KCC,5),CLEBF(KAA,5,5),COEF(KGG),CCDD(KJL) 
      COMMON/SPECFC/ALPHA(30,6)
      COMMON/DISW / DISWB(LXB*NXB),DISWA(LXA*NXA)
      COMPLEX       DISWA,DISWB
CCCCC
CCCCC               DEFINE BASIC VARIABLES AND CONSTANTS
CCCCC
      KCT1=KCT
      LAP1MI=LDWMIR(1)+1
      LAP1MX=LDWMXR(1)+1
      LASTEP=LDWSTR(1)
      LBP1MI=LDWMIR(2)+1
      LBP1MX=LDWMXR(2)+1
      LBSTEP=LDWSTR(2)
      NOLA=(LAP1MX-LAP1MI)/LASTEP+1
      NOLB=(LBP1MX-LBP1MI)/LBSTEP+1
      LAM1MX=LAMMXD(1)+1
      LAM2MX=LAMMXD(2)+1
      LAMP1X=LAMMXD(3)+1
      L1MAX=L1R(NLSMAX)
      LAM3MX=LAM1MX+LAM2MX-1
CCCCCC
      IF (KTLOUT(3).NE.0) THEN
      WRITE(8,825) LAM3MX,LAMP1X,LRP1MX,LLP1MX
  825 FORMAT(1H  ,5X,'(GFAC)LAM3MX,LAMP1X,LRP1MX,LLP1MX=',
     1       6I3)         
      ENDIF
      FPISQI=SQRT(1.0/(4.0*PI))
CCCCC
CCCCC                   CALCULATE G-FACTORS
CCCCC
      DO 10 I=1,KGG
        KAD(I)=0
   10 CONTINUE
      DO 11 I=1,KAA
      DO 11 J=1,5
      DO 11 K=1,5
        CLEBF(I,J,K)=0.
   11 CONTINUE
CCCCCC
      IF(KTLOUT(4).EQ.1) WRITE(8,830)
  830 FORMAT(//1H  ,'   KG  KC  KA  KB NLS LA1 LPP LA2  LL   K LT LAM',
     1              '    COEF    KA NLS    CLEBF')
      FAC1=PI*2.0
      KB=0
      KG=0
      KC=0
      KAMAX=1
      DO 120 NLS=1,NLSMAX
        L1=L1R(NLS)
        L1TW=L1*2
        IS=ISR(NLS)
        IT=ITR(NLS)
        KB=KBBSD(NLS)
        LLL1=0
        DO 100 LAM1P1=1,LAM1MX
          LLL1=LLL1+1
          LAM1=LAM1P1-1
          LAM1TW=LAM1*2
          LPP=L1-LAM1-1
          NLPPMX=2*LAM1+1
          DO 98 NLPP=1,NLPPMX
            LPP=LPP+1
            IF (LPP.LT.0)    GO TO 98
            LPPTW=LPP*2
            KB=KB+1
            DO 96 LAM2P1=1,LAM2MX
              LAM2=LAM2P1-1
              LAM2TW=LAM2*2
              LL=LPP-LAM2-2
              DO 94  NLL=1,LAM2P1
                LL=LL+2
                IF (LL.LT.0)   GO TO 94
	        IF (LL+1.GT.LLP1MX) GO TO 94               
                LLTW=LL*2
                KC=KC+1
                IF (KC.GT.KCC)  WRITE(8,840)
                KBC(KC)=KB
	          LLP1D(KC)=LL+1
                LAM2PC(KC)=LAM2P1
                LPP1C(KC)=LPP+1
                NINTC(KC)=IT*3+IS+1
                H1=LLTW+1
                H2=LAM2TW+1
                HAT=SQRT(H1)/H2
                ID=0
                T1=1.0
                IA=LLTW
                IB=LPPTW
                IC=LAM2TW
                DO 50 M2=1,LAM2P1
                  IE=(M2-1)*2
                  IFF=IE
                  CALL CLEB(IA,IB,IC,ID,IE,IFF,FACLOG,RAC) 
                  CLEBQ(KC,M2)=RAC*FAC1*HAT*T1
                  T1=2.0
   50           CONTINUE
CCCCC
                DO 86  LAMP1=1,LAMP1X
                  LAM=LAMP1-1
                  KAT=(LAMP1*LAM)/2+(LL+LAM-L1)/2+1
                  IF (LAM.GT.(LL+L1).OR.LAM.LT.IABS(LL-L1))    GO TO 86
                  IF (MOD(LL+LAM+KPARIT,2).NE.0)               GO TO 86
                  KA=KAT
                  IF (KA.GT.KAMAX) KAMAX=KA
                  LAMTW=LAM*2
                  KG=KG+1
                  COEF(KG)=0.0
                  KCD(KG)=KC
                  KAD(KG)=KA
                  NLSD(KG)=NLS
                  LAMD(KA)=LAM
                  IF (MOD(LAM1+LAM2+LAM,2).NE.0) GO TO 55
                  IF (LAM.GT.(LAM1+LAM2).OR.LAM.LT.IABS(LAM1-LAM2))
     1              GO TO 55
                  IA=LAM1TW
                  IB=LAM2TW
                  IC=LAMTW
                  CALL CLEBZ(IA,IB,IC,FACLOG,RAC) 
                  C1=RAC
                  H1=(LAM1TW+1)*(LAM2TW+1)
                  C3=SQRT(H1)
                  S1=1-MOD(LL,2)*2
                  IA=LLTW
                  IB=L1TW
                  IC=LAM2TW
                  ID=LAM1TW
                  IE=LAMTW
                  IFF=LPPTW
                  CALL RAC7(IA,IB,IC,ID,IE,IFF,FACLOG,RAC) 
                  C5=RAC
                  A1=1.0
                  COEF(KG)=COEF(KG)+C1*C3*C5*S1*FPISQI*A1
   55             CONTINUE
                  ID=0
                  IA=LLTW
                  IB=LAMTW
                  IC=L1TW
                  H1=LLTW+1
                  F1=FAC1*SQRT(H1)/FLOAT(L1TW+1)
                  T1=1.0
                  M1X=MIN0(L1MAX+1,LAMP1X)
                  DO 60 M1=1,M1X
                    IE=(M1-1)*2
                    IFF=IE
                    CALL CLEB(IA,IB,IC,ID,IE,IFF,FACLOG,RAC) 
                    T2=1-2*MOD(M1-1,2)
                    CLEBF(KA,M1,NLS)=RAC*F1*T1*T2
                    T1=2.0
   60             CONTINUE
CCCCCC
                  IF (KTLOUT(4).EQ.1) WRITE(8,841) KG,KC,KA,KBC(KC),
     1            NLS,LAM1,LPP,LAM2,LL,K,LT,LAM,COEF(KG),KA,LLP1D(KA),
     1            (CLEBF(KA,M1,NLS),M1=1,2),(CLEBQ(KC,M2),M2=1,LAM2MX)
CCCCCC
   86           CONTINUE
   94         CONTINUE
   96       CONTINUE
   98     CONTINUE
  100   CONTINUE
  120 CONTINUE
      KCMAX=KC
      KGMAX=KG
CCCCCC 
      IF (KTLOUT(3).NE.0) THEN
      WRITE(8,851) KGMAX,KCMAX,KAMAX
      ENDIF
      IF(KGMAX.LE.KGG.AND.KCMAX.LE.KCC.AND.KAMAX.LE.KAA)
     1                GO TO 130
      WRITE(8,852)
      STOP
  130 CONTINUE
  840 FORMAT(1H  ,'KC  IS OVER KCC')
  841 FORMAT(1H  ,12I4,F8.3,I6,I4,2F8.3,3X,5F8.3)
  851 FORMAT(1H  ,5X,'(GFAC)KGMAX,KCMAX,KAMAX=',3I3)
  852 FORMAT(1H ,'ONE OF KGMAX ETC IS OVER DIM. CALCULATION STOPS')
CCCCCC
      JLS=0
      NLTM1MX=NOLTR*MP1MAX
      DO 300 NLTM1=1,NLTM1MX
        NLT=(NLTM1-1)/MP1MAX+1
        M1P1=NLTM1-(NLT-1)*MP1MAX
        M1=M1P1-1
        LT=LTR(NLT)
        LTTW=LT*2
        M1TW=M1*2
        DO 244 NLS=1,NLSMAX
          L1=L1R(NLS)
          L1TW=L1*2
          DO 242 KA=1,KAMAX
            LAM=LAMD(KA)
            LL=KA*2-(LAM+2)*LAM+L1-2
            LAMTW=LAM*2
            LLTW=LL*2
            KT=(KCT-1)*2
            KTTW=KT*2
            DO 230 LRP1=1,LRP1MX
              LR=LRP1-1
              IF (LR.GT.KT+LAM)       GO TO 230
              IF (LR.LT.IABS(KT-LAM)) GO TO 230
              IF (MOD(KT+LAM+LR,2).NE.0) GO TO 230
              JLS=JLS+1
              LRTW=LR*2
              H1SQ=(KTTW+1)*(LAMTW+1)*(LRTW+1)*(L1TW+1)
              H1=SQRT(H1SQ)
              S1=(-1.0)**(KT+L1-LT)
              S2=(-1.0)**((LR+LL-KPARIT)/2)
              IA=KTTW
              IB=LAMTW
              IC=LRTW
              CALL CLEBZ(IA,IB,IC,FACLOG,RAC) 
              C1=H1*S1*RAC*S2
              IA=LLTW
              IB=LRTW
              IC=LTTW
              ID=M1TW
              IE=0
              IFF=M1TW
              CALL CLEB(IA,IB,IC,ID,IE,IFF,FACLOG,RAC) 
              C2=RAC
              IA=LLTW
              IB=LAMTW
              IC=LTTW
              ID=KTTW
              IE=L1TW
              IFF=LRTW
              CALL RAC7(IA,IB,IC,ID,IE,IFF,FACLOG,RAC) 
              C3=RAC
              CCDD(JLS)=C1*C2*C3
  230       CONTINUE
  242     CONTINUE
  244   CONTINUE
  300 CONTINUE
      JLSMX=JLS
CCCCCC
      IF (KTLOUT(3).NE.0) THEN
      WRITE(8,862) JLSMX,KCT
  862 FORMAT(1H ,5X,'(GFAC)JLSMX,KCT=',2(I8,2X),/)
      ENDIF
CCCCCC
      RETURN
      END
