      SUBROUTINE CROSS(FAC)
CCCCCC
      PARAMETER(NHS=9,NPS=9,LXA=150,LXB=150,LPH=10,LMM=4)
      PARAMETER(NXA=420,NXB=140,NIN=300)
      PARAMETER(KFM=10)
CCCCCC	
      COMMON/CNST  /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
      COMMON/CNTRL/ KTRLD(9),KTLOUT(24),NERUN,KRTYPE
      COMMON/DWCC/  LDWMIR(4),LDWMXR(4),WN(4),LDWSTR(4),CE(4),CFUNIR(4),
     1              ECM(4),XMESR(4),NXMXR(4),NXMIR(4),NOEXP,MXSTEP,
     2              TMASA,TMASB,TMASX,TMASY,PMASA,PMASB,PMASX,PMASY,
     3              TZA,TZB,TZX,TZY,
     4              PZA,PZB,PZX,PZY,XBARD(4),XBARID(4),WNXD(20),EXPD(20)
      COMMON/ANGLC/ NTHEB,NDTHE,THEB,DTHEB
      COMMON/SPSTAT/EET,SPECTA(LPH),ESP(LPH),ESH(LPH),NSP(LPH),LSP(LPH),
     1              JTWP(LPH),NSH(LPH),LSH(LPH),JTWH(LPH),IPAIR,NOSP,
     2              NOSH,NPST(LPH),NHST(LPH),NPAIRD,L1R(8),ISR(8),
     3              ITR(8),NLSMAX,JT,KPARIT,IST,NOLTR,N1,LTR(8),
     4              ITP(LPH),ITH(LPH),ITZP(LPH),ITZH(LPH)
      COMMON/RELKIN/RMASA,RMASB,TREL1,TREL2,OMEGR(20)
      COMMON/FFCC/  VRANG(6),VSTR(12,6),VV(NIN,24),XMESH,XMESC,
     1              XMESHD,XMESHE,RHED(NIN),
     2              NBCMI,NBCMX,NONB,NBSTEP,LASTEP,LBSTEP,MAMAX,MXMAX,
     3              MP1MAX,NOLA,NOLB,INTRAN,LMI,LMX,KCETN(2),
     4              NCMAX,NCMIN,LDMAX,JATRTY,NONA,NONAH,NONAR,NASTEP,
     5              NHDMX,NHEMX,NGAUSR,NBSTPD,NBSTPE,N1STEP,JATW,ISATW
      COMMON/CCKF/  L1KFD(KFM),ISKFD(KFM),KKFD(KFM),LTKFD(KFM),KFMAX,
     1              ITKFD(KFM),NLSKFD(KFM),NLTKFD(KFM),LTMIM,LTMAX
      COMMON/OVDE/  OVDD(20,20),OVED(20,20),OVDDP(10,10),OVDEP(10,10)   
      COMPLEX       OVED,OVDD,OVDDP,OVDEP                               
      COMMON/CTFAC/ TFAC(3,10,10),MJMAX                                 
      DIMENSION     CROSD(20,4)
      DIMENSION     LTOUT(89)
      COMPLEX       OV1,OV2,OV3,TTI,ZERO
CCCCC
CCCCC     NOTE THAT THIS PROGRAM ALLOWS ONLY ONE SET OF IS1
CCCCC
      ZERO=(0.0,0.0)
      TTI=(0.0,1.0)
      KOUT=KTLOUT(6)
      HBRSQ=HBAR*HBAR
      CONV=180.0/PI
      LAP1MI=LDWMIR(1)+1
      LAP1MX=LDWMXR(1)+1
      LASTEP=LDWSTR(1)
      IF(LASTEP.EQ.0) LASTEP=1
      LBP1MI=LDWMIR(2)+1
      LBP1MX=LDWMXR(2)+1
      LBSTEP=LDWSTR(2)
      IF(LBSTEP.EQ.0) LBSTEP=1
      NOLB=(LBP1MX-LBP1MI)/LBSTEP+1
      NOLA=(LAP1MX-LAP1MI)/LASTEP+1
      LAP1MX=LAP1MI+(NOLA-1)*LASTEP
      MXP1MX=MIN0(LTMAX,MXMAX)+1
      JTTW=JT*2
      ISTTW=IST*2
      KFMAX=NOLTR
CCCCC
CCCCC           TRANSITION DENSITIES
CCCCC
      IF(KTRLD(9).EQ.1) THEN
      REL1=(4.0*PI/(WN(1)*WN(2)))**2
      T1=10.*(WN(2)/WN(1))/(TREL1*TREL2*REL1)
      ELSE
      T1 =((AMAS/(2.0*PI*HBRSQ))**2)*RMASA*RMASA*10.0*WN(2)/WN(1)
      ENDIF
CCCCCC
      FLASTP=FLOAT(LASTEP*LASTEP)
      T2=1.0/FLOAT((JATW+1)*(ISATW+1))
      T3=FLOAT(JT*2+1)
      FAC=T1*FLASTP*T2*T3*0.8400
	TCROST=0.0
      TCROSD=0.0     
      TCROSE=0.0 
      DO 910 LTM=1,NOLTR
         CROSD(LTM,1)=0.0
         CROSD(LTM,2)=0.0
         CROSD(LTM,3)=0.0
 910  CONTINUE
      DO 915 LTM=1,NOLTR
         LTP1=LTR(LTM)+1
         LT1=LTP1-1
         LTOUT(LTM)=LTP1-1
         SCROS=0.0
         SCROSD=0.0
         SCROSE=0.0
         F1=1.0
         MP1XT=MIN0(LTP1,MP1MAX)
         DO 912 MP1=1,MP1XT
            OV3=OVED(LTM,MP1)
            OV2=OVDD(LTM,MP1)
            OV1=OV2+OV3
            OV1R=REAL(OV1)
            OV1I=AIMAG(OV1)
            OV2R=REAL(OV2)
            OV2I=AIMAG(OV2)
            OV3R=REAL(OV3)
            OV3I=AIMAG(OV3)
            T1=(OV1R**2+OV1I**2)*F1
            T2=(OV2R**2+OV2I**2)*F1
            T3=(OV3R**2+OV3I**2)*F1
            SCROS =SCROS +T1
            SCROSD=SCROSD+T2
            SCROSE=SCROSE+T3
            F1=2.0
  912    CONTINUE
         CROSD(LTM,1)=SCROS *FAC
         CROSD(LTM,2)=SCROSD*FAC
         CROSD(LTM,3)=SCROSE*FAC
         TCROST = TCROST + CROSD(LTM,1)
	 TCROSD = TCROSD + CROSD(LTM,2)
	 TCROSE = TCROSE + CROSD(LTM,3)
  915 CONTINUE
      WRITE(8,832) THEB
      WRITE(8,833)
      DO 920 LTM=1,NOLTR
      LT=LTR(LTM)
      WRITE(8,834) LT,CROSD(LTM,2),CROSD(LTM,3),
     1         CROSD(LTM,1)
  920 CONTINUE
      WRITE(8,835) TCROSD,TCROSE,TCROST
  832 FORMAT(///1H  ,'ANGLE=',F8.4,/)
  833 FORMAT(15X,'DIRECT',6X,'EXCHANGE',4X,'TOTAL',/)
  834 FORMAT(1H ,'LT=',I3,5X,3(2X,E8.3,2X))
  835 FORMAT(1H ,'TOTAL CS=',2X,3(2X,E8.3,2X))
CCCCCC
CCCCCC       CROSS SECTIONS FROM SOURCE FUNCTIONS
CCCCCC
      IF (KTRLD(6).NE.0) THEN
      WRITE(8,843)
  843 FORMAT(/,1H ,'CROSS SECTIONS FROM SOURCE FUNCTIONS',/)
      CROSD1=0.0                                                        
      CROSD2=0.0                                                      
      CROSD3=0.0                                                      
      MJP1MX=MJMAX+1                                                  
      DO 940 LXYZ=1,3                                                 
         SUM1=0.0                                                        
         SUM2=0.0                                                        
         SUM3=0.0                                                        
         DO 930 MJP1=1,MJP1MX                                            
            F1=2.0                                                          
            IF(MJP1.EQ.1) F1=1.0                                            
            OV2=OVDDP(LXYZ,MJP1)                                            
            OV3=OVDEP(LXYZ,MJP1)                                            
            OV1=OV2+OV3                                                     
            OV1R=REAL (OV1)                                                 
            OV1I=AIMAG(OV1)                                                 
            OV2R=REAL (OV2)                                                 
            OV2I=AIMAG(OV2)                                                 
            OV3R=REAL (OV3)                                                 
            OV3I=AIMAG(OV3)                                                 
            SUM1=SUM1+(OV1R**2+OV1I**2)*F1*FAC                              
            SUM2=SUM2+(OV2R**2+OV2I**2)*F1*FAC                              
            SUM3=SUM3+(OV3R**2+OV3I**2)*F1*FAC                              
  930    CONTINUE                                                        
         CROSD1=CROSD1+SUM1                                              
         CROSD2=CROSD2+SUM2                                              
         CROSD3=CROSD3+SUM3                                              
         WRITE(8,844) LXYZ,CROSD2,CROSD3,CROSD1           
  940 CONTINUE 
  844 FORMAT(1H ,'LXYZ=',I3,3X,3(2X,E8.3,2X))                         
CCCCCC
CCCCCC
      T0RD=REAL(OVDD(1,1))
      T0ID=AIMAG(OVDD(1,1))
      T2RD=REAL(OVDD(2,1))
      T2ID=AIMAG(OVDD(2,1))
      T0RE=REAL(OVED(1,1))
      T0IE=AIMAG(OVED(1,1))
      T2RE=REAL(OVED(2,1))
      T2IE=AIMAG(OVED(2,1))
      T0R =T0RD+T0RE
      T0I =T0ID+T0IE
      T2R =T2RD+T2RE
      T2I =T2ID+T2IE
      T0R2=T0R**2
      T0I2=T0I**2
      DOWN=T0R2+T0I2
      T02R=T0R*T2R
      T02I=T0I*T2I
      UPPER=T02R+T02I
      TERM2=2.*SQRT(2.)/3.*UPPER/DOWN
      DNN  =(-1./3.)+TERM2
      DNNR =-SUM1/CROSD1
      WRITE(8,851) T0RD,T0ID
      WRITE(8,852) T0RE,T0IE
      WRITE(8,853) T0R,T0I
      WRITE(8,854) T2RD,T2ID
      WRITE(8,855) T2RE,T2IE
      WRITE(8,856) T2R,T2I
      WRITE(8,857) T0R2,T0I2,DOWN
      WRITE(8,858) T02R,T02I,UPPER
      WRITE(8,859) TERM2
      WRITE(8,860) DNN
      WRITE(8,861) DNNR
  851 FORMAT(//1H ,'T00(D)         =',3F15.3)
  852 FORMAT(1H ,'T00(E)         =',3F15.3)
  853 FORMAT(1H ,'T00(Total)     =',3F15.3)
  854 FORMAT(/1H ,'T22(D)         =',3F15.3)
  855 FORMAT(1H ,'T22(E)         =',3F15.3)
  856 FORMAT(1H ,'T22(Total)     =',3F15.3)
  857 FORMAT(/1H ,'T0R2,T0I2,down =',3F15.3)
  858 FORMAT(1H ,'T02R,T02I,upper=',3F15.3)
  859 FORMAT(1H ,'Second Term    =',3F15.3)
  860 FORMAT(//1H ,'DNN values     =',3F15.3)
  861 FORMAT(1H ,'DNN real       =',3F15.3)
      ENDIF
CCCCCC
      RETURN
      END
