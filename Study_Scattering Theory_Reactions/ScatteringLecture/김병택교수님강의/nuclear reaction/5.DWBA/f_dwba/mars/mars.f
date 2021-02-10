C ABPBMARS-1-FOR-EFR-DWBA. EXACT FINITE RANGE DWBA CALCULATIONS FOR             
C 1   HEAVY-ION INDUCED NUCLEAR REACTION.  TAMURA,T., LOW,K.S.                  
C REF. IN COMP. PHYS. COMMUN. 8 (1974) 349                                      
C      PROGRAM MARS (INPUT,OUTPUT,TAPE5=INPUT,TAPE6=OUTPUT,TAPE11,TAPE12)        
      PROGRAM MARS         
      COMMON/CRAC/FACLOG(500),IA,IB,IC,ID,IE,IF,RAC,L9(10),U9                   
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/OUTP/ NOUT(20),IBSPR(10),BSPAR(20),OMPAR(20)                       
      COMMON COULST(1600),POTENT(3000),FORM(2500)                               
CCCCCC    /CRAC/   APPEARS IN MARS,JP1MAR,OVLAP,DWAMP.                          
CCCCCC    /DWBA/   APPEARS IN MARS,CBCTRL,INEXCM,OVLAP,INTGMS,INTPOL,           
CCCCCC                        URENOM,ELCROS,DWAMP,DWCSOS.                       
CCCCCC    /OUTP/   APPEARS IN MARS,JP1MAR,DWCROS.                               
CCCCCC    /JP1/    APPEARS IN JP1MAR,POTEMS,CBCTRL,INEXCM,OVLAP,INTGMS,         
CCCCCC                        INTPOL,URENOM,ELCROS,LBSELC,DWAMP,DWCROS.         
CCCCCC    /COU/    APPEARS IN JP1MAR,FLGLCH,OVLAP.                              
CCCCCC    /OVLC/   APPEARS IN OVLAP,INTGMS,INTPOL,URENOM.                       
CCCCCC    /ELAS/   APPEARS IN OVLAP,ELCROS.                                     
CCCCCC    /LEGENC/ APPEARS IN LEGNDR,DWCROS.                                    
      DIMENSION NAT(100),            KPARWR(2),KEVOD(2)                         
      DATA KPARWR,KEVOD/3H(+),3H(-),2H/2,1H /                                   
      DATA NEUTRN,NAT /2H N,                                                    
     *        2H H,2HHE,2HLI,2HBE,2H B,2H C,2H N,2H O,2H F,2HNE,2HNA,           
     *        2HMG,2HAL,2HSI,2H P,2H S,2HCL,2HAR,2H K,2HCA,2HSC,2HTI,           
     *        2H V,2HCR,2HMN,2HFE,2HCO,2HNI,2HCU,2HZN,2HGA,2HGE,2HAS,           
     *        2HSE,2HBR,2HKR,2HRB,2HSR,2H Y,2HZR,2HNB,2HMO,2HTC,2HRU,           
     *        2HRH,2HPD,2HAG,2HCD,2HIN,2HSN,2HSB,2HTE,2H I,2HXE,2HCS,           
     *        2HBA,2HLA,2HCE,2HPR,2HND,2HPR,2HSM,2HEU,2HGD,2HTB,2HDY,           
     *        2HHO,2HER,2HTM,2HYB,2HLU,2HHF,2HTA,2H W,2HRE,2HOS,2HIR,           
     *        2HPT,2HAU,2HHG,2HTL,2HPB,2HBI,2HPO,2HAT,2HRN,2HFR,2HRA,           
     *        2HAC,2HTH,2HPA,2H U,2HNP,2HPU,2HAM,2HCM,2HBK,2HCF,2HES,           
     *        2HFM/                                                             
CCCCCC      IF KTRLCB(5)=1, ONLY D.W. IN INCIDENT CHANNEL IS CALCULATED.        
CCCCCC      IF KTRLCB(7)=1 NO-RECOIL FINITE-RANGE CALCULATION IS MADE.          
CCCCCC      IF KTRLCB(8)=1 EXACT FINITE-RANGE CALCULATION (WITH RECOIL)         
CCCCCC                     IS MADE.                                             
CCCCCC      KEXTCB(1) TIMES MESH SIZE GIVES THE LOWER CUT-OFF RADIUS IN         
CCCCCC                THE OVERLAP INTEGRAL.                                     
CCCCCC      KEXTCB(2) IS, IF NON-ZERO, THE NUMBER OF MESH POINTS AT             
CCCCCC                WHICH THE FORM FACTOR IS READ IN.                         
CCCCCC      KEXTCB(3) IS, IF NON-ZERO, THE LOWER CUT-OFF LA AS DEFINED          
CCCCCC                BY KEXCST(1) IN SATURN.                                   
CCCCCC      KEXTCB(4) IS, IF NON-ZERO, EQUAL TO 2 AS DEFINED BY                 
CCCCCC                KEXCST(2) IN SATURN, IN WHICH CASE EVERY OTHER            
CCCCCC                LA IS CALCULATED ONLY.                                    
CCCCCC      KOUTCB(3) IS USED TO GENERATE OUTPUT IN DWAMP AND DWCROS.           
CCCCCC
      OPEN(UNIT=5, FILE='MARS.DAT',STATUS='OLD')
      OPEN(UNIT=6, FILE='MARS.OUT',STATUS='UNKNOWN')
CCCCCC
      REWIND 12
	REWIND 11                                      
      FACLOG(1)=0.0                                                             
      FACLOG(2)=0.0                                                             
      FN=1.0                                                                    
      DO 10 N=3,500                                                             
      FN=FN+1.0                                                                 
   10 FACLOG(N)=FACLOG(N-1)+ALOG(FN)                                            
  700 FORMAT(24I3)                                                              
  701 FORMAT(8I5,3F10.5)                                                        
  702 FORMAT(4I5,2F10.5)                                                        
  704 FORMAT(7F10.5)                                                            
  706 FORMAT(6I5,/,9F7.2,/,9F7.2)                                               
  100 DO 103 N=1,10                                                             
      KTRLCB(N)=0                                                               
      KOUTCB(N)=0                                                               
  103 KEXTCB(N)=0                                                               
CCCCCC  ******    INPUT CARDS IN MARS AND READING OF INFORMATION  ******        
CCCCCC  ******    FROM SATURN VIA I/O UNIT 12                     ******        
      READ(5,701) KTRLCB(5),KTRLCB(7),KTRLCB(8),KOUTCB(3),KTAP12,KRALSJ,        
     1          NTHETA,NONE  ,TETADG(1),DTETAD,ELABI                            
      IF(KTRLCB(5).EQ.1) GO TO 610                                              
      ALSJF(1)=1.0                                                              
      IF(KRALSJ.EQ.0) GO TO 107                                                 
      KEXTCB(5)=KRALSJ                                                          
      READ(5,704) (ALSJF(N),N=1,KRALSJ)                                         
  107 IF(KTRLCB(7).EQ.1) GO TO 111                                              
      READ(12)   MASCA,NZCA,KPCA,ICATW,MASCB,NZCB,KPCB,ICBTW,                   
     1           MASSA,NZSA,KPSA,ISATW,MASSB,NZSB,KPSB,ISBTW,                   
     2           JLSMAX,LDWMAX,NBREAD,LBSTTL,LACUT,LASTP                        
      READ(12)  (JJTWR(N),LLTWR(N),ISTWR (N),N=1,JLSMAX)                        
      READ(12)  NATOTL,NBMIN,NBMAX,MESFCA,MESFCB,XMESA,XMESB,QVALGR             
      READ(12) (IBSPR(N),N=1,6),(BSPAR(N),N=1,18)                               
      NBTOTL=NBMAX-NBMIN+1                                                      
      MSTOT=(NBTOTL-3)*MESFCB                                                   
      KEXTCB(1)=(NBMIN+1)*MESFCB                                                
      KEXTCB(3)=LACUT                                                           
      KEXTCB(4)=LASTP                                                           
      GO TO 171                                                                 
  111 IF(KTAP12.EQ.0) GO TO 114                                                 
      READ(12)   MASCA,NZCA,KPCA,ICATW,MASCB,NZCB,KPCB,ICBTW,                   
     1           MASSA,NZSA,KPSA,ISATW,MASSB,NZSB,KPSB,ISBTW,                   
     2           JLSMAX,LDWMAX,NBREAD,LBSTTL,LACUT,LASTP                        
      READ(12)  (JJTWR(N),LLTWR(N),ISTWR (N),N=1,JLSMAX)                        
      READ(12) KEXTCB(1),KEXTCB(2),MSTOT,MESFCA,XMESA,QVALGR                    
      READ(12) (IBSPR(N),N=1,6),(BSPAR(N),N=1,18)                               
      GO TO 121                                                                 
  114 READ(5,700)MASCA,NZCA,KPCA,ICATW,MASCB,NZCB,KPCB,ICBTW,                   
     1           MASSA,NZSA,KPSA,ISATW,MASSB,NZSB,KPSB,ISBTW,                   
     2           JLSMAX,LDWMAX                                                  
      READ(5,700)(JJTWR(N),LLTWR(N),ISTWR(N),N=1,JLSMAX)                        
      READ(5,702)KEXTCB(1),KEXTCB(2),MSTOT,MESFCA,XMESA,QVALGR                  
      READ(5,706)(IBSPR(N),N=1,6),(BSPAR(N),N=1,18)                             
CCCCCC  ******             READ IN NO-RECOIL FORM FACTOR          ******        
  121 DO 160 N1=1,JLSMAX                                                        
      NXFBS=(N1-1)*MSTOT-KEXTCB(1)                                              
      KEXCB1=KEXTCB(1)                                                          
      KEXCB2=KEXTCB(2)                                                          
      IF(KTAP12   .EQ.1 ) READ(12) KEXCB2,(COULST(NX),NX=1,KEXCB2)              
      IF(KTAP12   .EQ.1 ) GO TO 154                                             
      READ(5,153) (COULST(NX),NX=1,KEXCB2)                                      
  153 FORMAT(5E14.6)                                                            
  154 NXMIN=KEXCB1+1                                                            
      NXMAX=KEXCB1+MSTOT                                                        
      DO 155 NX=NXMIN,NXMAX                                                     
      NXF=NXFBS+NX                                                              
  155 FORM(NXF)=COULST(NX)                                                      
  160 CONTINUE                                                                  
  171 K1=IABS(KPCA+KPSA+KPCB+KPSB)                                              
      K1=MOD(K1,2)                                                              
      LLMAX=0                                                                   
      DO 177 JLS=1,JLSMAX                                                       
      LL=LLTWR(JLS)/2                                                           
      IF(MOD(LL,2).NE.K1) LL=LL-1                                               
      LLMAX=MAX0(LLMAX,LL)                                                      
  177 CONTINUE                                                                  
      LMXPAR=LLMAX                                                              
      DO 220 NA=2,NTHETA                                                        
  220 TETADG(NA)=TETADG(NA-1)+DTETAD                                            
      DO 222 NA=1,NTHETA                                                        
  222 TETARD(NA)=TETADG(NA)*0.0174532925                                        
CCCCCC  ******             PREPARATION FOR INITIAL OUTPUT         ******        
      NATWI=NAT(NZCA)                                                           
      NAPWI=NAT(NZSA)                                                           
      NATWE=NAT(NZCB)                                                           
      NAPWE=NAT(NZSB)                                                           
      IF(NZSA.EQ.0) NAPWI=NEUTRN                                                
      IF(NZSB.EQ.0) NAPWE=NEUTRN                                                
      KEOI=2-MOD(ICATW,2)                                                       
      IIWR=ICATW/KEOI                                                           
      KEOE=2-MOD(ICBTW,2)                                                       
      IEWR=ICBTW/KEOE                                                           
      KEOIS=2-MOD(ISATW,2)                                                      
      ISIWR=ISATW/KEOIS                                                         
      KEOES=2-MOD(ISBTW,2)                                                      
      ISEWR=ISBTW/KEOES                                                         
      KPARIS=KPSA+1                                                             
      KPARES=KPSB+1                                                             
      KPARI =KPCA+1                                                             
      KPARE =KPCB+1                                                             
CCCCCC  ******                      INITIAL OUTPUT                ******        
      WRITE(6,404)                                                              
  404 FORMAT(1H1,49X,20H********************)                                   
      IF(KTRLCB(8).EQ.1) WRITE(6,405)                                           
      IF(KTRLCB(7).EQ.1) WRITE(6,406)                                           
  405 FORMAT(//35X,50HEXACT FINITE RANGE DISTORTED WAVE BORN CALCULATION        
     1)                                                                         
  406 FORMAT(//33X,54HNO RECOIL FINITE RANGE DISTORTED WAVE BORN CALCULA        
     1TION)                                                                     
      WRITE(6,410)NATWI,MASCA,NAPWI,MASSA,NAPWE,MASSB,NATWE,MASCB,ELABI         
  410 FORMAT(//27X,4HFOR ,                                                      
     1  A2,1H(,I3,2H)(A2,1H(I3,2H),A2,1H(I3,2H)),A2,1H(,I3,10H)-REACTION        
     2   ,9H AT ELAB=F7.3,4H-MEV)                                               
      WRITE (6,411) IIWR,KEVOD(KEOI),KPARWR(KPARI),ISIWR,KEVOD(KEOIS),          
     *   KPARWR(KPARIS),ISEWR,KEVOD(KEOES),KPARWR(KPARES),IEWR,                 
     *   KEVOD(KEOE),KPARWR(KPARE),QVALGR                                       
  411 FORMAT (31X,I2,A2,A3,1H(,I2,A2,A3,1H,,I2,A2,A3,1H),I2,A2,A3,              
     1   10X,8HQ VALUE=,F7.3,4H MEV)                                            
      WRITE (11)   NATWI,MASCA ,NAPWI,MASSA ,NAPWE,MASSB ,NATWE,MASCB ,         
     1              IIWR,KEVOD(KEOI),KPARWR(KPARI),ISIWR,KEVOD(KEOIS),          
     2   KPARWR(KPARIS),ISEWR,KEVOD(KEOES),KPARWR(KPARES),IEWR,                 
     3   KEVOD(KEOE),KPARWR(KPARE)                                              
      REWIND 11                                                                 
      READ(11) (NOUT(N),N=1,20)                                                 
      WRITE(6,412)                                                              
  412 FORMAT(///,50X,20H********************///1H )                             
      WRITE(6,461) (JJTWR(N),LLTWR(N),ISTWR(N),N=1,JLSMAX)                      
  461 FORMAT(25X,29HFORM FACTORS HAVE (2J,2L,2S)=, 5(1H(,3I2,1H),1X))           
      WRITE(6,483) MSTOT                                                        
  483 FORMAT(25X,18HTOTAL MESH POINTS=,I3)                                      
      WRITE(6,491) (KTRLCB(N),N=1,10)                                           
  491 FORMAT(/25X,07HKTRLCB=,10I5)                                              
      WRITE(6,493) (KOUTCB(N),N=1,10)                                           
  493 FORMAT( 25X,07HKOUTCB=,10I5)                                              
      WRITE(6,495) (KEXTCB(N),N=1,10)                                           
  495 FORMAT( 25X,07HKEXTCB=,10I5)                                              
      IF(KTRLCB(5).EQ.1) GO TO 610                                              
      WRITE(6,496)                                                              
  496 FORMAT(//20X,5(1H*),22HBOUND STATE PARAMETERS,/8X,7HSYSTEM,3X,            
     1   23H(N,L,J)  BINDING ENERGY,7X,3HVSX,5X,4HVSOR,4X,4HDFNR,4X,            
     2   5HDFNSO,3X,3HRZR,5X,4HRZSO,4X,3HRZC,6X,4HXMES)                         
      I=1                                                                       
      WRITE(6,497)I,(IBSPR(N),N=1,3),(BSPAR(N),N=1,9)                           
      I=2                                                                       
      WRITE(6,497)I,(IBSPR(N),N=4,6),(BSPAR(N),N=10,18)                         
  497 FORMAT(9X,I2,5X,1H(,3I2,1H),4X,F8.3,7X,8F8.3)                             
      WRITE (6,502)                                                             
  502 FORMAT(//20X,5(1H*),22HFORM FACTOR PARAMETERS )                           
      IF(KTRLCB(8).EQ.0) GO TO 507                                              
      WRITE(6,504)NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,NBREAD                     
  504 FORMAT(/,8X,51HEXACT FINITE RANGE FORM FACTORS READ IN FROM TAPE12        
     1       /8X,42HNATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,NBREAD= ,6I5)            
      IF(LACUT+LASTP.EQ.0) GO TO 517                                            
      IF(LASTP.EQ.0) LASTP=1                                                    
      WRITE(6,505)LACUT,LASTP                                                   
  505 FORMAT(///25X,9HIN SATURN,//,8X,25HFORM FACTOR FOR THE FIRST,I5,          
     1 35H PARTIAL WAVES WERE NOT CALCULATED.,/8X,45HAND THE FORM FACTOR        
     2 IS CALCULATED AT STEPS OF,I3)                                            
      GO TO 517                                                                 
  507 R1=(KEXTCB(1)+1)*XMESA                                                    
      WRITE(6,508)R1,XMESA                                                      
  508 FORMAT(/8X,27HNO RECOIL FORM FACTOR FROM ,F8.3,9H FERMI AT,F8.3,          
     1    11H FERMI MESH)                                                       
      DO 515 N1=1,JLSMAX                                                        
      WRITE(6,514) N1                                                           
  514 FORMAT(1H ,06H  NFF=,I2)                                                  
      NXFBS=(N1-1)*MSTOT                                                        
      NXFMI=NXFBS+1                                                             
      NXFMX=NXFBS+MSTOT                                                         
      WRITE(6,516) (FORM(NXF),NXF=NXFMI,NXFMX)                                  
  515 CONTINUE                                                                  
  516 FORMAT(10E11.3)                                                           
CCCCCC  ******    CALL JP1MAR TO PREPARE FOR CALCULATING          ******        
CCCCCC  ******    DISTORTED WAVE IN THE INCIDENT  CHANNEL         ******        
  517 IF(NTHETA.GT.100) WRITE(6,521)NTHETA                                      
  521 FORMAT(///9H  NTHETA=,I5,7H GT 100)                                       
      IF(NTHETA.GT.100) STOP                                                    
      AOVBT=FLOAT(MASCA)/FLOAT(MASCB)                                           
      IF(MASSA.LT.MASSB)   AOVBT=1.0/AOVBT                                      
      IF(KTRLCB(8).EQ.1) AOVBT=1.0                                              
  610 KINEX=1                                                                   
      AOVERB=1.0                                                                
      IF(MASSA.LT.MASSB)   AOVERB=AOVBT                                         
      IF(KTRLCB(5).EQ.1) AOVERB=1.0                                             
      CALL JP1MAR                                                               
      IF(KTRLCB(5).EQ.1) GO TO 625                                              
CCCCCC  ******    CALL JP1MAR TO PREPARE FOR CALCULATING          ******        
CCCCCC  ******    DISTORTED WAVE IN THE EXIT CHANNEL              ******        
      KINEX=2                                                                   
      AOVERB=1.0                                                                
      IF(MASSA.GT.MASSB)   AOVERB=AOVBT                                         
      CALL JP1MAR                                                               
CCCCCC  ******    BEGINS TO CALCULATE THE DISTORTED WAVES         ******        
CCCCCC  ******    AND THE OVERLAP INTEGRALS                       ******        
  625 REWIND 11                                                                 
      CALL CBCTRL                                                               
CCCCCC  ******    BEGINS TO CALCULATE DWBA CROSS SECTIONS         ******        
      REWIND 11                                                                 
      CALL ELCROS                                                               
      IF(KTRLCB(5).EQ.1) GO TO 610                                              
      CALL DWAMP                                                                
      CALL DWCROS                                                               
      STOP                                                                    
      END                                                                       
      SUBROUTINE JP1MAR                                                         
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      COMMON/OUTP/ NOUT(20),IBSPR(10),BSPAR(20),OMPAR(20)                       
      COMMON/COU/ F(201),G(201),FD(200),GD(200),ETA,SIGMAZ,RHOMX,RD,            
     1            LMAX,KTOUT7                                                   
      COMMON F1(200),G1(200),FD1(200),GD1(200),EXSGR1(200),EXSGR2(200)          
      COMPLEX EXSGR1,EXSGR2                                                     
      COMMON      VCENR1(500),VCENI1(500),VCOUL1(500),                          
     1            VCENTR(500),VCENTI(500),VCOULM(500)                           
      DIMENSION OMCOM(10)                                                       
      EQUIVALENCE (OMCOM(1),VSX)                                                
      INTEGER PLUS,MINUS,EVEN,ODD,AMU,PMU,AMPMWR                                
      DATA PLUS,MINUS,EVEN,ODD,AMU,PMU/3H(+),3H(-),1H ,2H/2,3HAMU,3HPMU/        
      TTR=(1.0  ,0.0  )                                                         
      TTI=(0.0  ,1.0  )                                                         
      ZERO=(0.0  ,0.0  )                                                        
CCCCCC      IF KTRL(4)=1 R=RZERO*(TMAS**(1/3)+PMASS**(1/3))                     
CCCCCC      IF KTRL(9)=1 CALCULATION IS MADE ONLY FOR LA=KEXCOM(13)+1           
CCCCCC      KEXCOM(1)=MESH NUMBER UP TO WHICH LOWER CUT OFF IS MADE IN          
CCCCCC                SOLVING DIFFERENTIAL EQUATION FOR D.W.                    
CCCCCC      KEXCOM(11)=LA. THIS IS USED IN CBCTRL,DWBAMP AND LBSELC.            
CCCCCC      IF KEXCOM(12)=2 LA INTERPOLATION IS MADE.                           
CCCCCC      IF KEXCOM(13)=N IS NONZERO, FIRST N VALUES OF LA ARE IGNORED        
CCCCCC      IF KEXCOM(14)=0, LAP1MX=LDWMAX. OTHERWISE LAP1MX=KEXCOM(14).        
CCCCCC      KEXCOM(15)=LMXPAR.   KEXCOM(10)=LMIPAR.                             
CCCCCC      KTLOUT(3) IS USED IN CBCTRL TO GENERATE OUTPUT.                     
CCCCCC      KTLOUT(7) IS USED IN POTEMS AND IN FLGLCH.                          
CCCCCC      KTLOUT(10) IS USED IN OVLAP, INTGMS, INTPOL AND URENOM.             
CCCCCC                                                                          
  701 FORMAT(14I5)                                                              
  702 FORMAT(10F7.2)                                                            
      NDFMES=4                                                                  
      IF(KINEX.EQ.2) GO TO 135                                                  
      DO 115 N=1,15                                                             
      KTRL(N)=0                                                                 
      KEXCOM(N)=0                                                               
      KTLOUT(N)=0                                                               
  115 CONTINUE                                                                  
      READ(5,701) KTRL(4),KTRL(9),KEXCOM(1),(KEXCOM(N),N=12,14),                  
     1          KTLOUT(3),KTLOUT(7),KTLOUT(10)                                  
  135 READ(5,702) VSX,WSX,WSF,DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC            
      NXMAX=KEXTCB(1)+MSTOT                                                     
      NXPOT=MIN0(480,NXMAX)                                                     
CCCCCC   ******    NOTE THAT THE 480 ABOVE CORRESPONDS TO THE 500   ****        
CCCCCC   ******    LOCATIONS IN THE DIMENSION OF THE POTENTIALS     ****        
      IF(KEXTCB(4).NE.0) KEXCOM(12)=2                                           
      KEXCOM(13)=MAX0(KEXTCB(3),KEXCOM(13))                                     
      KPTOTL=KPCA+KPCB+KPSA+KPSB                                                
      IF(KEXCOM(13).NE.0) KEXCOM(13)=KEXCOM(13)+MOD(KPTOTL+KEXCOM(13),2)        
      IF(KEXCOM(14).EQ.0) KEXCOM(14)=LDWMAX                                     
      IF(KEXCOM(12).NE.0) KEXCOM(14)=KEXCOM(14)-MOD(KEXCOM(14)+KPTOTL+1,        
     1   2)                                                                     
      KEXCOM(15)=LMXPAR                                                         
      LMIPAR=LLTWR(JLSMAX)/2                                                    
      KEXCOM(10)=LMIPAR+MOD(LMXPAR+LMIPAR,2)                                    
      XMES2=XMESA/MESFCA*AOVERB                                                 
      XMES1=0.125*XMES2                                                         
      IF(KINEX.EQ.2) GO TO 146                                                  
CCCCCC  ******          BASIC PARAMETERS IN INCIDENT CHANNEL      ******        
      DO 142 N=1,10                                                             
  142 OMPAR(N)=OMCOM(N)                                                         
      ELAB=ELABI                                                                
      PMAS=MASSA                                                                
      TMAS=MASCA                                                                
      CHARGE=NZSA*NZCA                                                          
      GO TO 503                                                                 
CCCCCC  ******          BASIC PARAMETERS IN EXIT CHANNEL          ******        
  146 DO 148 N=1,10                                                             
  148 OMPAR(N+10)=OMCOM(N)                                                      
      ELABCI=ELABI* FLOAT(MASCA )/ FLOAT(MASSA+MASCA)                           
      ELABCE=ELABCI+QVALGR                                                      
      ELAB=ELABCE* FLOAT(MASCB+MASSB)/FLOAT(MASCB)                              
      PMAS=MASSB                                                                
      TMAS=MASCB                                                                
      CHARGE=NZSB*NZCB                                                          
      GO TO 507                                                                 
CCCCCC  ******                      INITIAL OUTPUT                ******        
  503 WRITE(6,506) ELAB                                                         
      GO TO 509                                                                 
  506 FORMAT(1H1,40X,18H(INCIDENT CHANNEL),10H     ELAB=F7.3)                   
  507 WRITE(6,508) ELAB                                                         
  508 FORMAT(1H1,40X,14H(EXIT CHANNEL),10H     ELAB=F7.3)                       
  509 WRITE(6,513)VSX,WSX,WSF,DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC          
  513 FORMAT(/21X,11HVSX,WSX,WSF,15X,1H=,3F8.3/21X,13HDFN,DFNW,DFNS,13X,        
     11H=,3F8.3/21X,27HRZERO,RZEROW,RZEROS,RZEROC=,4F8.3)                       
      WRITE(6,562) (KTRL(N),N=1,15)                                             
      WRITE(6,563) (KEXCOM(N),N=1,15)                                           
      WRITE(6,564) (KTLOUT(N),N=1,15)                                           
  558 WRITE(6,566) XMES1,XMES2                                                  
  562 FORMAT(  /7H KTRL  ,28I4)                                                 
  563 FORMAT(7H KEXCOM,28I4)                                                    
  564 FORMAT(7H KTLOUT,28I4)                                                    
  566 FORMAT(7H XMES1=F7.5,3X,7H XMES2=F7.5)                                    
CCCCCC  ******    C.M. ENERGY, WAVE NUMBER AND COULOMB PARAMETER  ******        
      RMAS=PMAS*TMAS/(PMAS+TMAS)                                                
      ECM   =ELAB*RMAS/PMAS                                                     
      ETUNIT=0.15745400                                                         
      WNUNIT=0.21870660                                                         
      E1= ABS(ECM)                                                              
      CE=ETUNIT*CHARGE*SQRT(RMAS/E1)                                            
      WN=WNUNIT*SQRT(RMAS*E1)                                                   
      WNINI=WNUNIT*SQRT(RMAS*(ECM+VSX))                                         
      ETA=CE                                                                    
      IF(ETA.GE.10.) GO TO 633                                                  
      ETA2=ETA*ETA                                                              
      ETA2A=2.0  *ETA                                                           
      ETA6=ETA2+16.0                                                            
      SIGMA0=-(ETA/(12.  *ETA6))*(1.  +(ETA2-48.  )/(30.  *(ETA6**2))           
     1       +((ETA2-160.  )*ETA2+1280.  )/(105.  *(ETA6**4)))                  
     2       -ETA+(ETA/2.  )* ALOG(ETA6)+3.5  * ATAN(0.25  *ETA)                
     4       -( ATAN(ETA)+ ATAN(0.5  *ETA)+ ATAN(ETA/3.0  ))                    
      GO TO 636                                                                 
  633 EINV1=1.0  /ETA                                                           
      EINV2=EINV1*EINV1                                                         
      EINV3=EINV1*EINV2                                                         
      EINV5=EINV3*EINV2                                                         
      EINV7=EINV5*EINV2                                                         
      EINV9=EINV7*EINV2                                                         
      SIGMA0=0.7853981634  +ETA*ALOG(ETA)-ETA                                   
     1      -(0.08333333333  *EINV1+0.00277777777  *EINV3                       
     2       +0.00079365079  *EINV5+0.00059523810  *EINV7                       
     3       +0.00084175084  *EINV9)                                            
  636 MODTPI=SIGMA0/6.2831853072                                                
      SGMAZZ    =SIGMA0-6.2831853072*MODTPI                                     
CCCCCC  ******                 FIXING RADIAL PARAMETERS           ******        
      A1=TMAS**0.333333333333                                                   
      IF(KTRL(4).NE.0) A1=A1+PMAS**0.333333333333                               
      XBAR=RZERO*A1                                                             
      XMAX=NXMAX*XMES2                                                          
      WRITE(6,663) ECM,WN,WNINI,CE,SGMAZZ                                       
  663 FORMAT(24H ECM,WN,WNINI,ETA,SIGMO= ,5E15.5)                               
      WRITE(6,673) KEXCOM(14),NXMAX,NXPOT,XMAX,XBAR                             
  673 FORMAT(19H KEX14,NXMAX,NXPOT=,3I5,11H XMAX,XBAR=,2E15.7)                  
CCCCCC  ******    CALCULATE COULOMB WAVE FUNCTION AT MATCHING RADIUS  **        
      LMAX=KEXCOM(14)                                                           
      IF(KINEX.EQ.2) LMAX=KEXCOM(14)+KEXCOM(15)                                 
      LMAXM1=LMAX-1                                                             
      X=XMAX                                                                    
      SIGMAZ=SGMAZZ                                                             
      RHOMX=X*WN                                                                
      RD=XMES2*WN                                                               
      KTOUT7=KTLOUT(7)                                                          
      CALL FLGLCH                                                               
      ETASQ=ETA*ETA                                                             
      SG=SGMAZZ                                                                 
      LI=1                                                                      
      EXSGR2(LI)=COS(SG)+TTI*SIN(SG)                                            
      FL=1.0                                                                    
      DO 710 L=1,LMAXM1                                                         
      DENOM= SQRT(1.0  /(ETASQ+FL*FL))                                          
      LI=LI+1                                                                   
      EXSGR2(LI    )=(FL+TTI*ETA)*EXSGR2(LI-1)*DENOM                            
  710 FL=FL+1.0                                                                 
      LSTEP=(LMAX+19)/20                                                        
      IF(KTLOUT(7).NE.0) WRITE(6,720)(EXSGR2(L),L=1,LMAX,LSTEP)                 
  720 FORMAT(3H EX,8E15.7)                                                      
      CALL INEXCM(1)                                                            
      CALL POTEMS                                                               
      IF(KINEX.EQ.2) GO TO 727                                                  
      IMAX=NXPOT+14                                                             
      IF(KEXCOM(1).NE.0)  IMAX=NXPOT+2-KEXCOM(1)                                
      DO 722 I=1,IMAX                                                           
      VCENR1(I)=VCENTR(I)                                                       
      VCENI1(I)=VCENTI(I)                                                       
      VCOUL1(I)=VCOULM(I)                                                       
  722 CONTINUE                                                                  
      KEX14=KEXCOM(14)                                                          
      DO 724 L=1,KEX14                                                          
      F1(L)=F(L)                                                                
      G1(L)=G(L)                                                                
      FD1(L)=FD(L)                                                              
      GD1(L)=GD(L)                                                              
      EXSGR1(L)=EXSGR2(L)                                                       
  724 CONTINUE                                                                  
      WNI=WN                                                                    
      GO TO 800                                                                 
  727 WNE=WN                                                                    
  800 RETURN                                                                    
      END                                                                       
      SUBROUTINE FLGLCH                                                         
      COMMON/COU/ F(201),G(201),FD(200),GD(200),ETA,SIGMAZ,RHOMX,RD,            
     1            LMAX,KTOUT7                                                   
      DIMENSION GT(2,2),WRONSK(200)                                             
CCCCCC      ******      CALCULATION OF G      ******                            
      ETATW=ETA+ETA                                                             
      ETASQ=ETA*ETA                                                             
      LMAX1=LMAX+1                                                              
      LCONV=LMAX+300                                                            
      NNN=0                                                                     
      DO 200 L=1,LMAX1                                                          
  200 WRONSK(L)=0.0                                                             
      EEPS1=1.E-10                                                              
      EEPS2=1.E-15                                                              
      SQ= SQRT(1.0  +ETASQ)                                                     
      M=1                                                                       
      RHOMXG=RHOMX+20.0  *RD                                                    
      NINC=0                                                                    
  250 SL1=1.0                                                                   
      TL1=0.                                                                    
      SC1=0.0                                                                   
      TC1=1.0  -ETA/RHOMXG                                                      
      SL=SL1                                                                    
      TL=TL1                                                                    
      SC=SC1                                                                    
      TC=TC1                                                                    
      IF( ABS(ETA)-EEPS1) 263,263,255                                           
  255 DO 260 N=1,25                                                             
      T1=N                                                                      
      T2=2.0  *T1-1.0                                                           
      T3=T1*(T1-1.0  )                                                          
      DENOM=2.0  *RHOMXG*T1                                                     
      C1=(ETA*T2)/DENOM                                                         
      C2=(ETASQ-T3)/DENOM                                                       
      SL2=C1*SL1-C2*TL1                                                         
      TL2=C1*TL1+C2*SL1                                                         
      SC2=C1*SC1-C2*TC1-SL2/RHOMXG                                              
      TC2=C1*TC1+C2*SC1-TL2/RHOMXG                                              
      SL=SL+SL2                                                                 
      TL=TL+TL2                                                                 
      SC=SC+SC2                                                                 
      TC=TC+TC2                                                                 
      SL1=SL2                                                                   
      TL1=TL2                                                                   
      SC1=SC2                                                                   
      TC1=TC2                                                                   
  260 CONTINUE                                                                  
  263 IF(M-1) 265,265,310                                                       
  265 IF( ABS(SL*TC-SC*TL-1.0  )-EEPS1)310,270,270                              
  270 NINC=NINC+1                                                               
      IF(NINC-50) 275,280,280                                                   
  275 RHOMXG=RHOMXG+20.0  *RD                                                   
      GO TO 250                                                                 
  280 NNN=1                                                                     
      GO TO 900                                                                 
  310 ARG=SIGMAZ+RHOMXG-ETA*ALOG(2.0  *RHOMXG)                                  
      SINE= SIN(ARG)                                                            
      COSI= COS(ARG)                                                            
      GT(1,M)=SL*COSI-TL*SINE                                                   
      GDT=SC*COSI-TC*SINE                                                       
      T1=SQ                                                                     
      GT(2,M)=((ETA+1.0  /RHOMXG)*GT(1,M)-GDT)/T1                               
      IF(M-1) 320,320,330                                                       
  320 M=M+1                                                                     
      RHOMXG=RHOMXG+RD                                                          
      GO TO 250                                                                 
CCCCCC      ******      INWARD SOLUTION OF DIFFERENTIAL EQUATION      **        
  330 RHOMXG=RHOMXG-RD                                                          
      IF( ABS(RHOMXG-RHOMX)-0.0001  ) 335,335,340                               
  335 G(1)=GT(1,1)                                                              
      G(2)=GT(2,1)                                                              
      GO TO 360                                                                 
  340 RDSQ=RD*RD                                                                
      RDSQ56=0.83333333  *RDSQ                                                  
      RDSQ12=0.1  *RDSQ56                                                       
      NMAX=(RHOMXG-RHOMX)/RD+0.1                                                
      DO 350 L=1,2                                                              
      RHOTM=RHOMXG                                                              
      FL=L*(L-1)                                                                
      PSI1=GT(L,2)                                                              
      PSI2=GT(L,1)                                                              
      R1=RHOTM+RD                                                               
      R2=RHOTM                                                                  
      R3=RHOTM-RD                                                               
      DO 345 N=1,NMAX                                                           
      DENOM=1.0  -RDSQ12*(FL/(R3*R3)+ETATW/R3-1.0  )                            
      FAC1=(1.0  -RDSQ12*(FL/(R1*R1)+ETATW/R1-1.0  ))/DENOM                     
      FAC2=(2.0  +RDSQ56*(FL/(R2*R2)+ETATW/R2-1.0  ))/DENOM                     
      PSI3=FAC2*PSI2-FAC1*PSI1                                                  
      PSI1=PSI2                                                                 
      PSI2=PSI3                                                                 
      R1=R2                                                                     
      R2=R3                                                                     
      R3=R3-RD                                                                  
  345 CONTINUE                                                                  
      G(L)=PSI3                                                                 
  350 CONTINUE                                                                  
  360 RO=RHOMX                                                                  
      ETAC=ETASQ                                                                
      GO=G(1)                                                                   
      GPO=(ETA+1.  /RO)*G(1)-SQ*G(2)                                            
      GD(1)=GPO                                                                 
      SUMP=0.                                                                   
      SUM=0.                                                                    
      L=0                                                                       
  410 L=L+1                                                                     
      L1=L+1                                                                    
      ZL=L                                                                      
      Z1= SQRT(ETAC+ZL*ZL)/ZL                                                   
      Z2=ETA/ZL+ZL/RO                                                           
      GN=(Z2*GO-GPO)/Z1                                                         
      GPN=Z1*GO-Z2*GN                                                           
      IF(L.GT.LMAX) GO TO 420                                                   
      GO=GN                                                                     
      G(L1)=GO                                                                  
      GPO=GPN                                                                   
      GD(L1)=GPO                                                                
      GO TO 410                                                                 
  420 S=1.  /(Z1*GN*GO)                                                         
      SP=(Z1*Z1-Z2*Z2)/(Z1*GPN*GPO)                                             
      SUM=SUM+S                                                                 
      SUMP=SUMP+SP                                                              
      GO=GN                                                                     
      GPO=GPN                                                                   
      LGCONV=L                                                                  
      IF( ABS(S/SUM).LT.EEPS1.AND. ABS(SP/SUMP).LT.EEPS2) GO TO 430             
      IF(L.LT.LCONV) GO TO 410                                                  
      WRITE(6,425)  L,GO,GPO,SUM,SUMP                                           
  425 FORMAT (18H NO CONV IN COULOM,I5,4E15.6)                                  
      NNN=5                                                                     
      GO TO 900                                                                 
  430 L=LMAX                                                                    
      L1=L+1                                                                    
      FO=G(L1)*SUM                                                              
      F(L1)=FO                                                                  
      FPO=GD(L1)*SUMP                                                           
      FD(L1)=FPO                                                                
  432 L=L-1                                                                     
      L1=L+1                                                                    
      ZL=L1                                                                     
      Z1= SQRT(ETAC+ZL*ZL)/ZL                                                   
      Z2=ETA/ZL+ZL/RO                                                           
      FN=(Z2*FO+FPO)/Z1                                                         
      FPN=Z2*FN-Z1*FO                                                           
      FO=FN                                                                     
      F(L1)=FO                                                                  
      FPO=FPN                                                                   
      FD(L1)=FPO                                                                
      IF(L.GE.1) GO TO 432                                                      
      T1=SQ                                                                     
      DO 590 L=1,LMAX                                                           
      FL=L                                                                      
      FL1=FL+1.0                                                                
      T2= SQRT(FL1*FL1+ETASQ)                                                   
      TS=FL/T1                                                                  
      TEST= ABS(F(L)*G(L+1)-F(L+1)*G(L)-TS)                                     
      WRONSK(L)=TEST                                                            
  580 T1=T2                                                                     
  590 CONTINUE                                                                  
      IF(KTOUT7   .EQ.0) GO TO 1000                                             
      WRITE(6,765)                                                              
  765 FORMAT(10H0IN FLGLCH)                                                     
      WRITE(6,770) ETA,SIGMAZ,RHOMX,RHOMXG,RD,LMAX,NINC,LGCONV                  
  770 FORMAT(28H ETA,SIGMAZ,RHOMX,RHOMXG,RD=5E14.6/18H LMAX,NINC,LGCONV=        
     1 ,3I5)                                                                    
      LSTEP=(LMAX+19)/20                                                        
      WRITE(6,781) (F (L),L=1,LMAX1,LSTEP)                                      
      WRITE(6,782) (G (L),L=1,LMAX1,LSTEP)                                      
      WRITE(6,783) (FD(L),L=1,LMAX ,LSTEP)                                      
      WRITE(6,784) (GD(L),L=1,LMAX ,LSTEP)                                      
      WRITE(6,785) (WRONSK(L),L=1,LMAX,LSTEP)                                   
  781 FORMAT(3H F ,7E15.7)                                                      
  782 FORMAT(3H G ,7E15.7)                                                      
  783 FORMAT(3H FD,7E15.7)                                                      
  784 FORMAT(3H GD,7E15.7)                                                      
  785 FORMAT(3H WR,7E15.7)                                                      
      GO TO 1000                                                                
  900 WRITE(6,920) NNN                                                          
  920 FORMAT(13H0IN FLGL NNN=,I1)                                               
 1000 RETURN                                                                    
      END                                                                       
      SUBROUTINE POTEMS                                                         
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      COMMON COULST(1600)                                                       
      COMMON      VCENR1(500),VCENI1(500),VCOUL1(500),                          
     1            VCENTR(500),VCENTI(500),VCOULM(500)                           
      COMMON FORM(2500)                                                         
      COMMON    PFORM1(4),PFORM2(4),PFORM3(4),PFORM4(4),PFORM5(4),              
     1          VLAMD1(3,3),VLAMD2(3,3),VTERM1(3,3),VTERM2(3,3)                 
     2         ,XMEM(500),EXPDX(4)                                              
      DIMENSION IVWR(7)                                                         
      DATA IVWR/4HX   ,4HVCEN,4HTRAL,4HWCEN,4HTRAL,4HVCOU,4HLOMB/               
      NTTLMS=NXPOT+12                                                           
      XBFAC=TMAS**0.333333333                                                   
      IF(KTRL(4).NE.0) XBFAC=XBFAC+PMAS**0.3333333333                           
      XBAR =RZERO *XBFAC                                                        
      XBARW=RZEROW*XBFAC                                                        
      XBARS=RZEROS*XBFAC                                                        
      XBARC=RZEROC*XBFAC                                                        
      VCLFC2=1.4398650  *CHARGE                                                 
      VCLFC1=VCLFC2*0.5  /XBARC                                                 
      DO 410 ND=1,NDFMES                                                        
C      IF(KEXCOM(1).EQ.0) GO TO 251                                              
C      ND=NDFMES                                                                 
C      NXIMIN=1                                                                  
C      NXMAX2=NXPOT+2                                                            
C      NXIMAX=NXMAX2-KEXCOM(1)                                                   
C      DX=XMES2                                                                  
C      X=DX*KEXCOM(1)                                                            
C      GO TO 265                                                                 
C  251 FND=2.0**(ND-1)                                                           
      FND=2.0**(ND-1)                                                           
      DX=XMES1*FND                                                              
      IF(ND-1) 263,259,263                                                      
  259 X=0.0                                                                     
      NXIMIN=1                                                                  
      NXIMAX=8                                                                  
      GO TO 265                                                                 
  263 NXIMIN=NXIMAX+1                                                           
      NXIMAX=NXIMAX+4                                                           
      IF(ND-NDFMES) 265,264,265                                                 
  264 NXMAX2=NTTLMS+2                                                           
      NXIMAX=NXMAX2                                                             
  265 PFORM1(1)= EXP((X-XBAR )/DFN )                                            
      PFORM1(2)= EXP((X-XBARW)/DFNW)                                            
      PFORM1(3)= EXP((X-XBARS)/DFNS)                                            
      EXPDX(1)=EXP(DX/DFN  )                                                    
      EXPDX(2)=EXP(DX/DFNW )                                                    
      EXPDX(3)=EXP(DX/DFNS )                                                    
      DO 400 NX=NXIMIN,NXIMAX                                                   
      X=X+DX                                                                    
      XMEM(NX)=X                                                                
      DO 267 N=1,3                                                              
      PFORM1(N)=EXPDX(N)*PFORM1(N)                                              
      PFORM2(N)=1.0  /(1.0  +PFORM1(N))                                         
      PFORM3(N)=PFORM1(N)*PFORM2(N)*PFORM2(N)                                   
  267 CONTINUE                                                                  
      VCENTR(NX)=-VSX*PFORM2(1)                                                 
      VCENTI(NX)=-WSX*PFORM2(2)-4.0  *WSF*PFORM3(3)                             
      IF(X-XBARC)271,271,272                                                    
  271 VCOULM(NX)=VCLFC1*(3.0  -((X/XBARC)**2))                                  
      GO TO 400                                                                 
  272 VCOULM(NX)=VCLFC2/X                                                       
  400 CONTINUE                                                                  
  410 CONTINUE                                                                  
      IF(KTLOUT(7).EQ.0) GO TO 600                                              
      NXMXW=NXMAX2                                                              
      IF(KEXCOM(1).NE.0) NXMXW=NXMAX2-KEXCOM(1)                                 
      NX1=1                                                                     
      KT7=1                                                                     
      IF(KTLOUT(7).LE.1) KT7=10                                                 
      IF(KTLOUT(7).LE.1.AND.KEXCOM(1).EQ.0) NX1=22                              
      WRITE(6,4520) IVWR(1)                                                     
      WRITE(6,4525)  (XMEM  (NX),NX=NX1,NXMXW,KT7)                              
      WRITE(6,4520) IVWR(2),IVWR(3)                                             
      WRITE(6,4525)  (VCENTR(NX),NX=NX1,NXMXW,KT7)                              
      WRITE(6,4520) IVWR(4),IVWR(5)                                             
      WRITE(6,4525)  (VCENTI(NX),NX=NX1,NXMXW,KT7)                              
      WRITE(6,4520) IVWR(6),IVWR(7)                                             
      WRITE(6,4525)  (VCOULM(NX),NX=NX1,NXMXW,KT7)                              
 4520 FORMAT(1X,2A4)                                                            
 4525 FORMAT(10E12.4)                                                           
  600 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE INEXCM(KWR)                                                    
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      KWRSEL=2*(KINEX-1)+KWR                                                    
      GO TO (101,151,201,251),KWRSEL                                            
CCCCCC  **********  STORING OF PARAMETERS IN THE INCIDENT CHANNEL  *****        
  101 CHARG1=CHARGE                                                             
      ELAB1=ELAB                                                                
      XMES11=XMES1                                                              
      XMES21=XMES2                                                              
      TMAS1=TMAS                                                                
      PMAS1=PMAS                                                                
      RMAS1=RMAS                                                                
      CE1=CE                                                                    
      ECM1=ECM                                                                  
      SGMAZ1=SGMAZZ                                                             
      WN1=WN                                                                    
      WNINI1=WNINI                                                              
      GO TO 500                                                                 
CCCCCC  **********  RECOVERING PARAMETERS IN THE INCIDENT CHANNEL  *****        
  151 CHARGE=CHARG1                                                             
      ELAB=ELAB1                                                                
      XMES1=XMES11                                                              
      XMES2=XMES21                                                              
      TMAS=TMAS1                                                                
      PMAS=PMAS1                                                                
      RMAS=RMAS1                                                                
      CE=CE1                                                                    
      ECM=ECM1                                                                  
      SGMAZZ=SGMAZ1                                                             
      WN=WN1                                                                    
      WNINI=WNINI1                                                              
      GO TO 500                                                                 
CCCCCC  **********  STORING OF PARAMETERS IN THE EXIT CHANNEL  *********        
  201 CHARG2=CHARGE                                                             
      ELAB2=ELAB                                                                
      XMES12=XMES1                                                              
      XMES22=XMES2                                                              
      TMAS2=TMAS                                                                
      PMAS2=PMAS                                                                
      RMAS2=RMAS                                                                
      CE2=CE                                                                    
      ECM2=ECM                                                                  
      SGMAZ2=SGMAZZ                                                             
      WN2=WN                                                                    
      WNINI2=WNINI                                                              
      GO TO 500                                                                 
CCCCCC  **********  RECOVERING PARAMETERS IN THE EXIT CHANNEL  *********        
  251 CHARGE=CHARG2                                                             
      ELAB=ELAB2                                                                
      XMES1=XMES12                                                              
      XMES2=XMES22                                                              
      TMAS=TMAS2                                                                
      PMAS=PMAS2                                                                
      RMAS=RMAS2                                                                
      CE=CE2                                                                    
      ECM=ECM2                                                                  
      SGMAZZ=SGMAZ2                                                             
      WN=WN2                                                                    
      WNINI=WNINI2                                                              
  500 RETURN                                                                    
      END                                                                       
      SUBROUTINE CBCTRL                                                         
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      COMMON COULST(1600),POTENT(3000),FORM(2500)                               
      DIMENSION NTPRD(50)                                                       
      KPSUM=MOD(KPCA+KPSA+KPCB+KPSB,2)                                          
      LAP1MI=KEXCOM(13)+1                                                       
      LAP1MX=KEXCOM(14)                                                         
      LAP1MX=(1-KTRL(9))*LAP1MX+KTRL(9)*(KEXCOM(13)+1)                          
      IF(KTRLCB(7).EQ.1) GO TO 271                                              
      NABLOK=NATOTL*LBSTTL                                                      
      ITPRD=1                                                                   
      NTPRD(1)=MIN0(NBREAD,NBTOTL)*NABLOK                                       
      IF(NBTOTL.LE.NBREAD) GO TO 258                                            
      ITPRD=2+(NBTOTL-NBREAD-1)/(NBREAD-3)                                      
      IF(ITPRD.GT.50) WRITE(6,245) ITPRD                                        
  245 FORMAT(1X,20(1H*),'IN CBCTRL, ITPRD=',I5,' GT.50.INCREASE NTPRD')        
      IF(ITPRD.GT.50) STOP                                                      
      I1=ITPRD-1                                                                
      DO 255 I=2,I1                                                             
  255 NTPRD(I)=NABLOK*(NBREAD-3)                                                
      NTPRD(ITPRD)=(MOD(NBTOTL-NBREAD-1,NBREAD-3)+1)*NABLOK                     
  258 IF(KEXCOM(13).EQ.0) GO TO  271                                            
      LM1=KEXCOM(13)                                                            
      DO 270 NT=1,LM1                                                           
      IF(NT.LE.KEXTCB(3)) GO TO 270                                             
      IF(KEXTCB(4).EQ.2.AND.MOD(NT,2).EQ.KPSUM) GO TO 270                       
      DO 260 IR=1,ITPRD                                                         
      NN=NTPRD(IR)                                                              
  260 READ(12) (FORM(N),N=1,NN)                                                 
  270 CONTINUE                                                                  
CCCCCC  ******           BEGINS TO CALCULATE DISTORTED WAVES      ******        
CCCCCC  ******           AND OVERLAP INTEGRAL                     ******        
  271 DO 900 LAP1=LAP1MI,LAP1MX                                                 
      LA=LAP1-1                                                                 
      KEXCOM(11)=LA                                                             
      IF(KEXCOM(12).EQ.0) GO TO 731                                             
      IF(MOD(LA,2).EQ.KPSUM) GO TO 731                                          
      IF(KTRLCB(7).EQ.1) GO TO 900                                              
      IF(KEXTCB(4).EQ.2) GO TO 900                                              
      DO 720 IR=1,ITPRD                                                         
      NN=NTPRD(IR)                                                              
  720 READ(12) (FORM(N),N=1,NN)                                                 
      GO TO 900                                                                 
CCCCCC  ******    CALCULATION OF DISTORTED WAVES IN INCIDENT CHANNEL  **        
  731 KINEX=1                                                                   
      CALL INEXCM(2)                                                            
      MXROW=1                                                                   
      LLROW(1)=2*LA                                                             
      IF(KTLOUT(3).EQ.0) GO TO 790                                              
      DO 750 N=1,MXROW                                                          
  750 LLROW(N)=LLROW(N)/2                                                       
      WRITE(6,762) (LLROW(N),N=1,MXROW)                                         
  762 FORMAT(3H L=,23I5)                                                        
      DO 775 N=1,MXROW                                                          
  775 LLROW(N)=LLROW(N)*2                                                       
  790 CALL OVLAP                                                                
      IF(KTRLCB(5).EQ.1) GO TO 900                                              
CCCCCC  ******    CALCULATION OF DISTORTED WAVES IN EXIT CHANNEL  ******        
CCCCCC  ******    AND THE OVERLAP INTEGRALS                       ******        
      KINEX=2                                                                   
      CALL INEXCM(2)                                                            
      CALL LBSELC                                                               
      IF(KTLOUT(3).EQ.0) GO TO 865                                              
      DO 852 N=1,MXROW                                                          
  852 LLROW(N)=LLROW(N)/2                                                       
      WRITE(6,762) (LLROW(N),N=1,MXROW)                                         
      DO 858 N=1,MXROW                                                          
  858 LLROW(N)=LLROW(N)*2                                                       
  865 CALL OVLAP                                                                
  900 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE LBSELC                                                         
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      LITR=KEXCOM(11)                                                           
      LBMX=LITR+KEXCOM(15)                                                      
      LBMI1=IABS(LITR-KEXCOM(15))                                               
      LBMI2=IABS(LITR-KEXCOM(10))                                               
      LBMI=MIN0(LBMI1,LBMI2)                                                    
      MXROW=(LBMX-LBMI)/2+1                                                     
      LBTW=2*LBMX                                                               
      DO 50 M1=1,MXROW                                                          
      LLROW(M1)=LBTW                                                            
   50 LBTW=LBTW-4                                                               
      RETURN                                                                    
      END                                                                       
      SUBROUTINE OVLAP                                                          
      COMMON/CRAC/FACLOG(500),IA,IB,IC,ID,IE,IF,RAC,L9(10),U9                   
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      COMMON/OVLC/NXIMIN,NXIMAX,NXMAX2,NXMAX3,NX,NNX,NXPRCH,NXPOTS,ND,          
     1            KSTORE,NOVE,NOVEMX,KUNDER,KFOMCK(40),CFUNIT,X,DX,DRSQ,        
     2            SIMPFC,PRE41,PRE42,PRE43,PRE44,PRE45                          
      COMMON/ELAS/ SCAMP(200,2)                                                 
      COMPLEX SCAMP                                                             
      COMMON/COU/ F(201),G(201),FD(200),GD(200)                                 
      COMMON F1(200),G1(200),FD1(200),GD1(200),EXSGR1(200),EXSGR2(200)          
      COMPLEX EXSGR1,EXSGR2                                                     
      COMMON              POTENT(3000),FORM(2500)                               
      COMMON     URFORD(15,5),UIFORD(15,5),SPFACT(15),                          
     1           FPRER(5,15),FPREI(5,15),FPRERM(4,15),FPREIM(4,15),             
     2           UCR1M(15),UCI1M(15),COEFST(15),                                
     3           UCR1(15),UCI1(15),UCR2(15),UCI2(15),                           
     4           URISV(1500),OVE(40)                                            
      COMPLEX    COEFST,URISV,OVE,URITTT,EX,EXRI,UFFD                           
      INTEGER PLUS,MINUS                                                        
      DATA PLUS,MINUS/3H(+),3H(-)/                                              
      DATA DD20,DDM20/ 1.E20, 1.0E-20/                                          
CCCCCC  ********     PREPARE FOR SKIPPING SOME OF THE OVERLAP  *********        
CCCCCC  ********     INTEGRALS FOR LOWER VALUES OF LA          *********        
      WNUTSQ=0.2187066**2                                                       
      CFUNIT=1.0/(RMAS*WNUTSQ)                                                  
      IF(KINEX.EQ.1) GO TO 41                                                   
      NOVRUN=0                                                                  
      LITW=2*KEXCOM(11)                                                         
      DO 10 M1=1,MXROW                                                          
      LETW=LLROW(M1)                                                            
      DO 10 JLS=1,JLSMAX                                                        
      LFOMTW=LLTWR(JLS)                                                         
      K1=LFOMTW-IABS(LITW-LETW)                                                 
      K2=LITW+LETW-LFOMTW                                                       
      K3=MIN0(K1,K2)                                                            
      NOVRUN=NOVRUN+1                                                           
      KFOMCK(NOVRUN)=1                                                          
      IF(K3.LT.0) KFOMCK(NOVRUN)=0                                              
   10 CONTINUE                                                                  
      DO 30 NOVE=1,NOVRUN                                                       
   30 OVE(NOVE  )=ZERO                                                          
CCCCCC  *****   SIMPSON FACTOR $ MULTIPLIER FOR STORMER INTEGRATION  ***        
   41 DM20LG=ALOG(DDM20)                                                        
      NXMAX2=NXMAX+2                                                            
      SIMPFC=1.                                                                 
      PRE41=19.  /240.                                                          
      PRE42=-2.  /5.                                                            
      PRE43=97.  /120.                                                          
      PRE44=-11.  /15.                                                          
      PRE45=299.  /240.                                                         
      KUNDER=0                                                                  
      IF(KEXCOM(1).NE.0) GO TO 716                                              
      WNS=WNINI                                                                 
      LL=LLROW(1)/2+1                                                           
      FL=LL                                                                     
      U1=FL*ALOG(2.  *XMES1*WNS)+FACLOG(LL)-FACLOG(2*LL)                        
      IF(U1.LE.DM20LG) KUNDER=1                                                 
  716 KSTORE=2-KINEX                                                            
CCCCCC    ********               ND DO-LOOP BEGINS HERE        *********        
      DO 1890 ND=1,NDFMES                                                       
      NXPOTS=4*ND-5                                                             
      FND=2.0**(ND-1)                                                           
      DX=XMES1*FND                                                              
C      IF(KEXCOM(1).EQ.0) GO TO 753                                              
C      ND=4                                                                      
C      DX=XMES2                                                                  
C      KEX1=KEXCOM(1)                                                            
C      NXPOTS=-1-KEXCOM(1)                                                       
C  753 RH=DX*WN                                                                  
      RH=DX*WN                                                                  
      DRSQ=(RH*RH)/ECM                                                          
      IF(KEXCOM(1).EQ.0.AND.ND.NE.1) GO TO 860                                  
CCCCCC  ******    DEFINE WAVE FUNCTION AT THE LOWEST MESH POINT(S)  ****        
      DO 805 M1=1,MXROW                                                         
      UCR1(M1)=0.0                                                              
      UCI1(M1)=0.0                                                              
      UCR2(M1)=0.0                                                              
      IF(KEXCOM(1).NE.0) UCR2(M1)=DDM20                                         
      UCI2(M1)=0.0                                                              
      DO 805 N=1,5                                                              
      FPRER(N,M1)=0.0                                                           
      FPREI(N,M1)=0.0
	IF (N.EQ.5) GO TO 805                                                           
      FPRERM(N,M1)=0.0                                                          
      FPREIM(N,M1)=0.0                                                          
  805 CONTINUE                                                                  
C      IF(KEXCOM(1).EQ.0) GO TO 826                                              
C      X=KEX1*DX                                                                 
C      NXIMIN=KEX1+2                                                             
C      NXIMAX=NXMAX+2                                                            
C      NXMAX3=NXIMAX-5                                                           
C      GO TO 1225                                                                
C  826 DO 835 M1=1,MXROW                                                         
      DO 835 M1=1,MXROW                                                         
      L1=LLROW(M1)                                                              
      LL=(L1/2)+1                                                               
      WNS=WNINI                                                                 
      IF(KUNDER.EQ.0) GO TO 831                                                 
      UCR2TT=DDM20                                                              
      GO TO 832                                                                 
  831 UCR2TT  =0.5  *((2.0  *DX*WNS)**LL)* EXP(FACLOG(LL)-FACLOG(2*LL))         
  832 UCR2(M1)=UCR2TT                                                           
      UCI2(M1)=0.0                                                              
      IF(LL.NE.2) GO TO 835                                                     
      FPRETT=CFUNIT*0.6666666666  *WNS*WNS*DRSQ                                 
      FPRERM(1,M1)=FPRETT                                                       
      FPRER(1,M1)=FPRERM(1,M1)                                                  
      FPREIM(1,M1)=0.0                                                          
      FPREI(1,M1)=FPREIM(1,M1)                                                  
  835 CONTINUE                                                                  
CCCCCC  *****   FIX THE RANGE OF MESH-POINTS FOR A GIVEN VALUE OF ND  **        
      X=0.0                                                                     
      NXIMIN=2                                                                  
      NXIMAX=8                                                                  
      NXPRCH=3                                                                  
      NNX=2                                                                     
      IF(NDFMES.NE.1) GO TO 1225                                                
      NXIMAX=NXMAX+2                                                            
      NXMAX3=NXIMAX-5                                                           
      GO TO 1225                                                                
  860 NXIMIN=5                                                                  
      X=3.0  *DX                                                                
      IF(ND.EQ.NDFMES) GO TO 863                                                
      NXIMAX=8                                                                  
      GO TO 865                                                                 
  863 NXIMAX=NXMAX+2                                                            
      NXMAX3=NXIMAX-5                                                           
  865 DO 875 M1=1,MXROW                                                         
      DO 870 NQ=1,4                                                             
      IF(KSTORE.EQ.1) URISV(NQ )=URISV(2*NQ)                                    
      FPRER(NQ,M1)=FPRERM(NQ,M1)*4.0                                            
  870 FPREI(NQ,M1)=FPREIM(NQ,M1)*4.0                                            
      FPRERM(1,M1)=FPRER(1,M1)                                                  
      FPRERM(2,M1)=FPRER(3,M1)                                                  
      FPREIM(1,M1)=FPREI(1,M1)                                                  
      FPREIM(2,M1)=FPREI(3,M1)                                                  
      UCR1(M1)=UCR1M(M1)                                                        
      UCI1(M1)=UCI1M(M1)                                                        
  875 CONTINUE                                                                  
      NXPRCH=5                                                                  
      NNX=3                                                                     
 1225 CALL INTGMS                                                               
      IF(ND-NDFMES) 1890,1805,1890                                              
CCCCCC  ****  CALCULATE DERIVATIVE OF WAVE FUNCTION AT MATCHING RADIUS**        
 1805 DO 1850 M1=1,MXROW                                                        
      UCR1(M1)=URFORD(M1,3)                                                     
      UCI1(M1)=UIFORD(M1,3)                                                     
      DENOM=1.0/(12.0*DX*WN)                                                    
      UCR2(M1)   =(8.0  *(URFORD(M1,4)-URFORD(M1,2))                            
     1                -(URFORD(M1,5)-URFORD(M1,1)))*DENOM                       
      UCI2(M1)   =(8.0  *(UIFORD(M1,4)-UIFORD(M1,2))                            
     1                -(UIFORD(M1,5)-UIFORD(M1,1)))*DENOM                       
 1850 CONTINUE                                                                  
 1890 CONTINUE                                                                  
 1900 CONTINUE                                                                  
CCCCCC  ******    MATCHING OF DISTORTED WAVES TO EXTERNAL SOLUTION  ****        
      DO 1959 M1=1,MXROW                                                        
      LP1=LLROW(M1)/2+1                                                         
      IF(KINEX.EQ.2) GO TO 1957                                                 
      EX=(UCR1(M1)  +TTI*UCI1(M1)  )*(GD1(LP1  )+TTI*FD1(LP1  ))                
     1  -(UCR2(M1)  +TTI*UCI2(M1)  )*(G1 (LP1  )+TTI*F1 (LP1  ))                
      COEFST(M1)  =-EXSGR1(LP1  )/EX                                            
      UFFD=(UCR2(M1)+TTI*UCI2(M1))*F1 (LP1)                                     
     1    -(UCR1(M1)+TTI*UCI1(M1))*FD1(LP1)                                     
      SCAMP(LP1,1)=UFFD/EX                                                      
      GO TO 1959                                                                
 1957 EX=(UCR1(M1)  +TTI*UCI1(M1)  )*(GD (LP1  )+TTI*FD (LP1  ))                
     1  -(UCR2(M1)  +TTI*UCI2(M1)  )*(G  (LP1  )+TTI*F  (LP1  ))                
      COEFST(M1)  =-EXSGR2(LP1  )/EX                                            
      UFFD=(UCR2(M1)+TTI*UCI2(M1))*F  (LP1)                                     
     1    -(UCR1(M1)+TTI*UCI1(M1))*FD (LP1)                                     
      SCAMP(LP1,2)=UFFD/EX                                                      
 1959 CONTINUE                                                                  
      IF(KTLOUT(10).EQ.0) GO TO 2105                                            
      WRITE(6,1963) (COEFST(M1)  ,M1=1,MXROW)                                   
 1963 FORMAT(/ 20H OVLAP-1963. COEFST=,3(E14.6,E15.6,3X))                       
 2105 IF(KINEX.EQ.2) GO TO 2210                                                 
CCCCCC  ****   RENORMALIZATION OF INCIDENT CHANNEL DISTORTED WAVE   ****        
      EXRI=COEFST(1)                                                            
      U1=URISV(MSTOT)*EXRI                                                      
      TEST= ABS(U1)*DDM20                                                       
      DO 2120 NX=1,MSTOT                                                        
      EX=URISV(NX )*EXRI                                                        
      EXR=EX                                                                    
      EXI=-TTI*EX                                                               
      IF( ABS(EXR).LT.TEST) EXR=0.                                              
      IF( ABS(EXI).LT.TEST) EXI=0.                                              
      URISV(NX )=EXR+TTI*EXI                                                    
 2120 CONTINUE                                                                  
      IF(KTLOUT(10).EQ.0)  GO TO 3000                                           
      WRITE(6,2173)KEXCOM(11)                                                   
 2173 FORMAT( 24H COUPMS-2120. URISV. LA=  ,I3)                                 
      WRITE(6,2175) (URISV(NX9),NX9=10,MSTOT,10)                                
 2175 FORMAT(5(E11.3,E11.3,2X))                                                 
      GO TO 3000                                                                
CCCCCC  ******    RENORMALIZATION OF THE OVERLAP INTEGRAL DUE TO  ******        
CCCCCC  ******    RENORMALIZATION OF THE DISTORTED WAVES IN EXIT  ******        
CCCCCC  ******    CHANNEL                                         ******        
 2210 NOVE=0                                                                    
      NOVRUN=0                                                                  
      DO 2265 M1=1,MXROW                                                        
      EX=COEFST(M1)                                                             
      DO 2263 JLS=1,JLSMAX                                                      
      NOVRUN=NOVRUN+1                                                           
      IF(KFOMCK(NOVRUN).EQ.0) GO TO 2263                                        
      NOVE=NOVE+1                                                               
      OVE(NOVE)=OVE(NOVE)*EX                                                    
 2263 CONTINUE                                                                  
 2265 CONTINUE                                                                  
      IF(NOVEMX.EQ.0)  NOVEMX=1                                                 
      WRITE (11)    NOVEMX,(OVE (NOVE),NOVE=1,NOVEMX)                           
      IF(KTLOUT(10).EQ.0) GO TO 3000                                            
      WRITE (6,2272)                                                            
 2272 FORMAT(18H   OVLAP-2272. OVE)                                             
      WRITE (6,2274)  (OVE (NOVE ),NOVE =1,NOVEMX)                              
 2274 FORMAT(5(E11.3,E11.3,2X))                                                 
      IF(KTLOUT(10).LT.2) GO TO 3000                                            
      WRITE(6,2276) (KFOMCK(N),N=1,NOVRUN)                                      
 2276 FORMAT(20H OVLAP-2276. KFOMCK=,20I5)                                      
 3000 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE INTGMS                                                         
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      COMMON/OVLC/NXIMIN,NXIMAX,NXMAX2,NXMAX3,NX,NNX,NXPRCH,NXPOTS,ND,          
     1            KSTORE,NOVE,NOVEMX,KUNDER,KFOMCK(40),CFUNIT,X,DX,DRSQ,        
     2            SIMPFC,PRE41,PRE42,PRE43,PRE44,PRE45                          
      COMMON COULST(1600)                                                       
      COMMON      VCENR1(500),VCENI1(500),VCOUL1(500),                          
     1            VCENTR(500),VCENTI(500),VCOULM(500)                           
      COMMON FORM(2500)                                                         
      COMMON     URFORD(15,5),UIFORD(15,5),SPFACT(15),                          
     1           FPRER(5,15),FPREI(5,15),FPRERM(4,15),FPREIM(4,15),             
     2           UCR1M(15),UCI1M(15),COEFST(15),                                
     3           UCR1(15),UCI1(15),UCR2(15),UCI2(15),                           
     4           URISV(1500),OVE(40)                                            
      COMPLEX    COEFST,URISV,OVE,URITTT,EX,EXRI                                
      DIMENSION CENTFG(15),URISVT(15),BR(15),BI(15)                             
      COMPLEX URISVT                                                            
      DO 1205 M1=1,MXROW                                                        
      LL=LLROW(M1)/2                                                            
      CENTFG(M1)=LL*(LL+1)*CFUNIT                                               
 1205 CONTINUE                                                                  
      KEXCB1=KEXTCB(1)                                                          
      KECB11=KEXCB1+1                                                           
      MSTTT=MSTOT+KEXCB1                                                        
      SIMPFC=XMES2/3.                                                           
      SIMP1=SIMPFC                                                              
      SIMP2=SIMPFC*2.                                                           
      SIMP4=SIMPFC*4.                                                           
CCCCCC  ******    SOLVE DIFFERENTIAL EQUATION FOR DISTORTED WAVE  ******        
      DO 2000 NX=NXIMIN,NXIMAX                                                  
      IF(KINEX.NE.2) GO TO 1233                                                 
      IF(KEXCB1.EQ.0) GO TO 1221                                                
      IF(NX-KECB11) 1233,1211,1230                                              
CCCCCC  ******    PREPARE SIMPSON MULTIPLIER FOR RB INTEGRATION   ******        
 1211 IIPP=0                                                                    
      FFF1=SIMP1                                                                
      GO TO 1233                                                                
 1221 IF(ND-2) 1226,1227,1228                                                   
 1226 FFF1=SIMP4                                                                
      IIPP=1                                                                    
      GO TO 1233                                                                
 1227 FFF1=SIMP2                                                                
      IIPP=0                                                                    
      GO TO 1233                                                                
 1228 IF(ND-3) 1233,1229,1230                                                   
 1229 IF(NX-6) 1233,1226,1227                                                   
 1230 IF(IIPP.NE.0) GO TO 1232                                                  
      IIPP=1                                                                    
      FFF1=SIMP4                                                                
      GO TO 1233                                                                
 1232 IIPP=0                                                                    
      FFF1=SIMP2                                                                
 1233 CONTINUE                                                                  
CCCCCC      ******      DIAGONAL ELEMENTS      *************************        
      N3=(NX+NXPOTS)                                                            
      X=X+DX                                                                    
      XX2=X*X                                                                   
      XX2INV=1.0  /XX2                                                          
      DO 1700 M1=1,MXROW                                                        
      CFFCLL=CENTFG(M1)                                                         
      IF(NX.GT.NXPOT+2) GO TO 1315                                              
CCCCCC  ******   MULTIPLICATION OF THE WAVE FUNCTION WITH POTENTIAL  ***        
      IF(KINEX.EQ.2) GO TO 1312                                                 
      AR2=VCENR1(N3)+VCOUL1(N3)+CFFCLL*XX2INV-ECM                               
      AI2=VCENI1(N3)                                                            
      GO TO 1320                                                                
 1312 AR2=VCENTR(N3)+VCOULM(N3)+CFFCLL*XX2INV-ECM                               
      AI2=VCENTI(N3)                                                            
      GO TO 1320                                                                
 1315 AR2=CFFCLL*XX2INV+1.4398650  *CHARGE/X-ECM                                
      AI2=0.                                                                    
 1320 BR2=AR2*UCR2(M1)-AI2*UCI2(M1)                                             
      BI2=AR2*UCI2(M1)+AI2*UCR2(M1)                                             
      BR(M1)=BR2*DRSQ                                                           
      BI(M1)=BI2*DRSQ                                                           
 1700 CONTINUE                                                                  
      NOVE=0                                                                    
      NOVRUN=0                                                                  
      DO 1790 M1=1,MXROW                                                        
      IF(ND.NE.1) GO TO 1704                                                    
      IF(NX.GT.4) GO TO 1704                                                    
      TERMR=BR(M1)                                                              
      FPRER(NX,M1)=TERMR                                                        
      TERMI=BI(M1)                                                              
      FPREI(NX,M1)=TERMI                                                        
      GO TO 1705                                                                
 1704 FPRER(5,M1)=BR(M1)                                                        
      FPREI(5,M1)=BI(M1)                                                        
      TERMR =PRE41*FPRER (1,M1)+PRE42*FPRER (2,M1)+PRE43*FPRER (3,M1)           
     1      +PRE44*FPRER (4,M1)+PRE45*FPRER (5,M1)                              
      TERMI =PRE41*FPREI (1,M1)+PRE42*FPREI (2,M1)+PRE43*FPREI (3,M1)           
     1      +PRE44*FPREI (4,M1)+PRE45*FPREI (5,M1)                              
 1705 ARR=2.0  *UCR2(M1)-UCR1(M1)+TERMR                                         
      AII=2.0  *UCI2(M1)-UCI1(M1)+TERMI                                         
      UCR1(M1)=UCR2(M1)                                                         
      UCR2(M1)=ARR                                                              
      UCI1(M1)=UCI2(M1)                                                         
      UCI2(M1)=AII                                                              
CCCCCC      ******      OVERLAP INTEGRAL      **************************        
      IF(KINEX.NE.2) GO TO 1771                                                 
CCCCCC  ******          PREPARATION FOR CALLING INTPOL        **********        
      NS=ND                                                                     
      IF(ND-3) 1713,1711,1714                                                   
 1711 IF(NX-6) 1771,1715,1712                                                   
 1712 NS=4                                                                      
 1713 IF(NX-8) 1771,1715,1771                                                   
 1714 NS=NX                                                                     
 1715 IF(NS.LE.KEXCB1) GO TO 1771                                               
      IF(NS.GT.MSTTT ) GO TO 1771                                               
      NT=NS-KEXCB1                                                              
      NTCOM=NT                                                                  
      IF(KTRLCB(8).EQ.0) GO TO 1736                                             
      IF(M1.NE.1) GO TO 1736                                                    
      IF(NT.NE.1) GO TO 1732                                                    
      NBRED3=NBREAD-3                                                           
      NTRDCK=NBRED3*MESFCB                                                      
      NBREDT=MIN0(NBREAD,NBTOTL)                                                
      GO TO 1734                                                                
 1732 MODNT=MOD(NT,NTRDCK)                                                      
      IF(MODNT.NE.1) GO TO 1736                                                 
      NBLEFT=NBTOTL-NBLMCM                                                      
      NBREDT=MIN0(NBRED3,NBLEFT)                                                
 1734 CALL INTPOL                                                               
CCCCCC  ******    IN INTPOL, INTEGRATION OVER RA IS CARRIED OUT   ******        
CCCCCC  ******    FOR A BLOCK OF RB VALUES (FOR EFR CASE ONLY)    ******        
 1736 URITTT=UCR2(M1)+TTI*UCI2(M1)                                              
      DO 1750 JLS=1,JLSMAX                                                      
      IF(KTRLCB(8).EQ.1) GO TO 1742                                             
CCCCCC  ********     OVERLAP INTEGRAL FOR NO-RECOIL DWBA       *********        
      NOVRUN=NOVRUN+1                                                           
      IF(KFOMCK(NOVRUN).EQ.0) GO TO 1750                                        
      LINEFM=(JLS -1)*MSTOT+NT                                                  
      F1=FFF1*FORM(LINEFM)                                                      
      EX=F1*URITTT*URISV(NT)                                                    
      GO TO 1744                                                                
CCCCCC  ******    INTEGRATION OVER RB FOR EFR-DWBA CALCULATION    ******        
 1742 NOVRUN=NOVRUN+1                                                           
      F1=FFF1*XMES2*NS                                                          
      IF(KFOMCK(NOVRUN).EQ.0) GO TO 1750                                        
      NX9=MSTOT+NOVE*NTRDCK+NT-NBSMCM                                           
      EX=URITTT*URISV(NX9)*F1                                                   
 1744 NOVE=NOVE+1                                                               
      NOVEMX=NOVE                                                               
      OVE(NOVE)=OVE(NOVE)+EX                                                    
 1750 CONTINUE                                                                  
 1771 IF(ND.EQ.NDFMES) GO TO 1774                                               
      IF(NX.NE.NXPRCH) GO TO 1774                                               
      FPRERM(NNX,M1)=BR(M1)                                                     
      FPREIM(NNX,M1)=BI(M1)                                                     
      IF(M1.NE.MXROW) GO TO 1773                                                
      NXPRCH=NXPRCH+2                                                           
      NNX=NNX+1                                                                 
 1773 IF(NX-NXIMAX+1.NE.0) GO TO 1774                                           
      UCR1M(M1)=UCR1(M1)                                                        
      UCI1M(M1)=UCI1(M1)                                                        
 1774 IF(ND.NE.1) GO TO 1775                                                    
      IF(NX.LE.4) GO TO 1790                                                    
 1775 DO 1778 NQ=1,4                                                            
      FPRER(NQ,M1)=FPRER(NQ+1,M1)                                               
 1778 FPREI(NQ,M1)=FPREI(NQ+1,M1)                                               
 1790 CONTINUE                                                                  
      IF(ND.NE.NDFMES) GO TO 1844                                               
CCCCCC  ******    PREPARATION FOR TAKING DERIVATIVE OF DISTORTED  ******        
CCCCCC  ******    WAVE AT THE MATCHING RADIUS                     ******        
      NXCH2=NX-NXMAX3                                                           
      IF(NXCH2.LE.0) GO TO 1844                                                 
      DO 1842 M1=1,MXROW                                                        
      URFORD(M1,NXCH2)=UCR2(M1)                                                 
 1842 UIFORD(M1,NXCH2)=UCI2(M1)                                                 
 1844 IF(KSTORE.EQ.0) GO TO 1847                                                
      IF(NX.LE.KEXCB1) GO TO 1847                                               
      IF(NX.GT.MSTTT ) GO TO 1847                                               
      NT=NX-KEXCB1                                                              
CCCCCC  ******    STORING OF DISTORTED WAVE IN THE INITIAL CHANNEL  ****        
      URITTT=UCR2( 1)+TTI*UCI2( 1)                                              
 1845 URISV(NT )=URITTT                                                         
CCCCCC  ******  OUTPUT OF D.W. AT EVERY MESH POINT  ********************        
 1847 IF(KTLOUT(10)+KINEX.LT.4) GO TO 1855                                      
      URISVT(    1)=UCR2(    1)+TTI*UCI2(    1)                                 
      URISVT(MXROW)=UCR2(MXROW)+TTI*UCI2(MXROW)                                 
      KWFOUT=0                                                                  
      IF(KEXCOM(11).EQ.0) KWFOUT=1                                              
      IF(KTRL(9).NE.0) KWFOUT=1                                                 
      IF(KWFOUT.EQ.0) GO TO 1855                                                
      IF(NX.LE.15) GO TO 1851                                                   
      IF(NX.GE.NXMAX3) GO TO 1851                                               
      NXMOD=MOD(NX,5)                                                           
      IF(NXMOD.NE.0) GO TO 1855                                                 
 1851 XX=X+DX                                                                   
      WRITE(6,1853)NX,XX,URISVT(1),URISVT(MXROW),OVE(1),OVE(NOVEMX)             
 1853 FORMAT(4H NX=,I3,3H X=,F9.4,4(E12.3,E12.3))                               
 1855 CONTINUE                                                                  
      CALL URENOM                                                               
 2000 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE INTPOL                                                         
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      COMMON/OVLC/NXIMIN,NXIMAX,NXMAX2,NXMAX3,NX,NNX,NXPRCH,NXPOTS,ND,          
     1            KSTORE,NOVE,NOVEMX,KUNDER,KFOMCK(40),CFUNIT,X,DX,DRSQ,        
     2            SIMPFC,PRE41,PRE42,PRE43,PRE44,PRE45                          
      COMMON COULST(1600),POTENT(3000),FORM(2500)                               
      COMMON     URFORD(15,5),UIFORD(15,5),SPFACT(15),                          
     1           FPRER(5,15),FPREI(5,15),FPRERM(4,15),FPREIM(4,15),             
     2           UCR1M(15),UCI1M(15),COEFST(15),                                
     3           UCR1(15),UCI1(15),UCR2(15),UCI2(15),                           
     4           URISV(1500),OVE(40)                                            
      COMPLEX    COEFST,URISV,OVE,URITTT,EX,EXRI                                
      COMMON FFORA(650),FFA(650),FA1(40),FA2(40),FA3(40),FA4(40),               
     1               FLA(30),FGR(30),FAN(30),FGE(30),                           
     2               AFLA(30),AFGR(30),AFAN(30),AFGE(30)                        
      DIMENSION FTAPE(2500)                                                     
      EQUIVALENCE (FTAPE(1),FORM(1))                                            
      DATA KINT/0/                                                              
      IF(KINT.EQ.1) GO TO 151                                                   
CCCCCC  ******    PREPARE FOR INTERPOLATION OF EFR-FORM FACTOR    ******        
      FMC=MESFCB                                                                
      DO 143 N=1,MESFCB                                                         
      FN=N                                                                      
      DELN=FN/FMC                                                               
      FLA(N) =-DELN*(DELN-1.  )*(DELN-2.  )/6.                                  
      FGR(N) =0.5  *(DELN+1.  )*(DELN-1.  )*(DELN-2.  )                         
      FAN(N) =-0.5  *(DELN+1.  )*DELN*(DELN-2.  )                               
      FGE(N) =(DELN+1.  )*DELN*(DELN-1.  )/6.                                   
  143 CONTINUE                                                                  
      FMC=MESFCA                                                                
      DO 145 N=1,MESFCA                                                         
      FN=N                                                                      
      DELN=FN/FMC                                                               
      AFLA(N)=-DELN*(DELN-1.  )*(DELN-2.  )/6.                                  
      AFGR(N)=0.5  *(DELN+1.  )*(DELN-1.  )*(DELN-2.  )                         
      AFAN(N)=-0.5  *(DELN+1.  )*DELN*(DELN-2.  )                               
      AFGE(N)=(DELN+1.  )*DELN*(DELN-1.  )/6.                                   
  145 CONTINUE                                                                  
      KINT=1                                                                    
  151 MFCATP=MESFCA                                                             
      IF(MESFCA.LT.0) MFCATP=1                                                  
      SIMPA1=XMES2 /3.                                                          
      SIMPA2=2.  *SIMPA1                                                        
      SIMPA4=4.  *SIMPA1                                                        
      NABLOK=NATOTL*LBSTTL                                                      
      NBBLOK=NBTOTL*NABLOK                                                      
      NAFMES=(NATOTL-1)*MFCATP+1                                                
      NBRED3=NBREAD-3                                                           
      NBSTBS=NBRED3*MESFCB                                                      
      IF(NTCOM.NE.1) GO TO 1831                                                 
CCCCCC  ******          PREPARE FOR READING FORM FACTOR FOR       ******        
CCCCCC  ******          THE LOWEST BLOCK OF RB POINTS             ******        
      NMI=1                                                                     
      NMX=NABLOK*NBREAD                                                         
      NBMINT=1                                                                  
      NBMAXT=NBRED3*MESFCB                                                      
      NBSMCM=0                                                                  
      NBLMCM=NBREAD                                                             
      GO TO 1845                                                                
CCCCCC  ********        PREPARE FOR READING FORM FACTOR FOR A   ********        
CCCCCC  ********        RB-BLOCK OTHER THAN THE LOWEST          ********        
 1831 NBSMCM=NBSMCM+(NBMAXT-NBMINT+1)                                           
      NBLMCM=NBLMCM+NBREDT                                                      
      NMI=NABLOK*3+1                                                            
      NMX=NABLOK*(NBREDT+3)                                                     
      NBMINT=NBSMCM+1                                                           
      NBMAXT=NBSMCM+NBREDT*MESFCB                                               
      NSHMX=NABLOK*3                                                            
      NSH=NABLOK*NBRED3                                                         
      DO 1840 N1=1,NSHMX                                                        
      N2=N1+NSH                                                                 
 1840 FTAPE(N1)=FTAPE(N2)                                                       
CCCCCC    ************           READING EFR FORM FACTOR     ***********        
 1845 READ(12) (FTAPE(N),N=NMI,NMX)                                             
      IF(KTLOUT(10).LT.5) GO TO 1907                                            
      IF(KEXCOM(11).GE.4) GO TO 1907                                            
      WRITE(6,1853) KEXCOM(11),NBMINT,NBMAXT,NBLMCM,NBSMCM                      
 1853 FORMAT(/55H INTPOL-1853. FTAPE. KEX11,NBMINT,NBMAXT,NBLMCM,NBSMCM=        
     1,5I5)                                                                     
      WRITE(6,1857) (FTAPE(N),N=1,NMX)                                          
 1857 FORMAT(1X,10E13.4)                                                        
CCCCCC  ********     REARRANGE FORM FACTOR TO PREPARE FOR     **********        
CCCCCC  ********     INTERPOLATION AND INTEGRATION OVER RA    **********        
 1907 FTBASE=0                                                                  
      NZ1=MSTOT+1                                                               
      NZ2=MSTOT+LBSTTL*NBSTBS                                                   
      DO 1911 NZ=NZ1,NZ2                                                        
 1911 URISV(NZ)=ZERO                                                            
      N2BASE=3*NABLOK                                                           
      DO 1980 NB=NBMINT,NBMAXT                                                  
      NBSTOR=NB-NBMINT+1                                                        
      IF(MESFCA.EQ.1.AND.MESFCB.EQ.1) GO TO 1937                                
      INTPCH=MOD(NB,MESFCB)                                                     
      IF(INTPCH.EQ.0) GO TO 1914                                                
      I=INTPCH                                                                  
CCCCCC  ******    INTERPOLATION OF FORM-FACTOR WITH RESPECT TO RB  *****        
      DO 1913 NA=1,NABLOK                                                       
      N1=FTBASE+NA                                                              
      N2=N1+NABLOK                                                              
      N3=N2+NABLOK                                                              
      N4=N3+NABLOK                                                              
      FFA(NA)=FLA(I)*FTAPE(N1)+FGR(I)*FTAPE(N2)                                 
     1       +FAN(I)*FTAPE(N3)+FGE(I)*FTAPE(N4)                                 
 1913 CONTINUE                                                                  
      GO TO 1917                                                                
 1914 FTBASE=FTBASE+NABLOK                                                      
      NNA1=FTBASE+NABLOK                                                        
      DO 1915 NA=1,NABLOK                                                       
      NNA1=NNA1+1                                                               
 1915 FFA(NA)=FTAPE(NNA1)                                                       
 1917 IF(MESFCA.GT.1) GO TO 1921                                                
CCCCCC  ******    ARRANGE FORM-FACTOR IN FFORA(N), WHEN MESFCA=1,  *****        
CCCCCC  ******    TO BE READY FOR INTEGRATION OVER RA              *****        
      DO 1920 N=1,NABLOK                                                        
 1920 FFORA(N)=FFA(N)                                                           
      GO TO 1939                                                                
CCCCCC  ********       INTERPOLATE FORM FACTOR WITH RESPECT    *********        
CCCCCC  ********       TO RA AND ARRANGE IT IN FFORA(N)        *********        
 1921 DO 1930 LBS=1,LBSTTL                                                      
      NABASE=2*LBSTTL                                                           
      FFORA(LBS)=FFA(LBS)                                                       
      FA1(LBS)=0.                                                               
 1926 FA2(LBS)=FFA(LBS)                                                         
      N3=LBS+LBSTTL                                                             
      FA3(LBS)=FFA(N3)                                                          
      N3=N3+LBSTTL                                                              
      FA4(LBS)=FFA(N3)                                                          
      DO 1930 NA=2,NAFMES                                                       
      NALBS=(NA-1)*LBSTTL+LBS                                                   
      INCHA=MOD(NA-1,MFCATP)                                                    
      IF(INCHA.EQ.0) GO TO 1927                                                 
      I=INCHA                                                                   
      FFORA( NALBS)=AFLA(I)*FA1(LBS)+AFGR(I)*FA2(LBS)                           
     1             +AFAN(I)*FA3(LBS)+AFGE(I)*FA4(LBS)                           
      GO TO 1930                                                                
 1927 NABASE=NABASE+LBSTTL                                                      
      NA4=NABASE+LBS                                                            
      FFORA(NALBS)=FA3(LBS)                                                     
      IF(NA.EQ.NAFMES) GO TO 1930                                               
      FA1(LBS)=FA2(LBS)                                                         
      FA2(LBS)=FA3(LBS)                                                         
      FA3(LBS)=FA4(LBS)                                                         
      FA4(LBS)=FFA(NA4)                                                         
      IF(NA+MFCATP.EQ.NAFMES) FA4(LBS)=0.                                       
 1930 CONTINUE                                                                  
      IF(KTLOUT(10).LT.5) GO TO 1939                                            
      IF(NB.GT.5) GO TO 1939                                                    
      WRITE(6,1932) (FFA(NA),NA=1,NABLOK)                                       
 1932 FORMAT(5H FFA ,10E12.4)                                                   
      DO 1934 LBS=1,LBSTTL                                                      
      NA1=LBS                                                                   
      NA2=(NAFMES-1)*LBSTTL+LBS                                                 
      WRITE(6,1933)(FFORA(N),N=NA1,NA2,LBSTTL)                                  
 1933 FORMAT(6H FFORA,10E12.4)                                                  
 1934 CONTINUE                                                                  
      GO TO 1939                                                                
CCCCCC  ******    ARRANGE FORM FACTOR IN FFORA(N) WHEN INTERPOLATION  **        
CCCCCC  ******    IS NOT MADE WITH RESPECT TO RA OR RB            ******        
 1937 NBASE=(NB-1)*NABLOK                                                       
      NALBS=0                                                                   
      DO 1938 NA=1,NATOTL                                                       
      DO 1938 LBS=1,LBSTTL                                                      
      N=NBASE+(NA-1)*LBSTTL+LBS                                                 
      NALBS=NALBS+1                                                             
 1938 FFORA( NALBS)=FTAPE(N)                                                    
CCCCCC  ********      PREPARATION FOR INTEGRATION OVER RA      *********        
 1939 NAINTG=NAFMES                                                             
      IF(MESFCA.GT.0) GO TO 1945                                                
      MSFCAB=-MESFCA                                                            
      NAINTG=((NAFMES-1)/MSFCAB)+1                                              
 1945 NAFHLF=(NAINTG+1)/2                                                       
      NSABAS=NB-NAFHLF                                                          
      NRABAS=NSABAS+KEXTCB(1)                                                   
CCCCCC  **********           INTEGRATION OVER RA BEGINS        *********        
      DO 1970 NA=1,NAINTG                                                       
      NN=NA                                                                     
      IF(MESFCA.LT.0) NN=(NA-1)*MSFCAB+1                                        
      NSA=NSABAS+NA                                                             
      IF(NSA.LT.1) GO TO 1970                                                   
      IF(NSA.GT.MSTOT) GO TO 1970                                               
      RA=XMES2 *(NRABAS+NA)                                                     
      NCK=MIN0(NA,NSA)                                                          
      IF(NCK.NE.1) GO TO 1954                                                   
      SIMPA=SIMPA1                                                              
      IIA=0                                                                     
      GO TO 1958                                                                
 1954 IF(IIA.EQ.1) GO TO 1956                                                   
      IIA=1                                                                     
      SIMPA=SIMPA4                                                              
      GO TO 1958                                                                
 1956 IIA=0                                                                     
      SIMPA=SIMPA2                                                              
 1958 EX=URISV(NSA)*RA*SIMPA                                                    
      NALBS=(NA-1)*LBSTTL                                                       
      DO 1965 LBS=1,LBSTTL                                                      
      NSLBS=MSTOT+(LBS-1)*NBSTBS+NBSTOR                                         
      NALBS=NALBS+1                                                             
      URISV(NSLBS)=URISV(NSLBS)+FFORA( NALBS)*EX                                
 1965 CONTINUE                                                                  
 1970 CONTINUE                                                                  
 1980 CONTINUE                                                                  
 1981 IF(KTLOUT(10).LT.5) GO TO 2000                                            
CCCCCC  ******    OUTPUT OF THE RESULT OF RA-INTEGRATION          ******        
      IF(KEXCOM(11).GT.3) GO TO 2000                                            
      WRITE(6,1992) KEXCOM(11)                                                  
 1992 FORMAT(/44H INTPOL-1992. FF TIMES URISV FOR KEXCOM(11)=,I2)               
      DO 1995 LBS=1,LBSTTL                                                      
      NZ1=MSTOT+(LBS-1)*NBSTBS+1                                                
      NZ2=NZ1+NBSTBS-1                                                          
      WRITE(6,1994) (URISV(NZ),NZ=NZ1,NZ2,10)                                   
 1994 FORMAT(5(E12.4,E12.4,1X))                                                 
 1995 CONTINUE                                                                  
 2000 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE URENOM                                                         
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      COMMON/OVLC/NXIMIN,NXIMAX,NXMAX2,NXMAX3,NX,NNX,NXPRCH,NXPOTS,ND,          
     1            KSTORE,NOVE,NOVEMX,KUNDER,KFOMCK(40),CFUNIT,X,DX,DRSQ,        
     2            SIMPFC,PRE41,PRE42,PRE43,PRE44,PRE45                          
      COMMON COULST(1600),POTENT(3000),FORM(2500)                               
      COMMON     URFORD(15,5),UIFORD(15,5),SPFACT(15),                          
     1           FPRER(5,15),FPREI(5,15),FPRERM(4,15),FPREIM(4,15),             
     2           UCR1M(15),UCI1M(15),COEFST(15),                                
     3           UCR1(15),UCI1(15),UCR2(15),UCI2(15),                           
     4           URISV(1500),OVE(40)                                            
      COMPLEX    COEFST,URISV,OVE,URITTT,EX,EXRI                                
      DATA DD20,DDM20/ 1.E20, 1.0E-20/                                          
      NXCH2=NX-NXMAX3                                                           
      KEXCB1=KEXTCB(1)                                                          
      IF(ND.EQ.NDFMES) GO TO 1811                                               
      IF(NX.EQ.NXIMAX) GO TO 1821                                               
      GO TO 2000                                                                
 1811 IF(NXCH2.LE.0) GO TO 1821                                                 
      IF(NXCH2.NE.5) GO TO 2000                                                 
      IF(KSTORE.EQ.0) GO TO 2000                                                
      IF(KUNDER.EQ.0) GO TO 1883                                                
      GO TO 1841                                                                
CCCCCC  ***************  RENORMALIZATION OF TEMPORAL WAVE  *************        
CCCCCC  ***************  FUNCTION AND/OR OVELAP INTEGRAL   *************        
 1821 IF(KUNDER.EQ.0) GO TO 2000                                                
      TEST=ABS(UCR2(1))+ABS(UCI2(1))                                            
      IF(TEST.LT.DD20) GO TO 2000                                               
      DO 1825 M1=1,MXROW                                                        
      UCR1(M1)=UCR1(M1)*DDM20                                                   
      UCI1(M1)=UCI1(M1)*DDM20                                                   
      UCR2(M1)=UCR2(M1)*DDM20                                                   
      UCI2(M1)=UCI2(M1)*DDM20                                                   
      DO 1825 NQ=1,4                                                            
      FPRER(NQ,M1)=FPRER(NQ,M1)*DDM20                                           
      FPREI(NQ,M1)=FPREI(NQ,M1)*DDM20                                           
 1825 CONTINUE                                                                  
      IF(ND.EQ.NDFMES) GO TO 1828                                               
      DO 1827 M1=1,MXROW                                                        
      UCR1M(M1)=UCR1M(M1)*DDM20                                                 
      UCI1M(M1)=UCI1M(M1)*DDM20                                                 
      DO 1827 NQ=1,4                                                            
      FPRERM(NQ,M1)=FPRERM(NQ,M1)*DDM20                                         
      FPREIM(NQ,M1)=FPREIM(NQ,M1)*DDM20                                         
 1827 CONTINUE                                                                  
 1828 IF(KINEX.EQ.1) GO TO 1841                                                 
      IF(NOVEMX.EQ.0) GO TO 2000                                                
      DO 1830 NOVE=1,NOVEMX                                                     
 1830 OVE(NOVE)=OVE(NOVE)*DDM20                                                 
      GO TO 2000                                                                
CCCCCC    ******    RENORMALIZATION OF STORED WAVE FUNCTIONS    ********        
 1841 IF(KSTORE.EQ.0) GO TO 2000                                                
      NX1MX=NX-KEXCB1                                                           
      TEST=1.0                                                                  
      IF(NX1MX.LT.1) GO TO 2000                                                 
      IF(NX.EQ.NXMAX2) GO TO 2000                                               
      DO 1875 NX1=1,NX1MX                                                       
      EX=URISV(NX1)                                                             
      EXR=EX                                                                    
      EXI=-TTI*EX                                                               
      IF( ABS(EXR).LT.TEST) EXR=0.0                                             
      IF( ABS(EXI).LT.TEST) EXI=0.0                                             
      URISV(NX1)=(EXR+TTI*EXI)*DDM20                                            
 1875 CONTINUE                                                                  
 1883 IF(KTLOUT(10).LT.2) GO TO 2000                                            
      IF(KTRL(9).NE.0) GO TO 1884                                               
      IF(KEXCOM(11).GT.3) GO TO 2000                                            
 1884 WRITE(6,1885) (URISV(N9),N9=10,MSTOT,10)                                  
 1885 FORMAT(5(E12.4,E12.4,1X))                                                 
 2000 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE ELCROS                                                         
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      COMMON/ELAS/ SCAMP(200,2)                                                 
      COMPLEX SCAMP                                                             
      COMMON F1(200),G1(200),FD1(200),GD1(200),EXSGR1(200),EXSGR2(200)          
      COMPLEX EXSGR1,EXSGR2                                                     
      COMMON XAMP(100),FCRI(100),CROSEL(100),SGRUTH(100),XRATIO(100),           
     1       P1(100),PL(3,100)                                                  
      COMPLEX XAMP,FCRI,ERI                                                     
      SQRT10=SQRT(10.0)                                                         
      LAP1MX=KEXCOM(14)                                                         
      IF(KEXTCB(3)+KEXTCB(4)+KEXCOM(12)+KEXCOM(13).NE.0) GO TO 510              
      KINM=2-KTRLCB(5)                                                          
      DO 500 KINEX=1,KINM                                                       
      CALL INEXCM(2)                                                            
      WNFAC=SQRT10/WN                                                           
      ETA=CE                                                                    
      DO 120 NA=1,NTHETA                                                        
      SN= SIN(0.5*TETARD (NA))                                                  
      SN=SN*SN                                                                  
      FLN=ETA*ALOG(SN)-2.  *SGMAZZ                                              
      FNO=ETA/(2.  *WN   *SN)                                                   
      ERI=-FNO*SQRT10*( COS(FLN)-TTI* SIN(FLN))                                 
      ER=ERI                                                                    
      EI=-TTI*ERI                                                               
      SGRUTH(NA)=ER*ER+EI*EI                                                    
      FCRI(NA)=ERI                                                              
      XAMP(NA)=ZERO                                                             
      PL(1,NA)=1.0                                                              
      PL(2,NA)=COS(TETARD(NA))                                                  
      P1(  NA)=COS(TETARD(NA))                                                  
  120 CONTINUE                                                                  
      DO 300 LP1=1,LAP1MX                                                       
      ETR=EXSGR1(LP1)                                                           
      ETI=-TTI*EXSGR1(LP1)                                                      
      IF(KINEX.EQ.2) ETR=EXSGR2(LP1)                                            
      IF(KINEX.EQ.2) ETI=-TTI*EXSGR2(LP1)                                       
      ERI=TTR*(ETR*ETR-ETI*ETI)+2.0*TTI*ETR*ETI                                 
      ERI=ERI*(2.*LP1-1.)*SCAMP(LP1,KINEX)*WNFAC                                
      DO 250 NA=1,NTHETA                                                        
      XAMP(NA)=XAMP(NA)+ERI*PL(1,NA)                                            
      PL(3,NA)=((2.0*LP1+1)*P1(NA)*PL(2,NA)-(LP1  )*PL(1,NA))/(LP1+1)           
      PL(1,NA)=PL(2,NA)                                                         
      PL(2,NA)=PL(3,NA)                                                         
  250 CONTINUE                                                                  
  300 CONTINUE                                                                  
      DO 350 NA=1,NTHETA                                                        
      ERI=XAMP(NA)+FCRI(NA)                                                     
      ER=ERI                                                                    
      EI=-TTI*ERI                                                               
      CROSEL(NA)=ER*ER+EI*EI                                                    
      XRATIO(NA)=CROSEL(NA)/SGRUTH(NA)                                          
  350 CONTINUE                                                                  
      M1=TMAS+0.1                                                               
      M2=PMAS+0.1                                                               
      I1=NZCA                                                                   
      I2=NZSA                                                                   
      IF(KINEX.EQ.2) I1=NZCB                                                    
      IF(KINEX.EQ.2) I2=NZSB                                                    
      WRITE(6,360) KINEX,M1,I1,M2,I2,ELAB                                       
  360 FORMAT(1H1,/ 17X,43HELASTIC SCATTERING CROSS SECTION IN CHANNEL,I3        
     1,/17X,14HWITH (A,Z) = (, I3, 1H, ,I3, 5H) + (, I3, 1H, ,I3, 1H),          
     2 8HAT ELAB=,F7.3,4H MEV,//)                                               
      WRITE(6,365)                                                              
  365 FORMAT(4X,2(5HTHETA,6X,5HSIGMA,6X,12HSIGMA/SRGUTH,6X)/)                   
      NHAF=NTHETA/2+MOD(NTHETA,2)                                               
      WRITE(6,370)(TETADG(N),CROSEL(N),XRATIO(N),TETADG(N+NHAF),                
     1      CROSEL(N+NHAF),XRATIO(N+NHAF),N=1,NHAF)                             
  370 FORMAT(2X,F7.2,2E14.4,5X,F7.2,2E14.4)                                     
      WRITE(6,390)                                                              
  390 FORMAT(///20X,21H  SCATTERING C-MATRIX)                                   
      WRITE(6,395) (SCAMP(L,KINEX),L=1,LAP1MX)                                  
  395 FORMAT(1X,5(1H(,F7.4,F7.4,2H) ))                                          
      IF(KTRLCB(5).EQ.1) GO TO 510                                              
  500 CONTINUE                                                                  
  510 RETURN                                                                    
      END                                                                       
      SUBROUTINE DWAMP                                                          
      COMMON/CRAC/FACLOG(500),IA,IB,IC,ID,IE,IF,RAC,L9(10),U9                   
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      COMMON HBAR(4000),OVE(40),ABAMP(200),AMPL(4),FASE(4),                     
     1       LLPRIN(50),MLPRIN(50)                                              
      COMPLEX HBAR,OVE,EX                                                       
      IF(KOUTCB(3).NE.0) WRITE(6,101)                                           
  101 FORMAT(1H1,20X,10(1H*),10X,21HDWAMP HAS BEEN CALLED,10X,10(1H*))          
      PAI=3.141592653589793                                                     
      TWOPAI=2.*PAI                                                             
      KPSUM=MOD(KPCA+KPSA+KPCB+KPSB,2)                                          
      LLCMMX=0                                                                  
      MLMM=0                                                                    
      DO 155 JLS=1,JLSMAX                                                       
      LL=LLTWR(JLS)/2                                                           
      MLMXP1=LL+1                                                               
      LLCMMX=LLCMMX+MLMXP1                                                      
      DO 150 MLM=1,MLMXP1                                                       
      MLMM=MLMM+1                                                               
      LLPRIN(MLMM)=LL                                                           
  150 MLPRIN(MLMM)=MLM-1                                                        
  155 CONTINUE                                                                  
      KINEX=2                                                                   
      CALL INEXCM(2)                                                            
      LEP1MX=KEXCOM(14)+LMXPAR                                                  
      NHBMAX=LEP1MX*LLCMMX                                                      
      IF(NHBMAX.GT.4000) WRITE(6,180) NHBMAX                                    
  180 FORMAT(///,20(1H*),7HNHBMAX=,I6,51H GT. 4000, INCREASE HBAR(4000)         
     1 IN DWAMP AND DWCROS)                                                     
      IF(NHBMAX.GT.4000) STOP                                                   
      DO 710 NHB=1,NHBMAX                                                       
  710 HBAR(NHB)=ZERO                                                            
      LAP1MI=KEXCOM(13)+1                                                       
      LAP1MX=KEXCOM(14)                                                         
      DO 800 LAP1=LAP1MI,LAP1MX                                                 
      LITR=LAP1-1                                                               
      IF(KEXCOM(12).NE.0.AND.MOD(LITR,2).NE.KPSUM) GO TO 800                    
      KEXCOM(11)=LITR                                                           
      LITW=2*LITR                                                               
      CALL LBSELC                                                               
      MXROWE=MXROW                                                              
      IF(MXROWE.EQ.0) GO TO 800                                                 
      READ(11) NOVEMX,(OVE(NOVE),NOVE=1,NOVEMX)                                 
      NOVE=0                                                                    
      DO 770 NLJE=1,MXROWE                                                      
      LETW=LLROW(NLJE)                                                          
      LETR=LETW/2                                                               
      HAT1=(LITW+1)*(LETW+1)                                                    
      HAT2=SQRT(HAT1)                                                           
      LLCUM=0                                                                   
      DO 760 JLS=1,JLSMAX                                                       
      LFOMTW=LLTWR(JLS)                                                         
      K1=LFOMTW-IABS(LITW-LETW)                                                 
      K2=LITW+LETW-LFOMTW                                                       
      K3=MIN0(K1,K2)                                                            
      IF(K3.GE.0) GO TO 719                                                     
      LLCUM=LLCUM+(LFOMTW/2)+1                                                  
      GO TO 760                                                                 
  719 NOVE=NOVE+1                                                               
      EX=OVE(NOVE)                                                              
      MLE1MX=LFOMTW/2+1                                                         
      DO 750 MLE1=1,MLE1MX                                                      
      LLCUM=LLCUM+1                                                             
      MLETW=2*(MLE1-1)                                                          
      IA=LITW                                                                   
      IB=LETW                                                                   
      IC=LFOMTW                                                                 
      ID=0                                                                      
      IE=MLETW                                                                  
      IF=MLETW                                                                  
      CALL CLEB                                                                 
      C4=RAC                                                                    
      IF(C4.EQ.0.  ) GO TO 750                                                  
      C5=1.                                                                     
      IF(KTRLCB(8).EQ.1) GO TO 733                                              
      IE=0                                                                      
      IF=0                                                                      
      CALL CLEB                                                                 
      C5=RAC                                                                    
      K1=IABS(LFOMTW-LITW-LETW)/4                                               
      S5=1-2*MOD(K1,2)                                                          
      HAT5=(LITW+1)*(LETW+1)/ FLOAT(LFOMTW+1)                                   
      HAT5= SQRT(HAT5)                                                          
      C5=S5*C5*HAT5                                                             
  733 MLETAB=IABS(MLETW)                                                        
      K1=(LETW-MLETAB)/2+1                                                      
      K2=(LETW+MLETAB)/2+1                                                      
      G1= EXP(0.5  *(FACLOG(K1)-FACLOG(K2)))                                    
      K1=(LETW+MLETW)/2                                                         
      S1=1-2*MOD(K1,2)                                                          
      F1=S1*C4*C5*G1*HAT2                                                       
      N1=LETR*LLCMMX+LLCUM                                                      
      HBAR(N1)=HBAR(N1)+F1*EX                                                   
  750 CONTINUE                                                                  
  760 CONTINUE                                                                  
  770 CONTINUE                                                                  
  800 CONTINUE                                                                  
      IF(KEXCOM(12).EQ.0) GO TO 816                                             
      LLCMX2=2*LLCMMX                                                           
      LE1MX=LEP1MX-1                                                            
      DO 815 LEP1=1,LEP1MX                                                      
      IF(MOD(LEP1,2).NE.0) GO TO 815                                            
      IF(LEP1.EQ.2) GO TO 805                                                   
      IF(LEP1.EQ.LE1MX) GO TO 805                                               
      N0BAS=(LEP1-4)*LLCMMX                                                     
      NMAX=4                                                                    
      GO TO 808                                                                 
  805 IF(LEP1.EQ.2) N0BAS=0                                                     
      IF(LEP1.EQ.LE1MX) N0BAS=(LEP1MX-3)*LLCMMX                                 
      NMAX=2                                                                    
  808 DO 814 M1=1,LLCMMX                                                        
      N0=N0BAS+M1                                                               
      N1=N0+LLCMMX*(NMAX-1)                                                     
      DO 811 N=1,NMAX                                                           
      FASE(N)=0.                                                                
      NT=N0+(N-1)*LLCMX2                                                        
      HREAL=HBAR(NT)                                                            
      HIMAG=-TTI*HBAR(NT)                                                       
      AMPL(N)= SQRT(HREAL**2+HIMAG**2)                                          
      IF(AMPL(N).EQ.0.0  ) GO TO 811                                            
      PHASE= ATAN2(HIMAG,HREAL)                                                 
      IF(PHASE.LT.0.0  ) PHASE=PHASE+TWOPAI                                     
      FASE(N)=PHASE                                                             
      IF(N.EQ.1) GO TO 811                                                      
  809 IF(FASE(N).GE.FASE(N-1)) GO TO 811                                        
      FASE(N)=FASE(N)+TWOPAI                                                    
      GO TO 809                                                                 
  811 CONTINUE                                                                  
      IF(NMAX.EQ.4) GO TO 812                                                   
      AMP1=(AMPL(1)+AMPL(2))/2.0                                                
      FAS1=(FASE(1)+FASE(2))/2.0                                                
      GO TO 813                                                                 
  812 AMP1=.0625  *(-(AMPL(1)+AMPL(4))+9.  *(AMPL(2)+AMPL(3)))                  
      FAS1=.0625  *(-(FASE(1)+FASE(4))+9.  *(FASE(2)+FASE(3)))                  
  813 HBAR(N1)=AMP1*( COS(FAS1)+TTI* SIN(FAS1))                                 
  814 CONTINUE                                                                  
  815 CONTINUE                                                                  
  816 IF(KOUTCB(3)-1) 900,823,817                                               
  817 WRITE(6,818)                                                              
  818 FORMAT(16H OVLAP4-818.HBAR)                                               
      DO 820  LE1=1,LEP1MX                                                      
      N1MIN=(LE1-1)*LLCMMX+1                                                    
      N1MAX=N1MIN+LLCMMX-1                                                      
      WRITE(6,819) (HBAR(N1),N1=N1MIN,N1MAX)                                    
  819 FORMAT(5(E11.3,E11.3,2X))                                                 
  820 CONTINUE                                                                  
  823 DO 850 LLC=1,LLCMMX                                                       
      DO 830 LE1=1,LEP1MX                                                       
      N1=(LE1-1)*LLCMMX+LLC                                                     
      EX=HBAR(N1)                                                               
      EXR=EX                                                                    
      EXI=-TTI*EX                                                               
      ABAMP(LE1)=EXR*EXR+EXI*EXI                                                
  830 CONTINUE                                                                  
      WRITE(6,835) LLPRIN(LLC),MLPRIN(LLC)                                      
  835 FORMAT(//17H DWAMP-835. ABAMP,5H  LL=,I3,5H  ML=,I3)                      
      WRITE(6,840) (ABAMP(LE1),LE1=1,LAP1MX)                                    
  840 FORMAT(1X,10E13.5)                                                        
  850 CONTINUE                                                                  
  900 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE DWCROS                                                         
      COMMON/DWBA/KTRLCB(10),KEXTCB(10),KOUTCB(10),JLSMAX,JJTWR(5),             
     1            LLTWR(5),ISTWR(5),ALSJF(5),NTHETA,TETADG(100),                
     2            TETARD(100),DTETAD,KINEX,MSTOT,AOVERB,WNI,WNE,ELABI,          
     3            QVALGR,NATOTL,NBTOTL,MESFCA,MESFCB,LBSTTL,LDWMAX,             
     4            KTAP12,LMXPAR,NBREAD,NBREDT,NBLMCM,NBSMCM,NTCOM,              
     5            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     6            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,XMESA,XMESB           
      COMMON/JP1/ KTRL(15),KEXCOM(15),KTLOUT(15),LLROW(15),PMAS,TMAS,           
     1            RMAS,CHARGE,CE,SGMAZZ,ECM,ELAB,WN,WNINI,DR,                   
     2            MXROW,NXMAX,NXPOT,NDFMES,XMES1,XMES2,VSX,WSX,WSF,             
     3            DFN,DFNW,DFNS,RZERO,RZEROW,RZEROS,RZEROC,TTR,TTI,ZERO         
      COMPLEX     TTR,TTI,ZERO                                                  
      COMMON/LEGENC/LCALTR,MMXTR,RADIAN,PLM20,PLM10,PL(15)                      
      COMMON/OUTP/ NOUT(20),IBSPR(10),BSPAR(20),OMPAR(20)                       
      COMMON HBAR(4000),XAMP(30,50),CROSEX(100,6),                              
     1       PM(30,20),PMST(2,30),LLOUT(7)                                      
      COMPLEX    HBAR,XAMP,EX,SUM,HBFAC                                         
      DATA NEFR,NNR,NZR,IINC,IEXIT / 3HEFR,3HNR ,3HZR ,3HIN ,3HOUT /            
      LLTWMX=0                                                                  
      LLCMMX=0                                                                  
      DO 155 JLS=1,JLSMAX                                                       
      LLTW=LLTWR(JLS)                                                           
      IF(LLTW.GT.LLTWMX) LLTWMX=LLTW                                            
      LLCMMX=LLCMMX+(LLTW/2)+1                                                  
      LLOUT(JLS)=LLTW/2                                                         
  155 CONTINUE                                                                  
      MMXTR1=LLTWMX/2+1                                                         
      MMXTR=MMXTR1-1                                                            
      LEP1MX=KEXCOM(14)+LMXPAR                                                  
      FOURPI=4.  *3.14159265                                                    
      SI1=(FOURPI/(WNI*WNE))**2                                                 
      SIMPSQ=SI1                                                                
      IF(KTRLCB(7).EQ.0.AND.KTRLCB(8).EQ.0) SIMPSQ=SIMPSQ/FOURPI                
      FAL=MASCA                                                                 
      FAS=MASSA                                                                 
      FBL=MASCB                                                                 
      FBS=MASSB                                                                 
      IF(KTRLCB(8).EQ.0) GO TO 703                                              
      FT=FAL+FAS                                                                
      FX=FAS-FBS                                                                
      SI2=(FAS*FBL/(FX*FT))**6                                                  
      IF(MASSA.LT.MASSB)   SI2=(FBS*FAL/(FX*FT))**6                             
      SIMPSQ=SI1*SI2                                                            
  703 RMASI=FAL*FAS/(FAL+FAS)                                                   
      RMASE=FBL*FBS/(FBL+FBS)                                                   
      SCALE1 = RMASI*RMASE*(WNE/WNI)*1.4489E-05                                 
      SCALE=SCALE1*       10.  *SIMPSQ/ FLOAT((ISATW+1)*(ICATW+1))              
      AOVBFC=FAL/FBL                                                            
      IF(MASCA.LT.MASCB)   AOVBFC=(1.0/AOVBFC)**2                               
      IF(MASSA.LT.MASSB)   SCALE3=(ISBTW+1)*(ICATW+1)                           
      IF(MASSA.GT.MASSB)   SCALE3=(ISATW+1)*(ICBTW+1)                           
      IF(KTRLCB(7).EQ.1) SCALE3=SCALE3*(AOVBFC**2)                              
  715 SCALE=SCALE*SCALE3                                                        
      NTHREP=(NTHETA-1)/30+1                                                    
      DO 900 NTHR=1,NTHREP                                                      
      NANGBS=30*(NTHR-1)                                                        
      N1=NTHETA-NANGBS                                                          
      NTHE=MIN0(30,N1)                                                          
      DO 805 NA=1,NTHE                                                          
      DO 805 LLC=1,LLCMMX                                                       
  805 XAMP(NA,LLC)=ZERO                                                         
      KEXC14=KEXCOM(14)                                                         
      DO 850 LEP1=1,KEXC14                                                      
      LETR=LEP1-1                                                               
      LCALTR=LETR                                                               
      IF(LETR-1) 808,811,813                                                    
  808 DO 809 IAN=1,NTHE                                                         
      PMST(1,IAN)=1.                                                            
  809 PM(IAN,1)=1.                                                              
      GO TO 825                                                                 
  811 DO 812 IAN=1,NTHE                                                         
      TH=TETARD(IAN+NANGBS)                                                     
      PM(IAN,1)= COS(TH)                                                        
      PMST(2,IAN)=PM(IAN,1)                                                     
  812 PM(IAN,2)= SIN(TH)                                                        
      GO TO 825                                                                 
  813 DO 820 IAN=1,NTHE                                                         
      RADIAN=TETARD(IAN+NANGBS)                                                 
      PLM20=PMST(1,IAN)                                                         
      PLM10=PMST(2,IAN)                                                         
      CALL LEGNDR                                                               
      DO 815 M=1,MMXTR1                                                         
  815 PM(IAN,M)=PL(M)                                                           
      PMST(1,IAN)=PMST(2,IAN)                                                   
      PMST(2,IAN)=PL(1)                                                         
  820 CONTINUE                                                                  
CCCCCC                                                                          
CCCCCC     **********     CALCULATION OF CROSS SECTION     **********           
CCCCCC                                                                          
  825 LLCUM=0                                                                   
      NHBBAS=(LEP1-1)*LLCMMX                                                    
      DO 840 JLS=1,JLSMAX                                                       
      LFOMTW=LLTWR(JLS)                                                         
      MLE1MX=LFOMTW/2+1                                                         
      DO 835 MLE1=1,MLE1MX                                                      
      LLCUM=LLCUM+1                                                             
      NHB=NHBBAS+LLCUM                                                          
      HBFAC=HBAR(NHB)                                                           
      IF(LEP1.LT.MLE1) GO TO 835                                                
      DO 830 NA=1,NTHE                                                          
      XAMP(NA,LLCUM)=XAMP(NA,LLCUM)+HBFAC*PM(NA,MLE1)                           
  830 CONTINUE                                                                  
  835 CONTINUE                                                                  
  840 CONTINUE                                                                  
  850 CONTINUE                                                                  
      DO 855 N=1,JLSMAX                                                         
      DO 855 NA=1,NTHE                                                          
      NB=NA+NANGBS                                                              
  855 CROSEX(NB,N)=0.0                                                          
      LLCUM=0                                                                   
      DO 870 JLS=1,JLSMAX                                                       
      MLE1MX=LLTWR(JLS)/2+1                                                     
      DO 860 MLE1=1,MLE1MX                                                      
      LLCUM=LLCUM+1                                                             
      DO 857 NA=1,NTHE                                                          
      EX=XAMP(NA,LLCUM)                                                         
      F1=EX* CONJG(EX)                                                          
      IF(MLE1.NE.1) F1=2.  *F1                                                  
      NB=NA+NANGBS                                                              
      CROSEX(NB,JLS)=CROSEX(NB,JLS)+F1                                          
  857 CONTINUE                                                                  
  860 CONTINUE                                                                  
  870 CONTINUE                                                                  
  900 CONTINUE                                                                  
      NCC1=JLSMAX+1                                                             
      KRALSJ=1                                                                  
      IF(KEXTCB(5).NE.0) KRALSJ=KEXTCB(5)                                       
      DO 990 KRA=1,KRALSJ                                                       
      IF(KRA.NE.1) GO TO 902                                                    
      SCALE=SCALE*ALSJF(1)                                                      
      GO TO 903                                                                 
  902 SCALE=SCALE*(ALSJF(KRA)/ALSJF(KRA-1))                                     
  903 DO 905 NA=1,NTHETA                                                        
  905 CROSEX(NA,NCC1)=0.0                                                       
      DO 950 JLS=1,JLSMAX                                                       
      DO 910 NA=1,NTHETA                                                        
      CROSEX(NA,JLS)=CROSEX(NA,JLS)*SCALE                                       
  910 CROSEX(NA,NCC1)=CROSEX(NA,NCC1)+CROSEX(NA,JLS)                            
  950 CONTINUE                                                                  
      NREFR=NZR                                                                 
      IF(KTRLCB(7).EQ.1) NREFR=NNR                                              
      IF(KTRLCB(8).EQ.1) NREFR=NEFR                                             
      WRITE(6,962) NREFR,(NOUT(N),N=1,8),ELABI                                  
  962 FORMAT(1H1,/1X,A3,9H DWBA FOR,2X,                                         
     1  A2,1H(,I3,2H)(A2,1H(I3,2H),A2,1H(I3,2H)),A2,1H(,I3,10H)-REACTION        
     2   ,9H AT ELAB=F7.3,4H-MEV)                                               
      WRITE(6,963)(NOUT(N),N=9,20),QVALGR                                       
  963 FORMAT (15X,I2,A2,A3,1H(,I2,A2,A3,1H,,I2,A2,A3,1H),I2,A2,A3,              
     1 10X,8HQ VALUE=,F7.3,4H MEV)                                              
      WRITE(6,964) (JJTWR(N),LLTWR(N),ISTWR(N),N=1,JLSMAX)                      
  964 FORMAT(/29H FORM FACTOR HAVE (2J,2L,2S)=,10(1H(,3I2,1H),1X))              
      LASTP=1                                                                   
      IF(KEXCOM(12).EQ.2) LASTP=2                                               
      XMES=XMESA/MESFCA                                                         
      WRITE(6,965)KEXCOM(14),KEXCOM(13),LASTP,NXMAX,KEXTCB(1),XMES              
  965 FORMAT( 9H LDW MAX=,I5,3X,12HLDW CUT OFF=,I5,3X,9HLDW STEP=,I5,           
     1       /9H MES MAX=,I5,3X,12HMES CUT OFF=,I5,3X,9HMES STEP=,F5.2)         
      IF(KTRLCB(8).EQ.1) WRITE(6,966)NATOTL,MESFCA,MESFCB                       
  966 FORMAT( 9H NATOTL =,I5,3X,12HMESFCA     =,I5,3X,9HMESFCB  =,I5)           
      WRITE(6,967)                                                              
  967 FORMAT(/1X,25HSYSTEM  (N,L,J)   B.E.   ,3HVSX,3X,4HVSOR,2X,4HDFNR,        
     12X,5HDFNSO, 2X,3HRZR,3X,4HRZSO,2X,3HRZC,3X,4HXMES)                        
      I=1                                                                       
      WRITE(6,968)I,(IBSPR(N),N=1,3),(BSPAR(N),N=1,9)                           
      I=2                                                                       
      WRITE(6,968)I,(IBSPR(N),N=4,6),(BSPAR(N),N=10,18)                         
  968 FORMAT(4X,I2,2X,1H(,3I2,1H),2F7.3,1X,7F6.3)                               
      WRITE(6,971)                                                              
  971 FORMAT(/52H CHANNEL  VSX    WSX    WSF    DFN  DFNW  DFNS  RZRO,          
     1  18H  RZRW  RZRS  RZRC)                                                  
      WRITE(6,973)IINC ,(OMPAR(N),N= 1,10)                                      
      WRITE(6,973)IEXIT,(OMPAR(N),N=11,20)                                      
  973 FORMAT(1X,A3,4X,3F7.3,8F6.3)                                              
      WRITE(6,977) ALSJF(KRA)                                                   
  977 FORMAT(15H ALSJ USED IS  ,8F10.5)                                         
      WRITE(6,982)(LLOUT(N),N=1,JLSMAX)                                         
  982 FORMAT(/2X,5HANGLE,5X,3HSUM,8X,2HL=,7(I1,8X,2HL=))                        
      DO 985 NA=1,NTHETA                                                        
      IF(MOD(NA-1,5).EQ.0.AND.NA.NE.1) WRITE(6,983)                             
  983 FORMAT( )                                                                 
      WRITE(6,984) TETADG(NA),CROSEX(NA,NCC1),(CROSEX(NA,N),N=1,JLSMAX)         
  984 FORMAT(1X,F6.2,5E11.3)                                                    
  985 CONTINUE                                                                  
  990 CONTINUE                                                                  
 1000 RETURN                                                                    
      END                                                                       
      SUBROUTINE LEGNDR                                                         
      COMMON/LEGENC/LCALTR,MMXTR,RADIAN,PLM20,PLM10,PL(15)                      
      FL=LCALTR                                                                 
      FLM1=FL-1                                                                 
      FLP1=FL+1                                                                 
      TWLM1=FL+FLM1                                                             
      TWLP1=FL+FLP1                                                             
      CO= COS(RADIAN)                                                           
      COSAB= ABS(CO)                                                            
      IF( ABS(COSAB-1.  ).LT.1.E-7) GO TO 221                                   
      SI=1.  / SIN(RADIAN)                                                      
      CT=2.  *CO*SI                                                             
      PL0  =(TWLM1*CO*PLM10-FLM1*PLM20)/FL                                      
      PLP10=(TWLP1*CO*PL0  -FL  *PLM10)/FLP1                                    
      PL(1)=PL0                                                                 
      PL(2)=FLP1*SI*(CO*PL0-PLP10)                                              
      MMAX=MIN0(LCALTR,MMXTR)+1                                                 
      IF(MMAX.LT.3) GO TO 200                                                   
      FM=0.                                                                     
      DO 150 M=3,MMAX                                                           
      FMP1=FM+1.                                                                
      PL(M)=CT*FMP1*PL(M-1)-(FL-FM)*(FL+FMP1)*PL(M-2)                           
      FM=FMP1                                                                   
  150 CONTINUE                                                                  
      GO TO 200                                                                 
  221 K1=MOD(LCALTR,2)                                                          
      PL(1)=CO**K1                                                              
      DO 225 M=2,MMAX                                                           
  225 PL(M)=0.                                                                  
  200 RETURN                                                                    
      END                                                                       
      SUBROUTINE CLEB                                                           
      COMMON / CRAC / FACLOG(500),IA,IB,IC,ID,IE,IF,RAC                         
      RAC=0.0                                                                   
      IF(ID+IE-IF) 1000,105,1000                                                
  105 K1=IA+IB+IC                                                               
      IF(K1-2*(K1/2)) 1000,110,1000                                             
  110 K1=IA+IB-IC                                                               
      K2=IC-IABS(IB-IA)                                                         
      K3=MIN0(K1,K2)                                                            
      IF(K3) 1000,130,130                                                       
  130 IF((-1)**(IB+IE)) 1000,1000,140                                           
  140 IF((-1)**(IC+IF)) 1000,1000,150                                           
  150 IF(IA-IABS(ID)) 1000,152,152                                              
  152 IF(IB-IABS(IE)) 1000,154,154                                              
  154 IF(IC-IABS(IF)) 1000,160,160                                              
  160 IF(IA) 1000,175,165                                                       
  165 IF(IB) 1000,175,170                                                       
  170 IF(IC) 1000,180,250                                                       
  175 RAC=1.0                                                                   
      GO TO 1000                                                                
  180 FB=IB+1                                                                   
      RAC=((-1.0  )**((IA-ID)/2))/ SQRT(FB)                                     
      GO TO 1000                                                                
  250 FC2=IC+1                                                                  
      IABCP=(IA+IB+IC)/2+1                                                      
      IABC=IABCP-IC                                                             
      ICAB=IABCP-IB                                                             
      IBCA=IABCP-IA                                                             
      IAPD=(IA+ID)/2+1                                                          
      IAMD=IAPD-ID                                                              
      IBPE=(IB+IE)/2+1                                                          
      IBME=IBPE-IE                                                              
      ICPF=(IC+IF)/2+1                                                          
      ICMF=ICPF-IF                                                              
      SQFCLG=0.5  *(ALOG(FC2)-FACLOG(IABCP+1)                                   
     1      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)                             
     2      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)                             
     3      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))                            
      NZMIC2=(IB-IC-ID)/2                                                       
      NZMIC3=(IA-IC+IE)/2                                                       
      NZMI=MAX0(0,NZMIC2,NZMIC3)+1                                              
      NZMX=MIN0(IABC,IAMD,IBPE)                                                 
      S1=1-2*MOD(NZMI-1,2)                                                      
      DO 400 NZ=NZMI,NZMX                                                       
      NZM1=NZ-1                                                                 
      NZT1=IABC-NZM1                                                            
      NZT2=IAMD-NZM1                                                            
      NZT3=IBPE-NZM1                                                            
      NZT4=NZ-NZMIC2                                                            
      NZT5=NZ-NZMIC3                                                            
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)                        
     1           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)                        
      SSTERM=S1* EXP(TERMLG)                                                    
      RAC=RAC+SSTERM                                                            
  400 S1=-S1                                                                    
      IF( ABS(RAC).LT.1.0E-10) RAC=0.0                                          
 1000 RETURN                                                                    
      END                                                                       
