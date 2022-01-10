C *******************************************************************
       PROGRAM FOLDER
C
C      Written by David P. Murdock
C      Last revision 1/15/90
C
C      Calculates relativistic optical potential from the RIA.
C      Uses either the RLF model of Horowitz (PRC 31,1340 (1985)) or
C      the parametrization of MRW (PRC 27,2123 (1983)).                 
C *******************************************************************
C
      COMMON/POTS/POT(500,2,2),POTDIR(500,2,2),POTEXC(500,2,2)
      COMMON/INPUT/IFOLD,IPICPL,IPBL,RHO(500,2,2)
      COMMON/CONST/PI,HBARC,PMASS,AMU
      COMMON/PB/PBCOEF(2,2)
      COMMON/RELLF/NMES,ILOR(32),IT(32),G2(2,32),G2O4P(2,32),
     1        AM(2,32),CUT(2,32),AW(2,32),AY(2,32),AZ(2,32)
      COMMON/KINEM/ELAB,PLAB,ECM,PCM,XKC,EKCFAC
      COMMON/GRIDS/DR,NR,RNUC,NRP,NX
C
      OPEN(UNIT=7,FILE='FOLDER.INP',STATUS='OLD')
      OPEN(UNIT=22,FILE='FOLDER.LOG',STATUS='UNKNOWN')
      OPEN(UNIT=11,FILE='TIMORA.DEN',STATUS='OLD')
      OPEN(UNIT=12,FILE='FOLDER.POT',STATUS='UNKNOWN')
C
      WRITE(6,'(''     PROGRAM FOLDER (1990)   '')')
      WRITE(22,'(''     PROGRAM FOLDER (1990)   '')')
      CALL SETUP
      IF(IFOLD.EQ.1) THEN
        CALL RLFSET
        CALL RLFFLD
        ENDIF
      IF(IFOLD.EQ.2) THEN
        CALL MRWSET
        CALL MRWFLD
        ENDIF
      CALL OUTPUT
C
      STOP
      END
C
C *******************************************************************
C
      BLOCK DATA
      COMMON/CONST/PI,HBARC,PMASS,AMU
      DATA PI,HBARC,PMASS,AMU/3.1415926535,197.328,939.,931.77/
      END
C
C ******************************************************************
C
      SUBROUTINE SETUP
      CHARACTER*72 ADUM
      COMMON/INPUT/IFOLD,IPICPL,IPBL,RHO(500,2,2)
      COMMON/CONST/PI,HBARC,PMASS,AMU
      COMMON/KINEM/ELAB,PLAB,ECM,PCM,XKC,EKCFAC
      COMMON/GRIDS/DR,NR,RNUC,NRP,NX
      COMMON/FALL/DAMP(2,2)
C
C
      READ(7,5) ADUM
5      FORMAT(A72)
      READ(7,*) IFOLD
      READ(7,5) ADUM
      READ(7,*) IPICPL,IPBL
      READ(7,5) ADUM
      READ(7,*) ELAB,AM2
      READ(7,5) ADUM
      READ(7,*) DR,NR,NRP,NX
C
      WRITE(6,34)
      WRITE(22,34)
      WRITE(6,35) ELAB,AM2
      WRITE(22,35) ELAB,AM2
      WRITE(6,36) DR,NR,NRP,NX
      WRITE(22,36) DR,NR,NRP,NX
      WRITE(6,37) IPICPL,IPBL
      WRITE(22,37) IPICPL,IPBL
34    FORMAT(/' Input:')
35    FORMAT(' Elab=',f7.2,'    target mass =',f7.2)
36    FORMAT(' dr=',f7.4,'   grid pts=',i4,
     1'   nodes (r)=',i3,'   nodes (x)=',i3)
37    FORMAT(' pion coupling=',i3,'    Pauli blocking? ',i3)
C
      READ(11,5) ADUM
      DO 70 IR=1,NR
70      READ(11,*) R,RHO(IR,1,2),RHO(IR,2,2),RHO(IR,1,1),
     1        RHO(IR,2,1)
C
      DO 80 IPN=1,2
      DO 80 ISV=1,2
      DAMP(IPN,ISV)=RHO(NR,IPN,ISV)/RHO(NR-1,IPN,ISV)
      DO 80 IR=NR+1,500
80      RHO(IR,IPN,ISV)=RHO(IR-1,IPN,ISV)*DAMP(IPN,ISV)
C
      RNUC=1.2*AM2**.333333
      ENLAB=ELAB+PMASS
      PLAB=SQRT((ELAB+PMASS)**2-PMASS**2)
      XKC=.5*SQRT(PLAB**2-ELAB**2)
      EKC=SQRT(XKC*XKC+PMASS*PMASS)
      EKCFAC=PMASS**2/EKC
      FM1=PMASS
      FM2=AM2*AMU
      FMT=FM1+FM2
      T1=SQRT(2.*ELAB*FM2+FMT*FMT)
      ECM=T1-FMT
      FK2=(ELAB**2+2.*ELAB*FM1)*(FM2/T1)**2
      PCM=SQRT(FK2)
      ENCM=SQRT(PMASS*PMASS+FK2)
      WRITE(22,82) elab,enlab,plab,encm,pcm,ekc,xkc
82    FORMAT(/' Proton energies and momenta (MeV) for Elab=',f7.2,':',
     *  /,'    (Proton energies include rest mass.)',
     *  /,17x,'E',10x,'p',/,
     *' Lab:',7x,f9.3,3x,f9.3,/,' N-Nucl CM',2x,f9.3,3x,f9.3,/,
     *' N-N CM:',4x,f9.3,3x,f9.3,/,)
C
      RETURN
      END
C
C ****************************************************************
C
      SUBROUTINE RLFSET
      COMMON/INPUT/IFOLD,IPICPL,IPBL,RHO(500,2,2)
      COMMON/KINEM/ELAB,PLAB,ECM,PCM,XKC,EKCFAC
      COMMON/CONST/PI,HBARC,PMASS,AMU
      COMMON/PB/PBCOEF(2,2)
      COMMON/RELLF/NMES,ILOR(32),IT(32),G2(2,32),G2O4P(2,32),
     1        AM(2,32),CUT(2,32),AW(2,32),AY(2,32),AZ(2,32)
      DIMENSION ILTAB(5,32),ITTAB(5,32),G2TAB(5,2,32),
     1        AMTAB(5,2,32),CUTTAB(5,2,32),PBTAB(5,2,2),ETAB(5),
     1        NMESON(5)
C
      WRITE(6,31)
      WRITE(22,31)
31    FORMAT(/' Folding with RLF model'/)
      WRITE(22,32)
      WRITE(6,32)
32    FORMAT(' Reading in RLF model parameters')
C
      OPEN(UNIT=19,FILE='RLF.DAT',STATUS='OLD')
      READ(19,*) NSET
      READ(19,*) (ETAB(I),I=1,NSET)
C
      DO 50 K=1,NSET
      READ(19,*) NMESON(K)
      READ(19,*) (ILTAB(K,I),I=1,NMESON(K))
      READ(19,*) (ITTAB(K,I),I=1,NMESON(K))
      DO 35 IRI=1,2
      READ(19,*) (G2TAB(K,IRI,I),I=1,NMESON(K))
      READ(19,*) (AMTAB(K,IRI,I),I=1,NMESON(K))
35      READ(19,*) (CUTTAB(K,IRI,I),I=1,NMESON(K))
      READ(19,*) PBTAB(K,1,1),PBTAB(K,1,2),PBTAB(K,2,1),
     1        PBTAB(K,2,2)
50      CONTINUE
C
      KSET=1
      ABSDIF=100.
      DO 100 ISET=1,NSET
      TDIF=ABS(ELAB-ETAB(ISET))
      IF(TDIF.LT.ABSDIF) THEN
        ABSDIF=TDIF
        KSET=ISET
        ENDIF
100      CONTINUE
      K1=KSET-1
      K2=KSET
      IF(ELAB.GT.ETAB(KSET) .OR. KSET.EQ.1) THEN
        K1=KSET
        K2=KSET+1
        ENDIF
C
      WRITE(22,105) ETAB(KSET)
105   FORMAT(' Using parameters for ',f6.2,' MeV:')
C
      NMES=NMESON(KSET)
      DO 80 IMES=1,NMES
      IT(IMES)=ITTAB(KSET,IMES)
      ILOR(IMES)=ILTAB(KSET,IMES)
      DO 80 IRI=1,2
      G2(IRI,IMES)=G2TAB(KSET,IRI,IMES)
      AM(IRI,IMES)=AMTAB(KSET,IRI,IMES)
      CUT(IRI,IMES)=CUTTAB(KSET,IRI,IMES)
      G2O4P(IRI,IMES)=G2(IRI,IMES)/4./PI
      AW(IRI,IMES)=CUT(IRI,IMES)**2/(CUT(IRI,IMES)**2
     1        -AM(IRI,IMES)**2)
      AY(IRI,IMES)=AM(IRI,IMES)**2/(CUT(IRI,IMES)**2
     1        -AM(IRI,IMES)**2)
      AZ(IRI,IMES)=CUT(IRI,IMES)**2/4./PMASS**2
      CUT(IRI,IMES)=CUT(IRI,IMES)/HBARC
      AM(IRI,IMES)=AM(IRI,IMES)/HBARC
80      CONTINUE
C
      WRITE(22,*) (ILOR(I),I=1,NMES)
      WRITE(22,*) (IT(I),I=1,NMES)
      do 85 iri=1,2
      WRITE(22,*) (G2(IRI,I),I=1,NMES)
      WRITE(22,*) (HBARC*AM(IRI,I),I=1,NMES)
85      WRITE(22,*) (HBARC*CUT(IRI,I),I=1,NMES)
C
      IF(IPBL.EQ.0) GO TO 160
      DO 150 ISV=1,2
      DO 150 IRI=1,2
150      PBCOEF(ISV,IRI)=PBTAB(K1,ISV,IRI)+(ELAB-ETAB(K1))*
     1        (PBTAB(K2,ISV,IRI)-PBTAB(K1,ISV,IRI))/
     1        (ETAB(K2)-ETAB(K1))
C
      WRITE(22,155) PBCOEF
155      FORMAT(' Use Pauli blocking, with coefficients:',
     1/,4f10.5)
C
160   RETURN
      END
C
C *******************************************************************
C
      SUBROUTINE RLFFLD
      COMMON/POTS/POT(500,2,2),POTDIR(500,2,2),POTEXC(500,2,2)
      COMMON/INPUT/IFOLD,IPICPL,IPBL,RHO(500,2,2)
      COMMON/CONST/PI,HBARC,PMASS,AMU
      COMMON/PB/PBCOEF(2,2)
      COMMON/KINEM/ELAB,PLAB,ECM,PCM,XKC,EKCFAC
      COMMON/GRIDS/DR,NR,RNUC,NRP,NX
      DIMENSION XRP(80),WRP(80),XX(80),WX(80)
C
      RMAX=2.5*RNUC
      CALL GAULEG(0.,RMAX,XRP,WRP,NRP)
      CALL GAULEG(-1.,1.,XX,WX,NX)
      RIAFAC=2.*PI*4.*PI*PCM*HBARC/939.
C
      write(6,98)
98    format('        R(FM)      SR(MEV)        SI       VR(MEV)',
     *'      VI')
      DO 1000 IR=1,NR
      R=IR*DR
      DO 1000 ISV=1,2
      DO 1000 IRI=1,2
      RISIGN=1.
      IF(IRI.EQ.2) RISIGN=-1.
C
      RPSUM=0.
      DO 500 IRP=1,NRP
      RP=XRP(IRP)
      CALL FDIR(RP,ISV,IRI,FPPD,FPND)
      XPSUM=0.
      XNSUM=0.
      DO 400 IX=1,NX
      X=XX(IX)
      RPLUS=SQRT(R*R+RP*RP+2.*R*RP*X)
      CALL RHOINT(ISV,RPLUS,RHOPD,RHOND)
      XPSUM=XPSUM+RHOPD*WX(IX)
      XNSUM=XNSUM+RHOND*WX(IX)
400      CONTINUE
500      RPSUM=RPSUM+(XPSUM*FPPD+XNSUM*FPND)*RP*RP*WRP(IRP)
      POTDIR(IR,ISV,IRI)=RPSUM*RIAFAC*RISIGN
C
      RPSUM=0.
      DO 700 IRP=1,NRP
      RP=XRP(IRP)
      RPK=RP*PCM/HBARC
      CALL FEXC(RP,ISV,IRI,FPPE,FPNE)
      XPSUM=0.
      XNSUM=0.
      DO 600 IX=1,NX
      X=XX(IX)
      RRP2=SQRT(RP*RP/4.+R*R-R*RP*X)
      CALL RHOINT(2,RRP2,RHOP,RHON)
      RPKF=RP*(1.5*PI*PI*(RHOP+RHON))**.33333333
      DENFAC=3.*SBESJ1(RPKF)/RPKF
      CALL RHOINT(ISV,RRP2,RHOPE,RHONE)
      XPSUM=XPSUM+RHOPE*WX(IX)*DENFAC
      XNSUM=XNSUM+RHONE*WX(IX)*DENFAC
600      CONTINUE
700      RPSUM=RPSUM+(XPSUM*FPPE+XNSUM*FPNE)*SBESJ0(RPK)
     1        *RP*RP*WRP(IRP)
      POTEXC(IR,ISV,IRI)=RPSUM*RIAFAC*RISIGN
C
      IF(IPBL.GE.1) THEN
        RHOTOT=RHO(IR,1,2)+RHO(IR,2,2)
        PBFACT=1.-PBCOEF(ISV,IRI)*(RHOTOT/0.1934)**.6666666
        POTDIR(IR,ISV,IRI)=PBFACT*POTDIR(IR,ISV,IRI)
        POTEXC(IR,ISV,IRI)=PBFACT*POTEXC(IR,ISV,IRI)
      ENDIF
      POT(IR,ISV,IRI)=POTDIR(IR,ISV,IRI)+POTEXC(IR,ISV,IRI)
      IF(20*NINT(IR/20.).EQ.IR .and. isv+iri.eq.4)
     1  WRITE(6,99)R,pot(ir,1,1),
     1  pot(ir,1,2),pot(ir,2,1),pot(ir,2,2)
99        FORMAT(1x,5f12.4)
1000      CONTINUE
C
      RETURN
      END
C
C ******************************************************************
C
      SUBROUTINE MRWSET
      COMMON/KINEM/ELAB,PLAB,ECM,PCM,XKC,EKCFAC
      COMMON/GRIDS/DR,NR,RNUC,NRP,NX
      COMMON/MRW/F0(2,2,2),BETA(2,2,2)
      DIMENSION F0TAB(10,2,2,2),BETAB(10,2,2,2),ETAB(10)
C
      WRITE(6,40)
      WRITE(22,40)
40    FORMAT(/' Folding with MRW parametrization',)
      IF(NRP+NX.LT.100) WRITE(6,42)
42    FORMAT(/' ==>Warning!<==  NRP and NX may be too small.',/,
     1' I suggest NRP and NX both >= 50.'/)
      WRITE(22,45)
45      FORMAT(' Reading in MRW parameters.')
      OPEN(UNIT=18,FILE='MRW.DAT',STATUS='OLD')
      DO 50 IE=1,10
      ETAB(IE)=100.*IE
      READ(18,*) F0TAB(IE,1,1,1),F0TAB(IE,1,1,2),
     1        F0TAB(IE,1,2,1),F0TAB(IE,1,2,2),F0TAB(IE,2,1,1),
     1        F0TAB(IE,2,1,2),F0TAB(IE,2,2,1),F0TAB(IE,2,2,2)
      READ(18,*) BETAB(IE,1,1,1),BETAB(IE,1,1,2),
     1        BETAB(IE,1,2,1),BETAB(IE,1,2,2),BETAB(IE,2,1,1),
     1        BETAB(IE,2,1,2),BETAB(IE,2,2,1),BETAB(IE,2,2,2)
50      CONTINUE
C
      KSET=1
      ABSDIF=100.
      DO 100 ISET=1,10
      TDIF=ABS(ELAB-ETAB(ISET))
      IF(TDIF.LT.ABSDIF) THEN
        ABSDIF=TDIF
        KSET=ISET
        ENDIF
100      CONTINUE
C
      WRITE(22,105) ETAB(KSET)
105      FORMAT(' Using parameters for', f10.5,' MeV')
C
      DO 150 IPN=1,2
      DO 150 ISV=1,2
      DO 150 IRI=1,2
      F0(IPN,ISV,IRI)=F0TAB(KSET,IPN,ISV,IRI)*1.E-06
      BETA(IPN,ISV,IRI)=BETAB(KSET,IPN,ISV,IRI)*1.E-06
150      CONTINUE
C
      WRITE(22,155)
155      FORMAT(' Using parameters:')
      WRITE(22,*) F0,BETA
      RETURN
      END
C
C *****************************************************************
C
      SUBROUTINE MRWFLD
      COMPLEX AISUM,AJSUM,F0PP,F0PN,BETAPP,BETAPN,T1,T2,
     1        RIAFAC
      COMMON/POTS/POT(500,2,2),POTDIR(500,2,2),POTEXC(500,2,2)
      COMMON/INPUT/IFOLD,IPICPL,IPBL,RHO(500,2,2)
      COMMON/CONST/PI,HBARC,PMASS,AMU
      COMMON/KINEM/ELAB,PLAB,ECM,PCM,XKC,EKCFAC
      COMMON/MRW/F0(2,2,2),BETA(2,2,2)
      COMMON/GRIDS/DR,NR,RNUC,NRP,NX
      DIMENSION XRP(80),WRP(80),XQ(80),WQ(80)
C
      RMAX=NR*DR*1.2
      QMAX=10.
      NQ=NX
      CALL GAULEG(0.,RMAX,XRP,WRP,NRP)
      CALL GAULEG(0.,QMAX,XQ,WQ,NQ)
C
      RIAFAC=-4.*PI*PCM*CMPLX(0.,1.)/PMASS*(2./PI)
     1        *HBARC**3
      write(6,98)
98    format('        R(FM)      SR(MEV)        SI       VR(MEV)',
     *'      VI')
      DO 1000 IR=1,NR
      R=IR*DR
      DO 800 ISV=1,2
      F0PP=CMPLX(F0(1,ISV,1),F0(1,ISV,2))
      F0PN=CMPLX(F0(2,ISV,1),F0(2,ISV,2))
      BETAPP=CMPLX(BETA(1,ISV,1),BETA(1,ISV,2))*HBARC**2
      BETAPN=CMPLX(BETA(2,ISV,1),BETA(2,ISV,2))*HBARC**2
      AJSUM=CMPLX(0.,0.)
      DO 300 J=1,NQ
      AISUM=CMPLX(0.,0.)
      DO 100 I=1,NRP
      CALL RHOINT(ISV,XRP(I),RHOP,RHON)
      T1=RHOP*F0PP*EXP(-BETAPP*XQ(J)*XQ(J))
      T2=RHON*F0PN*EXP(-BETAPN*XQ(J)*XQ(J))
100      AISUM=AISUM+SIN(XQ(J)*XRP(I))*WRP(I)*(T1+T2)*XRP(I)
300      AJSUM=AJSUM+SIN(XQ(J)*R)*WQ(J)*AISUM
      AJSUM=RIAFAC*AJSUM/R
      POT(IR,ISV,1)=REAL(AJSUM)
      POT(IR,ISV,2)=AIMAG(AJSUM)
800      CONTINUE
      IF(20*NINT(IR/20.).EQ.IR) WRITE(6,99)R,pot(ir,1,1),
     1  pot(ir,1,2),pot(ir,2,1),pot(ir,2,2)
99        FORMAT(1x,5f12.4)
1000      CONTINUE
C
      RETURN
      END
C
C ******************************************************************
C
      SUBROUTINE RHOINT(ISV,R,RHOP,RHON)
C
      COMMON/INPUT/IFOLD,IPICPL,IPBL,RHO(500,2,2)
      COMMON/GRIDS/DR,NR,RNUC,NRP,NX
      COMMON/FALL/DAMP(2,2)
C
      IR=R/DR
      IF(IR.GT.499) THEN
        RHOP=DAMP(1,ISV)*(IR-499)*RHO(499,1,ISV)
        RHON=DAMP(2,ISV)*(IR-499)*RHO(499,2,ISV)
        RETURN
      ENDIF
      IF(IR.LT.1) THEN
        RHOP=RHO(1,1,ISV)
        RHON=RHO(1,2,ISV)
        RETURN
      ENDIF
      IR1=IR+1
      RHOP=RHO(IR,1,ISV)+(R-IR*DR)
     1        *(RHO(IR1,1,ISV)-RHO(IR,1,ISV))/DR
      RHON=RHO(IR,2,ISV)+(R-IR*DR)
     1        *(RHO(IR1,2,ISV)-RHO(IR,2,ISV))/DR
      RETURN
      END
C
C ******************************************************************
C
      SUBROUTINE GAULEG(X1,X2,X,W,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL X1,X2,X(N),W(N)
      PARAMETER (EPS=3.D-14)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
C
C ******************************************************************
C
      SUBROUTINE FDIR(R,ISV,IRI,FPPD,FPND)
C
      COMMON/INPUT/IFOLD,IPICPL,IPBL,RHO(500,2,2)
      COMMON/CONST/PI,HBARC,PMASS,AMU
      COMMON/RELLF/NMES,ILOR(32),IT(32),G2(2,32),G2O4P(2,32),
     1        AM(2,32),CUT(2,32),AW(2,32),AY(2,32),AZ(2,32)
      COMMON/KINEM/ELAB,PLAB,ECM,PCM,XKC,EKCFAC
C
      SUMPP=0.
      DO 100 IMES=1,NMES
      IF(ILOR(IMES).NE.ISV) GO TO 100
      SUMPP=SUMPP+FFUN(G2O4P(IRI,IMES),AW(IRI,IMES),
     1        AM(IRI,IMES),CUT(IRI,IMES),R)
100      CONTINUE
      FPPD=SUMPP*EKCFAC*.5/XKC
C
      SUMPN=0.
      DO 200 IMES=1,NMES
      IF(ILOR(IMES).NE.ISV) GO TO 200
      FAC=1.-4.*IT(IMES)
      SUMPN=SUMPN+FFUN(G2O4P(IRI,IMES),AW(IRI,IMES),
     1        AM(IRI,IMES),CUT(IRI,IMES),R)*FAC
C
200      CONTINUE
      FPND=.5*(SUMPN*EKCFAC*.5/XKC+FPPD)
      RETURN
      END
C
C *******************************************************************
C
      SUBROUTINE FEXC(R,ISV,IRI,FPPE,FPNE)
C
      COMMON/INPUT/IFOLD,IPICPL,IPBL,RHO(500,2,2)
      COMMON/CONST/PI,HBARC,PMASS,AMU
      COMMON/RELLF/NMES,ILOR(32),IT(32),G2(2,32),G2O4P(2,32),
     1        AM(2,32),CUT(2,32),AW(2,32),AY(2,32),AZ(2,32)
      COMMON/KINEM/ELAB,PLAB,ECM,PCM,XKC,EKCFAC
C
      DIMENSION F(30),C(5,5)
      DATA C/.25,.25,.125,-.25,.25,
     1             1.,-.5,0.,-.5,-1.,
     2             3.,0.,-.5,0.,3.,
     3             -1.,-.5,0.,-.3,1.,
     4             .25,-.25,.125,.25,.25/
C
      DO 50 IMES=1,NMES
      IF(IPICPL.EQ.1 .OR. IRI.EQ.2 .OR.
     1        ILOR(IMES).NE.5) THEN
      F(IMES)=FFUN(G2O4P(IRI,IMES),AW(IRI,IMES),
     1        AM(IRI,IMES),CUT(IRI,IMES),R)
      GO TO 50
      ENDIF
C
      FHOLD=FFUNPV(G2O4P(IRI,IMES),AW(IRI,IMES),AY(IRI,IMES),
     1        AZ(IRI,IMES),AM(IRI,IMES),CUT(IRI,IMES),R)
      IF(ISV.EQ. 1) THEN
      F(IMES)=-FHOLD
      GO TO 50
      ENDIF
      IF(ISV.EQ.2) THEN
      F(IMES)=FHOLD
      GO TO 50
      ENDIF
50      CONTINUE
C
      SUMPP=0.
      XSIGN=-1.
      FAC=1.
      DO 100 IMES=1,NMES
100      SUMPP=SUMPP+FAC*XSIGN*C(ISV,ILOR(IMES))*F(IMES)
      FPPE=SUMPP*EKCFAC*.5/XKC
C
      SUMPN=0.
      XSIGN=1.
      DO 200 IMES=1,NMES
      FAC=1.-4.*IT(IMES)
200      SUMPN=SUMPN+FAC*XSIGN*C(ISV,ILOR(IMES))*F(IMES)
      FPNE=.5*(SUMPN*EKCFAC*.5/XKC+FPPE)
      RETURN
      END
C
C ******************************************************************
C
      FUNCTION FFUN(G2O4P,AW,AM,CUT,R)
      T1=EXPOK(CUT*R)
      T=AW*(EXPOK(AM*R)-T1)/R
      T=T-CUT*T1*.5
      FFUN=G2O4P*AW*T
      RETURN
      END
C
      FUNCTION FFUNPV(G2O4P,AW,AY,AZ,AM,CUT,R)
      T1=EXPOK(CUT*R)
      T=AY/R*(T1-EXPOK(AM*R))
      T=T+CUT*T1*.5
      FFUNPV=AZ*G2O4P*AW*T
      RETURN
      END
C
C *****************************************************************
C
      FUNCTION SBESJ0(X)
      SBESJ0=SIN(X)/X
      RETURN
      END
C
      FUNCTION SBESJ1(X)
      SBESJ1=(SIN(X)-X*COS(X))/X**2
      RETURN
      END
C
C *****************************************************************
C
      SUBROUTINE OUTPUT
      COMMON/POTS/POT(500,2,2),POTDIR(500,2,2),POTEXC(500,2,2)
      COMMON/GRIDS/DR,NR,RNUC,NRP,NX
C
      DO 100 IR=1,NR
      R=IR*DR
100      WRITE(12,10) R,POT(IR,1,1),POT(IR,1,2),POT(IR,2,1),
     1        POT(IR,2,2)
10      FORMAT(' ',F8.3,3X,6(E11.4,3X))
      WRITE(6,120)
      WRITE(22,120)
120   FORMAT(/,' Calculation complete.',/,
     1' Output written to file FOLDER.POT')
      RETURN
      END
C
C ******************************************************************
C
      FUNCTION EXPOK(X)
      IF(X.GT.72.)X=72.
      EXPOK=EXP(-X)
      RETURN
      END
C
C ******************************************************************
C ******************************************************************
