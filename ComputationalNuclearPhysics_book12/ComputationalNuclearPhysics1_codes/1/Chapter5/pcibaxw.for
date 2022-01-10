	program pcibXW
C		X : eXcitation energies
C		W : Wave functions
C	shortened version of PHINT
C
C      WRITTEN BY :  OLAF SCHOLTEN
C      ADDRESS :    Kernfysisch Versneller Instituut
C				Zernikelaan 25
C                        9747 AA Groningen
C                        NEDERLAND
C
      INTEGER EDM,EVOD
      LOGICAL PRINT,FIT,PRINTV,MULT,WRTAPE,SDEQSF,PRINTP
     &    ,NPLOG
      CHARACTER COMMNT*40, P*1
      DIMENSION CLN(3),CLP(3),CONP(6),COEF(6),COEL(8),COMUL(6)
     &          ,ENERGY(1),FPARM(23)
      DIMENSION XYZ(34567)
C		Array XYZ is the workingspace for calculation of hamiltonian and 
C			eigenvectors. In the progran it will be dynamically be 
C			assigned, according to the dimension of the problem and
C			number of eigenvectors to be calculated. When the problem 
C			requires too much memory space, the program prints an 
C			error message and goes on to the calculation for the next
C			spin and parity.
C     COMMON / HAM0   / PHEN(15, 8),GAM
      COMMON / HAM0   / PHEN(8,4),GAM
C		Coefficients to facilitate calc. of diag. matrix el.
      COMMON / ECSBK  / IECSBK(50,2,2)
C		Bookkeeping matrix for eigen-values and -vectors
C
      COMMON / TEXT   / PRINT,PRINTV,PRINTP,IWD
      COMMON / TEXT2  / P,COMMNT
C   PRINT     : =.T. : PRINT ALSO INITIAL MATRIX AND EIG. VECTORS
C   PRINTV    : =.F. : DO NOT PRINT INITIAL MATRIX
      COMMON / CONTR  / NPHMSU,NPHMAX,NEIG,IAI,IAM,IPPI,IPPM,EVOD
C		Specifies particular case to be calculated
C   NPHMSU    : NUMBER OF PH. INCLUDED IN THE CALC.
C   IAI       : INITIAL L VALUE TO BE COMP.
C   IAM       : FINAL L VALUE TO BE COMP.
C   IPPI      : =1 : COMP. POS. PAR. ; =2 : NO POS. PAR.
C   IPPM      : =1 : NO NEG. PAR. ; =2 : COMP. NEG. PAR.
C   NEIG      : MAX NUMBER OF EIGENVECTORS TO BE COMPUTED
C
      COMMON / FPAR   / HBAR3P,DP(5),F3P,EPSDP,D(5),F3,FELL,FQQ,FEX,
     &                  HBAR3,EPSD,RKAP3,CHO,CHON,CHOP,SDEQSF
C		Parameters for f-boson coupling
C   SDEQSF    : =.T. : TAKE S PHONON FOR D EQUAL TO S PHONON FOR F PH
C
      COMMON / H      / HBAR,C(3),F,G,CH1,CH2
C		Standard SU(5) form for H
      COMMON / MUL    / EPS,PAIR,ELL,QQ,OCT,HEX,MULT,CHQ
C		Multipole form for H (if requested)
      COMMON / NP     / ED,RKAP,CHN,CLN,CHP,CLP,NN,NP,NPLOG
C		IBA-2 form for H, if requested
C     COMMON / REDMAT / REDGS(375,5),BETFAC(56,2),IBDEL(2,5)
      COMMON / REDMAT / REDGS(57,5),BETFAC(16,2),IBDEL(2,3)
C		Values of C.F.P.
      EQUIVALENCE (ENERGY,XYZ)
      EQUIVALENCE (HBAR,COEL(1)),(EPS,COMUL(1))
      EQUIVALENCE (ED,CONP(1)),(HBAR3P,FPARM(1))
C
C     DATA NDM,IDIB,IDB,IDGS/14,5,56,375/
C		info on coefficients of fractional parentage
      DATA NDM,IDIB,IDB,IDGS/7,3,16,57/
C
C      DIMENSIONS FOR ARRAYS IN COMMON REDMAT
C      USED FOR THE COMPUTATION OF REDUCED MAT. ELEMENTS
C      OF D+
C      NDM=MAX NUMBER OF PHONONS USED IN THE PROGRAM
C         IDIB=NDM/3+1
C         IDB=(K+1)*(NDM-K)  ; K=(NDM+1)/2
C         IDGS=(K+1)*(NDM**2-NDM+1+3*K*(1-NDM+K))   ;  K=NDM/3
C      THESE DIMENSIONS MUST BE EXACTLY EQUAL TO THOSE IN
C       COMMON/REDMAT , USED IN MAIN PROG AND SUBROUTINE RED
C      DIMENSION OF PHEN = (NDM+1,NDM/2+1)
C      PHEN IS PLACED IN COMMON/HAM0/ USED IN MAIN PROG AND SUBR HAM023
C      NO TEST IS MADE ON THE DIMENSION OF PHEN
C
      DATA LENMIN,LENMAX/1,34567/
C		Dimension of working-space array XYZ
C
      OPEN(UNIT=2,FILE='PHINT.OUT')
      OPEN(UNIT=3,STATUS='OLD',FILE='PHINT.CFP',FORM='UNFORMATTED',
     & RECL=2048)
      OPEN(UNIT=10,STATUS='UNKNOWN',FILE='PHINT.WAV',FORM='UNFORMATTED')
C
      WRITE(2,200) 
      WRITE(*,200) 
  200 FORMAT(' --- Program PCIBAXW ,version JANUARY 1990 ---')
C
	CALL RDXW
      LEN=LENMAX-LENMIN
      IEDMA=SQRT(81.+LEN)-9
C
      NPLOG=NP.GT.0 .AND. NN.GT.0
      IF(NPLOG) NPHMAX=NN+NP
      K=IEDMA/2
      IF(NEIG.GT.K) NEIG=K
      IF(NPHMAX.LE.0) NPHMAX=NDM
      IF(NPHMSU.LE.0.OR.NPHMSU.GT.NPHMAX) NPHMSU=NPHMAX
      IF(NPHMSU.GT.NDM) NPHMSU=NDM
      EDM=2*NPHMSU+(IPPM-1)*3
      IF(IAM.LT.0) IAM=EDM
      IAI=IAI+1
      IAM=IAM+1
      NPHMSP=NPHMSU+1
      NPHMXP=NPHMAX+1
      SQR=SQRT(6./35.)
C		Convert input parameters for f-boson
      CHOT=(CHON+CHOP)/2 + CHO
      HBAR3P= HBAR3 - 2*RKAP3*(NPHMAX-1)
      F3P=F3 + 2*SQRT(35.)*FQQ + 2*RKAP3*CHOT
      EPSDP= EPSD + 2*RKAP3
      FEXP = FEX - RKAP3*(CHON*CHOP + CHO**2)*2/35.
      DP(1)=D(1) - FQQ*CHQ*SQR*4       - FELL*8 + FEXP*2
      DP(2)=D(2) - FQQ*CHQ*SQR         - FELL*6 - FEXP*7/2.
      DP(3)=D(3) + FQQ*CHQ*SQR*11./6.  - FELL*3 + FEXP*19/12.
      DP(4)=D(4) + FQQ*CHQ*SQR*5/2.    + FELL   + FEXP*35/12.
      DP(5)=D(5) - FQQ*CHQ*SQR*5/3.    + FELL*6 + FEXP*5/6.
C
C      CFP'S are read in from 'PHINT.CFP'
      CALL READRE(NDM,IDGS,IDB,IDIB,REDGS,BETFAC,IBDEL)
C
C       	Convert into parameters of COMMON /H/
      IF(NPLOG) CALL PROJ
      IF(MULT) CALL MLTP
C
      IF(.NOT. PRINT) CALL PRTFHD(NPHMSU,NPHMAX,COMMNT)
C
C      	Start real computation of energies
   	ZEROEN=0.
      ALPHA=(3*C(3)+4*C(2))/7
      GAM=(C(3)-C(2))/14
      TBETA=(C(1)-ALPHA+12*GAM)/5
      DO 150 N=1,NPHMSP
C      note : n= # of phonons+1
      NBPM=(N+1)/2
      PHEN(N,1)=(HBAR-6*GAM+ALPHA*(N-2)/2)*(N-1)
      DO 150 NBP=1,NBPM
      NB=NBP-1
      PHEN(N,NBP)=PHEN(N,1)+TBETA*NB*(2*N-2*NB+1)
  150 CONTINUE
	WRITE(*,*) 'EXCITATION ENERGIES'
      NSTTOT=1
C
C		Loop over parity
      DO 121 IPP=IPPI,IPPM
      P='+'
      IF(IPP.EQ.2) P='-'
C
C  		angular momentum loop
      DO 1000 IA=IAI,IAM
      IANG=IA-1
   80 K=NSTTOT+4*IEDMA+1
      IF(IEDMA.LE.0) GOTO 127
      LEN=LENMIN+1850+K
      IF(LEN.LE.LENMAX) GOTO 81
      IEDMA=(LENMAX-LENMIN-NSTTOT-1850)/4-1
      GOTO 80
   81 CALL XRFLX(LEN,IWD)
C		generate basis states
      CALL GENST(NPHMSU,IEDMA,IA,IPP,XYZ(K),XYZ(K+50),XYZ(NSTTOT),IED)
      IF(IED.LE.0) GOTO 127
      IED=IED+1
      MV=MIN0(IED,NEIG+1)
C		Calculate the positions of the different array's in 
C			the workspace array XYZ
      IBSRT=NSTTOT+4*IED
      IBSV=IBSRT+IED
      IBSEV=IBSV+IED*IED
      IBSI=IBSEV+MV*IED+1
      IBSR=IBSI+IED
      IDBLE=2
C      Set to use double precision in EIGSYM
      IBSR2=IBSR+IED*IDBLE
      IBSR3=IBSR2+IED*IDBLE
      IBSR4=IBSR3+IED*IDBLE
      IBSR5=IBSR4+IED*IDBLE
      IBSR6=IBSR5+IED*IDBLE
      WRITE(10) IA,IPP,IED-1,MV-1
      K=NSTTOT+4*IED-4
      WRITE(10) (XYZ(I),I=NSTTOT,K)
      LEN=IBSR6+IED*IDBLE+LENMIN
      IF(LEN.GT.LENMAX) GOTO 83
      CALL XRFLX(LEN,IWD)
C      Total needed length = lenmin+ibsr4+ied
      CALL HAMILT (IA,IPP,IED,MV,NSTTOT,ZEROEN 
     & ,XYZ(1),XYZ(NSTTOT),XYZ(IBSRT),XYZ(IBSV),XYZ(IBSEV)
     & ,XYZ(IBSI),XYZ(IBSR),XYZ(IBSR2),XYZ(IBSR3),XYZ(IBSR4),
     & XYZ(IBSR5),XYZ(IBSR6))
      K=IBSRT+IED-2
      WRITE(10) (XYZ(I),I=IBSRT,K)
      J=IBSV
      K=IBSEV
      DO 130 I=1,MV
      CALL MOVLEV(XYZ(K),XYZ(J),IED)
      K=K+IED
      J=J+IED-1
  130 CONTINUE
      J=IBSV+MV*(IED-1)-1
      WRITE(10) (XYZ(I),I=IBSV,J)
      NSTTOT=NSTTOT+IED-1
      GOTO 1000
  283 FORMAT(' REQUIRED LENGTH OF',I10,' FOR THIS PROBLEM IS TOO MUCH')
   83 WRITE(2,283) LEN
C      No allowed states
  127 WRITE(2,129)IANG,P
  129 FORMAT(//,2X,'NO STATE WITH L=',I2,' AND PARITY ',A1)
      IF(IA.GT.1) GOTO 1
      IF(IPP.EQ.2) GOTO 1
      IECSBK(2,1,1)=1
      IECSBK(2,2,1)=1
      GOTO 1000
    1 IECSBK(IA+1,1,IPP)=IECSBK(IA,1,IPP)
      IECSBK(IA+1,2,IPP)=IECSBK(IA,2,IPP)
      NZF=IECSBK(IA+1,1,IPP)
      NZSQ=IECSBK(IA+1,2,IPP)
      IF(IPP.EQ.2) GOTO 1000
      DO 2 I=1,IAI
      IECSBK(I,1,2)=NZF
      IECSBK(I,2,2)=NZSQ
    2 CONTINUE
 1000 CONTINUE
  121 CONTINUE
C
      WRITE(2,203) ZEROEN,HBAR+(NPHMAX-1)*(CH2-CH1)
  203 FORMAT(/' BINDING-ENERGY =',F8.4,' , EPS-EFF =',F8.4/)
C
C       ===== WRITEPHONONTAPE ========== WRITEPHONONTAPE =====
C
      IA=-1
      WRITE(10) IA,IPP,IED,MV
      WRITE(10) IECSBK
      WRITE(10) NPHMSU,NPHMAX,NEIG,IAI,IAM,IPPI,IPPM,IEDM
      WRITE(10) ZEROEN,COMMNT
      WRITE(10) COEL,COMUL,MULT,CHQ,CONP,NN,NP,FPARM,SDEQSF
      STOP ' full output on PHINT.OUT'
	END
	BLOCKDATA
      INTEGER EVOD
      LOGICAL PRINT,PRINTV,MULT,SDEQSF,PRINTP
     &         ,NPLOG
      DIMENSION CLN(3),CLP(3)
      CHARACTER COMMNT*40, P*1
      COMMON / TEXT   / PRINT,PRINTV,PRINTP,IWD
      COMMON / TEXT2  / P,COMMNT
      COMMON / CONTR  / NPHMSU,NPHMAX,NEIG,IAI,IAM,IPPI,IPPM,EVOD
      COMMON / FPAR   / HBAR3P,DP(5),F3P,EPSDP,D(5),F3,FELL,FQQ,FEX,
     &                  HBAR3,EPSD,RKAP3,CHO,CHON,CHOP,SDEQSF
      COMMON / H      / HBAR,C(3),F,G,CH1,CH2
      COMMON / MUL    / EPS,PAIR,ELL,QQ,OCT,HEX,MULT,CHQ
      COMMON / NP     / ED,RKAP,CHN,CLN,CHP,CLP,NN,NP,NPLOG
      COMMON/ECSBK/IECSBK(50,2,2)
      DATA IECSBK/200*1/
      DATA EVOD,NPHMAX,NPHMSU,IWD,IAI,IAM,IPPI,IPPM,NEIG,PRINT,PRINTV,
     &     PRINTP,MULT,NPLOG,SDEQSF/2,2,99,0,0,-1,1,1,4,
     &      5*.FALSE.,.TRUE./
      DATA CLN,CLP,RKAP,CHN,CHP,ED,NN,NP/10*0.,2*0/
      DATA FELL,FQQ,FEX,PAIR,ELL,QQ,OCT,HEX,CHQ,CH1,CH2,EPS,F,G,HBAR,C,
     &   EPSD,HBAR3,D,F3,RKAP3,CHO,CHON,CHOP/8*0.,-2.9580399,21*0./
      END
      SUBROUTINE GENST(NDN,IEDMAX,IA,IPP,IBK,ISTP,IST,NZ)
C
C		Calculate basis states for J=IA-1 and parity IPP
C      INPUT :
C          NDM=MAX. NO. OF PHONONS
C          IAI,IAM =INITIAL,FINAL L , ONLY USED FOR NEG. PAR.
C           IPPM=1 D ONLY POS PAR. , =2 OLSO NEG. PAR.
C          IEDMAX=DIMENSION OF V AND EIGV IN MAIN PROG. , A CHECK IS
C             DONE ON THIS
C      OUTPUT :
C           IST(1,K)=ND
C           IST(2,K)=NB
C           IST(3,K)=NC
C           IST(4,K)=LD
C          NEIG IS SET TO 0 IF DIMENSION OF EIG IS NOT LARGE ENOUGH
C                  (TESTED WITH IDEIG)
C          IECSBK(L+1,1,1)=FIRST POSITION IN IECST WITH L AND POS. PAR.
C          IECSBK(L+1,1,2)= SAME FOR NEG. PAR.
C          IECSBK(L+1,2,1)= SAME BUT IN EIG WITH POS  PAR.
C          IECSBK(L+1,2,2)= SAME FOR NEG. PAR.
C
      LOGICAL SDEQSF
      DIMENSION IBK(50),ICON(7,2),ISTP(3,600),IST(4,IEDMAX)
      COMMON/CONTR/NPHMSU,NPHMAX,NEIG,IAI,IAM,IPPI,IPPM ,IEVOD
      COMMON/ECSBK/IECSBK(50,2,2)
      COMMON / FPAR / SKP(23),SDEQSF
C
      NZ=1
C       Used only for print out of basis states statistics :
      NDM=NDN
      NDPM=NDM+1
      LPM=2*NDM+1
      LPI=IAI-3
      IF(LPI.LT.1) LPI=1
      IBK(LPI)=1
      LP=IAM+3
      IF(LP.LT.LPM) LPM=LP
      DO 1 LP=LPI,LPM
      L=LP-1
C
C      Find allowed states belonging to this L
C
      NDPI=LP/2+1
      DO 2 NDP=NDPI,NDPM
      ND=NDP-1
      NBPM=(NDP-NDPI)/2+1
      DO 3 NBP=1,NBPM
      NB=NBP-1
      NC3F=(NDP-NDPI-2*NB)/3
      NC3I=(ND-2*NB-L+2)/3
      IF(NC3I.LT.0) NC3I=0
      IF(NC3F.LT.NC3I) GOTO 3
      NCPI=NC3I+1
      NCPF=NC3F+1
      DO 4 NCP=NCPI,NCPF
      NC=NCP-1
      LAMBD=ND-2*NB-3*NC
      IF(L.EQ.(2*LAMBD-1)) GOTO 4
      ISTP(1,NZ)=ND
      ISTP(2,NZ)=NB
      ISTP(3,NZ)=NC
      NZ=NZ+1
    4 CONTINUE
    3 CONTINUE
    2 CONTINUE
C
      IBK(LP+1)=NZ
    1 CONTINUE
      LPM=LPM+1
      DO 6 I=LPM,50
    6 IBK(I)=NZ
C
      IF(IPP.NE.1) GOTO 199
      NZL=IBK(IA)-1
      NZ=IBK(IA+1)-NZL-1
      LD=IA-1
      IF(NZ.GT.IEDMAX .OR. NZ.LE.0) GOTO 99
      DO 10 I=1,NZ
      K=I+NZL
      IST(1,I)=ISTP(1,K)
      IST(2,I)=ISTP(2,K)
      IST(3,I)=ISTP(3,K)
      IST(4,I)=LD
   10 CONTINUE
      GOTO 98
C
C      NEGATIVE PARITY
C
  199 NZF=0
      IF(SDEQSF) NDM=NDN-1
      NDPM=NDM+1
      LTP=IA
      LT=IA-1
C
C      SELECT ALLOWED STATES BELONGING TO THIS LT
C
      LDPI=LTP-3
      LDPF=LTP+3
      IF(LDPI.LT.1) LDPI=1
      IF(LDPI.LT.(4-LT)) LDPI=4-LT
      IF(LDPF.GT.(2*NDM+1)) LDPF=2*NDM+1
      DO 102 LDP=LDPI,LDPF
      LDEL=LTP-LDP+4
      ICON(LDEL,1)=IBK(LDP)
      ICON(LDEL,2)=IBK(LDP+1)
      ICON(LDEL,2)=ICON(LDEL,2)-1
  102 CONTINUE
      DO 103 NP=1,NDPM
      N=NP-1
      DO 104 LDP=LDPI,LDPF
      LD=LDP-1
      LDEL=LTP-LDP+4
      NZ=ICON(LDEL,1)
      NR=ICON(LDEL,2)
      IF(NR-NZ) 104,105,105
  105 DO 106 I=NZ,NR
      ND=ISTP(1,I)
      IF(ND-N) 106,107,108
  107 CONTINUE
      NZF=NZF+1
      IF(NZF.GT.IEDMAX) GOTO 106
      IST(1,NZF)=ND
      IST(2,NZF)=ISTP(2,I)
      IST(3,NZF)=ISTP(3,I)
      IST(4,NZF)=LD
  106 CONTINUE
  108 ICON(LDEL,1)=I
  104 CONTINUE
  103 CONTINUE
      NZ=NZF
   99 IF(NZ.LT.IEDMAX) GOTO 98
      WRITE(2,201) NZ,IEDMAX,IA-1
  201 FORMAT(10X,'/// WARNING FROM GENST /// ',I6,
     & ' EXCEEDS MAXIMAL DIMENSION OF H =',I6,' ,AT L=',I3)
      NZ=0
   98 NZF=IECSBK(IA,1,IPP)+NZ
C      WRITE(2,210) ((IST(I,K),I=1,4),K=1,NZ)
C  210 FORMAT(1X,5(1X,4I5))
      NZSQ=MIN0(NZ,NEIG)*NZ+IECSBK(IA+1,1,IPP)
      IECSBK(IA+1,1,IPP)=NZF
      IECSBK(IA+1,2,IPP)=NZSQ
      IF(IPP.EQ.2) RETURN
      IAP=IA+1
      DO 97 I=1,IAP
      IECSBK(I,1,2)=NZF
   97 IECSBK(I,2,2)=NZSQ
      RETURN
      END
      SUBROUTINE HAMILT (IA,IPP,IED,MV,NSTTOT
     & ,ZEROEN
     & ,ENERGY,STATE,ROOT,V,EIGV
     & ,LIG,W,Q,WW,RM
     & ,AE,BE)
C		Calculate hamiltonian
      LOGICAL PRINT,PRINTV,PRINTP
      INTEGER ED,STATE
      DIMENSION ENERGY(NSTTOT),STATE(4,IED),ROOT(IED),V(IED,IED),
     &  EIGV(IED,MV),LIG(IED),W(IED),Q(IED),WW(IED),RM(IED)
     & ,AE(IED),BE(IED)
      CHARACTER COMMNT*40, P*1
      COMMON / CONTR  / NPHMSU,NPHMAX
      COMMON/ECSBK/IECSBK(50,2,2)
      COMMON/TEXT/PRINT,PRINTV,PRINTP
      COMMON / TEXT2  / P,COMMNT
C
      MV=MV-1
      ED=IED-1
      IANG=IA-1
      DO 125 I=1,IED
      DO 125 J=I,IED
  125 V(I,J)=0.
C
C ----- FILL IN FIRST HALF OF MATRIX 'V'
C
      CALL FVULVH(IANG,IPP,IED,V,LIG,STATE)
C
C     OUTPUT , LARGE PRINT OUT
C
      IF(.NOT.PRINT) GOTO 12
      IF(IA.NE.1) WRITE(2,200)
  200 FORMAT(1H1)
      CALL PRTFHD(NPHMSU,NPHMAX,COMMNT)
      WRITE(2,26)
   26 FORMAT(5(26H NR =|ND,NB,NC,LD,NF,LP> ;))
      NF=IPP-1
     	WRITE(2,201)(I,(STATE(J,I),J=1,4),NF,IANG,P,I=1,ED)
  201 FORMAT(5(1X,I2,' =|',4(I2,','),I1,',',I2,A1,'> ;'))
      IF(.NOT.PRINTV) GOTO 12
      WRITE (2, 330)
  330 FORMAT(//,15H INITIAL MATRIX,/)
      DO 91 I=1, ED
      WRITE(2,94) I,(V(J,I),J=1,I)
   91 CONTINUE
   12 CONTINUE
C
C  	fill other half of interaction matrix
C
      DO 2 I=1,ED
      DO 1 J=I,ED
    1 V(J,I)=V(I,J)
C
C     	Get energies in increasing order
C
    2 V(I,I)=V(I,I)+99.
C
C      ==== SUBROUTINE FROM 'MATHFTN' LIBRARY ====
C       SUBROUTINE 'EIGSYM' DIAGALISES MATRIX 'V(ED,ED)' AND
C        STORES 'MV' EIGENVECTORS IN 'EIGV(ED,MV)' . ALL THE
C        EIGENVALUES ARE STORED IN 'ROOT(ED)' , ORDERED IN INCREASING
C        ABSOLUTE MAGNITUDE WHEN 'MV' IS NEGATIVE .
C       'IEDMAX' IS THE DIMENSION OF 'V' AND 'EIGV' AS DECLARED IN
C        THE DIMENSION STATEMENT
C       ARRAY'S 'AE' TILL 'LIG' ARE WORK SPACE
C     LAST INDEX IN EIGV ENNUMERATES THE EIGENVECTORS
C
      NITER=0
      CALL EIGSYM(ED,IED,-MV,V,ROOT,EIGV,AE,BE,W,Q,WW,RM,LIG,NITER)
C
      DO 3 I=1,ED
    3 ROOT(I)=ROOT(I)-99.
      IF((IA+IPP).EQ.2) ZEROEN=ROOT(1)
      NK=IECSBK(IA,1,IPP)
      NKI=NK
      DO 193 J=1,ED
      LIG(J)=STATE(1,J)
C      STATE will be overwritten here, but nd necessary for prob.
      ENERGY(NK)=ROOT(J)-ZEROEN
      NK=NK+1
  193 CONTINUE
      NKF=NK-1
C     Last index in EIGV ennumerates the eigenvectors
C
C     OUTPUT
C
	J=MIN(MV,10)-1
	WRITE(*,80) IANG,P,(ENERGY(I),I=NKI,NKI+J)
80	FORMAT(I3,A1,':',10F7.3)
      IF(PRINT) GOTO 11
C      SMALL PRINT OUT
      WRITE(2,97) IANG,P,(ENERGY(I),I=NKI,NKF)
      GOTO 126
C      LARGE PRINT OUT
   11 WRITE(2,92) IANG,P,(ROOT(I),I=1,ED)
   92 FORMAT(//,' EIGENVALUES , L=',I2,A1/8(2X,16F8.4/))
      WRITE(2,97) IANG,P,(ENERGY(I),I=NKI,NKF)
   97 FORMAT(//,' ENERGIES , L=',I2,A1/8(2X,16F8.4/))
      WRITE(2,93)
   93 FORMAT(//,13H EIGENVECTORS,/)
      DO 95 I=1,ED
      WRITE(2,94)I,(EIGV(I,J),J=1,MV)
   94 FORMAT(' <',I2,'> ',17F7.3,4(/5X,17F7.3))
   95 CONTINUE
  126 CONTINUE
      IF(.NOT.PRINTP) RETURN
      NPHMSP=LIG(ED)+1
      DO 32 M=1,MV
      DO 33 I=1,NPHMSP
   33 AE(I)=0.
      DO 31 I=1,ED
      NP=LIG(I)+1
      E=EIGV(I,M)
   31 AE(NP)=AE(NP)+E*E
      WRITE(2,230) M,(AE(I),I=1,NPHMSP)
  230 FORMAT(1X,I3,18F7.4)
   32 CONTINUE
      RETURN
      END
*DECK FVULVH
      SUBROUTINE FVULVH(IBNG,IPP,IEDP1,V,ISTBK,STATE)
C
C		Calculate upper half of hamiltonian matrix
C      NPHMXP has here the meaning of NPHMSP
C      and determines the max. no. of phonons to be included
C      IPP=2*I+IIP WITH :
C                         IIP=1 => POS.PARITY
C                         IIP=2 => NEG.PARITY
C      OUTPUT :
C         ED = used dimension in 'V'
C         V = hamiltonian matrix
C
      LOGICAL DIAG,SDEQSF
      INTEGER ED,STATE
      DIMENSION STATE(4,IEDP1),ISTBK(IEDP1),V(IEDP1,IEDP1)
      COMMON/FPAR/ HBAR3,D(5),F3,EPSD(16),SDEQSF
      COMMON/STAB/N,NBA,NCA,LDA,NF,NR,NBB,NCB,LDB,NP,ODDSP,IANG
      COMMON/CONTR/NPHMSU,NPHMAX
C      Change nphmax only for negative parity states
      NPH=NPHMAX
      IF(SDEQSF .AND. (IPP.EQ.2)) NPHMAX=NPHMAX-1
      NPHMXP=NPHMAX+1
      IANG=IBNG
      ED=IEDP1-1
C
C      Set up array 'ISTBK' which keeps track of the places in
C       'STATE' where ND changes
C
      ISTBK(NPHMXP+2)=0
      ISTBK(NPHMXP+3)=0
      II=1
      ISTBK(1)=1
      DO 10 NP=1,NPHMXP
      N=NP-1
      IF(ED-II) 14,12,12
   12 DO 11 I=II,ED
      ND=STATE(1,I)
      IF(ND-N) 14,11,14
   11 CONTINUE
      I=ED+1
   14 ISTBK(NP+1)=I
      II=I
   10 CONTINUE
      NRA=0
      NF=IPP-1
C
C       set up for main loop
C
      I1=ISTBK(1)
      I2=ISTBK(2)-1
      I3=ISTBK(3)-1
      I4=ISTBK(4)-1
C
C      MAIN LOOP
C
      DO 100 NP=1,NPHMXP
      N=NP-1
      IF(I2-I1) 99,101,101
  101 I2P=I2+1
      IB1=I3-I2P
      I3P=I3+1
      IB2=I4-I3P
      DO 102 IA=I1,I2
      NRA=NRA+1
      NRB=NRA
      NBA=STATE(2,IA)
      NCA=STATE(3,IA)
      LDA=STATE(4,IA)
      DIAG=.TRUE.
      DO 110 IB=IA,I2
      NBB=STATE(2,IB)
      NCB=STATE(3,IB)
      LDB=STATE(4,IB)
      V(NRA,NRB)=HAM023(DIAG)
      DIAG=.FALSE.
  110 NRB=NRB+1
      IF(IB1) 103,104,104
  104 DO 120 IB=I2P,I3
      NBB=STATE(2,IB)
      NCB=STATE(3,IB)
      LDB=STATE(4,IB)
      H=0.
      IF(LDA.EQ.LDB) CALL HAM1(H)
      IF(NF) 122,122,121
C      F3*<N/S+*D*F+*F/ >
  121 LD=LDA-LDB+3
      IF(LD.LT.1.OR.LD.GT.5) GOTO 122
      H=H-F3*RED(N,NBA,NCA,LDA,LD,NBB-NBA,NCB-NCA)*
     *   RACAHI(IANG,3,LDA,2,LDB,3)*SQRT(FLOAT(NPHMAX-N))
  122 V(NRA,NRB)=H
  120 NRB=NRB+1
  103 IF(IB2) 102,106,106
  106 DO 130 IB=I3P,I4
      LDB=STATE(4,IB)
      IF(LDA.NE.LDB) GOTO 130
      NBB=STATE(2,IB)
      NCB=STATE(3,IB)
      CALL HAM2(V(NRA,NRB))
  130 NRB=NRB+1
  102 CONTINUE
C
C      keep track of bookkeeping
   99 I1=I2+1
      I2=I3
      I3=I4
      I4=ISTBK(NP+4)-1
  100 CONTINUE
      NPHMAX=NPH
      RETURN
      END
      FUNCTION HAM023(DIAG)
C
C      Computes matrix elements between states with the same phonon numbers,
C      matrix elements need not to be diagonal.
C      EPSD*< /(N-S+*S)*F+*F/ >
C
      LOGICAL DIAG
      COMMON/HAM0/PHEN(8 , 4),GAM
      COMMON/STAB/N,NBL,NCL,LL,N3M,NR,NBR,NCR,LR,NP,ODDSP,IANG
      COMMON/FPAR/ HBAR3,E(5),F3,EPSD
      COMMON/H/HBAR,C(3),F,G,CH1,CH2
      COMMON/CONTR/NPHMSU,NPHMAX
C
      HAM023=0.
      NS=N-1
C      N3M = number of f-bosons  ( 0 or 1 )
      IF(N3M) 9,5,6
C
C      one three minus phonon
    6 IF(.NOT.DIAG) GOTO 62
C
C      diagonal matrix element
C      LL=LR
      NBL=NBR
      NCL=NCR
C
   61 HAM023=HBAR3+EPSD*N
      DO 611 LSL=1,5
      LS=LL+LSL-3
      H=0.
      RD=0.
      DO 614 K=1,5
      R=RACAHI(LS,2,IANG,3,LL,K)
  614 H=H+(2*K+1)*E(K)*R*R
      DO 613 NBLSP=1,2
      NBLS=NBLSP-1
      NBS=NBL-NBLS
      DO 613 NBCS=1,2
      NCLS=NBCS-NBLSP
      NCS=NCL-NCLS
      R=RED(NS,NBS,NCS,LS,LSL,NBLS,NCLS)
  613 RD=RD+R*R
  611 HAM023=HAM023+H*RD
      GOTO 51
   62 CONTINUE
C      non diagonal part
C
      DO 621 LSL=1,5
      LS=LL+LSL-3
      LSR=LS-LR+3
      IF(LSR.LT.1) GOTO 621
      IF(LSR.GT.5) GOTO 9
      H=0.
      RD=0.
      DO 624 K=1,5
  624 H=H+(2*K+1)*RACAHI(LS,2,IANG,3,LL,K)*RACAHI(LS,2,IANG,3,LR,K)*E(K)
      DO 623 NBLSP=1,2
      NBLS=NBLSP-1
      NBS=NBL-NBLS
      NBRS=NBR-NBS
      DO 623 NBCS=1,2
      NCLS=NBCS-NBLSP
      NCS=NCL-NCLS
      NCRS=NCR-NCS
  623 RD=RD+RED(NS,NBS,NCS,LS,LSR,NBRS,NCRS)*
     *    RED(NS,NBS,NCS,LS,LSL,NBLS,NCLS)
      HAM023=HAM023+H*RD
  621 CONTINUE
      GOTO 9
    5 CONTINUE
C
C      no f-boson
      IF(.NOT.DIAG) RETURN
   51 HAM023=LL*(LL+1)*GAM+PHEN(N+1,NBL+1)+HAM023
     &   +CH2*(NPHMAX-N)*N+CH1*(NPHMAX-N)*(NPHMAX-N-1)/2.
    9 RETURN
      END
      SUBROUTINE RDXW
      INTEGER EVOD
	CHARACTER*1 ANS
	CHARACTER*10 TXT
      LOGICAL PRINT,PRINTV,MULT,WRTAPE,SDEQSF,PRINTP
     &    ,NPLOG
      CHARACTER COMMNT*40, P*1
      DIMENSION CLN(3),CLP(3),CONP(6),COEF(6),COEL(8),COMUL(6)
     &          ,ENERGY(1),FPARM(23)
      COMMON / TEXT   / PRINT,PRINTV,PRINTP,IWD
      COMMON / TEXT2  / P,COMMNT
      COMMON / CONTR  / NPHMSU,NPHMAX,NEIG,IAI,IAM,IPPI,IPPM,EVOD
      COMMON / FPAR   / HBAR3P,DP(5),F3P,EPSDP,D(5),F3,FELL,FQQ,FEX,
     &                  HBAR3,EPSD,RKAP3,CHO,CHON,CHOP,SDEQSF
      COMMON / H      / HBAR,C(3),F,G,CH1,CH2
      COMMON / MUL    / EPS,PAIR,ELL,QQ,OCT,HEX,MULT,CHQ
      COMMON / NP     / ED,RKAP,CHN,CLN,CHP,CLP,NN,NP,NPLOG
C
C   NPHMAX    : PLACE OF CUT OF
C   NPHMSU    : NUMBER OF PH. INCLUDED IN THE CALC.
C   SDEQSF    : =.T. : TAKE S PHONON FOR D EQUAL TO S PHONON FOR F PH
C   IAI       : INITIAL L VALUE TO BE COMP.
C   IAM       : FINAL L VALUE TO BE COMP.
C   IPPI      : =1 : COMP. POS. PAR. ; =2 : NO POS. PAR.
C   IPPM      : =1 : NO NEG. PAR. ; =2 : COMP. NEG. PAR.
C   PRINT     : =.T. : PRINT ALSO INITIAL MATRIX AND EIG. VECTORS
C   PRINTV    : =.F. : DO NOT PRINT INITIAL MATRIX
C   NEIG      : MAX NUMBER OF EIGENVECTORS TO BE COMPUTED
C   D     (5) : <(DF)L|H|(DF)L>
C   F3        : <(SF)3|H|(DF)3>*-SQRT(7)
C   EPSD      : <DF|ND*NF|DF>
C
      OPEN(UNIT=1,STATUS='SCRATCH')
	WRITE(*,200)
200	FORMAT('  Hamiltonian parameters, use IBA-2 projection? (n/y)')
	READ(*,1) ANS
    1 FORMAT(A1)
	IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') GOTO 10
	WRITE(*,201)
201	FORMAT(' Number of bosons =')
	READ(*,*) NPHMAX
	WRITE(*,202) 
202	FORMAT(' Hamiltonian parameters, use multipole form? (y/n)')
	READ(*,1) ANS
	IF(ANS.EQ.'N' .or. ANS.EQ.'n') GOTO 20
      MULT=.TRUE.
	WRITE(*,203) 
203	FORMAT(' Multipole form will be used.',/)
C   EPS       : PHONON ENERGY , ADS TO THE EFFECT OF THE MULTIPOLES
C                   ON HBAR
C   ELL       :   L.L   FORCE WITH STRENGTH PROP. TO 'ELL'
C   QQ        :   Q.Q   FORCE WITH STRENGTH PROP. TO 'QQ'
C   OCT       : OCTUPOLE PART
C   HEX       : HEXADECUPOLE PART
1002	CONTINUE
	CALL RDPAR('EPS ',EPS,1)
	CALL RDPAR('ELL ',ELL,1)
	CALL RDPAR(' QQ ',QQ,1)
	CALL RDPAR('CHQ ',CHQ,1)
	CALL RDPAR('OCT ',OCT,1)
	CALL RDPAR('HEX ',HEX,1)
	WRITE(*,204) EPS,ELL,QQ,CHQ,OCT,HEX
204 	FORMAT(2X,'EPS=',F7.4,' , ELL =',F7.4,' , QQ =',F7.4,
     &  ' , CHQ =',F7.4,/,2X,'OCT=',F7.4,' , HEX=',F7.4)
	WRITE(*,*) 'Are these parameters OK? y/n'
	READ(*,1) ANS
	IF(ANS.EQ.'N' .OR. ANS.EQ.'n') GOTO 1002
      GOTO 30
10	WRITE(*,210) 
210	FORMAT(' IBA-2 projection will be used')
	MULT=.TRUE.
11	WRITE(*,*) '# neutron bosons'
	READ(*,*) NN
	WRITE(*,212) 
212	FORMAT(' # proton bosons')
	READ(*,*) NP
	CALL RDPAR(' ED ',ED,1)
	CALL RDPAR('RKAP',RKAP,1)
	CALL RDPAR('CHN ',CHN,1)
	CALL RDPAR('CHP ',CHP,1)
	CALL RDPAR('CLN ',CLN,3)
	CALL RDPAR('CLP ',CLP,3)
      WRITE(*,213) NN,CHN,(CLN(I),I=1,3),NP,CHP,(CLP(I),I=1,3),ED,RKAP
  213 FORMAT(10X,'PROJECTION FROM :'/
     &        12X,'NN=',I2,' , CHN=',F5.2,' , CLN=',3(F5.2,1H,)/
     &        12X,'NP=',I2,' , CHP=',F5.2,' , CLP=',3(F5.2,1H,)/
     &       /15X,'ED=',F5.3,' , RKAP=',F7.4)
	WRITE(*,*) 'Are these parameters OK? y/n'
	READ(*,1) ANS
	IF(ANS.EQ.'N' .OR. ANS.EQ.'n') GOTO 11
	GOTO 30
20	WRITE(*,220) 
220	FORMAT(' Non-multipole form is used'/)
C   C     (3) : <(DD)L|H|(DD)L>
C   F         : <(SD)2|H|(DD)2>*SQRT(5/2)
C   G         : <(SS)0|H|(DD)0>/2
C   CH1       : <(SS)0|H|(SS)0>
C   CH2       : <(SD)2|H|(SD)2>
21	CALL RDPAR('HBAR',HBAR,1)
	CALL RDPAR('C(3)',C,3)
	CALL RDPAR(' F  ',F,1)
	CALL RDPAR(' G  ',G,1)
	CALL RDPAR('CH1 ',CH1,1)
	CALL RDPAR('CH2 ',CH2,1)
      WRITE(*,221) HBAR,C,F,G,CH1,CH2
221 	FORMAT(3X,'HBAR =',F8.5,' , C =',3F8.5,/,3X,'F =',F8.5,
     &          ' , G =',F8.5,' , CH1 =',F8.5,' , CH2 =',F8.5)
	WRITE(*,*) 'Are these parameters OK? y/n'
	READ(*,1) ANS
	IF(ANS.EQ.'N' .OR. ANS.EQ.'n') GOTO 21
30 	CONTINUE
	WRITE(*,230) 
230	FORMAT(' Should negative parity states be calculated? n/y')
	READ(*,1) ANS
	IF(ANS.NE.'Y' .AND. ANS.NE.'y') GOTO 50
	IPPM=2
      IF(.NOT. MULT) GOTO 40
	WRITE(*,231) 
231	FORMAT(' Multipole form will be used.'/)
C   HBAR3     : 3- PHONON ENERGY 
C   FELL      :   L.L   FORCE WITH STRENGTH PROP. TO 'ELL'
C   FQQ       :   Q.Q   FORCE WITH STRENGTH PROP. TO 'QQ'
C   RKAP3     : OCTUPOLE PART
C   CHOP&CHON :
31	CONTINUE
	CALL RDPAR('EPS3',HBAR3,1)
	CALL RDPAR('FELL',FELL,1)
	CALL RDPAR('FQQ ',FQQ,1)
	CALL RDPAR('KAP3',RKAP3,3)
	CALL RDPAR('CHO ',CHO,1)
      WRITE(*,232) HBAR3,FELL,FQQ,RKAP3,CHO
232 	FORMAT(3X,'EPS3 =',F8.5,' , FELL =',F8.5,' , FQQ =',F8.5,
     &          ' , RKAP3 =',F8.5,' , CHO =',F8.5)
	WRITE(*,*) 'Are these parameters OK? y/n'
	READ(*,1) ANS
	IF(ANS.EQ.'N' .OR. ANS.EQ.'n') GOTO 31
      GOTO 50
40	WRITE(*,240) 
240	FORMAT(' Non-multipole form is used.'/)
C   D     (5) : <(DF)L|H|(DF)L>
C   F3        : <(SF)3|H|(DF)3>*-SQRT(7)
C   EPSD      : <DF|ND*NF|DF>
41	CONTINUE
	CALL RDPAR('EPS3',HBAR3,1)
	WRITE(*,241) 
241	FORMAT(' D=')
	READ(*,3) DP
    3 FORMAT(5F10.2)
	CALL RDPAR(' F3 ',F3,1)
	CALL RDPAR('EPSD',EPSD,1)
      WRITE(*,242) HBAR3,DP,F3,EPSD
242 	FORMAT(3X,'EPS3 =',F8.5,' , D=',5F8.5,/,3X,'F3 =',F8.5,
     & ' , EPSD =',F8.5)
	WRITE(*,*) 'Are these parameters OK? y/n'
	READ(*,1) ANS
	IF(ANS.EQ.'N' .OR. ANS.EQ.'n') GOTO 41
50	WRITE(*,250) 
250	FORMAT(' Print eigenvectors? n/y')
	READ(*,1) ANS
	IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') PRINT=.TRUE.
  	WRITE(*,251) 
251	FORMAT(' Print hamiltonian matrix? n/y')
	READ(*,1) ANS
	IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') PRINTV=.TRUE.
    	WRITE(*,252) 
252	FORMAT(' Print probability distribution for eigenvectors? n/y')
	READ(*,1) ANS
	IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') PRINTP=.TRUE.
	WRITE(*,253) NEIG
253	FORMAT(' How many eigenvectors per spin, default =',I3,/)
	CALL RDPARI('NEIG',NEIG,1)
	RETURN
	END
