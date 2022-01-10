C   ================================================
C 
C 
C         DLIB library directory, routines for PC-version of PHINT
C          REVISED 20 JANUARY 1982 : FTC,AND CHQ INCLUDED
C         COLUMN 4 : 1 = NOT CALLED FROM LIBRARY ROUTINE
C 
C         -------------------------
C 
C 
C     SUBROUTINE    CALLS TO                COMMON BLOCKS
C 
C  1  DPD           RACAHI , RED
C  1  HAM1          RACAHI , RED            STAB  ,H     ,CONTR
C  1  HAM2                                  STAB  ,H     ,CONTR
C     MOVLEV
C  1  PROJ          MLTP                    NP    ,MUL   ,H
C     MLTP                                  MUL   ,H
C  1  NLINPT
C  1  PRTFHD                               FPAR  ,MUL   ,H
C  1  RDCRD
C  1  RDFTP                               FECSBK,FPAR,H,MUL,NP,FCONTR,CMT,CMT2
C  1  READRE
C     RED                                   REDMAT
C  1  XRFLX
C
C   ===============================================
C
      FUNCTION DPD(N,NBL,NCL,LL,NBR,NCR,LR,LC)
C 
C      DPD=<N,NBL,NCL,LL//[D+*D]LC//N,NBR,NCR,LR>
C                                        /SQRT(2*LC+1)
C 
      DPD=0.  
      NS=N-1
      DO 1 LSR=1,5  
      LS=LR+LSR-3  
      LSL=LS-LL+3
      IF(LSL.LT.1) GOTO 1
      IF(LSL.GT.5) GOTO 2
      RD=0.
      DO 3 NBRSP=1,2  
      NBRS=NBRSP-1  
      NBS=NBR-NBRS  
      NBLS=NBL-NBS
      DO 3 NBCS=1,2  
      NCRS=NBCS-NBRSP  
      NCS=NCR-NCRS  
      NCLS=NCL-NCS
      RD=RD+RED(NS,NBS,NCS,LS,LSR,NBRS,NCRS)*
     *      RED(NS,NBS,NCS,LS,LSL,NBLS,NCLS)
    3 CONTINUE
      DPD=DPD+RD*RACAHI(LL,LC,LS,2,LR,2)
    1 CONTINUE
    2 CONTINUE
      RETURN
      END
      SUBROUTINE HAM1(H)
C		delta Nd=1 hamiltonian matrixelements 
C      right side has (n+1) phonons
C      A*<LL,N/S+*D+*(D*D)2/LR,NR=N+1 >
C 
      COMMON/STAB/N,NBL,NCL,LL,NPL,NR,NBR,NCR,LR,NPR
      COMMON/H/E(4),F,G
      COMMON/CONTR/NPHMSU,NPHMAX,NEIG(6)
C      H1=SQRT(0.4)
      DATA H1/0.632455532/
      H=0.
      IF(N.GE.2) GOTO 1
      H=F*N*H1*SQRT(FLOAT(NPHMAX-N))
      RETURN
    1 CONTINUE
C 
C      L2 HAS N PHONONS
      L2PI=LR-1  
      IF(L2PI.LE.1) L2PI=1
      L2PF=LR+3
C 
C      L1 HAS (N-1) PHONONS
      L1PI=LL-1  
      IF(L1PI.LE.1) L1PI=1
      L1PF=LL+3
      DO 21 L2P=L2PI,L2PF  
      L2=L2P-1  
      L2R=L2-LR+3
      DO 22 NBD2P=1,2  
      NBD2=NBD2P-1  
      NB2=NBR-NBD2
      DO 22 NBC2=1,2  
      NCD2=NBC2-NBD2P  
      NC2=NCR-NCD2
      HAM=0
      DO 11 L1P=L1PI,L1PF  
      L1=L1P-1  
      L12=L1-L2+3
      IF(L12.LE.0) GOTO 11
      IF(L12.GT.5) GOTO 22
      L1L=L1-LL+3
      C=0.
      DO 12 NBD1P=1,2  
      NBD1=NBD1P-1  
      NB1=NBL-NBD1
      NB21=NB2-NB1
      DO 12 NBC1=1,2  
      NCD1=NBC1-NBD1P  
      NC1=NCL-NCD1
      NC21=NC2-NC1
      C=C+RED(N-1,NB1,NC1,L1,L1L,NBD1,NCD1)*RED(N-1,NB1,NC1,L1,L12,NB21,
     &         NC21)
   12 CONTINUE
      HAM=HAM+C*RACAHI(LL,2,L1,2,L2,2)
   11 CONTINUE
   22 H=H      +HAM*RED(N,NB2,NC2,L2,L2R,NBD2,NCD2)
   21 CONTINUE
      H=H*F*SQRT(FLOAT(NPHMAX-N))/(2*LL+1)
      RETURN
      END
      SUBROUTINE HAM2(H)
C 		Two d-boson changing term in H, use analytic formula
C      G*<LL,N/S+*S+*D*D/LR,NR=N+2>
C 
      COMMON/STAB/N,NBL,NCL,LL,NPL,NR,NBR,NCR,LR,NPR
      COMMON/CONTR/NPHMSU,NPHMAX,NEIG(6)
      COMMON/H/HBAR,C(3),F,G
      H=0.
      IF(NBL-NBR+1) 24,1,24
    1 IF(NCL-NCR) 24,2,24
    2 IF(LL-LR) 24,3,24
    3 K=2*NBL  
      H=G*SQRT(0.2*(K+2)*(2*N-K+5)*(NPHMAX-N)*(NPHMAX-1-N))
   24 RETURN
      END
      SUBROUTINE MOVLEV(A,B,N)
C		C copies N elements of array A to array B
      DIMENSION A(1),B(1)
      DO 1 I=1,N
      B(I)=A(I)
    1 CONTINUE
      RETURN
      END
      SUBROUTINE PROJ
C		Project IBA-2 parameters to IBA-1 
      COMMON/NP/ED,RKAP,CHN,CLN(3),CHP,CLP(3),NN,NP
      COMMON/MUL/EPS,PAIR,ELL,QQ,OCT,HEX,MULT,CHQ
      COMMON/H/HBAR,C(3),F,G,CH1,CH2
      M=(NN+NP)*(NN+NP-1)
      QQ=2.*RKAP*NP*NN/M
      CHQ=(CHN+CHP)*SQRT(5.)/2.
      PAIR=0.
      CHS=CHQ*CHQ/5.
      CH=CHS-CHN*CHP
      C(1)=NN*(NN-1)*CLN(1)/M+NP*(NP-1)*CLP(1)/M-  CH*QQ
      C(2)=NN*(NN-1)*CLN(2)/M+NP*(NP-1)*CLP(2)/M+3*CH*QQ/14.
      C(3)=NN*(NN-1)*CLN(3)/M+NP*(NP-1)*CLP(3)/M-2*CH*QQ/7.
      ELL=(C(1)-5*C(2)+54*C(3))/225.
      OCT=(7*C(2)-2*C(1)+9*ELL)/84.
      HEX=7*(C(3)-4*ELL-OCT)
      EPS=ED+(2-CHS/2.)*QQ-3*ELL-7*OCT-9*HEX
      CALL MLTP
      RETURN
      END
      SUBROUTINE MLTP
C		Convert multipole form to SU(5) form of hamiltonian
      COMMON/H/HBAR,C(3),F,G,CH1,CH2
      COMMON/MUL/EPS,PAIR,ELL,QQ,OCT,HEX,MULT,CHQ
      CHS=CHQ*CHQ/5.
      HBAR=EPS+(CHS/2.-2)*QQ+3*ELL+7*OCT+9*HEX
      CH1=PAIR
      CH2=QQ
      C(1)=CHS*QQ+5*PAIR-6*ELL-14*OCT+18*HEX
      C(2)=8*OCT-3*CHS*QQ/14.-3*ELL+36*HEX/7
      C(3)=OCT+2*CHS*QQ/7.+4*ELL+HEX/7
      F=CHQ*QQ
      G=SQRT(1.25)*(QQ-PAIR)
      RETURN
      END
      SUBROUTINE PRTFHD(NPHMSU,NPHMAX,COMMNT)
C		Print header line with parameter values
      LOGICAL MULT,SDEQSF,NPLOG
      CHARACTER*40 COMMNT
      COMMON/FPAR/HBAR3P,DP(5),F3P,EPSDP,D(5),F3,FELL,FQQ,FEX,
     &  HBAR3,EPSD,RKAP3,CHO,CHON,CHOP,SDEQSF
      COMMON/NP/ED,RKAP,CHN,CLN(3),CHP,CLP(3),NN,NP,NPLOG
      COMMON/MUL/EPS,PAIR,ELL,QQ,OCT,HEX,MULT,CHQ
      COMMON/H/HBAR,C(3),F,G,CH1,CH2
      WRITE(2,109) COMMNT
  109 FORMAT(1H0,1X,A40)
C      IF(SDEQSF) WRITE(2,202)
C  202 FORMAT(10X,'SD EQUALS SF')
      WRITE(2,221) NPHMAX,NPHMSU
  221 FORMAT(/10X,'TOTAL NUMBER OF BOSONS = ',I2
     &       /16X,'TRUNCATION AT ND = ',I2)
C 
      IF(.NOT.NPLOG) GOTO 50
      WRITE(2,222) NN,CHN,(CLN(I),I=1,3),NP,CHP,(CLP(I),I=1,3),ED,RKAP
  222 FORMAT(/10X,'PROJECTION FROM :'//
     &        12X,'NN=',I2,' , CHN=',F5.2,' , CLN=',3(F5.2,1H,)/
     &        12X,'NP=',I2,' , CHP=',F5.2,' , CLP=',3(F5.2,1H,)/
     &       /15X,'ED=',F5.3,' , RKAP=',F7.4//)
   50 IF(.NOT.MULT.AND. .NOT.NPLOG) GOTO 60
      WRITE(2,223) EPS,QQ,CHQ,ELL,OCT,HEX
  223 FORMAT(/10X,'MULTIPOLE EXPANTION :'//
     &        12X,'EPS=',F7.4,' , QQ  =',F7.4,' , CHQ=',F7.4/
     &        12X,'ELL=',F7.4,' , OCT =',F7.4,' , HEX=',F7.4)
   60 CONTINUE
      WRITE(2,101) CH1,CH2,EPSD,FELL,FQQ,RKAP3,CHO
  101 FORMAT(/, 3X,'CH1 =',F8.5,' , CH2 =',F8.5,' , EPSD =',F8.5,
     &          ' , FELL =',F8.5,' , FQQ =',F8.5,/
     &          34X,'KAP3 =',F8.5,' , CHO =',F8.5/)
      WRITE(2,102) HBAR,HBAR3,C(1),D(1),F,G,F3,C(2),D(2),C(3),D(3),D(4),
     &D(5)
  102 FORMAT(1H0,' 2+ ENERGY 3- ENERGY   I 2+_2+ INTER. I 2+_3- INTER.
     & ONE PHONON   TWO PHONON   F3 (S+F+DF)'/2X,F8.5,2X,F8.5,4X,1H0,2X,
     &F8.5,4X,1H1,2X,F8.5,6X,F8.5,5X,F8.5,6X,F8.5/24X,1H2,2X,F8.5,4X,1H2
     &,2X,F8.5/24X,1H4,2X,F8.5,4X,1H3,2X,F8.5/39X,1H4,2X,F8.5/39X,1H5,2X
     &,F8.5///)
      RETURN
      END
      SUBROUTINE RDFTP(IO,IA,IPP,IED,MV,IER,MVR,STATE,ENERG,EIG)
C		Read wavefunction file as written by PHINT
C		ENERG excitation energies
C		EIG   eigenvectors of lowest MVR states 
C			 for this particular (J=IA-1,IPP)
C		STATE basis states
      LOGICAL SDEQSF,MULT,SDEQ
      INTEGER STATE(4,IED),SCR
      DIMENSION ENERG(IED),EIG(IED,MV)
      COMMON/ FECSBK / IECSBK(50,2,2)
      COMMON/ FPAR   / FPARM(23),SDEQ
      COMMON/ H      / COEL(8)
      COMMON/ MUL    / COMUL(6),MULT,CHQ
      COMMON/ NP     / CONP(6),NN,NP
      COMMON/ FCONTR / IAI,IAM,IPPI,IPPM,SDEQSF,NPHMSU,
     &                    NPHMAX,NEIG
      CHARACTER*40 COMNT
      COMMON/ CMT    / ZEROEN
      COMMON /CMT2   / COMNT
      REWIND(IO)
      IER=0
      IF(IA.LE.0) GOTO 100
C 
    1 READ(IO) IAR,IPPR,IEDR,MVR
      IF(IAR.LE.0) GOTO 30
      IF(IAR.NE.IA) GOTO 10
      IF(IPPR.NE.IPP) GOTO 10
      IF(IEDR.NE.IED) GOTO 20
      MVR=MIN0(MV,MVR)
      READ(IO) ((STATE(I,J),I=1,4),J=1,IED)
      READ(IO) (ENERG(J),J=1,IED)
      READ(IO) ((EIG(J,I),J=1,IED),I=1,MVR)
      RETURN
   10 READ(IO) SCR
      READ(IO) SCR
      READ(IO) SCR
      GOTO 1
C 
   20 IER=1  
      RETURN
   30 IER=2  
      RETURN
C 
  100 READ(IO) IAR
      IF(IAR.LE.0 ) GOTO 110
      READ(IO) SCR
      READ(IO) SCR
      READ(IO) SCR
      GOTO 100
  110 READ(IO) IECSBK
      READ(IO) NPHMSU,NPHMAX,NEIG,IAI,IAM,IPPI,IPPM,IEDM
      READ(IO) ZEROEN,COMNT
      READ(IO) COEL,COMUL,MULT,CHQ,CONP,NN,NP,FPARM,SDEQSF
      SDEQ=SDEQSF
      RETURN
      END
	SUBROUTINE RDPAR(ALF,PAR,N)
	DIMENSION PAR(1),IPAR(1)
	CHARACTER*4 ALF
	CHARACTER*20 TXT
	WRITE(*,3) ALF
3	FORMAT(2X,A4,'=')
	READ(*,4) TXT
4	FORMAT(A20)
	IF(TXT.NE.'                    ') READ(TXT,2) (PAR(I),I=1,N)
2 	FORMAT(3F10.4)
	RETURN
	ENTRY RDPARI(ALF,IPAR,N)
	WRITE(*,3) ALF
	READ(*,4) TXT
	IF(TXT.NE.'                    ') READ(TXT,1) (IPAR(I),I=1,N)
1 	FORMAT(5I5)
	RETURN
	END
      SUBROUTINE READRE(NDM,IDGS,IDB,IDIB,REDGS,BETFAC,IBDEL)
C 
C      read cfp from file "PHINT.CFP"
C      THESE SHOULD BE WRITTEN WITH PROGRAM 'CFPGEN'
C 
      DIMENSION REDGS(IDGS,5),BETFAC(IDB,2),IBDEL(2,IDIB)
      CHARACTER*4 AND,AIDGS,AIDB,AIDIB
      DATA AND,AIDGS,AIDB,AIDIB/'ND  ','IDGS','IDB ','IDIB'/
      R E W I N D  3
C		First make shure that the dimensions of matrices used
C			in the main program match those used in CFPGEN
C			when making this particular PHINT.CFP file
      READ(3) NDMC,IDGSC,IDBC,IDIBC
      IF(NDM.EQ.NDMC) GOTO 1
      WRITE(2,8) AND   
      GOTO 9
    1 IF(IDGS.EQ.IDGSC) GOTO 2
      WRITE(2,8) AIDGS  
      GOTO 9
    2 IF(IDB.EQ.IDBC) GOTO 3
      WRITE(2,8) AIDB  
      GOTO 9
    3 IF(IDIB.EQ.IDIBC) GOTO 4
      WRITE(2,8) AIDIB  
      GOTO 9
    4 READ(3) REDGS,BETFAC,IBDEL
      RETURN
    9 CONTINUE
    8 FORMAT(10X,'/// ERROR ///  CFP-S READ IN WITH DIFFERENT ',A4)
      WRITE(2,10) NDMC,IDGSC,IDBC,IDIBC,NDM,IDGS,IDB,IDIB
   10 FORMAT(10X,'ON TAPE3-FILE : NDM=',I2,' IDGS=',I4,' IDB=',I3,
     & ' IDIB=',I1,/
     & 10X,'IN PROGRAM :',4I6)
      STOP
      END
      FUNCTION RED(ND,NB,NC,L,LD,NBD,NCD)
C 
C      RED=<ND+1,NB+NBD,NC+NCD,LPR//D+//ND,NB,NC,L> ; LD=L-LPR+3
C      LD needs to lie inbetween 1 and 5,  no check is made on this!!
C      value returned only non zero if :
C          NBD=0 or 1 and NBD+NCD=0 or 1
C       and K,KP<=L,LP<=2*K,2*KP ; K=ND-2*NB-3*NC
C 
C	7 BOSONS :
      COMMON/REDMAT/REDGS(57,5),BETFAC(16,2),IBDEL(2,3)
      RED=0.  
      LPR=L+3-LD
      IF(ND.LT.0) GOTO 999
      IF(NB.LT.0) GOTO 999
      IF(NBD) 999,100,101
  101 IF(NBD-1) 999,200,999
C 
  100 CONTINUE
C      difference in NBETA(=(ND-seniority)/2) equal to 0
C 
      IF(NC.LT.0) GOTO 999
      IF(NCD) 999,110,111
  111 IF(NCD-1) 999,120,999
C
C      CND=0
  110 CONTINUE
      ND3=ND-2*NB-3*NC
      IF(ND3) 999,112,112
  112 IF(L.LT.ND3) GOTO 999
      IF(L.GT.2*ND3) GOTO 999
      K=ND/2
      RED=BETFAC(NB+1+(K+1)*(ND-K),1)
      K=ND3*(ND3-1)/2+1+L+IBDEL(1,NC+1)
      RED=RED*REDGS(K,LD)
      RETURN
C 
C      NCD=1
C      NC = MIN
  120 CONTINUE
      ND3=ND-2*NB-3*NC-2
      IF(ND3) 999,122,122
  122 IF(LPR.LT.ND3) GOTO 999
      IF(LPR.GT.2*ND3) GOTO 999
      K=ND/2
      RED=BETFAC(NB+1+(K+1)*(ND-K),1)
      NCP=NC+1
      K=ND3*(ND3-1)/2+1+LPR+IBDEL(2,NCP)
      RED=RED*REDGS(K,LD)
      RETURN
C 
C      NBD=1
  200 CONTINUE
      IF(NCD) 211,210,999
  211 IF(NCD+1) 999,220,999
C 
C      NCD=0
  210 CONTINUE
      IF(NC.LT.0) GOTO 999
      ND3=ND-2*NB-3*NC-1
      IF(ND3) 999,212,212
  212 IF(LPR.LT.ND3) GOTO 999
      IF(LPR.GT.2*ND3) GOTO 999
      K=ND/2
      RED=BETFAC(NB+1+(K+1)*(ND-K),2)*(-1)**(L+LPR)
      K=ND3*(ND3-1)/2+1+LPR+IBDEL(1,NC+1)
      RED=RED*REDGS(K,6-LD)
      RETURN
C 
C      CND=-1
C      NC = MAX
  220 CONTINUE
      IF(NC.LT.1) GOTO 999
      ND3=ND-2*NB-3*NC
      IF(ND3) 999,222,222
  222 IF(L.LT.ND3) GOTO 999
      IF(L.GT.2*ND3) GOTO 999
      NCM=NC-1
      IF(NCM) 999,221,221
  221 CONTINUE
      K=ND/2
      RED=BETFAC(NB+1+(K+1)*(ND-K),2)*(-1)**(L+LPR)
      K=ND3*(ND3-1)/2+1+L+IBDEL(2,NC)
      RED=RED*REDGS(K,6-LD)
      RETURN
C 
  999 RED=0.
      RETURN
      END
      SUBROUTINE XRFLX(LEN,IWD)
      IF(IWD.EQ.0) RETURN
      WRITE(2,200) LEN
  200 FORMAT(11X,'ARRAY SPACE IN USE =',I6)
      RETURN
      END
