c$debug
c***********************************************************************
      PROGRAM DWUCK4
c     Main program
c
c     If you compile this program using the Microsoft optimizing compiler
c     you should compile it in 3 blocks in order to prevent code and data
c     segment overflows.  The first block contains the routines from Dwuck4
c     up to Bdwuck4, the second from routines Bdwuck4 up to Angmom and the
c     third from routines Angmom to the end.
c
c     This program contains two sets of routines for opening the disk files
c     and interfacing to the date and time routines.  One set, Dw4pc, Pcfile
c     Pcdate and Second, are for PC and the second set, Dw4vax, Vaxfile,
c     Vaxdat and Second, are for a VAX computer.  Delete the file you
c     don't need.  In the main routine Dwuck4 choose either the call to
c     subroutine Dw4pc or subroutine Dw4vax.
c
c     In addition the routine Pcfile allows for command line entry of the
c     input and output file names for compilation under the Microsoft
c     compiler version 5.0.
c***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array1/U(800,2),V(800,2),SPACE(4800)
      Common/array2/UB(800),ff(10)
C
C     STORE STANDARD ANGLE DATA
C
      ANGLE(1)=37.0
      ANGLE(2)=0.0
      ANGLE(3)=5.0
c     Store fundamental constants
      HBARC=197.327
      AMU=931.495
      FSCON=137.036
      CHSQ=HBARC/FSCON
      IBUFF=1600
      IFIRST=0
C
 1001 CONTINUE
c      CALL DW4PC (IDAT,IFIRST)
      CALL DW4VAX(IDAT,IFIRST)
      WRITE(6,9999)IFIRST
        call second(time)
        write(*,9900)time
        write(6,9900)time
      ifirst = 1
      ibf(3) = 0
C
C     TAPE 2 STORES RADIAL INTEGRALS
      REWIND 2
C
C     TAPE 4 STORES DISTORTED WAVES
      REWIND 4
C

C     READ CARD SET 1

      WRITE(6,9502)
      WRITE(6,9804)IDAT
c
c     icon    Control integers
c     alpha   Alphanumeric run identification
c
      READ (5,9802,END=2)ICON,ALPHA
      WRITE(6,9803)ICON,ALPHA

      IF(ICON(1).NE.9) GO TO 3
    2 STOP ' End of data'
    3 CONTINUE
C
C     DIVISION FOR OVERLAY MAY BE MADE  1-ADWUCK4, 2-BDWUCK4
C
        write(*,'('' '',15a4)')alpha
        call second(time)
        write(*,'(a)') ' Adwck4 entered'
      CALL ADWUCK4
C     CALL OVERLAY(FNAME,1,0)
      IF(IBF(3).EQ.0) then
        call second(time)
        write(*,'(a)') ' Bdwck4 entered'
        CALL BDWUCK4
C        CALL OVERLAY(FNAME,2,0)
      endif
      GO TO 1001
 9502 FORMAT('0CONTROL INTEGERS')
 9802 FORMAT(20I1,15A4)
 9803 FORMAT(1H ,20I2,4X,15A4)
 9804 FORMAT('  1 2 3 4 5 6 7 8 9 A B C D E F G H I J K'
     1,      '    RUN IDENTIFICATION  ',40X
     2       ,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
 9900 format(' Elapsed time in seconds =',f7.2)
 9999 FORMAT(I1,'DWUCK4-DISTORTED WAVES U.COLORADO- PC-VAX-VERSION'
     1,'  31/Jan /1990' )
      END


c***********************************************************************
c   The following five routines should be kept for compilation for a VAX
c   Choose the version of VAXFILE for either input/output routine names
c   to be entered from a command file or from within the program.
c
c   There are 2 versions of VAXFILE, one for use in a command file
c   for entering the I/O file names and the other for a prompt within
c   DWUCK4.  Be sure to delete the version of VAXFILE you do not use.
c***********************************************************************
      SUBROUTINE DW4VAX(IDAT,ifirst)
c
c     opens files - input, output and scratch files 2 and 4
c***********************************************************************
c
      DIMENSION IDAT(6)
C
      if(ifirst.eq.0) then
       open(unit=2,status='scratch',FORM='UNFORMATTED')
       open(unit=4,status='scratch',FORM='UNFORMATTED')
C
      CALL VAXFILE(5,6,'DWUCK4 AT YOUR SERVICE')
C
      endif
        call vaxdat(idat)
      RETURN
      END

c***********************************************************************
C
C     SUBROUTINE VAXFILE (INPUT, IOUTPUT, TITLE)
C
C                                        Dale E. Prull  03/02/83
C
C     This subroutine will open the following files from a command file:
C
C                    'INFILE'  -- logical unit INPUT
C                    'OUTFILE' -- logical unit IOUTPUT
C
C     The logical name translation for these files is obtained, and
C     a two line "title block" is written to 'OUTFILE' containing:
C
C                    TITLE , date, and time
C                    'INFILE', and 'OUTFILE' logical name translations
C
C
c***********************************************************************
      SUBROUTINE VAXFILE (INPUT, IOUTPUT, TITLE)
c***********************************************************************
C
      CHARACTER*(*) TITLE
C
      CHARACTER*9 DATEBUF
      CHARACTER*8 TIMEBUF
      CHARACTER*63 NAMEIN,NAMEOUT
      INTEGER*2 LENIN,LENOUT
      INTEGER*4 SYS$TRNLOG
C
C     open the input and output files
C     -------------------------------
C
      OPEN ( UNIT=INPUT,   FILE='INFILE',  STATUS='OLD', READONLY )
      OPEN ( UNIT=IOUTPUT, FILE='OUTFILE', STATUS='NEW',
     $                     CARRIAGECONTROL='FORTRAN' )
C
C     find the logical name translation
C     ---------------------------------
C
      ISTAT=SYS$TRNLOG ( 'INFILE',  LENIN,  NAMEIN ,,, )
      ISTAT=SYS$TRNLOG ( 'OUTFILE', LENOUT, NAMEOUT,,, )
C
C     get the date and time
C     ---------------------
C
      CALL DATE ( DATEBUF )
      CALL TIME ( TIMEBUF )
C
C     output two line title block
C     ---------------------------
C
      WRITE (IOUTPUT,101) TITLE , NAMEIN(1:LENIN) , DATEBUF, TIMEBUF
  101 FORMAT (' ', A25 ,' **** ',' INPUT  FILE NAME = ', A ,
     + T105, A, 2X, A )
C
      WRITE (IOUTPUT,102) NAMEOUT(1:LENOUT)
  102 FORMAT ( 33X, 'OUTPUT FILE NAME = ', A ,/)
      RETURN
      END

cc***********************************************************************
c      SUBROUTINE VAXFILE (INPUT, IOUTPUT, TITLE)
cc
cc     Opens input and output files at a prompt from within the program
cc***********************************************************************
c      IMPLICIT REAL*8(A-H,O-Z)
c      CHARACTER*25 TITLE,INFILE,OUTFILE
c      integer*2 istatus
cC
cC     open the input and output files
cC     -------------------------------
cC
c          write(*,'(a)') ' Enter INPUT  file name '
c          read (*,'(a)')  infile
c          write(*,'(a)') ' Enter OUTPUT file name '
c          read (*,'(a)') outfile
cc
c      OPEN ( UNIT=INPUT,   FILE=INFILE,  STATUS='OLD', READONLY )
c      OPEN ( UNIT=IOUTPUT, FILE=OUTFILE, STATUS='NEW',
c     $                     CARRIAGECONTROL='FORTRAN' )
cC
cC
cC     get the date and time
cC     ---------------------
cC
c      CALL DATE ( DATEBUF )
c      CALL TIME ( TIMEBUF )
cC
cC     output two line title block
cC     ---------------------------
cC
c      WRITE (IOUTPUT,101) TITLE , INFILE , DATEBUF, TIMEBUF
c  101 FORMAT (' ', A25 ,' **** ',' INPUT  FILE NAME = ', A ,
c     + T105, A, 2X, A )
cC
c      WRITE (IOUTPUT,102) OUTFILE
c  102 FORMAT ( 33X, 'OUTPUT FILE NAME = ', A ,/)
cC
c      RETURN
c      END

c***********************************************************************
      SUBROUTINE VAXDAT(IDAT)
c***********************************************************************
      character*1 STR(8)
      character*8 STX
      DIMENSION IDAT(6)
      equivalence (str,stx)
C
      CALL IDATE(IDAT(1),IDAT(2),IDAT(3))
C
      CALL TIME(STX)
      DO 160 I=1,3
      J=3*I-2
      IDAT(I+3)=(ICHAR(STR(J))-48)*10+ICHAR(STR(J+1))-48
  160 CONTINUE
      RETURN
      END

c***********************************************************************
      SUBROUTINE SECOND(TIME)
c
C     THIS SUBROUTINE INTERFACES THE VAX SECONDS ROUTINE
C     TO THE SECONDS CALL IN THE PROGRAMS
c***********************************************************************
      REAL*8 TIME
      DATA ICALL/0/
      IF(ICALL.EQ.0) CALL LIB$INIT_TIMER()
        ICALL=1
        CALL LIB$STAT_TIMER(2,LCSEC)
        TIME=0.01*FLOAT(LCSEC)
      RETURN
      END
c$nodebug
c***********************************************************************
      SUBROUTINE ADWUCK4
c
c     Control program for first portion of DWUCK4 which reads in data
c     and calculates the form factors
c***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array1/U(800,2),V(800,2),SPACE(4800)
      Common/array2/UB(800),ff(10)
      DIMENSION VB(800)
      EQUIVALENCE (space(2401),VB)
C
      IF(ICON(1).EQ.0) then
c
c       Use standard angle data
c
        II=(ANGLE(1)-ANGLE(2))/ANGLE(3)+1.0
        IF(II.LT.0) II=0
        ANGLE(1)=II
      else
C
C       READ CARD SET 2, Angle data
c
        READ (5,9001)ANGLE
      endif
      WRITE(6,9010)ANGLE
C
C
C     READ CARD SET 3, Angular momentum parameters
C
c     L      Maximum L for partial waves
c     Nltr   Number of angular momentum transfers
c     Ltrt   Orbital angular momentum transfer
c     Jtrt   Total   angular momentum transfer
c
      READ (5,9002)L,NLTR,(LTRT(I),I=1,NLTR),(JTRT(I),I=1,NLTR)
      WRITE(6,9008)L   ,(LTRT(I),I=1,NLTR)
      WRITE(6,9009)NLTR,(JTRT(I),I=1,NLTR)
      LC=0
      IF(L.LT.0) LC=IABS(L)
      L=IABS(L)
C
C
C     READ CARD SET 4
c
c     Drf    Radial increment
c     Rz     Lower radial cutoff
c     Rmaz   Upper Radial cutoff
c     Vce    Coulomb excitation scale factor for inelastic scattering reactions
c     Fnrng  Finite range parameter for stripping reactions
c     Amass  Mass for scaling radial increment, if zero, then scaling automatic
c
      READ (5,9001)DRF,RZ,RMAX,VCE,FNRNG,AMASS
      IF(ABS(RMAX).LT.ABS(RZ)) THEN
        temp=RZ
        RZ=RMAX
        RMAX=temp
      ENDIF
      WRITE(6,9011)DRF,RZ,RMAX,VCE,FNRNG
      IF(AMASS.NE.0.0) WRITE(6,9012)AMASS
C     default value for drf
      IF(DRF.EQ.0.0) DRF=0.1
      KZ=ABS(RZ  )/DRF+1.0E-08
      K =ABS(RMAX)/DRF+1.0E-08
      KC=2*((K+1)/2)
      LPLUS=L+1
      IND=0
      INCR=0
      IBF(5)=0
      IBF(6)=0
      WRITE(6,'(a)')' PARTICLE DATA '

C
C     READ IN DISTORTED WAVE INFORMATION
C
c     ibf(4) is flag to shut off spins if no spin orbit potential
c
      IBF(4)=0
c     Initial wave
      CALL FORMF(U(1,1),V(1,1),1,0,0,IVB)
      IF(IBF(4).NE.0) IBF(5)=1
      IBF(4)=0
c     Final wave
      CALL FORMF(U(1,2),V(1,2),2,0,0,IVB)
      IF(IBF(4).NE.0) IBF(6)=1
C
C     FORM FACTOR LOOP, up to 8 form factors may be used
C
      ILINE=0
      DO 35 II=1,NLTR
      IBF(4)=0
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
      IF(ICON(3).NE.0.or.II.EQ.1) then
      DO 10 I=1,800
      UB(I)=0.0
   10 CONTINUE
      ivb=0
      IF(ICON(2).EQ.0.or.ICON(2).EQ.3) then
C
C     Here for collective or single particle form factors
C
      WRITE(6,'(a)')'0FORM FACTOR DATA '
      IF(RMAX.LT.0.0) K=KC
C
      CALL FORMF(UB,VB,3,LTRT(II),0,IVB)
C
C     HERE FOR 2 NUCLEON FORM FACTORS
C
      elseif(icon(2).eq.2.or.icon(2).eq.3) then
      CALL CATHEN(LTRT(II),JTRT(II),space(   1),space(1601),UB)
      else
      write(6,'(a)')'0Warning, no form factor computed'
      endif
      endif
c
C     L.S FORM FACTOR KLUDGE
C     Move Spin-Orbit form factor to UB
      IF(IVB.EQ.1) THEN
        DO 14 M=1,800
        UB(M)=VB(M)
   14   CONTINUE
      ENDIF
      ISTRT(II)=IS(3)
      IF(ICON(4).EQ.0.or.icon(4).eq.3) then
        IWRITE=2
        call ffprint(iwrite,k,alpha,ub(1),dr(3)
     1              ,ltrt(ii),istrt(ii),jtrt(ii),idat)
      endif
C
      IF(ICON(3).NE.0.or.II.EQ.1) then
C       Apply non-local correction to form factor
        CALL FNLOC(U,V,UB,VB,LTRT(II),Z,FNRNG,K)
      endif
      IF(VCE.NE.0.0) then
c     Store end point for integration into complex R plane
      ivb = 4
        FF(  1)=3.0*VCE*CHSQ*Z(3)*ZA(3)
     1  *(RC(3)/(FLOAT(KC)*DR(3)))**LTRT(II)/FLOAT(2*LTRT(II)+1)
        FF(  2)=0.0
        FF(  3)=0.0
        FF(  4)=1.0
        if(ltrt(ii).eq.0) ff(  1)=ff(  1)/3.0
      ENDIF
      IF(IVB.EQ.2.OR.IVB.EQ.4) THEN
C
C     SET ERROR TRAPS FOR UNBOUND STRIPPING AND COULOMB EXCITATION CASES
C
        ILINE=-1
C       RESTRICT MAX L TRANSFER IN UNBOUND STRIPPING - STORAGE LIMITED
        IF(LTRT(II).GT.7) THEN
          WRITE(6,9903)LTRT(II)
          IBF(3)=1
        ENDIF
        IF(FM(1).LT.FM(2)) THEN
C         SET ERROR FLAG
          IBF(3)=1
          WRITE(6,'(a)')'0ERROR * UNBOUND STATE FOR PICK UP NOT ALLOWED'
        ENDIF
C
        IF(IVB.EQ.2) THEN
          T1=(SQRT(FLOAT(LTRT(II)*(LTRT(II)+1))+ETA(3)**2)+ETA(3))/FK(3)
          IF(ABS(RMAX).LT.T1) WRITE(6,9904)T1,RMAX
        ENDIF
        L = FK(1)*FLOAT(KC)*DR(1) + 12.5
        IF(L.LT.LPLUS) WRITE(6,9905)
C
      ENDIF
C
C     STORE FF PARAMETERS FOR USE IN UNBND STRIP. AND Coulomb excitation
C
      FF(  5)=FLOAT(IVB)
      FF(  6)=FK2(3)
      FF(  7)=ETA(3)
      FF(  8)=FLOAT(KC)
      IF(IVB.EQ.4) FF(  7)=0.0
c     write UB and FF onto disk
      WRITE(4) UB,FF
      IS(3)=ISTRT(II)
      IF(IS(3).EQ.0) JTRT(II)=LTRT(II)+LTRT(II)
      IND=MAX0(LTRT(II)+1,IND)
      INCR=MAX0(INCR,JTRT(II))
C
      IF(ICON(4).EQ.3) then
        IWRITE=1
        call ffprint(iwrite,k,alpha,ub(1),dr(3)
     1              ,ltrt(ii),istrt(ii),jtrt(ii),idat)
      endif
C
C     CHECK ON SPINS AND STATISTICS
C
      LTR=LTRT(II)
      JTR=JTRT(II)
      IF(PHASEF(IS(1)+IS(2)+IS(3)).LT.0.0) GO TO 33
      IF(PHASEF(        JTR+IS(3)).LT.0.0) GO TO 33
      GO TO 35
   33 CONTINUE
C
C     SET ERROR FLAG
C
      IBF(3)=1
      WRITE(6,9906)IS,LTR,JTR
   35 CONTINUE
c
      IF(NLTR.GT.8) IBF(3)=1
C
C     CHECK ON COMPATIBILITY OF LMAX AND STORAGE
C
      LPLUS=MIN0(LPLUS,400/(NS(1)+NS(2)),2800/(NS(1)*NS(2)*(IND+1)))
      L=LPLUS-1
      LPL2=LPLUS*2
      IF(RMAX.LT.0.0) K=KC
      IF(ILINE.NE.0) GO TO 40
      K=2*(K/2)
      L=FK(1)*ABS(RMAX)+4.0
      L=MIN0(LPLUS-1,L)
      IF(LC.NE.0) L=MIN0(L,LC)
      LPLUS=L+1
      LPL2=LPLUS*2
   40 WRITE(6,9506)L,K,NLTR
      IBUF=LPL2*IND+1
      INC=(INCR+IS(1)+IS(2))/2+1
      IBUF=MAX0(IBUF,INC*LPLUS+233)
      IBF(1)=IBUF
C
      IBF(7)=IS(1)
      IF(IBF(5).NE.0) GO TO 51
      IS(1)=0
      FS(1)=0.0
      NS(1)=1
   51 CONTINUE
      IBF(8)=IS(2)
      IF(IBF(6).NE.0) GO TO 52
      IS(2)=0
      FS(2)=0.0
      NS(2)=1
   52 CONTINUE
      IF(ICON(15).NE.0) THEN
C
C     PRINT OUT K(R)**2 FOR DISTORTED WAVES
C
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
      DO 55 I=1,2
      DO 55 N=1,2
      IF(N.EQ.1) THEN
      WRITE(6,9057)I
      ELSE
      WRITE(6,9058)I
      ENDIF
      IND=0
      R=DR(I)
      DO 54 M=1,K,5
      MK=M+M-1
      MK4=MIN0(MK +9,K+K)
      IF(N.EQ.1) THEN
      WRITE(6,9052)R,(U(J,I),J=MK,MK4)
      ELSE
      WRITE(6,9052)R,(V(J,I),J=MK,MK4)
      ENDIF
      R=R+5.0*DR(I)
   54 CONTINUE
   55 CONTINUE
      ENDIF
   56 CONTINUE
      RETURN
C
 9001 FORMAT(10F8.4)
 9002 FORMAT(18I3)
 9008 FORMAT(18H0CARD SET 3 DATA  ,9H   LMAX =,I4,14H          LTR=,8I4)
 9009 FORMAT(18H                  ,9H   NLTR =,I4,14H        2*JTR=,8I4)
 9010 FORMAT(18H0ANGLE DATA       ,9H   THETN=,F9.4,9H   THET1=,F9.4
     1,           9H   DTHET=,F9.4,9H   A-ANG=,F9.4,9H   B-ANG=,F9.4)
 9011 FORMAT(18H0CARD SET 4 DATA  ,9H   DRF  =,F9.4,9H   RZ   =,F9.4
     1,           9H   RMAX =,F9.4,9H   VCE  =,F9.4,9H   FNRNG=,F9.4)
 9012 FORMAT(18X,9H   AMASS=,F9.4)
 9052 FORMAT(1H ,F6.2,1P10E12.4)
 9057 FORMAT(38H0Central    K(R)**2 FOR DISTORTED WAVE,I2)
 9058 FORMAT(38H0Spin orbit K(R)**2 FOR DISTORTED WAVE,I2)
 9500 FORMAT(1P5E16.7)
 9506 FORMAT(1H0,17X,9H   LMAX =,I4,14H        NSTEP=,I4
     1,                             14H        NLTR =,I4)
c 9803 FORMAT(1H ,16I2,4X,15A4,I2,2(1H/,I2.2),I4,2(1H.,I2.2))
 9903 FORMAT('0ERROR ** LTR=',I3,'  is too large' )
 9904 FORMAT('0*****  Warning  *****  CLASSICAL TURNING POINT FOR BOUND'
     1,'STATE =',F8.4,'  IS GREATER THAN RMAX =',F8.4)
 9905 FORMAT('0*****  Warning *****  LMAX too large for reliable '
     1,'calculation of unbound stripping or coulomb excitation cases.'
     2,'  Increase RMAX ')
 9906 FORMAT(1H0,28HSPIN STATISTICS NOT CORRECT ,7H 2*IS1=,I3,7H 2*IS2=
     1,I3,7H 2*IS3=,I3,7H   LTR=,I3,7H 2*JTR=,I3)
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2),I4,2(1H.,I2))
      END

c***********************************************************************
      subroutine ffprint(ii,k,alpha,ub,drf,ltrt,istrt,jtrt,idat)
C
C     Prints form factor
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      dimension ub(800),idat(6),alpha(15)
c
      if(ii.eq.1) then
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
      endif
      R=DRF
      WRITE(6,9100)LTRT,ISTRT,JTRT
      WRITE(6,9101)
      DO 3006 M=1,K,5
      MK=M+M-1
      MK4=MIN0(MK+9,K+K)
      WRITE(6,9052)R ,(UB(N),N=MK,MK4)
      R=R+5.0*DRF
 3006 CONTINUE
      return
 9052 FORMAT(1H ,F6.2,1P10E12.4)
 9100 FORMAT(12H0FORM FACTOR,6X,9H    LTR =,I4,5X,9H   2*STR=,I4,5X
     1,                         9H   2*JTR=,I4)
 9101 FORMAT(55H0    R    RL,R        IM,R        RL,R+DR*1   IM,R+DR*1
     1,             48H   RL,R+DR*2   IM,R+DR*2   RL,R+DR*3   IM,R+DR*3
     2,             24H   RL,R+DR*4   IM,R+DR*4 )
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      end

c***********************************************************************
      SUBROUTINE FORMF(U,V,N,LAM,IK,IVB)
c
c     Processes kinematic input for waves and form factors and
c     calculates potentials or form factors
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC
      COMMON/POTTER/DRX,AFACT(2),VFACT,SFACT,E,RM,G(4),ETAX,FKX,ETAKX
     1             ,RCX,HBARC2,ABETA(3),FLDF(3)
     2             ,NX,LAMX,KMXX,KX,IBX,LPLUSX,ICON4,NSPC,IDIRAC,ICHK
      DIMENSION U(810),V(800),X(4)
      EQUIVALENCE (X,FN)
C
C     READ IN POTENTIAL CARDS FOR CARD SETS 5,6,OR 7
C
      ETA6=60.
      KMAX=400
      ICON4=ICON(4)
      NSPC=N
      IVB=0
      KB=K
      READ (5,9000)E,FM(N),Z(N),FMA(N),ZA(N),RY,AC(N),PNLOC(N),FS(N),QCD
      E=E+QCD
      IF(N.EQ.2) THEN
        Q=E
        E=(ECM(1)+Q)*(FM(2)+FMA(2))/FMA(2)
      ENDIF
      IS(N)=FS(N)
      NS(N)=IS(N)+1
      IF(AMASS.EQ.0.0) AMASS=FMA(1)
      IF(IK.EQ.0) DR(N)=DRF*AMASS/FMA(N)
c     photon does not have recoil factor
      if(is(2).eq.2.and.fm(2).eq.0.0) dr(n)=drf
      if(is(1).eq.1.and.fm(1).eq.0.0) dr(n)=drf
      KMXX=KMAX
      DRX=DR(N)
      AFACT(1)=FMA(N)**.333333333
      AFACT(2)=FM (N)**.333333333
      RC(N)=ABS(RY)*AFACT(1)
      IF(RY.LT.0.0) RC(N)=RC(N)+ABS(RY)*AFACT(2)
      DO 12 M=1,KMAX
      MK=M+M-1
      U(MK  )=0.0
      U(MK+1)=0.0
      V(MK  )=0.0
      V(MK+1)=0.0
   12 CONTINUE
      RM=0.0
      IF(E.EQ.0.0) GO TO 66
C
C        ICON(10).NE.0   GIVES RELATIVISTIC KINEMATICS
C
      IF(ICON(10).NE.0) THEN
        IF(N.NE.2) GO TO 26
        IF(QCD.GT.0.0) GO TO 26
        E=E+(ECM(1)+Q)**2/(2.0*FMA(2)*AMU)
   26   CONTINUE
        FM1=FM(N)*AMU
        FM2=FMA(N)*AMU
        FMT=FM1+FM2
C     KLUDGE FAKE KE-LAB FOR BOUND STATES
        IF(N.GE.3) E=((E+FMT)**2-FMT**2)/(2.0*FM2)
        WLAB=E+FMT
        WCM =SQRT(2.0*E*FM2+FMT**2)
        GAMM=WLAB/WCM
        W1=GAMM*(FMT*FM1+FM2*E)/WLAB
        W2=GAMM* FM2
        ECM(N)=WCM-FMT
        IF(ICON(10).EQ.1) THEN
C       OLD KINEMATICS  SELECTED BY ICON(10) = 1
          FMU(N)=W1*W2/(W1+W2)
          VFACT=2.0*FMU(N)/HBARC**2
          SFACT=2.0*fm(n)*fma(n)/(fm(n)+fma(n))/hbarc**2
          EFACT=VFACT
          ELSE
C       NEW KINEMATICS  SELECTED BY ICON(10) > 1
          FMU(N)=W1
          VFACT= 2.0*W1 /HBARC**2
          SFACT= 2.0*FM1/HBARC**2
          EFACT=VFACT
        ENDIF
        FMU(N)=FMU(N)/AMU
        FK2(N)=(W1**2-FM1**2)/HBARC**2
      ELSE
C     NON RELATIVISTIC KINEMATICS
        FMU(N)=FM(N)*(FMA(N)/(FM(N)+FMA(N)))
        IF(N.LT.3) ECM(N)=E*(FMA(N)/(FM(N)+FMA(N)))
        IF(N.GE.3) ECM(N)=E
        VFACT=2.0*FMU(N)*AMU/HBARC**2
        SFACT=VFACT
        EFACT=VFACT
        FK2(N)=SFACT*ECM(N)
      ENDIF
      FK(N)=SQRT(ABS(FK2(N)))
      ETAK=CHSQ*Z(N)*ZA(N)*EFACT
      ETA(N)=ETAK*0.5/FK(N)
      HBARC2=HBARC**2
C
C     ADD COULOMB AND KINETIC ENERGIES TO U
C
      RCX=RC(N)
      IF(RCX.EQ.0.0) RCX=DR(N)
      R=0.0
      FCOU=0.5*ETAK/RCX
      DO 42 M=1,400
      MK=M+M-1
      R=R+DR(N)
      IF(R.GT.RCX) GO TO 40
      FC=FCOU*(3.0-(R/RCX)**2)
      GO TO 41
   40 CONTINUE
      FC=ETAK/R
   41 CONTINUE
      IF(N.NE.3) U(MK  )=U(MK  )+FK2(N)-FC
      IF(N.EQ.3) U(MK+1)=U(MK+1)+FK2(N)-FC
   42 CONTINUE
      GO TO 67
   66 CONTINUE
      FK(N)=0.0
      ETA(N)=0.0
      FK2(N)=0.0
      ECM(N)=0.0
      FMU(N)=FM(N)*FMA(N)/(FM(N)+FMA(N))
      VFACT=2.0*FMU(N)*AMU/HBARC**2
      SFACT=VFACT
      ETAK=CHSQ*Z(N)*ZA(N)*VFACT
   67 CONTINUE
      IF(N.LT.3) GO TO 68
      IF(ICON(4).EQ.2) GO TO 69
   68 CONTINUE
      IF(N.GT.2) THEN
        Q=ECM(N)
      ELSE
        Q=ECM(N)-ECM(1)
      ENDIF
      WRITE(6,9010)N
      WRITE(6,9503)E,RY,AC(N),FS(N)
      WRITE(6,9504)FM(N),FMA(N),Q
      WRITE(6,9505)Z(N),ZA(N),PNLOC(N)
      WRITE(6,9500)
      RHO=FK(N)*RC(N)
      WRITE(6,9506)ECM(N),RC(N),RHO
      WRITE(6,9507)FK(N),ETA(N),DR(N)
      WRITE(6,9008)
   69 CONTINUE
      FS(N)=FS(N)/2.0
C     Set up variables for potential routine
      IBX=0
      ETAX=ETA(N)
      FKX=FK(N)
      ETAKX=ETAK
      RCX=RC(N)
      LAMX=LAM
      NX=N
      LPLUSX=LPLUS
      KX=K
      ICHK=0
      IDIRAC=0
C
      CALL POTS(U,V)
C
      DR(N)=DRX
      LPLUS=LPLUSX
      IBF(4)=IBX
C     Set nonlocality for Dirac-Darwin term
      IF(N.LE.2.AND.IDIRAC.NE.0.AND.PNLOC(N).EQ.0.0) PNLOC(N)=-1.0
      K=MIN0(MAX0(K,KX),KMAX)
      IF(N.LE.2) GO TO 3000
      IF(E.NE.0.0) GO TO 2000
      IF(IBF(4).NE.0) IVB=IBX
      IF(IVB.EQ.4) THEN
        U(801)=U(2*KC-1)*FLOAT(KC)*DR(3)
        U(802)=0.0
        U(803)=0.0
        U(804)=1.0
      ENDIF
      GO TO 3000
C
C     SINGLE PARTICLE ORBITAL
C
C
C     READ IN QUANTUM NUMBERS FOR SINGLE PARTICLE ORBITAL
C
 2000 CONTINUE
C     Set flags for unbound stripping case
      IF(E.GT.0.0) THEN
        IVB=2
        K=KC
      ENDIF
      READ(5,9000)G,VTRIAL,FISW
      FN =G(1)
      FL =G(2)
      FJ2=G(3)
      FSS=G(4)
      ISW=FISW
      WRITE(6,9500)
      FJ0=FJ2/2.0
      FS0=FSS/2.0
      IF(VTRIAL.EQ.0.0) VTRIAL=ETA6
      WRITE(6,9508)G,FISW
      FACT=(FJ0**2+FJ0-FL**2-FL-FS0**2-FS0)*0.5
      DO 2028 M=1,400
      MK=M+M-1
      U(MK  )=U(MK  )+V(MK  )*FACT
      U(MK+1)=U(MK+1)+V(MK+1)*FACT
 2028 CONTINUE
      WRITE(6,9500)
      CALL BIND(U,  V     ,DR(3),RM,FN,FL,K,FK(3),ETA(3),VTRIAL,ECM(3)
     1,FK2(3),ISW,IBF(3),U(801))
C
      IBF(2)=RM/DR(3)
      DO 2050 M=1,400
      MK=M+M-1
      V(M+400)=VTRIAL*U(MK  )+U(MK+1)
 2050 CONTINUE
      SUM=1.0
      FACT=PNLOC(3)**2/8.0
      IF(FACT.EQ.0.0) GO TO 2076
C
C     NON-LOCAL CORRECTION FOR SINGLE PARTICLE FUNCTION
C
      SUM=0.0
      R=0.0
      DO 2075 M=1,K
      MK=M+M-1
      R=R+DR(3)
      V(M)=V(M)*EXP(FACT*(FK2(3)-V(M+400)))
      SUM=SUM+(V(M)*R)**2
 2075 CONTINUE
      SUM=1.0/SQRT(SUM*DR(3))
C
C     DO NOT RENORMALIZE FOR POSITIVE ENERGY
C
      IF(FK2(3).GT.0.0) SUM=1.0
 2076 CONTINUE
      DO 2100 M=1,K
      MK=M+M-1
      IF(M.GT.K) V(M)=0.0
      U(MK  )=V(M)*SUM
      U(MK+1)=0.0
      V(M)=V(M+400)
 2100 CONTINUE
 3000 CONTINUE
      DO 3020 M=1,4
      X(M)=G(M)
 3020 CONTINUE
      RETURN
 9000 FORMAT(10F8.4)
 9008 FORMAT(21H0POTENTIAL PARAMETERS )
 9010 FORMAT( 9H0PARTICLE,I2,120(1H*))
 9500 FORMAT(1H ,3A6,5(3X,A6,F9.4))
 9503 FORMAT(18H INPUT DATA       ,9H   ELAB =,F9.4,9H   RC0  =,F9.4
     1           ,9H   AC   =,F9.4,9H   2*STR=,F9.4)
 9504 FORMAT(18X ,9H   MASSP=,F9.4,9H   MASST=,F9.4,9H   Q    =,F9.4)
 9505 FORMAT(18X ,9H   ZP   =,F9.4,9H   ZT   =,F9.4,9H   PNLOC=,F9.4)
 9506 FORMAT(18H DERIVED DATA     ,9H   ECM  =,F9.4,9H   RC   =,F9.4
     1           ,9H   RHO  =,F9.4)
 9507 FORMAT(18X ,9H   K    =,F9.4,9H   ETA  =,F9.4,9H   DR   =,F9.4)
 9508 FORMAT(18X ,9H   NODES=,F9.4,9H   L    =,F9.4,9H   2*J  =,F9.4
     1           ,9H   2*S  =,F9.4,9H   FISW =,F9.4)
 9511 FORMAT(18X ,9H   DAMP =,F9.4)
      END

c***********************************************************************
      SUBROUTINE FNLOC(U,V,W,VB,LTR,BKIN,FNRNG,KT)
c
c     Calculates non-local correction factor for distorted waves
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(800,2),V(800,2),W(800),VB(800),BKIN(3,14)
     1,RT(8,7),FM(3),FMA(3),RC(3),PNLOC(3),FMU(3),FK(3),FK2(3)
     1  ,ETA(3),DR(3)
      DATA  RT/  0.0   ,  4.4934,  7.7253, 10.9041, 14.0662, 17.2208,
     1 20.3713, 23.5195,  2.0816,  5.9404,  9.2058, 12.4044, 15.5792,
     2 18.7426, 21.8997, 25.0528,  3.3421,  7.2899, 10.6139, 13.8461,
     3 17.0429, 20.2219, 23.3905, 26.5526,  4.5141,  8.5838, 11.9727,
     4 15.2445, 18.4681, 21.6666, 24.8501, 28.0239,  5.6467,  9.8404,
     5 13.2956, 16.6093, 19.8624, 23.0828, 26.2833, 29.4706,  6.7565,
     6 11.0702, 14.5906, 17.9472, 21.2311, 24.4748, 27.6937, 30.8960,
     7  7.8511, 12.2793, 15.8632, 19.2627, 22.5781, 25.8461, 29.0843,
     8 32.3025/
C
C      FLOAT(III)=DFLOAT(III)
C
C     FINITE RANGE CORRECTION
C
C
C     POSITIVE FR PARAMETER CALCULATES A HULTHEN TYPE FR CORRECTION
C     NEGATIVE FR PARAMETER CALCULATES A EXPONENTIAL TYPE FR CORRECTION
C
      K=KT
      DO 20 I=1,3
      FM(I)   =BKIN(I,3)
      FMA(I)  =BKIN(I,4)
      RC(I)   =BKIN(I,5)
      PNLOC(I)=BKIN(I,7)
      FK(I)   =BKIN(I,10)
      FK2(I)  =BKIN(I,11)
      ETA(I)  =BKIN(I,12)
      DR(I)   =BKIN(I,13)
      FMU(I)  =BKIN(I,14)
   20 CONTINUE
      IF(FNRNG.EQ.0.0) GO TO 398
      IF(FM(1)-FM(2))143,398,142
  142 I1=1
      I2=2
      GO TO 144
  143 I1=2
      I2=1
C
C     I1 IS FOR LARGER MASS , I2 IS FOR SMALLER MASS
C
  144 CONTINUE
      FM1=FMU(I1)
      FM2=FMU(I2)
      TEMP=FM(I2)*FNRNG**2/FM(I1)
      FACT=TEMP*FMA(I1)/FMA(I2)
      FMX=FM(I1)-FM(I2)
      TEMP=TEMP*FMX
      IF(FMX.LT.1.9) GO TO 350
C
C     FORM FACTOR PART BY ROST-KUNZ METHOD
C     THIS OPTION CALCULATES AN EXPONENTIAL TYPE FR CORRECTION
C     FOR MULTI -NUCLEON TRANSFER REACTIONS
C
      DR2=DR(3)**2
      FLF=LTR*(LTR+1)
      R3=DR(3)*FLOAT(K+1)
      D3=0.0
      A3=0.0
      KM=K+K+1
      R2=R3-DR(3)
      D2=-FLF/R3
      A2=W (KM-2)*R2
      DO 220 M=1,K
      KM=KM-2
      R1=R2-DR(3)
      A1=W (KM-2)*R1
      D1=-FLF/R2**2
      IF(A3.EQ.0.0) GO TO 215
      D1=D1+(A3-2.0*A2+A1)/(A2*DR2)
      W (KM+2)=W (KM+2)*EXP(D2*FACT)
      DX2=(D3*A3-2.0*D2*A2+D1*A1)/(R2*DR2)-A2*D2**2/R2
      W (KM+2)=W (KM+2)+0.5*DX2*FACT**2
  215 CONTINUE
      IF(ABS(W (KM-2)).LT.ABS(W (KM  ))) GO TO 225
      R3=R2
      R2=R1
      D3=D2
      D2=D1
      A3=A2
      A2=A1
  220 CONTINUE
  225 CONTINUE
      KM=KM/2+1
      DO 227 M=1,KM
      MK=M+M-1
      VB(M)=W(MK)
      W(MK)=0.0
  227 CONTINUE
      RKM=FLOAT(KM)*DR(3)
      LL=MIN0(LTR,6)
      FLF=LL*(LL+1)
      DR12=DR(3)**2/12.0
      DO 250 N=1,8
      FKX=(RT(N,LL+1)/RKM)**2
      A1=0.0
      D1=0.0
      R=DR(3)
      A2=DR(3)
      D2=(1.0+DR12*(FKX-FLF/R**2))*A2
      FNORM=0.0
      AN=0.0
      DO 230 M=1,KM
      AN=AN+VB(M)*A2*R
      FNORM=FNORM+A2**2
      VB(M+400)=A2/R
      R=R+DR(3)
      D3=1.0+DR12*(FKX-FLF/R**2)
      A3=(12.0*A2-10.0*D2-D1)/D3
      D1=D2
      D2=D3*A3
      A1=A2
      A2=A3
  230 CONTINUE
      AN=AN-0.5*VB(KM)*A1*(R-DR(3))
      FNORM=FNORM-0.5*A1**2
      FNORM=AN*EXP(-FACT*FKX)/FNORM
      DO 240 M=1,KM
      MK=M+M-1
      W (MK)=W (MK)+VB(M+400)*FNORM
  240 CONTINUE
  250 CONTINUE
  340 CONTINUE
C
C     DISTORTED WAVE PART, GAUSSIAN FACTOR
C
      DO 345 M=1,K
      MK=M+M-1
      CTEMP1=TEMP*(U(MK  ,I2)/FM2-U(MK  ,I1)/FM1)
      CTEMP2=TEMP*(U(MK+1,I2)/FM2-U(MK+1,I1)/FM1)
      FACT=EXP(-CTEMP1)
      UT1= FACT*COS(CTEMP2)
      UT2=-FACT*SIN(CTEMP2)
      UF1   =UT1*W(MK  )-UT2*W(MK+1)
      UF2   =UT1*W(MK+1)+UT2*W(MK  )
      W(MK  )=UF1
      W(MK+1)=UF2
  345 CONTINUE
  348 CONTINUE
      GO TO 398
  350 CONTINUE
      IF(FNRNG.GT.0.0) GO TO 360
C
C     SINGLE PARTICLE PART, GAUSSIAN FACTOR
C
      DO 355 M=1,K
      MK=M+M-1
      FACT=EXP(-TEMP*VB(M)/FMU(3))
      W(MK  )=W(MK  )*FACT
      W(MK+1)=W(MK+1)*FACT
  355 CONTINUE
      GO TO 340
  360 CONTINUE
C
C     HULTHEN TYPE
C
      DO 365 M=1,K
      MK=M+M-1
      CTEMP1=VB(M)/FMU(3)+U(MK  ,I2)/FM2-U(MK  ,I1)/FM1
      CTEMP2=             U(MK+1,I2)/FM2-U(MK+1,I1)/FM1
      CTEMP1=1.0+TEMP*CTEMP1
      CTEMP2=    TEMP*CTEMP2
      DET=CTEMP1**2+CTEMP2**2
      UT1=(CTEMP1*W(MK  )+CTEMP2*W(MK+1))/DET
      UT2=(CTEMP1*W(MK+1)-CTEMP2*W(MK  ))/DET
      W(MK  )=UT1
      W(MK+1)=UT2
  365 CONTINUE
  398 CONTINUE
C
C     NON LOCALITY CORRECTION FACTOR FOR DISTORTED WAVES
C
      DO 410 I=1,2
      CALL FNLOC5(U(1,I),V(1,I),W(1),PNLOC(I),FK2(I),ETA(I)
     1,           FK(I),RC(I),DR(I),K)
  410 CONTINUE
      RETURN
      END

c***********************************************************************
      SUBROUTINE CATHEN(LOPT,JOPT,UB,VB,SI)
c
c     Two particle form factor subroutine, calculates microscopic
c     inelastic or two-nucleon transfer form factors
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC
      DIMENSION UB(800,2),VB(801),SI(800),C(10),QNUM(4,2),G(4),IQN(3)
     1         ,JI(2),JJ(2),FJT(2),CON(8)
      EQUIVALENCE (G,FN),(C(1),CNTROL),(C(2),QCODE),(C(3),FMUV)
     1,(C(4),OPT)
      DATA SQR4PI,SQRTEN,PI/3.54490780,3.162277660,3.141592/
C
      IQFLG=0
      KMFLG=0
      IQN(1)=LOPT
      IQN(3)=JOPT
      KT=399
      IBF(2)=1
      FK2(3)=0.0
  100 CONTINUE
      READ(5,9001)C
      CNTRL=ABS(CNTROL)
      IF(CNTRL.EQ.0.0) GO TO 3010
      IF(QCODE)2020,3000,1100
 1100 CONTINUE
      icode=abs(qcode)
C
C     LOOP FOR TWO ORBITALS
C
      ML=1
      DO 2015 I=1,2
      FJT(1)=C(5)
      FJT(2)=C(I+5)
      IF(I.EQ.2.AND.CNTRL.EQ.1.0) GO TO 2005
      IK=0
      IF(icode.EQ.7.OR.icode.EQ.8) ik=i-1
      G(1)=0.0
      G(2)=0.0
      G(3)=0.0
      G(4)=0.0
C
      CALL FORMF(UB(1,2),VB,3,LOPT,IK,IVB)
      IVB=0
C
      DRX=DR(3)
      IF(KMFLG.EQ.0) KM=IBF(2)
      KMFLG=1
      KT=MIN0(KT,K)
 2005 CONTINUE
      IF(FJT(1).EQ.0.0) FJT(2)=FJ2
C
C     PRINT OUT SINGLE PARTICLE INFORMATION
C
      WRITE(6,9501)
      WRITE(6,9100)FJ2,FJT
      IF(ICON(4).EQ.2) GO TO 2007
      IF(CNTRL.EQ.1.0.AND.I.EQ.2) GO TO 2007
      R=DRX
      K2=K+K
      DO 2006 M=1,K2,20
      MK4=MIN0(M+19,K2)
      WRITE(6,9052)R,(UB(N,2),N=M,MK4,2)
      R=R+10.0*DRX
 2006 CONTINUE
 2007 CONTINUE
      DO 2008 M=1,400
      MK=M+M-1
      UB(ML,1)=UB(MK,2)
      ML=ML+1
 2008 CONTINUE
      DO 2009 M=1,4
      QNUM(M,I)=G(M)
 2009 CONTINUE
      JI(I)=FJT(1)
      JJ(I)=FJT(2)
 2015 CONTINUE
C
 2020 ICODE=ABS(QCODE)
      KT=MIN0(KT,K)
      LTR1=IQN(1)+1
      IQN(2)=IS(3)
      IF(IS(3).EQ.0) IQN(3)=IQN(1)+IQN(1)
      WRITE(6,9051)QCODE,FMUV,OPT
C     ENTER TIME REVERSAL PHASE
      LVR=QNUM(2,1)+QNUM(2,2)
      OPT=OPT*PHASEF((LVR+IQN(1))/2)
      IFLAG=1
      IF(ICODE.LE.10) GO TO 2025
      IFLAG=2
      ICODE=ICODE-10
 2025 CONTINUE
C
      GO TO (2100,2100,2300,2990,2500,2100,2700,2700),ICODE
C
C     ICODE=1  YUKAWA POTENTIAL
C     ICODE=2  COULOMB POTENTIAL
C     ICODE=3  OPEP TENSOR POTENTIAL
C     ICODE=5  TWO PARTICLE TRANSFER
C     ICODE=6  ZERO RANGE KNOCK OUT
C     ICODE=7  NO RECOIL FORM FACTOR
C
 2100 LVR=0
      TOPT=1.0
      GO TO 2350
 2300 ICODE=1
      LVR=2
      TOPT=-SQRTEN
 2350 MNNL=IABS(LTR1-LVR-1)+1
      MXXL=LTR1+LVR
      LTR2=LTR1+LTR1-2
      LVR2=LVR+LVR
      LLX=QNUM(2,1)+QNUM(2,1)
      LPX=QNUM(2,2)+QNUM(2,2)
      JLX=QNUM(3,1)
      JPX=QNUM(3,2)
      IS1=QNUM(4,1)
      IS2=QNUM(4,2)
      OPT=OPT*PHASEF(LPX/2)
      DRX=DR(1)
      DO 2390 LAM =MNNL,MXXL,2
      LAM2 =LAM +LAM -2
      VOPT=TOPT*OPT*SQR4PI*SQRT(FLOAT(LAM2 +1)*FLOAT(IQN(2)+1))
     1 *VCC(LAM2 ,LVR2,LTR2 ,0,0)
     2 *RACAH(IQN(2),IQN(2),LAM2 ,LTR2 ,LVR2,IQN(3))
      RME=0.0
      IF(IQN(2).EQ.0) RME=1.0
      IF(IQN(2).EQ.2) RME=SQRT(FLOAT(IS(1)*(IS(1)+2))
     1*FLOAT(IS1*(IS1+2)))
      RME=RME*SQRT(FLOAT(IS1+1))*PHASEF((IS(1)-IS(2))/2)
     1 *SQRT(FLOAT((LLX+1)*(LAM2+1)))*VCC(LLX,LAM2,LPX,0,0)
     2 *SQRT(FLOAT(JLX+1)*FLOAT(JPX+1)*FLOAT(IQN(3)+1))
       rme =rme*WINEJ(LPX,LLX,LAM2 ,IS2,IS1,IQN(2),JPX,JLX,IQN(3))
     4 *PHASEF((JI(1)+IQN(3)-JLX-JJ(2))/2)
     5 *RACAH(JPX,JJ(2),JLX,JJ(1),JI(1),IQN(3))
     6 *SQRT(FLOAT(JJ(1)+1)*FLOAT(JJ(2)+1))
      VOPT=VOPT*RME*SQRT(FLOAT(LTR2+1)/FLOAT((JJ(1)+1)*(IQN(3)+1)))
     1 *(DR(3)/DR(1))**3
      SL=0.0
      IF(VOPT.EQ.0.0) GO TO 2390
      IF(ICODE.EQ.6 ) GO TO 2365
C
      KMAX=400
      CALL SLATR (KT,KMAX,DRX,VB(  1),LTR1,FMUV,ICODE)
      CALL SLATR (KT,KMAX,DRX,VB(801),LAM ,FMUV,ICODE)
C
C
      CALL RADIN
     1(KT,KMAX,DRX,FMUV,VB,UB(1,1),UB(401,1),SL,VOPT,SI(IFLAG),KM,SK)
C
      GO TO 2380
 2365 CONTINUE
      SCALE=VOPT/SQR4PI**2
      R=0.0
      MK=IFLAG
      DO 2370 M=1,KT
      R=R+DRX
      TEMP=UB(M    ,1)*UB(M+400,1)*SCALE
      SI(MK  )=SI(MK  )+TEMP
      SL      =SL      +TEMP*R**2
      IF(M.EQ.KM) SK=TEMP
      MK=MK+2
 2370 CONTINUE
      SL=SL*DRX
 2380 CONTINUE
 2385 I=LAM -1
      WRITE(6,9102)SL,IQN(1),I,RME
      WRITE(6,9002)KM,SK
 2390 CONTINUE
      GO TO 2990
C
C     HERE FOR TWO NUCLEON TRANSFER
C
 2500 CONTINUE
c
c     c1 = R1 scale
c     c2 = R2 scale
c     c3 = r1 scale
c     c4 = r2 scale
c     c5 = r0, integration scale length for relative coordinate
c     c6 = Pauli flag
c     c7 = order  of gaussian integration
c     c8 = number of integration points
c
      IF(FMUV.EQ.0.0) FMUV=1.7
      CON(1)= 1.0
      CON(2)= 1.0
      CON(3)= 0.5
      CON(4)=-0.5
      CON(5)= 2.0*fmuv
      con(6)= 0.0
      con(7)= 0.5
      con(8)= 0.0
      OPT=OPT*16./PI
      CALL DSTRIP
     1(IQN,DRX,KT,UB(1,1),UB(401,1),SI(IFLAG),QNUM,OPT,KM,SL,CON)
      WRITE(6,9002)KM,SL
      IQFLG=1
      GO TO 2990
 2700 CONTINUE
      T1= MIN (FM(1),FM(2))
      T2= MAX (FM(1),FM(2))
      CON(1)= 1.0
      CON(2)= 0.0
      CON(3)= T1/T2
      CON(4)= 1.0
      CON(5)= 1.0
      CON(7)= 0.0
      CON(8)= 0.0
      IF(ICODE.EQ.8) GO TO 2706
      T3=HBARC**2/(2.0*AMU*FMU(3))
      R=0.0
      DO 2705 M=1,KT
      R=R+DRX
      UB(M+400,1)=UB(M+400,1)*T3*  (FK2(3)-VB(M))
 2705 CONTINUE
 2706 CONTINUE
      M=(KT+KM)/2+400
      T1= LOG(UB(M-1,1)*FLOAT(M-401)/(UB(M  ,1)*FLOAT(M-400)))
      FMUV = DRX/T1
      CON(5)=fmuv
      DO 2710 M=1,KT
      UB(M+400,1)=UB(M+400,1)*EXP(T1*FLOAT(M))*DRX*FLOAT(M)/100.
 2710 CONTINUE
      OPT=OPT*FMUV*FMUV*SQR4PI/2.0
      CALL DSTRIP
     1(IQN,DRX,KT,UB(1,1),UB(401,1),SI(IFLAG),QNUM,OPT,KM,SL,CON)
      WRITE(6,9002)KM,SL
      IQFLG=1
      GO TO 2990
 2990 CONTINUE
 3000 CONTINUE
      IF(CNTROL.GT.0.0) GO TO 100
      IF(IQFLG .EQ.1) KT=FLOAT(KT)-1.50/DRX
      K=KT
 3010 RETURN
 9001 FORMAT(10F8.4)
 9002 FORMAT(15H FORM FACTOR,M=,I3,1PE18.6)
 9051 FORMAT(18H0PARAMETERS       ,9H   QCODE=,F9.4,9H   RANGE=,F9.4
     1 ,9H   VZERO=,F9.4)
 9052 FORMAT(1H ,F6.2,1P10E12.4)
 9100 FORMAT(27H0SINGLE PARTICLE FUNCTIONS ,8HCOUPLING ,6H 2*J1=,F3.0
     1 ,9H    2*J2=,F3.0,9H    2*JI=,F3.0)
 9102 FORMAT(19H VOLUME INTEGRAL = ,F10.4,8H    LTR=,I3,8H    LAM=,I3
     1,8H    RME=  ,F8.4)
 9501 FORMAT(1H ,18HFORM FACTOR DATA  )
      END

c***********************************************************************
      SUBROUTINE BDWUCK4
c
c     Control program for second part of DWUCK4 which integrates the
c     distorted waves, calculates radial integrals, inelastic
c     amplitudes and the crossections.
c
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      logical i_sym(2)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array1/U(800,2),V(800,2),SPACE(4800)
      Common/array2/UB(800),ff(10)

      DIMENSION DTEMP(8000),d(2400),sigplt(200,2),plm(800)
      Equivalence (DTEMP(1),U(1,1))
     1,           (dtemp(3201),d(1)),(dtemp(5601),plm(1))
     2,           (dtemp(7601),sigplt(1))
      DATA FOURPI/12.5663706/
C
      IBUF=IBF(1)
      write(*,'(a)')'   Subroutine INTEG  entered'
      CALL INTEG4
C
C
C     PRINT OUT ELASTIC CROSS-SECTIONS
C
      IF(ICON(6).ne.0.or.ICON(16).ne.0) then
        WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
        WRITE(6,9904)
C
        CALL ELSIG(dtemp,d,plm,sigplt,angle,fk,eta,rsig,alpha
     1            ,idat,is,icon,lplus)
C
      else
        WRITE(6,9002)(RSIG(N),N=1,2)
      endif
   66 CONTINUE
      IF(ICON(12).NE.0) then
c
c     print out radial wave functions
c
c     skip over form factors
      DO 72 M=1,NLTR
      READ(4)
   72 CONTINUE
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
C
      CALL TAPED
C
      endif
c
      do 40 i=1,2
      if(fm(i).eq.fma(i).and.z(i).eq.za(i)) then
        i_sym(i)=.true.
      else
        i_sym(i)=.false.
      endif
   40 continue
c
      RZ=FLOAT(KZ)*DRF
      RMAX=FLOAT(K)*DRF
      JR=NS(1)
      JS=NS(2)
      DO 101 II=1,NLTR
      LTR=LTRT(II)
      MPLUS=LTR+1
C
C     SPACE TO PROPER FORM FACTOR FOR THIS L-TRANSFER
C
      DO 70 IJ=1,NLTR
      IF(IJ.EQ.II) then
        READ (4)UB,FF
      else
        READ (4)
      endif
   70 continue
C
      write(*,'(a)')'   Subroutine RADINT entered'
      CALL RADINT(DTEMP(5601),DTEMP(   1),LTR)
C
      IF(ICON(7).ne.0) then
C
C     WRITE( RADIAL MATRIX ELEMENTS IF ICON(7).NE.0
C
      INCR=LPL2*LTR
      INC=1
      IS1=-IS(1)

      DO 100 I=1,JR
      IS2=-IS(2)
      DO  98 J=1,JS
      WRITE(6,9999)ALPHA,(IDAT(M),M=1,3)
      WRITE(6,9900)IS1,IS2
      WRITE(6,9901)
      DO 90 LL=1,LPLUS
      LM=LL-1
      IND=INC+INCR
      WRITE(6,9902)LM,(DTEMP(IT),DTEMP(IT+1),IT=INC,IND,LPL2)
      INC=INC+2
   90 CONTINUE
      INC=INC+INCR
      IS2=IS2+2
   98 CONTINUE
      IS1=IS1+2
  100 CONTINUE
      endif
  101 CONTINUE
      END FILE 2
      REWIND 2
C
C     RESTORE SPIN STORAGE
C
      IS(1)=IBF(7)
      FS(1)=FLOAT(IBF(7))/2.0
      NS(1)=IBF(7)+1
      JR=NS(1)
      IS(2)=IBF(8)
      FS(2)=FLOAT(IBF(8))/2.0
      NS(2)=IBF(8)+1
      JS=NS(2)
      DO 300 II=1,NLTR
      LTR=LTRT(II)
      JTR=JTRT(II)
      IS(3)=ISTRT(II)
C
C     CALCULATE NORMALIZATION FACTOR
C
c     fact normalizes the cross section
c     flfact normalizes the amplitudes
c
      flfact=1.0
C
      if(fm(2).eq.0.0.and.is(2).eq.2) then
c       (p,gamma) capture normalization
        fact=     2.0*fmu(1)*amu/(hbarc*fk(1))**2
     1            *chsq/fk(2)**2
        flfact=sqrt(float(2*ltr+1))
      else
          FACT=     2.0*FMU(1)/(HBARC*FK(1))**2
     1             *2.0*FMU(2)/(HBARC*FK(2))**2
     2             *AMU**2/FOURPI
        if(abs(fm(1)-fm(2)).GT.0.1) then
C         Stripping normalization factors
          flfact=100.0*SQRT(FLOAT(2*LTR+1)/FLOAT(JTR+1))
          fact=fact*float(jtr+1)
        endif
      endif
c
      FACT=FACT*FK(2)/FK(1)
      FN=(II-1)*ICON(3)*(ICON(3)-1)
C
      CALL BETAFN(DTEMP(   1),DTEMP(IBUF),LTR,JTR,flfact,i_sym)
C
      MPLUS=JTR/2+1
      IF(icon(3).le.1.or.(icon(3).eq.2.and.II.eq.NLTR)) then
      IF(ICON(8).ne.0) then
C
C     WRITE BETA TABLES IF ICON(8).NE.0
C
      I1=MPLUS+MPLUS
      IFACT=I1     *NS(1)*NS(2)
      KT=0
      IS2=-IS(2)
      DO 120 J=1,JS
      IS1=-IS(1)
      DO 116 I=1,JR
      WRITE(6,9999)ALPHA,(IDAT(M),M=1,3)
      WRITE(6,9905)IS2,IS1
      WRITE(6,9903)
      DO 114 LL=1,LPLUS
      LM=LL-1
      IND=KT+LM*IFACT+IBUF
      INDF=IND+I1-1
      WRITE(6,9902)LM,(DTEMP(INDEX),INDEX=IND,INDF)
  114 CONTINUE
      KT=KT+I1
      IS1=IS1+2
  116 CONTINUE
      IS2=IS2+2
  120 CONTINUE
      endif
C
C     SET HEADINGS FOR INELASTIC SIGMA
C
      WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
      WRITE(6,9507)DRF,RZ,RMAX,VCE,FS(1)
      WRITE(6,9508)FNRNG,PNLOC,FS(2)
      TEMP=ECM(2)-ECM(1)
      do 280 jj=1,nltr
      if((icon(3).le.1.and.jj.eq.ii).or.
     1   (icon(3).eq.2.and.ii.eq.nltr)) then
        WRITE(6,9505)LTRT(jj),JTRT(jj),IS(3),TEMP
      endif
  280 continue
C
        write(*,'(a)')'   Subroutine INSIG  entered'
        CALL INSIG(DTEMP(IBUF),DTEMP( 233),JTR,FACT)
C
      endif
  300 CONTINUE
      REWIND 2
      RETURN

 9002 FORMAT(1H0,'REACSIG 1',1PE14.5,45X,'REACSIG 2',1PE14.5)
 9505 FORMAT(18H ANG MOM TRANSFER ,9H   LTR  =,I4,14H.       2*JTR=,I4
     1,14H.       2*STR=,I4,14H.       Q    =,F9.4)
 9507 FORMAT(18H0BASIC DATA       ,9H   DR   =,F9.4,9H   RMIN =,F9.4
     1,                            9H   RMAX =,F9.4,9H   COUEX=,F9.4
     2,                            9H     IS1=,F9.4)
 9508 FORMAT(18X,9H   FNRNG=,F9.4,9H   PNLC1=,F9.4,9H   PNLC2=,F9.4
     1,                           9H   PNLC3=,F9.4,9H     IS2=,F9.4)
 9900 FORMAT(24H0 RADIAL MATRIX ELEMENTS ,9H,  J1=L1+,I2,2H/2
     1,                                   9H,  J2=L2+,I2,2H/2 )
 9901 FORMAT('0 L2  F(L2,/L2-LTR/  )    F(L2,/L2-LTR/+2)'
     1,        '    F(L2,/L2-LTR/+4)    F(L2,/L2-LTR/+6)')
 9902 FORMAT(I4,1P10E10.3/(4X,1P10E10.3))
 9903 FORMAT('0 L2  BETA(L2,0)          BETA(L2,1)          BETA(L2,2)'
     1,  '          BETA(L2,3)          BETA(L2,4)')
 9904 FORMAT(1H0,32HELASTIC SCATTERING CROSS-SECTION )
 9905 FORMAT(1H0,8H    MS2=,I2,2H/2,8H    MS1=,I2,2H/2)
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END

c***********************************************************************
      SUBROUTINE INTEG4
C
c     Subroutine to integrate radial differential equations for the
c     distorted waves.
c
c***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array1/U(800,2),V(800,2),SPACE(4800)
      Common/array2/UB(800),ff(10)
      DIMENSION F(400),FP(400),G(400),GP(400),S(400)
     1         ,F1(800),F2(800),Q1(800),Q2(800),C(800)
     2         ,D (800),X(800),DTEMP(3200)
     3         ,E (12),Q(4),A(2),B(2),CTEMP(2),LM(6)
     4         ,DRR2(2),DR2(2),R(2)
      EQUIVALENCE (U,DTEMP)
     1           ,(space(   1),Q1),(space( 801),Q2),(space(1601),F1)
     2           ,(space(2401),F2),(space( 801),X ),(space(   1),D )
     3           ,(DTEMP(   1),F ),(DTEMP( 401),FP),(DTEMP( 801),G )
     4           ,(DTEMP(1201),GP),(DTEMP(1601),S ),(DTEMP(2001),C )
      DATA ETA3/10.E+00/
C
      IWORD=0
      JT=NS(1)+NS(2)
      NP=LPL2*JT
      I=0
      DO 30 N=1,2
      DR2(N)=DR(N)**2/12.0
      R(N)=0.0
      JS=NS(N)
      DO 29 ISS=1,JS
      I=I+1
      LM(I)=0
   29 CONTINUE
   30 CONTINUE
      DO 40 IQ=1,NP
      F1(IQ)=0.0
      F2(IQ)=0.0
      Q1(IQ)=0.0
      Q2(IQ)=0.0
   40 CONTINUE
C
      DO 100 M=1,K
      MK=M+M-1
      IX=0
      I=0
      DO 90 N=1,2
      R(N)=R(N)+DR(N)
      DRR2(N)=DR2(N)/R(N)**2
      Q(1)=1.0+DR2(N)*U(MK  ,N)
      Q(2)=    DR2(N)*U(MK+1,N)
      LTEMP=2.0*FK(N)*R(N)+ETA3
      LTEMP=MIN0(LTEMP,LPLUS)
      FI=-FS(N)
      JS=NS(N)
      SFACT=FS(N)**2+FS(N)
      DO 89 ISS=1,JS
      I=I+1
      FL=0.0
      DO 80 LL=1,LPLUS
      FJ=FL+FI
      IX1=IX+LL+LL-1
      FLFACT=FL**2+FL
      FACT=DR2(N)*(FJ**2+FJ-FLFACT-SFACT)*0.5
      Q(3     )=Q(1)+FACT*V(MK  ,N)-DRR2(N)*FLFACT
      Q(4     )=Q(2)+FACT*V(MK+1,N)
      IF(LL.LE.LM(I)) GO TO 70
      IF(LTEMP.LT.LL) GO TO 72
      LM(I)=LM(I)+1
      IF(FJ-ABS(FL-FS(N)).LT.0.0) GO TO 72
c       calculate approximate starting value
        f2(ix1  )=1.0
        do 50 ii=1,ll
        f2(ix1  )=f2(ix1  )*(fk(n)*r(n))/float(2*ii-1)
   50   continue
      F2(IX1+1)=0.0
      Q2(IX1  )=Q(3)*f2(ix1  )
      Q2(IX1+1)=Q(4)*f2(ix1  )
C
C     EVALUATE Q AT ORIGIN FOR L=1
C
      IF(LL.EQ.2) Q1(IX+3)=-f2(ix1  )/6.0
      GO TO 72
   70 CONTINUE
c
c     Step equations forward by dr(n) via Numerov-Fox-Goodwin-Milne method
c
      CTEMP(1)=12.*F2(IX1  )-10.*Q2(IX1  )-Q1(IX1  )
      CTEMP(2)=12.*F2(IX1+1)-10.*Q2(IX1+1)-Q1(IX1+1)
      F1(IX1  )=F2(IX1  )
      F1(IX1+1)=F2(IX1+1)
      DET=Q(3)**2+Q(4)**2
      F2(IX1  )=(CTEMP(1)*Q(3     )+CTEMP(2)*Q(4     ))/DET
      F2(IX1+1)=(CTEMP(2)*Q(3     )-CTEMP(1)*Q(4     ))/DET
      Q1(IX1  )=Q2(IX1  )
      Q1(IX1+1)=Q2(IX1+1)
      Q2(IX1  )=CTEMP(1)
      Q2(IX1+1)=CTEMP(2)
   72 CONTINUE
      FL=FL+1.0
   80 CONTINUE
      FI=FI+1.0
      IX=IX+LPL2
   89 CONTINUE
   90 CONTINUE
C
C     WRITE RADIAL WAVE FUNCTIONS ON TAPE 4
C
      WRITE(4)(F2(J),J=1,NP)
  100 CONTINUE
      LX=1
      drrc = 0.1
      DO 120 N=1,2
      R2=FK(N)*R(N)
      R1=R2-DR(N)*FK(N)
      CALL COU(R1,R2,ETA(N),LPLUS,drrc,F(LX),FP(LX),G(LX),GP(LX),S(LX))
      RSIG(N)=0.0
      LX=LX+LPLUS
  120 CONTINUE
C
      IF(ICON(5).ne.0) then
        WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
        WRITE(6,9600)
        WRITE(6,9601)
      endif
c
c     Match solutions to asymptotic form
c
      DO 300 LL=1,LPLUS
      FL=FLOAT(LL-1)
      LX=LL
      I=0
      IX1=LL+LL-1
      DO 200 N=1,2
      JS=NS(N)
      FI=-FS(N)
      ARG=S(LX)-S(LX-LL+1)
      Q(1)=COS(ARG)
      Q(2)=SIN(ARG)
      Q(3)=Q(1)**2-Q(2)**2
      Q(4)=2.0*Q(1)*Q(2)
      DO 199 ISS=1,JS
      FJ=FL+FI
      I=I+1
      DET=F(LX)*GP(LX)-FP(LX)*G(LX)
      A(1)=(F1(IX1  )*GP(LX)-F2(IX1  )*G (LX))/DET
      A(2)=(F1(IX1+1)*GP(LX)-F2(IX1+1)*G (LX))/DET
      B(1)=(F2(IX1  )*F (LX)-F1(IX1  )*FP(LX))/DET
      B(2)=(F2(IX1+1)*F (LX)-F1(IX1+1)*FP(LX))/DET
      IF(LL.LE.LM(I).and.FJ-ABS(FL-FS(N)).ge.0.0) then
        DET=(A(1)+B(2))**2+(A(2)-B(1))**2
        CTEMP(1)=(A(1)+B(2))/DET
        CTEMP(2)=(B(1)-A(2))/DET
      else
          CTEMP(1)=0.0
          CTEMP(2)=0.0
      endif
C
C     C=NORMALIZATION CONSTANTS
C
      C(IX1  )=Q(1)*CTEMP(1)-Q(2)*CTEMP(2)
      C(IX1+1)=Q(1)*CTEMP(2)+Q(2)*CTEMP(1)
C
C     E=PARTIAL WAVE SCATTERING AMPLITUDES
C
      E(I+I-1)=B(1)*CTEMP(1)-B(2)*CTEMP(2)
      E(I+I  )=B(1)*CTEMP(2)+B(2)*CTEMP(1)
C
C     D=PARTIAL WAVE SCATTERING AMPLITUDES* COULOMB PHASE / WAVE NUMBER
C
      D(IX1  )=(Q(3)*E(I+I-1)-Q(4)*E(I+I  ))/FK(N)
      D(IX1+1)=(Q(3)*E(I+I  )+Q(4)*E(I+I-1))/FK(N)
      X(IX1  )=E(I+I-1)
      X(IX1+1)=E(I+I  )
C
C     CALCULATE REACTION SIGMA
C
      RSIG(N)=RSIG(N)+(2.0*FJ+1.0)*(E(I+I  )-E(I+I-1)**2-E(I+I  )**2)
      FI=FI+1.0
      IX1=IX1+LPL2
  199 CONTINUE
      LX=LX+LPLUS
  200 CONTINUE
      IF(ICON(5).ne.0) then
C
C     WRITE ELASTIC PARTIAL WAVE SCATTERING AMPLITUDES
C
        IX=LL-1
        I1=1
        I2=NS(1)*2
        WRITE(6,9602)IX,(E(INDEX  ),INDEX=I1,I2)
        I1=I2+1
        I2=I2+NS(2)*2
        WRITE(6,9603)IX,(E(INDEX  ),INDEX=I1,I2)
      endif
  300 CONTINUE
      DO 310 N=1,2
      RSIG(N)=RSIG(N)*12.566371 /((2.0*FS(N)+1.0)*FK(N)**2)
  310 CONTINUE
C
C     WRITE NORMALIZATION CONSTANTS ON TAPE 4
C
      WRITE(4)(C(I),I=1,NP)
C
C     WRITE END POINT FUNCTIONS FOR UNBOUND STATE STRIPPING
C
C     ILINE IS SET IN ADWUCK FORM FACTOR LOOP
      IF(ILINE.EQ.0) GO TO 400
      WRITE(4)(X(I),I=1,NP),(GP(I),FP(I),I=1,LPL2),(S(I),I=1,LPL2)
  400 CONTINUE
      END FILE 4
      REWIND 4
      RETURN
 9600 FORMAT(1H0,35HPARTIAL WAVE SCATTERING AMPLITUDES )
 9601 FORMAT(4H   L,20H   REAL D1   IMAG D1 ,20H   REAL D2   IMAG D2
     1,             20H   REAL D3   IMAG D3
     2,   4X,4H   L,20H   REAL D1   IMAG D1 ,20H   REAL D2   IMAG D2
     3,             20H   REAL D3   IMAG D3 )
 9602 FORMAT(1H ,    I3,6F10.6)
 9603 FORMAT(1H+,68X,I3,6F10.6)
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END

c***********************************************************************
      SUBROUTINE TAPED
C
c     Subroutine writes out radial wave functions,
c     spaced icon(12) points apart
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array1/U(800,2),V(800,2),SPACE(4800)

      DIMENSION C(800),FR(800),R(2),CTEMP(2)
      EQUIVALENCE (C(1),V(401,1)),(FR(1),U(1,2))
c
      IK=ICON(12)
      JT=NS(1)+NS(2)
      NP=LPL2*JT
      R(1)=0.0
      R(2)=0.0
      IKTEMP=0
      DO 50 M=1,K
      IKTEMP=IKTEMP+1
C
C     READ IN DISTORTED WAVE RADIAL FUNCTIONS
C
      READ (4)(FR(J),J=1,NP)
      IF(IKTEMP.NE.IK) GO TO 50
      IKTEMP=0
      IX=0
      DO 40 N=1,2
      R(N)=R(N)+DR(N)*FLOAT(IK)
      JX=NS(N)
      DO 39 J=1,JX
      DO 30 LL=1,LPLUS
      LK=LL+LL-1
      IX1=LK+IX
      IX2=IX1+1
C
C     NORMALIZE RADIAL FUNCTIONS
C
      CTEMP(1)=FR(IX1)*C(IX1)-FR(IX2)*C(IX2)
      CTEMP(2)=FR(IX1)*C(IX2)+FR(IX2)*C(IX1)
      FR(IX1)=CTEMP(1)
      FR(IX2)=CTEMP(2)
   30 CONTINUE
      IX=IX+LPL2
   39 CONTINUE
   40 CONTINUE
      WRITE(6,9001)R(1),R(2)
      WRITE(6,9601)
      DO 45 LL=1,LPLUS
      LM=LL-1
      LK1=LM+LL
      LK2=LPL2*(NS(1)-1)+LK1
      WRITE(6,9602)LM,(FR(LK),FR(LK+1),LK=LK1,LK2,LPL2)
      LK1=LK2+LPL2
      LK2=LPL2*(NS(2)-1)+LK1
      WRITE(6,9603)LM,(FR(LK),FR(LK+1),LK=LK1,LK2,LPL2)
   45 CONTINUE
   50 CONTINUE
      REWIND 4
      RETURN
 9001 FORMAT(1H0,3HR1=,F8.4,57X,3HR2=,F8.4)
 9601 FORMAT(4H   L,20H  REAL D1   IMAG D1 ,20H  REAL D2   IMAG D2 ,20H
     1 REAL D3   IMAG D3 ,4X,4H   L,20H  REAL D1   IMAG D1 ,20H  REAL D2
     2   IMAG D2 ,20H  REAL D3   IMAG D3 )
 9602 FORMAT(1H ,    I3,6F10.7)
 9603 FORMAT(1H+,68X,I3,6F10.7)
      END

c***********************************************************************
      SUBROUTINE RADINT(F,FLL,LTR)
c
c     Subroutine for computing radial integrals
c
c     f     Radial wave function storage
c     fll   Radial integral storage
c     ltr   Angular momentum transfer
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array1/U(800,2),V(800,2),SPACE(3200),BUFFR(1600)
      Common/array2/UB(800),ff(10)
      DIMENSION H(800), F(800),FLL(5600),C(800),CTEMP(2),VTEMP(2),DRX(2)
     1         ,D(800),FP(204),FA(204),F1(4),LPTBL(200),S(400)
     2         ,FFTEMP(2),WT(5)
      EQUIVALENCE (C,UB),(D(1),SPACE(2401)),(H(1),UB(1))
C
c     Xmax   maximum X = ikr off real axis
c     Nx     number of points off real axis, must be a multiple of 5
c     lpt0   starting position in buffr
C
      DATA XMAX,NX,LPT0/10.0, 50, 51/
C
      IF(RC(3).EQ.0.0) RC(3)=DR(3)
c     Coulomb excitation factor
      VFACT=3.0*VCE*CHSQ*Z(3)*ZA(3)*RC(3)**LTR/FLOAT(LTR+LTR+1)
      if(ltr.eq.0) vfact=vfact/3.0
      NMTR=LTR+1
      JR=NS(1)
      JS=NS(2)
      DRX(1)=2.0*DR(3)/3.0
      DRX(2)=2.0*DRX(1)
C     JACOBIAN FACTOR
      ANORM= MAX (DR(1),DR(2))**2/(DR(3)* MIN (DR(1),DR(2)))
      NP=LPL2*(JR+JS)
      IZZ=NS(1)*LPL2
      INCR=LPL2*NMTR
      INC=INCR*NS(1)*NS(2)
C     Unbound state and coulomb excitation end point data
C     ff(1) = FF REAL * RMAX
C     ff(2) = FF IMAG * RMAX
C     ff(3) = D  REAL
C     ff(4) = D  IMAG
      DO 2 M=1,4
    2 F1(M) =ff(M)
      IVB   =ff(5)
      FK2(3)=ff(6)
      ETA(3)=ff(7)
      if(ivb.eq.2) then
c       Unbound state
        FFTEMP(1)=UB(2*KC-1)
        FFTEMP(2)=UB(2*KC  )
      elseif(ivb.eq.4) then
c       Coulomb excitation
        fftemp(1)=f1(1)/(dr(3)*float(k))
        fftemp(2)=f1(2)/(dr(3)*float(k))
      endif
c     Spin orbit factor
      FSFACT=FS(1)*(FS(1)+1.)+FS(2)*(FS(2)+1.)
c     clear radial integral storage
      DO 3 II=1,INC
      FLL(II)=0.0
    3 CONTINUE
c
      MX=1
      R=0.0
      DO 100 M=1,K
      R=R+DRF
C
C     READ IN DISTORTED WAVE RADIAL FUNCTIONS
C
      READ (4)(F(J),J=1,NP)
      MK=M+M-1
      MX=3-MX
      VTEMP(1)=UB(MK  )*DRX(MX)
      VTEMP(2)=UB(MK+1)*DRX(MX)
      IF(VFACT.NE.0.0) then
C
C     COULOMB EXCITATION ADDITION TO FORM FACTOR
C
      IF(R.GT.RC(3)) then
      VTEMP(1)=VTEMP(1)+DRX(MX)*VFACT/R**NMTR
      endif
      endif
C
C**** READ IN NORMALIZATION CONSTANTS at last point
C
      IF(M.NE.K) GO TO 45
      READ (4)(C(I),I=1,NP)
      VTEMP(1)=VTEMP(1)*0.5
      VTEMP(2)=VTEMP(2)*0.5
   45 INC=0
      IY=0
      FS1=-FS(1)
      DO 96 I= 1,JR
      IZ=IZZ
      FS2=-FS(2)
      DO 95 J=1,JS
      DO 90 LL=1,LPLUS
      LK=LL+LL-1
      IZ1=LK+IZ
      LP1=IABS(LL-LTR-1)+1
      LP2=MIN0(LL+LTR,LPLUS)
      CTEMP(1)=VTEMP(1)*F(IZ1  )-VTEMP(2)*F(IZ1+1)
      CTEMP(2)=VTEMP(1)*F(IZ1+1)+VTEMP(2)*F(IZ1  )
      DO 80 LP=LP1,LP2,2
      IY1=LP+LP-1+IY
      IND=LK+INC
      IF(M.LT.KZ) GO TO 79
      FLL(IND  )=FLL(IND  )+CTEMP(1)*F(IY1  )-CTEMP(2)*F(IY1+1)
      FLL(IND+1)=FLL(IND+1)+CTEMP(1)*F(IY1+1)+CTEMP(2)*F(IY1  )
      IF(M.EQ.K) THEN
C
C     NORMALIZE RADIAL INTEGRALS
C
      CTEMP1=C(IZ1  )*FLL(IND  )-C(IZ1+1)*FLL(IND+1)
      CTEMP2=C(IZ1  )*FLL(IND+1)+C(IZ1+1)*FLL(IND  )
      FLL(IND  )=(CTEMP1*C(IY1  )-CTEMP2*C(IY1+1))*ANORM
      FLL(IND+1)=(CTEMP1*C(IY1+1)+CTEMP2*C(IY1  ))*ANORM
C     L*S FORM FACTOR,  use symmetrized form
      IF(IVB.EQ.1) THEN
        FJ1=FLOAT(LL)+FS1
        FJ2=FLOAT(LP)+FS2
        TEMP=(FJ1*(FJ1-1.)+FJ2*(FJ2-1.)-FLOAT(LL*(LL-1))
     1        -FLOAT(LP*(LP-1))-FSFACT)/4.
        FLL(IND  )=FLL(IND  )*TEMP
        FLL(IND+1)=FLL(IND+1)*TEMP
      ENDIF
      ENDIF
   79 LK=LK+LPL2
   80 CONTINUE
   90 CONTINUE
      INC=INC+INCR
      IZ=IZ+LPL2
      FS2=FS2+1.0
   95 CONTINUE
      IY=IY+LPL2
      FS1=FS1+1.0
   96 CONTINUE
  100 CONTINUE
      IF(ICON(20).NE.0) GO TO 1100
c     special flag to skip this
      IF(IVB.EQ.2.OR.IVB.EQ.4) THEN
C
C     UNBOUND STATE SECTION
C
C     read in scatt. amp. and asymptotic fcts - h(ODD) = g, h(EVEN) = f
c     and coulomb phases
      READ (4)(D(I),I=1,NP),(H(2*I-1),H(2*I),I=1,LPL2),(S(I),I=1,LPL2)
      R1=DR(1)*FLOAT(KC)
      R2=DR(2)*FLOAT(KC)
      R3=DR(3)*FLOAT(KC)
      DX=XMAX/FLOAT(NX)
      EZ=EXP(-DX)
      DX1=DX/(FK(1)-FK(2)*FMA(1)/FMA(2)-FK(3))
      IF(DX1.LE.0.0) THEN
        WRITE(6,'(a)')'0KINEMATICS INCONSISTENT - EXIT ON ERROR'
        STOP 'KINEMATICS ERROR'
      ENDIF
C
      NI=DX1/DR(1)/2.5+0.5
      NI=MAX0(NI,1)
      DX1=DX1/FLOAT(NI)
      DX2=DX1*FMA(1)/FMA(2)
C     CALCULATE FORM FACTOR extension into complex R plane  Rmax + iY
      IFLG=1
      CALL UNBIND(LTR,NX,NI,FA,F1(1),F1(3),ETA(3),FK2(3),DX1,R1,IFLG)
C
C     ADJUST FA FOR FINITE RANGE, NON LOCALITY, AND PHASE
C     FFTEMP IS FORM FACTOR AT RMAX
C     FA IS UNMODIFIED FORM FACTOR FROM BIND AT RMAX
C     (T1,T2) = FFTEMP/FA
      T5=FA(2*NX+1)**2+FA(2*NX+2)**2
      T1=(FFTEMP(1)*FA(2*NX+1)+FFTEMP(2)*FA(2*NX+2))/T5
      T2=(FFTEMP(2)*FA(2*NX+1)-FFTEMP(1)*FA(2*NX+2))/T5
      DX=DX1*FLOAT(NI)
      X =DX *FLOAT(NX)
      N1=NX+1
      M2=4*NX+5
C     DIVIDE FA BY R+iX AND MULTIPLY BY EXP(-i*DELTA)
      DO 120 M=1,N1
      M1=M+M-1
      M2=M2-2
      T5= R1*(R1*T1+X *T2)/(R1*R1+X*X)
      T6= R1*(R1*T2-X *T1)/(R1*R1+X*X)
      T3      =FA(M1  )*T5 - FA(M1+1)*T6
      FA(M1+1)=FA(M1+1)*T5 + FA(M1  )*T6
      FA(M1  )=T3
      FA(M2  )=T3
      FA(M2+1)=-FA(M1+1)
      X=X-DX
  120 CONTINUE
      M2=4*NX+4
C
      WT(1)=190.*DX/288.
      WT(2)=375.*DX/288.
      WT(3)=250.*DX/288.
      WT(4)=WT(3)
      WT(5)=WT(2)
      M2=4*NX+5
C     MULTIPLY FA BY INTEGRATION WEIGHTS
      DO 130 M=1,N1
      M1=M+M-1
      M2=M2-2
      MX=MOD(M-1,5)
      FA(M1+1)=FA(M1+1)*WT(MX+1)
      FA(M1  )=FA(M1  )*WT(MX+1)
      FA(M2  )=FA(M2  )*WT(MX+1)
      FA(M2+1)=FA(M2+1)*WT(MX+1)
  130 CONTINUE
C     ADJUST WT FOR END POINTS IN INTEGRATION RULE
      FA(M1  )=FA(M1  )*0.5
      FA(M1+1)=FA(M1+1)*0.5
      FA(M2  )=FA(M2  )*0.5
      FA(M2+1)=FA(M2+1)*0.5
C
      IF((2*LTR+1)*(2*NX+2).GT.IBUFF-LPT0) THEN
        WRITE(6,9900)LTR
 9900   FORMAT('0NOT ENOUGH CIRCULAR BUFFER SPACE FOR DEUTERON WAVE,'
     1  ,' LTR=' , I3, ' IS TOO LARGE')
        STOP 'BUFFER ERROR IN UNBIND'
      ENDIF
C
      INC=0
      IY1=0
      FS1=-FS(1)
      DO 1020 I=1,JR
      IY2=IZZ
      FS2=-FS(2)
      DO 1010 J=1,JS
      LPX=0
      LPT=LPT0
      DO 1000 LL=1,LPLUS
      LK=LL+LL-1
      LP1=IABS(LL-LTR-1)+1
      LP2=MIN0(LL+LTR,LPLUS)
      FJ2=FLOAT(LL-1)+FS2
C     CALCULATE PROTON WAVE extension into complex plane
      LM=LK+LPL2
      IZ2=IY2+LK
      IFLG=1
      CALL UNBIND(LL-1,NX,NI,FP(1),D(IZ2),H(LM)
     1,ETA(2),FK2(2),DX2,R2,IFLG)
C
      M2=4*NX+4
      LP=1
      IF(LL.EQ.1) LP=MIN0(LTR+1,LPLUS)
      IF(LPX.GE.LPLUS) GO TO 160
      DO 150 L=1,LP
C     CALCULATE DEUTERON WAVE AND PLACE IN CIRCULAR BUFFER
      LM=LPX+LPX+1
      IZ1=IY1+LM
      IFLG=0
      CALL UNBIND(LPX,NX,NI,BUFFR(LPT),D(IZ1),H(LM)
     1,           ETA(1),FK2(1),DX1,R1,IFLG)
      LPX=LPX+1
      LPTBL(LPX)=LPT-1
      LPT=LPT+(2*NX+2)
      IF(LPT+(2*NX+2).GT.IBUFF) LPT=LPT0
  150 CONTINUE
  160 CONTINUE
C
C     MULTIPLY FP BY FA
C
      LM=LL+LPLUS
      CP=COS(S(LM)-S(LPLUS+1))
      SP=SIN(S(LM)-S(LPLUS+1))
      NX1=2*NX+2
      DO 300 M=1,NX1
      M1=M+M-1
      T3      =FP(M1  )*FA(M1  )-FP(M1+1)*FA(M1+1)
      T4      =FP(M1  )*FA(M1+1)+FP(M1+1)*FA(M1  )
      FP(M1  )=CP*T3-SP*T4
      FP(M1+1)=CP*T4+SP*T3
  300 CONTINUE
      M2=4*NX+4
c
      LPT1=LPTBL(LP)+1
      LPT2=LPT1+(2*NX+2)-1
C
      DO 400 LP=LP1,LP2,2
      FJ1=FLOAT(LP-1)+FS1
      IND=LK+INC
      SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      SUM4=0.0
C     CONSTRUCT OUTGOING WAVE AMPLITUDE FROM D = ( EXP(2i*delta) - 1)/2i
      IZ1=IY1+LP+LP-1
      S1=1.0-2.0*D(IZ1+1)
      S2=    2.0*D(IZ1  )
      CD=COS(S(LP)-S(1))
      SD=SIN(S(LP)-S(1))
      M1=2*NX+3
      M2=2*NX+1
      NX1=NX+1
      IF(FJ1.LT.0.0.OR.FJ2.LT.0.0) GO TO 331
C     NUMERICAL INTEGRATION FOR RADIAL MATRIX ELEMENTS
C     INTEGRATION IS +i*DX FOR +X AND -i*DX FOR -X
C     WITH LIMITS ON X OF 0 --> inf HENCE THE + SIGN
C     WHEN COMBINING OUTGOING AND INGOING DEUTERON WAVES
C
      FE=1.0
      DO 330 M=1,NX1
      M1=M1-2
      M2=M2+2
      IX1=LPTBL(LP)+M1
C     MULTIPLY BY DEUTERON OUTGOING WAVE PART, UD(+)*FP(R+iX)
      SUM1=SUM1+(FP(M1  )*BUFFR(IX1  )-FP(M1+1)*BUFFR(IX1+1))*FE
      SUM2=SUM2+(FP(M1+1)*BUFFR(IX1  )+FP(M1  )*BUFFR(IX1+1))*FE
C     MULTIPLY BY DEUTERON INGOING  WAVE PART, UD(-)*FP(R-iX)
      SUM3=SUM3+(FP(M2  )*BUFFR(IX1  )+FP(M2+1)*BUFFR(IX1+1))*FE
      SUM4=SUM4+(FP(M2+1)*BUFFR(IX1  )-FP(M2  )*BUFFR(IX1+1))*FE
      FE=FE*EZ
  330 CONTINUE
  331 CONTINUE
      L1=LP+LP-1
      L2=LL+LL-1+LPL2
C     UD = (S*UD(+) + UD(-))/2
      T3=(S1*SUM1-S2*SUM2+SUM3)*0.5
      T4=(S1*SUM2+S2*SUM1+SUM4)*0.5
      T1=T3*CD-T4*SD
      T2=T3*SD+T4*CD
      FLL(IND  )=FLL(IND  )+T1*ANORM
      FLL(IND+1)=FLL(IND+1)+T2*ANORM
      LK=LK+LPL2
  400 CONTINUE
 1000 CONTINUE
      INC=INC+INCR
      IY2=IY2+LPL2
      FS1=FS1+1.0
 1010 CONTINUE
      IY1=IY1+LPL2
      FS2=FS2+1.0
 1020 CONTINUE
C
      ENDIF
C
 1100 CONTINUE
C
C
      REWIND 4
      JR=IBF(7)+1
      JS=IBF(8)+1
      IF(NLTR.NE.1) GO TO 1105
      IF(JR+JS.EQ.2) GO TO 2000
 1105 CONTINUE
C
C     WRITE RADIAL MATRIX ELEMENTS ON TAPE 2
C
      JSX=JS
      IF(IBF(6).EQ.0) JSX=1
      DO 1120 I=1,JR
      IY=I-1
      IF(IBF(5).EQ.0) IY=0
      DO 1110 J=1,JS
      IZ=J-1
      IF(IBF(6).EQ.0) IZ=0
      INC=INCR*(JSX*IY+IZ)+1
      INDEX=INC-1+INCR
      WRITE(2)(FLL(II),II=INC,INDEX)
 1110 CONTINUE
 1120 CONTINUE
 2000 CONTINUE
      RETURN
      END

c***********************************************************************
      SUBROUTINE UNBIND(L,NX,NI,Y,D,H,ETA,FK2,DX,R,IFLG)
C
c     Subroutine to extend wave functions into complex plane for unbound
c     stripping cases.
c
c     l      Orbital angular momentum
c     nx     Number of point in integration
c     ni     Number of intermediate points
c     y      Storage for resulting wave function
c     d      Scattering amplitude
c     h      Hankel functions at Rmax
c     eta    Coulomb parameter
c     fk2    Square of wave number
c     dx     Integration step size
c     r      Rmax on real axis
c     iflg   =0 for single side of real axis,  =1 for both sides
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      DIMENSION Y1(2),Z0(2),Z1(2),Y(204),D(2),H(2)
      DATA Y0/1.E+00/
C
      R2=R*R
      FL1=L*(L+1)
      DX12=DX*DX/12.
      X=DX*FLOAT(NX*NI+NI+1)
      Z=DX*SQRT(FK2)
      EZ=EXP(-Z)
c     set starting point
      S1=R2-X*X
      S2=2.*R*X
      ETAK2=2.0*ETA*SQRT(FK2)
      T0   =FL1/(S1*S1+S2+S2)
      T1   =ETAK2/(R2+X*X)
      IF(FK2.EQ.0.0) THEN
        Y0R = 1.0-FLOAT(L)*DX*X/(R2+X*X)
        Y0I =    -FLOAT(L)*DX*R/(R2+X*X)
      ELSE
        Y0R = EZ
        Y0I = 0.0
      ENDIF
      T3   =(( S1*T0+R*T1-FK2)*DX12+1.0)*Y0
      T4   =((-S2*T0-X*T1    )*DX12    )*Y0
      Z0(1)=T3*Y0R - T4*Y0I
      Z0(2)=T4*Y0R + T3*Y0I
C
      X=X-DX
      S1=R2-X*X
      S2=2.0*R*X
      T0=FL1/(S1*S1+S2*S2)
      T1=ETAK2/(R2+X*X)
      Y1(1)=Y0
      Y1(2)=0.0
      Z1(1)=Y1(1)*(( S1*T0+R*T1-FK2)*DX12+1.0)
      Z1(2)=Y1(1)* (-S2*T0-X*T1    )*DX12
C
C     INTEGRATE INWARDS FROM R + iX TO R - iX FOR IFLG = 1
C                       FROM R + iX TO R      FOR IFLG = 0
C     REMOVE EXP( X ) FACTOR SO THERE IS NO OVERFLOW
C
      N2=NX*(IFLG+1)+1
      IND=1
      DO 200 I=1,N2
      DO 190 J=1,NI
      X=X-DX
      S1=R2-X*X
      S2=2.0*R*X
      T0=  FL1/(S1*S1+S2*S2)
      T3=ETAK2/(R2+X*X)
      T1=(12.0*Y1(1)-10.0*Z1(1)-Z0(1)*EZ)*EZ
      T2=(12.0*Y1(2)-10.0*Z1(2)-Z0(2)*EZ)*EZ
      Z0(1)=Z1(1)
      Z0(2)=Z1(2)
      Z1(1)=T1
      Z1(2)=T2
      F1=( S1*T0+R*T3-FK2)*DX12+1.0
      F2=(-S2*T0-X*T3    )*DX12
      Y1(1)=(T1*F1+T2*F2)/(F1*F1+F2*F2)
      Y1(2)=(T2*F1-T1*F2)/(F1*F1+F2*F2)
  190 CONTINUE
      Y(IND  )=Y1(1)
      Y(IND+1)=Y1(2)
      IND=IND+2
  200 CONTINUE
C
C     E1,E2 ARE OUTGOING WAVE AMPLITUDES
C     Z0 IS NORMALIZATION CONSTANT TO AN
C     OUTGOING HANKEL FCT AT Z = RMAX + i0
      IF(IFLG.EQ.0) GO TO 500
C     CONSTRUCT AN OUTGOING WAVE AMPLITUDE FROM E
      F1=1.0
      F2=EXP(-2.0*Z*FLOAT(NI))
      T0=Y(N2)**2+Y(N2+1)**2
C     H(1) IS G
C     H(2) IS F
C     Z0   =  H(+)/Y(R + i0)
      Z0(1)= (H(1)*Y(N2)+H(2)*Y(N2+1))/T0
      Z0(2)= (H(2)*Y(N2)-H(1)*Y(N2+1))/T0
C
C     SPREAD FCT OUT INTO 2 TABLES, FIRST HALF
C     IS FOR OUTGOING WAVE AND SECOND HALF IS
C     INGOING WAVE
      IND=4*NX+3
      IK =2*NX+3
      NX1=NX+1
      DO 220 I=1,NX1
      IND=IND-2
      IK =IK -2
      Y(IK   )=Y(IK   )*F1
      Y(IK +1)=Y(IK +1)*F1
      Y(IND+2)=Y(IND  )
      Y(IND+3)=Y(IND+1)
      F1=F1*F2
  220 CONTINUE
C
      S1=1.0-2.0*D(2)
      S2=2.0*D(1)
      IND=2*NX+1
      IK =2*NX+3
      DO 250 I=1,NX1
      IND=IND+2
      IK =IK -2
C     T1, T2 = exp(-k*X)*H(R-iX)
C     T3, T4 = exp( k*X)*H(R+iX)
      T3=Y(IK   )*Z0(1)-Y(IK +1)*Z0(2)
      T4=Y(IK +1)*Z0(1)+Y(IK   )*Z0(2)
      T1=Y(IND  )*Z0(1)-Y(IND+1)*Z0(2)
      T2=Y(IND+1)*Z0(1)+Y(IND  )*Z0(2)
C     u(R-iX) = exp(-k*X)*(S*H(R-iX) - H(R+iX))/(2i)
C
      Y(IND  )=( T4+S1*T2+S2*T1)*0.5
      Y(IND+1)=( T3-S1*T1+S2*T2)*0.5
C     u(R+iX) = exp(-k*X)*(S*H(R+iX) - H(R-iX))/(2*i)
      Y(IK   )=( T2+S1*T4+S2*T3)*0.5
      Y(IK +1)=( T1-S1*T3+S2*T4)*0.5
  250 CONTINUE
      GO TO 1000
C
C     NORMALIZE 'DEUTERON' WAVE TO OUTGOING HANKEL FUNCTION
C
  500 CONTINUE
      N1=N2
      N2=N2+NX
      T0=Y(N2)**2+Y(N2+1)**2
C     Z0 = H(+)/Y( R + i0)
      Z0(1)=(H(1)*Y(N2)+H(2)*Y(N2+1))/T0
      Z0(2)=(H(2)*Y(N2)-H(1)*Y(N2+1))/T0
      DO 650 I=1,N1
      IK=I+I-1
      T1=Y(IK  )*Z0(1)-Y(IK+1)*Z0(2)
      T2=Y(IK  )*Z0(2)+Y(IK+1)*Z0(1)
      Y(IK  )=T1
      Y(IK+1)=T2
  650 CONTINUE
      N2=N1
 1000 CONTINUE
      RETURN
      END

c***********************************************************************
      SUBROUTINE BETAFN(D,E,LTR,JX,flfact,i_sym)
C
c     Subroutine to form scattering amplitude from radial matrix elements
c
c     D    Storage for radial matrix elements
c     E    Storage for scattering amplitudes
c     LTR  L  -transfer
c     JX   2*J-transfer
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      logical i_sym(2)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC
      DIMENSION D(918),E(5680)
      logical iflag1,iflag2
C
      ISB=IS(3)
      JR=NS(1)
      JS=NS(2)
c
c     Set up spin for photon
c
      if(fm(1).eq.0.0.and.is(1).eq.2) then
        isi=0
        iflag1=.true.
        z1= z(2)-fm(2)*(z(2)+za(2))/(fm(2)+fma(2))
        z2=(z(2)*fma(2)**2+za(2)*fm(2)**2)/(fm(2)+fma(2))**2
      else
        isi=is(1)
        iflag1=.false.
      endif
      if(fm(2).eq.0.0.and.is(2).eq.2) then
        isf=0
        iflag2=.true.
        z1= z(1)-fm(1)*(z(1)+za(1))/(fm(1)+fma(1))
        z2=(z(1)*fma(1)**2+za(1)*fm(1)**2)/(fm(1)+fma(1))**2
      else
        isf=is(2)
        iflag2=.false.
      endif
c
      MPLUS=JX/2+1
      INCR=LPL2*(LTR+1)
      IFACT=2*MPLUS*JR*JS
      J2K=(1.0+PHASEF(NS(2)))/2.0
      M2K=JX-MPLUS-MPLUS+2
      LX=LTR+LTR
      TEMP2=SQRT(FLOAT((JX+1)*(IS(1)+1)))*flfact
      IF(FN.EQ.0.0) THEN
c       clear amplitude storage unless for coherent sum
        IND=LPLUS*IFACT
        DO 10 M=1,IND
        E(M)=0.0
   10   CONTINUE
      ENDIF
      IS1=-IS(1)
      DO 95 I=1,JR
      IS2=-IS(2)
      DO 90 J=1,JS
      IF(NLTR.NE.1) GO TO 14
      IF(JR*JS.EQ.1) GO TO 15
      IF(IBF(5)+IBF(6).EQ.0) GO TO 15
   14 CONTINUE
C
C     READ RADIAL MATRIX ELEMENTS FROM TAPE 2
C
      READ (2)(D(INDEX),INDEX=1,INCR)
   15 continue
c     final L loop
      DO 80 LL=1,LPLUS
      lf=LL-1
      LLX=lf+lf
      JLX=LLX+IS2
      IF(JLX.LT.0) GO TO 75
      if(i_sym(2)) then
        if(phasef(lf).gt.0.0) then
          temp4=2.0*temp2
        else
          temp4=0.0
        endif
      else
        temp4=temp2
      endif
      TEMP4=temp4*SQRT(FLOAT(JLX+1))*float(llx+1)
      if(iflag2) then
         temp2=temp2*sqrt(float(ll)/(float(lf)+1.e-6))
         if(lf.eq.1) temp4=temp4*z1
         if(lf.eq.2) temp4=temp4*z2
      endif
      LSTOR=lf*IFACT
      LP1=IABS(LL-LTR-1)+1
      LP2=MIN0(LL+LTR,LPLUS)
      LK=0
c     initial L loop
      DO 60 LP=LP1,LP2,2
      li=lp-1
      LPX=LP+LP-2
      JPX=LPX+IS1
      IF(JPX.GE.0) then
      if(i_sym(1)) then
        if(phasef(li).gt.0.0) then
          temp3=2.0*temp4
        else
          temp3=0.0
        endif
      else
        temp3=temp4
      endif
      TEMP=temp3*SQRT(FLOAT(LPX+1))*PHASEF((LP-LL-LTR)/2)
     1    *VCC(LLX,LX,LPX,0,0)*WINEJ(LLX,isf,JLX,LX,ISB,JX,LPX,isi,JPX)
      INDEX=LK+LLX+1
      KT=0
      MS =-IS(2)
      DO 57 MS2=1,JS
      MSP=-IS(1)
      DO 55 MS1=1,JR
      VCP=VCC(LPX,IS(1),JPX,0,MSP)
      DO 50 M=1,MPLUS
      MK=M+M-1
      MX=MK-1+M2K
      ML2=MSP-MX-MS
      ML=IABS(ML2/2)
      IF(ML.LE.lf) then
        IND=LSTOR+KT+MK
        FACT=VCP*VCC(JLX,JX,JPX,MSP-MX,MX)*VCC(LLX,IS(2),JLX,ML2,MS)
     1      *SQRT(YXFCT(lf+ML,lf-ML))*TEMP
        E(IND  )=E(IND  )+D(INDEX  )*FACT
        E(IND+1)=E(IND+1)+D(INDEX+1)*FACT
      endif
   50 CONTINUE
      KT=KT+MPLUS+MPLUS
      MSP=MSP+2
   55 CONTINUE
      MS =MS +2
   57 CONTINUE
      endif
      LK=LK+LPL2
   60 CONTINUE
   75 CONTINUE
   80 CONTINUE
      IS2=IS2+2
   90 CONTINUE
      IS1=IS1+2
   95 CONTINUE
      RETURN
      END

c***********************************************************************
      SUBROUTINE INSIG(D,PLM,JTR,FACTR)
c
c     Calculate inelastic cross sections and spin observables.
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      COMMON ALPHA(15),IDAT(6),ICON(20),ANGLE(5),HBARC,AMU,AMASS,CHSQ,AA
     1,DRF,Z(3),ZA(3),FM(3),FMA(3),RC(3),AC(3),PNLOC(3),FS(3),ECM(3)
     2,FK(3),FK2(3),ETA(3),DR(3),FMU(3),FN,FL,FJ2,FSS,VCE,FNRNG,RSIG(2)
     3,K,KZ,LPLUS,LPL2,IS(3),NS(3),NLTR,LTRT(8),JTRT(8),ISTRT(8),IBF(8)
     4,KC,IBUFF,IWORD,ILINE,JLINE
      Common/array1/DTEMP(8000)
c
      DIMENSION D(5680),SUM(5,5,2),PLM(1230),POLZ(10),POL(10)
      DIMENSION SIGPLT(181),ASYPLT(181)
      EQUIVALENCE (DTEMP,SUM),(DTEMP(51),SIGPLT)
      EQUIVALENCE (ANGLE(1),THETAN),(ANGLE(2),THETA1),(ANGLE(3),DTHETA)
      DATA NX/181/
      data sqrt2,pi/1.41421356,3.14159265/
C
      NTHETA=THETAN
      asymax=0.0
      JR=NS(1)
      JS=NS(2)
      if(is(1).eq.2.and.fm(1).eq.0.0) then
        facta=factr/2.0
      else
        FACTA=FACTR/(2.0*FS(1)+1.0)
      endif
      FACTF=(FK(1)/FK(2))**2
      M2K=(1.0-PHASEF(IS(3)))/2.0
      NPLUS=(JTR+IS(1)+IS(2))/2+1
      MPLUS=JTR/2+1
      IFACT=2*MPLUS*JR*JS
      WRITE(6,9000)
      IF(NTHETA.EQ.0) GO TO 120
      THETA=THETA1
      DO 110 NTH=1,NTHETA
      CALL LGNDR(PLM,NPLUS,LPLUS,THETA)
      SIG=0.0
      DO 10 I=1,10
      POLZ(I)=0.0
   10 CONTINUE
      DO 100 M=1,MPLUS
      MK=M+M-1
      KT=MK
      IS2=-IS(2)
c     Loop for final   spin substates
      DO 70 J=1,JS
      IS1=-IS(1)
c     Loop for initial spin substates
      DO 60 I=1,JR
      ML=-(MK-1+M2K+IS2-IS1)/2
      PHAS=1.0
      IF(ML.GT.0) PHAS=PHASEF(ML)
      ML=IABS(ML)*LPLUS
      SUM1=0.0
      SUM2=0.0
      IND=KT
      DO 40 LL=1,LPLUS
      ML=ML+1
      SUM1=SUM1+D(IND  )*PLM(ML)
      SUM2=SUM2+D(IND+1)*PLM(ML)
      IND=IND+IFACT
   40 CONTINUE
      SUM(J,I,1)=SUM1*PHAS
      SUM(J,I,2)=SUM2*PHAS
      FACT=SUM(J,I,1)**2+SUM(J,I,2)**2
      IF(M+M2K.NE.1) FACT=FACT*2.0
      SIG=SIG+FACT
      KT=KT+MPLUS+MPLUS
      IS1=IS1+2
   60 CONTINUE
      IS2=IS2+2
   70 CONTINUE
      if(fm(1).eq.0.0.and.is(1).eq.2) then
c
c     Force gamma ray to be transverse in initial channel
c
      do 80 j = 1,3
      sum(j,2,1)=0.0
      sum(j,2,2)=0.0
   80 continue
      endif
c
      if(fm(2).eq.0.0.and.is(2).eq.2) then
c
c     Force gamma ray to be transverse in final  channel
c
c     Rotate spin components to outgoing direction and set A_0 = 0
c
      cs=cos(theta*(pi/180.))
      ss=sin(theta*(pi/180.))/sqrt2
      do 90 i=1,3
      t1=sum(1,i,1)*((1.+cs)/2.)-sum(2,i,1)*ss+sum(3,i,1)*((1-cs)/2.)
      t3=sum(1,i,1)*((1.-cs)/2.)+sum(2,i,1)*ss+sum(3,i,1)*((1+cs)/2.)
      sum(1,i,1)=t1
      sum(2,i,1)=0.
      sum(3,i,1)=t3
      t1=sum(1,i,2)*((1.+cs)/2.)-sum(2,i,2)*ss+sum(3,i,2)*((1-cs)/2.)
      t3=sum(1,i,2)*((1.-cs)/2.)+sum(2,i,2)*ss+sum(3,i,2)*((1+cs)/2.)
      sum(1,i,2)=t1
      sum(2,i,2)=0.
      sum(3,i,2)=t3
   90 continue
      endif
c
      CALL POLFCT(POL,SUM(1,1,1),SUM(1,1,2),JS,JR)
      DO 99 I=1,10
      IF(M+M2K.NE.1) POL(I)=POL(I)*2.0
      POLZ(I)=POLZ(I)+POL(I)
   99 CONTINUE
  100 CONTINUE
      IF(SIG.ne.0.0) then
        DO 105 I=1,10
        POLZ(I)=POLZ(I)/SIG
  105   CONTINUE
      endif
      SIG1=SIG *FACTA
      SIG2=SIG1*FACTF
      WRITE(6,9001)THETA,SIG1,(POLZ(I),I=1,4),SIG2,(POLZ(I),I=5,9),THETA
      IF(NTH.LE.NX) THEN
        SIGPLT(NTH)=SIG1
        ASYPLT(NTH)=POLZ(5)
        asymax=max(asymax,abs(asyplt(nth)))
      ENDIF
      THETA=THETA+DTHETA
  110 CONTINUE
  120 CONTINUE
C
C     CALCULATE TOTAL INELASTIC SIGMA
C
      SIG =0.0
      KT=0
      IS2=-IS(2)
      DO 170 J=1,JS
      IS1=-IS(1)
      DO 160 I=1,JR
      DO 150 M=1,MPLUS
      MK=M+M-1
      IND=KT+MK
      MK=MK-1+M2K
      ML=IABS(MK+IS2-IS1)/2
      DO 140 L=1,LPLUS
      LL=L-1
      IF(ML.le.LL) then
        FACT=D(IND  )**2+D(IND+1)**2
        FACT=FACT*YXFCT(LL-ML,LL+ML)/FLOAT(LL+LL+1)
        IF(MK.NE.0) FACT=FACT*2.0
        SIG=SIG+FACT
      endif
      IND=IND+IFACT
  140 CONTINUE
  150 CONTINUE
      IS1=IS1+2
      KT=KT+MPLUS+MPLUS
  160 CONTINUE
      IS2=IS2+2
  170 CONTINUE
      SIG1 =SIG *4.0*pi*FACTA
      SIG2 =SIG1 *FACTF
      WRITE(6,9002)SIG1,SIG2
      IF(NTHETA.ne.0) then
c
      NTH=MIN(NTHETA,NX)
      IF(ICON(13).ne.0) then
        WRITE(7,9011)ALPHA
        WRITE(7,9010)(SIGPLT(N),N=1,NTH)
      endif
      IF(ICON(9).NE.0) then
        WRITE(6,9011)ALPHA,(IDAT(I),I=1,3)
        CALL DWPLOT(NTH,ICON(9),SIGPLT,DTHETA,THETA1)
        NTEMP=-10
          IF(IS(1).NE.0.and.asymax.gt.0.05) THEN
            WRITE(6,9013)ALPHA,(IDAT(I),I=1,3)
            CALL DWPLOT(NTH,NTEMP  ,ASYPLT,DTHETA,THETA1)
          ENDIF
        endif
      endif
  200 RETURN
c
 9000 FORMAT('0 Theta',' Inelsig,fm**2','      Polz','       T22'
     1,   '       T21','       T20',' Inelsig,fm**2','      POLZ'
     2,   '       T22','       T21','       T20',    '       AYY'
     3,'  Theta')
 9001 FORMAT(0PF7.2,1PE14.4,0P4F10.4,1PE14.4,0P5F10.4,0PF7.2)
 9002 FORMAT('0Tot-sig',1PE14.4,40X,1PE14.4)
 9010 FORMAT(1PE10.4,5(4X,1PE10.4))
 9011 FORMAT('1Inelastic  ',15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
 9013 FORMAT('1Asymmetry  ',15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END
c***********************************************************************
      FUNCTION VCC(JX1,JX2,JX3,MX1,MX2)
c
c     Clebsch-Gordan cofficient routine
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FACT(32)
      EQUIVALENCE (FACT(1),F(2))
      EXTERNAL FACTOR
      COMMON/FACTRL/F(33)
C
      VCC=0.0
      J1=JX1
      J2=JX2
      J3=JX3
      M1=MX1
      M2=MX2
      IF(J1.LT.J2) GO TO 20
      IF(J3.LT.J2) GO TO 30
      ICNTR=0
      GO TO 40
   20 IF(J3.LT.J1) GO TO 30
      ICNTR=-1
      IT=J1
      J1=J2
      J2=IT
      IT=M1
      M1=M2
      M2=IT
      GO TO 40
   30 ICNTR=1
      IT=J2
      J2=J3
      J3=IT
      M2=-M1-M2
   40 CONTINUE
      JZ1=(J1+J2-J3)/2
      IF(JZ1.LT.0) GO TO 150
      JZ2=(J1+J3-J2)/2
      IF(JZ2.LT.0) GO TO 150
      JZ3=(J2+J3-J1)/2
      IF(JZ3.LT.0) GO TO 150
      IF(J1-IABS(M1).LT.0) GO TO 150
      IF(J2-IABS(M2).LT.0) GO TO 150
      IF(J3-IABS(M1+M2).LT.0) GO TO 150
      JT1=(J1-J3+M2)/2
      JT2=(J2-J3-M1)/2
      NUMIN=MAX0 (JT1,JT2,0)
      JT3=(J1-M1)/2
      JT4=(J2+M2)/2
      NUMAX=MIN0 (JT3,JT4,JZ1)
      JT5=(J2-M2)/2
      IF(NUMAX.ge.NUMIN) then
      J4=J1/2
      J5=J3/2
      PHAS=PHASEF(NUMIN)
      DO 100 NU=NUMIN,NUMAX
      VCC=VCC+PHAS      *(YXFCT(JT3-NU,J4)*YXFCT(NU-JT2,J5))
     1/(FACT(JT4-NU)*FACT(NU-JT1)*FACT(JZ1-NU)*FACT(NU))
      PHAS=-PHAS
  100 CONTINUE
      FCTOR=YXFCT(J4,(J1+M1)/2)*YXFCT(J4,JT3)*YXFCT((J1+J2+J3)/2+1,JZ2)
     1*YXFCT(J5,(J3+M1+M2)/2)*YXFCT(J5,(J3-M1-M2)/2)*FACT(JZ1)*FACT(JZ3)
     2*FACT(JT4)*FACT(JT5)*FLOAT(J3+1)
      VCC=SQRT(FCTOR)*VCC
      IF(ICNTR)120,150,110
  110 VCC=VCC*SQRT(FLOAT(J2+1)/FLOAT(J3+1))*PHASEF(JT3)
      GO TO 150
  120 VCC=VCC*PHASEF(JZ1)
      endif
  150 RETURN
      END

      BLOCK DATA FACTOR
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/FACTRL/F(33)
C
      DATA F   / 1.0000000000D+00, 1.0000000000D+00, 2.0000000000D+00
     1         , 6.0000000000D+00, 2.4000000000D+01, 1.2000000000D+02
     2         , 7.2000000000D+02, 5.0400000000D+03, 4.0320000000D+04
     3         , 3.6288000000D+05, 3.6288000000D+06, 3.9916800000D+07
     4         , 4.7900160000D+08, 6.2270208000D+09, 8.7178291200D+10
     5         , 1.3076743680D+12, 2.0922789888D+13, 3.5568742810D+14
     6         , 6.4023737057D+15, 1.2164510041D+17, 2.4329020082D+18
     7         , 5.1090942172D+19, 1.1240007278D+21, 2.5852016739D+22
     8         , 6.2044840173D+23, 1.5511210043D+25, 4.0329146113D+26
     9         , 1.0888869450D+28, 3.0488834461D+29, 8.8417619937D+30
     $         , 2.6525285981D+32, 8.2228386542D+33, 2.6313083693D+35/
C    $         , 8.6833176188D+36, 2.9523279904D+38, 1.0333147966D+40
C    $         , 3.7199332679D+41, 1.3763753091D+43, 5.2302261747D+44
C    $         , 2.0397882081D+46, 8.1591528325D+47, 3.3452526613D+49
C    $         , 1.4050061178D+51, 6.0415263063D+52, 2.6582715748D+54
C    $         , 1.1962222087D+56, 5.5026221598D+57, 2.5862324151D+59
C    $         , 1.2413915593D+61, 6.0828186403D+62, 3.0414093202D+64
C    $         , 1.5511187533D+66/
      END
C
c***********************************************************************
      FUNCTION PHASEF(N)
c
c     Calculates (-1)**N
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)

      PHASEF=FLOAT(1-2*IABS(N-2*(N/2)))
      RETURN
      END

c***********************************************************************
      FUNCTION YXFCT(M,N)
c
c     calculates N_factorial/M_factorial
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      logical ictrl
c
      YXFCT=1.0
      NUMAX=M-N
      if    (NUMAX.gt.0) then
        ictrl=.true.
        FCTOR=N
      elseif(NUMAX.lt.0) then
        ictrl=.false.
        NUMAX=-NUMAX
        FCTOR=M
      else
        return
      endif
      DO 50 NU=1,NUMAX
      FCTOR=FCTOR+1.0
      YXFCT=YXFCT*FCTOR
   50 CONTINUE
      IF(ictrl) YXFCT=1.0/YXFCT
      RETURN
      END

c***********************************************************************
      FUNCTION RACAH(J1,J2,J3,J4,J5,J6)
c
c     Calculates Racah coefficients
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION FACT(32)
      EXTERNAL FACTOR
      COMMON/FACTRL/F(33)
      EQUIVALENCE (FACT(1),F(2))
c
      RACAH=0.0
      Z1=DELR(J1,J2,J5)
      IF(Z1.EQ.0.0) GO TO 90
      Z1=DELR(J3,J4,J5)*Z1
      IF(Z1.EQ.0.0) GO TO 90
      Z2=DELR(J1,J3,J6)
      IF(Z2.EQ.0.0) GO TO 90
      Z2=DELR(J2,J4,J6)*Z2
      IF(Z2.EQ.0.0) GO TO 90
      Z1=SQRT(Z1/Z2)*Z2
      JT1=(J1+J2+J5)/2
      JT2=(J3+J4+J5)/2
      JT3=(J1+J3+J6)/2
      JT4=(J2+J4+J6)/2
      JZ1=(J1+J2+J3+J4)/2
      JZ2=(J1+J4+J5+J6)/2
      JZ3=(J2+J3+J5+J6)/2
c
      NUMIN=MAX0 (JT1,JT2,JT3,JT4)
      NUMAX=MIN0 (JZ1,JZ2,JZ3)
      IF(NUMAX.ge.NUMIN) then
      PHASE=PHASEF(NUMIN+JZ1)*Z1
      DO 80 NU=NUMIN,NUMAX
      JY1=NU-JT1
      JY2=NU-JT2
      JY3=NU-JT3
      JY4=JZ1-NU
      JY5=JZ2-NU
      JY6=JZ3-NU
      FCTOR=FACT(JY1)*FACT(JY2)*FACT(JY3)*YXFCT(NU+1,NU-JT4)
     1*FACT(JY4)*FACT(JY5)*FACT(JY6)
      RACAH=RACAH+PHASE/FCTOR
      PHASE=-PHASE
   80 CONTINUE
      endif
   90 RETURN
      END

c***********************************************************************
      FUNCTION DELR(J1,J2,J3)
c
c     Triangular function for RACAH
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION FACT(32)
      EXTERNAL FACTOR
      COMMON/FACTRL/F(33)
      EQUIVALENCE (FACT(1),F(2))
c
      JZ1=(J1+J2-J3)/2
      IF(JZ1.LT.0) GO TO 130
      JZ2=(J1-J2+J3)/2
      IF(JZ2.LT.0) GO TO 130
      JZ3=(J2+J3-J1)/2
      IF(JZ3.LT.0) GO TO 130
      JZ4=(J1+J2+J3)/2+1
      IF(JZ3.LT.JZ2) GO TO 80
      IF(JZ3.LT.JZ1) GO TO 70
      DELR=YXFCT(JZ4,JZ3)*FACT(JZ1)*FACT(JZ2)
      RETURN
   70 DELR=YXFCT(JZ4,JZ1)*FACT(JZ2)*FACT(JZ3)
      RETURN
   80 IF(JZ2.LT.JZ1) GO TO 70
      DELR=YXFCT(JZ4,JZ2)*FACT(JZ1)*FACT(JZ3)
      RETURN
  130 DELR=0.0
      RETURN
      END

c***********************************************************************
      FUNCTION WINEJ(J1,J2,J3,J4,J5,J6,J7,J8,J9)
c
c     Nine-J symbol function
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      WINEJ=0.0
      MUMIN=MAX0(IABS(J1-J9),IABS(J2-J6),IABS(J4-J8))
      MUMAX=MIN0(J1+J9,J2+J6,J4+J8)
      IF(MUMAX.ge.MUMIN) then
        DO 20 MU=MUMIN,MUMAX,2
        PROD=RACAH(J1,J4,J9,J8,J7,MU)*RACAH(J2,J5,MU,J4,J8,J6)*
     1     RACAH(J9,MU,J3,J2,J1,J6)*FLOAT(MU+1)
        WINEJ=WINEJ+PROD
   20   CONTINUE
        WINEJ=WINEJ*PHASEF((J1+J3+J5+J8)/2+J2+J4+J9)
      endif
      RETURN
      END

c***********************************************************************
      SUBROUTINE BIND(U,F,DR,RM,FNODE,FL,KMAX,FK,ETA,V,E,FK2,ISW,IERR,D)
c
c     Calculates a single particle wave function
c***********************************************************************
C
C     U       = POTENTIAL STORAGE
C        ODD  =             - SCALED PART OF V *(2M/HBAR**2)
C        EVEN = (ENERGY - NON SCALED PART OF V)*(2M/HBAR**2)
C     F       = OUTPUTTED WAVE FUNCTION
C     DR      = RADIAL INCREMENT
C     RM      = MATCHING RADIUS, INPUT IS DEFAULT, OUTPUT IS VALUE USED
C     FNODE   = NO. NODES IN WAVE FCT., EXCLUDING ORIGIN AND INFINITY
C     FL      = ORBITAL ANGULAR MOMENTUM
C     KMAX    = MAX NO. OF POINTS + 2048*(MIN POINT)
C     FK      = ]WAVE NUMBER] E.G. SIN(FK*R) OR EXP(-FK*R)
C     ETA     = COULOMB PARAMETER
C     V       = SCALE FACTOR FOR ODD POTENTIAL VALUES
C     E       = ENERGY
C     FK2     = E*(2*M/HBAR**2)
C     ISW     = SEARCH SWITCH
C        0    = SEARCH ON V
C        1    = SEARCH ON E
C        2    = NO SEARCH, FOR E>0 ONLY (PHASE SHIFT IS CALCULATED)
C        3    = Search on V, box b.c.
C        4    = Search on E, box b.c.
C     IERR    = ERROR FLAG, IS SET IF ERROR OCCURS
C     D       = LAST POINT DATA FOR UNBOUND STRIPPING
C     D(1),D(2) ARE REAL AND IMAG PARTS OF EXP(2I*DELTA)
C     D(3),D(4) ARE GL(RMAX) AND FL(RMAX)
c
      IMPLICIT REAL*8(A-H,O-Z)
      logical dplus,dminus,turnpt
      CHARACTER*2 A
      DIMENSION U(800),F(400),T(8),D(4)
      DATA ICMAX/32/
      DATA DLAMY,DETX/.0001,.001/
      DATA XFACT,ELOG10/0.707,2.302585/
C
C     CHANGES TO ALLOW INTEGRATION OF HARD CORE POTENTIALS
C     THEY ARE KLUDGED TO ALLOW BACKWARD COMPATIBILITY
C
      K1 = KMAX/2048
      KT = KMAX-K1*2048
      K1=K1+1
C
      CFACT=1.0
      dv=0.225*v
      v0=0.0
      v1=0.0
      d0 =-1.0e20
      d1 = 1.0e20
      dplus =.false.
      dminus=.false.
      IXNT=0
      IPRNT=0
      IF(IERR.LT.0) IPRNT=1
      IERR=0
      ICNT=ICMAX
      IF(ISW.EQ.2) ICNT=0
      FLP=FL*(FL+1.0)
      DR2=DR**2/12.0
      LL=FL+1.0
      ITER=0
C
C     CALCULATE OUTER BOUNDARY CONDITION
C
   10 CONTINUE
      RMAX=DR*FLOAT(KT+1)
      IF(FK2.GT.0.0) THEN
C     POSITIVE ENERGY BOUNDARY CONDITION
      IF(ISW.EQ.0.OR.ISW.EQ.1.OR.ISW.EQ.2) THEN
        M=MAX0(LL,3)
        T(2)=FK*RMAX-FK*DR
        T(1)=T(2)-FK*DR
        CALL COU(T(1),T(2),ETA,M ,DR,F( 1),F(21),F(41),F(61),F(81))
        D(4)=F(LL+20)
        D(3)=F(LL+60)
        T(8)=F(LL+60)
        T(7)=F(LL+40)
        T(6)=F(LL+20)
        T(5)=F(LL   )
        GP=D(3)
        M1=1
      ELSEIF(ISW.EQ.3.OR.ISW.EQ.4) THEN
        T(8)=0.0
        T(7)=1.0E-15
        GP  =0.0
        M1=3
        F(KT  )=0.0
        F(KT-1)=T(7)
        RMAX=FLOAT(KT-1)*DR
      ELSE
        GO TO 900
      ENDIF
      IF(ISW.EQ.0.OR.ISW.EQ.1) THEN
C
C     CALCULATE QUANTITIES FOR WIDTH
C
      SCALE=ABS(F(LL+60))
      DX=.01
      DK=T(2)*DX
      T1=1.0-DX/2.0
      T2=1.0+DX/2.0
      T(3)=T(2)*T1*T2
      T(4)=T(2)*T2*T2
      ETAP=ETA/T2
      CALL COU(T(3),T(4),ETAP,M,DR,F( 1),F(21),F(41),F(61),F(81))
      SURF =        2.0*FK*(F(LL+60)-F(LL+40))/(F(LL+60)+F(LL+40))/DK
      T(3)=T(2)*T1*T1
      T(4)=T(2)*T2*T1
      ETAP=ETA/T1
      CALL COU(T(3),T(4),ETAP,M,DR,F( 1),F(21),F(41),F(61),F(81))
      SURF = SURF - 2.0*FK*(F(LL+60)-F(LL+40))/(F(LL+60)+F(LL+40))/DK
      SURF=SURF/(DX*2.0*FK2)
      T(1)=FK*RMAX
      T(2)=T(1)+FK*DR
      CALL COU(T(1),T(2),ETA,M ,DR,F( 1),F(21),F(41),F(61),F(81))
C     SCALE FUNCTIONS TO PREVENT OVERFLOW
      T(8)=F(LL+60)/SCALE
      T(7)=F(LL+40)/SCALE
      T(6)=F(LL+20)/SCALE
      T(5)=F(LL   )/SCALE
      D(3)=D(3)    /SCALE
      D(4)=D(4)    /SCALE
      GP=D(3)
      ENDIF
      ELSE
C     NEGATIVE ENERGY BOUNDARY CONDITION
      IF(ISW.EQ.0.OR.ISW.EQ.1) THEN
        T(8)=1.0E-15
        T(7)=T(8)*EXP(FK*DR+ETA* LOG(1.0+DR/RMAX))
        GP  =T(7)**2/T(8)
        M1=1
C     BOX BOUNDARY CONDITIONS
      ELSEIF(ISW.EQ.3.OR.ISW.EQ.4) THEN
        T(8)=0.0
        T(7)=1.0E-15
        GP  =0.0
        M1=3
        F(KT  )=0.0
        F(KT-1)=T(7)
        RMAX=FLOAT(KT-1)*DR
      ELSE
        GO TO 900
      ENDIF
      ENDIF
   40 CONTINUE
      turnpt=.true.
      KTEMP=KT
      KM=KTemp
C
C     INTEGRATE FROM INFINITY TO TURNING POINT
C
   50 continue
      R=RMAX
      FNCT=0.0
      IF(ISW.NE.2) THEN
      G5=U(2*KT )-FLP/R**2
      IF(G5.LT.0.0) THEN
        IFLAG=1
      ELSE
        IFLAG=0
      ENDIF
      Q6=1.0+DR2*(U(2*KT )-FLP/(R+DR)**2)
      Q5=1.0+DR2*G5
      F6=T(8)
      F5=T(7)
      W2=0.0
      FNORM2=-0.5*GP**2
C
      DO 100 M=M1,KTEMP
      MM=KT-M+1
      MK=MM+MM-1
      R=R-DR
      G4=U(MK+1)+V*U(MK  )-FLP/R**2
      Q4=1.0+DR2*G4
      F4=((12.0-10.0*Q5)*F5-Q6*F6)/Q4
      Q6=Q5
      Q5=Q4
      G6=G5
      G5=G4
      if(f6*f5.lt.0.0) fnct=fnct+1.0
      F6=F5
      F5=F4
      F(MM)=F4
      IF(G6*G5.LT.0.0) IFLAG=IFLAG+1
      if((.not.turnpt.and.f5.lt.f6).or.IFLAG.GE.2) GO TO 110
      FNORM2=FNORM2+F5**2
      W2=W2+U(MK+2)*F5**2
  100 CONTINUE
C
C     If no classical turning point is found, then find first maximum,
C     then if none, use external matching radius.
C
      if(turnpt) then
        turnpt=.false.
        go to 50
      else
        KTEMP=KT-INT(RM/DR)
        MM=KT-KTEMP+1
      endif
  105 continue
      F5=F(MM  )
      F6=F(MM+1)
  110 CONTINUE
C
C     INTEGRATE FROM ORIGIN TO TURNING POINT
C
      if(isw.lt.2) fnct = 0.0
      KM=MM+1
  120 CONTINUE
      ENDIF
      KS=FL/3.3+2.0
      W1=0.0
      FNORM1=0.0
      F2=0.0
      Q2=0.0
      R =0.0
C
      DO 200 M=1,KM
      MK=M+M-1
      R=R+DR
      Q3=1.0+DR2*(U(MK+1)+V*U(MK  )-FLP/R**2)
      IF(M.GT.KS) THEN
        F3=((12.0-10.0*Q2)*F2-Q1*F1)/Q3
      ELSE
        F3=R**LL
      ENDIF
      Q1=Q2
      Q2=Q3
      F1=F2
      F2=F3
      IF(ICNT.EQ.0) F(M)=F3
      FNORM1=FNORM1+F1**2
      W1=W1+U(MK)*F1**2
C     NODE COUNT
      IF(F1*F2.LT.0.0) FNCT=FNCT+1.0
  200 CONTINUE
      F12=(F1+F2)*0.5
      IF(F1*F2.lt.0.0.and.F12*F2.lt.0.0) FNCT=FNCT-1.0
C
      FN=FNODE-FNCT
      IF(ISW.EQ.2) GO TO 500
      F56=(F5+F6)*0.5
      DET=(F1*F6-F5*F2)/(F12*F56*DR)
      FNORM1=FNORM1/F2**2
      FNORM2=FNORM2/F6**2
      FNORM=FNORM1+FNORM2
      ITER=ITER+1
      IF(ICNT.EQ.0) GO TO 500
      ICNT=ICNT-1
      IF(ICNT.EQ.0) THEN
        WRITE(6,9001)ICMAX
        IERR=1
      ENDIF
      IF(ICNT.EQ.18) CFACT=CFACT*XFACT
      IF(ICNT.EQ.12) CFACT=CFACT*XFACT
      RM=(FLOAT(KM)-0.5)*DR
      IF((ISW.EQ.3.OR.ISW.EQ.4).AND.E.GT.0.0) THEN
        RSCALE=RMAX
      ELSE
        RSCALE=RM*1.5
      ENDIF
      ESCALE= ABS(((fnode+1.0)*3.0+fl)/(FK2*RSCALE**2))
      DLAMX =DLAMY*ESCALE
        vold=v
        eold=e
      IF(ISW.EQ.0.OR.ISW.EQ.3) THEN
C
C     CHOOSE NEXT GUESS FOR WELL DEPTH
C
        IF(FN.NE.0.0) THEN
c       Node count incorrect
          dlam= 1.625*cfact*fn*escale
          if(dlam.lt.-abs(dv/v)) then
            dlam = -abs(dv/v)
          else
            dlam = min(dlam,1.0d0)
          endif
          a='Vn'
        ELSE
          if(det.gt.0.0) then
            dplus =.true.
            if(det.lt.d1) then
              v1=v
              d1=det
            endif
          else
            dminus=.true.
            if(det.gt.d0) then
              v0=v
              d0=det
            endif
          endif
          IF(dplus.and.dminus) THEN
c         Interpolation search
            DLAM= ((d1*v0-d0*v1)/((d1-d0)*v) - 1.0)
            a='Vi'
          ELSE
c         Variational search
            DLAM= -DET/(V*DR*(W1/F12**2+W2/F56**2))
            IF(ABS(DLAM).GT.0.3) DLAM=SIGN(0.3d0,DLAM)
            a='Vv'
          ENDIF
        IXNT=IXNT+1
        ENDIF
        flam=1.0+dlam
        V=V*FLAM
        DLAMX=DLAMY
C
C     CHOOSE NEXT GUESS FOR BINDING ENERGY
C
      ELSEIF(ISW.EQ.1.OR.ISW.EQ.4) THEN
        IF(FN.NE.0.0) THEN
c       Node count incorrect
          DLAM= 1.625*cfact*fn*sign(escale,e)
          a='En'
        ELSE
          if(det.gt.0.0) then
            dplus =.true.
            if(det.lt.d1) then
              v1=e
              d1=det
            endif
          else
            dminus=.true.
            if(det.gt.d0) then
              v0=e
              d0=det
            endif
          endif
          IF(dplus.and.dminus) THEN
c         Interpolation search
            DLAM= ((d1*v0-d0*v1)/((d1-d0)*e) - 1.0)
            a='Ei'
          ELSE
c         Variational search
            DLAM= -cfact*DET/(DR*FK2*FNORM)
            a='Ev'
          IF(ABS(DLAM).GT.ESCALE) DLAM = SIGN(ESCALE,DLAM)
          ENDIF
        ENDIF
        FLAM=1.0+DLAM
        TEMP=SQRT(ABS(FLAM))
        FK=FK*TEMP
        ETA=ETA/TEMP
        TEMP=FK2*FLAM-FK2
        FK2=FK2+TEMP
        E=E*FLAM
        DO 485 M=1,KT
        MK=M+M-1
        U(MK+1)=U(MK+1)+TEMP
  485   CONTINUE
      ENDIF
c     diagnostic printout
c      if(icnt+1.eq.icmax) write(6,8888)
c 8888 format('     Iter',' Nodes',5x,'  Det',9x,'Fnorm',9x,'Old E'
c     1,  9x,'New E',9x,'Old V',9x,'New V',9x,' Flam',9x,'   RM')
c      WRITE(6,8889)A,ITER,FNCT,DET,FNORM,eold,E
c     1  ,vold,V,1.0+DLAM,rm
c 8889 FORMAT(' ',A2,I6,F6.0,8(2X,1PE12.5))
c
      IF(ABS(DET).LT.DETX .AND. ABS(DLAM).LT.DLAMX) THEN
        ICNT=0
      ENDIF
      GO TO 10
C
  500 CONTINUE
C
C     NORMALIZE WAVE FUNCTION
C
      IF(FK2.LT.0.0.OR.ISW.EQ.3.OR.ISW.EQ.4) THEN
C
C     NEGATIVE ENERGY AND BOX B.C.
C
        FNORM=SQRT(FNORM*DR)
      ELSEIF(FK.GT.0.0.AND.(ISW.EQ.0.OR.ISW.EQ.1)) THEN
C     POSITIVE ENERGY AND PI/2 PHASE SHIFT
C
        D(1)=0.0
        D(2)=1.0
        VOL=(FNORM1+FNORM2)*DR*F6**2
        WIDTH=2.0*E/(FK*(VOL+SURF))
        IF(WIDTH.LT.0.0) WRITE(6,9002)
        WIDTH=ABS(WIDTH)
        TEMP=( LOG(WIDTH)-2.0* LOG(SCALE))/ELOG10
        I1=TEMP-1.0
        A1=EXP(ELOG10*(TEMP-FLOAT(I1)))
        WRITE(6,9502)A1,I1,VOL*SCALE**2,SURF*SCALE**2
        FNORM=FK/SQRT(WIDTH*FK2*FK/(2.0*E))/F6
C
      ELSEIF(FK.GT.0.0.AND.ISW.EQ.2) THEN
C
C     FIND B.C. FOR E > 0  AND NO SEARCH (ISW = 2)
C
        DET=T(5)*T(8)-T(6)*T(7)
        A1=(F1  *T(8)-F2  *T(7))/DET
        B1=(T(5)*F2  -T(6)*F1  )/DET
        DET=1.0/SQRT(A1**2+B1**2)
        A1=A1*DET
        B1=B1*DET
C     A1=COS(DELTA), B1=SIN(DELTA)
        WRITE(6,9501)A1,B1
        D(1)=B1*A1
        D(2)=B1*B1
        FNORM=FK
        F2=1.0/DET
        F6=1.0
        DET=0.0
      ELSE
        GO TO 900
      ENDIF
C
      TEMP=1.0/(F2*FNORM)
      R=0.0
      DO 510 M=1,KM
      R=R+DR
      F(M)=F(M)*TEMP/R
  510 CONTINUE
      IF(KM.LT.KT) THEN
        KM=KM+1
        TEMP=1.0/(F6*FNORM)
        DO 520 M=KM,KT
        R=R+DR
        F(M)=F(M)*TEMP/R
  520 CONTINUE
      ENDIF
      IF(IPRNT.EQ.0.OR.IERR.NE.0)  WRITE(6,9500)V,DET,FNCT,RM,E,ITER
      RETURN
  900 CONTINUE
        IERR=1
        WRITE(6,9003)
        RETURN
C
 9001 FORMAT(28H0BOUND STATE SEARCH FAILS IN,I3,12H ITERATIONS  )
 9002 FORMAT(29H0WARNING, WIDTH IS NEGATIVE  ,50(1H*))
 9003 FORMAT('0FISW parameter is out of range or is illegal')
 9500 FORMAT(21X,6HV    =,F9.4,3X,6HDET  =,F9.4,3X,6HNODES=,F9.4,3X,
     1 6HRM   =,F9.4,3X,6HE    =,F9.4,3X,6HITER.=,I4  )
 9501 FORMAT(21X,6HCOSD =,F9.4,9H   SIND =,F9.4)
 9502 FORMAT(24H0SINGLE PARTICLE WIDTH = ,F7.4,1HE,I3,'  MEV     '
     1,13HVOLUME TERM =  ,1PE12.4,5X,14HSURFACE TERM =  ,1PE12.4)
      END


c***********************************************************************
      SUBROUTINE COU(R,RP,ETA,L,H,F,FP,G,GP,S)
c
c     Coulomb function subroutine, calculates the functions at 2 points,
c     R and RP for convenience in matchig to the output of the
c     integration formulae.
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(51),FP(51),G(51),GP(51),S(51)
C
      LL=L
      IF(LL.LT.50) THEN
        ELP=50.
        J=50
      ELSE
        ELP=LL
        J=LL
      ENDIF
      A=ATAN (ETA/ELP)
      B=SQRT (ETA**2+ELP**2)
      Y=A*(ELP-0.5)+ETA*( LOG(B)-1.)-SIN (A)/(12.*B)
     1  +SIN (3.*A)/(360.*B**3)-SIN (5.*A)/(1260.*B**5)
     2  +SIN (7.*A)/(1680.*B**7)-SIN (9.*A)/(1188.*B**9)
C
      K=J-1
      DO 100 I=1,K
      IF(J.LE.LL) S(J)=Y
      ELP=ELP-1.
      J=J-1
      Y=Y-ATAN (ETA/ELP)
  100 CONTINUE
      S(1)=Y
      DEL1=R-2.0*ETA
      RMAX=.41666667*(ETA**2+4.*ETA+3.)
      IF(RMAX.LT.10.0) RMAX=10.0
      DEL=R-RMAX
      IF(ETA.GE.5.0.AND.ABS (DEL1).LT.ABS (DEL)) THEN
      DEL=DEL1
        I1=2
C     Transition line expansion
      X=2.0*ETA
      T1=ETA**2
      T2=T1**2
      T3=ETA** .666666667
      T4=T3**2
      T5=T4**2
      T6=T3*T5
      T7=T4*T6
      T8=T3*T7
      T9=ETA** .166666667
      Y=1.22340402*T9*(1.+.495957017E-1/T4-.888888889E-2/T1
     1  +.245519918E-2/T6-.910895806E-3/T2+.845362E-3/T8)
      Z=-.707881773/T9*(1.-.172826039/T3+.317460317E-3/T1
     1  -.358121485E-2/T5+.311782468E-3/T2-.907396643E-3/T7)
      ELSE
C
      IF(DEL.GE.0.0.OR.ETA.EQ.0.0) THEN
        X=R
        I1=1
      ELSE
        X=RMAX
        I1=2
      ENDIF
C     Asymptotic expansion
      T1=ETA**2
      T2=2.*X
      T4=ETA/T2
      SS=1.
      TS=0.
      SL=0.
      TL=1.-ETA/X
      SSS=1.
      STS=0.
      SSL=0.
      STL=TL
      EN=0.
      DO 700 K=1,25
      IF(ABS (SS/SSS).GT.1.0E-10) THEN
        T5=EN+1.
        T6=T5+EN
        T7=EN*T5
        T8=T6*T4/T5
        T9=(T1-T7)/(T2*T5)
        T5=T8*SS-T9*TS
        TS=T8*TS+T9*SS
        SS=T5
        T5=T8*SL-T9*TL-SS/X
        TL=T8*TL+T9*SL-TS/X
        SL=T5
        SSS=SSS+SS
        STS=STS+TS
        SSL=SSL+SL
        STL=STL+TL
        EN=EN+1.
      ENDIF
  700 CONTINUE
C
      T3=X-ETA* LOG(T2)+S(1)
      T8=SIN (T3)
      T9=COS (T3)
      Y=SSS*T9-STS*T8
      Z=SSL*T9-STL*T8
      ENDIF
C
      DO 810 I=1,I1
      IF(I.EQ.I1) THEN
        G(1) = Y
        W    = Z
        DEL = RP - R
      ENDIF
C     Runge-Kutta integration
      N=ABS (DEL/H)+1.0
      DX=DEL/FLOAT(N)
      T1=DX/2.
      T2=DX/8.
      T3=2.0*ETA
      DO 805 J=1,N
      T4=DX*(T3/X-1.)*Y
      X=X+T1
      T5=DX*(T3/X-1.)*(Y+T1*Z+T2*T4)
      X=X+T1
      T6=DX*(T3/X-1.)*(Y+DX*Z+T1*T5)
      Y=Y+DX*(Z+(T4+2.*T5)/6.)
      Z=Z+(T4+4.*T5+T6)/6.
  805 CONTINUE
  810 CONTINUE
      GP(1)=Y
C
      T1=ETA**2
      T8=SQRT (1.+T1)
      G(2)=((1./R+ETA)*G(1)-W)/T8
      GP(2)=((1./RP+ETA)*Y-Z)/T8
      T2=1.
      T3=2.
C     Recur irregular function upwards
      DO 910 I=3,LL
      T4=T2+T3
      T5=T2*T3
      T6=T3*SQRT (T2**2+T1)
      T7=T2*SQRT (T3**2+T1)
      G (I)=(T4*(ETA+T5/R )*G (I-1)-T6*G (I-2))/T7
      GP(I)=(T4*(ETA+T5/RP)*GP(I-1)-T6*GP(I-2))/T7
      T2=T2+1.
      T3=T3+1.
  910 CONTINUE
      N=MAX0(L+11, INT(2.0*R+11.0) )
      I1=N
      Y =1.0E-20
      YP=Y
      Z =0.
      ZP=Z
C     Recur regular function downwards
      DO 930 I=1,I1
      T2=N
      T3=T2+1.
      T4=T2+T3
      T5=T2*T3
      T6=T2*SQRT (T3**2+T1)
      T7=T3*SQRT (T2**2+T1)
      X =(T4*(ETA+T5/R )*Y -T6*Z )/T7
      XP=(T4*(ETA+T5/RP)*YP-T6*ZP)/T7
      IF(N.LE.LL) THEN
        F (N)=X
        FP(N)=XP
      ELSE
        IF(ABS (X).GT.1.0) THEN
          Y =Y *1.0E-20
          YP=YP*1.0E-20
          X =X *1.0E-20
          XP=XP*1.0E-20
        ENDIF
      ENDIF
      N=N-1
      Z =Y
      ZP=YP
      Y =X
      YP=XP
  930 CONTINUE
      Z =1./(T8*(F (1)*G (2)-F (2)*G (1)))
      ZP=1./(T8*(FP(1)*GP(2)-FP(2)*GP(1)))
      DO 950 I=1,LL
      F (I)=F (I)*Z
  950 FP(I)=FP(I)*ZP
      RETURN
       END


c***********************************************************************
      SUBROUTINE DSTRIP(IQ,DR,K,F1,F2,FR,QNUM,OPT,KM,SL,C)
c
c     Calculates two nucleon transfer form factor via the
c     Bayman and Kallio method.
c***********************************************************************
c
c     IQ    Input quantum numbers for form factor
c     DR    Step size
c     K     Number of steps
c     F1    First  orbital
c     F2    Second orbital
c     FR    Output form factor
c     QNUM  Quantum numbers for orbitals
c     OPT
c
      implicit real*8(a-h,o-z)
      parameter (maxg=10, maxr=12)
      DIMENSION F1(400),F2(400),FR(800),QNUM(4,2),TVCC(10),IQ(3)
     1         ,D1(10),D2(10),C(32),AG(maxg),WG(maxg),BG(maxg)
     2         ,AR(maxr),WR(maxr)
C
      DATA KR,KX/0,0/
C
c     c1 = R1 scale
c     c2 = R2 scale
c     c3 = r1 scale
c     c4 = r2 scale
c     c5 = r0, integration scale
c     c6 = Pauli flag
c     c7 = order of integration
c     c8 = number of integration points
c
      IPAULI=C(6)
      R0    =C(5)
      L =IQ(1)
      IS=IQ(2)
      JJ=IQ(3)
      N1=QNUM(1,1)
      N2=QNUM(1,2)
      L1=QNUM(2,1)
      L2=QNUM(2,2)
      J1=QNUM(3,1)
      J2=QNUM(3,2)
      IS1=QNUM(4,1)
      IS2=QNUM(4,2)
      c3=c(3)
      c4=c(4)
      dr1=dr*c(1)
      dr2=dr*c(2)
      ITEMP=N1+N1+N2+N2
      NX=(((ITEMP+L1+L2+L+2)/4+2)+1+2)/2+2*IPAULI
      if(c(8).eq.0.0) then
        NR=(((ITEMP+L1+L2-L+3)/4+8)+1+2)/2+2*IPAULI
        NR=MIN0(NR,maxr)
        IF(NR.NE.KR) then
          CUT=0.0
          IERR=0
          ALFA=C(7)
          CALL GAUSSR(NR,KS,ALFA,AR,WR,IERR,CUT)
          NR=KS
          kr=nr
        endif
      else
        nr = c(8)
        do 45 i=1,nr
        ar(i)=c(i+8)
        wr(i)=c(i+nr+8)
        if(wr(i).ne.0.0) kr=i
   45   continue
      nr=kr
      endif
      NX=MIN0(NX,maxg)
      WRITE(6,9000)NR,NX
 9000 FORMAT(15H0 NO. R STEPS =,I3,18H     NO. X STEPS =,I3)
      IF(NX.NE.KX) then
        CALL LEGAUS(2*NX,AG,WG)
        KX=NX
        do 47 i=1,nx
        BG(i) = sqrt(1.0-ag(i)**2)
   47   continue
      endif
      LL=L+L
      LL1=L1+L1
      LL2=L2+L2
      FACT1=YXFCT(L1,LL1)/2.0**L1
      FACT2=YXFCT(L2,LL2)/2.0**L2
      FNORM=OPT*SQRT(FLOAT((LL1+1)*(LL2+1))/FLOAT(LL+1))
      IF(IPAULI.NE.0) GO TO 60
      TEMP=2.0
      DO 50 I=1,4
      IF(QNUM(I,1).NE.QNUM(I,2)) GO TO 55
   50 CONTINUE
      TEMP=4.0
   55 CONTINUE
      FNORM=FNORM/SQRT(TEMP)
   60 CONTINUE
      FNORM=FNORM*SQRT(FLOAT((LL+1)*(IS+1)*(J1+1)*(J2+1)))
     1     *WINEJ(LL1,IS1,J1,LL2,IS2,J2,LL,IS,JJ)
      FM1=1.0
      FM2=1.0
      FL1=L1
      FL2=L2
      LPL=MIN0(L1,L2)+1
      DO 80 M=1,LPL
      M2=M+M-2
      FM=M-1
      TVCC(M)=VCC(LL1,LL2,LL,M2,-M2)*2.0/SQRT(FM1*FM2)
      FM1=FM1*(FL1+FM+1)*(FL1-FM)
      FM2=FM2*(FL2+FM+1)*(FL2-FM)
   80 CONTINUE
      TVCC(1)=TVCC(1)/2.0
C
C     RS=r
C
C     R1=] C1*R+C3*r ]
C     R2=] C2*R+C4*r ]
C
C
C     CENTER OF MASS R LOOP
C
      R=0.0
      S=0.0
      DO 500 M=1,K
      R=R+DR1
      S=S+DR2
      RSQ=R**2
      SSQ=S**2
      SUMR=0.0
C
C     RELATIVE R LOOP
C
      DO 400 MR=1,KR
      RS=AR(MR)*R0
      SS=RS*C3
      RS=RS*C4
      RSSQ=RS**2+RSQ
      RPROD=2.0*R*RS
      SSSQ=SS**2+SSQ
      SPROD=2.0*S*SS
C
C     RELATIVE R ANGLE LOOP
C
      SUMX=0.0
      DO 300 MX=1,KX
      X =AG(MX)
      y =BG(mx)
      IX=0
  110 CONTINUE
      R1=SQRT(RSSQ+RPROD*X)
      R2=SQRT(SSSQ+SPROD*X)
      FK1=R1/DR
      K1=FK1
      K1=MAX0(K1,2)
      FK1=FK1-FLOAT(K1)
      IF(K1.GT.K) GO TO 300
      FK2=R2/DR
      K2=FK2
      K2=MAX0(K2,2)
      FK2=FK2-FLOAT(K2)
      IF(K2.GT.K) GO TO 300
      COS1=(R+RS*X)/R1
      COS2=(S+SS*X)/R2
      SIN1=abs(rs*y/r1)
      SIN2=abs(ss*y/r2)
  120 CONTINUE
      DM1=0.0
      DL1=FACT1
      IF(L1.EQ.0) GO TO 141
      DO 140 LI=1,L1
      DL1=DL1*SIN1
  140 CONTINUE
  141 CONTINUE
      D1(L1+1)=DL1
      DM2=0.0
      DL2=FACT2
      IF(L2.EQ.0) GO TO 151
      DO 150 LI=1,L2
      DL2=DL2*SIN2
  150 CONTINUE
  151 CONTINUE
      D2(L2+1)=DL2
      IF(L1.EQ.0) GO TO 171
      FJ1=1.0
      FL1=LL1
      FM1=LL1
      DO 170 LI=1,L1
      DK1=(FM1*COS1*DL1/SIN1-DM1)/(FJ1*FL1)
      FJ1=FJ1+1.0
      FL1=FL1-1.0
      FM1=FM1-2.0
      DM1=DL1
      DL1=DK1
      INDX=L1+1-LI
      D1(INDX)=DL1
  170 CONTINUE
  171 CONTINUE
      IF(L2.EQ.0) GO TO 181
      FJ2=1.0
      FL2=LL2
      FM2=LL2
      DO 180 LI=1,L2
      DK2=(FM2*COS2*DL2/SIN2-DM2)/(FJ2*FL2)
      FJ2=FJ2+1.0
      FL2=FL2-1.0
      FM2=FM2-2.0
      DM2=DL2
      DL2=DK2
      INDX=L2+1-LI
      D2(INDX)=DL2
  180 CONTINUE
  181 CONTINUE
      PROD=0.0
      DO 185 LI=1,LPL
      PROD=PROD+D1(LI)*D2(LI)*TVCC(LI)
  185 CONTINUE
  280 CONTINUE
C     INTERPOLATE - 4 POINT FORMULA
      FT1=-FK1*(FK1-1.)*(FK1-2.)*F1(K1-1)/6.
     1  +   (FK1**2-1.)*(FK1-2.)*F1(K1  )/2.
     2  -  FK1*(FK1+1.)*(FK1-2.)*F1(K1+1)/2.
     3  +        FK1*(FK1**2-1.)*F1(K1+2)/6.
      FT2=-FK2*(FK2-1.)*(FK2-2.)*F2(K2-1)/6.
     1  +   (FK2**2-1.)*(FK2-2.)*F2(K2  )/2.
     2  -  FK2*(FK2+1.)*(FK2-2.)*F2(K2+1)/2.
     3  +        FK2*(FK2**2-1.)*F2(K2+2)/6.
C     INTERPOLATE - 2 POINT FORMULA
C     FT1=FK1*(F1(K1+1)-F1(K1))+F1(K1)
C     FT2=FK2*(F2(K2+1)-F2(K2))+F2(K2)
c       if(m.eq.21) write(20,'(2i4,1p10e12.4)')mr,mx,x,r1,r2,r1**2+r2**2
c     1    ,ft1,ft2,ft1*ft2
      SUMX=SUMX+WG(MX)*PROD*FT1*FT2
      IF(IX.NE.0) GO TO 300
      IX=1
      IF(IPAULI.EQ.0) then
        ITEMP=K1
        K1=K2
        K2=ITEMP
        ATEMP=FK1
        FK1=FK2
        FK2=ATEMP
        ATEMP=COS1
        COS1=COS2
        COS2=ATEMP
        ATEMP=SIN1
        SIN1=SIN2
        SIN2=ATEMP
        IF(L1.EQ.L2) GO TO 280
        GO TO 120
      endif
      X=-X
      GO TO 110
  300 CONTINUE
      SUMR=SUMR+WR(MR)*SUMX
  400 CONTINUE
      SUMR=SUMR*FNORM
      FR(2*M-1)=FR(2*M-1)+SUMR
      IF(M.EQ.KM) SL=SUMR
  500 CONTINUE
 1000 CONTINUE
      RETURN
      END


c***********************************************************************
      SUBROUTINE DWPLOT(NTH0,NLOG,SIGPLT,DTHETA,THETA1)
C
C     This is the new version of the printer plot routine.
c***********************************************************************
C
c     NTH0    Number of points to be plotted
c     NLOG    Number of decades in log plot
c             Positive value for log    scale
c             Negative value for linear scale
c     SIGPLT  Array containing numbers to be plotted
c     DTHETA  Increment in angle value
c     THETA1  Initial angle value
c
      implicit real*8(a-h,o-z)
      DIMENSION SIGPLT(60)
      CHARACTER*1 DSTORE(101),BLANK,ASTER
     1           ,APLUS ,FENCE,TEMP
      DATA LEN,pi/100,3.141592/
      data blank,aster,aplus,fence/' ','*','+',']'/
C
      nth = abs(nth0)
      IF(NTH .EQ.0.OR.NLOG.EQ.0) GO TO 2000
      if(nlog.gt.0) then
        MAXN=0
        DO 30 I=1,NTH
        IF(SIGPLT(I).LE.0.0) SIGPLT(I)=1.0E-20
        SIGPLT(I)= LOG10(SIGPLT(I))+100.0
        NDELT=SIGPLT(I)
        MAXN=MAX0(MAXN,NDELT+1)
   30   CONTINUE
        MAXN=MAXN-NLOG
        nlogp=nlog
      else
        do 35 i=1,nth
        sigplt(i)=5.0*sigplt(i)
   35   continue
        maxn=-5
        nlogp=iabs(nlog)
      endif
      FMIN=MAXN
      NDELT=LEN/NLOGP
      FDELT=NDELT
      IWRITE=1
      GO TO 1000
   50 CONTINUE
      tint = THETA1
      DO 100 I=1,NTH
C
      if(nth0.lt.0) then
        if    (tint.gt. 1.0) then
          tint = 1.0
        elseif(tint.lt.-1.0) then
          tint =-1.0
        endif
        theta = acos(tint)*180./pi
      else
        theta = tint
      endif
      cs = cos(theta*pi/180.)
      IF(I.EQ.1.OR.I.EQ.NTH) THEN
        TEMP=APLUS
      ELSE
        TEMP=BLANK
      ENDIF
      DSTORE(    1)=APLUS
      DSTORE(LEN+1)=APLUS
      DO 70 J=2,LEN
      DSTORE(J)=TEMP
   70 CONTINUE
      if(nlog.lt.0) dstore(len/2+1)=aplus
C
      IF(I.EQ.1.OR.I.EQ.NTH) THEN
        DSTORE(1)=FENCE
        DO 80 J=1,NLOGP
        INDX=J*NDELT
        DSTORE(INDX+1)=FENCE
   80   CONTINUE
      ENDIF
C
      INDX=(SIGPLT(I)-FMIN)*FDELT+0.5
      IF(INDX.GE.0) DSTORE(INDX+1)=ASTER
C
      WRITE(6,9002)THETA,DSTORE,THETA,cs
      IF(INDX.GE.0) DSTORE(INDX+1)=BLANK
      tint = tint + DTHETA
  100 CONTINUE
C
      IWRITE=0
 1000 CONTINUE
      IF(NLOG.GT.0) THEN
        WRITE(6,9000)
      ELSE
        WRITE(6,9001)
      ENDIF
      IF(IWRITE.EQ.1) GO TO 50
 2000 RETURN
 9000 FORMAT('    THETA  ',101X ,'   THETA  ','cos(theta)')
 9001 FORMAT('    THETA -1.',47X,'0.0',47X ,' 1.  THETA  ','cos(theta)')
 9002 FORMAT(' ',F8.2 ,'  ',101A1,   F8.2,  f10.3 )
      END


c***********************************************************************
      SUBROUTINE GAUSSR(NMAX,INDEX,ALFA,AG,WG,IERR,CUTOFF)
C
c     Gauss-Hermite and gauss-laguerre point and weight routine
c***********************************************************************
C
C     IF ALFA IS INTEGER       -- GAUSS-LAGUERRE
C     IF ALFA IS INTEGER + 1/2 -- GAUSS-HERMITE
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AG(100),WG(100)
      data eps/1.e-6/
      DATA PI,SQRPI/3.14159265,1.77245385/
c
      INDEX=0
      FI=NMAX
      FKI=4.0*(FI+(ALFA+1.0)*0.5)
      FLN= LOG(FI)
C
      K=ALFA
      JJ=(ALFA-FLOAT(K))*2.0
      DY=0.0
      IF(JJ.NE.0) GO TO 19
      FNORM=1.0
      IF(K.EQ.0) GO TO 11
      DO 10 J=1,K
      FNORM=FNORM*(ALFA+FI+1.0-FLOAT(J))
   10 CONTINUE
   11 CONTINUE
      GO TO 25
   19 CONTINUE
      FNORM=SQRPI/2.0
      DO 20 J=1,NMAX
      FNORM=FNORM*(1.0+ALFA/FLOAT(J))
   20 CONTINUE
      K=ALFA+1.0
      IF(K.EQ.0) GO TO 25
      DO 22 J=1,K
      FNORM=FNORM*(FLOAT(J)-0.5)
   22 CONTINUE
   25 CONTINUE
      Y=0.0
      Z1=0.0
      Z2=0.0
C
      DO 80 J=1,NMAX
      FJ=J
      FKJ=4.0*(FJ+(ALFA+1.0)*0.5)
      Z=((FJ+0.5*ALFA-0.25)*PI)**2/FKI
      Z=Z*(1.0+Z*(1.0+PI* LOG(FJ)*FLN/(8.0*(FI+FLN-FJ+eps)))/(3.0*FKI))
      DELZ=Z-Z1
      Z1=Z
      Z=Y+DELZ
      Z3=Z
      DO 74 M=1,20
      A1=0.0
      A2=1.0
      DO 70 K=1,NMAX
      FK=K
      A3=((2.0*FK+ALFA-1.0-Z)*A2-(FK+ALFA-1.0)*A1)/FK
      A1=A2
      A2=A3
      B2=(FK*A2-(FK+ALFA)*A1)/Z
   70 CONTINUE
      Y=Z-A2/B2
      IF(ABS((Z-Y)/Z).LT.3.E-14) GO TO 75
      Z=Y
   74 CONTINUE
      IERR=J
   75 CONTINUE
      DZ=Z-Z2
      Z2=Z
      DZ=((FNORM/Z)/B2)/B2
      IF(DZ.LT.CUTOFF.AND.DZ.LT.DY) GO TO 100
      DY=DZ
      INDEX=J
      AG(J)=Y
      WG(J)=DZ
      IF(JJ.NE.0) AG(J)=SQRT(Y)
   80 CONTINUE
  100 CONTINUE
      RETURN
      END


c***********************************************************************
      SUBROUTINE LEGAUS(LL,X,W)
C
c     Gauss-Legendre point and weight routine for an even no. of points.
c***********************************************************************
C
C       ll= order i.e. the total number of points  -1< x <1
C       x = points, calculates and stores only for  0< x <1
C       w = weights
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(100),W(100)
c
c     starting values for recursion
c
      Z3=-1.6/FLOAT(LL+1)
      Z2=3.0*Z3
      Z1=5.0*Z3
c
      NL=(LL+1)/2
      DO 200 L=1,NL
      ZOLD=0.0
      Z=Z1+3.0*(Z3-Z2)
      DO 50 J=1,10
      P1=0.0
      P2=1.0
      DO 30 I=1,LL
      P3=(FLOAT(I+I-1)*Z*P2-FLOAT(I-1)*P1)/FLOAT(I)
      P1=P2
      P2=P3
   30 CONTINUE
      DP=FLOAT(LL)*(P1-Z*P2)/(1.0-Z*Z)
      Z=Z-P2/DP
      IF(ABS(Z-ZOLD)/Z.LT.1.0E-10) GO TO 51
      ZOLD=Z
   50 CONTINUE
      WRITE(6,9100)L
 9100 FORMAT(28H0NO CONVERGENCE FOR ZERO NO.,I4)
   51 CONTINUE
      X(L)=Z
      W(L)=2.0/((1.0-Z*Z)*DP*DP)
      Z1=Z2
      Z2=Z3
      Z3=Z
  200 CONTINUE
      RETURN
      END


c***********************************************************************
      SUBROUTINE LGNDR(PLM,MPLUS,LPLUS,THET)
c
c     Calculates Legendre polynomials Plm
c
c     mplus   number of m's >0
c     lplus   number of l's >0
c     thet    angle in degrees
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION PLM(459)
c
      THETA=THET /57.295779
      Y=COS(THETA)
      Z=SIN(THETA)
      IX=0
      DO 100 M=1,MPLUS
      LX=M-1
      L2=0
      P3=1.0
      FL1=LX
      IF(LX.EQ.0) GO TO 41
      DO 40 LT=1,LX
      FL1=FL1+1.0
      P3=P3*FL1*Z/2.0
   40 CONTINUE
   41 P2=0.0
      FL2=FL1+1.0
      FL3=1.0
      DO 90 LT=1,LPLUS
      IX1=IX+LT
      IF(L2-LX)50,70,60
   50 PLM(IX1)=0.0
      GO TO 75
   60 P3=(FL2*Y*P2-FL1*P1)/FL3
      FL1=FL1+1.0
      FL2=FL2+2.0
      FL3=FL3+1.0
   70 PLM(IX1)=P3
      P1=P2
      P2=P3
   75 L2=L2+1
   90 CONTINUE
      IX=IX+LPLUS
  100 CONTINUE
      RETURN
      END


c***********************************************************************
      SUBROUTINE POLFCT(POL,SUMR,SUMI,JSX,JRX)
c
c     Calculates spin observables, polarization, asymmetry and tensor.
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SUMR(5,5),SUMI(5,5),B(2,2),C(3,3),TSUM(10),POL(10)
     1,         s22(3,3),s21(3,3),s20(3,3),syy(3,3)
C
C     DATA FOR SIGMA-Y
C
C     SY MATRIX FOR SPIN-1/2
      DATA B/0.0,-1.0,1.0,0.0/
C
C     SY MATRIX FOR SPIN-1
      DATA C/0.0,-.707106781,0.0,.707106781,0.0,-.707106781
     1      ,0.0, .707106781,0.0/
C
C     SYY = 3*SY*SY - 2  MATRIX FOR AYY
      DATA SYY/-0.5,0.0,-1.5,0.0,1.0,0.0,-1.5,0.0,-0.5/
C
C     S22 = (   S^*S^        )*SQRT(3)/4
C     S21 = (   S^*SZ + SZ*S^)*SQRT(3)/2
C     S20 = ( 3*SZ*SZ -  S*S )/SQRT(2)
C
      DATA S22/0.0,0.0,1.73205081,0.0,0.0,0.0,0.0,0.0,0.0/
      DATA S21/0.0,-1.2247449,0.0,0.0,0.0,1.2247449,0.0,0.0,0.0/
      DATA S20/0.70710678,0.0,0.0,0.0,-1.4124214,0.0,0.0,0.0,0.70710678/
C
      DO 2 I=1,10
      POL(I)=0.0
    2 CONTINUE
      m=0
    5 continue
      if    (m.eq.0) then
        JS=JSX
        JR=JRX
        LDEL=0
        L1 = 4
        M=1
      elseif(m.eq.1) then
        JS=JRX
        JR=JSX
        LDEL=4
        L1 = 6
        M=2
      endif
c
      DO 100 I=1,JR
      DO 90 J=1,JS
      DO 7 K=1,10
      TSUM(K)=0.0
    7 CONTINUE
      DO 80 K=1,JS
      IF(JS.GT.3) GO TO 80
      if    (js.eq.1) then
        TEMP=1.0
      elseif(js.eq.2) then
        TEMP=B(J,K)
      elseif(js.eq.3) then
        TEMP=C(J,K)
      endif
   50 CONTINUE
      if    (m.eq.1) then
        T1=SUMR(K,I)
        T2=SUMI(K,I)
      elseif(m.eq.2) then
        T1=SUMR(I,K)
        T2=SUMI(I,K)
      endif
c     polarization or asymmetry
      TSUM(1)=TSUM(1)-T2*TEMP
      TSUM(2)=TSUM(2)+T1*TEMP
c
      IF(JS.eq.3) then
c      DO 75 L=1,3
c      LL=L+L+1
c      TEMP=D(J,K,L)
c      TSUM(LL  )=TSUM(LL  )+T1*TEMP
c      TSUM(LL+1)=TSUM(LL+1)+T2*TEMP
c   75 CONTINUE
c     tensor polarizations and asymmetries
      tsum( 3)=tsum( 3)+t1*s22(j,k)
      tsum( 4)=tsum( 4)+t2*s22(j,k)
      tsum( 5)=tsum( 5)+t1*s21(j,k)
      tsum( 6)=tsum( 6)+t2*s21(j,k)
      tsum( 7)=tsum( 7)+t1*s20(j,k)
      tsum( 8)=tsum( 8)+t2*s20(j,k)
      tsum( 9)=tsum( 9)+t1*syy(j,k)
      tsum(10)=tsum(10)+t2*syy(j,k)
      endif
   80 continue
c
      if    (m.eq.1) then
        T1=SUMR(J,I)
        T2=SUMI(J,I)
      elseif(m.eq.2) then
        T1=SUMR(I,J)
        T2=SUMI(I,J)
      endif
c
      DO 85 L=1,L1
      LL=L+L-1
      K=L+LDEL
      POL(K)=POL(K)+TSUM(LL  )*T1+TSUM(LL+1)*T2
   85 CONTINUE
   90 CONTINUE
  100 CONTINUE
      GO TO (5,105),M
  105 POL(7)=-POL(7)
      POL(5)=-POL(5)
      RETURN
      END


c***********************************************************************
      SUBROUTINE POTS(U,V)
c
c     Calculates the potentials or form factors
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/POTTER/DRX,AFACT(2),VFACT,SFACT,ENG,RM,G(4),ETAX,FKX,ETAKX
     1             ,RCX,HBARC2,ABETA(3),FLDF(3)
     2             ,NX,LAMX,KMXX,KX,IBX,LPLUSX,ICON4,NSPC,IDIRAC,ICHK

      CHARACTER*18 B(21),C(22)
      DIMENSION UT(5),CN(16),CP(16),YLAM(16),PLM(8)
     1         ,XG(8),WG(8),U(800),V(800),LDFRM(3)
C      DIMENSION UD(800)
      EQUIVALENCE (B(1),C(2)),(YLAM(1),PLM(1))
C
      DATA C( 1) /' NX=0 No potential'/
      DATA B( 1) /' NX=1 VOLUME  W-S '/
      DATA B( 2) /' NX=2 SURFACE W-S '/
      DATA B( 3) /' NX=3 2ND DERIV   '/
      DATA B( 4) /' NX=4 L*S VOLUME  '/
      DATA B( 5) /' NX=5 L*S SURFACE '/
      DATA B( 6) /' NX=6 VOL *R**POWR'/
      DATA B( 7) /' NX=7 SURF*R**POWR'/
      DATA B( 8) /' NX=8 EXTERN FORMF'/
      DATA B( 9) /' NX=9 HARMONIC OSC'/
      DATA B(10) /' NX=10 GAUSSIAN   '/
      DATA B(11) /' NX=11 DEFORM VOL '/
      DATA B(12) /' NX=12 DEFORM SURF'/
      DATA B(13) /' NX=13 HULTHEN    '/
      DATA B(14) /' NX=14 NO OPTION  '/
      DATA B(15) /' NX=15 NO OPTION  '/
      DATA B(16) /' NX=16 NO OPTION  '/
      DATA B(17) /' NX=17 NO OPTION  '/
      DATA B(18) /' NX=18 Volume L*S '/
      DATA B(19) /' NX=19 COUL EXCIT '/
      DATA B(20) /' NX=20 VECTOR     '/
      DATA B(21) /' NX=21 SCALAR     '/

C
      DATA NG,NGX/8,0/
C
      ETA4 = 6.0
      ETA5 =10.0
      SQRPI= 1.772453851
      PI=3.141592654
C
      IVFLAG=0
      ISFLAG=0
      FACT=VFACT
   70 CONTINUE
C
C     READ IN CARD SET 5,6,OR 7   POTENTIAL CARDS
C
      READ (5,9000)FZ,VR,RY,AR,VSOR,VI,RZ,AI,VSOI,PWR
C
      NZ=ABS(FZ)
      RR=ABS(RY)*AFACT(1)
      IF(RY.LT.0.0) RR=RR+ABS(RY)*AFACT(2)
      RI=ABS(RZ)*AFACT(1)
      IF(RZ.LT.0.0) RI=RI+ABS(RZ)*AFACT(2)
      IF(ICON4.NE.2.OR.NSPC.LT.3) THEN
        WRITE(6,9509)B(NZ),VR,RY,AR,RR,VSOR
        WRITE(6,9510)                VI,RZ,AI,RI,VSOI,PWR
      ENDIF
C     NZ= 1 VOLUME WOODS-SAXON
C     NZ= 2 SURFACE WOODS-SAXON
C     NZ= 3 SECOND DERIVATIVE WOODS-SAXON
C     NZ= 4 L*S POTENTIAL FOR WOODS-SAXON VOLUME
C     NZ= 5 L*S POTENTIAL FOR WOODS-SAXON SURFACE
C     NZ= 6 WOOD-SAXON*R**PWR
C     NZ= 7 1ST-DER WOOD-SAXON*R**PWR
C     NZ= 8 FORMF8 EXTERNAL FORM FACTOR
C     NZ= 9 HARMONIC OSCILLATOR
C     NZ=10 GAUSSIAN*R**PWR
C     NZ=11 DEFORMED VOLUME
C     NZ=12 DEFORMED SURFACE
C     NZ=13 Hulthen
C     Nz=18 Volume L*S
C     NZ=19 COULOMB-EXCITATION
C     NZ=20 VECTOR DIRAC
C     NZ=21 SCALAR DIRAC
C
      KFLAG=0
      IF(NX.LT.3) THEN
        VR=VR*FACT
        VI=VI*FACT
        KT=FKX  * MAX (RR,RI)+ETA5
        LPLUSX=MAX0(LPLUSX,KT)
        KT=(2.3*ETA4* MAX (AR,AI)+ MAX (RR,RI))/DRX
      ELSE
        IF(ENG.EQ.0.0) THEN
          KT=(2.3*ETA4* MAX (AR,AI)+ MAX (RR,RI))/DRX
        ELSE
          RM= MAX (RM,RR)
          RM= MAX (RM,RI)
          IF(RM.EQ.0.0) RM=1.0
          VR=VR*FACT
          VI=VI*FACT
          KT=(2.3*ETA4/SQRT(FKX**2+2.0*ETAX*FKX/RM))/DRX
        ENDIF
      ENDIF
      KX=MIN0(MAX0(KX,KT),KMXX)
   83 CONTINUE
      IF(AR.EQ.0.0) THEN
        F1=0.0
        F2=0.0
      ELSE
        F2=EXP(-DRX  /AR)
        F1=EXP( RR/AR)
      ENDIF
      IF(AI.EQ.0.0) THEN
        F3=0.0
        F4=0.0
      ELSE
        F4=EXP(-DRX  /AI)
        F3=EXP( RI/AI)
      ENDIF
C
      IF(NX.GE.3.AND.ENG.EQ.0.0.AND.NZ.LE.5) THEN
        IF(AR.NE.0.0) VR=VR*(RR/AR)**(PWR+1)
        IF(AI.NE.0.0) VI=VI*(RI/AI)**(PWR+1)
      ENDIF
      IF(NZ.EQ.0) GO TO 6050
      GO TO (100,200,300,400,500,600,700,800,900,1000,1100,1200
     1,1300,1400,1500,1600,1700,1800,1900,2000,2100),NZ
C
C     VOLUME WOODS SAXON
C
  100 CONTINUE
      DO 160 M=1,KX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )-VR*F1/(1.0+F1)
      U(MK+1)=U(MK+1)-VI*F3/(1.0+F3)
  160 CONTINUE
      GO TO 6000
C
C     1ST DERIVATIVE WOODS SAXON
C
  200 CONTINUE
      DO 260 M=1,KX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )+VR*F1/(1.0+F1)**2
      U(MK+1)=U(MK+1)+VI*F3/(1.0+F3)**2
  260 CONTINUE
      GO TO 6000
C
C     2ND DERIVATIVE WOODS SAXON
C
  300 CONTINUE
      DO 360 M=1,KX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )-VR*F1*(1.0-F1)/(1.0+F1)**3
      U(MK+1)=U(MK+1)-VI*F3*(1.0-F3)/(1.0+F3)**3
  360 CONTINUE
      GO TO 6000
C
C     L.S VOLUME WOODS SAXON
C
  400 CONTINUE
      IBX=1
      IF(AR.NE.0.0) VR=VR/AR
      IF(AI.NE.0.0) VI=VI/AI
      R=0.0
      DO 460 M=1,KX
      R=R+DRX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      V(MK  )=V(MK  )-VR*F1/(R*(1.0+F1)**2)
      V(MK+1)=V(MK+1)-VI*F3/(R*(1.0+F3)**2)
  460 CONTINUE
      GO TO 6000
C
C     L.S 1ST DERIVATIVE WOODS SAXON
C
  500 CONTINUE
      IBX=1
      IF(AR.NE.0.0) VR=VR/AR
      IF(AI.NE.0.0) VI=VI/AI
      R=0.0
      DO 560 M=1,KX
      R=R+DRX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      V(MK  )=V(MK  )+VR*F1*(1.0-F1)/(R*(1.0+F1)**3)
      V(MK+1)=V(MK+1)+VI*F3*(1.0-F3)/(R*(1.0+F3)**3)
  560 CONTINUE
      GO TO 6000
C
C     WOOD SAXON*R**PWR
C
  600 CONTINUE
      R=0.0
      DO 660 M=1,KX
      R=R+DRX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      RPWR=R**PWR
      U(MK  )=U(MK  )-VR*F1*RPWR/(1.0+F1)
      U(MK+1)=U(MK+1)-VI*F3*RPWR/(1.0+F3)
  660 CONTINUE
      GO TO 6000
C
C      1ST DERIVATIVE WOOD SAXON*R**PWR
C
  700 CONTINUE
      R=0.0
      DO 760 M=1,KX
      MK=M+M-1
      R=R+DRX
      F1=F1*F2
      F3=F3*F4
      RPWR=R**PWR
      U(MK  )=U(MK  )+VR*F1*RPWR/(1.0+F1)**2
      U(MK+1)=U(MK+1)+VI*F3*RPWR/(1.0+F3)**2
  760 CONTINUE
      GO TO 6000
C
C     EXTERNAL FORM FACTOR
C
  800 CONTINUE
      IF(NX.GE.3) THEN
        READ (5,9000)G
        WRITE(6,9508)G
      ENDIF
      READ(5,9000)F1,F2
C
C     F2 = 0  REAL CENTRAL
C     F2 = 1  IMAG CENTRAL
C     F2 = 2  REAL SPIN ORBIT
C     F2 = 3  IMAG SPIN ORBIT
C
      IF(F2.EQ.0.0.OR.F2.EQ.2) THEN
        F3=VR
        MK=1
      ELSE
        F3=VI
        MK=2
      ENDIF
      IF(F3.EQ.0.0) F3=1.0
      KT=F1
      DO 820 M=1,KT,5
      READ(5,9100)UT
      IF(M.GT.KMXX) GO TO 820
      DO 810 I=1,5
      IF(F2.LT.2.0) THEN
        U(MK  )=U(MK  )+UT(I)*F3
      ELSE
        V(MK  )=V(MK  )+UT(I)*F3
      ENDIF
      MK=MK+2
  810 CONTINUE
  820 CONTINUE
C
      KX=min0(KT,kmxx)

      GO TO 6000
C
C     HARMONIC OSCILLATOR, NORMALIZED   N*EXP(-(R/RY)**2/2)
C
  900 CONTINUE
      READ (5,9000)G
      WRITE(6,9508)G
      F1=1.0/RY**2
      F2=F1/RY
      F3=0.5
      F4=SQRPI*0.5
      L=G(2)
      IF(L.NE.0) THEN
        DO 930 I=1,L
        F3=F3+1.0
        F4=F4*F3
        F2=F2*F1
  930   CONTINUE
      ENDIF
      NN=G(1)
      T1=1.0
      T2=F4
C     LAGUERRE POLYNOMIAL COEFFICIENTS = (Abramowitz-Stegun)*(-1)**nn
      CN(1)=(-1.0)**NN
      IF(NN.NE.0) THEN
        DO 940 I=1,NN
        F3=F3+1.0
        T1=T1*FLOAT(I)
        T2=T2*F3
        CN(I+1)=-CN(I)*F1*FLOAT(NN+1-I)/(FLOAT(I)*F3)
  940   CONTINUE
      ENDIF
      ANORM=SQRT(2.0*F2*T1/T2)*T2/(T1*F4)
      IF(VR.NE.0.0) ANORM=VR*ANORM
      KT=10.0*RY/DRX
      KT=MIN0(KT,KMXX)
      R=0.0
      F1=F1/2.0
      DO 960 M=1,KT
      MK=M+M-1
      R=R+DRX
      R2=R**2
      F2=CN(1)
      IF(NN.EQ.0) GO TO 951
      F3=1.0
      DO 950 I=1,NN
      F3=F3*R2
      F2=F2+CN(I+1)*F3
  950 CONTINUE
  951 CONTINUE
      U(MK  )=U(MK  )+F2*ANORM*EXP(-F1*R2)*R**L
  960 CONTINUE
      GO TO 6000
C
C     EXP(-(R/R0)**2)*R**POWR
C
 1000 CONTINUE
      IF(VR.NE.0.0) THEN
        R=0.0
        DO 1060 M=1,KX
        MK=M+M-1
        R=R+DRX
        U(MK  )=U(MK  )-VR*EXP(-(R/RY)**2)*R**PWR
 1060 CONTINUE
      endif
      IF(VI.NE.0.0) THEN
        R=0.0
        DO 1070 M=1,KX
        MK=M+M-1
        R=R+DRX
        U(MK+1)=U(MK+1)-VI*EXP(-(R/RZ)**2)*R**PWR
 1070 CONTINUE
      endif
      GO TO 6000
C
C     DEFORMED VOLUME OR SURFACE BY YLM EXPANSION
C
 1100 CONTINUE
      IF(NGX.NE.NG) THEN
        CALL LEGAUS(2*NG,XG,WG)
        NGX=NG
      ENDIF
      T2=(-1.0)**LAMX
      IF(ICHK.EQ.0) THEN
        READ (5,9000) (ABETA(J),FLDF(J) ,J=1,3)
        WRITE(6,9000)
        WRITE(6,9512) (ABETA(J),FLDF(J) ,J=1,3)
      ENDIF
      LMAX=LAMX+1
      DO 1101 J=1,3
      LDFRM(J)=FLDF(J)
      LMAX=MAX0(LMAX,LDFRM(J)+1)
 1101 CONTINUE
      T2=(-1.0)**LAMX
C
      DO 1108 I=1,NG
      CN(I  )=0.0
      CN(I+NG)=0.0
      P2=0.0
      P3=1.0
C
      DO 1106 M=1,LMAX
      L=M-1
      FL=L-1
      IF(L.EQ.0) GO TO 1102
      P3=((2.0*FL+1.0)*XG(I)*P2-FL*P1)/(FL+1.0)
 1102 CONTINUE
      DO 1103 J=1,3
      IF(ABETA(J).EQ.0.0) GO TO 1103
      IF(L.NE.LDFRM(J)) GO TO 1103
      FACTOR=P3*ABETA(J)*SQRT(FLOAT(L+L+1)/(4.0*PI))
      CN(I  )=CN(I  )+FACTOR
      CN(I+NG)=CN(I+8)+FACTOR*(-1.0)**LDFRM(J)
 1103 CONTINUE
      IF(L.NE.LAMX) GO TO 1104
      YLAM(I  )= P3*WG(I)*SQRT(FLOAT(L+L+1)*PI)
      YLAM(I+NG)=YLAM(I  )*T2
 1104 CONTINUE
      P1=P2
      P2=P3
 1106 CONTINUE
 1108 CONTINUE
      DO 1170 I=1,2
      IF(I.EQ.1) THEN
        IF(VR.EQ.0.0) GO TO 1170
        VX=VR
        RX=RR
        AX=AR
        F1=1.0
        IFL=-1
      ELSE
        IF(VI.EQ.0.0) GO TO 1170
        VX=VI
        RX=RI
        AX=AI
        F1=1.0
        F2=F4
        IFL=0
      ENDIF
      DO 1135 J=1,16
      CP(J)=EXP((1.0+CN(J))*RX/AX)
 1135 CONTINUE
      IF(LAMX.EQ.0) VX=VX/(SQRPI*2.0)
      J=NZ-10
      IF(J.EQ.2) GO TO 1151
      DO 1150 M=1,KX
      MK=M+M+IFL
      VTEMP=0.0
      F1=F1*F2
      DO 1145 J=1,16
      F3=F1*CP(J)
      VTEMP=VTEMP-YLAM(J)*VX*F3/(1.0+F3)
 1145 CONTINUE
      U(MK  )=U(MK  )+VTEMP
 1150 CONTINUE
      GO TO 1170
 1151 CONTINUE
      DO 1160 M=1,KX
      MK=M+M+IFL
      VTEMP=0.0
      F1=F1*F2
      DO 1155 J=1,16
      F3=F1*CP(J)
      VTEMP=VTEMP+YLAM(J)*VX*F3/(1.0+F3)**2
 1155 CONTINUE
      U(MK  )=U(MK  )+VTEMP
 1160 CONTINUE
 1170 CONTINUE
      GO TO 6000
 1200 CONTINUE
      GO TO 1100
C
C     HULTHEN WAVE FUNCTION= NORM*(EXP(-R/RY)-EXP(-R/RZ))/R
C
 1300 CONTINUE
      READ (5,9000)G
      WRITE(6,9508)G
      T1=1.0/RY
      T2=1.0/RZ
      T3=T2**2-T1**2
      T4=SQRT(2.0*T1*(T1+T2)*T2)/(T2-T1)
      IF(VR.NE.0.0) KT=16.0* MIN (RY,RZ)/DRX
      IF(VR.EQ.0.0) KT=16.0* MAX (RY,RZ)/DRX
      KX=MIN0(KX,KMXX)
      F1=1.0
      F2=EXP(-DRX  *T1)
      F3=1.0
      F4=EXP(-DRX  *T2)
      R=0.0
      DO 1260 M=1,KX
      MK=M+M-1
      R=R+DRX
      F1=F1*F2
      F3=F3*F4
      V(M)=T3*F3/(F1-F3)
      TEMP=1.0
      IF(VR.NE.0.0) TEMP=V(M)/FACT
      U(MK  )=U(MK  )+TEMP*T4*(F1-F3)/R
 1260 CONTINUE
      GO TO 6000
 1400 CONTINUE
 1500 CONTINUE
 1600 CONTINUE
 1700 CONTINUE
       go to 6000
 1800 CONTINUE
      DO 1860 M=1,KX
      MK=M+M-1
      F1=F1*F2
      F3=F3*F4
      V(MK  )=V(MK  )-VR*F1/(1.0+F1)
      V(MK+1)=V(MK+1)-VI*F3/(1.0+F3)
 1860 CONTINUE
      GO TO 6000
 1900 CONTINUE
c
c     Coulomb excitation for a deformed uniform charge distribution
c
      IF(ICHK.EQ.0) THEN
        READ (5,9000) (ABETA(J),FLDF(J) ,J=1,3)
        WRITE(6,9000)
        WRITE(6,9512) (ABETA(J),FLDF(J) ,J=1,3)
      ENDIF
      if(vr.eq.0.0) vr=1.0
      do 1990 i=1,3
      beta=abeta(i)
      if(beta.ne.0.0) then
c     set flag for unbound stripping evaluation of coulomb ex.
      ibx=4
      flam1=beta*sqrt(float(2*lamx+1)/(4.0*pi))
      flam2=flam1*phasef(lamx)
c
      if(ngx.ne.ng) then
        ngx=ng
        call legaus(2*ng,xg,wg)
      endif
c
      an=0.0
      do 1920 k=1,ng
      p1=0.0
      p2=1.0
      if(lamx.gt.0) then
        do 1910 j=1,lamx
        p3=(float(2*j-1)*xg(k)*p2 - float(j-1)*p1)/float(j)
        p1=p2
        p2=p3
 1910 continue
      endif
      plm(k)=p2
c
c     calculate normalization
c
      r1=rcx*(1.0+flam1*p2)
      r2=rcx*(1.0+flam2*p2)
      an=an+(r1**3+r2**3)*wg(k)
 1920 continue
c
      an =an/3.0
c
      r  =0.0
c
      do 1940 m=1,kx
      mk=m+m-1
      r=r+drx
      sum=0.0
      do 1930 k=1,ng
      r1=rcx*(1.0+flam1*plm(k))
      r2=rcx*(1.0+flam2*plm(k))
c
      if(r.gt.r1) then
        s1=r1**(lamx+3)/(r**(lamx+1)*float(lamx+3))
      else
        if(lamx.eq.2) then
          s1=( log(r1/r)+1.0/float(lamx+3))*r**2
        else
          s1=(r**lamx/r1**(lamx-2)
     1      -float(2*lamx+1)*r**2/float(lamx+3))/float(2-lamx)
        endif
      endif
c
      if(r.gt.r2) then
        s2=r2**(lamx+3)/(r**(lamx+1)*float(lamx+3))
      else
        if(lamx.eq.2) then
          s2=( log(r2/r)+1.0/float(lamx+3))*r**2
        else
          s2=(r**lamx/r2**(lamx-2)
     1      -float(2*lamx+1)*r**2/float(lamx+3))/float(2-lamx)
        endif
      endif
c
      sum=sum+(s1+s2*phasef(lamx))*wg(k)*plm(k)
 1930 continue
c
      sum=sum/an
      if(beta.ne.0.0.and.lamx.ne.0) sum=sum/flam1
      u(mk  )=u(mk  )+sum*(vr*etakx/fact)
 1940 continue
      endif
 1990 continue
      GO TO 6000
C
C     VECTOR POTENTIAL
C     VSOR, VSOI ARE THE THIRD PARAMETERS
C     IN THE 3 PARAMETER FERMI MODEL
C         Õ1 + VSO?*(R/R?)**2å
C
 2000 CONTINUE
      IVFLAG=-1
      R=0.0
      DO 2060 M=1,KX
      MK=M+M-1
      R=R+DRX
      F1=F1*F2
      F3=F3*F4
      U(MK  )=U(MK  )-VR*F1*(1.0+VSOR*(R/RR)**2)/(1.0+F1)
      U(MK+1)=U(MK+1)-VI*F3*(1.0+VSOI*(R/RI)**2)/(1.0+F3)
 2060 CONTINUE
      GO TO 6020
C
C     SCALAR POTENTIAL
C     VSOR, VSOI ARE THE FERMI THIRD PARAMETERS
C     IN THE 3 PARAMETER FERMI MODEL
C         Õ1 + VSO?*(R/R?)**2å
C
 2100 CONTINUE
      ISFLAG=-1
      VR=VR*(SFACT/VFACT)
      VI=VI*(SFACT/VFACT)
      R=0.0
      DO 2160 M=1,KX
      MK=M+M-1
      R=R+DRX
      F1=F1*F2
      F3=F3*F4
      V(MK  )=V(MK  )-VR*F1*(1.0+VSOR*(R/RR)**2)/(1.0+F1)
      V(MK+1)=V(MK+1)-VI*F3*(1.0+VSOI*(R/RI)**2)/(1.0+F3)
 2160 CONTINUE
      GO TO 6020
C
C      END OF POTENTIALS
C
 6000 CONTINUE
      IDIRAC=1
 6020 CONTINUE
      IF(KFLAG.NE.0.OR.NZ.GT.5) GO TO 6050
      IF(ABS(VSOR)+ABS(VSOI).EQ.0.0) GO TO 6050
      NZ=NZ+3
      KFLAG=1
      VR=VR*VSOR/45.2
      VI=VI*VSOI/45.2
      GO TO 83
 6050 CONTINUE
      IF(FZ.GT.0.0) GO TO 70
C
C     PROCESS DIRAC POTENTIALS
C
C     ENTRY WITH  U -> K**2 - VFACT*V   VFACT = 2.0*W1 /HBARC**2
C                 V ->      - SFACT*S   SFACT = 2.0*FM1/HBARC**2
C
      IF(IVFLAG.NE.0.AND.ISFLAG.NE.0) THEN
      IF(IDIRAC.EQ.1) THEN
        WRITE(6,9515)
      ENDIF
      IDIRAC=-1
C
C      KT2=KX+KX
C      WRITE(20,7777)' ENTRY  POTENTIALS'
C      WRITE(20,7778)(U(I),I=1,KT2)
C      WRITE(20,7777)
C      WRITE(20,7778)(V(I),I=1,KT2)
      W1M1=(VFACT+SFACT)*0.5*HBARC2

      DO 6100 M=1,Kx
      MK=M+M-1
      VVR=(U(MK  )-FKX**2)/VFACT
      VVI= U(MK+1)        /VFACT
      VSR= V(MK  )        /SFACT
      VSI= V(MK+1)        /SFACT
      U(MK  )=U(MK  )+V(MK  ) +(VVR**2-VVI**2 - VSR**2+VSI**2)/HBARC2
      U(MK+1)=U(MK+1)+V(MK+1) +(VVR*VVI - VSR*VSI)*2.0/HBARC2
      T1 = W1M1 + VVR-VSR
      T2 =        VVI-VSI
      V(MK  )=0.5* LOG (T1**2 + T2**2)
      V(MK+1)=ATAN2(T2,T1)
 6100 CONTINUE
C      WRITE(20,7777)' SECOND POTENTIALS'
C      WRITE(20,7778)(U(I),I=1,KT2)
C      WRITE(20,7777)
C      WRITE(20,7778)(V(I),I=1,KT2)
C
      R=FLOAT(KX+1)*DRX
      MK=KX+KX-1
      D2=V(MK  )
      D1=D2
      A2=V(MK+1)
      A1=A2
      IBX=1
      DO 6150 M=2,KX
      MK=MK-2
      R=R-DRX
      D3=D2
      D2=D1
      D1=V(MK  )
      A3=A2
      A2=A1
      A1=V(MK+1)
C     FIRST  DERIVATIVE TERMS
      DPR=(D3-D1)/(2.0*DRX)
      DPI=(A3-A1)/(2.0*DRX)
      V(MK+2)=2.0*DPR/R
      V(MK+3)=2.0*DPI/R
C     SECOND DERIVATIVE TERMS
      DPPR=(D3+D1-2.0*D2)/DRX**2
      DPPI=(A3+A1-2.0*A2)/DRX**2
      UDR =0.5*DPPR-0.25*(DPR**2-DPI**2)+DPR/R
      UDI =0.5*DPPI-0.25*(2.0*DPR*DPI  )+DPI/R
C      UD(MK+2)=UDR
C      UD(MK+3)=UDI
      U(MK+2)=U(MK+2)+UDR
      U(MK+3)=U(MK+3)+UDI
 6150 CONTINUE
      V(1   )=V(3   )*2.0
      V(2   )=V(4   )*2.0
      U(1   )=U(1   )+UDR
      U(2   )=U(2   )+UDI
C      WRITE(20,7777)'TERTIARY POTENTIALS'
C      WRITE(20,7777)'CENTRAL POTENTIAL'
C      WRITE(20,7778)(U(I),I=1,KT2)
C      WRITE(20,7777)'SPIN ORBIT POTENTIAL'
C      WRITE(20,7778)(V(I),I=1,KT2)
C      WRITE(20,7777)'UD - DARWIN TERM'
C      WRITE(20,7778)(UD(I),I=1,KT2)
C 7777 FORMAT(A30)
C 7778 FORMAT(' ',1P10E12.4)
      ENDIF
C
      IF(IDIRAC.EQ.1) IDIRAC=0
      RETURN
C
 9000 FORMAT(10F8.4)
 9100 FORMAT(5E16.7)
 9508 FORMAT(18X,9H   NODES=,F9.4,9H   L    =,F9.4,9H   2*J  =,F9.4
     1          ,9H   2*S  =,F9.4,9H   FISW =,F9.4)
 9509 FORMAT(A18,9H   V RL =,F9.4,9H   R0RL =,F9.4,9H   A RL =,F9.4
     1          ,9H   R RL =,F9.4,9H   VSOR =,F9.4)
 9510 FORMAT(18X,9H   V IM =,F9.4,9H   R0IM =,F9.4,9H   A IM =,F9.4
     1          ,9H   R IM =,F9.4,9H   VSOI =,F9.4,9H   POWR =,F9.4)
 9512 FORMAT(18X,9H   BETA1=,F9.4,9H   LDFR1=,F9.4,9H   BETA2=,F9.4
     1          ,9H   LDFR2=,F9.4,9H   BETA3=,F9.4,9H   LDFR3=,F9.4)
 9515 FORMAT('0',20('*'),'  WARNING, Mixing of Dirac and non-Dirac'
     1,          ' potentials may be hazardous to your calculation'
     2,          20('*'))
      END

c***********************************************************************
      SUBROUTINE ELSIG(dtemp,d,plm,sigplt,angle,fk,eta,rsig,alpha
     1                ,idat,is,icon,lplus)
c
c     Subroutine to print out elastic scattering cross sections
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION dtemp(5600),d(800),plm(1000),SIGPLT(200,2)
     1,angle(3),fk(2),eta(2),rsig(2),alpha(15),idat(6),is(2),icon(20)
     2,SIG(2),RATIO(2),POL(12),INCR(2),SUM(5,5,2),sumx(50)
      Equivalence (sum,sumx)

      DATA NX,pi/201,3.141592654/
C
      NTHETA=angle(1)
      theta1=angle(2)
      dtheta=angle(3)
      IF(NTHETA.EQ.0) return
      lpl2=2*lplus
C
C     CONSTRUCT TABLE OF SCATTERING AMPLITUDES WITH ANG-MOM COEFFICENTS
C
      ICNT=LPL2*((is(1)+1)**2+(is(2)+1)**2)
      DO 3 IP=1,ICNT
      DTEMP(IP)=0.0
    3 CONTINUE
      JR=1
      JS=0
      ICNT=0
      IC1=0
      DO 26 N=1,2
      ISS=is(N)+1
      JS=JS+ISS
      INCR(N)=LPLUS*ISS*(ISS+1)
      FL=1.0
      IS2=IS(N)
      DO 25 LL=1,LPLUS
      L=LL-1
      LK=LL+LL-1
      L2=LK-1
      J2=L2-IS2
      IQ=LK+IC1
      DO 20 JJ=JR,JS
      IP=LK+ICNT
      MS2=-IS2
      DO 15 I=1,ISS
      VCI=VCC(L2,IS2,J2,0,MS2)*FL
      MT2=MS2
      DO 14 J=I,ISS
      ML2=MS2-MT2
      ML=IABS(ML2/2)
      IF(ML.GT.L) GO TO 13
      VTEMP=VCI*VCC(L2,IS2,J2,ML2,MT2)*SQRT(YXFCT(L+ML,L-ML))
      DTEMP(IP  )=DTEMP(IP  )+D(IQ  )*VTEMP
      DTEMP(IP+1)=DTEMP(IP+1)+D(IQ+1)*VTEMP
   13 CONTINUE
      IP=IP+LPL2
      MT2=MT2+2
   14 CONTINUE
      MS2=MS2+2
   15 CONTINUE
      J2=J2+2
      IQ=IQ+LPL2
   20 CONTINUE
      FL=FL+2.0
   25 CONTINUE
      ICNT=ICNT+INCR(N)
      IC1=IC1+LPL2*ISS
      JR=JS+1
   26 CONTINUE
C
C     CALCULATE ELASTIC CROSS SECTION
C
      IMAX=MAX0(is(1),is(2))+1
      THETA=THETA1
      WRITE(6,9000)
      DO 100 NTH=1,NTHETA
      CALL LGNDR(PLM,IMAX,LPLUS,THETA)
      ICNT=0
      DO 80 N=1,2
      ISS=is(N)+1
      DO 30 I=1,50
      SUMX(I)=0.0
   30 CONTINUE
      DO 40 LL=1,LPLUS
      LK=LL+LL-1
      IP=LK+ICNT
      DO 35 I=1,ISS
      DO 34 J=I,ISS
      ML=(J-I)*LPLUS+LL
      VTEMP=PLM(ML)
      SUM(J,I,1)=SUM(J,I,1)+DTEMP(IP  )*VTEMP
      SUM(J,I,2)=SUM(J,I,2)+DTEMP(IP+1)*VTEMP
      IP=IP+LPL2
   34 CONTINUE
   35 CONTINUE
   40 CONTINUE
      SIG(N)=0.0
      POL(N)=0.0
C
C     CALCULATE COULOMB AMPLITUDE
C
      if(theta.ne.0.0) then
      ARG=THETA*(pi/360.)
      S=SIN(ARG)**2
      TL= LOG(S)
      ARG=ETA(N)*TL
      FACT=ETA(N)/(2.0*FK(N)*S)
      CL=COS(ARG)*FACT
      SL=SIN(ARG)*FACT
      else
      cl=0.0
      sl=0.0
      fact=1.0
      endif
c
      DO 60 I=1,ISS
      DO 59 J=1,ISS
      ML=I-J
      if    (ml.gt.0) then
        PHAS=PHASEF(ML)
        SUM(J,I,1)=SUM(I,J,1)*PHAS
        SUM(J,I,2)=SUM(I,J,2)*PHAS
      elseif(ml.eq.0) then
        SUM(J,I,1)=SUM(J,I,1)-CL
        SUM(J,I,2)=SUM(J,I,2)+SL
      endif
        SIG(N)=SIG(N)+SUM(J,I,1)**2+SUM(J,I,2)**2
   59 CONTINUE
   60 CONTINUE
      CALL POLFCT(POL(N),SUM(1,1,1),SUM(1,1,2),ISS,ISS)
      POL(N)=POL(N)/SIG(N)
      SIG(N)=SIG(N)/float(is(n)+1)
      if(eta(n).eq.0.0) then
        ratio(n)=1.0
      else
        if(theta.eq.0.0) then
          sig(n)  =0.0
          pol(n)  =0.0
          ratio(n)=1.0
        else
          ratio(n)=sig(n)/fact**2
        endif
      endif

      ICNT=ICNT+INCR(N)
      IF(NTH.LE.NX) SIGPLT(NTH,N)=RATIO(N)
   80 CONTINUE
      WRITE(6,9001)(THETA,SIG(N),RATIO(N),POL(N),N=1,2)
   90 CONTINUE
      THETA=THETA+DTHETA
  100 CONTINUE
      WRITE(6,9002)(RSIG(N),N=1,2)
      NTH=MIN0(NTHETA,NX)
      IF(ICON(15).ne.0) then
        WRITE(7,9011)
        WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
        theta=theta1
        do 120 n=1,nth
        WRITE(7,9010)theta,SIGPLT(N,1),SIGPLT(N,2)
        theta=theta+dtheta
  120   continue
      endif
      IF(ICON(6).ne.0) then
        do 140 n=1,2
        WRITE(6,9999)ALPHA,(IDAT(I),I=1,3)
        CALL DWPLOT(NTH,ICON(6),SIGPLT(1,n),DTHETA,THETA1)
  140   continue
      endif
      RETURN
c
 9000 FORMAT(1H0,7H  Theta,'   Elsig,fm**2','   Sigma/Coul '
     1,'       Pol    ',5X,'   Theta','   Elsig,fm**2'
     2,'   Sigma/Coul ','       Pol    ')
 9001 FORMAT(2(0PF8.2,1P2E14.5,0PF11.4,8X))
 9002 FORMAT('0REACSIG',1PE14.5,41X,1PE14.5)
 9010 FORMAT(1PE10.4,5(4X,1PE10.4))
 9011 FORMAT('  Elastic cross section ')
 9999 FORMAT(1H1,15A4,I4,2(1H/,I2.2),I4,2(1H.,I2.2))
      END


c***********************************************************************
      SUBROUTINE SLATR (KT,KMAX,DRF,VB,MINL,FMU,ICODE)
c
c     computes the slater integrals for the microscopic inelastic cases.
c***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VB(800)
C
C     YUKAWA SLATER EXPANSION
C
      KM =KMAX
      R=0.0
      IF(ICODE.EQ.2) GO TO 101
C
C     HANKEL FUNCTION*EXP(+FMU*R)
C
      F1=1.0
      F2=EXP(-FMU*DRF)
      DO 100 M=1,KT
      R=R+DRF
      X=FMU*R
      F1=F1*F2
      AZ=(1.0-F1*F1)/(2.0*X)
      B2=1.0/X
      B1=B2*(1.0+B2)
      FL=-1.0
      DO 50 LL=1,MINL
      B3=FL*B2/X+B1
      B1=B2
      B2=B3
      FL=FL+2.0
   50 CONTINUE
      VB(M+KM )=B2
C
C     BESSEL FUNCTION*EXP(-FMU*R)
C
      MAX=2.0*X+10.0
      MAX=MAX0(MAX,MINL+5)
      A3=0.0
      A2=1.0
      FL=MAX+MAX+3
      DO 70 LL=1,MAX
      N=MAX+1-LL
      FL=FL-2.0
      A1=A3+FL*A2/X
      IF(N.LT.MINL) GO TO 69
      IF(ABS(A1).LT.1.0E+20) GO TO 69
      A1=A1*1.0E-20
      A2=A2*1.0E-20
   69 CONTINUE
      IF(N.EQ.MINL) TEMP=A1
      A3=A2
      A2=A1
   70 CONTINUE
      VB(M    )=TEMP*AZ/A1
  100 CONTINUE
      RETURN
C
C     COULOMB SLATER EXPANSION
C
  101 CONTINUE
      FL=MINL+MINL-1
      DO 200 M=1,KT
      R=R+DRF
      A2=1.0
      DO 105 LL=1,MINL
      A2=A2*R
  105 CONTINUE
      VB(M    )=A2/(R*FL)
      VB(M+KM )=1.0/A2
  200 CONTINUE
      RETURN
      END
CDWK407
      SUBROUTINE RADIN(KT,KMAX,DRF,FMU,VB,UB,UC,SL,OPT,SI,KMT,SK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION UB(400),UC(400),VB(800),SI(800),DG(2)
      DATA XRHO/80./
      FLOAT(III)=DFLOAT(III)
C
      KM =KMAX
      KM2=KMAX*2
      KM3=KMAX*3
      DG(1)=2.0*DRF/3.0
      DG(2)=2.0*DG(1)
      F0=EXP( FMU*DRF)
      MMAX=XRHO/FMU/DRF
      R2=0.0
      SUMA=0.0
      SUMB=0.0
      SUMC=0.0
      SUMD=0.0
      DO 200 M=1,KT
      MK=M+M-1
      R2=R2+DRF
      F2=EXP(-FLOAT(MIN0(MMAX,M)-1)*FMU*DRF)
      SLL=0.0
      MX=1
      R1=0.0
      DO 100 MM=1,KT
      R1=R1+DRF
      MX=3-MX
      IF(IABS(M-MM).GT.MMAX) GO TO 100
      F1=DG(MX)*UB(MM)*UC(MM)*R1**2
      IF(M.GT.MM) GO TO 80
      IF(M.EQ.MM) GO TO 90
      F2=F2/F0
      TEMP=VB(M)*VB(MM+KM3)*F2
      GO TO 95
   80 CONTINUE
      TEMP=VB(MM+KM2)*VB(M+KM )*F2
      F2=F2*F0
      GO TO 95
   90 TEMP=0.5*(VB(M)*VB(M+KM3)+VB(M+KM2)*VB(M+KM ))
   95 CONTINUE
      SLL=SLL+F1*TEMP
  100 CONTINUE
      SLL=SLL*OPT
      R22=R2**2
      IF(M.EQ.KMT) SK=SLL
      SI(MK)=SI(MK)+SLL
      SUMA=SUMA+UB(M)*UC(M)*R22
      SUMB=SUMB+UB(M)*UC(M)*R22**2
      SUMC=SUMC+SLL*R22
      SUMD=SUMD+SLL*R22**2
  200 CONTINUE
      SUMA=SUMA*DRF
      SUMB=SUMB*DRF
      SUMC=SUMC*DRF
      SUMD=SUMD*DRF
      SL=SUMC
      WRITE(6,9100)SUMA,SUMB,SUMC,SUMD
      RETURN
 9100 FORMAT(13H0        J0 =,F11.4,7H   J1 =,F11.4,7H   K0 =,F11.4
     1 ,7H   K1 =,F11.4)
      END


c***********************************************************************
      SUBROUTINE FNLOC5(U,V,W,PNLOC,FK2,FK,ETA,RC,DR,KT)
c
c     Non-locality correction factor for the distorted waves
c
c***********************************************************************
c
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(800),V(800),W(800)
      FACP=PNLOC**2/8.0
      IF(FACP.EQ.0.0) GO TO 410
C
      WRITE(6,9900)PNLOC
      IF(PNLOC.GT.0.0) THEN
        R=0.0
        ETAK=2.0*ETA*FK
      ELSE
        R1=DR*FLOAT(KT+1)
      ENDIF
      CTEMP1=0.0
      CTEMP2=0.0
      T1=0.0
      T2=0.0
C
      DO 400 M=1,KT
      IF(PNLOC.GT.0.0) THEN
C       PNLOC POSITIVE GIVES USUAL NON-LOCAL FACTOR WITH GAUSSIAN SHAPE
        MK=M+M-1
        R=R+DR
        IF(R.GT.RC) THEN
          VCOUL=ETAK/R
        ELSE
          VCOUL=0.5*ETAK*(3.0-(R/RC)**2)/RC
        ENDIF
        CTEMP1=-FACP*U(MK  )+FACP*(FK2-VCOUL)
        CTEMP2=-FACP*U(MK+1)
        T1=U(MK  )
        T2=U(MK+1)
      ELSE
C       PNLOC NEGATIVE GIVES DIRAC TYPE DARWIN FACTOR BASED ON SPIN ORBI
        MK=KT+KT-M-M+1
        R2=R1
        R1=R1-DR
        CTEMP1=CTEMP1+PNLOC*(R2*T1+R1*V(MK  ))*(DR*0.125)
        CTEMP2=CTEMP2+PNLOC*(R2*T2+R1*V(MK+1))*(DR*0.125)
        T1=V(MK  )
        T2=V(MK+1)
      ENDIF
C
      FACT=EXP(CTEMP1)
      UT1=FACT*COS(CTEMP2)
      UT2=FACT*SIN(CTEMP2)
      UF1   =UT1*W(MK  )-UT2*W(MK+1)
      UF2   =UT1*W(MK+1)+UT2*W(MK  )
      W(MK  )=UF1
      W(MK+1)=UF2
c      WRITE(20,7777)R1,T1,T2,UT1,UT2
c 7777 FORMAT(' R =',F6.2,'  POT =',1P2E12.4,'  FACTOR =',1P2E12.4)
  400 CONTINUE
  410 CONTINUE
      RETURN
 9900 FORMAT('0Non-locality correction is applied, parameter =',f8.3)
      END


