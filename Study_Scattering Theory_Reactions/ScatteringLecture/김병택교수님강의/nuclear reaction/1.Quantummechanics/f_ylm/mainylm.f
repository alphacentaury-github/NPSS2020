      PROGRAM MAIN 
      DIMENSION  PKK(30,30),PK(30,30)
	DIMENSION  YL(500,30,30)
      DOUBLE PRECISION A 

      open(unit=7, file='ylmcal.dat',status='old')
      open(unit=8, file='ylmcal.out',status='unknown')

      READ(7,10) AMIN,ASTEP
	READ(7,11) MUANG,NMAX,LLL,MMM
10    FORMAT(10F7.3)
11    FORMAT(10I5)  

	PI=4.0*ATAN(1.)
	DO 100 N=1,NMAX
	AA=AMIN+ASTEP*FLOAT(N-1)
	IF (MUANG.EQ.1) AA=COS(AA*PI/180.0)
	A=AA
	CALL YLCAL(A,LLL,MMM,PKK)
	DO 90 L=1,LLL
	DO 90 M=1,MMM
	YL(N,L,M)=PKK(L,M)
90    CONTINUE
100   CONTINUE

      LMAX=LLL
      DO 200 L=1,LMAX
	LREAL=L-1
	WRITE(8,20) LREAL
	MMAX=L
	DO 190 N=1,NMAX
	X=AMIN+ASTEP*FLOAT(N-1)
	WRITE(8,21) X, (YL(N,L,M),M=1,MMAX)
20    FORMAT(//1H ,'ABS OF SPHERICAL HARMONICS FOR L=',I4)
21    FORMAT(F9.4,10E13.5)
190   CONTINUE
200   CONTINUE

      STOP
      END
                                                                        
