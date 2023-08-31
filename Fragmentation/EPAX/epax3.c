#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* EPAX V31 function by Klaus Sümmerer, March 2012 */
/* version of 15.03.2012, as in PRC with erratum   */
/* translated to C by Helmut Weick 09.01.2013      */
/* compile with: cc epax_v31-new.c -lm -oepax_v31  */

/* prototyping */
double epaxv3(double Ap, double Zp, double At, double Zt, double A, double Z);

main(int argc, char *argv[])
{
  double Ap, Zp, At, Zt, A, Z, result;

  printf("\n"); 
  printf("*****************************************\n");
  printf("*                                       *\n");
  printf("*            EPAX Version 3.1           *\n");
  printf("*                                       *\n");
  printf("*      An Empirical Parametrization     *\n");
  printf("*      of Projectile-Fragmentation      *\n");
  printf("*            Cross Sections             *\n");
  printf("*      Klaus Suemmerer, 15.03.2012      *\n");
  printf("*****************************************\n");
  printf("\n");

/* start program as epax_v31-new 58 28 9 4 48 28 */
  if(argc>5) {
    sscanf(argv[1],"%lf",&Ap);
    sscanf(argv[2],"%lf",&Zp);
    sscanf(argv[3],"%lf",&At);
    sscanf(argv[4],"%lf",&Zt);
    sscanf(argv[5],"%lf",&A);
    sscanf(argv[6],"%lf",&Z);

    result = epaxv3(Ap, Zp, At, Zt, A, Z);
    if(result<0){
      printf("error: Input parameters incorrect!");
    }

    printf("\n");
    printf("*******************************************\n");
    printf("Fragmentation cross section:\n");
    printf(" %f %f -> %f %f -> %f %f\n",Ap,Zp,At,Zt,A,Z); 
    printf(" sigma = %e b \n",result);
    printf("*******************************************\n");
  }
  else {
    printf("Error, no valid input !\n");
    printf("usage:   epax_v31 <A_p> <Z_p> <A_t> <Z_t> <A_p> <Z_p> \n");
    printf("example: ./epax_v31 58 28 9 4 48 28 \n");
    printf("result should be sigma = 1.4075E-14 b \n\n");
  }
}




/***************************************************************************
*
* Returns cross section (in barn) for fragment (A,Z) in projectile 
* fragmentation of projectile (Ap,Zp) on target (At,Zt)
* Version: 3.1 15-March-2012 K.Suemmerer 
*
****************************************************************************/

double epaxv3(double Ap, double Zp, double At, double Zt, double A, double Z) 
{
  double At13, Ap13;
  double R,r_0,p,f_pt,f_mod_y,zprob,zbeta,dzbeta_p,delta,dq;
  double f_mod,slope,z0,a2,expo_fac,offset,dzprob,u_p,u_n;
  double expo,fract,yield_a,zbeta_p,fwhm;
  double cd[7],cp[3],cr[8],cn[5],cq[3],un[2],up[2];
  double cs[2],cl[3],cy[3],cb[3],norm[2],result;

  result = -999.0 ;

  Ap13 = pow(Ap,0.33333);
  At13 = pow(At,0.33333);

  /************ Constants of EPAX V 3.1 *******************************/

  cs[1]  = 0.270 ;        /* scaling factor (barn) */
  cs[2]  = 1.80 ;

  cd[1] = -1.0870E+00 ;   /* centroid rel. to beta-stability */
  cd[2] = +3.0470E-02 ;
  cd[3] = +2.1353E-04 ; 
  cd[4] = +7.1350E+01 ;
  cd[5] = -25.0 ;         /* correction close to proj.: centroid */
  cd[6] = 0.80 ; 

  cp[1] = -1.731 ;        /* mass yield slope P */  
  cp[2] = -0.01399 ;

  cr[1] = 2.78 ;          /* width parameter R */
  cr[2] = -0.015 ;         
  cr[3] = 0.320E-04 ;
  cr[4] = 0.0412 ;
  cr[5] = 0.124 ;
  cr[6] = 30.*sqrt(Ap) ;   /* correction close to proj.: width R */
  cr[7] = 0.85 ;

  un[1] = 1.65 ;           /* slope par. n-rich ride of Z distr. */
  up[1] = 2.10 ;           /* slope par. p-rich ride of Z distr. */  

  cn[1] = 0.000 ;          /* memory effect n-rich projectiles */
  cn[2] = 0.400 ;   
  cn[3] = 0.000 ;  
  cn[4] = 0.600 ;  

  cq[1] = -10.25 ;         /* memory effect p-rich projectiles */
  cq[2] = +10.25 ;  

  cl[1] = 1.20 ;           /* lin. slope on p-rich side */
  cl[2] = 0.647 ;   

  cy[1]  = 0.75 ;          /* correction of mass yield */
  cy[2]  = 0.10 ;           

  cb[1]  = 2.3e-03 ;       /* "brute force" correction for very n-rich */
  cb[2]  = 2.4 ;

  /*****************************************************************************

  /* calculate mass yield: */

  /* slope parameter */
  p = exp(cp[1]+cp[2]*Ap) ;
  f_pt = cs[1]*p*(Ap13 +At13-cs[2]) ;
  yield_a = f_pt * exp(-p * (Ap - A)) ;

  /* steeper slope near projectile */
  f_mod_y=1.0 ;
  if (A/Ap > cy[1]) {
    f_mod_y=exp(cy[2]*Ap*pow(A/Ap-cy[1],2)) ;
  }
  yield_a = yield_a * f_mod_y ;

  /* ---- calculate zprob as for V2.1 -----*/
  zbeta = A/(1.98+0.0155*pow(A,2./3.)) ;
  zbeta_p = Ap/(1.98+0.0155*pow(Ap,2./3.)) ;
  dzbeta_p = Zp - zbeta_p ;
  if(A > cd[4]) {
     delta = cd[1]+cd[2]*A ;
  }
  else
  {
     delta = cd[3]*A*A ; 
  } 

  f_mod = 1.0 ;
  if(A/Ap > cd[6]) { 
    f_mod = cd[5] * pow(A/Ap-cd[6],2) +1.0 ;
  }
  delta = delta * f_mod ;
  zprob = zbeta + delta ;

  if((Zp-zbeta_p) > 0) {  /* ------------- p-rich */
     dq = exp(cq[1]+cq[2]*(A/Ap)) ;
  }
  else                /* ------------- n-rich */
  {
     dq = A/Ap*(cn[1]+A/Ap*(cn[2]+A/Ap*(cn[3]+A/Ap*cn[4]))) ;
  }

  dzprob = dq*(Zp-zbeta_p) ; 

  /* small corr. since Xe-129 and Pb-208 are not on Z_beta line */
  zprob = zprob +dzprob +0.0020*A ;

  /* calculate width parameter */

  if ((Zp-zbeta_p)<0) {           /* n-rich */
     r_0 = cr[1]* exp(cr[4]*dzbeta_p) ;
  }
  else  
  {                               /* p-rich */
     r_0 = cr[1]* exp(cr[5]*dzbeta_p) ;
  }   
  R = r_0 * exp(cr[2]*A + cr[3]*A*A) ;
  f_mod = 1.0 ;
  if (A/Ap > cr[7]) {
    f_mod = exp(cr[6] * pow(A/Ap-cr[7],3.)) ;
  }
  R = R * f_mod ;

  u_p = up[1] ;
  u_n = un[1] ;
  if((zprob-Z) > 0) {
  /*     neutron-rich */
    expo = -R * pow(fabs(zprob-Z),u_n) ;
    fract = exp(expo) * sqrt(R/3.14159) ;
  }
  else {
  /* old V2.1 parameterization */
    slope = cl[1] + cl[2] * pow(A/2.,0.3) ;
    z0 = zprob+( (cl[1]+cl[2] * pow(A/2.,0.3) ) *log(10.)/(2.*R)) ;
    expo = -R * pow(fabs(zprob-Z),u_p) ;
    fract   =  exp(expo) * sqrt(R/3.14159) ;
    if( Z > z0 ) {
      expo  = -R * pow(fabs(zprob-z0),u_p) ;
      a2    =  exp(expo) * sqrt(R/3.14159) ;
      /* fract =  a2 * exp(slope*(z0-Z)) */
      fract = a2/pow( pow(10,slope),(Z-z0) ) ;
    }
  }

  fwhm = pow(0.693/R,1.0/u_n) + pow(0.693/R,1.0/u_p) ;
  norm[1] = 1.0 ;

  /* "brute force" scaling factor for n-rich projectiles only */

  if((zbeta_p-Zp)>0) { 
    expo_fac = -cb[1] * fabs(Zp-zbeta_p) ;
    if ((zbeta-Z)>(cb[2]+Zp-zbeta_p)) { 
      offset = Zp - zbeta_p + cb[2] ;
      norm[1]= pow(10.0,expo_fac * pow(zbeta-Z+offset,3)) ;
    }
    else
    {
      norm[1]=1.0 ;
    }
  }

  fwhm = pow(0.693/R,1.0/u_n) + pow(0.693/R,1.0/u_p) ; 
  result =  norm[1] * fract * yield_a ; 

  return(result);
}
