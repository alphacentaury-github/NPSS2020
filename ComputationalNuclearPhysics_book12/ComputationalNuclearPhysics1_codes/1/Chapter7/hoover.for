C ****************************************************************
      program HOOVER
C
C     Written by David P. Murdock
C     Last revision: 1/15/90
C
C     Simple version intended for use with simple FOLDER.  
C     Calculates elastic scattering observables for proton scattering
C     from a spin-zero nucleus by solving the Dirac equation.
C
C     ==> When compiling with MS Fortran, it may be necessary to use
C     the /Od switch.  
C *******************************************************************
C
      implicit double precision (a-h,o-z)
      character*72 adum
      common/pwave/lmax,f(0:100,2),g(0:100,2),sigma(0:100),
     *        p(0:100),p1(0:100),uval(0:200,2,2,2)
      common/kin/elab,ecm,fk,fk2,eta,rhonuc,z1,z2,am1,am2,
     *        fm1,fm2,fmu,sefac
      common/math/hbarc,pi,amu
      common/output/icon(16),theta0,dtheta,ntheta,dsig(0:300),
     *        ay(0:300),q(0:300)
      common/rnuc/dr,rmax,n,rcou,rval(500)
      common/potls/ipotop,param(4,2,3),w(500,2),wso(500,2),
     *        wcou(500)
C
C
C
      OPEN(UNIT=8,FILE='HOOVER.INP',STATUS='OLD')
      OPEN(UNIT=14,FILE='HOOVER.OBS',STATUS='UNKNOWN')
      OPEN(UNIT=23,FILE='HOOVER.LOG',STATUS='UNKNOWN')
C
      read(8,12) adum
12    format(a72)
      read(8,*) theta0,thmax,dtheta
      read(8,12) adum
      read(8,*) almax,dr,rmax
      read(8,12) adum
      read(8,*) elab,am1,z1,am2,z2,r0cou
      read(8,12) adum
      read(8,*) ipotop
C
      write(6,15)
      write(23,15)
15    format(' *****HOOVER! (1990)*****',//,
     *        ' Elastic scattering for spin--1/2 particles')
C
      write(6,31)
      write(23,31)
      write(6,32) theta0,thmax,dtheta
      write(23,32) theta0,thmax,dtheta
      write(6,33) almax,dr,rmax
      write(23,33) almax,dr,rmax
      write(6,34) elab,am1,z1,am2,z2
      write(23,34) elab,am1,z1,am2,z2
      write(6,35) ipotop
      write(23,35) ipotop
31    format(/' Program Input:')
32    format(5x,'theta0=',f7.3,'  thmax=',f8.4,'   dtheta=',f6.3)
33    format(5x,'alMax=',f7.2,'  dr=',f8.4,'   rmax=',f7.3)
34    format(5x,'Elab=',f7.3,'   AM1=',f8.3,'   Z1=',f7.3,
     *'   AM2=',f8.3,'   Z2=',f7.3)
35    format(5x,'Potential Option=',i2)
C
      lmax=almax+.001
      n=rmax/dr+.001
      ntheta=(thmax-theta0)/dtheta+1
      rcou=r0cou*am2**.333333
C
      write(6,40) LMAX
40    format(' ',I3,' partial waves')
C
      call driver
C
      stop
      end
C ******************************************************************
C
      BLOCK DATA INIT
C
      implicit double precision (a-h,o-z)
      common/math/hbarc,pi,amu
C
      data hbarc,amu,pi/197.329,931.501,3.1415927/
C
      end
C ******************************************************************
C
      subroutine DRIVER
      implicit double precision (a-h,o-z)
C
      call setup
      call getpot
      call coulom
      call intup
      call xsect
C
      return
      end
C ******************************************************************
C
C
      subroutine SETUP
      implicit double precision (a-h,o-z)
      common/pwave/lmax,f(0:100,2),g(0:100,2),sigma(0:100),
     *        p(0:100),p1(0:100),uval(0:200,2,2,2)
      common/kin/elab,ecm,fk,fk2,eta,rhonuc,z1,z2,am1,am2,
     *        fm1,fm2,fmu,sefac
      common/math/hbarc,pi,amu
      common/rnuc/dr,rmax,n,rcou,rval(500)
      common/potls/ipotop,param(4,2,3),w(500,2),wso(500,2),
     *        wcou(500)
Cc
      do 10 i=1,n
10    rval(i)=i*dr
C
      fm1=am1*amu
      fm2=am2*amu
      fmt=fm1+fm2
      t1=sqrt(2.*elab*fm2+fmt**2)
      fk2=(elab**2+2.*elab*fm1)*(fm2/t1)**2/hbarc**2
      fk=sqrt(fk2)
      fkmev=fk*hbarc
      rhonuc=fk*(n-1)*dr
C
      fmu=fm1
      ecm=sqrt((fk*hbarc)**2+fm1**2)
C
      sefac=2.*fmu/hbarc**2
      eta=z1*z2*fmu/(137.036*fk*hbarc)
C
      write(6,90) ecm, fkmev
      write(23,90) ecm, fkmev
90    format(/' Projectile total energy and momentum, CM:',/,
     *' E cm, MeV',
     *       f12.4, '       p cm, MeV= ',f12.4)
C
      return
      end
C *******************************************************************
C
C
      subroutine COULOM
C
      implicit double precision (a-h,o-z)
      common/kin/elab,ecm,fk,fk2,eta,rhonuc,z1,z2,am1,am2,
     *        fm1,fm2,fmu,sefac
      common/output/icon(16),theta0,dtheta,ntheta,dsig(0:300),
     *        ay(0:300),q(0:300)
      common/rnuc/dr,rmax,n,rcou,rval(500)
C
      call getsig
      call couwf(1,rhonuc,icon(13))
      call couwf(2,rhonuc+fk*dr,icon(13))
      return
      end
C ******************************************************************
C
C
      subroutine GETSIG
C
      implicit double precision (a-h,o-z)
      common/pwave/lmax,f(0:100,2),g(0:100,2),sigma(0:100),
     *        p(0:100),p1(0:100),uval(0:200,2,2,2)
      common/kin/elab,ecm,fk,fk2,eta,rhonuc,z1,z2,am1,am2,
     *        fm1,fm2,fmu,sefac
C
      t1=-eta+(eta/2.)*dlog(eta**2+16.)+3.5*atan(eta/4.)
      t2=atan(eta)+atan(eta/2.)+atan(eta/3.)
      t3=eta/(12.*(eta**2+16.))
      t4=1.+(eta**2-48.)/(30.*(eta**2+16.)**2)
      t5=(eta**4-160.*eta**2+1280.)/(105.*(eta**2+16.)**4)
      sigma(0)=t1-t2-t3*(t4+t5)
C
      do 10 i=1,lmax
10    sigma(i)=sigma(i-1)+atan(eta/i)
      return
      end
C ******************************************************************
C
C
      subroutine COUWF(ipt,rho,iter)
C
      implicit double precision (a-h,o-z)
      common/pwave/lmax,f(0:100,2),g(0:100,2),sigma(0:100),
     *        p(0:100),p1(0:100),uval(0:200,2,2,2)
      common/kin/elab,ecm,fk,fk2,eta,rhonuc,z1,z2,am1,am2,
     *        fm1,fm2,fmu,sefac
C
      call couval(ipt,rho)
      do 50 i=1,lmax
         t1=(2.*i+1.)*(eta+i*(i+1.)/rho)*g(i,ipt)
         t2=(i+1.)*sqrt(i**2+eta**2)*g(i-1,ipt)
         t3=i*sqrt((i+1.)**2+eta**2)
50    g(i+1,ipt)=(t1-t2)/t3
C
      if(iter.eq.1) go to 58
      do 56 i=1,lmax
         t1=(2.*i+1.)*(eta+i*(i+1.)/rho)*f(i,ipt)
         t2=(i+1.)*sqrt(i**2+eta**2)*f(i-1,ipt)
         t3=i*sqrt((i+1.)**2+eta**2)
56    f(i+1,ipt)=(t1-t2)/t3
      return
C
58    f(lmax+11,ipt)=0.
      f(lmax+10,ipt)=0.001
      i=lmax+10
60    t1=(2.*i+1.)*(eta+i*(i+1.)/rho)*f(i,ipt)
      t2=i*sqrt((i+1.)**2+eta**2)*f(i+1,ipt)
      t3=(i+1.)*sqrt(i**2+eta**2)
      f(i-1,ipt)=(t1-t2)/t3
      if(abs(f(i-1,ipt)).gt.1.0e+10) then
         do 65 j=i-1,lmax+11
65       f(j,ipt)=f(j,ipt)*1.0e-10
      endif
      i=i-1
      if(i.eq.0) go to 70
      go to 60
C
70    alpha=(f(0,ipt)*g(1,ipt)-f(1,ipt)*g(0,ipt))
     1        *sqrt(1.+eta**2)
      do 80 i=0,lmax+1
80    f(i,ipt)=f(i,ipt)/alpha
C
      return
      end
C ******************************************************************
C
C
      subroutine COUVAL(ipt,rho)
C
      implicit double precision (a-h,o-z)
      complex fg0,fg1,a1,a2,b1,b2,term0,term1,x
C
      common/pwave/lmax,f(0:100,2),g(0:100,2),sigma(0:100),
     *        p(0:100),p1(0:100),uval(0:200,2,2,2)
      common/kin/elab,ecm,fk,fk2,eta,rhonuc,z1,z2,am1,am2,
     *        fm1,fm2,fmu,sefac
      common/math/hbarc,pi,amu
C
      fg0=cmplx(1.,0.)
      fg1=cmplx(1.,0.)
      x=2.*cmplx(0.,1.)*rho
      term0=cmplx(1.,0.)
      term1=cmplx(1.,0.)
      a1=cmplx(0.d0,eta)
      a2=cmplx(1.d0,eta)
      b1=cmplx(-1.d0,eta)
      b2=cmplx(2.d0,eta)
      do 40 i=1,20
      term0=term0*a1*a2/(i*x)
      term1=term1*b1*b2/(i*x)
      fg0=fg0+term0
      fg1=fg1+term1
      a1=a1+1.
      a2=a2+1.
      b1=b1+1.
      b2=b2+1.
40    continue
C
50    thet0=rho-eta*dlog(2.*rho)+sigma(0)
      thet1=rho-eta*dlog(2.*rho)+sigma(1)-pi/2.
      f(0,ipt)=aimag(fg0)*cos(thet0)+real(fg0)*sin(thet0)
      f(1,ipt)=aimag(fg1)*cos(thet1)+real(fg1)*sin(thet1)
      g(0,ipt)=real(fg0)*cos(thet0)-aimag(fg0)*sin(thet0)
      g(1,ipt)=real(fg1)*cos(thet1)-aimag(fg1)*sin(thet1)
C
      return
      end
C ******************************************************************
C
C
      subroutine INTUP
C
      implicit double precision (a-h,o-z)
      complex u,up,um,s,sm,sp,temp
C
      common/pwave/lmax,f(0:100,2),g(0:100,2),sigma(0:100),
     *        p(0:100),p1(0:100),uval(0:200,2,2,2)
      common/rnuc/dr,rmax,n,rcou,rval(500)
      common/potls/ipotop,param(4,2,3),w(500,2),wso(500,2),
     *        wcou(500)
C
      dimension lfac(2)
C
      snum=dr*dr/12.
C
      do 500 l=0,lmax
         lfac(1)=l
         lfac(2)=-l-1
C
         do 495 isp=1,2
             um=cmplx(0.,0.)
             u=cmplx(1.e-10,0.)
             temp=cmplx(w(1,1),w(1,2))-l*(l+1.)/dr**2
             s=temp+lfac(isp)*cmplx(wso(1,1),wso(1,2))
             i=2
100          temp=cmplx(w(i,1),w(i,2))-l*(l+1.)/rval(i)**2
             sp=temp+lfac(isp)*cmplx(wso(i,1),wso(i,2))
             up=(2.*(1.-5.*snum*s)*u-(1.+snum*sm)*um)/(1.+snum*sp)
             if(i.ge.n) go to 200
             if(abs(up).ge.1.e+10) then
                  u=u*1.d-10
                  up=up*1.d-10
             endif
             um=u
             u=up
             sm=s
             s=sp
             i=i+1
             go to 100
200          uval(l,isp,1,1)=real(u)
             uval(l,isp,1,2)=aimag(u)
             uval(l,isp,2,1)=real(up)
             uval(l,isp,2,2)=aimag(up)
C
495      continue
500   continue
      return
      end
C *******************************************************************
C
C
      subroutine XSECT
C
      implicit double precision (a-h,o-z)
      complex c(0:100,2),ff(0:100,2),a,b,ci,t1(2),t2(2),
     *        umatch(2,2),e2isig,a1,b1
C
      common/pwave/lmax,f(0:100,2),g(0:100,2),sigma(0:100),
     *        p(0:100),p1(0:100),uval(0:200,2,2,2)
      common/kin/elab,ecm,fk,fk2,eta,rhonuc,z1,z2,am1,am2,
     *        fm1,fm2,fmu,sefac
      common/math/hbarc,pi,amu
      common/output/icon(16),theta0,dtheta,ntheta,dsig(0:300),
     *        ay(0:300),q(0:300)
C
      ci=cmplx(0.,1.)
C
      do 100 l=0,lmax
          e2isig=exp(2.*ci*sigma(l))
C
          do 35 isp=1,2
             do 25 ipt=1,2
25           umatch(isp,ipt)=cmplx(uval(l,isp,ipt,1),
     1           uval(l,isp,ipt,2))
             t1(isp)=umatch(isp,1)*f(l,2)-umatch(isp,2)*f(l,1)
             t2(isp)=umatch(isp,2)*g(l,1)-umatch(isp,1)*g(l,2)
     1        +ci*(umatch(isp,2)*f(l,1)-umatch(isp,1)*f(l,2))
             c(l,isp)=t1(isp)/t2(isp)
             ff(l,isp)=c(l,isp)*e2isig/fk
35        continue
C
C         if(icon(12).gt.0) call wavout
C
100   continue
C
      do 500 i=0,ntheta
          theta=theta0+i*dtheta
          call leg(theta)
          sin2=sin(theta*pi/360.)**2
          a1=-ci*eta*dlog(sin2)+2.*ci*sigma(0)
          b1=cmplx(-1.,0.)*eta/(2.d0*fk*sin2)
          a=cmplx(0.,0.)
          b=cmplx(0.,0.)
          do 200 l=0,lmax
             a=a+((l+1.)*ff(l,1)+l*ff(l,2))*p(l)
200       b=b+(ff(l,1)-ff(l,2))*p1(l)
          a=b1*exp(a1)+a
          b=b/ci
          dsigg=abs(a)**2+abs(b)**2
          dsig(i)=10.*dsigg
          ay(i)=2.*real(conjg(a)*b)/dsigg
          q(i)=-2.*aimag(conjg(a)*b)/dsigg
500    continue
C
C
      write(14,549)
549   format(1x,' theta,  cross sec (mb/sr), analyzing power, ',
     *'  spin rotation')
      do 550 i=1,ntheta
         theta=theta0+i*dtheta
550      write(14,553) theta, dsig(i),ay(i),q(i)
553   format(1x,f7.3,3(1pe17.5))
      write(6,570)
      write(23,570)
570   format(/' Calculation complete.',/,
     *' Observables written to file HOOVER.OBS')
      return
      end
C *****************************************************************
C
C
      subroutine LEG(theta)
C
      implicit double precision (a-h,o-z)
      common/pwave/lmax,f(0:100,2),g(0:100,2),sigma(0:100),
     *        p(0:100),p1(0:100),uval(0:200,2,2,2)
      common/math/hbarc,pi,amu
C
      x=cos(theta*pi/180.)
      p(0)=1.
      p(1)=x
      p1(0)=0.
      p1(1)=-sqrt(1.-x**2)
C
      do 100 l=1,lmax-1
          p(l+1)=((2.*l+1.)*x*p(l)-l*p(l-1))/(l+1.)
100   p1(l+1)=((2.*l+1.)*x*p1(l)-(l+1.)*p1(l-1))/l
C
      return
      end
C *****************************************************************
C
C
      subroutine GETPOT
C
      implicit double precision (a-h,o-z)
      complex vs,vv,vsp,vvp,vspp,vvpp,b,bp,bpp,
     *        vnuc,vcoul,vdar,vc,vso
C
      common/kin/elab,ecm,fk,fk2,eta,rhonuc,z1,z2,am1,am2,
     *        fm1,fm2,fmu,sefac
      common/math/hbarc,pi,amu
      common/rnuc/dr,rmax,n,rcou,rval(500)
      common/potls/ipotop,param(4,2,3),w(500,2),wso(500,2),
     *        wcou(500)
C
      dimension pot1(500,2),pot2(500,2),pot4(500,2)
C
      cfac=z1*z2*hbarc/137.06
C
      if(ipotop.eq.1) go to 1100
      if(ipotop.eq.0) go to 1400
C
C    READ IN DIRAC POTENTIALS
1100  OPEN(UNIT=12,FILE='FOLDER.POT',STATUS='OLD')
      do 1110 i=1,n
          read(12,*) rdum, pot1(i,1),pot1(i,2),pot2(i,1),pot2(i,2)
          if(rval(i).lt.rcou) pot4(i,1)=.5*cfac*
     1    (3.-(rval(i)/rcou)**2)/rcou
          if(rval(i).ge.rcou) pot4(i,1)=cfac/rval(i)
1110  continue
      go to 1800
C
C      ===> Wood-Saxon parameters for S and V potentials
C
1400  read(8,*) param(1,1,1),param(1,1,2),param(1,1,3)
      read(8,*) param(1,2,1),param(1,2,2),param(1,2,3)
      read(8,*) param(2,1,1),param(2,1,2),param(2,1,3)
      read(8,*) param(2,2,1),param(2,2,2),param(2,2,3)
C
      do 1420 iri=1,2
          factr=exp(dr/param(1,iri,3))
          rexp=exp(-param(1,iri,2)/param(1,iri,3))
          do 1415 i=1,n
             rexp=rexp*factr
1415      pot1(i,iri)=param(1,iri,1)/(1.+rexp)
1420  continue
      do 1440 iri=1,2
          factr=exp(dr/param(2,iri,3))
          rexp=exp(-param(2,iri,2)/param(2,iri,3))
          do 1435 i=1,n
             rexp=rexp*factr
1435      pot2(i,iri)=param(2,iri,1)/(1.+rexp)
1440  continue
      do 1460 i=1,n
          if(rval(i).lt.rcou) pot4(i,1)=.5*cfac*
     1       (3.-(rval(i)/rcou)**2)/rcou
          if(rval(i).ge.rcou) pot4(i,1)=cfac/rval(i)
1460  continue
C
1800  ddr=2.*dr
      drdr=dr*dr
      eplusm=ecm+fmu
C
      do 1850 i=1,n
          im=i-1
          i0=i
          ip=i+1
          if(i.eq.1) then
               im=1
               i0=2
               ip=3
          endif
          if(i.eq.n) then
               im=n-2
               i0=n-1
               ip=n
               endif
          vs=cmplx(pot1(i,1),pot1(i,2))
          vv=cmplx(pot2(i,1),pot2(i,2))
          vsp=cmplx((pot1(ip,1)-pot1(im,1))/ddr,
     1        (pot1(ip,2)-pot1(im,2))/ddr)
          vvp=cmplx((pot2(ip,1)-pot2(im,1))/ddr,
     1        (pot2(ip,2)-pot2(im,2))/ddr)
          vspp=cmplx((pot1(ip,1)-2.*pot1(i0,1)+pot1(im,1))/drdr,
     1        (pot1(ip,2)-2.*pot1(i0,2)+pot1(im,2))/drdr)
          vvpp=cmplx((pot2(ip,1)-2.*pot2(i0,1)+pot2(im,1))/drdr,
     1        (pot2(ip,2)-2.*pot2(i0,2)+pot2(im,2))/drdr)
C
          b=1.+(vs-vv)/eplusm
          bp=(vsp-vvp)/eplusm
          bpp=(vspp-vvpp)/eplusm
          vdar=(hbarc**2/(2.*fmu))*
     1        (-bp/(rval(i)*b)-bpp/(2.*b)+.75*(bp/b)**2)
          vnuc=vs+vv*ecm/fmu+(vs**2-vv**2)/(2.*fmu)
          vcoul=pot4(i,1)-pot4(i,1)*vv/fmu
          vc=vnuc+vdar+vcoul
          vso=hbarc*(-hbarc*bp/b)/(2.*fmu*rval(i))
          w(i,1)=fk2-sefac*real(vc)
          w(i,2)=-sefac*aimag(vc)
          wso(i,1)=-sefac*real(vso)
          wso(i,2)=-sefac*aimag(vso)
1850  continue
C
      return
      end
C ******************************************************************
C ******************************************************************
