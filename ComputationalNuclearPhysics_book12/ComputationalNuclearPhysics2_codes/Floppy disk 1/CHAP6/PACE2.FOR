c    PACE2.FOR
c
c     VAX users please note: Compile with /G_FLOATING option.
c     Thus:  FOR/G_FLOATING pace2
c            LINK           pace2
c            RUN            pace2
c            .
c            .
c            . . . . . . input deck
c
c
c     this is a modified version of  julian - the hillman-eyal evapor.
c     monte-carlo code coupling angular momentum. this version is called
c
c      p a c e 2  - projection angular-momentum coupled evaporation.
c
c     many new features have been added and a major part of the code
c     rewritten. it has been modified to have convenient, simple
c     input and to run much more efficiently. it can now handle
c     high excitation energies and spins in short running times, thus
c     enabling its use for high energy heavy ion reactions.
c
c 1)   the input has been shortened. there are default options
c      for almost everything the program needs.
c
c 2)   masses are read off the wapstra-1983 mass table
c      if table does not contain mass, liquid drop mass (lysekil-1966)
c      is substituted.
c
c 3)   the light emitted particles all have default optical model
c      parameters.
c
c 4)   the storage of tl's has been completely rewritten -
c     the same tl's are used for all compound nuclei.the shift in
c     the coulomb barrier during deexcitation has been accounted for
c     by calculating the tl's at an effective energy determined by
c     the shift.
c
c 5)  a default level density is taken from gilbert and cameron-
c     the little-a parameter can alternatively be taken as mass/const.
c     the program decides during deexcitation what level densities
c     and masses it needs, and this does not have to be decided in
c     the preparation of the input
c    ** the yrast line of gilbert and cameron can be replaced by the
c     rotating diffused suface yrast line of sierk. this is close to
c     the rotating liquid drop line of cohen, plasil and swiatecki but
c     improved. outside the range of sierk's values this defaults to the
c     r.l.d.m. model yrast line.
c
c 6)  fission is considered as a decay  mode.
c
c     the default fission barrier is the a. j. sierk diffuse surface
c     nucleus fission barrier. outside the range of sierk's values the
c     default value is the r.l.d.m. results, courtesy of f. plasil, ornl
c     you are welcome however,
c     to input your own barrier.  the saddle point level density is
c     determined by a little-a parameter times aratio.
c     aratio is actually  'a-fiss' / 'a-gamma'
c     ( see below )
c
c 7)  an internal discriminator prevents wasting time calculating
c     probabilities for highly unlikely events..
c
c 8)  the projection of the angular momentum is also followed
c     *******************************************************
c
c     the purpose of this is to obtain correct angular distributions
c     of the residual nuclei and the emitted particles. fractional
c     (i.e. spin 1/2 ) angular momentum is neglected in this part
c     ( and in this part only ) of the calculations since it should
c     not have a significant effect on these distributions.
c
c     for compound nucleus reactions ( mdir=0, see following) angular
c     distibutions are determined as funtion of angle around recoil axis
c     for deep inelastic products (mdir=1, see following) angular distr.
c     are determined as function of angle from axis perpendicular to
c     reaction plane. note that this is only partial information.
c     should more ang. dist. information be required, note that at the
c     end the velocity vectors of all final products are specified -
c     vz - perpendicular to plane. vy - recoil direction. vx perpendicul
c     to both these.
c
c 8a) if angular distribution is requested (idist=1) with partial trace-
c     -back mode (itrac=0, see below) the absolute value of m is printed
c     out as function of (ex,spin). this is useful for fission angular
c     distribution determination.
c
c 9)   a routine to determine the input fusion (compound nucleus)
c      cross section via the bass method (p.r.l. 1977) is included.
c **** important modification 2/24/84
c      cross section is now truncated so that calculations are performed
c      only for those partial waves which have a finite fission barrier.
c      of course, higher partial waves fission immediately and a
c      calculation is superfluous
c
c 10)  a non-statistical yrast cascade gamma decay chain has been
c      artificially incorporated to simulate gamma multiplicity and
c      energy results. when the level density table cannot accomodate
c      further decay due to spin inhibition (delta-j .gt. 2 ) a cascade
c      of gamma rays with  delta-j = 2 and  egamma = ex *(4j-2)/ j(j+1)
c      brings the inhibited cascade to it final ground state.
c
c      note - subroutine chnprb contains a low probability cutoff
c             of the form
c
c     if(iprob.eq.iprobo.and.icut.eq.0.and.tl(1).ge.1.e-10)return
c
c             if you are interested in individual xn xsections this
c             statement probably should be commented out. then you
c             also should worry about exact details of the level
c             density table at low energies. print it out (subr.
c             gclvd  at end - look at tape4)
c
c 11)  a cascade dump and trace feature has been added to  the program
c
c      every emitted particle is dumped on disk (file 2) together
c      with all relevant emission parameters which enable a
c      complete traceback determination of the cascades which lead
c      to final specific nuclei.thus one can correlate the type of
c      emitted particle,excitation energy, particle energy, spin etc.
c      see detailed working in subroutine statis.
c
c   *** the traceback features enable one to calculate very low
c   *** fission probabilities through dumping fprob, the fission
c   *** probability at each stage of deexcitation in each cascade.
c   *** the sum of fission probabilities is
c   ***   sum ( gamma-f/(gamma-n + gamma-f) )
c   ***    ( sum over all cascades).
c   *** this value is different than the exact value which contains
c   ***        gamma-n / ( gamma-n + gamma-f ) of the previous evap-
c   *** -oration stage. this is somewhat compensated by the fact that
c   *** if the fission probability is large, the actual fission in
c   *** the monte-carlo calculation will decrease the calculated value.
c   ***
c   *** in any case , for high fission probability you should use the
c   *** monte-carlo value and not the calculated one.
c
c 12) execution speed of subroutine chnprb (where the program spends
c     30 percent of its time) has been improved substituting  do loops
c     in place of if's.
c
c     additional important modification, february, 1982
c     *************************************************
c
c     the internal loop over the running sum of tl-s in chnprob has
c     been replaced by a difference of sums in a way which is rigorous
c     for spins of 0 and 1/2. the old part has been deleted. if you
c     consider evaporation of particles with spin greater than 1/2
c     please request special version of subroutine chnprb from avigdor
c
c     this has reduced running time again by 30 percent
c
c     *************************************************
c
c 13) runtime/io-time has been improved by storing the entire mass
c     table in memory. when running on older machines with minimal
c     memory you can revert to reading masses from file each time
c     you need one.
c
c 14) subroutine tlll has been modified for optimal vectorization.
c     this should not have any effect on non vector machines.
c     the same comment applies to subroutine gclvd.
c
c 15) the level density table has been modified at low excitation
c     energies to have a finer granularity. this prevents emission
c     of spurious alphas at very low energies.
c     see subroutine gclvd for details.
c
c 16)  ***** major modification 1/24/84 ********
c
c     subroutine tlom (old optical model) found to give
c     inaccurate results for alphas
c     replaced by subroutine tccal package courtesy of e.d. arthur,
c     los alamos national lab. , group t-2.
c
c     tape4 (fileno. 4) contains all optical model calculation info
c     as well as the level density table parameters as they change
c     throughout the deexcitation cascade.
c
c *** significant change . oct 15 1981.
c     ********************************
c
c     saddle point level density now rigorously calculated.
c     previously, g.s. level density was exponentiated by sqrt(af/an)
c     fission probabilities may differ by 20 percent.
c
c 17)  ***** modification 2/11/84
c     barfac factor deleted. fission is now calculated at
c     any barrier height.
c
c 18)  ***** important modification 2/11/84
c     c.p.s. rotating liquid drop barriers removed. subroutine fisrot
c     replaced by subroutine barfit containing arnie sierks modified
c     rotating nucleus barriers and yrast lines.
c
c      end modification list ************************************
c
c   the results of the calculations are given as distributions
c          of the residual nuclei.  energy and angular distributions of
c          the emitted particles and recoil energy distributions
c          of the residual nuclei are also available by supplying
c          an appropriate control, and energy distributions of gammas
c          emitted are also available by supplying an appropriate
c          control.
c                 input cards
c*************************
c
c     card 1.read ncasc,idist,input,fyrst,fissbr,aratio,mdir,facla,itrac
c           , noshl
c     *******
c              format(3i5,3f5.0,i5,f5.,2i5)
c  ncasc - number of cascades. (events in monte carlo calculation)
c  idist -  =1  detailed angular and energy distribution of residual
c               nuclei and evaporated particles.
c           =0  brief, schematic results of particle spectra and list
c               of evaporated (residual) nuclei population.
c  input -  =1  projectile + target input . agraz parameter (see below)
c               determines diffuseness of partial wave distribution
c           =2  compound nucleus input for single spin.
c           =3  compound nucleus input. spin distribution read in.
c           =4  compound nucleus input. spin distribution calculated
c               taking spin-cutoff paramater at given excitation energy.
c               (from gilbert and cameron).
c           =5  triangular  (sigma = 2l+1 ) cross section between
c                lminn and maximum spin
c        input can also be equal to minus these values.
c        in this case, the cascade dump or detailed trace of
c        the emission points in the e-j plane of all particles
c         leading to a given final product is supressed.
c         for this,  file 2 should be defined accordingly.
c  fyrst  -  parameter determining yrast line to be used.
c  ******** default modification 3/16/84.
c           gilbert and cameron is no longer a default option.
c           since it has caused many strange results,you have to
c           request it specifically - facla lt.0 provides the g and c
c           little-a.
c                                     fyrst lt. 0 provides the g-c
c           yrast line.
c       l t. 0.  gilbert-cameron spin cutoff parameter.
c                erot = (spin)**2/(2.*sigsq)
c       n e. 0.  erot = rotating liquid drop rotational energy,
c                multiplied by factor of fyrst.
c           =0.  value changed to fyrst = 1.
c            in both cases level density calculated at e = ex-erot.
c  fissbr  -  the program assumes the a.j.sierk modified rotating
c             liquid drop barrier if this is left blank. if you provide
c             a fission barrier of your own, the sierk barrier will be
c             renormalized accordingly.
c             if fissbr is positive it will be taken as the desired
c             zero spin fission barrier.
c             if fissbr is negative, its absolute value will be taken
c             as a factor to multiply the sierk barrier.
c  aratio  -  ratio of the fermi gas level density parameter 'little a'
c             at the saddle point to the ground state value.leave blank
c             if you want aratio = 1.000. the saddle point level dens-
c             -ity is determined by g.s. 'little-a' x aratio
c  *** this is major change *** 10/15/81
c  *** previously fission level density calculated in more approximate
c  *** manner. fission probability change of 20 percent upward results.
c  mdir    -  = 0 . compound nucleus is initially in  m=0 states and
c                   the  z-axis is the recoil axis.
c             = 1 . compound nucleus is initially in  m=j states, the
c                   z axis is perpendicular to recoil direction (and
c                   reaction plane). this calculation is appropriate for
c                   deep inelastic fragment deexcitation.
c  facla  -        level density parameter = mass/facla if not zero.
c                  if zero  little-a = mass/7.5
c                  if  .lt. 0  gilbert and cameron  value used
c  itrac -  this card controls the degree of event traceback if you have
c           opted for it.( i.e. 'input' is positive)
c
c           = 0 produces compact traceback, summed over all residues
c           = 1 detailed traceback leading to each individual isotope
c               separately
c           *** remember, to supress traceback completely,use negative
c           *** value of 'input'
c  noshl - =0 uses regular (wapstra) mass table off disk
c          =1 uses lysekil masses with shell correction zeroed for
c              the compound nucleus and its decay products, but  n o t
c              for the projectile and target where experimental masses
c              are used.
c             (this is believed to be appropriate for high ex when
c             shell effects have washed out ).
c
c  end card 1
c
c  card 2
c  *********
c
c      there are 2 possibilities . input = 1, or input gt 1
c  input=1
c        read izp,iap,izt,iat,sp,st,qcn
c        format(4i5,2f5.0,2f10.0)
c  izp,iap  -  z and a of projectile
c  izt,iat  -  z and a of target
c  sp, st   -  ground state spins of projectile and target
c  qcn      -  q value of reaction.
c          -  if left zero it is calculated from mass tables.
c
c  for input =1, an additional card 2a is necessary
c  card 2a  (if(input.eq.1))
c  *******
c         read elab,expsig,agraz,eloss,lminn
c         format(2f10.5,5x,2f5.0,i5)
c  elab  -  beam energy
c  expsig-  experimental fusion cross section if known. tl-s from
c           optical model shifted to reproduce this value if inputted,
c           preserving the l-diffuseness.
c  ******** expsig = 0 results in bass model(p.r.l. 1977) fusion
c           cross section being used.
c  agraz -  to bypass input channel optical model routine (tlom)
c           specify l - diffuseness of fusion cross section
c           if agraz eq 0. diffuseness will be set to 0.5  which is
c           essentially sharp cutoff
c  eloss -  energy loss of beam thru full target width. ( total de )
c           energies will be distributed between ebeam & ebeam-eloss
c  lminn    lowest partial wave l in calculation.
c           partial waves from l=0 to lminn excluded, enabling low-l
c           non-fusion window in reaction calculation.
c
c     ***** note - if you are running at high bombarding energies
c      for which the grazing angular momentum is above 75 hbar,
c      it is recomended to input  agraz gt 0, and to specify
c      an arbitrary value for expsig (or 0 = bass) which corresponds to
c      a fusion cross section with a limiting l-value around 80.
c      this will give you all the evaporation residue data and the
c      fission probabilities you need. for j gt 80  all nuclei will
c      fission anyway, and you will run out of dimension if you try.
c
c  card 2  input=2, 3, 4, or 5
c  ******
c          *******
c          read iz,ia,eexcn,erec,ajnuc,lminn
c          format(2i5,3f10.0,i5)
c  iz,ia,eexcn - z, a, excitation energy of compound nucleus
c  erec,ajnuc  - recoil energy and maximum spin of compound nucleus
c  lminn       - lowest partial wave in calculation
c                ( see above )
c  if input=3, an additional card 2a is necessary
c  card 2a (if input=3)
c  *******
c          read maxc,(sigjc(i),i=1,maxc)
c          format(i5,5x,7f10.5/(8f10.5))
c  maxc  - maximum spin index (i.e. j+1 for even nucleus, j+0.5 for
c          odd nucleus
c  sigjc - partial cross section for specific j
c
c  end card 2  (and 2a for input=1 and input=3)
c
c   optical model potentials are  read in
c       3 potentials are read in -
c 1) neutron. 2) proton. 3) alpha.
c
c for neutrons, protons, or alphas you may leave a blank card in stead
c of input. the program will take optical potential(s) from systematics.
c
c   cards 3a,3b,3c
c   *****************
c   each card is a complete set of potential parameters.
c     read    ,(av(j),j=1,3),(aw(j),j=1,3),(awv(j),j=1,3),rr,ar,ri,ai,
c    1riv,aiv,rcx
c     format(16f5.0)
c   depth of real nuclear well is given by
c          v0=av(1)+av(2))*ecm+av(3)*ecm**2
c   depth of surface-imaginary nuclear well is given by
c          w0=aw(1)+aw(2)*ecm+aw(3)*ecm**2
c   depth of volume-imaginary nuclear well is given by
c          w0=awv(1)+awv(2)*ecm+awv(3)*ecm**2
c          if rr,iri.lt.1.8   radius of real well is rr*amt**(1./3.)
c                        radius of imaginary well is ri*amt**(1/3)
c          if rr,ri.gt.1.8    radius of real well is  rr
c                        radius of imaginary well is    ri
c
c       same goes for coulomb radius rcx
c       all radii should conform to same type
c
c    card 4
c    *************
c
c    read fge1,fgm1,fge2,fgm2
c    format(20x,4f5.0)
c
c    fge1,fgm1,fge2,fgm2, - gamma transition strengths in weisskopf
c    units for e1,m1,e2,m2 transitions respectively.
c     leave blank to obtain program defaults.
c
c   cards 5,6,...  (if you want discrete levels in some level densities
c   **********
c      read ndl ,izdl ,indl
c      format(3i5)
c      read (edl(i),jdl(i),i=1,ndl)
c      format(6(f5.0,i5))
c      ndl is number of discrete levels for nucleus z=izdl ,n=indl
c      edl,jdl are energies of the levels and the spin indices of them.
c      *****  jdl= spin + 1  for even a.  spin + 1/2  for odd a.
c      ****   jdl = 1,2,3,...  always (never 0).
c
c      *******  terminate by blank card. repeat both cards for each
c      nucleus , up to 6 nuclei. maximum levels per nucleus = 36.
c     dimensions need to be changed for high ex or high spin
c     this is arranged to be done with a single ''change'' of a
c     text editor
c      max lev. den. dimension (rlev) =( emax x (maxc+10))
c      if greater than 19696, change all 19696 numbers
c     throughout the program to the larger number
c     emax is maximum excitation energy of compound nucleus
c     maxc  is maximum spin used in program.
c     max  energy dimension =emax +ixpr (no. of exptl. discrete levels)
c     if greater than 286 change  a l l   286-s in program to larger
c     number.
c
c     **********  jcl (job control) comments (necessary files)
c
c
c      disk file 2
c
c      this is the file which contains the full event traceback.
c
c
c      disk file 3
c
c      this is the mass defect file.
c
c
c      disk file 4
c
c      this file contains optical model and level density parameters
c      which may be usefull for diagnostic purposes
c
c
c     this file is in card image format and can be read with text editor
c      **** no code is bug-free. much effort has been made to provide
c      **** a reliable code but no guarantees can be given. please
c      **** consider your output carefully and see whether it makes
c      **** sense.
c      ****
c     good luck. call avigdor about bugs,problems etc.
c     fts no. 843-5475   commercial  505-667-5475
c
c     **********************************
c
      external xdata
      common /mass/ flmass, izpart(4), inpart(4), izr(4), inr(4), nfile,
     1 barfac, aratio, fbarr(110), lfiss
      common /gam/ fge1, fgm1, fge2, fgm2, ergc, discpr, iadef, ndisc,
     1 ergmin, ndl(6), izdl(6), indl(6), jdl(36,6), edl(36,6)
      common /srch/ prob(2222), iprob, mj(2222), emax(4), emin(4), be(4)
     1 , imodl, ge1, ge2
      common /pot/ att(5), av(4), aw(4), awv(4), aso(4), aw1(4), epp(36)
     1 , irad, imt, ecut, ethrs, riv, aiv, rsot, asot, imesh, rcx, vqs(3
     2 ), vq0(3)
      common /sch/ maxj(286,4), ebin(286,4), maxjs(286,4), rmass(4),
     1 mebin(4)
      common /dat/ spart(4), pmass(4)
      common /out/ izcs(9996), incs(9996), jcs(9996), ergcs(9996), mjcs(
     1 9996), rlev(19696,4)
      common /xqsig/ sigma, ncasc
      common /qmmj/ maxc, fyrst, facla, erot(110)
      common /xqres/ iz, in, energy, vzc, itrac, idist
      common /xqang/ sum(300,10), mdir, ct
      common /xqshl/ noshl, alev(100), ajcn(150)
      dimension probm(4), exrsm(4), ixrsm(4), ixpr(4), ixmin(4), tl(30),
     1 nspc(4,51), ibuf2(9), imgz(9996)
      real*8 prob,sprob,tprob,rlev,probm
c      call errset(219,1,1,0,0,220)
      idim=19696
      istop=1
cibm  call jstime(jtime1)
c     call gettr (time1)
      call plm
      modes=4
      imodl=4*19696
      be(4)=0.
      ecloss=0.
      emin(1)=.01
      emin(4)=0.
      print 850
      ibound=0
      inotgs=0
      inotg1=0
      snotgs=0.
      enotgs=0.
      nodes=modes-1
      read 940, ncasc,idist,input,fyrst,fissbr,aratio,mdir,facla,itrac
     1 ,noshl
      idump=0
      if (input.lt.0) idump=1
      input=iabs(input)
      if (itrac.lt.0) idump=1
      itrac=iabs(itrac)
      if (ncasc.gt.9996) ncasc=9996
c      aratio is (little-a-saddle)/(little-a-g.s.)
c      default for no input is 1.000
      if (aratio.eq.0.) aratio=1.
      acasc=ncasc
      if(acasc.lt.100.)acasc=100.
c     discpr is discriminator on probability vector (prob) to prevent
c            search for highly unlikely events. change to 0.001 if you
c            want rare events and are using good statistics.
      discpr=1./acasc
      call compos (iz,in,energy,vzc,maxc,sigma,ecloss,input)
c     compos forms initial spin distribution depending on input
      iadef=2
      if (in.ge.77.and.in.le.110) iadef=1
      if (in.gt.138) iadef=1
      if (iadef.eq.2.and.facla.lt.0.) print 860
      if (iadef.eq.1.and.facla.lt.0.) print 870
c      iadef=2 is assumed spherical nuclei for level densities.
c      iadef=1 is considered a deformed nucleus.
c      the decision is made only for the z,n of the initial nucleus.
      izmax=iz
      inmax=in
      call masses (iz,in,flmass)
      barfac=1.
      ia=izmax+inmax
      a=ia
      iz=izmax
      call barfit (iz,ia,0,brr,segs,selmax)
      if (fissbr.eq.0.) go to 20
      barfac=fissbr/brr
      if (fissbr.gt.0.) go to 10
      barfac=-fissbr
c     negative fissbr means it is the factor in abs. value
      fissbr=brr*barfac
c     barfac multiplies the sierk barrier to get your own.
   10 print 880, fissbr,barfac,aratio
   20 if (fissbr.eq.0) print 890, aratio,brr
      if (facla.lt.0.) print 900
      if(facla.eq.0.)print 902
      if(facla.eq.0.)facla=7.5
      if (facla.gt.0.) print 910, facla
      if (fyrst.lt.0.) print 810
      if(fyrst.eq.0.)print 812
      if(fyrst.eq.0.)fyrst=1.
      if (fyrst.gt.0.) print 800, fyrst
      if (noshl.eq.1) print 920
      temp=sqrt(energy*10./a)
      a3=a**.3333
      rd=a3+1.
      vcoulp=izmax/rd
      rd=a3+1.587
      vcoula=2.*izmax/rd
      emax(1)=16.
      emax(4)=16.
c     ***************************
      emax(3)=vcoula*1.1+5.*temp
      emin(2)=.4*vcoulp
      emin(3)=.5*vcoula
      emax(2)=vcoulp*1.1+5.*temp
c  the emin are the smallest energies for which reliable cross sections
c          are calculable.  the formulae were obtained from experience.
c     emax are maximum emission energies so as not to waste time
c     computing unlikely events
      do 30 ikq=1,4
      do 30 jkq=1,51
   30 nspc(ikq,jkq)=0
      print 930, (emin(i),i=1,4),(emax(i),i=1,4)
      print 950, ncasc,discpr
      if (idist.eq.0) go to 40
      call outem (1,i,a,a)
c     outem is particle emission output subroutine
c
   40 print 960
      do 60 i=1,3
c     ***************************
      read 820, (av(j),j=1,3),(aw(j),j=1,3),(awv(j),j=1,3),rr,ar,ri,ai
     1 ,riv,aiv,rcx
      if (av(1).eq.0..and.aw(1).eq.0.) call optpot (i,in,iz,rr,ar,ri,ai)
      imesh=100
      irad=1
      if (rr.gt.1.8.or.ri.gt.1.8.or.riv.gt.1.8) irad=0
      ecut=0.
      ethrs=0.
      rsot=0.
      asot=99.
      if (ar.eq.0.) ar=99.
      if (ai.eq.0.) ai=99.
      if (aiv.eq.0.) aiv=99.
c     zero diffuseness is default for no potential, but program tries to
c     calculate anyway. so a-s have finite value to prevent blowup
      imt=1
c      imt=1  is saxon wood surface derivative
c      imt=2 (change program) is gaussian surface potential
      att(1)=rr
      att(2)=ar
      att(3)=ri
      att(4)=ai
      do 50 j=1,4
      aw1(j)=0.
      aso(j)=0.
   50 continue
      call tlll (edum,idum,izmax,inmax,tl,idum,tl)
c     initialize transmission coeff routine
   60 continue
      if(ncasc.eq.0)stop
c     tlll prepares all n, p, and alpha transmission coefficients
c     in advance.
c     the change in transmission coefficients for second and sub-
c     -sequent decaysis simulated by shifting the energy at which
c     the tl-s are calculated according to the shift of the coulomb
c     barrier of the secondary residual nuclei.
c     tl-s at this energy are calculated by interpolation
c
c     **************************
      read 980, fge1,fgm1,fge2,fgm2
      if (fge1.ne.0.) go to 100
      if (a.gt.100.) go to 70
      fge1=.00008
      fge2=4.8
      fgm1=.025
      fgm2=.0195
      go to 100
   70 if (a.gt.126.) go to 80
      fge1=.000014
      fge2=5.9
      fgm1=.01
      fgm2=.00088
      go to 100
   80 if (a.gt.150.) go to 90
      fge1=.000046
      fge2=7.7
      fgm1=.007
      fgm2=.058
      go to 100
   90 fge1=.000011
      fge2=9.0
      fgm1=.010
      fgm2=1.2
  100 print 970, fge1,fgm1,fge2,fgm2
      if(input.eq.2)go to 101
      print 830
      print 840, (ajcn(i),alev(i),i=1,maxc)
  101 ndisc=1
c      discrete levels (if desired) for various level densities
c      are now read in. the number of final nuclei
c      that have discrete levels is ndisc.
c      note *****  ndisc should not be more than 6.
c           *****  the maximum number of levels for each nucleus is 36
  110 read (5,990,end=120) ndl(ndisc),izdl(ndisc),indl(ndisc)
      nll=ndl(ndisc)
      if (nll.gt.0.and.nll.le.36) go to 130
  120 ndisc=ndisc-1
      go to 140
  130 read 1000, (edl(i,ndisc),jdl(i,ndisc),i=1,nll)
      print 1010, izdl(ndisc),indl(ndisc),(edl(i,ndisc),jdl(i,ndisc),i=1
     1 ,nll)
      ndisc=ndisc+1
      go to 110
  140 if (input.eq.2) go to 150
c     input = 1  - projectile + target
c           = 2  - comp. nuc. , ex and j fixed
c           = 3  - comp. nuc. , ex fixed, j distributed
c           = 4  - comp. nuc. ,j distibution calculated internally
  150 id=1
c     prepare tables for recoiling nucleus.
      if (idist.ne.0) call moment (0,0.,0.,0,vzc,ncasc,id)
c
      if (input.eq.2) go to 210
c   the spin distribution of the starting excited nuclide has been
c          determined.  the determination of the probabilities of
c          the paths of deexcitation will be determined in the following
c          way
c         1) using random numbers, ncasc channels are established for
c            the compound nucleus.
c         2) the highest a in the array is determined.
c         3) the first nucleus of this a in the array is determined.
c         4) the highest energy for this nucleus is determined.
c         5) the first spin of this nucleus and energy is found.
c         6) the probability matrix for this nucleus energy and spin
c            is determined.
c         7) by means of random numbers, the deexcitation path is chosen
c            for all members of the array of that nuclide, energy and
c            spin.  the new z, n, energy, and spin are substituted for
c            the old values in the array.
c         8) the next spin for the (old) nuclide and energy is located.
c            if it exists, go to step 7.
c         9) the next higher energy for that (old) nuclide is located.
c            if it exists, got to step 5.
c        10) the next nuclide for that (old) a is located.
c            if it exists, go to step 4.
c        11) go to step 2.
c
c   sprob   is unnormalized sum of probabilities
c   prob    (equivalent - prob) is unnormalized probabilities
c   modes   is number of different particles (and gammas)
c              considered for evaporation
c   izpart  and  inpart  are z and n of emitted particles
c     the vector associated with each event is as follows
c
c     ergcs(j)- the excitation energy
c     izcs(j)- the atomic number ( z )
c     incs(j) - the number of neutrons ( n )
c     jcs(j)  - the spin index  i.e.  1,2,3,   etc.
c               spin = spin index - 1  for even a
c               spin = spin index - 1/2 for odd a
c     mjcs(j) - the spin projection  =  0,1,2,   spin.
c              we neglect 1/2 integer values for m only *******
c
c   spin distribution of initial nuclide normalized to sum=1.
      if (maxc.gt.100) maxc=100
c     for j.gt.100 there is no compound nucleus in practice.
      do 160 i=2,maxc
c     alev  is partial xsection for spin index i
c           (from  compose  thru common)
  160 alev(i)=alev(i)+alev(i-1)
      slev=alev(maxc)
      do 170 i=1,maxc
  170 alev(i)=alev(i)/slev
c     running sum for j random selection formed
      jjl=1
      do 200 l=1,ncasc
      izcs(l)=iz
      incs(l)=in
c      l  is no. of cascade being calculated
      ergcs(l)=energy-ecloss*ranf(0.)
      njl=alev(jjl)*ncasc+.5
      if (l.le.njl) go to 190
  180 jjl=jjl+1
      if (jjl.ge.maxc) go to 190
      njl=alev(jjl)*ncasc+.5
      if (l.gt.njl) go to 180
  190 jcs(l)=jjl
      mjcs(l)=(jjl-1)*mdir
  200 continue
c     print 8881,(i,incs(i),jcs(i),i=1,ncasc)
c8881 format((4(3i5,5x)))
      go to 230
  210 do 220 l=1,ncasc
      izcs(l)=iz
      incs(l)=in
      ergcs(l)=energy
      mjcs(l)=(maxc-1)*mdir
  220 jcs(l)=maxc
c     compound nucleus preparation finished
  230 call track (lk,lk,lk,eu,eu,eu,energy,1,sprob,ecloss)
      lat=1
      mat=ncasc
      iac=0
      izc=0
c
c     loop to find highest ex for highest a
c
c     the advantage of this particular code is that all events with
c     same a, z, ex and j use the same decay probability table
c     calculated by chnprb. this saves much computing time.
c

  240 do 310 i=lat,mat
      if (ergcs(i).eq.0.) go to 310
      ia=izcs(i)+incs(i)
      if (ia-iac) 310,260,250
  250 iac=ia
      la=i
c     lowest limit for a
c        loop to find highest a for particular z
      izc=izcs(i)
      inc=incs(i)
      ergc=0.
      ergmin=9999.
c      ergmin is minimum energy for given a
      go to 270
  260 if (izc.eq.0) go to 250
  270 ma=i
c     highest limit of a
      izcsi=izcs(i)
      incsi=incs(i)
      if (izcsi.ne.izc.or.incsi.ne.inc) go to 310
      if (ergcs(i)-ergc) 300,290,280
  280 ergc=ergcs(i)
c     find highest ex for highest a
      lerg=i
  290 merg=i
c     for the above la , ma -  lerg , merg are indices of highest
c       energy limits
  300 if (ergmin.gt.ergcs(i)) ergmin=ergcs(i)
  310 continue
      if (iac.eq.0.and.lat.eq.1.and.mat.eq.ncasc) go to 710
      if (iac.lt.6) go to 710
c
c     when all energies are zero  -
c                              the whole job is finished
c
      if (izc.eq.0) go to 230
      sc=.5*mod(izc+inc,2)
c     sc is spin index of nucleus -  0 or 1/2
      do 320 i=1,4
      ii=5-i
      call amass (rmass,izc,inc,ii,mebin(ii),maxj(1,ii),maxjs(1,ii),ebin
     1 (1,ii),rlev(1,ii),ixpr(ii),be)
  320 continue
      vqs(1)=0.
      vqs(2)=1.44*izr(2)/((izr(2)+inr(2))**.3333)
      vqs(3)=2.88*izr(3)/((izr(3)+inr(3))**.3333)
c     get mass excesses of nuclei
c
c   rmass  are masses of residual nuclides
c  the emin are the smallest energies for which reliable cross sections
c          are calculable.  the formulae were obtained from experience.
c
      do 330 i=1,nodes
      exrsm(i)=ergc+be(i)-emin(i)
c     ex of rsid.nuc. - maximum   maximum excitation energy of
c     particular residual nucleus
  330 continue
      exrsm(4)=ergc
      do 340 k=1,modes
      if (exrsm(k).le.0.) go to 340
c     for  ex.lt.0.  mode not available
      if (exrsm(k).le.ebin(mebin(k),k)) go to 340
      print 1020, izr(k),inr(k),exrsm(k),ebin(mebin(k),k),mebin(k)
      stop
  340 continue
  350 jc=jcs(lerg)
c    jc  is cascade spin index
      do 360 i=1,nodes
      if (exrsm(i).gt.0.) go to 390
  360 continue
c
c  evaporation of particles of particular z,n with energy e or less no
c      longer possible
c
      if (exrsm(4).le.0.) go to 370
      ibound=1
      go to 390
  370 lat=la
      mat=ma
      do 380 i=lat,mat
c     zero the energy of the cascade to avoid further calculations
      izcsi=izcs(i)
      incsi=izcs(i)
      if (izcsi.eq.izc.and.incsi.eq.inc) ergcsi=0.
  380 continue
      izc=0
      go to 240
c
c                                     new z,n of max a needed
c
  390 do 410 k=1,modes
c     a and ex in predetermined limits. j fixed. particles
c     and gammas can be emitted.
      if (exrsm(k).gt.0.) go to 400
      ixrsm(k)=0
      go to 410
  400 call find (exrsm(k),ixrsm(k),maxj(1,k),maxjs(1,k),ebin(1,k),mebin
     1 (k),k)
c     the following is a proceedure to find the integer index
c     associated with a given excitation energy in the level
c     file. ebin is the energy of the bin ixrsm.
      exmin=exrsm(k)-emax(k)+emin(k)-1.
      call find (exmin,ixmin(k),maxj(1,k),maxjs(1,k),ebin(1,k),mebin(k)
     1 ,k)
  410 continue
  420 do 430 j=1,2222
c     zero the probabilities before calculation
  430 prob(j)=0.
      fprob=0.
      sprob=0.
      iprob=0
      iprobp=0
c     iprobp counts particle output
      do 440 k=1,modes
      probm(k)=0.
c     probm is maximum probability for mode (channel) k.
      if (exrsm(k).le.0.) go to 440
      amr=izc+inc
c     chnprb gives probability for all cells in mode k
      call chnprb (k,ixrsm(k),jc,sc,spart(k),rlev(1,k),probm(k),maxj(1,k
     1 ),maxjs(1,k),ebin(1,k),amr,ixmin(k),mebin(k),ixpr(k))
  440 continue
      if (iprob.eq.0) go to 480
c    an excited trap
c     delta-j too large to decay.
c
      fprob=prob(iprob)*lfiss
      do 450 i=2,iprob
      prob(i)=prob(i)+prob(i-1)
  450 continue
      sprob=prob(iprob)
      if (sprob.ne.0.) go to 460
      iprob=0
      go to 480
  460 do 470 i=1,iprob
  470 prob(i)=prob(i)/sprob
      fprob=fprob/sprob
  480 mz=lerg
      mg=merg
      lerg=0
      il=0
      do 500 i=mz,mg
      qif=iabs(incs(i)-inc)+iabs(izcs(i)-izc)+abs(ergcs(i)-ergc)
      if (qif.ne.0.) go to 500
      if (jc.ne.jcs(i)) go to 490
      il=il+1
      imgz(il)=i
      go to 500
  490 if (lerg.eq.0) lerg=i
  500 continue
c     above loop sets up following:
c     il   - number of cascades having right z,n,spin,ex to use same
c            probability distribution
c     imgz - the running index number of these cascades
c     lerg - first cascade having different spin (same z,n,ex)
      if (il.eq.0) go to 630
      do 620 icsc=1,il
      i=imgz(icsc)
c     scan all events looking for ones with right a,z,ex,j
      izcsi=izcs(i)
      incsi=incs(i)
      jcsi=jcs(i)
      mjcsi=mjcs(i)
c     following insert forces decay thru e2 gamma yrast cascade
      if (iprob.ne.0) go to 540
      ergnot=ergcs(i)
  510 if (jcs(i).gt.1) go to 520
      ergcs(i)=0.
      jcs(i)=0.
      jcsi=0
      go to 610
  520 ajn=jcs(i)-1.+sc
      mode=4
      etran=ergcs(i)*(4.*ajn-2.)/(ajn*(ajn+1.))
      ergc=ergcs(i)
      ergcs(i)=ergcs(i)-etran
      if (ergcs(i).lt.0.) ergcs(i)=0.
      if (jcs(i).ne.jc) go to 530
      inotgs=inotgs+1
      ergnot=ergc
  530 jcsi=jcs(i)
      jcs(i)=jcs(i)-2
c     not decayed to g.s. (''yrast traps''- no e3 or m3 in program)
      go to 560
c     end yrast cascade insert.  15 march 1979
  540 call search (izcs(i),incs(i),ergcs(i),jcs(i),izc,inc,mode)
      if (ergcs(i).ge.0.) go to 560
      ergcs(i)=0.
      ergcsi=ergcs(i)
      call track (5,jcsi,jcsi,ergc,ergcsi,0.,ergc,2,sprob,ecloss)
      if (idump.eq.1) go to 550
c     setup output buffer for traceback
      ibuf2(1)=i
      ibuf2(2)=6
c     mode=6 is fission traceback
      ibuf2(3)=jcsi
      ibuf2(4)=0
      ibuf2(5)=mjcsi
      ibuf2(6)=ergc+1.
      ibuf2(7)=0
      ibuf2(8)=izc
      ibuf2(9)=inc
      write (2) ibuf2,fprob
  550 izcs(i)=-izcs(i)
c     negative izcs(i) implies fission of izcs(i)
      go to 620
c     search determines final nucleus
c
c     proceed to catalogue distributions.
c
  560 ep=ergc-ergcs(i)+be(mode)
      if (iprob.eq.0) ep=etran
      jcsf=jcs(i)
      ergcsi=ergcs(i)
      call track (mode,jcsi,jcsf,ergc,ergcsi,ep,ergc,2,sprob,ecloss)
c     mainstream traceback of general decay features
      ac=izcsi+incsi
      ap=izpart(mode)+inpart(mode)
      ep=ep*ac/(ap+ac)
      if (idump.eq.1) go to 590
      ie12=0
      ibuf2(1)=i
      if (mode.ne.4) go to 580
      jddj=jcsf-jcsi
      jddj=iabs(jddj)
      if (jddj.ne.2) go to 570
      ie12=1
      go to 580
  570 gmsum=ge1+ge2*ep**2
      gmsum=ge1/gmsum
c     if delta-j is 2 decay is e2
c     otherwise, decay mode (for traceback purposes only)
c     is decided by size of decay probability (e1 or e2 only)
      x=ranf(0.)
      if (x.gt.gmsum) ie12=1
  580 ibuf2(2)=mode+ie12
c      ibuf2(2) now gives mode=4 for e1 gammas and mode=5 for e2 gammas
c      *** this is for the dump-traceback only
      ibuf2(3)=jcsi
      ibuf2(4)=jcsf
      ibuf2(5)=mjcsi
      ibuf2(7)=ep*10.+1.
      ibuf2(6)=ergc+1.
      ibuf2(8)=izc
      ibuf2(9)=inc
      write (2) ibuf2,fprob
  590 jkq=ep+1.
      if (jkq.gt.50) jkq=50
      nspc(mode,jkq)=nspc(mode,jkq)+1
      if (idist.eq.0) go to 610
      if (mode.ne.modes) go to 600
      jcsff=jcsf-1
      if (mjcs(i).lt.-jcsff) mjcs(i)=-jcsff
      if (mjcs(i).gt.jcsff) mjcs(i)=jcsff
      ig=2
      if (ibound.eq.1) ig=4
c     store gamma energy
      call outem (ig,mode,ep,dum)
      go to 610
  600 call mjran (jcsi,jcsf,mjcsi,mjcsf,mode,ep,izcsi,incsi)
      mjcs(i)=mjcsf
      id=2
      call moment (i,ac,ap,mode,ep,ncasc,id)
c     store particle energy. forms output distributions by
c     call ing outem.
  610 if (ergcs(i).gt.0..and.iprob.eq.0) go to 510
  620 continue
c 630 call timend (ncasc,istop)
c     if (istop.eq.0) go to 710
  630 if (iprob.ne.0) go to 640
c     print3470,izc,inc,iac,ergc,jc,ibound,inotgs
      enotgs=enotgs+ergnot*(inotgs-inotg1)
      snotgs=snotgs+jc*(inotgs-inotg1)
      inotg1=inotgs
c     inotgs=0
c     remove previous 2 c-s in order to obtain details about traps
  640 if (lerg.eq.0) go to 650
c
c     no further spin at this a and z
c                                     new energy or nuclide needed
c
      jc=jcs(lerg)
      go to 420
c
c                                     new spin needed
c
  650 ibound=0
      lat=la
      mat=ma
      ergc=0.
      do 680 i=lat,mat
      if (ergcs(i).eq.0.) go to 680
      if (inc.ne.incs(i).or.izc.ne.izcs(i)) go to 680
      if (ergcs(i)-ergc) 680,670,660
  660 ergc=ergcs(i)
      lerg=i
  670 merg=i
  680 continue
      if (ergc.eq.0.) go to 700
      do 690 k=1,nodes
  690 exrsm(k)=ergc+be(k)-emin(k)
      exrsm(4)=ergc
      go to 350
c
c                                     same z,n new energy needed
c
  700 izc=0
      go to 240
c
c                            new nuclide of same a needed
c
  710 continue
      if (idist.eq.0) go to 730
      do 720 i=1,ncasc
      ac=izcs(i)+incs(i)
  720 if (izcs(i).gt.0) call moment (i,ac,0.,0,0.,0,3)
c     last parameter in moment
c        = 1  initialization of tables
c        = 2  collecting data during run
c        = 3  termination and output.
  730 print 1030
      call prodct (idist,izmax,inmax,ncasc)
      if (idist.ne.0) call outem (3,idum,dum,dum)
      if (inotgs.eq.0) go to 740
      enotgs=enotgs/inotgs
      pr=(100.*inotgs)/ncasc
      snotgs=snotgs/inotgs
      print 1040, pr,enotgs,snotgs
  740 print 1050
      do 750 ist=1,4
  750 be(ist)=0.
      do 760 ikq=1,50
      ikq1=ikq-1
      ikz=ikq
      if (ikq.eq.50) ikz=99
      if (nspc(1,ikq).ne.0.or.nspc(2,ikq).ne.0.or.nspc(4,ikq).ne.0.or
     1 .nspc(3,ikq).ne.0) print 1060, ikq1,ikz,(nspc(jkq,ikq),jkq=1,4)
      epart=ikq-.5
      do 760 ist=1,4
      nspc(ist,51)=nspc(ist,51)+nspc(ist,ikq)
      be(ist)=be(ist)+epart*nspc(ist,ikq)
  760 continue
      do 770 ist=1,4
  770 be(ist)=be(ist)/(nspc(ist,51)+1.e-9)
      print 1070, (nspc(jkq,51),jkq=1,4)
      print 1080, be
      do 780 ist=1,4
      be(ist)=nspc(ist,51)
  780 be(ist)=be(ist)/ncasc
      print 1090, be
      call track (lk,lk,lk,ergc,ergc,ergc,ergc,3,sprob,ecloss)
      print 1100
c     if (istop.eq.0) stop
c     you can get times for different code segments by activating
c     the following
cibm  call jstime(jtime2)
c     call gettr (time2)
      t1=(time1-time2)*1.0
c     print 1110, t1
      if (idump.eq.1) stop
      call statis
cibm  call jstime(jtime3)
c     call gettr (time3)
      t2=(time2-time3)*1.0
c     print 1120, t2
      stop
c
  800 format (35x,'sierk  diffuse surface nucleus yrast line used', 9x,
     1'*****'/55x,' multiplied by factor of ',10x,f5.2)
  810 format (35x,'gilbert and cameron        yrast line used',13x,'****
     1*')
  812 format (35x,'fyrst = 0.   changed to default value     ',13x,'****
     1*')
  820 format (16f5.0)
  830 format (/,47x,'partial cross sections (mb)'//46x,'j',21x,'sig (j)'
     1 //)
  840 format ((6(1x,f5.1,2x,f6.2,6x)))
  850 format (///////40(1h*),'  p a c e  2  ',40(1h*),' modified julian
     1'///' projection  angular-momentum  coupled  evaporation ','  mont
     2e-carlo  code'//' vax - v.m.s. /  i.b.m.      version.  avigdor,
     331 dec 1985'//' angular distributions obtained using m-states of a
     4ngular momentum '//' sierk fission barriers'//92(1h*)//)
  860 format (35x,'spherical nucleus level density'18x,' **********')
  870 format (35x,'deformed nucleus level density'18x,' ***********')
  880 format (35x,'input fission barrier (mev)   ',24x,f6.2/35x,'** sier
     1k barrier will be multipled by ',17x,f5.2/35x,'little-a sub-f / li
     2ttle-a sub-gamma   ',16x,f6.3)
  890 format (35x,'little-a sub-f / little-a sub-gamma   ',14x,f6.3/35x,
     1 'zero spin fission barrier   ',18x,f6.2)
  900 format (35x,'gilbert cameron little-a selected',22x,' **** ')
  902 format (35x,'default value of facla substituted',21x,' **** ')
  910 format (35x,'little-a  ',39x,' mass /',f4.1/)
  920 format (35x,'lysekil masses with no shell correction used ',11x,'*
     1***')
  930 format (/34x,' ecm range for  neut    prot    alph   gamma (mev)'/
     1 46x,4f8.2,'  min'/46x,4f8.2,'  max')
  940 format (3i5,3f5.0,i5,f5.0,2i5)
  950 format (//35x,'number of cascades   ',33x,i6/35x,'internal probabi
     1lity discriminator of program set to',f8.4/)
  960 format (35x,'optical model parameters for light emitted particles'
     1 /35x,10(1h*),16x,'  outputted to tape 4 (fileno. 4) ')
  970 format (/35x,'e.m. transition strengths in weisskopf units '//35x,
     1 '      e1        m1        e2        m2  '/35x,4(1x,f9.5))
  980 format (20x,4f5.0)
  990 format (3i5)
 1000 format ((6(f5.0,i5)))
 1010 format (///' discrete levels.  z n  ',2i4/(8(f6.2,i4,4x)/))
 1020 format ('0level density table for z ',i3,' n ',i3,' does not go to
     1 high enough energy.  ',f7.3,' required.   .'/f10.3,i7)
 1030 format ('1output results for compound nucleus decay '/' **********
     1*******************************')
 1040 format (//f7.2,' percent of cascades trapped before reaching groun
     1d state due to spin inhibition'//' average energy at which cascade
     2s were trapped is ',f6.2,' mev, average spin   ',f5.1,' hbar'//' *
     3*** successive decays through single yrast cascade assumed')
 1050 format ('1',30x,'c.m.  spectra  of   emitted   particles'/27x,47(1
     1 h-)//27x,'en(mev)  [  neut  [   prot  [  alpha  [  gamma'/27x,47(
     2 1h-)//)
 1060 format (27x,i2,' - ',i2,'  [ ',i5,'  [  ',i5,'  [  ',i5,'  [  ',i5
     1 )
 1070 format (//27x,' total   [ ',i5,'  [  ',i5,'  [  ',i5,'  [  ',i5)
 1080 format (//25x,'av. energy [ ',f5.2,'  [  ',f5.2,'  [  ',f5.2,'  [
     1 ',f5.2)
 1090 format (//23x,'multiplicity [ ',f5.2,'  [  ',f5.2,'  [  ',f5.2,'
     1[  ',f5.2)
 1100 format (//' end evaporation calculation.   shalom '/5(120(1h*)/)/1
     1 h1)
 1110 format (//' time spent in evaporation calculation = ',f7.2//)
 1120 format (//' time spent in traceback subroutine = ',f7.2//)
      end
      block data xdata
      common /mass/ flmass, izpart(4), inpart(4), izr(4), inr(4), nfile,
     1 barfac, aratio, fbarr(110), lfiss
      common /dat/ spart(4), pmass(4)
      data izpart, inpart, spart, pmass /0,1,2,0,1,0,2,0,.5,.5,0.,0.,8.0
     1 7169,7.28922,2.42494,0./
      end
      subroutine tlll (er,np,iz,in,tl,llmax)
c     written by a. gavron
c     generate and extrapolate transmission coefficients
      dimension tlpn(2,16,30), tlal(30,30), tlt(36,80), de(3), ll(3,30),
     1 de1(3), tl(30), nmode(3), lc(36)
      data nmode /4hneut,4hprot,4halph/
c     transmission coefficients written to tape 4
      common /mass/ flmass, izpart(4), inpart(4), izr(4), inr(4), nfile,
     1 barfac, aratio, fbarr(110), lfiss
      common /srch/ prob(2222), iprob, mj(2222), emax(4), emin(4), be(4)
     1 , imodl
      real*8 prob
      common /pot/ att(5), av(4), aw(4), awv(4), aso(4), aw1(4), ep(36),
     1 irad, imt, ecut, ethrs, riv, aiv, rsot, asot, imesh, rcx, vqs(3),
     2 vq0(3)
      data k /0/, ip /1/
      tll=alog(1.e-28)
      if (ip-3) 10,50,100
   10 ia=iz+in
      k=1
      de(ip)=(emax(ip)-emin(ip))/29.
      do 20 ien=1,30
   20 ep(ien)=emin(ip)+(ien-1)*de(ip)
      call tccal (ip,iz,ia,30,15,imt,att,riv,aiv,av,aw,awv,rcx,ep,ethrs
     1 ,ecut,aso,rsot,asot,irad,imesh,tlt,lc)
      do 40 ien=1,30
      lc1=lc(ien)
      sum=0.
      do 30 ij=1,lc1
      aj=ij-1
      sum=sum+(2.*aj+1.)*tlt(ien,ij)
   30 continue
      xsec=679.58/(ep(ien)*ia/(ia+1.))*sum
      write (4,140) ep(ien),xsec,(tlt(ien,kp),kp=1,lc1)
      ll(ip,ien)=lc1-1
      if (ll(ip,ien).gt.15) ll(ip,ien)=15
      do 40 l=1,16
      tlpn(ip,l,ien)=alog(1.e-28)
      if (l.le.lc1.and.tlt(ien,l).gt.0.) tlpn(ip,l,ien)=alog(tlt(ien,l)+
     1 1.e-28)
   40 continue
      ip=ip+1
      return
   50 de(3)=(emax(3)-emin(3))/29.
      do 60 ien=1,30
   60 ep(ien)=emin(3)+(ien-1)*de(3)
      call tccal (ip,iz,ia,30,29,imt,att,riv,aiv,av,aw,awv,rcx,ep,ethrs
     1 ,ecut,aso,rsot,asot,irad,imesh,tlt,lc)
      do 80 ien=1,30
      lc1=lc(ien)
      sum=0.
      do 70 ij=1,lc1
      aj=ij-1
      sum=sum+(2.*aj+1.)*tlt(ien,ij)
   70 continue
      xsec=679.58/(ep(ien)*4.*ia/(ia+4.))*sum
      write (4,140) ep(ien),xsec,(tlt(ien,kp),kp=1,lc1)
      ll(3,ien)=lc(ien)-1
      if (ll(3,ien).gt.29) ll(3,ien)=29
      do 80 l=1,30
      tlal(l,ien)=alog(1.e-28)
      if (l.le.lc1.and.tlt(ien,l).gt.0.) tlal(l,ien)=alog(tlt(ien,l)+1.e
     1 -28)
   80 continue
      do 90 ip=1,3
   90 de1(ip)=1./de(ip)
      vq0(1)=0.
      vq0(2)=1.44*(iz-1.)/((ia-1.)**.3333)
      vq0(3)=2.88*(iz-2.)/((ia-4.)**.3333)
      ip=4
      return
  100 ee=er+vq0(np)-vqs(np)
c     this proceedure adjusts the t-ls to the changing coulomb
c     barrier during deexcitation.
      if (ee.ge.emax(np)) ee=emax(np)-.01
      if (ee.le.emin(np)) ee=emin(np)+.01
      ne=(ee-emin(np))*de1(np)
      dee=ee-emin(np)-de(np)*ne
      ne=ne+1
      llmax=ll(np,ne)
      ll1=llmax+1
      if (np.eq.3) go to 120
c     loops below rewritten for vectorization
      do 110 lp1=1,ll1
      tl(lp1)=0.
      dtt=0.
      t0=tlpn(np,lp1,ne)
      dtt=(tlpn(np,lp1,ne+1)-t0)*de1(np)
      tll=t0+dtt*dee
      tl(lp1)=exp(tll)
  110 continue
      return
  120 tl(1)=0.
      lmt=llmax+1
      do 130 lp1=1,lmt
      tl(lp1)=0.
      dtt=0.
      t0=tlal(lp1,ne)
      dtt=(tlal(lp1,ne+1)-t0)*de1(np)
      tll=t0+dtt*dee
      tl(lp1)=exp(tll)
  130 continue
      k=k+1
      return
c
  140 format (1x,f6.2,2x,f8.2,10e10.2/(17x,10e10.2))
      end
      subroutine search (iz,in,erg,j,izc,inc,mode)
c
c     original julian subroutine
c
c     ***********
c
c     subroutine  determines final residual nucleus based on the generated
c     random number
c
      common /mass/ flmass, izpart(4), inpart(4), izr(4), inr(4), nfile,
     1 barfac, aratio, fbarr(110), lfiss
      common /srch/ prob(2222), iprob, mj(2222), emax(4), emin(4), be(4)
     1 , imodl
      common /sch/ maxj(286,4), ebin(286,4), maxjs(286,4), rmass(4),
     1 mebin(4)
      real*8 prob
      r=ranf(0.)
c      the search algorithm is to keep halving the vector and
c      decreasing limits till index is found. the index of the
c      determines the decay mode
      if (r.le.prob(iprob)) go to 10
      i=iprob
      go to 70
   10 if (r.gt.prob(1)) go to 20
      i=1
      go to 70
   20 is=0
      jf=iprob
   30 jf=jf/2+1
      if (jf.eq.2) jf=1
      is=is+jf
      if (is.gt.iprob) is=iprob
   40 if (r.gt.prob(is)) go to 30
   50 jf=jf/2+1
      if (jf.eq.2) jf=1
      is=is-jf
      if (is.lt.1) is=1
      if (jf.eq.1) go to 60
      go to 40
   60 if (r.le.prob(is)) go to 50
      i=is+1
   70 mode=mj(i)/19696+1
      if (mode.lt.5) go to 80
      erg=-1.
      return
c      ibin is place in 19696 long vector of one specific mode.
   80 iz=izc-izpart(mode)
      in=inc-inpart(mode)
      mji=mj(i)
      ibin=mod(mji,19696)
      if (ibin.eq.0) ibin=19696
      meb=mebin(mode)
      do 90 i=1,meb
      if (maxjs(i,mode).ge.ibin) go to 100
   90 continue
      print 110, mode,meb,maxjs(meb,mode),ibin
      stop
  100 j=maxj(i,mode)-maxjs(i,mode)+ibin
      erg=ebin(i,mode)
      return
c
  110 format (' search error.'/' mode,mebin,maxjs,ibin'/1x,4i5)
      end
      subroutine amass (a,izc,inc,mode,mebin,maxj,maxjs,ebin,rlev,ixpr
     1 ,be)
c
c     modified julian subroutine
c      get binding energies, level densities and fission barriers
c
c     **************
c      determine masses.
      dimension maxj(286), maxjs(286), ebin(286), rlev(19696), be(4), a(
     1 4)
      common /mass/ flmass, izpart(4), inpart(4), izr(4), inr(4), nfile,
     1 barfac, aratio, fbarr(110), lfiss
      common /qmmj/ maxc, fyrst, facla, erot(110)
      common /dat/ spart(4), pmass(4)
      common /gam/ fge1, fgm1, fge2, fgm2, ergc, discpr, iadef
      real*8 rlev
      maxxl=maxc+9
      if (maxxl.gt.110) maxxl=110
      izr(mode)=izc-izpart(mode)
      inr(mode)=inc-inpart(mode)
      call masses (izr(mode),inr(mode),a(mode))
      be(mode)=a(4)-a(mode)-pmass(mode)
      if (mode.ne.4) go to 20
      write (4,30)
      ia=izr(mode)+inr(mode)
      iz=izr(mode)
      do 10 i=1,maxxl
      l=i-1
      call barfit (iz,ia,l,fbarr(i),erot(i),selmax)
      fbarr(i)=fbarr(i)*barfac
      erot(i)=erot(i)*fyrst
   10 continue
   20 call range (mode,imin,imax)
      call levden (izr(mode),inr(mode),imin,imax,maxj,maxjs,ebin,mebin
     1 ,rlev,ixpr,be(mode),1.)
      return
c
   30 format (1h )
      end
      subroutine prodct (idist,iz,in,ncasc)
c
c     original julian subroutine
c
c     ***********
c
c      sorting routine of residual nuclei
c
      common /out/ izcs(9996), incs(9996), jcs(9996), ergcs(9996), mjcs(
     1 9996), rlev(19696,4)
      common /xqsig/ sigma
      common /xqang/ sum(300,10), mdir, ct
      dimension knuc(9996), pres(990), nres(990), nnuc(990), iznuc(990),
     1 innuc(990)
      equivalence (knuc(1),mjcs(1))
      do 10 j=1,990
   10 nres(j)=0
      ires=0
      ifiss=0
      nc=ncasc
      do 50 i=1,nc
      if (izcs(i).ge.0) go to 20
      ifiss=ifiss+1
      go to 50
   20 izc=iz-izcs(i)
      inc=in-incs(i)
      ia=izc+inc
      ibin=ia*(ia+1)/2+izc+1
      if (ibin.le.990) go to 40
      if (ires.eq.1) go to 30
      print 180
      ires=1
   30 print 190, izcs(i),incs(i)
      go to 50
   40 nres(ibin)=nres(ibin)+1
   50 continue
      print 280
      print 200
      print 280
      sres=ncasc
      iat=iz+in
      k=1
      sigs=0.
      nress=0
      press=0.
      do 90 i=1,44
c     if dimension changed from 990, change this       do range.
      ibin1=i*(i-1)/2+1
      ibin2=i*(i+1)/2
      izr=iz
      iat=iz+in-i+1
      do 80 j=ibin1,ibin2
      if (nres(j).eq.0) go to 70
      pres(j)=nres(j)/sres*100.
      inr=iat-izr
      sig=sigma*pres(j)/100.
      call lmnt (izr,ma1)
      print 210, izr,inr,iat,ma1,nres(j),pres(j),sig
      if (pres(j).lt.3.) go to 60
      iznuc(k)=izr
      innuc(k)=inr
      nnuc(k)=nres(j)
      k=k+1
   60 sigs=sigs+sig
      press=press+pres(j)
      nress=nress+nres(j)
   70 izr=izr-1
   80 continue
   90 continue
      pr=ifiss/sres*100.
      sig=sigma*pr/100.
      press=press+pr
      nress=nress+ifiss
      if (ifiss.gt.0) print 220, ifiss,pr,sig
      sigs=sigs+sig
      print 230
      print 240, nress,press,sigs
      print 280
      if (idist.eq.0) go to 170
      id=1
      call outres (id,idum,idum)
      k=k-1
      do 130 i=1,nc
      if (izcs(i).lt.0) go to 130
      do 100 j=1,k
      if (izcs(i).eq.iznuc(j).and.incs(i).eq.innuc(j)) go to 110
  100 continue
      jk=16
      go to 120
  110 knuc(i)=j
      jk=j
  120 ik=i
      id=2
      call outres (id,ik,jk)
  130 continue
      k=k+1
      print 250
      if (k.gt.16) k=16
      if (mdir.eq.0) print 260
      if (mdir.eq.1) print 270
      do 160 i=1,k
      print 280
      if (i.eq.k) go to 140
      print 290, iznuc(i),innuc(i)
      go to 150
  140 print 300
  150 print 280
      id=3
      jk=j
      ik=i
      if (ik.eq.k) ik=16
      call outres (id,jk,ik)
  160 continue
  170 return
c
  180 format (1x,'the following products were not binned'/1x,'  z    n')
  190 format (1x,i3,i5)
  200 format (51x,'yields of residual nuclei'/35x,'z',9x,'n',9x,'a',6x,'
     1events',6x,'percent',6x,'x-section(mb)')
  210 format (34x,i3,7x,i3,7x,i3,1x,a2,2x,i6,6xf7.3,10x,f9.2)
  220 format (34x,'total fission',15x,i6,6x,f7.3,10x,f9.2)
  230 format (62x,38('-'))
  240 format (52x,'total',5x,i6,6x,f7.3,10x,f9.2)
  250 format (//'1******  angular distribution results ******'/'********
     1************************************'//)
  260 format (///' ***** spin alignment perpendicular to recoil axis
     1   - standard compound nucleus angular distribution'//)
  270 format (///'****** spin alignment perpendicular to reaction plane
     1   - angular distribution is around z axis perpendicular to'/'
     2***** reaction plane'//)
  280 format (1x,132('-'))
  290 format (30x,'energy and angular distribution of residual nucleus -
     1 z=',i3,' n=',i3)
  300 format (30x,'energy and angular distribution of all residual nucle
     1i')
      end
      subroutine compnd (iz,in,energy,maxc)
c
c     original julian subroutine
c
c     ***********
c     form j distribution of initial compound nuclei
      common /xqshl/ noshl, alev(100), ajcn(150)
      z=iz
      a=iz+in
      sr=mod(iz+in,2)*.5
      call gcldp (a,z,0,fl,al,sig,f1,d,aconst,sr,1.)
      uex=energy-d
      sig=1./(2.*sig*sqrt(uex/al))
      if (maxc.gt.100) maxc=100
      do 10 j=1,maxc
      aj=j-1+sr
   10 alev(j)=(2.*aj+1.)*exp(-(aj+.5)**2*sig)
      return
      end
      subroutine compos (iz,in,eexcn,vzc,maxc,sigcn,ecloss,input)
c
c     modified julian subroutine
c     form initial partial wave distribution
c
c     ***********
      common /xqshl/ noshl, sigjc(100), ajcn(150)
      dimension tlcn(150), csch(8)
      data pi /3.14159265/
      if (input.gt.1) go to 350
c   amp     is a of projectile
c   azp     is z of projectile
c   sp      is spin of projectile
c   st      is spin of target
c   qcn     is q of compound nucleus formation
c   ipot    is optical model potential code
c
c     **************************
      read 410, izp,iap,izt,iat,sp,st,qcn
      if (qcn.ne.0.) go to 10
      inp=iap-izp
      int=iat-izt
      izcm=izp+izt
      incm=inp+int
      ntemp=noshl
      noshl=0
      call masses (izp,inp,dfcp)
      call masses (izt,int,dfct)
      noshl=ntemp
c     modification 1/31/83  g.s. of projectile and target
c                           are always shell corrected.
      call masses (izcm,incm,dfccm)
      qcn=dfcp+dfct-dfccm
   10 ident=0
      if (izp.eq.izt.and.iap.eq.iat) ident=1
      ispin=0
      if (sp.ne.0..or.st.ne.0.) ispin=1
      if (ident.eq.1.and.ispin.eq.1) print 400
c     **************************
      read 420, elab,expsig,agraz,eloss,lminn
      elb1=elab-.5*eloss
      ixper=expsig
      if (expsig.eq.0.) call bass (iap,izp,iat,izt,elb1,expsig)
      if (expsig.lt.1.) expsig=1.
c     agraz is diffuseness of transmission coefficients.
c     if agraz.ne.0  tlom is bypasssed
c     expsig is correct experimental xsection if exists. it is used
c     to cut the tl distribution to provide the correct xsection.
      azp=izp
      amp=iap
      azt=izt
      amt=iat
      iz=izp+izt
      inp=iap-izp
      int=iat-izt
      in=inp+int
      ia=iz+in
      amc=ia
      ecm=elab*amt/(amp+amt)
      amu=amp*amt/(amp+amt)
      alamda=4.651/sqrt(amu*ecm)
      ecloss=eloss*amt/(amp+amt)
      call barfit (iz,ia,0,brr,segs,selmax)
      sigma=31.4159*alamda**2*(selmax+1.)**2
c     selmax is ang mom at which barrier goes to zero
c     sigma is maximum cross section for which statistical model is
c     valid. higher partial waves have zero fission barrier
      sigtmp=expsig
      ilimit=0
      if (expsig.le.sigma.or.selmax.lt.10.) go to 20
      expsig=sigma
      ilimit=1
c     ilimit=1 implies cross section  w a s  truncated
c     ilimit=0 implies cross section not modified
   20 ec1=ecm-ecloss*.5
      eexcn=qcn+ecm
      rot=34.548*amc**(-1.6667)
      lyrast=sqrt(eexcn/rot)+.5
c     lyrast is spin of yrast (rigid body) level at c.n. excitation
c     energy.
      if (agraz.le..3) agraz=0.3
      ll=100
      l1=100
      do 30 i=1,100
      ex=i-40
      ex=ex/agraz
      if (ex.gt.20.) tlcn(i)=0.
      if (ex.lt.-20.) tlcn(i)=1.
   30 if (ex.ge.-20..and.ex.le.20.) tlcn(i)=1./(1.+exp(ex))
   40 sigom=0.
      am=amp*amt/(amp+amt)
      ak=0.218898*sqrt(am*ec1)
      plsq=10.0*pi/ak**2
      do 60 i=1,l1
      fac=1.
      ai1=i-1
      if (ident.eq.0) go to 50
      jmod=mod(i,2)
      fac=0.
      if (jmod.eq.1) fac=2.
   50 sigp=(2.*ai1+1.)*tlcn(i)*fac
   60 sigom=sigom+sigp
      sigom=sigom*plsq
c
c   normalization of optical model and experimental cross sections
c
c   t = 1/(1+exp((l-l(g))/v)
c   m*l+b = y where y = ln(1/t-1), m = 1/v, and b = -m * l(g)
c         solve for l(g) and v from calculated transmission coefficients
c   s = c * sum((2 * l + 1) * t  where l(g) in t is now l(c)
c   (l(g)+1)**2/(l(c)+1)**2 = sigma(calculated)/sigma(experimental)
c         is used to solve for l(c)
c
      exps=expsig/plsq
      tlls=0.
      tls=0.
      ijump=0
      do 90 i=2,l1
      if (tlcn(i).gt..9*tlcn(1)) go to 90
      if (tlcn(i).ge..004) go to 70
      l2=i-1
      go to 100
   70 if (ijump.eq.1) go to 80
      ijump=1
      l0=i-1
   80 ai=i-1
      th=alog(1./tlcn(i)-1.)
      tls=tls+th
      tlls=tlls+th*ai
   90 continue
      l2=l1
  100 ls=(l2-1)*l2/2
  110 l2s=ls*(2*l2-1)/3
      if (l0.le.1) go to 120
      l0s=(l0-1)*l0/2
      ls=ls-l0s
      l2s=l2s-l0s*(2*l0-1)/3
  120 sl=ls
      sl2=l2s
      al=l2-l0
      den=sl*sl-al*sl2
      elam=den/(tls*sl-al*tlls)
      elg=elam*(tls*sl2-sl*tlls)/den
  130 elgom=elg
      iter=0
  140 iter=iter+1
      if (iter.eq.101) go to 190
      sf=0.
      do 170 i=1,150
      ai=i-1
      expon=(ai-elg)/elam
      if (expon.gt.40.) go to 180
      if (expon.gt.-40.) go to 150
      den=1.
      go to 160
  150 expl=exp(expon)
      den=1.+expl
  160 el21=2.*ai+1.
      func=el21/den
      if (func.lt..001) go to 180
  170 sf=sf+func
  180 rsig=exps/sf
      if (abs(rsig-1.).lt..0001) go to 190
      elg=(elg+1.)*sqrt(rsig)-1.
      go to 140
  190 continue
      do 200 i=1,150
      ai=i-1
      expon=(ai-elg)/elam
      if (expon.lt.-60.) expon=-60.
      if (expon.gt.60.) expon=60.
c     60 depends on computer in use.
      tlcn(i)=1./(1.+exp(expon))
      if (tlcn(i).lt..001) go to 210
  200 continue
      l1=150
      go to 220
  210 l1=i-1
  220 sigcn=0.
      if(l1.gt.100)l1=100
  230 schx=sp+st
      schn=abs(sp-st)
      sin=(2.0*sp+1.0)*(2.0*st+1.0)
      isc=sp+st+0.01
      asc=isc
      sc=sp+st-asc
      do 310 i=1,l1
      ac1=i-1
      ajc=ac1+sc
      sumtl=0.
      almax=ajc+schx
      lmax=almax+0.01
      if (schx.lt.ajc) go to 240
      if (schn.gt.ajc) go to 250
      lmin=0
      go to 260
  240 almin=ajc-schx
      lmin=almin+0.01
      go to 260
  250 almin=schn-ajc
      lmin=almin+0.01
  260 l=lmin
  270 if (l.gt.ll) go to 290
      if (l.gt.l1) go to 290
      dl=abs(ajc-l)
      sl=ajc+l
      s1=amax1(schn,dl)
      s2=amin1(schx,sl)
      ds=s2-s1+1.0
      if (ds.lt.1.0) go to 280
      sumtl=sumtl+ds*tlcn(l+1)
  280 l=l+1
      if (l.gt.lmax) go to 290
      go to 270
  290 ajcn(i)=ajc
      sigjc(i)=(2.0*ajc+1.0)*plsq*sumtl/sin
      if (ident.eq.0) go to 300
      jmod=mod(i,2)
      if (jmod.eq.1) sigjc(i)=sigjc(i)*2.
      if (jmod.eq.0) sigjc(i)=0.
c      partial xsection (is  alev  in main program)
  300 sigcn=sigcn+sigjc(i)
  310 continue
      erec=elab*amp/amc
      vzc=sqrt(2.0*erec/(931.5*amc))
      vzcc=vzc*30.
      vzb=sqrt(2.0*elab/(931.5*amp))
      vzbc=vzb*30.
      maxc=l1
      print 450
      print 460
      print 470
      print 480, izp,inp,iap,sp
      print 490, izt,int,iat,st
      print 500, iz,in,ia
      print 510, elab
      print 520, ecm
      print 530, eexcn
      if (ecloss.gt.0.) print 540, ecloss
      print 550, erec
      print 560, vzcc
      print 570, vzc
      print 580, vzbc
      print 590, vzb
      if (expsig.eq.0.) go to 320
      if (lminn.gt.0) print 600, lminn
      expsig=sigtmp
      if (ixper.gt.0) print 610, expsig
      if (ixper.eq.0) print 620, expsig
      if (ilimit.eq.1) print 390, sigma
      print 630, elg,elam
      print 640, lyrast
      print 650, sigcn
  320 sc=mod(ia,2)*0.5
      if (lminn.eq.0) go to 340
      do 330 ll=1,lminn
c    ang. mom. lminn less 1 is index lminn
      sigcn=sigcn-sigjc(ll)
  330 sigjc(ll)=0.
      print 660, sigcn
  340 go to 380
  350 read 670, iz,ia,eexcn,erec,ajnuc,lminn
c      ajnuc is real true spin. (not the index)
      in=ia-iz
      amc=ia
      vzc=sqrt(2.0*erec/(931.5*amc))
      vzcc=vzc*30.
      sc=mod(ia,2)*0.5
      amaxc=ajnuc-sc+1.1
      if(amaxc.lt.1)amaxc=1
      maxc=amaxc
      if (input.eq.3) read 680, maxc,(sigjc(i),i=1,maxc)
      if (input.eq.4) call compnd (iz,in,eexcn,maxc)
      print 690
      print 470
      print 500, iz,in,ia
      print 530, eexcn
      print 550, erec
      print 560, vzcc
      print 570, vzc
      if (input.gt.2) go to 360
      print 700, ajnuc
      sigcn=100.
      go to 380
  360 sigcn=0.
      lm1=lminn+1
      do 370 i=1,maxc
      if (input.eq.5) sigjc(i)=0.
      ac1=i-1
      ajcn(i)=ac1+sc
      if (input.eq.5.and.i.ge.lm1) sigjc(i)=2.*ajcn(i)+1.
      sigcn=sigcn+sigjc(i)
  370 continue
      print 650, sigcn
  380 return
c
  390 format (/35x,'cross section truncated for calculation to',8x,f10.3
     1 /29x,'***** higher partial waves have zero fission barrier'//)
  400 format (///20(1h*),' warning-non zero spins. tlom needs correction
     1'' for identical particles here'///)
  410 format (4i5,2f5.1,2f10.5)
  420 format (2f10.5,5x,2f5.0,i5)
  430 format (///20(1h*),' max. spin probably too small (maxc,ajnuc)
     1  '/' or reaction xsection too big . change agraz ***',20(1h*))
  440 format (///' ******** warning. input channel l above 100 '//' ****
     1**** optical model l-diffuseness ignored.'//' ******** value set i
     2nternally to ',f5.2//' **************************************'//)
  450 format (1x,132('-'))
  460 format (1x,'starting conditions')
  470 format (/,74x,'z',9x,'n',9x,'a',10x,'spin')
  480 format (35x,'projectile',20x,3(7x,i3),7x,f5.1)
  490 format (35x,'target',24x,3(7x,i3),7x,f5.1)
  500 format (35x,'compound nucleus',14x,3(7x,i3))
  510 format (/,35x,'bombarding energy (mev)',27x,f10.3)
  520 format (35x,'center of mass energy (mev)',23x,f10.3)
  530 format (/,35x,'compound nucleus excitation energy (mev)',10x,f10.3
     1 )
  540 format (35x,'excitation energy loss thru target (mev)',9x,f10.2)
  550 format (35x,'compound nucleus recoil energy (mev)',14x,f10.3)
  560 format (35x,'compound nucleus recoil velocity (cm/nanosec)',3x,e12
     1 .5)
  570 format (35x,'compound nucleus velocity/c',21x,e12.5)
  580 format (35x,'beam velocity (cm/nanosec)',22x,e12.5)
  590 format (35x,'beam velocity/c',33x,e12.5)
  600 format (35x,'partial waves below l = ',i2,' excluded')
  610 format (35x,'experimental fusion cross section (mb)',12x,f10.3)
  620 format (35x,'bass 1977    fusion cross section (mb)',12x,f10.3)
  630 format (35x,'fusion l-grazing and l-diffuseness',6x,2f10.3)
  640 format (35x,'yrast spin at maximum excitation energy',11x,i10)
  650 format (35x,'compound nucleus formation cross section (mb)',5x,f10
     1 .3//)
  660 format (35x,'comp. nuc. cross section in window (mb)',11x,f10.3)
  670 format (2i5,3f10.5,i5)
  680 format (i5,5x,7f10.5/(8f10.5))
  690 format (1x,' starting compound nucleus'/1x,26(1h*))
  700 format (/,35x,'compound nucleus spin',34x,f5.1)
      end
      subroutine chnprb (ipot,ixr,isc,sc,se,rlev,probm,maxj,maxjs,ebin
     1 ,amr,ixmin,mebin,ixpr)
c
c     modified julian subroutine
c
c     ***********
      common /gam/ fge1, fgm1, fge2, fgm2, ergc, discpr, iadef
      common /mass/ flmass, izpart(4), inpart(4), izr(4), inr(4), nfile,
     1 barfac, aratio, fbarr(110), lfiss
      common /srch/ prob(2222), iprob, mj(2222), emax(4), emin(4), be(4)
     1 , imodl, ge1, ge2
      dimension maxj(286), maxjs(286), ebin(286), rlev(19696), maxjt(286
     1 ), maxjst(286), ebint(286), wue(3), wum(3), tl(30), tlsum(30),
     2 tls0(25)
      equivalence (tls0(2),tlsum(1))
      common /rlvv/ rlevt(19696)
      real*8 prob,rlevt,rlev,probm,pro
c
c   weisskopf units for gamma deexcitation
c     calculate decay width for specific a, z, ex, j. search is limited
c     to specific final states to save time.
c
      data wue, wum /6.8e-08,4.9e-14,2.3e-20,2.1e-08,1.5e-14,6.8e-21/,
     1 izprev, inprev /0,0/, tls0 /25*0./
      data pi2, tothrd /6.28318,.666667/
      sr=abs(sc-se)
c      isc = spin of level
c      se = particle spin
c      sc  = basic spin unit  ( 0 or 1/2 )
c      sr  = minimum spin we can decay to
      ixrsm=ixr
c      starting from highest level ixrsm , we go left (igo=-1)
c      or right to look for maximum probability.
      icut=0
      if (ixmin.lt.ixpr) icut=1
      jsc=isc
      jmax=isc
      modul=(ipot-1)*19696
      ajc=isc-1+sc
      sep01=se+.01
      sem01=se-.01
      srm1=sr-1.
   10 if (ipot.eq.4) go to 80
      iprobo=iprob
      ecm=ergc+be(ipot)-ebin(ixrsm)
      izz=izr(ipot)
      inn=inr(ipot)
      call tlll (ecm,ipot,izz,inn,tl,llmax)
      if (tl(1).lt.1.e-10) go to 70
      tlsum(1)=tl(1)
      llm=llmax+1
      do 20 i=2,llm
   20 tlsum(i)=tlsum(i-1)+tl(i)
      igo=-1
c      jsr = j of residual nucleus.
      jsr=jsc
      if (jsr.gt.maxj(ixrsm)) jsr=maxj(ixrsm)
      mdelj=maxjs(ixrsm)-maxj(ixrsm)
   30 ajr=jsr+srm1
      sumtl=0.
      dj=abs(ajr-ajc)
      sj=ajr+ajc
      almax=sj+sep01
c      almax = sum of spins
c     lmax=almax+0.01
      lmax=almax
      if (lmax.gt.llmax) lmax=llmax
c      minimum l contributing to transitions
      almin=dj-sem01
      lmin=almin
      if (lmin.gt.lmax) go to 50
c  ++++++
c  ++++++ follwing is substitute code for channel spin
c  ++++++ of  0  and  1/2 .
      lm1=lmin
      lm2=lmax+1
      sumtl=tlsum(lm2)-tlsum(lm1)
c     for lm1=0 tlsum is 0 through equivalence to tls0
      sumtl=sumtl+2.*se*(sumtl-tl(lm2)-tl(lm1+1))
c
c     sumtl is being multiplied by  1 for se=0  all l-s
c                                   2 for se=1/2 l above lmin below lmax
c                                   1 for se=1/2 l=lmin or lmax
c  ++++++ end of substitute code.
      mjs=mdelj+jsr
c      ordinal no. of this energy level
cpar  pro=rlev(mjs)*sumtl*.5
c     ignore parity factor.does not effect results
      pro=rlev(mjs)*sumtl
      if (pro.le.discpr*probm) go to 50
c      see remark about discpr in main prog.
      iprob=iprob+1
c      iprob is total no. of probabilities calculated
c     if(iprob.gt.2222)go to 60
      prob(iprob)=pro
      mj(iprob)=mjs+modul
      if (pro.lt.probm) go to 40
      probm=pro
      jmax=jsr
c      search in all ex and j directions for prob.ge..001probmax
   40 if (jsr.eq.1) go to 60
      jsr=jsr+igo
      if (jsr.gt.maxj(ixrsm)) go to 70
      go to 30
   50 if (igo.eq.1) go to 70
      if (jsr.eq.jsc) go to 40
   60 igo=1
      jsr=jsc+1
      if (jsr.le.maxj(ixrsm)) go to 30
   70 if (ixrsm.eq.ixmin) return
      ixrsm=ixrsm-1
      if (iprob.eq.iprobo.and.icut.eq.0.and.tl(1).ge.1.e-10) return
      jsc=jmax
      go to 10
c     calculate gamma decay probability
   80 amrt=amr**tothrd
      amrtt=amrt**2
      pamrt=pi2*amrt
      ge1=fge1*pamrt*wue(1)
      gm1=fgm1*wum(1)*pi2
      ge2=fge2*pi2*amrtt*wue(2)
      gm2=fgm2*pamrt*wum(2)
      ecm=ergc-ebin(ixrsm)
      ecube=ecm**3
      efive=ecube*ecm*ecm
      gge1=ge1*ecube
      ggm1=gm1*ecube
      gge2=ge2*efive
      ggm2=gm2*efive
cpar  proba=(gge1+ggm1)*.5
cpar  probb=(gge2+ggm2)*.5
      proba=gge1+ggm1
      probb=gge2+ggm2
c     ignore parity
      probc=proba+probb
      l=jsc-3
      do 140 i=1,5
      jr=i+l
      if (jr.le.0) go to 140
      if (jr.gt.maxj(ixrsm)) go to 140
      if (i.eq.1.or.i.eq.5) go to 100
      if (i.eq.2.or.i.eq.4) go to 110
      if (ajc.lt..1) go to 140
      if (ajc.lt..6) go to 120
   90 prob1=probc
      go to 130
  100 prob1=probb
      go to 130
  110 ajr=jr-1+sr
      if (ajc.lt..1.and.ajr.lt..1) go to 120
      go to 90
  120 prob1=proba
  130 mjs=maxjs(ixrsm)+jr-maxj(ixrsm)
      pro=rlev(mjs)*prob1
      if (pro.le.discpr*probm) go to 140
      iprob=iprob+1
      if (iprob.gt.2222) go to 150
      mj(iprob)=mjs+modul
      prob(iprob)=pro
      if (pro.gt.probm) probm=pro
  140 continue
      if (ixrsm.eq.ixmin) go to 160
      ixrsm=ixrsm-1
      jsc=jmax
      go to 80
  150 print 200
      stop
  160 barr=fbarr(isc)
c     calculate fission probability from bohr-wheeler
      lfiss=0
      if (ergc.ge.barr+ebin(1)) go to 161
      iprob=iprob+1
      if(iprob.gt.2222)go to 150
      mj(iprob)=imodl
      prob(iprob)=0.
c     bug fix 12/31/85
      return
  161 if (izprev.eq.izr(4).and.inprev.eq.inr(4)) go to 170
c     dont need to recalcultae level density if same nucleus
      izprev=izr(4)
      inprev=inr(4)
      ixmax=ergc+2.5
      ix=1
      mx=0
      call levden (izr(4),inr(4),ix,ixmax,maxjt,maxjst,ebint,mebint
     1 ,rlevt,mx,be(4),aratio)
  170 exres=ergc-barr
      call find (exres,ixrsm,maxjt,maxjst,ebint,mebint,1)
      sum=0.
      do 180 i=1,ixrsm
      n=ixrsm-i+1
      if (isc.gt.maxjt(n)) go to 180
      nlev=maxjst(n)-maxjt(n)+isc
      r=rlevt(nlev)
      sum=sum+r
      if (r.lt..001*sum) go to 190
  180 continue
cp204 pro=sum*.5
  190 pro=sum
c     ignore parity
      lfiss=1
      iprob=iprob+1
      if (iprob.gt.2222) go to 150
      mj(iprob)=imodl
      prob(iprob)=pro
      return
c
  200 format ('0dimension of prob not high enough - ****************
     1 ')
      end
      subroutine moment (i,a,ap,mode,ep,ncasc,id)
c
c     original julian subroutine
c
c     ***********
c      monitoring of outgoing particles.
      common /mom/ vy(9996), vz(9996), vx(9996)
c      ep parm. comes in from vzc when id=1
      dimension eres(9996), fct(9996)
      common /xqang/ sum(300,10), mdir, costh
      equivalence (eres(1),vy(1)), (fct(1),vz(1))
      if (id-2) 10,30,40
c
c   initialization
c
   10 n=1
      nn=n+ncasc-1
      do 20 j=n,nn
      vy(j)=ep*mdir
      vx(j)=0.
   20 vz(j)=ep*(1-mdir)
      go to 60
c
c   calculation of moment
c
   30 rn4=ranf(0.)
      phi=6.28318*rn4
      sox=2.*ap*ep/(a**2*931.5)
c     bug fix  a*(ap+a) replaced by a**2 12/15/82
c     if(sox.lt.0.)print 943,ap,ep,a
      vt=sqrt(sox)
      sox=1.-costh**2
      if (sox.lt.0.) sox=0.
      sinth=sqrt(sox)
      vzse=vt*costh
      vyse=vt*sinth*sin(phi)
      vxse=vt*sinth*cos(phi)
      vz(i)=vz(i)+vzse
      vy(i)=vy(i)+vyse
      vx(i)=vx(i)+vxse
      vzp=vz(i)-vzse*(a/ap+1.)
      vyp=vy(i)-vyse*(a/ap+1.)
      vxp=vx(i)-vxse*(a/ap+1.)
      vpp2=vzp*vzp+vyp*vyp+vxp*vxp
      epart=ap*vpp2*469.
      sox=0.
      if (vpp2.gt.0.0) sox=vzp/sqrt(vpp2)
c      do not use quick functions here
      ang=acos(sox)*180./3.1415927
      call outem (2,mode,epart,ang)
      go to 60
c
c   end calculation
c
   40 vfts=vx(i)**2+vy(i)**2+vz(i)**2
      eres(i)=0.5*a*vfts*931.5
      if (vfts.ne.0.) go to 50
      fct(i)=0.0
      go to 60
   50 fct(i)=acos(vz(i)/sqrt(vfts))*180./3.1415927
   60 return
c
      end
      subroutine outem (ictl,mode,eeml,aeml)
c
c     original julian subroutine
c     bin particle angular distribution
c
c     ***********
c
c   ictl is entry control
c   mode is kind of particle emitted or gamma
c   eeml is laboratory energy of emitted particle or gamma
c   aeml is laboratory angle of emitted particle
c   nt is number of events of each mode
c   n(mode,energy,angle) is distribution of events
c       31st energy is for energy overflow
c       32nd energy is for total at each angle
c       there are 18 angles
c   ng(energy) is for each energy number of gammas
c   nw(mode,bin,angle) is distribution of events into four energy bins
c   ngw(bin) is distribution of gammas into four energy bins
c   ew is energy of separation between four energy bins
c   ew1 is lower energy of each of four energy bins
c   ew2 is upper energy of each of three energy bins
c                       the fourth is infinity
c
      common /em/ nt(4), n(3,32,18), ng(32), ew(3), ew1(4), ew2(3), nw(3
     1 ,4,18), ngw(4), ngb(32), ngbw(4), ngt(32), ngtw(32)
      common /xqsig/ sigma, ncasc, dsig
      dimension dsig(18)
      data iwe, dele /1,1./
      go to (10,50,130,50), ictl
c   dele is the energy difference for the 30 energy bins
c   iwe is control  0 is no four energy bins  1 is yes
c
c   initiation
c
c     **************************
   10 continue
c      read10,dele,iwe,(ew(i),i=1,3)
c     **************************
      ew(1)=5.
      ew(2)=10.
      ew(3)=20.
      do 20 i=1,4
   20 nt(i)=0
      do 30 k=1,32
      ng(k)=0
      ngb(k)=0
      ngt(k)=0
      do 30 i=1,3
      do 30 l=1,18
   30 n(i,k,l)=0
      if (iwe.eq.0) return
      do 40 kw=1,4
      ngw(kw)=0
      ngbw(kw)=0
      ngtw(kw)=0
      do 40 i=1,3
      do 40 l=1,18
   40 nw(i,kw,l)=0
      ew1(1)=0.
      ew2(1)=ew(1)
      ew1(2)=ew(1)
      ew2(2)=ew(2)
      ew1(3)=ew(2)
      ew2(3)=ew(3)
      ew1(4)=ew(3)
      return
c
c   monitoring
c
   50 nt(mode)=nt(mode)+1
      k=eeml/dele
      kgm=eeml+1.0001
      if (kgm.gt.30) kgm=31
c      energy bins
      k=k+1
      if (k.gt.30) k=31
      if (mode.eq.4) go to 60
      l=aeml/10.
      l=l+1
c      angle bins
      n(mode,k,l)=n(mode,k,l)+1
      n(mode,32,l)=n(mode,32,l)+1
      go to 80
   60 ngt(kgm)=ngt(kgm)+1
      ngt(32)=ngt(32)+1
      if (ictl.eq.4) go to 70
      ng(kgm)=ng(kgm)+1
      ng(32)=ng(32)+1
      go to 80
   70 ngb(kgm)=ngb(kgm)+1
      ngb(32)=ngb(32)+1
   80 if (iwe.eq.0) return
      if (mode.eq.4) return
      if (eeml.ge.ew(1)) go to 90
      kw=1
      go to 120
   90 if (eeml.ge.ew(2)) go to 100
      kw=2
      go to 120
  100 if (eeml.ge.ew(3)) go to 110
      kw=3
      go to 120
  110 kw=4
  120 nw(mode,kw,l)=nw(mode,kw,l)+1
      return
c
c   graphing
c
  130 do 210 mode=1,3
      print 230
      if (mode.eq.3) go to 150
      if (mode.eq.2) go to 140
      print 240, nt(1)
      go to 160
  140 if (nt(2).eq.0) go to 210
      print 250, nt(2)
      go to 160
  150 if (nt(3).eq.0) go to 210
      print 260, nt(3)
  160 print 230
      print 270
      print 280
      print 290
      print 300
      print 310
      do 180 k=1,30
      ak=k
      e1=(ak-1.)*dele
      e2=e1+dele
      isum=0
      do 170 iq=1,18
  170 isum=isum+n(mode,k,iq)
      if (isum.gt.0) print 320, e1,e2,(n(mode,k,l),l=1,18)
  180 continue
      e1=30.*dele
      print 330, e1,(n(mode,31,l),l=1,18)
      print 310
      e1=0.
      print 350, (n(mode,32,l),l=1,18)
      do 190 i=1,18
      tet=(i-.5)*.1745
      st=sin(tet)
  190 dsig(i)=n(mode,32,i)/(ncasc*st*.1745*6.2832)*sigma
      print 340, dsig
      print 360
      if (iwe.eq.0) go to 210
      do 200 kw=1,3
      print 320, ew1(kw),ew2(kw),(nw(mode,kw,l),l=1,18)
  200 continue
      print 330, ew1(4),(nw(mode,4,l),l=1,18)
      print 230
  210 continue
      mode=4
      print 230
      print 370, nt(4)
      print 230
      print 380
      print 390
      print 230
      print 400
      print 410
      do 220 k=1,30
      e1=k-1
      e2=e1+1.
      if (ngt(k).gt.0) print 420, e1,e2,ng(k),ngb(k),ngt(k)
  220 continue
      e1=30.*dele
      if (ngt(31).gt.0) print 430, e1,ng(31),ngb(31),ngt(31)
      print 410
      e1=0.
      print 440, ng(32),ngb(32),ngt(32)
      print 450
      return
c
  230 format (1x,132('-'))
  240 format (37x,'neutron spectra in laboratory coordinates (',i6,1x,'e
     1vents )')
  250 format (38x,'proton spectra in laboratory coordinates (',i6,1x,'ev
     1ents )')
  260 format (38x,'alpha spectra in laboratory coordinates (',i6,1x,'eve
     1nts )')
  270 format (/,1x,'energy range [',45x,'angular range (deg)')
  280 format (5x,'(mev)',4x,'[     0    10    20    30    40    50    60
     1    70    80    90   100   110   120   130   140   150   160   170
     2')
  290 format (14x,'[',18(5x,'['))
  300 format (14x,'[    10    20    30    40    50    60    70    80
     190   100   110   120   130   140   150   160   170   180')
  310 format (1x,13('-'),'[',109('-'))
  320 format (1x,f5.1,1x,'-',f5.1,1x,'[',18i6)
  330 format (1x,'above',2x,f5.1,1x,'[',18i6)
  340 format (1x,'dsig/domeg',3x,'[',18f6.1)
  350 format (1x,'total',8x,'[',18i6)
  360 format (1x,123('-'))
  370 format (48x,'gamma ray spectrum   (',i6,1x,'events )')
  380 format (1x,'emission from unbound and bound states(*), and total g
     1amma ray spectrum')
  390 format (1x,'(*) note that emission of a particle from an unbound s
     1tate is not allowed in the code if ecm is less than emin')
  400 format (26x,'energy range (mev)   [',13x,'unbound',15x,'bound',15x
     1 ,'total')
  410 format (26x,21('-'),'[',60('-'))
  420 format (28x,f5.1,'  - ',f5.1,'     [',3(14x,i6))
  430 format (28x,'above    ',f5.1,'     [',3(14x,i6))
  440 format (28x,'total              [',3(14x,i6))
  450 format (26x,82('-'))
      end
      subroutine outres (ictrl,lp,np)
c
c     original julian subroutine
c     bin residue angular distribution
c
c     ***********
c
c   ictrl is control parameter
c   eresl is kinetic energy of residual nucleus in laboratory
c   aresl is laboratory angle of residual nucleus
c
      common /xqsig/ sigma, ncasc, dsig
      common /xqang/ sum(300,10), mdir, ct
      common /out/ knuc(9996,6), nr(16,32,37), tl(94)
      common /srch/ nrw(16,4,37)
      common /xqres/ iz, in, ex, vz
      common /mom/ eres(9996), ares(9996), vx(9996)
      dimension ewr(3), ewr1(4), ewr2(3), dsig(36)
      dimension angle1(36), angle2(36)
c
c   nr(energy,angle) is distribution of residual nuclei
c   nrw(bin,angle) is distribution of residual nuclei in 4 energy bins
c   ewr  are separation energies between four energy bins
c   ewr1 low energies of four energy bins
c   ewr2 high energies of four energy bins. fourth is infinity
c
      if (ictrl-2) 10,60,160
c
c   initialization
c
c     **************************
   10 continue
c      read4,elow,dele,delang,iwr,(ewr(i),i=1,3)
      erec=460.*(iz+in)*vz**2
      ilow=(1.-.4*ex/(iz+in))*erec+1.001
      elow=ilow
      if (elow.le.0.) elow=0.
      idele=(erec-elow)/16.+.5
      dele=idele
      if (dele.lt.1.) dele=1.
      delang=1.+4.*mdir
      iwr=1
      ewr(1)=elow
      irec=erec
      ewr(2)=irec
      if (ewr(2).lt.ewr(1)+1.) ewr(2)=ewr(1)+1.
      ewr(3)=2.*ewr(2)-elow
      if (delang.eq.0..or.delang.gt.5.) delang=5.
c
c   energy interval of 1st of 30 bins is from 0. to elow
c   dele is energy interval in other 29 bins
c      31st bin for overflow
c      32nd bin is for total
c   delang is angle interval in 18 angle bins
c   iwr is control - 0 is no four bins, 1 is yes
c
      do 20 k=1,32
      do 20 i=1,16
      do 20 l=1,37
   20 nr(i,k,l)=0
      do 30 m=1,36
      am=m
      angle1(m)=(am-1.)*delang
   30 angle2(m)=am*delang
      if (iwr.eq.0) go to 50
      do 40 i=1,16
      do 40 kw=1,4
      do 40 l=1,37
   40 nrw(i,kw,l)=0
      ewr1(1)=0.
      ewr2(1)=ewr(1)
      ewr1(2)=ewr(1)
      ewr2(2)=ewr(2)
      ewr1(3)=ewr(2)
      ewr2(3)=ewr(3)
      ewr1(4)=ewr(3)
   50 return
c
c   monitoring
c
   60 if (eres(lp).ge.elow) go to 70
      k=1
      go to 80
   70 k=(eres(lp)-elow)/dele
      k=k+2
      if (k.le.30) go to 80
      k=31
   80 continue
      l=ares(lp)/delang
      l=l+1
      if (l.le.36) go to 90
      l=37
   90 continue
      if (np.gt.15) go to 100
      nr(np,k,l)=nr(np,k,l)+1
      nr(np,32,l)=nr(np,32,l)+1
  100 nr(16,k,l)=nr(16,k,l)+1
      nr(16,32,l)=nr(16,32,l)+1
      if (iwr.eq.0) go to 50
      if (eres(lp).gt.ewr(1)) go to 110
      kw=1
      go to 140
  110 if (eres(lp).gt.ewr(2)) go to 120
      kw=2
      go to 140
  120 if (eres(lp).gt.ewr(3)) go to 130
      kw=3
      go to 140
  130 kw=4
  140 if (np.gt.15) go to 150
      nrw(np,kw,l)=nrw(np,kw,l)+1
  150 nrw(16,kw,l)=nrw(16,kw,l)+1
      go to 50
c
c   graphing
c
  160 continue
      print 250
      print 260, (angle1(m),m=1,18)
      print 270
      print 280, (angle2(m),m=1,18)
      print 290
      if (elow.lt.0.01) go to 170
      e2=elow
      print 300, e2,(nr(np,1,l),l=1,18)
  170 do 180 k=2,30
      ak=k
      e1=(ak-2.)*dele+elow
      e2=e1+dele
      if (isum(nr,np,k,1,18).gt.0) print 310, e1,e2,(nr(np,k,l),l=1,18)
  180 continue
      print 320, e2,(nr(np,31,l),l=1,18)
      print 290
      fac=sigma/(6.2832*ncasc)
      do 190 i=1,36
      tet=(i-.5)*.01745
  190 dsig(i)=fac*nr(np,32,i)/(sin(tet)*.01745)
      print 340, (nr(np,32,l),l=1,18)
      print 330, (dsig(l),l=1,18)
      print 350
      if (iwr.eq.0) go to 50
      do 200 i=1,3
      print 310, ewr1(i),ewr2(i),(nrw(np,i,l),l=1,18)
  200 continue
      print 320, ewr1(4),(nrw(np,4,l),l=1,18)
      print 360
      is=0
      do 210 ij=1,32
      do 210 lq=19,37
  210 is=is+nr(np,ij,lq)
      if (is.eq.0) return
      print 250
      print 370, (angle1(m),m=19,36)
      print 270
      print 380, (angle2(m),m=19,36),angle2(36)
      print 390
      if (elow.lt.0.01) go to 220
      e2=elow
      print 400, e2,(nr(np,1,l),l=19,37)
  220 do 230 k=2,30
      ak=k
      e1=(ak-2.)*dele+elow
      e2=e1+dele
      if (isum(nr,np,k,19,37).gt.0) print 410, e1,e2,(nr(np,k,l),l=19,37
     1 )
  230 continue
      print 420, e2,(nr(np,31,l),l=19,37)
      print 390
      print 430, (nr(np,32,l),l=19,37)
      print 330, (dsig(l),l=19,36)
      print 360
      if (iwr.eq.0) go to 50
      do 240 i=1,3
      print 410, ewr1(i),ewr2(i),(nrw(np,i,l),l=19,37)
  240 continue
      print 420, ewr1(4),(nrw(np,4,l),l=19,37)
      print 360
      go to 50
c
  250 format (/,' energy range [',48x,'angular range (deg)')
  260 format (5x,'(mev)    [',18(1x,f5.1),' [')
  270 format (14x,'[',18(4x,'[',1x),' [')
  280 format (14x,'[',18(1x,f5.1),' [')
  290 format (1x,13('-'),'[',109('-'),'[')
  300 format (1x,'below',2x,f5.1,' [',18i6,' [')
  310 format (1x,f5.1,' -',f5.1,' [',18i6,' [')
  320 format (1x,'above',2x,f5.1,' [',18i6,' [')
  330 format (1x,'dsig/domeg   ','[',18f6.0,'[')
  340 format (1x,'total',8x,'[',18i6,' [')
  350 format (1x,124('-'))
  360 format (1x,130('-'))
  370 format (5x,'(mev)    [',18(1x,f5.1),' [ above')
  380 format (14x,'[',18(1x,f5.1),' [',1x,f5.1)
  390 format (1x,13('-'),'[',109('-'),'[',6('-'))
  400 format (1x,'below',2x,f5.1,' [',18i6,' [',i6)
  410 format (1x,f5.1,' -',f5.1,' [',18i6,' [',i6)
  420 format (1x,'above',2x,f5.1,' [',18i6,' [',i6)
  430 format (1x,'total',8x,'[',18i6,' [',i6)
      end
      subroutine levden (iz,in,ifirst,imax,maxj,maxjs,ebin,mebin,rlev
     1 ,ixpr,be,aratio)
c
c     original julian subroutine
c     generate level density table
c
c     ***********
c      gilbert and cameron default level density subroutine
      dimension maxj(286), maxjs(286), ebin(286), rlev(19696), jbin(286)
      common /gam/ fge1, fgm1, fge2, fgm2, ergc, discpr, iadef, nd,
     1 ergmin, ndl(6), izdl(6), indl(6), jdl(36,6), edl(36,6)
      common /qmmj/ maxc, fyrst, facla, erot(110)
      real*8 rlev
      data iprnt /0/
      kk=0
      iprnt=iprnt+1
      ibin=1
      if (ifirst.lt.1) ifirst=1
      ixpr=0
      if (ifirst.gt.1) go to 30
      ixpr=1
      jbin(1)=1
      ebin(1)=0.
      if (nd.eq.0) go to 30
      do 20 i=1,nd
      if (iz.ne.izdl(i).or.in.ne.indl(i)) go to 20
      ixpr=ndl(i)
      do 10 j=1,ixpr
      jbin(j)=jdl(j,i)
   10 ebin(j)=edl(j,i)
      go to 30
   20 continue
   30 ia=iz+in
      if (ixpr.eq.0) go to 60
      do 50 i=1,ixpr
      max=jbin(i)
      maxj(i)=max
      do 40 j=1,max
      kk=kk+1
      rlev(kk)=0.
c     if (j.eq.max) rlev(kk)=1.e-25
c
c     this is where you need to renormalize the level density for
c     computers that have a limited exponent, as well as in GCLVD
c
      if (j.eq.max) rlev(kk)=1.0000
   40 continue
      maxjs(i)=kk
   50 continue
   60 amr=ia
      azr=iz
      delex=ibin
      del2=delex*.5
      if (ixpr.gt.0) ifirst=(ebin(ixpr)+.5+del2*.5)
      efirst=del2+ifirst
      emax=imax
      sr=mod(ia,2)*.5
      if (iprnt.le.5) write (4,80)
      if (iprnt.le.5) write (4,100) iz,ia,ixpr,ifirst,imax
      call gcldp (amr,azr,iadef,falit,alit,sigsq,fact1,delta,aconst,sr
     1 ,aratio)
      call gclvd (alit,sigsq,fact1,delta,efirst,emax,delex,sr,kk,mebin
     1 ,ixpr,maxj,maxjs,ebin,rlev,jbin)
      if (alit.gt.0.) go to 70
      print 90, iz,ia,ifirst,imax
      stop
   70 maxxl=maxc+9
      if (iprnt.le.5) write (4,80)
      write (4,100) iz,ia,ixpr,ifirst,imax,alit,maxxl,erot(maxxl),delta
     1 ,be
c     if(ixpr.gt.4)return
c     print 494,(ebin(i),jbin(i),i=1,ixpr)
c 494 format((8(1x,f8.3,1x,i3)))
      return
c
   80 format (//' level density generated for following nuclei'//'   iz
     1  ia   ixpr emin emax alit  lmax erot   del     b.e.'/)
   90 format (' level density generation error at z , a, ifirst, imax'4i
     1 5//)
  100 format (5i5,f7.1,i5,f7.2,f6.2,f9.3)
      end
      subroutine gcldp (amr,azr,iadef,falit,alit,sigsq,fact1,delta
     1 ,aconst,sr,aratio)
      dimension pzgc(98), pngc(150), szgc(98), sngc(150)
      common /xqtt/ ux, ex, tt
      common /qmmj/ maxc, fyrst, facla, erot(110)
      common /xqshl/ noshl, alev(100), ajcn(150)
c***  level density parameters were taken from a. gilbert and a.g.w.
c***  cameron, canadian journal of physics, volume 43 (1965) 1446.
c***  data  do not include level density parameters for z or n lt.11
      data pzgc /0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.46,0.,2.09,0.,1.62,0
     1 .,1.62,0.,1.83,0.,1.73,0.,1.35,0.,1.54,0.,1.20,0.,1.06,0.,1.36,0.
     2 ,1.43,0.,1.17,0.,1.24,0.,1.20,0.,1.28,0.,1.28,0.,1.35,0.,1.36,0.,
     3 1.19,0.,1.14,0.,1.12,0.,1.58,0.,1.17,0.,1.18,0.,1.22,0.,0.97,0.,0
     4 .92,0.,0.62,0.,0.68,0.,0.64,0.,0.72,0.,0.75,0.,0.71,0.,0.87,0.,0.
     5 83,0.,0.89,0.,0.79,0.,0.89,0.,0.78,0.,0.69,0.,0.61,0.,0.72,0.,0.7
     6 7/
      data pngc /0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.67,0.,1.80,0.,1.67,0
     1 .,1.86,0.,2.04,0.,1.64,0.,1.44,0.,1.54,0.,1.30,0.,1.27,0.,1.29,0.
     2 ,1.41,0.,1.50,0.,1.50,0.,1.43,0.,1.88,0.,1.47,0.,1.57,0.,1.46,0.,
     3 0.93,0.,0.72,0.,1.12,0.,1.29,0.,0.94,0.,1.24,0.,1.25,0.,1.14,0.,1
     4 .32,0.,1.15,0.,1.24,0.,1.43,0.,1.09,0.,1.20,0.,1.04,0.,0.70,0.,0.
     5 85,0.,0.76,0.,0.92,0.,0.99,0.,1.10,0.,0.92,0.,0.73,0.,0.70,0.,0.8
     6 7,0.,0.61,0.,0.69,0.,0.55,0.,0.40,0.,0.73,0.,0.58,0.,0.86,0.,1.13
     7 ,0.,0.84,0.,0.79,0.,0.82,0.,0.71,0.,0.41,0.,0.38,0.,0.67,0.,0.61,
     8 0.,0.78,0.,0.67,0.,0.67,0.,0.79,0.,0.60,0.,0.57,0.,0.49,0.,0.43,0
     9 .,0.50,0.,0.39/
      data szgc /0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-2.91,-4.17,-5.72,-7.80,-
     1 8.97,-9.70,-10.10,-10.70,-11.38,-12.07,-12.55,-13.24,-13.93,-14.7
     2 1,-15.53,-16.37,-17.36,-18.52,-18.44,-18.19,-17.68,-17.09,-16.65,
     3 -16.66,-16.59,-16.35,-16.18,-16.41,-16.60,-16.54,-16.42,-16.84,-1
     4 7.22,-17.42,-17.52,-17.82,-18.19,-18.58,-19.11,-19.83,-19.14,-18.
     5 35,-17.40,-16.54,-15.68,-14.75,-13.71,-12.87,-12.18,-11.61,-11.09
     6 ,-10.78,-10.53,-10.41,-10.21,-9.85,-9.36,-8.97,-8.56,-8.13,-7.68,
     7 -7.33,-7.11,-7.16,-7.05,-6.81,-6.56,-6.95,-7.52,-8.03,-8.41,-8.86
     8 ,-7.71,-6.38,-5.47,-4.78,-4.37,-4.17,-4.12,-4.29,-4.61,-5.04,-5.4
     9 8,-5.96,-6.40,-6.87,-7.20,-7.74/
      data sngc /0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,6.80,7.53,7.55,7.21,7.44,
     1 8.07,8.94,9.81,10.60,11.39,12.54,13.68,14.34,14.19,13.83,13.50,13
     2 .00,12.13,12.60,13.26,14.13,14.92,15.60,16.38,17.08,17.55,17.98,1
     3 8.33,18.56,18.71,18.65,18.55,18.52,18.34,18.01,17.38,16.56,15.62,
     4 14.38,12.88,13.24,13.71,14.40,15.16,15.89,16.43,16.97,17.59,18.08
     5 ,18.72,19.22,19.51,19.73,19.91,20.06,20.16,20.09,19.83,19.41,19.0
     6 6,18.66,17.73,17.03,16.44,16.00,15.33,14.49,13.42,12.28,11.14,10.
     7 10,9.09,10.00,10.64,11.18,11.70,12.22,12.71,13.05,12.99,12.62,12.
     8 11,11.66,11.21,10.81,10.38,10.03,09.65,09.38,08.99,8.62,8.33,8.10
     9 ,7.82,7.56,7.33,7.15,6.83,6.69,6.55,6.53,6.49,6.39,5.82,5.26,4.53
     $ ,3.83,3.08,2.37,1.72,1.05,0.27,-0.69,-1.69,-2.58,-3.16,-1.72,-0.4
     $ 1,0.71,1.66,2.62,3.22,3.76,4.10,4.46,4.83,5.09,5.18,5.17,5.10,5.0
     $ 5,5.04,5.03,4.99,4.98,5.11,5.27,5.39,5.37,5.30/
c
      izr=azr+0.01
      cnr=amr-azr+0.5
      inr=cnr
      if (facla.lt.0.) go to 10
      alit=amr/facla*aratio
c  gilbert & cameron bypassed for little-a
      go to 30
   10 if (izr.ge.11.and.inr.ge.11) go to 20
      print 60, izr,inr
      aconst=0.
      return
   20 continue
      aconst=0.120
      if (iadef.eq.2) aconst=0.142
      falit=0.00917*(szgc(izr)+sngc(inr))+aconst
      alit=falit*amr*aratio
   30 sigsq=0.0888*alit*amr**0.66667
      fact1=1.0/(12.0*8.0**0.5*alit**0.25*sigsq**1.5)
c      obtain rotation energies from c.p.s.
      maxxl=maxc+9
      if (maxxl.gt.110) maxxl=110
   50 delta=(pzgc(izr)+pngc(inr))*(1-noshl)
      ux=2.5+150./amr
      ex=ux+pzgc(izr)+pngc(inr)
      tt=sqrt(alit/ux)-1.5/ux
      tt=1./tt
      return
c
   60 format (1x,'level density parameters for z=',i3,1x,'and/or n=',i3,
     1 1x,'are not included as data in subroutine gcldp')
      end
      subroutine gclvd (alit,sigsq,fact1,delta,efirst,emax,delex,sr,kk
     1 ,mebin,ixpr,maxj,maxjs,ebin,rlev,jbin)
c     this should not affect other compilers.
      dimension maxj(286), maxjs(286), ebin(286), rlev(19696), jbin(286)
     1 , twosr(110)
      common /qmmj/ maxc, fyrst, facla, erot(110)
      common /xqtt/ ux, ex, tt
      common /gam/ fge1, fgm1, fge2, fgm2, ergc, discpr, iadef
       real*8 rlev,tmpexp
      data iprnt /0/
      iprnt=iprnt+1
      ii=286-ixpr
      maxxl=maxc+9
      if (maxxl.gt.110) maxxl=110
c *****   change 2/7/84 no backbending assumed in level density
c
      eexrst=efirst
      do 70 iex=1,ii
      ie=iex+ixpr
       ebn=ebin0(eexrst)
      eexrs=eexrst+ebn
      eexrst=eexrs
      if (eexrst.lt.ex) eexrs=ex
      if (eexrst.gt.emax) go to 80
      ebin(ie)=eexrst
      uex1=eexrst-delta
      uex=eexrs-delta
      if (uex.gt.0.) go to 10
      maxj(ie)=1
      kk=kk+1
      rlev(kk)=0.
      maxjs(ie)=kk
      go to 70
   10 temp=sqrt(uex/alit)
      fact2=fact1/(uex**1.25*temp**1.5)
      factt=1.
      if (eexrst.le.ex) factt=exp((eexrst-ex)/tt)
      fact2=fact2*factt
      if (maxxl.gt.110) maxxl=110
c    following loop replaced for optimization enhancement
      jmx=maxxl
      do 20 j=1,maxxl
c     calculate Gilber and Cameron EROT if necessary
      if (fyrst.lt.0.) then
        al=j-1+sr
        erot(j)=(al+.5)**2/(2.*sigsq*temp)
c     sigsq*temp is spin cutoff squared in G and C formalism.
      endif
      jj=j
      if (erot(j).gt.uex) go to 30
      if (erot(j).gt.uex1) go to 30
      twosr(j)=2*j-1+2.*sr
   20 continue
      go to 40
   30 jmx=jj-1
      if (jmx.gt.0) go to 40
      kk=kk+1
      rlev(kk)=0
      maxj(ie)=1
      maxjs(ie)=kk
      go to 70
   40 do 50 j=1,jmx
c
c     This is where you need to renormalize level density for computers
c     that have a limited range of exponents in addition to subroutine LEVDEN
c     10**-25 = e**(-57.5647)
c
c     rlev(kk+j)=twosr(j)*fact2*exp(2.*sqrt(alit*(uex-erot(j)))-57.5647)
      tmpexp=2.*sqrt(alit*(uex-erot(j)))
      rlev(kk+j)=twosr(j)*fact2*dexp(tmpexp)
     1 * ebn
   50 continue
c     print 105,uex,alit,fact2,(erot(jj),jj=1,100,10)
c 105 format(2f6.1,e10.1,10f7.3)
      kk=kk+jmx
      if (kk.le.19696) go to 60
      print 130
      alit=0.
      return
   60 maxjs(ie)=kk
      maxj(ie)=jmx
   70 continue
      print 100, eexrst,emax
      alit=0.
      return
   80 mebin=ie-1
      if (iprnt.gt.4) return
      write (4,110)
      do 90 ii=1,40
      j2=maxjs(ii)
      j1=j2-maxj(ii)+1
      jj1=0
      jj2=j2-j1
      write (4,120) ii,jj2,ebin(ii),(rlev(jj),jj=j1,j2)
   90 continue
      return
c
  100 format (' error. level density table max energy = ',f7.1/'   maxim
     1um demanded is ',f7.1/' edit source code increasing all 286-s'//)
  110 format (//' level density table'/' no.     jmax  exc.     level de
     1nsity'//)
  120 format (1x,i4,4x,i4,f6.1,1p10e10.1/((19x,10e10.1)))
  130 format (' ***** level density exceeded 19696 ')
      end
      function ebin0 (x)
      i=x*.25
      if (i-1) 10,20,30
   10 ebin0=.1
      return
   20 ebin0=.3
      return
   30 ebin0=1.
      return
      end
      subroutine masses (iz,in,defec)
c     written by a. gavron
      common /xqshl/ noshl
      dimension dmt(262,15), izz1(262), izz2(262)
      data ae, m /2ha=,0/
      ia=iz+in
      if (noshl.eq.1) go to 40
      if (m.ne.0) go to 30
      open(3,file='mass3.dat',status='old',form='formatted')
      do 20 im=1,262
   10 read (3,70,end=40) a4,am,z1,z2
      if (a4.ne.ae) go to 10
      m=am+.1
      if (m.ne.im) go to 50
      iz1=z1+.01
      iz2=z2+.01
      izz1(im)=iz1
      izz2(im)=iz2
      iz12=iz2-iz1+1
   20 read (3,80) (dmt(im,i),i=1,iz12)
      close(3)
   30 iz2=izz2(ia)
      iz1=izz1(ia)
      if (iz.lt.iz1.or.iz.gt.iz2) go to 40
      defec=dmt(ia,iz-iz1+1)*.001
      return
c  60 print 100,iz,in
   40 continue
      aaa=iz+in
      zzz=iz
      defec=ymass(zzz,aaa)
      if (defec.eq.999.) print 60
      return
   50 print 90, ia,im,m
      stop
c
   60 format (/////60(1h*)//' warning - no physical mass defect availabl
     1e. '//' decay mode to this nucleus  and following chain suppressed
     2'/' edit mass table to correct for missing mass'//60(1h*))
   70 format (a2,f3.0,4x,f3.0,4x,f3.0)
   80 format (8f10.0)
   90 format (' mass sequence error'/' ia im m ',3i6)
      end
      function isum (n,i,j,l1,l2)
      dimension n(16,32,27)
      isum=0
      do 10 l=l1,l2
   10 isum=isum+n(i,j,l)
      return
      end
      subroutine optpot (i,in,iz,rr,ar,ri,ai)
c     written by a. gavron
      common /pot/ att(5), av(4), aw(4), awv(4), aso(4), aw1(4), epp(36)
     1 , irad, imt, ecut, ethrs, riv, aiv, rsot, asot, imesh, rcx, vqs(3
     2 ), vq0(3)
c   - subroutine to provide default optical model parameters
      a=in+iz
      if (i.gt.1) go to 10
c   - neutrons
      av(1)=47.01
      av(2)=-.267
      av(3)=-.0018
c
      rr=1.322-a*7.6e-4+a**2*4.e-6-a**3*8.e-9
      ar=.66
      rcx=0.
      aw(1)=9.52
      aw(2)=-.053
      aw(3)=0.
      ri=1.266-a*3.7e-4+a**2*4.e-6-a**3*4.e-9
      ai=.48
      rmatch=0.
      awv(1)=0.
      awv(2)=0.
      awv(3)=0.
      ri=1.25
      ai=0.47
      return
   10 if (i.gt.2) go to 20
c     protons
c     protons and neutrons from perey and perey
      av(1)=53.3+27*(in-iz)/a+.4*iz/a**.3333
      av(2)=-.55
      av(3)=0.
      rr=1.25
      ar=.65
      rcx=1.25
      aw(1)=13.5
      aw(2)=0.
      aw(3)=0.
      ri=1.25
      ai=0.47
      riv=0.
      aiv=.5
      awv(1)=0.
      awv(2)=0.
      awv(3)=0.
      return
   20 av(1)=50.
c     alphas from igo and huizenga
      av(2)=0.
      av(3)=0.
      rr=1.17*a**.3333+1.77
      ar=.576
      rcx=rr-1.77
      awv(1)=3.+.105*a
      awv(2)=0.
      awv(3)=0.
      riv=rr
      aiv=ar
      ri=0.
      ai=.49
      aw(1)=0.
      aw(2)=0.
      aw(3)=0.
      return
      end
      subroutine range (mode,imin,imax)
c     written by a. gavron
c     this subroutine determines the range in ex of the level density
c     table. (relic from days when memory was scarce...)
      common /xqres/ iz, in, energy, vzc
      common /mass/ flmass, izpart(4), inpart(4), izq(4), inq(4), nfile,
     1 barfac, aratio, fbarr(110), lfiss
      common /dat/ spart(4), pmass(4)
      common /srch/ prob(2222), iprob, mj(2222), qmax(4), qmin(4), be(4)
     1 , imodl
      common /gam/ fge1, fgm1, fge2, fgm2, ergc, discpr, iadef, ndisc,
     1 ergmin
      real*8 prob
      izr=izq(mode)
      inr=inq(mode)
      imm=286
      inn=imm-1
      a=izr+inr
      z=izr
      emax=ergc+be(mode)-qmin(mode)
      emin=ergmin+be(mode)-qmax(mode)-2.*qmax(4)
      if (mode.ne.4) go to 10
      emin=emin-qmax(4)
c     the range for mode=4 allows for 4 gamma decays. if you specify
c     a large gamma decay width, find-errors may result. increase the
c     coefficient of qmax(4) to increase level density range.
      br=fbarr(1)
      br=br+1.
      emin1=ergmin-br-15.
c     this assures enough level density for fission calculation
      if (emin1.lt.emin) emin=emin1
   10 imin=emin-2.5
      if (imin.lt.0) imin=0
      imax=emax+2.5
      if (imax.gt.imm) imax=imm
      if (imin.gt.inn) imin=inn
      if (imax.gt.0) go to 20
      imax=2
   20 return
      end
      subroutine find (e,ixrsm,maxj,maxjs,ebin,mebin,mode)
c     written by a. gavron
c     sub. to find index ixrsm corresponding to nearest energy below e
      dimension maxj(286), maxjs(286), ebin(286)
      common /mass/ flmass, izpart(4), inpart(4), izr(4), inr(4), nfile,
     1 barfac, aratio, fbarr(110), lfiss
      data kk /0/
      ixrsm=1
      if (e.lt.0.) return
      if (e.ge.ebin(1).and.e.le.ebin(mebin)) go to 10
      kk=kk+1
      if (kk.le.10) print 50, e,ebin(1),ebin(mebin),izr(mode),inr(mode)
      if (kk.eq.10) print 40
      ixrsm=mebin
      if (e.lt.ebin(1)) ixrsm=1
c     find limits error means that the energy range of the level
c     density table is too small to accomodate the possible decay
c     modes. either this is an input error (if you are inputting
c     a specific range for a level density ) or subroutine range
c     needs modification. consult additional comment in subroutine
c     range.
      return
   10 ixrsm=e-ebin(1)+1.05
      if (ixrsm.gt.mebin) ixrsm=mebin
   20 if (ebin(ixrsm).gt.e) go to 30
      ixrsm=ixrsm+1
      go to 20
   30 ixrsm=ixrsm-1
      if (ebin(ixrsm).gt.e) go to 30
      return
c
   40 format (20(1h*),' tenth warning. further warnings suppressed'/' be
     1tter look into it  ************************')
   50 format (' find limits error. e=',f7.2,' not between ',2f7.2,' z='
     1 ,i3,'  n=',i3/' consult comments in subroutine range and check in
     2put data for level density ')
      end
      subroutine timend (nc,istop)
c     written by a. gavron
c     if time runs out, call this to see how far you got. (needs to be
c     activated in main program)
      common /out/ izcs(9996), incs(9996), jcs(9996), ergcs(9996), mjcs(
     1 9996), rlev(19696,4)
      real*8 rlev
cibm  call jstime(itime)
c     call gettr (time)
      if (time.gt.6.0) return
      s=0.
      sf=0.
      sz=0.
      sn=0.
      ss=0.
      se=0.
      do 20 i=1,nc
      s=s+1.
      if (izcs(i).lt.0) go to 10
      sz=sz+izcs(i)
      sn=sn+incs(i)
      se=se+ergcs(i)
      ss=ss+jcs(i)
      go to 20
   10 sf=sf+1
   20 continue
      sz=sz/(s-sf)
      sn=sn/(s-sf)
      ss=ss/(s-sf)
      se=se/(s-sf)
      sf=sf/s*100.
      print 30, sz,sn,ss,se,sf
      print 40
      istop=0
      return
c
   30 format (' time end diagnostics'/' average z    n     j    e   %fis
     1sion'/4x,5f6.1)
   40 format (' your job is about to run out of time ********* '/' incre
     1ase time parameter, considering the excitation energy you  are lef
     2t at ')
      end
      subroutine fisrot (a,z,al,bf,barfac)
c     cohen-plasil-swiatecki rotating liquid drop fission barrier
c     provided by courtesy of dr. f. plasil, oak ridge national lab.
c
c     barrier=delsp-delr
      dimension x1b(6,11), x2b(6,11), x3b(10,20)
      data x1b /.28,.243,.221,.208,.195,.18,.211,.186,.17,.1506,.136,.12
     1 ,.152,.131,.1155,.096,.0795,.0625,.09725,.0795,.065,.0506,.0375,.
     2 0253,.05771,.0455,.03414,.0235,.014,.0065,.03325,.0235,.0153,.008
     3 1,.001,.0,.01625,.009,.0032,.0,.0,.0,.0071,.0,.0,.0,.0,.0,.0,.0,.
     4 0,.0,.0,.0,.0,.0,.0,.0,.0,.0,0.,0.,0.,0.,0.,0./
      data x2b /.18,.1695,.1515,.133,.1155,.0949,.1495,.1363,.1165,.099,
     1 .0815,.0594,.12,.1032,.0864,.0678,.0469,.028,.09,.0725,.0556,.037
     2 ,.019,.0057,.0625,.045,.0304,.016,.005,0.,.0406,.0264,.0151,.0052
     3 ,0.,0.,.0253,.0144,.0027,0.,0.,0.,.0141,.006,0.,0.,0.,0.,.0065,.0
     4 008,0.,0.,0.,0.,.002,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
      data x3b /.0949,.0755,.0564,.0382,.0223,.0121,.00588,.00242,.00069
     1 ,.0001,.0873,.0684,.049,.0306,.0162,.0074,.00267,.00055,0.,0.,.08
     2 01,.061,.0418,.0235,.0108,.00373,.00071,0.,0.,0.,.073,.054,.035,.
     3 0178,.0062,.00125,0.,0.,0.,0.,.0661,.047,.0284,.012,.0025,0.,0.,0
     4 .,0.,0.,.0594,.0404,.022,.0065,0.,0.,0.,0.,0.,0.,.0528,.034,.0159
     5 ,.002,0.,0.,0.,0.,0.,0.,.0465,.0277,.01,0.,0.,0.,0.,0.,0.,0.,.040
     6 1,.0217,.0044,0.,0.,0.,0.,0.,0.,0.,.0339,.0158,.00024,0.,0.,0.,0.
     7 ,0.,0.,0.,.028,.0106,0.,0.,0.,0.,0.,0.,0.,0.,.0219,.0064,0.,0.,0.
     8 ,0.,0.,0.,0.,0.,.0164,.0025,0.,0.,0.,0.,0.,0.,0.,0.,.0122,0.,0.,0
     9 .,0.,0.,0.,0.,0.,0.,.0085,0.,0.,0.,0.,0.,0.,0.,0.,0.,.0057,0.,0.,
     $ 0.,0.,0.,0.,0.,0.,0.,.0035,0.,0.,0.,0.,0.,0.,0.,0.,0.,.0016,0.,0.
     $ ,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0
     $ .,0.,0.,0.,0.,0./
      an=a-z
      paren=1.-1.7826*((an-z)/a)**2
      eso=17.9439*paren*a**.66667
      x=0.019655*z*(z/a)/paren
      y=1.9254*al**2/(paren*a**2.3333)
      ix=20.*x+.999
      cx=ix
      bx=20.*x+.999
      dx=bx-cx
      if (x.gt..25) go to 10
      by=10.*y+.999
      if (by.gt.9.) by=9.
      if (by.lt.1.) by=1.
      iy=by
      cy=iy
      dy=by-cy
      b2=(x1b(ix+1,iy+1)-x1b(ix,iy+1))*dx+x1b(ix,iy+1)
      b1=(x1b(ix+1,iy)-x1b(ix,iy))*dx+x1b(ix,iy)
      bf=(b2-b1)*dy+b1
      go to 30
   10 if (x.gt..5) go to 20
      by=20.*y+.999
      if (by.gt.11.) by=11.
      if (by.lt.1.) by=1.
      ix=ix-5
      iy=by
      cy=iy
      dy=by-cy
      b1=(x2b(ix+1,iy)-x2b(ix,iy))*dx+x2b(ix,iy)
      b2=(x2b(ix+1,iy+1)-x2b(ix,iy+1))*dx+x2b(ix,iy+1)
      bf=(b2-b1)*dy+b1
      go to 30
   20 if (x.gt..95) x=.95
      ix=20.*x+.999
      ix=ix-10
      by=100.*y+.999
      if (by.gt.19.) by=19.
      if (by.lt.1.) by=1.
      iy=by
      cy=iy
      dy=by-cy
      b1=(x3b(ix+1,iy)-x3b(ix,iy))*dx+x3b(ix,iy)
      b2=(x3b(ix+1,iy+1)-x3b(ix,iy+1))*dx+x3b(ix,iy+1)
      bf=(b2-b1)*dy+b1
   30 bf=bf*barfac*eso
      return
      end
      subroutine barfit (iz,ia,il,bfis,segs,selmax)
c
csrk This subroutine returns the barrier height bfis, the ground-state
csrk energy segs, in MeV, and the angular momentum at which the fission
csrk barrier disappears,  Lmax,  in units of h-bar,
csrk when called with integer arguments iz, the atomic number,
csrk ia, the atomic mass number, and il, the angular momentum in units
csrk of h-bar, (Planck's constant divided by 2*pi).
c
csrk      The fission barrier for  il = 0  is calculated from a 7th orde
csrk fit in two variables to 638 calculated fission barriers for z value
csrk from 20 to 110.  These  638 barriers are fit with an rms deviation
csrk 0.10 MeV by this 49-parameter function.
csrk If  barfit  is called with (iz,ia) values outside the range of the
csrk the barrier height is set to 0.0, and a message is
csrk printed on the default output file.
c
csrk      For il values not equal to zero, the values of
csrk L at which the barrier is  80%  and  20%  of the L=0 value are
csrk respectively fit to 20-parameter functions of  Z  and  A, over a mo
csrk restricted range of  A  values, than is the case for  L = 0.
csrk The value of L where the barrier disappears, Lmax
csrk is fit to a 24-parameter function of Z and A,  with the same range
csrk Z  and  A  values as  l-80  and  l-20.
csrk      Once again, if an  (iz,ia) pair is outside of the range of
csrk validity of the fit, the barrier value is set to  0.0  and a messag
csrk is printed.  These three values  (Bfis(L=0),L-80, and L-20) and the
csrk constraints  of  Bfis = 0 and  d(Bfis)/dL = 0 at L = Lmax and L = 0
csrk lead to a fifth-order fit to Bfis(L) for L> L-20.  The first three
csrk constraints lead to a third-order fit for the region L < L-20.
c
csrk      The ground-state energies are calculated from a 120-parameter
csrk fit in Z, A, and L to 214 ground-state energies for 36 different
csrk Z  and  A  values.
csrk (The range of Z and A is the same as for L-80, L-20, and L-max)
c
csrk      The calculated barriers from which the fits were
csrk made were calculated in 1983-1984 by A. J. Sierk of Los Alamos
csrk National Laboratory   Group T-9, using  Yukawa-plus-exponential dou
csrk folded nuclear energy, exact Couloub diffuseness corrections,
csrk and diffuse-matter moments of inertia. The parameters of the model
csrk are those derived by Moller and Nix in 1979:
csrk r-0 = 1.16 fm, as = 21.13 MeV, kappa-s = 2.3  a = 0.68 fm.
csrk The diffuseness of the matter and charge distributions used
csrk corresponds to a surface diffuseness parameter (defined by Myers)
csrk of 0.99 fm.  The calculated barriers for L = 0 are
csrk accurate to a little less than 0.1 MeV;  the output from this
csrk subroutine is a little less accurate.  Worst errors may be as large
csrk as 0.5 MeV; characteristic uncertainty is in the range of 0.1-0.2
csrk MeV.   The rms deviation of the ground-state fit from the 214 input
csrk values is 0.20 MeV.  The maximum error occurs for light nuclei
csrk in the region where the ground state is prolate, and may be greater
csrk than 1.0 MeV for very neutron-deficient nuclei, with L near Lmax.
csrk most nuclei likely to be encountered in real experiments, the maxim
csrk error is closer to 0.5 MeV, again for light nuclei and L near Lmax.
c
csrk      Lmax may be in error by up to 1.0 for light nuclei,
csrk but is probably good to about 0.2-0.3 for heavy nuclei which do
csrk not have prolate ground states  (Z >75 or so).
c
csrk written by A. J. Sierk,  LANL  T-9
csrk Version 1.0   February, 1984
c
csrk the following is NECESSARY for 32-bit machines like DEC VAX, IBM,et
c
      real*8 elzcof,pa,pz,aa,zz,el20,el80,bfis0,elmax,pl,ell
      real*8 egs
c
      dimension elzcof(7,7), elmcof(5,4), emncof(5,4), pa(7), pz(7), pl(
     1 10)
      dimension emxcof(6,4), egscof(5,6,4), egs1(5,6), egs2(5,6), egs3(5
     1 ,6), egs4(5,6)
      equivalence (egs1,egscof), (egs2,egscof(1,1,2)), (egs3,egscof(1,1,
     1 3)), (egs4,egscof(1,1,4))
      data iprnt /0/
      data emncof /-9.01100e+2,-1.40818e+3,2.77000e+3,-7.06695e+2,8.8986
     1 7e+2,1.35355e+4,-2.03847e+4,1.09384e+4,-4.86297e+3,-6.18603e+2,-3
     2 .26367e+3,1.62447e+3,1.36856e+3,1.31731e+3,1.53372e+2,7.48863e+3,
     3 -1.21581e+4,5.50281e+3,-1.33630e+3,5.05367e-02/
      data elmcof /1.84542e+3,-5.64002e+3,5.66730e+3,-3.15150e+3,9.54160
     1 e+2,-2.24577e+3,8.56133e+3,-9.67348e+3,5.81744e+3,-1.86997e+3,2.7
     2 9772e+3,-8.73073e+3,9.19706e+3,-4.91900e+3,1.37283e+3,-3.01866e+1
     3 ,1.41161e+3,-2.85919e+3,2.13016e+3,-6.49072e+2/
      data emxcof /9.43596e4,-2.241997e5,2.223237e5,-1.324408e5,4.68922e
     1 4,-8.83568e3,-1.655827e5,4.062365e5,-4.236128e5,2.66837e5,-9.9324
     2 2e4,1.90644e4,1.705447e5,-4.032e5,3.970312e5,-2.313704e5,7.81147e
     3 4,-1.322775e4,-9.274555e4,2.278093e5,-2.422225e5,1.55431e5,-5.787
     4 42e4,9.97505e3/
      data elzcof /5.11819909e+5,-1.30303186e+6,1.90119870e+6,-1.2062824
     1 2e+6,5.68208488e+5,5.48346483e+4,-2.45883052e+4,-1.13269453e+6,2.
     2 97764590e+6,-4.54326326e+6,3.00464870e+6,-1.44989274e+6,-1.020266
     3 10e+5,6.27959815e+4,1.37543304e+6,-3.65808988e+6,5.47798999e+6,-3
     4 .78109283e+6,1.84131765e+6,1.53669695e+4,-6.96817834e+4,-8.565598
     5 35e+5,2.48872266e+6,-4.07349128e+6,3.12835899e+6,-1.62394090e+6,1
     6 .19797378e+5,4.25737058e+4,3.28723311e+5,-1.09892175e+6,2.0399726
     7 9e+6,-1.77185718e+6,9.96051545e+5,-1.53305699e+5,-1.12982954e+4,4
     8 .15850238e+4,7.29653408e+4,-4.93776346e+5,6.01254680e+5,-4.013082
     9 92e+5,9.65968391e+4,-3.49596027e+3,-1.82751044e+5,3.91386300e+5,-
     $ 3.03639248e+5,1.15782417e+5,-4.24399280e+3,-6.11477247e+3,3.66982
     $ 647e+2/
      data egs1 /1.927813e5,7.666859e5,6.628436e5,1.586504e5,-7.786476e3
     1 ,-4.499687e5,-1.784644e6,-1.546968e6,-4.020658e5,-3.929522e3,4.66
     2 7741e5,1.849838e6,1.641313e6,5.229787e5,5.928137e4,-3.017927e5,-1
     3 .206483e6,-1.124685e6,-4.478641e5,-8.682323e4,1.226517e5,5.015667
     4 e5,5.032605e5,2.404477e5,5.603301e4,-1.752824e4,-7.411621e4,-7.98
     5 9019e4,-4.175486e4,-1.024194e4/
      data egs2 /-6.459162e5,-2.903581e6,-3.048551e6,-1.004411e6,-6.5582
     1 20e4,1.469853e6,6.564615e6,6.843078e6,2.280839e6,1.802023e5,-1.43
     2 5116e6,-6.322470e6,-6.531834e6,-2.298744e6,-2.639612e5,8.665296e5
     3 ,3.769159e6,3.899685e6,1.520520e6,2.498728e5,-3.302885e5,-1.42931
     4 3e6,-1.512075e6,-6.744828e5,-1.398771e5,4.958167e4,2.178202e5,2.4
     5 00617e5,1.167815e5,2.663901e4/
      data egs3 /3.117030e5,1.195474e6,9.036289e5,6.876190e4,-6.814556e4
     1 ,-7.394913e5,-2.826468e6,-2.152757e6,-2.459553e5,1.101414e5,7.918
     2 994e5,3.030439e6,2.412611e6,5.228065e5,8.542465e3,-5.421004e5,-2.
     3 102672e6,-1.813959e6,-6.251700e5,-1.184348e5,2.370771e5,9.459043e
     4 5,9.026235e5,4.116799e5,1.001348e5,-4.227664e4,-1.738756e5,-1.795
     5 906e5,-9.292141e4,-2.397528e4/
      data egs4 /-1.072763e5,-5.973532e5,-6.151814e5,7.371898e4,1.255490
     1 e5,2.298769e5,1.265001e6,1.252798e6,-2.306276e5,-2.845824e5,-2.09
     2 3664e5,-1.100874e6,-1.009313e6,2.705945e5,2.506562e5,1.274613e5,6
     3 .190307e5,5.262822e5,-1.336039e5,-1.115865e5,-5.715764e4,-2.56098
     4 9e5,-2.228781e5,-3.222789e3,1.575670e4,1.189447e4,5.161815e4,4.87
     5 0290e4,1.266808e4,2.069603e3/
c
csrk The program starts here
c
      if (iz.lt.19.or.iz.gt.111) go to 70
      if (iz.gt.102.and.il.gt.0) go to 80
      z=float(iz)
      a=float(ia)
      el=float(il)
      amin=1.2*z+0.01*z*z
      amax=5.8*z-0.024*z*z
      if (a.lt.amin.or.a.gt.amax) go to 90
      aa = dble(2.5d-3*a)
      zz = dble(1.d-2*z)
      ell = dble(1.d-2*el)
      bfis0 = 0.d0
      call lpoly (zz,7,pz)
      call lpoly (aa,7,pa)
      do 10 i=1,7
      do 10 j=1,7
      bfis0=bfis0+elzcof(j,i)*pz(j)*pa(i)
   10 continue
      egs = 0.d0
      segs = sngl(egs)
      bfis = sngl(bfis0)
      amin2=1.4*z+0.009*z*z
      amax2=20.+3.0*z
      if ((a.lt.amin2-5..or.a.gt.amax2+10.).and.il.gt.0) go to 100
      call lpoly (zz,5,pz)
      call lpoly (aa,4,pa)
      el80 = 0.d0
      el20 = 0.d0
      elmax = 0.d0
      do 20 i=1,4
      do 20 j=1,5
        el80 = el80 + dble(elmcof(j,i))*pz(j)*pa(i)
        el20 = el20 + dble(emncof(j,i))*pz(j)*pa(i)
   20 continue
      sel80 = sngl(el80)
      sel20 = sngl(el20)
      call lpoly (zz,6,pz)
      call lpoly (ell,9,pl)
      do 30 i=1,4
      do 30 j=1,6
        elmax = elmax + dble(emxcof(j,i))*pz(j)*pa(i)
   30 continue
      selmax = sngl(elmax)
      if (il.lt.1) return
      x=sel20/selmax
      y=sel80/selmax
      if (el.gt.sel20) go to 40
      q=0.2/(sel20**2*sel80**2*(sel20-sel80))
      qa=q*(4.*sel80**3-sel20**3)
      qb=-q*(4.*sel80**2-sel20**2)
      bfis=bfis*(1.+qa*el**2+qb*el**3)
      go to 50
   40 aj=(-20.*x**5+25.*x**4-4.)*(y-1.)**2*y*y
      ak=(-20.*y**5+25.*y**4-1.)*(x-1.)**2*x*x
      q=0.2/((y-x)*((1.-x)*(1.-y)*x*y)**2)
      qa=q*(aj*y-ak*x)
      qb=-q*(aj*(2.*y+1.)-ak*(2.*x+1.))
      z=el/selmax
      a1=4.*z**5-5.*z**4+1.
      a2=qa*(2.*z+1.)
      bfis=bfis*(a1+(z-1.)*(a2+qb*z)*z*z*(z-1.))
   50 if (bfis.le.0.0) bfis=0.0
      if (el.gt.selmax) bfis=0.0
c
csrk Now calculate rotating ground-state energy
c
      if (el.gt.selmax) go to 110
      do 60 k=1,4
      do 60 l=1,6
      do 60 m=1,5
            egs=egs+dble(egscof(m,l,k))*pz(l)*pa(k)*pl(2*m-1)
   60 continue
      segs=sngl(egs)
      if (segs.lt.0.0) segs=0.0
      return
   70 if(iprnt.eq.0)print 130
      go to 90
   80 if(iprnt.eq.0)print 140
      go to 90
   90 if(iprnt.eq.0)print 150, ia
      aa=ia
      az=iz
      al=il
      call fisrot (aa,az,al,bfis,1.)
      go to 110
  100 if(iprnt.eq.0)print 160, ia,il
  110 aa=ia
      az=iz
      al=il
      segs=yrast(aa,az,al)
      if (iprnt.eq.0) print 120,ia,iz,il
      iprnt=1
c
  120 format (35x,'c.p.s. values substituted for out-of-range sierk valu
     1es'/35x,'a=',i4,4x,'z=',i3,4x,'l=',i3)
  121 format (35x,'c.p.s. values substituted for out-of-range sierk valu
     1e of fission barrier'/35x,'   l=',f6.1,4x,'barr=',f6.1)
  130 format (/10x,'*  *  *  *  barfit called with  z  less than 19 or '
     1 ,' greater than 111.  bfis is set to c.p.s.   *  *  *')
  140 format (/10x,'*  *  *  *  barfit called with  z  greater than 102'
     1 ,' and  L  not equal to zero.  bfis is set to c.p.s.   *  *  *')
  150 format (/10x,'*  *  *  *  barfit called with  a =',i3,', outside '
     1 ,'the allowed values for z = ',i3,' *  *  *  *')
  160 format (/10x,'*  *  *  *  barfit called with  a  =',i3,', outside'
     1 ,' the allowed values for z = ',i3/26x,'for nonzero  L =',i3,'  *
     2  *  *  *')
      end
      subroutine lpoly (x,n,pl)
c
csrk this subroutine calculates the ordinary Legendre Polynomials of
csrk order 0 to n-1 of argument  x  and stores them in the vector
csrk pl.  They are calculated by recursion relation from the first two
csrk polynomials.
c
csrk written by A. J. Sierk   LANL  T-9  February,1984
c
csrk NOTE:  pl and x must be real*8 on 32-bit computers!
c
      real*8 pl,x
c
      dimension pl(20)
      pl(1)=1.0
      pl(2)=x
      do 10 i=3,n
      pl(i)=((2*i-3)*x*pl(i-1)-(i-2)*pl(i-2))/(i-1)
   10 continue
      return
      end
      subroutine plm
c     written by a. gavron
      dimension bin(3), pl(24,24)
c     calculate associated legendre funtions and form table of running
c     sums for later selection of angle.
c     to save space the running sum  sum(n,iang)  is stored for each
c     p(l,m) at  n=l*(l-1)/2+m . the polynomial is of order (l-1,m-1).
c     iang is a cos(angle) index from .05 to .95 in ten steps.
c ******** warning. for  m = l  and l above 20 results are inaccurate.
c ******** real*8 should be changed to real*16 and dsqrt by sqrt in
c ******** fortran h extended, if angular distribution at these spins
c ******** is important.
      real*8 pl,z,temp
      common /xqang/ sum(300,10), mdir, ct
      data bin /.01667,.05000,.08333/
      call faclog
c     this forms factorials for later  c.g.  calculations
      do 10 ik=1,300
      do 10 jk=1,10
   10 sum(ik,jk)=0.
      do 40 iang=1,10
      do 40 iiter=1,3
      z=iang*.1-bin(iiter)
      pl(1,1)=1.
      pl(2,1)=z
      temp=1.-z**2
      pl(2,2)=sqrt(temp)
      pl(3,1)=(3.*z**2-1.)*.5
      pl(3,2)=pl(2,2)*3.*z
      pl(3,3)=3.*temp
      do 20 l=4,24
      ll=l-1
      pl(l,1)=((2*ll-1)*z*pl(l-1,1)-(ll-1)*pl(l-2,1))/ll
      pl(l,2)=((2*ll-1)*z*pl(l-1,2)-ll*pl(l-2,2))/(ll-1)
      do 20 m=3,l
      mm=m-1
      pl(l,m)=(ll+mm-1)*pl(l-1,m-1)-(ll-mm+1)*z*pl(l,m-1)
      pl(l,m)=pl(l,m)/pl(2,2)
c     index of (l,m) = n
   20 continue
c           do 3 l=1,6
c   3 if(iiter.eq.2)print 4,iang,ll,(pl(l,m),m=1,l)
c   4 format(2i10,6f10.4)
      do 30 l=1,24
      do 30 m=1,l
      plmm=pl(l,m)*1.e-16
      n=l*(l-1)/2+m
   30 sum(n,iang)=sum(n,iang)+plmm**2
   40 continue
c
      do 60 n=1,300
      do 50 iang=2,10
   50 sum(n,iang)=sum(n,iang)+sum(n,iang-1)
      do 60 iang=1,10
   60 sum(n,iang)=sum(n,iang)/sum(n,10)
c           do 10 l=1,10
c           do 10 m=1,l
c     n=l*(l-1)/2+m
c  10 print 11,l,m,n,(sum(n,i),i=1,10)
c  11 format(1x,3i4,10f7.4)
      return
      end
      function c3j (fj1,fj2,fj3,fm1,fm2,fm3)
c     obtained from avri fraenkel.
      double precision fa, fact, a1, a2, a3, sign
      common /fac/ fa(203)
      dimension fact(202)
      equivalence (fa(2),fact(1))
      a3=0.
      g1=fj1+fj2
      g2=fj1-fj2
      g3=fj3-g2
      g4=g1-fj3+.01
      g5=g2+fj3
      g6=g1+fj3
      a=fj3-abs(g2)+.01
      if (a.lt.0.or.g4.lt.0.) go to 40
      m1=g4
      m2=g5+.01
      m3=g3+.01
      m4=g6+1.01
      a1=fact(m1)+fact(m2)+fact(m3)-fact(m4)
      h1=fj1+fm1
      h2=fj2+fm2+.01
      h3=fj3+fm3
      p1=fj1-fm1+.01
      p2=fj2-fm2
      p3=fj3-fm3
      p4=g5-p1+.02
      p5=g3-h2+.02
      m1=h1+.01
      m2=h2
      m3=h3+.01
      m4=p1
      m5=p2+.01
      m6=p3+.01
      m7=m1-m5
      mmax=max (m1,m2,m3,m4,m5,m6)
      if (mmax.lt.200) go to 60
      print 50, mmax
   60 a2=(a1+fact(m1)+fact(m2)+fact(m3)+fact(m4)+fact(m5)+fact(m6))/2.d0
      sign=-1.
      nj=g4+1.
      do 10 ii=1,nj
      iz=ii-1
      sign=-sign
      b2=g4-iz
      b3=p1-iz
      if (b3.lt.0.) go to 10
      b4=p4+iz
      if (b4.lt.0.) go to 10
      b5=p5+iz
      if (b5.lt.0.) go to 10
      b6=h2-iz
      if (b6.lt.0.) go to 10
      m1=iz
      m2=b2
      m3=b3
      m4=b4
      m5=b5
      m6=b6
      a3=a3+sign/(dexp(fact(m1)+fact(m2)+fact(m3)+fact(m4)+fact(m5)+fact
     1 (m6)-a2))
   10 continue
      k=m7/2
      a3=a3*sqrt(2.*fj3+1.)
      if (m7-2*k) 20,30,20
   20 c3j=-a3
      go to 40
   30 c3j=a3
   40 continue
      return
c
   50 format (' c3j routine acessed beyond dimension.  mmax = ',i4/'  **
     1**** angular distribution may be in error '///)
      end
      subroutine faclog
      double precision fa, a
      common /fac/ fa(203)
      fa(1)=0.d0
      fa(2)=0.d0
      do 10 i=3,203
      a=i-1
   10 fa(i)=fa(i-1)+dlog(a)
      return
      end
      subroutine mjran (ji,jf,mi,mf,np,ep,iz,in)
c     written by a. gavron
c     this subroutine determines the m-state distribution of the
c     angular momentum i.e., the orientation in space of the compound
c     nucleus.
c     ********* warning -  this is done ignoring the spin of the
c     emitted paricle since it is not judged to have a significant
c     effect on the results.
c     the reason for including m-state distribution is to obtain
c     correct angular distributions for particles and residual nuclei.
c
c     the initial m-state distibution is assumed to be m=0 for all
c     events produced by a compound nucleus reaction.
c     for a fragment produced in a deep inelastic collision, m=j
c     is assumed; the z axis is then taken perpendicular to the
c     reaction plane and the direction of motion of the fragment is
c     the x axis.
c     mjcs(i) in the main program is actual m value - not an index
c             as used for jcs(i)  (which is j+1 or j+0.5)
c             again note that we use integer m's and neglect spin.
      dimension tl(30), tm(200)
      common /xqang/ sum(300,10), mdir, ct
c     sum is table of integrals for chosing angle for a given spherical
c     harmonic
      call tlll (ep,np,iz,in,tl,lmax)
c     from transmission coeff vector we select the l value that was
c     responsible for the transition from  ji  to  jf.
c - lmax is actual l  (not  l+1)
      l1=iabs(ji-jf)
      l2=ji+jf-2
      if (l2.gt.lmax) l2=lmax
      l1=l1+1
      l2=l2+1
      la=l1+1
      ll=l1
      if (la.ge.l2) go to 30
      do 10 i=la,l2
   10 tl(i)=tl(i)+tl(i-1)
      fac=1./tl(l2)
      x=ranf(0.)
      do 20 i=l1,l2
      ll=i
      if (x.le.tl(i)*fac) go to 30
   20 continue
c     m-s are actual values and not indices ( i.e. m+1)
   30 l=ll-1
c     l is now known- we proceed to choose  m.
      ma=jf-1+mi
      if (ma.gt.l) ma=l
      ma=-ma
      mb=jf-1-mi
      if (mb.gt.l) mb=l
      nm=mb-ma+1
      fji=ji-1
      fjf=jf-1
      fl=l
      fmi=mi
      m=1
      if (nm.le.1) go to 70
      do 40 i=1,nm
      m=ma+i-1
      fm=m
      fmf=fmi+fm
   40 tm(i)=c3j(fji,fl,fjf,fmi,fm,fmf)**2
c     m is chosen according to weight of  3-j  (c.g.) symbol squared
      do 50 i=2,nm
   50 tm(i)=tm(i)+tm(i-1)
      fac=1./tm(nm)
      x=ranf(0.)
      do 60 i=1,nm
      m=i
      if (x.lt.tm(i)*fac) go to 70
   60 continue
   70 m=ma+m-1
      mf=mi+m
      m=iabs(m)
      n=ll*(ll-1)/2+m+1
c     find sum table for this l-m combination
      x=ranf(0.)
      do 80 ict=1,10
      ct=ict
      if (x.le.sum(n,ict)) go to 90
c      choose angle cosine
   80 continue
   90 ct=(ct-ranf(0.))*.1
      x=ranf(0.)
      if (x.gt..5) ct=-ct
      return
c
      end
      subroutine track (k,ji,jf,ei,ef,ep,energy,ic,sprob,ecloss)
c     written by a. gavron
c     provide brief summary of decay modes at various ex.
      dimension np(5,3), jd(5,3), jd2(5,3), de(5,3), part(5), e(3),
     1 gammv(3), tauv(3), ntot(3), javr(3), jrms(3)
      real*8 sprob
      common /sch/ maxj(286,4), ebin(286,4), maxjs(286,4), rmass(4),
     1 mebin(4)
      common /out/ izcs(9996), incs(9996), jcs(9996), ergcs(9996), mjcs(
     1 9996), rlev(19696,4)
      real*8 rlev
      data np, jd, jd2 /45*0/, de /15*0./, gammv, tauv /6*0./, part /4hn
     1eut,4hprot,4halph,4hgamm,4hfiss/, ntot /3*0/, javr, jrms /6*0/
      if (ic-2) 10,20,50
   10 e(1)=energy-ecloss-.1
      e(2)=energy*.5+.5
      e(3)=10.
      return
   20 j=0
      if (ei.gt.e(1).and.ef.le.e(1)) j=1
      if (j.eq.0.and.ei.gt.e(2).and.ef.le.e(2)) j=2
      if (j.eq.0.and.ei.gt.e(3).and.ef.le.e(3)) j=3
      if (j*sprob.eq.0.) return
      meb=mebin(4)
      do 30 i=1,meb
      mi=meb-i+1
      if (abs(ebin(mi,4)-energy).lt..501) go to 40
   30 continue
      print 100, ebin(1,4),ebin(meb,4),energy,j
   40 indx=maxjs(mi,4)-maxj(mi,4)+ji
      gamma=sprob/(rlev(indx,4)*6.2832+1.e-36)+1.e-36
      tau=6.6e-22/gamma
      gammv(j)=gammv(j)+gamma
      tauv(j)=tauv(j)+tau
      ntot(j)=ntot(j)+1
      javr(j)=javr(j)+ji
      jrms(j)=jrms(j)+ji**2
      np(k,j)=np(k,j)+1
      jd(k,j)=jd(k,j)+jf-ji
      jd2(k,j)=jd2(k,j)+(jf-ji)**2
      de(k,j)=de(k,j)+ei-ef
      return
   50 print 70, e(1),e(2),e(3)
      do 60 i=1,3
      atot=ntot(i)+1.e-9
      gammv(i)=gammv(i)/atot
      tauv(i)=tauv(i)/atot
      ajavr=javr(i)
      ajrms=jrms(i)
      ajavr=ajavr/atot
      ajrms=sqrt(ajrms/atot-ajavr**2)
      print 80, e(i),gammv(i),tauv(i),ajavr,ajrms
      do 60 k=1,5
      an=np(k,i)
      if (an.le.0.) go to 60
      a=jd(k,i)/(an+1.e-9)
      b=jd2(k,i)
      b=sqrt(b/(an+1.e-9))
      c=de(k,i)/(an+1.e-9)
      if (k.lt.5) print 90, part(k),an,a,b,c
      if (k.eq.5) print 90, part(5),an
   60 continue
      return
c
   70 format ('1 track down of decay modes at ',f5.0,',',f5.0,' and',f5.
     1 0,' mev excitation '/2x,59(1h*)/)
   80 format (///'    ex = ',f5.0,'    gamma = ',1pe11.2,' mev    lifeti
     1me = ',e11.2,' sec '//'    average j = ',0pf5.1,5x,' stand. dev. =
     2 ',f5.1//20x,' part  num    delj   rms-dj  dex '/)
   90 format (21x,a4,f6.0,2f7.1,f8.1)
  100 format (' track error. ld range = ',2f10.4,'  energy = ',f10.4,' j
     1 = ',i5//)
      end
      function ymass (z,a)
c     from program alice (f. plasil & m. blann)
c     lysekil masses.
      common /xqshl/ noshl
      dimension em(10), xk(10), y(2), f(2), emp(10)
      data iprint/0/
      zee=z
      ymass=999.
c     if (z.lt.21.) go to 170
      em(1)=0.0
      em(2)=2.0
      em(3)=8.0
      em(4)=14.0
      em(5)=28.0
      em(6)=50.0
      em(7)=82.0
      em(8)=126.0
      em(9)=184.0
      em(10)=258.0
      cay1=1.15303
      cay2=0.0
      cay3=200.0
      cay4=11.0
      cay5=8.07144
      cay6=7.28899
      gamma=1.7826
      a1=15.4941
      a2=17.9439
      a3=0.7053
      d=0.444
      c=5.8
      smalc=0.325
      do 10 i=1,10
      emp(i)=em(i)**(5.0/3.0)
   10 continue
      do 20 i=1,9
      xk(i)=0.6*(emp(i+1)-emp(i))/(em(i+1)-em(i))
   20 continue
      rz=.863987/a3
      l=0
      un=a-z
      a3rt=a**(1.0/3.0)
      a2rt=sqrt(a)
      a3rt2=a3rt**2.0
      zsq=z**2.0
      sym=((un-z)/a)**2
      acor=1.0-gamma*sym
      parmas=cay5*un+cay6*z
      volnuc=-1.0*a1*acor*a
      sufnuc=a2*acor*a3rt2
      coulmb=a3*zsq/a3rt
      fuzsur=-1.0*cay1*zsq/a
      n=un+.01
      iz=z+.01
      oddev=-1.0*(1.0+2.0*(n/2)-un+2.*(iz/2)-z)/sqrt(a)*cay4
     1 *(1-noshl)
c  correction 6/5/84
     1 *(1-noshl)
c correction 6/4/84
      if (sym.gt.0.4) wterm=0.
      if (sym.gt.0.4) go to 30
      wterm=-1.*cay2*a3rt2*exp(-1.*cay3*sym)
   30 continue
      wotnuc=parmas+coulmb+fuzsur+oddev+wterm
      smass=wotnuc+volnuc+sufnuc
      c2=(sufnuc+wterm)/(a**(2.0/3.0))
      x=coulmb/(2.0*(sufnuc+wterm))
      barr=0.0
      y(1)=un
      y(2)=z
      do 60 j=1,2
      do 40 i=1,9
      if (y(j)-em(i+1)) 50,50,40
   40 continue
      stop
   50 f(j)=xk(i)*(y(j)-em(i))-.6*(y(j)**(5./3.)-emp(i))
   60 continue
      s=(2.0/a)**(2.0/3.0)*(f(1)+f(2))-smalc*a**(1./3.)
      ee=2.*c2*d**2*(1.0-x)
      ff=.42591771*c2*d**3*(1.+2.*x)/a3rt
      sshell=c*s
      v=sshell/ee
      eps=1.5*ff/ee
      if (ee*(1.-3.*v).le.0.0) go to 70
      qcalc=0.0
      theta=0.0
      shll=sshell
      go to 150
   70 to=1.0
   80 do 90 ipq=1,10
      t=to-(1.-eps*to-v*(3.-2.*to**2)*exp(-to**2))/(-eps+v*(10.*to-4.*to
     1 **3)*exp(-to**2))
      if (t.le.0.0) go to 110
      if (abs(t-to).lt.0.0001) go to 100
      to=t
   90 continue
      go to 130
  100 if (2.*ee*(1.-2.*eps*t-v*(3.-12.*t**2+4.*t**4)*exp(-t**2)).gt.0.0)
     1 go to 140
  110 do 120 i=1,20
      to=float(i)/10.
      gl=ee*(1.-eps*to-v*(3.-2.*to**2)*exp(-to**2))
      if (gl.ge.0.0) go to 80
  120 continue
  130 continue
      go to 160
  140 theta=t
      alpha0=d*sqrt(5.)/a**(1./3.)
      alpha=alpha0*theta
      sigma=alpha*(1.+alpha/14.)
      qcalc=.004*z*(rz*a3rt)**2*(exp(2.*sigma)-exp(-sigma))
      shll=ee*t**2-ff*t**3+sshell*(1.-2.*t**2)*exp(-t**2)
  150 cmass=smass+shll*(1-noshl)
      if (noshl.eq.0.and.iprint.eq.0) print 180, a,z,cmass
      iprint=1
      ymass=cmass
  160 continue
      if(z.gt.20)return
  170 if (a.eq.01..and.z.eq.01.) ymass=7.29
      if (a.eq.01..and.z.eq.00.) ymass=8.07
      if (a.eq.02..and.z.eq.01.) ymass=13.14
      if (a.eq.03..and.z.eq.01.) ymass=14.95
      if (a.eq.03..and.z.eq.02.) ymass=14.93
      if (a.eq.04..and.z.eq.02.) ymass=2.43
      if (a.eq.06..and.z.eq.03.) ymass=14.09
      if (a.eq.07..and.z.eq.03.) ymass=14.91
      if (a.eq.09..and.z.eq.04.) ymass=11.35
      if (a.eq.10..and.z.eq.05.) ymass=12.05
      if (a.eq.11..and.z.eq.05.) ymass=8.67
      if (a.eq.12..and.z.eq.06.) ymass=0.
      if (a.eq.13..and.z.eq.06.) ymass=3.13
      if (a.eq.14..and.z.eq.07.) ymass=2.86
      if (a.eq.15..and.z.eq.07.) ymass=0.10
      if (a.eq.16..and.z.eq.08.) ymass=-4.74
      if (a.eq.17..and.z.eq.08.) ymass=-.81
      if (a.eq.18..and.z.eq.08.) ymass=-.78
      if (a.eq.19..and.z.eq.09.) ymass=-3.43
      if (a.eq.20..and.z.eq.10.) ymass=-7.04
      if (a.eq.21..and.z.eq.10.) ymass=-5.73
      if (a.eq.22..and.z.eq.10.) ymass=-8.03
      if (a.eq.23..and.z.eq.11.) ymass=-9.53
      if (a.eq.24..and.z.eq.12.) ymass=-13.93
      if (a.eq.25..and.z.eq.12.) ymass=-13.19
      if (a.eq.26..and.z.eq.12.) ymass=-16.21
      if (a.eq.27..and.z.eq.13.) ymass=-17.2
      if (a.eq.28..and.z.eq.14.) ymass=-21.49
      if (a.eq.29..and.z.eq.14.) ymass=-21.89
      if (a.eq.30..and.z.eq.14.) ymass=-24.44
      if (a.eq.31..and.z.eq.15.) ymass=-24.44
      if (a.eq.32..and.z.eq.16.) ymass=-26.01
      if (a.eq.33..and.z.eq.16.) ymass=-26.58
      if (a.eq.34..and.z.eq.16.) ymass=-29.93
      if (a.eq.36..and.z.eq.16.) ymass=-30.66
      if (a.eq.35..and.z.eq.17.) ymass=-29.01
      if (a.eq.37..and.z.eq.17.) ymass=-31.77
      if (a.eq.36..and.z.eq.18.) ymass=-30.23
      if (a.eq.38..and.z.eq.18.) ymass=-34.72
      if (a.eq.40..and.z.eq.18.) ymass=-35.04
      if (a.eq.39..and.z.eq.19.) ymass=-33.80
      if (a.eq.40..and.z.eq.19.) ymass=-33.53
      if (a.eq.41..and.z.eq.19.) ymass=-35.55
      if (a.eq.40..and.z.eq.20.) ymass=-34.85
      if (a.eq.42..and.z.eq.20.) ymass=-38.54
      if (a.eq.43..and.z.eq.20.) ymass=-38.40
      if (a.eq.44..and.z.eq.20.) ymass=-41.46
      if (a.eq.46..and.z.eq.20.) ymass=-43.14
      if (a.eq.48..and.z.eq.20.) ymass=-44.22
      if(iprint.eq.0)print 180
      iprint=1
      return
c
  180 format(/35x,'lysekil masses substituted for nuclei not in table'/
     1 35x,'a = ',f6.0,'   z = ',f5.0,'   defect = ',f9.2//)
      end
      subroutine bass (ia1,iz1,ia2,iz2,e,sigf0)
c     written by a. gavron
c     bass cross section for fusion
c     print 60
      z1=iz1
      z2=iz2
      a1=ia1
      a2=ia2
      r10=a1**.3333*1.16
      r20=a2**.3333*1.16
      r1=r10-1.6124/r10
      r2=r20-1.6124/r20
      ecm=e*a2/(a1+a2)
      sigf0=0.
      sigf1=0.
      i0=10.*(r1+r2)
      do 40 i=i0,300
      r=.1*i
      s=r-r1-r2
      if (s.le.0.) go to 40
      g=1./(.03*exp(s/3.30)+.0061*exp(s/.65))
      vnuc=-r1*r2*g/(r1+r2)
      vcoul=1.44*z1*z2/r
      v=vnuc+vcoul
      sigf=31.416*r**2*(1.-v/ecm)
      if (sigf0.gt.0.) go to 10
      sigf0=sigf
      ra=r
      va=v
   10 if (sigf.lt.sigf1) go to 20
      go to 30
   20 s=r-r1-r2+.1
      g=1./(.03*exp(s/3.30)+.0061*exp(s/.65))
      vnuc=-r1*r2*g/(r1+r2)
      vcoul=1.44*z1*z2/(r+.1)
      v2=vnuc+vcoul
      sigf2=31.416*(r+.1)**2*(1.-v2/ecm)
      if (sigf2.lt.sigf) go to 30
      sigf0=sigf
      ra=r
      va=v
      go to 50
   30 sigf1=sigf
   40 continue
   50 continue
      if (sigf0.lt.0.) sigf0=0.
c     print 70, e,sigf0,ra,va
      return
c
      end
      function yrast (a,z,al)
c     rotating liquid drop yrast line.
c     ***** courtesy of dr. f. plasil, o.r.n.l.
      dimension x1h(6,11), x2h(6,11), x3h(10,20)
      data x1h /.0,.0,.0,.0,.0,.0,-.0057,-.0058,-.006,-.0061,-.0062,-.00
     1 63,-.0193,-.0203,-.0211,-.022,-.023,-.0245,-.0402,-.0427,-.0456,-
     2 .0497,-.054,-.0616,-.0755,-.0812,-.0899,-.0988,-.109,-.12,-.1273,
     3 -.1356,-.147,-.1592,-.1745,-.1897,-.1755,-.1986,-.2128,-.2296,-.2
     4 51,-.26,-.255,-.271,-.291,-.301,-.327,-.335,-.354,-.36,-.365,-.37
     5 2,-.403,-.42,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
      data x2h /0.,0.,0.,0.,0.,0.,-.0018,-.0019,-.00215,-.0024,-.0025,-.
     1 003,-.0063,-.00705,-.0076,-.0083,-.0091,-.0095,-.015,-.0158,-.016
     2 6,-.0192,-.0217,-.025,-.0245,-.0254,-.029,-.0351,-.0478,-.0613,-.
     3 0387,-.0438,-.0532,-.0622,-.0845,-.0962,-.0616,-.0717,-.0821,-.09
     4 72,-.1123,-.1274,-.0793,-.1014,-.1138,-.1262,-.1394,-.1526,-.12,-
     5 .134,-.1503,-.1666,-.1829,-.1992,-.1528,-.171,-.1907,-.2104,-.230
     6 1,-.2498,0.,0.,0.,0.,0.,0./
      data x3h /0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-.00012,-.00014,-.00016,-.
     1 00018,-.0002,-.00024,-.00029,-.00036,-.00065,-.00089,-.00047,-.00
     2 05,-.00058,-.00065,-.00074,-.00085,-.00101,-.00124,-.00138,-.0017
     3 8,-.001,-.00105,-.00124,-.00138,-.00156,-.00179,-.00275,-.00292,-
     4 .003,-.003,-.00176,-.0019,-.00211,-.00235,-.00263,-.00298,-.00449
     5 ,-.0053,-.0053,-.0053,-.003,-.00308,-.00318,-.00352,-.00392,-.004
     6 17,-.0062,-.0062,-.0062,-.0062,-.00374,-.0041,-.00444,-.00488,-.0
     7 0521,-.00545,-.0066,-.0066,-.0066,-.0066,-.0053,-.0055,-.00585,-.
     8 0064,-.00695,-.007,-.007,-.007,-.007,-.007,-.00632,-.007,-.00742,
     9 -.00792,-.00856,-.009,-.009,-.009,-.009,-.009,-.0079,-.0085,-.010
     $ 22,-.0119,-.012,-.012,-.012,-.012,-.012,-.012,-.00944,-.0102,-.01
     $ 42,-.0182,-.019,-.019,-.019,-.019,-.019,-.019,-.0112,-.0133,-.018
     $ 2,-.0238,-.024,-.024,-.024,-.024,-.024,-.024,-.01303,-.0178,-.022
     $ 6,-.0274,-.028,-.028,-.028,-.028,-.028,-.028,-.0165,-.0254,-.0343
     $ ,-.0343,-.034,-.034,-.034,-.034,-.034,-.034,-.0203,-.033,-.04,-.0
     $ 4,-.04,-.04,-.04,-.04,-.04,-.04,-.025,-.0406,-.046,-.047,-.047,-.
     $ 047,-.047,-.047,-.047,-.047,-.03036,-.0482,-.048,-.048,-.048,-.04
     $ 8,-.048,-.048,-.048,-.048,-.0363,-.0558,-.056,-.056,-.056,-.056,-
     $ .056,-.056,-.056,-.056,-.04234,-.0634,-.064,-.064,-.064,-.064,-.0
     $ 64,-.064,-.064,-.064,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
      an=a-z
      paren=1.-1.7826*((an-z)/a)**2
      eso=17.9439*paren*a**.66667
      x=0.019655*z*(z/a)/paren
      ero=34.548*al**2/a**1.6667
      y=1.9254*al**2/(paren*a**2.3333)
      ix=20.*x+1.
      cx=ix
      bx=20.*x+1.
      dx=bx-cx
      if (x.gt..25) go to 10
      by=10.*y+1.
      if (by.gt.9.) by=9.
      if (by.lt.1.) by=1.
      iy=by
      cy=iy
      dy=by-cy
      h1=(x1h(ix+1,iy)-x1h(ix,iy))*dx+x1h(ix,iy)
      h2=(x1h(ix+1,iy+1)-x1h(ix,iy+1))*dx+x1h(ix,iy+1)
      hf=(h2-h1)*dy+h1
      go to 30
   10 if (x.gt..5) go to 20
      by=20.*y+1.
      if (by.gt.11.) by=11.
      if (by.lt.1.) by=1.
      ix=ix-5
      iy=by
      cy=iy
      dy=by-cy
      h1=(x2h(ix+1,iy)-x2h(ix,iy))*dx+x2h(ix,iy)
      h2=(x2h(ix+1,iy+1)-x2h(ix,iy+1))*dx+x2h(ix,iy+1)
      hf=(h2-h1)*dy+h1
      go to 30
   20 if (x.gt..95) x=.95
      ix=20.*x+1.
      ix=ix-10
      by=100.*y+1.
      if (by.gt.19.) by=19.
      if (by.lt.1.) by=1.
      iy=by
      cy=iy
      dy=by-cy
      h1=(x3h(ix+1,iy)-x3h(ix,iy))*dx+x3h(ix,iy)
      h2=(x3h(ix+1,iy+1)-x3h(ix,iy+1))*dx+x3h(ix,iy+1)
      hf=(h2-h1)*dy+h1
   30 yrast=ero+hf*eso
      return
      end
      subroutine statis
c     written by a. gavron
c     statistical summary of decay modes. user can substitute his.
      dimension ntitl(6), ecm(5), dspin(5), nspc(30,2), pf(30), npf(30),
     1 ecmf(5), dspinf(5), ntotf(5)
      common /qmmj/ maxc, fyrst, facla, erot(110)
      common /xqsig/ sigma, ncasc
      common /xqres/ ir, iq, energy, vzc, itrac, idist
      common /out/ izcs(9996), incs(9996), izq(200), inq(200), ndist(30,
     1 20,6), iebin(30), jbin(20), jsum(16), ntot(5), yjm(30,20), ntfiss
     2 (30,20,6)
      data ntitl /4hneut,4hprot,4halph,4hg-e1,4hg-e2,4hfiss/, nspc /60*0
     1 /
      print 470
      sfe=0.
      tfe=0.
      sfs=0.
      tfs=0.
      gfs=0.
      sfm=0.
      tfm=0.
      do 20 i=1,30
      pf(i)=0.
      npf(i)=0
      do 20 j=1,20
      yjm(i,j)=0.
      do 10 k=1,6
   10 ntfiss(i,j,k)=0
   20 ndist(i,j,6)=0
      do 30 i=1,31
   30 erot(i)=.1*i-.05
      rewind 2
   40 read (2,end=50) indx,mode,jc,jf,mp,ie,ip,izc,inc
      if (mode.lt.4) go to 40
      jx=mode-3
      if (ip.gt.30) ip=30
      nspc(ip,jx)=nspc(ip,jx)+1
      go to 40
   50 print 480, (erot(i),erot(i+1),nspc(i,1),nspc(i,2),i=1,29)
      print 490, erot(31),nspc(30,1),nspc(30,2)
      ib=energy*.03333
      ib=ib+1
      do 60 i=1,30
      if (i.le.20) jbin(i)=-(5*i-1)
   60 iebin(i)=i*ib
      inuc=1
      izq(inuc)=izcs(1)
      inq(inuc)=incs(1)
   70 rewind 2
      do 80 k=1,5
      ecm(k)=0.
      dspin(k)=0.
      ntot(k)=0
      ecmf(k)=0.
      dspinf(k)=0.
      ntotf(k)=0.
      do 80 i=1,30
      do 80 j=1,20
   80 ndist(i,j,k)=0
      if (itrac.eq.0) go to 210
   90 read (2,end=100) indx,mode,jc,jf,mp,ie,ip,izc,inc,fprob
c      the above can provide detailed information on all cascades
c      leading to any specific final nucleus.
c      indx -  izcs(indx) and incs(indx) are the z and n of the final
c       nucleus in the chain of decay.
c      mode - type of particle emitted. 1-n, 2-p, 3-alpha, 4-gamma
c      jc, jf - initial and final spin indices for this particle
c       emission. (each read statement corresponds to 1 emitted part.)
c      mp - the projection of jc on the z axis. ( fractional spins
c       are neglected for the projection. for jc,jf, the actual
c       spin = jc,jf-1   in even mass nucleus
c            = jc,jf -1/2   in odd mass nucleus )
c      ie - excitation energy at emitting level
c      ip - particle energy (c.m.) multiplied by factor of 2
c      izc, inc - z n , of emitting nucleus
c      if izc is negative,it means that the nucleus fissioned
c
      if (izq(inuc).ne.izcs(indx).or.inq(inuc).ne.incs(indx)) go to 90
      ip=(ip+9)/10
c     change scale to 1 mev bins
      gfs=gfs+fprob
      sfe=sfe+ie*fprob
      tfe=tfe+ie**2*fprob
      sfs=sfs+jc*fprob
      tfs=tfs+jc**2*fprob
      sfm=sfm+fprob*mp
      tfm=tfm+fprob*mp**2
      ie1=ie/ib+1
      if (ie1.gt.30) ie1=30
      jc1=jc/5+1
      if (jc1.gt.20) jc1=20
      ndist(ie1,jc1,mode)=ndist(ie1,jc1,mode)+1
      if (mode.gt.5) go to 90
      ntot(mode)=ntot(mode)+1
      ecm(mode)=ecm(mode)+ip
      dspin(mode)=dspin(mode)+jc-jf
      go to 90
  100 print 500
      nzn=0
      do 110 i=1,ncasc
      if (izq(inuc).ne.izcs(i).or.inq(inuc).ne.incs(i)) go to 110
      nzn=nzn+1
  110 continue
      do 160 k=1,5
      if (ntot(k).eq.0) go to 160
      ifiss=0
      izpr=izq(inuc)
      if (izq(inuc).gt.0) go to 120
      ifiss=1
      izpr=-izq(inuc)
  120 am=ntot(k)
      am=am/nzn
      ec=ecm(k)/ntot(k)-.5
      dsp=dspin(k)/ntot(k)
      print 510, izpr,inq(inuc),ntitl(k),ntot(k),nzn,am,ec,dsp
      if (ifiss.eq.1) print 520
      print 90001
90001 format('    ex / spin = ',
     1  '0-4  5-9  10-14 15-19 20-24 25-29 30-34 35-39 40-44 45-49 ',
     1'50-54 55-59 60-64 65-69 70-74 75-79 sum  avrg  stdv '/)
      do 140 i=1,30
      isum=0
      iavr=0
      isgm=0
      do 130 j=1,20
      iavr=iavr+ndist(i,j,k)*(5*j-3)
      isgm=isgm+ndist(i,j,k)*(5*j-3)**2
      isum=isum+ndist(i,j,k)
  130 continue
      if (isum.eq.0) go to 140
      avr=iavr
      avr=avr/isum
      std=isgm
      std=zsqrt(std/isum-avr**2)
      iee1=iebin(i)-ib
      print 530, iee1,iebin(i),(ndist(i,j,k),j=1,16),isum,avr,std
  140 continue
      do 150 j=1,16
      jsum(j)=0
      do 150 i=1,30
  150 jsum(j)=jsum(j)+ndist(i,j,k)
      print 550, jsum
  160 continue
      do 180 i=1,ncasc
      do 170 j=1,inuc
      if (izq(j).eq.izcs(i).and.inq(j).eq.incs(i)) go to 180
  170 continue
      inuc=inuc+1
      izq(inuc)=izcs(i)
      inq(inuc)=incs(i)
      go to 70
  180 continue
      do 190 k=1,5
      ecm(k)=0.
      dspin(k)=0.
      ntot(k)=0
      do 190 i=1,30
      do 190 j=1,20
  190 ndist(i,j,k)=0
      do 200 i=1,30
      do 200 j=1,20
  200 ndist(i,j,6)=0
      rewind 2
  210 read (2,end=230) indx,mode,jc,jf,mp,ie,ip,izc,inc,fprob
      ie1=ie/ib+1
      ip=(ip+9)/10
      gfs=gfs+fprob
      sfe=sfe+ie*fprob
      tfe=tfe+ie**2*fprob
      sfs=sfs+jc*fprob
      tfs=tfs+jc**2*fprob
      sfm=sfm+fprob*mp
      tfm=tfm+fprob*mp**2
c      change scale to 1 mev bins
      if (ie1.gt.30) ie1=30
      jc1=jc/5+1
      if (jc1.gt.20) jc1=20
      yjm(ie1,jc1)=yjm(ie1,jc1)+abs(mp)
      pf(ie1)=pf(ie1)+fprob
      npf(ie1)=npf(ie1)+1
      if (mode.eq.6) go to 220
      if (izcs(indx).lt.0) go to 220
      ndist(ie1,jc1,mode)=ndist(ie1,jc1,mode)+1
c     traceback ndist is residues only
      ntot(mode)=ntot(mode)+1
      ecm(mode)=ecm(mode)+ip
      dspin(mode)=dspin(mode)+jc-jf
      go to 210
  220 ntfiss(ie1,jc1,mode)=ntfiss(ie1,jc1,mode)+1
      if (mode.eq.6) go to 210
      ntotf(mode)=ntotf(mode)+1
      ecmf(mode)=ecmf(mode)+ip
      dspinf(mode)=dspinf(mode)+jc-jf
      go to 210
  230 print 570
      nzn=0
      do 240 i=1,ncasc
      if (izcs(i).ge.0) nzn=nzn+1
  240 continue
      do 260 i=1,30
      do 260 j=1,20
      sum=ntfiss(i,j,6)+1.e-19
      do 250 k=1,5
      sum=sum+ndist(i,j,k)+ntfiss(i,j,k)
  250 continue
      yjm(i,j)=yjm(i,j)/sum
  260 continue
      if (idist.eq.0) go to 280
      print 440
      do 270 i=1,30
      iee1=iebin(i)-ib
      print 450, iee1,iebin(i),(yjm(i,j),j=1,17)
  270 continue
      print 560, (jbin(i),i=1,16)
  280 print 460
      do 320 k=1,5
      if (ntot(k).eq.0) go to 320
      am=ntot(k)
      am=am/nzn
      ec=ecm(k)/ntot(k)-.5
      dsp=dspin(k)/ntot(k)
      print 580, ntitl(k),ntot(k),nzn,am,ec,dsp
      print 90001
      do 300 i=1,30
      isum=0
      iavr=0
      isgm=0
      do 290 j=1,20
      iavr=iavr+ndist(i,j,k)*(5*j-3)
      isgm=isgm+ndist(i,j,k)*(5*j-3)**2
      isum=isum+ndist(i,j,k)
  290 continue
      if (isum.eq.0) go to 300
      avr=iavr
      avr=avr/isum
      std=isgm
      std=zsqrt(std/isum-avr**2)
      iee1=iebin(i)-ib
      print 530, iee1,iebin(i),(ndist(i,j,k),j=1,16),isum,avr,std
  300 continue
      do 310 j=1,16
      jsum(j)=0
      do 310 i=1,30
  310 jsum(j)=jsum(j)+ndist(i,j,k)
      print 550, jsum
  320 continue
      print 590
      sumpf=0.
      do 330 i=1,30
      ni=31-i
      if (npf(ni).le.0) go to 330
      sumpf=sumpf+pf(ni)
      pf(ni)=pf(ni)/npf(ni)
      iee1=iebin(ni)-ib
      print 600, iee1,iebin(ni),pf(ni)
  330 continue
      sumpf=sumpf/ncasc
      print 610, sumpf
      if (sumpf.lt.1.e-6) return
      sfe=sfe/gfs
      tfe=zsqrt(tfe/gfs-sfe**2)*2.36
      sfs=sfs/gfs
      tfs=zsqrt(tfs/gfs-sfs**2)*2.36
      sfm=sfm/gfs
      tfm=sqrt(tfm/gfs)
      print 640, sfe,tfe,sfs,tfs
      if (idist.gt.0) print 650, sfm,tfm
      ntt=0
      do 340 i=1,30
      do 340 j=1,20
      do 340 k=1,6
  340 ndist(i,j,k)=ntfiss(i,j,k)
      do 350 i=1,30
      do 350 j=1,20
  350 ntt=ntt+ndist(i,j,6)
      if (ntt.eq.0) return
      print 620, ntt
      do 420 k=1,6
      if (k.eq.6) go to 360
      if (ntotf(k).eq.0) go to 420
      am=ntotf(k)
      am=am/ntt
      ec=ecmf(k)/ntotf(k)-.5
      dsp=dspinf(k)/ntotf(k)
      print 430, ntitl(k),ntotf(k),am,ec,dsp
      go to 380
  360 print 630, ntitl(k),ntt
      se=0.
      te=0.
      ss=0.
      ts=0.
      do 370 i=1,30
      do 370 j=1,20
      spin=5*j-3
      exen=iebin(i)-.5*ib
      se=se+exen*ndist(i,j,6)
      te=te+exen**2*ndist(i,j,6)
      ss=ss+spin*ndist(i,j,6)
      ts=ts+spin**2*ndist(i,j,6)
  370 continue
      se=se/ntt
      te=zsqrt(te/ntt-se**2)*2.36
      ss=ss/ntt
      ts=zsqrt(ts/ntt-ss**2)*2.36
      print 640, se,te,ss,ts
  380 continue
      print 90001
      do 400 i=1,30
      isum=0
      iavr=0
      isgm=0
      do 390 j=1,20
      iavr=iavr+ndist(i,j,k)*(5*j-3)
      isgm=isgm+ndist(i,j,k)*(5*j-3)**2
  390 isum=isum+ndist(i,j,k)
      if (isum.eq.0) go to 400
      avr=iavr
      avr=avr/isum
      std=isgm
      std=std/isum-avr**2
      std=zsqrt(std)
      iee1=iebin(i)-ib
      print 530, iee1,iebin(i),(ndist(i,j,k),j=1,16),isum,avr,std
  400 continue
      do 410 j=1,16
      jsum(j)=0
      do 410 i=1,30
  410 jsum(j)=jsum(j)+ndist(i,j,k)
      print 550, jsum
  420 continue
      return
c
  430 format (//' pre-fiss decay.  mode = ',a4,4x,'total = ',i4,4x,'mult
     1iplicity =',f5.2,5x,'average ecm =',f6.2,5x,'average spin removed
     2=',f5.1/1x,35(1h*)//)
  440 format ('1 abs value of m states at intermediate j vs ex '//)
  450 format (1x,i3,'-',i3,3x,17f6.1)
  460 format (' evaporation residue decay summary'/' *******************
     1**************'//)
  470 format ('1',/////' ********* complete traceback diagnostic of part
     1icle and gamma emission *********'///)
  480 format (/////20x,'components  of  gamma  spectrum *** '//16x,' ene
     1rgy        e1-spec  e2-spec'//(14x,f5.2,' -',f5.2,5x,i5,5x,i5))
  490 format (14x,'above  ',f5.2,5x,i5,5x,i5)
  500 format (1h1)
  510 format (///' z n =',2i4,'   mode = ',a4,'   total = ',i6,'  out of
     1 ',i4,' events.    multiplicity = ',f5.2,'  average ecm = ',f5.2,'
     2 mev'/1x,13(1h*),3x,'average spin removed = ',f5.1//)
  520 format (' *************  fission follows  particle ','emission cas
     1cades  ************* '//)
  530 format (1x,i3,'-',i3,3x,17i6,2f6.1)
  540 format (/'  ex / j',i9,15i6,'   sum  avrg  stdv'/)
  550 format (' sum',7x,16i6//)
  560 format (/'  ex / j',i9,15i6/)
  570 format (1h1)
  580 format (///' decay to e.r.-s   mode = ',a4,'   total = ',i6,'  out
     1of ',i4,' events.    multiplicity = ',f5.2,'  average ecm = ',f5.2
     2 ,' mev'/1x,13(1h*),3x,'average spin removed = ',f5.1//)
  590 format ('1 fission probability as function of excitation'//'
     1              ex. en.    probability     '/)
  600 format (17x,i3,'-',i3,6x,1pe10.2)
  610 format (' ***** total sum of fission probabilities ',1pe11.2)
  620 format (///'1fission summary. total events = ',i4/' **************
     1**********************'///' paricles preceeding fission'//' ******
     2********************* '//)
  630 format (///' mode = ',a4,4x,' total number = ',i4/1x,35(1h*)//)
  640 format (//' excitation energy window - average = ',f6.1,'
     1  fwhm = ',f5.1//' spin window              - average = ',f6.1,'
     2         fwhm = ',f5.1//)
  650 format (' average  spin projection ',f5.1,'       average rms  spi
     1n projection ',f5.1//)
      end
      subroutine lmnt (iz,mmm)
      dimension ielem(102)
c
      data ielem /2hh ,2hhe,2hli,2hbe,2hb ,2hc ,2hn ,2ho ,2hf ,2hne,2hna
     1 ,2hmg,2hal,2hsi,2hp ,2hs ,2hcl,2har,2hk ,2hca,2hsc,2hti,2hv ,2hcr
     2 ,2hmn,2hfe,2hco,2hni,2hcu,2hzn,2hga,2hge,2has,2hse,2hbr,2hkr,2hrb
     3 ,2hsr,2hy ,2hzr,2hnb,2hmo,2htc,2hru,2hrh,2hpd,2hag,2hcd,2hin,2hsn
     4 ,2hsb,2hte,2hi ,2hxe,2hcs,2hba,2hla,2hce,2hpr,2hnd,2hpm,2hsm,2heu
     5 ,2hgd,2htb,2hdy,2hho,2her,2htm,2hyb,2hlu,2hhf,2hta,2hw ,2hre,2hos
     6 ,2hir,2hpt,2hau,2hhg,2htl,2hpb,2hbi,2hpo,2hat,2hrn,2hfr,2hra,2hac
     7 ,2hth,2hpa,2hu ,2hnp,2hpu,2ham,2hcm,2hbk,2hcf,2hes,2hfm,2hmd,2hno
     8 /
      mmm=ielem(iz)
      return
      end
      function zsqrt (x)
      if (x.lt.0.) x=0.
      zsqrt=sqrt(x)
      return
      end
      function icoupl (sc,se)
      icoupl=-1
      if (se+sc.eq.0.) icoupl=-2
      return
      end
      subroutine tccal (jah,jz,ja,ne,kuc,jmt,btt,siv,bmv,bv,bw,bwv,scx
     1 ,ep,fthrs,fcut,bso,ssot,bsot,irad,jmesh,tp,lmx)
c
c      subroutine tccal package , e.d. arthur, t-2, l.a.n.l. contains
c      following subroutines -   jmtran    rpot    outc    faze2
c                                neuts     cool    zeros   series
c                                int1      int2    int3    int4
c                                prparm
c
c     adapted for vax by a. gavron, 3/28/84
c
c      subroutine package calculates transmission coefficients
c
c     calculates probabilities particle   emission channels
c
      dimension p(80), att(5), av(4), aw(4), lmx(36), ep(36), tp(36,80)
      dimension bv(4), bw(4), bso(4), bwv(4), btt(5)
      common /q55/ aso(4)
      common /q1/ en, ec, iz, ia, is, imdl, si, sii, a, fact, fact1,
     1 ipun
      common /q3/ iah
      common /q10/ pc, pm, tc, tm, iim, imesh, juc, v0, w0, rcx, wvolf,
     1 at(5), rmatch
      common /q30/ aw1(4), w1, ethrs, ecut
      common /q99/ rsot, asot, riv, aiv, awv(4), wvo
      dimension xpart1(6), pc1(6), pm1(6)
      data xpart1 /4hneut,4hprot,4halph,4hdeut,4htrit
     1 ,4hhe-3   /
      data pc1 /0.,1.,2.,1.,1.,2./
      data pm1 /1.00866522,1.007825,4.002635,2.0141,3.01605,3.016035/
c
      iz=jz
      imt=jmt
      riv=siv
      iah=jah
      juc=kuc
      ethrs=fthrs
      imesh=jmesh
      rsot=ssot
      ia=ja
      rcx=scx
      ecut=fcut
      aiv=bmv
      asot=bsot
      do 10 i=1,4
      av(i)=bv(i)
      aw(i)=bw(i)
      aso(i)=bso(i)
      awv(i)=bwv(i)
   10 continue
      do 20 i=1,5
      att(i)=btt(i)
   20 continue
c
c     particle    iah
c     ********    ***
c     neutrons    1
c      protons    2
c     alphas      3
c     deuterons   4
c     tritons     5
c     he3         6
c
c     calculate continuum channel probabilities
      tc=iz
      iim=imt
      iza=za
      tm=ia
      xpart=xpart1(iah)
      write (4,80) xpart
      pc=pc1(iah)
      pm=pm1(iah)
      nep=0
      jup2=2*juc
c
c      ecut eq 0 ignore
c     first set of imag param for imag surface term
      if (ne.eq.0) go to 70
c
c     set up continuum calculation
      nop=jup2
c
c     store zero transmission coefficients for zero incident energy
c
c     set up optical model parameters
      tc=(1.0-pcs)*tc
      do 30 i=1,4
   30 at(i)=att(i)
      call prparm (av,aw)
c
      wvolf=0.
c     store transmission coefficients
      do 60 n=1,ne
c
c     calculate coefficients
      v0=av(1)+av(2)*ep(n)+av(3)*ep(n)**2+av(4)*alog(ep(n))
      w0=aw(1)+aw(2)*ep(n)+aw(3)*ep(n)**2+aw(4)*alog(ep(n))
      att(5)=aso(1)+aso(2)*ep(n)+aso(3)*ep(n)**2+aso(4)*alog(ep(n))
      at(5)=att(5)
      wvo=awv(1)+awv(2)*ep(n)+awv(3)*(ep(n)**2)+awv(4)*alog(ep(n))
      wvo=amax1(wvo,0.)
      v0=amax1(v0,0.)
      w0=amax1(w0,0.)
      if (ecut.gt.0.) w1=aw1(1)+aw1(2)*(ep(n)-ethrs)+aw1(3)*(ep(n)**2)
     1 +aw1(4)*alog(ep(n))
      if (ecut.gt.0..and.ep(n).gt.ecut) w0=amax1(w1,0.)
      vs=at(5)
      vs=amax1(vs,0.)
      call jmtran (p,ep(n),lmax,irad)
      do 40 i=1,nop
      tp(n,i)=0.
      if (i.gt.lmax) go to 40
      tp(n,i)=p(i)
   40 continue
      amu=pm*tm/(pm+tm)
      ecm=ep(n)
      alamda=4.651/sqrt(amu*ecm)
      sum=0.
      do 50 i=1,lmax
      al=i-1
      sum=sum+(2.*al+1)*p(i)
   50 continue
      sum=sum*31.4159*alamda**2
      lmx(n)=lmax
c     write(4,70) ep(n),alamda,sum,lmx(n)
ctc   write (4,80) (tp(n,i),i=1,lmax)
   60 continue
c
c     punch coefficints if requested
c
   70 continue
      return
c
ctc80 format (6e16.5)
c
   80 format (10x,a10' particle emission')
      end
      subroutine jmtran (p,ele,lmaxo,irad)
c
c     calculates spherical transmission coefficients
c     supplied by jim ferguson, 7/7/72
c
c     p - array in which the  transmission coefficients will be returned
c     ele - incident particle energy
c
      dimension p(80), f(44), g(44), fp(44), gp(44)
      common /q99/ rsot, asot, riv, aiv, awv(4), wvo
      common /q10/ pc, pm, tc, tm, iim, imesh, juc, v0, w0, rcx, wvolf,
     1 rs, ar, ri, am, vs, rmatch
      common /q20/ eta(2,80)
      common /jf1/ vr(305,80), wr(305), psi(305), psii(305), psip, psipi
      common /jf2/ h, ekin, atot
      do 10 i=1,80
   10 p(i)=0.
c
      itm=pm+.001
      itm=mod(itm,2)
      izin=itm+1
c
      atot=tm+pm
      cons=pm*tm/atot*0.0478326
      v0=-v0*cons
      w0=-w0*cons
      vx=-vs*cons
      wvo=-wvo*cons
      ekin=ele*cons
      call rpot (vx,irad)
c
      x=imesh
      jmax=2*juc
      lmax=juc+1
      cay=sqrt(ekin)
      rho=x*h*cay
      ee=ekin*atot/(pm*tm*0.0478326)
      etta=.5*pc*tc*2.818*.511*cay/ee
      if (pc.eq.0.) go to 20
      call cool (rho,etta,lmax,f,g,fp,gp)
      go to 40
c
   20 call neuts (rho,lmax,f,g,fp,gp)
      do 30 l=1,lmax
      g(l)=-1.*g(l)
   30 gp(l)=-1.*gp(l)
   40 continue
c     alphas and deuterons yields t-sub-l
c     nskip=1
c     if(izin.eq.1)nskip=2
c  remove two above c-s to get spin dependence,and delete 1 below
      nskip=2
      do 60 j=1,jmax,nskip
      jj=j
c
      call outc (ekin,jj)
      l=j/2+1
c
      call faze2 (psi(imesh),psii(imesh),psip,psipi,g(l),f(l),gp(l),fp(l
     1 ),etar,etai)
      p(j)=1.-etar**2-etai**2
      if (p(j).lt.5.e-14) go to 50
      if(p(j).lt.1.e-4)p(j)=-4.*(f(l)/psi(imesh))**2*psi(imesh)*
     1 (psipi-psip*(psii(imesh)/psi(imesh)))/
     1 (1.+(psii(imesh)/psi(imesh))**2)
c
   50 eta(1,j)=etar
      eta(2,j)=etai
      if (p(j).le.5.0e-14) p(j)=0.
   60 continue
c
      lmaxo=0
      do 70 j=1,80,2
      l=j/2+1
      p(l)=p(j)
      if (p(l).le.0.) go to 80
   70 lmaxo=lmaxo+1
   80 continue
      return
      end
      subroutine rpot (vx,irad)
c
      dimension xtemp1(80), xtemp2(80)
      common /q99/ rsot, asot, riv, aiv, awv(4), wvo
      common /q10/ pc, pm, tc, tm, iim, imesh, juc, v0, w0, rcx, wvolf,
     1 rs, ar, ri, am, vs, rmatch
      common /jf1/ vr(305,80), wr(305), psi(305), psii(305), psip, psipi
      common /jf2/ h, ekin, atot
c
      if (irad.eq.0) go to 10
c     irad=0  implies radii are full and not r0
c     irad=1  implies radii are radius parameters
c     calculate radii
      r=rs*tm**.333333
      rso=rsot*tm**.333333
      rimv=riv*tm**.333333
      rim=ri*tm**.333333
      if (rcx.le.0.) rcx=rs
      rc=rcx*tm**.333333
      go to 20
   10 r=rs
      rso=rsot
      rimv=riv
      rim=ri
      rc=rcx
      if (rc.lt..01) rc=r
c
c     matching radius
   20 x=imesh
      rcmplx=0.
      if (rim.gt.0.) rcmplx=rim+12.*am
      if (rimv.gt.0.) rcmplx=rimv+10.*aiv
      rmatch=amax1(r+10.*ar,rcmplx)
ctc   write(4,909)rmatch,rc
ctc   change ctc to blanks for additional optical model info
ct909 format(' matching radius = ',f7.2,'  rc = ',f7.2)
      h=rmatch/x
c     rmatch could be read in here if necessary
c
      cons2=0.0478326*tm*pm/atot*pc*tc*.511*2.818
      cons=.5/rc
      cons=cons*cons2
      vzer=v0
      wim=w0
      vso=vx
      vol=wvolf
      qq=1.
      jmax=2*juc
      do 30 k=1,jmax
      l=k/2
      j=k-mod(k-1,2)
      xtemp1(k)=(l+1)*l
      xtemp2(k)=(j*(j+2)-4*l*(l+1)-3)*.25
   30 continue
c
      q=0.
      i1=imesh+5
      do 70 i=1,i1
      q=q+h
      es=exp((q-r)/ar)
c     unique rso and aso spin-orbit param
      eso=exp((q-rso)/asot)
      eso=2.*eso/(q*asot*(1.+eso)**2)
      do 40 k=1,jmax
      esp=vzer/(1.+es)+xtemp1(k)/q**2+vso*xtemp2(k)*eso
      vr(i,k)=esp
   40 continue
c
      if (iim.eq.2) go to 50
      esp=(q-rim)/am
      esp1=exp(esp)
      esp=4.*esp1/(1.+esp1)**2
      esp2=(q-rimv)/aiv
      esp3=exp(esp2)
      esp=wim*esp+wvo/(1.+esp3)
      go to 60
c
   50 esp=((q-rim)/am)**2
      if (esp.gt.50.) esp=50.
      esp1=exp(esp)
      esp=wim*(exp(-esp)*(1.-vol)+vol/(1.+esp1))
   60 wr(i)=esp
   70 continue
c
      if (pc.eq.0.) go to 130
      q=0.
      if (pc.eq.2.) go to 90
      do 80 i=1,i1
      q=q+h
      es=cons2/q
      if (q.lt.rc) es=cons*(3.-(q/rc)**2)
      do 80 k=1,jmax
   80 vr(i,k)=es+vr(i,k)
      go to 130
c
c     use the ford-hill coulomb potential for alpha particles
   90 con3=1.627*tm**0.333333
      con32=con3*con3
      do 120 i=1,i1
      q=q+h
      xrad=q/rc
      if (xrad.lt.1.) go to 100
      brak=((rc/q+con3/2.)*exp(con3*(rc-q)/rc))/(con3*con32/3.+2.*con3
     1 +exp(-con3))
      es=(cons2/q)*(1.-xrad*brak)
      go to 110
  100 cx=con3*xrad
      brak1=((1.-exp(cx))/cx+0.5*exp(cx))*(exp(-con3))/con32
      brak2=(1./3.+2./con32+(exp(-con3))/(con3*con32))
      es=(cons2/q)*xrad*(1./con32+0.5-(xrad*xrad)/6.+brak1)/brak2
  110 do 120 k=1,jmax
      vr(i,k)=es+vr(i,k)
  120 continue
  130 return
      end
      subroutine outc (e,j)
      common /q10/ pc, pm, tc, tm, iim, imesh, juc, v0, w0, rcx, wvolf,
     1 rs, ar, ri, am, vs, rmatch
      common /jf1/ vr(305,80), wr(305), psi(305), psii(305), psip, psipi
      common /jf2/ h, ekin, atot
      i1=imesh
      k=j
      i2=i1+2
      ist=j/2-2
      if (ist.lt.1) ist=1
      e2=.0001
      e2i=0.
      c=e2
      ci=0.
      qq=h*h/12.
      qb=16./(21.*h*sqrt(e))
      do 70 i=ist,i2
      zz=(vr(i,k)-e)*qq
      zzi=wr(i)*qq
      c=c+12*(e2*(zz*(1.+zz)-zzi**2)-e2i*zzi*(1.+2*zz))
      ci=ci+12*(e2i*zz*(1.+zz)+zzi*(e2*(1.+2*zz)-e2i*zzi))
      e3=e2+c
      e3i=e2i+ci
      ij=i-i1+3
      bot=(1.-zz)**2+zzi**2
      x=(e2*(1.-zz)-e2i*zzi)/bot
      xi=(e2*zzi+e2i*(1.-zz))/bot
      if (i.lt.i1-2) go to 30
      go to (10,20,30,40,50), ij
   10 psip=-.578125*x+.425*(zz*x-zzi*xi)
      psipi=-.578125*xi+.425*(zz*xi+zzi*x)
      go to 60
   20 psip=.5*x+7.4*(x*zz-xi*zzi)+psip
      psipi=.5*xi+7.4*(zz*xi+zzi*x)+psipi
      go to 60
   30 psi(i)=x
      psii(i)=xi
      go to 60
   40 psip=-.5*x-7.4*(x*zz-xi*zzi)+psip
      psipi=-.5*xi-7.4*(zz*xi+zzi*x)+psipi
      go to 60
   50 psip=.578125*x-.425*(zz*x-zzi*xi)+psip
      psipi=.578125*xi-.425*(zz*xi+zzi*x)+psipi
      psip=qb*psip
      psipi=qb*psipi
   60 continue
      e2=e3
   70 e2i=e3i
c     write(4,150) j,psip,psipi
c 150 format(/5x,'j = ',i6,'  derivatives are psip = ',1pe12.5,
c    ''  and psipi = ',1pe12.5)
c     write(4,155) psi(imesh),psii(imesh)
c 155 format(5x,'at the matching radius psi = ',1pe12.5,
c    ''  and psii = ',1pe12.5/)
      return
      end
      subroutine neuts (bound,lmax,bes,pneum,besp,pneump)
      dimension bes(44), pneum(44), besp(44), pneump(44), w(44)
      si=sin(bound)
      co=cos(bound)
      bes(1)=si/bound
      bes(2)=si/bound**2-co/bound
      pneum(1)=-co/bound
      pneum(2)=-co/bound**2-si/bound
      lll=lmax+1
      do 10 l=3,lll
      el=l-2
      bes(l)=(2.*el+1.)*bes(l-1)/bound-bes(l-2)
   10 pneum(l)=(2.*el+1.)*pneum(l-1)/bound-pneum(l-2)
      do 20 l=1,lmax
      el=l-1
      besp(l)=el*bes(l)/bound-bes(l+1)
      pneump(l)=el*pneum(l)/bound-pneum(l+1)
      besp(l)=bes(l)+besp(l)*bound
   20 pneump(l)=pneum(l)+pneump(l)*bound
      do 30 l=1,lmax
      bes(l)=bound*bes(l)
   30 pneum(l)=bound*pneum(l)
c     write(4,200) lmax,bound
c     do 40 l=1,lmax
c     w(l)=-besp(l)*pneum(l)+bes(l)*pneump(l)
c     write(4,300) bes(l),besp(l),pneum(l),pneump(l),w(l)
c  40 continue
      return
c
      end
      subroutine cool (rho,eta,lmax,f,g,fp,gp)
      dimension g(44), gp(44), f(44), fp(44), fb(90), w(44)
      call zeros (eta,rho,go,gop)
      g(1)=((1./rho+eta)*go-gop)/sqrt(1.+eta**2)
      gp(1)=sqrt(1.+eta**2)*go-(1./rho+eta)*g(1)
      xi=1.
      g(2)=((2.*xi+1.)*(eta+xi*(xi+1.)/rho)*g(1)-(xi+1.)*sqrt(xi**2+eta*
     1 *2)*go)/(xi*sqrt((xi+1.)**2+eta**2))
      l=lmax-1
      do 10 i=2,l
      xi=i-1
      xx=i
      if (i.eq.2) go to 10
      g(i)=((2.*xi+1.)*(eta+xi*(xi+1.)/rho)*g(i-1)-(xi+1.)*sqrt(xi**2
     1 +eta**2)*g(i-2))/(xi*sqrt((xi+1.)**2+eta**2))
   10 gp(i)=(sqrt(xx**2+eta**2)*g(i-1)-((xx**2)/rho+eta)*g(i))/xx
      q=lmax-1
      if (q.lt.rho+10.) q=rho+10.
      imax=q+.0001
      fb(imax+1)=0.
      fb(imax)=1.
      l=imax-1
      do 20 j=1,l
      i=imax-j
      xi=i
   20 fb(i)=((2.*xi+3.)*(eta+(xi+1.)*(xi+2.)/rho)*fb(i+1)-(xi+1.)*sqrt(
     1 (xi+2.)**2+eta**2)*fb(i+2))/((xi+2.)*sqrt((xi+1.)**2+eta**2))
      fbo=(3.*(eta+2./rho)*fb(1)-sqrt(4.+eta**2)*fb(2))/(2.*sqrt(eta**2+
     1 1.))
      alpha=(fbo*g(1)-fb(1)*go)*sqrt(1.+eta*eta)
      do 30 l=2,lmax
      j=lmax+2-l
      xl=j-1
      g(j)=g(j-1)
      gp(j)=gp(j-1)
   30 f(j)=fb(j-1)/alpha
      f(1)=fbo/alpha
      gp(1)=gop
      g(1)=go
c
c     try different recurrence relation
c     write(4,400) lmax,eta,rho
      lmm=lmax-1
      do 40 l=1,lmm
      xl=l
      t1=((xl**2)/rho+eta)*f(l)
      t2=f(l+1)*sqrt(xl**2+eta**2)
      fp(l)=(t1-t2)/xl
      w(l)=fp(l)*g(l)-f(l)*gp(l)
   40 continue
      xlmax=lmax-1
      etalm=eta/xlmax
      fp(lmax)=f(lmax-1)*sqrt(1.+etalm**2)-(xlmax/rho+etalm)*f(lmax)
      w(lmax)=fp(lmax)*g(lmax)-f(lmax)*gp(lmax)
c     write(4,45) (f(l),fp(l),g(l),gp(l),w(l), l=1,lmax)
      return
c
      end
      subroutine zeros (eta,rho,go,gop)
      dimension sigma(400)
      common /bint/ xg(16), ag(16)
      data xg /.9894009349,.9445750230,.8656312023,.7554044083,.61787624
     1 44,.4580167776,.2816035507,.0950125098,-.0950125098,-.2816035507,
     2 -.4580167776,-.6178762444,-.7554044083,-.8656312023,-.9445750230,
     3 -.9894009349/
      data ag /.0271524594,.0622535239,.0951585117,.1246289712,.14959598
     1 88,.1691565193,.1826034150,.1894506104,.1894506104,.1826034150,.1
     2 691565193,.1495959888,.1246289712,.0951585117,.0622535239,.027152
     3 4594/
      if (eta-4.) 20,20,10
   10 rzero=2.*eta
      eta3=exp(alog(eta)/3.)
      sigma(1)=1.223404016*sqrt(eta3)*(1.+.049595701650/(eta*eta3)-.0088
     1 888888889/(eta**2)+.0024551991810/(eta3*eta**3)-.0009108958061/
     2 (eta**4)+.0002534684115/(eta3*eta**5))
      sigma(2)=-.7078817734*(1.-.172826039/(eta3**2)+.0003174603174/(eta
     1 **2)-.00358121485/((eta*eta3)**2)+.000311782468/(eta**4)-.0009073
     2 966427/((eta**4)*eta3**2))/sqrt(eta3)
      call series (sigma,rzero,rho,eta,gp,gpp)
      go=gp
      gop=gpp
      return
   20 rhop=4.
      if (eta.gt.2.) rhop=2.*eta
      up=(22.+2.*eta)/rhop
      down=0.
      call int1 (down,up,rhop,eta,ans1)
      up=12.
      down=0.
      call int2 (down,up,rhop,eta,ans2)
      down=0.
      up=(22.+2.*eta)/rhop
      call int3 (down,up,rhop,eta,ans3)
      up=12.
      down=0.
      call int4 (down,up,rhop,eta,ans4)
      go=ans1-ans2
      gop=-ans3-ans4
      ao=sqrt((1.-exp(-6.283184*eta))/(6.283184*eta))
      sigma(1)=go*ao*rhop
      sigma(2)=ao*(gop*rhop+go)
      call series (sigma,rhop,rho,eta,gp,gpp)
      go=gp
      gop=gpp
      return
      end
      subroutine series (sigma,rzero,rho,eta,an,anp)
      dimension sigma(400)
   10 del=abs(rho-rzero)
      n=0
   20 if (del-.2) 40,40,30
   30 n=n+1
      del=del/5.
      go to 20
   40 del=(rho-rzero)/(5.**n)
      sigma(2)=sigma(2)*del
      sigma(3)=-.5*(del**2)*(1.-2*eta/rzero)*sigma(1)
      sum=sigma(1)+sigma(2)+sigma(3)
      super=sigma(2)/del+2.*sigma(3)/del
      do 70 i=4,200
      x=i-2
      sigma(i)=(del*x*(x-1.)*sigma(i-1)+(rzero-2.*eta)*(del**2)*sigma(i-
     1 2)+(del**3)*sigma(i-3))/(-1.*rzero*x*(x+1.))
      sum=sum+sigma(i)
      super=super+(x+1.)*sigma(i)/del
      if (abs(sigma(i)/sum)-.00000001) 50,50,70
   50 if (abs(x*sigma(i)/(del*super))-.00000001) 60,60,70
   60 imax=i
      go to 80
   70 continue
      if (imax.gt.198) write (4,110)
   80 if (n) 90,100,90
   90 rzero=rzero+del
      sigma(1)=sum
      sigma(2)=super
      go to 10
  100 an=sum
      anp=super
      return
c
  110 format (/'lookout,man')
      end
      subroutine int1 (down,up,rho,eta,sum)
      common /bint/ xg(16), ag(16)
      f(q,eta,rho)=exp(-rho*q+2.*eta*atan(q))
      sum=0.
      del=.2*(up-down)
      do 20 j=1,5
      up=down+del
      do 10 i=1,16
   10 sum=sum+.5*del*ag(i)*f(.5*del*xg(i)+.5*(up+down),eta,rho)
   20 down=down+del
      return
      end
      subroutine int2 (down,up,rho,eta,sum)
      common /bint/ xg(16), ag(16)
      f(q,eta,rho)=(1.-(tanh(q)**2))*sin(rho*tanh(q)-2.*eta*q)
      sum=0.
      del=.2*(up-down)
      do 20 j=1,5
      up=down+del
      do 10 i=1,16
   10 sum=sum+.5*del*ag(i)*f(.5*del*xg(i)+.5*(up+down),eta,rho)
   20 down=down+del
      return
      end
      subroutine int3 (down,up,rho,eta,sum)
      common /bint/ xg(16), ag(16)
      f(q,eta,rho)=q*exp(-rho*q+2.*eta*atan(q))
      sum=0.
      del=.2*(up-down)
      do 20 j=1,5
      up=down+del
      do 10 i=1,16
   10 sum=sum+.5*del*ag(i)*f(.5*del*xg(i)+.5*(up+down),eta,rho)
   20 down=down+del
      return
      end
      subroutine int4 (down,up,rho,eta,sum)
      common /bint/ xg(16), ag(16)
      f(q,eta,rho)=tanh(q)*(1.-(tanh(q)**2))*cos(rho*tanh(q)-2.*eta*q)
      sum=0.
      del=.2*(up-down)
      do 20 j=1,5
      up=down+del
      do 10 i=1,16
   10 sum=sum+.5*del*ag(i)*f(.5*del*xg(i)+.5*(up+down),eta,rho)
   20 down=down+del
      return
      end
      subroutine prparm (av,aw)
c
c     prints spherical optical model parameters
c
c              av - arithmetic expression for the real well depth
c              aw - arithmetic expression for the imaginary well depth
c
      dimension av(4), aw(4), sog(2,3)
      common /q10/ pc, pm, tc, tm, iim, imesh, juc, v0, w0, rcx, wvolf,
     1 at(5), rmatch
      common /q30/ aw1(4), w1, ethrs, ecut
      common /q55/ aso(4)
      common /q99/ rsot, asot, riv, aiv, awv(4), wvo
c     data sog /60hsaxon derivative            gaussian        volume sa
c    1xon    /
c
      write (4,50)
      write (4,60)
      write (4,70) pc,pm,tc,tm
      if (iim.ne.4)go to 20
   10 write (4,80) av,aw,aw1
      go to 30
   20 continue
c     write (4,70) (sog(i,iim),i=1,2)
      write (4,40)
      write (4,90) av
      write (4,100) aw
      if (ecut.gt.0.) write (4,110) ecut,aw1
   30 continue
      write (4,120) (at(k),k=1,4),aso
      write (4,130) rsot,asot
      write (4,140) (awv(k),k=1,4),riv,aiv
c
      return
c
   40 format (20x,'saxon real well',26x,'imaginary well'/)
   50 format (/)
   60 format (10x,'transmission coefficients calculated from the ','foll
     1owing'/10x,'optical model parameters'/)
   70 format (20x,'projectile     charge',f5.1,10x,'mass',1pe14.6/20x,'t
     1arget       charge',0pf5.1,10x,'mass',1pe14.6/)
   80 format (20x' saxon real well and combination vol and saxon derivat
     1ive imag well'/20x' v = '4f8.2/20x' w(vol) = '4f8.2/20x' w(d) = '
     2 4f8.2)
   90 format (20x' v = 'f8.3' + ('f8.3') x e + ('f8.3') x e2 + ('f8.3')
     1x ln(e)'/)
  100 format (20x' w = 'f8.3' + ('f8.3') x e + ('f8.3') x e2 + ('f8.3')
     1x ln(e)'/)
  110 format (20x'above e =  'f10.3'  mev,the following surf deriv pot i
     1s used'/20x' w = 'f8.3' + ('f8.3' )xe + ('f8.3' )xe2 +( 'f8.3' )xl
     2n(e)'/)
  120 format (20x,'rr = ',1pe14.6/20x,'ar = ',1pe14.6/20x,'ri = ',1pe14.
     1 6/20x,'ai = ',1pe14.6//20x,'vso = ',0p,f8.3' + ('f8.3') x e + ('f
     2 8.3')x e2 + ('f8.3') x ln(e)'/)
  130 format (20x'rso = ',1pe14.6/20x'aso = '1pe14.6//)
  140 format (/20x' the following saxon vol param are used '/20x'wvo= 'f
     1 8.3' + ('f8.3')x e +('f8.3')x e2 + ('f8.3') x ln(e)'/20x'riv = '1
     2 pe14.6/20x'aiv = '1pe14.6//)
      end
      subroutine faze2 (zwr,zwi,zwpr,zwpi,zg,zf,zgp,zfp,cr,ci)
      pr(ar,ai,br,bi)=ar*br-ai*bi
      pi(ar,ai,br,bi)=ar*bi+br*ai
      qr(tr,ti,br,bi)=(tr/br+(ti/br)*(bi/br))/(1.+(bi/br)**2)
      qi(tr,ti,br,bi)=(ti/br-(tr/br)*(bi/br))/(1.+(bi/br)**2)
      amx=amax1(abs(zwr),abs(zwi),abs(zwpr),abs(zwpi),abs(zg)
     1,abs(zf),abs(zgp),abs(zfp))
      amn=amin1(abs(zwr),abs(zwi),abs(zwpr),abs(zwpi),abs(zg)
     1,abs(zf),abs(zgp),abs(zfp))
      fac=sqrt(amn*amx)
      wr=zwr/fac
      wi=zwi/fac
      wpr=zwpr/fac
      wpi=zwpi/fac
      g=zg/fac
      f=zf/fac
      gp=zgp/fac
      fp=zfp/fac
      d=-1.*f
      dp=-1.*fp
      tr=-pr(wr,wi,gp,dp)+pr(wpr,wpi,g,d)
      ti=-pi(wr,wi,gp,dp)+pi(wpr,wpi,g,d)
      br=pr(wpr,wpi,g,f)-pr(wr,wi,gp,fp)
      bi=pi(wpr,wpi,g,f)-pi(wr,wi,gp,fp)
      cr=qr(tr,ti,br,bi)
      ci=qi(tr,ti,br,bi)
      return
      end
      function ranf(i)
      data irn/35791/
      ranf=ran(irn)
      return
      end
