c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c     where M_B is the Born amplitude and 
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
      subroutine setvirtual(p,vflav,res_virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      integer nleg
      parameter (nleg=nlegborn)
      real * 8 p(0:3,nleg)
      integer vflav(nleg)
      real * 8 res_virtual

      logical ini
      data ini/.true./
      save ini

      integer ileg,vflav_loc(nleg)

      integer nf
      real *8 alphas,gs,gw_loc,s,m,m2t,q2,t1,u1,t,u,lns,lnt,lnu,musq,mu2
      real *8 dotp

      double complex qlI1,qlI2,qlI3,qlI4

      real *8 as1f,bs1f,bs2f,bs3f,bs4f,bs5f,
     $cs1f,cs2f,cs3f,cs4f,cs5f,cs6f,cs7f,cs8f,ds1f,ds2f,ds3f

      integer ep

cccccccccccccccccccccccccccccccc    
c     common bl. originally present in lh_readin, needed
c     by my_setpara
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb
ccccccccccccccccccccccccccccccccc

      real *8 ckmtmp
      real *8 vep,vep2,born_pow,born,ewcoupl,virtual_ep0,mt
      real *8 A0sq
      external A0sq

c     check
      if (abs(vflav(4)).ne.6) then
         write(*,*) 'setvirtual: ERROR in flavor assignement'
         call exit(1)
      endif

ccccccccccccccccccccccccccccccccccccccc
c     charge conjugation
c     if ttype=-1, then bflav has been filled with tbar-production flavours.
c     Subroutines here work for t-production flavour assignment.
c     Therefore, invert the sign of local flavours.
      do ileg=1,nleg
         vflav_loc(ileg)= ttype *vflav(ileg)
      enddo
ccccccccccccccccccccccccccccccccccccccc




c     initialise QCDloop libraries
      if(ini) then
         call qlinit
         ini=.false.
      endif

c     define symbols used below
c     ca, cf, nc, pi already defined in pwhg_math.h

      nf=5

      alphas=st_alpha
      gs=sqrt(4d0*pi*alphas)
      gw_loc=sqrt(4d0*pi*alphaem_pow/sthw2_pow)

      
      s=2*dotp(p(0,1),p(0,2))
      m=topmass_pow
      m2t=m**2
      q2=wmass_pow**2
      mu2=st_muren2


      if(vflav_loc(1).eq.0) then
c     the computation assumes g b -> t w and t1 and u1 are
c     defined as
c     t1 = t-m2t
c     u1 = u-m2t
c     Here, in the code, momenta are ordered such that 
c     g b -> w t
         t1=-2*dotp(p(0,1),p(0,4))
         u1=-2*dotp(p(0,1),p(0,3)) + wmass_pow**2 - m2t 
         
         t=t1+m2t
         u=u1+m2t

c     CKM coupling
         if(vflav_loc(2).eq.1) then
            ckmtmp=Vtd**2
         elseif(vflav_loc(2).eq.3) then
            ckmtmp=Vts**2
         elseif(vflav_loc(2).eq.5) then
            ckmtmp=Vtb**2
         else
            write(*,*) 'Error 1 in setvirtual'
            call exit(1)
         endif

      else
         t1=-2*dotp(p(0,1),p(0,3)) + wmass_pow**2 - m2t 
         u1=-2*dotp(p(0,1),p(0,4))

         t=t1+m2t
         u=u1+m2t

c     CKM coupling
         if(vflav_loc(1).eq.1) then
            ckmtmp=Vtd**2
         elseif(vflav_loc(1).eq.3) then
            ckmtmp=Vts**2
         elseif(vflav_loc(1).eq.5) then
            ckmtmp=Vtb**2
         else
            write(*,*) 'Error 2 in setvirtual'
            call exit(1)
         endif

      endif

      lns=log(s/m2t)
      lnt=log(-t1/m2t)
      lnu=log(-u1/m2t)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Since a factor (mu/topmass)^2ep is included explicitely in CW result,
c     I evaluate QCDloop integrals using musq=topmass^2, since
c     (mu)^2ep * Int = (mu/topmass)^2ep * (topmass^2ep * Int)
      musq=m2t

c     Moreover, in doing this, we have the prefactor (mu/topmass)^(2ep) = (mu/Q)^2ep * (Q/topmass)^2ep

c     The POWHEG-BOX expression has to be obtained by setting Q=mu. In CW result, the extra-terms
c     arising from the expansion of (mu/topmass)^2ep have to be included. In the last version, this
c     expansion was not done with mathematica.
c     Added extra pieces by hand in the fortran code, see later.

c     I want finite parts !
      ep=0

c     Evaluation of EZ scalar integrals. These have been checked with CW ones
c     (see old revision of this file).

      as1f=dreal(qlI1(m2t,               musq,ep))

      bs1f=dreal(qlI2(s,0d0,0d0,         musq,ep))
      bs2f=dreal(qlI2(m2t,0d0,m2t,       musq,ep))
      bs3f=dreal(qlI2(u,0d0,m2t,         musq,ep))       
      bs4f=dreal(qlI2(q2,0d0,m2t,        musq,ep))
      bs5f=dreal(qlI2(t,0d0,m2t,         musq,ep))

      cs1f=dreal(qlI3(s  ,0d0,0d0,0d0,0d0,0d0,      musq,ep))
      cs2f=dreal(qlI3(u  ,m2t,0d0,0d0,m2t,0d0,      musq,ep))
      cs3f=dreal(qlI3(t  ,m2t,0d0,0d0,m2t,0d0,      musq,ep))
      cs4f=dreal(qlI3(s  ,q2 ,m2t,0d0,0d0,m2t,      musq,ep))
      cs5f=dreal(qlI3(0d0,q2 ,u  ,0d0,0d0,m2t,      musq,ep))  
      cs6f=dreal(qlI3(0d0,q2 ,t  ,0d0,0d0,m2t,      musq,ep))
      cs7f=dreal(qlI3(m2t,0d0,t  ,0d0,m2t,m2t,      musq,ep))
      cs8f=dreal(qlI3(u  ,0d0,q2 ,0d0,m2t,m2t,      musq,ep))

      ds1f=dreal(qlI4(0d0,0d0,m2t ,q2 ,s   ,u   ,0d0,0d0,0d0,m2t,
     $     musq,ep))
      ds2f=dreal(qlI4(0d0,0d0,m2t ,q2 ,s   ,t   ,0d0,0d0,0d0,m2t,
     $     musq,ep))
      ds3f=dreal(qlI4(0d0,m2t,0d0 ,q2 ,t   ,u   ,0d0,0d0,m2t,m2t,
     $     musq,ep))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     CW result, from virtual_re2.nb
c     (gw->gw_loc, deleted '\' character, m**->mt**)

      mt=topmass_pow

c     double pole
c$$$      virtual_ep2=-gs**4*gw**2*(mu2/mt**2)**ep*
c$$$     -     (Born*CF*(CA + CF)*NC)/(8.*ep**2*Pi**2)/
c$$$     -     (64.*CF*NC**2)


c     single pole
c$$$      virtual_ep1=-gs**4*gw**2*(mu2/mt**2)**ep*
c$$$     -     ( - 
c$$$     -     (Born*CF*NC*(CA*(-11 + 6*lns + 6*lnt - 6*lnu) + 
c$$$     -     3*CF*(-5 + 4*lnu) + 2*nf))/(48.*ep*Pi**2)/
c$$$     -     (64.*CF*NC**2))


c     From pole structure, it follows that (in this version) of virtual
c     corrections, we have
c     born_pow=1/64./nc *gs**2 *gw**2 *Born.

c     Crossing if gluon is leg 2 has been already performed by exchanging
c     t1 and u1. The following lines have been taken from Born.f
c     ew coupling
      ewcoupl=4d0*pi*alphaem_pow/sthw2_pow
      born_pow=(4*pi*st_alpha)*(ewcoupl/8.)/4./3./8.
      born_pow=born_pow *A0sq(s,u1,t1,m2t,wmass_pow**2)

c     Therefore:
      Born=born_pow*(1/64./nc *gs**2 *gw_loc**2)**(-1)
      
c     In the following line I delete from the fortran code
c     the factor (mu2/mt**2)**ep. Effects of its expansion on
c     the double and single pole is taken into account later, at the
c     end.
c$$$      virtual_ep0=-gs**4*gw_loc**2*(mu2/mt**2)**ep*
      virtual_ep0=-gs**4*gw_loc**2 *
     -     ( - 
     -       (CA*CF*CS3f*NC*(2*mt**6 + 2*mt**4*(s + t1) + 
     -            2*q2*(2*q2**2 + s**2 + t1**2 - 2*q2*(s + t1)) + 
     -            mt**2*(-6*q2**2 + 2*q2*(s + t1) + (s + t1)**2)))/
     -        (2.*Pi**2*q2*s) + 
     -       ((CA - 2*CF)*CF*CS5f*NC*
     -          (2*mt**6*(s + t1) + mt**4*(s**2 + 3*s*t1 + 2*t1**2) + 
     -            2*q2*(s + t1)*
     -             (2*q2**2 + 2*s**2 + 2*s*t1 + t1**2 - 2*q2*(2*s + t1))
     -              + mt**2*(s + t1)*
     -             (-6*q2**2 + (s + t1)**2 + 2*q2*(2*s + t1))))/
     -        (2.*Pi**2*q2*s*t1) - 
     -       (2*BS3f*(CA - 2*CF)*CF*NC*
     -          (2*mt**4*q2*(s + t1) + 
     -            2*t1*(-q2 + s + t1)**2*(q2 + s + t1) + 
     -            mt**2*(2*q2*(s + t1)**2 + (s + t1)**3 - 
     -               2*q2**2*(s + 2*t1))))/
     -        (Pi**2*t1*(s + t1)**2*(mt**2 - q2 + s + t1)) - 
     -       ((CA - 2*CF)*CF*CS7f*NC*
     -          (2*mt**6*(s**2 + s*t1 + t1**2) + 
     -            mt**4*(t1**3 + 2*q2*s*(2*s + t1)) + 
     -            2*q2*t1**2*
     -             (2*q2**2 + s**2 + 2*s*t1 + 2*t1**2 - 2*q2*(s + 2*t1))
     -              + mt**2*t1*
     -             (t1*(s + t1)**2 + 2*q2*t1*(3*s + 2*t1) - 
     -               2*q2**2*(2*s + 3*t1))))/(2.*Pi**2*q2*s*t1**2) - 
     -       (CA*CF*DS2f*NC*(2*mt**6*(2*s + t1) + 
     -            2*q2*t1*(2*q2**2 + s**2 + t1**2 - 2*q2*(s + t1)) + 
     -            2*mt**4*(2*q2*s + t1*(s + t1)) + 
     -            mt**2*(2*q2*t1*(s + t1) + t1*(s + t1)**2 - 
     -               2*q2**2*(4*s + 3*t1))))/(2.*Pi**2*q2*t1) - 
     -       ((CA - 2*CF)*CF*CS2f*NC*
     -          (mt**6*(s + t1) - mt**4*(3*q2 - s - t1)*(s + t1) + 
     -            2*mt**2*q2*(3*q2*(s + t1) - 
     -               2*(s**2 + 3*s*t1 + t1**2)) - 
     -            2*q2*(s**3 + 5*s**2*t1 + 5*s*t1**2 + t1**3 + 
     -               2*q2**2*(s + t1) - q2*(3*s**2 + 8*s*t1 + 3*t1**2)))
     -          )/(2.*Pi**2*q2*s*t1) - 
     -       ((CA - 2*CF)*CF*DS1f*NC*
     -          (2*mt**8 + mt**6*(-2*q2 + 3*s + 4*t1) + 
     -            mt**4*(-6*q2**2 + 3*q2*s + 2*s**2 + 5*s*t1 + 
     -               3*t1**2) + 
     -            2*q2*(-2*q2**3 + 2*s**3 + 4*s**2*t1 + 3*s*t1**2 + 
     -               t1**3 + q2**2*(6*s + 4*t1) - 
     -               q2*(6*s**2 + 8*s*t1 + 3*t1**2)) + 
     -            mt**2*(10*q2**3 + (s + t1)**3 - 6*q2**2*(3*s + 2*t1) + 
     -               q2*(7*s**2 + 8*s*t1 + 3*t1**2))))/(2.*Pi**2*q2*t1)!!!!!!!!:\\
     -        + ((CA - 2*CF)*CF*CS8f*NC*
     -          (2*mt**6*(s + t1)**2 + 
     -            mt**4*(s + t1)*(8*q2*s + t1*(s + t1)) + 
     -            2*q2*(s + t1)**2*
     -             (2*q2**2 + s**2 + 2*s*t1 + 2*t1**2 - 2*q2*(s + 2*t1))
     -              + mt**2*((s + t1)**4 + 
     -               2*q2*(s + t1)**2*(3*s + 2*t1) - 
     -               2*q2**2*(3*s**2 + 10*s*t1 + 3*t1**2))))/
     -        (2.*Pi**2*q2*s*t1*(s + t1)) - 
     -       ((CA - 2*CF)*CF*DS3f*NC*
     -          (2*mt**8*(2*s + t1) + 
     -            mt**6*(4*s**2 + 8*s*t1 + t1*(-2*q2 + 3*t1)) + 
     -            mt**4*(-6*q2**2*(2*s + t1) + 
     -               t1*(3*s**2 + 5*s*t1 + 2*t1**2) + 
     -               q2*(4*s**2 + 4*s*t1 + 3*t1**2)) + 
     -            2*q2*t1*(-2*q2**3 + s**3 + 3*s**2*t1 + 4*s*t1**2 + 
     -               2*t1**3 + q2**2*(4*s + 6*t1) - 
     -               q2*(3*s**2 + 8*s*t1 + 6*t1**2)) + 
     -            mt**2*(t1*(s + t1)**3 + 2*q2**3*(4*s + 5*t1) + 
     -               q2*t1*(3*s**2 + 8*s*t1 + 7*t1**2) - 
     -               2*q2**2*(4*s**2 + 10*s*t1 + 9*t1**2))))/
     -        (2.*Pi**2*q2*s*t1) - 
     -       (CF*CS1f*NC*(-2*CF*
     -             (2*mt**6 + mt**4*(s + 2*t1) + 
     -               2*q2*(2*q2**2 + 2*s**2 + 2*s*t1 + t1**2 - 
     -                  2*q2*(2*s + t1)) + 
     -               mt**2*(-6*q2**2 + (s + t1)**2 + 2*q2*(2*s + t1))) + 
     -            CA*(4*mt**6 + mt**4*(3*s + 4*t1) + 
     -               2*q2*(4*q2**2 - 6*q2*s + 3*s**2 - 4*q2*t1 + 
     -                  2*s*t1 + 2*t1**2) + 
     -              2*mt**2*(-6*q2**2 + (s + t1)**2 + q2*(3*s + 2*t1))))
     -          )/(2.*Pi**2*q2*t1) + 
     -       (CF*CS6f*NC*(mt**2 - q2 + t1)*
     -         (CA*(mt**4 - 2*mt**2*q2 + 2*q2*(2*q2 - 2*s - t1))*t1**2 + 
     -            2*CF*(2*mt**6*(2*s + t1) + 
     -               mt**4*(4*q2*s + t1*(2*s + t1)) + 
     -               2*q2*t1*
     -                (2*q2**2 + s**2 + 2*s*t1 + 2*t1**2 - 
     -                  2*q2*(s + 2*t1)) + 
     -               mt**2*(t1*(s + t1)**2 + 2*q2*t1*(s + 2*t1) - 
     -                  2*q2**2*(4*s + 3*t1)))))/(2.*Pi**2*q2*s*t1**2)!!!!!!!:\\
     -        - (BS5f*CF*NC*(CA*(mt**2 + t1)*(mt**2 - q2 + t1)*
     -             (-2*q2*t1**4 + mt**8*(s + t1) + 
     -               2*mt**2*q2*t1*(q2**2 - q2*s + (s - 2*t1)*t1) + 
     -               mt**6*(q2*s + 2*t1*(s + t1)) + 
     -               mt**4*(q2**2*(-2*s - 3*t1) + q2*(3*s - 2*t1)*t1 + 
     -                  t1**2*(s + t1))) - 
     -            CF*(2*mt**12*(s + t1) - 
     -               2*q2*t1**4*
     -                (4*q2**2 - 4*q2*s - 6*q2*t1 + 3*s*t1 + 2*t1**2) + 
     -               mt**10*t1*(-2*q2 + 7*(s + t1)) + 
     -               mt**8*(q2*(4*s - t1)*t1 - 6*q2**2*(s + t1) + 
     -                  9*t1**2*(s + t1)) - 
     -               2*mt**2*q2*t1**2*
     -                (3*q2**3 + q2**2*(-s + t1) + t1**2*(9*s + 7*t1) - 
     -                  q2*t1*(6*s + 13*t1)) + 
     -               mt**4*t1*
     -                (-4*q2**4 + q2**2*t1*(-7*s + t1) + 
     -                  t1**3*(s + t1) - q2*t1**2*(14*s + 13*t1) + 
     -                  q2**3*(6*s + 19*t1)) + 
     -               mt**6*(2*q2*(s - t1)*t1**2 + 5*t1**3*(s + t1) + 
     -                  2*q2**3*(2*s + 5*t1) - q2**2*t1*(17*s + 19*t1)))
     -           ))/(Pi**2*q2*t1**3*(mt**2 + t1)*(mt**2 - q2 + t1)**2) - 
     -       (2*BS1f*CF*NC*(-(CA*q2*s*
     -               (mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))*
     -               ((-q2 + s)*t1 + mt**2*(2*s + t1))) + 
     -            CF*(-2*mt**10*(s + t1) + 
     -               q2*(q2 - s)**2*t1*
     -                (4*q2**2 - 6*q2*s + 2*s**2 - 4*q2*t1 + 3*s*t1) + 
     -               2*mt**8*
     -                (2*q2*s + 2*s**2 + 3*q2*t1 + s*t1 - t1**2) + 
     -               mt**6*(-2*q2**2*(3*s + 5*t1) + 
     -                  2*s*(-s**2 + s*t1 + 2*t1**2) + 
     -                  q2*(3*s**2 + 6*s*t1 + 4*t1**2)) + 
     -               mt**4*(-2*s**2*t1*(s + t1) + 2*q2**3*(4*s + 7*t1) + 
     -                  2*q2**2*(5*s**2 - 4*s*t1 - 3*t1**2) + 
     -                  q2*s*(-10*s**2 - 8*s*t1 + 7*t1**2)) - 
     -               mt**2*q2*
     -                (4*q2**3*(s + 3*t1) + 
     -                  2*q2*s*(s**2 - 6*s*t1 - t1**2) - 
     -                  q2**2*(7*s**2 + 14*s*t1 + 8*t1**2) + 
     -                  s**2*(s**2 + 14*s*t1 + 18*t1**2)))))/
     -       (Pi**2*q2*(mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))**2*t1) + 
     -       (CF*CS4f*NC*(CA*s*(mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))*
     -             (mt**10 - 2*q2*(q2 - s)**3*(2*q2 - s - 2*t1) - 
     -               mt**8*(5*q2 + s - 2*t1) + 
     -               mt**6*(13*q2**2 + q2*(6*s - 4*t1) - s*(s + 4*t1)) + 
     -               2*mt**2*q2*
     -                (7*q2**3 - 7*q2*s**2 - 2*q2**2*(3*s + 2*t1) + 
     -                  2*s**2*(3*s + 4*t1)) + 
     -               mt**4*(-19*q2**3 + s**2*(s + 2*t1) + 
     -                  q2**2*(-7*s + 6*t1) - q2*s*(11*s + 8*t1))) + 
     -            2*CF*(2*mt**16 + mt**14*(-10*q2 - 9*s + 2*t1) + 
     -               mt**12*(14*q2**2 + 19*q2*s + 20*s**2 - 8*q2*t1 - 
     -                  4*s*t1 + t1**2) - 
     -               2*q2*(q2 - s)**5*
     -                (2*q2**2 + 2*s**2 + 2*s*t1 + t1**2 - 
     -                  2*q2*(2*s + t1)) - 
     -               2*mt**8*(35*q2**4 - 6*s**4 + q2**2*s*(s - 16*t1) + 
     -                  5*q2**3*(3*s - 2*t1) + 4*s**3*t1 + 
     -                  q2*s*(9*s**2 + 6*s*t1 + 2*t1**2)) + 
     -               mt**10*(14*q2**3 + 6*q2**2*(2*s + t1) - 
     -                  s*(23*s**2 - 6*s*t1 + t1**2) - 
     -                  q2*(15*s**2 + 4*s*t1 + 3*t1**2)) + 
     -               mt**6*(98*q2**5 - s**4*(s - 6*t1) - 
     -                  5*q2**4*(13*s + 10*t1) + 
     -                  2*q2*s**2*(13*s**2 + 4*s*t1 - 9*t1**2) - 
     -                  2*q2**3*(10*s**2 + 4*s*t1 - 5*t1**2) + 
     -                  2*q2**2*s*(11*s**2 + 30*s*t1 + 4*t1**2)) + 
     -               mt**2*(q2 - s)*
     -                (26*q2**6 - s**4*(s + t1)**2 - 
     -                  22*q2**5*(4*s + t1) - 
     -                  4*q2**2*s**2*(9*s**2 + 9*s*t1 + 4*t1**2) + 
     -                  2*q2*s**3*(8*s**2 + 9*s*t1 + 5*t1**2) - 
     -                  2*q2**3*s*(4*s**2 + 6*s*t1 + 7*t1**2) + 
     -                  q2**4*(91*s**2 + 54*s*t1 + 9*t1**2)) - 
     -               mt**4*(70*q2**6 - 3*q2**5*(53*s + 16*t1) + 
     -                  2*q2**2*s**2*(3*s**2 - 4*s*t1 - 20*t1**2) + 
     -                  2*q2**3*s*(11*s**2 + 24*s*t1 - 5*t1**2) + 
     -                  s**4*(2*s**2 + 4*s*t1 + t1**2) + 
     -                  q2**4*(78*s**2 + 68*s*t1 + 15*t1**2) - 
     -                  q2*s**3*(19*s**2 + 48*s*t1 + 42*t1**2)))))/
     -       (2.*Pi**2*q2*s*(mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))**2*
     -          t1) + (BS2f*CF*NC*
     -          (CA*mt**2*s*(mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))*
     -             (mt**2 + t1)*(mt**2 - q2 + t1)*
     -             (mt**10*(s + t1) - 
     -               mt**8*(s**2 - s*t1 - 2*t1**2 + q2*(s + 2*t1)) - 
     -               mt**6*(s**3 + 4*s**2*t1 + 2*s*t1**2 - t1**3 + 
     -                  q2**2*(3*s + 2*t1) + 
     -                  q2*(4*s**2 + 4*s*t1 + 6*t1**2)) + 
     -               2*q2*t1*
     -                (q2**4 - 4*q2**3*s + q2**2*(6*s**2 - 2*t1**2) + 
     -                  q2*(-4*s**3 + 7*s*t1**2) + 
     -                  s*(s**3 - 5*s*t1**2 - 2*t1**3)) + 
     -               mt**4*(q2**3*(5*s + 8*t1) - 
     -                  q2**2*(s**2 + 3*s*t1 - 6*t1**2) + 
     -                  s*(s**3 + s**2*t1 - 2*s*t1**2 - 2*t1**3) - 
     -                  q2*(5*s**3 + 14*s**2*t1 + 4*s*t1**2 + 6*t1**3))!!!!!!!!!!!!!:\\
     -                + mt**2*
     -                (q2**4*(-2*s - 7*t1) + s**2*t1*(s + t1)**2 + 
     -                  2*q2*s**2*(s**2 - s*t1 - 7*t1**2) + 
     -                  2*q2**3*(3*s**2 + 7*s*t1 - t1**2) + 
     -                  q2**2*
     -                   (-6*s**3 - 6*s**2*t1 + 6*s*t1**2 + 9*t1**3)))!!!!!!!!!!!!!:\\
     -             + CF*(6*mt**22*(2*s + t1) + 
     -               mt**20*(-38*s**2 + 35*s*t1 + 25*t1**2 - 
     -                  12*q2*(5*s + 3*t1)) + 
     -               mt**18*(30*s**3 - 163*s**2*t1 + 15*s*t1**2 + 
     -                  40*t1**3 + 12*q2**2*(7*s + 6*t1) + 
     -                  q2*(56*s**2 - 205*s*t1 - 138*t1**2)) + 
     -               4*q2*(q2 - s)**4*(q2 - t1)*t1**3*
     -                (2*q2**3 - s**3 - s**2*t1 - s*t1**2 - t1**3 - 
     -                  4*q2**2*(s + t1) + 
     -                  q2*(3*s**2 + 4*s*t1 + 3*t1**2)) + 
     -               mt**16*(84*q2**3*s + 20*s**4 + 180*s**3*t1 - 
     -                  257*s**2*t1**2 - 52*s*t1**3 + 29*t1**4 + 
     -                  q2**2*(174*s**2 + 389*s*t1 + 246*t1**2) + 
     -                  q2*(22*s**3 + 206*s**2*t1 - 270*s*t1**2 - 
     -                     207*t1**3)) + 
     -               mt**14*(-40*s**5 + 370*s**3*t1**2 - 
     -                  174*s**2*t1**3 - 77*s*t1**4 + 7*t1**5 - 
     -                  84*q2**4*(5*s + 3*t1) + 
     -                  q2**3*(-380*s**2 + 11*s*t1 + 42*t1**2) + 
     -                  q2**2*
     -                   (-32*s**3 + 729*s**2*t1 + 883*s*t1**2 + 
     -                     339*t1**3) + 
     -                  q2*(40*s**4 + 158*s**3*t1 + 352*s**2*t1**2 - 
     -                     149*s*t1**3 - 149*t1**4)) + 
     -               mt**12*(18*s**6 - 101*s**5*t1 - 105*s**4*t1**2 + 
     -                  358*s**3*t1**3 - 31*s**2*t1**4 - 41*s*t1**5 - 
     -                  2*t1**6 + 84*q2**5*(7*s + 6*t1) - 
     -                  5*q2**4*(2*s**2 + 233*s*t1 + 168*t1**2) - 
     -                  q2**3*
     -                   (116*s**3 + 1366*s**2*t1 + 762*s*t1**2 - 
     -                     81*t1**3) + 
     -                  q2**2*
     -                   (40*s**4 - 132*s**3*t1 + 1051*s**2*t1**2 + 
     -                     1055*s*t1**3 + 237*t1**4) + 
     -                  q2*(-72*s**5 + 84*s**4*t1 + 386*s**3*t1**2 + 
     -                     396*s**2*t1**3 + 21*s*t1**4 - 47*t1**5)) + 
     -               mt**2*(q2 - s)*t1**2*
     -                (20*q2**8 - s**3*(2*s - t1)*t1**2*(s + t1)**2 - 
     -                  3*q2**7*(36*s + 37*t1) + 
     -                  q2**6*(240*s**2 + 425*s*t1 + 193*t1**2) - 
     -                  2*q2**5*
     -                   (140*s**3 + 339*s**2*t1 + 255*s*t1**2 + 
     -                     76*t1**3) + 
     -                  q2*s**2*t1*
     -                   (-12*s**4 + 16*s**3*t1 + 26*s**2*t1**2 + 
     -                     17*s*t1**3 + 19*t1**4) + 
     -                  q2**4*
     -                   (180*s**4 + 574*s**3*t1 + 518*s**2*t1**2 + 
     -                     243*s*t1**3 + 61*t1**4) + 
     -                  q2**2*s*
     -                   (8*s**5 + 73*s**4*t1 + 19*s**3*t1**2 - 
     -                     36*s**2*t1**3 - 41*s*t1**4 + 3*t1**5) - 
     -                  q2**3*
     -                   (60*s**5 + 271*s**4*t1 + 234*s**3*t1**2 + 
     -                     78*s**2*t1**3 + 49*s*t1**4 + 11*t1**5)) - 
     -               mt**10*(2*s**7 - 45*s**6*t1 + 93*s**5*t1**2 + 
     -                  152*s**4*t1**3 - 184*s**3*t1**4 - 
     -                  21*s**2*t1**5 + 8*s*t1**6 + t1**7 + 
     -                  84*q2**6*(5*s + 6*t1) - 
     -                  q2**5*(528*s**2 + 1985*s*t1 + 1470*t1**2) + 
     -                  q2**4*
     -                   (130*s**3 + 241*s**2*t1 + 1035*s*t1**2 + 
     -                     1053*t1**3) + 
     -                  q2**3*
     -                   (80*s**4 + 524*s**3*t1 + 1924*s**2*t1**2 + 
     -                     1346*s*t1**3 - 59*t1**4) - 
     -                  q2**2*
     -                   (88*s**5 + 460*s**4*t1 + 68*s**3*t1**2 + 
     -                     643*s**2*t1**3 + 609*s*t1**4 + 78*t1**5) + 
     -                  q2*(-16*s**6 + 229*s**5*t1 + 74*s**4*t1**2 - 
     -                     336*s**3*t1**3 - 276*s**2*t1**4 - 
     -                     79*s*t1**5 + t1**6)) + 
     -               mt**4*t1*
     -                (12*q2**9 - q2**8*(106*s + 177*t1) + 
     -                  q2**7*(350*s**2 + 882*s*t1 + 531*t1**2) + 
     -                  s**4*t1**2*
     -                   (6*s**3 + 11*s**2*t1 - 5*s*t1**2 - 10*t1**3) - 
     -                  q2**6*
     -                   (580*s**3 + 1815*s**2*t1 + 1725*s*t1**2 + 
     -                     601*t1**3) + 
     -                  q2**5*
     -                   (520*s**4 + 1986*s**3*t1 + 2162*s**2*t1**2 + 
     -                     1131*s*t1**3 + 309*t1**4) + 
     -                  q2**2*s*t1*
     -                   (-81*s**5 + 143*s**4*t1 + 463*s**3*t1**2 + 
     -                     370*s**2*t1**3 + 34*s*t1**4 - 4*t1**5) + 
     -                  q2*s**2*t1*
     -                   (6*s**5 - 76*s**4*t1 - 111*s**3*t1**2 - 
     -                     91*s**2*t1**3 - 76*s*t1**4 - 2*t1**5) - 
     -                  q2**4*
     -                   (242*s**5 + 1239*s**4*t1 + 1200*s**3*t1**2 + 
     -                     485*s**2*t1**3 + 229*s*t1**4 + 80*t1**5) + 
     -                  q2**3*
     -                   (46*s**6 + 438*s**5*t1 + 159*s**4*t1**2 - 
     -                     472*s**3*t1**3 - 178*s**2*t1**4 + 
     -                     16*s*t1**5 + 10*t1**6)) - 
     -               mt**6*(6*q2**8*(4*s + 15*t1) - 
     -                  q2**7*(116*s**2 + 641*s*t1 + 654*t1**2) + 
     -                  q2**6*
     -                   (220*s**3 + 1477*s**2*t1 + 2295*s*t1**2 + 
     -                     1263*t1**3) - 
     -                  q2**5*
     -                   (200*s**4 + 1566*s**3*t1 + 2792*s**2*t1**2 + 
     -                     2143*s*t1**3 + 921*t1**4) - 
     -                  s**2*t1**2*
     -                   (4*s**5 + 22*s**4*t1 - 21*s**3*t1**2 - 
     -                     39*s**2*t1**3 + 8*s*t1**4 + 2*t1**5) + 
     -                  q2**4*
     -                   (80*s**5 + 892*s**4*t1 + 1230*s**3*t1**2 + 
     -                     598*s**2*t1**3 + 245*s*t1**4 + 269*t1**5) + 
     -                  q2**3*
     -                   (-4*s**6 - 325*s**5*t1 + 42*s**4*t1**2 + 
     -                     1080*s**3*t1**3 + 676*s**2*t1**4 + 
     -                     202*s*t1**5 - 36*t1**6) + 
     -                  q2*s*t1*
     -                   (4*s**6 + 44*s**5*t1 + 155*s**4*t1**2 + 
     -                     181*s**3*t1**3 + 114*s**2*t1**4 - 
     -                     2*s*t1**5 - 6*t1**6) - 
     -                  q2**2*
     -                   (4*s**7 - 69*s**6*t1 + 161*s**5*t1**2 + 
     -                     931*s**4*t1**3 + 726*s**3*t1**4 + 
     -                     44*s**2*t1**5 - 6*s*t1**6 - 4*t1**7)) + 
     -               mt**8*(12*q2**7*(13*s + 24*t1) - 
     -                  q2**6*(446*s**2 + 1585*s*t1 + 1302*t1**2) + 
     -                  q2**5*
     -                   (446*s**3 + 1962*s**2*t1 + 2710*s*t1**2 + 
     -                     1635*t1**3) - 
     -                  q2**4*
     -                   (220*s**4 + 668*s**3*t1 + 547*s**2*t1**2 + 
     -                     414*s*t1**3 + 665*t1**4) + 
     -                  s**2*t1*
     -                   (-2*s**5 + 41*s**4*t1 - 48*s**3*t1**2 - 
     -                     97*s**2*t1**3 + 54*s*t1**4 + 12*t1**5) + 
     -                  2*q2**3*
     -                   (52*s**5 - 86*s**4*t1 - 530*s**3*t1**2 - 
     -                     725*s**2*t1**3 - 432*s*t1**4 + 23*t1**5) + 
     -                  q2**2*
     -                   (-38*s**6 + 151*s**5*t1 + 1000*s**4*t1**2 + 
     -                     610*s**3*t1**3 + 176*s**2*t1**4 + 
     -                     130*s*t1**5 + 2*t1**6) + 
     -                  q2*(-2*s**7 + 26*s**6*t1 - 234*s**5*t1**2 - 
     -                     229*s**4*t1**3 + 30*s**3*t1**4 + 
     -                     86*s**2*t1**5 + 40*s*t1**6 + 2*t1**7)))))/
     -        (Pi**2*q2*s*(mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))**2*
     -       t1**3*(mt**2 + t1)*(mt**2 - q2 + t1)*(mt**2 - q2 + s + t1))
     -         + (2*BS4f*CF*NC*
     -          (CA*s*(mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))*t1*
     -             (mt**2 - q2 + t1)*
     -             (2*mt**8*q2*(s + t1) + 
     -               mt**6*((s + t1)**3 - 2*q2**2*(3*s + 4*t1) - 
     -                  2*q2*(2*s**2 + s*t1 - t1**2)) + 
     -               q2*(q2 - s)*t1*
     -                (2*q2**3 - 2*q2**2*(s + t1) + (s + t1)**3 - 
     -                  2*q2*(s**2 + s*t1 + t1**2)) + 
     -               mt**4*(-6*q2**2*t1**2 - (s - t1)*(s + t1)**3 + 
     -                  6*q2**3*(s + 2*t1) - 
     -                  2*q2*t1*(5*s**2 + 7*s*t1 + 2*t1**2)) - 
     -               mt**2*(s*t1*(s + t1)**3 + 2*q2**4*(s + 4*t1) - 
     -                  2*q2**3*(2*s**2 + 3*s*t1 + 3*t1**2) + 
     -                  q2**2*(s**3 + s**2*t1 - 11*s*t1**2 - 5*t1**3) - 
     -                  q2*(s**4 + 5*s**3*t1 + 3*s**2*t1**2 - 
     -                     3*s*t1**3 - 2*t1**4))) + 
     -            CF*(2*mt**14*q2*
     -                (3*s**3 + 7*s**2*t1 + 7*s*t1**2 + 3*t1**3) + 
     -               mt**12*(5*s*t1*(s + t1)**3 - 
     -                  2*q2**2*
     -                   (18*s**3 + 45*s**2*t1 + 46*s*t1**2 + 21*t1**3)!!!!!!!!!:\\
     -                   - 2*q2*
     -                   (12*s**4 + 19*s**3*t1 + 5*s**2*t1**2 - 
     -                     11*s*t1**3 - 9*t1**4)) + 
     -               mt**10*(-6*s*(s - 2*t1)*t1*(s + t1)**3 + 
     -                  6*q2**3*
     -                   (15*s**3 + 41*s**2*t1 + 43*s*t1**2 + 21*t1**3)!!!!!!!!!!!!:\\
     -                   + 12*q2**2*
     -                   (6*s**4 + 8*s**3*t1 - 5*s**2*t1**2 - 
     -                     14*s*t1**3 - 9*t1**4) + 
     -                  q2*(36*s**5 + s**4*t1 - 118*s**3*t1**2 - 
     -                     134*s**2*t1**3 - 30*s*t1**4 + 21*t1**5)) - 
     -               q2*(q2 - s)**2*t1*
     -                (-2*s*t1*(s + t1)**3*(s**2 + s*t1 + t1**2) + 
     -                  q2**5*(6*s**2 + 8*s*t1 + 6*t1**2) - 
     -                  2*q2**4*
     -                   (9*s**3 + 23*s**2*t1 + 23*s*t1**2 + 9*t1**3) + 
     -                  q2**3*
     -                   (21*s**4 + 82*s**3*t1 + 114*s**2*t1**2 + 
     -                     82*s*t1**3 + 21*t1**4) - 
     -                  2*q2**2*
     -                   (6*s**5 + 31*s**4*t1 + 59*s**3*t1**2 + 
     -                     59*s**2*t1**3 + 31*s*t1**4 + 6*t1**5) + 
     -                  q2*(3*s**6 + 20*s**5*t1 + 53*s**4*t1**2 + 
     -                     68*s**3*t1**3 + 53*s**2*t1**4 + 20*s*t1**5 + 
     -                     3*t1**6)) + 
     -               mt**8*(s*t1*(s + t1)**3*
     -                   (3*s**2 - 9*s*t1 + 10*t1**2) - 
     -                  10*q2**4*
     -                   (12*s**3 + 37*s**2*t1 + 40*s*t1**2 + 21*t1**3)!!!!!!!!!!!!!!!:\\
     -                   + q2**3*
     -                   (-48*s**4 + 26*s**3*t1 + 402*s**2*t1**2 + 
     -                     510*s*t1**3 + 270*t1**4) + 
     -                  q2**2*
     -                   (-48*s**5 + 23*s**4*t1 + 260*s**3*t1**2 + 
     -                     186*s**2*t1**3 - 48*s*t1**4 - 105*t1**5) + 
     -                  q2*(-24*s**6 + 25*s**5*t1 + 44*s**4*t1**2 - 
     -                     100*s**3*t1**3 - 180*s**2*t1**4 - 
     -                     73*s*t1**5 + 12*t1**6)) + 
     -               mt**6*(-(s*t1*(s + t1)**3*
     -                     (2*s**3 + 3*s*t1**2 - 3*t1**3)) + 
     -                  10*q2**5*
     -                   (9*s**3 + 33*s**2*t1 + 37*s*t1**2 + 21*t1**3)!!!!!!!!!!!:\\
     -                   - 8*q2**4*
     -                   (6*s**4 + 38*s**3*t1 + 101*s**2*t1**2 + 
     -                     100*s*t1**3 + 45*t1**4) + 
     -                  q2**3*
     -                   (24*s**5 + 66*s**4*t1 + 56*s**3*t1**2 + 
     -                     308*s**2*t1**3 + 392*s*t1**4 + 210*t1**5) + 
     -                  q2**2*
     -                   (24*s**6 + 45*s**5*t1 + 88*s**4*t1**2 + 
     -                     308*s**3*t1**3 + 252*s**2*t1**4 + 
     -                     43*s*t1**5 - 48*t1**6) + 
     -                  q2*(6*s**7 - 30*s**6*t1 - 46*s**5*t1**2 - 
     -                     47*s**4*t1**3 - 98*s**3*t1**4 - 
     -                     110*s**2*t1**5 - 46*s*t1**6 + 3*t1**7)) - 
     -               mt**4*(3*s**3*t1**2*(s + t1)**3*(s + 2*t1) + 
     -                  6*q2**6*
     -                   (6*s**3 + 29*s**2*t1 + 34*s*t1**2 + 21*t1**3)!!!!!!!!!!!!:\\
     -                   - 6*q2**5*
     -                   (12*s**4 + 61*s**3*t1 + 127*s**2*t1**2 + 
     -                     115*s*t1**3 + 45*t1**4) + 
     -                  q2**4*
     -                   (48*s**5 + 293*s**4*t1 + 677*s**3*t1**2 + 
     -                     925*s**2*t1**3 + 663*s*t1**4 + 210*t1**5) + 
     -                  q2*s*t1*
     -                   (-7*s**6 + 8*s**5*t1 + 46*s**4*t1**2 + 
     -                     92*s**3*t1**3 + 99*s**2*t1**4 + 48*s*t1**5 + 
     -                     10*t1**6) - 
     -                  q2**3*
     -                   (24*s**6 + 75*s**5*t1 + 148*s**4*t1**2 + 
     -                     148*s**3*t1**3 + 216*s**2*t1**4 + 
     -                     189*s*t1**5 + 72*t1**6) + 
     -                  q2**2*
     -                   (12*s**7 + 21*s**6*t1 - 66*s**5*t1**2 - 
     -                     175*s**4*t1**3 - 272*s**3*t1**4 - 
     -                     139*s**2*t1**5 - 18*s*t1**6 + 9*t1**7)) + 
     -               mt**2*(-(s**3*t1**3*(s + t1)**3*(s + 3*t1)) + 
     -                  q2**7*
     -                   (6*s**3 + 50*s**2*t1 + 62*s*t1**2 + 42*t1**3)!!!!!!!!!!!!!:\\
     -                   - 4*q2**6*
     -                   (6*s**4 + 44*s**3*t1 + 87*s**2*t1**2 + 
     -                     78*s*t1**3 + 27*t1**4) + 
     -                  q2**5*
     -                   (36*s**5 + 261*s**4*t1 + 646*s**3*t1**2 + 
     -                     762*s**2*t1**3 + 462*s*t1**4 + 105*t1**5) + 
     -                  q2*s**2*t1*
     -                   (3*s**6 + 14*s**5*t1 + 21*s**4*t1**2 + 
     -                     8*s**3*t1**3 - 23*s**2*t1**4 - 34*s*t1**5 - 
     -                     13*t1**6) - 
     -                  q2**2*s*t1*
     -                   (23*s**6 + 82*s**5*t1 + 69*s**4*t1**2 - 
     -                     8*s**3*t1**3 - 73*s**2*t1**4 - 30*s*t1**5 + 
     -                     t1**6) - 
     -                  q2**4*
     -                   (24*s**6 + 211*s**5*t1 + 546*s**4*t1**2 + 
     -                     766*s**3*t1**3 + 618*s**2*t1**4 + 
     -                     275*s*t1**5 + 48*t1**6) + 
     -                  q2**3*
     -                   (6*s**7 + 96*s**6*t1 + 254*s**5*t1**2 + 
     -                     315*s**4*t1**3 + 254*s**3*t1**4 + 
     -                     148*s**2*t1**5 + 62*s*t1**6 + 9*t1**7)))))/
     -       (Pi**2*s*(mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))**2*t1**2*
     -          (mt**2 - q2 + t1)**2*(s + t1)**2) - 
     -       (CF*NC*(-3*(-4*CA*(mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))*
     -                t1*(mt**2 + t1)*(mt**2 - q2 + t1)*
     -                (mt**4*(s + t1)**3 - 2*q2**2*t1*(s**2 + t1**2) + 
     -                  mt**2*q2*(2*s**3 + s**2*t1 + t1**3)) - 
     -               CF*(48*mt**16*(2*s**2 + 3*s*t1 + t1**2) - 
     -                  8*mt**14*(s + t1)*
     -                   (24*s**2 - 35*s*t1 - 23*t1**2 + 
     -                     6*q2*(4*s + 3*t1)) - 
     -                  4*mt**12*(s + t1)*
     -                   (48*q2**2*s - 24*s**3 + 172*s**2*t1 - 
     -                     49*s*t1**2 - 65*t1**3 + 
     -                     2*q2*(24*s**2 + 88*s*t1 + 63*t1**2)) + 
     -                  4*mt**10*(s + t1)*
     -                   (24*q2**3*(8*s + 5*t1) + 
     -                     4*q2**2*(36*s**2 - 8*s*t1 + 3*t1**2) - 
     -                     q2*t1*(228*s**2 + 283*s*t1 + 164*t1**2) + 
     -                     2*t1*
     -                      (43*s**3 - 110*s**2*t1 - 12*s*t1**2 + 
     -                        21*t1**3)) - 
     -                  q2*(q2 - s)*t1**3*
     -                   (32*q2**4*(s + t1) - 
     -                     8*q2**3*(7*s**2 + 18*s*t1 + 7*t1**2) + 
     -                     s*t1*(s + t1)*
     -                      (8*s**2 - Born*s*t1 + 8*t1**2) + 
     -                     q2**2*(s + t1)*
     -                      (40*s**2 - (-112 + Born)*s*t1 + 40*t1**2) + 
     -                     q2*
     -                      (-16*s**4 + (-48 + Born)*s**3*t1 + 
     -                        2*(-48 + Born)*s**2*t1**2 + 
     -                        (-48 + Born)*s*t1**3 - 16*t1**4)) + 
     -                  mt**2*t1**2*
     -                   (-128*q2**6*(s + t1) + 
     -                     24*s**3*t1**2*(s + t1)**2 + 
     -                     16*q2**5*(28*s**2 + 53*s*t1 + 23*t1**2) + 
     -                     q2**4*
     -                      (-560*s**3 + (-1452 + Born)*s**2*t1 + 
     -                        (-1240 + Born)*s*t1**2 - 284*t1**3) + 
     -                     2*q2*s*t1*(s + t1)*
     -                      (16*s**3 - (38 + Born)*s**2*t1 + 
     -                        (-46 + Born)*s*t1**2 - 32*t1**3) + 
     -                     2*q2**3*
     -                      (144*s**4 - (-552 + Born)*s**3*t1 + 
     -                        (550 - 3*Born)*s**2*t1**2 - 
     -                        2*(-102 + Born)*s*t1**3 + 46*t1**4) + 
     -                     q2**2*
     -                      (-48*s**5 + (-308 + Born)*s**4*t1 + 
     -                        (-368 + 3*Born)*s**3*t1**2 + 
     -                        4*(37 + Born)*s**2*t1**3 + 
     -                        2*(52 + Born)*s*t1**4 - 24*t1**5)) - 
     -                  mt**8*
     -                   (48*q2**4*(14*s**2 + 29*s*t1 + 15*t1**2) - 
     -                     16*q2**3*
     -                      (12*s**3 + 125*s**2*t1 + 186*s*t1**2 + 
     -                        73*t1**3) + 
     -                     q2*t1*(s + t1)*
     -                      (16*s**3 + 1544*s**2*t1 + 
     -                        (1072 + Born)*s*t1**2 + 388*t1**3) + 
     -                     4*q2**2*
     -                      (72*s**4 - 252*s**3*t1 - 515*s**2*t1**2 - 
     -                        217*s*t1**3 - 34*t1**4) + 
     -                     4*t1*
     -                      (-4*s**5 - 115*s**4*t1 + 12*s**3*t1**2 + 
     -                        160*s**2*t1**3 + 24*s*t1**4 - 13*t1**5))!!!!!!!!!!!!:\\
     -                   + mt**6*
     -                   (48*q2**5*(4*s**2 + 13*s*t1 + 9*t1**2) - 
     -                     8*q2**4*
     -                      (48*s**3 + 275*s**2*t1 + 416*s*t1**2 + 
     -                        189*t1**3) - 
     -                     2*q2*t1*(s + t1)*
     -                      (-8*s**4 + 50*s**3*t1 - 
     -                        (-600 + Born)*s**2*t1**2 + 
     -                        (318 + Born)*s*t1**3 + 46*t1**4) + 
     -                     q2**2*t1*
     -                      (-760*s**4 + 8*s**3*t1 + 
     -                        (2036 + 3*Born)*s**2*t1**2 + 
     -                        (1240 + 3*Born)*s*t1**3 + 100*t1**4) + 
     -                     4*q2**3*
     -                      (48*s**4 + 340*s**3*t1 + 675*s**2*t1**2 + 
     -                        647*s*t1**3 + 240*t1**4) + 
     -                     8*t1**2*
     -                      (7*s**5 + 40*s**4*t1 + 16*s**3*t1**2 - 
     -                        23*s**2*t1**3 - 5*s*t1**4 + t1**5)) - 
     -                  mt**4*t1*
     -                   (96*q2**6*(s + t1) - 
     -                     8*q2**5*(70*s**2 + 163*s*t1 + 93*t1**2) + 
     -                     4*q2**4*
     -                      (216*s**3 + 668*s**2*t1 + 735*s*t1**2 + 
     -                        259*t1**3) - 
     -                     q2*s*t1*(s + t1)*
     -                      (40*s**3 - (160 + Born)*s**2*t1 + 
     -                        4*(-119 + Born)*s*t1**2 - 
     -                        (248 + Born)*t1**3) - 
     -                     4*s*t1**2*
     -                      (16*s**4 + 39*s**3*t1 + 16*s**2*t1**2 - 
     -                        9*s*t1**3 - 2*t1**4) - 
     -                     q2**3*
     -                      (432*s**4 + 2072*s**3*t1 + 
     -                        (2344 - 3*Born)*s**2*t1**2 - 
     -                        3*(-404 + Born)*s*t1**3 + 348*t1**4) + 
     -                     q2**2*
     -                      (32*s**5 + 712*s**4*t1 + 652*s**3*t1**2 - 
     -                        (856 + 5*Born)*s**2*t1**3 - 
     -                        (624 + 5*Born)*s*t1**4 + 12*t1**5)))) - 
     -            2*(11*CA - 2*nf)*
     -             (mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))*t1*
     -             (mt**2 + t1)*(mt**2 - q2 + t1)*(s + t1)*
     -             (2*mt**6*(s + t1) + 
     -               2*q2*t1*(2*q2**2 + s**2 + t1**2 - 2*q2*(s + t1)) + 
     -               2*mt**4*(q2*s + t1*(s + t1)) + 
     -               mt**2*(2*q2*t1*(s + t1) + t1*(s + t1)**2 - 
     -                  2*q2**2*(2*s + 3*t1)))*Log(mu2/mt**2)))/
     -        (12.*Pi**2*q2*s*(mt**4 + (q2 - s)**2 - 2*mt**2*(q2 + s))*
     -          t1**3*(mt**2 + t1)*(mt**2 - q2 + t1)*(s + t1)))/
     -  (64.*CF*NC**2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Now, I need to add finite parts arising from different prefactor:
      Vep2=-(CF+CA)
      Vep =-(5./2.*CF + (11*CA-2*nf)/6. +(CA-2*CF)*lnu -CA*lns -CA*lnt)

      virtual_ep0=virtual_ep0
     $     + Vep2 *(log(mu2/m2t))**2 /2.  *born_pow *(st_alpha/(2*pi))
     $     + Vep  *log(mu2/m2t)           *born_pow *(st_alpha/(2*pi))

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccc
c     assign output
      res_virtual=virtual_ep0 /(st_alpha/(2*pi)) *ckmtmp
ccccccccccccccccccccccccccccccc


      end


