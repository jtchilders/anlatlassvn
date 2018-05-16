cccccccccccccccccccccccccccccccccccccccc
c     Total was checked using bambo PS generator.
c     Agreement level: (bambo-mine)/bambo = 0.00014
c     Error: err/bambo=err/mine=0.00013
c     using psgen=0 (and cthgen=0)
cccccccccccccccccccccccccccccccccccccccc
      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      double precision xborn(ndiminteg-3)

      integer psgen,cthgen
      logical ini
      data ini/.true./
      save ini,psgen,cthgen

      double precision xx(ndiminteg-3),tmp,m2t,m2b,mt,mb,xplus,xminus
      double precision tau,tau_min,tau_max
      double precision ycm,ycm_min,ycm_max
      double precision s45,s45min,s45max
      double precision cth2,cth3,phi2,z
      double precision jacb,shat,s_had
      double precision beta,vec(3),pb(0:3),pj(0:3),tmom,bmom,ptot(0:3)

      integer k,ileg,mu

      integer ixx_tau,ixx_ycm,ixx_s45,ixx_cth3,ixx_cth2,ixx_phi2
      parameter(
     $     ixx_tau  =1,
     $     ixx_ycm  =2,
     $     ixx_s45  =3,
     $     ixx_cth3 =4,
     $     ixx_cth2 =5,
     $     ixx_phi2 =6)
      double precision lambda,dotp
      external lambda,dotp

ccccccccccccccccccccccccccccccc
c     to check with bambo PS generator
      integer n
      real * 8 et,xm(3),p(4,3),xpar(5),wt
      logical check_bambo
      parameter (check_bambo=.false.)
ccccccccccccccccccccccccccccccc

      if(ini) then
c     sanity check
         if(check_bambo) then
            write(*,*) 'INTEGRATING USING BAMBO PS GENERATOR'
            if((ndiminteg-3).ne.7) then
               write(*,*) 'Problem in Born_phsp, bambo'
               call exit(1)
            endif
         else
            if((ndiminteg-3).ne.6) then
               write(*,*) 'Problem in Born_phsp'
               call exit(1)
            endif
         endif
c     Parameter to generate phase space
c     psgen=0:     flat in 1/tau
c     psgen=1:     flat in tau
         psgen=0
         if(psgen.eq.0) then
            write(*,*) '********************'
            write(*,*) ' 1/tau Imp sampling '
            write(*,*) '********************'
         endif
c     Parameter to generate phase space
c     cthgen=0: cosines flat (used during all tests)
c     cthgen=1: imp. sampled cosines
         cthgen=1
         if(cthgen.eq.1) then
            write(*,*) '*********************'
            write(*,*) ' cthgen Imp sampling '
            write(*,*) '*********************'
         endif
         do ileg=1,nlegborn
            kn_masses(ileg)=0
         enddo
         kn_masses(nlegreal)=0
         ini=.false.
      endif

c     t is massive
      kn_masses(3)=topmass_pow
c     b is massive
      kn_masses(4)=bmass_pow

c     local copy of variables
c     xx(1) -> tau
c     xx(2) -> ycm
      do k=1,ndiminteg-3
         xx(k)=xborn(k)
      enddo

      s_had=kn_sbeams
      jacb=1d0

      m2b=bmass_pow**2
      mb=bmass_pow
      m2t=topmass_pow**2
      mt=topmass_pow

      tau_min=(mt+mb)**2/s_had
      tau_max=1d0

cccccccccccc
c     tau
cccccccccccc
      if(psgen.eq.0)then
c     imp. sampling (flat in 1/tau )
         tmp=1d0/tau_max + (1d0/tau_min-1d0/tau_max)*xx(ixx_tau)
         tau=1d0/tmp
         jacb=jacb * tau**2 * (1d0/tau_min-1d0/tau_max)
      elseif(psgen.eq.1) then
c     uniform generation
         tau=tau_min + (tau_max-tau_min)*xx(ixx_tau)
         jacb=jacb * (tau_max-tau_min)
      else
         write(*,*) 'Wrong psgen in gen_born_vars'
         call exit(1)
      endif

cccccccccccc
c     ycm
cccccccccccc
      ycm_min=  log(tau)/2
      ycm_max= -log(tau)/2
      ycm = ycm_min + xx(ixx_ycm)*(ycm_max-ycm_min)
      jacb=jacb * (ycm_max-ycm_min)

      shat=tau*s_had

      if(check_bambo) goto 1111 ! now I have shat and ycm, so I can call bambo


cccccccccccc
c     b-j virtuality
cccccccccccc
      s45min=(0.+mb)**2
      s45max=(sqrt(shat) - mt)**2
      s45= s45min + xx(ixx_s45)*(s45max-s45min)
      jacb=jacb * (s45max-s45min)

cccccccccccc
c     top direction in CM frame
c     chose top momentum along x axis
c     azimuthal rotation at the very end, before event
c     generation
cccccccccccc
      if(cthgen.eq.0) then
         cth3= -1. + xx(ixx_cth3) *2.
         jacb=jacb*2.
      elseif(cthgen.eq.1) then
         z=1-2*xx(ixx_cth3)
         jacb=jacb*2
         cth3=1.5d0*(z-z**3/3)
         jacb=jacb*1.5d0*(1-z**2)
      else
         write(*,*) 'Wrong cthgen in gen_born_vars'
         call exit(1)
      endif
         
cccccccccccc
c     angles of bottom in frame where bottom and jet are
c     back-to-back
cccccccccccc
      if(cthgen.eq.0) then
         cth2= -1. + xx(ixx_cth2) *2.
         jacb=jacb*2.
      elseif(cthgen.eq.1) then
         z=1-2*xx(ixx_cth2)
         jacb=jacb*2
         cth2=1.5d0*(z-z**3/3)
         jacb=jacb*1.5d0*(1-z**2)
      endif
      phi2= 2. * pi *xx(ixx_phi2)
      jacb=jacb*2.*pi

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     physical phase space a la Byckling-Kajantie
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     all the pi factors
      jacb=jacb * (2.*pi)**4 / (2.*pi)**9  
      jacb=jacb * (2.*pi) ! *( d_phi/(2d0*pi) )
c     jacobian 1
      jacb=jacb*sqrt(lambda(shat,s45,m2t))/8./shat
c     jacobian 2
      jacb=jacb*sqrt(lambda(s45,0d0,m2b))/8./s45
      
      kn_jacborn=jacb
      kn_sborn=shat

ccccccccccccccccccccccccccccccccc
c     Here I check with bambo
c     Notice that I need one more random number with bambo.
c     Therefore I need to add 1 to ndiminteg when I do this check.
      if(.not.check_bambo) goto 3333
 1111 continue
      xm(1)=topmass_pow
      xm(2)=bmass_pow
      xm(3)=0.
      do ileg=1,5
         xpar(ileg)=xborn(ileg+2)
      enddo
      call bambo(3,sqrt(shat),xm,xpar,p,wt)
      kn_jacborn=wt*jacb
      kn_sborn=shat
c     Feynman x's
      xplus=sqrt(tau) * exp(ycm)
      xminus=sqrt(tau) * exp(-ycm)
      kn_xb1=xplus
      kn_xb2=xminus
c     initial state particles
      kn_cmpborn(0,1)=sqrt(shat)/2
      kn_cmpborn(0,2)=kn_cmpborn(0,1)
      kn_cmpborn(3,1)=kn_cmpborn(0,1)
      kn_cmpborn(3,2)=-kn_cmpborn(0,2)
      kn_cmpborn(1,1)=0
      kn_cmpborn(1,2)=0
      kn_cmpborn(2,1)=0
      kn_cmpborn(2,2)=0  
      do ileg=3,5
c!         print*, sqrt(dabs(p(4,ileg-2)**2-p(1,ileg-2)**2-p(2,ileg-2)**2-p(3,ileg-2)**2))
         kn_cmpborn(0,ileg)=p(4,ileg-2)
         do mu=1,3
            kn_cmpborn(mu,ileg)=p(mu,ileg-2)
         enddo
      enddo
      goto 2222 !now I can boost all momenta back to lab frame
ccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccc
c     compute momenta
cccccccccccccccccccccccccccccc
 3333 continue

c     Feynman x's
      xplus=sqrt(tau) * exp(ycm)
      xminus=sqrt(tau) * exp(-ycm)
      kn_xb1=xplus
      kn_xb2=xminus

c     initial state particles
      kn_cmpborn(0,1)=sqrt(shat)/2
      kn_cmpborn(0,2)=kn_cmpborn(0,1)
      kn_cmpborn(3,1)=kn_cmpborn(0,1)
      kn_cmpborn(3,2)=-kn_cmpborn(0,2)
      kn_cmpborn(1,1)=0
      kn_cmpborn(1,2)=0
      kn_cmpborn(2,1)=0
      kn_cmpborn(2,2)=0  


c     top in CM frame
      tmom=sqrt(lambda(shat,s45,m2t)/4./shat)
      kn_cmpborn(1,3)= tmom * sqrt(1.-cth3**2)
      kn_cmpborn(2,3)= 0.
      kn_cmpborn(3,3)= tmom * cth3
      kn_cmpborn(0,3)= sqrt(tmom**2 + m2t)

c     b-j kinematics in b-j rest frame
      bmom=sqrt(lambda(s45,0d0,m2b)/4./s45)
      pb(1)= bmom * sqrt(1.-cth2**2) * sin(phi2)
      pb(2)= bmom * sqrt(1.-cth2**2) * cos(phi2)
      pb(3)= bmom * cth2
      pb(0)= sqrt(bmom**2 + m2b)
      do k=1,3
         pj(k)=-pb(k)
      enddo
      pj(0)=bmom

c$$$      do k=0,3
c$$$         print*, (pb(k)+pj(k))/sqrt(s45)
c$$$      enddo

c     now boost pb and pj momenta in CM frame
c     (pb+pj)_cm = (sqrt(s)-E_t , -pt_x , -pt_y , -pt_z)
      beta=tmom/(sqrt(shat)-kn_cmpborn(0,3))
      do k=1,3
         vec(k)=-kn_cmpborn(k,3)/tmom
      enddo
      call mboost(1,vec,beta,pb,kn_cmpborn(0,4))
      call mboost(1,vec,beta,pj,kn_cmpborn(0,5))

 2222 continue

c$$$      do mu=0,3
c$$$         ptot(mu)=0.
c$$$         do ileg=3,5
c$$$            ptot(mu)=ptot(mu)+kn_cmpborn(mu,ileg)
c$$$         enddo
c$$$         print*,ptot(mu)/sqrt(shat)
c$$$      enddo

      call checkmomzero(nlegborn,kn_cmpborn)

c now boost everything along 3
      beta=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn-2,vec,beta,kn_cmpborn(0,3),kn_pborn(0,3))
      do mu=0,3
         kn_pborn(mu,1)=kn_xb1*kn_beams(mu,1)
         kn_pborn(mu,2)=kn_xb2*kn_beams(mu,2)
      enddo
      call checkmomzero(nlegborn,kn_pborn)

c      print*, dotp(kn_pborn(0,3),kn_pborn(0,3))

      kn_minmass=topmass_pow+bmass_pow
      end


      function lambda(x,y,z)
      double precision lambda,x,y,z
      lambda = x**2 + y**2 + z**2 -2.*x*y -2.*x*z -2.*y*z
      end

      subroutine born_suppression(fact)
      implicit none
      logical ini
      data ini/.true./
      real * 8 fact
      if (ini) then
         write(*,*) '**************************'
         write(*,*) 'No Born suppression factor'
         write(*,*) '**************************'
         ini=.false.
      endif
      fact=1d0
      end


      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      real * 8 muf,mur
      logical ini
      data ini/.true./
      real *8 muref
      real *8 dotp
      real *8 pt2b,mtb
      external dotp
      if (ini) then
c$$$         write(*,*) '*************************************'
c$$$         write(*,*) '    Factorization and renormalization '
c$$$         write(*,*) '    scales set to top mass         '
c$$$         write(*,*) '    (renscfact,facscfact) = ',st_renfact,st_facfact
c$$$         write(*,*) '*************************************'
         write(*,*) '*************************************'
         write(*,*) '    Factorization and renormalization '
         write(*,*) '    scales set to 4*(mb^2+ptb^2)'
         write(*,*) '    (renscfact,facscfact) = ',st_renfact,st_facfact
         write(*,*) '*************************************'
         ini=.false.
      endif
c$$$      muf=topmass_pow
c$$$      mur=topmass_pow
      pt2b=kn_pborn(1,4)**2+kn_pborn(2,4)**2
      mtb=sqrt( bmass_pow**2 + pt2b )
      muf=4d0*mtb
      mur=4d0*mtb
      end



ccccccccccccccccccccccccccccccccccccccccccccccccc
c     In the following, routines from FNO paper
ccccccccccccccccccccccccccccccccccccccccccccccccc

c Maps 3n-4 dimensional unit cube into n-particle phase space,
c and computes jacobian.
c Phase space is normalized in standard way (as (2.8) in fno2006)
c
c n: number of particles
c et: CM energy
c xm(n): masses
c xpar(3n-4): the 3n-4 parameters between zero and 1 that parametrize the phase space
c p(4,n): output 5 vectors: px py pz energy
c wt: output weight, wt d xpar_1 ... d xpar_{3n-4} is the phase space element
      subroutine bambo(n,et,xm,xpar,p,wt)
      implicit none
      integer n
      real * 8 et,xm(n),p(4,n),xpar(3*n-4),wt
      integer m,ipar,j,nn
      real * 8 esys,beta,vec(3),recm(4),mrec2,mrec,mrmax,mrmin,mr,
     #   mrecmin,x,xmin,xmax,phi,cth,sth,km,em,xjac,tmp
      real * 8 pi
      parameter (pi=3.141592653589793d0)

      if(n.lt.2) then
         write(*,*) ' error: cannot do 1 body phase space'
         stop
      endif

      xjac=1
c initial energy of subsystem
      esys=et
c initial boost parameters of subsystem
      beta=0
      vec(1)=1
      vec(2)=0
      vec(3)=0

c phasespace parametrized by xpar(ipar)
      ipar=1
c minimum mass of subsystem
      mrecmin=0
      do j=1,n
         mrecmin=mrecmin+xm(j)
      enddo

      do m=1,n-1
c mass of recoil system
         mrecmin=mrecmin - xm(m)
c     Recoil sistem must have mass >= sum j=2,m-1 m_j
         if(et.lt.xm(m)+mrecmin) then
            wt=0
            return
         endif
c     generate uniform 2-body phase space as if recoil system was a single
c     particle of mass from mrec to maximum, scaling as if massless
c     n-m phase space: mrec^(2*nn-2)/(2*nn-2) - mrec^(2*nn) Esys^(-2) /(2*nn)
         if(m.lt.n-1) then
            nn=(n-m)
c mr is mrec^2/esys^2
            mrmax = ((esys-xm(m))/esys)**2

c 2 body massless phase space is proportional to
c (E^2-mrec^2)*mrec^(2nn-4) d mrec^2
c      \propto Esys^2 mrec^(2*nn-2)/(2*nn-2) - mrec^(2*nn) /(2*nn)

            xmax = mrmax**(nn-1)/(nn-1)-mrmax**nn/nn
            mrmin = (mrecmin/esys)**2
            xmin = mrmin**(nn-1)/(nn-1)-mrmin**nn/nn
            x = xmax*xpar(ipar)+(1-xpar(ipar))*xmin
            xjac=xjac*(xmax-xmin)
c solves m^(n-1)/(n-1)-m^(n) / n =x
            call solvespec(nn,x,mr)
c jacobian from dx to d mr:
            xjac = xjac/(mr**(nn-2) - mr**(nn-1))
            mrec2=mr*esys**2
            xjac=xjac*esys**2
c divide by 2 pi, when introducing 2 pi delta(prec^2-mrec^2) d mrec^2/(2 pi)
            xjac=xjac/(2 * pi)
            mrec=sqrt(mrec2)
            ipar=ipar+1
         else
            mrec=mrecmin
         endif
         phi=xpar(ipar)*2*pi
         xjac=xjac*2*pi
         ipar=ipar+1
         cth=1-xpar(ipar)*2
         ipar=ipar+1
         xjac=xjac*2
         sth=sqrt(1-cth**2)
c energy of mth particle
         em=(esys**2-mrec**2+xm(m)**2)/(2*esys)
c its momentum
         km=sqrt(em**2-xm(m)**2)
c     include two body phase space factor
c     
         xjac=xjac*km/esys /(16*pi**2)
c     
         p(3,m)=km*cth
         p(1,m)=km*sth*cos(phi)
         p(2,m)=km*sth*sin(phi)
         p(4,m)=em
c     recoil momentum
         recm(1)=-p(1,m)
         recm(2)=-p(2,m)
         recm(3)=-p(3,m)
         recm(4)=esys-p(4,m)
c     Boost p_m and recm according to beta of current system
         call mboost_14(1,vec,beta,p(1,m),p(1,m))
         call mboost_14(1,vec,beta,recm,recm)
         if(m.eq.n-1) then
c nothing else to be done! last momentum is recoil momentum
            p(1,n)=recm(1)
            p(2,n)=recm(2)
            p(3,n)=recm(3)
            p(4,n)=recm(4)
         else
c     Compute velocity of recoil system; at next iteration
c     we compute the phase space of recoil system, so esys=mrec.
            tmp=sqrt(recm(1)**2+recm(2)**2+recm(3)**2)
            vec(1)=recm(1)/tmp
            vec(2)=recm(2)/tmp
            vec(3)=recm(3)/tmp
            beta=tmp/recm(4)
            esys=mrec
         endif
      enddo
      if(ipar.ne.3*n-3) then
         write(*,*) '?'
      endif
      wt=xjac
      end


      subroutine solvespec(n,x,m)
      implicit none
      real * 8 x,m,m0,dm0,vm0
      integer n
c solves the equation
c     m^(n-1)/(n-1)-m^(n) / n =x
c for integer n
      if(x.lt.0.or.x.gt.1d0/(n-1)-1d0/n) then
         m=-1
         return
      endif
      if(n.eq.2) then
         m=1-sqrt(1-2*x)
      else
         m0=0.5d0
 1       continue
         dm0=m0**(n-2)-m0**(n-1)
         vm0=m0**(n-1)/(n-1)-m0**n/n
c solve vm0+(m-m0)*dm0=x
         m=(x-vm0)/dm0+m0
         if(m.gt.1) then
            m=1
         elseif(m.lt.0) then
            m=0
         endif
         if(abs(m-m0)/(abs(m)+abs(m0)).lt.1d-13) then
            return
         else
            m0=m
            goto 1
         endif
      endif
      end



      subroutine mboost_14(m,vec,beta,vin,vout)
c     boosts the m vectors vin(4,m) into the vectors vout(4,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (x,y,z,t).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(4,m),vout(4,m)
      real * 8 betav,gamma
      real * 8 vdotb
      integer ipart,idim
      gamma=1/sqrt(1-beta**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*vec(1)
     #         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
         do idim=1,3
            vout(idim,ipart)=vin(idim,ipart)
     #           +vec(idim)*((gamma-1)*vdotb
     #           +gamma*beta*vin(4,ipart))
         enddo
         vout(4,ipart)=gamma*(vin(4,ipart)+vdotb*beta)
      enddo
      end



