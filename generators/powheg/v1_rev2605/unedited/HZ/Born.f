c     computation of the Born amplitude (result from
c     maple+average over initial spin and colour)

      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      real *8 p(0:3,nlegborn)
      integer bflav(nlegborn)
      real *8 bornjk(nlegborn,nlegborn),bmunu(0:3,0:3,nlegborn)
      integer i,j,k,mu,nu,n
      real * 8 born,amp2,gw,props,couplz,T3L,T3Q,chargeQ,chargeL
      real * 8 p13,p14,p23,p24,p34,p12,p45,VQ,AQ,VL,AL,t0
      real * 8 dotp
      real * 8 pp(0:3,nlegborn)
      external dotp
c
      amp2 = 0d0
      n=3d0 !number of colours
      gw=ph_unit_e/ph_sthw
      ph_cthw = sqrt(1-ph_sthw**2)
      couplz = ph_unit_e/(2*ph_sthw*ph_cthw)
c
c     vectorial and axial couplings to Z boson
      if (mod(abs(bflav(4)),2).eq.1) then
c     LEPTON
         chargeL = -1
         T3L = -1d0/2d0
      elseif (mod(abs(bflav(4)),2).eq.0) then
c     NEUTRINO
         chargeL = 0
         T3L = 1d0/2d0
      endif

      if (mod(abs(bflav(1)),2).eq.0) then
c     UP TYPE QUARK
         chargeQ = 2d0/3d0
         T3Q = 1d0/2d0
      elseif (mod(abs(bflav(1)),2).eq.1) then
c     DOWN TYPE QUARK
         chargeQ = -1d0/3d0
         T3Q = -1d0/2d0
      endif

      VL = T3L - 2*chargeL*ph_sthw**2
      AL = -T3L
      VQ = T3Q - 2*chargeQ*ph_sthw**2
      AQ = -T3Q
c
      pp = p

      if (bflav(4).lt.0)then
         call swap_momenta(pp(:,4),pp(:,5))
      endif

      if (bflav(1).lt.0.and.bflav(2).gt.0) then
         call swap_momenta(pp(:,1),pp(:,2))
      endif

      p12=dotp(pp(0,1),pp(0,2))
      p13=dotp(pp(0,1),pp(0,3))
      p14=dotp(pp(0,1),pp(0,4))
      p23=dotp(pp(0,2),pp(0,3))
      p24=dotp(pp(0,2),pp(0,4))
      p34=dotp(pp(0,3),pp(0,4))
      p45=dotp(pp(0,4),pp(0,5))
c
      t0 = ((-64*p24*p14-32*p24*p13+32*p14*p12-32*p14*p23+32*p24*p12)*VQ
     #**2+(-64*p24*p14-32*p24*p13+32*p14*p12-32*p14*p23+32*p24*p12)*AQ**
     #2)*VL**2+(128*p24*p12-128*p24*p13-128*p14*p12+128*p14*p23)*AL*AQ*V
     #Q*VL+(-64*p24*p14-32*p24*p13+32*p14*p12-32*p14*p23+32*p24*p12)*AL*
     #*2*VQ**2+(-64*p24*p14-32*p24*p13+32*p14*p12-32*p14*p23+32*p24*p12)
     #*AL**2*AQ**2
c
c     average over initial spins and colours
      amp2=t0*n/4d0/n**2
c
c     Z propagators
      props = 1/((2*p12-ph_Zmass2)**2+ph_ZmZw**2)/
     $     ((2*p45-ph_Zmass2)**2+ph_ZmZw**2)

      born=amp2*props
c
c     coupling constants and Z mass
c     factor couplz^2 from each weak vertex: two vertices, couplz^4
c     factor 4mz^4/v^2 from Higgs vertex; but v^2 = 4mw^2/gw^2
c     from Higgs vertex: mz^4 gw^2/mw^2
      born=born*(couplz)**2*(couplz)**2*(ph_Zmass2*gw)**2/ph_Wmass2
c     
c     initialization of bornjk
      do j=1,nlegborn
        do k=1,nlegborn
           bornjk(j,k)=0d0
        enddo
      enddo
c     colour correlated born amplitude:
      bornjk(1,1)=-CF*born
      bornjk(2,2)=-CF*born
      bornjk(1,2)=CF*born
      bornjk(2,1)=CF*born
c
c     spin correlated born amplitude
      do j=1,nlegborn
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,j)=0d0
            enddo
         enddo
      enddo
      end




      subroutine borncolour_lh     
c Sets up the colour for the given flavour configuration
c already filled in the Les Houches interface.
c In case there are several colour structure, one
c should pick one with a probability proportional to
c the value of the corresponding cross section, for the
c kinematics defined in the Les Houches interface
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
c     neutral particles
      icolup(1,3)=0
      icolup(2,3)=0
      icolup(1,4)=0
      icolup(2,4)=0
      icolup(1,5)=0
      icolup(2,5)=0
c     colored particles
      if((idup(1).gt.0).and.(idup(2).lt.0)) then
         icolup(1,1)=501
         icolup(2,1)=0
         icolup(1,2)=0
         icolup(2,2)=501
      elseif((idup(1).lt.0).and.(idup(2).gt.0)) then
         icolup(1,1)=0
         icolup(2,1)=501
         icolup(1,2)=501
         icolup(2,2)=0
      else
         write(*,*) ' invalid flavour'
         call pwhg_exit(-1)
      endif
      end


      subroutine finalize_lh
c     Set up the resonances whose mass must be preserved
c     on the Les Houches interface.
c     
c     vector boson id and decay
      integer idvecbos,vdecaymode
      common/cvecbos/idvecbos,vdecaymode
c     lepton masses
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass

      call add_resonance(idvecbos,4,5)
c     The following routine also performs the reshuffling of momenta if
c     a massive decay is chosen
      call momenta_reshuffle(4,5,6,decmass,decmass)

c     fix here the W decay mode
      id5=vdecaymode
      id6=-vdecaymode
      call change_id_particles(5,6,id5,id6)

      end




      subroutine change_id_particles(i1,i2,id1,id2)
      implicit none
      include 'LesHouches.h'
      integer i1,i2,id1,id2
      idup(i1)=id1
      idup(i2)=id2
      end



c     i1<i2
      subroutine momenta_reshuffle(ires,i1,i2,m1,m2)
      implicit none
      include 'LesHouches.h'
      integer ires,i1,i2
      real * 8 m1,m2
      real * 8 ptemp(0:3),pfin(0:3),beta(3),betainv(3),modbeta,m
      real * 8 mod_pfin,m0
      integer j,id,dec
      if (i1.ge.i2) then
         write(*,*) 'wrong sequence in momenta_reshuffle'
         stop
      endif
cccccccccccccccccccccccccccccc
c construct boosts from/to vector boson rest frame 
      do j=1,3
         beta(j)=-pup(j,ires)/pup(4,ires)
      enddo
      modbeta=sqrt(beta(1)**2+beta(2)**2+beta(3)**2)
      do j=1,3
         beta(j)=beta(j)/modbeta
         betainv(j)=-beta(j)
      enddo

      m0 = pup(5,ires)
      mod_pfin=
     $     1/(2*m0)*sqrt(abs((m0**2-m1**2-m2**2)**2 - 4*m1**2*m2**2))
               
cccccccccccccccccccccccccccccccccccccccc
c     loop of the two decay products
      
      do dec=1,2
         if(dec.eq.1) then
            id=i1
            m=m1
         else
            id=i2
            m=m2
         endif
         ptemp(0)=pup(4,id)
         do j=1,3
            ptemp(j)=pup(j,id)
         enddo
         call mboost(1,beta,modbeta,ptemp,ptemp)
         pfin(0)=sqrt(mod_pfin**2 + m**2)
         do j=1,3
            pfin(j)=ptemp(j)*mod_pfin/ptemp(0)
         enddo
         call mboost(1,betainv,modbeta,pfin,ptemp)
         do j=1,3
            pup(j,id)=ptemp(j)
         enddo
         pup(4,id)=ptemp(0)
         pup(5,id)=sqrt(abs(pup(4,id)**2-pup(1,id)**2
     $        -pup(2,id)**2-pup(3,id)**2))
         
      enddo

      end
