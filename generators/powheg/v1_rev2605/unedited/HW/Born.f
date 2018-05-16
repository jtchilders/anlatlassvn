c     computation of the born amplitude (result from
c     maple+average over initial spin and colour)
      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'      
      real *8 p(0:3,nlegborn),born,bornjk(nlegborn,nlegborn)
      real *8 bmunu(0:3,0:3,nlegborn)
      integer bflav(nlegborn)
      integer i,j,k,mu,nu,n
      real * 8 amp2,bv,ba,gw,mw,props
      real *8 p13,p14,p23,p24,p34,p12,p33
      real *8 dotp
      external dotp

      amp2 = 0d0
      n=3d0 !number of initial colour
c
      gw=ph_unit_e/ph_sthw
      mw=ph_wmass
c
      p12=dotp(p(0,1),p(0,2))
      p13=dotp(p(0,1),p(0,3))
      p14=dotp(p(0,1),p(0,4))
      p23=dotp(p(0,2),p(0,3))
      p24=dotp(p(0,2),p(0,4))
      p33=dotp(p(0,3),p(0,3))
      p34=dotp(p(0,3),p(0,4))
c
c     vectorial and axial parts
      bv = -256*p24*p14-128*p24*p13+128*p24*p12-128*p14*p23+128*p14*p12
      ba = -128*p24*p13+128*p24*p12+128*p14*p23-128*p14*p12
c
c     W- PRODUCTION
      if (bflav(4).gt.0) then
c     case one: q(p1),qbar(p2) -> H W- 
         if (bflav(1).gt.0) then
            amp2=bv+ba          
c     case two: qbar(p1),q(p2) -> H W- 
         elseif (bflav(1).lt.0) then
            amp2=bv-ba        
         endif
c     W+ PRODUCTION
      elseif (bflav(4).lt.0) then
c     case one: q(p1),qbar(p2) -> H W+ 
         if (bflav(1).gt.0) then
            amp2=bv-ba         
c     case two: qbar(p1),q(p2) -> H  W+ 
         elseif (bflav(1).lt.0) then
            amp2=bv+ba     
         endif
      endif

c     average over initial spins and colours
      amp2=amp2*n/4d0/n**2
c
c     W propagators
      props = 1/((2*p12-mw**2)**2+ph_Wwidth**2*mw**2)/((p33-2*p13-2*p23+
     #2*p12-mw**2)**2+ph_Wwidth**2*mw**2)

      born=amp2*props
c     CKM matrix
      i=bflav(1)
      j=bflav(2)
      if (mod(i,2).eq.0) then
         born=born*ph_CKM(abs(i)/2,(abs(j)+1)/2)**2
      elseif (mod(j,2).eq.0) then
         born=born*ph_CKM((abs(i)+1)/2,abs(j)/2)**2
      endif
c
c     coupling constants and W mass
c     factor gw^2/8 from each weak vertex: two vertices, gw^4/64
c     factor 4mw^4/v^2 from Higgs vertex; but v^2 = 4mw^2/gw^2
c     from Higgs vertex: mw^2 gw^2
      born=born* (gw**2/8d0) * (gw**2/8d0) * (mw**2*gw**2)
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

c     spin correlated born amplitude
      do j=1,nlegborn
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,j)=0d0
            enddo
         enddo
      enddo
      end



c$$$      subroutine pconj(p,n)
c$$$      implicit none
c$$$      real * 8 p(0:3,n)
c$$$      integer n,j,mu
c$$$      do j=1,n
c$$$         do mu=1,3
c$$$            p(mu,j)=-p(mu,j)
c$$$         enddo
c$$$      enddo
c$$$      end
c$$$
c$$$      subroutine borncolour_lh
c$$$c Wrapper subroutine to call the MadGraph code to associate
c$$$c a (leading) color structure to an event.
c$$$      implicit none
c$$$      include 'nlegborn.h'
c$$$      include 'LesHouches.h'
c$$$      integer bflav(nlegborn),color(2,nlegborn)
c$$$      integer i,j,itmp
c$$$
c$$$      integer idvecbos,vdecaymode
c$$$      common/cvecbos/idvecbos,vdecaymode
c$$$
c$$$      do i=1,nlegborn
c$$$         bflav(i)=idup(i)
c$$$         if (bflav(i).eq.21) bflav(i)=0
c$$$      enddo
c$$$
c$$$      if(idvecbos.eq.-24) then
c$$$         call cconj(bflav,nlegborn)
c$$$      endif
c$$$
c$$$      call born_color(bflav,color)
c$$$
c$$$c      write(*,*) 'color'
c$$$c      write(*,*) (color(1,i),i=1,nlegborn)
c$$$c      write(*,*) (color(2,i),i=1,nlegborn)
c$$$
c$$$      if(idvecbos.eq.-24) then
c$$$         call cconj(bflav,nlegborn)
c$$$         do j=1,nlegborn
c$$$            itmp = color(1,j)
c$$$            color(1,j) = color(2,j)
c$$$            color(2,j) = itmp
c$$$         enddo
c$$$      endif
c$$$
c$$$      do i=1,2
c$$$         do j=1,nlegborn
c$$$            icolup(i,j)=color(i,j)
c$$$         enddo
c$$$      enddo
c$$$
c$$$      call borncolour_lh_tmp      
c$$$
c$$$      do i=1,2
c$$$         do j=1,nlegborn     
c$$$            if (color(i,j).ne.0) then
c$$$               write(*,*) icolup(i,j)/color(i,j)
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$
c$$$
c$$$      end
c$$$


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
       implicit none
       integer id5,id6
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
      call momenta_reshuffle(4,5,6,decmass)

c     fix here the W decay mode
      id5=vdecaymode
      id6=-vdecaymode + sign(1,idvecbos) 
      call change_id_particles(5,6,id5,id6)

      end



c     i1<i2
      subroutine momenta_reshuffle(ires,i1,i2,decmass)
      implicit none
      include 'LesHouches.h'
      integer ires,i1,i2,j
      real * 8 ptemp(0:3),ptemp1(0:3),beta(3),betainv(3),modbeta,decmass
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
cccccccccccccccccccccccccccccccccccccccc
c first decay product (massive)
      ptemp(0)=pup(4,i1)
      do j=1,3
         ptemp(j)=pup(j,i1)
      enddo
      call mboost(1,beta,modbeta,ptemp,ptemp)
      ptemp1(0)=0.5d0*(pup(5,ires)+(decmass**2)/pup(5,ires))
      do j=1,3
         ptemp1(j)=ptemp(j)/ptemp(0)*sqrt(ptemp1(0)**2 -decmass**2)
      enddo
      call mboost(1,betainv,modbeta,ptemp1,ptemp)
      do j=1,3
         pup(j,i1)=ptemp(j)
      enddo
      pup(4,i1)=ptemp(0)
      pup(5,i1)=sqrt(pup(4,i1)**2-pup(1,i1)**2
     $     -pup(2,i1)**2-pup(3,i1)**2)
      
c second decay product (massless)

      ptemp(0)=pup(4,i2)
      do j=1,3
         ptemp(j)=pup(j,i2)
      enddo
      call mboost(1,beta,modbeta,ptemp,ptemp)
      ptemp1(0)=0.5d0*(pup(5,ires)-(decmass**2)/pup(5,ires))
      do j=1,3
         ptemp1(j)=ptemp(j)/ptemp(0)*ptemp1(0)
      enddo
      call mboost(1,betainv,modbeta,ptemp1,ptemp)
      do j=1,3
         pup(j,i2)=ptemp(j)
      enddo
      pup(4,i2)=ptemp(0)
c abs to avoid tiny negative values
      pup(5,i2)=sqrt(abs(pup(4,i2)**2-pup(1,i2)**2
     $     -pup(2,i2)**2-pup(3,i2)**2))
cccccccccccccccccccccccccccccccccccccccc
      end



      subroutine change_id_particles(i1,i2,id1,id2)
      implicit none
      include 'LesHouches.h'
      integer i1,i2,id1,id2
      idup(i1)=id1
      idup(i2)=id2
      end


