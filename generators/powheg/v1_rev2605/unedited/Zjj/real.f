      subroutine setreal(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include '../include/pwhg_flst.h'
      include '../include/pwhg_math.h'
      include '../include/pwhg_st.h'
      include '../include/pwhg_kn.h' !:
      double precision p(0:3,nlegreal),amp2
      integer rflav(nlegreal)

      double precision kr_mad(0:3,nlegreal),amp2mad
      integer ileg,mu
      logical ini
      data ini/.true./
      save ini

      double precision amp2mcfm
      integer ngluons
      external ngluons
      integer gluonpos(3)


      if(ini) then 
c     setting madgraph parameters (needed for madgraph subroutines)
         call set_madgraph_parameters
         ini=.false.
      endif
c     set madgraph parameters that can change on an event-by-event basis
      call mad_setparam
      
      do ileg=1,nlegreal
         do mu=0,3
            kr_mad(mu,ileg)=p(mu,ileg)
         enddo
      enddo
c     to avoid bugs in HELAS, restore exact masslessness of  incoming partons 
      kr_mad(0,1)=dabs(kr_mad(3,1))
      kr_mad(0,2)=dabs(kr_mad(3,2))
      call compreal(kr_mad,rflav,amp2mad)


ccccccccccccccccccccccccccccccccccc
c$$$      call setrealmcfm(p,rflav,amp2mcfm)
c$$$
c$$$      if(ngluons(rflav,gluonpos).eq.3) then
c$$$
c$$$      if(rflav(1).eq.0.and.rflav(2).eq.0.and.rflav(5)*rflav(6)*rflav(7).eq.0) then
c$$$         write(*,*) '>>> SHOULD BE EQUAL >>> ',amp2mcfm/amp2mad
c$$$      endif
c$$$      if(rflav(5).eq.0.and.rflav(6).eq.0.and.rflav(7).eq.0) then
c$$$         write(*,*) '>>> SHOULD BE EQUAL >>> ',amp2mcfm/amp2mad
c$$$      endif
c$$$      if(rflav(2).eq.0.and.rflav(5)*rflav(6)*rflav(7).eq.0) then
c$$$         write(*,*) '>>> SHOULD BE EQUAL >>> ',amp2mcfm/amp2mad
c$$$      endif
c$$$      if(rflav(1).eq.0.and.rflav(5)*rflav(6)*rflav(7).eq.0) then
c$$$         write(*,*) '>>> SHOULD BE EQUAL >>> ',amp2mcfm/amp2mad
c$$$      endif
c$$$
c$$$         if(dabs(amp2mcfm/amp2mad-1.).gt.0.01) then
c$$$            write(*,*) 'Should be 1 ',amp2mcfm/amp2mad, rflav,kn_csitilde,kn_y
c$$$c            stop
c$$$         endif
c$$$
c$$$
c$$$      elseif(ngluons(rflav,gluonpos).eq.1) then

c$$$         if(rflav(1).gt.0.and.rflav(2).lt.0) then
c$$$            write(*,*) '>>> SHOULD BE EQUAL >>> ',amp2mcfm/amp2mad
c$$$         endif
c$$$         if(rflav(1).lt.0.and.rflav(2).gt.0) then
c$$$            write(*,*) '>>> SHOULD BE EQUAL >>> ',amp2mcfm/amp2mad
c$$$         endif
c$$$         if(rflav(1).ne.0.and.rflav(2).eq.0) then
c$$$            write(*,*) '>>> SHOULD BE EQUAL >>> ',amp2mcfm/amp2mad
c$$$         endif
c$$$         if(rflav(1).eq.0.and.rflav(2).ne.0) then
c$$$            write(*,*) '>>> SHOULD BE EQUAL >>> ',amp2mcfm/amp2mad
c$$$         endif
c$$$         if(rflav(1).gt.0.and.rflav(2).gt.0) then
c$$$            write(*,*) '>>> SHOULD BE EQUAL >>> ',amp2mcfm/amp2mad
c$$$         endif
c$$$         if(rflav(1).lt.0.and.rflav(2).lt.0) then
c$$$            write(*,*) '>>> SHOULD BE EQUAL >>> ',amp2mcfm/amp2mad
c$$$         endif

c$$$         if(dabs(amp2mcfm/amp2mad-1.).gt.0.01) then
c$$$            write(*,*) 'Should be 1 ',amp2mcfm/amp2mad, rflav,kn_csitilde,kn_y
c$$$c            stop
c$$$         endif
c$$$
c$$$      else
c$$$         write(*,*) 'Error in setreal, ngluons'
c$$$         call exit(1)
c$$$      endif



ccccccccccccccccccccccccccccccccccc



      amp2=amp2mad

c      if(amp2.lt.0) amp2=amp2mad

c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2*pi))



      end



