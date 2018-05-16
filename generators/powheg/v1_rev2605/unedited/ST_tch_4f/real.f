      subroutine setreal(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include '../include/pwhg_flst.h'
c      include '../include/pwhg_math.h'
      include '../include/pwhg_st.h'
      include 'PhysPars.h'
      double precision p(0:3,nlegreal),amp2
      integer rflav(nlegreal)
      double precision kr_mad(0:3,nlegreal),amp2mad,kr_mad_cm(0:3,nlegreal)
      integer ileg,imu
      double precision vec(3),beta
      logical ggtype,qqtype,gqtype
ccccccccccccccccccccccccccccccccccc
      include 'MCFM_include/constants.f'
      include 'MCFM_include/noglue.f'
      include 'MCFM_include/scale.f'
      include 'MCFM_include/qcdcouple.f'
      include 'MCFM_include/stopscales.f'
      include 'MCFM_include/masses.f'
      include 'MCFM_include/ckm.f'
      integer j,k
      double precision pmcfm(mxpart,1:4),rtmp
      double precision realmcfm(-nf:nf, -nf:nf)
      logical ini
      save ini
      data ini/.true./
      logical checkreal
      parameter (checkreal=.false.)
      if(checkreal) then
         noglue=.false.
         ggonly=.false.
         gqonly=.false.
c     put all ckm to 0 in mcfm
         do ileg=-nf,nf
            vsum(ileg)=0.
            do imu=-nf,nf
               Vsq(ileg,imu)=0.
            enddo
         enddo
c     as in ckmfill.f, when nwz=1
         Vsq(2,-1)=ph_CKM(1,1)**2
         Vsq(2,-3)=ph_CKM(1,2)**2
         Vsq(2,-5)=ph_CKM(1,3)**2
         Vsq(4,-1)=ph_CKM(2,1)**2
         Vsq(4,-3)=ph_CKM(2,2)**2
         Vsq(4,-5)=ph_CKM(2,3)**2
         do j=-nf,nf
            do k=-nf,nf
               Vsq(j,k)=Vsq(k,j)
            enddo
         enddo
         do j=1,5
            Vsum(+j)=Vsq(+j,-1)+Vsq(+j,-2)+Vsq(+j,-3)+Vsq(+j,-4)+Vsq(+j,-5)
            Vsum(-j)=Vsq(-j,+1)+Vsq(-j,+2)+Vsq(-j,+3)+Vsq(-j,+4)+Vsq(-j,+5)
         enddo
         Vsum(0)=0
c$$$  print*, 'VSUM ',vsum(1),vsum(2),vsum(3),vsum(4),vsum(5)
c$$$  print*, 'VSUM ',vsum(-1),vsum(-2),vsum(-3),vsum(-4),vsum(-5)
c$$$  stop
         if(ttype.eq.-1) then
            write(*,*) 'ttype=-1 while checking real'
            write(*,*) 'this was not tested'
            stop
         endif
      endif
ccccccccccccccccccccccccccccccccccc


c     set madgraph parameters that can change on an event-by-event basis
      call mad_setparam
      do ileg=1,nlegreal
         do imu=0,3
            kr_mad(imu,ileg)=p(imu,ileg)
         enddo
      enddo
c     to avoid bugs in HELAS, restore exact masslessness of  incoming partons 
      kr_mad(0,1)=dabs(kr_mad(3,1))
      kr_mad(0,2)=dabs(kr_mad(3,2))


      call checkmomzero(nlegreal,kr_mad)


c$$$c     boost momenta in cm, it may improve madgraph
c$$$      beta=-(kr_mad(3,1)+kr_mad(3,2))/(kr_mad(0,1)+kr_mad(0,2))
c$$$      vec(1)=0
c$$$      vec(2)=0
c$$$      vec(3)=1
c$$$      call mboost(nlegreal,vec,beta,kr_mad,kr_mad_cm)
c$$$c     to avoid bugs in HELAS, restore exact masslessness of  incoming partons 
c$$$      kr_mad_cm(0,1)=dabs(kr_mad_cm(3,1))
c$$$      kr_mad_cm(0,2)=dabs(kr_mad_cm(3,2))
c$$$c      call checkmomzero(nlegreal,kr_mad_cm)
c$$$c      do mu=0,3
c$$$c         write(*,*) kr_mad_cm(mu,1) , kr_mad_cm(mu,2)
c$$$c      enddo

      if(ttype.eq.1) then
         call compreal(kr_mad,rflav,amp2mad)
      elseif(ttype.eq.-1) then
         call compreal_tb(kr_mad,rflav,amp2mad)
      else
         write(*,*) 'wrong ttype in real.f'
         call exit(-1)
      endif
      amp2=amp2mad

      if(checkreal) then
c     Notice that some checks don't work because of
c     how MCFM includes the ckm matrix (cfr Vsum and Vsq).
c     For many cases where there is a disagreement, I was able 
c     to go into the (j,k) loop of qg_tbq_g and check
c     exactly the values before the multiplication
c     for vsum is performed.
cccccccccccccccccccccccccccccccccccccccccccccccccccc
         write(*,*) 'Cannot do these checks.'
         write(*,*) 'Code needed kept only in branches'
         call exit(-1)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
         call mom_to_MCFM_real(p,pmcfm)

cccccccccccccccccccc
c     COPY&PASTE FROM virtual.f
         if (ini) then
            call virtual_initialize_MCFM
            ini=.false.
         endif
ccccccccccccccccccccccccccccccccccccc
c     !ER:
c     relevant for st 4 flavour
         scale = sqrt(st_muren2)
         renscale_H=sqrt(st_muren2)
         renscale_L=sqrt(st_muren2)
         facscale_H=sqrt(st_muren2)
         facscale_L=sqrt(st_muren2)
         
         as_H = st_alpha
         as_L = st_alpha
         mt=topmass_pow
         mb=bmass_pow
cccccccccccccccccccccccccccccccccccccc
c$$$c         print*, '-------------------'
c$$$         call qg_tbq_g(pmcfm,realmcfm)
c$$$
c$$$c         ckm_tb=CKM(abs(rflav(3)),abs(rflav(4)))**2
c$$$
c$$$         
c$$$         if(rflav(1).eq.0.and.rflav(2).eq.0) then
c$$$c     GG
c$$$            if(rflav(5).gt.0) then
c$$$c     normal ordering
c$$$               print*, 'GG ',rflav(5),rflav(6),realmcfm(rflav(1),rflav(2))/amp2
c$$$            elseif(rflav(5).lt.0) then
c$$$c     exchange momenta 5 and 6
c$$$               rtmp=pmcfm(5,4)
c$$$               pmcfm(5,4)=pmcfm(6,4)
c$$$               pmcfm(6,4)=rtmp
c$$$               do imu=1,3
c$$$                  rtmp=pmcfm(5,imu)
c$$$                  pmcfm(5,imu)=pmcfm(6,imu)
c$$$                  pmcfm(6,imu)=rtmp
c$$$               enddo
c$$$               call qg_tbq_g(pmcfm,realmcfm)
c$$$               print*, 'GG ',rflav(5),rflav(6),realmcfm(rflav(1),rflav(2))/amp2
c$$$            endif
c$$$         elseif(
c$$$     $           (rflav(1).eq.0.and.rflav(2).ne.0).or.
c$$$     $           (rflav(2).eq.0.and.rflav(1).ne.0)) then
c$$$c     QG / QG
c$$$            print*, 'QG ',rflav(1),rflav(2),realmcfm(rflav(1),rflav(2))/amp2
c$$$         elseif(
c$$$     $           rflav(1).lt.0.and.rflav(2).gt.0) then
c$$$c     QBARQ
c$$$            print*, 'QBARQ ',rflav(1),rflav(2),rflav(5),rflav(6),realmcfm(rflav(1),rflav(2))/amp2
c$$$         elseif(
c$$$     $           rflav(1).gt.0.and.rflav(2).lt.0) then
c$$$c     QQBAR
c$$$            print*, 'QQBAR ',rflav(1),rflav(2),rflav(5),rflav(6),realmcfm(rflav(1),rflav(2))/amp2
c$$$
c$$$
c$$$         elseif(
c$$$     $           rflav(1).gt.0.and.rflav(2).gt.0) then
c$$$c     QQ
c$$$            print*, 'QQ ',rflav(1),rflav(2),rflav(5),rflav(6),realmcfm(rflav(1),rflav(2))/amp2
c$$$         elseif(
c$$$     $           rflav(1).lt.0.and.rflav(2).lt.0) then
c$$$c     QBARQBAR
c$$$            print*, 'QBARQBAR ',rflav(1),rflav(2),rflav(5),rflav(6),realmcfm(rflav(1),rflav(2))/amp2
c$$$         else
c$$$            print*, 'unrecognized flavour string in setreal'
c$$$            stop
c$$$         endif
      endif

      


c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2*pi))


      end


      subroutine mom_to_MCFM_real(cmpborn,pmcfm)
      implicit none
      include 'nlegborn.h'
      include 'MCFM_include/constants.f'
      double precision cmpborn(0:3,6)
      double precision pmcfm(mxpart,1:4)
      integer ileg,mu
      do ileg=1,2
         pmcfm(ileg,4)=-cmpborn(0,ileg)
         do mu=1,3
            pmcfm(ileg,mu)=-cmpborn(mu,ileg)
         enddo
      enddo
      do ileg=3,6
         pmcfm(ileg,4)=cmpborn(0,ileg)
         do mu=1,3
            pmcfm(ileg,mu)=cmpborn(mu,ileg)
         enddo    
      enddo
c$$$      do ileg=1,6
c$$$         write(*,*) (pmcfm(ileg,mu), mu=1,4)
c$$$      enddo
      end



ccccccccccccccccccccccccccccccccccccccccccccccc
c     From here on, to check with MCFM reals
ccccccccccccccccccccccccccccccccccccccccccccccc
