      subroutine setreal(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'

      real * 8 p(0:3,nlegreal)
      integer rflav(nlegreal)
      real * 8 amp2,amp2_mad

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


cccccccccccccccccc
      integer nleg
      parameter (nleg=nlegreal)
      real *8 kr(0:3,nleg)
      integer rflav_ME(nleg)
      real *8 ewcoupl
      integer mu,ileg

      integer ftemp,mflav(nleg)
      real *8 ktemp,kr_mad(0:3,nleg)

      real *8 dotp
      external dotp

      logical ini
      data ini/.true./
      save ini
      
ccccccccccccccccccccccccccccccccccccccc
c     charge conjugation
c     if ttype=-1, then rflav has been filled with tbar
c     production flavours. Subroutines here work for t flavour.
c     Therefore, invert the sign of local flavours.
      do ileg=1,nleg
         rflav_ME(ileg)= ttype *rflav(ileg)
      enddo
ccccccccccccccccccccccccccccccccccccccc

c     local copy of input variables (p->kr)
      do ileg=1,nleg
         do mu=0,3
            kr(mu,ileg)=p(mu,ileg)
         enddo
      enddo

c     check
      if ((abs(rflav(3)).ne.24).or.(abs(rflav(4)).ne.6)) then
         write(*,*) 'real_ampsq: ERROR in flavor assignement'
         call exit(1)
      endif

c     ew coupling
      ewcoupl=4d0*pi*alphaem_pow/sthw2_pow

ccccccccccccccccccccccccccccccccccccccccccc
c     >>> WT CHANNEL <<<
ccccccccccccccccccccccccccccccccccccccccccc

c     USING MADGRAPH SUBROUTINES
      do ileg=1,5
         mflav(ileg)=rflav_ME(ileg)
         do mu=0,3
            kr_mad(mu,ileg)=kr(mu,ileg)
         enddo
      enddo
c     to avoid bugs in HELAS, restore exact masslessness of incoming partons 
      kr_mad(0,1)=dabs(kr_mad(3,1))
      kr_mad(0,2)=dabs(kr_mad(3,2))
c     reassign here helas couplings and parameters that 
c     can change on an event-by-event basis
      alfas=st_alpha
      mtMS=sqrt(dotp(kr_mad(0,4),kr_mad(0,4)))
      tmass=mtMS
      twidth=0d0
      wwidth=0d0
      call my_setpara
c     invert 3rd and 4th particles before passing the array to
c     madgraph.
      ftemp=mflav(4)
      mflav(4)=mflav(3)
      mflav(3)=ftemp
      do mu=0,3
         ktemp=kr_mad(mu,4)
         kr_mad(mu,4)=kr_mad(mu,3)
         kr_mad(mu,3)=ktemp
      enddo


c     protect from extreme soft region: if one of the 2 incoming momenta has
c     an energy much larger than the FKS parton, a problem can occur in
c     madgraph subroutines (division by zero due to roundings)
c$$$      if((kr_mad(0,5)/kr_mad(0,1).lt.1d-9).or.
c$$$     $(kr_mad(0,5)/kr_mad(0,2).lt.1d-9)) then
c$$$         write(*,*) 'Critical kinematical point. Real contrib dropped'
c$$$         amp2_mad=1d-20
c$$$      else
         call choose_real_process(kr_mad,mflav,amp2_mad)
c$$$      endif

ccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccc
c     assign output
      amp2=amp2_mad/(st_alpha/2./pi)
cccccccccccccccccccccc
      end
