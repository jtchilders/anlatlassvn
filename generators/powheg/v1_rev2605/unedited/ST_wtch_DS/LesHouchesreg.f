      subroutine gen_leshouches_reg
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'LesHouches.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      integer ileg,fl,tmp
      integer tcol,col2
      data tcol/501/
      data col2/502/
      integer igluon
      data igluon/21/
c id of the event
      idprup=lprup(1)
      xwgtup=+1
c     aqedup must can be set here
c     since the process-dependent PhysPars.h 
c     is well defined
      aqedup=alphaem_pow
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      aqcdup=st_alpha
      nup=nlegreal

c     assign istup, mothup and idup
      do ileg=1,nup
         fl=flst_regular(ileg,rad_realreg)
         if(fl.eq.0) fl=igluon
         idup(ileg)=fl
         if(ileg.gt.2) then
            istup(ileg)=1
            mothup(1,ileg)=1
            mothup(2,ileg)=2
         else
            istup(ileg)=-1
            mothup(1,ileg)=0
            mothup(2,ileg)=0
         endif
         spinup(ileg)=9
         vtimup(ileg)=0
      enddo

c     assign momenta (3rd is W, 4th is top)
      call momenta_lh(kn_preal,nup)

c     cross check
      if((abs(flst_regular(3,rad_realreg)).ne.24).or.
     $     (abs(flst_regular(4,rad_realreg)).ne.6)) then
         write(*,*) 'Error in gen_leshouches_reg'
         call exit(1)
      endif

c     Set color connections for all particles.
c     The only flavour configurations that can be
c     called here is u ubar -> w t bbar.
c     In fact, this is the only fully regular subprocess.
c     It is also a pure double-resonant process.
c     Therefore, we deal only with this case. All the other
c     double-resonant subprocesses have also at least one
c     divergent graph, i.e. they have a valid underlying-Born
c     configuration. Therefore, their color structure is
c     builded automatically from the Born color flow, which
c     for Wt-channel is unambiguos.
      if((idup(1).eq.igluon).or.
     $(idup(2).eq.igluon).or.
     $(idup(5).eq.igluon).or.
     $(abs(idup(1)).eq.abs(idup(5))).or.
     $(abs(idup(2)).eq.abs(idup(5)))) then
         write(*,*) 'Error 1 in gen_leshouches_reg ',
     $        idup(1),idup(2),idup(5)
         call exit(1)
      endif

      do ileg=1,nup
         if(abs(idup(ileg)).eq.24) then
            icolup(1,ileg)=0
            icolup(2,ileg)=0
         elseif(abs(idup(ileg)).eq.6) then
            icolup(1,ileg)=tcol
            icolup(2,ileg)=0
         elseif(ttype*idup(ileg).gt.0) then
            icolup(1,ileg)=tcol
            icolup(2,ileg)=0
         elseif(ttype*idup(ileg).lt.0) then
            icolup(1,ileg)=0
            icolup(2,ileg)=col2
         else
            write(*,*) 'Error 2 in gen_leshouches_reg ',
     $           idup(1),idup(2),idup(5)
            call exit(1)
         endif
      enddo

c     if charge-conjugated process, exchange color
c     with anticolors
      if(ttype.lt.0) then
         do ileg=1,nup
            tmp=icolup(1,ileg)
            icolup(1,ileg)=icolup(2,ileg)
            icolup(2,ileg)=tmp
         enddo
      endif

c     add resonance 
      call finalize_lh
c     Don't forget to set scale for scalup equal to the pt of the 
c     radiation (whatever it is now!)
      scalup=sqrt(rad_pt2max)
      end
