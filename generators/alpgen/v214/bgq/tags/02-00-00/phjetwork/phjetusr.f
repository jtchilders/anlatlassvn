c      data resc/1d3/
c-----------------------------------------------------------------
      subroutine usrcut(lnot,wusr)
c-----------------------------------------------------------------
c PRIMARY CUTS ALREADY APPLIED TO PHASE-SPACE GENERATION:
c     ptjmin < pt(jet) < ptjmax for all light jets
c     -etajmax < eta(jet) < etajmax for all light jets
c     delta R(jj) > drjmin for all (light jet, light jet) pairs
c USE THIS ROUTINE TO ENFORCE OTHER CUTS    
c
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
      double precision wusr
      integer lnot,i
c initialize output parameters
      lnot=0
      wusr=1d0
c
c apply cut (negative threshold: cut off) 
      if (ptphmax.gt.0) then
         do i= njets+3,njets+nph+2
c           write(6,*) "DEBUG",pt(i),ptphmax
            if (pt(i).gt.ptphmax)           goto 10
         enddo
      endif
c cut passed
      return
c if(cut-not-passed) goto 10
 10   lnot= 1
      return
      end

c-------------------------------------------------------------------
      subroutine alshis
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
      real ptbin,ptmax,xmbin,xmmax
      ptbin=2e0
      ptmax=200e0
      xmbin=4e0
      xmmax=400e0
c
      call mbook(1,'pt_1',2.*ptbin,0e0,2.*ptmax)
      call mbook(2,'pt_2',ptbin,0e0,ptmax)
      call mbook(3,'pt_3',ptbin,0e0,ptmax)
      call mbook(4,'pt_4',ptbin,0e0,ptmax)
      call mbook(5,'pt_5',ptbin,0e0,ptmax)
      call mbook(6,'pt_6',ptbin,0e0,ptmax)
      end

c-------------------------------------------------------------------
      subroutine alfhis
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
      integer i
      real  xnorm
      character *1 jet(9)
      data jet/'1','2','3','4','5','6','7','8','9'/
c
      open(unit=99,file=topfile,err=101,status='unknown')
      if(imode.le.1) then
         xnorm=sngl(avgwgt/totwgt)
      elseif(imode.eq.2) then
         xnorm=1e0/real(unwev)
      else
         write(*,*) 'imode type not allowed, stop'
         stop
      endif
c
      do i=1,200
         if(i.ne.61) call mopera(i,'F',i,i,xnorm,1.)
         call mfinal(i)
      enddo 
c
      do i=1,njets
         call mtop(i,99,'pt'//jet(i),' ','LOG')
      enddo

 100  close(99)
 101  return
      end

      subroutine monitor(n,mon_fname)
c     This routine is called by default every 100K events.
c     The user can use it to get regular updates on the run
c     while this is progressing. Textual output can be written to file
c     fname, where partial cross-sections and and generation
c     efficiencies have already been printed by default
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
      integer n
      character *50 mon_fname
c
      if(evgen) then
         if(mod(n,100000).eq.0) then
c     save histograms' contents
            call msave
c     print out histograms
            call alfhis
c     restore original contents, to proceed with analysis
            call mrestore
         endif 
      endif 
      end
c-------------------------------------------------------------------
      subroutine aleana(jproc,wgt)
c     analyse event, fill histograms
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
      real*8 wgt,tmp,ptlep,etmiss,xmz,mll
      real rwgt
      integer i,j,jproc,ord(10)
c
      rwgt=real(wgt)
      if(rwgt.lt.0e0) then
         write(*,*) 'negative wgt=',wgt
         return
      elseif (rwgt.eq.0e0) then
         return
      endif 
c
c     reordering according to pt
      call alusor(ptj,njets,ord,2)              
      do i=1,njets
         call mfill(i,real(ptj(ord(njets+1-i))),rwgt)
      enddo
c
      end
















