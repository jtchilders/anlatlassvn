c-------------------------------------------------------------------
      subroutine usrcut(lnot,wusr)
c-------------------------------------------------------------------
c PRIMARY CUTS ALREADY APPLIED TO PHASE-SPACE GENERATION:
c     ptjmin < pt(jet) < ptjmax for all light jets
c     -etajmax < eta(jet) < etajmax for all light jets
c     delta R(jj) > drjmin for all (light jet, light jet) pairs
c     pt(lept)>ptlmin  etmiss > minetmiss 
c     abs(eta(lept)) < etalmax 
c     lepton/jet isolation 
c     pt(ph) > ptphmin for all photons
c     -etaphmax < eta(ph) < etaphmax for all light photons
c     delta R(j-ph) > drphjmin for all (light jet, photon) pairs
c     delta R(l-ph) > drphlmin for all (lepton, photon) pairs
c
c USE THIS ROUTINE TO ENFORCE OTHER CUTS    
      implicit none
      include 'alpgen.inc'
      include 'wphjet.inc'
      integer lnot
      double precision wusr
c
      lnot=0
      wusr=1d0
c
c  USR will add possible extra cuts at this point. 
c     if(cut-not-passed) goto 10
      return
 10   lnot= 1
      wusr= 0.d0 
      end

c-------------------------------------------------------------------
      subroutine alshis
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wphjet.inc'
      real ptbin,ptmax,xmbin,xmmax
      character*1 ijet(6)
      integer i
      data ijet/'1','2','3','4','5','6'/
      ptbin=2.5e0
      ptmax=100*ptbin
      xmbin=4e0
      xmmax=400e0
c
      do i=1,min(5,njets)
        call mbook(i,'pt j'//ijet(i),ptbin,0e0,ptmax)
        call mbook(5+i,'eta j'//ijet(i),0.1,-3e0,3e0)
      enddo
      call mbook(12,'ptlept',2.,0e0,200.)
      call mbook(13,'mW',0.5,70.,110.)
      call mbook(14,'etal',0.2,-5.,5.)
      end
c-------------------------------------------------------------------
      subroutine alfhis
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wphjet.inc'
      integer i
      real  xnorm
      character *1 jet(9)
      data jet/'1','2','3','4','5','6','7','8','9'/
c debug
      integer idbg
      double precision fcount
      common/fldbg/fcount(16),idbg
      data idbg/0/
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
      do i=1,min(5,njets)
        call mtop(i,99,'pt j'//jet(i),' ','LOG')
      enddo
      do i=1,min(5,njets)
        call mtop(5+i,99,'eta j'//jet(i),' ','LIN')
      enddo
c
      call mtop(12,99,'ptl',' ','LIN')
      call mtop(13,99,'mW',' ','LIN')
      call mtop(14,99,'etal',' ','LIN')
c
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
      include 'wphjet.inc'
      integer n
      character *50 mon_fname
c
      if(evgen) then
         if(mod(n,1000000).eq.0) then
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
      include 'wphjet.inc'
      real*8 wgt,xmw
      real rwgt
      integer i,jproc,ord(10)
c
      rwgt=real(wgt)
      if(rwgt.lt.0e0) then
         write(*,*) 'negative wgt=',wgt
         return
      elseif (rwgt.eq.0e0) then
         return
      endif 
c
      call mfill(12,real(ptlep),rwgt)
      if(njets.eq.0) return
      call alusor(ptj,njets,ord,2)
      do i=1,min(5,njets)
        call mfill(i,real(ptj(ord(njets-i+1))),rwgt)
        call mfill(5+i,real(etaj(ord(njets-i+1))),rwgt)
      enddo
      xmw=sqrt(pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2)
      call mfill(13,real(xmw),rwgt)
      call mfill(14,real(etalep),rwgt)
      end
      
