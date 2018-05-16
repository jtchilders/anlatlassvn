c-------------------------------------------------------------------
      subroutine usrcut(lnot,wusr)
c-------------------------------------------------------------------
c PRIMARY CUTS ALREADY APPLIED TO PHASE-SPACE GENERATION:
c     ptjmin < pt(jet) < ptjmax for all light jets
c     -etajmax < eta(jet) < etajmax for all light jets
c     delta R(jj) > drjmin for all (light jet, light jet) pairs
c     delta R(bj) > drjmin for all (heavy quark, light jet) pairs
c     ptbmin < pt(b) < ptbmax for all heavy quarks
c     -etabmax < eta(b) < etabmax for all heavy quarks 
c     delta R(bb) > drjmin for all (heavy quark, heavy quark) pairs
c     lepton/jet isolation 
c     pt(ph) > ptphmin for all photons
c     -etaphmax < eta(ph) < etaphmax for all light photons
c     delta R(j-ph) > drphjmin for all (light jet, photon) pairs
c     delta R(l-ph) > drphlmin for all (lepton, photon) pairs
c
c USE THIS ROUTINE TO ENFORCE OTHER CUTS    
      implicit none
      include 'alpgen.inc'
      include 'wphqq.inc'
      integer lnot,i
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
      ptbin=2e0
      ptmax=200e0
      xmbin=4e0
      xmmax=400e0
      call mbook(1,'pt_1',2.*ptbin,0e0,2.*ptmax)
      call mbook(2,'pt_2',ptbin,0e0,ptmax)
      call mbook(3,'pt_3',ptbin,0e0,ptmax)
      call mbook(4,'pt_4',ptbin,0e0,ptmax)
      call mbook(5,'pt_5',ptbin,0e0,ptmax)
      call mbook(6,'pt_6',ptbin,0e0,ptmax)
      call mbook(7,'ptel',ptbin,0e0,ptmax)
      call mbook(8,'pt_b',ptbin,0e0,ptmax)
      call mbook(9,'dR(b-bbar)',0.2,0e0,6e0)
      call mbook(10,'mbb ',xmbin,0e0,xmmax)
      call mbook(11,'dR(b-jetmax)',0.2,0e0,6e0)
      call mbook(12,'dR(b-jetmin)',0.2,0e0,6e0)
      end

c-------------------------------------------------------------------
      subroutine alfhis
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wphqq.inc'
      integer i
      real  xnorm
      character *1 jet(9)
      data jet/'1','2','3','4','5','6','7','8','9'/
c debug
      integer idbg
      double precision fcount
      common/fldbg/fcount(16),idbg
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
c dbg:
c      if(idbg.eq.1) then
c        write(*,*) 'quarks: ',(fcount(i)*xnorm,i=1,4)
c        write(*,*) 'leptons: ',(fcount(i)*xnorm,i=11,16)
c        call mfinal(100)
c        call mtop(100,99,'rn',' ','LIN')
c        stop
c      endif

      do i=1,200
         if(i.ne.61) call mopera(i,'F',i,i,xnorm,1.)
         call mfinal(i)
      enddo 
      do i=1,njets+2
         call mtop(i,99,'pt'//jet(i),' ','LOG')
      enddo
      call mtop(7,99,'pt(lept)',' ','LOG')
      call mtop(8,99,'pt(b,bbar)',' ','LOG')
      call mtop(9,99,'dR(b-bbar)',' ','LIN')
      call mtop(10,99,'m(b-bbar)',' ','LOG')
      call mtop(11,99,'dR(b-jetmax)',' ','LIN')
      call mtop(12,99,'dR(b-jetmin)',' ','LIN')
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
      include 'wphqq.inc'
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
      include 'wphqq.inc'
      real*8 mbb,wgt,tmp
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
c     reordering according to pt
      call alusor(ptj,njets+2,ord,2)              
      do i=1,njets+2
         call mfill(i,real(ptj(ord(njets+3-i))),rwgt)
      enddo
      call mfill(7,real(ptlep),rwgt)
      call mfill(8,real(ptb),rwgt)
      call mfill(8,real(ptbb),rwgt)

      call mfill(9,real(drbb),rwgt)
c     dR(b-jetmax)
      tmp=drbj(ord(njets+2))
      if(tmp.eq.0) tmp=drbj(ord(njets+1))
      call mfill(11,real(tmp),rwgt)
c     dR(bbar-jetmax)
      tmp=drbbj(ord(njets+2))
      if(tmp.eq.0) tmp=drbbj(ord(njets+1))
      call mfill(11,real(tmp),rwgt)
c     dR(b-jetmin)
      tmp=drbj(ord(1))
      if(tmp.eq.0) tmp=drbj(ord(2))
      call mfill(12,real(tmp),rwgt)
c     dR(bbar-jetmin)
      tmp=drbbj(ord(1))
      if(tmp.eq.0) tmp=drbbj(ord(2))
      call mfill(12,real(tmp),rwgt)
c     b-bbar mass
      mbb=(pbott(4)+pbbar(4))**2
      do i=1,3
         mbb=mbb-(pbott(i)+pbbar(i))**2
      enddo
      mbb=sqrt(mbb)
      call mfill(10,real(mbb),rwgt)
      end


