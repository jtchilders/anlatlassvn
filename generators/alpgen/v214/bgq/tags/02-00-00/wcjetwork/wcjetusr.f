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
c
c USE THIS ROUTINE TO ENFORCE OTHER CUTS    
      implicit none
      include 'alpgen.inc'
      include 'wcjet.inc'
      integer lnot,i
      double precision wusr,ptmax
c
      lnot=0
      wusr=1d0
      ptmax=-100
      do i=1,njets+1
         ptmax=max(ptmax,ptj(i))
c         write(6,*) " pti: " , i ,njets, ptj(i)
      enddo
c      write(6,*) " ptmaxi: " , ptmax ,ptj1min,njets
      if(ptj1min.gt.0.and.ptmax.lt.ptj1min) goto 10
      if(ptj1max.gt.0.and.ptmax.gt.ptj1max) goto 10
c      write(6,*)"accepted "
c
c$$$      do i=1,nfspart
c$$$        if(abs(ifl(i+2)).eq.4) then
c$$$          if(ptj(i).lt.ptcmin) goto 10
c$$$          if(abs(etaj(i)).gt.etacmax) goto 10
c$$$        else
c$$$          if(ptj(i).lt.20d0) goto 10
c$$$          if(abs(etaj(i)).gt.2.5d0) goto 10
c$$$        endif
c$$$      enddo
c
c  USR will add possible extra cuts at this point. 
      return
 10   lnot= 1
      end

c-------------------------------------------------------------------
      subroutine alshis
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wcjet.inc'
      real ptbin,ptmax,xmbin,xmmax
      ptbin=2e0
      ptmax=200e0
      xmbin=4e0
      xmmax=400e0
c
      call mbook(1,'ptj1',ptbin,0e0,ptmax)
      call mbook(2,'ptj2',ptbin,0e0,ptmax)
      call mbook(3,'ptj3',ptbin,0e0,ptmax)
      call mbook(4,'ptj4',ptbin,0e0,ptmax)
      call mbook(5,'ptj5',ptbin,0e0,ptmax)
      call mbook(7,'ptlept',2.,0e0,200.)
      call mbook(11,'etaj1',0.5,-5.,5.)
      call mbook(12,'etaj2',0.5,-5.,5.)
      call mbook(13,'etaj3',0.5,-5.,5.)
      call mbook(14,'etaj4',0.5,-5.,5.)
      call mbook(15,'etaj5',0.5,-5.,5.)
      call mbook(17,'etal',0.5,-5.,5.)
      call mbook(20,'ht',5e0,0e0,500e0)
c
      end
c-------------------------------------------------------------------
      subroutine alfhis
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wcjet.inc'
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
      call multitop(1,99,1,1,'ptj1',' ','LOG')
      do i=2,njets
         call mtfill(i,99,1,1,'ptj'//jet(i),' ','LOG')
      enddo
      call newplot
      call mtop(7,99,'ptlep',' ','LOG')
c
      call multitop(11,99,1,1,'etaj1',' ','LIN')
      do i=2,njets
         call mtfill(10+i,99,1,1,'etaj'//jet(i),' ','LIN')
      enddo
      call newplot
      call mtop(17,99,'etalep',' ','LIN')
c
      call mtop(20,99,'Ht',' ','LIN')

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
      include 'wcjet.inc'
      integer n
      character *15 mon_fname
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
      include 'wcjet.inc'
      real*8 wgt,xmw,ht
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
      call mfill(7,real(ptlep),rwgt)
      call mfill(17,real(etalep),rwgt)
      if(njets.eq.0) return
c
      call alusor(ptj,njets,ord,2)
      ht=ptlep+ptmiss
      do i=1,njets
        ht=ht+ptj(i)
        call mfill(i,real(ptj(ord(njets+1-i))),rwgt)
        call mfill(10+i,real(etaj(ord(njets+1-i))),rwgt)
      enddo
      call mfill(20,real(ht),rwgt)
      end
      
