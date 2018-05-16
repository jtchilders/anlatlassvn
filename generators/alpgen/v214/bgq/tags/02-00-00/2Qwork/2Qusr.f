c-------------------------------------------------------------------
      subroutine alshis
c-------------------------------------------------------------------
      include 'alpgen.inc'
      include '2Q.inc'
      real*4 ctbin,ctmin,ctmax
      ptbin=2e0
      ptmax=200e0
      xmbin=4e0
      xmmax=400e0
      call mbook(1,'pt_1',ptbin,0e0,ptmax)
      call mbook(2,'pt_2',ptbin,0e0,ptmax)
      call mbook(3,'pt_3',ptbin,0e0,ptmax)
      call mbook(4,'pt_4',ptbin,0e0,ptmax)
      call mbook(5,'pt_5',ptbin,0e0,ptmax)
      call mbook(6,'pt_6',ptbin,0e0,ptmax)
      call mbook(8,'pt_Q',2.*ptbin,0e0,2.*ptmax)
      call mbook(9,'dR(Q-Qbar)',0.2,0e0,6e0)
      if(ihvy.eq.5) call mbook(10,'mQQ ',xmbin,0e0,xmmax)
      if(ihvy.eq.6) call mbook(10,'mQQ ',3.*xmbin,0e0,3.*xmmax)
      call mbook(13,'ptjet(incl)',ptbin,0e0,ptmax)
      end

      subroutine usrcut(lnot,weight)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     Applies additional kinematical cuts to the final state 
c     during the phase-space generation                          c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'alpgen.inc'
      include '2Q.inc'
      integer lnot
      double precision weight

      weight= 1.d0
      lnot= 0
c
c     if(cut-not-passed) goto 10
c
 20   return
 10   lnot= 1
      return
      end

c-------------------------------------------------------------------
      subroutine alfhis
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include '2Q.inc'
c debug
      double precision rate(16)
      common/dbg/rate
c
      integer i
      real  xnorm
      character *1 jet(9)
      data jet/'1','2','3','4','5','6','7','8','9'/
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
c
      do i=1,200
         call mopera(i,'F',i,i,xnorm,1.)
         call mfinal(i)
      enddo 
      do i=1,njets
         call mtop(i,99,'pt'//jet(i),' ','LOG')
      enddo
      call mtop(8,99,'pt(Q,Qbar)',' ','LOG')
      call mtop(9,99,'dR(Q-Qbar)',' ','LIN')
      call mtop(10,99,'m(Q-Qbar)',' ','LOG')
      call mtop(13,99,'pt(jet)',' ','LOG')
c
      close(99)
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
      include '2Q.inc'
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
      include '2Q.inc'
      real*8 mQQ,wgt
      real rwgt
      integer i,jproc,ord(10)
c
      integer nmax        !maximum number of external particles, change here
      parameter (nmax=10)          

      real*8 mbb,mnn,mmt,gi(4),ptrb,ptrbbar,ptt,ptm,ptnt,ptnm
      real*8 ctnn,ctbb,ctmt,ctbts,ctnts,ctmts
      real*8 qboost(4),pbstar(4),qtop(4),pmustar(4),pnustar(4)
      data gi/1.d0,-1.d0,-1.d0,-1.d0/

*
      rwgt=real(wgt)

      if(rwgt.lt.0e0) then
         write(*,*) 'negative wgt=',wgt
         return
      elseif (rwgt.eq.0e0) then
         return
      endif
c     reordering according to pt
      if(ihvy.eq.5) call alusor(ptj,njets+2,ord,2)              
      if(njets.gt.0.and.ihvy.eq.6) call alusor(ptj,njets,ord,2)              
      do i=1,njets
         call mfill(i,real(ptj(ord(njets+1-i))),rwgt)
         call mfill(13,real(ptj(ord(njets+1-i))),rwgt)
      enddo
      call mfill(8,real(pt(3)),rwgt)
      call mfill(8,real(pt(4)),rwgt)
c     dR(Q,Qbar)
      call mfill(9,real(dr(3,4)),rwgt)
c     Q-Qbar mass
      if(ihvy.eq.5) then
         mQQ=(pbott(4)+pbbar(4))**2
         do i=1,3
            mQQ=mQQ-(pbott(i)+pbbar(i))**2
         enddo
         if(mQQ.le.0)then
            write(*,*) 'mQQ^2=',mQQ,' <0'
            stop
         endif 
      elseif(ihvy.eq.6) then
         mQQ=(ptop(4)+ptbar(4))**2
         do i=1,3
            mQQ=mQQ-(ptop(i)+ptbar(i))**2
         enddo
         if(mQQ.le.0)then
            write(*,*) 'mQQ^2=',mQQ,' <0'
            stop
         endif 
      endif
      mQQ=sqrt(mQQ)
      call mfill(10,real(mQQ),rwgt)

      end


