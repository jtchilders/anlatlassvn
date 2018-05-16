c-------------------------------------------------------------------
      subroutine alshis
c-------------------------------------------------------------------
      include 'alpgen.inc'
      include 'QQh.inc'
      ptbin=4e0
      ptmax=400e0
      xmbin=4e0
      xmmax=400e0
      call mbook(1,'pt_1',ptbin,0e0,ptmax)
      call mbook(2,'pt_2',ptbin,0e0,ptmax)
      call mbook(3,'pt_3',ptbin,0e0,ptmax)
      call mbook(4,'pt_4',ptbin,0e0,ptmax)
      if(ihvy.eq.5) then 
         call mbook(8,'pt_b_min',2.*ptbin,0e0,2.*ptmax)
         call mbook(9,'dR(b-bbar)',0.2,0e0,6e0)
         call mbook(10,'mbb',xmbin,0e0,xmmax)
      endif
      if(ihvy.eq.6) then 
         call mbook(10,'mtt',3.*xmbin,0e0,3.*xmmax)
      endif
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
      include 'QQh.inc'
      integer lnot
      double precision weight

      weight= 1.d0
      lnot= 0
c      if(cuts not passed) goto 10

 20   return
 10   lnot= 1
      return
      end

c-------------------------------------------------------------------
      subroutine alfhis
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'QQh.inc'
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
c debug
c      if(imode.eq.2) then
c        write(*,*) (rate(i)*xnorm,i=1,6)
c        write(*,*) (rate(i)*xnorm,i=11,16)
c      endif

      do i=1,200
         call mopera(i,'F',i,i,xnorm,1.)
         call mfinal(i)
      enddo 
      do i=1,njets
         call mtop(i,99,'pt'//jet(i),' ','LOG')
      enddo

      if(ihvy.eq.5) then
         call mtop(8,99,'pt_b_min',' ','LOG')
         call mtop(9,99,'dR(b-bbar)',' ','LIN')
         call mtop(10,99,'mbb',' ','LOG')
      endif
      if(ihvy.eq.6) then
         call mtop(10,99,'mtt',' ','LOG')
      endif
      call mtop(13,99,'pt(jet)',' ','LOG')
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
      include 'QQh.inc'
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
      include 'QQh.inc'
      double precision mQQ,wgt
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
c     reordering according to pt
      if(njets.gt.0) then
        call alusor(ptj,njets,ord,2)
        do i=1,njets
          call mfill(i,real(ptj(ord(njets+1-i))),rwgt)
          call mfill(13,real(ptj(ord(njets+1-i))),rwgt)
        enddo
      endif
      call mfill(8,real(pt(4)),rwgt)
      call mfill(8,real(pt(5)),rwgt)
c     dR(Q,Qbar)
      call mfill(9,real(dr(4,5)),rwgt)
c     Q-Qbar mass
      if(ihvy.eq.5) then
         mQQ=(pbott(4)+pbbar(4))**2
         do i=1,3
            mQQ=mQQ-(pbott(i)+pbbar(i))**2
         enddo
      elseif(ihvy.eq.6) then
         mQQ=(ptop(4)+ptbar(4))**2
         do i=1,3
            mQQ=mQQ-(ptop(i)+ptbar(i))**2
         enddo
      endif
      mQQ=sqrt(mQQ)
      call mfill(10,real(mQQ),rwgt)
      end
