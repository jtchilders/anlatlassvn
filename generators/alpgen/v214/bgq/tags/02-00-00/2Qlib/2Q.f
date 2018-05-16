c-------------------------------------------------------------------
      subroutine alsprc
c     assigns the hard process code
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      ihrd=6
      end

c-------------------------------------------------------------------
      subroutine alhset
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include '2Q.inc'
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c     ngrid is the total number of grids allowd for P.S. variables.
c     jgrid(jproc) labels the grid associated to a given jproc.
c       jgrid= 1   -->   g g
c       jgrid= 2   -->   q qbar and qbar q
c       jgrid= 3   -->   g q    and g qbar
c       jgrid= 4   -->   q g    and qbar g
c       jgrid= 5   -->   q q     and qbar qbar
      data jgrid/1,2, 3,4, 1,2*2, 5, 2,5, 3,4,3,4,1,1, 84*0/
      integer i
      character*2 Q(4:6)
      character*5 Qbar(4:6)
      data Q/' c',' b',' t'/,Qbar/' cbar',' bbar',' tbar'/
c     parameters for the gauge invariance prescription:
      winsize  = 2.d0/pi
      resonance= 'n'
      wmode    = 'yy'
c     process input parameters
      npart=njets+4
      nprtns=njets+2
      if(ihvy.eq.4) then
        idecay='n'
        mq=mc
        ptQmin=ptcmin
        etaQmax=etacmax
        drQmin=drcmin
        paruse(21,6)=0
        paruse(31,6)=0
        paruse(41,6)=0
        paruse(51,6)=0
      elseif(ihvy.eq.5) then
        idecay='n'
        mq=mb
        ptQmin=ptbmin
        etaQmax=etabmax
        drQmin=drbmin
        paruse(22,6)=0
        paruse(32,6)=0
        paruse(42,6)=0
        paruse(52,6)=0
      elseif(ihvy.eq.6) then
        if(itdec.eq.1) then
          idecay='y'
        else
          idecay='n'
        endif
        mq=mt
        ptQmin=0
        etaQmax=100
        drQmin=0
        paruse(21,6)=0
        paruse(31,6)=0
        paruse(41,6)=0
        paruse(51,6)=0
        paruse(22,6)=0
        paruse(32,6)=0
        paruse(42,6)=0
        paruse(52,6)=0
      endif
c masses
      do i=1,njets+4
         p(5,i)=0
      enddo
      p(5,3)=mq
      p(5,4)=mq
c
      if(njets.eq.0) then
         jprocmax=2
         ngrid   = 2
      elseif(njets.eq.1) then
         jprocmax=4
         ngrid   = 4
      elseif(njets.eq.2) then
         jprocmax=10
         ngrid   = 5
      elseif(njets.eq.3) then
         jprocmax=14
         ngrid   = 5
      elseif(njets.eq.4) then
         jprocmax=16
         ngrid   = 5
      elseif(njets.eq.5) then
         jprocmax=16
         ngrid   = 5
      elseif(njets.eq.6) then
         jprocmax=16
         ngrid   = 5
      else 
         print*,'njets=',njets,' not yet available'
         stop
      endif 
      write(*,*) ' '
c run parameters:
c write process code, code=.. for QQbar +jets
      if(njets.gt.0) then
        write(niosta,*) Q(ihvy),Qbar(ihvy)
     $       ,'+',njets,' jets'
      else
        write(niosta,*) Q(ihvy),Qbar(ihvy)
      endif 
      write(niosta,*)
     $     '======================================='
      write(niosta,*) 'Heavy quark mass: m(',Q(ihvy),')=',mq
      write(niosta,*)
     $     'Generation cuts for the partonic event sample:' 
      write(niosta,*) '     Heavy quarks:'
      if(ihvy.eq.5.or.ihvy.eq.4) then
        write(niosta,*) '    Heavy quark jets:'
        write(niosta,*) 'ptmin=',ptQmin,' |etamax|=',etaQmax
     $       ,' dR(bb)>',drQmin 
      elseif(ihvy.eq.6) then
        write(niosta,*) 'No generation cuts on t tbar'
      endif
      if(njets.gt.0) then
        write(niosta,*) '     Light jets:'
        write(niosta,*) 'ptmin=',ptjmin,' |etamax|=',etajmax
     $       ,' dR(j-j),dR(Q-j)>',drjmin 
      write(niosta,*) 'ptj1min=',ptj1min,' ptj1max=',ptj1max
      endif
      end
c-------------------------------------------------------------------
      subroutine selflav(jproc,xlum,afl)
c     evaluates parton luminosities, and assigns PDG-code flavours
c     gluon=21 cbar=-4 sbar=-3 ubar=-2 dbar=-1 d=1 u=2 s=3 c=4
c     jproc
c
c---Q Qbar final state (+ up to 4 jets)
c
c  1  g g    -> Q Qbar (+ gluons)
c  2  q qbar -> Q Qbar (+ gluons)   q=u,d,c,s
c     qbar q -> Q Qbar (+ gluons)   q=u,d,c,s
c
c---Q Qbar + 1 q jet final state
c
c  3  g q    -> Q Qbar q    (+ gluons)
c     g qbar -> Q Qbar qbar (+ gluons)
c  4  q g    -> Q Qbar q    (+ gluons)
c     qbar g -> Q Qbar qbar (+ gluons)
c
c---Q Qbar + 2 q jets final state
c
c  5  g g    -> Q Qbar q qbar   (+ gluons)
c  6  q qbar -> Q Qbar q qbar   (+ gluons)   q=u,d,c,s
c     qbar q -> Q Qbar q qbar   (+ gluons)   q=u,d,c,s
c  7  q qbar -> Q Qbar q' qbar' (+ gluons)   q=u,d,c,s
c     qbar q -> Q Qbar q' qbar' (+ gluons)   q=u,d,c,s
c  8  q q'       -> Q Qbar q q'           (+ gluons)    q=u,d,c,s
c     qbar q'bar -> Q Qbar qbar q'bar     (+ gluons)    q=u,d,c,s
c  9  q qbar'    -> Q Qbar q qbar'        (+ gluons)    q=u,d,c,s
c     qbar q'    -> Q Qbar qbar q'        (+ gluons)    q=u,d,c,s
c 10  q q        -> Q Qbar q q            (+ gluons)    q=u,d,c,s
c     qbar qbar  -> Q Qbar qbar qbar      (+ gluons)    q=u,d,c,s
c
c---Q Qbar + 3 q jets final state
c
c 11  g q    -> Q Qbar q q' qbar'    (+ gluon)
c     g qbar -> Q Qbar qbar q' qbar' (+ gluon)
c 12  q g    -> Q Qbar q q' qbar'    (+ gluon)
c     qbar g -> Q Qbar qbar q' qbar' (+ gluon)
c 13  g q    -> Q Qbar q q qbar      (+ gluon)
c     g qbar -> Q Qbar qbar q qbar   (+ gluon)
c 14  q g    -> Q Qbar q q qbar      (+ gluon)
c     qbar g -> Q Qbar qbar q qbar   (+ gluon)
c
c---Q Qbar + 4 q jets final state
c
c 15  g g    -> Q Qbar q qbar q qbar   (+ gluons)  q=u,d,c,s
c 16  g g    -> Q Qbar q qbar q' qbar' (+ gluons)
c-----------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include '2Q.inc'
c commons
      integer icconj
      common/hwconv/icconj
      integer afl(maxpar)
      real tmp(100),slum,cwgt,swgt,rn,tmptot
c
      double precision xlum,xrn                                 
      integer i,k,itmp,icount,init,jproc,ng,j
      real cfac(0:6),ifact(0:6)
      data cfac/8e0,6*3e0/,init/0/
      integer imap(-4:4)
      data imap/-2,-1,-2,-1,0,1,2,1,2/
      integer ib,ibb,iflaux
c
c   for top decay
c
      complex*16 vtp(4),utpb(4) !array containing the individual field 
      real*8 ptp(4),ptpb(4)
      common/tspinors/vtp,utpb
      real*8 pf(4),pfb(4),pb(4)
      common/pdec/pf,pfb,pb
      integer jt,jtb
      integer exch(4)
      data exch/2,3,4,1/
c
c     overall efficiency for extraction of colour states:
c     >(non-zero color states) / 3**nq*8**ng
c     ccoef(i,j) for njets=j and i=#(light quark pairs)
      double precision effco
      double precision ccoef(3,0:8) !ccoeff(# of q-qb pairs,  # of gluons)
      data ccoef/0.333333333,0.185185185,0.127572016,
     >           0.1666666667,0.12037037,0.0936213992,
     >           0.114583333,0.0902777778,0.0743312757,
     >           0.087239583,0.072337963,0.06171232,
     >           0.070475260,0.0603841146,0.05278461,
     >           0.059122721,0.051834672,-1.d0,
     >           0.050923665,0.045412134,-1.d0,
     >           0.042060375,-1.d0,-1.d0,
     >           0.037214041,-1.d0,-1.d0/
      integer nlq  ! number of light quarks
      save init,cfac,ifact,nlq

c
      if(init.eq.0) then
         njets=npart-4
         init=1
         ifact(0)=1e0
         do i=1,6
            ifact(i)=ifact(i-1)/real(i)
         enddo
         nlq=4
         if(ihvy.eq.4) nlq=3
      endif
c
      do i=1,npart
         ifl(i)=0
      enddo 
 1    call randa(xrn)
      rn=real(xrn)
      if(1e0-rn.lt.1e-7) goto 1

      icount=0
      slum=0e0
      tmptot=0e0
c
c     jproc=1   njets>=0
c     g g -> Q Qbar + (njets) g
      if(jproc.eq.1) then
         slum = f1(0)*f2(0)
         ifl(1) = 0
         ifl(2) = 0
         do k=5,njets+4
            ifl(k)=0
         enddo
         itmp= 1
         effco = ccoef(1,njets+2)
         goto 100
c
c     jproc=2   njets>=0
c  2  q qbar -> Q Qbar + (njets) g
c     and qbar q -> Q Qbar + (njets) g   q=u,d,c,s
      else if(jproc.eq.2) then
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(-i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=-i
                  do k=5,njets+4
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(2,njets)
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=3   njets>=1
c  3  g q    -> Q Qbar q (njets-1) g 
c     and g qbar  -> Q Qbar qbar (njets-1) g
      else if(jproc.eq.3) then
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(0)*f2(i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=i
                  ifl(5)=i
                  do k=6,njets+4
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(2))
                  effco = ccoef(2,njets)
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=4   njets>=1
c  4  q g    -> Q Qbar (njets-1) g
c     and qbar g  -> Q Qbar qbar (njets-1) g
      else if(jproc.eq.4) then
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(0)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=0
                  ifl(5)=i
                  do k=6,njets+4
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(2,njets)
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=5   njets>=2
c  5  g g -> Q Qbar q qbar (njets-2 g) q=u,d,c,s
      else if(jproc.eq.5) then
         slum=f1(0)*f2(0)
c
c  slum x number of possible f.s. q qbar pairs
c
         slum=slum*dble(nlq)
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=1.
            endif
         enddo 
         icount=0
         rn=rn*2d0*dble(nlq)
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=0
                  ifl(5)=i
                  ifl(6)=-i
                  do k=7,njets+4
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(5))
                  effco = ccoef(2,njets)
                  goto 100
               endif 
            endif
         enddo 
c     jproc=6   njets>=2
c  6  q qbar -> Q Qbar q qbar (njets-2) g
c     and qbar q -> Q Qbar q qbar (njets-2) g   q=u,d,c,s
      else if(jproc.eq.6) then
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(-i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=-i
                  ifl(5)=i
                  ifl(6)=-i
                  do k=7,njets+4
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(3,njets-2)
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=7   njets>=2
c  7  q qbar -> Q Qbar q' qbar' (njets-2) g
c     and qbar q -> Q Qbar q' qbar' (njets-2) g   q=u,d,c,s
      else if(jproc.eq.7) then
         do i=-nlq,nlq
            if(i.ne.0) then
               do j=-nlq,nlq
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(-i)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               do j=-nlq,nlq
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=-i
                        ifl(5)=j
                        ifl(6)=-j
                        do k=7,njets+4
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,njets-2)
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
cc
c
c     jproc=8   njets>=2
c  8  q q' -> Q Qbar q q' (njets-2) g
c     and qbar qbar' -> Q Qbar qbar qbar' (njets-2) g   q=u,d,c,s
      else if(jproc.eq.8) then
         do i=-nlq,nlq
            if(i.ne.0) then
               do j=-nlq,nlq
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(j)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               do j=-nlq,nlq
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=j
                        ifl(5)=i
                        ifl(6)=j
                        do k=7,njets+4
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,njets-2)
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c
c     jproc=9   njets>=2
c  9  q qbar' -> Q Qbar q qbar' (njets-2) g
c     and qbar q' -> Q Qbar qbar q' (njets-2) g   q=u,d,c,s
      else if(jproc.eq.9) then
         do i=-nlq,nlq
            if(i.ne.0) then
               do j=-nlq,nlq
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.lt.0) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(j)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               do j=-nlq,nlq
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.lt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=j
                        ifl(5)=i
                        ifl(6)=j
                        do k=7,njets+4
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,njets-2)
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c
c     jproc=10   njets>=2
c 10  q q -> Q Qbar q q (njets-2) g
c     and qbar qbar -> Q Qbar qbar qbar (njets-2) g   q=u,d,c,s
      else if(jproc.eq.10) then
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=i
                  ifl(5)=i
                  ifl(6)=i
                  do k=7,njets+4
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(3,njets-2)
c statistical factor to account for indistinguishable quarks
                  slum= slum/2.d0
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=11   njets>=3
c  11 g q    -> Q Qbar q q' qbar' (njets-3) g 
c     and g qbar  -> Q Qbar qbar q' qbar' (njets-3) g
      else if(jproc.eq.11) then
         do i=-nlq,nlq
            if(i.ne.0) then
               do j=-nlq,nlq
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(0)*f2(i)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               do j=-nlq,nlq
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=0
                        ifl(2)=i
                        ifl(5)=j
                        ifl(6)=-j
                        ifl(7)=i
                        do k=8,njets+4
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(2))
                        effco = ccoef(3,njets-2)
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c
c     jproc=12   njets>=3
c  12 q g    -> Q Qbar q q' qbar' (njets-3) g 
c     and qbar g   -> Q Qbar qbar q' qbar' (njets-3) g
      else if(jproc.eq.12) then
         do i=-nlq,nlq
            if(i.ne.0) then
               do j=-nlq,nlq
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(0)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               do j=-nlq,nlq
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=0
                        ifl(5)=j
                        ifl(6)=-j
                        ifl(7)=i
                        do k=8,njets+4
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,njets-2)
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c
c     jproc=13   njets>=3
c  13 g q    -> Q Qbar q q qbar (njets-3) g 
c     and g qbar  -> Q Qbar qbar q qbar (njets-3) g
      else if(jproc.eq.13) then
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(0)*f2(i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=i
                  ifl(5)=i
                  ifl(6)=-i
                  ifl(7)=i
                  do k=8,njets+4
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(2))
                  effco = ccoef(3,njets-2)
c statistical factor to account for indistinguishable quarks
                  slum= slum/2.d0
                  goto 100
               endif
            endif
         enddo 
c  14 q g    -> Q Qbar q q qbar (njets-3) g 
c     and qbar g  -> Q Qbar qbar q qbar (njets-3) g
      else if(jproc.eq.14) then
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(0)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlq,nlq
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=0
                  ifl(5)=i
                  ifl(6)=-i
                  ifl(7)=i
                  do k=8,njets+4
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(3,njets-2)
c statistical factor to account for indistinguishable quarks
                  slum= slum/2.d0
                  goto 100
               endif
            endif
         enddo 
c     jproc=15   njets>=4
c  15 g g    -> Q Qbar q qbar q qbar (njets-4) g 
c     and qbar g  -> Q Qbar q qbar q qbar (njets-4) g
      else if(jproc.eq.15) then
         do i=1,4
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(0)*f2(0)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=1,4
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
               ifl(1)=0
               ifl(2)=0
               ifl(5)=i
               ifl(6)=-i
               ifl(7)=i
               ifl(8)=-i
               do k=9,njets+4
                  ifl(k)=0
               enddo
               itmp= 1
               effco = ccoef(3,njets-2)
c statistical factor to account for indistinguishable quarks
               slum= slum/4.d0
               goto 100
            endif
         enddo 
c
c     jproc=16   njets>=4
c  16 g g    -> Q Qbar q qbar q' qbar' (njets-4) g 
      else if(jproc.eq.16) then
         do i=1,nlq
            do j=1,nlq
               if(i.ne.j) then
                  icount=icount+1
                  tmp(icount)=f1(0)*f2(0)
                  slum=slum+tmp(icount)
               endif
            enddo
         enddo 
         icount=0
         rn=rn*slum
         do i=1,nlq
            do j=1,nlq
               if(i.ne.j) then
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=0
                     ifl(2)=0
                     ifl(5)=i
                     ifl(6)=-i
                     ifl(7)=j
                     ifl(8)=-j
                     do k=9,njets+4
                        ifl(k)=0
                     enddo
                     itmp= 1
                     effco = ccoef(3,njets-2)
                     goto 100
                  endif
               endif
            enddo
         enddo 
      else
         write(*,*) 'jproc not defined, slum=0'
         slum=0
         stop
      endif
c
      xlum=-1d0
      return
 100  continue

      ib=3
      ibb=4
      ifl(3)=ihvy
      ifl(4)=-ifl(3)
      if(itmp.lt.0) then
         iflaux=ifl(3)
         ifl(3)=ifl(4)
         ifl(4)=iflaux
         ib=4
         ibb=3 
      endif
c
      do i=1,maxpar
         afl(i)=ifl(i)
      enddo
c
c     evaluate colour weight factors
      do i=1,2
         if(ifl(i).eq.0) ifl(i)=21
      enddo
      ng=0
      cwgt=1e0
      do i=3,njets+4
         cwgt=cwgt*cfac(abs(ifl(i)))
         if(ifl(i).eq.0) then
            ifl(i)=21
            ng=ng+1
         endif
      enddo
c     evaluate spin weight factors
      swgt=2e0
      swgt=swgt**(njets+2)
*
      xlum=dble(slum*cwgt*swgt*ifact(ng)) /resc**(njets+2)
      xlum=xlum*effco
c
c     decaying top quarks
c
      if(idecay.eq.'y'.and.ihvy.eq.6) then 
         do i=1,4
            do j=1,3
               do k=1,2
                  idec(i,j,k)=0.d0
               enddo
            enddo
         enddo

         jtl= 1
         jtbl= 2
         jt= 3
         jtb= 4
         if(itmp.lt.0) then
            jt= 4
            jtb= 3
            jtl= 2
            jtbl= 1
         endif
         ptpb(1)= p(4,jtb)
         ptp(1)= p(4,jt)
         do i=1,3
            ptpb(i+1)= p(i,jtb)
            ptp(i+1)= p(i,jt)
         enddo
         call vtop(ptpb,mq,mb,mw,vtp)
         do i=1,4
            idec(i,1,jtbl)=pb(exch(i))
            idec(i,2,jtbl)=pf(exch(i))
            idec(i,3,jtbl)=pfb(exch(i))
         enddo
         call utopb(ptp,mq,mb,mw,utpb)
         do i=1,4
            idec(i,1,jtl)=pb(exch(i))
            idec(i,2,jtl)=pf(exch(i))
            idec(i,3,jtl)=pfb(exch(i))
         enddo
      endif
      call usrfll(ib,ibb)
      end
c
c
      subroutine usrfll(ib,ibb)
      implicit none
      include 'alpgen.inc'
      include '2Q.inc'
      integer i,j,ib,ibb
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c
      do i=1,4
         do j=1,2
            pin(i,j)=p(i,j)
         enddo
         do j=1,npart-2
            pout(i,j)=p(i,j+2)
         enddo
      enddo

      if(ihvy.le.5) then
         do i=1,4
            do j=1,njets+2
               pjet(i,j)=p(i,j+2)
               ptj(j)=pt(j+2)
               etaj(j)=eta(j+2)
            enddo
         enddo
         do j=1,njets+2
            do i=j+1,njets+2
               drjj(i,j)=dr(i+2,j+2)
               drjj(j,i)=drjj(i,j)
            enddo
         enddo
         drbb=drjj(1,2)
         do i=1,4
            pbott(i)=p(i,ib)
            pbbar(i)=p(i,ibb)
         enddo
         ptb=pt(ib)
         ptbb=pt(ibb)
         etab=eta(ib)
         etabb=eta(ibb)
         do i=1,njets
            drbj(i)=dr(ib,4+i)
            drbbj(i)=dr(ibb,4+i)
         enddo
      elseif(ihvy.eq.6) then
         do i=1,4
            ptop(i)=p(i,ib)
            ptbar(i)=p(i,ibb)
         enddo
         do i=1,4
            do j=1,njets
               pjet(i,j)=p(i,j+4)
               ptj(j)=pt(j+4)
               etaj(j)=eta(j+4)
            enddo
         enddo
         do j=1,njets
            do i=j+1,njets
               drjj(i,j)=dr(i+4,j+4)
               drjj(j,i)=drjj(i,j)
            enddo
         enddo
      endif
      end
c
      subroutine setdec(nwrt,iflwrt,icuwrt,pwrt,decbr)
      implicit none
      include 'alpgen.inc'
      include '2Q.inc'
c debug
      double precision rate(16)
      common/dbg/rate
      data rate/16*0d0/
c locals
      integer maxdec
      parameter (maxdec=40)
      integer ip, ic, il, irn,id1,id2,i,itmp,ich,idecmode
     $     ,itdec0
      integer iwfl(2,2)
      integer init
      data init/0/
      integer icab(-4:4,2) 
      double precision xrn,rangen2
      double precision pdec(5),m1,m2
c     cabibbo partner: (*,1)=cabibbo allowed, (*,2)-cabibbo suppressed
c     where 1=d 2=u 3=s 4=c 0=g and negatives are antiparticles
      data icab/-3,-4,-1,-2,0,2,1,4,3, 
     +          -1,-2,-3,-4,0,4,3,2,1/
      integer ilepc(3),inu(3),iup(2),idn(2)
      data ilepc/11,13,15/,inu/12,14,16/,iup/2,4/,idn/1,3/
c BR's
      double precision br(7)
      data br/3*0.148148148,0.444444444,0.111111111,0.444444444,1d0/
      double precision dmass(16)
      data dmass/16*0d0/
c arguments
      integer nwrt,iflwrt(maxdec),icuwrt(2,maxdec)
      double precision pwrt(5,maxdec),decbr
c
      save idecmode,dmass
c
      do ip=1,npart
        iflwrt(ip)=ifl(ip)
        do ic=1,2
          icuwrt(ic,ip)=icu(ic,ip)
        enddo
        do il=1,5
          pwrt(il,ip)=p(il,ip)
        enddo
      enddo
      nwrt=npart
      decbr=1d0
      if(abs(ifl(3)).ne.6) return
      if(idecay.ne.'y') return
c
c     add the W decay products; the W itself and the b will be
C     reconstructed  from momentum conservation by the Les Houches
C     interface
c
      if(init.eq.0) then
        idecmode=itdecmode
c        write(*,*) 'select top decay modes:'
c        write(*,*) '1: e nu b bbar + 2 jets'
c        write(*,*) '2: mu nu b bbar + 2 jets'
c        write(*,*) '3: tau nu b bbar + 2 jets'
c        write(*,*) '4: e/mu/tau nu b bbar + 2 jets'
c        write(*,*) '5: l nu l'' nu  b bbar (l,l''=e/mu/tau)'
c        write(*,*) '6: b bbar + 4 jets'
c        write(*,*) '7: fully inclusive'
c        read(*,*) idecmode
        dmass(4)=1.5
        dmass(5)=mb
        dmass(11)=0.5d-3
        dmass(13)=0.10566d0
        dmass(15)=1.777d0
        init=1
      endif
      decbr=br(idecmode)
c     randomize decays between t and tbar:
      xrn=rangen2(1)
      if(2*xrn.lt.1) then
        id1=jtl
        id2=jtbl
        ich=1
      else
        id1=jtbl
        id2=jtl
        ich=0
      endif
c freeze top dec mode:
      itdec0=idecmode
c     select flavours of W decay products
 1    if(idecmode.le.3) then
        iwfl(1,id1)=ilepc(idecmode)+ich   ! tbar->l- or t->nu
        iwfl(2,id1)=-inu(idecmode)+ich  !  tbar->nu_lbar or t->l+
 11     xrn=rangen2(1)
        irn=1+int(2*xrn)
        if(irn.gt.2) goto 11
        itmp=iup(irn)-ich 
        iwfl(1,id2)=itmp      ! t->u/c or tbar->d/s
        xrn=rangen2(1)
        if(xrn.lt.scab2) then
          iwfl(2,id2)=-icab(itmp,2)
        else
          iwfl(2,id2)=-icab(itmp,1)
        endif
      elseif(idecmode.eq.4) then
 12     xrn=rangen2(1)
        irn=1+int(3*xrn)
        if(irn.gt.3) goto 12
        iwfl(1,id1)=ilepc(irn)+ich   ! tbar->l- or t->nu
        iwfl(2,id1)=-inu(irn)+ich  !  tbar->nu_lbar or t->l+
 13     xrn=rangen2(1)
        irn=1+int(2*xrn)
        if(irn.gt.2) goto 13
        itmp=iup(irn)-ich 
        iwfl(1,id2)=itmp      ! t->u/c or tbar->d/s
        xrn=rangen2(1)
        if(xrn.lt.scab2) then
          iwfl(2,id2)=-icab(itmp,2)
        else
          iwfl(2,id2)=-icab(itmp,1)
        endif
      elseif(idecmode.eq.5) then
 14     xrn=rangen2(1)
        irn=1+int(3*xrn)
        if(irn.gt.3) goto 14
        iwfl(1,id1)=ilepc(irn)+ich   ! tbar->l- or t->nu
        iwfl(2,id1)=-inu(irn)+ich  !  tbar->nu_lbar or t->l+
 15     xrn=rangen2(1)
        irn=1+int(3*xrn)
        if(irn.gt.3) goto 15
        iwfl(1,id2)=inu(irn)-ich   ! t->nu or tbar->l-
        iwfl(2,id2)=-ilepc(irn)-ich  !   t->l+ or tbar->nubar
      elseif(idecmode.eq.6) then
 16     xrn=rangen2(1)
        irn=1+int(2*xrn)   !irn=1,2
        if(irn.gt.2) goto 16
        itmp=idn(irn)+ich 
        iwfl(1,id1)=itmp      ! tbar->d/s or t->u/c
        xrn=rangen2(1)
        if(xrn.lt.scab2) then
          iwfl(2,id1)=-icab(itmp,2)
        else
          iwfl(2,id1)=-icab(itmp,1)
        endif
 17     xrn=rangen2(1)
        irn=1+int(2*xrn)   !irn=1,2
        if(irn.gt.2) goto 17
        itmp=iup(irn)-ich 
        iwfl(1,id2)=itmp      ! t->u/c or tbar->d/s
        xrn=rangen2(1)
        if(xrn.lt.scab2) then
          iwfl(2,id2)=-icab(itmp,2)
        else
          iwfl(2,id2)=-icab(itmp,1)
        endif
      elseif(idecmode.eq.7) then
        xrn=rangen2(1)
        if(xrn.lt.0.44444444444) then
          idecmode=4
        elseif(xrn.lt.0.555555555) then 
          idecmode=5
        else
          idecmode=6
        endif
        goto 1
      endif
c     restore dec mode label
      idecmode=itdec0
c
      do i=1,2
        if(i.eq.1) then
          itmp=jtl
        else 
          itmp=jtbl
        endif
        nwrt=nwrt+2
        do il=1,4
          pwrt(il,nwrt-1)=idec(il,2,itmp)
          pwrt(il,nwrt)=idec(il,3,itmp)
          pdec(il)=pwrt(il,nwrt-1)+pwrt(il,nwrt)
        enddo
        iflwrt(nwrt-1)=iwfl(1,itmp)
        iflwrt(nwrt)=iwfl(2,itmp)
        m1=dmass(abs(iflwrt(nwrt-1)))
        m2=dmass(abs(iflwrt(nwrt)))
        pwrt(5,nwrt-1)=0
        pwrt(5,nwrt)=0
        pdec(5)=mw
c rescale momenta to incorporate mass effects in the decays
        call rescms(pdec,pwrt(1,nwrt-1),pwrt(1,nwrt),m1,m2)
c
        if(abs(iflwrt(nwrt-1)).le.6) then
          icuwrt(1,nwrt-1)=nwrt-1
          icuwrt(2,nwrt)=nwrt-1
        else
          icuwrt(1,nwrt-1)=0
          icuwrt(2,nwrt)=0
        endif
        icuwrt(2,nwrt-1)=0
        icuwrt(1,nwrt)=0
      enddo  
c$$$c  debug
c$$$      do i=nwrt-3,nwrt
c$$$        do il=1,16
c$$$          if(abs(iflwrt(i)).eq.il) rate(il)=rate(il)+1d0
c$$$        enddo
c$$$      enddo

      end
*           
      subroutine phspace(lnot,pswgt,djpd,djg)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c  (njets)-body phase space for the process:                     c
c                                                                c
c  h(1)  h(2) -> q(3) anti_q(4) h(5) ... h(njets)                c
c                                                                c
c  where h is any light-flavoured quark or gluon                 c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'alpgen.inc'
      include '2Q.inc'
      real*8 dummy,djg
      real *8 pswgt
      real *8 djpd,factor
      real *8 cutkin(10)
      real *8 wgt
      common/loccut/cutkin
      real *8 pl(maxpar),y(maxpar)
      real *8 pcm(0:3,maxpar)
*
c-    debugging variables
*
      real *8 totpt
      integer nx,ninit,i,l,m,lnot,ndummy
      parameter (nx= 20)
      data ninit/0/
      save ninit
      if(ninit.eq.0) then
*
c-       setup local generation cuts
*
         cutkin(1)=ptjmin
         cutkin(2)=ptQmin
         cutkin(3)=etajmax
         cutkin(4)=etaQmax
         cutkin(5)=drjmin
         cutkin(6)=drQmin
         cutkin(7)=ptj1min
         cutkin(8)=ptj1max
         cutkin(9)=0.d0
         cutkin(10)=0.d0
c-
         ninit=1
      endif
*
c-    The generation starts
*
      pswgt=0.d0
      if(ihvy.eq.6) then
         call momgen_ttj(njets,mq,roots,x1,x2,pcm,
     +                   wgt,lnot)
      else if(ihvy.le.5) then
         call momgen_bbj(njets,mq,roots,x1,x2,pcm,
     +                   wgt,lnot)
      endif
      djg= 1.d0 ! dummy variable
      if (lnot.eq.1) then
         pswgt= 0.d0
         goto 100
      endif
*
c-    will write factor=factor0/(x1*x2), with factor0 function of njets etc.
*
      factor= 1d0/(2.d0*pi)**(3*(njets+2)-4)/2.d0/s/x1/x2
*
c-    initial state momenta in the LAB frame:
*
      pcm(0,1)= roots/2.d0*x1   
      pcm(1,1)= 0.d0   
      pcm(2,1)= 0.d0   
      pcm(3,1)= roots/2.d0*x1   

      pcm(0,2)= roots/2.d0*x2   
      pcm(1,2)= 0.d0   
      pcm(2,2)= 0.d0   
      pcm(3,2)=-roots/2.d0*x2   

c-    Rapidities and pseudo-rapidities (in the LAB system), pt's, deltar's:

      do l= 3,njets+4
          pt(l) = sqrt(pcm(1,l)**2+pcm(2,l)**2)
          pl(l) = pcm(3,l)
          eta(l)= -log(tan(0.5d0*atan2(pt(l),pl(l))))
          y(l)  = 0.5d0*log((pcm(0,l)+pcm(3,l))/(pcm(0,l)-pcm(3,l)))
      enddo
*
c-    Calculates jet-jet distances:
*
      do l= 3,njets+3
         do m= l+1,njets+4
            if(min(pt(l),pt(m)).gt.0d0) then
               dphi(l,m)= 
     +              (pcm(1,l)*pcm(1,m)+pcm(2,l)*pcm(2,m))
     +              /pt(l)/pt(m)
               if(dabs(dphi(l,m)).gt.1.d0) then
c                  write(*,*) 'partons',l,m,', cos(Dphi)=', dphi(l,m)
c     +                 ,', set to +-1'
c                  write(*,*) 'pt(',l,')=',pt(l),'p(',l,')=',  (pcm(i,l)
c     +                 ,i=0,3)
c                  write(*,*) 'pt(',m,')=',pt(m),'p(',m,')=',  (pcm(i,m)
c     +                 ,i=0,3)
                  if (dphi(l,m).gt.0.d0) dphi(l,m)= 1.d0
                  if (dphi(l,m).lt.0.d0) dphi(l,m)=-1.d0
               endif
               dphi(l,m)=acos(dphi(l,m))
               dphi(m,l)=dphi(l,m)
            else
c***fix***
c                dphi(l,m)=pi
c                dphi(m,l)=pi
               goto 100
            endif
         enddo
      enddo
*
      do l= 3,njets+3
         do m= l+1,njets+4
            dr(l,m)= sqrt(dphi(l,m)**2+(eta(l)-eta(m))**2)
            dr(m,l)=dr(l,m)
         enddo
      enddo
*
c-    Redefine the momenta:           
*                                     
      do l= 1,njets+4                 
         do m=1,2                     
            p(m,l)=pcm(m,l)           
         enddo                        
         p(4,l)= pcm(0,l)             
         p(3,l)= pcm(3,l)             
      enddo                           
comment
c      do l= 1,npart
c        print*,'osh',l,'=',p(4,l)**2-p(1,l)**2-p(2,l)**2-p(3,l)**2
c      enddo     
c
c      do l= 1,4
c        dummy= -p(l,1)-p(l,2)
c        do m= 3,npart
c          dummy= dummy+p(l,m)
c        enddo
c        print*,'summ',l,'=',dummy
c      enddo
c      print*,'   '
comment
*                                     
c-    Initial state pt, eta, and y:         
*                                     
      do l= 1,2                       
         pt(l) = 0.d0                 
         y(l)  = 1.d6                 
         eta(l)= y(l)                 
      enddo                           
*
c-    Call to the the cut routine:
*
      if(ihvy.eq.6) call chkcut_ttj(lnot,pt,p,eta,dr,njets)
      if(ihvy.le.5) call chkcut_bbj(lnot,pt,p,eta,dr,njets)

      if (lnot.eq.1) then
         pswgt= 0.d0
         goto 100
      endif
*
      pswgt = factor*wgt
*
c-    Evaluate q2: will include several possible options:
*
      totpt= 0.d0
*
c-    Total et^2 of jets
*
      do i=3,njets+4
        totpt=totpt+pt(i)**2+p(5,i)**2
      enddo 
      if(iqopt.eq.0) then 
        qsq=1d0
      elseif(iqopt.eq.1) then
        qsq=totpt
      elseif(iqopt.eq.2) then
        qsq=roots*roots*x1*x2
      endif
      qsq=qfac**2*qsq
      
 100  continue
      end
*
      subroutine chkcut_ttj(lnot,pt,p,eta,dr,njets)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     Applies kinematical cuts to the final state during the phase
c     -space generation                                          c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 ptjmin,drQmin,drjmin,etaQmax,etajmax,ptQmin
      real*8 ptj1min,ptj1max,highjet
      integer maxpar,ninit,njets,j,lnot,i
      real*8 cutkin(10)
      common/loccut/cutkin
      data ninit/0/
      parameter (maxpar=10)
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),p(5,maxpar)
      save ninit,ptjmin,ptQmin,etajmax,etaQmax,drjmin,drQmin
      save ptj1min,ptj1max
      if(ninit.eq.0) then
         ninit=1
         ptjmin=cutkin(1)
         ptQmin=cutkin(2)
         etajmax=cutkin(3)
         etaQmax=cutkin(4)
         drjmin=cutkin(5)
         drQmin=cutkin(6)
         ptj1min=cutkin(7)
         ptj1max=cutkin(8)
      endif

      lnot= 0

c     impose leading jet cut
      highjet = pt(3)
      do i=4,njets+4
         if(pt(i).gt.highjet) highjet=pt(i)
      enddo
      if (ptj1max.gt.0.and.highjet.gt.ptj1max) goto 10
      if (ptj1min.gt.0.and.highjet.lt.ptj1min) goto 10

      if(njets.gt.0) then
c     impose minimum pt and require eta within allowed range, for jets
         do i=5,njets+4
            if (pt(i).lt.ptjmin)           goto 10
            if (abs(eta(i)).gt.etajmax)    goto 10
         enddo 
c     require dR(jet-jet)<drjmin
         do i=5,njets+3
            do j=i+1,njets+4
               if(dr(i,j).lt.drjmin)  goto 10
            enddo
         enddo
      endif 

 5    return

 10   lnot= 1
      return
      end

      subroutine chkcut_bbj(lnot,pt,p,eta,dr,njets)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     Applies kinematical cuts to the final state during the phase
c     -space generation                                          c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 ptjmin,drQmin,drjmin,etaQmax,etajmax,ptQmin
      integer maxpar,ninit,njets,j,lnot,i
      real*8 cutkin(10),highjet,ptj1min,ptj1max
      common/loccut/cutkin
      data ninit/0/
      parameter (maxpar=10)
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),p(5,maxpar)
      save ninit,ptjmin,ptQmin,etajmax,etaQmax,drjmin,drQmin
      save ptj1min,ptj1max
      if(ninit.eq.0) then
         ninit=1
         ptjmin=cutkin(1)
         ptQmin=cutkin(2)
         etajmax=cutkin(3)
         etaQmax=cutkin(4)
         drjmin=cutkin(5)
         drQmin=cutkin(6)
         ptj1min=cutkin(7)
         ptj1max=cutkin(8)
      endif

      lnot= 0

c     impose leading jet cut
      highjet = pt(3)
      do i=4,njets+4
         if(pt(i).gt.highjet) highjet=pt(i)
      enddo
      if (ptj1max.gt.0.and.highjet.gt.ptj1max) goto 10
      if (ptj1min.gt.0.and.highjet.lt.ptj1min) goto 10

c     impose minimum pt and require eta within allowed range, for b's
      do i=3,4
         if (pt(i).lt.ptQmin)           goto 10
         if (abs(eta(i)).gt.etaQmax)    goto 10
      enddo 
c     require dR(b-bbar)>drQmin
      if(dr(3,4).lt.drQmin)      go to 10

      if(njets.gt.0) then
c     impose minimum pt and require eta within allowed range, for jets
         do i=5,njets+4
            if (pt(i).lt.ptjmin)           goto 10
            if (abs(eta(i)).gt.etajmax)    goto 10
         enddo 
c     require dR(b-jet)<drjmin, dR(jet-jet)<drjmin
         do i=5,njets+4
            if(dr(3,i).lt.drjmin)        goto 10
            if(dr(4,i).lt.drjmin)        goto 10
         enddo
         do i=5,njets+3
            do j=i+1,njets+4
               if(dr(i,j).lt.drjmin)  goto 10
            enddo
         enddo
      endif 
 5    return

 10   lnot= 1
      return
      end
*
      subroutine momgen_bbj(njets,qm,roots,x1,x2,p,wgt,lw)
*
c-    Generator of np particles in the LAB frame.
*
      implicit none
      real*8 roots,x1,x2,wgt,wgt1,wgt2,zero,etacut
      real*8 ptjmin,ptQmin,etajmax,etaQmax,drjmin,drQmin
      real*8 qm2,qm,cutkin
      double precision s,xmsum0,pt0lmax,pt0lsum,eta0min
      double precision wtau,wjr,en,pz
      double precision ag,tau0,tau
      real*8 rootsh,ranram,sq,y0,yr,cxmb
      real*8 cbmin,rcnmin,rcxml,cxml
      real*8 xmsum,ptlim
      integer lw,npar,mpar,np,njets,lflag,lim,j,l,k,lw1
      parameter (npar= 20)
      double precision pm(0:4,npar),pt0(npar),pt1(npar),xm(npar),
     +       xmr(npar),eta0(npar),pt0_l(npar),pt1_l(npar),p(0:3,npar)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/loccut/cutkin(10)
      dimension ranram(111)
      data mpar/0/
      real*8 djb0,djb1,apw(1:2),ran0
      real *8 djbin,djbintot,dummy
      integer nvar,nbin,nv,m,ndummy,md
      common/psopt/nvar,nv
      integer nct,nx1,nx2,maxn,lmin,lmax
      parameter (maxn= 100)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      integer mask,mmask
      real*8 peropt
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      save
*
      if ((njets.gt.6).or.(njets.lt.0)) then
        write (6,*) 'WRONG NJETS VALUE'
        stop
      endif
      if (mpar.eq.0) then
        mpar   = 1
        np     = njets+2
        lflag  = 0
        lim    = 1
        zero   = 0.d0
        qm2    = qm*qm
        etacut = 40.d0
        s      = roots*roots
*
        ptjmin=cutkin(1)
        ptQmin=cutkin(2)
        etajmax=cutkin(3)
        etaQmax=cutkin(4)
        drjmin=cutkin(5)
        drQmin=cutkin(6)
*
       do j= 1,2
         xm(j)    = qm        
         pt0_l(j) = ptQmin
         pt1_l(j) = roots/2.d0
         eta0(j)  = etaQmax
       enddo
*
        do j= 3,np
          xm(j)    = zero        
          pt0_l(j) = ptjmin
          pt1_l(j) = roots/2.d0
          eta0(j)  = etajmax
        enddo
*
c-      Parameters for the x1,x2 integration:
*
        ag  = 0.98d0
*
        cbmin= (2.d0*ptQmin*sin(drQmin/2.d0))**2
        rcnmin= 2.d0*ptjmin*sin(drjmin/2.d0)
*
c-      l is the number of light jets
*
        l= njets                               
        rcxml= int(l/2)*rcnmin
        cxml= rcxml*rcxml
        cxmb = max(4.d0*qm2,cbmin)
        xmsum0 = 0.d0
        pt0lmax= 0.d0
        pt0lsum= 0.d0
        eta0min= 100.d0
        do j= 1,np
           xmsum0 = xmsum0+xm(j)
           pt0lmax= max(pt0lmax,pt0_l(j))
           pt0lsum= pt0lsum+pt0_l(j)
           eta0min= min(eta0min,eta0(j))
        enddo
        xmsum= xmsum0
*
c-      tau0 is the lower cut on tau= x1*x2
*
        tau0= 1.d0/s*(xmsum)**2
        tau0= max(tau0,pt0lsum**2/s)
        tau0= max(tau0,1.d0/s*
     +          (2d0*qm+float(njets)*ptjmin)**2)
*
        tau0= max(tau0,cxmb/s)
        tau0= max(tau0,cxml/s)
*
c-      sets the ratios of calls mom/momr
*
        apw(1)= 0.9d0
        apw(2)= 1.d0-apw(1)
      endif
*
      lw= 0
      lmin= (jgrid(jproc)-1)*nv+2
      lmax= jgrid(jproc)*nv+1
*
      call rans(ran0)
      if (ran0.le.apw(1)) then
         md= 1
      else
         md= 2 
      endif
*
      djb1= 1.d0
      djb0= 1.d0
      m= 0
      do l= lmin,lmax
        m= m+1
        if (md.eq.1) then
          mmask(l)= 1
        else
          mmask(l)= mask(l)
        endif
        if (mmask(l).eq.1) then
          call onedimbin(1,nbin,djbin,l,ndummy,dummy)
          call grans(nbin,nx1(l),ranram(m))
          if (mask(l).eq.1) then
            djb1= djb1*djbin
          elseif (mask(l).eq.0) then
            djb0= djb0*djbin
          endif
        else
          ranram(m)= 0.d0
        endif
      enddo
      djbintot= djb1
      if (md.eq.1) djbintot= djbintot*djb0 
*
c-    Generates x1 and x2: 
*
      call ppeaka(0,tau0,ag,tau0,1.d0,tau,wtau,ranram(1),lw1)
      sq = sqrt(tau)
      y0 =-0.5d0*log(tau) 
      call ttriangle(0,-y0,y0,-y0,y0,yr,wjr,ranram(2),lw1)
*
      x1 = sq*exp(yr) 
      x2 = sq*exp(-yr) 
      rootsh= roots*sq
*
c-    Protection:
*
      if (rootsh.lt.(xmsum+pt0lsum)) goto 100
*
c-    Rescalings and transformations to feed mom and etas:
*
      do j= 1,np
        ptlim = s*s
     +           +(xm(j)**2-(xmsum-xm(j))**2)**2
     +        -2.d0*s*(xm(j)**2+(xmsum-xm(j))**2)   
        if (ptlim.le.0.d0) goto 100
        ptlim =  0.5d0/roots*sqrt(ptlim)
        pt0(j)= pt0_l(j)/roots
        pt1(j)= min(pt1_l(j),ptlim)/roots
*
c-      Protection:
*
        if (pt0(j).gt.pt1(j)) goto 100
        xmr(j)= xm(j)/roots
      enddo
*
c-    momenta:
*
      en    = 0.5d0*(x1+x2)      ! Rescaled initial total energy
      pz    = 0.5d0*(x1-x2)      ! Rescaled initial longitudinal momentum
*
      if (md.eq.1) then
        call mom(0,np,pt0,pt1,eta0,xmr,en,pz,pm,wgt1,ranram,lw)
        if (lw.eq.0) then
          call momr(1,np,xmr,en,pz,pm,wgt2,lw)
          if (lw.eq.0) then
            if (wgt1.eq.0.d0.and.wgt2.eq.0.d0) goto 100
            wgt= wgt1*wgt2/
     +          (wgt2*apw(1)*djbintot+wgt1*apw(2)*djb1)
          else
            wgt= wgt1/(apw(1)*djbintot)
          endif
        else
          goto 100
        endif
      else
        call momr(0,np,xmr,en,pz,pm,wgt2,lw)
        if (lw.eq.0) then
          call mom(1,np,pt0,pt1,eta0,xmr,en,pz,pm,wgt1,ranram,lw)
          if (lw.eq.0) then
            m= 0
            do l= lmin,lmax
              m= m+1
              if (mask(l).eq.0) then
                nbin= min(nx1(l),int(1.d0+dfloat(nx1(l))*ranram(m)))
                call onedimbin(3,nbin,djbin,l,ndummy,dummy)
                djb0= djb0*djbin
              endif
            enddo
            if (wgt1.eq.0.d0.and.wgt2.eq.0.d0) goto 100
            wgt= wgt1*wgt2/
     +          (wgt1*apw(2)*djbintot+wgt2*apw(1)*djb0*djb1)
          else
            wgt= wgt2/(apw(2)*djbintot)
          endif
        else
          goto 100
        endif
      endif
*
c-    from the x1,x2 integration
*
      wgt= wgt*wtau*wjr
*
c-    Rescaling of momenta and weight:
*
      do j= 1,np
        do k= 0,3
           p(k,j+2)= pm(k,j)*roots
        enddo
      enddo
*
      wgt= wgt*s**(np-2)
*
      return
 100  lw= 1
      wgt= 0.d0
      return
      end
*
      subroutine momgen_ttj(njets,qm,roots,x1,x2,p,wgt,lw)
*
c-    Generator of np particles in the LAB frame.
*
      implicit none
      real*8 qm,roots,x1,x2,wgt,wgt1,wgt2,zero,qm2,etacut
      real*8 ptjmin,cutkin,etajmax,drjmin,s,sw
      real*8 ag,tau0,tau,xmsum,xmsum0,ptlim,cntt
      real*8 rootsh,sq,y0,yr,cxmw,cxpw,dj1,dj2,wtau,wjr
      real*8 pt0lmax,pt0lsum,en,pz,eta0min,ran0,ran(1:2)
      integer lw,lw1,npar,mpar,np,njets,nw,j,k
      parameter (npar= 20)
      real*8 pm(0:4,npar),pt0(npar),pt1(npar),xm(npar),xmr(npar),
     +       eta0(npar),pt0_l(npar),pt1_l(npar),p(0:3,npar)
      real*8 ranram(111)
      common/loccut/cutkin(10)
      integer nvar,nv,m,l,nbin,ndummy,mask,mmask,md
      integer nct,nx1,nx2,maxn,lmin,lmax
      parameter (maxn= 100)
      real*8 dj,djbintot,djbin,dummy,apw(1:2)
      real*8 djb0,djb1
      common/psopt/nvar,nv
      real*8 peropt
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      data mpar/0/
      save
*
      if (mpar.eq.0) then
        mpar   = 1
        np     = njets+2
        nw     = np-1
        zero   = 0.d0
        qm2    = qm*qm
        etacut = 40.d0
        s      = roots*roots
*
        ptjmin = cutkin(1)
        etajmax= cutkin(3)
        drjmin = cutkin(5)
*
c-      Parameters for the light jets:
*
        do j= 1,nw-1
          pt0_l(j) = ptjmin
          pt1_l(j) = roots/2.d0
          eta0(j)  = etajmax
          xm(j)    = zero        
        enddo
*
c-      Parameters for the the ttbar system:
*
        cntt      = 1.d0
        pt0_l(nw) = zero
        pt1_l(nw) = roots/2.d0
        eta0(nw)  = etacut
        xm(nw)    = zero  ! provisional ttbar mass!!!
*
        cxmw      = 4.d0*qm2
        cxpw      = s
*
c-      Parameters for the x1,x2 integration:
*
        ag  = 0.98d0
*
        xmsum0 = 0.d0
        pt0lmax= 0.d0
        pt0lsum= 0.d0
        eta0min= 100.d0
        do j= 1,nw
          xmsum0 = xmsum0+xm(j)
          pt0lmax= max(pt0lmax,pt0_l(j))
          pt0lsum= pt0lsum+pt0_l(j)
          eta0min= min(eta0min,eta0(j))
        enddo
*
c-      sets the ratios of calls mom/momr
*
        apw(1)= 0.9d0
        apw(2)= 1.d0-apw(1)
      endif
*
      lw  = 0
      lmin= (jgrid(jproc)-1)*nv+2
      lmax= jgrid(jproc)*nv+1
*
      dj= 0.d0
      call rans(ran0)
      if (ran0.le.apw(1)) then
         md= 1
      else
         md= 2 
      endif
*
      djb1= 1.d0
      djb0= 1.d0
      m= 0
      do l= lmin,lmax
        m= m+1
        if (md.eq.1) then
          mmask(l)= 1
        else
          mmask(l)= mask(l)
        endif
        if (mmask(l).eq.1) then
          call onedimbin(1,nbin,djbin,l,ndummy,dummy)
          call grans(nbin,nx1(l),ranram(m))
          if (mask(l).eq.1) then
            djb1= djb1*djbin
          elseif (mask(l).eq.0) then
            djb0= djb0*djbin
          endif
        else
          ranram(m)= 0.d0
        endif
      enddo
      djbintot= djb1
      if (md.eq.1) djbintot= djbintot*djb0 
*
c-    njets= 0 :
*
      if (njets.eq.0) then
        call photsm(0,cntt,cxmw,cxpw,sw,dj1,ranram(3))
        tau= sw/s
        y0 =-0.5d0*log(tau) 
        sq = sqrt(tau)
        call ttriangle(0,-y0,y0,-y0,y0,yr,wjr,ranram(1),lw1)
        x1 = sq*exp(yr) 
        x2 = sq*exp(-yr) 
*
        pm(0,nw)= roots/2.d0*(x1+x2)
        pm(1,nw)= 0.d0
        pm(2,nw)= 0.d0
        pm(3,nw)= roots/2.d0*(x1-x2)
*
        ran(1)= ranram(2)
        call rans(ran0)
        ran(2)= ran0
        call dec2fm(0,sw,pm(0,nw),
     +              qm2,qm2,p(0,nw+2),p(0,nw+3),dj2,ran)
*
        wgt= wjr/s/dj1/dj2/djbintot        
        return
      endif
*
c-    Generates sw=stt:
*
      call photsm(0,cntt,cxmw,cxpw,sw,dj1,ranram(2*nw+1))
*
c-    The mass of the last particle and the upgrading of xmsum
*
      xm(nw)= sqrt(sw)
      xmsum = xmsum0+xm(nw)
*
c-    tau0 is the lower cut on tau= x1*x2
*
      tau0   = 1.d0/s*(xmsum)**2
      tau0   = max(tau0,4.d0*pt0lmax**2/s)  
      tau0   = max(tau0,pt0lsum**2/s)
      tau0   = max(tau0,1.d0/s*(2.d0*qm
     +                  +float(njets)*ptjmin)**2)
*
c-    Generates x1 and x2: 
*
      call ppeaka(0,tau0,ag,tau0,1.d0,tau,wtau,ranram(1),lw1)
      sq = sqrt(tau)
      y0 =-0.5d0*log(tau) 
      call ttriangle(0,-y0,y0,-y0,y0,yr,wjr,ranram(2),lw1)
*
      x1 = sq*exp(yr) 
      x2 = sq*exp(-yr) 
      rootsh= roots*sq
*
c-    Protection:
*
      if (rootsh.lt.xmsum) goto 100
*
c-    Rescalings and transformations to feed mom:
*
      do j= 1,nw
         ptlim = s*s
     +               +(xm(j)**2-(xmsum-xm(j))**2)**2
     +        -2.d0*s*(xm(j)**2+(xmsum-xm(j))**2)   
         if (ptlim.le.0.d0) goto 100
         ptlim =  0.5d0/roots*sqrt(ptlim)
         pt0(j)= pt0_l(j)/roots
         pt1(j)= min(pt1_l(j),ptlim)/roots
*
c-      Protection:
*
        if (pt0(j).gt.pt1(j)) goto 100
        xmr(j)= xm(j)/roots
      enddo
*
c-    momenta:
*
      en    = 0.5d0*(x1+x2)      ! Rescaled initial total energy
      pz    = 0.5d0*(x1-x2)      ! Rescaled initial longitudinal momentum
*
      if (md.eq.1) then
        call mom(0,nw,pt0,pt1,eta0,xmr,en,pz,pm,wgt1,ranram,lw)
        if (lw.eq.0) then
          call momr(1,nw,xmr,en,pz,pm,wgt2,lw)
          if (lw.eq.0) then
            if (wgt1.eq.0.d0.and.wgt2.eq.0.d0) goto 100
            wgt= wgt1*wgt2/
     +          (wgt2*apw(1)*djbintot+wgt1*apw(2)*djb1)
          else
            wgt= wgt1/(apw(1)*djbintot)
          endif
        else
          goto 100
        endif
      else
        call momr(0,nw,xmr,en,pz,pm,wgt2,lw)
        if (lw.eq.0) then
          call mom(1,nw,pt0,pt1,eta0,xmr,en,pz,pm,wgt1,ranram,lw)
          if (lw.eq.0) then
            m= 0
            do l= lmin,lmax
              m= m+1
              if (mask(l).eq.0) then
                nbin= min(nx1(l),int(1.d0+dfloat(nx1(l))*ranram(m)))
                call onedimbin(3,nbin,djbin,l,ndummy,dummy)
                djb0= djb0*djbin
              endif
            enddo
            if (wgt1.eq.0.d0.and.wgt2.eq.0.d0) goto 100
            wgt= wgt1*wgt2/
     +          (wgt1*apw(2)*djbintot+wgt2*apw(1)*djb0*djb1)
          else
            wgt= wgt2/(apw(2)*djbintot)
          endif
        else
          goto 100
        endif
      endif
*
c-    from the x1,x2 integration
*
      wgt= wgt*wtau*wjr
*
c-    Rescaling of momenta and weight:
*
      do j= 1,nw
        do k= 0,3
          if (j.ne.nw) then
            p(k,j+4)= pm(k,j)*roots
          else
            pm(k,j) = pm(k,j)*roots
          endif
        enddo
      enddo
*
      wgt= wgt*s**(nw-2)     
*
      ran(1)= ranram(2*nw)
      call rans(ran0)
      ran(2)= ran0
      call dec2fm(0,sw,pm(0,nw),
     +            qm2,qm2,p(0,3),p(0,4),dj2,ran)
*
      wgt= wgt/dj1/dj2
*
      return
 100  lw= 1
      wgt= 0.d0
      return
      end
*
c-------------------------------------------------------------------
      subroutine alsgrd
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include '2Q.inc'
      integer nct1
      integer nvar,nch,n,nv
      integer init,j,k
      real*8 v1,ni
      real*8 al1,bet1
      integer nct,nx1,nx2
      integer maxn,ncmax
      parameter (maxn= 100,ncmax= 1000)  
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      common/rdandwrt/al1(ncmax,maxn),nct1(maxn)
      common/ausil/init(ncmax),bet1(0:ncmax,maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/psopt/nvar,nv
      common/book/v1(ncmax),ni(maxn)
      data v1/ncmax*0d0/,ni/maxn*0/
      integer mask,mmask 
      real*8 peropt
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      data mask/maxn*0/,mmask/maxn*1/,peropt/maxn*0.9d0/
* 
c-    initialise size of grids:
*
c-    for MC over jproc:

      nct(1)   = 0              ! first bin of variable 1
      nx1(1)   = jprocmax       ! number of jprocs
*
c-    for the Phase-space reweighting:
*
      if(ihvy.le.5) then
         nv  = 2*(njets+2)-1
      elseif(ihvy.eq.6) then
         nv  = 2*(njets+2)-1
      else
         print*,'wrong value of ihvy in setgrid'
         stop
      endif
      nvar= nv*ngrid
      do n= 2,nvar+1
        nch= 10
*
        nct(n)   =  nct(n-1)+nx1(n-1)   
        nx1(n)   =  nch 
        if (ihvy.le.5) then
*
c-        set the channels shared by mom and momr: 
*  
          if (mod(n-1,nv).eq.1)    mask(n)= 1
          if (mod(n-1,nv).eq.2)    mask(n)= 1
        else
*
c-        set the channels shared by mom and momr: 
*  
          if (mod(n-1,nv).eq.1)    mask(n)= 1
          if (mod(n-1,nv).eq.2)    mask(n)= 1
          if (mod(n-1,nv).eq.nv-1) mask(n)= 1
          if (mod(n-1,nv).eq.0)    mask(n)= 1
        endif
*
c-      set the percentages of optimization: 
*
        if (mask(n).eq.1) peropt(n)= 0.9d0 
        if (mask(n).eq.0) peropt(n)= 1.d0-dfloat(nv-1)/200.d0
      enddo
c-
      do n= 1,nvar+1
        nct1(n)= nx1(n)
        do j= 1,nct1(n)                                         

          init(nct(n)+j)= 0                                     

c-        the weights are initialized to be all equal 

          al1(j,n)= 1.d0/nct1(n)                                  
        enddo

        do j= 1,nct1(n)                                         
          bet1(j,n)= 0                                         
          do k= 1,j                                          
            bet1(j,n)= bet1(j,n)+al1(k,n)                          
          enddo                                              
        enddo                                                
      enddo
*
c-    protection:
*
      if (nct(nvar+1)+nct1(nvar+1).gt.ncmax) then
        print*,'INCREASE NCMAX'
        stop
      endif
*
      if ((nvar+1).gt.maxn) then
        print*,'INCREASE MAXN'
        stop
      endif
      end
*
c-------------------------------------------------------------------
      subroutine dumpwgt
c-------------------------------------------------------------------
c     routine for the debugging of individual large-weight events
c     Not functional to the rnning of the code
      implicit none
      include 'alpgen.inc'
      include '2Q.inc'
      real *8 tmpmtt,ptmp(4)
      integer i
      write(niosta,*) 'x1,x2=',x1,x2
      write(niosta,*) 'p(Q1)=',(p(i,3),i=1,3)
      write(niosta,*) 'p(Q2)=',(p(i,4),i=1,3)
      write(niosta,*) 'p(j1)=',(p(i,5),i=1,3)
      write(niosta,*) 'p(j2)=',(p(i,6),i=1,3)
      write(niosta,*) 'p(j3)=',(p(i,7),i=1,3)
      write(niosta,*) 'p(j4)=',(p(i,8),i=1,3)
c
c b-bbar inv mass
      ptmp(4)=p(4,3)+p(4,4)
      tmpmtt=ptmp(4)**2
      do i=1,3
         ptmp(i)=p(i,3)+p(i,4)
         tmpmtt=tmpmtt-ptmp(i)**2
      enddo
      tmpmtt=sqrt(tmpmtt)
      write(niosta,*) 'm(b-bbar)=',tmpmtt
c
c j-j inv mass
      ptmp(4)=p(4,5)+p(4,6)
      tmpmtt=ptmp(4)**2
      do i=1,3
         ptmp(i)=p(i,5)+p(i,6)
         tmpmtt=tmpmtt-ptmp(i)**2
      enddo
      tmpmtt=sqrt(tmpmtt)
      write(niosta,*) 'm(jj)=',tmpmtt
      end
c
c***** alpha-specific routines
c***** alpha-specific routines
C***********************************************************************
       subroutine matrix(flvmlm,posi,impul,hel,rst,labcol
     >                    ,wgcol,colfl)
C***********************************************************************
C
C Input parameters: FLVMLM(NMAX) particles according
C                   to michelangelo convention, POSI(NMAX) positions of
C                   incoming particles IMPUL(4,8) particles four momenta,
C                   HEL, array of particle helicities 
C                   LABCOL to choose whether su(3) mode only (LABCOL=0) or
C                   color flow is required (LABCOL=1) or dual mode only
C                   (LABCOL=2), WGCOL random number for color flow unweighting
C Output parameters:   RESULT squared matrix element (modulus), COLFLOW 
C                       string representing the selected coulor flow 
C
C
C At present:
C
C     the color flow is returned in the following way: each particle is
C     represented by a pair of integers (n1,n2); two particle (n1,n2), (n3,n4)
C     are colour connected if either n2=n3 or n4=n1 (in this case the coulor
C     flows from 2 to 1). COLFLOW will contain a string of pairs of integers 
C     ordered as follows: d (or dbar) nu_ebar bbar ubar (or u) e^- b glu glu
C     (with gluons ordered according to the ordering of momenta).
C     (COLFLOW=(n1,....,nn,....) only the first nn=2*particles_number
C       elements to be used)   
C
C
C   IMPUL(J,NPART) ===   J=1 Energy of the NPART-th particle
C   IMPUL(J,NPART) ===   J=2,3,4 x,y,z components  of the NPART-th particle 
C                                                         three momentum
C    the order of the momenta is: 
C        dbar nubar bbar u e- b ( glu glu glu)  (processes 1,2,3,4)
C        dbar nubar cbar bbar u e- c b (glu)  (processes > 5 )
C    where all particles are assumed outcoming and incoming gluons
C    must be before outcoming ones
C
         implicit none
C
         integer nmax        !maximum number of external particles, 
         parameter (nmax=10)
         integer nlb
         parameter (nlb=4)    
         real*8 impul(4,nmax),wgcol    !the contribution to the integral
         real*8 rst
c         complex*16 rst
         integer labcol,colfl(2*nmax),flvmlm(nmax)
c
         integer posi(2)
         integer hel(nmax)             !particles helicities
C
         integer flgdual          !dual (0) or su3 (1) amplitudes
         common/dual/flgdual
         integer color(2*nmax)    !coulor string
         integer colst(2*nmax)
         common/colore/color
         integer ncls,ant4(4,2),ant6a(6,6),ant6b(6,6),ant2(2,1)
         integer nant,nantl,j3
         parameter (ncls=2)       !different class of processes
         integer nx,proc(nmax),antns(6),ndual,j1,j2
         parameter (ndual=40400)
         real*8 dualamp(ndual),colfac,damp,dampref,avgspin
         integer colaux(ndual,2*nmax),ndl,ndla(0:10,1:3),class
         integer rep(nmax),nqrk,nprt,nglu,nlep,ngb,nphot,nh
         common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot,nh
         integer fq(nmax),pq(nmax)
         data ndla/0,0,1,5,23,119,719,5039,40319,0,0,
     >          0,1,5,23,119,719,5039,0,0,0,0,
     >          0,2,11,59,359,6*0/

         save ndla,ant4,ant6a,ant2
         complex*16 result
         integer inpint(1000)
         common/initinter/inpint  
C
         data ant2  /1,2/
         data ant4  /1,4,2,3,
     >               1,3,2,4/
         data ant6a  /1,5,2,6,3,4,
     >               1,6,2,4,3,5,
     >               1,6,2,5,3,4,
     >               1,4,2,6,3,5,
     >               1,4,2,5,3,6,
     >               1,5,2,4,3,6/         
C
         data inpint/ 6,
***     >               11,1,1,1,1,  11,1,2,1,2,  11,3,2,3,2,   2,1,1,1,2,  ! guu, gdd, gbb, w+ud
     >           11,1,1,1,1,    11,3,1,3,1,  11,1,2,1,2,  11,3,2,3,2,    ! guu, gtt, gdd, gbb
     >           11,2,1,2,1,  11,2,2,2,2,                                !gcc,gss
***     >                3,1,4,1,3,  11,2,1,2,1,                            ! w-en, gcc
     >                0,                                                 ! N of yukawa
     >                2,                                                 ! N self-gauge
     >                11,11,11,  12,11,11,                               ! ggg Auxgg
     >                961*-100/
***     >                971*-100/
c
         avgspin=0.d0
C
         if (labcol.eq.0) then
C
          flgdual=1
          call matrix0(flvmlm,posi,impul,hel,result)
c
         elseif (labcol.eq.1) then
c
c     
          j1=0
          do j3=1,nmax
           if (abs(flvmlm(j3)).gt.0.and.abs(flvmlm(j3)).le.6) then
            j1=j1+1
            fq(j1)=flvmlm(j3)
            pq(j1)=j3
           endif
          enddo
          if(j1.ne.nqrk) then
           write(*,*)'something wrong in matrix, NQRK not ok',nqrk,j3
           stop
          endif
c
          flgdual=0
          do j1=1,2*nmax
           colst(j1)=color(j1)
          enddo 
c
          do j1=1,nmax
           proc(j1)=abs(rep(j1))
          enddo
c
          ndl=0
          if(nqrk.eq.2) then
           nantl=1
          elseif(nqrk.eq.4) then
           nantl=2
          elseif(nqrk.eq.6) then
           nantl=6
          else
           write(*,*)'Cannot deal with',nqrk,'quarks'
           stop
          endif
          do j3=1,nantl
           do j1=1,nqrk
            if(nqrk.eq.2) then
             antns(j1)=pq(ant2(j1,j3))
            elseif(nqrk.eq.4) then
             antns(j1)=pq(ant4(j1,j3))
            elseif(nqrk.eq.6) then
             antns(j1)=pq(ant6a(j1,j3))
            endif
           enddo
           j2=nqrk/2
           do j1=0,ndla(nglu,j2)
            nx=j1
            if(nqrk.eq.2) then
             call gendual2q(proc,antns,nglu,colst,nx,colfac) 
            elseif(nqrk.eq.4) then
             call gendual4q(proc,antns,nglu,colst,nx,colfac) 
            elseif(nqrk.eq.6) then
             call gendual6q(proc,antns,nglu,colst,nx,colfac) 
            endif 
            if(abs(colfac).gt.1.d-20) then
             ndl=ndl+1
             call matrix0(flvmlm,posi,impul,hel,result)
             dualamp(ndl)=abs(result*colfac)**2
             do j2=1,2*nmax
              colaux(ndl,j2)=color(j2)
             enddo
            endif
           enddo
          enddo
c
          damp=0
          do j1=1,ndl
           damp=damp+dualamp(j1)
          enddo
          dampref=damp*wgcol
          j1=0
          damp=0.
          do while(damp.lt.dampref.and.j1.lt.ndl)
           j1=j1+1
           damp=damp+dualamp(j1)
          enddo
          ndl=j1
          if(ndl.eq.0) then
             do j1=1,2*nmax
                colfl(j1)=colst(j1)
             enddo
          else 
             do j1=1,2*nmax
                colfl(j1)=colaux(ndl,j1)
             enddo
          endif
C
         else
          write(6,*)'wrong LABCOL in matrix',labcol
          stop
         endif
C
         rst=abs(result)**2
c         rst=result
C
         return
         end
c
c
c***** alpha-specific routines
*
      subroutine mom(lflag,np,pt0,pt1,eta0,xmr,en0,pz0,
     .               pm,wgt,ranram,lw)
      implicit none
      real*8 pi,wgt,rpr,phr,ran0,wj,en0,pz0,prx,pry
      real*8 etai1,etai2,etaj1,etaj2,xmri2,bxmri2,wjeta
      real*8 cn,wjaci,wjacj,phi,ranj,ref,alimp,alimm
      real*8 p0m,p0p 
      real*8 phj,alim,phip,phim,etap,etam
      integer npar,n,nri,nr,iter,i,j,k,lw,icont,lflag
      parameter (pi= 3.14159265358979323846264338327950d0)
      parameter (npar= 20)
      real*8 pt0(npar),pt1(npar),eta0(npar),
     .       xmr(npar),bxmr(npar)
      real*8 pm(0:4,npar),pt(npar),eta(npar)
      integer jp(npar),init,np
      real*8 p0r(npar),p1r(npar),qcut(npar),pcut(npar),bpt0(npar)
      real*8 bet0(npar),ph(npar)
      real*8 ranram(111)
      real*8 en,pz,ga,de,p,q,a,b,cut,ausp,ausm,v,vmr,vpr,v2
      real*8 det,rx,x1m,x2m,pmod,qmod,den
      real*8 argsinh,aus,betp,alpp
      real*8 etalim,csi
      data init/0/
      save
*
      if (init.eq.0) then
        init= 1
        do iter= 1,np
          jp(iter)= iter
        enddo
        n  = np-1
        nri= 3
        etalim= 100.d0
*
        do iter= 1,n
          bpt0(jp(iter))= 0.d0
          pcut(jp(iter))= 0.d0
          qcut(jp(iter))= 0.d0
          bet0(jp(iter))= etalim
          do k= iter+1,np
              bpt0(jp(iter))= bpt0(jp(iter))+pt0(jp(k))
              pcut(jp(iter))= max(pcut(jp(iter)),pt0(jp(k)))
              bet0(jp(iter))= min(bet0(jp(iter)),eta0(jp(k)))
          enddo
        enddo
        qcut(jp(n))= pt0(jp(np))
*
      endif
*
      do iter= 1,n
        bxmr(jp(iter))= 0.d0
        do k= iter+1,np
            bxmr(jp(iter))= bxmr(jp(iter))+xmr(jp(k))
        enddo
      enddo
*
      lw  = 0
      nr  = nri
      rpr = 0.d0
      wgt = 1.d0
      en  = en0
      pz  = pz0
      if (lflag.eq.0) then
        phr = 0.d0
      elseif (lflag.eq.1) then
        prx = 0.d0
        pry = 0.d0
      else
        goto 101
      endif
*
      do iter= 1,n
        i= jp(iter)
        j= jp(iter+1)
*
c-      Generation of transverse momenta when lflag= 0
*
        ga = en+pz
        de = en-pz 
*
        if (ga.lt.0.d0) goto 100
        if (de.lt.0.d0) goto 100
        v  = sqrt(ga*de)
*
c-      Protection:
*
        vmr= v-rpr
        if (vmr.lt.0.d0) then
          vmr= dabs(vmr)
          v  = vmr+2.d0*rpr
          de = ga/v/v
        endif
        vpr= vmr+2.d0*rpr
        v2 = v*v
*
        if (ga.ge.de) then
          alpp= 2.d0*ga*exp(-eta0(i))
        else
          alpp= 2.d0*de*exp(-eta0(i))
        endif
        betp= ((ga+de)*cosh(eta0(i))-dabs(ga-de)*sinh(eta0(i)))
        csi = min(alpp,betp)
*
        if (lflag.eq.1) then
          pt(i)= sqrt(pm(1,i)**2+pm(2,i)**2)
          if (iter.ne.1) then
            phi= ((pm(1,i)*prx+pm(2,i)*pry)/pt(i)/rpr)
            if (phi.gt. 1.d0) phi= 1.d0
            if (phi.lt.-1.d0) phi=-1.d0
            phi= acos(phi)
          endif
          prx  = prx-pm(1,i)
          pry  = pry-pm(2,i)
          pt(j)= sqrt(prx*prx+pry*pry)
*
          eta(i)= -log(tan(0.5d0*atan2(pt(i),pm(3,i))))
        endif
*
        xmri2 = xmr(i)**2
        bxmri2= bxmr(i)**2
*
        aus   =  vpr*vmr 
        ausp  = (aus+xmri2-bxmri2)
        det   = (ausp**2-4.d0*xmri2*aus)
        if (det.lt.0.d0) then
          det= abs(det)
          if (det.gt.1.d-12) goto 100
        endif
        det   = sqrt(det)
        if (ausp.gt.0.d0) then
           p0p= (rpr*ausp+v*det)/2.d0/aus
           p0m=-(ausp**2-4.d0*v2*xmri2)/4.d0/aus/p0p
        else
           p0m= (rpr*ausp-v*det)/2.d0/aus
           p0p=-(ausp**2-4.d0*v2*xmri2)/4.d0/aus/p0m
        endif
        p0r(i)= max(pt0(i),p0m,qcut(i)-rpr)
        p1r(i)= min(p0p,v-bpt0(i))
        if (csi.gt.2.d0*v) then
           p1r(i)= min(p1r(i),(vpr*vmr)/(csi-2.d0*rpr))
        endif
        cn= 1.d0
        if (p0r(i).eq.0.d0) cn= 0.98d0
        if (p0r(i).gt.p1r(i)) goto 100
        call ppeaka(lflag,max(pt0(i),xmr(i),1.d-2),
     .              cn,p0r(i),p1r(i),pt(i),wjaci,ranram(nr),lw)
        nr= nr+1
        if (lw.eq.1) goto 100
*
        if (lflag.eq.0) call rans(ran0)
*
        if (iter.eq.1) then
          if (lflag.eq.0) then
            phi= 2.d0*pi*ran0
            pt(j)= pt(i)
            phj= phi+pi
          endif
          wgt= wgt*2.d0*pi*pt(i)*wjaci
        else
          if (lflag.eq.0) then
            if (ran0.le.0.5d0) then
              ref = 1.d0
            else
              ref =-1.d0
            endif
            call rans(ranj)
          endif
*
          p0r(j)= max(dabs(pt(i)-rpr),qcut(i),pt(i)-v+2.d0*pcut(i))
          det= (v-sqrt(pt(i)**2+xmri2))**2-bxmri2
          if (det.lt.0.d0) then
             det= abs(det)
             if (det.gt.1.d-12) goto 100
          endif
          det   = sqrt(det)
          p1r(j)= min(pt(i)+rpr,det)
          if (csi.gt.2.d0*v) then
            p1r(j)= min(p1r(j),sqrt(pt(i)*(pt(i)-csi)+v2))
          endif
          if (p0r(j).gt.p1r(j)) goto 100
*
          alimp= (rpr**2+pt(i)**2-p1r(j)**2)/2.d0/rpr/pt(i)
          alimm= (rpr**2+pt(i)**2-p0r(j)**2)/2.d0/rpr/pt(i)
*
c-        Protections:
*
          if (alimp.gt. 1.d0) alimp= 1.d0
          if (alimp.lt.-1.d0) alimp=-1.d0 
          if (alimm.gt. 1.d0) alimm= 1.d0
          if (alimm.lt.-1.d0) alimm=-1.d0
*          
          phip= acos(alimp)
          phim= acos(alimm)
          if (phim.gt.phip) goto 100
*
          call fflat(lflag,phi,phim,phip,wjacj,ranj,lw)
          if (lw.eq.1) goto 100
          wgt= wgt*wjaci*wjacj*pt(i)*2.d0
*
          if (lflag.eq.0) then
            phi= phi*ref
            pt(j)= dabs(rpr**2+pt(i)**2-2.d0*rpr*pt(i)*cos(phi))
            pt(j)= sqrt(pt(j))
*
c-          Protections:
*
            if (rpr.eq.0.d0)   rpr  = 1.d-25
            if (pt(j).eq.0.d0) pt(j)= 1.d-25
            alim= (rpr**2+pt(j)**2-pt(i)**2)/2.d0/rpr/pt(j)
*
c-          Protections:
*
            if (alim.gt. 1.d0) alim= 1.d0
            if (alim.lt.-1.d0) alim=-1.d0
*
            phj=-acos(alim)*ref
*
          endif
        endif  
*
        if (lflag.eq.0) then
          phi= phi+phr
          phj= phj+phr
          pm(1,i)= pt(i)*sin(phi)
          pm(2,i)= pt(i)*cos(phi)
          pm(1,j)= pt(j)*sin(phj)
          pm(2,j)= pt(j)*cos(phj)
          ph(i)  = phi
          if (iter.eq.n) ph(j)  = phj
*
          phr= phj
        else
          ph(i)= acos(pm(2,i)/pt(i)) 
          if (pm(1,i).lt.0.d0) ph(i)= 2.d0*pi-ph(i)
          if (iter.eq.n) then
            ph(j)= acos(pm(2,j)/pt(j)) 
            if (pm(1,j).lt.0.d0) ph(j)= 2.d0*pi-ph(j)
          endif
        endif
*
        rpr= pt(j)
*
        p  = pt(i)
        q  = pt(j)
        a  = 1.d0+xmri2/p/p
        b  = max(1.d0+bxmri2/q/q,1.d0+4.d0/q/q*pcut(i)*
     .          (pcut(i)-q),(bpt0(i)/q)**2)
*
        cut= (sqrt(a)*p+sqrt(b)*q)**2-v2
        if (cut.gt.1.d-12) goto 100
        ausm= a*p*p-b*q*q
        det = (v2+ausm)**2-4.d0*v2*a*p*p
        if (det.lt.0.d0) then
           det= abs(det)
           if (det.gt.1.d-12) goto 100
        endif
        det = sqrt(det)
        rx  = v2+ausm
        x1m= (rx+det)/2.d0/p/de  
        x2m= ga*a/de/x1m                  ! x2m < x1m
        etai1= argsinh(0.5d0*(x1m-a/x1m))
        etai2= argsinh(0.5d0*(x2m-a/x2m))
        if (iter.le.n-1) then
*
c-        Generation of etas when lflag= 0
*
          etam= max(etai2,-eta0(i))        
          etap= min(etai1, eta0(i))        
          if (etam.ge.etap) goto 100
*
comment 
c          call ttriangle(lflag,etai2,etai1,etam,etap,eta(i),
c     .                  wjeta,ranram(nr),lw)
          if (iter.ge.2) then
            k    = jp(iter-1)
            alimp= 4.d0*(xmr(k)+xmr(i))**2   
            call spk(lflag,alimp,0.98d0,pt(k),pt(i),ph(k),ph(i),
     .               xmr(k),xmr(i),eta(k),
     .               etam,etap,eta(i),wjeta,ranram(nr),lw)
c            call etagen(lflag,ph(k),ph(i),eta(k),
c     .                  etam,etap,eta(i),wjeta,ranram(nr),lw)
          else
            call ttriangle(lflag,etai2,etai1,etam,etap,eta(i),
     .                     wjeta,ranram(nr),lw)
          endif
comment
          nr= nr+1
          if (lw.eq.1) goto 100
          wgt   = wgt*wjeta
          if (lflag.eq.0) then
            pm(0,i)= sqrt(pt(i)**2*cosh(eta(i))**2+xmri2)   
            pm(3,i)= pt(i)*sinh(eta(i))                         
            pm(0,j)= en-pm(0,i)
            pm(3,j)= pz-pm(3,i)
          endif
          en     = en-pm(0,i)
          pz     = pz-pm(3,i)
          wgt    = wgt*sqrt(1.d0-(xmr(i)/pm(0,i))**2)
        else
*
c-        Cuts on the last 2 etas
*
          etaj1= argsinh((pz-pt(i)*sinh(etai1))/q)
          etaj2= argsinh((pz-pt(i)*sinh(etai2))/q)
*
c-        Look whether 0,1 or 2 solutions may contribute
*
          icont= 0                        
          if ((dabs(etai1).lt.eta0(i)).and.
     .        (dabs(etaj1).lt.eta0(j))) icont= icont+1
          if ((dabs(etai2).lt.eta0(i)).and.
     .        (dabs(etaj2).lt.eta0(j))) icont= icont+2
*
          if (icont.eq.0) then
             goto 100
          elseif (icont.eq.1) then
             if (lflag.eq.0) then
               eta(i)= etai1
               eta(j)= etaj1
             endif
             wj= 1.d0
          elseif (icont.eq.2) then
             if (lflag.eq.0) then
               eta(i)= etai2
               eta(j)= etaj2
             endif
             wj= 1.d0
          else
             if (lflag.eq.0) then
               call rans(ran0)
               if (ran0.lt.0.5d0) then
                 eta(i)= etai1
                 eta(j)= etaj1
               else
                 eta(i)= etai2
                 eta(j)= etaj2
               endif
             endif
             wj= 2.d0
          endif
          if (lflag.eq.0) then
             pm(0,i)= sqrt(pt(i)**2*cosh(eta(i))**2+xmri2)   
             pm(3,i)= pt(i)*sinh(eta(i))                         
             pm(0,j)= en-pm(0,i)
             pm(3,j)= pz-pm(3,i)                                  
          elseif (lflag.eq.1) then
             en     = en-pm(0,i)-pm(0,j)
             pz     = pz-pm(3,i)-pm(3,j)
             eta(j) =-log(tan(0.5d0*atan2(pt(j),pm(3,j))))
          else
             goto 101
          endif
          pmod= pm(0,i)**2-xmri2
          qmod= pm(0,j)**2-xmr(j)*xmr(j)
          if (pmod.lt.0.d0) goto 100    
          if (qmod.lt.0.d0) goto 100
          pmod= sqrt(pmod)
          qmod= sqrt(qmod)
          den = dabs(pm(0,j)*pm(3,i)-pm(0,i)*pm(3,j))
          wgt = wgt*wj*
     .          pmod*qmod/p/q/cosh(eta(j))/cosh(eta(i))/den
        endif
      enddo     
      wgt= wgt/2.d0**np      
*
      return
 100  wgt= 0.d0
      lw= 1
      return
 101  print*,'ERROR IN SUBROUTINE MOM'
      stop
      end
*
      subroutine momr(lflag,np,xmr,en0,pz0,
     .                pm,wgt,lw)
      implicit none
      integer lflag,np,npar,lw,j
      parameter (npar= 100)
      real*8 xmr(npar),en0,pz0,pm(0:4,npar),wgt,dj
      real*8 x1,x2,bvel,gvel,rootshr,pr(4,npar)
*
      lw= 0
      if (lflag.eq.1) then
        en0= 0.d0
        pz0= 0.d0
        do j= 1,np
          en0= en0+pm(0,j)          
          pz0= pz0+pm(3,j)
        enddo
      endif
      x1  = dabs(en0+pz0)
      x2  = dabs(en0-pz0)
      bvel= (x1-x2)/(x1+x2)
      gvel= 1.d0/sqrt(1.d0-bvel*bvel)
      rootshr= sqrt(x1*x2)
*
      if (lflag.eq.1) then
        do j= 1,np
          pr(4,j)= gvel*(pm(0,j)-bvel*pm(3,j))
          pr(1,j)= pm(1,j)
          pr(2,j)= pm(2,j)
          pr(3,j)= gvel*(pm(3,j)-bvel*pm(0,j))
        enddo 
      endif
      call rambo(lflag,np,rootshr,xmr,pr,dj)
      wgt= 1.d0/dj
      if (lflag.eq.1) return
*
      do j= 1,np
        pm(0,j)= gvel*(pr(4,j)+bvel*pr(3,j))
        pm(1,j)= pr(1,j)
        pm(2,j)= pr(2,j)
        pm(3,j)= gvel*(pr(3,j)+bvel*pr(4,j))
      enddo
*
      return
      end

