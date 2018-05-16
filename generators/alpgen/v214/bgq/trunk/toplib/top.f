c-------------------------------------------------------------------
      subroutine alsprc
c     assigns the hard process code
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      ihrd=13
      end

c-------------------------------------------------------------------
      subroutine alhset
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'top.inc'
      real*8 m1aux,m2aux
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c     ngrid is the total number of grids allowd for P.S. variables.
c     jgrid(jproc) labels the grid associated to a given jproc.
      integer jgridp(100,4)
      data jgridp/1,2,3,4,96*0,
     +            1,2,3,4,96*0,
     +            1,2,1,2,3,4,5,6,3,4,5,6,88*0,
     +            1,1,2,3,2,3,4,5,6,7,4,5,6,7,1,1,84*0/
      integer jpmx(0:6,4) !jpmx(njets,process)
      data jpmx/-1,2,4,4,4,4,4,
     +          -1,2,4,4,4,4,4,
     +           4,12,12,12,12,12,12,
     +          -1,6,14,16,16,16,16/
      integer njmx(0:6,4) !njmx(njets,process)
      data njmx/-1,2,4,4,4,4,4,
     +          -1,2,4,4,4,4,4,
     +           2,6,6,6,6,6,6,
     +          -1,3,7,7,7,7,7/
      integer i
      character*2 Q(4:6)
      character*5 Qbar(4:6)
      data Q/' c',' b',' t'/,Qbar/' cbar',' bbar',' tbar'/
c
c
c     parameters for the gauge invariance prescription:
      winsize  = 2.d0/pi
      resonance= 'n'
      wmode    = 'nn'
c     process input parameters
      if(itopprc.ne.3) njets=njets+1
      npart=njets+3
      if(itopprc.ge.3) npart=njets+4
      jprocmax=jpmx(njets,itopprc)
      if(njets.gt.1) jprocmax= -1
      if(jprocmax.lt.0) then
         write(*,*) 'jets number not allowed ',njets
         stop
      endif
      ngrid=njmx(njets,itopprc)
      do i=1,100
         jgrid(i)=jgridp(i,itopprc)
      enddo
      mq= mt
      mq2= mb
      idecay='y'
c masses
      do i=1,njets+4
         p(5,i)=0
      enddo
      p(5,3)=mt
c run parameters:
      if(itopprc.eq.1) then
         write(niosta,*) 'top + ',njets,' jets'
      elseif(itopprc.eq.2) then
         write(niosta,*) 'top + bbar + ',njets-1,' jets'
      elseif(itopprc.eq.3) then
         write(niosta,*) 'top + ',njets,' jets + W'
      elseif(itopprc.eq.4) then
         write(niosta,*) 'top + bbar + ',njets-1,' jets + W'
      endif
      write(niosta,*)
     $     '======================================='
      write(niosta,*) 'Heavy quark masses: mt and mb=',mt,mb
      write(niosta,*)
     $     'Generation cuts for the partonic event sample:' 
      write(niosta,*) 'No generation cuts on t'
      if(itopprc.eq.1) then
         write(niosta,*) '     Light jets:'
         write(niosta,*) 'ptmin=',ptjmin,' |etamax|=',etajmax
     $        ,' dR(j-j),dR(b-j)>',drjmin 
      elseif(itopprc.eq.2) then
         write(niosta,*) ' b jet: '
         write(niosta,*) 'ptbmin=',ptbmin,' |etabmax|=',etabmax
         if(njets-1.gt.0) then
            write(niosta,*) '     Light jets:'
            write(niosta,*) 'ptmin=',ptjmin,' |etamax|=',etajmax
     $           ,' dR(j-j),dR(b-j)>',drjmin 
         endif
      elseif(itopprc.eq.3) then
         write(niosta,*) '     Light jets:'
         write(niosta,*) 'ptmin=',ptjmin,' |etamax|=',etajmax
     $        ,' dR(j-j),dR(b-j)>',drjmin 
      elseif(itopprc.eq.4) then
         write(niosta,*) ' b jet: '
         write(niosta,*) 'ptbmin=',ptbmin,' |etabmax|=',etabmax
         if(njets-1.gt.0) then
            write(niosta,*) '     Light jets:'
            write(niosta,*) 'ptmin=',ptjmin,' |etamax|=',etajmax
     $           ,' dR(j-j),dR(b-j)>',drjmin 
         endif
      endif
      end
c-------------------------------------------------------------------
      subroutine selflav(jproc,xlum,afl)
      implicit none
      include 'alpgen.inc'
      include 'top.inc'
      integer afl(maxpar)
      integer jproc
      double precision xlum
c
c   for top decay
c
      complex*16 vtp(4),utpb(4) !array containing the individual field 
      real*8 ptp(4),ptpb(4)
      common/tspinors/vtp,utpb
      real*8 pf(4),pfb(4),pb(4)
      common/pdec/pf,pfb,pb
      integer exch(4),i,j,k
      data exch/2,3,4,1/
c
c   for W decay
c
      complex*16 eps(4,maxpar),aux(4)
      common/vdec/eps
      real*8 qw(4),qf(4),qfb(4)
c
      if(itopprc.eq.1) call selflav1(jproc,xlum,afl)
      if(itopprc.eq.2) call selflav2(jproc,xlum,afl)
      if(itopprc.eq.3) call selflav3(jproc,xlum,afl)
      if(itopprc.eq.4) call selflav4(jproc,xlum,afl)
c
c     decaying top quark
c
      if(idecay.eq.'y') then 
         do i=1,4
            do j=1,3
               do k=1,2
                  idec(i,j,k)=0.d0
               enddo
            enddo
         enddo

         if(ifl(3).eq.-6) then
            ptpb(1)= p(4,3)
            do i=1,3
               ptpb(i+1)= p(i,3)
            enddo
            call vtop(ptpb,mq,mq2,mw,vtp)
            do i=1,4
               idec(i,1,1)=pb(exch(i))
               idec(i,2,1)=pf(exch(i))
               idec(i,3,1)=pfb(exch(i))
            enddo
         elseif(ifl(3).eq.6) then
            ptp(1)= p(4,3)
            do i=1,3
               ptp(i+1)= p(i,3)
            enddo
            call utopb(ptp,mq,mq2,mw,utpb)
            do i=1,4
               idec(i,1,1)=pb(exch(i))
               idec(i,2,1)=pf(exch(i))
               idec(i,3,1)=pfb(exch(i))
            enddo
         else
            write(*,*) 'something wrong in top assignment in selflav'
            stop
         endif
*
*  W decay
*
         if(itopprc.ge.3) then
            qw(1)=p(4,njets+4)
            do j=2,4
               qw(j)=p(j-1,njets+4)
            enddo
            call vpol(mw,qw,-1.d0,1.d0,qf,qfb,aux)
            do i=1,4
               wdec(i,1)=qf(exch(i))
               wdec(i,2)=qfb(exch(i))
            enddo
            do j=1,4
               eps(j,1)=aux(j)
            enddo
         endif
      endif
      call usrfll(3)
      return
      end
c
      subroutine usrfll(it)
      implicit none
      include 'alpgen.inc'
      include 'top.inc'
      integer i,j,it,ib,ibb
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
      do i=1,4
         ptop(i)=p(i,it)
      enddo
      if(itopprc.eq.1.or.itopprc.eq.3) then
         do i=1,4
            do j=1,njets
               pjet(i,j)=p(i,j+3)
               ptj(j)=pt(j+3)
               etaj(j)=eta(j+3)
            enddo
         enddo
         do j=1,njets
            do i=j+1,njets
               drjj(i,j)=dr(i+3,j+3)
               drjj(j,i)=drjj(i,j)
            enddo
         enddo
      elseif(itopprc.eq.2.or.itopprc.eq.4) then
         do i=1,4
            do j=1,njets-1
               pjet(i,j)=p(i,j+3)
               ptj(j)=pt(j+3)
               etaj(j)=eta(j+3)
            enddo
         enddo
         do j=1,njets-1
            do i=j+1,njets-1
               drjj(i,j)=dr(i+3,j+3)
               drjj(j,i)=drjj(i,j)
            enddo
         enddo
      endif
      end
c-------------------------------------------------------------------
      subroutine selflav1(jproc,xlum,afl)
c     evaluates parton luminosities, and assigns PDG-code flavours
c     gluon=21 cbar=-4 sbar=-3 ubar=-2 dbar=-1 d=1 u=2 s=3 c=4
c     jproc
c
c   at present CKM assumed diagonal
c
c---t q final state (+ up to 4 gluons)
c
c  1  b q    -> t q (+ gluons)
c  2  q b    -> t q (+ gluons)
c
c---t + 2q final state (+ gluons)
c
c  3  g b    -> t q qbar' (+ gluons)
c  4  b g    -> t q qbar' (+ gluons)
c

c-----------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'top.inc'
c commons
      integer icconj
      common/hwconv/icconj
      integer afl(maxpar)
      real tmp(100),slum,cwgt,swgt,rn,tmptot
      real fcab(2)
      integer icab(-4:4,2) 
c     cabibbo partner: (*,1)=cabibbo allowed, (*,2)-cabibbo suppressed
c     where 1=d 2=u 3=s 4=c 0=g and negatives are antiparticles
      data icab/-3,-4,-1,-2,0,2,1,4,3, 
     +          -1,-2,-3,-4,0,4,3,2,1/
      integer iqpp(2,2,4)
c     iqpp(k,j,i):  if q=i, q'=cabibbo partner of i (allowed for j=1,
c     suppressed for j=2) then iqpp are the two (k=1,2) flavours .ne. (q.or.q')
      data iqpp/
c       du   dc   ud   us   sc   su   cs   cd   
     +  3,4, 2,3, 3,4, 1,4, 1,2, 1,4, 1,2, 2,3/
c
      double precision xlum,xrn 
      integer i,k,itmp,icount,init,jproc,ng,j,k1,j1
      real cfac(0:6),ifact(0:6)
      data cfac/8e0,6*3e0/,init/0/
      integer imap(-4:4)
      data imap/-2,-1,-2,-1,0,1,2,1,2/
      integer ibb,iflaux
      integer nlb
      parameter (nlb=13)
      integer i1,i2,nflv1,nflv2,i3,nf1,nf2
      integer nflpr41(20),nflpr51(20)
      integer labfl4(4,nlb),labfl5(5,nlb)
      character*2 lbflcr41(4,8),lbflcr51(5,4),lbc
c
c     overall efficiency for extraction of colour states:
c     >(non-zero color states) / 3**nq*8**ng
c     ccoef(i,j) for njets=j and i=#(light quark pairs)
      double precision effco
      double precision ccoef(3,0:8) !ccoeff(# of q-qb pairs,  # of gluons)
      data ccoef/0.333333333,0.185185185,0.127572016,
     +           0.1666666667,0.12037037,0.0936213992,
     +           0.114583333,0.0902777778,0.0743312757,
     +           0.087239583,0.072337963,0.06171232,
     +           0.070475260,0.0603841146,0.05278461,
     +           0.059122721,0.051834672,-1.d0,
     +           0.050923665,0.045412134,-1.d0,
     +           0.042060375,-1.d0,-1.d0,
     +           0.037214041,-1.d0,-1.d0/
      integer nlq  ! number of light quarks
      save init,cfac,ifact,nlq,iqpp,icab,fcab
      save labfl4,nflpr41,labfl5,nflpr51
c
      data lbflcr41 /'bq','uq','tq','dq', 'bq','cq','tq','sq', !jproc=1,2 
     +               'bq','db','tq','ub', 'bq','sb','tq','cb', !jproc=1,2
     +               'bb','ub','tb','db', 'bb','cb','tb','sb', !jproc=1,2
     +               'bb','dq','tb','uq', 'bb','sq','tb','cq'/ !jproc=1,2
      data lbflcr51/ 'gg','bq','tq','ub','dq', 'gg','bq','tq','cb','sq', !jproc=3,4 
     +               'gg','bb','tb','uq','db', 'gg','bb','tb','cq','sb'/ !jproc=3,4
      data nflpr41 /8,19*-100/
      data nflpr51 /4,19*-100/
c
      if(init.eq.0) then
         fcab(1)=real(ccab2)
         fcab(2)=real(scab2)
         njets=npart-3
         init=1
         ifact(0)=1e0
         do i=1,6
            ifact(i)=ifact(i-1)/real(i)
         enddo
c
c  initializing flavours for the various classes of processes
         do j1=1,2
            if(j1.eq.1) nflv1=nflpr41(1)
            if(j1.eq.2) nflv1=nflpr51(1)
c
            do j=1,nflv1
               do k=1,4+(j1-1)
                  if(j1.eq.1) lbc=lbflcr41(k,j)
                  if(j1.eq.2) lbc=lbflcr51(k,j)
                  if(lbc.eq.'gg') then
                     i3=0
                  elseif(lbc.eq.'dq') then
                     i3=1
                  elseif(lbc.eq.'uq') then
                     i3=2
                  elseif(lbc.eq.'sq') then
                     i3=3
                  elseif(lbc.eq.'cq') then
                     i3=4
                  elseif(lbc.eq.'bq') then
                     i3=5
                  elseif(lbc.eq.'tq') then
                     i3=6
                  elseif(lbc.eq.'db') then
                     i3=-1
                  elseif(lbc.eq.'ub') then
                     i3=-2
                  elseif(lbc.eq.'sb') then
                     i3=-3
                  elseif(lbc.eq.'cb') then
                     i3=-4
                  elseif(lbc.eq.'bb') then
                     i3=-5
                  elseif(lbc.eq.'tb') then
                     i3=-6
                  else
                     write(*,*)'something wrong in lbflcr4'
                     stop
                  endif
                  if(j1.eq.1) labfl4(k,j)=i3
                  if(j1.eq.2) labfl5(k,j)=i3
               enddo
            enddo 
         enddo
c
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
c  1  b q -> t q' + (njets-1) g        q=u,d,c,s    njets>0
c  2  q b -> t q' + (njets-1) g
c
      if(jproc.le.2) then
         nflv1=1
         nflv2=nflpr41(1)
         if(jproc.eq.1) then
            i1=1
            i2=2
         else
            i1=2
            i2=1
         endif
         i=0
         do j=nflv1,nflv2       !! FCAB????
            i=i+1
            nf1=labfl4(i1,j)
            nf2=labfl4(i2,j)
            tmp(i)=f1(nf1)*f2(nf2) ! *fcab(k)
            slum=slum+tmp(i)
         enddo
         rn=rn*slum
         i=0
         do j=nflv1,nflv2
            i=i+1
            tmptot=tmptot+tmp(i)
            if(tmptot.ge.rn) then 
               ifl(1)=labfl4(i1,j)
               ifl(2)=labfl4(i2,j)
               ifl(3)=labfl4(3,j)
               ifl(4)=labfl4(4,j)
               do k1=5,njets+3
                  ifl(k1)=0
               enddo
               effco = ccoef(2,njets-1)
               goto 100
            endif
         enddo
c
c  3  g b -> t qbar q' + (njets-2) g     q=u,d,c,s    njets>1
c  4  b g -> t qbar q' + (njets-2) g
c
      elseif(jproc.ge.3.and.jproc.le.4) then
         nflv1=1
         nflv2=nflpr51(1)
         if(jproc.eq.3) then
           i1=1
           i2=2
         else
           i1=2
           i2=1
         endif
         i=0
         do j=nflv1,nflv2                       !! FCAB????
           i=i+1
           nf1=labfl5(i1,j)
           nf2=labfl5(i2,j)
           tmp(i)=f1(nf1)*f2(nf2)            ! *fcab(k)
           slum=slum+tmp(i)
         enddo
         rn=rn*slum
         i=0
         do j=nflv1,nflv2
           i=i+1
           tmptot=tmptot+tmp(i)
           if(tmptot.ge.rn) then 
             ifl(1)=labfl5(i1,j)
             ifl(2)=labfl5(i2,j)
             ifl(3)=labfl5(3,j)
             ifl(4)=labfl5(4,j)
             ifl(5)=labfl5(5,j)
             do k1=6,njets+3
                ifl(k1)=0
             enddo
c             itmp=1
c             ibl=ifl(i1)
             effco = ccoef(2,njets-1)
             goto 100
           endif
         enddo
c
      else
         write(*,*) 'jproc not defined, slum=0'
         slum=0
         stop
      endif
c
      xlum=-1d0
      return
 100  continue
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
      do i=3,njets+3
         cwgt=cwgt*cfac(abs(ifl(i)))
         if(ifl(i).eq.0) then
            ifl(i)=21
            ng=ng+1
         endif
      enddo
c     evaluate spin weight factors
      swgt=2e0
      swgt=swgt**(njets+1)
*
      xlum=dble(slum*cwgt*swgt*ifact(ng)) /resc**(njets-1)
      xlum=xlum*effco
*
      end
c
c-------------------------------------------------------------------
      subroutine selflav2(jproc,xlum,afl)
c     evaluates parton luminosities, and assigns PDG-code flavours
c     gluon=21 cbar=-4 sbar=-3 ubar=-2 dbar=-1 d=1 u=2 s=3 c=4
c     jproc
c
c   at present CKM assumed diagonal
c
c---t q final state (+ up to 4 gluons)
c
c  1  q qbar' -> t bbar (+ gluons)                        1
c  2  qbar q' -> t bbar (+ gluons)                        2
c
c---t + 2q final state (+ gluons)
c
c  3  g q    -> t b q' (+ gluons)                         3
c  4  q g    -> t b q' (+ gluons)                         4
c
c-----------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'top.inc'
c commons
      integer icconj
      common/hwconv/icconj
      integer afl(maxpar)
      real tmp(100),slum,cwgt,swgt,rn,tmptot
      real fcab(2)
      integer icab(-4:4,2) 
c     cabibbo partner: (*,1)=cabibbo allowed, (*,2)-cabibbo suppressed
c     where 1=d 2=u 3=s 4=c 0=g and negatives are antiparticles
      data icab/-3,-4,-1,-2,0,2,1,4,3, 
     +          -1,-2,-3,-4,0,4,3,2,1/
      integer iqpp(2,2,4)
c     iqpp(k,j,i):  if q=i, q'=cabibbo partner of i (allowed for j=1,
c     suppressed for j=2) then iqpp are the two (k=1,2) flavours .ne. (q.or.q')
      data iqpp/
c       du   dc   ud   us   sc   su   cs   cd   
     +  3,4, 2,3, 3,4, 1,4, 1,2, 1,4, 1,2, 2,3/
c
      double precision xlum,xrn 
      integer i,k,itmp,icount,init,jproc,ng,j,k1,j1
      real cfac(0:6),ifact(0:6)
      data cfac/8e0,6*3e0/,init/0/
      integer imap(-4:4)
      data imap/-2,-1,-2,-1,0,1,2,1,2/
      integer ibb,iflaux
      integer nlb
      parameter (nlb=13)
      integer i1,i2,nflv1,nflv2,i3,nf1,nf2
      integer nflpr42(20),nflpr52(20),nflpr62(20)
      integer labfl42(4,nlb),labfl52(5,nlb),labfl62(6,nlb)
      character*2 lbflcr42(4,4),lbflcr52(5,8),lbflcr62(6,4),lbc
c
c     overall efficiency for extraction of colour states:
c     >(non-zero color states) / 3**nq*8**ng
c     ccoef(i,j) for njets=j and i=#(light quark pairs)
      double precision effco
      double precision ccoef(3,0:8) !ccoeff(# of q-qb pairs,  # of gluons)
      data ccoef/0.333333333,0.185185185,0.127572016,
     +           0.1666666667,0.12037037,0.0936213992,
     +           0.114583333,0.0902777778,0.0743312757,
     +           0.087239583,0.072337963,0.06171232,
     +           0.070475260,0.0603841146,0.05278461,
     +           0.059122721,0.051834672,-1.d0,
     +           0.050923665,0.045412134,-1.d0,
     +           0.042060375,-1.d0,-1.d0,
     +           0.037214041,-1.d0,-1.d0/
      integer nlq  ! number of light quarks
      save init,cfac,ifact,nlq,iqpp,icab,fcab
      save labfl42,nflpr42,labfl52,nflpr52,labfl62,nflpr62
c
      data lbflcr42 /'uq','db','tq','bb', 'cq','sb','tq','bb', !jproc=1,2 
     +               'ub','dq','tb','bq', 'cb','sq','tb','bq'/ !jproc=1,2
      data lbflcr52/ 'gg','uq','tq','bb','dq', 'gg','cq','tq','bb','sq', !jproc=3,4 
     +               'gg','ub','tb','bq','db', 'gg','cb','tb','bq','sb', !jproc=3,4
     +               'gg','dq','tb','bq','uq', 'gg','sq','tb','bq','cq', !jproc=3,4 
     +               'gg','db','tq','bb','ub', 'gg','sb','tq','bb','cb'/ !jproc=3,4

      data nflpr42 /4,19*-100/
      data nflpr52 /8,19*-100/
c
      if(init.eq.0) then
         fcab(1)=real(ccab2)
         fcab(2)=real(scab2)
         njets=npart-3
         init=1
         ifact(0)=1e0
         do i=1,6
            ifact(i)=ifact(i-1)/real(i)
         enddo
c
c  initializing flavours for the various classes of processes
         do j1=1,2
            if(j1.eq.1) nflv1=nflpr42(1)
            if(j1.eq.2) nflv1=nflpr52(1)
c
            do j=1,nflv1
               do k=1,4+(j1-1)
                  if(j1.eq.1) lbc=lbflcr42(k,j)
                  if(j1.eq.2) lbc=lbflcr52(k,j)
                  if(lbc.eq.'gg') then
                     i3=0
                  elseif(lbc.eq.'dq') then
                     i3=1
                  elseif(lbc.eq.'uq') then
                     i3=2
                  elseif(lbc.eq.'sq') then
                     i3=3
                  elseif(lbc.eq.'cq') then
                     i3=4
                  elseif(lbc.eq.'bq') then
                     i3=5
                  elseif(lbc.eq.'tq') then
                     i3=6
                  elseif(lbc.eq.'db') then
                     i3=-1
                  elseif(lbc.eq.'ub') then
                     i3=-2
                  elseif(lbc.eq.'sb') then
                     i3=-3
                  elseif(lbc.eq.'cb') then
                     i3=-4
                  elseif(lbc.eq.'bb') then
                     i3=-5
                  elseif(lbc.eq.'tb') then
                     i3=-6
                  else
                     write(*,*)'something wrong in lbflcr4'
                     stop
                  endif
                  if(j1.eq.1) labfl42(k,j)=i3
                  if(j1.eq.2) labfl52(k,j)=i3
               enddo
            enddo 
         enddo
c
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
c  1  q qbar'-> t bb + (njets-1) g     q=u,d,c,s    njets>0
c  2  qbar q'-> tb b + (njets-1) g
c
      if(jproc.le.2) then
         nflv1=1
         nflv2=nflpr42(1)
         if(jproc.eq.1) then
           i1=1
           i2=2
         else
           i1=2
           i2=1
         endif
         i=0
         do j=nflv1,nflv2                       !! FCAB????
           i=i+1
           nf1=labfl42(i1,j)
           nf2=labfl42(i2,j)
           tmp(i)=f1(nf1)*f2(nf2)            ! *fcab(k)
           slum=slum+tmp(i)
         enddo
         rn=rn*slum
         i=0
         do j=nflv1,nflv2
           i=i+1
           tmptot=tmptot+tmp(i)
           if(tmptot.ge.rn) then 
             ifl(1)=labfl42(i1,j)
             ifl(2)=labfl42(i2,j)
             ifl(3)=labfl42(3,j)
             ifl(4)=labfl42(4,j)
             do k1=5,njets+3
                ifl(k1)=0
             enddo
             effco = ccoef(2,njets-1)
c             endif
             goto 100
           endif
         enddo
c
c  3  g q -> t bb q'+ (njets-2) g        q=u,d,c,s    njets>1
c  4  q g -> t bb q'+ (njets-2) g
c
      elseif(jproc.ge.3) then
         nflv1=1
         nflv2=nflpr52(1)
         if(jproc.eq.3) then
           i1=1
           i2=2
         else
           i1=2
           i2=1
         endif
         i=0
         do j=nflv1,nflv2                       !! FCAB????
           i=i+1
           nf1=labfl52(i1,j)
           nf2=labfl52(i2,j)
           tmp(i)=f1(nf1)*f2(nf2)            ! *fcab(k)
           slum=slum+tmp(i)
         enddo
         rn=rn*slum
         i=0
         do j=nflv1,nflv2
           i=i+1
           tmptot=tmptot+tmp(i)
           if(tmptot.ge.rn) then 
             ifl(1)=labfl52(i1,j)
             ifl(2)=labfl52(i2,j)
             ifl(3)=labfl52(3,j)
             ifl(4)=labfl52(4,j)
             ifl(5)=labfl52(5,j)
             do k1=6,njets+3
                ifl(k1)=0
             enddo
c             itmp=1
c             ibl=ifl(i1)
             effco = ccoef(2,njets-1)
             goto 100
           endif
         enddo
c
      else
         write(*,*) 'jproc not defined, slum=0'
         slum=0
         stop
      endif
c
      xlum=-1d0
      return
 100  continue
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
      do i=3,njets+3
         cwgt=cwgt*cfac(abs(ifl(i)))
         if(ifl(i).eq.0) then
            ifl(i)=21
            ng=ng+1
         endif
      enddo
c     evaluate spin weight factors
      swgt=2e0
      swgt=swgt**(njets+1)
*
      xlum=dble(slum*cwgt*swgt*ifact(ng)) /resc**(njets-1)
      xlum=xlum*effco
*
      end
c
c-------------------------------------------------------------------
      subroutine selflav3(jproc,xlum,afl)
c     evaluates parton luminosities, and assigns PDG-code flavours
c     gluon=21 cbar=-4 sbar=-3 ubar=-2 dbar=-1 d=1 u=2 s=3 c=4
c     jproc
c
c   at present CKM assumed diagonal
c
c---t + W final state (+ gluons)
c
c  1  g b    -> t W- (+ gluons)                                      1
c  2  b g    -> t W- (+ gluons)                                      2
c  3  g bbar -> tbar W+ (+ gluons)                                   1
c  4  bbar g -> tbar W+ (+ gluons)                                   2
c
c---t + W + q final state (+ gluons)
c
c  5  q b    -> t W- q (+ gluons)   (q=u,d,c,s,ubar,dbar,cbar,sbar)  3
c  6  b q    -> t W- q (+ gluons)                                    4
c  7  qbar b -> t W- qbar + (njets-1) g                              5
c  8  b qbar -> t W- qbar + (njets-1) g                              6
c  9  q bbar -> tbar W+ q + (njets-1) g                              3
c 10  bbar q -> tbar W+ q + (njets-1) g                              4
c 11  qbar bbar -> tbar W+ qbar + (njets-1) g                        5
c 12  bbar qbar -> tbar W+ qbar + (njets-1) g                        6
c
c-----------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'top.inc'
c commons
      integer icconj
      common/hwconv/icconj
      integer afl(maxpar)
      real tmp(100),slum,cwgt,swgt,rn,tmptot
      real fcab(2)
      integer icab(-4:4,2) 
c     cabibbo partner: (*,1)=cabibbo allowed, (*,2)-cabibbo suppressed
c     where 1=d 2=u 3=s 4=c 0=g and negatives are antiparticles
      data icab/-3,-4,-1,-2,0,2,1,4,3, 
     +          -1,-2,-3,-4,0,4,3,2,1/
      integer iqpp(2,2,4)
c     iqpp(k,j,i):  if q=i, q'=cabibbo partner of i (allowed for j=1,
c     suppressed for j=2) then iqpp are the two (k=1,2) flavours .ne. (q.or.q')
      data iqpp/
c       du   dc   ud   us   sc   su   cs   cd   
     +  3,4, 2,3, 3,4, 1,4, 1,2, 1,4, 1,2, 2,3/
c
      double precision xlum,xrn 
      integer i,k,itmp,icount,init,jproc,ng,j,k1,j1
      real cfac(0:30),ifact(0:6)
      data cfac/8e0,6*3e0,24*1e0/,init/0/
      integer imap(-4:4)
      data imap/-2,-1,-2,-1,0,1,2,1,2/
      integer ibb,iflaux

      integer nlb
      parameter (nlb=16)
      integer i1,i2,nflv1,nflv2,i3,nf1,nf2
      integer nflpr43(20),nflpr53(20)
      integer labfl43(4,nlb),labfl53(5,nlb)
      character*2 lbflcr43(4,2),lbflcr53(5,16),lbc
c
c     overall efficiency for extraction of colour states:
c     >(non-zero color states) / 3**nq*8**ng
c     ccoef(i,j) for njets=j and i=#(light quark pairs)
      double precision effco
      double precision ccoef(3,0:8) !ccoeff(# of q-qb pairs,  # of gluons)
      data ccoef/0.333333333,0.185185185,0.127572016,
     +           0.1666666667,0.12037037,0.0936213992,
     +           0.114583333,0.0902777778,0.0743312757,
     +           0.087239583,0.072337963,0.06171232,
     +           0.070475260,0.0603841146,0.05278461,
     +           0.059122721,0.051834672,-1.d0,
     +           0.050923665,0.045412134,-1.d0,
     +           0.042060375,-1.d0,-1.d0,
     +           0.037214041,-1.d0,-1.d0/
      integer nlq  ! number of light quarks
      save init,cfac,ifact,nlq,iqpp,icab,fcab
      save labfl43,nflpr43,labfl53,nflpr53
c
      data lbflcr43/ 'gg','bq','tq','wm', 'gg','bb','tb','wp'/ !jproc=1,2,3,4 
      data lbflcr53/ 'uq','bq','tq','uq','wm', 'dq','bq','tq','dq','wm',  !jproc=5,6 
     +               'cq','bq','tq','cq','wm', 'sq','bq','tq','sq','wm',
     +               'ub','bq','tq','ub','wm', 'db','bq','tq','db','wm',  !jproc=7,8
     +               'cb','bq','tq','cb','wm', 'sb','bq','tq','sb','wm',
     +               'uq','bb','tb','uq','wp', 'dq','bb','tb','dq','wp',  !jproc=9,10
     +               'cq','bb','tb','cq','wp', 'sq','bb','tb','sq','wp',
     +               'ub','bb','tb','ub','wp', 'db','bb','tb','db','wp',  !jproc=11,12
     +               'cb','bb','tb','cb','wp', 'sb','bb','tb','sb','wp'/
      data nflpr43 /1,2,18*-100/
      data nflpr53 /0,4,8,12,16,15*-100/
c
      if(init.eq.0) then
         fcab(1)=real(ccab2)
         fcab(2)=real(scab2)
         njets=npart-4
         init=1
         ifact(0)=1e0
         do i=1,6
            ifact(i)=ifact(i-1)/real(i)
         enddo
c
c  initializing flavours for the various classes of processes
         do j1=1,2
            if(j1.eq.1) nflv1=nflpr43(2)
            if(j1.eq.2) nflv1=nflpr53(5)
c
            do j=1,nflv1
               do k=1,4+(j1-1)
                  if(j1.eq.1) lbc=lbflcr43(k,j)
                  if(j1.eq.2) lbc=lbflcr53(k,j)
                  if(lbc.eq.'gg') then
                     i3=0
                  elseif(lbc.eq.'dq') then
                     i3=1
                  elseif(lbc.eq.'uq') then
                     i3=2
                  elseif(lbc.eq.'sq') then
                     i3=3
                  elseif(lbc.eq.'cq') then
                     i3=4
                  elseif(lbc.eq.'bq') then
                     i3=5
                  elseif(lbc.eq.'tq') then
                     i3=6
                  elseif(lbc.eq.'db') then
                     i3=-1
                  elseif(lbc.eq.'ub') then
                     i3=-2
                  elseif(lbc.eq.'sb') then
                     i3=-3
                  elseif(lbc.eq.'cb') then
                     i3=-4
                  elseif(lbc.eq.'bb') then
                     i3=-5
                  elseif(lbc.eq.'tb') then
                     i3=-6
                  elseif(lbc.eq.'wm') then
                     i3=-24
                  elseif(lbc.eq.'wp') then
                     i3=+24
                  else
                     write(*,*)'something wrong in lbflcr4'
                     stop
                  endif
                  if(j1.eq.1) labfl43(k,j)=i3
                  if(j1.eq.2) labfl53(k,j)=i3
               enddo
            enddo 
         enddo
c
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
c  1  g b -> t W + njets g        q=u,d,c,s    njets>=0
c  2  b g -> t W + njets g
c
      if(jproc.le.4) then
         if(mod(jproc,2).eq.1) then
           i1=1
           i2=2
         else
           i1=2
           i2=1
         endif
c         i=0
c         do j=nflv1,nflv2                       !! FCAB????
c           i=i+1
         if(jproc.le.2) then
           j=1
         else
           j=2
         endif
         nf1=labfl43(i1,j)
         nf2=labfl43(i2,j)
         slum=f1(nf1)*f2(nf2)  
c         enddo
c         rn=rn*slum
c         i=0
c         do j=nflv1,nflv2
c           i=i+1
c           tmptot=tmptot+tmp(i)
c           if(tmptot.ge.rn) then 
         ifl(1)=labfl43(i1,j)
         ifl(2)=labfl43(i2,j)
         ifl(3)=labfl43(3,j)
         ifl(njets+4)=labfl43(4,j)
         do k1=4,njets+3
            ifl(k1)=0
         enddo
         effco = ccoef(1,njets+1)
         goto 100
c      enddo
c
c  5  q b -> t W- q + (njets-1) g        q=u,d,c,s,ub,db,cb,sb    njets>=1
c  6  b q -> t W- q + (njets-1) g
c  7  qbar b -> t W- qbar + (njets-1) g
c  8  b qbar -> t W- qbar + (njets-1) g
c  9  q bbar -> tbar W+ q + (njets-1) g
c 10  bbar q -> tbar W+ q + (njets-1) g
c 11  qbar bbar -> tbar W+ qbar + (njets-1) g
c 12  bbar qbar -> tbar W+ qbar + (njets-1) g
c
      elseif(jproc.ge.5.and.jproc.le.12) then
         i1=(jproc-3)/2
         i2=mod(jproc,2)
         nflv1=nflpr53(i1)+1
         nflv2=nflpr53(i1+1)
         if(i2.eq.1) then
           i1=1
           i2=2
         else
           i1=2
           i2=1
         endif
         i=0
         do j=nflv1,nflv2                       !! FCAB????
           i=i+1
           nf1=labfl53(i1,j)
           nf2=labfl53(i2,j)
           tmp(i)=f1(nf1)*f2(nf2)            ! *fcab(k)
           slum=slum+tmp(i)
         enddo
         rn=rn*slum
         i=0
         do j=nflv1,nflv2
           i=i+1
           tmptot=tmptot+tmp(i)
           if(tmptot.ge.rn) then 
             ifl(1)=labfl53(i1,j)
             ifl(2)=labfl53(i2,j)
             ifl(3)=labfl53(3,j)
             ifl(4)=labfl53(4,j)
             ifl(njets+4)=labfl53(njets+4,j)
             do k1=5,njets+3
                ifl(k1)=0
             enddo
c             itmp=1
c             ibl=ifl(i1)
             effco = ccoef(2,njets-1)
             goto 100
           endif
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
      do i=3,njets+3
         cwgt=cwgt*cfac(abs(ifl(i)))
         if(ifl(i).eq.0) then
            ifl(i)=21
            ng=ng+1
         endif
      enddo
c     evaluate spin weight factors
      swgt=2e0
      swgt=swgt**(njets+1)
      swgt=swgt*3e0              !for W
*
      xlum=dble(slum*cwgt*swgt*ifact(ng)) /resc**(njets+1)
      xlum=xlum*effco
*
      end
c
c-------------------------------------------------------------------
      subroutine selflav4(jproc,xlum,afl)
c     evaluates parton luminosities, and assigns PDG-code flavours
c     gluon=21 cbar=-4 sbar=-3 ubar=-2 dbar=-1 d=1 u=2 s=3 c=4
c     jproc
c
c   at present CKM assumed diagonal
c
c---t + b + W final states
c
c  1  g g    -> t bbar W- (+ gluons)                         1
c  2  g g    -> tbar b W+ (+ gluons)                         1
c  3  q qbar -> t bbar W- (+ gluons)                         2
c  4  qbar q -> t bbar W- (+ gluons)                         3
c  5  q qbar -> tbar b W+ (+ gluons)                         2
c  6  qbar q -> tbar b W+ (+ gluons)                         3
c
c---t + b + W + q final state (+ gluons)
c
c  7  g q    -> t bbar W- q (+ gluons)   (q=u,d,c,s)         4
c  8  q g    -> t bbar W- q (+ gluons)                       5
c  9  g qbar -> t bbar W- qbar (+ gluons)                    6
c  10 qbar g -> t bbar W- qbar (+ gluons)                    7
c  11 g q    -> tbar b W+ q (+ gluons)                       4
c  12 q g    -> tbar b W+ q (+ gluons)                       5
c  13 g qbar -> tbar b W+ qbar (+ gluons)                    6
c  14 qbar g -> tbar b W+ qbar (+ gluons)                    7
c
c---t + b + W + q qbar final state
c
c  15 g g    -> t bbar W- q qbar (+ gluons) (q=u,d,c,s)      1
c  16 g g    -> tbar b W+ q qbar (+ gluons) (q=u,d,c,s)      1
c
c-----------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'top.inc'
c commons
      integer icconj
      common/hwconv/icconj
      integer afl(maxpar)
      real tmp(100),slum,cwgt,swgt,rn,tmptot
      real fcab(2)
      integer icab(-4:4,2) 
c     cabibbo partner: (*,1)=cabibbo allowed, (*,2)-cabibbo suppressed
c     where 1=d 2=u 3=s 4=c 0=g and negatives are antiparticles
      data icab/-3,-4,-1,-2,0,2,1,4,3, 
     +          -1,-2,-3,-4,0,4,3,2,1/
      integer iqpp(2,2,4)
c     iqpp(k,j,i):  if q=i, q'=cabibbo partner of i (allowed for j=1,
c     suppressed for j=2) then iqpp are the two (k=1,2) flavours .ne. (q.or.q')
      data iqpp/
c       du   dc   ud   us   sc   su   cs   cd   
     +  3,4, 2,3, 3,4, 1,4, 1,2, 1,4, 1,2, 2,3/
c
      double precision xlum,xrn 
      integer i,k,itmp,icount,init,jproc,ng,j,k1,j1
      real cfac(0:30),ifact(0:6)
      data cfac/8e0,6*3e0,24*1e0/,init/0/
      integer imap(-4:4)
      data imap/-2,-1,-2,-1,0,1,2,1,2/
      integer ibb,iflaux

      integer nlb
      parameter (nlb=16)
      integer i1,i2,nflv1,nflv2,i3,nf1,nf2
      integer nflpr54(20),nflpr64(20),nflpr74(20)
      integer labfl54(5,nlb),labfl64(6,nlb),labfl74(7,nlb)
      character*2 lbflcr54(5,10),lbflcr64(6,16),lbflcr74(7,8),lbc
c
c     overall efficiency for extraction of colour states:
c     >(non-zero color states) / 3**nq*8**ng
c     ccoef(i,j) for njets=j and i=#(light quark pairs)
      double precision effco
      double precision ccoef(3,0:8) !ccoeff(# of q-qb pairs,  # of gluons)
      data ccoef/0.333333333,0.185185185,0.127572016,
     +           0.1666666667,0.12037037,0.0936213992,
     +           0.114583333,0.0902777778,0.0743312757,
     +           0.087239583,0.072337963,0.06171232,
     +           0.070475260,0.0603841146,0.05278461,
     +           0.059122721,0.051834672,-1.d0,
     +           0.050923665,0.045412134,-1.d0,
     +           0.042060375,-1.d0,-1.d0,
     +           0.037214041,-1.d0,-1.d0/
      integer nlq  ! number of light quarks
      save init,cfac,ifact,nlq,iqpp,icab,fcab
      save labfl54,nflpr54,labfl64,nflpr64,labfl74,nflpr74
c
      data lbflcr54/ 'gg','gg','tq','bb','wm', !jproc=1
     +               'gg','gg','tb','bq','wp', !jproc=2 
     +               'uq','ub','tq','bb','wm', 'dq','db','tq','bb','wm', !jproc=3,4
     +               'cq','cb','tq','bb','wm', 'sq','sb','tq','bb','wm', 
     +               'uq','ub','tb','bq','wp', 'dq','db','tb','bq','wp', !jproc=5,6
     +               'cq','cb','tb','bq','wp', 'sq','sb','tb','bq','wp'/ 
      data lbflcr64/ 'gg','uq','tq','bb','uq','wm',  !jproc=7,8
     +               'gg','dq','tq','bb','dq','wm', 
     +               'gg','cq','tq','bb','cq','wm', 
     +               'gg','sq','tq','bb','sq','wm', 
     +               'gg','ub','tq','bb','ub','wm',  !jproc=9,10
     +               'gg','db','tq','bb','db','wm', 
     +               'gg','cb','tq','bb','cb','wm', 
     +               'gg','sb','tq','bb','sb','wm', 
     +               'gg','uq','tb','bq','uq','wp',  !jproc=11,12
     +               'gg','dq','tb','bq','dq','wp', 
     +               'gg','cq','tb','bq','cq','wp', 
     +               'gg','sq','tb','bq','sq','wp', 
     +               'gg','ub','tb','bq','ub','wp',  !jproc=13,14
     +               'gg','db','tb','bq','db','wp', 
     +               'gg','cb','tb','bq','cb','wp', 
     +               'gg','sb','tb','bq','sb','wp'/
      data lbflcr74/ 'gg','gg','tq','bb','uq','ub','wm',  !jproc=15
     +               'gg','gg','tq','bb','dq','db','wm',
     +               'gg','gg','tq','bb','cq','cb','wm',
     +               'gg','gg','tq','bb','sq','sb','wm',
     +               'gg','gg','tb','bq','uq','ub','wp',  !jproc=16
     +               'gg','gg','tb','bq','dq','db','wp',
     +               'gg','gg','tb','bq','cq','cb','wp',
     +               'gg','gg','tb','bq','sq','sb','wp'/

      data nflpr54 /0,1,2,6,10,15*-100/
      data nflpr64 /0,4,8,12,16,15*-100/
      data nflpr74 /0,4,8,17*-100/
c
      if(init.eq.0) then
         fcab(1)=real(ccab2)
         fcab(2)=real(scab2)
         njets=npart-4
         init=1
         ifact(0)=1e0
         do i=1,6
            ifact(i)=ifact(i-1)/real(i)
         enddo
c
c  initializing flavours for the various classes of processes
         do j1=1,3
            if(j1.eq.1) nflv1=nflpr54(5)
            if(j1.eq.2) nflv1=nflpr64(5)
            if(j1.eq.3) nflv1=nflpr74(3)
c
            do j=1,nflv1
               do k=1,5+(j1-1)
                  if(j1.eq.1) lbc=lbflcr54(k,j)
                  if(j1.eq.2) lbc=lbflcr64(k,j)
                  if(j1.eq.3) lbc=lbflcr74(k,j)
                  if(lbc.eq.'gg') then
                     i3=0
                  elseif(lbc.eq.'dq') then
                     i3=1
                  elseif(lbc.eq.'uq') then
                     i3=2
                  elseif(lbc.eq.'sq') then
                     i3=3
                  elseif(lbc.eq.'cq') then
                     i3=4
                  elseif(lbc.eq.'bq') then
                     i3=5
                  elseif(lbc.eq.'tq') then
                     i3=6
                  elseif(lbc.eq.'db') then
                     i3=-1
                  elseif(lbc.eq.'ub') then
                     i3=-2
                  elseif(lbc.eq.'sb') then
                     i3=-3
                  elseif(lbc.eq.'cb') then
                     i3=-4
                  elseif(lbc.eq.'bb') then
                     i3=-5
                  elseif(lbc.eq.'tb') then
                     i3=-6
                  elseif(lbc.eq.'wm') then
                     i3=-24
                  elseif(lbc.eq.'wp') then
                     i3=+24
                  else
                     write(*,*)'something wrong in lbflcr4'
                     stop
                  endif
                  if(j1.eq.1) labfl54(k,j)=i3
                  if(j1.eq.2) labfl64(k,j)=i3
                  if(j1.eq.3) labfl74(k,j)=i3
               enddo
            enddo 
         enddo
c
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
c  1  g g    -> t bbar W- (+ gluons)
c  2  g g    -> tbar b W+ (+ gluons)
c  3  q qbar -> t bbar W- (+ gluons)
c  4  qbar q -> t bbar W- (+ gluons)
c  5  q qbar -> tbar b W+ (+ gluons)
c  6  qbar q -> tbar b W+ (+ gluons)
c
      if(jproc.le.6) then
         if(jproc.le.2) then
            i1=1
            i2=2
            itmp=jproc
         else
            if(mod(jproc,2).eq.0) then
               i2=1
               i1=2
            else
               i1=1
               i2=2
            endif
            itmp=2+(jproc-1)/2
         endif
         nflv1=nflpr54(itmp)+1
         nflv2=nflpr54(itmp+1)
         i=0
         do j=nflv1,nflv2                       !! FCAB????
           i=i+1
           nf1=labfl54(i1,j)
           nf2=labfl54(i2,j)
           tmp(i)=f1(nf1)*f2(nf2)            ! *fcab(k)
           slum=slum+tmp(i)
         enddo
         rn=rn*slum
         i=0
         do j=nflv1,nflv2
           i=i+1
           tmptot=tmptot+tmp(i)
           if(tmptot.ge.rn) then 
             ifl(1)=labfl54(i1,j)
             ifl(2)=labfl54(i2,j)
             ifl(3)=labfl54(3,j)
             ifl(4)=labfl54(4,j)
             ifl(njets+4)=labfl54(5,j)
             do k1=5,njets+3
                ifl(k1)=0
             enddo
             if(jproc.le.2) then
                effco = ccoef(1,njets+1)
             else
                effco = ccoef(2,njets-1)
             endif
             goto 100
           endif
         enddo
c
c  7  g q    -> t bbar W- q (+ gluons)   (q=u,d,c,s)
c  8  q g    -> t bbar W- q (+ gluons)
c  9  g qbar -> t bbar W- qbar (+ gluons)
c  10 qbar g -> t bbar W- qbar (+ gluons)
c  11 g q    -> tbar b W+ q (+ gluons)
c  12 q g    -> tbar b W+ q (+ gluons)
c  13 g qbar -> tbar b W+ qbar (+ gluons)
c  14 qbar g -> tbar b W+ qbar (+ gluons)
c
      elseif(jproc.ge.7.and.jproc.le.14) then
         i1=(jproc-5)/2
         i2=mod(jproc,2)
         nflv1=nflpr64(i1)+1
         nflv2=nflpr64(i1+1)
         if(i2.eq.1) then
           i1=1
           i2=2
         else
           i1=2
           i2=1
         endif
         i=0
         do j=nflv1,nflv2                       !! FCAB????
           i=i+1
           nf1=labfl64(i1,j)
           nf2=labfl64(i2,j)
           tmp(i)=f1(nf1)*f2(nf2)            ! *fcab(k)
           slum=slum+tmp(i)
         enddo
         rn=rn*slum
         i=0
         do j=nflv1,nflv2
           i=i+1
           tmptot=tmptot+tmp(i)
           if(tmptot.ge.rn) then 
             ifl(1)=labfl64(i1,j)
             ifl(2)=labfl64(i2,j)
             ifl(3)=labfl64(3,j)
             ifl(4)=labfl64(4,j)
             ifl(5)=labfl64(5,j)
             ifl(njets+4)=labfl64(6,j)
             do k1=6,njets+3
                ifl(k1)=0
             enddo
             effco = ccoef(2,njets-1)
             goto 100
           endif
         enddo
c  15 g g    -> t bbar W- q qbar (+ gluons)     (q=u,d,c,s)
c  16 g g    -> tbar b W+ q qbar (+ gluons)
c
      elseif(jproc.ge.15.and.jproc.le.16) then
         i1=1
         i2=2
         itmp=jproc-14
         nflv1=nflpr74(itmp)+1
         nflv2=nflpr74(itmp+1)
         i=0
         do j=nflv1,nflv2                       !! FCAB????
           i=i+1
           nf1=labfl74(i1,j)
           nf2=labfl74(i2,j)
           tmp(i)=f1(nf1)*f2(nf2)            ! *fcab(k)
           slum=slum+tmp(i)
         enddo
         rn=rn*slum
         i=0
         do j=nflv1,nflv2
           i=i+1
           tmptot=tmptot+tmp(i)
           if(tmptot.ge.rn) then 
             ifl(1)=labfl74(i1,j)
             ifl(2)=labfl74(i2,j)
             ifl(3)=labfl74(3,j)
             ifl(4)=labfl74(4,j)
             ifl(5)=labfl74(5,j)
             ifl(6)=labfl74(6,j)
             ifl(njets+4)=labfl74(7,j)
             do k1=7,njets+3
                ifl(k1)=0
             enddo
             effco = ccoef(2,njets-1)
             goto 100
           endif
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
      do i=3,njets+3
         cwgt=cwgt*cfac(abs(ifl(i)))
         if(ifl(i).eq.0) then
            ifl(i)=21
            ng=ng+1
         endif
      enddo
c     evaluate spin weight factors
      swgt=2e0
      swgt=swgt**(njets+1)
      swgt=swgt*3e0              !for W
*
      xlum=dble(slum*cwgt*swgt*ifact(ng)) /resc**(njets+1)
      xlum=xlum*effco
*
      end
*
      subroutine setdec(nwrt,iflwrt,icuwrt,pwrt,decbr)
      implicit none
      include 'alpgen.inc'
      include 'top.inc'
c debug
      double precision rate(16)
      common/dbg/rate
      data rate/16*0d0/
c locals
      integer maxdec
      parameter (maxdec=40)
      integer ip, ic, il, irn,id1,id2,i,itmp,ich,itdmode
     $     ,itdec0,iwdmode,iwdec0
      integer iwfl(3)
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
      double precision br(6)
      data br/3*0.111111111,0.333333333,0.66666666666,1d0/
      double precision dmass(16)
      data dmass/16*0d0/
c arguments
      real*8 xmw
c
      integer nwrt,iflwrt(maxdec),icuwrt(2,maxdec)
      double precision pwrt(5,maxdec),decbr
c
      real*8 ptopu(4),pbtopu(4),pftopu(4),pfbtopu(4)
c
      save itdmode,iwdmode,dmass
c
      integer ievnt
      data ievnt/0/
      ievnt=ievnt+1
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
        itdmode=itdecmode
c        write(*,*) 'select top decay modes:'
c        write(*,*) '1: e nu b'
c        write(*,*) '2: mu nu b'
c        write(*,*) '3: tau nu b'
c        write(*,*) '4: e/mu/tau nu b'
c        write(*,*) '5: b + 2 jets'
c        write(*,*) '6: fully inclusive'
c        read(*,*) itdmode
        if(itopprc.ge.3) then
          iwdmode=iwdecmode
c           write(*,*) 'select W decay modes:'
c           write(*,*) '1: e nu'
c           write(*,*) '2: mu nu'
c           write(*,*) '3: tau nu'
c           write(*,*) '4: e/mu/tau nu'
c           write(*,*) '5: 2 jets'
c           write(*,*) '6: fully inclusive'
c           read(*,*) iwdmode
        endif
        dmass(4)=1.5
        dmass(5)=mb
        dmass(11)=0.5d-3
        dmass(13)=0.10566d0
        dmass(15)=1.777d0
        init=1
      endif
*
* top decay
*
      decbr=br(itdmode)
      if(ifl(3).eq.6) then
         ich=1
         iwfl(3)=5
      elseif(ifl(3).eq.-6) then
         ich=0
         iwfl(3)=-5
      else
         write(*,*) 'problems in setdec'
         stop
      endif
c freeze top dec mode:
      itdec0=itdmode
c     select flavours of W decay products
 1    if(itdmode.le.3) then
        iwfl(1)=ilepc(itdmode)+ich   ! tbar->l- or t->nu
        iwfl(2)=-inu(itdmode)+ich  !  tbar->nu_lbar or t->l+
      elseif(itdmode.eq.4) then
        xrn=rangen2(1)
        irn=min(1+int(3*xrn),3) 
        iwfl(1)=ilepc(irn)+ich   ! tbar->l- or t->nu
        iwfl(2)=-inu(irn)+ich  !  tbar->nu_lbar or t->l+
      elseif(itdmode.eq.5) then
        xrn=rangen2(1)
        irn=min(1+int(2*xrn),2)   !irn=1,2
        itmp=idn(irn)+ich 
        iwfl(1)=itmp      ! tbar->d/s or t->u/c
        xrn=rangen2(1)
        if(xrn.lt.scab2) then
          iwfl(2)=-icab(itmp,2)
        else
          iwfl(2)=-icab(itmp,1)
        endif
      elseif(itdmode.eq.6) then
        xrn=rangen2(1)
        if(xrn.lt.0.333333333) then
          itdmode=4
        else 
          itdmode=5
        endif
        goto 1
      endif
c     restore dec mode label
      itdmode=itdec0
c
      nwrt=nwrt+2
      do il=1,4
         pwrt(il,nwrt-1)=idec(il,2,1)
         pwrt(il,nwrt)=idec(il,3,1)
         pdec(il)=pwrt(il,nwrt-1)+pwrt(il,nwrt)
      enddo
      iflwrt(nwrt-1)=iwfl(1)
      iflwrt(nwrt)=iwfl(2)
      m1=dmass(abs(iflwrt(nwrt-1)))
      m2=dmass(abs(iflwrt(nwrt)))
      pwrt(5,nwrt-1)=0
      pwrt(5,nwrt)=0
c      print*,'new event ',ievnt, ' mw=',sqrt(abs(pdec(4)**2-pdec(1)**2
c     $     -pdec(2)**2-pdec(3)**2))
c      print*,sqrt(abs(idec(4,1,1)**2-idec(1,1,1)**2-idec(2,1,1)**2
c     $     -idec(3,1,1)**2))
      pdec(5)=mw
c   rescale momenta to incorporate mass effects in the decays
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

c
      if(itopprc.le.2) return
*
* W decay
*
      decbr=decbr*br(iwdmode)
      if(ifl(npart).eq.24) then
         ich=1
      elseif(ifl(npart).eq.-24) then
         ich=0
      else
         write(*,*) 'problems in setdec'
         stop
      endif
c freeze top dec mode:
      iwdec0=iwdmode
c     select flavours of W decay products
 2    if(iwdmode.le.3) then
        iwfl(1)=ilepc(iwdmode)+ich   ! W- ->l- or W+ ->nu
        iwfl(2)=-inu(iwdmode)+ich  !  W- ->nu_lbar or W+ ->l+
      elseif(iwdmode.eq.4) then
        xrn=rangen2(1)
        irn=min(1+int(3*xrn),3) 
        iwfl(1)=ilepc(irn)+ich   ! tbar->l- or t->nu
        iwfl(2)=-inu(irn)+ich  !  tbar->nu_lbar or t->l+
      elseif(iwdmode.eq.5) then
        xrn=rangen2(1)
        irn=min(1+int(2*xrn),2)   !irn=1,2
        itmp=idn(irn)+ich 
        iwfl(1)=itmp      ! tbar->d/s or t->u/c
        xrn=rangen2(1)
        if(xrn.lt.scab2) then
          iwfl(2)=-icab(itmp,2)
        else
          iwfl(2)=-icab(itmp,1)
        endif
      elseif(iwdmode.eq.6) then
        xrn=rangen2(1)
        if(xrn.lt.0.333333333) then
          iwdmode=4
        else 
          iwdmode=5
        endif
        goto 2
      endif
c     restore dec mode label
      iwdmode=iwdec0
c
      nwrt=nwrt+2
      do il=1,4
         pwrt(il,nwrt-1)=wdec(il,1)
         pwrt(il,nwrt)=wdec(il,2)
         pdec(il)=p(il,npart)
      enddo
*
      iflwrt(nwrt-1)=iwfl(1)
      iflwrt(nwrt)=iwfl(2)
      m1=dmass(abs(iflwrt(nwrt-1)))
      m2=dmass(abs(iflwrt(nwrt)))
      pwrt(5,nwrt-1)=0
      pwrt(5,nwrt)=0
      pdec(5)=mw
c   rescale momenta to incorporate mass effects in the decays
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

c$$$c  debug
c$$$      do i=nwrt-3,nwrt
c$$$        do il=1,16
c$$$          if(abs(iflwrt(i)).eq.il) rate(il)=rate(il)+1d0
c$$$        enddo
c$$$      enddo

      do i=1,4
         ptopu(i)= pwrt(i,3)
         pftopu(i)=pwrt(i,nwrt-3)
         pfbtopu(i)=pwrt(i,nwrt-2)
         pbtopu(i)=idec(i,1,1)
      enddo
      
c      print*,'ptop   = ',(ptopu(i),i=1,4)
c      print*,'pb     = ',(pbtopu(i),i=1,4)
c      print*,'pftop  = ',(pftopu(i),i=1,4)
c      print*,'pfbtop = ',(pfbtopu(i),i=1,4)
c      
      
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
      include 'top.inc'
      real*8 dummy,djg
      real *8 pswgt
      real *8 djpd,factor
      real *8 cutkin(10)
      real *8 wgt
      common/loccut/cutkin
      real *8 pl(maxpar),y(maxpar)
      real *8 pcm(0:3,maxpar)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
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
         if(itopprc.eq.1.or.itopprc.eq.3) then
            cutkin(1)=ptjmin
            cutkin(3)=etajmax
            cutkin(5)=drjmin
            cutkin(6)=0.d0
            cutkin(7)=ptj1min
            cutkin(8)=ptj1max
            cutkin(9)=0.d0
            cutkin(10)=0.d0
         elseif(itopprc.eq.2.or.itopprc.eq.4) then
            cutkin(1)=ptjmin
            cutkin(2)=ptbmin
            cutkin(3)=etajmax
            cutkin(4)=etabmax
            cutkin(5)=drjmin
            cutkin(6)=0.d0
            cutkin(7)=ptj1min
            cutkin(8)=ptj1max
            cutkin(9)=0.d0
            cutkin(10)=0.d0
         else
            print*,'itopprc ',itopprc,' not yet available'
         endif
c-
         ninit=1
      endif
*
c-    The generation starts
*

      pswgt=0.d0
      if(itopprc.eq.1) call momgen1(njets,mq,mq2,roots,x1,x2,pcm,
     +                              wgt,lnot)
      if(itopprc.eq.2) call momgen2(njets,mq,mq2,roots,x1,x2,pcm,
     +                              wgt,lnot)
      if(itopprc.eq.3) call momgen3(njets,mq,mq2,roots,x1,x2,pcm,
     +                              wgt,lnot)
      if(itopprc.eq.4) call momgen4(njets,mq,mq2,roots,x1,x2,pcm,
     +                              wgt,lnot)
      djg= 1.d0 ! dummy variable
      if (lnot.eq.1) then
         pswgt= 0.d0
         goto 100
      endif
*
c-    will write factor=factor0/(x1*x2), with factor0 function of njets etc.
*
      if(itopprc.le.2) then
         factor= 1d0/(2.d0*pi)**(3*(njets+1)-4)/2.d0/s/x1/x2
      elseif(itopprc.ge.3) then
         factor= 1d0/(2.d0*pi)**(3*(njets+1+1)-4)/2.d0/s/x1/x2
      endif
*
c-    initial state momenta in the LAB frame:
*
c      pcm(0,1)= roots/2.d0*x1   
c      pcm(1,1)= 0.d0   
c      pcm(2,1)= 0.d0   
c      pcm(3,1)= roots/2.d0*x1   
c
c      pcm(0,2)= roots/2.d0*x2   
c      pcm(1,2)= 0.d0   
c      pcm(2,2)= 0.d0   
c      pcm(3,2)=-roots/2.d0*x2   

c-    Rapidities and pseudo-rapidities (in the LAB system), pt's, deltar's:

      do l= 3,njets+3
          pt(l) = sqrt(pcm(1,l)**2+pcm(2,l)**2)
          pl(l) = pcm(3,l)
          eta(l)= -log(tan(0.5d0*atan2(pt(l),pl(l))))
          y(l)  = 0.5d0*log((pcm(0,l)+pcm(3,l))/(pcm(0,l)-pcm(3,l)))
      enddo
*
c-    Calculates jet-jet distances:
*
      do l= 3,njets+2
         do m= l+1,njets+3
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
      do l= 3,njets+2
         do m= l+1,njets+3
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
*
      if(itopprc.eq.1) then
         if(jproc.eq.1.or.jproc.eq.4) then
            p(5,1)= mq2
            p(5,2)= 0.d0
            p(5,4)= 0.d0
         elseif(jproc.eq.2.or.jproc.eq.3) then
            p(5,2)= mq2
            p(5,1)= 0.d0
            p(5,4)= 0.d0
         endif
      elseif(itopprc.eq.3) then
         p(5,npart)=mw
         if(mod(jproc,2).eq.0) then
            p(5,1)= mq2
            p(5,2)= 0.d0
         else
            p(5,2)= mq2
            p(5,1)= 0.d0
         endif
      elseif(itopprc.eq.2.or.itopprc.eq.4) then
         if(itopprc.eq.4) p(5,npart)=mw
         p(5,1)= 0.d0
         p(5,2)= 0.d0
         p(5,4)= mq2
      else
         print*,'jproc not yet defined',jproc
         stop
      endif
*
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
      call chkcut(lnot,pt,p,eta,dr,njets)
*
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
      do i=3,njets+3
        totpt=totpt+pt(i)**2+p(5,i)**2
      enddo 
      if(itopprc.ge.3) then
        totpt=totpt+mw**2
      endif
      if(iqopt.eq.0) then 
         qsq=1d0
      elseif(iqopt.eq.1) then
         qsq=totpt
      elseif(iqopt.eq.2) then
         qsq=x1*x2*roots**2
      endif
      qsq=qfac**2*qsq
      
 100  continue
      end
*
      subroutine chkcut(lnot,pt,p,eta,dr,njets)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     Applies kinematical cuts to the final state during the phase
c     -space generation                                          c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 ptjmin,drjmin,etabmax,etajmax,ptbmin
      integer maxpar,ninit,njets,j,lnot,i
      real*8 cutkin(10)
      real*8 ptj1min,ptj1max,highjet
      common/loccut/cutkin
      data ninit/0/
      parameter (maxpar=10)
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),p(5,maxpar)
      integer itopprc
      common/singletop/itopprc
      real*8 mtt
c particle parameters
      real *8 amass(24),mlep(3),mc,mb,mt,mw,wwid,mz,zwid,mh,hwid
      common/pparam/amass,mlep,mc,mb,mt,mw,wwid,mz,zwid,mh,hwid
      save ninit,ptjmin,ptbmin,etajmax,etabmax,drjmin
      save ptj1min,ptj1max
      if(ninit.eq.0) then
         ninit=1
         ptjmin=cutkin(1)
         ptbmin=cutkin(2)
         etajmax=cutkin(3)
         etabmax=cutkin(4)
         drjmin=cutkin(5)
         ptj1min=cutkin(7)
         ptj1max=cutkin(8)
      endif

      lnot= 0
c     impose leading jet cut
      highjet = pt(3)
      do i=4,njets+3
         if(pt(i).gt.highjet) highjet=pt(i)
      enddo
      if (ptj1max.gt.0.and.highjet.gt.ptj1max) goto 10
      if (ptj1min.gt.0.and.highjet.lt.ptj1min) goto 10

      if(itopprc.eq.2.or.itopprc.eq.4) then
         if (pt(4).lt.ptbmin)           goto 10
         if (abs(eta(4)).gt.etabmax)    goto 10
      endif
c
      if(itopprc.eq.1) then
c     impose minimum pt and require eta within allowed range, for jets
         do i=4,njets+3
            if (pt(i).lt.ptjmin)           goto 10
            if (abs(eta(i)).gt.etajmax)    goto 10
         enddo
c     require dR(jet-jet)<drjmin
         do i=4,njets+2
            do j=i+1,njets+3
               if(dr(i,j).lt.drjmin)  goto 10
            enddo
         enddo
      elseif(itopprc.eq.2) then
         do i=5,njets+3
            if (pt(i).lt.ptjmin)           goto 10
            if (abs(eta(i)).gt.etajmax)    goto 10
         enddo
c     require dR(jet-jet)<drjmin
         do i=4,njets+2
            do j=i+1,njets+3
               if(dr(i,j).lt.drjmin)  goto 10
            enddo
         enddo
      elseif(itopprc.eq.3) then
         do i=4,njets+3
            if (pt(i).lt.ptjmin)           goto 10
            if (abs(eta(i)).gt.etajmax)    goto 10
         enddo
c     require dR(jet-jet)<drjmin
         do i=4,njets+2
            do j=i+1,njets+3
               if(dr(i,j).lt.drjmin)  goto 10
            enddo
         enddo
      elseif(itopprc.eq.4) then
         do i=5,njets+3
            if (pt(i).lt.ptjmin)           goto 10
            if (abs(eta(i)).gt.etajmax)    goto 10
         enddo
c     require dR(jet-jet)<drjmin
         do i=4,njets+2
            do j=i+1,njets+3
               if(dr(i,j).lt.drjmin)  goto 10
            enddo
         enddo

         mtt= (p(4,4)+p(4,njets+4))**2-(p(1,4)+p(1,njets+4))**2
     +       -(p(2,4)+p(2,njets+4))**2-(p(3,4)+p(3,njets+4))**2
         mtt= sqrt(mtt)

         if(abs(mtt-mt).lt.5.d0) goto 10

      endif

 5    return

 10   lnot= 1
      return
      end
*
      subroutine momgen1(njets,qm1,qm2,roots,x1,x2,p,wgt,lw)
*
c-    Generator of np particles in the LAB frame.
*
      implicit none
      real*8 qm1,qm2,roots,x1,x2,wgt,wgt1,wgt2,zero,etacut
      real*8 qm12,qm22
      real*8 ptjmin,cutkin,etajmax,drjmin,s,stw,ptbmin,etabmax
      real*8 ag,tau0,tau,xmsum,xmsum0,ptlim
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
      real*8 pin(0:3,2)
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
*
      data mpar/0/
      save
*
      if (mpar.eq.0) then
        mpar   = 1
        np     = njets+1
c        nw     = np-1
        zero   = 0.d0
        qm12   = qm1*qm1
        qm22   = qm2*qm2
        etacut = 40.d0
        s      = roots*roots
*
        ptjmin = cutkin(1)
        etajmax= cutkin(3)
        drjmin = cutkin(5)
*
c-      Parameters for the light jets:
*
c        pt0_l(2) = ptbmin
c        pt1_l(2) = roots/2.d0
c        eta0(2)  = etabmax
c        xm(2)    = qm2          !b mass
        do j= 2,np
          pt0_l(j) = ptjmin
          pt1_l(j) = roots/2.d0
          eta0(j)  = etajmax
          xm(j)    = zero
        enddo
*
c-      Parameters for the the top system:
*
        pt0_l(1) = zero
        pt1_l(1) = roots/2.d0
        eta0(1)  = etacut
        xm(1)    = qm1               !t mass
*
c-      Parameters for the x1,x2 integration:
*
        ag  = 0.98d0
*
c-      sets the ratios of calls mom/momr
*
        apw(1)= 0.9d0
        apw(2)= 1.d0-apw(1)
      endif
*
      xm(2)    = zero
*
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
c-    Generates sw=stt:
*
c      call photsm(0,cntt,cxmw,cxpw,sw,dj1,ranram(2*nw+1))
*
c-    The upgrading of xmsum
*
c      xmsum = xmsum0+xm(nw)
      xmsum = xmsum0
*
c-    tau0 is the lower cut on tau= x1*x2
*
      tau0   = 1.d0/s*(xmsum)**2
      tau0   = max(tau0,4.d0*pt0lmax**2/s)  
      tau0   = max(tau0,pt0lsum**2/s)
      tau0   = max(tau0,1.d0/s*(qm1
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
*
      pin(0,1)= roots/2.d0*x1   
      pin(1,1)= 0.d0   
      pin(2,1)= 0.d0   
      pin(3,1)= roots/2.d0*x1   

      pin(0,2)= roots/2.d0*x2   
      pin(1,2)= 0.d0   
      pin(2,2)= 0.d0   
      pin(3,2)=-roots/2.d0*x2   
      if(jproc.eq.1.or.jproc.eq.4) then
         if(pin(0,1).lt.qm2) goto 100
         pin(3,1)= sqrt((roots/2.d0*x1)**2-qm2*qm2)
      else if(jproc.eq.2.or.jproc.eq.3) then
         if(pin(0,2).lt.qm2) goto 100
         pin(3,2)= -sqrt((roots/2.d0*x2)**2-qm2*qm2)
      endif


c      rootsh= roots*sq
      rootsh= sqrt((pin(0,1)+pin(0,2))**2-(pin(3,1)+pin(3,2))**2)
*
c-    Protection:
*
      if (rootsh.lt.xmsum) goto 100
*
c-    Rescalings and transformations to feed mom:
*
      do j= 1,np
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
      pz    = (pin(3,1)+pin(3,2))/roots
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
c          if (j.ne.nw) then
            p(k,j+2)= pm(k,j)*roots
c          else
c            pm(k,j) = pm(k,j)*roots
c          endif
        enddo
      enddo
*
      do j=1,2
         do k= 0,3
            p(k,j)= pin(k,j)
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
      subroutine momgen2(njets,qm1,qm2,roots,x1,x2,p,wgt,lw)
*
c-    Generator of np particles in the LAB frame.
*
      implicit none
      real*8 qm1,qm2,roots,x1,x2,wgt,wgt1,wgt2,zero,etacut
      real*8 qm12,qm22
      real*8 ptjmin,cutkin,etajmax,drjmin,s,stw,ptbmin,etabmax
      real*8 ag,tau0,tau,xmsum,xmsum0,ptlim
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
      real*8 pin(0:3,2)
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
*
      data mpar/0/
      save
*
      if (mpar.eq.0) then
        mpar   = 1
        np     = njets+1
c        nw     = np-1
        zero   = 0.d0
        qm12   = qm1*qm1
        qm22   = qm2*qm2
        etacut = 40.d0
        s      = roots*roots
*
        ptjmin = cutkin(1)
        ptbmin = cutkin(2)
        etajmax= cutkin(3)
        etabmax= cutkin(4)
        drjmin = cutkin(5)
*
c-      Parameters for the light and b jets:
*
        pt0_l(2) = ptbmin
        pt1_l(2) = roots/2.d0
        eta0(2)  = etabmax
        xm(2)    = qm2          !b mass
        do j= 3,np
          pt0_l(j) = ptjmin
          pt1_l(j) = roots/2.d0
          eta0(j)  = etajmax
          xm(j)    = zero
        enddo
*
c-      Parameters for the the top system:
*
        pt0_l(1) = zero
        pt1_l(1) = roots/2.d0
        eta0(1)  = etacut
        xm(1)    = qm1               !t mass
*
c-      Parameters for the x1,x2 integration:
*
        ag  = 0.98d0
*
c-      sets the ratios of calls mom/momr
*
        apw(1)= 0.9d0
        apw(2)= 1.d0-apw(1)
      endif
*
      xm(3)    = zero
*
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
c-    Generates sw=stt:
*
c      call photsm(0,cntt,cxmw,cxpw,sw,dj1,ranram(2*nw+1))
*
c-    The upgrading of xmsum
*
c      xmsum = xmsum0+xm(nw)
      xmsum = xmsum0
*
c-    tau0 is the lower cut on tau= x1*x2
*
      tau0   = 1.d0/s*(xmsum)**2
      tau0   = max(tau0,4.d0*pt0lmax**2/s)  
      tau0   = max(tau0,pt0lsum**2/s)
      tau0   = max(tau0,1.d0/s*(qm1+qm2+ptbmin
     +                  +float(njets-1)*ptjmin)**2)
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
*
      pin(0,1)= roots/2.d0*x1   
      pin(1,1)= 0.d0   
      pin(2,1)= 0.d0   
      pin(3,1)= roots/2.d0*x1   

      pin(0,2)= roots/2.d0*x2   
      pin(1,2)= 0.d0   
      pin(2,2)= 0.d0   
      pin(3,2)=-roots/2.d0*x2   
      rootsh= roots*sq
c      rootsh= sqrt((pin(0,1)+pin(0,2))**2-(pin(3,1)+pin(3,2))**2)
*
c-    Protection:
*
      if (rootsh.lt.xmsum) goto 100
*
c-    Rescalings and transformations to feed mom:
*
      do j= 1,np
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
      pz = 0.5d0*(x1-x2)
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
c          if (j.ne.nw) then
            p(k,j+2)= pm(k,j)*roots
c          else
c            pm(k,j) = pm(k,j)*roots
c          endif
        enddo
      enddo
*
      do j=1,2
         do k= 0,3
            p(k,j)= pin(k,j)
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
      subroutine momgen3(njets,qm1,qm2,roots,x1,x2,p,wgt,lw)
*
c-    Generator of np particles in the LAB frame.
*
      implicit none
c particle parameters
      real *8 amass(24),mlep(3),mc,mb,mt,mw,wwid,mz,zwid,mh,hwid
      common/pparam/amass,mlep,mc,mb,mt,mw,wwid,mz,zwid,mh,hwid
      real*8 qm1,qm2,roots,x1,x2,wgt,wgt1,wgt2,zero,etacut
      real*8 qm12,qm22,mw2
      real*8 ptjmin,cutkin,etajmax,drjmin,s,stw,ptbmin,etabmax,cntw
      real*8 ag,tau0,tau,xmsum,xmsum0,ptlim
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
      real*8 pin(0:3,2),paux(0:3)
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      integer i
*
      data mpar/0/
      save
*
      if (mpar.eq.0) then
        mpar   = 1
        np     = njets+2
        nw     = np-1
        zero   = 0.d0
        qm12   = qm1*qm1
        qm22   = qm2*qm2
        mw2    = mw*mw
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
c-      Parameters for the the top+W system:
*
        cntw= 1.d0
        pt0_l(nw) = zero
        pt1_l(nw) = roots/2.d0
        eta0(nw)  = etacut
        xm(nw)    = zero               !provisional (t+W) mass
        cxmw      = (qm1+mw)*(qm1+mw)
        cxpw      = s
*
c-      Parameters for the x1,x2 integration:
*
        ag  = 0.98d0
*
c-      sets the ratios of calls mom/momr
*
        apw(1)= 0.9d0
        apw(2)= 1.d0-apw(1)
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

        call photsm(0,cntw,cxmw,cxpw,stw,dj1,ranram(3))
        tau= stw/s
        y0 =-0.5d0*log(tau) 
        sq = sqrt(tau)
        call ttriangle(0,-y0,y0,-y0,y0,yr,wjr,ranram(1),lw1)
        x1 = sq*exp(yr) 
        x2 = sq*exp(-yr) 
*
        pin(0,1)= roots/2.d0*x1   
        pin(1,1)= 0.d0   
        pin(2,1)= 0.d0   
        pin(3,1)= roots/2.d0*x1   
        
        pin(0,2)= roots/2.d0*x2   
        pin(1,2)= 0.d0   
        pin(2,2)= 0.d0   
        pin(3,2)=-roots/2.d0*x2   
        if(jproc.eq.1.or.jproc.eq.3) then
           if(pin(0,2).lt.qm2) goto 100
           pin(3,2)= -sqrt((roots/2.d0*x2)**2-qm2*qm2)
        else if(jproc.eq.2.or.jproc.eq.4) then
           if(pin(0,1).lt.qm2) goto 100
           pin(3,1)= sqrt((roots/2.d0*x1)**2-qm2*qm2)
        endif
*
        pm(0,nw)= pin(0,1)+pin(0,2)
        pm(1,nw)= 0.d0
        pm(2,nw)= 0.d0
        pm(3,nw)= pin(3,1)+pin(3,2)
*
        ran(1)= ranram(2)
        call rans(ran0)
        ran(2)= ran0

        stw= pm(0,nw)**2-pm(3,nw)**2
        if(stw.lt.(qm1+mw)*(qm1+mw)) goto 100

        call dec2fm(0,stw,pm(0,nw),
     +              qm12,mw2,p(0,nw+2),p(0,nw+3),dj2,ran)
*
        wgt= wjr/s/dj1/dj2/djbintot        
        do j=1,2
           do k= 0,3
              p(k,j)= pin(k,j)
           enddo
        enddo

        return
      endif
*
c-    Generates sw=stw:
*
      call photsm(0,cntw,cxmw,cxpw,stw,dj1,ranram(2*nw+1))
*
c-    The upgrading of xmsum
*
      xm(nw)= sqrt(stw)
      xmsum = xmsum0+xm(nw)
*
c-    tau0 is the lower cut on tau= x1*x2
*
      tau0   = 1.d0/s*(xmsum)**2
      tau0   = max(tau0,4.d0*pt0lmax**2/s)  
      tau0   = max(tau0,pt0lsum**2/s)
      tau0   = max(tau0,1.d0/s*(qm1+mw
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
*
      pin(0,1)= roots/2.d0*x1   
      pin(1,1)= 0.d0   
      pin(2,1)= 0.d0   
      pin(3,1)= roots/2.d0*x1   

      pin(0,2)= roots/2.d0*x2   
      pin(1,2)= 0.d0   
      pin(2,2)= 0.d0   
      pin(3,2)=-roots/2.d0*x2   
      if(mod(jproc,2).eq.0) then
         if(pin(0,1).lt.qm2) goto 100
         pin(3,1)= sqrt((roots/2.d0*x1)**2-qm2*qm2)
      else 
         if(pin(0,2).lt.qm2) goto 100
         pin(3,2)= -sqrt((roots/2.d0*x2)**2-qm2*qm2)
      endif
*c      rootsh= roots*sq
      rootsh= sqrt((pin(0,1)+pin(0,2))**2-(pin(3,1)+pin(3,2))**2)
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
      pz    =(pin(3,1)+pin(3,2))/roots
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
      do j=1,2
         do k= 0,3
            p(k,j)= pin(k,j)
         enddo
      enddo
*
      wgt= wgt*s**(nw-2)     
*
      ran(1)= ranram(2*nw)
      call rans(ran0)
      ran(2)= ran0
      call dec2fm(0,stw,pm(0,nw),
     +            qm12,mw2,p(0,3),p(0,4),dj2,ran)
*
      wgt= wgt/dj1/dj2
*
* momenta are stored in the order: t,b,jets,W
*
      do i=0,3
         paux(i)=p(i,4)
      enddo
      do i=0,3
         do j=4,njets+3
            p(i,j)=p(i,j+1)
         enddo
         p(i,np+2)= paux(i)
      enddo
*
      return
 100  lw= 1
      wgt= 0.d0
      return
      end
*
      subroutine momgen4(njets,qm1,qm2,roots,x1,x2,p,wgt,lw)
*
c-    Generator of np particles in the LAB frame.
*
      implicit none
c particle parameters
      real *8 amass(24),mlep(3),mc,mb,mt,mw,wwid,mz,zwid,mh,hwid
      common/pparam/amass,mlep,mc,mb,mt,mw,wwid,mz,zwid,mh,hwid
      real*8 qm1,qm2,roots,x1,x2,wgt,wgt1,wgt2,zero,etacut
      real*8 qm12,qm22,mw2
      real*8 ptjmin,cutkin,etajmax,drjmin,s,stw,ptbmin,etabmax,cntw
      real*8 ag,tau0,tau,xmsum,xmsum0,ptlim
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
      real*8 pin(0:3,2),paux(0:3)
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      integer i
*
      data mpar/0/
      save
*
      if (mpar.eq.0) then
        mpar   = 1
        np     = njets+2
        nw     = np-1
        zero   = 0.d0
        qm12   = qm1*qm1
        qm22   = qm2*qm2
        mw2    = mw*mw
        etacut = 40.d0
        s      = roots*roots
*
        ptjmin = cutkin(1)
        ptbmin = cutkin(2)
        etajmax= cutkin(3)
        etabmax= cutkin(4)
        drjmin = cutkin(5)
*
c-      Parameters for the b jets:
*
        pt0_l(1) = ptbmin
        pt1_l(1) = roots/2.d0
        eta0(1)  = etabmax
        xm(1)    = qm2
*
c-      Parameters for the light jets:
*
        do j= 2,nw-1
          pt0_l(j) = ptjmin
          pt1_l(j) = roots/2.d0
          eta0(j)  = etajmax
          xm(j)    = zero
        enddo
*
c-      Parameters for the the top+W system:
*
        cntw= 1.d0
        pt0_l(nw) = zero
        pt1_l(nw) = roots/2.d0
        eta0(nw)  = etacut
        xm(nw)    = zero               !provisional (t+W) mass
        cxmw      = (qm1+mw)*(qm1+mw)
        cxpw      = (roots-qm2)**2
*
c-      Parameters for the x1,x2 integration:
*
        ag  = 0.98d0
*
c-      sets the ratios of calls mom/momr
*
        apw(1)= 0.9d0
        apw(2)= 1.d0-apw(1)
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
c-    Generates sw=stw:
*
      call photsm(0,cntw,cxmw,cxpw,stw,dj1,ranram(2*nw+1))
*
c-    The upgrading of xmsum
*
      xm(nw)= sqrt(stw)
      xmsum = xmsum0+xm(nw)
*
c-    tau0 is the lower cut on tau= x1*x2
*
      tau0   = 1.d0/s*(xmsum)**2
      tau0   = max(tau0,4.d0*pt0lmax**2/s)  
      tau0   = max(tau0,pt0lsum**2/s)
      tau0   = max(tau0,1.d0/s*(qm1+mw+ptbmin
     +                  +float(njets-1)*ptjmin)**2)
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
*
      pin(0,1)= roots/2.d0*x1   
      pin(1,1)= 0.d0   
      pin(2,1)= 0.d0   
      pin(3,1)= roots/2.d0*x1   

      pin(0,2)= roots/2.d0*x2   
      pin(1,2)= 0.d0   
      pin(2,2)= 0.d0   
      pin(3,2)=-roots/2.d0*x2   
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
      do j=1,2
         do k= 0,3
            p(k,j)= pin(k,j)
         enddo
      enddo
*
      wgt= wgt*s**(nw-2)     
*
      ran(1)= ranram(2*nw)
      call rans(ran0)
      ran(2)= ran0
      call dec2fm(0,stw,pm(0,nw),
     +            qm12,mw2,p(0,3),p(0,4),dj2,ran)
*
      wgt= wgt/dj1/dj2
*
* momenta are stored in the order: t,jets,W
*
      do i=0,3
         paux(i)=p(i,4)
      enddo
      do i=0,3
         do j=4,njets+3
            p(i,j)=p(i,j+1)
         enddo
         p(i,np+2)= paux(i)
      enddo
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
      include 'top.inc'
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
      if(itopprc.le.2) then
         nv  = 2*(njets+1)-1
      else
         nv  = 2*(njets+2)-1
      endif
      nvar= nv*ngrid
      do n= 2,nvar+1
        nch= 10
*
        nct(n)   =  nct(n-1)+nx1(n-1)   
        nx1(n)   =  nch 
        if(itopprc.ne.3) then
           if (mod(n-1,nv).eq.1)    mask(n)= 1
           if (mod(n-1,nv).eq.2)    mask(n)= 1
        elseif(itopprc.eq.3) then
           if (mod(n-1,nv).eq.1)    mask(n)= 1
           if (mod(n-1,nv).eq.2)    mask(n)= 1
           if (mod(n-1,nv).eq.nv-1) mask(n)= 1
           if (mod(n-1,nv).eq.0)    mask(n)= 1
        else
           write(*,*) 'problems with itopprc'
           stop
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
      include 'top.inc'
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
      end
c
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
      data inpint/ 24,          ! N V-A 
     + 1,1,1,1,1,   1,1,2,1,2,   2,1,1,1,2,   3,1,2,1,1, ! zuu, zdd, w+ud, w-ud,
     + 10,1,1,1,1,  10,1,2,1,2,   1,2,1,2,1,   1,2,2,2,2, ! Auu, Add, zcc, zss
     + 2,2,1,2,2,   3,2,2,2,1,  10,2,1,2,1,  10,2,2,2,2, ! W+cs, W-sc, Acc, Ass
     + 11,1,1,1,1,  11,1,2,1,2,  11,2,1,2,1,  11,2,2,2,2, ! guu, gdd, gcc, gss
     + 11,3,1,3,1,  11,3,2,3,2,  10,3,1,3,1, 10,3,2,3,2,
     + 1,3,1,3,1,    1,3,2,3,2,   
     + 2,3,1,3,2,    3,3,2,3,1,
     + 1,                   ! N of yukawa
     + 3,1,3,1,
     + 15,                  ! N self-gauge
     + 1, 2, 3,  10, 2, 3,   4, 2, 3,   5, 1, 1, ! ZWW, AWW, auxiliary
     + 5, 1,10,   5,10,10,   5, 2, 3,   6, 1, 2, ! auxiliary
     + 6,10, 2,   7, 1, 3,   7,10, 3,   8, 2, 2, ! auxiliary
     + 9, 3, 3,  11,11,11,  12,11,11, ! auxiliary, ggg, Xgg
     + 4,                   ! N H-GAUGE
     + 1, 1, 1,   1, 2, 3,   3,1,1,   3,2,3, ! HZZ, HWW, HHZZ, HHWW
     + 3,                   ! N of self higgs
     + 1, 2, 3,             ! HHH, HHHH e aux per higgs quartici
     + 811*-100       /     ! EOF signal
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

