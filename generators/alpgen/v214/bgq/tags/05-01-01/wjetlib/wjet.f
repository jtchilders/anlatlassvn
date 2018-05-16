c-------------------------------------------------------------------
      subroutine alsprc
c     assigns the hard process code
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      ihrd=3
      end


c-------------------------------------------------------------------
      subroutine alhset
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wjet.inc'
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c     ngrid is the total number of grids allowd for P.S. variables.
c     jgrid(jproc) labels the grid associated to a given jproc.
c       jgrid= 1   -->   q qbar and qbar q (no quarks in final state)
c       jgrid= 2   -->   g q    and g qbar
c       jgrid= 3   -->   q g    and qbar g
c       jgrid= 4   -->   g g
c       jgrid= 5   -->   qq->qq
      data jgrid/1,2,3,4,1,2*5,5*1,3*5,3,2,2*3,2*2,3,2,2*4,75*0/
      integer i
c
c     parameters for the gauge invariance prescription:
      winsize  = 2.d0/pi
      resonance= 'n'
      wmode    = 'yy'
c
c     process input parameters
      npart=njets+4
      nprtns=njets
      nw=1
      if(njets.eq.0) then
         jprocmax=1
         ngrid   = 1
      elseif(njets.eq.1) then
         jprocmax=3
         ngrid   = 3
      elseif(njets.eq.2) then
         jprocmax=15
         ngrid   = 5
      elseif(njets.eq.3) then
         jprocmax=23
         ngrid   = 5
      elseif(njets.eq.4) then
         jprocmax=25
         ngrid   = 5
      elseif(njets.eq.5) then
         jprocmax=25
         ngrid   = 5
      elseif(njets.eq.6) then
         jprocmax=25
         ngrid   = 5
      endif 
c
      do i=1,njets+4
         p(5,i)=0
      enddo
c
c run parameters:
c write process code, code=3 for W+jets
      if(njets.eq.0) then
        write(niosta,*) 'W production'
      else
        write(niosta,*) 'W +',njets,' jets'
      endif
      write(niosta,*) 'W-> ell nu'
      write(niosta,*) '======================================='
      write(niosta,*)
     $     'Generation cuts for the partonic event sample:' 
      write(niosta,*) '     Light jets:'
      if(njets.gt.0) then
        write(niosta,*) 'ptmin=',ptjmin,' |etamax|=',etajmax
     $       ,' dR(j-j)>',drjmin 
      endif
      write(niosta,*) '     Leptons:'
      write(niosta,*) 'ptmin(lep)=',ptlmin,' |etamax|=',etalmax
     $     ,' Et(miss)>',metmin,' dR(l-j)>',drlmin  
      write(niosta,*) 'pthrmin=',pthrmin,' pthrmax=',pthrmax
      end
                              
c-------------------------------------------------------------------
      subroutine selflav(jproc,xlum,afl)
c     evaluates parton luminosities, and assigns PDG-code flavours
c     gluon=21 cbar=-4 sbar=-3 ubar=-2 dbar=-1 d=1 u=2 s=3 c=4
c     jproc
c     1  q qbar' -> W  ,   qbar q' -> W
c     2  g q -> q' W    and    g qbar -> q'bar W
c     3  q g -> q' W    and  qbar g -> q'bar W
c     4  gg -> q qbar' W              
c     5  q qbar' -> W q'' qbar'' ,   qbar q' -> W q'' qbar''
c     6  q(bar) q(bar)'' -> W q'(bar) q(bar)''  
c     7  q''(bar) q(bar) -> W q'(bar) q''(bar)
c     8  q'' q''bar -> W q q'bar
c     9  q qbar' -> W q qbar ,   qbar q' -> W qbar q 
c     10 qbar' q -> W q qbar ,   q' qbar -> W qbar q
c     11 q qbar -> W q qbar'  + c.c.
c     12 q qbar -> W q' qbar
c     13 q q -> W q q'
c     14 q q' -> W q q
c     15 q q' -> W q' q'
c     16 q g -> W q' q'' qbar'' ,   qbar g -> W qbar' q'' qbar''
c     17 g q -> W q' q'' qbar'' ,   g qbar -> W qbar' q'' qbar''
c     18 q g -> W q q qbar'
c     19 q g -> W q' q qbar
c     20 g q -> W q q qbar' 
c     21 g q -> W q' q qbar
c     22 q g -> W q' q' qbar'
c     23 g q -> W q' q' qbar'
c     24 gg -> q qbar' q'' qbar'' W              
c     25 gg -> q qbar q qbar' W  + c.c.
c-----------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wjet.inc'
c commons
      integer icconj
      common/hwconv/icconj
c
      integer il1,il2,afl(maxpar)
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
c     parton charges
      real dch,uch,dbch,ubch,chrg(-6:6)
      parameter (dch=-0.333333333e0,uch=0.666666666e0)
      parameter (dbch=0.333333333e0,ubch=-0.666666666e0)
      data chrg/ubch,dbch,ubch,dbch,ubch,dbch,0e0,  
     +           dch, uch, dch, uch, dch, uch/
c
      real *8 xlum,xrn                                 
      integer i,j,k,l,ik,is,itmp,icount,init,jproc,ng,nlqp
      real cfac(0:6),ifact(0:6)
      data cfac/8e0,6*3e0/,init/0/
      integer imap(-4:4)
      data imap/-2,-1,-2,-1,0,1,2,1,2/
c
c     overall efficiency for extraction of colour states:
c     #(non-zero color states) / 3**nq*8**ng
c     ccoef(i,j) for njets=j and i=#(light quark pairs)
      real*8 ccoef(3,0:6)    !ccoeff(number of q-pb pairs, number of gluons)
      data ccoef/  0.333333333,0.185185185,0.127572016,
     $             0.166666667,0.12037037,0.0936213992,
     $             0.114583333,0.0902777778,0.0743312757,
     $             0.0872395833,0.072337963,-1.d0,
     $             0.0704752604,0.0603841146,-1.d0,
     $             0.0591227214,-1d0,-1d0,
     $             0.0509236654,-1d0,-1d0/
c 2 quark pairs, ng=0,1,2,3:
c    15/81,0d0, 78/648,0d0, 468/5184, 93/729, 
c 3 quark pairs, ng=0,1,2,3:
c    3000/41472, 546/5832, 20034/331776, 3468/46656
c 1 quark pair, ng=0,1,2,3,4,5,6:
c  3/9, 12/72, 66/576, 402/4608, 2598/36864, 17436/294912, 120144/2359296
c
      integer il,inu,iflag
      save init,chrg,cfac,ifact,iqpp,icab,fcab
c
      if(init.eq.0) then
         njets=npart-4
         init=1
         ifact(0)=1e0
         do i=1,6
            ifact(i)=ifact(i-1)/real(i)
         enddo
         fcab(1)=real(ccab2)
         fcab(2)=real(scab2)
      endif
c
      do i=1,npart
         ifl(i)=0
         afl(i)=0
      enddo 
 1    call randa(xrn)
      rn=real(xrn)
      if(1e0-rn.lt.1e-7) goto 1

      icount=0
      slum=0e0
      tmptot=0e0
      do i=1,npart
         afl(i)=0
      enddo
c
c     start two light-quark processes
c
c     jproc=1   njets>=0
c     q qbar' -> W  ,   qbar q' -> W
      if(jproc.eq.1) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(-icab(i,j))*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=-icab(i,j)
                     do k=3,njets+2
                        ifl(k)=0
                     enddo
c     d ubar -> W-
                     itmp=imap(ifl(1))
                     afl(1)=imap(ifl(1))
                     afl(2)=imap(ifl(2))
                     goto 100
                  endif 
               enddo
            endif
         enddo 
c
c     jproc=2   njets>=1
c     g q -> q' W    and    g qbar -> q'bar W
      elseif(jproc.eq.2) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(0)*f2(i)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=0
                     ifl(2)=i
                     ifl(3)=icab(i,j)
                     do k=4,njets+2
                        ifl(k)=0
                     enddo 
c     g d -> u W-
                     itmp=imap(ifl(2))
                     afl(2)=imap(ifl(2))
                     afl(3)=imap(ifl(3))
                     goto 100
                  endif 
               enddo
            endif
         enddo 
c
c     jproc=3   njets>=1
c     q g -> q' W    and  qbar g -> q'bar W
      elseif(jproc.eq.3) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(0)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=0
                     ifl(3)=icab(i,j)
                     do k=4,njets+2
                        ifl(k)=0
                     enddo 
c     d g -> u W-
                     itmp=imap(ifl(1))
                     afl(1)=imap(ifl(1))
                     afl(3)=imap(ifl(3))
                     goto 100
                  endif 
               enddo
            endif
         enddo 
c     jproc=4   njets>=2
c     gg -> q qbar' W              
      elseif(jproc.eq.4) then
         do i=1,4
            do j=1,2
               icount=icount+1
               tmp(icount)=fcab(j)
               slum=slum+tmp(icount)
            enddo
         enddo 
         icount=0
         rn=rn*slum
         slum=slum*f1(0)*f2(0)
         do i=1,4
            do j=1,2
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=0
                  ifl(3)=i
                  ifl(4)=-icab(i,j)
                  do k=5,njets+2
                     ifl(k)=0
                  enddo 
c     gg -> u dbar W-
                  itmp=-imap(ifl(4))
                  afl(3)= imap(ifl(3))
                  afl(4)= imap(ifl(4))
                  goto 100
               endif
            enddo
         enddo 
c
c     start 4 light quark processes  (in order of increasing number of
c     minimum njets, different flavour first, like-flavours after
c
c     jproc=5  1001  njets>=2
c     q qbar' -> W q'' qbar'' ,   qbar q' -> W q'' qbar''
      elseif(jproc.eq.5) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(-icab(i,j))*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
c     include factor of two for the two possible q'' flavours .ne. q or
c     q'
         slum=2e0*slum
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  do k=1,2
c     consider the two possible q'' flavours
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=-icab(i,j)
                        ifl(3)=iqpp(k,j,abs(i))
                        ifl(4)=-ifl(3)
                        do l=5,njets+2
                           ifl(l)=0
                        enddo 
c     d ubar -> W- c cbar
                        itmp=imap(ifl(1))
                        afl(1)=imap(ifl(1))
                        afl(2)=imap(ifl(2))
                        afl(3)=sign(4,ifl(3))
                        afl(4)=-afl(3)
                        goto 100
                     endif 
                  enddo
               enddo 
            endif
         enddo
c
c     jproc=6 1002  njets>=2
c     q(bar) q(bar)'' -> W q'(bar) q(bar)''  
      elseif(jproc.eq.6) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  is=1
                  do k=1,4
                     ik=k
                     if(k.gt.2) then
                        ik=k-2
                        is=-1
                     endif
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(is*iqpp(ik,j,abs(i)))*fcab(j)
                     slum=slum+tmp(icount)
                  enddo
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  is=1
                  do k=1,4
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ik=k
                        if(k.gt.2) then
                           ik=k-2
                           is=-1
                        endif
                        ifl(1)=i
                        ifl(2)=is*iqpp(ik,j,abs(i))
                        ifl(3)=icab(i,j)
                        ifl(4)=ifl(2)
                        do l=5,njets+2
                           ifl(l)=0
                        enddo 
c     d c(bar) -> W- u c(bar)
                        itmp=imap(ifl(1))
                        afl(1)=imap(ifl(1))
                        afl(3)=imap(ifl(3))
                        afl(2)=sign(4,ifl(2))
                        afl(4)=afl(2)
                        goto 100
                     endif 
                  enddo
               enddo
            endif
         enddo 
c
c     jproc=7 1002 njets>=2
c      q''(bar) q(bar) -> W q'(bar) q''(bar)
      elseif(jproc.eq.7) then
         do i=-4,4
             if(i.ne.0) then
               do j=1,2
                  is=1
                  do k=1,4
                     ik=k
                     if(k.gt.2) then
                        ik=k-2
                        is=-1
                     endif
                     icount=icount+1
                     tmp(icount)=f1(is*iqpp(ik,j,abs(i)))*f2(i)*fcab(j)
                     slum=slum+tmp(icount)
                  enddo
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  is=1
                  do k=1,4
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ik=k
                        if(k.gt.2) then
                           ik=k-2
                           is=-1
                        endif
                        ifl(1)=is*iqpp(ik,j,abs(i))
                        ifl(2)=i
                        ifl(3)=icab(i,j)
                        ifl(4)=ifl(1)
                        do l=5,njets+2
                           ifl(l)=0
                        enddo 
c      c(bar) d -> W u c(bar)
                        itmp=imap(ifl(2))
                        afl(2)=imap(ifl(2))
                        afl(3)=imap(ifl(3))
                        afl(1)=sign(4,ifl(1))
                        afl(4)=afl(1)
                        goto 100
                     endif 
                  enddo
               enddo
            endif
         enddo 
c
c     jproc=8 njets>=2
c      q'' q''bar -> W q q'bar
      elseif(jproc.eq.8) then
         do i=-4,4
            if(i.ne.0) then
               is=1
               do k=1,4
                  if(k.ne.abs(i)) then
                     do j=1,2
                        if(icab(k,j).ne.abs(i)) then
                           icount=icount+1
                           tmp(icount)=f1(i)*f2(-i)*fcab(j)
                           slum=slum+tmp(icount)
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do k=1,4
                  if(k.ne.abs(i)) then
                     do j=1,2
                        if(icab(k,j).ne.abs(i)) then
                           icount=icount+1
                           tmptot=tmptot+tmp(icount)
                           if(tmptot.ge.rn) then 
                              ifl(1)=i
                              ifl(2)=-i
                              ifl(3)=k
                              ifl(4)=-icab(k,j)
                              do l=5,njets+2
                                 ifl(l)=0
                              enddo 
c      c cbar -> W- u dbar 
                              itmp=-imap(ifl(4))
                              afl(1)=sign(4,ifl(1))
                              afl(2)=-afl(1)
                              afl(3)=imap(ifl(3))
                              afl(4)=imap(ifl(4))
                              goto 100
                           endif
                        endif
                     enddo
                  endif
               enddo 
            endif
         enddo 
c
c    then processes with like flavours
c
c     jproc=9  njets>=2
c     q qbar' -> W q qbar ,   qbar q' -> W qbar q 
      elseif(jproc.eq.9) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(-icab(i,j))*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=-icab(i,j)
                     ifl(3)=i
                     ifl(4)=-ifl(3)
                     do l=5,njets+2
                        ifl(l)=0
                     enddo 
c     d ubar -> W- d dbar
                     itmp=imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     goto 100
                  endif 
               enddo
            endif
         enddo
c     jproc=10 njets>=2
c     qbar' q -> W q qbar ,   q' qbar -> W qbar q
      elseif(jproc.eq.10) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(-icab(i,j))*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(2)=i
                     ifl(1)=-icab(i,j)
                     ifl(3)=i
                     ifl(4)=-ifl(3)
                     do l=5,njets+2
                        ifl(l)=0
                     enddo 
c     ubar d -> W- d dbar
                     itmp=imap(ifl(2))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     goto 100
                  endif 
               enddo 
            endif
         enddo
c
c     jproc=11 njets>=2
c     q qbar -> W q qbar'  + c.c.
c
      elseif(jproc.eq.11) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(-i)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=-i
                     ifl(3)=i
                     ifl(4)=-icab(i,j)
                     do l=5,njets+2
                        ifl(l)=0
                     enddo 
c     dbar d -> W- dbar u 
                     itmp=imap(ifl(2))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     goto 100
                  endif
               enddo 
            endif
         enddo
c
c     jproc=12 nejts>=4
c     q qbar -> W q' qbar
c
      elseif(jproc.eq.12) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(-i)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=-i
                     ifl(3)=icab(i,j)
                     ifl(4)=-i
                     do l=5,njets+2
                        ifl(l)=0
                     enddo 
c     d dbar -> W- u dbar
                     itmp=imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     goto 100
                  endif
               enddo 
            endif
         enddo
c
c     jproc=13  njets>=2
c     q q -> W q q'
c
      elseif(jproc.eq.13) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(i)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=i
                     ifl(3)=i
                     ifl(4)=icab(i,j)
                     do l=5,njets+2
                        ifl(l)=0
                     enddo 
c     d d -> W- d u
                     itmp=imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     goto 100
                  endif
               enddo 
            endif
         enddo
c
c     jproc=14 njets>=2
c     q q' -> W q q
c
      elseif(jproc.eq.14) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(icab(i,j))*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=icab(i,j)
                     ifl(3)=i
                     ifl(4)=i
                     do l=5,njets+2
                        ifl(l)=0
                     enddo 
c     dbar ubar -> W- dbar dbar
                     itmp=-imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     goto 100
                  endif
               enddo 
            endif
         enddo
c
c     jproc=15 njets>=2
c     q q' -> W q' q'
c
      elseif(jproc.eq.15) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(icab(i,j))*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=icab(i,j)
                     ifl(3)=ifl(2)
                     ifl(4)=ifl(2)
                     do l=5,njets+2
                        ifl(l)=0
                     enddo 
c     ubar dbar -> W dbar dbar 
                     itmp=-imap(ifl(2))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     goto 100
                  endif
               enddo 
            endif
         enddo
c
c     jproc=16 njets>=3
c     q g -> W q' q'' qbar'' ,   qbar g -> W qbar' q'' qbar''
      elseif(jproc.eq.16) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(0)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
c     include factor of two for the two possible q'' flavours .ne. q or
c     q'
         slum=2e0*slum
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  do k=1,2
c     consider the two possible q'' flavours
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=0
                        ifl(3)=icab(i,j)
                        ifl(4)=iqpp(k,j,abs(i))
                        ifl(5)=-ifl(4)
                        do l=6,njets+2
                           ifl(l)=0
                        enddo 
c     d g -> W- u c cbar
                        itmp=imap(ifl(1))
                        afl(1)=imap(ifl(1))
                        afl(3)=imap(ifl(3))
                        afl(4)=sign(4,ifl(4))
                        afl(5)=-afl(4)
                        goto 100
                     endif 
                  enddo
               enddo
            endif
         enddo
c
c     jproc=17  njets>=3
c     g q -> W q' q'' qbar'' ,   g qbar -> W qbar' q'' qbar''
      elseif(jproc.eq.17) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(0)*f2(i)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
c     include factor of two for the two possible q'' flavours .ne. q or
c     q'
         slum=2e0*slum
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  do k=1,2
c     consider the two possible q'' flavours
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=0
                        ifl(2)=i
                        ifl(3)=icab(i,j)
                        ifl(4)=iqpp(k,j,abs(i))
                        ifl(5)=-ifl(4)
                        do l=6,njets+2
                           ifl(l)=0
                        enddo 
c     g d -> W- u c cbar
                        itmp=imap(ifl(2))
                        afl(2)=imap(ifl(2))
                        afl(3)=imap(ifl(3))
                        afl(4)=sign(4,ifl(4))
                        afl(5)=-afl(4)
                        goto 100
                     endif 
                  enddo
               enddo 
            endif
         enddo
c
c     jproc=18  njets>=3
c     q g -> W q q qbar'
      elseif(jproc.eq.18) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(0)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=0
                     ifl(3)=i
                     ifl(4)=i
                     ifl(5)=-icab(i,j)
                     do l=6,njets+2
                        ifl(l)=0
                     enddo 
c     dbar g -> W dbar dbar u
                     itmp=-imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     afl(5)= imap(ifl(5))
                     goto 100
                  endif 
               enddo
            endif
         enddo
c
c     jproc=19  njets>=3
c     q g -> W q' q qbar
      elseif(jproc.eq.19) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(0)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=0
                     ifl(3)=icab(i,j)
                     ifl(4)=i
                     ifl(5)=-i
                     do l=6,njets+2
                        ifl(l)=0
                     enddo 
c     d g -> W u d dbar
                     itmp=imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     afl(5)= imap(ifl(5))
                     goto 100
                  endif 
               enddo
            endif
         enddo
c
c     jproc=20  njets>=3
c     g q -> W q q qbar' 
      elseif(jproc.eq.20) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(0)*f2(i)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=0
                     ifl(2)=i
                     ifl(3)=i
                     ifl(4)=i
                     ifl(5)=-icab(i,j)
                     do l=6,njets+2
                        ifl(l)=0
                     enddo 
c     g dbar -> W dbar dbar u
                     itmp=-imap(ifl(2))
                     afl(2)= imap(ifl(2))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     afl(5)= imap(ifl(5))
                     goto 100
                  endif 
               enddo
            endif
         enddo
c
c     jproc=21  njets>=3
c     g q -> W q' q qbar
      elseif(jproc.eq.21) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(0)*f2(i)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=0
                     ifl(2)=i
                     ifl(3)=icab(i,j)
                     ifl(4)=i
                     ifl(5)=-i
                     do l=6,njets+2
                        ifl(l)=0
                     enddo 
c     g d -> W u d dbar
                     itmp=imap(ifl(2))
                     afl(2)= imap(ifl(2))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     afl(5)= imap(ifl(5))
                     goto 100
                  endif 
               enddo
            endif
         enddo
c
c     jproc=22 njets>=3
c     q g -> W q' q'' qbar'' ,   qbar g -> W qbar' q'' qbar''
      elseif(jproc.eq.22) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(i)*f2(0)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=i
                     ifl(2)=0
                     ifl(3)=icab(i,j)
                     ifl(4)=ifl(3)
                     ifl(5)=-ifl(4)
                     do l=6,njets+2
                        ifl(l)=0
                     enddo 
c     d g -> W- u u ubar
                     itmp=imap(ifl(1))
                     afl(1)=imap(ifl(1))
                     afl(3)=imap(ifl(3))
                     afl(4)=sign(4,ifl(4))
                     afl(5)=-afl(4)
                     goto 100
                  endif 
               enddo
            endif
         enddo
c
c     jproc=23  njets>=3
c     g q -> W q' q' qbar' ,   g qbar -> W qbar q' qbar'
      elseif(jproc.eq.23) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=f1(0)*f2(i)*fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=0
                     ifl(2)=i
                     ifl(3)=icab(i,j)
                     ifl(4)=ifl(3)
                     ifl(5)=-ifl(4)
                     do l=6,njets+2
                        ifl(l)=0
                     enddo 
c     g d -> W- u u ubar
                     itmp=imap(ifl(2))
                     afl(2)=imap(ifl(2))
                     afl(3)=imap(ifl(3))
                     afl(4)=sign(4,ifl(4))
                     afl(5)=-afl(4)
                     goto 100
                  endif 
               enddo 
            endif
         enddo
c
c     jproc=24  njets>=4
c     gg -> q qbar' q'' qbar'' W              
      elseif(jproc.eq.24) then
         do i=1,4
            do j=1,2
               icount=icount+1
               tmp(icount)=fcab(j)
               slum=slum+tmp(icount)
            enddo
         enddo 
         icount=0
c     include factor of two for the two possible q'' flavours .ne. q or
c     q'
         slum=2e0*slum
         rn=rn*slum
         slum=slum*f1(0)*f2(0)
         do i=1,4
            do j=1,2
               icount=icount+1
               do k=1,2
c     consider the two possible q'' flavours
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=0
                     ifl(2)=0
                     ifl(3)=i
                     ifl(4)=-icab(i,j)
                     ifl(5)=iqpp(k,j,abs(i))
                     ifl(6)=-ifl(5)
                     do l=7,njets+2
                        ifl(l)=0
                     enddo 
c     gg -> u dbar c cbar  W-              
                     itmp=-imap(ifl(4))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     afl(5)= sign(4,ifl(5))
                     afl(6)=-afl(5)
                     goto 100
                  endif 
               enddo
            enddo
         enddo
c
c     jproc=25  njets>=4
c     gg -> q qbar q qbar' W  + c.c.
      elseif(jproc.eq.25) then
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmp(icount)=fcab(j)
                  slum=slum+tmp(icount)
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         slum=slum*f1(0)*f2(0)
         do i=-4,4
            if(i.ne.0) then
               do j=1,2
                  icount=icount+1
                  tmptot=tmptot+tmp(icount)
                  if(tmptot.ge.rn) then 
                     ifl(1)=0
                     ifl(2)=0
                     ifl(3)=i
                     ifl(4)=-i
                     ifl(5)=i
                     ifl(6)=-icab(i,j)
                     do k=7,njets+2
                        ifl(k)=0
                     enddo 
c     gg -> dbar d dbar u W-  + c.c.
                     itmp=-imap(ifl(3))
                     afl(3)= imap(ifl(3))
                     afl(4)= imap(ifl(4))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     goto 100
                  endif
               enddo
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
      il1=njets+3
      il2=njets+4
c
c     label type of flavour configuration:
c
c     d ubar -> e- nubar:
c     input reference leptons
      ifl(il1)=11
      ifl(il2)=-12
      il=il1
      inu=il2
c     u dbar -> nu e+    ISOSPIN flip:
      if(itmp.eq.2) then
         ifl(il1)= 12
         ifl(il2)=-11
         il=il2
         inu=il1
c     dbar u -> e+ nu    CHARGE conjugation:
      elseif(itmp.eq.-1) then
         ifl(il1)=-11
         ifl(il2)= 12
c     ubar d -> nubar e-   ISOSPIN flip and CHARGE conjugation
      elseif(itmp.eq.-2) then
         ifl(il1)=-12
         ifl(il2)= 11
         il=il2
         inu=il1
      endif 
c
      afl(il1)=ifl(il1)
      afl(il2)=ifl(il2)
c     check lepton charges
      wchrg=chrg(ifl(1))+chrg(ifl(2))
      do i=3,njets+2
         wchrg=wchrg-chrg(ifl(i))
      enddo
      if(wchrg*ifl(il).gt.0.or.abs(ifl(il)).ne.11) then
         write(*,*) 'inconsistent electron charge assignement:'
         write(*,*) 'charge:',wchrg,' type:',ifl(il)
      endif
c
c if W+, then charge conjugate the flavours for ALPHA to process
c and store information in icconj variable (will be used by SETCOL
      icconj=1
      if(wchrg.gt.0) then
         icconj=-1
         do i=1,npart
            afl(i)=-afl(i)
         enddo
      endif
c
c     complete flavour-dependent momenta assignements in the user common
c     block, and apply remaining flavour-depepdent cuts
      call usrfll(il,inu,iflag)
      if(iflag.eq.1) then
         xlum=-1d0
         return
      endif
c
c     evaluate colour weight factors
      do i=1,2
         if(ifl(i).eq.0) ifl(i)=21
      enddo
      ng=0
      cwgt=1e0
      do i=3,njets+2
         cwgt=cwgt*cfac(abs(ifl(i)))
         if(ifl(i).eq.0) then
            ifl(i)=21
            ng=ng+1
         endif
      enddo
c     evaluate spin weight factors
      swgt=2e0
      swgt=swgt**(njets+2)
c
      xlum=dble(slum*cwgt*swgt*ifact(ng))/resc**njets
      if(jproc.eq.14.or.jproc.eq.25.or.jproc.eq.20.or.jproc.eq.18.or
     $     .jproc.eq.15.or.jproc.eq.22.or.jproc.eq.23)  xlum=0.5d0*xlum ! identical quarks in final state
c     effco=ccoef(#quark pairs,#gluons)
      if(jproc.le.4) then
         nlqp=1
      elseif(jproc.le.25) then
         nlqp=2
      else
         write(*,*) 'jproc not valid'
         stop
      endif 
      xlum=xlum*ccoef(nlqp,2+njets-2*nlqp)
c
      end                                                
c     end selflav


      subroutine usrfll(il,inu,iflag)
      implicit none
      include 'alpgen.inc'
      include 'wjet.inc'
      integer i,j,il,inu,iflag
c
      iflag=0
c assign momenta to the usr common block
c
c jets first
      do i=1,4
         do j=1,2
            pin(i,j)=p(i,j)
         enddo
         do j=1,npart-2
            pout(i,j)=p(i,j+2)
         enddo
         do j=1,njets
            pjet(i,j)=p(i,j+2)
         enddo
      enddo
      do j=1,njets
         ptj(j)=pt(j+2)
         etaj(j)=eta(j+2)
         do i=j+1,njets
            drjj(i,j)=dr(i+2,j+2)
            drjj(j,i)=drjj(i,j)
         enddo
      enddo
c leptons
      do i=1,4
         plep(i)=p(i,il)
         pnu(i)=p(i,inu)
         pw(i)=plep(i)+pnu(i)
      enddo
      ptlep=pt(il)
      ptmiss=pt(inu)
      etalep=eta(il)
      do i=1,njets
         drlj(i)=dr(il,2+i)
      enddo
c
c     impose minimum pt and require eta within allowed range, for leptons
      if (ptlep.lt.ptlmin)           goto 10
      if (ptmiss.lt.metmin)           goto 10
      if (abs(etalep).gt.etalmax)    goto 10
c
c     isolate lepton from light jets
      do i=1,njets
         if(drlj(i).lt.drlmin)        goto 10
      enddo 
      return
 10   iflag=1
      end


      subroutine setdec(nwrt,iflwrt,icuwrt,pwrt,decbr)
      implicit none
      include 'alpgen.inc'
      include 'wjet.inc'
c locals
      integer maxdec
      parameter (maxdec=40)
      integer ip, ic, il, irn,i,itmp,idecmode
     $     ,idec0
      integer init
      data init/0/
      integer icab(-4:4,2) 
      double precision xrn,xrn1,xrn2,tmp,rangen2
      double precision pdec(5),m1,m2
c     cabibbo partner: (*,1)=cabibbo allowed, (*,2)-cabibbo suppressed
c     where 1=d 2=u 3=s 4=c 0=g and negatives are antiparticles
      data icab/-3,-4,-1,-2,0,2,1,4,3, 
     +          -1,-2,-3,-4,0,4,3,2,1/
c     BR's: ev wgt already include BR(W->e), br given here is correction
C     factor for other decay modes
      double precision br(6)
      data br/3*1d0,3d0,6d0,9d0/
      double precision dmass(16)
      data dmass/16*0d0/
c arguments
      integer nwrt,iflwrt(maxdec),icuwrt(2,maxdec)
      double precision pwrt(5,maxdec),decbr
c debug
      integer idbg
      double precision fcount
      common/fldbg/fcount(16),idbg
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
c
c     add the W decay products; the W itself and the b will be
C     reconstructed  from momentum conservation by the Les Houches
C     interface
c
      if(init.eq.0) then
        idecmode=iWdecmode
c        write(*,*) 'select W decay modes:'
c        write(*,*) '1: e nu'
c        write(*,*) '2: mu nu'
c        write(*,*) '3: tau nu'
c        write(*,*) '4: e/mu/tau nu'
c        write(*,*) '5: q q''bar'
c        write(*,*) '6: fully inclusive'
c        read(*,*) idecmode
        dmass(4)=1.5
        dmass(11)=0.5d-3
        dmass(13)=0.10566d0
        dmass(15)=1.777d0
        init=1
        do i=1,16
          fcount(i)=0d0
        enddo
      endif
      decbr=br(idecmode)
c freeze W dec mode:
      idec0=idecmode
c     select flavours of W decay products

 1    if(idecmode.le.3) then
        itmp=iflwrt(nwrt-1)
        iflwrt(nwrt-1)=itmp+2*sign(idecmode-1,itmp)
        itmp=iflwrt(nwrt)
        iflwrt(nwrt)=itmp+2*sign(idecmode-1,itmp)
      elseif(idecmode.eq.4) then
 12     xrn=rangen2(1)
        irn=int(3d0*xrn)
        if(irn.gt.2) goto 12
        itmp=iflwrt(nwrt-1)
        iflwrt(nwrt-1)=itmp+2*sign(irn,itmp)
        itmp=iflwrt(nwrt)
        iflwrt(nwrt)=itmp+2*sign(irn,itmp)
      elseif(idecmode.eq.5) then
 16     xrn=rangen2(1)
        irn=int(2d0*xrn)   !irn=0,1
        if(irn.gt.1) goto 16
        itmp=iflwrt(nwrt-1)
        iflwrt(nwrt-1)=itmp-sign(10-2*irn,itmp)
        if(itmp.gt.0) then
          icuwrt(1,nwrt-1)=nwrt-1
          icuwrt(2,nwrt)=nwrt-1
        else
          icuwrt(2,nwrt-1)=nwrt-1
          icuwrt(1,nwrt)=nwrt-1
        endif
        xrn=rangen2(1)
        itmp=iflwrt(nwrt-1)
        if(xrn.lt.scab2) then
          iflwrt(nwrt)=-icab(itmp,2)
        else
          iflwrt(nwrt)=-icab(itmp,1)
        endif
      elseif(idecmode.eq.6) then
        xrn=rangen2(1)
        if(xrn.lt.0.33333333333) then
          idecmode=4
        else
          idecmode=5
        endif
        goto 1
      endif
c     restore dec mode label
      idecmode=idec0
c
      m1=dmass(abs(iflwrt(nwrt-1)))
      m2=dmass(abs(iflwrt(nwrt)))
      pwrt(5,nwrt-1)=0
      pwrt(5,nwrt)=0
      tmp=pw(4)
      pdec(4)=tmp
      tmp=tmp*tmp
      do il=1,3
        pdec(il)=pw(il)
        tmp=tmp-pdec(il)**2
      enddo
      pdec(5)=sqrt(tmp)
c     rescale momenta to incorporate mass effects in case of W->c
C     and W->tau decays
      call rescms(pdec,pwrt(1,nwrt-1),pwrt(1,nwrt),m1,m2)
c
c$$$c     debug
c$$$      do i=nwrt-1,nwrt
c$$$        do il=1,16
c$$$          if(abs(iflwrt(i)).eq.il) fcount(il)=fcount(il)+1d0
c$$$        enddo
c$$$      enddo

      end


      subroutine phspace(lnot,pswgt,djpd,djg)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c  (njets+2)-body phase space for the process:                   c
c                                                                c
c  h(1)  h(2) -> h(3) ...  h(njets) e nu                         c
c                                                                c
c  where h is any light-flavoured quark or gluon                 c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'alpgen.inc'
      include 'wjet.inc'
      real *8 dummy,djg
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
      integer nx,ninit,i,l,m,lnot
      parameter (nx= 20)
      data ninit/0/
      save ninit
      if(ninit.eq.0) then
*
c-       setup local generation cuts
*
         cutkin(1)=ptjmin
         cutkin(3)=etajmax
         cutkin(5)=drjmin
         cutkin(7)=ptlmin
         cutkin(8)=etalmax
         cutkin(9)=drlmin
         cutkin(10)=metmin
c-
c-         vol0= (pi/2.d0)**5*s**4/120.d0/24.d0*(1.d0/(5.d0-ag))**2
c-         volm= 1.d0+(tau0**(5.d0-ag))*((5.d0-ag)*log(tau0)-1.d0)
c-         vol = vol0*volm 
c-
         ninit=1
      endif
*
c-    The generation starts
*
      pswgt=0.d0
*
      call momgen(njets,mw,wwid,roots,x1,x2,pcm,wgt,lnot)
      djg= 1.d0 ! dummy variable
      if (lnot.eq.1) then
         pswgt= 0.d0
         goto 100
      endif


*
c-    will write factor=factor0/(x1*x2), with factor0 function of njets etc.
*
      factor= 1d0/(2.d0*pi)**(3*njets+2)/2.d0/s/x1/x2
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
     .              (pcm(1,l)*pcm(1,m)+pcm(2,l)*pcm(2,m))
     .              /pt(l)/pt(m)
               if(dabs(dphi(l,m)).gt.1.d0) then
c$$$                  write(*,*) 'partons',l,m,', cos(Dphi)=', dphi(l,m)
c$$$     .                 ,', set to +-1'
c$$$                  write(*,*) 'pt(',l,')=',pt(l),'p(',l,')=',  (pcm(i,l)
c$$$     .                 ,i=0,3)
c$$$                  write(*,*) 'pt(',m,')=',pt(m),'p(',m,')=',  (pcm(i,m)
c$$$     .                 ,i=0,3)
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
      call chkcut(lnot,pt,p,eta,dr,njets)
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
c-    Sum of trasnverse mass mt^2 of jets
*
      do i=3,njets+2
         totpt=totpt+pt(i)**2
      enddo 
      if(iqopt.eq.0) then 
         qsq=1d0
      elseif(iqopt.eq.1) then
         qsq=mw**2+totpt
      elseif(iqopt.eq.2) then
         qsq=mw**2
      elseif(iqopt.eq.3) then
         qsq=mw**2+(p(1,njets+3)+p(1,njets+4))**2+(p(2,njets+3)+p(2
     $        ,njets+4))**2
      elseif(iqopt.eq.4) then
         qsq=totpt
      elseif(iqopt.eq.5) then
        totpt=0
        do i=3,njets+4
          totpt=totpt+pt(i)
        enddo 
         qsq=totpt**2
      endif
      qsq=qfac**2*qsq
 100  continue
      end

      subroutine chkcut(lnot,pt,p,eta,dr,njets)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     Applies kinematical cuts to the final state during the phase
c     -space generation
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 ptjmin,drjmin,etajmax
      real*8 ptlmin,metmin
      integer maxpar,ninit,njets,j,lnot,i
      real*8 cutkin(10)
      common/loccut/cutkin
      data ninit/0/
      parameter (maxpar=10)
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),p(5,maxpar)
      save ninit,ptjmin,etajmax,drjmin,ptlmin
     $     ,metmin 
      if(ninit.eq.0) then
         ninit=1
         ptjmin=cutkin(1)
         etajmax=cutkin(3)
         drjmin=cutkin(5)
         ptlmin=cutkin(7)
         metmin=cutkin(10)
      endif
      lnot= 0

c     impose minimum pt and require eta within allowed range, for jets
      do i=3,njets+2
         if (pt(i).lt.ptjmin)           goto 10
         if (abs(eta(i)).gt.etajmax)    goto 10
      enddo 
c     require dR(jet-jet)<drjmin
      do i=3,njets+1
         do j=i+1,njets+2
            if(dr(i,j).lt.drjmin)  goto 10
         enddo
      enddo
c
c     lepton cuts enforced in usrfll
c
      return
c     
 10   lnot= 1
      return
      end

      subroutine momgen(njets,rm,ga,roots,x1,x2,p,wgt,lw)
*
c-    Generator of np particles in the LAB frame.
*
      implicit none
      real*8 rm,ga,roots,x1,x2,wgt,zero,rm2,etacut
      real *8 wgt1, wgt2
      real*8 ptjmin,cutkin,etajmax,drjmin,s,sw
      real*8 ag,tau0,tau,xmsum,xmsum0,ptlim
      real*8 rootsh,sq,y0,yr,cxmw,cxpw,dj1,dj2,wtau,wjr
      real*8 pt0lmax,pt0lsum,en,pz,eta0min,ran0
      integer lw,lw1,npar,mpar,np,njets,lim,nw,j,k
      parameter (npar= 20)
      real*8 pm(0:4,npar),pt0(npar),pt1(npar),xm(npar),xmr(npar),
     .       eta0(npar),pt0_l(npar),pt1_l(npar),p(0:3,npar)
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
        lim    = 1
        zero   = 0.d0
        rm2    = rm*rm
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
c-      Parameters for the W:
*
        pt0_l(nw) = zero
        pt1_l(nw) = roots/2.d0
        eta0(nw)  = etacut
        xm(nw)    = zero  ! provisional W mass!!!
*
        cxmw      = zero
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
        cxpw= s
        call rans(ran0)
        call resonm(0,rm,ga,lim,cxmw,cxpw,sw,dj1,ran0)
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
        call dec2fm(0,sw,pm(0,nw),
     .              zero,zero,p(0,nw+2),p(0,nw+3),dj2,ranram(2))
*
        wgt= wjr/s/dj1/dj2/djbintot        
        return
      endif
*
c-    Generates the W Breit-Wigner:
*
      cxpw= s
      call rans(ran0)
      call resonm(0,rm,ga,lim,cxmw,cxpw,sw,dj1,ran0)
*
c-    The mass of the last particle and the upgrading of xmsum
*
      xm(nw)= sqrt(sw)
      xmsum = xmsum0+xm(nw)
*
c-    tau0 is the lower cut on tau= x1*x2
c-    In the following we allow an extra 10*Gamma(W) to accomodate BW tails
*
      tau0   = 1.d0/s*(xmsum)**2
      tau0   = max(tau0,4.d0*pt0lmax**2/s)  
      tau0   = max(tau0,pt0lsum**2/s)
      tau0   = max(tau0,1.d0/s*((rm-10.d0*ga)
     .                  +float(njets)*ptjmin)**2)
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
     .               +(xm(j)**2-(xmsum-xm(j))**2)**2
     .        -2.d0*s*(xm(j)**2+(xmsum-xm(j))**2)   
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
            p(k,j+2)= pm(k,j)*roots
          else
            pm(k,j) = pm(k,j)*roots
          endif
        enddo
      enddo
*
      wgt= wgt*s**(nw-2)     
*
      call dec2fm(0,sw,pm(0,nw),
     .            zero,zero,p(0,nw+2),p(0,nw+3),dj2,ranram(2*nw))
*
      wgt= wgt/dj1/dj2
*
      return
 100  lw= 1
      wgt= 0.d0
      return
      end
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
      parameter (npar= 20)
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
*
c-------------------------------------------------------------------
      subroutine alsgrd
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wjet.inc'
      integer nct1
      integer nvar,nch,n,nv
      integer init,j,k
      real*8 ni
      real*8 v1
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
      data v1/ncmax*0d0/,ni/maxn*0d0/
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
      nv  = 2*(njets+2)-1
      nvar= nv*ngrid
      do n= 2,nvar+1
        nch= 10
*
        nct(n)   =  nct(n-1)+nx1(n-1)   
        nx1(n)   =  nch 
*
c-      set the channels shared by mom and momr: 
*
        if (mod(n-1,nv).eq.1)    mask(n)= 1
        if (mod(n-1,nv).eq.2)    mask(n)= 1
        if (mod(n-1,nv).eq.nv-1) mask(n)= 1
        if (mod(n-1,nv).eq.0)    mask(n)= 1
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
         integer posi(2)
         integer hel(nmax)             !particles helicities
C
         integer flgdual          !dual (0) or su3 (1) amplitudes
         common/dual/flgdual
         integer color(2*nmax)    !coulor string
         integer colst(2*nmax)
         common/colore/color
         integer ncls,ant4(4,2),ant6(6,6),ant2(2,1)
         integer ant4eq(4,2)
         integer nant,nantl,j3
         parameter (ncls=2)       !different class of processes
         integer nx,proc(nmax),antns(6),ndual,j1,j2
         parameter (ndual=40320)
         real*8 dualamp(ndual),colfac,damp,dampref,avgspin
         integer colaux(ndual,2*nmax),ndl,ndla(0:10,1:3),class
         integer rep(nmax),nqrk,nprt,nglu,nlep,ngb,nphot,nh
         common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot,nh
         integer fq(nmax),pq(nmax)
         data ndla/0,0,1,5,23,119,719,5039,40319,0,0,
     >          0,1,5,23,119,719,5039,0,0,0,0,
     >          0,2,11,59,359,6*0/
         save ndla,ant2,ant4,ant6
         complex*16 result
         integer inpint(1000)
         common/initinter/inpint  
c
         integer jp
         common/jpr/jp
         integer cls(100)
         data cls/4*0           ! 2 quarks processes
     >        ,4*1              ! 4 quarks unequal flavours 
     >        ,7*2              ! 4 quarks equal flavours 
     >        ,2*1              ! 4 quarks unequal flavours 
     >        ,4*2              ! 4 quarks equal flavours 
     >        ,1                ! 4 quarks unequal flavours 
     >        ,2                ! 4 quarks equal flavours 
     >        ,77*-99/
C
         data ant2  /1,2/
         data ant4  /1,4,2,3,
     >               1,3,2,4/
         data ant6  /1,5,2,6,3,4,
     >               1,6,2,4,3,5,
     >               1,6,2,5,3,4,
     >               1,4,2,6,3,5,
     >               1,4,2,5,3,6,
     >               1,5,2,4,3,6/         
C
c
c         data inpint/ 6,                                                ! N V-A 
c     >                1,1,1,1,1,   1,1,2,1,2,   2,1,1,1,2,   3,1,2,1,1, ! zuu, zdd, w+ud, w-ud,
c     >               10,1,1,1,1,  10,1,2,1,2,                           ! Auu, Add
c     >                0,                                                ! N of yukawa
c     >               13,                                                ! N self-gauge
c     >                1, 2, 3,  10, 2, 3,   4, 2, 3,   5, 1, 1,         ! ZWW, AWW, auxiliary
c     >                5, 1,10,   5,10,10,   5, 2, 3,   6, 1, 2,         ! auxiliary
c     >                6,10, 2,   7, 1, 3,   7,10, 3,   8, 2, 2,         ! auxiliary
c     >                9, 3, 3,                                          ! auxiliary
c     >                2,                                                ! N H-GAUGE
c     >                1, 1, 1,   1, 2, 3,                               ! HZZ, HWW
c     >                921*-100       /                                 ! EOF signal
c
         data inpint/ 7,
     >               11,1,1,1,1,  11,1,2,1,2,  11,3,2,3,2,   2,1,1,1,2,  ! guu, gdd, gbb, w+ud
     >                3,1,4,1,3,  11,2,1,2,1,  11,3,1,3,1,               ! w-en, gcc
     >                0,                                                 ! N of yukawa
     >                2,                                                 ! N self-gauge
     >                11,11,11,  12,11,11,                               ! ggg Auxgg
     >                956*-100/
c
         if (nqrk.eq.4.or.nqrk.eq.2) then
         else
          write(*,*)'nqrk=',nqrk ,'not yet supported'
         endif
         avgspin=0.d0
c
         if (labcol.eq.0) then
C
          flgdual=1
          call matrix0(flvmlm,posi,impul,hel,result)
c
         elseif (labcol.eq.1) then
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
             antns(j1)=pq(ant6(j1,j3))
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


