c-------------------------------------------------------------------
      subroutine alsprc
c     assigns the hard process code
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      ihrd=1
      end
c-------------------------------------------------------------------
      subroutine alhset
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wqq.inc'
      integer i,jprocmax,jproc,ngrid,jgrid,maxsubproc
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
      character*2 Q(4:6)
      character*5 Qbar(4:6)
      data Q/' c',' b',' t'/,Qbar/' cbar',' bbar',' tbar'/
c
c     parameters for the gauge invariance prescription:
      winsize  = 2.d0/pi
      resonance= 'n'
      wmode    = 'yy'
c     
c     process input parameters
c
      nfspart=njets+2
      npart=nfspart+4
      nprtns=nfspart
      nw=1
      if(ihvy.eq.4) then
        mq=mc
        ptQmin=ptcmin
        etaQmax=etacmax
        drQmin=drcmin
        paruse(21,1)=0
        paruse(31,1)=0
        paruse(41,1)=0
        paruse(51,1)=0
      elseif(ihvy.eq.5) then
        mq=mb
        ptQmin=ptbmin
        etaQmax=etabmax
        drQmin=drbmin
        paruse(22,1)=0
        paruse(32,1)=0
        paruse(42,1)=0
        paruse(52,1)=0
      elseif(ihvy.eq.6) then
        mq=mt
        ptQmin=0
        etaQmax=100
        drQmin=0
        paruse(21,1)=0
        paruse(31,1)=0
        paruse(41,1)=0
        paruse(51,1)=0
        paruse(22,1)=0
        paruse(32,1)=0
        paruse(42,1)=0
        paruse(52,1)=0
      endif
c     set masses
      do i=1,npart
         p(5,i)=0
      enddo
      p(5,3)=mq
      p(5,4)=mq
c
      if(nfspart.eq.2) then
         jprocmax=1
         ngrid   = 1
         navg=1
      elseif(nfspart.eq.3) then
         jprocmax=3
         ngrid=3
         navg=1
      elseif(nfspart.eq.4) then
         jprocmax=15
         ngrid=5
         navg=1
      elseif(nfspart.eq.5) then
         jprocmax=23
         ngrid=5
         navg=1
      elseif(nfspart.eq.6) then
         jprocmax=25
         ngrid=5
         navg=1
      endif 
c     
      if(nfspart.eq.2) then
        write(niosta,*) 'W',Q(ihvy),Qbar(ihvy)
      else
        write(niosta,*) 'W',Q(ihvy),Qbar(ihvy),
     $       ' +',nfspart-2,' jets'
      endif
      write(niosta,*) 'W-> ell nu'
      write(niosta,*) '======================================='
      write(niosta,*) Q(ihvy),' mass:',mq
      write(niosta,*)
     $     'Generation cuts for the partonic event sample:' 
      if(nfspart.gt.2) then
        write(niosta,*) '     Light jets:'
        write(niosta,*) 'ptmin=',ptjmin,' |etamax|=',etajmax
     $       ,' dR(j-j),dR(Q-j)>',drjmin 
      endif
      write(niosta,*) Q(ihvy),' quarks:'
      if(ihvy.eq.6) then
        write(niosta,*) 'no generation cuts for top'
      else
        write(niosta,*) 'ptmin=',ptQmin,' |etamax|=',etaQmax
     $       ,' dR(QQ)>',drQmin 
      endif
      write(niosta,*) '     Leptons:'
      write(niosta,*) 'ptmin(lep)=',ptlmin,' |etamax|=',etalmax
     $     ,' Et(miss)>',metmin,' dR(l-j)>',drlmin  
      write(niosta,*) 'pthrmin=',pthrmin,' pthrmax=',pthrmax
      end
                              
c-------------------------------------------------------------------
      subroutine selflav(jproc,xlum,afl)
c was:     subroutine selflav(jproc,npart,as,f1,f2,xlum,ifl,afl)
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
c     22 q g -> W q' q' qbar' ,   qbar g -> W qbar' q' qbar'
c     23 g q -> W q' q' qbar' ,   g qbar -> W qbar' q' qbar'
c     24 gg -> q qbar' q'' qbar'' W              
c     25 gg -> q qbar q qbar' W  + c.c.
c-----------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wqq.inc'
c commons
      integer icconj
      common/hwconv/icconj
c
c      double precision as,pi
c      parameter (pi=3.14159265358979324d0)
c      integer maxpar
c      parameter (maxpar=10)
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
c     ccoef(i,j) for nfspart=j and i=#(light quark pairs)
      real*8 ccoef(3,0:4)    !ccoeff(number of q-pb pairs, number of gluons)
      data ccoef/-1.d0,0.185185185,0.127572016,-1.d0,0.12037037,
     >            0.0936213992,-1.d0,0.0902777778,0.0743312757,
     >           -1.d0,0.072337963,-1.d0,-1.d0,0.0603841146,-1.d0/
c old:
c      real *8 ccoef(2,6)
c      data ccoef/ 2*0, 0.185185185,0.,0.12037037,0.,0.0902777778,0
c     $ .127572016,0.072337963,0.0936213992,0.0603841146,0.0743312757/
c
c     $  2*0, 15/81,0d0, 78/648,0d0, 468/5184, 93/729, 
c     $  3000/41472, 546/5832, 20034/331776, 3468/46656 /
c
      integer ib,ibb,il,inu,iflag
      save init,chrg,cfac,ifact,iqpp,icab,fcab
c
      if(init.eq.0) then
         fcab(1)=real(ccab2)
         fcab(2)=real(scab2)
         nfspart=npart-4
         init=1
         ifact(0)=1e0
         do i=1,6
            ifact(i)=ifact(i-1)/real(i)
         enddo
      endif
c
      do i=1,npart
         ifl(i)=0
         afl(i)=0
      enddo 
c
 1    call randa(xrn)
      rn=real(xrn)
      if(1e0-rn.lt.1e-7) goto 1

      icount=0
      slum=0e0
      tmptot=0e0
c
c     start two light-quark processes
c
c     jproc=1   nfspart>=2
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
                     do k=5,nfspart+2
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
c     jproc=2   nfspart>=3
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
                     ifl(5)=icab(i,j)
                     do k=6,nfspart+2
                        ifl(k)=0
                     enddo 
c     g d -> u W-
                     itmp=imap(ifl(2))
                     afl(2)=imap(ifl(2))
                     afl(5)=imap(ifl(5))
                     goto 100
                  endif 
               enddo
            endif
         enddo 
c
c     jproc=3   nfspart>=3
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
                     ifl(5)=icab(i,j)
                     do k=6,nfspart+2
                        ifl(k)=0
                     enddo 
c     d g -> u W-
                     itmp=imap(ifl(1))
                     afl(1)=imap(ifl(1))
                     afl(5)=imap(ifl(5))
                     goto 100
                  endif 
               enddo
            endif
         enddo 
c     jproc=4   nfspart>=4
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
                  ifl(5)=i
                  ifl(6)=-icab(i,j)
                  do k=7,nfspart+2
                     ifl(k)=0
                  enddo 
c     gg -> u dbar W-
                  itmp=-imap(ifl(6))
                  afl(5)= imap(ifl(5))
                  afl(6)= imap(ifl(6))
                  goto 100
               endif
            enddo
         enddo 
c
c     start 4 light quark processes  (in order of increasing number of
c     minimum nfspart, different flavour first, like-flavours after
c
c     jproc=5  1001  nfspart>=4
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
                        ifl(5)=iqpp(k,j,abs(i))
                        ifl(6)=-ifl(5)
                        do l=7,nfspart+2
                           ifl(l)=0
                        enddo 
c     d ubar -> W- c cbar
                        itmp=imap(ifl(1))
                        afl(1)=imap(ifl(1))
                        afl(2)=imap(ifl(2))
                        afl(5)=sign(4,ifl(5))
                        afl(6)=-afl(5)
                        goto 100
                     endif
                  enddo
               enddo 
            endif
         enddo
c
c     jproc=6 1002  nfspart>=4
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
                        ifl(5)=icab(i,j)
                        ifl(6)=ifl(2)
                        do l=7,nfspart+2
                           ifl(l)=0
                        enddo 
c     d c(bar) -> W- u c(bar)
                        itmp=imap(ifl(1))
                        afl(1)=imap(ifl(1))
                        afl(5)=imap(ifl(5))
                        afl(2)=sign(4,ifl(2))
                        afl(6)=afl(2)
                        goto 100
                     endif 
                  enddo
               enddo
            endif
         enddo 
c
c     jproc=7 1002 nfspart>=4
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
                        ifl(5)=icab(i,j)
                        ifl(6)=ifl(1)
                        do l=7,nfspart+2
                           ifl(l)=0
                        enddo 
c      c(bar) d -> W u c(bar)
                        itmp=imap(ifl(2))
                        afl(2)=imap(ifl(2))
                        afl(5)=imap(ifl(5))
                        afl(1)=sign(4,ifl(1))
                        afl(6)=afl(1)
                        goto 100
                     endif 
                  enddo
               enddo
            endif
         enddo 
c
c     jproc=8 nfspart>=4
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
                              ifl(5)=k
                              ifl(6)=-icab(k,j)
                              do l=7,nfspart+2
                                 ifl(l)=0
                              enddo 
c      c cbar -> W- u dbar 
                              itmp=-imap(ifl(6))
                              afl(1)=sign(4,ifl(1))
                              afl(2)=-afl(1)
                              afl(5)=imap(ifl(5))
                              afl(6)=imap(ifl(6))
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
c     jproc=9  nfspart>=4
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
                     ifl(5)=i
                     ifl(6)=-ifl(5)
                     do l=7,nfspart+2
                        ifl(l)=0
                     enddo 
c     d ubar -> W- d dbar
                     itmp=imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     goto 100
                  endif 
               enddo
            endif
         enddo
c     jproc=10 nfspart>=4
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
                     ifl(5)=i
                     ifl(6)=-ifl(5)
                     do l=7,nfspart+2
                        ifl(l)=0
                     enddo 
c     ubar d -> W- d dbar
                     itmp=imap(ifl(2))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     goto 100
                  endif 
               enddo 
            endif
         enddo
c
c     jproc=11 nfspart>=4
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
                     ifl(5)=i
                     ifl(6)=-icab(i,j)
                     do l=7,nfspart+2
                        ifl(l)=0
                     enddo 
c     dbar d -> W- dbar u 
                     itmp=imap(ifl(2))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
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
                     ifl(5)=icab(i,j)
                     ifl(6)=-i
                     do l=7,nfspart+2
                        ifl(l)=0
                     enddo 
c     d dbar -> W- u dbar
                     itmp=imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     goto 100
                  endif
               enddo 
            endif
         enddo
c
c     jproc=13  nfspart>=4
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
                     ifl(5)=i
                     ifl(6)=icab(i,j)
                     do l=7,nfspart+2
                        ifl(l)=0
                     enddo 
c     d d -> W- d u
                     itmp=imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     goto 100
                  endif
               enddo 
            endif
         enddo
c
c     jproc=14 nfspart>=4
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
                     ifl(5)=i
                     ifl(6)=i
                     do l=7,nfspart+2
                        ifl(l)=0
                     enddo 
c     dbar ubar -> W- dbar dbar
                     itmp=-imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     goto 100
                  endif
               enddo 
            endif
         enddo
c
c     jproc=15 nfspart>=4
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
                     ifl(5)=ifl(2)
                     ifl(6)=ifl(2)
                     do l=7,nfspart+2
                        ifl(l)=0
                     enddo 
c     ubar dbar -> W dbar dbar 
                     itmp=-imap(ifl(2))
                     afl(1)= imap(ifl(1))
                     afl(2)= imap(ifl(2))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     goto 100
                  endif
               enddo 
            endif
         enddo
c
c     jproc=16 nfspart>=5
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
                        ifl(5)=icab(i,j)
                        ifl(6)=iqpp(k,j,abs(i))
                        ifl(7)=-ifl(6)
                        do l=8,nfspart+2
                           ifl(l)=0
                        enddo 
c     d g -> W- u c cbar
                        itmp=imap(ifl(1))
                        afl(1)=imap(ifl(1))
                        afl(5)=imap(ifl(5))
                        afl(6)=sign(4,ifl(6))
                        afl(7)=-afl(6)
                        goto 100
                     endif 
                  enddo
               enddo
            endif
         enddo
c
c     jproc=17  nfspart>=5
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
                        ifl(5)=icab(i,j)
                        ifl(6)=iqpp(k,j,abs(i))
                        ifl(7)=-ifl(6)
                        do l=8,nfspart+2
                           ifl(l)=0
                        enddo 
c     g d -> W- u c cbar
                        itmp=imap(ifl(2))
                        afl(2)=imap(ifl(2))
                        afl(5)=imap(ifl(5))
                        afl(6)=sign(4,ifl(6))
                        afl(7)=-afl(6)
                        goto 100
                     endif 
                  enddo
               enddo 
            endif
         enddo
c
c     jproc=18  nfspart>=5
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
                     ifl(5)=i
                     ifl(6)=i
                     ifl(7)=-icab(i,j)
                     do l=8,nfspart+2
                        ifl(l)=0
                     enddo 
c     dbar g -> W dbar dbar u
                     itmp=-imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     afl(7)= imap(ifl(7))
                     goto 100
                  endif 
               enddo
            endif
         enddo
c
c     jproc=19  nfspart>=5
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
                     ifl(5)=icab(i,j)
                     ifl(6)=i
                     ifl(7)=-i
                     do l=8,nfspart+2
                        ifl(l)=0
                     enddo 
c     d g -> W u d dbar
                     itmp=imap(ifl(1))
                     afl(1)= imap(ifl(1))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     afl(7)= imap(ifl(7))
                     goto 100
                  endif 
               enddo
            endif
         enddo
c
c     jproc=20  nfspart>=5
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
                     ifl(5)=i
                     ifl(6)=i
                     ifl(7)=-icab(i,j)
                     do l=8,nfspart+2
                        ifl(l)=0
                     enddo 
c     g dbar -> W dbar dbar u
                     itmp=-imap(ifl(2))
                     afl(2)= imap(ifl(2))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     afl(7)= imap(ifl(7))
                     goto 100
                  endif 
               enddo
            endif
         enddo
c
c     jproc=21  nfspart>=5
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
                     ifl(5)=icab(i,j)
                     ifl(6)=i
                     ifl(7)=-i
                     do l=8,nfspart+2
                        ifl(l)=0
                     enddo 
c     g d -> W u d dbar
                     itmp=imap(ifl(2))
                     afl(2)= imap(ifl(2))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     afl(7)= imap(ifl(7))
                     goto 100
                  endif 
               enddo
            endif
         enddo
c
c     jproc=22 nfspart>=5
c     q g -> W q' q' qbar' ,   qbar g -> W qbar' q' qbar'
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
                     ifl(5)=icab(i,j)
                     ifl(6)=ifl(5)
                     ifl(7)=-ifl(6)
                     do l=8,nfspart+2
                        ifl(l)=0
                     enddo 
c     d g -> W- u c cbar
                     itmp=imap(ifl(1))
                     afl(1)=imap(ifl(1))
                     afl(5)=imap(ifl(5))
                     afl(6)=afl(5)
                     afl(7)=-afl(6)
                     goto 100
                  endif 
               enddo
            endif
         enddo
c
c     jproc=23  nfspart>=5
c     g q -> W q' q' qbar' ,   g qbar -> W qbar' q' qbar'
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
                     ifl(5)=icab(i,j)
                     ifl(6)=ifl(5)
                     ifl(7)=-ifl(6)
                     do l=8,nfspart+2
                        ifl(l)=0
                     enddo 
c     g d -> W- u c cbar
                     itmp=imap(ifl(2))
                     afl(2)=imap(ifl(2))
                     afl(5)=imap(ifl(5))
                     afl(6)=afl(5)
                     afl(7)=-afl(6)
                     goto 100
                  endif 
               enddo 
            endif
         enddo
c
c     jproc=24  nfspart>=6
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
                     ifl(5)=i
                     ifl(6)=-icab(i,j)
                     ifl(7)=iqpp(k,j,abs(i))
                     ifl(8)=-ifl(7)
                     do l=9,nfspart+2
                        ifl(l)=0
                     enddo 
c     gg -> u dbar c cbar  W-              
                     itmp=-imap(ifl(6))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     afl(7)= sign(4,ifl(7))
                     afl(8)=-afl(7)
                     goto 100
                  endif 
               enddo
            enddo
         enddo
c
c     jproc=25  nfspart>=6
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
                     ifl(5)=i
                     ifl(6)=-i
                     ifl(7)=i
                     ifl(8)=-icab(i,j)
                     do k=9,nfspart+2
                        ifl(k)=0
                     enddo 
c     gg -> dbar d dbar u W-  + c.c.
                     itmp=-imap(ifl(5))
                     afl(5)= imap(ifl(5))
                     afl(6)= imap(ifl(6))
                     afl(7)= imap(ifl(7))
                     afl(8)= imap(ifl(8))
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
      il1=nfspart+3
      il2=nfspart+4
c
c     label type of flavour configuration:
c
c     d ubar -> b bbar e- nubar:
c     input reference b bbar and leptons
      ifl(il1)=11
      ifl(il2)=-12
      ib=3
      ibb=4
      il=il1
      inu=il2
c     u dbar -> b bbar nu e+    ISOSPIN flip:
      if(itmp.eq.2) then
         ifl(il1)= 12
         ifl(il2)=-11
         il=il2
         inu=il1
c     dbar u -> bbar b e+ nu    CHARGE conjugation:
      elseif(itmp.eq.-1) then
         ifl(il1)=-11
         ifl(il2)= 12
         ib=4
         ibb=3
c     ubar d -> bbar b nubar e-   ISOSPIN flip and CHARGE conjugation
      elseif(itmp.eq.-2) then
         ifl(il1)=-12
         ifl(il2)= 11
         ib=4
         ibb=3
         il=il2
         inu=il1
      endif 
c
      ifl(ib)=ihvy
      ifl(ibb)=-ihvy
      afl(il1)=ifl(il1)
      afl(il2)=ifl(il2)
      afl(ib)=ifl(ib)
      afl(ibb)=ifl(ibb)
c     check lepton charges
      wchrg=chrg(ifl(1))+chrg(ifl(2))
      do i=5,nfspart+2
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
      call usrfll(ib,ibb,il,inu,iflag)
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
      do i=3,nfspart+2
         cwgt=cwgt*cfac(abs(ifl(i)))
         if(ifl(i).eq.0) then
            ifl(i)=21
            ng=ng+1
         endif
      enddo
c     evaluate spin weight factors
      swgt=2e0
      swgt=swgt**(nfspart+2)
c
      xlum=dble(slum*cwgt*swgt*ifact(ng))/resc**nfspart
      if(jproc.eq.14.or.jproc.eq.25.or.jproc.eq.20.or.jproc.eq.18.or
     $   .jproc.eq.15.or.jproc.eq.22.or.jproc.eq.23)  xlum=0.5d0*xlum ! identical quarks in final state
c     effco=ccoef(#quark pairs,#gluons)
      if(jproc.le.4) then
         nlqp=1
      elseif(jproc.le.25) then
         nlqp=2
      else
         write(*,*) 'jproc not valid'
         stop
      endif 
      xlum=xlum*ccoef(nlqp+1,nfspart-2*nlqp)
c
      end
c     end selflav

      subroutine usrfll(ib,ibb,il,inu,iflag)
      implicit none
      include 'alpgen.inc'
      include 'wqq.inc'
      integer i,j,il,inu,ib,ibb,iflag
c
      iflag=0
c assign momenta to the usr common block
      do i=1,4
         do j=1,2
            pin(i,j)=p(i,j)
         enddo
         do j=1,npart-2
            pout(i,j)=p(i,j+2)
         enddo
         do j=1,nfspart
            pjet(i,j)=p(i,j+2)
         enddo
      enddo
      do j=1,nfspart
         ptj(j)=pt(j+2)
         etaj(j)=eta(j+2)
         do i=j+1,nfspart
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
      do i=1,nfspart
         drbj(i)=dr(ib,2+i)
         drbbj(i)=dr(ibb,2+i)
      enddo
      do i=1,4
         plep(i)=p(i,il)
         pnu(i)=p(i,inu)
         pw(i)=plep(i)+pnu(i)
      enddo
      ptlep=pt(il)
      ptmiss=pt(inu)
      etalep=eta(il)
      do i=1,nfspart
         drlj(i)=dr(il,2+i)
      enddo
      drlb=dr(il,ib)
      drlbb=dr(il,ibb)
c
c     impose minimum pt and require eta within allowed range, for leptons
      if (ptlep.lt.ptlmin)           goto 10
      if (ptmiss.lt.metmin)           goto 10
      if (abs(etalep).gt.etalmax)    goto 10
c
c     isolate lepton from b's and light jets
      if(ihvy.ne.6) then
        if(drlb.lt.drlmin)        goto 10
        if(drlbb.lt.drlmin)        goto 10
      endif
      do i=1,nfspart
         if(drlj(i).lt.drlmin)        goto 10
      enddo 
c
      ptbjet=0d0
      return
 10   iflag=1
      end

      subroutine setdec(nwrt,iflwrt,icuwrt,pwrt,decbr)
      implicit none
      include 'alpgen.inc'
      include 'wqq.inc'
c locals
      integer maxdec
      parameter (maxdec=40)
      integer ip, ic, il, irn,i,itmp,idecmode
     $     ,idec0
      integer init
      data init/0/
      integer icab(-4:4,2) 
      double precision xrn,rangen2,tmp
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
      data idbg/0/
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
      tmp=pw(4)
      pdec(4)=tmp
      tmp=tmp*tmp
      do il=1,3
        pdec(il)=pw(il)
        tmp=tmp-pdec(il)**2
      enddo
      pdec(5)=sqrt(tmp)
c     rescale momenta to incorporate mass effects in case of W->c
C     decays
      m1=dmass(abs(iflwrt(nwrt-1)))
      m2=dmass(abs(iflwrt(nwrt)))
      call rescms(pdec,pwrt(1,nwrt-1),pwrt(1,nwrt),m1,m2)
c
c     debug
c      do i=nwrt-1,nwrt
c        do il=1,16
c          if(abs(iflwrt(i)).eq.il) fcount(il)=fcount(il)+1d0
c        enddo
c      enddo

      end
*
      subroutine phspace(lnot,pswgt,djpd,djg)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c  (nfspart+2)-body phase space for the process:                   c
c                                                                c
c  h(1)  h(2) -> b(3) anti_b(4) h(5) ...  h(nfspart) e nu          c
c                                                                c
c  where h is any light-flavoured quark or gluon                 c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'alpgen.inc'
      include 'wqq.inc'
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
         cutkin(2)=ptQmin
         cutkin(3)=etajmax
         cutkin(4)=etaQmax
         cutkin(5)=drjmin
         cutkin(6)=drQmin
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
      call momgen(nfspart,mw,wwid,mq,roots,x1,x2,pcm,wgt,lnot)
      djg= 1.d0 ! dummy variable
      if (lnot.eq.1) then
         pswgt= 0.d0
         goto 100
      endif
* 
c-    will write factor=factor0/(x1*x2), with factor0 function of nfspart etc.
*
      factor= 1d0/(2.d0*pi)**(3*nfspart+2)/2.d0/s/x1/x2
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

      do l= 3,nfspart+4
          pt(l) = sqrt(pcm(1,l)**2+pcm(2,l)**2)
          pl(l) = pcm(3,l)
          eta(l)= -log(tan(0.5d0*atan2(pt(l),pl(l))))
          y(l)  = 0.5d0*log((pcm(0,l)+pcm(3,l))/(pcm(0,l)-pcm(3,l)))
      enddo
*
c-    Calculates jet-jet distances:
*
      do l= 3,nfspart+3
         do m= l+1,nfspart+4
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
      do l= 3,nfspart+3
         do m= l+1,nfspart+4
            dr(l,m)= sqrt(dphi(l,m)**2+(eta(l)-eta(m))**2)
            dr(m,l)=dr(l,m)
         enddo
      enddo
*
c-    Redefine the momenta:           
*                                     
      do l= 1,nfspart+4                 
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
      call chkcut(lnot,pt,p,eta,dr,nfspart,ihvy)
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
c-    Sum of transverse masses of jets
*
      do i=3,nfspart+2
         totpt=totpt+pt(i)**2+p(5,i)**2
      enddo 
      if(iqopt.eq.0) then 
         qsq=1d0
      elseif(iqopt.eq.1) then
         qsq=mw**2+totpt
      elseif(iqopt.eq.2) then
         qsq=mw**2
      elseif(iqopt.eq.3) then
         qsq=mw**2+(p(1,nfspart+3)+p(1,nfspart+4))**2+(p(2,nfspart+3)
     $        +p(2,nfspart+4))**2
      elseif(iqopt.eq.4) then
         qsq=totpt
      elseif(iqopt.eq.5) then
        totpt=0
        do i=3,nfspart+4
          totpt=totpt+pt(i)
        enddo 
         qsq=totpt**2
      endif
      qsq=qfac**2*qsq
 100  continue
      end
*     
      subroutine chkcut(lnot,pt,p,eta,dr,nfspart,ihvy)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     c
c     Applies generic kinematical cuts to the final state 
c     during the phase-space generation
c     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 ptjmin,drQmin,drjmin,etaQmax,etajmax,ptQmin
      real*8 ptlmin,metmin
      integer i,j,maxpar,ninit,nfspart,lnot,ihvy
      real*8 cutkin(10)
      common/loccut/cutkin
      data ninit/0/
      parameter (maxpar=10)
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),p(5,maxpar)
      save ninit,ptjmin,ptQmin,etajmax,etaQmax,drjmin,drQmin,ptlmin
     $     ,metmin 
      if(ninit.eq.0) then
        ninit=1
        ptjmin=cutkin(1)
        ptQmin=cutkin(2)
        etajmax=cutkin(3)
        etaQmax=cutkin(4)
        drjmin=cutkin(5)
        drQmin=cutkin(6)
        ptlmin=cutkin(7)
        metmin=cutkin(10)
        ptlmin=min(ptlmin,metmin)
      endif

      lnot= 0

c     impose minimum pt and require eta within allowed range, for b's
      if(ihvy.ne.6) then
        do i=3,4
          if (pt(i).lt.ptQmin)           goto 10
          if (abs(eta(i)).gt.etaQmax)    goto 10
        enddo 
c     require dR(b-bbar)<drQmin
        if(dr(3,4).lt.drQmin)           goto 10
      endif
c     
      if(nfspart.gt.2) then
c     impose minimum pt and require eta within allowed range, for jets
        do i=5,nfspart+2
          if (pt(i).lt.ptjmin)           goto 10
          if (abs(eta(i)).gt.etajmax)    goto 10
        enddo 
c     require dR(b-jet)<drjmin, dR(jet-jet)<drjmin
        do i=5,nfspart+2
          if(ihvy.ne.6) then
            if(dr(3,i).lt.drjmin)        goto 10
            if(dr(4,i).lt.drjmin)        goto 10
          endif
          do j=i+1,nfspart+2
            if(dr(i,j).lt.drjmin)  goto 10
          enddo
        enddo
      endif 
c     
c     lepton cuts enforced in usrfll
c     
      return
 10   lnot= 1
      return
      end
*
*
      subroutine momgen(nfspart,rm,ga,mQ,roots,x1,x2,p,wgt,lw)
*
c-    Generator of np particles in the LAB frame.
*
      implicit none
      real *8 wgt1, wgt2
      real*8 rm,ga,mQ,roots,x1,x2,wgt,zero,rm2,mQ2,etacut
      real*8 ptjmin,ptQmin,cutkin,etajmax,drmQin,drjmin,s,sw
      real*8 ag,tau0,tau,xmsum,xmsum0,ptlim
      real*8 rootsh,sq,y0,yr,cxmb,cxmw,cxpw,dj1,dj2,wtau,wjr
      real*8 etamQax,pt0lmax,pt0lsum,en,pz,eta0min,ran0
      integer lw,lw1,npar,mpar,np,nfspart,lim,nw,j,k
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
        np     = nfspart+2
        nw     = np-1
        lim    = 1
        zero   = 0.d0
        rm2    = rm*rm
        mQ2    = mQ*mQ
        etacut = 40.d0
        s      = roots*roots
*
        ptjmin = cutkin(1)
        ptQmin = cutkin(2)
        etajmax= cutkin(3)
        etamQax= cutkin(4)
        drjmin = cutkin(5)
        drmQin = cutkin(6)
*
c-      Parameters for b and bbar:
*
        do j= 1,2
          pt0_l(j) = ptQmin
          pt1_l(j) = roots/2.d0
          eta0(j)  = etamQax
          xm(j)    = mQ
        enddo
*
        cxmb = (2.d0*ptQmin*sin(drmQin/2.d0))**2
        cxmb = max(4.d0*mQ2,cxmb)
*
c-      Parameters for the light jets:
*
        do j= 3,nw-1
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
c-    Generates the W Breit-Wigner:
*
      cxpw= (roots-sqrt(cxmb))**2
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
      tau0   = max(tau0,1.d0/s*((rm-10.d0*ga)+2.d0*sqrt(ptQmin**2+mQ
     .       **2)+float(nfspart-2)*ptjmin)**2)
      tau0   = max(tau0,cxmb/s)
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
c-------------------------------------------------------------------
      subroutine alsgrd
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'wqq.inc'
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
      nv  = 2*(nfspart+2)-1
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
*
c-------------------------------------------------------------------
      subroutine dumpwgt
c-------------------------------------------------------------------
c     routine for the debugging of individual large-weight events
c     Not functional to the rnning of the code
      implicit none
      include 'alpgen.inc'
      include 'wqq.inc'
      real *8 tmpmw,ptmp(4)
      integer i
c      write(*,*) 'x1,x2=',x1,x2
c      write(*,*) 'p(b)=',(p(i,3),i=1,3)
c      write(*,*) 'p(bbar)=',(p(i,4),i=1,3)
c      write(*,*) 'p(j1)=',(p(i,5),i=1,3)
c      write(*,*) 'p(j2)=',(p(i,6),i=1,3)
c      write(*,*) 'p(e)=',(p(i,7),i=1,3)
c      write(*,*) 'p(nu)=',(p(i,8),i=1,3)
      write(niosta,*) 'x1,x2=',x1,x2
      write(niosta,*) 'p(b)=',(p(i,3),i=1,3)
      write(niosta,*) 'p(bbar)=',(p(i,4),i=1,3)
      write(niosta,*) 'p(j1)=',(p(i,5),i=1,3)
      write(niosta,*) 'p(j2)=',(p(i,6),i=1,3)
      write(niosta,*) 'p(e)=',(p(i,7),i=1,3)
      write(niosta,*) 'p(nu)=',(p(i,8),i=1,3)
c
      pw(4)=p(4,7)+p(4,8)
      tmpmw=pw(4)**2
      do i=1,3
         pw(i)=p(i,7)+p(i,8)
         tmpmw=tmpmw-pw(i)**2
      enddo
      tmpmw=sqrt(tmpmw)
c      write(*,*) 'm(e-nu)=',tmpmw
      write(niosta,*) 'm(e-nu)=',tmpmw
c
c b-bbar inv mass
      ptmp(4)=p(4,3)+p(4,4)
      tmpmw=ptmp(4)**2
      do i=1,3
         ptmp(i)=p(i,3)+p(i,4)
         tmpmw=tmpmw-ptmp(i)**2
      enddo
      tmpmw=sqrt(tmpmw)
c      write(*,*) 'm(b-bbar)=',tmpmw
      write(niosta,*) 'm(b-bbar)=',tmpmw
c
c j-j inv mass
      ptmp(4)=p(4,5)+p(4,6)
      tmpmw=ptmp(4)**2
      do i=1,3
         ptmp(i)=p(i,5)+p(i,6)
         tmpmw=tmpmw-ptmp(i)**2
      enddo
      tmpmw=sqrt(tmpmw)
c      write(*,*) 'm(jj)=',tmpmw
      write(niosta,*) 'm(jj)=',tmpmw
      end


c     
c*****alpha-specific routines
C***********************************************************************
      subroutine matrix(flvmlm,posi,impul,hel,rst,labcol
     >     ,wgcol,colfl)
C***********************************************************************
C     
C     Input parameters: FLVMLM(NMAX) particles according
C     to michelangelo convention, POSI(NMAX) positions of
C     incoming particles IMPUL(4,8) particles four momenta,
C     HEL, array of particle helicities 
C     LABCOL to choose whether su(3) mode only (LABCOL=0) or
C     color flow is required (LABCOL=1) or dual mode only
C     (LABCOL=2), WGCOL random number for color flow unweighting
C     Output parameters:   RESULT squared matrix element (modulus),
c     COLFLOW 
C     string representing the selected coulor flow 
C     
C     
C     At present:
C     
C     the color flow is returned in the following way: each particle is
C     represented by a pair of integers (n1,n2); two particle (n1,n2),
c     (n3,n4)
C     are colour connected if either n2=n3 or n4=n1 (in this case the
c     coulor
C     flows from 2 to 1). COLFLOW will contain a string of pairs of
c     integers 
C     ordered as follows: d (or dbar) nu_ebar bbar ubar (or u) e^- b glu
c     glu
C     (with gluons ordered according to the ordering of momenta).
C     (COLFLOW=(n1,....,nn,....) only the first nn=2*particles_number
C     elements to be used)   
C     
C     
C     IMPUL(J,NPART) ===   J=1 Energy of the NPART-th particle
C     IMPUL(J,NPART) ===   J=2,3,4 x,y,z components  of the NPART-th
c     particle 
C     three momentum
C     the order of the momenta is: 
C     dbar nubar bbar u e- b ( glu glu glu)  (processes 1,2,3,4)
C     dbar nubar cbar bbar u e- c b (glu)  (processes > 5 )
C     where all particles are assumed outcoming and incoming gluons
C     must be before outcoming ones
C     
      implicit none
C     
      integer nmax              !maximum number of external particles, 
      parameter (nmax=10)
      integer nlb
      parameter (nlb=4)    
      real*8 impul(4,nmax),wgcol !the contribution to the integral
      real*8 rst
c     complex*16 rst
      integer labcol,colfl(2*nmax),flvmlm(nmax)
      integer posi(2)
      integer hel(nmax)         !particles helicities
C     
      integer flgdual           !dual (0) or su3 (1) amplitudes
      common/dual/flgdual
      integer color(2*nmax)     !coulor string
      integer colst(2*nmax)
      common/colore/color
      integer ncls,nprc,ant4(4,2),ant6a(6,6),ant2(2,1)
      integer nant,nantl,j3
      parameter (ncls=2)        !different class of processes
      integer nx,proc(nmax),antns(6),ndual,j1,j2
      parameter (ndual=41000)
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
c     
c     data inpint/ 6,                                                ! N V-A 
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
     >     11,1,1,1,1,  11,1,2,1,2,  11,3,2,3,2,   2,1,1,1,2, ! guu, gdd, gbb, w+ud
     >     3,1,4,1,3,  11,2,1,2,1, 11,3,1,3,1, ! w-en, gcc
     >     0,                   ! N of yukawa
     >     2,                   ! N self-gauge
     >     11,11,11,  12,11,11, ! ggg Auxgg
     >     956*-100/
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
c     rst=result
C     
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



