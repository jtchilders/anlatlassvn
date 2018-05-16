c-------------------------------------------------------------------
      subroutine alsprc
c     assigns the hard process code
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      ihrd=11
      end

c-------------------------------------------------------------------
      subroutine alhset
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c     ngrid is the total number of grids allowd for P.S. variables.
c     jgrid(jproc) labels the grid associated to a given jproc.
c       jgrid= 1   -->   q qbar and qbar q (no quarks in final state)
c       jgrid= 2   -->   g g
c       jgrid= 3   -->   g q    and g qbar
c       jgrid= 4   -->   q g    and qbar g
c       jgrid= 5   -->   qq->qq
      data jgrid/2*1,2*3,2*4,2*2,2*1,2*5,2*1,5,1,5,3*1,2*5,2*1,
     >           3,4,3,4,3,4,3,4,3,4,3,4,64*0/
      integer i
c     parameters for the gauge invariance prescription:
      winsize  = 2.d0/pi
      resonance= 'n'
      wmode    = 'yy'
c
c     process input parameters
c
      npart =njets+nph+2
      nprtns=njets
      if (nph.lt.1) then
         print*,' at least 1 photon is required!'
         stop
      endif
      if(njets.eq.0) then
         jprocmax=  2
         ngrid   =  1
         navg    =  1
      elseif(njets.eq.1) then
         jprocmax=  6
         ngrid   =  4
         navg    =  1 
      elseif(njets.eq.2) then
         jprocmax= 24
         ngrid   =  5
         navg    =  1
      elseif(njets.eq.3.or.njets.le.6) then
         jprocmax= 36
         ngrid   =  5
         navg    =  1
      else 
         print*,'njets=',njets,' not yet available'
         stop
      endif 
c     masses
      do i=1,npart
         p(5,i)=0
      enddo
c
      if(njets.eq.0) then
        write(niosta,*) nph,' photons'
      else
        write(niosta,*) nph,' photons +',njets,' jets'
      endif
      write(niosta,*) '======================================='
      write(niosta,*)
     +     'Generation cuts for the partonic event sample:' 
      if(njets.gt.0) then
        write(niosta,*) '     Light jets:'
        write(niosta,*) 'ptmin=',ptjmin,' |etamax|=',etajmax
     +       ,' dR(j-j)>',drjmin 
      endif
      write(niosta,*) '     Photons:'
      write(niosta,*) 'ptphmin=',ptphmin,' ptphmax=',ptphmax,
     +     ' |etaphax|=',etaphmax
     +     ,' dR(j-photon)>',drphjmin,' dR(ph-ph)>',drphmin  
c      write(niosta,*) 'pthrmin=',pthrmin,' pthrmax=',pthrmax
      end

c-------------------------------------------------------------------
      subroutine selflav(jproc,xlum,afl)
c     evaluates parton luminosities, and assigns PDG-code flavours
c     gluon=21 cbar=-4 sbar=-3 ubar=-2 dbar=-1 d=1 u=2 s=3 c=4
c     we keep separated processes in which Z is coupled to qu or qd
c     
c     Z stands for nph photons
c
c     jproc
c 0 light jets
c     1   qu qubar -> Z   and    qubar qu -> Z
c     2   qd qdbar -> Z   and    qdbar qd -> Z
c 1 light jets
c     3   g qu -> qu Z   and    g qubar -> qubar Z
c     4   g qd -> qd Z   and    g qdbar -> qdbar Z
c     5   qu g -> qu Z   and    qubar g -> qubar Z
c     6   qd g -> qd Z   and    qdbar g -> qdbar Z
c 2 light jets
c     7   g g -> qu qubar Z
c     8   g g -> qd qdbar Z  
c     9   qu qubar -> qu qubar Z   and   qubar qu -> qu qubar Z
c     10  qd qdbar -> qd qdbar Z   and   qdbar qd -> qd qdbar Z 
c     11  qu qu -> qu qu Z and qubar qubar -> qubar qubar Z
c     12  qd qd -> qd qd Z and qdbar qdbar -> qdbar qdbar Z   
c     13  qu qubar -> qu' qubar' Z   and   qubar qu -> qu' qubar' Z
c     14  qd qdbar -> qd' qdbar' Z   and   qdbar qd -> qd' qdbar' Z
c     15  qu qu' -> qu qu' Z and qubar qubar' -> qubar qubar' Z
c     16  qu qubar' -> qu qubar' Z and qubar qu' -> qubar qu'
c     17  qd qd' -> qd qd' Z and qdbar qdbar' -> qdbar qdbar' Z
c     18  qd qdbar' -> qd qdbar' Z and qdbar qd' -> qdbar qd' Z
c     19  qu qubar -> qd qdbar Z   and   qubar qu -> qd qdbar Z (generic d & u)
c     20  qd qdbar -> qu qubar Z   and   qdbar qd -> qu qubar Z (generic d & u)
c     21  qu qd -> qu qd Z and qubar qdbar -> qubar qdbar Z (generic d & u)
c     22  qd qu -> qu qd Z and qdbar qubar -> qubar qdbar Z (generic d & u)
c     23  qu qdbar -> qu qdbar Z and qubar qd -> qubar qd Z (generic d & u)
c     24  qd qubar -> qd qubar Z and qdbar qu -> qdbar qu Z (generic d & u)
c     25  g qu     -> qu qu qubar Z   and   g qubar -> qubar qu qubar Z
c     26  qu g     -> qu qu qubar Z   and   g qubar -> qubar qu qubar Z
c     27  g qu     -> qu qu' qubar' Z   and   g qubar -> qubar qu' qubar' Z
c     28  qu g     -> qu qu' qubar' Z   and   g qubar -> qubar qu' qubar' Z
c     29  g qu     -> qu qd qdbar Z   and   g qubar -> qubar qd qdbar Z (generic d & u)
c     30  qu g     -> qu qd qdbar Z   and   g qubar -> qubar qd qdbar Z (generic d & u)
c     31  g qd     -> qd qd qdbar Z   and   g qdbar -> qdbar qd qdbar Z
c     32  qd g     -> qd qd qdbar Z   and   g qdbar -> qdbar qd qdbar Z
c     33  g qd     -> qd qd' qdbar' Z   and   g qdbar -> qdbar qd' qdbar' Z
c     34  qd g     -> qd qd' qdbar' Z   and   g qdbar -> qdbar qd' qdbar' Z
c     35  g qd     -> qd qu qubar Z   and   g qdbar -> qdbar qu qubar Z (generic d & u)
c     36  qd g     -> qd qu qubar Z   and   g qdbar -> qdbar qu qubar Z (generic d & u)

c-------------------------------------------------------------------   
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
      integer afl(maxpar),iflag
      real tmp(100),slum,cwgt,swgt,rn,tmptot
c
c     N.B. 1=d 2=u 3=s 4=c 5=b and negatives are antiparticles
c
c     definition of qu
c
      integer qu(-2:2)
      data qu/-4,-2,0,2,4/
c
c     definition of qd
c
      integer qd(-2:2)
      data qd/-3,-1,0,1,3/
c
c     color
c
      real *8 xlum,xrn                                 
      integer i,j,k,l,itmp,icount,init,jproc,ng,nlqp
      real cfac(0:5),ifact(0:8)
      data cfac/8e0,5*3e0/,init/0/
c
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
      save init,cfac,ifact
c
      if(init.eq.0) then
         njets=npart-2-nph
         ifact(0)=1e0
         do i=1,8
            ifact(i)=ifact(i-1)/real(i)
         enddo
      endif
c     
      do i=1,npart
         ifl(i)=0
      enddo
c
 1    call randa(xrn)
      rn=real(xrn)
      if(1e0-rn.lt.1e-7) goto 1
c
      icount=0
      slum=0e0
      tmptot=0e0
c
c     -----------------------------------
c       start two light-quark processes   
c     -----------------------------------
c     ------------------------------------
c     jproc=1   njets>=0                  done
c     qu qubar -> Z   and    qubar qu -> Z 
c     ------------------------------------
c
      if(jproc.eq.1) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(qu(i))*f2(-qu(i))
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=qu(i)
                  ifl(2)=-qu(i)
                  do k=3,njets+2
                     ifl(k)=0
                  enddo
c
                  itmp=qu(i)
                  goto 100
               endif 
            endif
         enddo
c 
c     ------------------------------------
c     jproc=2   njets>=0                  done
c     qd qdbar -> Z   and    qdbar qd -> Z 
c     ------------------------------------
c
      elseif(jproc.eq.2) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(qd(i))*f2(-qd(i))
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=qd(i)
                  ifl(2)=-qd(i)
                  do k=3,njets+2
                     ifl(k)=0
                  enddo
c
                  itmp=qd(i)
                  goto 100
               endif 
            endif
         enddo
c 
c     -----------------------------------
c     jproc=3   njets>=1                 done
c     g qu -> qu Z   and   g qubar -> qubar Z
c     -----------------------------------
c
      elseif(jproc.eq.3) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(0)*f2(qu(i))
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=qu(i)
                  ifl(3)=qu(i)
                  do k=4,njets+2
                     ifl(k)=0
                  enddo
c
                  itmp=qu(i)
                  goto 100   
               endif 
            endif
         enddo
c 
c     -----------------------------------
c     jproc=4   njets>=1                 done
c     g qd -> qd Z   and   g qdbar -> qdbar Z
c     -----------------------------------
c
      elseif(jproc.eq.4) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(0)*f2(qd(i))
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=qd(i)
                  ifl(3)=qd(i)
                  do k=4,njets+2
                     ifl(k)=0
                  enddo
c
                  itmp=qd(i)
                  goto 100   
               endif 
            endif
         enddo 
c
c     -----------------------------------
c     jproc=5   njets>=1                 done
c     qu g -> qu Z   and   qubar g -> qubar Z
c     -----------------------------------
c
      elseif(jproc.eq.5) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(qu(i))*f2(0)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=qu(i)
                  ifl(2)=0
                  ifl(3)=qu(i)
                  do k=4,njets+2
                     ifl(k)=0
                  enddo
c
                  itmp=qu(i)
                  goto 100   
               endif 
            endif
         enddo 
c
c     -----------------------------------
c     jproc=6   njets>=1                 done          
c     qd g -> qd Z   and   qdbar g -> qdbar Z
c     -----------------------------------
c
      elseif(jproc.eq.6) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(qd(i))*f2(0)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=qd(i)
                  ifl(2)=0
                  ifl(3)=qd(i)
                  do k=4,njets+2
                     ifl(k)=0
                  enddo
c
                  itmp=qd(i)
                  goto 100   
               endif 
            endif
         enddo
c
c     --------------------------
c     jproc=7   njets>=2       done
c     gg -> qu qubar Z
c     --------------------------
c             
      elseif(jproc.eq.7) then
         do i=1,2
            icount=icount+1
            tmp(icount)=f1(0)*f2(0)
            slum=slum+tmp(icount)
         enddo 
         icount=0
         rn=rn*slum
         do i=1,2
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
               ifl(1)=0
               ifl(2)=0
               ifl(3)=qu(i)
               ifl(4)=-qu(i)
               do k=5,njets+2
                  ifl(k)=0
               enddo 
c
               itmp=1
               goto 100
            endif
         enddo
c 
c     --------------------------
c     jproc=8   njets>=2        done       
c     gg -> qd qdbar Z
c     --------------------------
c             
      elseif(jproc.eq.8) then
         do i=1,2
            icount=icount+1
            tmp(icount)=f1(0)*f2(0)
            slum=slum+tmp(icount)
         enddo 
         icount=0
         rn=rn*slum
         do i=1,2
            icount=icount+1
            tmptot=tmptot+tmp(icount)
            if(tmptot.ge.rn) then 
               ifl(1)=0
               ifl(2)=0
               ifl(3)=qd(i)
               ifl(4)=-qd(i)
               do k=5,njets+2
                  ifl(k)=0
               enddo 
c
               itmp=1
              goto 100
            endif
         enddo 
c
c     ------------------------------
c     start 4 light quark processes 
c     ------------------------------
c     ------------------------------
c     jproc=9  njets>=2              done
c     qu qubar -> qu qubar Z   and   qubar qu -> qu qubar Z
c     ------------------------------
c
      elseif(jproc.eq.9) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(qu(i))*f2(-qu(i))
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=qu(i)
                  ifl(2)=-qu(i)
                  ifl(3)=qu(i)
                  ifl(4)=-qu(i)
                  do l=5,njets+2
                     ifl(l)=0
                  enddo 
c
                  itmp=qu(i)
                  goto 100
               endif
            endif
         enddo
c
c     ------------------------------
c     jproc=10 njets>=2             done
c     qd qdbar -> qd qdbar Z   and   qdbar qd -> qd qdbar Z
c     ------------------------------
c
      elseif(jproc.eq.10) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(qd(i))*f2(-qd(i))
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=qd(i)
                  ifl(2)=-qd(i)
                  ifl(3)=qd(i)
                  ifl(4)=-qd(i)
                  do l=5,njets+2
                     ifl(l)=0
                  enddo 
c
                  itmp=qd(i)
                  goto 100
               endif
            endif
         enddo
c
c     ------------------------------
c     jproc=11 njets>=2             done
c     qu qu -> qu qu Z   and   qubar qubar -> qubar qubar Z
c     ------------------------------
c
      elseif(jproc.eq.11) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(qu(i))*f2(qu(i))
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=qu(i)
                  ifl(2)=qu(i)
                  ifl(3)=qu(i)
                  ifl(4)=qu(i)
                  do l=5,njets+2
                     ifl(l)=0
                  enddo 
c
                  itmp=qu(i)
c     statistical factor to account for the presence of indist. particles
                  slum= slum/2.d0
                  goto 100
               endif
            endif
         enddo
c
c     ------------------------------
c     jproc=12 njets>=2             done
c     qd qd -> qd qd Z   and   qdbar qdbar -> qdbar qdbar Z
c     ------------------------------
c
      elseif(jproc.eq.12) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(qd(i))*f2(qd(i))
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=qd(i)
                  ifl(2)=qd(i)
                  ifl(3)=qd(i)
                  ifl(4)=qd(i)
                  do l=5,njets+2
                     ifl(l)=0
                  enddo 
c
                  itmp=qd(i)
c     statistical factor to account for the presence of indist. particles
                  slum= slum/2.d0
                  goto 100
               endif
            endif
         enddo
c
c     -------------------------------------------------
c     jproc=13  njets>=2                               done
c     qu qubar -> qu' qubar' Z   and   qubar qu -> qu' qubar' Z
c     -------------------------------------------------
c
      elseif(jproc.eq.13) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then 
                        icount=icount+1
                        tmp(icount)=f1(qu(i))*f2(-qu(i))
                        slum=slum+tmp(icount)
                     endif
                  endif 
                enddo
            endif
         enddo 
         icount=0
c
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                    if(abs(i-j).eq.1) then
                     icount=icount + 1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=qu(i)
                        ifl(2)=-qu(i)
                        ifl(3)=qu(j)
                        ifl(4)=-qu(j)
                        do l=5,njets+2
                           ifl(l)=0
                        enddo 
c
                        itmp=qu(i)
                        goto 100
                      endif 
                   endif
                endif 
             enddo
           endif
         enddo

c
c     -------------------------------------------------
c     jproc=14  njets>=2                            done
c     qd qdbar -> qd' qdbar' Z   and   qdbar qd -> qd' qdbar' Z
c     -------------------------------------------------
c
      elseif(jproc.eq.14) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then 
                        icount=icount+1
                        tmp(icount)=f1(qd(i))*f2(-qd(i))
                        slum=slum+tmp(icount)
                     endif
                  endif 
                enddo
            endif
         enddo 
         icount=0
c
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then
                        icount=icount + 1
                        tmptot=tmptot+tmp(icount)
                        if(tmptot.ge.rn) then 
                           ifl(1)=qd(i)
                           ifl(2)=-qd(i)
                           ifl(3)=qd(j)
                           ifl(4)=-qd(j)
                           do l=5,njets+2
                              ifl(l)=0
                           enddo 
c
                           itmp=qd(i)
                           goto 100
                        endif 
                     endif
                  endif 
               enddo
            endif
         enddo
c
c     ------------------------------
c     jproc=15 njets>=2              done
c     qu qu' -> qu qu' Z   and qubar qubar' -> qubar qubar' Z
c     ------------------------------
c
      elseif(jproc.eq.15) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(qu(i)).ne.abs(qu(j)).and.
     #                  i*j.gt.0) then 
                        icount=icount+1
                        tmp(icount)=f1(qu(i))*f2(qu(j))
                        slum=slum+tmp(icount)
                     endif
                  endif 
                enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then      
                     if(abs(qu(i)).ne.abs(qu(j)).and.
     #                  i*j.gt.0) then 
                        icount=icount+1
                        tmptot=tmptot+tmp(icount)
                        if(tmptot.ge.rn) then 
                           ifl(1)=qu(i)
                           ifl(2)=qu(j)
                           ifl(3)=qu(i)
                           ifl(4)=qu(j)
                           do l=5,njets+2
                              ifl(l)=0
                           enddo 
c
                           itmp=qu(i)
                           goto 100
                        endif
                     endif
                  endif
               enddo  
            endif
         enddo
c
c     ------------------------------
c     jproc=16 njets>=2               done
c     qu qubar' -> qu qubar' Z and qubar qu' -> qubar qu'
c     ------------------------------
c
      elseif(jproc.eq.16) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(qu(i)).ne.abs(qu(j)).and.
     #                  i*j.lt.0) then 
                        icount=icount+1
                        tmp(icount)=f1(qu(i))*f2(qu(j))
                        slum=slum+tmp(icount)
                     endif
                  endif 
                enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then      
                     if(abs(qu(i)).ne.abs(qu(j)).and.
     #                  i*j.lt.0) then 
                        icount=icount+1
                        tmptot=tmptot+tmp(icount)
                        if(tmptot.ge.rn) then 
                           ifl(1)=qu(i)
                           ifl(2)=qu(j)
                           ifl(3)=qu(i)
                           ifl(4)=qu(j)
                           do l=5,njets+2
                              ifl(l)=0
                           enddo 
c
                           itmp=qu(i)
                           goto 100
                        endif
                     endif
                  endif
               enddo  
            endif
         enddo
c
c     ------------------------------
c     jproc=17 njets>=2             done
c     qd qd' -> qd qd' Z   and qdbar qdbar' -> qdbar qdbar' Z
c     ------------------------------
c
      elseif(jproc.eq.17) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(qd(i)).ne.abs(qd(j)).and.
     #                  i*j.gt.0) then 
                        icount=icount+1
                        tmp(icount)=f1(qd(i))*f2(qd(j))
                        slum=slum+tmp(icount)
                     endif
                  endif 
                enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then      
                     if(abs(qd(i)).ne.abs(qd(j)).and.
     #                  i*j.gt.0) then 
                        icount=icount+1
                        tmptot=tmptot+tmp(icount)
                        if(tmptot.ge.rn) then 
                           ifl(1)=qd(i)
                           ifl(2)=qd(j)
                           ifl(3)=qd(i)
                           ifl(4)=qd(j)
                           do l=5,njets+2
                              ifl(l)=0
                           enddo 
c
                           itmp=qd(i)
                           goto 100
                        endif
                     endif
                  endif
               enddo  
            endif
         enddo
c
c     ------------------------------
c     jproc=18 njets>=2             done
c     qd qdbar' -> qd qdbar' Z and qdbar qd' -> qdbar qd' Z
c     ------------------------------
c
      elseif(jproc.eq.18) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(qd(i)).ne.abs(qd(j)).and.
     #                  i*j.lt.0) then 
                        icount=icount+1
                        tmp(icount)=f1(qd(i))*f2(qd(j))
                        slum=slum+tmp(icount)
                     endif
                  endif 
                enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then      
                     if(abs(qd(i)).ne.abs(qd(j)).and.
     #                  i*j.lt.0) then 
                        icount=icount+1
                        tmptot=tmptot+tmp(icount)
                        if(tmptot.ge.rn) then 
                           ifl(1)=qd(i)
                           ifl(2)=qd(j)
                           ifl(3)=qd(i)
                           ifl(4)=qd(j)
                           do l=5,njets+2
                              ifl(l)=0
                           enddo 
c
                           itmp=qd(i)
                           goto 100
                        endif
                     endif
                  endif
               enddo  
            endif
         enddo
c
c     ------------------------------
c     jproc=19 njets>=2             done
c     qu qubar -> qd qdbar Z   and   qubar qu -> qd qdbar Z (generic qu and qd)
c     ------------------------------
c
      elseif(jproc.eq.19) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(qu(i))*f2(-qu(i))
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=qu(i)
                        ifl(2)=-qu(i)
                        ifl(3)=qd(j)             
                        ifl(4)=-qd(j)            
                        do l=5,njets+2
                           ifl(l)=0
                        enddo 
c
                        itmp=qu(i)
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo
c
c     ------------------------------
c     jproc=20 njets>=2            done
c     qd qdbar -> qu qubar Z   and   qdbar qd -> qu qubar Z (generic qu and qd)
c     ------------------------------
c
      elseif(jproc.eq.20) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(qd(i))*f2(-qd(i))
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then
                        ifl(1)=qd(i)
                        ifl(2)=-qd(i)
                        ifl(3)=qu(j)      
                        ifl(4)=-qu(j)     
                        do l=5,njets+2
                           ifl(l)=0
                        enddo 
C
                        itmp=qd(i)
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo
c
c     ------------------------------
c     jproc=21 njets>=2             done
c     qu qd -> qu qd Z and qubar qdbar -> qubar qdbar Z (generic qu and qd)
c     ------------------------------
c
      elseif(jproc.eq.21) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(qu(i))*f2(qd(j))
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=qu(i)
                        ifl(2)=qd(j)
                        ifl(3)=qu(i)
                        ifl(4)=qd(j)
                        do l=5,njets+2
                           ifl(l)=0
                        enddo 
c
                        itmp=qu(i)
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo
c
c     ------------------------------
c     jproc=22 njets>=2              done
c     qd qu -> qd qu Z and the qdbar qubar -> qdbar qubar Z (genric qu and qd)
c     ------------------------------
c
      elseif(jproc.eq.22) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(qd(i))*f2(qu(j))
                     slum=slum+tmp(icount)
                  endif
               enddo 
            endif
         enddo
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=qd(i)
                        ifl(2)=qu(j)
                        ifl(3)=qd(i)
                        ifl(4)=qu(j)
                        do l=5,njets+2
                           ifl(l)=0
                        enddo 
c
                        itmp=qd(i)
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo
c
c     ------------------------------
c     jproc=23 njets>=2             done
c     qu qdbar -> qu qdbar Z and qubar qd -> qubar qd Z (generic qu and qd)
c     ------------------------------
c
      elseif(jproc.eq.23) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.lt.0) then
                     icount=icount+1
                     tmp(icount)=f1(qu(i))*f2(qd(j))
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.lt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=qu(i)
                        ifl(2)=qd(j)
                        ifl(3)=qu(i)
                        ifl(4)=qd(j)
                        do l=5,njets+2
                           ifl(l)=0
                        enddo 
c
                        itmp=qu(i)
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo
c
c     ------------------------------
c     jproc=24 njets>=2             done
c     qd qubar -> qd qubar Z and the qdbar qu -> qdbar qu Z (generic qu and qd)
c     ------------------------------
c
      elseif(jproc.eq.24) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.lt.0) then
                     icount=icount+1
                     tmp(icount)=f1(qd(i))*f2(qu(j))
                     slum=slum+tmp(icount)
                  endif
               enddo 
            endif
         enddo
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0.and.i*j.lt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=qd(i)
                        ifl(2)=qu(j)
                        ifl(3)=qd(i)
                        ifl(4)=qu(j)
                        do l=5,njets+2
                           ifl(l)=0
                        enddo 
c
                        itmp=qd(i)
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo
c 
c     -----------------------------------
c     jproc=25   njets>=3
c     g qu -> qu qu qubar Z   and   g qubar -> qubar qu qubar Z
c     -----------------------------------
c
      elseif(jproc.eq.25) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(0)*f2(qu(i))
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=qu(i)
                  ifl(3)=qu(i)
                  ifl(4)=qu(i)
                  ifl(5)=-qu(i)
                  do k=6,njets+2
                     ifl(k)=0
                  enddo
c     to account for 2 id. particles
                  slum= slum/2.d0
c
                  itmp=qu(i)
                  goto 100   
               endif 
            endif
         enddo
c 
c     -----------------------------------
c     jproc=26   njets>=3
c     qu g -> qu qu qubar Z   and   qubar g -> qubar qu qubar Z
c     -----------------------------------
c
      elseif(jproc.eq.26) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(qu(i))*f2(0)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=qu(i)
                  ifl(2)=0
                  ifl(3)=qu(i)
                  ifl(4)=qu(i)
                  ifl(5)=-qu(i)
                  do k=6,njets+2
                     ifl(k)=0
                  enddo
c     to account for 2 id. particles
                  slum= slum/2.d0
c
                  itmp=qu(i)
                  goto 100   
               endif 
            endif
         enddo
c 
c     -----------------------------------
c     jproc=27   njets>=3
c     g qu -> qu qu' qubar' Z   and   g qubar -> qubar qu' qubar' Z
c     -----------------------------------
c
      elseif(jproc.eq.27) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then 
                        icount=icount+1
                        tmp(icount)=f1(0)*f2(qu(i))
                        slum=slum+tmp(icount)
                     endif
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then 
                        icount=icount+1
                        tmptot=tmptot+tmp(icount)
                        if(tmptot.ge.rn) then 
                           ifl(1)=0
                           ifl(2)=qu(i)
                           ifl(3)=qu(i)
                           ifl(4)=qu(j)
                           ifl(5)=-qu(j)
                           do k=6,njets+2
                              ifl(k)=0
                           enddo
c     to account for 2 id. particles
                           slum= slum/2.d0
c
                           itmp=qu(i)
                           goto 100   
                        endif 
                     endif
                  endif
               enddo
            endif
         enddo
c 
c     -----------------------------------
c     jproc=28   njets>=3
c     qu g -> qu qu' qubar' Z   and   qubar g -> qubar qu' qubar' Z
c     -----------------------------------
c
      elseif(jproc.eq.28) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then 
                        icount=icount+1
                        tmp(icount)=f1(qu(i))*f2(0)
                        slum=slum+tmp(icount)
                     endif
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then 
                        icount=icount+1
                        tmptot=tmptot+tmp(icount)
                        if(tmptot.ge.rn) then 
                           ifl(1)=qu(i)
                           ifl(2)=0
                           ifl(3)=qu(i)
                           ifl(4)=qu(j)
                           ifl(5)=-qu(j)
                           do k=6,njets+2
                              ifl(k)=0
                           enddo
c     to account for 2 id. particles
                           slum= slum/2.d0
c
                           itmp=qu(i)
                           goto 100   
                        endif 
                     endif
                  endif
               enddo
            endif
         enddo
c 
c     -----------------------------------
c     jproc=29   njets>=3
c     g qu -> qu qd qdbar Z   and   g qubar -> qubar qd qdbar Z (generic u & d)
c     -----------------------------------
c
      elseif(jproc.eq.29) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     icount=icount+1
                     tmp(icount)=f1(0)*f2(qu(i))
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=0
                        ifl(2)=qu(i)
                        ifl(3)=qu(i)
                        ifl(4)=qd(j)
                        ifl(5)=-qd(j)
                        do k=6,njets+2
                           ifl(k)=0
                        enddo
c     to account for 2 id. particles
                        slum= slum/2.d0
c
                        itmp=qu(i)
                        goto 100   
                     endif
                  endif
               enddo
            endif
         enddo
c 
c     -----------------------------------
c     jproc=30   njets>=3
c     qu g -> qu qd qdbar Z   and   qubar g -> qubar qd qdbar Z (generic u & d)
c     -----------------------------------
c
      elseif(jproc.eq.30) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     icount=icount+1
                     tmp(icount)=f1(qu(i))*f2(0)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=qu(i)
                        ifl(2)=0
                        ifl(3)=qu(i)
                        ifl(4)=qd(j)
                        ifl(5)=-qd(j)
                        do k=6,njets+2
                           ifl(k)=0
                        enddo
c     to account for 2 id. particles
                        slum= slum/2.d0
c
                        itmp=qu(i)
                        goto 100   
                     endif
                  endif
               enddo
            endif
         enddo
c 
c     -----------------------------------
c     jproc=31   njets>=3
c     g qd -> qd qd qdbar Z   and   g qdbar -> qdbar qd qdbar Z
c     -----------------------------------
c
      elseif(jproc.eq.31) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(0)*f2(qd(i))
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=qd(i)
                  ifl(3)=qd(i)
                  ifl(4)=qd(i)
                  ifl(5)=-qd(i)
                  do k=6,njets+2
                     ifl(k)=0
                  enddo
c     to account for 2 id. particles
                  slum= slum/2.d0
c
                  itmp=qd(i)
                  goto 100   
               endif 
            endif
         enddo
c 
c     -----------------------------------
c     jproc=32   njets>=3
c     qd g -> qd qd qdbar Z   and   qdbar g -> qdbar qd qdbar Z
c     -----------------------------------
c
      elseif(jproc.eq.32) then
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(qd(i))*f2(0)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=qd(i)
                  ifl(2)=0
                  ifl(3)=qd(i)
                  ifl(4)=qd(i)
                  ifl(5)=-qd(i)
                  do k=6,njets+2
                     ifl(k)=0
                  enddo
c     to account for 2 id. particles
                  slum= slum/2.d0
c
                  itmp=qd(i)
                  goto 100   
               endif 
            endif
         enddo
c 
c     -----------------------------------
c     jproc=33   njets>=3
c     g qd -> qd qd' qdbar' Z   and   g qdbar -> qdbar qd' qdbar' Z
c     -----------------------------------
c
      elseif(jproc.eq.33) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then 
                        icount=icount+1
                        tmp(icount)=f1(0)*f2(qd(i))
                        slum=slum+tmp(icount)
                     endif
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then 
                        icount=icount+1
                        tmptot=tmptot+tmp(icount)
                        if(tmptot.ge.rn) then 
                           ifl(1)=0
                           ifl(2)=qd(i)
                           ifl(3)=qd(i)
                           ifl(4)=qd(j)
                           ifl(5)=-qd(j)
                           do k=6,njets+2
                              ifl(k)=0
                           enddo
c     to account for 2 id. particles
                           slum= slum/2.d0
c
                           itmp=qd(i)
                           goto 100   
                        endif 
                     endif
                  endif
               enddo
            endif
         enddo
c 
c     -----------------------------------
c     jproc=34   njets>=3
c     qd g -> qd qd' qdbar' Z   and   qdbar g -> qdbar qd' qdbar' Z
c     -----------------------------------
c
      elseif(jproc.eq.34) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then 
                        icount=icount+1
                        tmp(icount)=f1(qd(i))*f2(0)
                        slum=slum+tmp(icount)
                     endif
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     if(abs(i-j).eq.1) then 
                        icount=icount+1
                        tmptot=tmptot+tmp(icount)
                        if(tmptot.ge.rn) then 
                           ifl(1)=qd(i)
                           ifl(2)=0
                           ifl(3)=qd(i)
                           ifl(4)=qd(j)
                           ifl(5)=-qd(j)
                           do k=6,njets+2
                              ifl(k)=0
                           enddo
c     to account for 2 id. particles
                           slum= slum/2.d0
c
                           itmp=qd(i)
                           goto 100   
                        endif 
                     endif
                  endif
               enddo
            endif
         enddo
c 
c     -----------------------------------
c     jproc=35   njets>=3
c     g qd -> qd qu qubar Z   and   g qdbar -> qdbar qu qubar Z (generic u & d)
c     -----------------------------------
c
      elseif(jproc.eq.35) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     icount=icount+1
                     tmp(icount)=f1(0)*f2(qd(i))
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=0
                        ifl(2)=qd(i)
                        ifl(3)=qd(i)
                        ifl(4)=qu(j)
                        ifl(5)=-qu(j)
                        do k=6,njets+2
                           ifl(k)=0
                        enddo
c     to account for 2 id. particles
                        slum= slum/2.d0
c
                        itmp=qd(i)
                        goto 100   
                     endif
                  endif
               enddo
            endif
         enddo
c 
c     -----------------------------------
c     jproc=36   njets>=3
c     qd g -> qd qu qubar Z   and   qdbar g -> qdbar qu qubar Z (generic u & d)
c     -----------------------------------
c
      elseif(jproc.eq.36) then
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     icount=icount+1
                     tmp(icount)=f1(qd(i))*f2(0)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-2,2
            if(i.ne.0) then
               do j=-2,2
                  if(j.ne.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=qd(i)
                        ifl(2)=0
                        ifl(3)=qd(i)
                        ifl(4)=qu(j)
                        ifl(5)=-qu(j)
                        do k=6,njets+2
                           ifl(k)=0
                        enddo
c     to account for 2 id. particles
                        slum= slum/2.d0
c
                        itmp=qd(i)
                        goto 100   
                     endif
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
c
 100  continue 
c
      do i= njets+3,npart
         ifl(i)= 22
      enddo

      do i=1,10
         afl(i)=ifl(i)
      enddo

c     complete flavour-dependent momenta assignements in the user common block:
      call usrfll(iflag)
      if(iflag.eq.1) then
         xlum=-1d0
         return
      endif
c
c     evaluate colour weight factors; the gluons are renamed 21
c
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
c
c     evaluate spin weight factors
c
      swgt=2e0
      swgt=swgt**(njets+nph)

      xlum=dble(slum*cwgt*swgt*ifact(ng)*ifact(nph))/resc**njets
c     
c         effco=ccoef(#quark pairs,#gluons)
c         nlqp=number of light quark pairs in the process.
c         2+njets-2*nlqp
c         is the number of external gluons
c
      if(jproc.le.8) then
         nlqp= 1
      elseif(jproc.gt.8.and.jproc.le.36) then
         nlqp= 2
      else
         write(*,*) 'jproc not valid'
         stop
      endif 
      xlum=xlum*ccoef(nlqp,2+njets-2*nlqp)
c     
      end         

      subroutine usrfll(iflag)
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
      integer iflag,i,j
c     
      iflag=0
c assign momenta to the usr common block: CHECK USE OF EQUIVALENCE STATEMENT
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
         do j=1,nph
            pphot(i,j)=p(i,j+2+njets)
         enddo
      enddo
*
      do j=1,njets
         ptj(j) = pt(j+2)
         etaj(j)= eta(j+2)
         do i= j+1,njets
            drjj(i,j)= dr(i+2,j+2)
            drjj(j,i)= drjj(i,j)
         enddo
      enddo
*
      do j=1,nph
         ptp(j) = pt(j+2+njets)
         etap(j)= eta(j+2+njets)
         do i= j+1,nph
            drpp(i,j)= dr(i+2+njets,j+2+njets)
            drpp(j,i)= drpp(i,j)
         enddo
      enddo
*     
      return
 10   iflag=1
      end
c
      subroutine setdec(nwrt,iflwrt,icuwrt,pwrt,decbr)
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
c debug
      double precision rate(16)
      common/dbg/rate
      data rate/16*0d0/
c locals
      integer maxdec
      parameter (maxdec=40)
      integer ip, ic, il
c arguments
      integer nwrt,iflwrt(maxdec),icuwrt(2,maxdec)
      double precision pwrt(5,maxdec),decbr
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
      end
c
      subroutine phspace(lnot,pswgt,djpd,djg)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c  (njets+nph)-body phase space for the process:                 c
c                                                                c
c  h(1) h(2) -> h(3)..h(njets+2) ph(njets+3)... ph(njets+nph+2)  c
c                                                                c
c  where h is any light-flavoured quark or gluon                 c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
      real *8 djg,dummy
      real *8 pswgt
      real *8 djpd,factor
      real *8 cutkin(12)
      real *8 wgt
      common/loccut/cutkin
      real *8 pl(maxpar),y(maxpar)
      real *8 pcm(0:3,maxpar)
      real *8 totpt,q0sq
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
         cutkin(7)=ptphmin
         cutkin(8)=etaphmax
         cutkin(9)=drphjmin
         cutkin(10)=drphmin
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
      call momgen(njets,nph,roots,x1,x2,pcm,wgt,lnot)
      djg= 1.d0 ! dummy variable
      if (lnot.eq.1) then
         pswgt= 0.d0
         goto 100
      endif
*
c-    will write factor=factor0/(x1*x2), with factor0 function of njets etc.
*
      factor= 1d0/(2.d0*pi)**(3*(njets+nph)-4)/2.d0/s/x1/x2
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

      do l= 3,njets+nph+2
          pt(l) = sqrt(pcm(1,l)**2+pcm(2,l)**2)
          pl(l) = pcm(3,l)
          eta(l)= -log(tan(0.5d0*atan2(pt(l),pl(l))))
          y(l)  = 0.5d0*log((pcm(0,l)+pcm(3,l))/(pcm(0,l)-pcm(3,l)))
      enddo
*
c-    Calculates jet-jet distances:
*
      do l= 3,njets+nph+1
         do m= l+1,njets+nph+2
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
      do l= 3,njets+nph+1
         do m= l+1,njets+nph+2
            dr(l,m)= sqrt(dphi(l,m)**2+(eta(l)-eta(m))**2)
            dr(m,l)=dr(l,m)
         enddo
      enddo
*
c-    Redefine the momenta:           
*                                     
      do l= 1,njets+nph+2                
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
      call chkcut(lnot,pt,p,eta,dr,njets,nph)
      if (lnot.eq.1) then
         pswgt= 0.d0
         goto 100
      endif
*
      pswgt = factor*wgt
*
c-    Evaluate q2: will include several possible options:
*
c-    Total et^2 of jets
*
      totpt= 0.d0
      do i=3,njets+2
         totpt=totpt+pt(i)**2
      enddo 
*
      if(iqopt.eq.0) then 
         qsq=1d0
      elseif(iqopt.eq.1) then
*
c-       Total et^2 of photons
*
         q0sq= 0.d0
         do i=njets+3,njets+2+nph
           q0sq=q0sq+pt(i)**2
         enddo 
         qsq=q0sq+totpt
      elseif(iqopt.eq.2) then
         qsq=totpt
      endif
      qsq=qfac**2*qsq
 100  continue
      end
*     
      subroutine chkcut(lnot,pt,p,eta,dr,njets,nph)
c=======================================================
c=======================================================
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     Applies kinematical cuts to the final state during the phase
c     -space generation
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 ptjmin,drjmin,etajmax
      real*8 ptphmin,drphjmin,drphmin,etaphmax
      integer maxpar,ninit,njets,nph,j,lnot,i
      real*8 cutkin(12)
      common/loccut/cutkin
      data ninit/0/
      parameter (maxpar=10)
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),p(5,maxpar)
      save ninit,ptjmin,etajmax,drjmin,ptphmin,etaphmax,drphjmin,drphmin
      if(ninit.eq.0) then
         ninit=1
         ptjmin  = cutkin(1)
         etajmax = cutkin(3)
         drjmin  = cutkin(5)
         ptphmin = cutkin(7)
         etaphmax= cutkin(8)
         drphjmin = cutkin(9)
         drphmin = cutkin(10)
      endif
      lnot= 0
*
c     Cuts involving jets, if any...
*
      if(njets.gt.0) then
         do i= 3,njets+2
            if (pt(i).lt.ptjmin)           goto 10
            if (abs(eta(i)).gt.etajmax)    goto 10
         enddo 
*
c        dR(jet-jet)<drjmin
*
         do i= 3,njets+1
           do j= i+1,njets+2
             if(dr(i,j).lt.drjmin)         goto 10
           enddo
         enddo
      endif 
*
c-    cuts involving photons 
*
      do i= njets+3,njets+nph+2
        if (pt(i).lt.ptphmin)           goto 10
        if (abs(eta(i)).gt.etaphmax)    goto 10
      enddo
*
c     dR(jet-photon)<drphjmin
*
      do i= njets+3,njets+nph+2
        do j= 3,njets+2
          if(dr(i,j).lt.drphjmin)         goto 10
        enddo
      enddo
*
c     dR(phot-phot)<drphmin
*
      do i= njets+3,njets+nph+1
        do j= i+1,njets+nph+2
          if(dr(i,j).lt.drphmin)         goto 10
        enddo
      enddo
*
c     stop here the cuts, for the generic user:
*
      return
 10   lnot= 1
      return
      end
*
      subroutine momgen(njets,nph,roots,x1,x2,p,wgt,lw)
*
c-    Generator of np particles in the LAB frame.
*
      implicit none
      real*8 wgt1,wgt2
      real*8 roots,x1,x2,wgt,zero,etacut
      real*8 cutkin,s
      real*8 ag,tau0,tau,xmsum,xmsum0,ptlim
      real*8 rootsh,sq,y0,yr,wtau,wjr
      real*8 pt0lmax,pt0lsum,en,pz,eta0min,ran0
      integer lw,lw1,npar,mpar,np,njets,nph,j,k
      real *8 ptjmin,ptjmax,etajmin,etajmax,drjmin,
     +        ptphmin,etaphmax,drphjmin
      parameter (npar= 20)
      real*8 pm(0:4,npar),pt0(npar),pt1(npar),xm(npar),xmr(npar),
     .       eta0(npar),pt0_l(npar),pt1_l(npar),p(0:3,npar)
      real*8 ranram(111)
      common/loccut/cutkin(12)
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
        np     = njets+nph
        zero   = 0.d0
        etacut = 40.d0
        s      = roots*roots
*
        ptjmin = cutkin(1)
        etajmax= cutkin(3)
        drjmin = cutkin(5)
*
        ptphmin = cutkin(7) 
        etaphmax= cutkin(8)
        drphjmin = cutkin(9)
*
c-      Parameters for the light jets:
*
        do j= 1,njets
          pt0_l(j) = ptjmin
          pt1_l(j) = roots/2.d0
          eta0(j)  = etajmax
          xm(j)    = zero        
        enddo
*
c-      Parameters for the photons:
*
        do j= njets+1,np
          pt0_l(j) = ptphmin
          pt1_l(j) = roots/2.d0
          eta0(j)  = etaphmax
          xm(j)    = zero        
        enddo
*
c-      Parameters for the x1,x2 integration:
*
        ag  = 0.98d0
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
      xmsum = xmsum0
*
c-    tau0 is the lower cut on tau= x1*x2
*
      tau0   = 1.d0/s*(xmsum)**2
      tau0   = max(tau0,4.d0*pt0lmax**2/s)  
      tau0   = max(tau0,pt0lsum**2/s)
      tau0   = max(tau0,1.d0/s*(float(njets)*ptjmin
     +                         +float(nph)  *ptphmin)**2)

      if(tau0.gt.1.d0) goto 100
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
c-------------------------------------------------------------------
      subroutine alsgrd
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'phjet.inc'
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
      nv  = 2*(njets+nph)-1
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
      include 'phjet.inc'
      end

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
         integer ncls,ant2(2,1),ant4(4,2),ant6a(6,6)
         integer nant,nantl,j3
         parameter (ncls=2)       !different class of processes
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
         data inpint/12,                                                ! N V-A 
c     >        1,1,1,1,1,   1,1,2,1,2,   1,2,1,2,1,   1,2,2,2,2, ! zuu, zdd, zcc, zss 
     >       10,1,1,1,1,  10,1,2,1,2,  10,2,1,2,1,  10,2,2,2,2, ! Auu, Add, Acc, Ass
     >       11,1,1,1,1,  11,1,2,1,2,  11,2,1,2,1,  11,2,2,2,2, ! guu, gdd, gcc, gss
     >       11,3,1,3,1,  11,3,2,3,2,  10,3,1,3,1,  10,3,2,3,2, ! gtt, gbb, Att, Abb
c     >        1,3,1,3,1,   1,3,2,3,2,   1,1,4,1,4,  10,1,4,1,4, ! Zee, Aee,Ztt, Zbb,
c     >        1,1,3,1,3,                                         ! Znene
     >                0,                                          ! N of yukawa
     >                2,                                          ! N self-gauge
     >                11,11,11,  12,11,11,                        !  ggg, Xgg
     >                0,                                          ! N H-GAUGE
     >                0,                                          ! N of self higgs
c     >                839*-100 /                                  ! EOF signal
     >                929*-100 /                                  ! EOF signal

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
c         rst=result
C
         return
         end

      subroutine mom(lflag,np,pt0,pt1,eta0,xmr,en0,pz0,
     .               pm,wgt,ranram,lw)
      implicit none
      real*8 pi,wgt,rpr,phr,ran0,wj,en0,pz0,prx,pry
      real*8 etai1,etai2,etaj1,etaj2,xmri2,bxmri2,wjeta
      real*8 cn,wjaci,wjacj,phi,ranj,ref,alimp,alimm
      real*8 p0m,p0p 
      real*8 phj,alim,phip,phim,etap,etam,ares
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
      real*8 asinh,aus,betp,alpp
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
          ares = pm(2,i)/pt(i)
          if (ares.gt. 1.d0) ares= 1.d0
          if (ares.lt.-1.d0) ares=-1.d0
          ph(i)= acos(ares) 
          if (pm(1,i).lt.0.d0) ph(i)= 2.d0*pi-ph(i)
          if (iter.eq.n) then
            ares =  pm(2,j)/pt(j)
            if (ares.gt. 1.d0) ares= 1.d0
            if (ares.lt.-1.d0) ares=-1.d0
            ph(j)= acos(ares) 
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
        etai1= asinh(0.5d0*(x1m-a/x1m))
        etai2= asinh(0.5d0*(x2m-a/x2m))
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
          etaj1= asinh((pz-pt(i)*sinh(etai1))/q)
          etaj2= asinh((pz-pt(i)*sinh(etai2))/q)
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


