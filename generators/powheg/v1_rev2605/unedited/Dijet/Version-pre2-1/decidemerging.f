      subroutine getbornfac(nlegs,flav,pin,bornfac)
c this is a remodelled setlocalscales0 routine, that only
c returns the Born factor.
      implicit none
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      integer nlegs
      integer flav(nlegs)
      real * 8 pin(0:3,nlegs),bornfac,nlofac
      integer onem
      parameter (onem=1000000)
      real * 8 scales(nlegs),p(0:3,nlegs),ptot(0:3),
     1     lscalej,lscalek
      integer j,k,l,lflav(nlegs),jmerge,kmerge,inlofac
      integer mergedfl
      real * 8 q2merge,q2merge0,renfac2,facfact2,alphas,mu2,muf2
      real * 8 sudakov,expsudakov_mrg,pwhg_alphas,b0,powheginput
      external sudakov,expsudakov_mrg,pwhg_alphas,powheginput
      real * 8 q2mergeMAX
      logical raisingscales,ini
      save raisingscales,ini
      data ini/.true./
      if(ini) then
         if(powheginput("#raisingscales").eq.1) then
            raisingscales = .true.
         else
            raisingscales = .false.
         endif
         ini = .false.
      endif
      lflav=flav
      p=pin
      scales=0
      q2mergeMAX=-1d10
      do l=1,nlegs
         call findNearestNeighbours(nlegs,p,lflav,jmerge,kmerge,
     $        mergedfl,q2merge)
         if(q2merge.lt.1d10) then
c     perform the merging
            if(q2merge.gt.q2mergeMAX) q2mergeMAX=q2merge
            lscalej=scales(jmerge)
            lscalek=scales(kmerge)
            scales(jmerge)=q2merge
            if(lscalej.eq.0) then
c     This is the first merge; it sets the low scale for
c     all partons; no Sudakov factor or reweighting is introduced
               do j=1,nlegs
                  scales(j)=q2merge
               enddo
c save this scale; it is the Q_0 scale that appears in all Sudakovs
               q2merge0=q2merge
               bornfac=0
c     Provide alpha_S reweighting for the first merge
               alphas=pwhg_alphas(max(q2merge,1d0),
     1              st_lambda5MSB,st_nlight)
               nlofac=alphas
               inlofac=1
            else
               bornfac=bornfac+
     1            expsudakov_mrg(q2merge0,q2merge,lscalej,lflav(jmerge))
               bornfac=bornfac+
     1            expsudakov_mrg(q2merge0,q2merge,lscalek,lflav(kmerge))
c provide alpha_S reweighting
               alphas=pwhg_alphas(max(q2merge,1d0),
     1              st_lambda5MSB,st_nlight)
               nlofac=nlofac+alphas
               inlofac=inlofac+1
            endif
            if(jmerge.gt.2) then
               p(:,jmerge)=p(:,jmerge)+p(:,kmerge)
            else
               p(3,jmerge)=p(3,jmerge)-p(3,kmerge)
               p(0,jmerge)=abs(p(3,jmerge))
            endif
            lflav(kmerge)=onem
            lflav(jmerge)=mergedfl
         else
            goto 99
         endif
      enddo
 99   continue
c     No more merging is possible. Define the initial scale as
c     the invariant mass of the remaining system
      ptot=0
      do j=3,nlegs
         if(lflav(j).ne.onem) then
            ptot=ptot+p(:,j)
         endif
      enddo
      if(raisingscales) then
        q2merge=max(q2mergeMAX,
     $              ptot(0)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2)
      else
        q2merge=ptot(0)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2
      endif
      if(scales(1).gt.0) then
         do j=1,nlegs
            if(abs(lflav(j)).le.st_nlight) then
               bornfac=bornfac+
     1              expsudakov_mrg(q2merge0,q2merge,scales(j),lflav(j))
            endif
         enddo
      else
c If scales(1)=0 no merge has taken place: no sudakovs.
         inlofac=0
         bornfac=0
         nlofac=0
      endif
      if(st_bornorder.gt.inlofac) then
         alphas=pwhg_alphas(max(q2merge,1d0),
     1           st_lambda5MSB,st_nlight)
         nlofac=nlofac+alphas*(st_bornorder-inlofac)
         inlofac=st_bornorder
      endif
      if (inlofac /= 0) nlofac=nlofac/inlofac
      bornfac=1+nlofac*bornfac
      end

      subroutine getbornfac_1(iflag,nlegs,flav,pin,bornfac)
c this is a remodelled setlocalscales0 routine, that only
c returns the Born factor.
      implicit none
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      integer iflag,nlegs
      integer flav(nlegs)
      real * 8 pin(0:3,nlegs),bornfac,nlofac
      integer onem
      parameter (onem=1000000)
      real * 8 scales(nlegs),p(0:3,nlegs),ptot(0:3),
     1     lscalej,lscalek
      integer j,k,l,lflav(nlegs),jmerge,kmerge,inlofac
      integer mergedfl
      real * 8 q2merge,q2merge0,renfac2,facfact2,alphas,mu2,muf2
      real * 8 sudakov,expsudakov_mrg,pwhg_alphas,b0,powheginput
      external sudakov,expsudakov_mrg,pwhg_alphas,powheginput
      real * 8 q2mergeMAX
      logical raisingscales,ini
      save raisingscales,ini
      data ini/.true./
      if(ini) then
         if(powheginput("#raisingscales").eq.1) then
            raisingscales = .true.
         else
            raisingscales = .false.
         endif
         ini = .false.
      endif
      lflav=flav
      p=pin
      scales=0
      q2mergeMAX=-1d10
      do l=1,nlegs
         call findNearestNeighbours(nlegs,p,lflav,jmerge,kmerge,
     $        mergedfl,q2merge)
         if(q2merge.lt.1d10) then
c if iflag=1 and l=1, perform the first merging and cycle
            if(l.le.iflag) then
               if(jmerge.gt.2) then
                  p(:,jmerge)=p(:,jmerge)+p(:,kmerge)
               else
                  p(3,jmerge)=p(3,jmerge)-p(3,kmerge)
                  p(0,jmerge)=abs(p(3,jmerge))
               endif
               lflav(kmerge)=onem
               lflav(jmerge)=mergedfl
               cycle
            endif
c     perform the merging
            if(q2merge.gt.q2mergeMAX) q2mergeMAX=q2merge
            lscalej=scales(jmerge)
            lscalek=scales(kmerge)
            scales(jmerge)=q2merge
            if(lscalej.eq.0) then
c     This is the first merge; it sets the low scale for
c     all partons; no Sudakov factor or reweighting is introduced
               do j=1,nlegs
                  scales(j)=q2merge
               enddo
c save this scale; it is the Q_0 scale that appears in all Sudakovs
               q2merge0=q2merge
               bornfac=0
c     Provide alpha_S reweighting for the first merge
               alphas=pwhg_alphas(max(q2merge,1d0),
     1              st_lambda5MSB,st_nlight)
               nlofac=alphas
               inlofac=1

            else
               bornfac=bornfac+
     1            expsudakov_mrg(q2merge0,q2merge,lscalej,lflav(jmerge))
               bornfac=bornfac+
     1            expsudakov_mrg(q2merge0,q2merge,lscalek,lflav(kmerge))
c provide alpha_S reweighting
               alphas=pwhg_alphas(max(q2merge,1d0),
     1              st_lambda5MSB,st_nlight)
               nlofac=nlofac+alphas
               inlofac=inlofac+1
               bornfac=1 + bornfac*nlofac/inlofac
               return

            endif

            if(jmerge.gt.2) then
               p(:,jmerge)=p(:,jmerge)+p(:,kmerge)
            else
               p(3,jmerge)=p(3,jmerge)-p(3,kmerge)
               p(0,jmerge)=abs(p(3,jmerge))
            endif
            lflav(kmerge)=onem
            lflav(jmerge)=mergedfl
         else
            goto 99
         endif
      enddo
 99   continue
c     No more merging is possible. Define the initial scale as
c     the invariant mass of the remaining system
      ptot=0
      do j=3,nlegs
         if(lflav(j).ne.onem) then
            ptot=ptot+p(:,j)
         endif
      enddo
      if(raisingscales) then
        q2merge=max(q2mergeMAX,
     $              ptot(0)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2)
      else
        q2merge=ptot(0)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2
      endif
      if(scales(1).gt.0) then
         do j=1,nlegs
            if(abs(lflav(j)).le.st_nlight) then
               bornfac=bornfac+
     1              expsudakov_mrg(q2merge0,q2merge,scales(j),lflav(j))
            endif
         enddo
      else
c If scales(1)=0 no merge has taken place: no sudakovs.
         write(*,*) ' WARNING: getbornfac_1, no merging', iflag 
         write(*,*) 'lflav',lflav
         write(*,*) 'flav',flav
         write(*,*) 'q2merge',q2merge
         inlofac=0
         bornfac=1
         nlofac=0
      endif
      if (inlofac /= 0) nlofac=nlofac/inlofac
      bornfac=1+nlofac*bornfac
      end


C ------------------------------------------------ C
C - Inputs:                                      - C
C - *******                                      - C
C - q2h  - Upper node scale / bound on Sudakov   - C
C - q2l  - Lower node scale / bound on Sudakov   - C
C - flav - flavour index for the evolving parton - C
C -                                              - C
C - Outputs:                                     - C
C - ********                                     - C
C - expsudakov - The Sudakov form factor's expon - C
C -              -ent MODULO a factor of minus   - C
C -              alphaS, integrated with alphaS  - C
C -              fixed. Summed over with the     - C
C -              relevant alphaS factors this is - C
C -              used in compensating the NLO    - C
C -              correction induced when the     - C
C -              Sudakov multiplies the Born.    - C
C -                                              - C
C ------------------------------------------------ C
      function expsudakov_mrg(q20,q2h,q2l,flav)
      implicit none
      real * 8 expsudakov_mrg,q2h,q2l,q20
      integer flav
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 b0,c,b,lam2
      real * 8 powheginput
      logical flagllsudmerg,ini
      save flagllsudmerg,ini
      data ini/.true./
      if(ini) then
         if(powheginput("#flagllsudmerg").eq.1) then
            flagllsudmerg = .true.
            write(*,*) ' Using only LL sudakov for merging'
         else
            write(*,*) ' Using NLL sudakov for merging'
            flagllsudmerg = .false.
         endif
         ini = .false.
      endif
      lam2=st_lambda5MSB**2
      if(q20.le.lam2.or.q2l.lt.lam2.or.q2h.lt.lam2) then
c in this case everything is zero, irrelevant
         expsudakov_mrg=0
         return
      endif
      if(q2l.ge.q2h.or.q2h.le.q20.or.flg_bornonly) then
         expsudakov_mrg=0
         return
      endif
      b0=(33-2*st_nlight)/12d0
      if(flav.eq.0) then
         c=3
         b=b0/3
      else
         c=4d0/3
         b=3d0/4
      endif
      if(flagllsudmerg) b = 0d0
      if(q2l.le.q20) then
         expsudakov_mrg=
     1        c/pi*(0.25d0*log(q2h/q20)**2 - b*log(q2h/q20))
      else
         expsudakov_mrg=
     1        c/pi*(0.25d0*log(q2h/q20)**2 - b*log(q2h/q20))
     2      - c/pi*(0.25d0*log(q2l/q20)**2 - b*log(q2l/q20))
      endif
      end

