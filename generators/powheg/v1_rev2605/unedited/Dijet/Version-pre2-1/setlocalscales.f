      subroutine setlocalscales(iuborn,imode,rescfac)
c returns the rescaling factor including sudakov form factors and
c coupling rescaling, for born (imode=1) and NLO corrections (imode=2)
      implicit none
      integer iuborn,imode
      real * 8 rescfac
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      include 'pwhg_st.h'
      integer j,k,mu
      logical samekin,sameflv
      integer flav(nlegborn)
      real * 8 op(0:3,nlegborn),basicfac,bornfac,nlofac
      real * 8 savbasicfac(maxprocborn),savbornfac(maxprocborn),
     1     savnlofac(maxprocborn),savmuf2(maxprocborn)
      logical valid(maxprocborn)
      save savbasicfac,savbornfac,savnlofac,savmuf2,valid
      save op
      data valid/maxprocborn*.false./
      do j=1,nlegborn
         do mu=0,3
            if(op(mu,j).ne.kn_cmpborn(mu,j)) then
               do k=1,flst_nborn
                  valid(k)=.false.
               enddo
               op=kn_cmpborn
            endif
         enddo
      enddo
 33   continue
      if(valid(iuborn)) then
         basicfac=savbasicfac(iuborn)
         nlofac=savnlofac(iuborn)
         bornfac=savbornfac(iuborn)
         st_mufact2=savmuf2(iuborn)
      else
         flav=flst_born(:,iuborn)
         call setlocalscales0(flav,kn_cmpborn,
     1        basicfac,bornfac,nlofac)
         savbasicfac(iuborn)=basicfac
         savnlofac(iuborn)=nlofac
         savbornfac(iuborn)=bornfac
         savmuf2(iuborn)=st_mufact2
         valid(iuborn)=.true.
      endif
      if(imode.eq.1) then
         rescfac=basicfac*bornfac
      elseif(imode.eq.2) then
         rescfac=basicfac*nlofac
      endif
      end

      subroutine setlocalscales0
     1     (flav,pin,basicfac,bornfac,nlofac)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      integer flav(nlegborn)
      real * 8 pin(0:3,nlegborn),basicfac,bornfac,nlofac
      integer onem
      parameter (onem=1000000)
      real * 8 scales(nlegborn),p(0:3,nlegborn),ptot(0:3),
     1     lscalej,lscalek
      integer j,k,l,lflav(nlegborn),jmerge,kmerge,inlofac
      integer mergedfl
      real * 8 q2merge,q2merge0,renfac2,facfact2,alphas,mu2,muf2
      real * 8 sudakov,expsudakov,pwhg_alphas,b0,powheginput
      external sudakov,expsudakov,pwhg_alphas,powheginput
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
      renfac2=st_renfact**2
      facfact2=st_facfact**2
      lflav=flav
      p=pin
      scales=0
      q2mergeMAX=-1d10
      do l=1,nlegborn
         call findNearestNeighbours(p,lflav,jmerge,kmerge,mergedfl,
     $                              q2merge)
         if(q2merge.lt.1d10) then
c     perform the merging
            if(q2merge.gt.q2mergeMAX) q2mergeMAX=q2merge
            lscalej=scales(jmerge)
            lscalek=scales(kmerge)
            scales(jmerge)=q2merge
            if(lscalej.eq.0) then
c     This is the first merge; it sets the low scale for
c     all partons; no Sudakov factor or reweighting is introduced
               do j=1,nlegborn
                  scales(j)=q2merge
               enddo
c save this scale; it is the Q_0 scale that appears in all Sudakovs
               q2merge0=q2merge
               bornfac=0
c     Provide alpha_S reweighting for the first merge
               alphas=pwhg_alphas(max(q2merge*renfac2,1d0),
     1              st_lambda5MSB,st_nlight)
               basicfac=alphas/st_alpha
               nlofac=basicfac
               mu2=max(q2merge*renfac2,1d0)
               muf2=max(q2merge*facfact2,1d0)
               inlofac=1
            else
c provide Sudakov
               basicfac=basicfac*
     1              sudakov(q2merge0,q2merge,lscalej,lflav(jmerge))
               basicfac=basicfac*
     1              sudakov(q2merge0,q2merge,lscalek,lflav(kmerge))
               bornfac=bornfac+
     1              expsudakov(q2merge0,q2merge,lscalej,lflav(jmerge))
               bornfac=bornfac+
     1              expsudakov(q2merge0,q2merge,lscalek,lflav(kmerge))
c provide alpha_S reweighting
               alphas=pwhg_alphas(max(q2merge*renfac2,1d0),
     1              st_lambda5MSB,st_nlight)
               basicfac=basicfac*alphas/st_alpha
               mu2=mu2*max(q2merge*renfac2,1d0)
               nlofac=nlofac+alphas/st_alpha
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
      do j=3,nlegborn
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
         do j=1,nlegborn
            if(abs(lflav(j)).le.st_nlight) then
               basicfac=basicfac*
     1              sudakov(q2merge0,q2merge,scales(j),lflav(j))
               bornfac=bornfac+
     1              expsudakov(q2merge0,q2merge,scales(j),lflav(j))
            endif
         enddo
      else
c If scales(1)=0 no merge has taken place: no sudakovs.
         mu2=1
         muf2=max(q2merge*facfact2,1d0)
         inlofac=0
         bornfac=0
         basicfac=1
         nlofac=0
      endif
      if(st_bornorder.gt.inlofac) then
         alphas=pwhg_alphas(max(q2merge*renfac2,1d0),
     1           st_lambda5MSB,st_nlight)
c         do j=inlofac+1,st_bornorder
c            mu2=mu2*max(q2merge*renfac2,1d0)
c            nlofac=nlofac+alphas/st_alpha
c            basicfac=basicfac*alphas/st_alpha
c         enddo
         mu2=mu2*max(q2merge*renfac2,1d0)**(st_bornorder-inlofac)
         nlofac=nlofac+alphas/st_alpha*(st_bornorder-inlofac)
         basicfac=basicfac*(alphas/st_alpha)**(st_bornorder-inlofac)
         inlofac=st_bornorder
      endif
      nlofac=nlofac/inlofac
      mu2=mu2**(1d0/inlofac)
      b0=(33-2*st_nlight)/(12*pi)
      bornfac=1+st_alpha*nlofac*
     1     (bornfac+st_bornorder*b0*log(mu2/st_muren2))
      st_mufact2=muf2
      end

C ------------------------------------------------ C
C - Inputs:                                      - C
C - *******                                      - C
C - p        - Underlying born momenta           - C
C - lflav    - Flavour list derived from         - C
C -            flst_born by subjecting it to     - C
C -            repeated QCD clusterings.         - C
C -                                              - C
C - Outputs:                                     - C
C - ********                                     - C
C - jmerge   - Index in lflav of one of the two  - C
C -            closest partons.                  - C
C - kmerge   - Index in lflav of the             - C
C -            corresponding parton.             - C
C - mergedfl - Flavour of parton resulting from  - C 
C -            combination.                      - C
C - q2merge  - pT^2 scale associated to the      - C
C -            merging of jmerge and kmerge.     - C
C -                                              - C
C - checked 24/03/12                             - C
C ------------------------------------------------ C
      subroutine findNearestNeighbours(p,lflav,jmerge,kmerge,mergedfl,
     $                                 q2merge)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
C - Input / output:
      real * 8 p(0:3,nlegborn)
      integer  lflav(nlegborn)
      integer  jmerge,kmerge,mergedfl
      real * 8 q2merge
C - Local variables:
      real * 8 ycm
      integer  onem
      parameter (onem=1000000)
      integer  j,k
      integer  fl1,fl2,fl
      real * 8 yj,phij,q2j
      real * 8 yk,phik,q2k
      real * 8 dphi
      real * 8 q2

      q2merge=1d10
      ycm=log(p(0,1)/p(0,2))/2
      mergedfl=onem
      do j=3,nlegborn
         if(abs(lflav(j)).gt.st_nlight) goto 11
         yj=0.5d0*log((p(0,j)+p(3,j))/(p(0,j)-p(3,j)))
         if(yj.gt.ycm) then
            call validmergeisr(lflav,1,j,fl1)
            if(fl1.ne.onem) then
               q2j = p(1,j)**2+p(2,j)**2
               if(q2j.lt.q2merge) then
                  q2merge=q2j
                  jmerge=1
                  kmerge=j
                  mergedfl=fl1
               endif
            endif
         else
            call validmergeisr(lflav,2,j,fl2)
            if(fl2.ne.onem) then
               q2j = p(1,j)**2+p(2,j)**2
               if(q2j.lt.q2merge) then
                  q2merge=q2j
                  jmerge=2
                  kmerge=j
                  mergedfl=fl2
               endif
            endif
         endif
         do k=j+1,nlegborn
            if(abs(lflav(k)).gt.st_nlight) goto 12
            call validmergefsr(lflav,j,k,fl)
            if(fl.ne.onem) then
               yk=0.5d0*log((p(0,k)+p(3,k))/(p(0,k)-p(3,k)))
               call phipt2(p(:,k),phik,q2k)
               call phipt2(p(:,j),phij,q2j)
               dphi=abs(phik-phij)
               if(dphi.gt.2*pi) dphi=dphi-2*pi
               if(dphi.gt.pi) dphi=2*pi-dphi
               q2=((yk-yj)**2+dphi**2)*min(q2k,q2j)
               if(q2.lt.q2merge) then
                  q2merge=q2
                  jmerge=j
                  kmerge=k
                  mergedfl=fl
               endif
            endif
 12         continue
         enddo
 11      continue
      enddo
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
C - sudakov - The Sudakov form factor.           - C
C -                                              - C
C ------------------------------------------------ C
      function sudakov(q20,q2h,q2l,flav)
      implicit none
      real * 8 sudakov,q2h,q2l,q20
      integer flav
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 lam2
      logical isQuark
      real * 8 theExponentN,theExponentD
      logical ini
      data ini/.true./
      save ini
      if(ini) then
         call sudakov_plotter
         ini=.false.
      endif
      lam2=st_lambda5MSB**2
      if(q20.le.lam2.or.q2l.lt.lam2.or.q2h.lt.lam2) then
         sudakov=0
         goto 999
      endif
      if(q2l.ge.q2h.or.q2h.le.q20) then
         sudakov=1
         goto 999
      endif
      if(flav.eq.0) then
         isQuark=.false.
      else
         isQuark=.true.
      endif
      if(q2l.le.q20) then
        call sudakov_exponent(q20,q2h,q2h,theExponentN,
     $                         isQuark,2,.true.)
        sudakov=exp(theExponentN)
      else
         call sudakov_exponent(q20,q2h,q2h,theExponentN,
     $                         isQuark,2,.true.)
         call sudakov_exponent(q20,q2l,q2l,theExponentD,
     $                         isQuark,2,.true.)
         sudakov=exp(theExponentN-theExponentD)
      endif
 999  continue
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
C - sudakov - The Sudakov form factor.           - C
C -           Uses 1-loop alpha, no CMW scheme   - C
C -           etc.                               - C
C -                                              - C
C ------------------------------------------------ C
      function LO_sudakov(q20,q2h,q2l,flav)
      implicit none
      real * 8 LO_sudakov,q2h,q2l,q20
      integer flav
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 b0,c,b,lam2
      lam2=st_lambda5MSB**2
      if(q20.le.lam2.or.q2l.lt.lam2.or.q2h.lt.lam2) then
         LO_sudakov=0
         goto 999
      endif
      if(q2l.ge.q2h.or.q2h.le.q20) then
         LO_sudakov=1
         goto 999
      endif
      b0=(33-2*st_nlight)/12d0
      if(flav.eq.0) then
         c=3
         b=b0/3
      else
         c=4d0/3
         b=3d0/4
      endif
      if(q2l.le.q20) then
         LO_sudakov= exp(
     1        -c/b0*(  log(log(q2h/lam2)/log(q20/lam2))
     2        *(0.5d0*log(q2h/lam2)-b)
     3        -0.5d0*log(q2h/q20)
     4        ))
      else
         LO_sudakov= exp(
     1        -c/b0*(  log(log(q2h/lam2)/log(q20/lam2))
     2        *(0.5d0*log(q2h/lam2)-b)
     3        -0.5d0*log(q2h/q20)
     4        )
     5        +c/b0*(  log(log(q2l/lam2)/log(q20/lam2))
     6        *(0.5d0*log(q2l/lam2)-b)
     7        -0.5d0*log(q2l/q20)
     8        ))
      endif
 999  continue
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
      function expsudakov(q20,q2h,q2l,flav)
      implicit none
      real * 8 expsudakov,q2h,q2l,q20
      integer flav
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 b0,c,b,lam2
      lam2=st_lambda5MSB**2
      if(q20.le.lam2.or.q2l.lt.lam2.or.q2h.lt.lam2) then
c in this case everything is zero, irrelevant
         expsudakov=0
         return
      endif
      if(q2l.ge.q2h.or.q2h.le.q20.or.flg_bornonly) then
         expsudakov=0
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
      if(q2l.le.q20) then
         expsudakov=
     1        c/pi*(0.25d0*log(q2h/q20)**2 - b*log(q2h/q20))
      else
         expsudakov=
     1        c/pi*(0.25d0*log(q2h/q20)**2 - b*log(q2h/q20))
     2      - c/pi*(0.25d0*log(q2l/q20)**2 - b*log(q2l/q20))
      endif
      end

C ---------------------------------------------- C
C - Inputs:                                    - C
C - *******                                    - C
C - flav - flavour list derived from flst_born - C
C -        by subjecting it to repeated QCD    - C
C -        compatible clusterings.             - C
C -  i   - index of i-th initial-state partON  - C
C -        in flav: hence i = 1 or 2 only.     - C
C -  j   - index of j-th final-state partICLE  - C
C -        particle in flav.                   - C
C -                                            - C
C - Outputs:                                   - C
C - ********                                   - C
C - fl   - Would-be PDG code of spacelike      - C
C -        "mother" parton obtained by merging - C
C -        (~on-shell) incoming parton i with  - C
C -        outgoing particle j:                - C
C -          i -> fl + j                       - C
C -        Note gluons have id=0 in Powheg-Box - C
C -        instead of 21. If the splitting is  - C
C -        not possible in QCD, fl=1000000 ;   - C
C -        this setting signals to the rest of - C
C -        the algorithm that this is not a    - C
C -        candidate pair for combination.     - C
C -                                            - C
C ---------------------------------------------- C
      subroutine validmergeisr(flav,i,j,fl)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      integer onem
      parameter (onem=1000000)
      integer flav(nlegborn),i,j,fl
      integer lflav(nlegborn)
      logical validflav
      external validflav
      if(i.gt.2.or.j.le.2) then  ! Remove when development is finished.
         write(*,*) 'validmergeisr: fatal error'
         write(*,*) 'Routine demands an i.s. and f.s. particle'
         write(*,*) 'index for the 2nd and 3rd input values   '
         write(*,*) 'respectively. Quitting.'
         call exit(-1)
      endif
      if(abs(flav(i)).gt.st_nlight.or.abs(flav(j)).gt.st_nlight) then
         fl=onem
         return
      endif
      if(flav(i).eq.flav(j)) then
c g -> g g or q -> g q
         fl=0
         goto 999
      endif
      if(flav(j).eq.0) then
c q -> q g
         fl=flav(i)
         goto 999
      endif
      if(flav(i).eq.0) then
c g -> qbar q
         fl=-flav(j)
         goto 999
      endif
      fl=onem
      return
 999  continue
C - Check that the flavour list that results from the merging
C - is acceptable e.g. check that for HJJ you don't get back to
C - qqbar->H; if you do then set fl to 1000000, as if the 
C - branching were not possible in QCD s.t. it will be neglected
C - as a candidate for clustering.
      lflav=flav
      lflav(j)=onem
      lflav(i)=fl
      if(.not.validflav(lflav)) then
         fl=onem
      endif
      end


C ---------------------------------------------- C
C - Inputs:                                    - C
C - *******                                    - C
C - flav - flavour list derived from flst_born - C
C -        by subjecting it to repeated QCD    - C
C -        compatible clusterings.             - C
C -  i   - index of i-th final-state partICLE  - C
C -        in flav: hence i = 1 or 2 only.     - C
C -  j   - index of j-th final-state partICLE  - C
C -        particle in flav.                   - C
C -                                            - C
C - Outputs:                                   - C
C - ********                                   - C
C - fl   - Would-be PDG code of timelike       - C
C -        "mother" parton obtained by merging - C
C -        outgoing particles i and j:         - C
C -          fl -> i + j                       - C
C -        Note gluons have id=0 in Powheg-Box - C
C -        instead of 21. If the splitting is  - C
C -        not possible in QCD, fl=1000000 ;   - C
C -        this setting signals to the rest of - C
C -        the algorithm that this is not a    - C
C -        candidate pair for combination.     - C
C -                                            - C
C ---------------------------------------------- C
      subroutine validmergefsr(flav,i,j,fl)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      integer onem
      parameter (onem=1000000)
      integer flav(nlegborn),i,j,fl
      integer lflav(nlegborn)
      logical validflav
      external validflav
      if(i.le.2.or.j.le.2) then  ! Remove when development is finished.
         write(*,*) 'validmergefsr: fatal error'
         write(*,*) 'Routine demands an f.s. and f.s. particle'
         write(*,*) 'index for the 2nd and 3rd input values   '
         write(*,*) 'respectively. Quitting.'
         call exit(-1)
      endif
      if(abs(flav(i)).gt.st_nlight.or.abs(flav(j)).gt.st_nlight) then
         fl=onem
         return
      endif
      if(flav(i).eq.-flav(j)) then
c g -> g g or g -> q qbar
         fl=0
         goto 999
      endif
      if(flav(j).eq.0) then
c q -> q g
         fl=flav(i)
         goto 999
      endif
      if(flav(i).eq.0) then
c q -> g q
         fl=flav(j)
         goto 999
      endif
      fl=onem
      return
 999  continue
C - Check that the flavour list that results from the merging
C - is acceptable e.g. check that for HJJ you don't get back to
C - qqbar->H; if you do then set fl to 1000000, as if the 
C - branching were not possible in QCD s.t. it will be neglected
C - as a candidate for clustering.
      lflav=flav
      lflav(j)=onem
      lflav(i)=fl
      if(.not.validflav(lflav)) then
         fl=onem
      endif
      end


C ---------------------------------------------- C
C - Inputs:                                    - C
C - *******                                    - C
C - p    - p(0) = Energy, p(3) = p_Z           - C
C -                                            - C
C - Outputs:                                   - C
C - ********                                   - C
C - y    - Rapidity                            - C
C - phi  - phi                                 - C
C - q2   - pT^2 w.r.t the beam                 - C
C -                                            - C
C ---------------------------------------------- C
      subroutine phipt2(p,phi,q2)
      implicit none
      real * 8 p(0:3),phi,q2
      q2=p(1)**2+p(2)**2
      phi=atan2(p(2),p(1))
      end

C ********* DDT / Ellis-Veseli / Nason-Ridolfi Sudakov ************ C
C -                                                               - C
C - Output:                                                       - C
C - ========                                                      - C
C - The value of the Sudakov exponent defined as the integral,    - C
C - from  Log [ ql^2/Lambda^2 ]  up to Log [  qh^2/Lambda^2 ], of - C
C -                                                               - C
C -    d Log[ q^2/Lambda^2 ]                                      - C
C -  - {                                                          - C
C -      aSBar*A1*Log[m^2/q^2] + aSBar^2*A2*Log[m^2/q^2]          - C
C -    + aSBar*B1              + aSBar^2*B2                       - C
C -    }                                                          - C
C -                                                               - C
C - where aSBar = aS/2/Pi.                                        - C
C -                                                               - C
C - For m2=qh2, except for an overall factor of two this is the   - C
C - Sudakov form factor of eq. 32 in the Ellis-Veseli paper - in  - C
C - that paper they have two quark lines to consider while here   - C
C - we only want to consider one line at a time. The factor of    - C
C - two is manifest in the code below as our A1, A2, B1, B2       - C
C - coefficients defined to be HALF of the Ellis-Veseli ones.     - C
C -                                                               - C
C - The m2 dependence is a relic of the Nason & Ridolfi form of   - C
C - the Sudakov form factor, which has the numerator in the large - C
C - log equal to mZZ but the upper bound on the Sudakov integral  - C
C - is Q^2. There doesn't seem to be any problem arising when     - C
C - you just call the routine with m^2=Q^2, but maybe if this     - C
C - gets resolved I can re-do the mathematica integral.           - C
C -                                                               - C
C - The analytic integral was done in Mathematica assuming no     - C
C - flavour thresholds. When q^2 is below the b or c quark        - C 
C - flavour thresholds a numerical integration is done instead    - C
C - using dgauss. The numerical integration and analytic results  - C
C - agree very well above these thresholds - try resetting        - C
C - debuggingEpsilon below.                                       - C
C -                                                               - C 
C - For the default values of A1, A2, B1, B2 in the code below    - C
C - the Sudakov should correspond to that of Nason and Ridolfi,   - C
C - which has an effective B2 term by virtue of the fact that     - C
C - the CMW alpha_S is used to multiply the leading & subleading  - C
C - term. At least with the calculation done in the way it is we  - C
C - can easily play around with the coefficients.                 - C 
C -                                                               - C 
C - To Use:                                                       - C
C - =======                                                       - C
C - q2l = The scale of the lower  clustered node.                 - C
C - q2h = The scale of the higher clustered node.                 - C
C - m2  = The boson mass squared (argument of the                 - C
C -       log in the exponent of N.R. eq 4.8).                    - C
C - theExponent                                                   - C 
C -     = The value of the curly brackets in N.R. eq 4.8.         - C
C - isQuark = .true. for a quark propagator                       - C 
C - theAccuracy = 0 for 1-loop alphaS and A2=B2=0,                - C 
C -             = 1 for 2-loop alphaS and Powheg A & B coeffs     - C 
C -             = 2 for 2-loop alphaS and NLL A & B coeffs        - C 
C -                                                               - C 
C - Notes:                                                        - C
C - ======                                                        - C
C - Details for the integration in the sudakov exponent can be    - C
C - found in the Mathematica notebook: menlops/DDT_exponent.nb .  - C 
C -                                                               - C 
C - The Mathematica notebook shows plots in which the 5-flavour   - C
C - and 4-flavour alphaS differ by <2% at pT=3 GeV and 4% at      - C
C - pT=2 GeV. Using the C.M.W. alphaS (aS -> aS*(1+aS*K/2*pi))    - C
C - increases these differences but they remain small: 2% at pT=4 - C
C - GeV and 6% at pT=2 GeV. Note well that since the program      - C
C - matches alphaS at flavour thresholds, not                     - C
C - alphaS*(1+alphaS*K/2*pi), since K too actually depends on     - C
C - the number of flavours, this means the 3,4 and 5 flavour      - C
C - C.M.W. alphaS*(1+alphaS*K/2*pi) DO NOT match at the flavour   - C
C - thresholds in pT! Whereas alphaS 4 and 5 flavour couplings    - C
C - match at 5 GeV, the nf dependence of K means that the 4 and 5 - C
C - flavour alphaS*(1+alphaS*K/2*pi) actually meet at about 9 GeV - C
C - instead.                                                      - C
C -                                                               - C
C ***************************************************************** C
      subroutine sudakov_exponent(q2l,q2h,m2,theExponent,isQuark,
     $                            theAccuracy,fixed_nf)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_math.h'
      real * 8 q2l,q2h,m2,theExponent
      logical  isQuark
      integer  theAccuracy
      integer  nf
      real * 8 bnf,bpnf,K
      real * 8 A1,B1,A2,B2,zeta3
      real * 8 Lq2l,Lq2h,Lm2
      real * 8 pwhg_alphas
      external pwhg_alphas
      real * 8 aSbar
      real * 8 theA1coeff,theA2coeff,theB1coeff,theB2coeff
      real * 8 eps
      logical  isQuarkLine
      integer  accuracy
      real * 8 m2_common
      common/sudakov_integral/isQuarkLine,accuracy,m2_common
      real * 8 dgauss,sudakov_exponent_integrand
      external dgauss,sudakov_exponent_integrand
      logical  pwhg_isfinite
      external pwhg_isfinite
      real * 8 debuggingEpsilon,tmp

      logical  fixed_nf ! Assume the number of active flavours in 
                        ! is always st_nlight when evaluating the
                        ! Sudakov.

C - Fractional difference (%) between analytic and numerical Sudakov
C - exponent integration which leads to an error on the screen if
C - it is exceeded. Making it negative (recommended) deactivates
C - this debugging.
      debuggingEpsilon = -999d0 !  1d-9

      isQuarkLine  = isQuark

      accuracy = theAccuracy

      if(fixed_nf) then
         nf=st_nlight
      else
         if(q2l.lt.rad_charmthr2) then
            nf=3
         elseif(q2l.lt.rad_bottomthr2) then
            nf=4
         else
            nf=5
         endif
      endif

      bnf  = (11d0*CA-2d0*nf)/12/Pi
      bpnf = (153 - 19d0*nf) / Pi / 2 / (33 - 2*nf)

      K  =  (67d0/18-Pi**2/6)*CA-5d0/9*nf
      if(isQuarkLine) then
         A1 =   Cf
         A2 =   Cf*K
         B1 =  -3d0/2*Cf
C - Powheg spurious B2:
         B2 =  -3d0/2*Cf*K
C - Davies and Stirling B2 (divided by two):
C - TO DO: adjust by Ellis-Veseli qT-space prescription for B2.
         zeta3 = 1.20206
c$$$         B2    = (     Cf*Cf * (Pi*Pi-3d0/4-12*zeta3)
c$$$     $               + Cf*CA * (11*Pi*Pi/9-193d0/12+6*zeta3)
c$$$     $           + Cf*nf*0.5 * (17d0/3-4d0/9*Pi*Pi)
c$$$     $           ) / 2d0
      else
         A1 =  CA
         A2 =  CA*K
         B1 = -2*Pi*bnf
C - Powheg spurious B2:
         B2 = -2*Pi*bnf*K
C - Grazzini and De Florian B2 (divided by two):
C - TO DO: adjust by Ellis-Veseli qT-space prescription for B2.
         zeta3 = 1.20206
c$$$         B2    = (       CA*CA * (23d0/6+22*Pi*Pi/9-6*zeta3)
c$$$     $               + 4*Cf*nf*0.5d0
c$$$     $               - CA*nf*0.5*(2d0/3+8d0*Pi*Pi/9)
c$$$     $               - 11d0/2*CA*Cf
c$$$     $           ) / 2d0
      endif

      Lq2l = Log(q2l/st_lambda5MSB/st_lambda5MSB)
      Lq2h = Log(q2h/st_lambda5MSB/st_lambda5MSB)
      Lm2  = Log(m2/st_lambda5MSB/st_lambda5MSB)

      aSbar = pwhg_alphas(q2l,st_lambda5MSB,-1)/2/Pi

c$$$C - We insist the DDT/Ellis-Veseli Sudakov exponent is positive
c$$$C - everywhere. This effectively defines the upper bound of the
c$$$C - Sudakov exponent integral.
c$$$      if( ( (A1+A2*aSbar)*(Lm2-Lq2l) + (B1+B2*aSbar) ) .lt.0) then
c$$$         theExponent=0d0
c$$$         return
c$$$      endif

      if(accuracy.eq.0) then
         A2   = 0d0 ! NLL coefficient
         B2   = 0d0 ! NNLL coefficient
         bpnf = 0d0 ! Reduce 2-loop to 1-loop alpha in the calculation
      else if(accuracy.eq.2) then
         B2   = 0d0 ! NNLL coefficient
      else if(accuracy.eq.3) then
         B1   = 0d0 ! NLL coefficient
         A2   = 0d0 ! NLL coefficient
         B2   = 0d0 ! NNLL coefficient
         bpnf = 0d0 ! Reduce 2-loop to 1-loop alpha in the calculation
      endif

      if(q2l.ge.rad_bottomthr2.or.fixed_nf) then

         theA1coeff =
     $        ( (Lq2h - Lq2l)
     $           - Lm2*Log(Lq2h/Lq2l)
     $        )/(2*bnf*Pi)
     $      + bpnf*( 2*Lm2*(Lq2h - Lq2l)
     $             + 2*Lm2*Lq2h*Log(Lq2l)
     $             + Lq2l*Lq2h*Log(Lq2l)**2 
     $             - 2*Lm2*Lq2l*Log(Lq2h)
     $             - Lq2l*Lq2h*Log(Lq2h)**2
     $             )/(4*bnf**2*Lq2l*Lq2h*Pi)

         theA2coeff = 
     $         ( Lm2*(Lq2l - Lq2h)
     $         - Lq2l*Lq2h*Log(Lq2l)
     $         + Lq2l*Lq2h*Log(Lq2h)
     $         )/(4d0*bnf**2*Lq2l*Lq2h*Pi**2)
     $       + bpnf*( 0.5*(Lq2h-Lq2l)*( Lm2*(Lq2l + Lq2h)
     $                                - 4*Lq2l*Lq2h
     $                                )
     $              + (Lm2 - 2*Lq2l)*Lq2h**2*Log(Lq2l)
     $              - Lq2l**2*(Lm2 - 2*Lq2h)*Log(Lq2h)
     $              )/(4d0*bnf**3*Lq2l**2*Lq2h**2*Pi**2)
     $       + bpnf**2*(  (Lq2l - Lq2h)*( 8*Lm2*( Lq2l**2
     $                                          + Lq2l*Lq2h
     $                                          + Lq2h**2
     $                                          )
     $                                  - 27*Lq2l*Lq2h*(Lq2l + Lq2h)
     $                                  )
     $                 - 6*(4*Lm2 - 9*Lq2l)*Lq2h**3*Log(Lq2l)
     $                 + 18*(-2*Lm2 + 3*Lq2l)*Lq2h**3*Log(Lq2l)**2
     $                 + 6*Lq2l**3*(4*Lm2 - 9*Lq2h)*Log(Lq2h)
     $                 + 18*Lq2l**3*(2*Lm2 - 3*Lq2h)*Log(Lq2h)**2
     $                 )/(432.*bnf**4*Lq2l**3*Lq2h**3*Pi**2)

         theB1coeff = 
     $       - Log(Lq2h/Lq2l)/(2.*bnf*Pi)
     $       - bpnf*( Lq2l - Lq2h - Lq2h*Log(Lq2l)
     $              + Lq2l*Log(Lq2h)
     $              )/(2.*bnf**2*Lq2l*Lq2h*Pi) 

         theB2coeff =
     $       - (1/Lq2l - 1/Lq2h)/(4.*bnf**2*Pi**2)
     $       - bpnf*(   Lq2l**2 - Lq2h**2 - 2*Lq2h**2*Log(Lq2l)
     $              + 2*Lq2l**2*Log(Lq2h)
     $              )/(8.*bnf**3*Lq2l**2*Lq2h**2*Pi**2)
     $       - bpnf**2*( 2*Lq2h**3 - 2*Lq2l**3
     $                 + 3*Lq2h**3*Log(Lq2l)*(2 + 3*Log(Lq2l)) 
     $                 - 3*Lq2l**3*Log(Lq2h)*(2 + 3*Log(Lq2h))
     $                 )/(108.*bnf**4*Lq2l**3*Lq2h**3*Pi**2)

         theExponent = A1*theA1coeff + A2*theA2coeff
     $               + B1*theB1coeff + B2*theB2coeff

         if(debuggingEpsilon.gt.0d0) then
            tmp=dgauss(sudakov_exponent_integrand,Lq2l,Lq2h,eps)
            tmp=100d0*abs((tmp-theExponent)/(tmp+theExponent))
            if(tmp.gt.debuggingEpsilon) then
               write(*,*) ''
               write(*,*) 'sudakov_exponent: debug mode'
               write(*,*) '============================'
               write(*,*) 'sqrt(q2l) (GeV) = ',sqrt(q2l)
               write(*,*) 'sqrt(q2h) (GeV) = ',sqrt(q2h)
               write(*,*) 'sqrt(m2) (GeV)  = ',sqrt(m2)
               write(*,*) 'analytic - numerical exponent (%) = ',tmp
               write(*,*) 'analytic exponent = ',theExponent
               write(*,*) 'A1 term = ',A1*theA1coeff
               write(*,*) 'A2 term = ',A1*theA2coeff
               write(*,*) 'B1 term = ',A1*theB1coeff
               write(*,*) 'B2 term = ',A1*theB2coeff
            endif
         endif
      else
         eps = 1d-6
         m2_common=m2
         theExponent=dgauss(sudakov_exponent_integrand,Lq2l,Lq2h,eps)
      endif

c$$$      if(theExponent.gt.0d0) theExponent=0d0;

      if(.not.pwhg_isfinite(theExponent)) then
         write(6,*) ' '
         write(6,*) 'Warning: sudakov_exponent is weird.'
         write(6,*) 'theExponent = ',theExponent
         write(6,*) 'exp(theExponent) = ',exp(theExponent)
         write(6,*) 'q_low   = ',sqrt(q2l)
         write(6,*) 'q_hi    = ',sqrt(q2h)
         write(6,*) 'm       = ',sqrt(m2)
      endif

      end


C ***************************************************************** C
C - The integrand in Sudakov exponent times q^2 i.e. we effectiv  - C
C - ely integrate in log(q2/lambda^2) (Lq2) to better sample      - C
C - the Sudakov peak region. This can be used to replace the      - C
C - above analytic result in the region below the b- and c-quark  - C
C - thresholds when a variable flavour number should come into    - C
C - effect (not taken account of in the analytic result!). We see - C
C - excellent agreement above the b-quark threshold between       - C
C - numerical and analytic results (basically exact agreement),   - C
C - while below pT=5 GeV we see differences of at most 5%, rising - C
C - to 10% and 50% as pT goes below the charm threshold (~2 GeV). - C
C - I suggest using this routine AS AN ENHANCEMENT of the one     - C
C - above, i.e. to call the one here with a low value of EPS for  - C
C - events in which pT < b-quark threshold.                       - C
C ***************************************************************** C
      function sudakov_exponent_integrand(Lq2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_math.h'
      real * 8 sudakov_exponent_integrand
      real * 8 Lq2
      logical  isQuarkLine
      integer  accuracy
      real * 8 m2_common
      common/sudakov_integral/isQuarkLine,accuracy,m2_common
      real * 8 q2,aSbar
      integer  nf
      real * 8 bnf,bpnf,K
      real * 8 A1,B1,A2,B2,zeta3
      real * 8 pwhg_alphas,pwhg_alphas0
      external pwhg_alphas,pwhg_alphas0

      q2=st_lambda5MSB**2*exp(Lq2)

      if(q2.lt.rad_charmthr2) then
         nf=3
      elseif(q2.lt.rad_bottomthr2) then
         nf=4
      else
         nf=5
      endif

      if(accuracy.eq.0) then    ! Do 1-loop running AND keep nf=5
         nf = 5
         aSbar = pwhg_alphas0(q2,st_lambda5MSB,nf)/2/Pi
      elseif(accuracy.eq.3) then! Do 1-loop running AND keep nf=5
         nf = 5
         aSbar = pwhg_alphas0(q2,st_lambda5MSB,nf)/2/Pi
      else                     ! Do 2-loop running
         aSbar = pwhg_alphas(q2,st_lambda5MSB,-1)/2/Pi
      endif

      bnf  = (11d0*CA-2d0*nf)/12/Pi
      bpnf = (153 - 19d0*nf) / Pi / 2 / (33 - 2*nf)

      K  =  (67d0/18-Pi**2/6)*CA-5d0/9*nf
      if(isQuarkLine) then
         A1 =   Cf
         A2 =   Cf*K
         B1 =  -3d0/2*Cf
C - Powheg spruious B2:
         B2 =  -3d0/2*Cf*K
C - Davies and Stirling B2 (divided by two):
C - TO DO: adjust by Ellis-Veseli qT-space prescription for B2.
         zeta3 = 1.20206
c$$$         B2    = (     Cf*Cf * (Pi*Pi-3d0/4-12*zeta3)
c$$$     $               + Cf*CA * (11*Pi*Pi/9-193d0/12+6*zeta3)
c$$$     $           + Cf*nf*0.5 * (17d0/3-4d0/9*Pi*Pi)
c$$$     $           ) / 2d0
      else
         A1 =  CA
         A2 =  CA*K
         B1 = -2*Pi*bnf
C - Powheg spruious B2:
         B2 = -2*Pi*bnf*K
C - Grazzini and De Florian B2 (divided by two):
C - TO DO: adjust by Ellis-Veseli qT-space prescription for B2.
         zeta3 = 1.20206
c$$$         B2    = (       CA*CA * (23d0/6+22*Pi*Pi/9-6*zeta3)
c$$$     $               + 4*Cf*nf*0.5d0
c$$$     $               - CA*nf*0.5*(2d0/3+8d0*Pi*Pi/9)
c$$$     $               - 11d0/2*CA*Cf
c$$$     $           ) / 2d0
      endif

      if(accuracy.eq.0) then
         A2  = 0d0 ! NLL coefficient
         B2  = 0d0 ! NNLL coefficient
      else if(accuracy.eq.2) then
         B2  = 0d0 ! NNLL coefficient
      else if(accuracy.eq.3) then
         B1  = 0d0 ! NLL coefficient
         A2  = 0d0 ! NLL coefficient
         B2  = 0d0 ! NNLL coefficient
      endif

      sudakov_exponent_integrand=
     $     - ( (A1*aSbar+A2*aSbar*aSbar)*log(m2_common/q2)
     $       + (B1*aSbar+B2*aSbar*aSbar)
     $       )

      end

      subroutine sudakov_plotter()
      implicit none
      include 'pwhg_st.h'
      real * 8 the_q,the_Lq2,the_qh
      real * 8 the_Sudakov,one_minus_alpha
      real * 8 the_step,the_start,the_end
      integer  no_steps,ixx,iun
      logical  isQuark

      no_steps=1000

      the_qh=120d0

      call newunit(iun)
      open(unit=iun,file='Sudakov.top')

      write(iun,*) 'SET FONT DUPLEX'
      write(iun,*) 'SET WINDOW X 1 12 Y 2.5 9'
      write(iun,*) 'SET INTENSITY 5'
      write(iun,*) 'SET LABELS SIZE 2.5'
      write(iun,*) 'SET TITLE TOP SCALE 1.0'
      write(iun,*) 'SET TITLE SIZE 2.0'
      write(iun,*) 'TITLE TEXT 7.0 6.0 "quark A1,B1"'
      write(iun,*) 'CASE               "           "'
      write(iun,*) '5.5 6.0'
      write(iun,*) '6.5 6.0'
      write(iun,*) 'JOIN TEXT RED'
      write(iun,*) 'TITLE TEXT 7.0 5.5 "quark A1,B1,A2"'
      write(iun,*) 'CASE               "              "'
      write(iun,*) '5.5 5.5'
      write(iun,*) '6.5 5.5'
      write(iun,*) 'JOIN TEXT GREEN'
      write(iun,*) 'TITLE TEXT 7.0 5.0 "quark Powheg"'
      write(iun,*) 'CASE               "            "'
      write(iun,*) '5.5 5.0'
      write(iun,*) '6.5 5.0'
      write(iun,*) 'JOIN TEXT BLUE'
      write(iun,*) 'TITLE TEXT 7.0 4.5 "gluon A1,B1"'
      write(iun,*) 'CASE               "           "'
      write(iun,*) '5.5 4.5'
      write(iun,*) '6.5 4.5'
      write(iun,*) 'JOIN TEXT RED DASHES'
      write(iun,*) 'TITLE TEXT 7.0 4.0 "gluon A1,B1,A2"'
      write(iun,*) 'CASE               "              "'
      write(iun,*) '5.5 4.0'
      write(iun,*) '6.5 4.0'
      write(iun,*) 'JOIN TEXT GREEN DASHES'
      write(iun,*) 'TITLE TEXT 7.0 3.5 "gluon Powheg"'
      write(iun,*) 'CASE               "            "'
      write(iun,*) '5.5 3.5'
      write(iun,*) '6.5 3.5'
      write(iun,*) 'JOIN TEXT GREEN DASHES'
      write(iun,*) 'TITLE BOTTOM "p0T1 (GeV)"'
      write(iun,*) 'CASE         " X X S   S"'
      write(iun,*) 'SET ORDER X Y DY'

C - Plot Sudakov vs kT
      the_start= 0.5d0
      the_end  = 
     $     st_lambda5MSB
     $   * exp(0.5d0*(Log(the_qh**2/st_lambda5MSB**2)-1.5d0))
      the_step = (the_end-the_start)/dble(no_steps)

      isQuark = .true.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,0,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN RED'
      isQuark = .true.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,2,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN GREEN'
      isQuark = .true.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,1,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN BLUE'
      isQuark = .false.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,0,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN RED DASHES'
      isQuark = .false.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,2,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN GREEN DASHES'
      isQuark = .false.
      do ixx=1,no_steps
         the_q = the_start + ixx*the_step
         call sudakov_exponent(the_q**2,the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,1,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_q,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN BLUE DASHES'

      write(iun,*) ''
      write(iun,*) 'NEW FRAME'
      write(iun,*) 'SET FONT DUPLEX'
      write(iun,*) 'SET WINDOW X 1 12 Y 2.5 9'
      write(iun,*) 'SET INTENSITY 5'
      write(iun,*) 'SET LABELS SIZE 2.5'
      write(iun,*) 'SET TITLE TOP SCALE 1.0'
      write(iun,*) 'SET TITLE SIZE 2.0'
      write(iun,*) 'TITLE TEXT 3.0 8.0 "quark A1,B1"'
      write(iun,*) 'CASE               "           "'
      write(iun,*) '1.5 8.0'
      write(iun,*) '2.5 8.0'
      write(iun,*) 'JOIN TEXT RED'
      write(iun,*) 'TITLE TEXT 3.0 7.5 "quark A1,B1,A2"'
      write(iun,*) 'CASE               "              "'
      write(iun,*) '1.5 7.5'
      write(iun,*) '2.5 7.5'
      write(iun,*) 'JOIN TEXT GREEN'
      write(iun,*) 'TITLE TEXT 3.0 7.0 "quark Powheg"'
      write(iun,*) 'CASE               "            "'
      write(iun,*) '1.5 7.0'
      write(iun,*) '2.5 7.0'
      write(iun,*) 'JOIN TEXT BLUE'
      write(iun,*) 'TITLE TEXT 3.0 6.5 "gluon A1,B1"'
      write(iun,*) 'CASE               "           "'
      write(iun,*) '1.5 6.5'
      write(iun,*) '2.5 6.5'
      write(iun,*) 'JOIN TEXT RED DASHES'
      write(iun,*) 'TITLE TEXT 3.0 6.0 "gluon A1,B1,A2"'
      write(iun,*) 'CASE               "              "'
      write(iun,*) '1.5 6.0'
      write(iun,*) '2.5 6.0'
      write(iun,*) 'JOIN TEXT GREEN DASHES'
      write(iun,*) 'TITLE TEXT 3.0 5.5 "gluon Powheg"'
      write(iun,*) 'CASE               "            "'
      write(iun,*) '1.5 5.5'
      write(iun,*) '2.5 5.5'
      write(iun,*) 'JOIN TEXT BLUE DASHES'
      write(iun,*) 'TITLE BOTTOM "Log(p0T1223/L051223)"'
      write(iun,*) 'CASE         "   S X XX X FX XX XS"'
      write(iun,*) 'SET ORDER X Y DY'

C - Plot Sudakov vs Log(kT/Lambda)
      the_start= Log(0.5d0**2/st_lambda5MSB**2)
      the_end  = Log(the_qh**2/st_lambda5MSB**2)-1.5d0
      the_step = (the_end-the_start)/dble(no_steps)

      isQuark = .true.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                         the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,0,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN RED'
      isQuark = .true.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                         the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,2,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN GREEN'
      isQuark = .true.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                         the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,1,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN BLUE'
      isQuark = .false.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                         the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,0,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN RED DASHES'
      isQuark = .false.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                         the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,2,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN GREEN DASHES'
      isQuark = .false.
      do ixx=1,no_steps
         the_Lq2 = the_start + ixx*the_step
         call sudakov_exponent(exp(the_Lq2)*st_lambda5MSB**2,
     $                         the_qh**2,the_qh**2,
     $                         the_Sudakov,isQuark,1,.true.)
         the_Sudakov=exp(the_Sudakov)
         write(iun,'(F8.5,A,F8.6)') the_Lq2,'  ',the_Sudakov
      enddo
      write(iun,*) 'JOIN BLUE DASHES'

      close(iun)

      end 
