      subroutine setreal(p,rflav,result)
      implicit none
      include 'nlegborn.h'
c      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'constants.f'
      include 'plabel.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'mcfmtopwhg.f'
      integer nlegs
      parameter (nlegs=nlegreal)
c For t tbar with decay, nlegborn=12
      real * 8 p(0:3,nlegs), result
      integer rflav(nlegs),dkflag
c maxpart is in constants.f
      real * 8 mcfmp(mxpart,4)
      double precision msq(-nf:nf,-nf:nf)
      integer iem,iep,inu,inub,ib,ibb,ig,itop,iatop,iwp,iwm,fl,j,mu,
     & ro,myin(nlegs)
      logical isodd,isaquark
      real * 8 psav(0:3,nlegs)
      save psav,msq,dkflag
      real * 8 brcorrect
      external brcorrect
      if(dkflag.ne.kn_resemitter) goto 11
      do mu=0,3
         do j=1,nlegs
            if(p(mu,j).ne.psav(mu,j)) goto 11
         enddo
      enddo
      goto 22

 11   continue
      psav = p
      dkflag = kn_resemitter
c set scale in mcfm blocks
      scale=sqrt(st_muren2)
      musq=st_muren2
      as=st_alpha
      gsq=4d0*pi*as
      ason2pi=as/(2d0*pi)
c identify semileptonic decays
c Loop over all fs particles, excluding the radiated one
      do j=3,nlegs-1
c W decay products are non-b, non-t fermions, odd are down type (e or d,s)
c even are up type (nu or u,c)
         fl=rflav(j)
c first non-b products
         if(abs(fl).ne.5) then
            if(abs(fl).ne.6.and.abs(fl).lt.20) then
c now only non b t dec. products are allowed
               if(isodd(fl)) then
c e mu tau, d,s
                  if(fl.gt.0) then
                     iem=j
                  else
                     iep=j
                  endif
               else
c nu, u,c
                  if(fl.gt.0) then
                     inu=j
                  elseif (fl .lt. 0) then
                     inub=j
                  else 
                     ig=j
                  endif
               endif
            elseif(abs(fl) .eq. 24) then
               if (fl .eq. 24) then
                   iwp=j
               else
                   iwm=j
               endif
            elseif(abs(fl) .eq. 6) then
               if (fl .eq. 6) then
                   itop=j
               else
                   iatop=j
               endif
            endif
         else
c b
            if(fl.gt.0) then
               ib=j
            else
               ibb=j
            endif
         endif
      enddo
C     setup plabels for MCFM
      plabel(:)='ig'
      if (isaquark(rflav(inu))) plabel(3)='pp'
      if (isaquark(rflav(iem))) plabel(7)='pp'


c Pedantic e nu example
c if(isodd(rflav(j)).and.rflav(j).gt.0.and.rflav(j).lt.20) iem=j
c if(iseven(rflav(j)).and.rflav(j).gt.0.and.rflav(j).lt.20.and.rflav(j).ne.6) inu=j
c         if(rflav(j).eq.11) iem=j
c         if(rflav(j).eq.-12) inub=j
c         if(rflav(j).eq.-11) iep=j
c         if(rflav(j).eq.12) inu=j
c         if(rflav(j).eq.5) ib=j
c         if(rflav(j).eq.-5) ibb=j
c
c iem, inub, ibb come from tb      
c iep, inu, ib come from t
c incoming partons

      ig=13


C     vector for redirection of powheg vector onto mcfm
      myin(:)=-1
      myin(1)=1
      myin(2)=2
      myin(inu)=3
      myin(iep)=4
      myin(ib)=5
      myin(ibb)=6
      myin(iem)=7
      myin(inub)=8
      myin(ig)=9


      do j=1,nlegs
      do ro=1,4
      if (myin(j) .gt. 0) then
      if (j .le. 2) then
c incoming partons
      mcfmp(myin(j),ro)=-p(pwhg(ro),j)
      else
      mcfmp(myin(j),ro)=+p(pwhg(ro),j)
      endif
      endif
      enddo
      enddo

      if (dkflag .eq. 0) then
C---top case
         call qqb_QQbdk_g(mcfmp,msq)
      elseif (dkflag .eq. itop) then
C---top case
         call dk1qqb_QQb_g(mcfmp,msq)
      elseif (dkflag .eq. iatop) then
C---anti-top case
C---top case
         call dk2qqb_QQb_g(mcfmp,msq)
      elseif (dkflag .eq. iwp) then
C---W boson 1
         call dkW1qqb_QQb_g(mcfmp,msq)
      elseif (dkflag .eq. iwm) then
C---W boson 2
         call dkW2qqb_QQb_g(mcfmp,msq)
      endif

 22   continue

      result=msq(rflav(1),rflav(2))/ason2pi

c supply strong corrections to branching ratios, if needed
      result=result*brcorrect(p)

      end
