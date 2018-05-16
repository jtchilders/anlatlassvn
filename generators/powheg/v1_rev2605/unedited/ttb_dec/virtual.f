c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c     where M_B is the Born amplitude and 
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
      subroutine setvirtual(p,vflav,vres)
      implicit none
      include 'nlegborn.h'
c      include 'pwhg_math.h'
c      include 'pwhg_st.h'
c      include 'pwhg_kn.h'
      include 'constants.f'
      include 'plabel.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'process.f'
      include 'mcfmtopwhg.f'
      integer nlegs
      parameter (nlegs=nlegborn)
c For t tbar with decay, nlegborn=12
      real * 8 p(0:3,nlegs),vres
      integer vflav(nlegs),ro
c mxpart is in constants.f
      real * 8 mcfmp(mxpart,4)
      double precision msqp(-nf:nf,-nf:nf),
     & msqsl(-nf:nf,-nf:nf),msqh(-nf:nf,-nf:nf)
      integer iem,iep,inu,inub,ib,ibb,fl,j,mu,myin(nlegs)
      logical isaquark, isodd
      real * 8 psav(0:3,nlegs)
      save psav, msqp,msqsl,msqh
      logical ini
      integer jnlowhich
      common/cjnlowhich/jnlowhich
      data ini/.true./
      save ini
      real * 8 powheginput,brcorrect
      external powheginput,brcorrect
      if(ini) then
         scheme='dred'
         epinv=0
         epinv2=0
         ini=.false.
      endif
      do mu=0,3
         do j=1,nlegs
            if(p(mu,j).ne.psav(mu,j)) goto 11
         enddo
      enddo
      goto 22

 11   continue
      psav=p
c identify semileptonic decays
      do j=3,nlegs
c W decay products are non-b, non-t fermions, odd are down type (e or d,s)
c even are up type (nu or u,c)
         fl=vflav(j)
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
                  else
                     inub=j
                  endif
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
      if (isaquark(vflav(inu))) then
      plabel(3)='pp'
      case = 'tt_bbh'
      endif
      if (isaquark(vflav(iem))) then
      plabel(7)='pp'
      case = 'tt_bbh'
      endif

c Pedantic e nu example
c if(isodd(vflav(j)).and.vflav(j).gt.0.and.vflav(j).lt.20) iem=j
c if(iseven(vflav(j)).and.vflav(j).gt.0.and.vflav(j).lt.20.and.vflav(j).ne.6) inu=j
c         if(vflav(j).eq.11) iem=j
c         if(vflav(j).eq.-12) inub=j
c         if(vflav(j).eq.-11) iep=j
c         if(vflav(j).eq.12) inu=j
c         if(vflav(j).eq.5) ib=j
c         if(vflav (j).eq.-5) ibb=j
c
c iem, inub, ibb come from tb      
c iep, inu, ib come from t


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

C----Translate powheg momenta into MCFM notation
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


C----initialize matrix elements
      msqp(:,:)=0d0
      msqsl(:,:)=0d0
      msqh(:,:)=0d0

C-----Radiation in production only
      if(jnlowhich.eq.0.or.jnlowhich.eq.1)
     1     call qqb_QQbdkBSY_v(mcfmp,msqp) 

C-----Radiation in topdecay only
      if(jnlowhich.eq.0.or.jnlowhich.eq.2)
     1  call dkqqb_QQb_v(mcfmp,msqsl)

C-----Radiation in topdecay only
      if(jnlowhich.eq.4) then
c I am not sure if this is correct; both t and anti-t virtual
c corrections are included, but they should not affect the total
        call dkqqb_QQb_v(mcfmp,msqsl)
      endif

      if(jnlowhich.eq.0.or.jnlowhich.eq.3 .and.
     &  ((plabel(3) .eq. 'pp') .or. (plabel(7) .eq. 'pp')) ) then
c Check consistency of MCFM result
c         ll=log((p(0,6)**2-p(1,6)**2-p(2,6)**2-p(3,6)**2)/st_muren2)
c         call  qqb_QQbdk(mcfmp,msqh)
c         msqh=msqh*(-ll**2+3*ll+pi**2-7-3d0/2)*ason2pi*cf
c         write(*,*) msqh(0,0), 'aaaaa'
         call dkWqqb_QQb_v(mcfmp,msqh)
c         write(*,*) msqh(0,0)
      endif


 22   continue

      vres = msqp(vflav(1),vflav(2))
     &      +msqsl(vflav(1),vflav(2))
     &       +msqh(vflav(1),vflav(2))
      vres = vres/ason2pi

c Supply strong correction to branching ratios, if needed
      vres = vres * brcorrect(p)

      end

      function isaquark(i)
      implicit none
      logical isaquark
      integer i
      isaquark=(abs(i).gt.0.and.abs(i).le.4)
      end
