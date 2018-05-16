      subroutine qqb_w2jet_g_pwhg_nc(p,rflav,jperm,res)
************************************************************************
*     Author: R.K.Ellis                                                *
*     June, 2012.                                                      *
c---  Matrix element squared averaged over initial colors and spins    *
c     q(-p1)+qbar(-p2) -->  g*  + W^+ + g(p7)                          *
c                           |     |                                    *
c                           |      --> nu(p3)+e^+(p4)                  *
c                           |   or --> e^-(p3)+nu~(p4)                 *
c                           |                                          *
c                            ---> f(p5)+f(p6)                          *
c     This routine only works:-                                        *
c     1) for the (idiosyncratic) ordering in rflav dictated by mcfm    *
c     2) for diagonal CKM                                              *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'flags.f'
      include 'nwz.f'
      include 'ckmallowed.f'
      integer nq,k,j,i5,i6,i7,rflav(7),up,dn,i,nu
      integer maxperm,jperm
      parameter (maxperm=6)      
      double precision p(mxpart,4),res
      double precision, save ::
     &mmsq_gg(maxperm),mmsq_qqb(maxperm),mmsq_qbq(maxperm),
     &mmsq_qg(maxperm),mmsq_gq(maxperm),mmsq_gqb(maxperm),
     &mmsq_qbg(maxperm),
     &QQ_ud_dd(maxperm),QQ_us_ds(maxperm),QQ_uu_du(maxperm),
     &QQ_du_dd(maxperm),QQ_su_sd(maxperm),QQ_uu_ud(maxperm),
     &QbQb_ds_us(maxperm),QbQb_du_uu(maxperm),QbQb_dd_ud(maxperm),
     &QbQb_sd_su(maxperm),QbQb_ud_uu(maxperm),QbQb_dd_du(maxperm),
     &RRb_dd_du(maxperm),RRb_ss_du(maxperm),RRb_uu_du(maxperm),
     &RRb_ud_dd(maxperm),RRb_ud_ss(maxperm),RRb_ud_uu(maxperm),
     &QQb_ud_dd(maxperm),QQb_us_ds(maxperm),QQb_ud_uu(maxperm),
     &QQb_du_dd(maxperm),QQb_su_sd(maxperm),QQb_uu_ud(maxperm),
     &RbR_dd_du(maxperm),RbR_ss_du(maxperm),RbR_uu_du(maxperm),
     &RbR_ud_dd(maxperm),RbR_ud_ss(maxperm),RbR_ud_uu(maxperm),
     &QbQ_ud_dd(maxperm),QbQ_us_ds(maxperm),QbQ_ud_uu(maxperm),
     &QbQ_du_dd(maxperm),QbQ_su_sd(maxperm),QbQ_uu_ud(maxperm),
     &QG_d_ddu(maxperm),QG_d_dsc(maxperm),QG_u_ddd(maxperm),
     &QG_u_dcc(maxperm),QG_u_duu(maxperm),QG_u_udu(maxperm),
     &GQ_d_ddu(maxperm),GQ_d_dsc(maxperm),GQ_u_ddd(maxperm),
     &GQ_u_dcc(maxperm),GQ_u_duu(maxperm),GQ_u_udu(maxperm),
     &QbG_d_ddu(maxperm),QbG_u_ucs(maxperm),QbG_u_uud(maxperm),
     &QbG_d_uuu(maxperm),QbG_d_ucc(maxperm),QbG_d_udd(maxperm),
     &GQb_d_ddu(maxperm),GQb_u_ucs(maxperm),GQb_u_uud(maxperm),
     &GQb_d_udd(maxperm),GQb_d_ucc(maxperm),GQb_d_uuu(maxperm),
     & propsq
      logical first,qa,aq,qq,aa,gq,ga,qg,ag,gg,diagonal
      integer jj(-nf:nf),kk(-nf:nf)
      data jj/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      data kk/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      data first/.true./
      save first,up,dn
C     variables needed to avoid recalculating same stuff
      integer, parameter :: nchn = 9 
      logical ::  recalc_g(nchn), recalc_q(nchn)
      double precision, save :: opin_g(mxpart,4,nchn,maxperm)
      double precision, save :: opin_q(mxpart,4,nchn,maxperm)
      integer, parameter :: iqa = 1 
      integer, parameter :: iaq = 2 
      integer, parameter :: iqg = 3 
      integer, parameter :: igq = 4 
      integer, parameter :: iag = 5 
      integer, parameter :: iga = 6 
      integer, parameter :: iqq = 7 
      integer, parameter :: iaa = 8 
      integer, parameter :: igg = 9 
      integer :: ichn 
      if (first) then
         first=.false.
         if (abs(nwz) == 1) then
C----setup allowed values of ckm (diagonal)
C----This routine only works with a diagonal CKM
            diagonal=.true.
            call setupckmallowed(nwz,diagonal)
            else 
            write(6,*) 'Unacceptable value of nwz, nwz=',nwz
            stop
         endif 
         up=(3+nwz)/2
         dn=(3-nwz)/2
      endif

C------first determine if it is a two- or four-quark process 
      nq = 0 
      do j=1,7
         if (j < 3 .or. j > 4) then 
            if (abs(rflav(j)) > 0) then 
               nq = nq+1
            endif
         endif
      enddo
      if (nq == 2) then 
         Gflag = .true. 
         Qflag = .false. 
      elseif (nq == 4) then 
         Qflag = .true. 
         Gflag = .false. 
      else 
         write(*,*) 'nq',nq 
         write(*,*) 'rflav',rflav
         stop 'nq out of range' 
      endif
      j=rflav(1)
      k=rflav(2)
      i5=rflav(5)
      i6=rflav(6)
      i7=rflav(7)

      qa=(j > 0) .and. (k < 0)
      aq=(j < 0) .and. (k > 0) 
      qq=(j > 0) .and. (k > 0) 
      aa=(j < 0) .and. (k < 0) 
      qg=(j > 0) .and. (k == 0) 
      ag=(j < 0) .and. (k == 0) 
      gq=(j == 0) .and. (k > 0) 
      ga=(j == 0) .and. (k < 0) 
      gg=(j == 0) .and. (k == 0) .and. ckmallowed(-i5,-i6)

      if (qa) then 
         ichn = iqa
      elseif (aq) then 
         ichn = iaq
      elseif (qg) then 
         ichn = iqg
      elseif (gq) then 
         ichn = igq
      elseif (ag) then 
         ichn = iag
      elseif (ga) then 
         ichn = iga
      elseif (qq) then 
         ichn = iqq
      elseif (aa) then 
         ichn = iaa
      elseif (gg) then 
         ichn = igg
      else
         stop 'qqb_w2jet_g_pwhg: channel not recognized' 
      endif
!     decide if recalc
      recalc_q = .false. 
      recalc_g = .false. 

!      call spinoru(7,p,za,zb)

      s(3,4)=two*(p(3,4)*p(4,4)-p(3,1)*p(4,1)
     &     -p(3,2)*p(4,2)-p(3,3)*p(4,3))

      propsq=s(3,4)**2/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

      if (Gflag) then
         do nu=1,4
            do i=1,7
               if(opin_g(i,nu,ichn,jperm).ne.p(i,nu)) then
                  recalc_g(ichn) = .true.
                  goto 10
               endif
            enddo
         enddo
      else 
         do nu=1,4
            do i=1,7
               if(opin_q(i,nu,ichn,jperm).ne.p(i,nu)) then
                  recalc_q(ichn) = .true.
                  goto 10
               endif
            enddo
         enddo
         
      endif

 10    continue 

       if (recalc_g(ichn)) then
          opin_g(:,:,ichn,jperm) = p
       elseif (recalc_q(ichn)) then
          opin_q(:,:,ichn,jperm) = p
       endif  
      

       if (Gflag) then
************************************************************************
*     Calculate contributions from the QQGG matrix elements            *
************************************************************************
          
          if (recalc_g(ichn)) then
             call spinoru(7,p,za,zb)
             if (recalc_g(igg)) then 
                call xwqqggg(5,1,2,7,6,3,4,mmsq_gg(jperm))
             elseif (recalc_g(iqa)) then 
                call xwqqggg(1,5,6,7,2,3,4,mmsq_qqb(jperm))
             elseif (recalc_g(iaq)) then 
                call xwqqggg(2,5,6,7,1,3,4,mmsq_qbq(jperm))
             elseif (recalc_g(iqg)) then 
                call xwqqggg(1,2,6,7,5,3,4,mmsq_qg(jperm))
             elseif (recalc_g(igq)) then 
                call xwqqggg(2,1,6,7,5,3,4,mmsq_gq(jperm))
             elseif (recalc_g(iga)) then 
                call xwqqggg(5,1,6,7,2,3,4,mmsq_gqb(jperm))
             elseif (recalc_g(iag)) then 
                call xwqqggg(5,2,6,7,1,3,4,mmsq_qbg(jperm))
             endif
          endif
       endif

      if (Qflag) then
************************************************************************
*     Calculate contributions from the QQBQQB matrix elements          *
************************************************************************

         if (recalc_q(ichn)) then
            call spinoru(7,p,za,zb)
            if (recalc_q(iqq)) then
c---  basic Q-Q amplitudes 
            call addhel(1,5,2,6,7,3,4,QQ_ud_dd(jperm),QQ_us_ds(jperm),
     &              QQ_uu_du(jperm))
            call addhel(2,6,1,5,7,3,4,QQ_du_dd(jperm),QQ_su_sd(jperm),
     &           QQ_uu_ud(jperm))
            elseif (recalc_q(iaa)) then
c---  basic Qb-Qb amplitudes 
            call addhel(5,1,6,2,7,3,4,QbQb_dd_ud(jperm),
     &              QbQb_ds_us(jperm),QbQb_du_uu(jperm))
            call addhel(6,2,5,1,7,3,4,QbQb_dd_du(jperm),
     &           QbQb_sd_su(jperm),QbQb_ud_uu(jperm))
            elseif (recalc_q(iqa)) then
c---  basic Q-Qb amplitudes
c---  annihilation 
            call addhel(6,5,1,2,7,3,4,RRb_dd_du(jperm),
     &              RRb_ss_du(jperm),RRb_uu_du(jperm))
            call addhel(1,2,6,5,7,3,4,RRb_ud_dd(jperm),
     &           RRb_ud_ss(jperm),RRb_ud_uu(jperm))
c---  scattering 
            call addhel(1,5,6,2,7,3,4,QQb_ud_dd(jperm),
     &           QQb_us_ds(jperm),QQb_ud_uu(jperm))
            call addhel(6,2,1,5,7,3,4,QQb_du_dd(jperm),
     &           QQb_su_sd(jperm),QQb_uu_ud(jperm))
            elseif (recalc_q(iaq)) then
c---  basic Qb-Q amplitudes
c---  annihilation 
            call addhel(6,5,2,1,7,3,4,RbR_dd_du(jperm),
     &              RbR_ss_du(jperm),RbR_uu_du(jperm))
            call addhel(2,1,6,5,7,3,4,RbR_ud_dd(jperm),
     &           RbR_ud_ss(jperm),RbR_ud_uu(jperm))
c---  scattering 
            call addhel(2,5,6,1,7,3,4,QbQ_ud_dd(jperm),
     &           QbQ_us_ds(jperm),QbQ_ud_uu(jperm))
            call addhel(6,1,2,5,7,3,4,QbQ_du_dd(jperm),
     &           QbQ_su_sd(jperm),QbQ_uu_ud(jperm))
            elseif (recalc_q(iqg)) then
c---  basic Q-G amplitudes
            call addhel(7,6,1,5,2,3,4,QG_d_ddu(jperm),
     &              QG_d_dsc(jperm),QG_u_udu(jperm))
            call addhel(1,5,7,6,2,3,4,QG_u_ddd(jperm),
     &           QG_u_dcc(jperm),QG_u_duu(jperm))
            elseif (recalc_q(iag)) then
c---  basic QB-G amplitudes
            call addhel(6,7,5,1,2,3,4,QbG_d_ddu(jperm),
     &              QbG_u_ucs(jperm),QbG_u_uud(jperm))
            call addhel(5,1,6,7,2,3,4,QbG_d_udd(jperm),
     &           QbG_d_ucc(jperm),QbG_d_uuu(jperm))
            elseif (recalc_q(igq)) then
c---  basic G-Q amplitudes
            call addhel(7,5,2,6,1,3,4,GQ_d_ddu(jperm),
     &              GQ_d_dsc(jperm),GQ_u_udu(jperm))
            call addhel(2,6,7,5,1,3,4,GQ_u_ddd(jperm),
     &           GQ_u_dcc(jperm),GQ_u_duu(jperm))
            elseif (recalc_q(iga)) then
c---  basic G-QB amplitudes
            call addhel(5,7,6,2,1,3,4,GQb_d_udd(jperm),
     &              GQb_u_ucs(jperm),GQb_u_uud(jperm))
            call addhel(6,2,5,7,1,3,4,GQb_d_ddu(jperm),
     &           GQb_d_ucc(jperm),GQb_d_uuu(jperm))
            endif
         endif
      endif

      if (Gflag) then
************************************************************************
*     Sum the contributions from the QQGG matrix elements              *
************************************************************************

c---  note the identical particle factor of 1/6 for the
c---  q-qb initial states, due to 3 gluons in the final state     
         if     (gg) then
            res=propsq*mmsq_gg(jperm)*(gwsq**2/4d0/esq**2)
         elseif (qa) then
            res=propsq*mmsq_qqb(jperm)
     .           *(aveqq/avegg)*(gwsq**2/4d0/esq**2)/6d0
         elseif (aq) then
            res=propsq*mmsq_qbq(jperm)
     .           *(aveqq/avegg)*(gwsq**2/4d0/esq**2)/6d0
         elseif (qg) then
            res=half*propsq*mmsq_qg(jperm)
     .           *(aveqg/avegg)*(gwsq**2/4d0/esq**2) 
         elseif (ag) then
            res=half*propsq*mmsq_qbg(jperm)
     .           *(aveqg/avegg)*(gwsq**2/4d0/esq**2) 
         elseif (gq) then
            res=half*propsq*mmsq_gq(jperm)
     .           *(aveqg/avegg)*(gwsq**2/4d0/esq**2) 
         elseif (ga) then
            res=half*propsq*mmsq_gqb(jperm)
     .           *(aveqg/avegg)*(gwsq**2/4d0/esq**2) 
         endif
      endif
      

      if (Qflag) then

************************************************************************
*     Sum the contributions from the QQBQQB matrix elements            *
************************************************************************

         if (qq) then
c---  Q Q --> Q Q
            if ((jj(j)==up).and.(kk(k)==dn).and.ckmallowed(j,-i5))then
               if (i5==i6) then
                  res=QQ_ud_dd(jperm)*0.5d0
               elseif (i5/=i6) then
                  res=QQ_us_ds(jperm)
               endif
           elseif((jj(j)==dn).and.(kk(k)==up).and.ckmallowed(k,-i6))then
               if (i5==i6) then
                  res=QQ_du_dd(jperm)*0.5d0
               elseif (i5/=i6) then
                  res=QQ_su_sd(jperm)
               endif
            elseif ((jj(j)==up) .and. (kk(k)==up)) then
               if (j == k) then
                  res=QQ_uu_du(jperm)
               else
                  if (ckmallowed(j,-i5)) res=QQ_us_ds(jperm)
                  if (ckmallowed(k,-i6)) res=QQ_su_sd(jperm)
               endif
            endif

         elseif (aa) then
c---  Qb Qb --> Qb Qb
            if((jj(j)==-up).and.(kk(k)==-dn).and.ckmallowed(k,-i6)) then
               if (i5==i6) then
                  res=QbQb_ud_uu(jperm)*0.5d0
               else
                  res=QbQb_sd_su(jperm)
               endif
         elseif((jj(j)==-dn).and.(kk(k)==-up).and.ckmallowed(j,-i5))then
               if (i5==i6) then
                  res=QbQb_du_uu(jperm)*0.5d0
               else
                  res=QbQb_ds_us(jperm)
               endif
            elseif ((jj(j) == -dn) .and. (kk(k) == -dn)) then
               if (j == k) then
                  res=QbQb_dd_ud(jperm)
               else
                  if (ckmallowed(j,-i5)) res=QbQb_ds_us(jperm)
                  if (ckmallowed(k,-i6)) res=QbQb_sd_su(jperm)
               endif
            endif

         elseif (qa) then
c---  Q Qb --> Q Qb
            if ((jj(j)==up) .and. (kk(k)==-dn)) then
               if (ckmallowed(j,k)) then
                  if (ckmallowed(j,-i5)) then
                     res=RRb_ud_dd(jperm)
                  elseif (ckmallowed(k,-i6)) then
                     res=RRb_ud_uu(jperm)
                  else
                     res=RRb_ud_ss(jperm)
                  endif
               else
                  if (ckmallowed(j,-i5)) res=QQb_us_ds(jperm)
                  if (ckmallowed(k,-i6)) res=QQb_su_sd(jperm)
               endif
            elseif ((jj(j)==up) .and. (kk(k)==-up)) then
               if (j==-k) then
                  if (ckmallowed(j,-i5).and.ckmallowed(-i5,-i6)) then
                     res=RRb_uu_du(jperm)
                  elseif (ckmallowed(-i5,-i6)) then
                     res=RRb_ss_du(jperm)
                  endif
               else
                  res=QQb_us_ds(jperm)
               endif
            elseif ((jj(j)==dn) .and. (kk(k)==-dn)) then
               if (j==-k) then
                  if (ckmallowed(k,-i6).and.ckmallowed(-i5,-i6)) then
                     res=RRb_dd_du(jperm)
                  elseif (ckmallowed(-i5,-i6)) then
                     res=RRb_ss_du(jperm)
                  endif
               else
                  res=QQb_su_sd(jperm)
               endif
            endif
         elseif (aq) then
c---  Qb Q --> Q Qb
            if ((jj(j)==-dn) .and. (kk(k)==up)) then
               if (ckmallowed(j,k)) then
                  if (ckmallowed(k,-i5)) then
                     res=RbR_ud_dd(jperm)
                  elseif (ckmallowed(j,-i6)) then
                     res=RbR_ud_uu(jperm)
                  else
                     res=RbR_ud_ss(jperm)
                  endif
               else
                  if (ckmallowed(k,-i5)) res=QbQ_us_ds(jperm)
                  if (ckmallowed(j,-i6)) res=QbQ_su_sd(jperm)
               endif
            elseif ((jj(j)==-up) .and. (kk(k)==up)) then
               if (j==-k) then
                  if (ckmallowed(k,-i5).and.ckmallowed(-i5,-i6)) then
                     res=RbR_uu_du(jperm)
                  elseif (ckmallowed(-i5,-i6)) then
                     res=RbR_ss_du(jperm)
                  endif
               else
                  res=QbQ_us_ds(jperm)
               endif
            elseif ((jj(j)==-dn) .and. (kk(k)==dn)) then
               if (j==-k) then
                  if (ckmallowed(j,-i6).and.ckmallowed(-i5,-i6)) then
                     res=RbR_dd_du(jperm)
                  elseif (ckmallowed(-i5,-i6)) then
                     res=RbR_ss_du(jperm)
                  endif
               else
                  res=QbQ_su_sd(jperm)
               endif
            endif
         elseif (qg) then
c---  Q G --> Q q qb
            if( jj(j)==up) then
               if ((i6==-i7) .and. ckmallowed(j,-i5)) then
                  if  (i5==i6) then
                     res=(aveqg/aveqq)*QG_u_ddd(jperm)*0.5d0
                  else
                     res=(aveqg/aveqq)*QG_u_dcc(jperm)
                  endif

               elseif (ckmallowed(-i6,-i7)) then 
                  if (ckmallowed(j,-i6)) then
                     res=(aveqg/aveqq)*QG_u_udu(jperm)
                  else
                     res=(aveqg/aveqq)*QG_d_dsc(jperm)
                  endif
               endif
            elseif( jj(j) == dn) then
               if (ckmallowed(-i6,-i7)) then
                  if (i5 == i6) then
                     res=(aveqg/aveqq)*QG_d_ddu(jperm)*0.5d0
                  else
                     res=(aveqg/aveqq)*QG_d_dsc(jperm)
                  endif
               endif
            endif

         elseif (ag) then
c---  QB G --> QB qb q
            if( jj(j)==-dn) then
               if ((i6==-i7) .and. ckmallowed(j,-i5)) then
                  if  (i5==i6) then
                     res=(aveqg/aveqq)*QbG_d_uuu(jperm)*0.5d0
                  else
                     res=(aveqg/aveqq)*QbG_d_ucc(jperm)
                  endif

               elseif (ckmallowed(-i6,-i7)) then 
                  if (ckmallowed(j,-i6)) then
                     res=(aveqg/aveqq)*QbG_d_ddu(jperm)
                  else
                     res=(aveqg/aveqq)*QbG_u_ucs(jperm)
                  endif
               endif
            elseif( jj(j)==-up) then
               if (ckmallowed(-i6,-i7)) then
                  if (i5 == i6) then
                     res=(aveqg/aveqq)*QbG_u_uud(jperm)*0.5d0
                  else
                     res=(aveqg/aveqq)*QbG_u_ucs(jperm)
                  endif
               endif
            endif
         elseif (gq) then
c---  G Q --> Q q qb
            if( kk(k)==up) then
               if ((i5==-i7) .and. ckmallowed(k,-i6)) then
                  if  (i5==i6) then
                     res=(aveqg/aveqq)*GQ_u_ddd(jperm)*0.5d0
                  else
                     res=(aveqg/aveqq)*GQ_u_dcc(jperm)
                  endif
               elseif (ckmallowed(-i5,-i7)) then 
                  if (ckmallowed(k,-i5)) then
                     res=(aveqg/aveqq)*GQ_u_udu(jperm)
                  else
                     res=(aveqg/aveqq)*GQ_d_dsc(jperm)
                  endif
               endif
            elseif( kk(k)==dn) then
               if (ckmallowed(-i5,-i7)) then
                  if (i5 == i6) then
                     res=aveqg/aveqq*GQ_d_ddu(jperm)*0.5d0
                  else
                     res=aveqg/aveqq*GQ_d_dsc(jperm)
                  endif
               endif
            endif

         elseif (ga) then
c---  G QB --> QB qb q
            if( kk(k)==-dn) then
               if ((i5==-i7) .and. ckmallowed(k,-i6)) then
                  if (i5==i6) then
                     res=(aveqg/aveqq)*GQb_d_uuu(jperm)*0.5d0
                  else
                     res=(aveqg/aveqq)*GQb_d_ucc(jperm)
                  endif
               elseif (ckmallowed(-i5,-i7)) then 
                  if (ckmallowed(k,-i5)) then
                     res=(aveqg/aveqq)*GQb_d_udd(jperm)
                  else
                     res=(aveqg/aveqq)*GQb_u_ucs(jperm)
                  endif
               endif
            elseif( kk(k)==-up) then
               if (ckmallowed(-i5,-i7)) then
                  if (i5 == i6) then
                     res=(aveqg/aveqq)*GQb_u_uud(jperm)*0.5d0
                  else
                     res=(aveqg/aveqq)*GQb_u_ucs(jperm)
                  endif
               endif
            endif
         endif

      endif

      return
      end


