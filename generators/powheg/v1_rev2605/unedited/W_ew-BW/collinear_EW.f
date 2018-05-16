c*******************************************************************
      subroutine sigcollremn_EW_CC(xjac1,xjac2,z1,z2,coll) !WZGRAD EDIT
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_wzgrad.h' 
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
      integer jb,fl1,fl2
      real*8 CKM(12),fcollz1,fcollz2,splitz1,splitz2,dotp, 
     &       o_fac,sinv(1,2),shat,sig012(2),rescoll_EW_CC(12),
     &       rescoll_QCD(flst_nborn),sig_coll,z1,z2
      real*8 pdfb1(-6:6),pdfb2(-6:6),pdfs1(-6:6),pdfs2(-6:6)
      real*8 coll(12),xjac1,xjac2
      complex*16 prop,prop2
      external dotp
c     begin WZGRAD EDIT------------------      
      character*20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      real*8 powheginput
      external powheginput
c     end WZGRAD EDIT--------------------      
      
      if(abs(powheginput('idvecbos')).eq.24)then

         shat = 2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,2)) 
         sinv(1,2) = shat

         fcollz1=(1d0+z1**2)/(1d0-z1)*dlog((1d0-z1)/z1)
     &          -3d0/2d0/(1d0-z1)+2d0*z1+3d0
         fcollz2=(1d0+z2**2)/(1d0-z2)*dlog((1d0-z2)/z2)
     &          -3d0/2d0/(1d0-z2)+2d0*z2+3d0

         splitz1=alpha0/pi/2d0*((1d0+z1**2)/(1d0-z1)*dlog(sinv(1,2)/
     &           mu_f**2*(1d0-z1)**2/z1*deltac/2d0)+1d0-z1-lfc*fcollz1)
         splitz2=alpha0/pi/2d0*((1d0+z2**2)/(1d0-z2)*dlog(sinv(1,2)/
     &           mu_f**2*(1d0-z2)**2/z2*deltac/2d0)+1d0-z2-lfc*fcollz2)

         if(powheginput('idvecbos').eq.24)then
         sig012(1) = 4d0*(2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,3)))**2
         sig012(2) = 4d0*(2d0*dotp(kn_cmpborn(0,2),kn_cmpborn(0,3)))**2
         elseif(powheginput('idvecbos').eq.-24)then
         sig012(1) = 4d0*(2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,4)))**2
         sig012(2) = 4d0*(2d0*dotp(kn_cmpborn(0,2),kn_cmpborn(0,4)))**2
         endif

         call CKMconvert(CKM) 

         if(wopt.eq.1)then
         prop = 1d0/dcmplx(shat-xmw,ph_WmWw)
         elseif(wopt.eq.2)then
         prop =1d0/dcmplx(shat-xmw,ph_Wwidth*shat/mw)
         endif
         prop2 = prop*dconjg(prop)
         if(rep.eq.1)then
         o_fac = (pi*ph_alphaem/ph_sthw2)**2/3d0/2d0/kn_sborn !has flux
         elseif(rep.eq.2)then
         o_fac = (w2*xmw*gfermi)**2/3d0/2d0/kn_sborn !has flux
         endif 

c        get pdfs at underlying born x values
         call pdfcall(1,kn_xb1,pdfb1)
         call pdfcall(2,kn_xb2,pdfb2)
c        get pdfs at underlying born x/z value values
         call pdfcall(1,kn_xb1/z1,pdfs1)
         call pdfcall(2,kn_xb2/z2,pdfs2)

      do jb=1,12
      fl1=flst_born(1,jb)
      fl2=flst_born(2,jb)
         if(jb.le.6)then
         sig_coll=o_fac*prop2*sig012(2)
         else
         sig_coll=o_fac*prop2*sig012(1)
         endif
      if(z1.gt.1d0-deltas)splitz1=0d0
      if(z2.gt.1d0-deltas)splitz2=0d0
         if(jb.le.6)then
         coll(jb)=CKM(jb)**2*sig_coll*kn_jacborn*
     &   ((1d0/9d0)*pdfs1(fl1)*pdfb2(fl2)*splitz1/z1*xjac1
     &   +(4d0/9d0)*pdfb1(fl1)*pdfs2(fl2)*splitz2/z2*xjac2)
         else
         coll(jb)=CKM(jb)**2*sig_coll*kn_jacborn*
     &   ((4d0/9d0)*pdfs1(fl1)*pdfb2(fl2)*splitz1/z1*xjac1
     &    +(1d0/9d0)*pdfb1(fl1)*pdfs2(fl2)*splitz2/z2*xjac2)
         endif
      enddo
      else
         do jb=1,12
         coll(jb)=0d0
         enddo
      endif
      end
c******************************************************************
