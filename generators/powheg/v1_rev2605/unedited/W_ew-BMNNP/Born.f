      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs)
      integer bflav(nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),bbmunu(0:3,0:3),born,colcf
      integer j,k,mu,nu
c Colour factors for colour-correlated Born amplitudes;
c Rule from 2.98 in FNO2007, leads to B_ij=Cj * B,
c where i#j
      call compborn(p,bflav,born,bbmunu)
      do j=1,nlegs
         if(abs(bflav(j)).le.6) then
            if(bflav(j).eq.0) then
               do mu=0,3
                  do nu=0,3
                     bmunu(mu,nu,j)=bbmunu(mu,nu)
                  enddo
               enddo
            endif
            do k=j+1,nlegs
               if(abs(bflav(k)).le.6) then
                  colcf=cf
               endif
               bornjk(j,k)=born*colcf
               bornjk(k,j)=bornjk(j,k)
            enddo
         endif
      enddo
      end


c     Example
c     q q'-> e+ ve~
      subroutine compborn(p,bflav,born,bmunu)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
c -*- Fortran -*-
c      character *2 flav(-5:5)
      real * 8 charge(-5:5)
c      data (charge(ijkh),ijkh=-5,5) 
c      data (flav(ijkh),ijkh=-5,5) 
c      data flav
c     #     /'b~','c~','s~','u~','d~','g','d','u','s','c','b'/
      data charge
     #     / 0.33333333333333333333d0, !   1d0/3
     #      -0.66666666666666666667d0, !  -2d0/3
     #       0.33333333333333333333d0, !   1d0/3 
     #      -0.66666666666666666667d0, !   -2d0/3
     #       0.33333333333333333333d0, !   1d0/3 
     #       0d0,                      !   0d0   
     #      -0.33333333333333333333d0, !   -1d0/3
     #       0.66666666666666666667d0, !   2d0/3   
     #      -0.33333333333333333333d0, !   -1d0/3
     #       0.66666666666666666667d0, !   2d0/3 
     #      -0.33333333333333333333d0/ !   -1d0/3
c      include 'QuarkFlavs.h'
      include 'PhysPars.h'
      integer nleg
      parameter (nleg=nlegborn)
      real * 8 p(0:3,nleg)
      integer bflav(nleg)
      real * 8 amp2,born,bmunu(0:3,0:3)
      integer ferm_type(nleg)
      integer i,j
      real * 8 ferm_charge(nleg)
c     vector boson id and decay
      integer idvecbos,vdecaymode
      common/cvecbos/idvecbos,vdecaymode   

      if (abs(bflav(3)).le.6.or.abs(bflav(4)).le.6) then
         write(*,*) 'born_ampsq: ERROR in flavor assignement'
         stop
      endif
 
c     i is the flavour index of first incoming parton
c     j is the flavour index of second incoming parton
c     with the convention:
c     
c      -6  -5  -4  -3  -2  -1  0  1  2  3  4  5  6                    
c      t~  b~  c~  s~  u~  d~  g  d  u  s  c  b  t                    
      
      i = bflav(1)
      j = bflav(2)
      ferm_charge(1) = charge(i)
      ferm_charge(2) = charge(j)
      ferm_type(1) = i/abs(i)
      ferm_type(2) = j/abs(j)


c     antilepton-neutrino from W decay
      ferm_type(3) = bflav(3)/abs(bflav(3))
      ferm_charge(3) = ferm_type(3)*(-1d0)
      ferm_type(4) = -ferm_type(3)
      ferm_charge(4) = 0d0

      
      if(idvecbos.eq.24) then
         call q_aqp_to_al_vl(p,ferm_type,ferm_charge,
     $        amp2)
      elseif(idvecbos.eq.-24) then
         call q_aqp_to_l_avl(p,ferm_type,ferm_charge,
     $        amp2)
      else
         write(*,*) 'ERROR: this subroutine deals only with W+ or W- '
         call exit(1)
      endif
      
      
      if(mod(abs(i),2).eq.0) then
         born=amp2*ph_CKM(abs(i)/2,(abs(j)+1)/2)**2
      elseif(mod(abs(i),2).eq.1) then   
         born=amp2*ph_CKM(abs(j)/2,(abs(i)+1)/2)**2
      endif

      do i=0,3
         do j=0,3
            bmunu(i,j)=0d0
         enddo
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine q_aqp_to_l_avl(pphy,fermion_type,fermion_charge,
     #     amp2)
   
      implicit none
c the 5 4-momentum vectors
c p(i,1) is the i-th component of vector p1...   
      integer nleg
      parameter (nleg=4)
      integer fermion_type(nleg),i
      real * 8 fermion_charge(nleg)
      real * 8 pphy(0:3,nleg)
      real * 8 amp2
      real * 8 ferm_charge(nleg)
      integer ferm_type(nleg)

       if ((fermion_type(3).ne.1).and.(fermion_type(4).ne.-1)) then
         write(*,*) 'ERROR: this subroutine deals only with W- decay'
         stop
      endif

      do i=1,nleg
      
         ferm_charge(i) = -fermion_charge(i)
         ferm_type(i) = -fermion_type(i)
      enddo
            
      
      call q_aqp_to_al_vl(pphy,ferm_type,ferm_charge,
     #     amp2)

      end

*
** function to calculate bornme of dy process
** mu = md = 0
**
** (pu) u  \            / e+ (pl)
**          \          /
**           \________/
**           /   W+   \ 
**          /          \
** (pd) d~ /            \ nu_e (pn)
*
      subroutine q_aqp_to_al_vl(pphy,fermion_type,fermion_charge,
     $     amp2)
      implicit none
c the nleg 4-momentum vectors
c p(i,1) is the i-th component of vector p1...   
      integer nleg
      parameter (nleg=4)
      integer fermion_type(nleg)
      real * 8 fermion_charge(nleg)
      real * 8 pphy(0:3,nleg)
      real * 8 amp2
      include 'pwhg_math.h'
      include 'mathx.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      real*8 dotp
      external dotp
*
      real*8 pu(0:3),pd(0:3),pn(0:3),pl(0:3),ptmp(0:3)
      complex*16 den
      real*8 pupl,pdpn,s
      real*8 tmp
      real*8 p(0:3,nleg)
      real*8 ferm_charge(nleg)
      integer ferm_type(nleg)
      integer i,nu
*

      if ((fermion_type(3).ne.-1).and.(fermion_type(4).ne.1)) then
         write(*,*) 'ERROR: this subroutine deals only with W+ decay'
         stop
      endif
     

c  copy to local variables
      do i=1,nleg
         do nu=0,3
            p(nu,i) = pphy(nu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo
      
c     exchance particle 1 and 2
      if (ferm_type(1).eq.-1) then
         if (ferm_type(2).eq.1) then
            call exchange_momenta(p(0,1),p(0,2))
            tmp = ferm_charge(1)
            ferm_charge(1)=ferm_charge(2)
            ferm_charge(2)=tmp
            tmp = ferm_type(1)
            ferm_type(1)=ferm_type(2)
            ferm_type(2)=tmp
         else
            write(*,*) 'Error in the type of the quark 1-2'
            stop
         endif
      endif

      if (ferm_charge(1)*ferm_type(1).lt.0d0) then
          do nu=0,3
              pu(nu) = p(nu,2)
              pd(nu) = p(nu,1)
          enddo
      else
          do nu=0,3
              pu(nu) = p(nu,1)
              pd(nu) = p(nu,2)
          enddo
      endif

      do nu=0,3
          pl(nu) = p(nu,3)
          pn(nu) = p(nu,4)
          ptmp(nu) = pu(nu) + pd(nu)
      enddo

      pupl = dotp(pu,pl)
      pdpn = dotp(pd,pn)
      s = dotp(ptmp,ptmp)

      den  = 1.d0/( s - ph_Wmass2 + ii*ph_WmWw )

      amp2 = 16  *pupl*pdpn
     +       *den*conjg(den)
     +       *g2*conjg(g2)/4d0 
     +        /4d0/nc

      end
