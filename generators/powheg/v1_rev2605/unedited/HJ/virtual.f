      subroutine setvirtual(p,vflav,virtual)
c Virtual needs to be provided by the user and put here
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 p(0:3,nlegborn)
      integer vflav(nlegborn)
      double precision res(0:3),c1,c2,c4
      real * 8 bornjk(nlegborn,nlegborn),bmunu(0:3,0:3,nlegborn),born
      real * 8 virtual,virtual_DR
      real * 8 s,t,u,virtgg,virtqa,virtaq,virtqg,virtgq
      real * 8 dotp
      external dotp

      call i2MCFM_2_POWHEG(p,vflav,res) 

      virtual_DR=res(0)/(st_alpha/(2*pi))

c from dimensional reduction to dimensional regularization
      call from_DR_to_CDR(p,vflav,virtual_DR,virtual)


cc The following code is used to check the double pole
cc in the virtual amplitude. Must set in
cc set_interface_MCFM.f:      ret_poles=.true.
c
c      call setborn(p,vflav,born,bornjk,bmunu)
c      if(vflav(1).eq.0) then
c         c1=3
c      else
c         c1=4d0/3
c      endif
c      if(vflav(2).eq.0) then
c         c2=3
c      else
c         c2=4d0/3
c      endif
c      if(vflav(4).eq.0) then
c         c4=3
c      else
c         c4=4d0/3
c      endif
c      write(*,*) ' Pole check: ',vflav,
c     1 (res(2)/(st_alpha/(2*pi)))/(born*(-(c1+c2+c4))),virtual_DR

      end

c     virtual_DR : finite part of the virtual in Dimensional Reduction, stripped off by the factor as/(2pi)
c     returns virtual_CDR, i.e. the virtual in conventional dimensional regularization.
      subroutine from_DR_to_CDR(p,vflav,virtual_DR,virtual_CDR)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      real * 8 p(0:3,nlegborn)
      integer vflav(nlegborn)
      real * 8 virtual_DR,virtual_CDR
      real * 8 bornjk(nlegborn,nlegborn),bmunu(0:3,0:3,nlegborn),born
      integer i
      real * 8 gammag,gammaq
      parameter (gammag=(nc*1d0)/6d0,gammaq=(nc**2-1d0)/(4d0*nc)) 

      call setborn(p,vflav,born,bornjk,bmunu)
      virtual_CDR=virtual_DR
c     initial-state partons
      do i=1,2
         if (vflav(i).eq.0) then
            virtual_CDR=virtual_CDR-born*gammag
         else
            virtual_CDR=virtual_CDR-born*gammaq
         endif
      enddo
      do i=flst_lightpart,nlegborn
         if (vflav(i).eq.0) then
            virtual_CDR=virtual_CDR-born*gammag
         else
            virtual_CDR=virtual_CDR-born*gammaq
         endif         
      enddo
      end
