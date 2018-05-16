      subroutine setvirtual(p,vflav,virtual)
c Virtual needs to be provided by the user and put here
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 p(0:3,nlegborn),virtual
      integer vflav(nlegborn)
      real * 8 plocal(0:3,nlegborn)
      integer vlocal(nlegborn)
      real * 8 born,bornjk(nlegborn,nlegborn),bmunu(0:3,0:3,nlegborn)
      real * 8 virtual_DR,res(0:3),c1,c2,c4,c5,tmp
      logical polecheck,ini
      data ini/.true./
      save polecheck,ini
      real * 8 powheginput
      external powheginput
      if(ini) then
         if(powheginput("#polecheck").gt.0) then
            polecheck=.true. 
         else
            polecheck=.false.
         endif
         ini=.false.
      endif
      plocal=p
      vlocal=vflav
      call i2MCFM_2_POWHEG(plocal,vlocal,res)
      if(res(0).eq.0.and.vlocal(4).ne.vlocal(5)) then
         plocal(:,4)=p(:,5)
         plocal(:,5)=p(:,4)
         vlocal(4)=vflav(5)
         vlocal(5)=vflav(4)
         call i2MCFM_2_POWHEG(plocal,vlocal,res)
      endif
      virtual_DR=res(0)/(st_alpha/(2*pi))

c from dimensional reduction to dimensional regularization
      call from_DR_to_CDR(p,vflav,virtual_DR,virtual)

c The following code is used to check the double pole
c in the virtual amplitude. Must set in
c set_interface_MCFM.f:      ret_poles=.true.
      if(polecheck) then
         call setborn(p,vflav,born,bornjk,bmunu)
         if(vflav(1).eq.0) then
            c1=3
         else
            c1=4d0/3
         endif
         if(vflav(2).eq.0) then
            c2=3
         else
            c2=4d0/3
         endif
         if(vflav(4).eq.0) then
            c4=3
         else
            c4=4d0/3
         endif
         if(vflav(5).eq.0) then
            c5=3
         else
            c5=4d0/3
         endif
         tmp=(res(2)/(st_alpha/(2*pi)))/(born*(-(c1+c2+c4+c5)))
         if(abs(tmp-1).gt.1d-4) then
            write(*,*) ' Pole check: ',vflav,
     1           tmp,virtual_DR
         endif
      endif

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
