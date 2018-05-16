      subroutine gen_born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 xborn(ndiminteg-3)
      call born_phsp(xborn)
      call compute_csimax_fsr
      end


      subroutine compute_csimax_fsr
      implicit none
c Compute csimax for all possible final state emitters;
c for initial state emitters it is not possible, since
c csimax depends upon y in this case.
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      integer j
      real * 8 q0,mrec2,m
      logical valid_emitter
      external valid_emitter
      do j=0,nlegborn
         if(valid_emitter(j)) then
            if(j.gt.2) then
               m=kn_masses(j)
               q0=2*kn_cmpborn(0,1)
               mrec2=abs((q0-kn_cmpborn(0,j))**2
     #               -kn_cmpborn(1,j)**2-kn_cmpborn(2,j)**2
     #               -kn_cmpborn(3,j)**2)
               if(m.eq.0) then
                  kn_csimax_arr(j)=1-mrec2/q0**2
               else
                  kn_csimax_arr(j)=1-(sqrt(mrec2)+m)**2/q0**2
               endif
            endif
         else
            kn_csimax_arr(j)=-1
         endif
      enddo
      end

