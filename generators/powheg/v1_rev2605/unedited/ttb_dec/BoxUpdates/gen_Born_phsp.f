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
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      integer j,kres
      real * 8 q0,mrec2,pj(0:3)
      real * 8 dotp
      logical valid_emitter
      external valid_emitter

      do j=0,nlegborn
         if(valid_emitter(j)) then
            if(j.gt.2) then
               if(flg_withresrad) then
                  kres=flst_bornres(j,1)
               else
                  kres=-1
               endif
               if(kres.gt.0) then
                  call boost2reson(kn_cmpborn(:,kres),1,
     1                 kn_cmpborn(:,j),pj)
                  q0=sqrt(dotp(kn_cmpborn(:,kres),kn_cmpborn(:,kres)))
               else
                  pj=kn_cmpborn(:,j)
                  q0=2*kn_cmpborn(0,1)
               endif
               mrec2=(q0-pj(0))**2-pj(1)**2-pj(2)**2-pj(3)**2
               kn_csimax_arr(j)=1-mrec2/q0**2
            endif
         else
            kn_csimax_arr(j)=-1
         endif
      enddo
      end

