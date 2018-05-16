!---- Driver routine for H(->bb/tautau) + j in i2MCFM 
!---- C. Williams August 2011

      subroutine gg_hg_vi(p,res) 
      implicit none 
      include 'MCFM_Include/interface_settings.f' 
      include 'MCFM_Include/constants.f' 
      include 'MCFM_Include/scheme.f'
      double precision p(mxpart,4),res(0:3)
      double precision eps_0,eps_1,eps_m1
      double precision Born(-nf:nf,-nf:nf)
      integer i

      eps_0=0d0
      eps_1=0d0
      eps_m1=0d0
      
      do i=0,3
         res(i)=0d0
      enddo
      
      if(ret_poles) then 
         call set_epinv(0d0)                 
         call gg_hg_eval_v(p,eps_0)
!-----Calculate whole matrix element with eps = 1      
         call set_epinv(1d0) 
         call gg_hg_eval_v(p,eps_1) 
!-----Calculate whole matrix element with eps = -1      
         call set_epinv(-1d0) 
         call gg_hg_eval_v(p,eps_m1) 
              
       
         
!------Fill res 
!---- Finite piece is simply eps=0 
         res(0)=eps_0
!---- 1/eps pole piece is 1/2(eps_1 - eps_m1)
         res(1)=half*(eps_1-eps_m1)
!-----1/eps^2 pole piece is 1/2(eps_1+eps_m1)-eps_0
         res(2)=half*(eps_1+eps_m1)-eps_0
      
      else
         call set_epinv(0d0) 
         call gg_hg_eval_v(p,eps_0) 
         res(0)=eps_0
         res(1)=0d0
         res(2)=0d0
      endif
  
!  call gg_hg(p,Born)
!  res(3)=Born(ret1,ret2)
      res(3)=0d0 

      return 
      end subroutine 
