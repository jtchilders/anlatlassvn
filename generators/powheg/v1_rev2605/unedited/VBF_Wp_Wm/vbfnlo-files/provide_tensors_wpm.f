
      subroutine provide_tensors_wpm(ttype)
      implicit none
      include 'tensor.inc'
      include 'tensor_born.inc'
      include 'tensor_real.inc'
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer i,j,jj,mu,nu
      integer ttype

      if (ttype.eq.1) then

      do i = 1,6   
         wep(i) = wep_born(i) 
         wve(i) = wve_born(i) 
         wmu(i) = wmu_born(i) 
         wvm(i) = wvm_born(i) 
         wp(i)  = wp_born(i) 
         wm(i)  = wm_born(i) 
         wpp(i) = wpp_born(i) 
         wmp(i) = wmp_born(i)          
      enddo !i

      xp = xp_born
      xm = xm_born

      do mu = 0,4
         qp(mu) = qp_born(mu)
         qm(mu) = qm_born(mu)
      enddo  !mu 

      do mu = 0,5
         aww(mu) = aww_born(mu)
         zww(mu) = zww_born(mu)
      enddo   !mu

      do mu = 0,3
         do nu = 0,3
            do j = 1,3
               aaww(mu,nu,j) = aaww_born(mu,nu,j)
               azww(mu,nu,j) = azww_born(mu,nu,j)
               zaww(mu,nu,j) = zaww_born(mu,nu,j)
               zzww(mu,nu,j) = zzww_born(mu,nu,j)
               wwww5(mu,nu,j) = wwww5_born(mu,nu,j)
               wwww6(mu,nu,j) = wwww6_born(mu,nu,j)
               do jj = 1,2
                  NCwpa(mu,nu,jj,j) = NCwpa_born(mu,nu,jj,j) 
                  NCwpz(mu,nu,jj,j) = NCwpz_born(mu,nu,jj,j)
                  CCwpa(mu,nu,jj,j) = CCwpa_born(mu,nu,jj,j) 
                  CCwpz(mu,nu,jj,j) = CCwpz_born(mu,nu,jj,j)
                  NCwma(mu,nu,jj,j) = NCwma_born(mu,nu,jj,j) 
                  NCwmz(mu,nu,jj,j) = NCwmz_born(mu,nu,jj,j)
                  CCwma(mu,nu,jj,j) = CCwma_born(mu,nu,jj,j) 
                  CCwmz(mu,nu,jj,j) = CCwmz_born(mu,nu,jj,j) 
               enddo !jj
            enddo !j
         enddo !nu
      enddo !mu

      elseif (ttype.eq.3) then !real-emission type

      do i = 1,6   
         wep(i) = wep_real(i) 
         wve(i) = wve_real(i) 
         wmu(i) = wmu_real(i) 
         wvm(i) = wvm_real(i) 
         wp(i)  = wp_real(i) 
         wm(i)  = wm_real(i) 
         wpp(i) = wpp_real(i) 
         wmp(i) = wmp_real(i)          
      enddo !i

      xp = xp_real
      xm = xm_real

      do mu = 0,4
         qp(mu) = qp_real(mu)
         qm(mu) = qm_real(mu)
      enddo  !mu 

      do mu = 0,5
         aww(mu) = aww_real(mu)
         zww(mu) = zww_real(mu)
      enddo   !mu

      do mu = 0,3
         do nu = 0,3
            do j = 1,3
               aaww(mu,nu,j) = aaww_real(mu,nu,j)
               azww(mu,nu,j) = azww_real(mu,nu,j)
               zaww(mu,nu,j) = zaww_real(mu,nu,j)
               zzww(mu,nu,j) = zzww_real(mu,nu,j)
               wwww5(mu,nu,j) = wwww5_real(mu,nu,j)
               wwww6(mu,nu,j) = wwww6_real(mu,nu,j)
               do jj = 1,2
                  NCwpa(mu,nu,jj,j) = NCwpa_real(mu,nu,jj,j) 
                  NCwpz(mu,nu,jj,j) = NCwpz_real(mu,nu,jj,j)
                  CCwpa(mu,nu,jj,j) = CCwpa_real(mu,nu,jj,j) 
                  CCwpz(mu,nu,jj,j) = CCwpz_real(mu,nu,jj,j)
                  NCwma(mu,nu,jj,j) = NCwma_real(mu,nu,jj,j) 
                  NCwmz(mu,nu,jj,j) = NCwmz_real(mu,nu,jj,j)
                  CCwma(mu,nu,jj,j) = CCwma_real(mu,nu,jj,j) 
                  CCwmz(mu,nu,jj,j) = CCwmz_real(mu,nu,jj,j) 
               enddo !jj
            enddo !j
         enddo !nu
      enddo !mu

      else
         stop 'bad ttype in provide_tensors'
      endif   

      return
      
      end
