
      subroutine vtoww_born_reset
      implicit none
      include 'tensor_born.inc'
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer i,j,jj,mu,nu


      do i = 1,6   
         wep_born(i) = zero 
         wve_born(i) = zero 
         wmu_born(i) = zero 
         wvm_born(i) = zero 
         wp_born(i)  = zero 
         wm_born(i)  = zero 
         wpp_born(i) = zero 
         wmp_born(i) = zero 
      enddo !i

      xp_born = zero 
      xm_born = zero 

      do mu = 0,4
         qp_born(mu) = zero 
         qm_born(mu) = zero 
      enddo  !mu 

      do mu = 0,5
         aww_born(mu) = zero
         zww_born(mu) = zero
      enddo   

      do mu = 0,3
         do nu = 0,3
            do j = 1,3
               aaww_born(mu,nu,j) = zero
               azww_born(mu,nu,j) = zero
               zaww_born(mu,nu,j) = zero
               zzww_born(mu,nu,j) = zero
               wwww5_born(mu,nu,j) = zero
               wwww6_born(mu,nu,j) = zero
               do jj = 1,2
                  NCwpa_born(mu,nu,jj,j) = zero 
                  NCwpz_born(mu,nu,jj,j) = zero
                  CCwpa_born(mu,nu,jj,j) = zero 
                  CCwpz_born(mu,nu,jj,j) = zero
                  NCwma_born(mu,nu,jj,j) = zero 
                  NCwmz_born(mu,nu,jj,j) = zero
                  CCwma_born(mu,nu,jj,j) = zero 
                  CCwmz_born(mu,nu,jj,j) = zero 
               enddo
            enddo
         enddo
      enddo
      return
      
      end

c==============================

      subroutine vtoww_real_reset
      implicit none
      include 'tensor_real.inc'
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer i,j,jj,mu,nu

      do i = 1,6   
         wep_real(i) = zero 
         wve_real(i) = zero 
         wmu_real(i) = zero 
         wvm_real(i) = zero 
         wp_real(i)  = zero 
         wm_real(i)  = zero 
         wpp_real(i) = zero 
         wmp_real(i) = zero 
      enddo !i

      xp_real = zero 
      xm_real = zero 

      do mu = 0,4
         qp_real(mu) = zero 
         qm_real(mu) = zero 
      enddo  !mu 

      do mu = 0,5
         aww_real(mu) = zero
         zww_real(mu) = zero
      enddo   
      do mu = 0,3
         do nu = 0,3
            do j = 1,3
               aaww_real(mu,nu,j) = zero
               azww_real(mu,nu,j) = zero
               zaww_real(mu,nu,j) = zero
               zzww_real(mu,nu,j) = zero
               wwww5_real(mu,nu,j) = zero
               wwww6_real(mu,nu,j) = zero
               do jj = 1,2
                  NCwpa_real(mu,nu,jj,j) = zero 
                  NCwpz_real(mu,nu,jj,j) = zero
                  CCwpa_real(mu,nu,jj,j) = zero 
                  CCwpz_real(mu,nu,jj,j) = zero
                  NCwma_real(mu,nu,jj,j) = zero 
                  NCwmz_real(mu,nu,jj,j) = zero
                  CCwma_real(mu,nu,jj,j) = zero 
                  CCwmz_real(mu,nu,jj,j) = zero 
               enddo
            enddo
         enddo
      enddo
      return
      
      end
