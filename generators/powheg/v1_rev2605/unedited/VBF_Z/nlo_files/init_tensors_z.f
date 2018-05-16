      
      subroutine vtoll_reset
      implicit none
      include 'tensorl-hel.inc'
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,h,mu,nu
      
      do mu = 0,3	 
         do nu = 0,3
            do h = -1,1
	       do j = 1,3
	       
                  aaee(mu,nu,j,h) = zero
                  azee(mu,nu,j,h) = zero
                  zaee(mu,nu,j,h) = zero
                  zzee(mu,nu,j,h) = zero
                  CCee(mu,nu,j,h) = zero
                  CCee6(mu,nu,j,h) = zero 
	       
	       enddo            !j
            enddo               !h
         enddo                  !nu
      enddo                     !mu
      
c      print*,' vv->ee tensors initialized to ',zero 
      return
      end
