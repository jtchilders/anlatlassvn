      subroutine qqb_wpwp_qqb_vpolesonly(p,msqv,ch)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     . p(mxpart,4),msqB_str(-nf:nf,-nf:nf,3),dot
      double complex virt_q_bq(3), virt_bq_q(3) 
      double complex virt_q_q(3), virt_bq_bq(3) 
      double complex virtcom
      double complex xl12,xl27,xl78,xl18,xl17,xl28
      double complex lnrat 
      double precision, parameter :: Nc = 3d0 
      character ch*3
      integer j,k 
      
      xl12=-Lnrat(musq,abs(two*dot(p,1,2))) ! t 
      xl78=-Lnrat(musq,abs(two*dot(p,7,8))) ! t 
      xl18=-Lnrat(musq,abs(two*dot(p,1,8))) ! s 
      xl27=-Lnrat(musq,abs(two*dot(p,2,7))) ! s 
      xl28=-Lnrat(musq,abs(two*dot(p,2,8))) ! u 
      xl17=-Lnrat(musq,abs(two*dot(p,1,7))) ! u 

      virtcom = -4d0*Cf*epinv2*epinv
     .     + epinv*(-Cf*6d0 + Nc*(11d0/3d0))

      call qqb_wpwp_qqb(p,msq,ch)
      call qqb_wpwp_qqb_str(p,msqB_str,ch)

!     -- chn qqb 
      virt_q_bq(1) = ason2pi * (virtcom 
     .     + epinv*((-4d0*Cf+2d0*Nc)*(xl17+xl28)
     .     +        ( 4d0*Cf-1d0*Nc)*(xl12+xl78) 
     .     +        ( 2d0*Cf-1d0*Nc)*(xl18+xl27)
     .     ))

      virt_q_bq(2) = ason2pi * (virtcom 
     .     + epinv*((-4d0*Cf+2d0*Nc)*(xl17+xl28)
     .     +        ( 4d0*Cf-1d0*Nc)*(xl18+xl27)
     .     +        ( 2d0*Cf-1d0*Nc)*(xl12+xl78) 
     .     ))


      virt_q_bq(3) = ason2pi * (virtcom 
     .     + epinv*((-2d0*Cf+2d0*Nc)*(xl17+xl28)
     .     +        ( 2d0*Cf-1d0*Nc)*(xl18+xl27)
     .     +        ( 2d0*Cf-1d0*Nc)*(xl12+xl78) 
     .     ))


!     -- chn qbq 
      virt_bq_q(1) = ason2pi * (virtcom 
     .     + epinv*((-4d0*Cf+2d0*Nc)*(xl27+xl18)
     .     +        ( 4d0*Cf-1d0*Nc)*(xl12+xl78)
     .     +        ( 2d0*Cf-1d0*Nc)*(xl28+xl17)
     .     ))

      virt_bq_q(2) = ason2pi * (virtcom 
     .     + epinv*((-4d0*Cf+2d0*Nc)*(xl27+xl18)
     .     +        ( 4d0*Cf-1d0*Nc)*(xl28+xl17)
     .     +        ( 2d0*Cf-1d0*Nc)*(xl12+xl78) 
     .     ))


      virt_bq_q(3) = ason2pi * (virtcom 
     .     + epinv*((-2d0*Cf+2d0*Nc)*(xl27+xl18)
     .     +        ( 2d0*Cf-1d0*Nc)*(xl28+xl17)
     .     +        ( 2d0*Cf-1d0*Nc)*(xl12+xl78) 
     .     ))

!     -- chn qqq 
      virt_q_q(1) = ason2pi * (virtcom 
     .     + epinv*((-4d0*Cf+2d0*Nc)*(xl12+xl78)
     .     +        ( 4d0*Cf-1d0*Nc)*(xl18+xl27)
     .     +        ( 2d0*Cf-1d0*Nc)*(xl17+xl28)
     .     ))

      virt_q_q(2) = ason2pi * (virtcom 
     .     + epinv*((-4d0*Cf+2d0*Nc)*(xl12+xl78)  
     .     +        ( 4d0*Cf-1d0*Nc)*(xl17+xl28)  
     .     +        ( 2d0*Cf-1d0*Nc)*(xl27+xl18)  
     .     ))

      virt_q_q(3) = ason2pi * (virtcom 
     .     + epinv*((-2d0*Cf+2d0*Nc)*(xl12+xl78)  
     .     +        ( 2d0*Cf-1d0*Nc)*(xl18+xl27)  
     .     +        ( 2d0*Cf-1d0*Nc)*(xl17+xl28)  
     .     ))


!     -- chn qbb 
      virt_bq_bq(1) = ason2pi * (virtcom
     .     + epinv*((-4d0*Cf+2d0*Nc)*(xl78+xl12)
     .     +        ( 4d0*Cf-1d0*Nc)*(xl27+xl18)
     .     +        ( 2d0*Cf-1d0*Nc)*(xl17+xl28)
     .     ))

      virt_bq_bq(2) = ason2pi * (virtcom 
     .     + epinv*((-4d0*Cf+2d0*Nc)*(xl12+xl78)  
     .     +        ( 4d0*Cf-1d0*Nc)*(xl17+xl28)  
     .     +        ( 2d0*Cf-1d0*Nc)*(xl27+xl18)  
     .     ))

      virt_bq_bq(3) = ason2pi * (virtcom 
     .     + epinv*((-2d0*Cf+2d0*Nc)*(xl12+xl78)  
     .     +        ( 2d0*Cf-1d0*Nc)*(xl18+xl27)  
     .     +        ( 2d0*Cf-1d0*Nc)*(xl17+xl28)  
     .     ))


c-----------------------------------------------
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf

      if (j.gt.0.and.k.lt.0) then 
         if ((j == 2 .and. k == -3) .or. (j == 4 .and. k == -1)) then 
            msqv(j,k) = msqB_str(j,k,1)*virt_q_bq(1) ! s
         else
            msqv(j,k) = 
     .           + msqB_str(j,k,2)*virt_q_bq(1) ! t 
     .           + msqB_str(j,k,1)*virt_q_bq(2) ! s
     .           + msqB_str(j,k,3)*virt_q_bq(3) ! s*t 
         endif

      endif 

      if (j.lt.0.and.k.gt.0) then 
         if ((j == -3 .and. k == 2) .or. (j == -1 .and. k == 4)) then  
            msqv(j,k) = msqB_str(j,k,1)*virt_bq_q(1) ! s
         else
            msqv(j,k) = 
     .           + msqB_str(j,k,2)*virt_bq_q(1)
     .           + msqB_str(j,k,1)*virt_bq_q(2)
     .           + msqB_str(j,k,3)*virt_bq_q(3)
         endif
      endif 


      if (j.gt.0.and.k.gt.0) then 
         if (j == k) then 
            msqv(j,k) = 
     .           + msqB_str(j,k,1)*virt_q_q(1)
     .           + msqB_str(j,k,2)*virt_q_q(2)
     .           + msqB_str(j,k,3)*virt_q_q(3)
         else
            msqv(j,k) = msqB_str(j,k,1)*virt_q_q(1)
         endif

      endif 

      if (j.lt.0.and.k.lt.0) then 
         if (j == k) then 
            msqv(j,k) = 
     .           + msqB_str(j,k,1)*virt_bq_bq(1)
     .           + msqB_str(j,k,2)*virt_bq_bq(2)
     .           + msqB_str(j,k,3)*virt_bq_bq(3)
         else
            msqv(j,k) = msqB_str(j,k,1)*virt_bq_bq(1)
         endif

      endif 

      enddo
      enddo


      return
      end


