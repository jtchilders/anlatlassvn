c
      subroutine setreal(p,fermion_flav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'cvecbos.h'
      integer nleg
      parameter (nleg=nlegreal)
      real * 8 p(0:3,nleg)
      integer fermion_flav(nleg)
      real * 8 amp2
      real*8 col_dec,had_sum
      integer ttype

      integer j
      real*8 pwsum
      logical, save :: firstreal = .true. 
      real*8 pidr 
      common /ps_idr/ pidr    

c----------------------------------------------------

      if (firstreal) then
         pidr = 0d0
         firstreal = .false.
      endif     
      
      pwsum = 0d0
      do j = 1,nleg
            pwsum = pwsum+real(j)*p(0,j)
      enddo        

      had_sum = 1d0
      col_dec = 1d0
      if (decmode_slp) then
         col_dec = 3d0
         if (abs(vdecaymodeWp).eq.107) had_sum = 2d0
      elseif (decmode_slm) then
         col_dec = 3d0
         if (abs(vdecaymodeWm).eq.107) had_sum = 2d0
      endif   

      if ( pwsum.ne.pidr) then ! new PS point -> compute tensors
         if (decmode_slp) then
           call compute_tensors_wpm_real_slp(p) 
         elseif (decmode_slm) then
           call compute_tensors_wpm_real_slm(p) 
         else
           call compute_tensors_wpm_real(p) 
         endif   
         pidr = 0d0
         do j = 1,nleg
            pidr = pidr+real(j)*p(0,j)
         enddo  
      endif
      ttype = 3
      call provide_tensors_wpm(ttype)
      call compreal_wpm_ew(p,fermion_flav,amp2)

      amp2 = amp2*col_dec*had_sum

c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2d0*pi))

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compreal_wpm_ew(pin,bflav,amp2)
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
      include 'cvecbos.h'
c
      integer nlegs
      parameter (nlegs=nlegreal)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs)
      real*8 amp2 
c
c vbfnlo stuff:
      include 'global.inc'
      integer bos
      real*8 p(0:3,np), v(0:3,nv)
      real*8 pbar(0:3,4+nv),qbar(0:4)
      real*8 polcol,polcolq,polcolg
      real*8 res(3)

      real*8 N ! color factors
      parameter(N=3d0)

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare local variables
c
      integer i,mu
      integer FSIGN(4+nv),gsign,physToDiag(5)
      
      real*8 nans(2,2,3),cans(2,3)
      logical nc_type
      integer k,icc

      integer ftype(1:9)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      bos = 34

      polcol = 0d0
      polcolq = 1d0/(4d0*N**2)
      polcolg = 1d0/(4d0*N*(N**2-1))

      ftype(1:9) = 0
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         do i = 1,4
            v(mu,i) = pin(mu,i+2)
         enddo !i  

         if (bflav(1)*bflav(2).ne.0) then ! final-state gluon 
            if (bflav(9).eq.0) then
               p(mu,3) = pin(mu,7)
               p(mu,4) = pin(mu,8) 
               p(mu,5) = pin(mu,9)     
            elseif (bflav(7).eq.0) then
               p(mu,3) = pin(mu,8)
               p(mu,4) = pin(mu,9) 
               p(mu,5) = pin(mu,7)  
            elseif (bflav(8).eq.0) then
               p(mu,3) = pin(mu,7)
               p(mu,4) = pin(mu,9) 
               p(mu,5) = pin(mu,8)  
            endif
         else   ! initial-state gluon 
            p(mu,3) = pin(mu,7)
            p(mu,4) = pin(mu,8) 
            p(mu,5) = pin(mu,9) 
         endif ! fin/in state gluon
      enddo ! mu

      if (bflav(1).gt.0.and.bflav(2).gt.0) then

C*******************  q1 q3 ---> q2 q4 g W W   **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1

      polcol = polcolq

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)

      if ((ftype(1).eq.ftype(7))) then
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = -3*(bflav(1)/abs(bflav(1)))+4
         k = 7-ftype(icc)
      endif   

      elseif (bflav(1).gt.0.and.bflav(2).lt.0) then
c            
C******************* q1 qb4 ---> q2 qb3 g W W   **********************
      
      physToDiag(1)=1    
      physToDiag(2)=4
      physToDiag(3)=2    
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign    =  1
      
      polcol = polcolq
      
c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)

      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = -3*(bflav(1)/abs(bflav(1)))+4
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).lt.0.and.bflav(2).gt.0) then
      
C******************* qbar2 q3 ---> qbar1 q4 g W W   **********************
      
      physToDiag(1)=2    
      physToDiag(2)=3
      physToDiag(3)=1    
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1
      gsign    =  1
      
      polcol = polcolq

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)

      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = -3*(bflav(1)/abs(bflav(1)))+4
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).lt.0.and.bflav(2).lt.0) then

C*******************  qbar2 qb4 ---> qbar1 qb3 g W W   **********************

      physToDiag(1)=2    
      physToDiag(2)=4
      physToDiag(3)=1    
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) = -1
      gsign    =  1
      
      polcol = polcolq

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)

      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = -3*(bflav(1)/abs(bflav(1)))+4
         k = 7-ftype(icc)
      endif !nc/cc  
     
      elseif (bflav(1).eq.0.and.bflav(2).gt.0.and.bflav(7).gt.0) then
         
c*******************  g q ---> q q qb W W   **********************

      physToDiag(1)=5          
      physToDiag(2)=3           
      physToDiag(3)=2           
      physToDiag(4)=4
      physToDiag(5)=1
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign    = -1
      
      polcol = polcolg  

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(9)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)

      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc
         icc = 3*(bflav(2)/abs(bflav(2)))+5
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).eq.0.and.bflav(2).gt.0.and.bflav(7).lt.0) then
         
c*******************  g q ---> qb q q W W   **********************
c
      physToDiag(1)=5          
      physToDiag(2)=3           
      physToDiag(3)=1           
      physToDiag(4)=4
      physToDiag(5)=2
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) =  1
      fsign(4) =  1
      gsign    = -1
      
      polcol = polcolg

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(9)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)     

      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = 3*(bflav(2)/abs(bflav(2)))+5
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).gt.0.and.bflav(2).eq.0.
     &    and.bflav(7).gt.0.and.bflav(8).gt.0) then
      
C*******************  q g ---> q q qb W W   **********************
      
      physToDiag(1)=1             
      physToDiag(5)=3            
      physToDiag(3)=2             
      physToDiag(2)=5
      physToDiag(4)=4
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol = polcolg
c
c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(9)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)
      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = -3*(bflav(1)/abs(bflav(1)))+4
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).gt.0.and.bflav(2).eq.0.and.
     &        bflav(7).gt.0.and.bflav(8).lt.0) then
      
C*******************  q g ---> q qb q W W   **********************

      physToDiag(1)=1             
      physToDiag(4)=3            
      physToDiag(3)=2             
      physToDiag(2)=5
      physToDiag(5)=4
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol = polcolg

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(9)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)

      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = -3*(bflav(1)/abs(bflav(1)))+4
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).eq.0.and.bflav(2).lt.0.and.bflav(7).gt.0) then
        
C*******************  g qbar ---> q qb qb W W  **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=2
      physToDiag(4)=3
      physToDiag(5)=1
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign    = -1

      polcol = polcolg  

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(9)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)

      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = 3*(bflav(2)/abs(bflav(2)))+5
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).eq.0.and.bflav(2).lt.0.and.bflav(7).lt.0) then
        
C*******************  g qbar ---> qb qb q W W  **********************

      physToDiag(1)=5
      physToDiag(2)=4
      physToDiag(3)=1
      physToDiag(4)=3
      physToDiag(5)=2
c
      fsign(1) = -1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1
      gsign    = -1

      polcol = polcolg

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(9)),2)         
      ftype(2) = 2-mod(abs(bflav(2)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)

      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = 3*(bflav(2)/abs(bflav(2)))+5
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).lt.0.and.bflav(2).eq.0.and.bflav(8).lt.0) then
 
C*******************  qbar2 g ---> qbar1 qb3 q4 W W   **********************
c
      physToDiag(1)=2             
      physToDiag(5)=4!3            
      physToDiag(3)=1             
      physToDiag(2)=5
      physToDiag(4)=3!4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      
      polcol = polcolg

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(9)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)

      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = -3*(bflav(1)/abs(bflav(1)))+4
         k = 7-ftype(icc)
      endif !nc/cc  

      elseif (bflav(1).lt.0.and.bflav(2).eq.0.and.bflav(8).gt.0) then
 
C*******************  qbar2 g ---> qbar1 q4 q3bar W W   **********************
c
      physToDiag(1)=2             
      physToDiag(5)=3            
      physToDiag(3)=1             
      physToDiag(2)=5
      physToDiag(4)=4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      
      polcol = polcolg

c up- or down-type quark:      
      ftype(1) = 2-mod(abs(bflav(1)),2)         
      ftype(2) = 2-mod(abs(bflav(9)),2)      
      ftype(7) = 2-mod(abs(bflav(7)),2)   
      ftype(8) = 2-mod(abs(bflav(8)),2)

      if ((ftype(1).eq.ftype(7))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = -3*(bflav(1)/abs(bflav(1)))+4
         k = 7-ftype(icc)
      endif !nc/cc  

      else
         
         print*,'this flav combination is not implemented'
         print*,'bflav=',bflav

      endif
         
C*****************  end of process evaluation  **********************

      do mu = 0,3
         do i = 1,5
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
         qbar(mu) = pbar(mu,5)
      enddo 
      qbar(4) = 0d0

c lepton momenta:
      do mu = 0,3             ! kinematics for 4-lepton final state
         pbar(mu,5) = v(mu,2) ! l+
         pbar(mu,6) = v(mu,1) ! nu
         pbar(mu,7) = v(mu,3) ! nubar
         pbar(mu,8) = v(mu,4) ! l-
      enddo
      fsign(5) = -1
      fsign(6) = 1
      fsign(7) = -1
      fsign(8) = 1

c k-ordering is: uucc,uuss,ddcc,ddss,udsc,ducs

      call qqwwqqj_channel(pbar,fsign,qbar,gsign,k,res(1))
      
      amp2 = res(1)*polcol
c symmetry factor for leptons:
      amp2 = amp2*wsymfact

      return
      end
