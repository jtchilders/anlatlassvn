ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     compute leptonic tensors for VBF pp->W+W-jj 
c
      subroutine compute_tensors_wpm_real(pin)

c     fully leptonic decays
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
c
      integer nlegs
      parameter (nlegs=nlegreal)
      real*8 pin(0:3,nlegs)
c
c vbfnlo stuff:
      include 'global.inc'
      include 'coupl.inc'
      include 'tensor_real.inc'
      include 'higgs_graphs.h'
      real*8 p(0:3,np), v(0:3,nv)
      real*8 q12(0:4,3),q34(0:4,3),qww(0:3)

      complex*16 zero
      parameter (zero=(0d0,0d0))

      double complex dotrc
      external dotrc
c
c declare local variables
      integer i,j,mu
      integer ttype
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ttype = 3 ! tensors are real-type

      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         do i = 1,4
            v(mu,i) = pin(mu,i+2)
         enddo !i     
         p(mu,3) = pin(mu,7)
         p(mu,4) = pin(mu,8) 
         p(mu,5) = pin(mu,9)  
      enddo ! mu

c compute invariants:
      call calc_invariants(p,v,q12,q34,1)
c
c reset leptonic tensors:
      call vtoww_real_reset
c compute leptonic tensors:
      CALL IXXXXX(v(0,2),ZERO ,1,-1,wep_real)  !e+       
      CALL OXXXXX(v(0,1),ZERO ,-1,1,wve_real)  !ve 
      CALL OXXXXX(v(0,4),ZERO ,-1,1,wmu_real)  !mu-      
      CALL IXXXXX(v(0,3),ZERO ,1,-1,wvm_real)  !vm~
      CALL JIOXXX(wep_real,wve_real,GWF ,WMASS,WWIDTH,wp_real) !W+
      CALL JIOXXX(wvm_real,wmu_real,GWF ,WMASS,WWIDTH,wm_real) !W-

      do mu = 0,3
         qp_real(mu) = v(mu,1)+v(mu,2)
         qm_real(mu) = v(mu,3)+v(mu,4)
         qww(mu) = qp_real(mu) + qm_real(mu)
      enddo
      qp_real(4)=qp_real(0)**2-qp_real(1)**2-
     &           qp_real(2)**2-qp_real(3)**2
      qm_real(4)=qm_real(0)**2-qm_real(1)**2-
     &           qm_real(2)**2-qm_real(3)**2

c determine "reduced" W+- polarization vectors 
c wpp(mu) = wp(mu) - xp*qp(mu) and wmp(mu) = wm(mu) - xm*qm(mu)
            xp_real = 0d0
            xm_real = 0d0
            do mu = 0,3
               wpp_real(mu+1) = wp_real(mu+1) - xp_real*qp_real(mu)
               wmp_real(mu+1) = wm_real(mu+1) - xm_real*qm_real(mu)
            enddo
            do mu = 5,6
               wpp_real(mu) = wp_real(mu)
               wmp_real(mu) = wm_real(mu)
            enddo

c leptonic tensors
C for W+W-
      if (.not.higgs_only) then
         call atoww(v,aww_real,ttype)
         call ztoww(v,zww_real,ttype)
      endif

      do j = 2,3

      if (.not.higgs_only) then
         call aatoww(q12(0,j),q34(0,j),v,aaww_real(0,0,j),ttype)
         call aztoww(q12(0,j),q34(0,j),v,azww_real(0,0,j),ttype)
         call aztoww(q34(0,j),q12(0,j),v,zaww_real(0,0,j),ttype)
      endif
      call zztoww(q12(0,j),q34(0,j),v,zzww_real(0,0,j),ttype)
      call wwtoww(q12(0,j),q34(0,j),v,wwww6_real(0,0,j),ttype)!q12=W-
      call wwtoww(q34(0,j),q12(0,j),v,wwww5_real(0,0,j),ttype)!q12=W+

      if (.not.higgs_only) then
C for WV --> e+ nu_e
c NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
      call WVtoWP(2,q34(0,j),v,NCwpa_real(0,0,1,j),
     &        NCwpz_real(0,0,1,j),ttype) !emit W- on upper
      call WVtoWP(2,q12(0,j),v,NCwpa_real(0,0,2,j),
     &        NCwpz_real(0,0,2,j),ttype)!emit W- on lower
      call WVtoWP(1,q34(0,j),v,CCwpa_real(0,0,1,j),
     &        CCwpz_real(0,0,1,j),ttype)!emit W- on upper
      call WVtoWP(1,q12(0,j),v,CCwpa_real(0,0,2,j),
     &        CCwpz_real(0,0,2,j),ttype)!emit W- on lower
C for WV --> mu- nu_mu
      call WVtoWM(2,q34(0,j),v,NCwma_real(0,0,1,j),
     &        NCwmz_real(0,0,1,j),ttype)!emit W+ on upper
      call WVtoWM(2,q12(0,j),v,NCwma_real(0,0,2,j),
     &        NCwmz_real(0,0,2,j),ttype)!emit W+ on lower
      call WVtoWM(1,q34(0,j),v,CCwma_real(0,0,1,j),
     &        CCwmz_real(0,0,1,j),ttype)!emit W+ on upper
      call WVtoWM(1,q12(0,j),v,CCwma_real(0,0,2,j),
     &        CCwmz_real(0,0,2,j),ttype)!emit W+ on lower
         
      endif ! H->WW only

      enddo !j

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     compute tensors for VBF pp->W+W-jj with semi-leptonic decay
c
      subroutine compute_tensors_wpm_real_slp(pin)
c
c     "slp": (e+ve) is replaced with (ud~) in final state
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
c
      integer nlegs
      parameter (nlegs=nlegreal)
      real*8 pin(0:3,nlegs)
c
c vbfnlo stuff:
      include 'global.inc'
      include 'coupl.inc'
      include 'tensor_real.inc'
      include 'higgs_graphs.h'
      real*8 p(0:3,np), v(0:3,nv)
      real*8 q12(0:4,3),q34(0:4,3),qww(0:3)

      complex*16 zero
      parameter (zero=(0d0,0d0))

      double complex dotrc
      external dotrc
c
c declare local variables
      integer i,j,mu
      integer ttype
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ttype = 3 ! tensors are real-emission type
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         do i = 1,4
            v(mu,i) = pin(mu,i+2)
         enddo !i     
         p(mu,3) = pin(mu,7)
         p(mu,4) = pin(mu,8) 
         p(mu,5) = pin(mu,9)  
      enddo ! mu

c compute invariants:
      call calc_invariants(p,v,q12,q34,1)
c
c reset (semi-)leptonic tensors:
      call vtoww_real_reset
c compute (semi-)leptonic tensors:
      CALL IXXXXX(v(0,2),ZERO ,1,-1,wve_real) !b~               
      CALL OXXXXX(v(0,1),ZERO ,-1,1,wep_real) !t        
      CALL OXXXXX(v(0,4),ZERO ,-1,1,wmu_real) !mu-      
      CALL IXXXXX(v(0,3),ZERO ,1,-1,wvm_real) !vm~

      CALL JIOXXX(wve_real,wep_real,GWF ,WMASS,WWIDTH,wp_real) !W+
      CALL JIOXXX(wvm_real,wmu_real,GWF ,WMASS,WWIDTH,wm_real) !W-

      do mu = 0,3
         qp_real(mu) = v(mu,1)+v(mu,2)
         qm_real(mu) = v(mu,3)+v(mu,4)
         qww(mu) = qp_real(mu) + qm_real(mu)
      enddo
      qp_real(4)=qp_real(0)**2-qp_real(1)**2-
     &           qp_real(2)**2-qp_real(3)**2
      qm_real(4)=qm_real(0)**2-qm_real(1)**2-
     &           qm_real(2)**2-qm_real(3)**2

c determine "reduced" W+- polarization vectors 
c wpp(mu) = wp(mu) - xp*qp(mu) and wmp(mu) = wm(mu) - xm*qm(mu)
      xp_real = 0d0
      xm_real = 0d0
      do mu = 0,3
         wpp_real(mu+1) = wp_real(mu+1) - xp_real*qp_real(mu)
         wmp_real(mu+1) = wm_real(mu+1) - xm_real*qm_real(mu)
      enddo
      do mu = 5,6
         wpp_real(mu) = wp_real(mu)
         wmp_real(mu) = wm_real(mu)
      enddo

c (semi-)leptonic tensors
C for W+W-
        if (.not.higgs_only) then
           call atoww_slp(v,aww_real,ttype)
           call ztoww_slp(v,zww_real,ttype)
        endif

        do j = 2,3

        if (.not.higgs_only) then
           call aatoww_slp(q12(0,j),q34(0,j),v,aaww_real(0,0,j),ttype)
           call aztoww_slp(q12(0,j),q34(0,j),v,azww_real(0,0,j),ttype)
           call aztoww_slp(q34(0,j),q12(0,j),v,zaww_real(0,0,j),ttype)
        endif
        call zztoww_slp(q12(0,j),q34(0,j),v,zzww_real(0,0,j),ttype)
        call wwtoww_slp(q12(0,j),q34(0,j),v,wwww6_real(0,0,j),ttype)!q12=W-
        call wwtoww_slp(q34(0,j),q12(0,j),v,wwww5_real(0,0,j),ttype)!q12=W+

        if (.not.higgs_only) then
cC for WV --> tb
cc NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
        call WVtoWP_slp(2,q34(0,j),v,NCwpa_real(0,0,1,j),
     &                  NCwpz_real(0,0,1,j),ttype)     !emit W- on upper
        call WVtoWP_slp(2,q12(0,j),v,NCwpa_real(0,0,2,j),
     &                  NCwpz_real(0,0,2,j),ttype)     !emit W- on lower
        call WVtoWP_slp(1,q34(0,j),v,CCwpa_real(0,0,1,j),
     &                  CCwpz_real(0,0,1,j),ttype)     !emit W- on upper
        call WVtoWP_slp(1,q12(0,j),v,CCwpa_real(0,0,2,j),
     &                  CCwpz_real(0,0,2,j),ttype)     !emit W- on lower
C for WV --> mu- nu_mu
        call WVtoWM(2,q34(0,j),v,NCwma_real(0,0,1,j),
     &        NCwmz_real(0,0,1,j),ttype)     !emit W+ on upper
        call WVtoWM(2,q12(0,j),v,NCwma_real(0,0,2,j),
     &        NCwmz_real(0,0,2,j),ttype)     !emit W+ on lower
        call WVtoWM(1,q34(0,j),v,CCwma_real(0,0,1,j),
     &        CCwmz_real(0,0,1,j),ttype)     !emit W+ on upper
        call WVtoWM(1,q12(0,j),v,CCwma_real(0,0,2,j),
     &        CCwmz_real(0,0,2,j),ttype)     !emit W+ on lower
         
        endif ! H->WW only

      enddo !j  

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     compute tensors for VBF pp->W+W-jj with semi-leptonic decay
c
      subroutine compute_tensors_wpm_real_slm(pin)
c
c     "slm": (vm~ mu) is replaced with (u~ d) in final state
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
c
      integer nlegs
      parameter (nlegs=nlegreal)
      real*8 pin(0:3,nlegs)
c
c vbfnlo stuff:
      include 'global.inc'
      include 'coupl.inc'
      include 'tensor_real.inc'
      include 'higgs_graphs.h'
      real*8 p(0:3,np), v(0:3,nv)
      real*8 q12(0:4,3),q34(0:4,3),qww(0:3)

      complex*16 zero
      parameter (zero=(0d0,0d0))

      double complex dotrc
      external dotrc
c
c declare local variables
      integer i,j,mu
      integer ttype
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ttype = 3 ! tensors are real.emission type
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         do i = 1,4
            v(mu,i) = pin(mu,i+2)
         enddo !i     
         p(mu,3) = pin(mu,7)
         p(mu,4) = pin(mu,8) 
         p(mu,5) = pin(mu,9) 
      enddo ! mu

c compute invariants:
      call calc_invariants(p,v,q12,q34,1)
c
c reset (semi-)leptonic tensors:
      call vtoww_real_reset
c compute (semi-)leptonic tensors:
      CALL IXXXXX(v(0,2),ZERO ,1,-1,wep_real)!e+       
      CALL OXXXXX(v(0,1),ZERO ,-1,1,wve_real)!ve 
      CALL OXXXXX(v(0,4),ZERO ,-1,1,wmu_real)!b    
      CALL IXXXXX(v(0,3),ZERO ,1,-1,wvm_real)!t~
      CALL JIOXXX(wep_real,wve_real,GWF ,WMASS,WWIDTH,wp_real)!W+
      CALL JIOXXX(wvm_real,wmu_real,GWF ,WMASS,WWIDTH,wm_real)!W-

      do mu = 0,3
         qp_real(mu) = v(mu,1)+v(mu,2)
         qm_real(mu) = v(mu,3)+v(mu,4)
         qww(mu) = qp_real(mu) + qm_real(mu)
      enddo
      qp_real(4)=qp_real(0)**2-qp_real(1)**2-
     &           qp_real(2)**2-qp_real(3)**2
      qm_real(4)=qm_real(0)**2-qm_real(1)**2-
     &           qm_real(2)**2-qm_real(3)**2

c determine "reduced" W+- polarization vectors 
c wpp(mu) = wp(mu) - xp*qp(mu) and wmp(mu) = wm(mu) - xm*qm(mu)
      xp_real = 0d0
      xm_real = 0d0
      do mu = 0,3
         wpp_real(mu+1) = wp_real(mu+1) - xp_real*qp_real(mu)
         wmp_real(mu+1) = wm_real(mu+1) - xm_real*qm_real(mu)
      enddo
      do mu = 5,6
         wpp_real(mu) = wp_real(mu)
         wmp_real(mu) = wm_real(mu)
      enddo

c (semi-)leptonic tensors
C for W+W-
         if (.not.higgs_only) then
            call atoww_slm(v,aww_real,ttype)
            call ztoww_slm(v,zww_real,ttype)
         endif

         do j = 2,3

         if (.not.higgs_only) then
            call aatoww_slm(q12(0,j),q34(0,j),v,aaww_real(0,0,j),ttype)
            call aztoww_slm(q12(0,j),q34(0,j),v,azww_real(0,0,j),ttype)
            call aztoww_slm(q34(0,j),q12(0,j),v,zaww_real(0,0,j),ttype)
         endif
         call zztoww_slm(q12(0,j),q34(0,j),v,zzww_real(0,0,j),ttype)
         call wwtoww_slm(q12(0,j),q34(0,j),v,wwww6_real(0,0,j),ttype)!q12=W-
         call wwtoww_slm(q34(0,j),q12(0,j),v,wwww5_real(0,0,j),ttype)!q12=W+

         if (.not.higgs_only) then
cC for WV --> tb
cc NCw tensors for NC process (k=1...4), CCw for CC (k=5,6)
         call WVtoWP(2,q34(0,j),v,NCwpa_real(0,0,1,j),
     &               NCwpz_real(0,0,1,j),ttype)     !emit W- on upper
         call WVtoWP(2,q12(0,j),v,NCwpa_real(0,0,2,j),
     &               NCwpz_real(0,0,2,j),ttype)     !emit W- on lower
         call WVtoWP(1,q34(0,j),v,CCwpa_real(0,0,1,j),
     &               CCwpz_real(0,0,1,j),ttype)     !emit W- on upper
         call WVtoWP(1,q12(0,j),v,CCwpa_real(0,0,2,j),
     &               CCwpz_real(0,0,2,j),ttype)     !emit W- on lower
C for WV --> mu- nu_mu
         call WVtoWM_slm(2,q34(0,j),v,NCwma_real(0,0,1,j),
     &               NCwmz_real(0,0,1,j),ttype)     !emit W+ on upper
         call WVtoWM_slm(2,q12(0,j),v,NCwma_real(0,0,2,j),
     &               NCwmz_real(0,0,2,j),ttype)     !emit W+ on lower
         call WVtoWM_slm(1,q34(0,j),v,CCwma_real(0,0,1,j),
     &               CCwmz_real(0,0,1,j),ttype)     !emit W+ on upper
         call WVtoWM_slm(1,q12(0,j),v,CCwma_real(0,0,2,j),
     &               CCwmz_real(0,0,2,j),ttype)     !emit W+ on lower
         
         endif ! H->WW only
         

      enddo !j
   
      return
      end
