      subroutine setvirtual(p,vflav,virtual)

      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'

      real * 8 p(0:3,nlegborn)
      integer vflav(nlegborn)
      real * 8 virtual

      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 bmunu(0:3,0:3,nlegs),bbmunu(0:3,0:3),born,colcf
      real *8 powheginput
      external powheginput 

      logical, save :: firsttime = .true. 

      integer fakevirt
      save fakevirt 

c================================================

      virtual = 0d0

      if (firsttime) then
        fakevirt=powheginput("#fakevirt")
        if (fakevirt == 1) write(*,*) 'WARNING: Using fakevirt !'
        firsttime = .false.
      endif

      if(fakevirt.eq.1) then  

      call compute_tensors_z(p)     
      call compborn_z_ew(p,vflav,born,bbmunu) 
 
      virtual = 0.2d0*born
      
      else

c numbering of momenta is q(1) q(2) -> l+(3)l-(4)q(5)q(6)
c
      call compute_tensors_z(p)    
      call compvirt_z_ew(p,vflav,virtual) 

c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      virtual = virtual/(st_alpha/(2d0*pi))

      endif

      return 
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compvirt_z_ew(pin,bflav,virtual)
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
      include 'cvecbos.h'
c
      integer nlegs,nf
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs)
      real*8 virtual
c
c vbfnlo stuff:
      include 'global.inc'
      integer bos,nlo
      real*8 p(0:3,np), v(0:3,nv)
      real*8 pbar(0:3,4+nv), polcol
      real*8 res,tri,box,pent,resv

      real*8 N ! color factors
      parameter(N=3d0)

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare local variables
c
      integer i,j,mu,nu
      integer FSIGN(4+nv),physToDiag(4)
      
      integer ftype(1:6)
      integer icc

      integer k

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      bos = 2
      nlo = 1 ! NLO
      polcol =1d0/(4d0*N**2)

      virtual = 0d0
      
      ftype(1:6) = 1
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         do i = 1,2
            v(mu,i) = pin(mu,i+2)
         enddo !i     
         p(mu,3) = pin(mu,5)
         p(mu,4) = pin(mu,6) 
         p(mu,5) = 0d0   
      enddo ! mu

c lepton momenta:
      do mu = 0,3             ! kinematics for 4-lepton final state
         pbar(mu,5) = v(mu,1) ! l+
         pbar(mu,6) = v(mu,2) ! l-
      enddo
      fsign(5) = -1
      fsign(6) = 1
      
      if (bflav(1).gt.0.and.bflav(2).gt.0) then

C*******************  q1 q3 ---> q2 q4 Z  **********************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1

      elseif (bflav(1).gt.0.and.bflav(2).lt.0) then

C******************* q1 qb4 ---> q2 qb3 Z  **********************
c      
      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) = -1      
      

      elseif (bflav(1).lt.0.and.bflav(2).gt.0) then

C******************* qbar2 q3 ---> qbar1 q4 Z  **********************
      
      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) =  1
      fsign(4) =  1


      elseif (bflav(1).lt.0.and.bflav(2).lt.0) then

C*******************  qbar2 qb4 ---> qbar1 qb3 Z  **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
c
      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) = -1      

      elseif (bflav(1).lt.0.and.bflav(2).lt.0) then
         
         write(*,*) 'wrong value of bflav(1) and bflav(2)'
         write(*,*) 'bflav(1) = ',bflav(1)
         write(*,*) 'bflav(2) = ',bflav(2)

      endif
            
C*****************  end of process evaluation  **********************

c get the amplitude squared:
      do mu = 0,3
         do i = 1,4
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
      enddo

c use value of bflav to select appropriate uucc, uuss etc.
      if(mod(abs(bflav(1)),2).eq.0) ftype(1) = 2 ! fermion1 = up-type 
      if(mod(abs(bflav(2)),2).eq.0) ftype(2) = 2 ! fermion2 = up-type 
      if(mod(abs(bflav(5)),2).eq.0) ftype(5) = 2 ! fermion7 = up-type 
      if(mod(abs(bflav(6)),2).eq.0) ftype(6) = 2 ! fermion8 = up-type 

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = -2*(bflav(1)/abs(bflav(1)))+3
         k = 7-ftype(icc)
      endif !nc/cc  

c all NLO pieces in one go:
c
      call qqzqq(pbar,fsign,1,1,k,resv)
	 
c triangles:
c      call qqzqq(pbar,fsign,4,1,k,tri)
c boxes:
c      call qqzqq(pbar,fsign,-4,1,k,box)
c
c      res = tri+box

      res = resv
c
      virtual = res*polcol  
c symmetry factor for leptons:
      virtual = virtual*wsymfact

      return
      end
