c
      subroutine setborn(p,bflav,born,bornjk,bmunu)

      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),bbmunu(0:3,0:3),born,colcf
      integer bflav(nlegs)
      integer j,k,mu,nu

c----------------------------------------------------
c
c numbering of momenta is q(1) q(2) -> l+(3)l-(4)q(5)q(6)
c
      call compute_tensors_z(p) 
      call compborn_z_ew(p,bflav,born,bbmunu) 

      do j=1,nlegs
         if(abs(bflav(j)).le.6) then
             if(bflav(j).eq.0) then
               do mu=0,3
                  do nu=0,3
                     bmunu(mu,nu,j)=bbmunu(mu,nu)
                  enddo
               enddo
             endif                 
            do k=j+1,nlegs
               if (((j.eq.1).and.(k.eq.5)).or.
     #             ((j.eq.2).and.(k.eq.6))) then   
                  colcf = cf    
               else
                  colcf = 0
               endif
               bornjk(j,k)=born*colcf
               bornjk(k,j)=bornjk(j,k)

            enddo !k
         endif !abs(bflav)
      enddo !j

      end
c
c==================================================
c
      subroutine compborn_z_ew(pin,bflav,born,bmunu)
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
      real*8 bmunu(0:3,0:3),born
c
c vbfnlo stuff:
      include 'global.inc'
      integer bos,nlo
      real*8 p(0:3,np), v(0:3,nv)
      real*8 pbar(0:3,4+nv), polcol
      real*8 res

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
      nlo = 0 ! LO
      polcol =1d0/(4d0*N**2)

      ftype(1:6) = 1
      
      born = 0d0
      bmunu(0:3,0:3) = 0d0

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
      do mu = 0,3             ! kinematics for 2-lepton final state
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

      else
         
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
      if(mod(abs(bflav(5)),2).eq.0) ftype(5) = 2 ! fermion5 = up-type 
      if(mod(abs(bflav(6)),2).eq.0) ftype(6) = 2 ! fermion6 = up-type 

      if ((ftype(1).eq.ftype(5))) then !nc
         k = -2*ftype(1)-ftype(2)+7
      else !cc 
         icc = -2*(bflav(1)/abs(bflav(1)))+3
         k = 7-ftype(icc)
      endif !nc/cc  

      call qqzqq(pbar,fsign,0,1,k,res)

      born = res*polcol
c symmetry factor for leptons:
      born = born*wsymfact  

      do mu=0,3
         do nu=0,3
            bmunu(mu,nu)=0d0
         enddo
      enddo

      return
      end
c
c==================================================
c
      subroutine borncolour_lh
      implicit none
c Sets up the colour for the given flavour configuration
c already filled in the Les Houches interface.
c In case there are several colour structures, one
c should pick one with a probability proportional to
c the value of the corresponding cross section, for the
c kinematics defined in the Les Houches interface

      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer i

C     -- neutral particles
      icolup(1,3)=0
      icolup(2,3)=0
      icolup(1,4)=0
      icolup(2,4)=0

c     -- colored particles
      icolup(1,1)=0
      icolup(2,1)=0
      icolup(1,2)=0
      icolup(2,2)=0
      icolup(1,5)=0
      icolup(2,5)=0
      icolup(1,6)=0
      icolup(2,6)=0

      if(idup(1).gt.0) then
         icolup(1,1)=501
         icolup(2,1)=0
      else
         icolup(1,1)=0
         icolup(2,1)=501
      endif
      
      if(idup(2).gt.0) then
         icolup(1,2)=502
         icolup(2,2)=0
      else
         icolup(1,2)=0
         icolup(2,2)=502
      endif
      
      do i=1,2
         icolup(i,5)=icolup(i,1)
         icolup(i,6)=icolup(i,2)
      enddo
      end
c
c==================================================
c
      subroutine finalize_lh
c     Set up the resonances whose mass must be preserved
c     on the Les Houches interface. 

      implicit none 

      include 'pwhg_physpar.h'   
      include 'cvecbos.h'
      logical ini
      data ini/.true./
      save ini
c     lepton masses
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass

c Z id:
      idvecbos=+23
      call add_resonance(idvecbos,3,4)
      if(ini) then
         write(*,*) 'Adding resoncance'
      endif

c  now add finite masses:
         call lhefinitemasses
         if(ini) then
            write(*,*) 'Adding finite electron/muon masses'
            ini=.false.
         endif
      
      end
