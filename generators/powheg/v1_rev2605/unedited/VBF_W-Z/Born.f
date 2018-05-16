      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'process.inc'
      include 'vbfnlo-files/global.inc'  
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),bbmunu(0:3,0:3),born,colcf
      integer bflav(nlegs)
      integer k,j, mu, nu

        select case (procID)
        case (Zjj_l)
           call compborn_ewZ(p,bflav,born,0) 
        case(wpjj,wmjj)
           call compborn_ewWp(p,bflav,born,0)         
        end select

      do j=1,nlegs
         if(abs(bflav(j)).le.6) then                
            do k=j+1,nlegs
               if (((j.eq.1).and.(k.eq.5)).or.
     #                 ((j.eq.2).and.(k.eq.6))) then   
               bornjk(j,k)=born*cf
               bornjk(k,j)=born*cf                 
               else
               bornjk(j,k)=0d0
               bornjk(k,j)=0d0
               endif

            enddo
         endif
      enddo

      do j=1,nlegs
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,j)=0d0 ! no gluons-spin-correlated=0
            enddo
         enddo
      enddo
      
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compborn_ewZ(pin,bflav,born,with_NLO)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'

      include 'pwhg_st.h'
      include 'vbfnlo-files/global.inc'
      integer jsig, jsig1, jsig3, jsig5,flavcombi, with_NLO
      common /chelsum/ jsig,jsig1,jsig3,jsig5

      double precision random
      external random

      logical ini
      data ini/.true./
      save ini
c
      integer nlegs
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs),bflavtwofl(nlegs)
      real * 8 born,uucc,uuss,ddcc,ddss,udsc,ducs
c
c vbfnlo stuff:

      integer bos,lflavr(5:6) !decay to leptons rahter than neutrinos, call qqbqqi
      real*8 p(0:3,4), v(0:3,nv)
      real*8 pbar(0:3,4+nv), polcol
      real*8 res(2),resb(2)

      real*8 N ! color factors
      parameter(N=3d0)

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare local variables
c
      integer i,mu
      integer FSIGN(4+nv),physToDiag(4) 

      bos=2
C Select helicity combination by setting jsig externally. This program will
c then fill jsig1,jsig3,jsig5 to correspond to the chiralities of the 1-2,
c 3-4 and 5-6 fermion lines for external use, e.g. in calculating polarized
c cross sections. This is option 2 and it requires that jsig lies within
c the allowed sig range [1,sigmax] for the process. If jsig is outside this
c range the sum, NOT average!, over all polarization states will be calculated
c 

!
      if (bos.eq.2) then
!          jsig = min(8*random()+1d0,8.01d0)
!   for above choice: W-exchange sometimes 0 -> /born in compare_vecsb
         jsig = 0
      else
         jsig = 0
      endif


      polcol = 1d0/(4d0*N**2)   

      if( ini) then
      lflavr(5)=2
      lflavr(6)=2
      call qqbqqi(2,lflavr) ! to set up couplings for decay Z->l+l-
      ini=.false.
      endif

      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         v(mu,2) = pin(mu,4)
         v(mu,1) = pin(mu,3)            
         p(mu,3) = pin(mu,5)
         p(mu,4) = pin(mu,6) 
         do i=1,6
            pbar(mu,i)=0d0
         enddo
      enddo ! mu
      uucc=0d0
      uuss=0d0
      ddcc=0d0
      ddss=0d0
      udsc=0d0
      ducs=0d0
      if (with_NLO.eq.1) then
         call BCD_fill(p(0,1),p(0,2),                       !beam
     1                 p(0,3),p(0,4),                       !jet
     2                 v(0,1),v(0,2))                       !decay   
      endif

c
C*******************  q1 q3 ---> q2 q4  V  **********************
      if (bflav(1).gt.0.and.bflav(2).gt.0) then          
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
C NOTE: for call of wbf_zh it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      
      elseif (bflav(1).gt.0.and.bflav(2).lt.0) then 
C*******************  q1 qb4 ---> q2 qb3 V  **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = -1
      fsign(4) = -1

      elseif (bflav(1).lt.0.and.bflav(2).gt.0) then 

C*******************  qbar2 q3 ---> qbar1 q4 V   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      

c      if (ldebug) call resprint2(nmin,nmax,res)
      elseif (bflav(1).lt.0.and.bflav(2).lt.0) then 

C*******************  qbar2 qb4 ---> qbar1 qb3 V  ******************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3

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
         pbar(mu,5) = v(mu,1) ! ebar
         pbar(mu,6) = v(mu,2) ! e
      enddo

         fsign(5) = -1
         fsign(6) = 1

!neutral currents

      do i=1, nlegs
      if( abs(bflav(i)).gt.2) then
        bflavtwofl(i)=sign(abs(bflav(i))-2,bflav(i))
      else
        bflavtwofl(i)= bflav(i)
      endif
      enddo
      
      
      if((abs(bflavtwofl(1))-abs(bflavtwofl(5)).eq.0).and.(abs(bflavtwofl(2))-abs(bflavtwofl(6)).eq.0)) then
         if (abs(bflavtwofl(1)).eq.abs(bflavtwofl(2))) then
             if(abs(bflavtwofl(1)).eq.1) then
               flavcombi=4
               call qqzqq(pbar,fsign,with_NLO,flavcombi,
     1              uucc,uuss,ddcc,ddss,udsc,ducs)                 
               born=ddss*polcol

             elseif (abs(bflavtwofl(1)).eq.2) then
               flavcombi=1
      call qqzqq(pbar,fsign, with_NLO,flavcombi,
     1              uucc,uuss,ddcc,ddss,udsc,ducs)                 
               born=uucc*polcol

             endif

         elseif(abs(bflavtwofl(1)).ne.abs(bflavtwofl(2))) then !ud, du
             if(abs(bflavtwofl(1)).eq.1) then
               flavcombi=3
                     call qqzqq(pbar,fsign, with_NLO,flavcombi,
     1              uucc,uuss,ddcc,ddss,udsc,ducs)  
               born=ddcc*polcol

             elseif (abs(bflavtwofl(1)).eq.2) then
               flavcombi=2
               
                     call qqzqq(pbar,fsign, with_NLO,flavcombi,
     1              uucc,uuss,ddcc,ddss,udsc,ducs)  
               born=uuss*polcol

             endif
         endif

      elseif((abs(bflavtwofl(1))-abs(bflavtwofl(5)).ne.0).and.(abs(bflavtwofl(2))-abs(bflavtwofl(6)).ne.0)) then  !charged currents
             if((bflavtwofl(1).eq.1).or.(bflavtwofl(1).eq.-2)) then
               flavcombi=6
               
                     call qqzqq(pbar,fsign, with_NLO,flavcombi,
     1              uucc,uuss,ddcc,ddss,udsc,ducs)  
               born=ducs*polcol

             elseif ((bflavtwofl(1).eq.2).or.(bflavtwofl(1).eq.-1)) then
               flavcombi=5
                     call qqzqq(pbar,fsign, with_NLO,flavcombi,
     1              uucc,uuss,ddcc,ddss,udsc,ducs)  
               born=udsc*polcol
             endif

       else
        write(*,*) "flavour combination not found, stop"
        write(*,*) bflav
        stop

      endif
         

      return
      end

      subroutine compborn_ewWp(pin,bflav,born,with_NLO)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'vbfnlo-files/global.inc'
      include 'process.inc'
      integer jsig, jsig1, jsig3, jsig5,flavcombi,with_NLO
      common /chelsum/ jsig,jsig1,jsig3,jsig5

      double precision random
      external random

      logical ini
      data ini/.true./
      save ini
c
      integer nlegs
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs),bflavtwofl(nlegs)
      real * 8 born,uucs,ddcs,udcc,udss
c
c vbfnlo stuff:

      integer bos,lflavr(5:6) !decay to leptons rahter than neutrinos, call qqbqqi
      real*8 p(0:3,4), v(0:3,nv)
      real*8 pbar(0:3,4+nv), polcol
      real*8 res(2),resb(2)

      real*8 N ! color factors
      parameter(N=3d0)

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare local variables
c
      integer i,mu, temp
      integer FSIGN(4+nv),physToDiag(4) 
      select case (procID)
      case(wpjj)
      bos=3
      case(wmjj)
      bos=4
      end select
C Select helicity combination by setting jsig externally. This program will
c then fill jsig1,jsig3,jsig5 to correspond to the chiralities of the 1-2,
c 3-4 and 5-6 fermion lines for external use, e.g. in calculating polarized
c cross sections. This is option 2 and it requires that jsig lies within
c the allowed sig range [1,sigmax] for the process. If jsig is outside this
c range the sum, NOT average!, over all polarization states will be calculated
c 


         jsig = 0

      uucs = 0d0
      ddcs = 0d0
      udcc = 0d0
      udss = 0d0

      polcol = 1d0/(4d0*N**2)   

      if( ini) then
      if(bos.eq.3) then
         lflavr(5) = 2          
         lflavr(6) = 1 
         call qqbqqi(3,lflavr)
      elseif (bos .eq.4) then
         lflavr(5) = 1
         lflavr(6) = 2
         call qqbqqi(4,lflavr)
      endif
      ini=.false.
      endif

      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
!          do i = 1,2
            v(mu,1) = pin(mu,3)
            v(mu,2) = pin(mu,4)            
!          enddo !i    
         p(mu,3) = pin(mu,5)
         p(mu,4) = pin(mu,6) 
      enddo ! mu
      flavcombi=0
      
      if (with_NLO.eq.1) then
         call BCD_fill(p(0,1),p(0,2),                       !beam
     1                 p(0,3),p(0,4),                       !jet
     2                 v(0,1),v(0,2))                       !decay   
      endif

c
C*******************  q1 q3 ---> q2 q4  V  **********************
      if (bflav(1).gt.0.and.bflav(2).gt.0) then          
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
C NOTE: for call of wbf_zh it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      
      elseif (bflav(1).gt.0.and.bflav(2).lt.0) then 
C*******************  q1 qb4 ---> q2 qb3 V  **********************

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = -1
      fsign(4) = -1

      elseif (bflav(1).lt.0.and.bflav(2).gt.0) then 

C*******************  qbar2 q3 ---> qbar1 q4 V   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      

c      if (ldebug) call resprint2(nmin,nmax,res)
      elseif (bflav(1).lt.0.and.bflav(2).lt.0) then 

C*******************  qbar2 qb4 ---> qbar1 qb3 V  ******************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
         

c      if (ldebug) call resprint2(nmin,nmax,res)
 
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
            pbar(mu,5) = v(mu,1)   ! lbar   for decay
            pbar(mu,6) = v(mu,2)   ! l
      enddo

         fsign(5) = -1
         fsign(6) = 1





      do i=1, nlegs
      if( abs(bflav(i)).gt.2) then
        bflavtwofl(i)=sign(abs(bflav(i))-2,bflav(i))
      else
        bflavtwofl(i)= bflav(i)
      endif

      enddo

      if(bflavtwofl(1).lt.0) then
         temp=bflavtwofl(1)
         bflavtwofl(1)=-bflavtwofl(5)
         bflavtwofl(5)=-temp
      endif

      if(bflavtwofl(2).lt.0) then
         temp=bflavtwofl(2)
         bflavtwofl(2)=-bflavtwofl(6)
         bflavtwofl(6)=-temp
      endif
      
      select case (procID)
      
      case(Wpjj)
      
      
       if(bflavtwofl(1)-bflavtwofl(5).eq.0) then
       
         if (bflavtwofl(1).eq.bflavtwofl(2)) then
               flavcombi=1
               call qqwpqq(pbar,fsign, with_NLO,flavcombi,
     1              uucs,ddcs,udcc,udss)                 
               born=uucs*polcol
!          
          elseif(abs(bflavtwofl(1)).ne.abs(bflavtwofl(2))) then ! du
                     flavcombi=2
                     call qqwpqq(pbar,fsign, with_NLO,flavcombi,
     1              uucs,ddcs,udcc,udss)  
               born=ddcs*polcol
         endif

      elseif(bflavtwofl(1)-bflavtwofl(5).ne.0) then  
           if (bflavtwofl(1).eq.bflavtwofl(2)) then
               flavcombi=3
              call qqwpqq(pbar,fsign, with_NLO,flavcombi,
     1              uucs,ddcs,udcc,udss)  
               born=udcc*polcol
             elseif (bflavtwofl(1).ne.bflavtwofl(2)) then
               flavcombi=4
                     call qqwpqq(pbar,fsign, with_NLO,flavcombi,
     1              uucs,ddcs,udcc,udss)  
               born=udss*polcol
             endif

        else
        write(*,*) "flavour combination not found, stop"
        write(*,*) bflav
        stop

       endif
       
      case(Wmjj)
      
             if(bflavtwofl(1)-bflavtwofl(5).eq.0) then
       
         if (bflavtwofl(1).eq.bflavtwofl(2)) then
               flavcombi=2
               call qqwmqq(pbar,fsign, with_NLO,flavcombi,
     1              uucs,ddcs,udcc,udss)                 
               born=ddcs*polcol
          
          elseif(abs(bflavtwofl(1)).ne.abs(bflavtwofl(2))) then ! du
                     flavcombi=1
                     call qqwmqq(pbar,fsign, with_NLO,flavcombi,
     1              uucs,ddcs,udcc,udss)  
               born=uucs*polcol
         endif

      elseif(bflavtwofl(1)-bflavtwofl(5).ne.0) then  
           if (bflavtwofl(1).eq.bflavtwofl(2)) then
               flavcombi=4
              call qqwmqq(pbar,fsign, with_NLO,flavcombi,
     1              uucs,ddcs,udcc,udss)  
               born=udss*polcol
             elseif (bflavtwofl(1).ne.bflavtwofl(2)) then
               flavcombi=3
                     call qqwmqq(pbar,fsign, with_NLO,flavcombi,
     1              uucs,ddcs,udcc,udss)  
               born=udcc*polcol
             endif

        else
        write(*,*) "flavour combination not found, stop"
        write(*,*) bflav
        stop

       endif
       
       end select
!       

      return
      end
c
      subroutine borncolour_lh
c Sets up the colour for the given flavour configuration
c already filled in the Les Houches interface.
c In case there are several colour structure, one
c should pick one with a probability proportional to
c the value of the corresponding cross section, for the
c kinematics defined in the Les Houches interface
      implicit none
      integer i
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
c     neutral particles
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
      icolup(1,3)=0
      icolup(2,3)=0
      icolup(1,4)=0
      icolup(2,4)=0

      end

      subroutine finalize_lh
      implicit none
      include 'LesHouches.h'
      include 'process.inc'     
      include 'pwhg_math.h' 
      include 'vbfnlo-files/global.inc' 
c     Set up the resonances whose mass must be preserved
c     on the Les Houches interface.

      select case(procID)
      case(Zjj_l)
         call add_resonance(23,3,4)
      case(Wpjj)
         call add_resonance(24,3,4)
      case(Wmjj)
         call add_resonance(-24,3,4)
      end select



c     give masses to final-state light particles
      call lhefinitemasses
      end
