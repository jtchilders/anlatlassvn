c
      subroutine setreal(p,fermion_flav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'process.inc'
      include 'vbfnlo-files/global.inc'        
      integer nleg
      parameter (nleg=nlegreal)
      real * 8 p(0:3,nleg)
      integer fermion_flav(nleg)
      real * 8 amp2


      call calc_invariants
        
      select case (procID)
        case (Zjj_l)
           call compreal_ewZ(p,fermion_flav,amp2) 
        case(wpjj,wmjj)
           call compreal_ewWp(p,fermion_flav,amp2)              
      end select      

      
c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2d0*pi))
      
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine compreal_ewZ(pin,bflav,amp2)
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
      include 'pwhg_br.h'
      include 'pwhg_kn.h'     
!       include 'cvecbos.h'
c
      integer nlegs
      parameter (nlegs=nlegreal)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs), flavor,bflavborn(nlegs),bflavborntwofl(nlegs)
      real*8 amp2!, res2(2) 
      integer jsig, jsig1, jsig3, jsig5,flavcombi
      common /chelsum/ jsig,jsig1,jsig3,jsig5
c
c vbfnlo stuff:
      include 'vbfnlo-files/global.inc'
      integer bos
      real*8 p(0:3,5), v(0:3,2),qbar(0:4)
      real*8 pbar(0:3,nlegs-1)
      real*8 polcol,polcolq,polcolg
      real*8 res(2),temp

      real*8 N ! color factors
      parameter(N=3d0)
      logical idferm

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare local variables
c
      logical ini
      data ini/.true./
      save ini
      integer i,mu
      integer FSIGN(6),gsign,physToDiag(5)
      real * 8 uucc, uuss, ddcc,
     &       ddss, udsc, ducs

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  
    
      bos = 2
      idferm=.false.
      jsig=0

      uucc=0d0
      uuss=0d0 
      ddcc=0d0
      ddss=0d0 
      udsc=0d0 
      ducs=0d0

      polcol = 0d0
      polcolq = 1d0/(4d0*N**2)
      polcolg = 1d0/(4d0*N*(N**2-1))
 
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2) 
         v(mu,1) = pin(mu,3)   
         v(mu,2) = pin(mu,4)           
         p(mu,3) = pin(mu,5)
         p(mu,4) = pin(mu,6) 
         p(mu,5) = pin(mu,7)
      enddo   



C*******************  q1 q3 ---> q2 q4 g V V   **********************
      if (bflav(1).gt.0.and.bflav(2).gt.0) then            
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
C NOTE: for call of wbf_zh3j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1

      polcol=polcolq
      
      do i=1,nlegreal
        bflavborn(i)=bflav(i)
      enddo      
C*******************  q1 qb4 ---> q2 qb3 g V V   **********************
      elseif (bflav(1).gt.0.and.bflav(2).lt.0) then 
      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = -1
      fsign(4) = -1
      gsign    = 1

      polcol=polcolq


      do i=1,nlegreal
        bflavborn(i)=bflav(i)
      enddo
      
C*******************  qbar2 q3 ---> qbar1 q4 g V V   **********************
      elseif (bflav(1).lt.0.and.bflav(2).gt.0) then 

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1
      
      polcol=polcolq
      
      do i=1,nlegreal
        bflavborn(i)=bflav(i)
      enddo      
C*******************  qbar2 qb4 ---> qbar1 qb3 g V V  ******************
      elseif (bflav(1).lt.0.and.bflav(2).lt.0) then 
      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
      gsign    = 1         
      polcol=polcolq
      
      do i=1,nlegreal
        bflavborn(i)=bflav(i)
      enddo    
      
            
      elseif (bflav(1).eq.0.and.bflav(2).gt.0.and.bflav(5).gt.0) then
         
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
      
      polcol=polcolg
       bflavborn(1)=-bflav(7)      
      do i=2,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo     

      elseif (bflav(1).eq.0.and.bflav(2).gt.0.and.bflav(5).lt.0) then
         
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
      
      polcol=polcolg
       bflavborn(1)=-bflav(7)      
      do i=2,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo             


      elseif (bflav(1).gt.0.and.bflav(2).eq.0.
     &    and.bflav(5).gt.0.and.bflav(6).gt.0) then
      
C*******************  q g ---> q q qb W W   **********************
      
      physToDiag(1)=1  
      physToDiag(2)=5      
      physToDiag(3)=2             
      physToDiag(4)=4
      physToDiag(5)=3       
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol=polcolg
      
       bflavborn(1)=bflav(1)      
       bflavborn(2)=-bflav(7)    
      do i=3,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo   

      elseif (bflav(1).gt.0.and.bflav(2).eq.0.and.
     &        bflav(5).gt.0.and.bflav(6).lt.0) then

C*******************  q g ---> q qb q W W   **********************
     
      physToDiag(1)=1 
      physToDiag(2)=5
      physToDiag(3)=2             
      physToDiag(4)=3
      physToDiag(5)=4
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol=polcolg
      
       bflavborn(1)=bflav(1)      
       bflavborn(2)=-bflav(7)    
      do i=3,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo   



      elseif (bflav(1).eq.0.and.bflav(2).lt.0.and.bflav(5).gt.0) then
        
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

      polcol=polcolg

       bflavborn(1)=-bflav(7)      
      do i=2,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo            

      elseif (bflav(1).eq.0.and.bflav(2).lt.0.and.bflav(5).lt.0) then
        
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

      polcol=polcolg

       bflavborn(1)=-bflav(7)      
      do i=2,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo

      elseif (bflav(1).lt.0.and.bflav(2).eq.0.and.bflav(5).lt.0.and.bflav(6).lt.0) then
 
C*******************  qbar2 g ---> qbar1 qb3 q4 W W   **********************
c
      physToDiag(1)=2!1
      physToDiag(2)=5        
      physToDiag(3)=1!2             
      physToDiag(4)=3 
      physToDiag(5)=4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      polcol=polcolg
      
       bflavborn(1)=bflav(1)      
       bflavborn(2)=-bflav(7)    
      do i=3,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo         

      elseif (bflav(1).lt.0.and.bflav(2).eq.0.and.bflav(5).lt.0.and.bflav(6).gt.0) then
! C*******************  qbar2 g ---> qbar1 q4 qb3 W W   **********************
! c
      physToDiag(1)=2!1
      physToDiag(2)=5        
      physToDiag(3)=1!2             
      physToDiag(4)=4 
      physToDiag(5)=3

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      polcol=polcolg
      
       bflavborn(1)=bflav(1)      
       bflavborn(2)=-bflav(7)    
      do i=3,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo         

      
      else
         
         write(*,*) 'wrong value of bflav(1) and bflav(2)'
         write(*,*) 'bflav(1) = ',bflav(1)
         write(*,*) 'bflav(2) = ',bflav(2)
         write(*,*) 'bflav', bflav
         stop

      endif
            
C*****************  end of process evaluation  **********************

      do mu = 0,3
         do i = 1,5
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
         qbar(mu) = pbar(mu,5)
!          print*, pbar(mu,5), qbar(mu)
      enddo
      qbar(4)=0
        do mu = 0,3             ! kinematics for Z-->l+l- decay
            pbar(mu,5) = v(mu,1) ! ebar
            pbar(mu,6) = v(mu,2) ! e
         enddo
         fsign(5) = -1
         fsign(6) = 1



!neutral currents
      do i=1, nlegs
      if( abs(bflavborn(i)).gt.2) then
        bflavborntwofl(i)=sign(abs(bflavborn(i))-2,bflavborn(i))
      else
        bflavborntwofl(i)= bflavborn(i)
      endif
      enddo



      if((abs(bflavborntwofl(1))-abs(bflavborntwofl(5)).eq.0).and.(abs(bflavborntwofl(2))-abs(bflavborntwofl(6)).eq.0)) then
         if (abs(bflavborntwofl(1)).eq.abs(bflavborntwofl(2))) then
             if(abs(bflavborntwofl(1)).eq.1) then
               flavcombi=4
      call qqzqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucc,uuss,ddcc,ddss,udsc,ducs)           
               amp2=ddss*polcol

             elseif (abs(bflavborntwofl(1)).eq.2) then
               flavcombi=1   
      call qqzqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucc,uuss,ddcc,ddss,udsc,ducs)         
               amp2=uucc*polcol
             endif

         elseif(abs(bflavborntwofl(1)).ne.abs(bflavborntwofl(2))) then !ud, du
             if(abs(bflavborntwofl(1)).eq.1) then
               flavcombi=3
      call qqzqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucc,uuss,ddcc,ddss,udsc,ducs)   
               amp2=ddcc*polcol

             elseif (abs(bflavborntwofl(1)).eq.2) then
               flavcombi=2
      call qqzqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucc,uuss,ddcc,ddss,udsc,ducs)   
               amp2=uuss*polcol
             endif
         endif

      elseif((abs(bflavborntwofl(1))-abs(bflavborntwofl(5)).ne.0).and.(abs(bflavborntwofl(2))-abs(bflavborntwofl(6)).ne.0)) then  !charged currents
             if((bflavborntwofl(1).eq.1).or.(bflavborntwofl(1).eq.-2)) then
               flavcombi=6
      call qqzqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucc,uuss,ddcc,ddss,udsc,ducs)   
               amp2=ducs*polcol
             elseif ((bflavborntwofl(1).eq.2).or.(bflavborntwofl(1).eq.-1)) then
               flavcombi=5

      call qqzqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucc,uuss,ddcc,ddss,udsc,ducs)   
               amp2=udsc*polcol
             endif

      else
        write(*,*) "flavour combination not found, stop"
        write(*,*) bflav
        write(*,*) 'uborn', bflavborn
        write(*,*) 'uborn2flav',bflavborntwofl
        stop

      endif


      return
      end

      subroutine compreal_ewWp(pin,bflav,amp2)
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
      include 'pwhg_br.h'
      include 'pwhg_kn.h'     
c
      integer nlegs,lflavr(5:6)
      parameter (nlegs=nlegreal)
      real*8 pin(0:3,nlegs)  
      integer bflav(nlegs), flavor,bflavborn(nlegs),bflavborntwofl(nlegs)
      real*8 amp2!, res2(2) 
      integer jsig, jsig1, jsig3, jsig5,flavcombi
      common /chelsum/ jsig,jsig1,jsig3,jsig5
c
c vbfnlo stuff:
      include 'vbfnlo-files/global.inc'
      include 'process.inc'
      integer bos
      real*8 p(0:3,5), v(0:3,2),qbar(0:4)
      real*8 pbar(0:3,nlegs-1)
      real*8 polcol,polcolq,polcolg
      real*8 res(2),temp

      real*8 N ! color factors
      parameter(N=3d0)
      logical idferm

      complex*16 zero
      parameter (zero=(0d0,0d0))
c
c declare local variables
c
      logical ini
      data ini/.true./
      save ini
      integer i,mu
      integer FSIGN(6),gsign,physToDiag(5)
      real * 8 uucc(3), uuss(3), ddcc(3),
     &       ddss(3), udsc(3), ducs(3)
     &       ,uucs,ddcs,udcc,udss

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  

      select case (procID)
      case(wpjj)
      bos=3
      case(wmjj)
      bos=4
      end select
      idferm=.false.
      jsig=0

      polcol = 0d0
         polcolq = 1d0/(4d0*N**2d0)
         polcolg = 1d0/(4d0*N*(N**2d0-1d0))
 
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         v(mu,1) = pin(mu,3)
         v(mu,2) = pin(mu,4)    
         p(mu,3) = pin(mu,5)
         p(mu,4) = pin(mu,6) 
         p(mu,5) = pin(mu,7)
      enddo   

C*******************  q1 q3 ---> q2 q4 g V V   **********************
      if ((bflav(1).gt.0).and.(bflav(2).gt.0)) then            
c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon
C NOTE: for call of wbf_zh3j it is important that p(*,1,*) and p(*,3,*)
c correspond to 1-2 fermion line ALWAYS, i.e physToDiag(1/2)={1,3} and 
c similarly physToDiag(3/4)={2,4} for the 3-4 fermion line
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1

      polcol=polcolq
      
      do i=1,nlegreal
        bflavborn(i)=bflav(i)
      enddo      
C*******************  q1 qb4 ---> q2 qb3 g V V   **********************
      elseif ((bflav(1).gt.0).and.(bflav(2).lt.0)) then 
      physToDiag(1)=1    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=2    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = -1
      fsign(4) = -1
      gsign    = 1

      polcol=polcolq


      do i=1,nlegreal
        bflavborn(i)=bflav(i)
      enddo
      
C*******************  qbar2 q3 ---> qbar1 q4 g V V   **********************
      elseif ((bflav(1).lt.0).and.(bflav(2).gt.0)) then 

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=3
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=4
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1
      gsign    = 1
      
      polcol=polcolq
      
      do i=1,nlegreal
        bflavborn(i)=bflav(i)
      enddo      
C*******************  qbar2 qb4 ---> qbar1 qb3 g V V  ******************
      elseif ((bflav(1).lt.0).and.(bflav(2).lt.0)) then 
      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=4
      physToDiag(3)=1    !physToDiag(3/4) are labels of outgoing quarks.
      physToDiag(4)=3
      physToDiag(5)=5   ! gluon

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1 
      fsign(4) = -1
      gsign    = 1         
      polcol=polcolq
      
      do i=1,nlegreal
        bflavborn(i)=bflav(i)
      enddo    
      
            
      elseif ((bflav(1).eq.0).and.(bflav(2).gt.0).and.(bflav(5).gt.0)) then
         
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
      
      polcol=polcolg
       bflavborn(1)=-bflav(7)      
      do i=2,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo     

      elseif ((bflav(1).eq.0).and.(bflav(2).gt.0).and.(bflav(5).lt.0)) then
         
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
      
      polcol=polcolg
       bflavborn(1)=-bflav(7)      
      do i=2,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo             


      elseif (bflav(1).gt.0.and.bflav(2).eq.0.
     &    and.bflav(5).gt.0.and.bflav(6).gt.0) then
      
C*******************  q g ---> q q qb W W   **********************
      
      physToDiag(1)=1  
      physToDiag(2)=5      
      physToDiag(3)=2             
      physToDiag(4)=4
      physToDiag(5)=3       
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol=polcolg
      
       bflavborn(1)=bflav(1)      
       bflavborn(2)=-bflav(7)    
      do i=3,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo   

      elseif (bflav(1).gt.0.and.bflav(2).eq.0.and.
     &        bflav(5).gt.0.and.bflav(6).lt.0) then

C*******************  q g ---> q qb q W W   **********************
     
      physToDiag(1)=1 
      physToDiag(2)=5
      physToDiag(3)=2             
      physToDiag(4)=3
      physToDiag(5)=4
c
      fsign(1) =  1
      fsign(2) =  1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1

      polcol=polcolg
      
       bflavborn(1)=bflav(1)      
       bflavborn(2)=-bflav(7)    
      do i=3,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo   



      elseif (bflav(1).eq.0.and.bflav(2).lt.0.and.bflav(5).gt.0) then
        
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

      polcol=polcolg

       bflavborn(1)=-bflav(7)      
      do i=2,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo            

      elseif (bflav(1).eq.0.and.bflav(2).lt.0.and.bflav(5).lt.0) then
        
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

      polcol=polcolg

       bflavborn(1)=-bflav(7)      
      do i=2,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo

      elseif (bflav(1).lt.0.and.bflav(2).eq.0.and.bflav(5).lt.0.and.bflav(6).lt.0) then
 
C*******************  qbar2 g ---> qbar1 qb3 q4 W W   **********************
c
      physToDiag(1)=2!1
      physToDiag(2)=5        
      physToDiag(3)=1!2             
      physToDiag(4)=3 
      physToDiag(5)=4

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      polcol=polcolg
      
       bflavborn(1)=bflav(1)      
       bflavborn(2)=-bflav(7)    
      do i=3,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo         

      elseif (bflav(1).lt.0.and.bflav(2).eq.0.and.bflav(5).lt.0.and.bflav(6).gt.0) then
C*******************  qbar2 g ---> qbar1 q4 qb3 W W   **********************
c
      physToDiag(1)=2!1
      physToDiag(2)=5        
      physToDiag(3)=1!2             
      physToDiag(4)=4 
      physToDiag(5)=3

      fsign(1) = -1
      fsign(2) = -1
      fsign(3) = -1
      fsign(4) =  1
      gsign    = -1
      polcol=polcolg
      
       bflavborn(1)=bflav(1)      
       bflavborn(2)=-bflav(7)    
      do i=3,nlegreal-1
        bflavborn(i)=bflav(i)
      enddo         

      
      else
         
         write(*,*) 'wrong value of bflav(1) and bflav(2)'
         write(*,*) 'bflav(1) = ',bflav(1)
         write(*,*) 'bflav(2) = ',bflav(2)
         write(*,*) 'bflav', bflav
         stop

      endif
            
C*****************  end of process evaluation  **********************

      do mu = 0,3
         do i = 1,5
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
         qbar(mu) = pbar(mu,5)
      enddo
      qbar(4)=0
        do mu = 0,3             ! kinematics for Z-->l+l- decay
            pbar(mu,5) = v(mu,1) ! ebar
            pbar(mu,6) = v(mu,2) ! e
         enddo
         fsign(5) = -1
         fsign(6) = 1


!neutral currents
      do i=1, nlegs
      if( abs(bflavborn(i)).gt.2) then
        bflavborntwofl(i)=sign(abs(bflavborn(i))-2,bflavborn(i))
      else
        bflavborntwofl(i)= bflavborn(i)
      endif
      enddo



      if(bflavborntwofl(1).lt.0) then
         temp=bflavborntwofl(1)
         bflavborntwofl(1)=-bflavborntwofl(5)
         bflavborntwofl(5)=-temp
      endif

      if(bflavborntwofl(2).lt.0) then
         temp=bflavborntwofl(2)
         bflavborntwofl(2)=-bflavborntwofl(6)
         bflavborntwofl(6)=-temp
      endif
      
      select case (procID)
      
      case(Wpjj)
      
      if(bflavborntwofl(1)-bflavborntwofl(5).eq.0) then
       
         if (bflavborntwofl(1).eq.bflavborntwofl(2)) then
               flavcombi=1

      call qqWPqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucs,ddcs,udcc,udss)                 
               amp2=uucs*polcol
          
          elseif(abs(bflavborntwofl(1)).ne.abs(bflavborntwofl(2))) then ! du
                 flavcombi=2
      call qqWPqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucs,ddcs,udcc,udss)     
               amp2=ddcs*polcol
         endif

      elseif(bflavborntwofl(1)-bflavborntwofl(5).ne.0) then  
           if (bflavborntwofl(1).eq.bflavborntwofl(2)) then
               flavcombi=3
      call qqWPqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucs,ddcs,udcc,udss)     
               amp2=udcc*polcol
             elseif (bflavborntwofl(1).ne.bflavborntwofl(2)) then
               flavcombi=4
      call qqWPqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucs,ddcs,udcc,udss)     
               amp2=udss*polcol
             endif

      else
        write(*,*) "flavour combination not found, stop"
        write(*,*) bflav
        write(*,*) 'uborn', bflavborn
        write(*,*) 'uborn2flav',bflavborntwofl
        stop

      endif
      case(Wmjj)
       if(bflavborntwofl(1)-bflavborntwofl(5).eq.0) then
       
         if (bflavborntwofl(1).eq.bflavborntwofl(2)) then
               flavcombi=2
      call qqWMqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucs,ddcs,udcc,udss)                 
               amp2=ddcs*polcol

          elseif(abs(bflavborntwofl(1)).ne.abs(bflavborntwofl(2))) then ! du

                 flavcombi=1
      call qqWMqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucs,ddcs,udcc,udss)     
               amp2=uucs*polcol
         endif

      elseif(bflavborntwofl(1)-bflavborntwofl(5).ne.0) then  
           if (bflavborntwofl(1).eq.bflavborntwofl(2)) then
               flavcombi=4
      call qqWMqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucs,ddcs,udcc,udss)     
               amp2=udss*polcol
             elseif (bflavborntwofl(1).ne.bflavborntwofl(2)) then
               flavcombi=3
      call qqWMqqj(pbar,fsign,qbar,gsign,flavcombi,
     1            uucs,ddcs,udcc,udss)     
               amp2=udcc*polcol
             endif

      else
        write(*,*) "flavour combination not found, stop"
        write(*,*) bflav
        write(*,*) 'uborn', bflavborn
        write(*,*) 'uborn2flav',bflavborntwofl
        stop

      endif      

      end select


      return
      end
