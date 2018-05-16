      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'global.inc'
      include 'process.inc'

      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs)
      integer vflav(nlegs)
      real * 8 virtual
      include 'PhysPars.h'
      real *8 powheginput
      external powheginput 

      logical firsttime
      save firsttime
      data firsttime/.true./

      integer fakevirt
      save fakevirt 


      if (firsttime) then
         fakevirt=powheginput("#fakevirt")
         if (fakevirt.eq.1) write(*,*) 'WARNING: Using fakevirt !'
         firsttime = .false.
      endif

      if(kn_jacborn.eq.0d0) then
        virtual=1d-20
      else
 
     

        if(fakevirt.eq.1) then    
          select case(procID)     
          case(Zjj_l)
             call compborn_ewZ(p,vflav,virtual,0)       
          case(Wpjj,Wmjj)
             call compborn_ewWp(p,vflav,virtual,0)               
          end select
           virtual=0.2d0*virtual
        else
          call calc_invariants 
        
          select case(procID)
        
          case(Zjj_l)
             call compborn_ewZ(p,vflav,virtual,1)       
          case(Wpjj,Wmjj)
             call compborn_ewWp(p,vflav,virtual,1)               
          end select
             
          virtual = virtual/(st_alpha/(2d0*pi))
        endif
      endif
      return

      end


