      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 masswindow
      logical verbose
      parameter(verbose=.true.)
      integer i,j
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
C     EW input as in WpWp paper 
      ph_Zmass   = 91.1876d0     
      ph_Zwidth  =  2.4952d0
      ph_Wmass   = 80.419d0     
      ph_Wwidth  =  2.141d0
      ph_alphaem = 1d0/128.802d0
      ph_sthw2   = 0.2222d0 

c     DIAGONAL CKM 
      ph_CKM(1,1)=1d0 
      ph_CKM(1,2)=0d0 
      ph_CKM(1,3)=0d0
      ph_CKM(2,1)=0d0 
      ph_CKM(2,2)=1d0 
      ph_CKM(2,3)=0d0
      ph_CKM(3,1)=0d0
      ph_CKM(3,2)=0d0
      ph_CKM(3,3)=1d0

c     number of light flavors
      st_nlight = 5


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      ph_Zmass2 = ph_Zmass**2
      ph_Wmass2 = ph_Wmass**2

c     set mass windows around Z-mass peak in unit of ph_Zwidth
c     It is used in the generation of the Born phase space
CAVEAT : masswindow should be a parameter passed by the user
      masswindow = 30
      ph_Zmass2low=(ph_Zmass-masswindow*ph_Zwidth)**2
      ph_Zmass2high=(ph_Zmass+masswindow*ph_Zwidth)**2
      ph_ZmZw = ph_Zmass * ph_Zwidth

c     set mass window around W-mass peak in unit of ph_Wwidth
c     It is used in the generation of the Born phase space
CAVEAT : masswindow should be a parameter passed by the user
C      masswindow = 30
      masswindow = 0 ! GZ 
      ph_Wmass2low=(ph_Wmass-masswindow*ph_Wwidth)**2
      ph_Wmass2high=(ph_Wmass+masswindow*ph_Wwidth)**2
      ph_WmWw = ph_Wmass * ph_Wwidth

      ph_unit_e = sqrt(4*pi*ph_alphaem)

      if(verbose) then
      write(*,*) '*************************************'
      write(*,*) 'Z mass = ',ph_Zmass
      write(*,*) 'Z width = ',ph_Zwidth
      write(*,*) 'W mass = ',ph_Wmass
      write(*,*) 'W width = ',ph_Wwidth
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) '(unit_e)^2 = ',ph_unit_e**2   
      write(*,*) '(g_w)^2 = ',ph_unit_e*ph_unit_e/ph_sthw2   
      write(*,*) 'CKM matrix' 
      do i=1,3
         write(*,*) (ph_CKM(i,j),j=1,3)
      enddo
      write(*,*) '*************************************'
      write(*,*)
      write(*,*) '*************************************'
      write(*,*) sqrt(ph_Wmass2low),'< M_W <',sqrt(ph_Wmass2high)
      write(*,*) '*************************************'
      endif
      end

