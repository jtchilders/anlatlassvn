c--- this is a set of switchyard functions that points
c--- to two sets of functions for fermion strings:
c---  one compatible with the numerical evaluation of top decay amplitudes
c---  and the other written for the implementation of Korner-Merebashvili

      double complex function string0(zp1,zp2)
      implicit none
      include 'spinorsw.f'
      double complex zp1(4),zp2(4),string0_KM!,string0_w
       
      if (spinorsw .eq. 'KM') then
        string0 = string0_KM(zp1,zp2)
      else
        write(6,*) 'Error: spinorsw not set correctly.'
	stop
c        string0 = string0_w(zp1,zp2)
      endif
      
      return 
      end

      double complex function string1(zp1,zp2,zp3)
      implicit none
      include 'spinorsw.f'
      double complex zp1(4),zp2(4),zp3(4),string1_KM!,string1_w
       
      if (spinorsw .eq. 'KM') then
        string1 = string1_KM(zp1,zp2,zp3)
      else
        write(6,*) 'Error: spinorsw not set correctly.'
	stop
c        string1 = string1_w(zp1,zp2,zp3)
      endif
      
      return 
      end

      double complex function string2(zp1,zp2,zp3,zp4)
      implicit none
      include 'spinorsw.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4),string2_KM!,string2_w
       
      if (spinorsw .eq. 'KM') then
        string2 = string2_KM(zp1,zp2,zp3,zp4)
      else
        write(6,*) 'Error: spinorsw not set correctly.'
	stop
c        string2 = string2_w(zp1,zp2,zp3,zp4)
      endif
      
      return 
      end

      double complex function string3(zp1,zp2,zp3,zp4,zp5)
      implicit none
      include 'spinorsw.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),
     . string3_KM!,string3_w
       
      if (spinorsw .eq. 'KM') then
        string3 = string3_KM(zp1,zp2,zp3,zp4,zp5)
      else
        write(6,*) 'Error: spinorsw not set correctly.'
	stop
c        string3 = string3_w(zp1,zp2,zp3,zp4,zp5)
      endif
      
      return 
      end

      double complex function string4(zp1,zp2,zp3,zp4,zp5,zp6)
      implicit none
      include 'spinorsw.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4),
     & string4_KM!,string4_w
       
      if (spinorsw .eq. 'KM') then
        string4 = string4_KM(zp1,zp2,zp3,zp4,zp5,zp6)
      else
        write(6,*) 'Error: spinorsw not set correctly.'
	stop
c        string4 = string4_w(zp1,zp2,zp3,zp4,zp5,zp6)
      endif
      
      return 
      end

      double complex function string5(zp1,zp2,zp3,zp4,zp5,zp6,zp7)
      implicit none
      include 'spinorsw.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4),
     & zp7(4),string5_KM!,string5_w
       
      if (spinorsw .eq. 'KM') then
        string5 = string5_KM(zp1,zp2,zp3,zp4,zp5,zp6,zp7)
      else
        write(6,*) 'Error: spinorsw not set correctly.'
	stop
c        string5 = string5_w(zp1,zp2,zp3,zp4,zp5,zp6,zp7)
      endif
      
      return 
      end

      double complex function string6(zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8)
      implicit none
      include 'spinorsw.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4),
     & zp7(4),zp8(4),string6_KM!,string6_w
       
      if (spinorsw .eq. 'KM') then
        string6 = string6_KM(zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8)
      else
        write(6,*) 'Error: spinorsw not set correctly.'
	stop
c        string6 = string6_w(zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8)
      endif
      
      return 
      end

      double complex function string7(zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8,
     & zp9)
      implicit none
      include 'spinorsw.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4),
     & zp7(4),zp8(4),zp9(4),string7_KM!,string7_w
       
      if (spinorsw .eq. 'KM') then
        string7 = string7_KM(zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8,
     & zp9)
      else
        write(6,*) 'Error: spinorsw not set correctly.'
	stop
c        string7 = string7_w(zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8,
c     & zp9)
      endif
      
      return 
      end

      double complex function string8(zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8,
     & zp9,zp10)
      implicit none
      include 'spinorsw.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4),
     & zp7(4),zp8(4),zp9(4),zp10(4),string8_KM
       
      if (spinorsw .eq. 'KM') then
        string8 = string8_KM(zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8,
     & zp9,zp10)
      else
        write(6,*) 'string8 not defined if spinorsw not KM!'
	stop
      endif
      
      return 
      end

      double complex function string9(zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8,
     & zp9,zp10,zp11)
      implicit none
      include 'spinorsw.f'
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4),
     & zp7(4),zp8(4),zp9(4),zp10(4),zp11(4),string9_KM
       
      if (spinorsw .eq. 'KM') then
        string9 = string9_KM(zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8,
     & zp9,zp10,zp11)
      else
        write(6,*) 'string8 not defined if spinorsw not KM!'
	stop
      endif
      
      return 
      end



