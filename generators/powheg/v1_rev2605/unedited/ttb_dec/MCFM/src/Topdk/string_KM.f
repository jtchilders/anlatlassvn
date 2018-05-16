c-----functions for fermion strings

       double complex function string0_KM(zp1,zp2)
       implicit none
       double complex zp1(4),zp2(4),res
       call psp(zp2,zp1,res)
       string0_KM = res
       return 
       end


       double complex function string1_KM(zp1,zp2,zp3)
       implicit none
       double complex zp1(4),zp2(4),zp3(4),sp2a(4), res
       call spb(zp1,zp2,sp2a)
       call psp(sp2a,zp3,res)
       string1_KM = res
       return 
       end

      double complex function string2_KM(zp1,zp2,zp3,zp4)
      implicit none
      double complex zp1(4),zp2(4),zp3(4),zp4(4),sp2a(4),sp2b(4),res

      call spb(zp1,zp2,sp2a)
      call spb(sp2a,zp3,sp2b)
      call psp(sp2b,zp4,res)

      string2_KM = res

      return 
      end

      double complex function string3_KM(zp1,zp2,zp3,zp4,zp5)
      implicit none
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4)
      double complex sp2a(4),sp2b(4),res

      call spb(zp1,zp2,sp2a)
      call spb(sp2a,zp3,sp2b)
      call spb(sp2b,zp4,sp2a)
      call psp(sp2a,zp5,res)

      string3_KM = res

      return 
      end

      double complex function string4_KM(zp1,zp2,zp3,zp4,zp5,zp6)
      implicit none
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4)
      double complex sp2a(4),sp2b(4),res

      call spb(zp1,zp2,sp2a)
      call spb(sp2a,zp3,sp2b)
      call spb(sp2b,zp4,sp2a)
      call spb(sp2a,zp5,sp2b)
      call psp(sp2b,zp6,res)

      string4_KM = res

      return 
      end


      double complex function string5_KM(zp1,zp2,zp3,zp4,zp5,zp6,zp7)
      implicit none
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4),zp7(4)
      double complex sp2a(4),sp2b(4),res
      call spb(zp1,zp2,sp2a)
      call spb(sp2a,zp3,sp2b)
      call spb(sp2b,zp4,sp2a)
      call spb(sp2a,zp5,sp2b)
      call spb(sp2b,zp6,sp2a)
      call psp(sp2a,zp7,res)
      string5_KM = res
      return 
      end



      double complex function string6_KM(zp1,zp2,zp3,zp4,zp5,zp6,zp7,
     & zp8)
      implicit none
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4)
      double complex zp7(4),zp8(4)
      double complex sp2a(4),sp2b(4),res
      call spb(zp1,zp2,sp2a)
      call spb(sp2a,zp3,sp2b)
      call spb(sp2b,zp4,sp2a)
      call spb(sp2a,zp5,sp2b)
      call spb(sp2b,zp6,sp2a)
      call spb(sp2a,zp7,sp2b)
      call psp(sp2b,zp8,res)
      string6_KM = res
      return 
      end


      double complex function string7_KM(zp1,zp2,zp3,zp4,zp5,zp6,zp7,
     & zp8,zp9)
      implicit none
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4)
      double complex zp7(4),zp8(4),zp9(4)
      double complex sp2a(4),sp2b(4),res

      call spb(zp1,zp2,sp2a)
      call spb(sp2a,zp3,sp2b)
      call spb(sp2b,zp4,sp2a)
      call spb(sp2a,zp5,sp2b)
      call spb(sp2b,zp6,sp2a)
      call spb(sp2a,zp7,sp2b)
      call spb(sp2b,zp8,sp2a)
      call psp(sp2a,zp9,res)
      string7_KM = res

      return 
      end



      double complex function string8_KM(zp1,zp2,zp3,zp4,zp5,zp6,zp7,
     & zp8,zp9,zp10)
      implicit none
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4)
      double complex zp7(4),zp8(4),zp9(4),zp10(4)
      double complex sp2a(4),sp2b(4),res

      call spb(zp1,zp2,sp2a)
      call spb(sp2a,zp3,sp2b)
      call spb(sp2b,zp4,sp2a)
      call spb(sp2a,zp5,sp2b)
      call spb(sp2b,zp6,sp2a)
      call spb(sp2a,zp7,sp2b)
      call spb(sp2b,zp8,sp2a)
      call spb(sp2a,zp9,sp2b)
      call psp(sp2b,zp10,res)
      string8_KM = res

      return 
      end



      double complex function string9_KM(zp1,zp2,zp3,zp4,zp5,zp6,zp7,
     & zp8,zp9,zp10,zp11)
      implicit none
      double complex zp1(4),zp2(4),zp3(4),zp4(4),zp5(4),zp6(4)
      double complex zp7(4),zp8(4),zp9(4),zp10(4), zp11(4)
      double complex sp2a(4),sp2b(4),res


      call spb(zp1,zp2,sp2a)
      call spb(sp2a,zp3,sp2b)
      call spb(sp2b,zp4,sp2a)
      call spb(sp2a,zp5,sp2b)
      call spb(sp2b,zp6,sp2a)
      call spb(sp2a,zp7,sp2b)
      call spb(sp2b,zp8,sp2a)
      call spb(sp2a,zp9,sp2b)
      call spb(sp2b,zp10,sp2a)
      call psp(sp2a,zp11,res)

      string9_KM = res

      return 
      end


