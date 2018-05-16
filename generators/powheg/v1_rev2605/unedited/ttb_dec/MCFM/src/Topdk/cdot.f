      double complex function cdot(zk1,zk2)
      implicit none
      include 'spinorsw.f'
      double complex zk1(4),zk2(4)

c--- two variants here, to handle different 4-vector conventions      
      if (spinorsw .eq. 'KM') then
        cdot=zk1(1)*zk2(1)-zk1(2)*zk2(2)-zk1(3)*zk2(3)-zk1(4)*zk2(4)
      else
        write(6,*) 'Error: spinorsw not set correctly.'
	stop
c        cdot=zk1(4)*zk2(4)-zk1(1)*zk2(1)-zk1(2)*zk2(2)-zk1(3)*zk2(3)
      endif
      
      return
      end
