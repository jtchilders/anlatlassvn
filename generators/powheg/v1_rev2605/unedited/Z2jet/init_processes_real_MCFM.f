      subroutine init_processes_real_MCFM
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "nwz.f"
      include "montecarlorpp.f"
      double precision powheginput
      external powheginput 
 
      integer i1,i2,i5,i6,i7,awx,ii5,ii6,ii7,lwx,
     & i5in,i6in,i7in,nlf,vdecaymodeZ

      integer fno(-5:5),bno(-5:5),cno(-5:5),sno(-5:5),dno(-5:5),
     & uno(-5:5),Qx3(-5:5),qsum
      data fno/-1,-1,-1,-1,-1,+0,+1,+1,+1,+1,+1/
      data bno/-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1/
      data cno/ 0,-1, 0, 0, 0, 0, 0, 0, 0,+1, 0/
      data sno/ 0, 0,-1, 0, 0, 0, 0, 0,+1, 0, 0/
      data uno/ 0, 0, 0,-1, 0, 0, 0,+1, 0, 0, 0/
      data dno/ 0, 0, 0, 0,-1, 0,+1, 0, 0, 0, 0/
      data Qx3  /1,-2,1,-2,1,0,-1,2,-1,2,-1/
      logical fermionno,bcon,chargecon,ucon,dcon,ccon,scon,allglue

      
      nwz=0
      vdecaymodeZ = powheginput('vdecaymodeZ') 
      if (vdecaymodeZ .ge. 11 .and. vdecaymodeZ .le. 16) then 
         lwx=vdecaymodeZ	 
         awx=-vdecaymodeZ
      else
         write(*,*) ' vdecaymodeZ=',vdecaymodeZ,' not supported'
         call pwhg_exit(-1)
      endif    


      flst_nreal=0
      nlf=5

      do i1=-nlf,nlf
      do i2=-nlf,nlf
      do i5=-nlf,nlf
      do i6=i5,nlf
      do i7=i6,nlf

      allglue=(i1==0).and.(i2==0).and.(i5==0).and.(i6==0).and.(i7==0)
      ucon=(UNO(i1)+UNO(i2)-UNO(i5)-UNO(i6)-UNO(i7) .eq. 0)
      dcon=(dno(i1)+dno(i2)-dno(i5)-dno(i6)-dno(i7) .eq. 0)
      ccon=(CNO(i1)+CNO(i2)-CNO(i5)-CNO(i6)-CNO(i7) .eq. 0)
      scon=(sNO(i1)+sNO(i2)-sNO(i5)-sNO(i6)-sNO(i7) .eq. 0)
      bcon=(bno(i1)+bno(i2)-bno(i5)-bno(i6)-bno(i7) .eq. 0)
      qsum=Qx3(i1)+Qx3(i2)-Qx3(i5)-Qx3(i6)-Qx3(i7)

      chargecon=(qsum .eq. 0)
      fermionno=(fno(i1)+fno(i2)-fno(i5)-fno(i6)-fno(i7) .eq. 0)

      if ((fermionno) .and. (chargecon).and. (.not.allglue)
     & .and.dcon.and.ucon.and.scon.and.ccon.and.bcon
     & ) then
         i5in=i5
         i6in=i6
         i7in=i7
         flst_nreal=flst_nreal+1
         call reorderZ0real(i1,i2,i5in,i6in,i7in,ii5,ii6,ii7)
         flst_real(1,flst_nreal) = i1
         flst_real(2,flst_nreal) = i2
         flst_real(3,flst_nreal) = lwx
         flst_real(4,flst_nreal) = awx
         flst_real(5,flst_nreal) = ii5
         flst_real(6,flst_nreal) = ii6
         flst_real(7,flst_nreal) = ii7
      endif

      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
 
