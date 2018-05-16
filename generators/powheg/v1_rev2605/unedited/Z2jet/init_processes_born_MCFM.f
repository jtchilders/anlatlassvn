      subroutine init_processes_born_MCFM
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "nwz.f"
      include "montecarlorpp.f"
  
      integer j,i1,i2,i5,i6,awx,ii5,ii6,lwx,nlf,
     & kswap(175),kmad(175),n,nproc,vdecaymodeZ 
      logical swap56(175)
      double precision powheginput
      external powheginput 
      integer fno(-5:5),bno(-5:5),cno(-5:5),sno(-5:5),dno(-5:5),
     & uno(-5:5),Qx3(-5:5),qsum
      data fno/-1,-1,-1,-1,-1,+0,+1,+1,+1,+1,+1/
      data bno/-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1/
      data cno/ 0,-1, 0, 0, 0, 0, 0, 0, 0,+1, 0/
      data sno/ 0, 0,-1, 0, 0, 0, 0, 0,+1, 0, 0/
      data uno/ 0, 0, 0,-1, 0, 0, 0,+1, 0, 0, 0/
      data dno/ 0, 0, 0, 0,-1, 0,+1, 0, 0, 0, 0/
      data Qx3  /1,-2,1,-2,1,0,-1,2,-1,2,-1/
      logical fermionno,bcon,chargecon,ucon,dcon,ccon,scon,udb,csb,
     & allglue

      nwz=0
      vdecaymodeZ = powheginput('vdecaymodeZ') 
      if (vdecaymodeZ .ge. 11 .and. vdecaymodeZ .le. 16) then 
         lwx=vdecaymodeZ	 
         awx=-vdecaymodeZ
      else
         write(*,*) ' vdecaymodeZ=',vdecaymodeZ,' not supported'
         call pwhg_exit(-1)
      endif    

      flst_nborn=0
      nlf=5
      do i1=-nlf,nlf
      do i2=-nlf,nlf
      do i5=-nlf,nlf
      do i6=i5,nlf

      allglue=(i1==0).and.(i2==0).and.(i5==0).and.(i6==0)
      ucon=(uno(i1)+uno(i2)-uno(i5)-uno(i6) .eq. 0)
      dcon=(dno(i1)+dno(i2)-dno(i5)-dno(i6) .eq. 0)
      ccon=(CNO(i1)+CNO(i2)-CNO(i5)-CNO(i6) .eq. 0)
      scon=(sNO(i1)+sNO(i2)-sNO(i5)-sNO(i6) .eq. 0)
      bcon=(bno(i1)+bno(i2)-bno(i5)-bno(i6) .eq. 0)
      qsum=Qx3(i1)+Qx3(i2)-Qx3(i5)-Qx3(i6)
      chargecon=(qsum .eq. 0)
      fermionno=(fno(i1)+fno(i2)-fno(i5)-fno(i6) .eq. 0)

      if ((fermionno) .and. (chargecon) .and.(.not.allglue) 
     & .and.(dcon).and.(ucon).and.(scon).and.(ccon).and.(bcon)
     & ) then
         call reorderZ0lord(i1,i2,i5,i6,ii5,ii6)
         flst_nborn=flst_nborn+1
         flst_born(1,flst_nborn) = i1
         flst_born(2,flst_nborn) = i2
         flst_born(3,flst_nborn) = lwx
         flst_born(4,flst_nborn) = awx
         flst_born(5,flst_nborn) = ii5
         flst_born(6,flst_nborn) = ii6
      endif

      enddo
      enddo
      enddo
      enddo

      end
 
 

