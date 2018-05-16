      subroutine init_processes_born_MCFM
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "nwz.f"
      include "montecarlorpp.f"
  
      integer j,i1,i2,i5,i6,awx,ii5,ii6,lwx,nlf,
     & kswap(114),kmad(114),temp(6,114),n,nproc
      logical swap56(114),diagonal
      double precision powheginput

      integer fno(-5:5),bno(-5:5),cno(-5:5),sno(-5:5),dno(-5:5),
     & uno(-5:5),Qx3(-5:5),qsum
      data fno/-1,-1,-1,-1,-1,+0,+1,+1,+1,+1,+1/
      data bno/-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1/
      data cno/ 0,-1, 0, 0, 0, 0, 0, 0, 0,+1, 0/
      data sno/ 0, 0,-1, 0, 0, 0, 0, 0,+1, 0, 0/
      data uno/ 0, 0, 0,-1, 0, 0, 0,+1, 0, 0, 0/
      data dno/ 0, 0, 0, 0,-1, 0,+1, 0, 0, 0, 0/
      data Qx3  /1,-2,1,-2,1,0,-1,2,-1,2,-1/
      logical fermionno,bquarkno,chargecon,ucon,dcon,ccon,scon,udb,csb
      integer vdecaymodeW

      vdecaymodeW = powheginput('vdecaymodeW')
      if (abs(vdecaymodeW) .ne. 11 .and. abs(vdecaymodeW) .ne. 13) then 
         write(*,*) 'vdecaymodeW set to', vdecaymodeW 
         stop 'vdecaymodeW invalid: W-decay only to e or mu'
      endif
      if (vdecaymodeW .gt. 0) then 
         nwz = -1 ! ( 11 -> e-,  13 -> mu-) 
      else 
         nwz = 1  ! (-11 -> e+, -13 -> mu+) 
      endif
C      nwz=int(powheginput("wtype"))
C-----This process only implemented for diagonal CKM
      diagonal=.true.
      call setupckmallowed(nwz,diagonal)

      if (nwz .eq. 1) then
         if (vdecaymodeW .eq. -11) then 
            awx=ea_pdg
            lwx=nel_pdg	 
         elseif (vdecaymodeW .eq. -13) then 
            awx=-13
            lwx=14	 
         endif
c         open (unit=66,file='MCFMlordWp.dat',status='unknown')
      elseif (nwz .eq. -1) then
         if (vdecaymodeW .eq. 11) then 
            awx=nea_pdg
            lwx=el_pdg	 
         elseif (vdecaymodeW .eq. 13) then 
            awx=13
            lwx=-14 
         endif
c         open (unit=66,file='MCFMlordWm.dat',status='unknown')
      else
         write(6,*) 'Error in init_processes.f'
         write(6,*) 'Set jwnz to -1 or +1 in input, current value=',nwz
         stop
      endif

      flst_nborn=0
      nlf=5
      do i1=-nlf,nlf
      do i2=-nlf,nlf
      do i5=-nlf,nlf
      do i6=i5,nlf

      ucon=(UNO(i1)+UNO(i2)-UNO(i5)-UNO(i6)-nwz .eq. 0)
      dcon=(dno(i1)+dno(i2)-dno(i5)-dno(i6)+nwz .eq. 0)
      ccon=(CNO(i1)+CNO(i2)-CNO(i5)-CNO(i6)-nwz .eq. 0)
      scon=(sNO(i1)+sNO(i2)-sNO(i5)-sNO(i6)+nwz .eq. 0)
      udb=(ucon .and. dcon)
      csb=(ccon .and. scon)
      qsum=Qx3(i1)+Qx3(i2)-Qx3(i5)-Qx3(i6)-3*nwz
      chargecon=(qsum .eq. 0)
      fermionno=(fno(i1)+fno(i2)-fno(i5)-fno(i6) .eq. 0)
      bquarkno=(bno(i1)+bno(i2)-bno(i5)-bno(i6) .eq. 0)

      if ((fermionno) .and. (bquarkno) .and. (chargecon)
     & ) then
      if (((udb .eqv. .true.) .and. (csb .eqv. .false.))
     &.or.((csb .eqv. .true.) .and. (udb .eqv. .false.))) then

         call reorderlord(nwz,i1,i2,i5,i6,ii5,ii6)
         flst_nborn=flst_nborn+1
         temp(1,flst_nborn) = i1
         temp(2,flst_nborn) = i2
         temp(3,flst_nborn) = lwx
         temp(4,flst_nborn) = awx
         temp(5,flst_nborn) = ii5
         temp(6,flst_nborn) = ii6
c         write(66,55) flst_nborn,i1,i2,lwx,awx,ii5,ii6 

      endif
      endif

      enddo
      enddo
      enddo
      enddo

C 55   format(7(1x,i3))
C      close(unit=66)

      do n=1,6
      do nproc=1,114      
      flst_born(n,nproc)=temp(n,nproc)
      enddo
      enddo

CC----reorder for comparison with Madgraph
C      if (nwz == 1) then
C      open(unit=88,file='swaplordWp.dat',status='unknown')
C      elseif (nwz == -1) then
C      open(unit=88,file='swaplordWm.dat',status='unknown')
C      endif
C
C      do j=1,114
C      read(88,*) kmad(j),kswap(j),swap56(j)
C      write(6,*) kmad(j),kswap(j),swap56(j)
C      enddo
C      do n=1,6
C      do nproc=1,114      
C      flst_born(n,kmad(nproc))=temp(n,nproc)
C      enddo
C      enddo
C
C
C      do nproc=1,114      
C      write(6,*) 'aftermadswap:flst_nborn:',nproc,
C     & flst_born(1,nproc),flst_born(2,nproc),flst_born(3,nproc),
C     & flst_born(4,nproc),flst_born(5,nproc),flst_born(6,nproc)
C      enddo

      return
      end
 
 

