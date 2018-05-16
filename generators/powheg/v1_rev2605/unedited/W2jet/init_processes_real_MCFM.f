      subroutine init_processes_real_MCFM
      implicit none
 
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "nwz.f"
      include "montecarlorpp.f"
 
      double precision powheginput
      integer i1,i2,i5,i6,i7,awx,ii5,ii6,ii7,lwx,j,
     & i5in,i6in,i7in,nlf,kmad(186),kswap(186),temp(7,186),n,nproc,
     & fno(-5:5),bno(-5:5),cno(-5:5),sno(-5:5),dno(-5:5),
     & uno(-5:5),Qx3(-5:5),qsum
      data fno/-1,-1,-1,-1,-1,+0,+1,+1,+1,+1,+1/
      data bno/-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1/
      data cno/ 0,-1, 0, 0, 0, 0, 0, 0, 0,+1, 0/
      data sno/ 0, 0,-1, 0, 0, 0, 0, 0,+1, 0, 0/
      data uno/ 0, 0, 0,-1, 0, 0, 0,+1, 0, 0, 0/
      data dno/ 0, 0, 0, 0,-1, 0,+1, 0, 0, 0, 0/
      data Qx3  /1,-2,1,-2,1,0,-1,2,-1,2,-1/
      logical fermionno,bquarkno,chargecon,ucon,dcon,ccon,scon,udb,csb,
     & swap56(186),diagonal
      integer vdecaymodeW
      
      vdecaymodeW = powheginput('vdecaymodeW')
      if (abs(vdecaymodeW) .ne. 11 .and. abs(vdecaymodeW) .ne. 13) then 
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
c         open (unit=66,file='MCFMrealWp.dat',status='unknown')
      elseif (nwz .eq. -1) then
         if (vdecaymodeW .eq. 11) then 
            awx=nea_pdg
            lwx=el_pdg	 
         elseif (vdecaymodeW .eq. 13) then 
            awx=13
            lwx=-14 
         endif
c         open (unit=66,file='MCFMrealWm.dat',status='unknown')
      else
         write(6,*) 'Error in init_processes.f'
         write(6,*) 'Set jwnz to -1 or +1 in input, current value=',nwz
         stop
      endif


      
      flst_nreal=0
      nlf=5

      do i1=-nlf,nlf
      do i2=-nlf,nlf
      do i5=-nlf,nlf
      do i6=i5,nlf
      do i7=i6,nlf

      ucon=(UNO(i1)+UNO(i2)-UNO(i5)-UNO(i6)-UNO(i7)-nwz .eq. 0)
      dcon=(dno(i1)+dno(i2)-dno(i5)-dno(i6)-dno(i7)+nwz .eq. 0)
      ccon=(CNO(i1)+CNO(i2)-CNO(i5)-CNO(i6)-CNO(i7)-nwz .eq. 0)
      scon=(sNO(i1)+sNO(i2)-sNO(i5)-sNO(i6)-sNO(i7)+nwz .eq. 0)
      udb=(ucon .and. dcon)
      csb=(ccon .and. scon)
      qsum=Qx3(i1)+Qx3(i2)-Qx3(i5)-Qx3(i6)-Qx3(i7)-3*nwz
      chargecon=(qsum .eq. 0)
      fermionno=(fno(i1)+fno(i2)-fno(i5)-fno(i6)-fno(i7) .eq. 0)
      bquarkno=(bno(i1)+bno(i2)-bno(i5)-bno(i6)-bno(i7) .eq. 0)
      if ((fermionno) .and. (bquarkno) .and. (chargecon)
     & ) then
      if (((udb .eqv. .true.) .and. (csb .eqv. .false.))
     &.or.((csb .eqv. .true.) .and. (udb .eqv. .false.))) then


         i5in=i5
         i6in=i6
         i7in=i7
         call reorderreal(nwz,i1,i2,i5in,i6in,i7in,ii5,ii6,ii7)
         flst_nreal=flst_nreal+1
         temp(1,flst_nreal) = i1
         temp(2,flst_nreal) = i2
         temp(3,flst_nreal) = lwx
         temp(4,flst_nreal) = awx
         temp(5,flst_nreal) = ii5
         temp(6,flst_nreal) = ii6
         temp(7,flst_nreal) = ii7
C         write(66,55) 
C     &               flst_nreal,i1,i2,lwx,awx,ii5,ii6,ii7


      endif
      endif

      enddo
      enddo
      enddo
      enddo
      enddo
c 55   format(8(1x,i3))
c      close(unit=66)


C----reorder for comparison with Madgraph
c      if (nwz == 1) then
c      open(unit=88,file='swaprealWp.dat',status='unknown')
c      elseif (nwz == -1) then
c      open(unit=88,file='swaprealWm.dat',status='unknown')
c      endif

c      do j=1,186
c      write(6,*) j
c      read(88,*) kmad(j),kswap(j),swap56(j)
c      write(6,*) kmad(j),kswap(j),swap56(j)
c      enddo
      do n=1,7
      do nproc=1,186
      flst_real(n,nproc)=temp(n,nproc)
c      flst_real(n,kmad(nproc))=temp(n,nproc)
      enddo
      enddo


      do nproc=1,186      
      write(6,*) 'flst_nreal:',nproc,
     & flst_real(1,nproc),flst_real(2,nproc),flst_real(3,nproc),
     & flst_real(4,nproc),flst_real(5,nproc),flst_real(6,nproc),
     & flst_real(7,nproc)

      enddo

      return
      end
 
