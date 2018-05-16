      subroutine setupBmunu(myp,myin,fxn_gvec,nlegs,bflav,Bmunu)
      implicit none
C-----Author R.K.Ellis, April 2012
C     sets up the Bmunu in the form wanted by powheg using gvec routines.
C     myp(mxpart,4)           ! momenta in MCFM notation
C     myin(nlegs)             ! a redirection vector projecting 
C                             ! powheg momentum labels onto mcfm labels
C     fxn_gvec                ! the appropriate MCFM gvec type routine
C     bflav(nlegs)            ! the powheg born parton flavor descriptor
C     nlegs                   ! the powheg number of legs
C     Bmunu(0:3,0:3,nlegs)    ! the powheg Bmunu array
      include 'constants.f'
      include 'mcfmtopwhg.f'
      integer nlegs,bflav(nlegs),j,fl1,fl2,ro,si,myin(nlegs),in
      double precision bmunu(0:3,0:3,nlegs),myp(mxpart,4)
      double precision e1(4),e2(4),e3(4),
     & msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf),msq3(-nf:nf,-nf:nf)

      fl1=bflav(1)
      fl2=bflav(2)

      do j=1,nlegs
      if (bflav(j) .eq. 0) then
      in=myin(j)
      call findperp(myp,in,e1,e2,e3)
      call fxn_gvec(myp,e1,in,msq1)
      call fxn_gvec(myp,e2,in,msq2)
      call fxn_gvec(myp,e3,in,msq3)

      do ro=1,4
      do si=1,4
      Bmunu(pwhg(ro),pwhg(si),j)=
     & msq1(fl1,fl2)*e1(ro)*e1(si)+msq2(fl1,fl2)*e2(ro)*e2(si)
     & +0.5d0*(msq3(fl1,fl2)-msq1(fl1,fl2)-msq2(fl1,fl2))
     & *(e3(ro)*e3(si)-e1(ro)*e1(si)-e2(ro)*e2(si))
      enddo
      enddo

c      write(6,*) 'e1',e1
c      write(6,*) 'e2',e2
c      write(6,*) 'e3',e3
c      write(6,*) '3-1-2',msq3(fl1,fl2)-msq1(fl1,fl2)-msq2(fl1,fl2)
c      write(6,*) 'e3.x^2',
c     & +e3(4)*e3(4)*(e3(4)*e3(4)-e1(4)*e1(4)-e2(4)*e2(4))
c     & -e3(1)*e3(1)*(e3(1)*e3(1)-e1(1)*e1(1)-e2(1)*e2(1))
c     & -e3(2)*e3(2)*(e3(2)*e3(2)-e1(2)*e1(2)-e2(2)*e2(2))
c     & -e3(3)*e3(3)*(e3(3)*e3(3)-e1(3)*e1(3)-e2(3)*e2(3))
c       write(6,*) 'd',dotpp(e3,e3)**2-dotpp(e3,e1)**2-dotpp(e3,e2)**2
c      write(6,*) '1',msq1(fl1,fl2),
c     & +e1(1)*e1(1)*Bmunu(1,1,in)
c     & +e1(1)*e1(2)*Bmunu(1,2,in)
c     & +e1(2)*e1(1)*Bmunu(2,1,in)
c     & +e1(2)*e1(2)*Bmunu(2,2,in)
c      write(6,*) '2',msq2(fl1,fl2),
c     & +e2(1)*e2(1)*Bmunu(1,1,in)
c     & +e2(2)*e2(1)*Bmunu(2,1,in)
c     & +e2(1)*e2(2)*Bmunu(1,2,in)
c     & +e2(2)*e2(2)*Bmunu(2,2,in)

c      write(6,*) '3',msq3(fl1,fl2),
c     & +e3(1)*e3(1)*Bmunu(1,1,in)
c     & +e3(1)*e3(2)*Bmunu(1,2,in)
c     & +e3(2)*e3(1)*Bmunu(2,1,in)
c     & +e3(2)*e3(2)*Bmunu(2,2,in)

      endif
      enddo
      return
      end


      subroutine rotate3vec(dir,sinphi,cosphi,vec)
c----Rotates 3-vector vec counterclockwise around the direction 
c----dir (|dir|=1) with angle phi, given sin phi and cos phi. 
      implicit none
      real * 8 sinphi,cosphi,dir(3),vec(3)
      real * 8 dircrossvec(3),dirdotvec
      integer i
      dircrossvec(1)=dir(2)*vec(3)-dir(3)*vec(2)
      dircrossvec(2)=dir(3)*vec(1)-dir(1)*vec(3)
      dircrossvec(3)=dir(1)*vec(2)-dir(2)*vec(1)
      dirdotvec=dir(1)*vec(1)+dir(2)*vec(2)+dir(3)*vec(3)
      do i=1,3
      vec(i)=vec(i)+sinphi*dircrossvec(i)
     &  -(1d0-cosphi)*(vec(i)-dir(i)*dirdotvec)
      enddo
      end


      subroutine findperp(p,j,e1,e2,e3)
      implicit none
      include 'constants.f'
      integer j
      double precision p(mxpart,4),e1(4),e2(4),dir(3),norm,
     & cosphi,sinphi,e3(4),n1(4),n2(4)
      data n1/1d0,0d0,0d0,0d0/
      data n2/0d0,1d0,0d0,0d0/
      sinphi=1d0
      cosphi=0d0
      if ((p(j,2) .eq. 0d0) .and. (p(j,1) .eq. 0d0)) then
      e1(:)=n1(:)
      e2(:)=n2(:)
      e3(:)=n1(:)+n2(:)
      return
      else
c     Construct e; First construct a vector in the plane of p_j and
c     the third axis, orthogonal to p_j. 
c     Then rotate it counterclockwise around the j direction 
      e1(1)=p(j,1)
      e1(2)=p(j,2)
      e1(3)=-(p(j,1)**2+p(j,2)**2)/p(j,3)
      norm=sqrt(e1(1)**2+e1(2)**2+e1(3)**2)
      e1(4)=0
      e1(:)=e1(:)/norm
      e2(:)=e1(:)
      dir(1)=p(j,1)/p(j,4)
      dir(2)=p(j,2)/p(j,4)
      dir(3)=p(j,3)/p(j,4)
      call rotate3vec(dir,sinphi,cosphi,e2(1))
      e3(:)=e1(:)+e2(:)
      endif

      end

