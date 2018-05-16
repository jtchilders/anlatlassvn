module madds
   use precision, only: ki, ki_ql, ki_lt
!AC!use precision_golem, only: ki_gol => ki
    use avh_olo_kinds, only: ki_avh => kindr2
   use constants
   use options
   use mfunctions
   use notfirst
   implicit none

   private

   interface add4
      module procedure add4_rm
      module procedure add4_cm
   end interface add4

   interface add3
      module procedure add3_rm
      module procedure add3_cm
   end interface add3

   interface add2
      module procedure add2_rm
      module procedure add2_cm
   end interface add2

   interface add1
      module procedure add1_rm
      module procedure add1_cm
   end interface add1

   ! If s_mat is allocated the addX routines will read their invariants
   ! from the matrix rather than to recompute them.
   !
   ! The matrix should be initialized as follows:
   !
   ! s_mat(i, j) = Vi(i-1).Vi(j-1) - msq(i) - msq(j)
   !
   !
   ! Example:
   !
   !   box diagram in gg>tt~
   !
   !   g(k1) ~~~~~~*~~~~~*====== t~(k4)
   !               S     I
   !               S     I
   !   g(k2) ~~~~~~*~~~~~*====== t(k3)
   !
   !   allocate(s_mat(4,4))
   !   s_mat(:,:) = 0.0_ki
   !   s_mat(1,3) = s      - 0.0_ki - 0.0_ki  ! = s
   !   s_mat(1,4) = mT**2  - 0.0_ki - mT**2   ! = 0.0_ki
   !   s_mat(2,4) = t      - 0.0_ki - mT**2   ! = t - mT**2
   !   s_mat(3,4) = mT**2  - 0.0_ki - mT**2   ! = 0.0_ki
   !   s_mat(4,4) = 0.0_ki - mT**2  - mT**2   ! = - 2.0_ki * mT**2
   !   call samurai( .... )
   !   deallocate(s_mat)

   complex(ki), dimension(:,:), allocatable, public :: s_mat

     interface
        function qlI4(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
           use precision, only: ki_ql
           implicit none
           real(ki_ql), intent(in) :: p1,p2,p3,p4,s12,s23
           real(ki_ql), intent(in) :: m1,m2,m3,m4,mu2
           integer, intent(in) :: ep
           complex(ki_ql) :: qlI4
        end function qlI4
     end interface
     interface
        function qlI3(p1,p2,p3,m1,m2,m3,mu2,ep)
           use precision, only: ki_ql
           implicit none
           real(ki_ql), intent(in) :: p1,p2,p3
           real(ki_ql), intent(in) :: m1,m2,m3,mu2
           integer, intent(in) :: ep
           complex(ki_ql) :: qlI3
        end function qlI3
     end interface
     interface
        function qlI2(p1,m1,m2,mu2,ep)
           use precision, only: ki_ql
           implicit none
           real(ki_ql), intent(in) :: p1
           real(ki_ql), intent(in) :: m1,m2,mu2
           integer, intent(in) :: ep
           complex(ki_ql) :: qlI2
        end function qlI2
     end interface
     interface
        function qlI1(m1,mu2,ep)
           use precision, only: ki_ql
           implicit none
           real(ki_ql), intent(in) :: m1,mu2
           integer, intent(in) :: ep
           complex(ki_ql) :: qlI1
        end function qlI1
     end interface

!AC!     interface
!AC!        function gD0(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
!AC!           use precision_golem, only: ki
!AC!           implicit none
!AC!           real(ki), intent(in) :: p1,p2,p3,p4,s12,s23
!AC!           real(ki), intent(in) :: m1,m2,m3,m4
!AC!           real(ki), intent(in) :: mu2
!AC!           integer, intent(in) :: ep
!AC!           complex(ki) :: gD0
!AC!        end function gD0
!AC!     end interface
!AC!     interface
!AC!        function gC0(p1,p2,p3,m1,m2,m3,mu2,ep)
!AC!           use precision_golem, only: ki
!AC!           implicit none
!AC!           real(ki), intent(in) :: p1,p2,p3
!AC!           real(ki), intent(in) :: m1,m2,m3
!AC!           real(ki), intent(in) :: mu2
!AC!           integer, intent(in) :: ep
!AC!           complex(ki) :: gC0
!AC!        end function gC0
!AC!     end interface
!AC!     interface
!AC!        function gB0(p1,m1,m2,mu2,ep)
!AC!           use precision_golem, only: ki
!AC!           implicit none
!AC!           real(ki), intent(in) :: p1
!AC!           real(ki), intent(in) :: m1,m2
!AC!           real(ki), intent(in) :: mu2
!AC!           integer, intent(in) :: ep
!AC!           complex(ki) :: gB0
!AC!        end function gB0
!AC!     end interface
!AC!     interface
!AC!        function gD0C(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
!AC!           use precision_golem, only: ki
!AC!           implicit none
!AC!           complex(ki), intent(in) :: p1,p2,p3,p4,s12,s23
!AC!           complex(ki), intent(in) :: m1,m2,m3,m4
!AC!           real(ki), intent(in) :: mu2
!AC!           integer, intent(in) :: ep
!AC!           complex(ki) :: gD0C
!AC!        end function gD0C
!AC!     end interface
!AC!     interface
!AC!        function gC0C(p1,p2,p3,m1,m2,m3,mu2,ep)
!AC!           use precision_golem, only: ki
!AC!           implicit none
!AC!           complex(ki), intent(in) :: p1,p2,p3
!AC!           complex(ki), intent(in) :: m1,m2,m3
!AC!           real(ki), intent(in) :: mu2
!AC!           integer, intent(in) :: ep
!AC!           complex(ki) :: gC0C
!AC!        end function gC0C
!AC!     end interface
!AC!     interface
!AC!        function gB0C(p1,m1,m2,mu2,ep)
!AC!           use precision_golem, only: ki
!AC!           implicit none
!AC!           complex(ki), intent(in) :: p1
!AC!           complex(ki), intent(in) :: m1,m2
!AC!           real(ki), intent(in) :: mu2
!AC!           integer, intent(in) :: ep
!AC!           complex(ki) :: gB0C
!AC!        end function gB0C
!AC!     end interface

!AC!     interface
!AC!        function D0(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           real(ki_lt), intent(in) :: p1,p2,p3,p4,s12,s23
!AC!           real(ki_lt), intent(in) :: m1,m2,m3,m4
!AC!           complex(ki_lt) :: D0
!AC!        end function D0
!AC!     end interface
!AC!     interface
!AC!        function C0(p1,p2,p3,m1,m2,m3)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           real(ki_lt), intent(in) :: p1,p2,p3
!AC!           real(ki_lt), intent(in) :: m1,m2,m3
!AC!           complex(ki_lt) :: C0
!AC!        end function C0
!AC!     end interface
!AC!     interface
!AC!        function B0(p1,m1,m2)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           real(ki_lt), intent(in) :: p1
!AC!           real(ki_lt), intent(in) :: m1,m2
!AC!           complex(ki_lt) :: B0
!AC!        end function B0
!AC!     end interface
!AC!     interface
!AC!        function B1(p1,m1,m2)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           real(ki_lt), intent(in) :: p1
!AC!           real(ki_lt), intent(in) :: m1,m2
!AC!           complex(ki_lt) :: B1
!AC!        end function B1
!AC!     end interface
!AC!     interface
!AC!        function B00(p1,m1,m2)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           real(ki_lt), intent(in) :: p1
!AC!           real(ki_lt), intent(in) :: m1,m2
!AC!           complex(ki_lt) :: B00
!AC!        end function B00
!AC!     end interface
!AC!     interface
!AC!        function B11(p1,m1,m2)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           real(ki_lt), intent(in) :: p1
!AC!           real(ki_lt), intent(in) :: m1,m2
!AC!           complex(ki_lt) :: B11
!AC!        end function B11
!AC!     end interface
!AC!     interface
!AC!        function A0(m1)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           real(ki_lt), intent(in) :: m1
!AC!           complex(ki_lt) :: A0
!AC!        end function A0
!AC!     end interface
!AC!     interface
!AC!        function D0C(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           complex(ki_lt), intent(in) :: p1,p2,p3,p4,s12,s23
!AC!           complex(ki_lt), intent(in) :: m1,m2,m3,m4
!AC!           complex(ki_lt) :: D0C
!AC!        end function D0C
!AC!     end interface
!AC!     interface
!AC!        function C0C(p1,p2,p3,m1,m2,m3)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           complex(ki_lt), intent(in) :: p1,p2,p3
!AC!           complex(ki_lt), intent(in) :: m1,m2,m3
!AC!           complex(ki_lt) :: C0C
!AC!        end function C0C
!AC!     end interface
!AC!     interface
!AC!        function B0C(p1,m1,m2)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           complex(ki_lt), intent(in) :: p1
!AC!           complex(ki_lt), intent(in) :: m1,m2
!AC!           complex(ki_lt) :: B0C
!AC!        end function B0C
!AC!     end interface
!AC!     interface
!AC!        function B1C(p1,m1,m2)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           complex(ki_lt), intent(in) :: p1
!AC!           complex(ki_lt), intent(in) :: m1,m2
!AC!           complex(ki_lt) :: B1C
!AC!        end function B1C
!AC!     end interface
!AC!     interface
!AC!        function B00C(p1,m1,m2)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           complex(ki_lt), intent(in) :: p1
!AC!           complex(ki_lt), intent(in) :: m1,m2
!AC!           complex(ki_lt) :: B00C
!AC!        end function B00C
!AC!     end interface
!AC!     interface
!AC!        function B11C(p1,m1,m2)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           complex(ki_lt), intent(in) :: p1
!AC!           complex(ki_lt), intent(in) :: m1,m2
!AC!           complex(ki_lt) :: B11C
!AC!        end function B11C
!AC!     end interface
!AC!     interface
!AC!        function A0C(m1)
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           complex(ki_lt), intent(in) :: m1
!AC!           complex(ki_lt) :: A0C
!AC!        end function A0C
!AC!     end interface
!AC!     interface
!AC!        function getlambda()
!AC!           use precision, only: ki_lt
!AC!           implicit none
!AC!           real(ki_lt) :: getlambda
!AC!        end function getlambda
!AC!     end interface

   ! cache sizes depending on nleg and istop: cachedim<nleg>(<istop>)
   ! These numbers are exactly the same as returned in calls
   ! to cachedim
   integer, dimension(1:1), parameter, public :: cachedim1 = (/1/)
   integer, dimension(1:2), parameter, public :: cachedim2 = (/7,5/)
   integer, dimension(1:3), parameter, public :: cachedim3 = (/19,16,1/)
   integer, dimension(1:4), parameter, public :: cachedim4 = (/39,35,5,1/)
   integer, dimension(1:4), parameter, public :: cachedim5 = (/70,65,15,5/)
   integer, dimension(1:4), parameter, public :: cachedim6 = (/116,110,35,15/)
   integer, dimension(1:4), parameter, public :: cachedim7 = (/182,175,70,35/)
   integer, dimension(1:4), parameter, public :: cachedim8 = (/274,266,126,70/)
   integer, dimension(1:4), parameter, public :: cachedim9 = &
  & (/399,390,210,126/)
   integer, dimension(1:4), parameter, public :: cachedim10 = &
  & (/565,555,330,210/)
   integer, dimension(1:4), parameter, public :: cachedim11 = &
  & (/781,770,495,330/)
   integer, dimension(1:4), parameter, public :: cachedim12 = &
  & (/1057,1045,715,495/)

   public :: add4, add3, add2, add1
   public :: add4_rm, add3_rm, add2_rm, add1_rm
   public :: add4_cm, add3_cm, add2_cm, add1_cm

   public :: cachedim

contains

   pure subroutine cachedim(dim,nleg,istop)
      implicit none
      integer, intent(in) :: nleg,istop
      integer, intent(out) :: dim
      integer :: n4, n3, n2, n1, j1, j2, j3, j4, icut1, icut2, icut3, icut4

      n1 = 0
      n2 = 0
      n3 = 0
      n4 = 0

      if (nleg.ge.4) then
         goto 20
      elseif (nleg.eq.3) then
         goto 30
      elseif (nleg.eq.2) then
         goto 40
      elseif ((nleg.eq.1).or.(nleg.le.0)) then
         goto 50
      endif

 20   continue

      if (istop.ge.5) goto 99

      n4 = nleg*(nleg-1)*(nleg-2)*(nleg-3)/24
      if (istop.ge.4) goto 99

 30   continue

      n3 = nleg*(nleg-1)*(nleg-2)/6
      if (istop.ge.3) goto 99

 40   continue

      n2 = nleg*(nleg-1)/2
      if (istop.ge.2) goto 99

 50    continue

      n1 = nleg

 99   continue

      dim=n4+n3+5*n2+n1

  end subroutine cachedim

  subroutine add4_rm(nleg,c4,cut4,Vi,msq,tot4,tot4r,scale2,&
        cache_flag, cache_offset, scalar_cache)
          use avh_olo, only: olo_d0
      implicit none
      integer, intent(in) :: cut4,nleg
      complex(ki), dimension(0:4), intent(in) :: c4
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi
      real(ki), dimension(0:nleg-1), intent(in) :: msq
      real(ki), intent(in) :: scale2
      complex(ki), dimension(-2:0), intent(out) :: tot4
      complex(ki), intent(out) :: tot4r

      logical, intent(in), optional :: cache_flag
      integer, intent(inout), optional :: cache_offset
      complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache

      real(ki), dimension(4):: Vi1, Vi2, Vi3, Vi21, Vi31, Vi32
      real(ki) :: m0, m1, m2, m3, V1, V2, V3, V21, V31, V32
      integer :: i,j1,j2,j3,j4
      complex(ki) :: c40
          complex(ki_avh), dimension(0:2) :: vald0
      complex(ki) :: ctmp
!AC!      complex(ki), dimension(-2:0) :: d0t
      integer :: ep, cache_index

      if (notfirsti.eqv.(.false.)) then
         if (isca .eq. 2) then
                call avh_olo_mu_set(real(sqrt(scale2),ki_avh))
         elseif (isca .eq. 4) then
!AC!             call setmudim(real(scale2, ki_lt))
         endif
         notfirsti=.true.
      endif


      j4=cut4/1000
      j3=(cut4-j4*1000)/100
      j2=(cut4-j4*1000-j3*100)/10
      j1=cut4-j4*1000-j3*100-j2*10

      m0=msq(j1)
      m1=msq(j2)
      m2=msq(j3)
      m3=msq(j4)

      if (allocated(s_mat)) then
         ! s_mat(i+1, j+1) = (Vi(i,:) - Vi(j,:))**2 - msq(i) - msq(j)
         V1  = s_mat(j2+1, j1+1) + msq(j2) + msq(j1)
         V2  = s_mat(j3+1, j1+1) + msq(j3) + msq(j1)
         V3  = s_mat(j4+1, j1+1) + msq(j4) + msq(j1)
         V21 = s_mat(j3+1, j2+1) + msq(j3) + msq(j2)
         V31 = s_mat(j4+1, j2+1) + msq(j4) + msq(j2)
         V32 = s_mat(j4+1, j3+1) + msq(j4) + msq(j3)
      else
         Vi1(:)=Vi(j2,:)-Vi(j1,:)
         Vi2(:)=Vi(j3,:)-Vi(j1,:)
         Vi3(:)=Vi(j1,:)-Vi(j4,:)
         Vi21(:)=Vi(j3,:)-Vi(j2,:)
         Vi31(:)=Vi(j4,:)-Vi(j2,:)
         Vi32(:)=Vi(j4,:)-Vi(j3,:)

         V1=sdot(Vi1,Vi1)
         V2=sdot(Vi2,Vi2)
         V3=sdot(Vi3,Vi3)
         V21=sdot(Vi21,Vi21)
         V31=sdot(Vi31,Vi31)
         V32=sdot(Vi32,Vi32)
      end if

      c40=c4(0)
      tot4r=-c4(4)/six

  1   Format(A3,I4,A1,I2,A5,D24.15,A1,D24.15,A3)

      if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
      if     (isca.eq.1) then
          do ep=-2,0
             if (present(cache_flag)) then
                if (cache_flag) then
                   ctmp = scalar_cache(ep,cache_index)
                else
                   ctmp=qlI4(&
                   & real(V1,ki_ql),real(V21,ki_ql),real(V32,ki_ql),&
                   & real(V3,ki_ql),real(V2,ki_ql),real(V31,ki_ql),&
                   & real(m0,ki_ql),real(m1,ki_ql),real(m2,ki_ql),&
                   & real(m3,ki_ql),real(scale2,ki_ql),ep)
                   scalar_cache(ep,cache_index) = ctmp
                end if
             else
                ctmp=qlI4(&
                & real(V1,ki_ql),real(V21,ki_ql),real(V32,ki_ql),&
                & real(V3,ki_ql),real(V2,ki_ql),real(V31,ki_ql),&
                & real(m0,ki_ql),real(m1,ki_ql),real(m2,ki_ql),&
                & real(m3,ki_ql),real(scale2,ki_ql),ep)
             end if
             tot4(ep)=c40*ctmp
             if (verbosity.ge.2) write(iout,1) &
         & 'I4(',cut4,',',ep,') = (',real(ctmp),',',aimag(ctmp),'  )' 
          enddo
!AC!      print*, "isca=1: QCDLoop not available"
!AC!      stop
      elseif (isca.eq.2) then
             if (present(cache_flag)) then
                if (cache_flag) then
                   vald0(0) = scalar_cache( 0,cache_index)
                   vald0(1) = scalar_cache(-1,cache_index)
                   vald0(2) = scalar_cache(-2,cache_index)
                else
                   call olo_d0(vald0,&
                & real(V1,ki_avh),real(V21,ki_avh),real(V32,ki_avh),&
                & real(V3,ki_avh),real(V2,ki_avh),real(V31,ki_avh),&
                & real(m0,ki_avh),real(m1,ki_avh),real(m2,ki_avh),&
                & real(m3,ki_avh))
                   scalar_cache( 0,cache_index) = vald0(0)
                   scalar_cache(-1,cache_index) = vald0(1)
                   scalar_cache(-2,cache_index) = vald0(2)
                end if
             else
                call olo_d0(vald0,&
                & real(V1,ki_avh),real(V21,ki_avh),real(V32,ki_avh),&
                & real(V3,ki_avh),real(V2,ki_avh),real(V31,ki_avh),&
                & real(m0,ki_avh),real(m1,ki_avh),real(m2,ki_avh),&
                & real(m3,ki_avh))
             end if
             do ep=-2,0
                tot4(ep)= c40*vald0(-ep) 
                if (verbosity.ge.2) write(iout,1) &
                & 'I4(',cut4,',',ep,') = (',real(vald0(-ep)),',',aimag(vald0(-ep)),'  )' 
             enddo
!AC!      print*, "isca=2: OneLOop not available"
!AC!      stop
      elseif (isca.eq.3) then
!AC!        call gtrunc_rm(abs(V32)+abs(V31), &
!AC!        &   V1,V2,V3,V21,V32,V31,m0,m1,m2,m3)
!AC!        do ep=-2,0
!AC!           if (present(cache_flag)) then
!AC!              if (cache_flag) then
!AC!                 d0t(ep) = scalar_cache(ep,cache_index)
!AC!              else
!AC!                 d0t(ep)=gD0(real(V1,ki_gol),real(V21,ki_gol),&
!AC!                 &        real(V32,ki_gol),real(V3,ki_gol),&
!AC!                 &        real(V2,ki_gol),real(V31,ki_gol),&
!AC!                 &        real(m0,ki_gol),real(m1,ki_gol),&
!AC!                 &        real(m2,ki_gol),real(m3,ki_gol),&
!AC!                 &        real(scale2,ki_gol),ep)
!AC!                 scalar_cache(ep,cache_index) = d0t(ep)
!AC!              end if
!AC!           else
!AC!              d0t(ep)=gD0(real(V1,ki_gol),real(V21,ki_gol),&
!AC!              &        real(V32,ki_gol),real(V3,ki_gol),&
!AC!              &        real(V2,ki_gol),real(V31,ki_gol),&
!AC!              &        real(m0,ki_gol),real(m1,ki_gol),&
!AC!              &        real(m2,ki_gol),real(m3,ki_gol),&
!AC!              &        real(scale2,ki_gol),ep)
!AC!           end if
!AC!        end do
!AC!        !d0t( 0) = d0t(0) + log(scale2) * (d0t(-1) &
!AC!        !      & + 0.5_ki * log(scale2) * d0t(-2))
!AC!        !d0t(-1) = d0t(-1) + log(scale2) * d0t(-2)
!AC!        if (verbosity.ge.2) then
!AC!           do ep=-2,0
!AC!              write(iout,1) 'I4(',cut4,',',ep,&
!AC!              & ') = (',real(d0t(ep)),',',aimag(d0t(ep)),'  )' 
!AC!           end do
!AC!        end if
!AC!        tot4(:) = d0t(:) * c40
            print*, "isca=3: Golem95 not available"
            stop
      elseif (isca.eq.4) then
!AC!        tot4(-2) = 0
!AC!        tot4(-1) = 0
!AC!        tot4(0)  = 0
!AC!        ep = -dim(0, int(getlambda()))
!AC!        if (present(cache_flag)) then
!AC!           if (cache_flag) then
!AC!              ctmp = scalar_cache(ep,cache_index)
!AC!           else
!AC!              call gtrunc_rm(abs(V32)+abs(V31), &
!AC!              &   V1,V2,V3,V21,V32,V31,m0,m1,m2,m3)
!AC!              ctmp=D0(&
!AC!              & real(V1,ki_lt),real(V21,ki_lt),real(V32,ki_lt),&
!AC!              & real(V3,ki_lt),real(V2,ki_lt),real(V31,ki_lt),&
!AC!              & real(m0,ki_lt),real(m1,ki_lt),real(m2,ki_lt),&
!AC!              & real(m3,ki_lt))
!AC!              scalar_cache(ep,cache_index) = ctmp
!AC!           end if
!AC!        else
!AC!           call gtrunc_rm(abs(V32)+abs(V31), &
!AC!           &   V1,V2,V3,V21,V32,V31,m0,m1,m2,m3)
!AC!           ctmp=D0(&
!AC!           & real(V1,ki_lt),real(V21,ki_lt),real(V32,ki_lt),&
!AC!           & real(V3,ki_lt),real(V2,ki_lt),real(V31,ki_lt),&
!AC!           & real(m0,ki_lt),real(m1,ki_lt),real(m2,ki_lt),&
!AC!           & real(m3,ki_lt))
!AC!        end if
!AC!        tot4(ep)=c40*ctmp
!AC!        if (verbosity.ge.2) write(iout,1) &
!AC!     & 'I4(',cut4,',',ep,') = (',real(ctmp),',',aimag(ctmp),'  )' 
             print*, "isca=4: LoopTools not available"
             stop
      else
         print*, 'error in add4'
         stop
      endif
      if (present(cache_flag)) cache_offset = cache_offset + 1
      tot4(0)=tot4(0) + tot4r
  end subroutine add4_rm

  subroutine add3_rm(nleg,c3,cut3,Vi,msq,tot3,tot3r,scale2,&
        cache_flag, cache_offset, scalar_cache)
          use avh_olo, only: olo_c0
      implicit none

      integer, intent(in) :: nleg, cut3
      complex(ki), dimension(0:9), intent(in) :: c3
      real(ki), dimension(0:nleg-1,4) ::Vi
      real(ki), dimension(0:nleg-1):: msq
      complex(ki), dimension(-2:0), intent(out) :: tot3
      complex(ki), intent(out) :: tot3r
      real(ki), intent(in) :: scale2

      logical, intent(in), optional :: cache_flag
      integer, intent(inout), optional :: cache_offset
      complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache

      integer :: j1,j2,j3
      real(ki) :: m0, m1, m2, V1, V2, V3
      real(ki), dimension(4):: Vi1, Vi2, Vi3
      complex(ki) :: c30
          complex(ki_avh), dimension(0:2) :: valc0
      complex(ki) :: ctmp
!AC!      complex(ki), dimension(-2:0) :: c0t
      integer :: ep, cache_index

      if (notfirsti.eqv.(.false.)) then
         if (isca .eq. 2) then
                call avh_olo_mu_set(real(sqrt(scale2),ki_avh))
         elseif (isca .eq. 4) then
!AC!             call setmudim(real(scale2, ki_lt))
         endif
         notfirsti=.true.
      endif


      j3=cut3/100
      j2=(cut3-j3*100)/10
      j1=cut3-j3*100-j2*10

      m0=msq(j1)
      m1=msq(j2)
      m2=msq(j3)

      if (allocated(s_mat)) then
         ! s_mat(i+1, j+1) = (Vi(i,:) - Vi(j,:))**2 - msq(i) - msq(j)
         V1 = s_mat(j2+1, j1+1) + msq(j2) + msq(j1)
         V2 = s_mat(j3+1, j2+1) + msq(j3) + msq(j2)
         V3 = s_mat(j3+1, j1+1) + msq(j3) + msq(j1)
      else
         Vi1(:)=Vi(j2,:)-Vi(j1,:)
         Vi2(:)=Vi(j3,:)-Vi(j2,:)
         Vi3(:)=Vi(j1,:)-Vi(j3,:)

         V1=sdot(Vi1,Vi1)
         V2=sdot(Vi2,Vi2)
         V3=sdot(Vi3,Vi3)
      end if

      c30=c3(0)

      tot3r=+c3(7)/two

  1  Format(A3,I3,A1,I2,A5,D24.15,A1,D24.15,A3)

      if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
      if (isca.eq.1) then
          do ep=-2,0
             if (present(cache_flag)) then
                if (cache_flag) then
                   ctmp = scalar_cache(ep,cache_index)
                else
                   ctmp=qlI3(&
                      real(V1,ki_ql),real(V2,ki_ql),real(V3,ki_ql),&
                      & real(m0,ki_ql),real(m1,ki_ql),real(m2,ki_ql),&
                      & real(scale2,ki_ql),ep)
                   scalar_cache(ep,cache_index) = ctmp
                end if
             else
                ctmp=qlI3(&
                   real(V1,ki_ql),real(V2,ki_ql),real(V3,ki_ql),&
                   & real(m0,ki_ql),real(m1,ki_ql),real(m2,ki_ql),&
                   & real(scale2,ki_ql),ep)
             end if
             tot3(ep)=c30*ctmp 
             if (verbosity.ge.2) write(iout,1) &
             &'I3(',cut3,',',ep,') = (',real(ctmp),&
             &',',aimag(ctmp),'  )' 
          enddo
!AC!      print*, "isca=1: QCDLoop not available"
!AC!      stop
      elseif  (isca.eq.2) then
              if (present(cache_flag)) then
                 if (cache_flag) then
                    valc0(0) = scalar_cache( 0,cache_index)
                    valc0(1) = scalar_cache(-1,cache_index)
                    valc0(2) = scalar_cache(-2,cache_index)
                 else
                    call olo_c0(valc0,&
                 & real(V1,ki_avh),real(V2,ki_avh),real(V3,ki_avh),&
                 & real(m0,ki_avh),real(m1,ki_avh),real(m2,ki_avh))
                    scalar_cache( 0,cache_index) = valc0(0)
                    scalar_cache(-1,cache_index) = valc0(1)
                    scalar_cache(-2,cache_index) = valc0(2)
                 end if
              else
                 call olo_c0(valc0,&
                 & real(V1,ki_avh),real(V2,ki_avh),real(V3,ki_avh),&
                 & real(m0,ki_avh),real(m1,ki_avh),real(m2,ki_avh))
              end if
              do ep=-2,0
                 tot3(ep)= c30*valc0(-ep) 
                 if (verbosity.ge.2) write(iout,1) &
         &'I3(',cut3,',',ep,') = (',real(valc0(-ep)),',',aimag(valc0(-ep)),'  )' 
              enddo
!AC!      print*, "isca=2: OneLOop not available"
!AC!      stop
      elseif (isca.eq.3) then
!AC!        call gtrunc_rm(abs(V1)+abs(V2)+abs(V3), &
!AC!        &                V1,V2,V3,m0,m1,m2)
!AC!        do ep=-2,0
!AC!           if (present(cache_flag)) then
!AC!              if (cache_flag) then
!AC!                 c0t(ep) = scalar_cache(ep,cache_index)
!AC!              else
!AC!                 c0t(ep)=gC0(real(V1,ki_gol),real(V2,ki_gol),&
!AC!                 &        real(V3,ki_gol),real(m0,ki_gol),&
!AC!                 &        real(m1,ki_gol),real(m2,ki_gol),&
!AC!                 &        real(scale2,ki_gol),ep)
!AC!                 scalar_cache(ep,cache_index) = c0t(ep)
!AC!              end if
!AC!           else
!AC!              c0t(ep)=gC0(real(V1,ki_gol),real(V2,ki_gol),&
!AC!              &        real(V3,ki_gol),real(m0,ki_gol),&
!AC!              &        real(m1,ki_gol),real(m2,ki_gol),&
!AC!              &        real(scale2,ki_gol),ep)
!AC!           end if
!AC!        end do
!AC!        !c0t( 0) = c0t(0) + log(scale2) * (c0t(-1) &
!AC!        !      & + 0.5_ki * log(scale2) * c0t(-2))
!AC!        !c0t(-1) = c0t(-1) + log(scale2) * c0t(-2)
!AC!        if (verbosity.ge.2) then
!AC!           do ep=-2,0
!AC!              write(iout,1) 'I3(',cut3,',',ep,&
!AC!             & ') = (',real(c0t(ep)),',',aimag(c0t(ep)),'  )' 
!AC!           end do
!AC!        end if
!AC!        tot3(:) = c0t(:) * c30
            print*, "isca=3: Golem95 not available"
            stop
      elseif (isca.eq.4) then
!AC!        tot3(-2) = 0
!AC!        tot3(-1) = 0
!AC!        tot3( 0) = 0
!AC!        call gtrunc_rm(abs(V1)+abs(V2)+abs(V3), &
!AC!        &  V1,V2,V3,m0,m1,m2)
!AC!        ep = -dim(0, int(getlambda()))
!AC!        if (present(cache_flag)) then
!AC!           if (cache_flag) then
!AC!              ctmp = scalar_cache(ep,cache_index)
!AC!           else
!AC!              ctmp=C0(&
!AC!              &  real(V1,ki_lt),real(V2,ki_lt),real(V3,ki_lt),&
!AC!              &  real(m0,ki_lt),real(m1,ki_lt),real(m2,ki_lt))
!AC!              scalar_cache(ep,cache_index) = ctmp
!AC!           end if
!AC!        else
!AC!           ctmp=C0(&
!AC!           &  real(V1,ki_lt),real(V2,ki_lt),real(V3,ki_lt),&
!AC!           &  real(m0,ki_lt),real(m1,ki_lt),real(m2,ki_lt))
!AC!        end if
!AC!        tot3(ep)=c30*ctmp 
!AC!        if (verbosity.ge.2) write(iout,1) &
!AC!       &'I3(',cut3,',',ep,') = (',real(ctmp),',',aimag(ctmp),'  )'
             print*, "isca=4: LoopTools not available"
             stop
      else
         print*, 'error in add3'
         stop
      endif
      tot3(0)=tot3(0) + tot3r
      if (present(cache_flag)) cache_offset = cache_offset + 1
  end subroutine add3_rm

  subroutine add2_rm(nleg,c2,cut2,k1,k2,msq,tot2,tot2r,scale2, &
        cache_flag, cache_offset, scalar_cache)
          use avh_olo, only: olo_b11
      implicit none

      integer, intent(in) :: nleg, cut2
      complex(ki), dimension(0:9), intent(in) :: c2
      real(ki), dimension(4), intent(in) :: k1, k2
      real(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(-2:0), intent(out) :: tot2
      complex(ki), intent(out) :: tot2r
      real(ki), intent(in) :: scale2

      logical, intent(in), optional :: cache_flag
      integer, intent(inout), optional :: cache_offset
      complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache

      real(ki) :: m0, m1, K11, K12, B06
      integer :: ep, cache_index
      integer :: i1,i2
      complex(ki), dimension(-2:0) :: J0, J1, J00, J01, J11 
          complex(ki_avh), dimension(0:2) :: scf2, scf1, scf0, scf 
          complex(ki) :: xbb, xb0, xb00
          real(ki) :: bkv
          integer :: i

      if (notfirsti.eqv.(.false.)) then
         if (isca .eq. 2) then
                call avh_olo_mu_set(real(sqrt(scale2),ki_avh))
         elseif (isca .eq. 4) then
!AC!             call setmudim(real(scale2, ki_lt))
         endif
         notfirsti=.true.
      endif

      i2=cut2/10
      i1=cut2-i2*10

      m0=msq(i1)
      m1=msq(i2)

      if (allocated(s_mat)) then
         K11 = s_mat(i2+1, i1+1) + msq(i1) + msq(i2)
      else
         K11 = sdot(K1,K1)
      end if
      K12=sdot(K1,K2)


      B06=-(K11-three*(m0+m1))/six

      tot2r= + B06*c2(9)

      if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
      if (isca.eq.1) then
          if (present(cache_flag)) then
             if (cache_flag) then
                J0(:)  = scalar_cache(:,cache_index+0)
                J1(:)  = scalar_cache(:,cache_index+1)
                J01(:) = scalar_cache(:,cache_index+2)
                J11(:) = scalar_cache(:,cache_index+3)
                J00(:) = scalar_cache(:,cache_index+4)
             else
                do ep=-2,0
                   J00(ep) = qlI2(&
               &     real(K11,ki_ql),real(m0,ki_ql),real(m0,ki_ql),&
               &     real(scale2,ki_ql),ep)
                   J11(ep) = qlI2(&
               &     real(K11,ki_ql),real(m1,ki_ql),real(m1,ki_ql),&
               &     real(scale2,ki_ql),ep)
                   J01(ep)=qlI2(&
               &     real(K11,ki_ql),real(m0,ki_ql),real(m1,ki_ql),&
               &     real(scale2,ki_ql),ep)
                   J0(ep) = qlI1(real(m0,ki_ql),real(scale2,ki_ql),ep)
                   J1(ep) = qlI1(real(m1,ki_ql),real(scale2,ki_ql),ep)
                end do
                scalar_cache(:,cache_index+0) = J0(:)
                scalar_cache(:,cache_index+1) = J1(:)
                scalar_cache(:,cache_index+2) = J01(:)
                scalar_cache(:,cache_index+3) = J11(:)
                scalar_cache(:,cache_index+4) = J00(:)
             end if
          else
             do ep=-2,0
                J00(ep) = qlI2(&
               &     real(K11,ki_ql),real(m0,ki_ql),real(m0,ki_ql),&
               &     real(scale2,ki_ql),ep)
                J11(ep) = qlI2(&
               &     real(K11,ki_ql),real(m1,ki_ql),real(m1,ki_ql),&
               &     real(scale2,ki_ql),ep)
                J01(ep)=qlI2(&
               &     real(K11,ki_ql),real(m0,ki_ql),real(m1,ki_ql),&
               &     real(scale2,ki_ql),ep)
                J0(ep) = qlI1(real(m0,ki_ql),real(scale2,ki_ql),ep)
                J1(ep) = qlI1(real(m1,ki_ql),real(scale2,ki_ql),ep)
             end do
          end if
          ! The remaining steps come in another if-statement
!AC!      print*, "isca=1: QCDLoop not available"
!AC!      stop
       elseif  (isca.eq.2) then
             xbb=c2(0)
             xb0=c2(1)
             xb00=c2(2)
             bkv=K12
             if (present(cache_flag)) then
                if (cache_flag) then
                   scf(:)  = scalar_cache(:,cache_index+0)
                   scf0(:) = scalar_cache(:,cache_index+1)
                   scf1(:) = scalar_cache(:,cache_index+2)
                   scf2(:) = scalar_cache(:,cache_index+3)
                else
                   call olo_b11(scf2,scf0,scf1,scf,&
                  & real(K11,ki_avh),real(m0,ki_avh),real(m1,ki_avh))
                   scalar_cache(:,cache_index+0) = scf(:)
                   scalar_cache(:,cache_index+1) = scf0(:)
                   scalar_cache(:,cache_index+2) = scf1(:)
                   scalar_cache(:,cache_index+3) = scf2(:)
                end if
             else
                call olo_b11(scf2,scf0,scf1,scf,&
               & real(K11,ki_avh),real(m0,ki_avh),real(m1,ki_avh))
             end if
             tot2(0)=xbb*scf(0)+xb0*bkv*scf1(0)+xb00*bkv*bkv*scf2(0)
             tot2(0)=tot2(0)+ B06*c2(9)
             tot2(-1)=xbb*scf(1)+xb0*bkv*scf1(1)+xb00*bkv*bkv*scf2(1)
             tot2(-2)=xbb*scf(2)+xb0*bkv*scf1(2)+xb00*bkv*bkv*scf2(2)
              if (verbosity.ge.2) then
                 do i=0,2
                    write(iout,903) 'B_0 (',cut2,',',-i,') = ',scf(i)
                    write(iout,903) 'B_1 (',cut2,',',-i,') =',scf1(i)
                    write(iout,903) 'B_11(',cut2,',',-i,') =',scf2(i)
                 enddo
              endif
!AC!      print*, "isca=2: OneLOop not available"
!AC!      stop
      elseif (isca.eq.3) then
!AC!        if (present(cache_flag)) then
!AC!           if (cache_flag) then
!AC!              J0(:)  = scalar_cache(:,cache_index+0)
!AC!              J1(:)  = scalar_cache(:,cache_index+1)
!AC!              J01(:) = scalar_cache(:,cache_index+2)
!AC!              J11(:) = scalar_cache(:,cache_index+3)
!AC!              J00(:) = scalar_cache(:,cache_index+4)
!AC!           else
!AC!              call gtrunc_rm(abs(K11)+1.0_ki, K11,m0,m1)
!AC!              do ep=-2,0
!AC!                 J00(ep)= gB0(real(K11,ki_gol),&
!AC!             &           real(m0,ki_gol),real(m0,ki_gol),&
!AC!             &           real(scale2,ki_gol),ep)
!AC!                 J11(ep)= gB0(real(K11,ki_gol),&
!AC!             &           real(m1,ki_gol),real(m1,ki_gol),&
!AC!             &           real(scale2,ki_gol),ep)
!AC!                 J01(ep)= gB0(real(K11,ki_gol),&
!AC!             &           real(m0,ki_gol),real(m1,ki_gol),&
!AC!             &           real(scale2,ki_gol),ep)
!AC!                 J0(ep) = gA0(real(m0,ki_gol),&
!AC!             &           real(scale2,ki_gol),ep)
!AC!                 J1(ep) = gA0(real(m1,ki_gol),&
!AC!             &           real(scale2,ki_gol),ep)
!AC!              end do
!AC!              scalar_cache(:,cache_index+0) = J0(:)
!AC!              scalar_cache(:,cache_index+1) = J1(:)
!AC!              scalar_cache(:,cache_index+2) = J01(:)
!AC!              scalar_cache(:,cache_index+3) = J11(:)
!AC!              scalar_cache(:,cache_index+4) = J00(:)
!AC!           end if
!AC!        else
!AC!           call gtrunc_rm(abs(K11)+1.0_ki, K11,m0,m1)
!AC!           do ep=-2,0
!AC!              J00(ep)= gB0(real(K11,ki_gol),&
!AC!          &           real(m0,ki_gol),real(m0,ki_gol),&
!AC!          &           real(scale2,ki_gol),ep)
!AC!              J11(ep)= gB0(real(K11,ki_gol),&
!AC!          &           real(m1,ki_gol),real(m1,ki_gol),&
!AC!          &           real(scale2,ki_gol),ep)
!AC!              J01(ep)= gB0(real(K11,ki_gol),&
!AC!          &           real(m0,ki_gol),real(m1,ki_gol),&
!AC!          &           real(scale2,ki_gol),ep)
!AC!              J0(ep) = gA0(real(m0,ki_gol),&
!AC!          &           real(scale2,ki_gol),ep)
!AC!              J1(ep) = gA0(real(m1,ki_gol),&
!AC!          &           real(scale2,ki_gol),ep)
!AC!           end do
!AC!        end if
            print*, "isca=3: Golem95 not available"
            stop
      elseif (isca.eq.4) then
!AC!        ep = -dim(0, int(getlambda()))
!AC!        J00(:) = 0.0_ki_lt
!AC!        J11(:) = 0.0_ki_lt
!AC!        J0(:) = 0.0_ki_lt
!AC!        J1(:) = 0.0_ki_lt
!AC!        if (present(cache_flag)) then
!AC!           if (cache_flag) then
!AC!              J0(:)  = scalar_cache(:,cache_index+0)
!AC!              J1(:)  = scalar_cache(:,cache_index+1)
!AC!              J01(:) = scalar_cache(:,cache_index+2)
!AC!              J11(:) = scalar_cache(:,cache_index+3)
!AC!              J00(:) = scalar_cache(:,cache_index+4)
!AC!           else
!AC!              call gtrunc_rm(abs(K11)+1.0_ki, K11,m0,m1)
!AC!              J00(ep) = B0(real(K11,ki_lt),&
!AC!             &             real(m0,ki_lt),real(m0,ki_lt))
!AC!              J11(ep) = B0(real(K11,ki_lt),&
!AC!             &             real(m1,ki_lt),real(m1,ki_lt))
!AC!              J01(ep) = B0(real(K11,ki_lt),&
!AC!             &             real(m0,ki_lt),real(m1,ki_lt))
!AC!              J0(ep)  = A0(real(m0,ki_lt))
!AC!              J1(ep)  = A0(real(m1,ki_lt))
!AC!              scalar_cache(:,cache_index+0) = J0(:)
!AC!              scalar_cache(:,cache_index+1) = J1(:)
!AC!              scalar_cache(:,cache_index+2) = J01(:)
!AC!              scalar_cache(:,cache_index+3) = J11(:)
!AC!              scalar_cache(:,cache_index+4) = J00(:)
!AC!           end if
!AC!        else
!AC!           call gtrunc_rm(abs(K11)+1.0_ki, K11,m0,m1)
!AC!           J00(ep) = B0(real(K11,ki_lt),&
!AC!          &             real(m0,ki_lt),real(m0,ki_lt))
!AC!           J11(ep) = B0(real(K11,ki_lt),&
!AC!          &             real(m1,ki_lt),real(m1,ki_lt))
!AC!           J01(ep) = B0(real(K11,ki_lt),&
!AC!          &             real(m0,ki_lt),real(m1,ki_lt))
!AC!           J0(ep)  = A0(real(m0,ki_lt))
!AC!           J1(ep)  = A0(real(m1,ki_lt))
!AC!        end if
            print*, "isca=4: LoopTools not available"
            stop
       else
          print*, 'error in add2'
          stop
       endif

       if (isca.eq.1 .or. isca.eq.3 .or. isca.eq.4) then
          if (abs(K11).gt.zip1) then
             do ep=-2,0
                tot2(ep)=-(K12*(two*K12*(m0 - m1)*c2(2) +  &
        &         K11*(-three*c2(1) + two*K12*c2(2)))*J0(ep))/(six*K11**2) &
        &         + ((two*K12**2*(m0 - m1)**2*c2(2) +  &
        &         K11*K12*(-three*m0*c2(1) + three*m1*c2(1) +  &
        &         two*K12*m0*c2(2) - four*K12*m1*c2(2)) +  &
        &         K11**2*(6*c2(0) + K12*(-three*c2(1) + two*K12*c2(2))))* &
        &         J01(ep))/(six*K11**2) +  &
        &         (K12*(two*K12*(m0 - m1)*c2(2) +  &
        &         K11*(-three*c2(1) + four*K12*c2(2)))*J1(ep))/(six*K11**2)
             enddo
             tot2(0)=tot2(0)+(K12**2*c2(2))/18.0_ki &
        &       - (K12**2*m0*c2(2))/(six*K11) -  &
        &       (K12**2*m1*c2(2))/(six*K11) + B06*c2(9)
          else
             if (m1.eq.m0) then
                do ep=-2,0
                   tot2(ep)=(c2(0) + (K12*(-three*c2(1) &
        &                    + two*K12*c2(2)))/six)*J00(ep)
                enddo
                tot2(0)=tot2(0)+B06*c2(9)
             else
                do ep=-2,0
                   tot2(ep)=(K12*m0**2*(-three*m0*c2(1) + three*m1*c2(1) &
        &              + two*K12*m0*c2(2))* &
        &              J00(ep))/(six*(m0 - m1)**3) + c2(0)*J01(ep) -  &
        &              (K12*m1*(m0*m1*(9*c2(1) - 6*K12*c2(2)) -  &
        &              6*m0**2*(c2(1) - K12*c2(2)) +  &
        &              m1**2*(-three*c2(1) + two*K12*c2(2)))*J11(ep))/ &
        &              (six*(m0 - m1)**3)
                enddo
                tot2(0)=tot2(0) + (-three*K12*m0*c2(1))/(four*(m0 - m1)) +  &
        &           (K12*m1*c2(1))/(four*m0 - four*m1) +  &
        &           (11.0_ki*K12**2*m0**2*c2(2))/(18.0_ki*(m0 - m1)**2) -  &
        &           (7.0_ki*K12**2*m0*m1*c2(2))/(18.0_ki*(m0 - m1)**2) +  &
        &           (K12**2*m1**2*c2(2))/(9.0_ki*(m0 - m1)**2) + B06*c2(9)
             endif
          endif
          if (verbosity.ge.2) then
             do ep=0,2
                write(iout,903) 'B0 (',cut2,',',-ep,') = ',J0(-ep)
                write(iout,903) 'B1 (',cut2,',',-ep,') =',J1(-ep)
                write(iout,903) 'B00(',cut2,',',-ep,') =',J00(-ep)
                write(iout,903) 'B11(',cut2,',',-ep,') =',J11(-ep)
             end do
          endif
       end if
       if (present(cache_flag)) cache_offset = cache_offset + 5

 903   format(a5,I2,a1,I2,a4,2(D24.15))

  end subroutine add2_rm

  subroutine add1_rm(nleg,c1,cut1,msq,tot1,scale2, &
        cache_flag, cache_offset, scalar_cache)
          use avh_olo, only: olo_a0
      implicit none 
      integer, intent(in) :: nleg, cut1
      complex(ki), dimension(0:4), intent(in) :: c1
      real(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(-2:0), intent(out) :: tot1
      real(ki), intent(in) :: scale2

      logical, intent(in), optional :: cache_flag
      integer, intent(inout), optional :: cache_offset
      complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache

      integer ::j1
      real(ki) :: m0
          complex(ki_avh), dimension(0:2) :: vala0
      complex(ki) :: ctmp
      integer :: ep, cache_index

      if (notfirsti.eqv.(.false.)) then
         if (isca .eq. 2) then
                call avh_olo_mu_set(real(sqrt(scale2),ki_avh))
         elseif (isca .eq. 4) then
!AC!             call setmudim(real(scale2, ki_lt))
         endif
         notfirsti=.true.
      endif

      j1=cut1

      m0=msq(j1)

  1   Format(A3,I2,A1,I2,A5,D24.15,A1,D24.15,A3)

      if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
      if      (isca.eq.1) then
          do ep=-2,0
             if (present(cache_flag)) then
                if (cache_flag) then
                   ctmp = scalar_cache(ep,cache_index)
                else
                   ctmp = qlI1(real(m0,ki_ql),real(scale2,ki_ql),ep)
                   scalar_cache(ep,cache_index) = ctmp
                end if
             else
                ctmp = qlI1(real(m0,ki_ql),real(scale2,ki_ql),ep)
             end if
             tot1(ep)=c1(0)*ctmp
             if (verbosity.ge.2) write(iout,1) &
         & 'I1(',cut1,',',ep,') = (',real(ctmp),',',aimag(ctmp),'  )' 
          enddo
!AC!      print*, "isca=1: QCDLoop not available"
!AC!      stop
      elseif  (isca.eq.2) then         
          if (present(cache_flag)) then
             if (cache_flag) then
                vala0(0) = scalar_cache( 0,cache_index)
                vala0(1) = scalar_cache(-1,cache_index)
                vala0(2) = scalar_cache(-2,cache_index)
             else
                call olo_a0(vala0,real(m0,ki_avh))
                scalar_cache( 0,cache_index) = vala0(0)
                scalar_cache(-1,cache_index) = vala0(1)
                scalar_cache(-2,cache_index) = vala0(2)
             end if
          else
             call olo_a0(vala0,real(m0,ki_avh))
          end if
          do ep=-2,0
             tot1(ep)= c1(0)*vala0(-ep) 
             if (verbosity.ge.2) write(iout,1) &
              & 'I1(',cut1,',',ep,') = (',&
              & real(vala0(-ep)),',',aimag(vala0(-ep)),'  )' 
          enddo
!AC!      print*, "isca=2: OneLOop not available"
!AC!      stop
      elseif (isca.eq.3) then
!AC!      if (verbosity.ge.2) write(iout,1) &
!AC!       & 'I1(',cut1,',',-2,') = (',0.0_ki,',',0.0_ki,'  )' 
!AC!      tot1(-2) = (0.0_ki,0.0_ki)
!AC!      if (present(cache_flag)) then
!AC!         if (cache_flag) then
!AC!            do ep=-1,0
!AC!               ctmp = scalar_cache(ep,cache_index)
!AC!               tot1(ep) = c1(0)*ctmp
!AC!               if (verbosity.ge.2) write(iout,1) &
!AC!                 & 'I1(',cut1,',',ep,') = (',&
!AC!                 & real(ctmp),',',aimag(ctmp),'  )' 
!AC!            enddo
!AC!         else
!AC!            scalar_cache(-2,cache_index) = czip
!AC!            call gtrunc_rm(1.0_ki, m0)
!AC!            do ep=-1,0
!AC!               ctmp = gA0(real(m0,ki_gol),&
!AC!              &        real(scale2,ki_gol),ep)
!AC!               scalar_cache(ep,cache_index) = ctmp
!AC!               tot1(ep) = c1(0)*ctmp
!AC!               if (verbosity.ge.2) write(iout,1) &
!AC!                 & 'I1(',cut1,',',ep,') = (',&
!AC!                 & real(ctmp),',',aimag(ctmp),'  )' 
!AC!            enddo
!AC!         end if
!AC!      else
!AC!         call gtrunc_rm(1.0_ki, m0)
!AC!         do ep=-1,0
!AC!            ctmp = gA0(real(m0,ki_gol),&
!AC!           &        real(scale2,ki_gol),ep)
!AC!            tot1(ep) = c1(0)*ctmp
!AC!            if (verbosity.ge.2) write(iout,1) &
!AC!              & 'I1(',cut1,',',ep,') = (',&
!AC!              & real(ctmp),',',aimag(ctmp),'  )' 
!AC!         enddo
!AC!      end if
          print*, "isca=3: Golem95 not available"
          stop
      elseif (isca.eq.4) then
!AC!      tot1(-2) = 0
!AC!      tot1(-1) = 0
!AC!      tot1( 0) = 0
!AC!      ep = -dim(0, int(getlambda()))
!AC!      if (present(cache_flag)) then
!AC!         if (cache_flag) then
!AC!            ctmp = scalar_cache(ep,cache_index)
!AC!         else
!AC!            ctmp = A0(real(m0,ki_lt))
!AC!            scalar_cache(ep,cache_index) = ctmp
!AC!         end if
!AC!      else
!AC!         ctmp = A0(real(m0,ki_lt))
!AC!      end if
!AC!      tot1(ep)=c1(0)*ctmp
!AC!      if (verbosity.ge.2) write(iout,1) &
!AC!     & 'I1(',cut1,',',ep,') = (',real(ctmp),',',aimag(ctmp),'  )' 
             print*, "isca=4: LoopTools not available"
             stop
      else
         print*, 'error in add1'
         stop
      endif

      if (present(cache_flag)) cache_offset = cache_offset + 1
  end subroutine add1_rm

  subroutine add4_cm(nleg,c4,cut4,Vi,msq,tot4,tot4r,scale2,&
        cache_flag, cache_offset, scalar_cache)
          use avh_olo, only: olo_d0
      implicit none
      integer, intent(in) :: cut4,nleg
      complex(ki), dimension(0:4), intent(in) :: c4
      real(ki), dimension(0:nleg-1,4), intent(in) :: Vi
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      real(ki), intent(in) :: scale2
      complex(ki), dimension(-2:0), intent(out) :: tot4
      complex(ki), intent(out) :: tot4r

      logical, intent(in), optional :: cache_flag
      integer, intent(inout), optional :: cache_offset
      complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache

      integer :: i,j1,j2,j3,j4

      real(ki) :: V1, V2, V3, V21, V31, V32
      complex(ki) :: m0, m1, m2, m3
      real(ki), dimension(4):: Vi1, Vi2, Vi3, Vi21, Vi31, Vi32
      complex(ki) :: c40
          complex(ki_avh), dimension(0:2) :: vald0
      complex(ki) :: ctmp
!AC!      complex(ki), dimension(-2:0) :: d0t
      integer :: ep, cache_index

      if (notfirsti.eqv.(.false.)) then
         if (isca .eq. 2) then
                call avh_olo_mu_set(real(sqrt(scale2),ki_avh))
         elseif (isca .eq. 4) then
!AC!             call setmudim(real(scale2, ki_lt))
         endif
         notfirsti=.true.
      endif

      j4=cut4/1000
      j3=(cut4-j4*1000)/100
      j2=(cut4-j4*1000-j3*100)/10
      j1=cut4-j4*1000-j3*100-j2*10

      m0=msq(j1)
      m1=msq(j2)
      m2=msq(j3)
      m3=msq(j4)

      if (allocated(s_mat)) then
         ! s_mat(i+1, j+1) = (Vi(i,:) - Vi(j,:))**2 - msq(i) - msq(j)
         V1  = s_mat(j2+1, j1+1) + msq(j2) + msq(j1)
         V2  = s_mat(j3+1, j1+1) + msq(j3) + msq(j1)
         V3  = s_mat(j4+1, j1+1) + msq(j4) + msq(j1)
         V21 = s_mat(j3+1, j2+1) + msq(j3) + msq(j2)
         V31 = s_mat(j4+1, j2+1) + msq(j4) + msq(j2)
         V32 = s_mat(j4+1, j3+1) + msq(j4) + msq(j3)
      else
         Vi1(:)=Vi(j2,:)-Vi(j1,:)
         Vi2(:)=Vi(j3,:)-Vi(j1,:)
         Vi3(:)=Vi(j1,:)-Vi(j4,:)
         Vi21(:)=Vi(j3,:)-Vi(j2,:)
         Vi31(:)=Vi(j4,:)-Vi(j2,:)
         Vi32(:)=Vi(j4,:)-Vi(j3,:)

         V1=sdot(Vi1,Vi1)
         V2=sdot(Vi2,Vi2)
         V3=sdot(Vi3,Vi3)
         V21=sdot(Vi21,Vi21)
         V31=sdot(Vi31,Vi31)
         V32=sdot(Vi32,Vi32)
      end if

      c40=c4(0)
      tot4r=-c4(4)/six

  1  Format(A3,I4,A1,I2,A5,D24.15,A1,D24.15,A3)

      if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
      if     (isca.eq.1) then
          print*, "isca=1: QCDLoop does not support complex masses."
!AC!      print*, "isca=1: QCDLoop not available"
         stop
      elseif (isca.eq.2) then
             if (present(cache_flag)) then
                if (cache_flag) then
                   vald0(0) = scalar_cache( 0,cache_index)
                   vald0(1) = scalar_cache(-1,cache_index)
                   vald0(2) = scalar_cache(-2,cache_index)
                else
                   call olo_d0(vald0,&
                  &      cmplx(V1,0.0_ki_avh,ki_avh), &
                  &      cmplx(V21,0.0_ki_avh,ki_avh), &
                  &      cmplx(V32,0.0_ki_avh,ki_avh), &
                  &      cmplx(V3,0.0_ki_avh,ki_avh), &
                  &      cmplx(V2,0.0_ki_avh,ki_avh), &
                  &      cmplx(V31,0.0_ki_avh,ki_avh), &
                  &      cmplx(real(m0,ki_avh),aimag(m0),ki_avh),&
                  &      cmplx(real(m1,ki_avh),aimag(m1),ki_avh),&
                  &      cmplx(real(m2,ki_avh),aimag(m2),ki_avh),&
                  &      cmplx(real(m3,ki_avh),aimag(m3),ki_avh))
                   scalar_cache( 0,cache_index) = vald0(0)
                   scalar_cache(-1,cache_index) = vald0(1)
                   scalar_cache(-2,cache_index) = vald0(2)
                end if
             else
                call olo_d0(vald0,&
               &      cmplx(V1,0.0_ki_avh,ki_avh), &
               &      cmplx(V21,0.0_ki_avh,ki_avh), &
               &      cmplx(V32,0.0_ki_avh,ki_avh), &
               &      cmplx(V3,0.0_ki_avh,ki_avh), &
               &      cmplx(V2,0.0_ki_avh,ki_avh), &
               &      cmplx(V31,0.0_ki_avh,ki_avh), &
               &      cmplx(real(m0,ki_avh),aimag(m0),ki_avh),&
               &      cmplx(real(m1,ki_avh),aimag(m1),ki_avh),&
               &      cmplx(real(m2,ki_avh),aimag(m2),ki_avh),&
               &      cmplx(real(m3,ki_avh),aimag(m3),ki_avh))
             end if
             do ep=-2,0
                tot4(ep)= c40*vald0(-ep) 
                if (verbosity.ge.2) write(iout,1) &
                & 'I4(',cut4,',',ep,') = (',real(vald0(-ep)),',',&
                & aimag(vald0(-ep)),'  )' 
             enddo
!AC!      print*, "isca=2: OneLOop not available"
!AC!      stop
      elseif (isca.eq.3) then
!AC!        call gtrunc_rm(abs(V32)+abs(V31), V1,V2,V3,V21,V32,V31)
!AC!        call gtrunc_cm(abs(V32)+abs(V31), m0,m1,m2,m3)
!AC!        do ep=-2,0
!AC!           if (present(cache_flag)) then
!AC!              if (cache_flag) then
!AC!                 d0t(ep) = scalar_cache(ep,cache_index)
!AC!              else
!AC!                 d0t(ep)=gD0C(&
!AC!                &        cmplx(V1,0.0_ki_gol,ki_gol), &
!AC!                &        cmplx(V21,0.0_ki_gol,ki_gol), &
!AC!                &        cmplx(V32,0.0_ki_gol,ki_gol), &
!AC!                &        cmplx(V3,0.0_ki_gol,ki_gol), &
!AC!                &        cmplx(V2,0.0_ki_gol,ki_gol), &
!AC!                &        cmplx(V31,0.0_ki_gol,ki_gol), &
!AC!                &        cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!                &        cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!                &        cmplx(real(m2,ki_gol),aimag(m2),ki_gol),&
!AC!                &        cmplx(real(m3,ki_gol),aimag(m3),ki_gol),&
!AC!                &        real(scale2,ki_gol),ep)
!AC!                 scalar_cache(ep,cache_index) = d0t(ep)
!AC!              end if
!AC!           else
!AC!              d0t(ep)=gD0C(&
!AC!             &        cmplx(V1,0.0_ki_gol,ki_gol), &
!AC!             &        cmplx(V21,0.0_ki_gol,ki_gol), &
!AC!             &        cmplx(V32,0.0_ki_gol,ki_gol), &
!AC!             &        cmplx(V3,0.0_ki_gol,ki_gol), &
!AC!             &        cmplx(V2,0.0_ki_gol,ki_gol), &
!AC!             &        cmplx(V31,0.0_ki_gol,ki_gol), &
!AC!             &        cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!             &        cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!             &        cmplx(real(m2,ki_gol),aimag(m2),ki_gol),&
!AC!             &        cmplx(real(m3,ki_gol),aimag(m3),ki_gol),&
!AC!             &        real(scale2,ki_gol),ep)
!AC!           end if
!AC!        end do
!AC!        !d0t( 0) = d0t(0) + log(scale2) * (d0t(-1) &
!AC!        !      & + 0.5_ki * log(scale2) * d0t(-2))
!AC!        !d0t(-1) = d0t(-1) + log(scale2) * d0t(-2)
!AC!        if (verbosity.ge.2) then
!AC!           do ep=-2,0
!AC!              write(iout,1) 'I4(',cut4,',',ep,&
!AC!              & ') = (',real(d0t(ep)),',',aimag(d0t(ep)),'  )' 
!AC!           end do
!AC!        end if
!AC!        tot4(:) = d0t(:) * c40
            print*, "isca=3: Golem95 not available"
            stop
      elseif (isca.eq.4) then
!AC!        tot4(-2) = 0
!AC!        tot4(-1) = 0
!AC!        tot4(0) = 0
!AC!        ep = -dim(0, int(getlambda()))
!AC!        if (present(cache_flag)) then
!AC!           if (cache_flag) then
!AC!              ctmp = scalar_cache(ep,cache_index)
!AC!           else
!AC!              call gtrunc_rm(abs(V32)+abs(V31), V1,V2,V3,V21,V32,V31)
!AC!              call gtrunc_cm(abs(V32)+abs(V31), m0,m1,m2,m3)
!AC!              ctmp=D0C(&
!AC!             &       cmplx(V1,0.0_ki_lt,ki_lt), &
!AC!             &       cmplx(V21,0.0_ki_lt,ki_lt), &
!AC!             &       cmplx(V32,0.0_ki_lt,ki_lt), &
!AC!             &       cmplx(V3,0.0_ki_lt,ki_lt), &
!AC!             &       cmplx(V2,0.0_ki_lt,ki_lt), &
!AC!             &       cmplx(V31,0.0_ki_lt,ki_lt), &
!AC!             &       cmplx(real(m0,ki_lt),aimag(m0),ki_lt),&
!AC!             &       cmplx(real(m1,ki_lt),aimag(m1),ki_lt),&
!AC!             &       cmplx(real(m2,ki_lt),aimag(m2),ki_lt),&
!AC!             &       cmplx(real(m3,ki_lt),aimag(m3),ki_lt))
!AC!              scalar_cache(ep,cache_index) = ctmp
!AC!           end if
!AC!        else
!AC!           call gtrunc_rm(abs(V32)+abs(V31), V1,V2,V3,V21,V32,V31)
!AC!           call gtrunc_cm(abs(V32)+abs(V31), m0,m1,m2,m3)
!AC!           ctmp=D0C(&
!AC!          &       cmplx(V1,0.0_ki_lt,ki_lt), &
!AC!          &       cmplx(V21,0.0_ki_lt,ki_lt), &
!AC!          &       cmplx(V32,0.0_ki_lt,ki_lt), &
!AC!          &       cmplx(V3,0.0_ki_lt,ki_lt), &
!AC!          &       cmplx(V2,0.0_ki_lt,ki_lt), &
!AC!          &       cmplx(V31,0.0_ki_lt,ki_lt), &
!AC!          &       cmplx(real(m0,ki_lt),aimag(m0),ki_lt),&
!AC!          &       cmplx(real(m1,ki_lt),aimag(m1),ki_lt),&
!AC!          &       cmplx(real(m2,ki_lt),aimag(m2),ki_lt),&
!AC!          &       cmplx(real(m3,ki_lt),aimag(m3),ki_lt))
!AC!        end if
!AC!        tot4(ep)=c40*ctmp
!AC!        if (verbosity.ge.2) write(iout,1) &
!AC!     & 'I4(',cut4,',',ep,') = (',real(ctmp),',',aimag(ctmp),'  )' 
             print*, "isca=4: LoopTools not available"
             stop
      else
         print*, 'error in add4'
         stop
      endif
      if (present(cache_flag)) cache_offset = cache_offset + 1
      tot4(0)=tot4(0) + tot4r
  end subroutine add4_cm

  subroutine add3_cm(nleg,c3,cut3,Vi,msq,tot3,tot3r,scale2,&
        cache_flag, cache_offset, scalar_cache)
          use avh_olo, only: olo_c0
      implicit none

      integer, intent(in) :: nleg, cut3
      complex(ki), dimension(0:9), intent(in) :: c3
      real(ki), dimension(0:nleg-1,4) ::Vi
      complex(ki), dimension(0:nleg-1):: msq
      complex(ki), dimension(-2:0), intent(out) :: tot3
      complex(ki), intent(out) :: tot3r
      real(ki), intent(in) :: scale2

      logical, intent(in), optional :: cache_flag
      integer, intent(inout), optional :: cache_offset
      complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache

      integer :: i,j1,j2,j3
      real(ki) :: V1, V2, V3
      complex(ki) :: m0, m1, m2
      real(ki), dimension(4):: Vi1, Vi2, Vi3
      complex(ki) :: c30
          complex(ki_avh), dimension(0:2) :: valc0
      complex(ki) :: ctmp
!AC!      complex(ki), dimension(-2:0) :: c0t
      integer :: ep, cache_index

      if (notfirsti.eqv.(.false.)) then
         if (isca .eq. 2) then
                call avh_olo_mu_set(real(sqrt(scale2),ki_avh))
         elseif (isca .eq. 4) then
!AC!             call setmudim(real(scale2, ki_lt))
         endif
         notfirsti=.true.
      endif

      j3=cut3/100
      j2=(cut3-j3*100)/10
      j1=cut3-j3*100-j2*10

      m0=msq(j1)
      m1=msq(j2)
      m2=msq(j3)

      if (allocated(s_mat)) then
         ! s_mat(i+1, j+1) = (Vi(i,:) - Vi(j,:))**2 - msq(i) - msq(j)
         V1 = s_mat(j2+1, j1+1) + msq(j2) + msq(j1)
         V2 = s_mat(j3+1, j2+1) + msq(j3) + msq(j2)
         V3 = s_mat(j3+1, j1+1) + msq(j3) + msq(j1)
      else
         Vi1(:)=Vi(j2,:)-Vi(j1,:)
         Vi2(:)=Vi(j3,:)-Vi(j2,:)
         Vi3(:)=Vi(j1,:)-Vi(j3,:)

         V1=sdot(Vi1,Vi1)
         V2=sdot(Vi2,Vi2)
         V3=sdot(Vi3,Vi3)
      end if

      c30=c3(0)

      tot3r=+c3(7)/two

  1   Format(A3,I3,A1,I2,A5,D24.15,A1,D24.15,A3)

      if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
      if (isca.eq.1) then
          print*, "isca=1: QCDLoop does not support complex masses."
!AC!      print*, "isca=1: QCDLoop not available"
         stop
      elseif  (isca.eq.2) then
              if (present(cache_flag)) then
                 if (cache_flag) then
                    valc0(0) = scalar_cache( 0,cache_index)
                    valc0(1) = scalar_cache(-1,cache_index)
                    valc0(2) = scalar_cache(-2,cache_index)
                 else
                    call olo_c0(valc0,&
                   &      cmplx(V1,0.0_ki_avh,ki_avh), &
                   &      cmplx(V2,0.0_ki_avh,ki_avh), &
                   &      cmplx(V3,0.0_ki_avh,ki_avh), &
                   &      cmplx(real(m0,ki_avh),aimag(m0),ki_avh),&
                   &      cmplx(real(m1,ki_avh),aimag(m1),ki_avh),&
                   &      cmplx(real(m2,ki_avh),aimag(m2),ki_avh))
                    scalar_cache( 0,cache_index) = valc0(0)
                    scalar_cache(-1,cache_index) = valc0(1)
                    scalar_cache(-2,cache_index) = valc0(2)
                 end if
              else
                 call olo_c0(valc0,&
                &      cmplx(V1,0.0_ki_avh,ki_avh), &
                &      cmplx(V2,0.0_ki_avh,ki_avh), &
                &      cmplx(V3,0.0_ki_avh,ki_avh), &
                &      cmplx(real(m0,ki_avh),aimag(m0),ki_avh),&
                &      cmplx(real(m1,ki_avh),aimag(m1),ki_avh),&
                &      cmplx(real(m2,ki_avh),aimag(m2),ki_avh))
              end if
              do ep=-2,0
                 tot3(ep)= c30*valc0(-ep) 
                 if (verbosity.ge.2) write(iout,1) &
         &'I3(',cut3,',',ep,') = (',real(valc0(-ep)),',',aimag(valc0(-ep)),'  )' 
              enddo
!AC!      print*, "isca=2: OneLOop not available"
!AC!      stop
      elseif (isca.eq.3) then
!AC!        call gtrunc_rm(abs(V1)+abs(V2)+abs(V3), V1,V2,V3)
!AC!        call gtrunc_cm(abs(V1)+abs(V2)+abs(V3), m0,m1,m2)
!AC!        do ep=-2,0
!AC!           if (present(cache_flag)) then
!AC!              if (cache_flag) then
!AC!                 c0t(ep) = scalar_cache(ep,cache_index)
!AC!              else
!AC!                 c0t(ep)=gC0C(&
!AC!                &        cmplx(V1,0.0_ki_gol,ki_gol), &
!AC!                &        cmplx(V2,0.0_ki_gol,ki_gol), &
!AC!                &        cmplx(V3,0.0_ki_gol,ki_gol), &
!AC!                &        cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!                &        cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!                &        cmplx(real(m2,ki_gol),aimag(m2),ki_gol),&
!AC!                &        real(scale2,ki_gol),ep)
!AC!                 scalar_cache(ep,cache_index) = c0t(ep)
!AC!              end if
!AC!           else
!AC!              c0t(ep)=gC0C(&
!AC!             &        cmplx(V1,0.0_ki_gol,ki_gol), &
!AC!             &        cmplx(V2,0.0_ki_gol,ki_gol), &
!AC!             &        cmplx(V3,0.0_ki_gol,ki_gol), &
!AC!             &        cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!             &        cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!             &        cmplx(real(m2,ki_gol),aimag(m2),ki_gol),&
!AC!             &        real(scale2,ki_gol),ep)
!AC!           end if
!AC!        end do
!AC!        !c0t( 0) = c0t(0) + log(scale2) * (c0t(-1) &
!AC!        !      & + 0.5_ki * log(scale2) * c0t(-2))
!AC!        !c0t(-1) = c0t(-1) + log(scale2) * c0t(-2)
!AC!        do ep=-2,0
!AC!           if (verbosity.ge.2) write(iout,1) 'I3(',cut3,',',ep,&
!AC!             & ') = (',real(c0t(ep)),',',aimag(c0t(ep)),'  )' 
!AC!        end do
!AC!        tot3(:) = c0t(:) * c30
            print*, "isca=3: Golem95 not available"
            stop
      elseif (isca.eq.4) then
!AC!        tot3(-2) = 0
!AC!        tot3(-1) = 0
!AC!        tot3(0) = 0
!AC!        call gtrunc_rm(abs(V1)+abs(V2)+abs(V3), V1,V2,V3)
!AC!        call gtrunc_cm(abs(V1)+abs(V2)+abs(V3), m0,m1,m2)
!AC!        ep = -dim(0, int(getlambda()))
!AC!        if (present(cache_flag)) then
!AC!           if (cache_flag) then
!AC!              ctmp = scalar_cache(ep,cache_index)
!AC!           else
!AC!              ctmp=C0C(&
!AC!             &       cmplx(V1,0.0_ki_lt,ki_lt), &
!AC!             &       cmplx(V2,0.0_ki_lt,ki_lt), &
!AC!             &       cmplx(V3,0.0_ki_lt,ki_lt), &
!AC!             &       cmplx(real(m0,ki_lt),aimag(m0),ki_lt),&
!AC!             &       cmplx(real(m1,ki_lt),aimag(m1),ki_lt),&
!AC!             &       cmplx(real(m2,ki_lt),aimag(m2),ki_lt))
!AC!              scalar_cache(ep,cache_index) = ctmp
!AC!           end if
!AC!        else
!AC!           ctmp=C0C(&
!AC!          &       cmplx(V1,0.0_ki_lt,ki_lt), &
!AC!          &       cmplx(V2,0.0_ki_lt,ki_lt), &
!AC!          &       cmplx(V3,0.0_ki_lt,ki_lt), &
!AC!          &       cmplx(real(m0,ki_lt),aimag(m0),ki_lt),&
!AC!          &       cmplx(real(m1,ki_lt),aimag(m1),ki_lt),&
!AC!          &       cmplx(real(m2,ki_lt),aimag(m2),ki_lt))
!AC!        end if
!AC!        tot3(ep)=c30*ctmp 
!AC!        if (verbosity.ge.2) write(iout,1) &
!AC!       &'I3(',cut3,',',ep,') = (',real(ctmp),',',aimag(ctmp),'  )'
             print*, "isca=4: LoopTools not available"
             stop
      else
         print*, 'error in add3'
         stop
      endif
      tot3(0)=tot3(0) + tot3r
      if (present(cache_flag)) cache_offset = cache_offset + 1
  end subroutine add3_cm

  subroutine add2_cm(nleg,c2,cut2,k1,k2,msq,tot2,tot2r,scale2,&
        cache_flag, cache_offset, scalar_cache)
          use avh_olo, only: olo_b11
      implicit none

      integer, intent(in) :: nleg, cut2
      complex(ki), dimension(0:9), intent(in) :: c2
      real(ki), dimension(4), intent(in) :: k1, k2
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(-2:0), intent(out) :: tot2
      complex(ki), intent(out) :: tot2r
      real(ki), intent(in) :: scale2

      logical, intent(in), optional :: cache_flag
      integer, intent(inout), optional :: cache_offset
      complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache

      integer :: i1,i2
      real(ki) :: K11, K12
      complex(ki) :: B06
      complex(ki) :: m0, m1
      integer :: ep, cache_index

      complex(ki), dimension(-2:0) :: J0, J1, J00, J01, J11 
          complex(ki_avh), dimension(0:2) :: scf2, scf1, scf0, scf 
          complex(ki) :: xbb, xb0, xb00
          real(ki) :: bkv
          integer :: i

      if (notfirsti.eqv.(.false.)) then
         if (isca .eq. 2) then
                call avh_olo_mu_set(real(sqrt(scale2),ki_avh))
         elseif (isca .eq. 4) then
!AC!             call setmudim(real(scale2, ki_lt))
         endif
         notfirsti=.true.
      endif

      i2=cut2/10
      i1=cut2-i2*10

      m0=msq(i1)
      m1=msq(i2)

      if (allocated(s_mat)) then
         K11 = s_mat(i2+1, i1+1) + msq(i1) + msq(i2)
      else
         K11 = sdot(K1,K1)
      end if
      K12=sdot(K1,K2)

      B06=-(K11-three*(m0+m1))/six

      tot2r= + B06*c2(9)

      if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
      if (isca.eq.1) then
          print*, "isca=1: QCDLoop does not support complex masses."
!AC!      print*, "isca=1: QCDLoop not available"
          stop
       elseif  (isca.eq.2) then
             xbb=c2(0)
             xb0=c2(1)
             xb00=c2(2)
             bkv=K12
             if (present(cache_flag)) then
                if (cache_flag) then
                   scf(:)  = scalar_cache(:,cache_index+0)
                   scf0(:) = scalar_cache(:,cache_index+1)
                   scf1(:) = scalar_cache(:,cache_index+2)
                   scf2(:) = scalar_cache(:,cache_index+3)
                else
                   call olo_b11(scf2,scf0,scf1,scf,&
                  &             cmplx(K11,0.0_ki_avh,ki_avh), &
                  &             cmplx(real(m0,ki_avh),aimag(m0),ki_avh),&
                  &             cmplx(real(m1,ki_avh),aimag(m1),ki_avh))
                   scalar_cache(:,cache_index+0) = scf(:)
                   scalar_cache(:,cache_index+1) = scf0(:)
                   scalar_cache(:,cache_index+2) = scf1(:)
                   scalar_cache(:,cache_index+3) = scf2(:)
                end if
             else
                call olo_b11(scf2,scf0,scf1,scf,&
               &             cmplx(K11,0.0_ki_avh,ki_avh), &
               &             cmplx(real(m0,ki_avh),aimag(m0),ki_avh),&
               &             cmplx(real(m1,ki_avh),aimag(m1),ki_avh))
             end if
             tot2(0)=xbb*scf(0)+xb0*bkv*scf1(0)+xb00*bkv*bkv*scf2(0)
             tot2(0)=tot2(0)+ B06*c2(9)
             tot2(-1)=xbb*scf(1)+xb0*bkv*scf1(1)+xb00*bkv*bkv*scf2(1)
             tot2(-2)=xbb*scf(2)+xb0*bkv*scf1(2)+xb00*bkv*bkv*scf2(2)
              if (verbosity.ge.2) then
                 do i=0,2
                    write(iout,903) 'B_0 (',cut2,',',-i,') = ',scf(i)
                    write(iout,903) 'B_1 (',cut2,',',-i,') =',scf1(i)
                    write(iout,903) 'B_11(',cut2,',',-i,') =',scf2(i)
                 enddo
              endif
!AC!      print*, "isca=2: OneLOop not available"
!AC!      stop
      elseif (isca.eq.3) then
!AC!        if (present(cache_flag)) then
!AC!           if (cache_flag) then
!AC!              J0(:)  = scalar_cache(:,cache_index+0)
!AC!              J1(:)  = scalar_cache(:,cache_index+1)
!AC!              J01(:) = scalar_cache(:,cache_index+2)
!AC!              J11(:) = scalar_cache(:,cache_index+3)
!AC!              J00(:) = scalar_cache(:,cache_index+4)
!AC!           else
!AC!              call gtrunc_rm(abs(K11)+1.0_ki, K11)
!AC!              call gtrunc_cm(abs(K11)+1.0_ki, m0,m1)
!AC!              do ep=-2,0
!AC!                 J00(ep)= gB0C(&
!AC!                &   cmplx(K11,0.0_ki_gol,ki_gol), &
!AC!                &   cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!                &   cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!                &   real(scale2,ki_gol),ep)
!AC!                 J11(ep)= gB0C(&
!AC!                &   cmplx(K11,0.0_ki_gol,ki_gol), &
!AC!                &   cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!                &   cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!                &   real(scale2,ki_gol),ep)
!AC!                 J01(ep)= gB0C(&
!AC!                &   cmplx(K11,0.0_ki_gol,ki_gol), &
!AC!                &   cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!                &   cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!                &   real(scale2,ki_gol),ep)
!AC!                 J0(ep) = gA0C(&
!AC!                &   cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!                &   real(scale2,ki_gol),ep)
!AC!                 J1(ep) = gA0C(&
!AC!                &   cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!                &   real(scale2,ki_gol),ep)
!AC!              end do
!AC!              scalar_cache(:,cache_index+0) = J0(:)
!AC!              scalar_cache(:,cache_index+1) = J1(:)
!AC!              scalar_cache(:,cache_index+2) = J01(:)
!AC!              scalar_cache(:,cache_index+3) = J11(:)
!AC!              scalar_cache(:,cache_index+4) = J00(:)
!AC!           end if
!AC!        else
!AC!           call gtrunc_rm(abs(K11)+1.0_ki, K11)
!AC!           call gtrunc_cm(abs(K11)+1.0_ki, m0,m1)
!AC!           do ep=-2,0
!AC!              J00(ep)= gB0C(&
!AC!             &   cmplx(K11,0.0_ki_gol,ki_gol), &
!AC!             &   cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!             &   cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!             &   real(scale2,ki_gol),ep)
!AC!              J11(ep)= gB0C(&
!AC!             &   cmplx(K11,0.0_ki_gol,ki_gol), &
!AC!             &   cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!             &   cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!             &   real(scale2,ki_gol),ep)
!AC!              J01(ep)= gB0C(&
!AC!             &   cmplx(K11,0.0_ki_gol,ki_gol), &
!AC!             &   cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!             &   cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!             &   real(scale2,ki_gol),ep)
!AC!              J0(ep) = gA0C(&
!AC!             &   cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!             &   real(scale2,ki_gol),ep)
!AC!              J1(ep) = gA0C(&
!AC!             &   cmplx(real(m1,ki_gol),aimag(m1),ki_gol),&
!AC!             &   real(scale2,ki_gol),ep)
!AC!           end do
!AC!        end if
            print*, "isca=3: Golem95 not available"
            stop
      elseif (isca.eq.4) then
!AC!        J00(:) = 0.0_ki_lt
!AC!        J11(:) = 0.0_ki_lt
!AC!        J0(:) = 0.0_ki_lt
!AC!        J1(:) = 0.0_ki_lt
!AC!        ep = -dim(0, int(getlambda()))
!AC!        if (present(cache_flag)) then
!AC!           if (cache_flag) then
!AC!              J0(:)  = scalar_cache(:,cache_index+0)
!AC!              J1(:)  = scalar_cache(:,cache_index+1)
!AC!              J01(:) = scalar_cache(:,cache_index+2)
!AC!              J11(:) = scalar_cache(:,cache_index+3)
!AC!              J00(:) = scalar_cache(:,cache_index+4)
!AC!           else
!AC!              call gtrunc_rm(abs(K11)+1.0_ki, K11)
!AC!              call gtrunc_cm(abs(K11)+1.0_ki, m0,m1)
!AC!              J00(ep) = B0C(&
!AC!             &   cmplx(K11,0.0_ki_lt,ki_lt), &
!AC!             &   cmplx(real(m0,ki_lt),aimag(m0),ki_lt),&
!AC!             &   cmplx(real(m0,ki_lt),aimag(m0),ki_lt))
!AC!              J11(ep) = B0C(&
!AC!             &   cmplx(K11,0.0_ki_lt,ki_lt), &
!AC!             &   cmplx(real(m1,ki_lt),aimag(m1),ki_lt),&
!AC!             &   cmplx(real(m1,ki_lt),aimag(m1),ki_lt))
!AC!              J01(ep) = B0C(&
!AC!             &   cmplx(K11,0.0_ki_lt,ki_lt), &
!AC!             &   cmplx(real(m0,ki_lt),aimag(m0),ki_lt),&
!AC!             &   cmplx(real(m1,ki_lt),aimag(m1),ki_lt))
!AC!              J0(ep)  = A0C(&
!AC!             &   cmplx(real(m0,ki_lt),aimag(m0),ki_lt))
!AC!              J1(ep)  = A0C(&
!AC!             &   cmplx(real(m1,ki_lt),aimag(m1),ki_lt))
!AC!              scalar_cache(:,cache_index+0) = J0(:)
!AC!              scalar_cache(:,cache_index+1) = J1(:)
!AC!              scalar_cache(:,cache_index+2) = J01(:)
!AC!              scalar_cache(:,cache_index+3) = J11(:)
!AC!              scalar_cache(:,cache_index+4) = J00(:)
!AC!           end if
!AC!        else
!AC!           call gtrunc_rm(abs(K11)+1.0_ki, K11)
!AC!           call gtrunc_cm(abs(K11)+1.0_ki, m0,m1)
!AC!           J00(ep) = B0C(&
!AC!          &           cmplx(K11,0.0_ki_lt,ki_lt), &
!AC!          &           cmplx(real(m0,ki_lt),aimag(m0),ki_lt),&
!AC!          &           cmplx(real(m0,ki_lt),aimag(m0),ki_lt))
!AC!           J11(ep) = B0C(&
!AC!          &           cmplx(K11,0.0_ki_lt,ki_lt), &
!AC!          &           cmplx(real(m1,ki_lt),aimag(m1),ki_lt),&
!AC!          &           cmplx(real(m1,ki_lt),aimag(m1),ki_lt))
!AC!           J01(ep) = B0C(&
!AC!          &           cmplx(K11,0.0_ki_lt,ki_lt), &
!AC!          &           cmplx(real(m0,ki_lt),aimag(m0),ki_lt),&
!AC!          &           cmplx(real(m1,ki_lt),aimag(m1),ki_lt))
!AC!           J0(ep)  = A0C(&
!AC!          &           cmplx(real(m0,ki_lt),aimag(m0),ki_lt))
!AC!           J1(ep)  = A0C(&
!AC!          &           cmplx(real(m1,ki_lt),aimag(m1),ki_lt))
!AC!        end if
            print*, "isca=4: LoopTools not available"
            stop
       else
          print*, 'error in add2'
          stop
       endif

       if (isca.eq.1 .or. isca.eq.3 .or. isca.eq.4) then
          if (abs(K11).gt.zip1) then
             do ep=-2,0
                tot2(ep)=-(K12*(two*K12*(m0 - m1)*c2(2) +  &
        &         K11*(-three*c2(1) + two*K12*c2(2)))*J0(ep))/(six*K11**2) &
        &         + ((two*K12**2*(m0 - m1)**2*c2(2) +  &
        &         K11*K12*(-three*m0*c2(1) + three*m1*c2(1) +  &
        &         two*K12*m0*c2(2) - four*K12*m1*c2(2)) +  &
        &         K11**2*(6*c2(0) + K12*(-three*c2(1) + two*K12*c2(2))))* &
        &         J01(ep))/(six*K11**2) +  &
        &         (K12*(two*K12*(m0 - m1)*c2(2) +  &
        &         K11*(-three*c2(1) + four*K12*c2(2)))*J1(ep))/(six*K11**2)
             enddo
             tot2(0)=tot2(0)+(K12**2*c2(2))/18.0_ki &
        &       - (K12**2*m0*c2(2))/(six*K11) -  &
        &       (K12**2*m1*c2(2))/(six*K11) + B06*c2(9)
          else
             if (m1.eq.m0) then
                do ep=-2,0
                   tot2(ep)=(c2(0) + (K12*(-three*c2(1) &
        &                    + two*K12*c2(2)))/six)*J00(ep)
                enddo
                tot2(0)=tot2(0)+B06*c2(9)
             else
                do ep=-2,0
                   tot2(ep)=(K12*m0**2*(-three*m0*c2(1) + three*m1*c2(1) &
        &              + two*K12*m0*c2(2))* &
        &              J00(ep))/(six*(m0 - m1)**3) + c2(0)*J01(ep) -  &
        &              (K12*m1*(m0*m1*(9*c2(1) - 6*K12*c2(2)) -  &
        &              6*m0**2*(c2(1) - K12*c2(2)) +  &
        &              m1**2*(-three*c2(1) + two*K12*c2(2)))*J11(ep))/ &
        &              (six*(m0 - m1)**3)
                enddo
                tot2(0)=tot2(0) + (-three*K12*m0*c2(1))/(four*(m0 - m1)) +  &
        &           (K12*m1*c2(1))/(four*m0 - four*m1) +  &
        &           (11.0_ki*K12**2*m0**2*c2(2))/(18.0_ki*(m0 - m1)**2) -  &
        &           (7.0_ki*K12**2*m0*m1*c2(2))/(18.0_ki*(m0 - m1)**2) +  &
        &           (K12**2*m1**2*c2(2))/(9.0_ki*(m0 - m1)**2) + B06*c2(9)
             endif
          endif
          if (verbosity.ge.2) then
             do ep=0,2
                write(iout,903) 'B0 (',cut2,',',-ep,') = ',J0(-ep)
                write(iout,903) 'B1 (',cut2,',',-ep,') =',J1(-ep)
                write(iout,903) 'B00(',cut2,',',-ep,') =',J00(-ep)
                write(iout,903) 'B11(',cut2,',',-ep,') =',J11(-ep)
             end do
          endif
       end if
       if (present(cache_flag)) cache_offset = cache_offset + 5

 903   format(a5,I2,a1,I2,a4,2(D24.15))
  end subroutine add2_cm

  subroutine add1_cm(nleg,c1,cut1,msq,tot1,scale2,&
        cache_flag, cache_offset, scalar_cache)
          use avh_olo, only: olo_a0
      implicit none
      integer, intent(in) :: nleg, cut1
      complex(ki), dimension(0:4), intent(in) :: c1
      complex(ki), dimension(0:nleg-1), intent(in) :: msq
      complex(ki), dimension(-2:0), intent(out) :: tot1
      real(ki), intent(in) :: scale2

      logical, intent(in), optional :: cache_flag
      integer, intent(inout), optional :: cache_offset
      complex(ki), intent(inout), optional, dimension(-2:0,*) :: scalar_cache


      integer ::j1
      complex(ki) :: m0
          complex(ki_avh), dimension(0:2) :: vala0
      complex(ki) :: ctmp
      integer :: ep, cache_index

      if (notfirsti.eqv.(.false.)) then
         if (isca .eq. 2) then
                call avh_olo_mu_set(real(sqrt(scale2),ki_avh))
         elseif (isca .eq. 4) then
!AC!             call setmudim(real(scale2, ki_lt))
         endif
         notfirsti=.true.
      endif

      j1=cut1

      m0=msq(j1)

  1   Format(A3,I2,A1,I2,A5,D24.15,A1,D24.15,A3)

      if (present(cache_flag)) cache_index = lbound(scalar_cache,2)+cache_offset
      if      (isca.eq.1) then
          print*, "isca=1: QCDLoop does not support complex masses."
!AC!      print*, "isca=1: QCDLoop not available"
         stop
      elseif  (isca.eq.2) then         
          if (present(cache_flag)) then
             if (cache_flag) then
                vala0(0) = scalar_cache( 0,cache_index)
                vala0(1) = scalar_cache(-1,cache_index)
                vala0(2) = scalar_cache(-2,cache_index)
             else
                call olo_a0(vala0,&
               &     cmplx(real(m0,ki_avh),aimag(m0),ki_avh))
                scalar_cache( 0,cache_index) = vala0(0)
                scalar_cache(-1,cache_index) = vala0(1)
                scalar_cache(-2,cache_index) = vala0(2)
             end if
          else
             call olo_a0(vala0,&
               & cmplx(real(m0,ki_avh),aimag(m0),ki_avh))
          end if
          do ep=-2,0
             tot1(ep)= c1(0)*vala0(-ep) 
             if (verbosity.ge.2) write(iout,1) &
            & 'I1(',cut1,',',ep,') = (',real(vala0(-ep)),',',&
            & aimag(vala0(-ep)),'  )' 
          enddo
!AC!      print*, "isca=2: OneLOop not available"
!AC!      stop
      elseif (isca.eq.3) then
!AC!      tot1(-2) = (0.0_ki,0.0_ki)
!AC!      if (verbosity.ge.2) write(iout,1) &
!AC!     & 'I1(',cut1,',',-2,') = (',0.0_ki,',',0.0_ki,'  )' 
!AC!      if (present(cache_flag)) then
!AC!         if (cache_flag) then
!AC!            do ep=-1,0
!AC!               ctmp = scalar_cache(ep,cache_index)
!AC!               tot1(ep) = c1(0)*ctmp
!AC!               if (verbosity.ge.2) write(iout,1) &
!AC!                 & 'I1(',cut1,',',ep,') = (',&
!AC!                 & real(ctmp),',',aimag(ctmp),'  )' 
!AC!            enddo
!AC!         else
!AC!            scalar_cache(-2,cache_index) = czip
!AC!            call gtrunc_cm(1.0_ki, m0)
!AC!            do ep=-1,0
!AC!               ctmp = gA0C(&
!AC!              & cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!              &        real(scale2,ki_gol),ep)
!AC!               scalar_cache(ep,cache_index) = ctmp
!AC!               tot1(ep) = c1(0)*ctmp
!AC!               if (verbosity.ge.2) write(iout,1) &
!AC!                 & 'I1(',cut1,',',ep,') = (',&
!AC!                 & real(ctmp),',',aimag(ctmp),'  )' 
!AC!            enddo
!AC!         end if
!AC!      else
!AC!         call gtrunc_cm(1.0_ki, m0)
!AC!         do ep=-1,0
!AC!            ctmp = gA0C(&
!AC!           & cmplx(real(m0,ki_gol),aimag(m0),ki_gol),&
!AC!           &        real(scale2,ki_gol),ep)
!AC!            tot1(ep) = c1(0)*ctmp
!AC!            if (verbosity.ge.2) write(iout,1) &
!AC!              & 'I1(',cut1,',',ep,') = (',&
!AC!              & real(ctmp),',',aimag(ctmp),'  )' 
!AC!         enddo
!AC!      end if
          print*, "isca=3: Golem95 not available"
          stop
      elseif (isca.eq.4) then
!AC!      tot1(-2) = 0
!AC!      tot1(-1) = 0
!AC!      tot1(0) = 0
!AC!      ep = -dim(0, int(getlambda()))
!AC!      if (present(cache_flag)) then
!AC!         if (cache_flag) then
!AC!            ctmp = scalar_cache(ep,cache_index)
!AC!         else
!AC!            ctmp = A0C(cmplx(real(m0,ki_lt),aimag(m0),ki_lt))
!AC!            scalar_cache(ep,cache_index) = ctmp
!AC!         end if
!AC!      else
!AC!         ctmp = A0C(cmplx(real(m0,ki_lt),aimag(m0),ki_lt))
!AC!      end if
!AC!      tot1(ep)=c1(0)*ctmp
!AC!      if (verbosity.ge.2) write(iout,1) &
!AC!     & 'I1(',cut1,',',ep,') = (',real(ctmp),',',aimag(ctmp),'  )' 
          print*, "isca=4: LoopTools not available"
          stop
      else
         print*, 'error in add1'
         stop
      endif
      if (present(cache_flag)) cache_offset = cache_offset + 1
  end subroutine add1_cm

!AC!function gA0(m0,mu2,ep)
!AC!   implicit none
!AC!   real(ki_gol), intent(in) :: m0, mu2
!AC!   integer, intent(in) :: ep
!AC!   complex(ki_gol) :: gA0
!AC!   if(ep.eq.(-2) .or. m0.eq.0.0_ki_gol) then
!AC!      gA0 = (0.0_ki_gol, 0.0_ki_gol)
!AC!   elseif(ep.eq.(-1)) then
!AC!      gA0 = m0 * gB0(0.0_ki_gol,m0,m0,mu2,-1)
!AC!   else
!AC!      gA0 = m0 * (gB0(0.0_ki_gol,m0,m0,mu2,0) &
!AC!          &    +  gB0(0.0_ki_gol,m0,m0,mu2,-1))
!AC!   end if
!AC!end function gA0
!AC!function gA0C(m0,mu2,ep)
!AC!   implicit none
!AC!   complex(ki_gol), intent(in) :: m0
!AC!   real(ki_gol), intent(in) :: mu2
!AC!   integer, intent(in) :: ep
!AC!   complex(ki_gol) :: gA0C
!AC!   if(ep.eq.(-2) .or. m0.eq.(0.0_ki_gol,0.0_ki_gol)) then
!AC!      gA0C = (0.0_ki_gol, 0.0_ki_gol)
!AC!   elseif(ep.eq.(-1)) then
!AC!      gA0C = m0 * gB0C((0.0_ki_gol,0.0_ki_gol),m0,m0,mu2,-1)
!AC!   else
!AC!      gA0C = m0 * (gB0C((0.0_ki_gol,0.0_ki_gol),m0,m0,mu2,0) &
!AC!          &    +  gB0C((0.0_ki_gol,0.0_ki_gol),m0,m0,mu2,-1))
!AC!   end if
!AC!end function gA0C

!AC!pure subroutine gtrunc_rm(ref,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
!AC!   implicit none
!AC!   real(ki), intent(in) :: ref
!AC!   real(ki), intent(inout), optional :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
!AC!   real(ki), parameter :: small = 1.0E-08_ki
!AC!   if(present(s1)) then
!AC!      if(abs(s1/ref) .lt. small) s1 = 0.0_ki
!AC!   end if
!AC!   if(present(s2)) then
!AC!      if(abs(s2/ref) .lt. small) s2 = 0.0_ki
!AC!   end if
!AC!   if(present(s3)) then
!AC!      if(abs(s3/ref) .lt. small) s3 = 0.0_ki
!AC!   end if
!AC!   if(present(s4)) then
!AC!      if(abs(s4/ref) .lt. small) s4 = 0.0_ki
!AC!   end if
!AC!   if(present(s5)) then
!AC!      if(abs(s5/ref) .lt. small) s5 = 0.0_ki
!AC!   end if
!AC!   if(present(s6)) then
!AC!      if(abs(s6/ref) .lt. small) s6 = 0.0_ki
!AC!   end if
!AC!   if(present(s7)) then
!AC!      if(abs(s7/ref) .lt. small) s7 = 0.0_ki
!AC!   end if
!AC!   if(present(s8)) then
!AC!      if(abs(s8/ref) .lt. small) s8 = 0.0_ki
!AC!   end if
!AC!   if(present(s9)) then
!AC!      if(abs(s9/ref) .lt. small) s9 = 0.0_ki
!AC!   end if
!AC!   if(present(s10)) then
!AC!      if(abs(s10/ref) .lt. small) s10 = 0.0_ki
!AC!   end if
!AC!end  subroutine gtrunc_rm
!AC!pure subroutine gtrunc_cm(ref,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
!AC!   implicit none
!AC!   real(ki), intent(in) :: ref
!AC!   complex(ki), intent(inout), optional :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
!AC!   real(ki), parameter :: small = 1.0E-08_ki
!AC!   if(present(s1)) then
!AC!      if(abs(s1/ref) .lt. small) s1 = (0.0_ki, 0.0_ki)
!AC!   end if
!AC!   if(present(s2)) then
!AC!      if(abs(s2/ref) .lt. small) s2 = (0.0_ki, 0.0_ki)
!AC!   end if
!AC!   if(present(s3)) then
!AC!      if(abs(s3/ref) .lt. small) s3 = (0.0_ki, 0.0_ki)
!AC!   end if
!AC!   if(present(s4)) then
!AC!      if(abs(s4/ref) .lt. small) s4 = (0.0_ki, 0.0_ki)
!AC!   end if
!AC!   if(present(s5)) then
!AC!      if(abs(s5/ref) .lt. small) s5 = (0.0_ki, 0.0_ki)
!AC!   end if
!AC!   if(present(s6)) then
!AC!      if(abs(s6/ref) .lt. small) s6 = (0.0_ki, 0.0_ki)
!AC!   end if
!AC!   if(present(s7)) then
!AC!      if(abs(s7/ref) .lt. small) s7 = (0.0_ki, 0.0_ki)
!AC!   end if
!AC!   if(present(s8)) then
!AC!      if(abs(s8/ref) .lt. small) s8 = (0.0_ki, 0.0_ki)
!AC!   end if
!AC!   if(present(s9)) then
!AC!      if(abs(s9/ref) .lt. small) s9 = (0.0_ki, 0.0_ki)
!AC!   end if
!AC!   if(present(s10)) then
!AC!      if(abs(s10/ref) .lt. small) s10 = (0.0_ki, 0.0_ki)
!AC!   end if
!AC!end  subroutine gtrunc_cm

end module madds

