c     slepton pair production
c     Born contribution:
c     parton parton -> slepton slepton
c     2012-04 Andreas v. Manteuffel

c      -6  -5  -4  -3  -2  -1  0  1  2  3  4  5  6
c      t~  b~  c~  s~  u~  d~  g  d  u  s  c  b  t


c     sets Born cross section
      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      integer bflav(nlegborn)
      real *8 p(0:3,nlegborn),bornjk(nlegborn,nlegborn)
      real *8 bmunu(0:3,0:3,nlegborn),born,colcf
      integer j,k,mu,nu

      if ( bflav(2).ne.-bflav(1) ) then
         print*,'the current Born subprocess is ',bflav
         print*,'check your settings in init_processes'
         stop
      endif
      ! born without color correlations
      call compborn(p,bflav,born)
      ! spin and color correlations
      do j=1,nlegborn
         if(abs(bflav(j)).le.6) then
      ! gluon spin correlations as in (2.8) of 1002.2581v1
      ! no contribution for external gluon: Bmunu=0
            if(bflav(j).eq.0) then
               do mu=0,3
                  do nu=0,3
                     bmunu(mu,nu,j) = 0d0
                  enddo
               enddo
            endif
      ! color correlations as in (2.6) of 1002.2581v1
      ! Bij = Cf B for i != j incoming q/qbar
            do k=j+1,nlegborn
               if(abs(bflav(k)).le.6) then
                  colcf = cf
               else
                  colcf = 0d0
               endif
               bornjk(j,k)=born*colcf
               bornjk(k,j)=bornjk(j,k)
            enddo
         endif
      enddo
      end



c     compute helicity averaged Born cross section
c     neutral current slepton pairs
c     q(1) qbar(2) -> Z/gamma -> slepton^-(3) slepton^+(4)
c     average over incoming helicites
      subroutine compborn(p,flav,born)
      implicit none
      include 'nlegborn.h'
      integer flav(nlegborn)
      real *8 p(0:3,nlegborn),bornlr,bornrl,born
      call compbornhel(p,flav,bornlr,bornrl)
      ! incoming quarks assumed massless
      born = (1d0/2)**2 * (bornlr + bornrl)
      end



c     compute helicity dependent Born cross sections
c     neutral current slepton pairs
c     q(1) qbar(2) -> Z/gamma -> slepton^-(3) slepton^+(4)
c     helicity amplitudes for incoming quark helicities
      subroutine compbornhel(p,flav,bornlr,bornrl)
      implicit none
#include "PhysPars.h"
      include 'nlegborn.h'
      include 'pwhg_math.h'
      integer flav(nlegborn)
      real *8 p(0:3,nlegborn),bornlr,bornrl
      real *8 propfac_a,s,t,u,p3sq,p4sq,pre
      complex *16 propfac_z
      real *8 aapre,aalr,aarl,zapre,zalr,zarl,zzpre,zzlr,zzrl
      real *8 delta_ij,e4
      real *8 e_q, t3_q, gl_q, gr_q
      real *8 e_l, t3_l, gl_l, gr_l
      complex *16 cgl_sl, cgr_sl
      real *8 sprod, kronecker ! functions from funcbasic.f
      real *8 alemeff
      logical verbose
c      parameter (verbose=.false.)
      logical init_born
      data init_born/.true./

      if (init_born) then
         verbose = .true.
         init_born = .false.
      else
         verbose = .false.
      endif

      if (.not.( mod(flav(3),100).eq.mod(-flav(4),100) .and.
     &      (flav(3).eq.1000011 .or. flav(3).eq.2000011 .or.
     &       flav(3).eq.1000013 .or. flav(3).eq.2000013 .or.
     &       flav(3).eq.1000015 .or. flav(3).eq.2000015))) then
        write(*,*) 'compbornhel: final state not supported'
        stop
      endif

      if ((abs(flav(1)).ge.1).and.(abs(flav(1)).le.6)
     &    .and.(flav(1).eq.(-flav(2)))) then
         ! q qbar in initial state
         ! kinematics
         p3sq = sprod(p(0,3),p(0,3))
         p4sq = sprod(p(0,4),p(0,4))
         s =  2*sprod(p(0,1),p(0,2))
         t = -2*sprod(p(0,1),p(0,3)) + p3sq
         u = -2*sprod(p(0,2),p(0,3)) + p3sq
         ! photon and Z propagator factors
         propfac_a  = 1d0
         propfac_z = dcmplx(s)/dcmplx(s-ph_Zmass2,ph_Zmass*ph_Zwidth)
         ! quark couplings
         if (mod(abs(flav(1)),2).eq.0) then
           ! up-type quark
           if (verbose) write(*,*) "up-type quark selected"
           e_q = 2d0 / 3d0
           t3_q = 1d0 / 2d0
         else
           ! down-type quark
           if (verbose) write(*,*) "down-type quark selected"
           e_q = -1d0 / 3d0
           t3_q = -1d0 / 2d0
         endif
         gl_q = -(t3_q - e_q*ph_sthw2)/(ph_sthw*ph_cthw)
         gr_q = -(0d0  - e_q*ph_sthw2)/(ph_sthw*ph_cthw)
         ! lepton couplings (charged lepton l-)
         e_l  = -1d0
         t3_l = -1d0 / 2d0
         gl_l = -(t3_l - e_l*ph_sthw2)/(ph_sthw*ph_cthw)
         gr_l = -(0d0  - e_l*ph_sthw2)/(ph_sthw*ph_cthw)
         ! slepton couplings
         delta_ij = kronecker(abs(flav(3)),abs(flav(4)))
         cgl_sl = gl_l*ph_c_U3l*conjg(ph_c_U4l)
         cgr_sl = gr_l*ph_c_U3r*conjg(ph_c_U4r)
         ! |gamma|^2       (for specific q hel)
         aapre = delta_ij *e_q**2 *e_l**2 *propfac_a**2
         aalr = aapre
         aarl = aapre
         ! 2 Re(Z gamma^*) (for specific q hel)
         zapre = 2d0*dreal(delta_ij*e_q*e_l*propfac_a*propfac_z
     &           *(cgl_sl+cgr_sl))
         zalr = zapre*gl_q
         zarl = zapre*gr_q
         ! |Z|^2           (for specific q hel)
         zzpre = abs(propfac_z)**2 *(abs(cgl_sl+cgr_sl))**2
         zzlr = zzpre*gl_q**2
         zzrl = zzpre*gr_q**2
         ! electromagnetic unit charge to power 4
         e4 = 16d0 * pi**2 * alemeff(s)**2
         ! averaged |M|^2  (still specific q hel)
         pre = (1d0/NC)**2 ! color average
     &         *4d0*NC*e4 *(t*u - p3sq*p4sq)/(s**2)
         bornlr = pre * (aalr + zalr + zzlr)
         bornrl = pre * (aarl + zarl + zzrl)
      else ! no q qbar initial state
         bornlr = 0d0;
         bornrl = 0d0;
      endif

      if (verbose) then
        write(*,'(A, F7.0)') "sqrt(s)=",sqrt(s)
        write(*,'(A, F7.0)') "m3=",sqrt(p3sq)
        write(*,'(A, F7.0)') "m4=",sqrt(p4sq)
        write(*,*) "Born helicity amps:"
        write(*,*) bornlr, bornrl
      endif

      end
