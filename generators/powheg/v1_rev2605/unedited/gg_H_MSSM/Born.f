      subroutine setborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs)
      integer bflav(nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),born
      real * 8 amp2

      amp2 = 0d0
      call compborn(p,bflav,amp2,bornjk,bmunu)      
      if(amp2.eq.0d0) then
        write(*,*) 'WARNING setborn: returning 0 amplitude!'
        write(*,*) bflav(1),' ',bflav(2),'->',bflav(3)
      endif
      born=amp2
      end

      subroutine compborn(p,bflav,born,bornjk,bmunu)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'Flags.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs),bornjk(nlegs,nlegs)
      integer bflav(nlegs)
      real * 8 bmunu(0:3,0:3,nlegs),born
      integer mu,nu,i,j,k
      real * 8 colcf
      real * 8 gtens(0:3,0:3)
      data gtens/1d0, 0d0, 0d0, 0d0,
     #           0d0,-1d0, 0d0, 0d0,
     #           0d0, 0d0,-1d0, 0d0,
     #           0d0, 0d0, 0d0,-1d0/
      save gtens
      real * 8 prevborn
      save prevborn
      data prevborn /-1d0/
      integer initborn
      save initborn
      data initborn /0/

      call M2_gg_h(p,born)

c     Born amplitude and the color- and spin-correlated ones, needed as
c     counterterms for real contributions, must be evaluated with the
c     approximation, that is full mass dependence.
c
      do j=1,nlegs
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,j)=0d0
            enddo
         enddo
         do k=j+1,nlegs
            bornjk(j,k)=0d0
            bornjk(k,j)=bornjk(j,k)
         enddo
         if(abs(bflav(j)).le.6) then
C     Spin-correlated Born amplitudes
            if(bflav(j).eq.0) then
               if(j.eq.1) then 
                  i=2
               elseif(j.eq.2) then
                  i=1
               else
                  write(*,*) 'Error in bmunu'
                  call exit(1)
               endif
               do mu=0,3
                  do nu=0,3
                     bmunu(mu,nu,j)=born*((p(mu,j)*p(nu,i)+p(nu
     $                    ,j)*p(mu,i))/kn_sborn-gtens(mu,nu)/2d0) 
                  enddo
               enddo
              endif
c     Colour factors for colour-correlated Born amplitudes; Rule from
c     2.98 in FNO2007, leads to B_ij=Cj * B, where i#j
            do k=j+1,nlegs
               if(bflav(k).eq.0) then
                  colcf=ca
                  bornjk(j,k)=born*colcf
                  bornjk(k,j)=bornjk(j,k)
               endif
            enddo
         endif
      enddo
      end
      
      subroutine M2_gg_h(pphy,amp2)
c     Born matrix element times normalizations and averages.
c     IMPORTANT: the flux factor 1/2s is intentionally missing 
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'
      include 'Flags.h'
      real * 8 s,v2,tiny,tmp,xnorm
      parameter (v2=0.70710678118654757d0)
      parameter (tiny=1.d-8)
      real * 8 amp2
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 pphy(0:3,nlegs)
      integer i
      real * 8 m12,m0,y12,y0,bornqcd
      complex * 16 ampl, reduced,x12,x0,aux
      common /bornampl/ampl,bornqcd
      external reduced
      real * 8 getdeltaew,ewcorr
      external getdeltaew
      common /ew/ewcorr

      s=(pphy(0,1)+pphy(0,2))**2-(pphy(1,1)+pphy(1,2))**2
     $     -(pphy(2,1)+pphy(2,2))**2 -(pphy(3,1)+pphy(3,2))**2
      xnorm=1/(2.d0*pi)*1/(pi)*s/(256*v2)

c     1/(2pi) comes from the 2*pi*delta(s-m^2) of phase space
      ampl = dcmplx(0d0)
      if ((ih.eq.1).or.(ih.eq.2)) then

      do i=1,afer
         m12=mfer(i)
         y12=m12**2/mh2
         x12 = reduced(1d0/y12)
         aux = lambdafer(i)*trfer(i)*
     &        (-4d0)*y12*(2d0-(1d0-4d0*y12)*0.5d0*log(x12)**2)
         ampl = ampl+aux
c         write(*,*) 'fermion ', i , ' ampl ', aux
      end do
c     Scalars with full mass dependence
      if (flg_lhscalars.eq.0) then
         do i=1,asca
            m0 = msca(i)
            y0 = m0**2/mh2
            x0 = reduced(1d0/y0)
            aux = lambdasca(i)*trsca(i)*(mmaa/m0)**2*
     &           4d0*y0*(1d0+2d0*y0*0.5d0*log(x0)**2)
            ampl = ampl + aux
c        write(*,*) 'scalar ', afer+i , ' ampl ', aux
         end do
      else
c     Scalars in the light Higgs mass limit
         do i=1,asca
            m0 = msca(i)
            y0 = m0**2/mh2
            aux = lambdasca(i)*trsca(i)*(mmaa/m0)**2*
     &             (-1d0/3d0)
            ampl = ampl + aux
c            write(*,*) 'scalar ', afer+i , ' ampl ', aux
         end do
      endif
c     pseudoscalar
      else
         do i=1,afer
            m12=mfer(i)
            y12=m12**2/mh2
            x12 = reduced(1d0/y12)
            aux = lambdafer(i)*trfer(i)*4d0*y12*0.5d0*log(x12)**2
            ampl = ampl+aux
c          write(*,*) 'fermion ', i , ' ampl ', aux
        end do
c        write(*,*) ampl
      endif

c     If enabled, here we add the EW corrections

      tmp = ampl * dconjg(ampl)

      if (flg_ew.eq.1) then
         ewcorr = getdeltaew(ampl)
      else
         ewcorr = 0d0
      endif

      bornqcd = xnorm*st_alpha*st_alpha*ph_GF*2d0*s*tmp
      amp2 = bornqcd*(1d0+ewcorr)
c     the multiplication for 2s is needed to remove the flux factor
       end

      subroutine borncolour_lh
c     Sets up the colour for the given flavour configuration already
c     filled in the Les Houches interface.  In case there are several
c     colour structure, one should pick one with a probability
c     proportional to the value of the corresponding cross section, for
c     the kinematics defined in the Les Houches interface
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
c colours of incoming gluons
      integer icolgi(2)
      data icolgi/ 502, 501 /
      save icolgi

      icolup(1,1)=icolgi(1)
      icolup(2,1)=icolgi(2)
      icolup(1,2)=icolgi(2)
      icolup(2,2)=icolgi(1)
c neutral particles
      icolup(1,3)=0
      icolup(2,3)=0

      end

      subroutine finalize_lh
c     Set up the resonances whose mass must be preserved
c     on the Les Houches interface.
c
c     Since the Higgs decay is always handled by the SMC
c     no resonance is present here. Higgs boson is a final state
c     particle. 
c      call add_resonance(25,3,4)
      end
