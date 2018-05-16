c Find the underlying Born momenta from the real momenta and the
c emitter-readiated pair
      subroutine compuborn(em,rad,cmppborn)
      implicit none
      integer em,rad
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 cmppborn(0:3,nlegreal)
      if(em.lt.3) then
         call findubisr(rad,cmppborn)
      else
         call findubfsr(em,rad,cmppborn)
      endif
      end
      
      subroutine findubfsr(i,j,cmppborn)
      implicit none
      integer i,j
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 cmppborn(0:3,nlegreal)
      include 'pwhg_kn.h'
      real * 8 krecv(3),q0,q2,krec,beta,
     1 k0rec,k,vec(3)
      cmppborn(0:3,1)=kn_cmpreal(0:3,1)
      cmppborn(0:3,2)=kn_cmpreal(0:3,2)
      q0=2*cmppborn(0,1)
      q2=q0**2
c recoil system momentum 
      k0rec=q0-kn_cmpreal(0,i)-kn_cmpreal(0,j)
      krecv=-kn_cmpreal(1:3,i)-kn_cmpreal(1:3,j)
      krec=sqrt(krecv(1)**2+krecv(2)**2+krecv(3)**2)
      beta=(q2-(k0rec+krec)**2)/(q2+(k0rec+krec)**2)
      vec=krecv/krec
      call mboost(nlegreal-2,vec,beta,
     1     kn_cmpreal(:,3:nlegreal),cmppborn(:,3:nlegreal))
      k=(q2-(k0rec**2-krec**2))/(2*q0)
      cmppborn(0,i)=k
      cmppborn(1:3,i)=-vec*k
      cmppborn(:,j)=0
      end

      subroutine findubisr(j,cmppborn)
      implicit none
      integer j
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 cmppborn(0:3,nlegreal)
      include 'pwhg_kn.h'
      real * 8 krecv(3),q0,q2,krec,k0rec,
     1     krecperp,mrec2,beta,vec(3)
      cmppborn(0:3,1)=kn_cmpreal(0:3,1)
      cmppborn(0:3,2)=kn_cmpreal(0:3,2)
      q0=2*cmppborn(0,1)
      q2=q0**2
c recoil system momentum 
      k0rec=q0-kn_cmpreal(0,j)
      krecv=-kn_cmpreal(1:3,j)
      krec=sqrt(krecv(1)**2+krecv(2)**2+krecv(3)**2)
      mrec2=(k0rec**2-krec**2)
      beta=-krecv(3)/k0rec
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegreal-2,vec,beta,
     1     kn_cmpreal(:,3:nlegreal),cmppborn(:,3:nlegreal))
c Now the transverse boost
      krecperp=sqrt(krecv(1)**2+krecv(2)**2)
      vec(3)=0
      vec(1:2)=krecv(1:2)/krecperp
      beta=-krecperp/sqrt(mrec2+krecperp**2)
      call mboost(nlegreal-2,vec,beta,
     1     cmppborn(:,3:nlegreal),cmppborn(:,3:nlegreal))
      cmppborn(:,j)=0
      end

      function dijterm(em,rad,alr)
      implicit none
      real * 8 dijterm
      integer em,rad,alr
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_par.h'
      integer rflav(nlegreal)
      real * 8 cmppborn(0:3,nlegreal)
      real * 8 avub,getdistance1,getdistance2,getdistance,dalr
      integer nub,mergeisr,mergefsr,i,j,k,ifl1,ifl2,onem
      parameter (onem=1000000)
      logical ini,olddij
      data ini/.true./
      save ini,olddij
      real * 8 powheginput
      external powheginput
      if(ini) then
         if(powheginput("#olddij").eq.1) then
            olddij=.true.
         else
            olddij=.false.
         endif
         ini=.false.
      endif
      if(olddij) then
         dijterm=kn_dijterm(em,rad)
         return
      endif
      rflav(:)=flst_alr(:,alr)
c find UB flavour
      if(em.eq.0) then
         continue
      elseif(em.eq.1) then
         rflav(1)=mergeisr(rflav(1),rflav(rad))
      elseif(em.eq.2) then
         rflav(2)=mergeisr(rflav(2),rflav(rad))
      else
         rflav(em)=mergefsr(rflav(em),rflav(rad))
      endif
c invalidater rad parton with impossible pdg code
      rflav(rad)=onem
      call compuborn(em,rad,cmppborn)
c looop over all possible singularities
      dalr=getdistance(em,rad,kn_cmpreal)
      if(abs(dalr/kn_dijterm(em,rad)-1).gt.1d-6) then
         write(*,*) 'dalr', dalr/kn_dijterm(em,rad)
      endif
c get average singularity from underlying Born
      avub=1
      nub=0
      do j=3,nlegreal
         ifl1=mergeisr(rflav(1),rflav(j))
         ifl2=mergeisr(rflav(2),rflav(j))
         if(ifl1.eq.rflav(1).and.ifl2.eq.rflav(2)) then
            avub=avub*(1+dalr/getdistance(0,j,cmppborn))
            nub=nub+1
         else
            if(ifl1.ne.onem) then
               avub=avub*(1+dalr/getdistance(1,j,cmppborn))
               nub=nub+1
            endif
            if(ifl2.ne.onem) then
               avub=avub*(1+dalr/getdistance(2,j,cmppborn))
               nub=nub+1
            endif
         endif
         do k=j+1,nlegreal
            if(mergefsr(rflav(j),rflav(k)).ne.onem) then
               avub=avub*(1+dalr/getdistance(j,k,cmppborn))
               nub=nub+1
            endif
         enddo
      enddo
      dijterm=dalr*avub
      end

      function getdistance(em,rad,cmp)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_par.h'
      real * 8 getdistance
      integer em,rad
      real * 8 cmp(0:3,nlegreal),y
      real * 8 dotp
      external dotp
      if(em.lt.3) then
         y=1-dotp(cmp(0,1),cmp(0,rad))/(cmp(0,1)*cmp(0,rad))
         if(em.eq.0) then
            getdistance=(cmp(0,rad)**2*(1-y**2))**par_diexp
         elseif(em.eq.1) then
            getdistance=(cmp(0,rad)**2*2*(1-y))**par_diexp
         elseif(em.eq.2) then
            getdistance=(cmp(0,rad)**2*2*(1+y))**par_diexp
         endif
      else
         getdistance=(2*dotp(cmp(0,em),cmp(0,rad))*
     1        cmp(0,em)*cmp(0,rad)/(cmp(0,em)+cmp(0,rad))**2
     2        )**par_dijexp
      endif
      end

      function mergeisr(i,j)
      implicit none
      integer mergeisr,i,j,onem
      parameter (onem = 1000000)
      if(abs(i).gt.5.or.abs(j).gt.5) then
         mergeisr = onem
      elseif(j.eq.0) then
         mergeisr = i
      elseif(i.eq.0) then
         mergeisr = -j
      elseif(i.eq.j) then
         mergeisr = 0
      else
         mergeisr = onem
      endif
      end

      function mergefsr(i,j)
      implicit none
      integer mergefsr,i,j,onem
      parameter (onem = 1000000)
      if(abs(i).gt.5.or.abs(j).gt.5) then
         mergefsr = onem
      elseif(j.eq.0) then
         mergefsr = i
      elseif(i.eq.0) then
         mergefsr = j
      elseif(i.eq.-j) then
         mergefsr = 0
      else
         mergefsr = onem
      endif
      end
