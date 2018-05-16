      subroutine setreal(pin,rflav,realres)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'mcfmtopwhg.f'
      integer nlegs,ro
      parameter (nlegs=nlegreal)
      double precision pin(0:3,nlegs),mcfmp(mxpart,4),realres,tmp
      integer rflav(nlegs),perm(nlegs),rreal(nlegs),j,jperm
      logical flavperm      
c set mcfm strong coupling and scales
      call st_mcfm
c ---

c find the rflav structure in the list of real diagrams, and find the appropriate                                 
c permutation                                                                                                     
      do j=1,flst_nreal
         if(flavperm(nlegreal,flst_real(:,j),rflav,perm)) go to 10
      enddo
      write(*,*) ' setreal:',rflav,' not found'
      call exit(-1)
 10      continue
      rreal=flst_real(:,j)
      if(perm(5).eq.5) then
         if(perm(6).eq.6) then
            jperm=1
         else
            jperm=2
         endif
      elseif(perm(5).eq.6) then
         if(perm(6).eq.5) then
            jperm=3
         else
            jperm=4
         endif
      elseif(perm(5).eq.7) then
         if(perm(6).eq.5) then
            jperm=5
         else
            jperm=6
         endif
      endif

*  MCFM notation
*     q(-p1) +qbar(-p2) -> nu(p3)+e+(p4)+j(p5)+j(p6)+j(p7)
*     q(-p1) +qbar(-p2) -> e-(p3)+nu~(p4)+j(p5)+j(p6)+j(p7)

      do j=1,nlegs
         do ro=1,4
            if (j .le. 2) then
c     incoming partons                                                                                                
               mcfmp(j,ro)=-pin(pwhg(ro),perm(j))
            else
               mcfmp(j,ro)=+pin(pwhg(ro),perm(j))
            endif
         enddo
      enddo

      call qqb_z2jet_g_pwhg_nc(mcfmp,rreal,jperm,realres)
cc Code to test teh caching mechanism
c      call qqb_z2jet_g_pwhg_nocache(mcfmp,rreal,realres)
c      if(realres.ne.0) then
c         if(abs(tmp/realres-1).gt.1d-15) then
c            write(*,*) ' _nc does not work',tmp,realres
c         endif
c      else
c         if(tmp.ne.0) then
c            write(*,*) ' _nc does not work',tmp,realres
c         endif
c      endif

c Cancel as/(2pi) It will be put back by real_ampsq
      realres = realres/ason2pi
       
      return
      end


      function flavperm(n,aflav,bflav,perm)
c returns true if the flavour structures aflav and bflav are                                                          
c equivalent up to a permutation of the final state lines,                                                            
c with aflav(j)=bflav(perm(j))                                                                                        
c false otherwise.                                                                                                    
      implicit none
      logical flavperm
      integer n, aflav(n),bflav(n),perm(n)
      integer j,k,itmp
      do j=1,n
         perm(j)=j
      enddo
      do j=1,n
         if(aflav(j).ne.bflav(perm(j))) then
            if(j.le.2) then
               flavperm=.false.
               return
            endif
            do k=j+1,n
               if(aflav(j).eq.bflav(perm(k))) then
                  itmp=perm(j)
                  perm(j)=perm(k)
                  perm(k)=itmp
                  goto 10
               endif
            enddo
            flavperm=.false.
            return
         endif
 10            continue
      enddo
      flavperm=.true.
      end

