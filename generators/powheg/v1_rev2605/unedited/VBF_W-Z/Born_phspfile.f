      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'process.inc'
      include 'vbfnlo-files/global.inc'
      logical ini
      data ini/.true./
      save ini
      real * 8 xborn(ndiminteg-3)      
      real * 8 powheginput, PS, bornxsec
      common /PS/ PS, bornxsec
      external powheginput
      if(ini) then
         PS=1d0
         PS=powheginput("#Phasespace")
      endif
      if (PS.eq.2d0) then
          if(ini) then
             bornxsec=powheginput("#bornxsec")
             print*, "---------------------------"
             print*, "--     Events as PS      --" 
             print*, "---------------------------"       
             ini=.false.
          endif
          call born_phsp_file(xborn)
      else
          call born_phsp2(xborn)
      endif
      end


      subroutine born_phsp_file(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-3)
      real * 8 m2,xjac,taumin,tau,y,beta,betaCM,vec(3),cth,s,
     #     z,zhigh,zlow, borntemp(0:4,7)
      integer mu,k,j,a,i,uni
      logical ini
      data ini/.true./
      save ini
      real * 8 Vmass2,Vmass2low,Vmass2high,VmVw  
      real * 8 m2jj,pV(0:3),pVmod,pVmod2,pJ(0:3,2),cthj,phij,pcmjj(0:3),
     #     pcmjjmod,ptmp(0:3,2)
      real * 8 mass, dotp ,test   
      external mass , dotp     
      logical check
      parameter(check=.false.)
      logical BW
      parameter (BW=.true.)
      real * 8 epsilon
      parameter (epsilon=1d-10)
      real * 8 pt1cut,pt2cut,pt1,pt2,m2jjmin
      real * 8 BW_fixed,BW_running
      character*16 event
      logical higgsfixedwidth
      save higgsfixedwidth
      real * 8 powheginput
      external powheginput
 
      real*8 K1(0:3),K2(0:3),V(0:3,4),px,py,pz,e


      real * 8 PS, bornxsec
      common /PS/ PS, bornxsec      
      real * 8 sumdps
      integer counter
      common /sumdps/ sumdps, counter
      
      
      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0d0
         enddo
         kn_masses(nlegreal)=0d0
         ini=.false.
         pt1cut = 0d0
         pt2cut = 0d0
         sumdps=0d0
         counter=0
!          call newunit(uni)
         uni=44
         OPEN (uni, FILE='event.total.lhe',err=999)
         m2jjmin = 0d0
         if ((pt1cut.ne.0d0).or.(pt2cut.ne.0d0).or.(m2jjmin.ne.0)) then
            write(*,*) '**************************************'
            write(*,*) '****       CUTS IN PLACE!!!      *****' 
            write(*,*) '**************************************'
         endif
      endif

      
      do i=1,7
         do mu=0,4
            borntemp(mu,i)=0d0
         enddo
      enddo
       E=0d0
      px=0d0
      py=0d0
      pz=0d0

      read (uni,*,END=999) event
      do while (event(1:6).ne.'<event')
           read(uni,fmt='(a)',END=999) event    
      enddo

         read(uni,*,END=999) a,a,kn_jacborn
         do i=1,7 !events

             read(uni,*,END=999) a,a,a,a,a,a,borntemp(1,i),borntemp(2,i),
     &             borntemp(3,i),borntemp(0,i),borntemp(4,i)
         enddo

         read (uni,*,END=999) event

         do i=1,4

             borntemp(0,i)= borntemp(1,i)**2+borntemp(2,i)**2
     1                     + borntemp(3,i)**2
             borntemp(0,i)=sqrt(borntemp(0,i))
            if(i.le.2) then !last momentum from momentum conservation, therefore subtract all p,E of initial state particles
               px=px-borntemp(1,i)
               py=py-borntemp(2,i)
               pz=pz-borntemp(3,i)   
               E=E-borntemp(0,i)                         
            elseif(i.gt.2) then !last momentum from momentum conservation, therefore add all p,E of final state particles
               px=px+borntemp(1,i)
               py=py+borntemp(2,i)
               pz=pz+borntemp(3,i)      
               E=E+borntemp(0,i)                
            endif             
         enddo
          i=6
             borntemp(0,i)= borntemp(1,i)**2+borntemp(2,i)**2
     1                     + borntemp(3,i)**2
             borntemp(0,i)=sqrt(borntemp(0,i))
                       
 !last momentum from momentum conservation, therefore add all p,E of final state particles
               px=px+borntemp(1,i)
               py=py+borntemp(2,i)
               pz=pz+borntemp(3,i)      
               E=E+borntemp(0,i)                

         borntemp(1,7)= -px
         borntemp(2,7)= -py
         borntemp(3,7)= -pz
         
             borntemp(0,7)= borntemp(1,7)**2+borntemp(2,7)**2
     1                     + borntemp(3,7)**2
             borntemp(0,7)=sqrt(borntemp(0,7))         
         E=E+borntemp(0,7) 
         borntemp(0,1)=borntemp(0,1)+E/2.  !make sure that energy & momenta are conserved
         borntemp(0,2)=borntemp(0,2)+E/2.
         borntemp(3,1)=borntemp(3,1)+E/2.
         borntemp(3,2)=borntemp(3,2)-E/2.

         do mu=0,3
         do i=1,2
         kn_pborn(mu,i)=borntemp(mu,i)
         enddo
         
         do i=3,4
         kn_pborn(mu,i+2)=borntemp(mu,i)
         enddo
         do i=5,6
         kn_pborn(mu,i-2)=borntemp(mu,i+1)
         enddo         
         enddo


      kn_xb1=2d0*kn_pborn(0,1)/sqrt(kn_sbeams)
      kn_xb2=2d0*kn_pborn(0,2)/sqrt(kn_sbeams)

      betaCM=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
c     boost in the CM frame
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn,vec,-betaCM,kn_pborn(0,1),kn_cmpborn(0,1))     
      

      kn_minmass=0.1d0 !min. mjj from vbfnlo
      
      kn_sborn=2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,2))
      
      kn_jacborn=kn_jacborn*bornxsec
      kn_jacborn=kn_jacborn*(2d0*kn_sborn) !flux  
      kn_jacborn=kn_jacborn/3.89379323d8 !conversion into pb
      counter=counter+1

      
      if (check) then
         write(*,*) ''
         write(*,*) 'new set'
         do j=1,nlegborn
            write(*,*) 'mom ',j,(kn_pborn(mu,j),mu=0,3)
            write(*,*) 'mass ',j,mass(kn_pborn(0,j))
         enddo
         call checkmomzero(nlegborn,kn_pborn)
      endif

      
      
      if (check) then
         print*, "--------- new event -----------"
         print*, ""
         print*, "x1    = ", kn_xb1
         print*, "x2    = ", kn_xb2
         print*, "wgt   = ", kn_jacborn
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"

         do mu=0,3 
            k1(mu)=kn_pborn(mu,1)
            k2(mu)=kn_pborn(mu,2)
            v(mu,1)=kn_pborn(mu,3)
            v(mu,2)=kn_pborn(mu,4)
            v(mu,3)=kn_pborn(mu,5)
            v(mu,4)=kn_pborn(mu,6)            
            test = k1(mu)+k2(mu)-v(mu,1)-v(mu,2)
     &               -v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
         enddo         
         
         print*, "------ momenta ------------"
         print*, k1
         print*, k2
         print*, v(0,1),v(1,1),v(2,1),v(3,1)
         print*, v(0,2),v(1,2),v(2,2),v(3,2)
         print*, v(0,3),v(1,3),v(2,3),v(3,3)
         print*, v(0,4),v(1,4),v(2,4),v(3,4)         
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotp(k1(0),k1(0))
         print*,"p2^2 = ",dotp(k2(0),k2(0))
         print*,"v1^2 = ",dotp(v(0,1),v(0,1))
         print*,"v2^2 = ",dotp(v(0,2),v(0,2))
         print*,"v3^2 = ",dotp(v(0,3),v(0,3))
         print*,"v4^2 = ",dotp(v(0,4),v(0,4))         
         print*, "-----------------------"
      write(*,*) '----> CM FRAME' 
         print*, "--- momenta conservation ? ---"
         do mu=0,3 
            k1(mu)=kn_cmpborn(mu,1)
            k2(mu)=kn_cmpborn(mu,2)
            v(mu,1)=kn_cmpborn(mu,3)
            v(mu,2)=kn_cmpborn(mu,4)
            v(mu,3)=kn_cmpborn(mu,5)
            v(mu,4)=kn_cmpborn(mu,6)            
            test = k1(mu)+k2(mu)-v(mu,1)-v(mu,2)
     &               -v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
         enddo
         print*, "------ momenta ------------"
         print*, k1
         print*, k2
         print*, v(0,1),v(1,1),v(2,1),v(3,1)
         print*, v(0,2),v(1,2),v(2,2),v(3,2)
         print*, v(0,3),v(1,3),v(2,3),v(3,3)
         print*, v(0,4),v(1,4),v(2,4),v(3,4)         
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotp(k1(0),k1(0))
         print*,"p2^2 = ",dotp(k2(0),k2(0))
         print*,"v1^2 = ",dotp(v(0,1),v(0,1))
         print*,"v2^2 = ",dotp(v(0,2),v(0,2))
         print*,"v3^2 = ",dotp(v(0,3),v(0,3))
         print*,"v4^2 = ",dotp(v(0,4),v(0,4))    
         print*, "-----------------------"
      endif
      return
999   write(*,*) "reopening file"
      CLOSE (uni)
      OPEN (uni, FILE='event.total.lhe',err=999)      
      
      do i=1,7
         do mu=0,4
            borntemp(mu,i)=0d0
         enddo
      enddo
       E=0d0
      px=0d0
      py=0d0
      pz=0d0

      read (uni,*) event
      do while (event(1:6).ne.'<event')
           read(uni,fmt='(a)') event    
      enddo

         read(uni,*) a,a,kn_jacborn
         do i=1,7 !events

         read(uni,*) a,a,a,a,a,a,borntemp(1,i),borntemp(2,i),
     &             borntemp(3,i),borntemp(0,i),borntemp(4,i)
         enddo

         read (uni,*) event

         do i=1,4
             borntemp(0,i)= borntemp(1,i)**2+borntemp(2,i)**2
     1                     + borntemp(3,i)**2
             borntemp(0,i)=sqrt(borntemp(0,i))
            if(i.le.2) then !last momentum from momentum conservation, therefore subtract all p,E of initial state particles
               px=px-borntemp(1,i)
               py=py-borntemp(2,i)
               pz=pz-borntemp(3,i)   
               E=E-borntemp(0,i)                         
            elseif(i.gt.2) then !last momentum from momentum conservation, therefore add all p,E of final state particles
               px=px+borntemp(1,i)
               py=py+borntemp(2,i)
               pz=pz+borntemp(3,i)      
               E=E+borntemp(0,i)                
            endif             
         enddo
          i=6
             borntemp(0,i)= borntemp(1,i)**2+borntemp(2,i)**2
     1                     + borntemp(3,i)**2
             borntemp(0,i)=sqrt(borntemp(0,i))
                       
 !last momentum from momentum conservation, therefore add all p,E of final state particles
               px=px+borntemp(1,i)
               py=py+borntemp(2,i)
               pz=pz+borntemp(3,i)      
               E=E+borntemp(0,i)                

         borntemp(1,7)= -px
         borntemp(2,7)= -py
         borntemp(3,7)= -pz
         
             borntemp(0,7)= borntemp(1,7)**2+borntemp(2,7)**2
     1                     + borntemp(3,7)**2
             borntemp(0,7)=sqrt(borntemp(0,7))         
         E=E+borntemp(0,7) 
         borntemp(0,1)=borntemp(0,1)+E/2.  !make sure that energy & momenta are conserved
         borntemp(0,2)=borntemp(0,2)+E/2.
         borntemp(3,1)=borntemp(3,1)+E/2.
         borntemp(3,2)=borntemp(3,2)-E/2.

         do mu=0,3
         do i=1,2
         kn_pborn(mu,i)=borntemp(mu,i)
         enddo
         
         do i=3,4
         kn_pborn(mu,i+2)=borntemp(mu,i)
         enddo
         do i=5,6
         kn_pborn(mu,i-2)=borntemp(mu,i+1)
         enddo         
         enddo


         kn_xb1=2d0*kn_pborn(0,1)/sqrt(kn_sbeams)
         kn_xb2=2d0*kn_pborn(0,2)/sqrt(kn_sbeams)

      betaCM=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
c     boost in the CM frame
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn,vec,-betaCM,kn_pborn(0,1),kn_cmpborn(0,1))     
      

      kn_minmass=0.1d0 !min. mjj from vbfnlo
      
      kn_sborn=2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,2))
      
      kn_jacborn=kn_jacborn*bornxsec
      kn_jacborn=kn_jacborn*(2d0*kn_sborn) !flux  
      kn_jacborn=kn_jacborn/3.89379323d8 !conversion into pb
      counter=counter+1      
      
      end



