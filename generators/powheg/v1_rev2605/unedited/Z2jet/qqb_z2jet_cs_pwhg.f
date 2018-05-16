      subroutine qqb_z2jet_cs_pwhg(bflav,bornjk)
      implicit none
      include 'constants.f'
      include 'flags.f'
      include 'msq_cs.f'
      integer nlegs,bflav(6),bf1,bf2,bf5,bf6,j,k,ig,ih,ia,iq
      parameter (nlegs=6)
      logical ggaq,qagg,aqgg,qgqg,agag,gqqg,gaag,aaaa,aqqa,qqqq,qaqa,
     & test
      double precision bornjk(6,nlegs),lo

C--initialize to zero
      bornjk(:,:)=zip

      bf1=bflav(1)
      bf2=bflav(2)
      bf5=bflav(5)
      bf6=bflav(6)
      qagg=(bf1>0).and.(bf2<0).and.(bf5==0).and.(bf6==0)
      aqgg=(bf1<0).and.(bf2>0).and.(bf5==0).and.(bf6==0)
      qgqg=(bf1>0).and.(bf2==0).and.(bf5>0).and.(bf6==0)
      agag=(bf1<0).and.(bf2==0).and.(bf5<0).and.(bf6==0)
      gqqg=(bf2>0).and.(bf1==0).and.(bf5>0).and.(bf6==0)
      gaag=(bf2<0).and.(bf1==0).and.(bf5<0).and.(bf6==0)
      ggaq=(bf1==0).and.(bf2==0).and.(bf5<0).and.(bf6>0)
      aaaa=(bf1<0).and.(bf2<0).and.(bf5<0).and.(bf6<0)
      aqqa=(bf1<0).and.(bf2>0).and.(bf5>0).and.(bf6<0)
      qqqq=(bf1>0).and.(bf2>0).and.(bf5>0).and.(bf6>0)
      qaqa=(bf1>0).and.(bf2<0).and.(bf5>0).and.(bf6<0)

      IF (Gflag) then
C-----Remember total cross section is              
C-----msq_cs(1,j,k)+msq_cs(2,j,k)+msq_cs(0,j,k)    
C-----where      
C-----msq_cs(1,j,k) propto A1^2   
C-----msq_cs(2,j,k) propto A2^2   
C-----msq_cs(0,j,k) propto -1/xn^2*|(A1+A2)|^2     

      if (ggaq) then
      ig=1
      ih=2
      ia=5
      iq=6
      elseif (qagg) then
C-----5<-->1,2<-->6 of above
      ig=5
      ih=6
      ia=1
      iq=2
      elseif (aqgg) then
c--- (1<-->2) of above
      ig=5
      ih=6
      ia=2
      iq=1
      elseif (gqqg) then
C----(2<--5) of qagg
      ig=1
      ih=6
      ia=2
      iq=5
      elseif (qgqg) then
C---(1<-->2) of above
      ig=2
      ih=6
      ia=1
      iq=5
      elseif (gaag) then
      ig=1
      ih=6
      ia=5
      iq=2
      elseif (agag) then
C---(1<-->2) of above
      ig=2
      ih=6
      ia=5
      iq=1
      else 
      write(6,*) ' qqb_z2jet_cs_pwhg:Unimplemented GFLAG option'
      stop
      endif
      lo=msq_cs(1,bf1,bf2)+msq_cs(2,bf1,bf2)+msq_cs(0,bf1,bf2)
      bornjk(ig,ig)=-xn*lo
      bornjk(ih,ih)=-xn*lo
      bornjk(ia,ia)=-cf*lo
      bornjk(iq,iq)=-cf*lo
      bornjk(ig,ih)=xn/2d0*(msq_cs(1,bf1,bf2)+msq_cs(2,bf1,bf2))
      bornjk(ig,ia)=xn/2d0*(msq_cs(1,bf1,bf2)+msq_cs(0,bf1,bf2))
      bornjk(ig,iq)=xn/2d0*(msq_cs(2,bf1,bf2)+msq_cs(0,bf1,bf2))
      bornjk(ia,iq)=-0.5d0/xn*(lo+xn**2*msq_cs(0,bf1,bf2))
      bornjk(ih,ia)=bornjk(ig,iq)
      bornjk(ih,iq)=bornjk(ig,ia)
C----symmetrization
      bornjk(ih,ig)=bornjk(ig,ih)
      bornjk(ia,ig)=bornjk(ig,ia)
      bornjk(iq,ig)=bornjk(ig,iq)
      bornjk(iq,ia)=bornjk(ia,iq)
      bornjk(ia,ih)=bornjk(iq,ig)
      bornjk(iq,ih)=bornjk(ia,ig)

      elseif (Qflag) then
      lo=msq_cs(0,bf1,bf2)+msq_cs(1,bf1,bf2)+msq_cs(2,bf1,bf2)
      bornjk(1,1)=-CF*lo;
      bornjk(2,2)=-CF*lo;
      bornjk(5,5)=-CF*lo;
      bornjk(6,6)=-CF*lo;

      if (aaaa .or. qqqq) then
      bornjk(1,2)=CF*msq_cs(0,bf1,bf2)+lo/xn;
      bornjk(1,5)=-lo/2/xn+CF*msq_cs(2,bf1,bf2);
      bornjk(1,6)=-lo/2/xn+CF*msq_cs(1,bf1,bf2);
      bornjk(2,5)=-lo/2/xn+CF*msq_cs(1,bf1,bf2);
      bornjk(2,6)=-lo/2/xn+CF*msq_cs(2,bf1,bf2);
      bornjk(5,6)=CF*msq_cs(0,bf1,bf2)+lo/xn;
      elseif (aqqa) then
      bornjk(1,2)=-lo/2/xn+CF*msq_cs(2,bf1,bf2);
      bornjk(1,5)=CF*msq_cs(0,bf1,bf2)+lo/xn;
      bornjk(1,6)=-lo/2/xn+CF*msq_cs(1,bf1,bf2);
      bornjk(2,5)=-lo/2/xn+CF*msq_cs(1,bf1,bf2);
      bornjk(2,6)=CF*msq_cs(0,bf1,bf2)+lo/xn;
      bornjk(5,6)=-lo/2/xn+CF*msq_cs(2,bf1,bf2);
      elseif (qaqa) then
      bornjk(1,2)=-lo/2/xn+CF*msq_cs(2,bf1,bf2);
      bornjk(2,5)=CF*msq_cs(0,bf1,bf2)+lo/xn;
      bornjk(2,6)=-lo/2/xn+CF*msq_cs(1,bf1,bf2);
      bornjk(1,5)=-lo/2/xn+CF*msq_cs(1,bf1,bf2);
      bornjk(1,6)=CF*msq_cs(0,bf1,bf2)+lo/xn;
      bornjk(5,6)=-lo/2/xn+CF*msq_cs(2,bf1,bf2);
      else
      write(6,*) 'Unimpemented ordering in qqb_w2jet_cs_pwhg'
      stop
      endif

C fill other non-zero values
      do j=1,nlegs-1
         do k=j+1,nlegs
            bornjk(k,j)=bornjk(j,k)
         enddo
      enddo

      endif

C      test= (abs(sum(bornjk(1,:))/lo) > 1d-11)
C     & .or. (abs(sum(bornjk(2,:))/lo) > 1d-11)
C     & .or. (abs(sum(bornjk(5,:))/lo) > 1d-11)
C     & .or. (abs(sum(bornjk(6,:))/lo) > 1d-11)
C      if (test) then
C      write(6,*) bf1,bf2,bf5,bf6
C      write(6,*) 'sum1',sum(bornjk(1,:))
C      write(6,*) 'sum2',sum(bornjk(2,:))
C      write(6,*) 'sum5',sum(bornjk(5,:))
C      write(6,*) 'sum6',sum(bornjk(6,:))
C      stop 'qqb_z2jet_cs_pwhg'
C      endif


       return
       end
