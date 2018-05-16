      subroutine pwhgrescale_single(i,scale)
C********************************************************************
C The following function rescale the content of the histogram i by the
C scale factor.  IMPORTANT: must be called AFTER the statistical
C analysis has been performed and before histograms are written to the
C TOPDRAWER file. For example : 
C    call pwhgsetout ( or pwhgaddout)
C    call pwhgrescale_single(scale)
C    call pwhgtopout 
C**********************************************************************
      implicit none
      include 'pwhg_book.h'
      real * 8 scale,dummy_real
      integer i,dummy_int
      character * 3 tag
      call pwhggettag(i,tag)
      if(tag.eq.'YST') then
C PWHGOPERA(I,'F',J,K,X,Y) multiplies hist I by the factor X, and puts
C the result in hist K;
         call pwhgopera(i+nmh3,'F',dummy_int,i+nmh3,scale,dummy_real)
C i+nmh3 contains the values 
         call pwhgopera(i+nmh4,'F',dummy_int,i+nmh4,scale,dummy_real)
C i+nmh4 contains the errors
      endif
      end           
      

      subroutine pwhgdiscardup
      implicit none
      include 'pwhg_book.h'
      integer j,l
      character * 3 tag
C
c     nmb and nxm defined in  pwhg_book.h 
c     nmb is the number of histograms, nxm the number of bins   
c
c     Empty the original histogram j without accumulating the values
c     Useful when an error occurs and the points have to be discarded
      
      do j=1,nmh
         call pwhggettag(j,tag)
         if(tag.eq.'YST') then
            do l=1,nbin(j)
               hist(j,l)=0
            enddo
            ient(j)=0
            uscore(j)=0
            oscore(j)=0
         endif
      enddo
      call increasecnt("discarded points")
      return
      end


