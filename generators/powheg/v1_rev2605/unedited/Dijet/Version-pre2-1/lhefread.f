c...lhefheader(nlf)
c...reads initialization information from a les houches events file on unit nlf. 
      subroutine lhefreadhdr(nlf)
      implicit none
      integer nlf
      character * 200 string
      integer ipr
      include 'LesHouches.h'
 1    read(nlf,fmt='(a)',err=998,end=998) string
      if(string(1:5).eq.'<init') then
         read(nlf,*) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &        pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
         do ipr=1,nprup
            read(nlf,*) xsecup(ipr),xerrup(ipr),xmaxup(ipr),
     &           lprup(ipr)
         enddo
         goto 999
      else
         goto 1
      endif
 998  write(*,*) 'lhefreadhdr: could not find <init> data'
      call exit(1)
 999  end


c...reads event information from a les houches events file on unit nlf. 
      subroutine lhefreadev(nlf)
      implicit none
      integer nlf
      character * 200 string
      include 'LesHouches.h'
      integer i,j
 1    continue
      string=' '
      read(nlf,fmt='(a)',err=777,end=666) string
      if(string.eq.'</LesHouchesEvents>') then
         goto 998
      endif
      if(string(1:6).eq.'<event') then
c on error try next event. The error may be cause by merging
c truncated event files. On EOF return with no event found
         read(nlf,*,end=998,err=1)nup,idprup,xwgtup,scalup,aqedup,aqcdup
         do i=1,nup
            read(nlf,*,end=998,err=1) idup(i),istup(i),mothup(1,i),
     &           mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &           vtimup(i),spinup(i)
         enddo
         call lhefreadextra(nlf)
         goto 999
      else
         goto 1
      endif
c no event found:
 777  continue
      print *,"Error in reading"
      print *,string
      nup=0
      return
 666  continue
      print *,"reached EOF"
      print *,string
      nup=0
      return
 998  continue
      print *,"read </LesHouchesEvents>"
      nup=0      
 999  end


      subroutine lhefreadextra(nlf)
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_weights.h'
      character * 200 string
      integer nlf
      weights_num = 0
 1    continue
      string=' '
      read(unit=nlf,fmt='(a)',end=998) string
      string=adjustl(string)
      if(string.eq.'<event>') then
         backspace nlf
         return
      endif
      if(string(1:11).eq.'#new weight') then
         if(weights_num.eq.weights_max) then
            write(*,*) ' too many weights!'
            write(*,*) ' increase weights_max'
            call exit(-1)
         endif
         weights_num = weights_num + 1
         read(string(38:),*) weights_val(weights_num),
     1                       weights_renfac(weights_num),
     2                       weights_facfac(weights_num),
     3                       weights_npdf1(weights_num),
     4                       weights_npdf2(weights_num),
     5                       weights_whichpdf(weights_num)
      endif
      if(string.eq.'# Start extra-info-previous-event') then
         read(nlf,'(a)') string
         read(string(3:),*) rad_kinreg
         read(nlf,'(a)') string
         read(string(3:),*) rad_type
      endif
      goto 1
      return
 998  continue
      end
