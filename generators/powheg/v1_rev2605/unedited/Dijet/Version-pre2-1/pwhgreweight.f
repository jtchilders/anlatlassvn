      subroutine pwhgnewweight
      implicit none
      integer iret
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      include 'LesHouches.h'
      integer maxev
      integer gen_seed,gen_n1,gen_n2
      common/cgenrand/gen_seed,gen_n1,gen_n2
      real * 8 newweight
      logical pwhg_isfinite
      external pwhg_isfinite
      character * 3 whichpdfpk
      character * 200 string
      call lhefreadevnew(97,99,iret)
      if(iret.lt.0) then
         write(*,*) ' End of event file! Aborting ...'
         call exit(-1)
      endif
      call setrandom(gen_seed,gen_n1,gen_n2)
      if(rad_type.eq.1) then
         call gen_btilderw
         newweight=rad_btilde_arr(rad_ubornidx)*
     1        rad_btilde_sign(rad_ubornidx)
      elseif(rad_type.eq.2) then
         call gen_sigremnantrw
         newweight=rad_damp_rem_arr(rad_realalr)
      elseif(rad_type.eq.3) then
         call gen_sigremnantrw
         newweight=rad_reg_arr(rad_realreg)
      else
         write(*,*) 'Error in pwhgnewweight, invalid rad_type: ',
     $        rad_type
         call exit(-1)
      endif

      if(.not.pwhg_isfinite(newweight)) newweight=0d0
      write(string,*) '#new weight,renfact,facfact,pdf1,pdf2',
     1        xwgtup*newweight/rad_currentweight,st_renfact,
     2        st_facfact,pdf_ndns1,pdf_ndns2,' ',whichpdfpk()
      write(99,'(a)') trim(adjustl(string))
      write(99,'(a)') '</event>'
      end


      subroutine gen_btilderw
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer mcalls,icalls
      include 'cgengrids.h'
      real * 8 xx(ndiminteg)      
      real * 8 btilde
      external btilde
      mcalls=0
      icalls=0
      call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,2,
     #    mcalls,icalls,xx)
      end

      subroutine gen_sigremnantrw
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'cgengrids.h'
      real * 8 xx(ndiminteg)
      integer mcalls,icalls
      real * 8 sigremnant
      external sigremnant
      mcalls=0
      icalls=0
      call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1 xmmmrm,ifoldrm,2,mcalls,icalls,xx)
      end



      subroutine openoutputrw
      implicit none
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      include 'pwhg_rnd.h'
      if(rnd_cwhichseed.ne.'none') then
         open(unit=99,file=pwgprefix(1:lprefix)//'events-rwgt-'
     1        //rnd_cwhichseed//'.lhe'
     2     ,status='unknown')
      else
         open(unit=99,file=pwgprefix(1:lprefix)//'events-rwgt.lhe'
     1     ,status='unknown')
      endif
      end

c...reads event information from a les houches events file on unit nlf. 
      subroutine lhefreadevnew(nlf,nuo,iret)
      implicit none
      integer nlf,nuo,iret
      character * 500 string
      include 'LesHouches.h'
      integer i,j
      iret=0
 1    continue
      string=' '
      read(nlf,fmt='(a)',err=777,end=666) string
      write(nuo,'(a)') trim(string)
      if(string.eq.'</LesHouchesEvents>') then
         goto 998
      endif
      if(string(1:6).eq.'<event') then
c on error try next event. The error may be cause by merging
c truncated event files. On EOF return with no event found
         read(nlf,'(a)') string
         write(nuo,'(a)') trim(string)
         read(string,fmt=*,end=998,err=1)
     1        nup,idprup,xwgtup,scalup,aqedup,aqcdup
         do i=1,nup
            read(nlf,'(a)') string
            write(nuo,'(a)') trim(string)
            read(string,fmt=*,end=998,err=1)
     1           idup(i),istup(i),mothup(1,i),
     &           mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &           vtimup(i),spinup(i)
         enddo
         call lhefreadextrarw(nlf,nuo,iret)
         goto 999
      else
         goto 1
      endif
c no event found:
 777  continue
      write(*,*) "Error in reading"
      write(*,*) string
      call exit(-1)
 666  continue
      iret=-1
      return
 998  continue
      print *,"read </LesHouchesEvents>"
      iret=-1
      nup=0      
 999  end


      subroutine lhefreadextrarw(nlf,nou,iret)
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      character * 500 string
      integer nlf,nou,iret
      logical readrw
      integer gen_seed,gen_n1,gen_n2
      common/cgenrand/gen_seed,gen_n1,gen_n2
      readrw = .false.
 1    continue
      read(unit=nlf,fmt='(a)',end=998) string
      if(string.eq.'</event>') then
c Don't write the end event record; first we must output the new weight
         return
      endif
      if(string.eq.'<event>') then
         if(.not.readrw) then
            write(*,*) 
     $ 'Error in lhefreadextra, while reading rwg informations'
            write(*,*)'Abort run'
            call exit(-1)
         endif
         backspace nlf
         return
      endif
      if(string.eq.'# Start extra-info-previous-event') then
         write(nou,'(a)') trim(string)
         read(nlf,'(a)') string
         write(nou,'(a)') trim(string)
         read(string(3:),*) rad_kinreg
         read(nlf,'(a)') string
         read(string(3:),*) rad_type
      endif
      write(nou,'(a)') trim(string)
      if(flg_newweight) then
c read a string; if it starts with #rwgt, read first rad_type from the
c string, then all other information, depending upon rad_type.
c set readrw to true
         string=adjustl(string)
         if(string(1:5).eq.'#rwgt') then
            string(1:5)=' '
c     do things
c            print*, 'FOUND'
            read(string,*) rad_type
            if(rad_type.eq.1) then
c     btilde
               read(string,*)rad_type,
     $              rad_ubornidx,rad_currentweight,
     $              gen_seed,gen_n1,gen_n2
            elseif(rad_type.eq.2) then
c     remnant
               read(string,*)rad_type,
     $              rad_realalr,rad_currentweight,
     $              gen_seed,gen_n1,gen_n2
            elseif(rad_type.eq.3) then
c     regular
               read(string,*)rad_type,
     $              rad_realreg,rad_currentweight,
     $              gen_seed,gen_n1,gen_n2
            else
               write(*,*) 'Invalid rad_type in lhefwriteevrw: ',rad_type
               call exit(-1)
            endif
c     if all went ok, set readrw to true
            readrw=.true.
         endif
      endif
      goto 1
 998  continue
      iret=-1
      end
