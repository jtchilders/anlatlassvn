c...lhefheader(nlf)
c...writes initialization information to a les houches events file on unit nlf. 
      subroutine lhefwritehdr(nlf)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      integer nlf
      real * 8 version
      common/pwghvq/version
      integer ipr,iran,n1ran,n2ran
      include 'LesHouches.h'
      write(nlf,'(a)') '<LesHouchesEvents version="1.0">'
      write(nlf,'(a)') '<!--'
      write(nlf,'(a,f3.1)') 'file generated with POWHEG-BOX version ',
     #     version
      write(nlf,'(a)') 'Input file powheg.input contained:'
      call wrtpowheginput(nlf)
      write(nlf,'(a)') 'End of powheg.input content'
      call rm48ut(iran,n1ran,n2ran)
      write(nlf,*) 'Random number generator initialized with: ',
     # iran,' ',n1ran,' ',n2ran
      write(nlf,'(a)') '-->'
      write(nlf,'(a)') '<init>'
      write(nlf,110) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     &pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
      do 100 ipr=1,nprup
         write(nlf,120) xsecup(ipr),xerrup(ipr),xmaxup(ipr),
     &        lprup(ipr)
 100  continue
      write(nlf,'(a)') '</init>'
 110  format(1p,2(1x,i8),2(1x,e12.5),6(1x,i6))
 120  format(1p,3(1x,e12.5),1x,i6)
      end


      subroutine printleshouches
c useful for debugging
      call lhefwritev(6)
      end

c...lhefeader(nlf)
c...writes event information to a les houches events file on unit nlf. 
      subroutine lhefwritev(nlf)
      implicit none
      integer nlf
      include 'LesHouches.h'
      include 'pwhg_flg.h'
      integer i,j
      write(nlf,'(a)')'<event>'
      write(nlf,210) nup,idprup,xwgtup,scalup,aqedup,aqcdup
      do 200 i=1,nup
         write(nlf,220) idup(i),istup(i),mothup(1,i),
     & mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     & vtimup(i),spinup(i)
 200  continue
      if(flg_pdfreweight) call lhefwritepdfrw(nlf)
      if(flg_debug) call lhefwritextra(nlf)
      write(nlf,'(a)')'</event>'      
 210  format(1p,2(1x,i6),4(1x,e12.5))
 220  format(1p,i8,5(1x,i5),5(1x,e16.9),1x,e12.5,1x,e10.3)
      end

c...lheftrailer(nlf)
c...writes trailer to a les houches events file on unit nlf. 
      subroutine lhefwritetrailer(nlf)
      implicit none
      integer nlf,iran,n1ran,n2ran
      write(nlf,'(a)') '</LesHouchesEvents>'
c     save last random number
      call rm48ut(iran,n1ran,n2ran)
      write(nlf,*) '#Random number generator exit values: ',
     # iran,' ',n1ran,' ',n2ran
      end

      subroutine lhefwritextra(nlf)
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      integer nlf
      integer iran,n1ran,n2ran
      write(nlf,'(a)') '# Start extra-info-previous-event'
      write(nlf,*) '# ',rad_kinreg,'       rad_kinreg'
      write(nlf,*) '# ',rad_type,'         rad_type'
      write(nlf,*) '# ',rad_ubornidx,'         rad_ubornidx'
      write(nlf,*) '# ',rad_realalr,'         rad_realalr'
      call getcurrentrandom(iran,n1ran,n2ran)
      write(nlf,*) '# ', iran,' ',n1ran,' ',n2ran,
     #     "    previous event's random seeds "
      write(nlf,'(a)') '# End extra-info-previous-event'
      end

      subroutine lhefwritepdfrw(nlf)
      implicit none
      integer nlf
      integer id1,id2
      real * 8 x1,x2,xf1,xf2,xmufact
      call pdfreweightinfo(id1,id2,x1,x2,xmufact,xf1,xf2)
      write(nlf,111)'#pdf ',id1,id2,x1,x2,xmufact,xf1,xf2
 111  format(a,2(1x,i2),5(1x,d14.8))
      end

