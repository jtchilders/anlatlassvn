      subroutine pdfset(PARM,VALUE)
      call pdfdummieserr('PDFSET')
      end

      subroutine structp
      call pdfdummieserr('STRUCTP')
      end

      subroutine structm
      call pdfdummieserr('STRUCTM')
      end

      subroutine pdfdummieserr(str)
      character *(*) str
      write(*,*) 'Dummy ',str,' interface; link the true one'
      write(*,*) 'if you need it.'
      call exit(-1)
      end

