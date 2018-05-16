      subroutine OLP_Start(contract_file_name,ierr)
      implicit none
      character * 50 contract_file_name
      integer ierr
      end
      
      subroutine OLP_EvalSubProcess(label, momenta, mu, parameters, res)
      implicit none
      integer label
      real * 8 mu
      real * 8 momenta(50)
      real * 8 parameters(10)
      real * 8 res(4)
      end

      subroutine OLP_Finalize()
      implicit none
      end

      subroutine OLP_Option(line,stat)
      implicit none
      character * 50 line
      integer stat
      end
