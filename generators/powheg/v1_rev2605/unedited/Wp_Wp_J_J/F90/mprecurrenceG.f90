!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECrecurrenceG.f90.

module mprecurrence
  use mpmodule; use mpconverter
  use mpsimpleoperations; use mpadvancedoperations
  use mpsimpleoperations_c; use mpadvancedoperations_c
  use types
  use mprecurrencebitsone
  use mprecurrencebitstwo
  use mprecurrencebitsthree
  use mprecurrencebitsfour
  implicit none 
  public  


  !====================================================================
  !                  Notation FOR CURRENTS
  !      the generic notation is a_bcd ... etc;
  !
  !      ``a''  refers to an off-shell line (e.g. g, f, bf), 
  !
  !      if the fermion direction of the off-shell line is incoming 
  !      it is refered to as ``f'', if it is outgoing, as ``bf''; 
  !
  !      all symbols after _ refer to the fermionic content of 
  !      on-shell lines; the rule is that outgoing fermion on-shell 
  !      line is ``bf'' and the incoming on-shell fermion line is ``f''
  !      for simplicity, f_bf current is called  f,
  !      and bf_f current is called bf
  !
  !     for currents with the W, notations are as above but there 
  !     is additional letter indicating presence of a corresponding 
  !     vector boson (e.g. fW ... )
  !
  !     for all currents with outgoing fermion lines, last few arguments
  !     of the current refer to the number of gluon lines to the left of 
  !     the first fermion line, the number of gluons between first and 
  !     second fermion lines,  the number of gluons between second and 
  !     third fermion lines, etc. For example, g_bff(...,ng1,ng2) 
  !     refers to the gluon current, where (clockwise) the first 
  !     on-shell fermion line is outgoing, the second fermion line is 
  !     incoming and the number of gluons between off-shell gluon and 
  !     the first fermion line is ng1 and the number of gluons between 
  !     the first and the second fermion lines is ng2
  !
  !     we also pass the information about fermion flavors to 
  !     the currents
  !=====================================================================

end module mprecurrence

