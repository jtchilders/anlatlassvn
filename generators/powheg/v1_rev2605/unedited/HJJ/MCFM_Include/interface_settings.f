
      
      integer ret1,ret2,ret3,ret4
      logical efficient,params 
      integer n_particles
      logical read_ew,read_mv,read_mq,read_ml,read_ckm,read_as
      logical read_scheme
      logical inherit_as,ret_poles
      double precision inh_as 
      integer pid
      logical alphas_eq1,alphaew_eq1
      character*4 inscheme 
      
      common/inscheme/inscheme
      common/efficient/efficient
      common/pid/pid
      common/ME_ij/ret1,ret2,ret3,ret4
      common/n_particles/n_particles
      common/params/params
      common/read_p/read_ew,read_mv,read_mq,read_ml,read_ckm,read_as
     &     ,read_scheme
      common/alpha_eq1/alphas_eq1,alphaew_eq1
      common/inh_as/inherit_as,inh_as
      common/ret_poles/ret_poles
