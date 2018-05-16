c          -*- fortran -*-
      real * 8 par_csicut
      parameter (par_csicut=1)
      integer par_maxseeds,par_maxxgriditerations
      parameter (par_maxseeds=200,par_maxxgriditerations=10)
      
      real * 8 par_diexp,par_dijexp,par_2gsupp,
     1         par_fsrtinycsi,par_fsrtinyy,
     2         par_isrtinycsi,par_isrtinyy
      common/pwhg_par/par_diexp,par_dijexp,par_2gsupp,
     1         par_fsrtinycsi,par_fsrtinyy,
     2         par_isrtinycsi,par_isrtinyy
      save /pwhg_par/
