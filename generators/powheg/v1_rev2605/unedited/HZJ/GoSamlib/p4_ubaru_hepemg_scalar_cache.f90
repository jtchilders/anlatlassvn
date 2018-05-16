module     p4_ubaru_hepemg_scalar_cache
   use precision, only: ki_sam => ki
   use madds
   implicit none
   save

   private
!---#[ scalar integral cache for samurai:
    logical, public  :: samurai_cache_flag_g0
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_g0
    logical, public :: samurai_cache_flag_d541
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_d541
    logical, public  :: samurai_cache_flag_g1
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_g1
    logical, public :: samurai_cache_flag_d279
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d279
    logical, public :: samurai_cache_flag_d333
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d333
    logical, public :: samurai_cache_flag_d421
    complex(ki_sam), dimension(-2:0,cachedim2(1)), public :: samurai_cache_d421
    logical, public :: samurai_cache_flag_d537
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_d537
    logical, public :: samurai_cache_flag_d605
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d605
    logical, public  :: samurai_cache_flag_g2
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_g2
    logical, public :: samurai_cache_flag_d323
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d323
    logical, public :: samurai_cache_flag_d413
    complex(ki_sam), dimension(-2:0,cachedim2(1)), public :: samurai_cache_d413
    logical, public :: samurai_cache_flag_d533
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_d533
    logical, public :: samurai_cache_flag_d613
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d613
    logical, public  :: samurai_cache_flag_g3
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_g3
    logical, public :: samurai_cache_flag_d101
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_d101
    logical, public :: samurai_cache_flag_d553
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d553
    logical, public  :: samurai_cache_flag_g4
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_g4
    logical, public :: samurai_cache_flag_d77
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_d77
    logical, public  :: samurai_cache_flag_g5
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_g5
    logical, public :: samurai_cache_flag_d45
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_d45
    logical, public :: samurai_cache_flag_d277
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d277
!---#] scalar integral cache for samurai:

   public :: invalidate_cache
contains
   subroutine invalidate_cache()
      implicit none
      samurai_cache_flag_g0 = .false.
      samurai_cache_flag_d541 = .false.
      samurai_cache_flag_g1 = .false.
      samurai_cache_flag_d279 = .false.
      samurai_cache_flag_d333 = .false.
      samurai_cache_flag_d421 = .false.
      samurai_cache_flag_d537 = .false.
      samurai_cache_flag_d605 = .false.
      samurai_cache_flag_g2 = .false.
      samurai_cache_flag_d323 = .false.
      samurai_cache_flag_d413 = .false.
      samurai_cache_flag_d533 = .false.
      samurai_cache_flag_d613 = .false.
      samurai_cache_flag_g3 = .false.
      samurai_cache_flag_d101 = .false.
      samurai_cache_flag_d553 = .false.
      samurai_cache_flag_g4 = .false.
      samurai_cache_flag_d77 = .false.
      samurai_cache_flag_g5 = .false.
      samurai_cache_flag_d45 = .false.
      samurai_cache_flag_d277 = .false.
   end subroutine invalidate_cache
end module p4_ubaru_hepemg_scalar_cache