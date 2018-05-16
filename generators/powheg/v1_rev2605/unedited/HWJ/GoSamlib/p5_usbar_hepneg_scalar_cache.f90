module     p5_usbar_hepneg_scalar_cache
   use precision, only: ki_sam => ki
   use madds
   implicit none
   save

   private
!---#[ scalar integral cache for samurai:
    logical, public  :: samurai_cache_flag_g0
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_g0
    logical, public :: samurai_cache_flag_d70
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d70
    logical, public  :: samurai_cache_flag_g1
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_g1
    logical, public :: samurai_cache_flag_d61
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_d61
    logical, public  :: samurai_cache_flag_g2
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_g2
    logical, public :: samurai_cache_flag_d34
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d34
    logical, public :: samurai_cache_flag_d40
    complex(ki_sam), dimension(-2:0,cachedim2(1)), public :: samurai_cache_d40
    logical, public :: samurai_cache_flag_d60
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_d60
    logical, public :: samurai_cache_flag_d92
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d92
    logical, public  :: samurai_cache_flag_g3
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_g3
    logical, public :: samurai_cache_flag_d30
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d30
    logical, public :: samurai_cache_flag_d37
    complex(ki_sam), dimension(-2:0,cachedim2(1)), public :: samurai_cache_d37
    logical, public :: samurai_cache_flag_d59
    complex(ki_sam), dimension(-2:0,cachedim4(1)), public :: samurai_cache_d59
    logical, public :: samurai_cache_flag_d94
    complex(ki_sam), dimension(-2:0,cachedim3(1)), public :: samurai_cache_d94
!---#] scalar integral cache for samurai:

   public :: invalidate_cache
contains
   subroutine invalidate_cache()
      implicit none
      samurai_cache_flag_g0 = .false.
      samurai_cache_flag_d70 = .false.
      samurai_cache_flag_g1 = .false.
      samurai_cache_flag_d61 = .false.
      samurai_cache_flag_g2 = .false.
      samurai_cache_flag_d34 = .false.
      samurai_cache_flag_d40 = .false.
      samurai_cache_flag_d60 = .false.
      samurai_cache_flag_d92 = .false.
      samurai_cache_flag_g3 = .false.
      samurai_cache_flag_d30 = .false.
      samurai_cache_flag_d37 = .false.
      samurai_cache_flag_d59 = .false.
      samurai_cache_flag_d94 = .false.
   end subroutine invalidate_cache
end module p5_usbar_hepneg_scalar_cache