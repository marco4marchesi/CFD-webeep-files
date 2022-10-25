!
      module grid_3d
!
!     +--------------------------------------------------------------------------------+
      use dynamic_vector
      use grid_1d

      implicit none

      type section_grid

        ! nodes
        integer :: nj_d, nj_b
        real(kind=8), dimeNsion(:,:), pointer :: rr
        integer,      dimeNsion(:),   pointer :: jd_jb, bound_p
        ! elements
        integer :: nm_d, nm_b
        integer,      dimeNsion(:),   pointer :: ele_type_d, ele_type_b, bound_m
        type (d_i_v), dimeNsion(:),   pointer :: j_m_d, j_m_b
        type (d_i_v), dimeNsion(:),   pointer :: ma_m_d, ma_m_b
         
      end type
!
      integer, public :: Ntetra, Nhexas, Nprism, Npyram
!
      real*8, dimension(:), pointer :: extrusion_path
!     +--------------------------------------------------------------------------------+
!
!
      contains
! 
! 
      subroutine extrude_path ( name )
!     +--------------------------------------------------------------------------------+
      use init_data, only : idf_cfg, read_record
      implicit none
!
      character(len=64), intent(in) :: name
      integer :: i
!     +--------------------------------------------------------------------------------+
!
      write (*,*) '   Extruding...'
!      call read_path_1D
      call read_record ( 'G_EXTRUSION_PARAMETERS' )
      call build_grid_1D( name, 1, extrusion_path )
!
      end subroutine extrude_path



       
      subroutine save_grid_3d(idf, name)
      !----------------------------------------------------------------------
      implicit none

      integer,           intent(in) :: idf
      character(len=64), intent(in) :: name
      !----------------------------------------------------------------------
       
      end subroutine save_grid_3d



 
      subroutine plot_grid_3d(idf, name)
      !----------------------------------------------------------------------
      implicit none

      integer,           intent(in) :: idf
      character(len=64), intent(in) :: name
      !----------------------------------------------------------------------
       
      end subroutine plot_grid_3d


      end module grid_3d
