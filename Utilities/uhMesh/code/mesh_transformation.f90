!
      module mesh_transformation
!
!     +--------------------------------------------------------------------------------+
      use init_data
      use nodes
!
      implicit none
!     +--------------------------------------------------------------------------------+
!
!
      contains
!
!
      subroutine scale_mesh
!     +-----------------------------------------------+
      implicit none
      integer :: i
!     +-----------------------------------------------+
!
      do i = 1, Np_d
        rr(1,i) = rr(1,i)*xfactor
        rr(2,i) = rr(2,i)*yfactor
        if ( size(rr,1) .eq. 3 ) rr(3,i) = rr(3,i)*zfactor
      end do
!
      return      
      end subroutine scale_mesh
!
      end module mesh_transformation
