!
      module smooth_bline
      implicit none
!
!
      contains
!
!
      subroutine compute_kernel( ksize, g )
!     +----------------------------------------------------------------+
!     Compute the Gaussian kernel of any size
!
      implicit none
!
      integer, intent(in) :: ksize
      real*8, dimension(-ksize:ksize), intent(out) :: g
!
      integer :: j
!
      real*8, dimension(-ksize:ksize) :: g0
      real*8 :: C
!     +----------------------------------------------------------------+
!
      g = 0.d0
      g0 = 0.d0
!
      if ( ksize .eq. 0 ) then 
        g0 = 1.d0
      else
        do j = -ksize, ksize
          g0(j) = exp( -real(j)*real(j) / (2*ksize**2) )
        enddo
      end if
!
      C = sum(g0)
      C = 1.d0 / C
!
      do j = -ksize, ksize
        g(j) = C * g0(j)
      enddo
!
      end subroutine compute_kernel
!
!
!
      subroutine convolution( kdGMax, nsl, signal )
!     +----------------------------------------------------------------+
      implicit none
!
      integer*4, intent(in) :: kdGMax, nsl
      real*8, dimension(:,:), intent(inout) :: signal
!
      real*8, dimension(:,:), allocatable :: g
      real*8, dimension(size(signal,2)) :: f
!
      real*8 :: C, gf
!
      integer*4 :: i, j, k, q, r
      integer*4 :: npnt, ksize
!     +----------------------------------------------------------------+
!
      npnt = size(signal,2)
!
      allocate (g(-kdGMax:kdGMax,npnt))
      g = 0.d0
!
!     +----------------------------------------------------------------+
!     | Compute the convolution                                        |
!     +----------------------------------------------------------------+
!
      do k = 1, nsl
        do q = 1, size(signal,1)
!
          f = signal(q,:)
          
          do i = 1, npnt
          
            if ( i .le. kdGMax ) then
              ksize = i - 1
            else if ( i .ge. npnt-kdGMax+1 ) then
              ksize = npnt - i
            else
              ksize = kdGMax
            end if
             
            call compute_kernel( ksize, g(-ksize:ksize,i) )
!
            gf = 0.d0
            do j = -ksize, ksize
              gf = gf + g(j,i) * f(i+j)
            enddo
            signal(q,i) = gf
!
          enddo
!
        enddo
      enddo
!
      end subroutine convolution
!
!
      end module smooth_bline
