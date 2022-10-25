!
      module export_module
!
!     +--------------------------------------------------------------------------------+
      implicit none
      private
!      
      integer, dimension(:), allocatable :: BC
      integer ::  extr_dir
      integer :: ziBC, zeBC
      real*8 :: zi, ze
      logical :: tbc
!
      integer, parameter :: XCRD = 1
      integer, parameter :: YCRD = 2
      integer, parameter :: ZCRD = 3      
!
      public :: save_npc_mesh, & 
                save_su2_mesh, &
                save_fns_mesh, &
		save_vtk_mesh, &
		save_foam_mesh
!     +--------------------------------------------------------------------------------+
!
!
      contains
!
!
      subroutine save_npc_mesh ( grid_name )
!     +-------------------------------------------------------------------------+
      use grid_1d
      use init_data
      use import_module
      use mesh_structure
      use nodes 
      use node_pair_structure
      use np_topology_gen
!
      implicit none
!
      character (len=64), intent(in)  :: grid_name
      integer :: name_length
!     +-------------------------------------------------------------------------+
!
      name_length = len_trim(grid_name)
      open (unit=11, file='nodes.'//trim(grid_name))
      call save_nodes ( 11, trim(grid_name), name_length )
      close (11)
!
      open (unit=11, file='grid.'//trim(grid_name))
      call save_mesh ( 11, trim(grid_name), name_length )
      close (11)
!
      return
      end subroutine save_npc_mesh
!
!
!
!
!
      subroutine save_fns_mesh ( grid_name )
!     +-------------------------------------------------------------------------+
      use grid_1d
      use init_data
      use import_module
      use mesh_structure
      use nodes 
      use node_pair_structure
      use np_topology_gen
!
      implicit none
!
      character(len=64), intent(in)  :: grid_name
!
      integer, dimension(:,:), allocatable :: fbound
      integer, dimension(:,:), allocatable :: jb_jd
      integer, dimension(:,:), allocatable :: bfaces
      integer, dimension(:), allocatable :: prd
      integer, dimension(:), allocatable :: etype, nds, bid
      integer, dimension(:), allocatable :: bpoints2D, bound_p2D
      integer, dimension(:), allocatable :: bpoints2Dt, bound_p2Dt
      integer, dimension(:), allocatable :: bpoints3D, bound_p3D
      integer, dimension(:), allocatable :: periodic
!
      integer, dimension(4,5) :: prism_faces
      integer, dimension(4,6) :: brick_faces
      integer, dimension(8) :: connect3D
      integer, dimension(4) :: connect
      integer, dimension(4) :: elgeom
      integer, dimension(3) :: adj_elm_i
!
      integer :: c, i, j, k, s, f, cpy, counter, Nbf, ndpel, bp, idx
      integer :: offset_prev, offset_next, NnodeP, index, cface, bcnt
      integer :: n1, n2, np, prev, Npext, pBC, Nbpnts, cntD, Nbl
!
      real*8, dimension(:,:), allocatable :: dx
      real*8, dimension(:), allocatable :: heigh
      real*8, dimension(:), pointer :: extrusion_path
      real*8, dimension(2) :: prv, nxt, avg
      real*8 :: off, extC
!
      character(len=64) :: fnsp_name
      logical, dimension(:), allocatable :: flag
!     +-------------------------------------------------------------------------+
!
      fnsp_name = trim(grid_name)//'.fns'
!
      Nbl = b_grid % number_of_edges          
      allocate ( BC(Nbl) )                    
      do i = 1, Nbl
        read (b_grid % edges(i) % edgename, '(i4)') BC(i)
      end do
!
      call set_path_1D
      call build_grid_1D( grid_name, 1, extrusion_path )
      call read_npc_mesh ( grid_name ) 
!
!     +------------------+
!     | Inverting jd_jb  |
!     +------------------+
!
      allocate ( jb_jd(2,Np_d) )
      jb_jd = 0
!
      do i = 1, Np_d
        bcnt = 1
        do j = 1, Np_b
!
          if ( jd_jb(j) .eq. i ) then
            jb_jd(bcnt,i) = j
            bcnt = bcnt + 1
            if ( bcnt .gt. 2 ) exit
          end if
!
        end do
      end do
!
!     +-------------------+
!     |  Mesh dimensions  |
!     +-------------------+
!
      Npext = size(extrusion_path)
      Nnodes = Np_d * Npext
      NnodeP = Np_d * (Npext-1)
      Nbpnts = Np_b/2 * Npext
      Nbf = 2*Ne_d + Ne_b*(Npext - 1)
!
      elgeom = 0
      elgeom(1) = count( ele_type_d == 3 )
      elgeom(3) = count( ele_type_d == 2 )
!
      write (*,*) '  - Nodes     ', Nnodes
      if ( count( ele_type_d == 3 ) .gt. 0 ) write (*,*) '  - Hexahedra ', elgeom(1) * (Npext-1)
      if ( count( ele_type_d == 2 ) .gt. 0 ) write (*,*) '  - Prisms    ', elgeom(3) * (Npext-1)
!
      open ( unit=11, file=trim(fnsp_name) )
!
      write (11,*) Nnodes, Nbf, 1
      write (11,*) count( elgeom .ne. 0 )
      do i = 1, 4
        if ( elgeom(i) .ne. 0 ) then
          write (11,*) i, elgeom(i) * ( Npext - 1 )
        end if
      end do

      write (11,*) 1.0, 1.0
      write (11,*) 'Converted UH-Mesh grid'        
!
!  write (11,*) 'variables = "x", "y", "z", "ib"'
!  write (11,*) 'zone N=',Nnodes,', E=', elgeom(3)*( Npext - 1 ),', ZONETYPE=FEBRICK, DATAPACKING=POINT'
!
!     +-----------------+
!     |  Writing nodes  |
!     +-----------------+
!
      allocate ( periodic(Np_d) )
      cntD = 0
      do s = 1, Npext
!
        extC = extrusion_path(s)
        do i = 1, Np_d
!
          select case ( extr_dir )
!
          case ( XCRD )
          
          pBC = 4100            
          if ( all(jb_jd(:,i) .ne. 0) ) then
            write (11,'(1x,3(2x,e22.16),2x,i9)') extC, rr(:,i), minval(BC(bound_p(jb_jd(:,i))))
          else
            if ( s .eq. 1 .or. s .eq. Npext ) then
              write (11,'(1x,3(2x,e22.16),2x,i9)') extC, rr(:,i), pBC
            else
              write (11,'(1x,3(2x,e22.16),2x,i9)') extC, rr(:,i), 0
            end if
          end if
!
          case ( YCRD )
!
          pBC = 4200            
          if ( all(jb_jd(:,i) .ne. 0) ) then
            write (11,'(1x,3(2x,e22.16),2x,i9)') rr(1,i), extC, rr(2,i), minval(BC(bound_p(jb_jd(:,i))))
          else
            if ( s .eq. 1 .or. s .eq. Npext ) then
              write (11,'(1x,3(2x,e22.16),2x,i9)') rr(1,i), extC, rr(2,i), pBC
            else
              write (11,'(1x,3(2x,e22.16),2x,i9)') rr(1,i), extC, rr(2,i), 0
            end if
          end if
!
          case ( ZCRD )
!
          pBC = 4300            
          if ( all(jb_jd(:,i) .ne. 0) ) then
!
            if ( s .eq. 1 ) then
              cntD = cntD + 1
              write (11,'(1x,3(2x,e22.16),2x,i9)') rr(:,i), extC, minval(BC(bound_p(jb_jd(:,i))))
              periodic(i) = cntD
            else if ( s .eq. Npext ) then
              cntD = cntD + 1
              write (11,'(1x,3(2x,e22.16),2x,i9)') rr(:,i), extC, -periodic(i)              
            else
              cntD = cntD + 1
              write (11,'(1x,3(2x,e22.16),2x,i9)') rr(:,i), extC, minval(BC(bound_p(jb_jd(:,i))))              
            end if
!
          else
!
            if ( s .eq. 1 ) then
              cntD = cntD + 1
              write (11,'(1x,3(2x,e22.16),2x,i9)') rr(:,i), extC, pBC
              periodic(i) = cntD
            else if ( s .eq. Npext ) then
              cntD = cntD + 1
              write (11,'(1x,3(2x,e22.16),2x,i9)') rr(:,i), extC, -periodic(i)
            else
              cntD = cntD + 1
              write (11,'(1x,3(2x,e22.16),2x,i9)') rr(:,i), extC, 0
            end if
!
          end if
!
          end select
!
        end do
      end do
!
!     +-------------------------------+
!     |  Boundary points for jbounds  |
!     +-------------------------------+
!
!        allocate ( bpoints2D(Np_b/2) )
!        bpoints2D = 0
!        allocate ( bound_p2D(Np_b/2) )
!        bound_p2D = 0
!
!        index = 0
!        do i = 1, Np_b
!          if ( .not.(any(bpoints2D .eq. jd_jb(i))) ) then
!            index = index + 1
!            bpoints2D(index) = jd_jb(i)
!            bound_p2D(index) = bound_p(i)
!          end if
!        end do
!
      allocate ( bpoints2Dt(Np_d) )
      bpoints2Dt = 0
      allocate ( bound_p2Dt(Np_d) )
      bound_p2Dt = 0
      
      do i = 1, Np_b
        bpoints2Dt(jd_jb(i)) = jd_jb(i)
        bound_p2Dt(jd_jb(i)) = max( BC(bound_p(i)), bound_p2Dt(jd_jb(i)) )
      end do

      allocate ( bpoints2D(Np_b/2) )
      bpoints2D = 0
      allocate ( bound_p2D(Np_b/2) )
      bound_p2D = 0

      index = 0
      do i = 1, Np_d
        if ( bpoints2Dt(i) .ne. 0 ) then
          index = index + 1
          bpoints2D(index) = bpoints2Dt(i)
          bound_p2D(index) = bound_p2Dt(i)
        end if
      end do

      deallocate ( bpoints2Dt )
      deallocate ( bound_p2Dt )

!
      open (unit=15, file='boundary_layer.dat')
!
!     +-----------------------------------+
!     |  Elements and Faces connectivity  |
!     +-----------------------------------+
!
      allocate ( fbound(3,Nbf) )
      fbound = 0
!
      brick_faces(:,1) = (/ 1,4,3,2 /)
      brick_faces(:,2) = (/ 1,5,8,4 /)
      brick_faces(:,3) = (/ 4,8,7,3 /)
      brick_faces(:,4) = (/ 2,3,7,6 /)
      brick_faces(:,5) = (/ 1,2,6,5 /)
      brick_faces(:,6) = (/ 5,6,7,8 /)
!
      prism_faces(:,1) = (/ 2,3,6,5 /)
      prism_faces(:,2) = (/ 3,1,4,6 /)
      prism_faces(:,3) = (/ 1,2,5,4 /)
      prism_faces(:,4) = (/ 1,3,2,2 /)
      prism_faces(:,5) = (/ 4,5,6,6 /)
!
      allocate ( etype(Ne_d*(Npext-1)) )
      etype = 0
      allocate ( bpoints3D(Nbpnts) )
      bpoints3D = 0
      allocate ( bound_p3D(Nnodes) )
      bound_p3D = 0
!
      offset_prev = 0
      offset_next = 0
      do s = 1, Npext
        bpoints3D(1+offset_next:Np_b/2+offset_next) = bpoints2D + offset_prev
        bound_p3D(bpoints3D(1+offset_next:Np_b/2+offset_next)) = bound_p2D
        offset_prev = offset_prev + Np_d
        offset_next = offset_next + Np_b/2
      end do
!
!     +-------------+
!     |  Hexahedra  |
!     +-------------+
!
      counter = 0
      cface = 0
!
      offset_prev = 0
      offset_next = Np_d
!
      do s = 1, Npext - 1
!
        extC = extrusion_path(s)
!
        if ( allocated(nds) ) deallocate (nds)
        allocate ( nds(4) )
        nds = 0
!
        do i = 1, Ne_d
          if ( ele_type_d(i) .eq. 3 ) then
!
            counter = counter + 1
            connect = j_m_d(i) % vec
!
!    +------------------------------------+
!    |  Sorting nodes in clockwise order  |
!    +------------------------------------+
!
            cpy = connect(1)
            connect(1) = connect(4)
            connect(4) = cpy
            cpy = connect(3)
            connect(3) = connect(2)
            connect(2) = cpy
!
!     +-------------------------+
!     |  Elements connectivity  |
!     +-------------------------+
!
            connect3D(1:4) = connect + offset_prev
            connect3D(5:8) = connect + offset_next
!
            write (11,'(1x,8(1x,i8))') connect3D
            write (15,*) counter
!
!     +------------------------+
!     |  Jbounds connectivity  |
!     +------------------------+
!
            do f = 1, 6
!
              nds = connect3D(brick_faces(:,f))
              if ( all(nds .le. Np_d) .or. all(nds .gt. NnodeP) ) then
!
                cface = cface + 1
                fbound(1,cface) = pBC
                fbound(2,cface) = f
                fbound(3,cface) = counter
!
              else
!
                bp = 0
                do k = 1, 4
                  if ( any(bpoints3D .eq. nds(k)) ) bp = bp + 1
                end do
!
                if ( bp .eq. 4 ) then
                  cface = cface + 1
                  fbound(1,cface) = minval( bound_p3D(nds) )
                  fbound(2,cface) = f
                  fbound(3,cface) = counter
                end if
!
              end if
            end do
          end if
        end do
!
        offset_prev = offset_prev + Np_d
        offset_next = offset_next + Np_d
!
      end do
!
!     +----------+
!     |  Prisms  |
!     +----------+
!
      offset_prev = 0
      offset_next = Np_d
!
      do s = 1, Npext - 1
!
        extC = extrusion_path(s)
        prev = counter
!
        do i = 1, Ne_d
          if ( ele_type_d(i) .eq. 2 ) then

            counter = counter + 1
            connect(1:3) = j_m_d(i) % vec              
!
!    +------------------------------------+
!    |  Sorting nodes in clockwise order  |
!    +------------------------------------+
!
           cpy = connect(1)
           connect(1) = connect(3)
           connect(3) = cpy
           connect(4) = cpy
!
           connect3D(1:3) = connect(1:3) + offset_prev
           connect3D(4:6) = connect(1:3) + offset_next
           connect3D(7) = 0
           connect3D(8) = 0
!
           write (11,'(1x,8(1x,i9))') connect3D(1:6)
!
!     +------------------------+
!     |  Jbounds connectivity  |
!     +------------------------+
!
            do f = 1, 5
!
              if ( f .le. 3 ) ndpel = 4
              if ( f .gt. 3 ) ndpel = 3
!               
              if ( allocated(nds) ) deallocate (nds)
              allocate ( nds(ndpel) )
              nds = connect3D(prism_faces(1:ndpel,f))                
!
              if ( all(nds .le. Np_d) .or. all(nds .gt. NnodeP) ) then
!
                cface = cface + 1
                fbound(1,cface) = pBC
                fbound(2,cface) = f
                fbound(3,cface) = counter
!
              else 
!
                  allocate ( bfaces(ndpel,2) )
                  bfaces = 0
                  bp = 0
!
                  idx = 1
                  do k = 1, ndpel
                    if ( any(bpoints3D .eq. nds(k)) ) then
!
                      bp = bp + 1
                      bfaces(k,1) = bound_p3D(nds(k))
!
                    end if
                  end do
!
                  do k = 1, size(bfaces,1)
                    bfaces(k,2) = count( bfaces(:,1) .eq. bfaces(k,1) )
                  end do
!
                  if ( bp .eq. ndpel ) then
                    cface = cface + 1
                    fbound(1,cface) = minval( bound_p3D(nds) )
                    fbound(2,cface) = f
                    fbound(3,cface) = counter
                  end if
!
                  deallocate ( bfaces )
!
              end if
            end do
!
!     +------------------------------------------------------+
!     |  Check if the triangle is inside the boundary layer  |
!     +------------------------------------------------------+
!
 !            adj_elm_i = ma_m_d(i) % vec
!              do j = 1, 3
!  !
!                if ( adj_elm_i(j) .eq. 0 ) cycle 
!                if ( ele_type_d(adj_elm_i(j)) .eq. 3 ) then
!  !
!                  do k = 1, 3
!                    if ( ele_type_d(adj_elm_i(k)) .eq. 2 ) then
!                    
!                      if ( any( ele_type_d(ma_m_d(adj_elm_i(k)) % vec) .eq. 3 ) ) then
!                        if ( counter .eq. prev+1 ) then
!                          write (15,*) counter
!                          prev = counter
!                          exit
!                        end if
!                      end if
!                    end if
!                  end do
!  !
!                end if
!              end do

          end if
        end do
!
        offset_prev = offset_prev + Np_d
        offset_next = offset_next + Np_d
!
      end do        
!
      do i = 1, Nbf
        write (11,*) fbound(:,i)
      end do
!
      return
      end subroutine save_fns_mesh
!
!
!
!
!
      subroutine save_su2_mesh ( grid_name, mtype )
!     +-------------------------------------------------------------------------+
      use import_module
      use init_data
      use nodes 
      use mesh_structure
      use node_pair_structure
      use np_topology_gen
      use grid_1d
      use grid_3d      
      use grid_types
!
      implicit none
!
      character(len=64), intent(in)  :: grid_name
      integer, intent(in) :: mtype
      character(len=64) :: su2_name
!
      integer, dimension(:,:), allocatable :: jb_jd
      integer, dimension(2) :: elgeom
!
      integer :: c, i, j, k, s, f, Nbm, boundm_idx
      integer :: bcnt, offset, Nbf, NnodeP
      integer :: Nbl, Nbpnts, cntD, Npext
!
      integer, dimension(3), parameter :: etyp2D = (/ 3, 5, 9 /)
      integer, dimension(3), parameter :: etyp3D = (/ 0, 13, 12 /)
      integer, parameter :: XCRD = 1
      integer, parameter :: YCRD = 2
      integer, parameter :: ZCRD = 3
!
      real*8 :: extC      
!     +-------------------------------------------------------------------------+
!
      su2_name = trim(grid_name)//'.su2'
      Nbl = b_grid % number_of_edges
!
      if ( .not. allocated(jd_jb) ) then
        call read_npc_mesh ( grid_name )
      end if
!
!     +------------------+
!     | Inverting jd_jb  |
!     +------------------+
!
      allocate ( jb_jd(2,Np_d) )
      jb_jd = 0      
!
      do i = 1, Np_d      
        bcnt = 1
        do j = 1, Np_b
!
          if ( jd_jb(j) .eq. i ) then
            jb_jd(bcnt,i) = j
            bcnt = bcnt + 1
            if ( bcnt .gt. 2 ) exit
          end if
!
        end do
      end do
!
      do j = 1, Ne_b
        allocate ( j_m_b(j)%ivec(2) )
        j_m_b(j)%ivec(1) = j_m_b(j)%vec(2)
        j_m_b(j)%ivec(2) = j_m_b(j)%vec(1)
      end do      
!
!
      open ( unit=11, file=trim(su2_name) )
!        
      select case (mtype)      
!
      case (2)
!
!     +-------------------+
!     |  Mesh dimensions  |
!     +-------------------+
!
        Nnodes = Np_d
        Nbpnts = Np_b/2
        Nbm = maxval(bound_m)
!
        elgeom = 0
        elgeom(1) = count( ele_type_d == 2 )
        elgeom(2) = count( ele_type_d == 3 )
!
        write (*,*) '   - Nodes          ', Nnodes
        if ( count( ele_type_d == 3 ) .gt. 0 ) write (*,*) '   - Quadrilaterals ', elgeom(2)
        if ( count( ele_type_d == 2 ) .gt. 0 ) write (*,*) '   - Triangles      ', elgeom(1)
!
        write (11,'(a6,1x,i1)') 'NDIME=', meshtype
!
!     +------------------------+
!     |  Writing connectivity  |
!     +------------------------+
!
        write (11,'(a6,1x,i10)') 'NELEM=', Ne_d
        do i = 1, Ne_d
          write (11,*) etyp2D(ele_type_d(i)), j_m_d(i)%vec-1, i-1 
        end do
!
!     +-----------------+
!     |  Writing nodes  |
!     +-----------------+
!
        write (11,'(a6,1x,i10)') 'NPOIN=', Np_d
        do i = 1, Np_d
          write (11,*) rr(:,i), i-1
        end do
!
!     +----------------------+
!     |  Writing boundaries  |
!     +----------------------+
!
        write (11,'(a6,1x,i10)') 'NMARK=', Nbm
!
        do i = 1, Nbm
!
          if ( associated( b_grid%edges )) then
            write (11,'(a11,1x,a)') 'MARKER_TAG=', trim(b_grid%edges(i)%edgename)
          else
            write (11,'(a11,1x,i2)') 'MARKER_TAG=', i
          end if
          write (11,'(a13,1x,i10)') 'MARKER_ELEMS=', count(bound_m .eq. i)
!
          boundm_idx = 0
          do j = 1, Ne_b
            if ( bound_m(j) .eq. i ) then
              write (11,*) 3, jd_jb(j_m_b(j)%vec)-1
              boundm_idx = boundm_idx + 1
            end if
            if ( boundm_idx .eq. count(bound_m .eq. i) ) exit
          end do
!
        end do
!
!
      case (3)
!      
        Npext = size(extrusion_path)
        Nnodes = Np_d*Npext
        NnodeP = Np_d * (Npext-1)
        Nbpnts = Np_b/2 * Npext
        Nbf = 2*Ne_d + Ne_b*(Npext - 1) 
        Nbm = maxval(bound_m)
!
        elgeom = 0
        elgeom(1) = count( ele_type_d == 2 )
        elgeom(2) = count( ele_type_d == 3 )
!
        write (*,*) '   - Nodes     ', Nnodes
        if ( count( ele_type_d == 3 ) .gt. 0 ) write (*,*) '   - Hexahedra ', elgeom(2)*(Npext-1)
        if ( count( ele_type_d == 2 ) .gt. 0 ) write (*,*) '   - Prisms    ', elgeom(1)*(Npext-1)
!
        write (11,'(a6,1x,i1)') 'NDIME=', meshtype
!
!     +------------------------+
!     |  Writing connectivity  |
!     +------------------------+
!
        write (11,'(a6,1x,i10)') 'NELEM=', Ne_d*(Npext-1)
        cntD = 0
        offset = 0        
        do s = 1, Npext-1
          do i = 1, Ne_d
            write (11,*) etyp3D(ele_type_d(i)), offset+j_m_d(i)%vec-1, offset+j_m_d(i)%vec-1+Np_d, cntD
            cntD = cntD + 1 
          end do
          offset = offset + Np_d
        end do
!
!     +-----------------+
!     |  Writing nodes  |
!     +-----------------+
!
        write (11,'(a6,1x,i10)') 'NPOIN=', Nnodes
        cntD = 0
        do s = 1, Npext
          extC = extrusion_path(s)
          do i = 1, Np_d            
            write (11,'(1x,3(2x,e22.16),2x,i9)') rr(:,i), extC, cntD
            cntD = cntD + 1
          end do
        end do
!
!     +----------------------+
!     |  Writing boundaries  |
!     +----------------------+
!
        write (11,'(a6,1x,i10)') 'NMARK=', Nbm+2
!
        do i = 1, Nbm
!
          if ( associated( b_grid%edges )) then
            write (11,'(a11,1x,a)') 'MARKER_TAG=', trim(b_grid%edges(i)%edgename)
          else
            write (11,'(a11,1x,i2)') 'MARKER_TAG=', i
          end if
          write (11,'(a13,1x,i10)') 'MARKER_ELEMS=', count(bound_m .eq. i)*(Npext-1)
!
          boundm_idx = 0
          do j = 1, Ne_b
            if ( bound_m(j) .eq. i ) then
              offset = 0
              do s = 1, Npext-1
                write (11,*) 9, offset+jd_jb(j_m_b(j)%vec)-1, offset+jd_jb(j_m_b(j)%ivec)-1+Np_d
                offset = offset + Np_d
              end do
              boundm_idx = boundm_idx + 1
            end if
            if ( boundm_idx .eq. count(bound_m .eq. i) ) exit
          end do
!
        end do
!
        write (11,'(a11,1x,a)') 'MARKER_TAG=', trim(floor_tag)
        write (11,'(a13,1x,i10)') 'MARKER_ELEMS=', Ne_d
        do j = 1, Ne_d
          if ( ele_type_d(j) == 2 ) write(11,*) 5, j_m_d(j)%vec-1
          if ( ele_type_d(j) == 3 ) write(11,*) 9, j_m_d(j)%vec-1          
        end do
!
        write (11,'(a11,1x,a)') 'MARKER_TAG=', trim(roof_tag)
        write (11,'(a13,1x,i10)') 'MARKER_ELEMS=', Ne_d
        offset = (Npext-1)*Np_d
        do j = 1, Ne_d
          if ( ele_type_d(j) == 2 ) write(11,*) 5, offset+j_m_d(j)%vec-1
          if ( ele_type_d(j) == 3 ) write(11,*) 9, offset+j_m_d(j)%vec-1          
        end do
!
      end select
      return
      end subroutine save_su2_mesh
!
!
!
!
!
      subroutine save_vtk_mesh ( grid_name )
!     +-------------------------------------------------------------------------+
      use import_module
      use init_data
      use nodes 
      use mesh_structure
      use node_pair_structure
      use np_topology_gen
      use grid_1d
      use grid_types
      use arrays, only: b_jjdir
!
      implicit none
!
      character(len=64), intent(in)  :: grid_name
      character(len=64) :: vtk_name
!
      integer, dimension(:,:), allocatable :: jb_jd
      integer, dimension(2) :: elgeom
!
      integer :: c, i, j, k, s, f, Nbm, boundm_idx
      integer :: bcnt
      integer :: Nbl, Nbpnts, Ncel, Ncsi
!
      integer, dimension(3), parameter :: etyp = (/ 3, 5, 9 /)
      integer, parameter :: XCRD = 1
      integer, parameter :: YCRD = 2
      integer, parameter :: ZCRD = 3
!     +-------------------------------------------------------------------------+
!
      vtk_name = trim(grid_name)//'.vtk'
      Nbl = b_grid % number_of_edges
!
      if ( .not. allocated(jd_jb) ) then
        call read_npc_mesh ( grid_name )
      end if
!
!     +------------------+
!     | Inverting jd_jb  |
!     +------------------+
!
      allocate ( jb_jd(2,Np_d) )
      jb_jd = 0      
!
      do i = 1, Np_d      
        bcnt = 1
        do j = 1, Np_b
!
          if ( jd_jb(j) .eq. i ) then
            jb_jd(bcnt,i) = j
            bcnt = bcnt + 1
            if ( bcnt .gt. 2 ) exit
          end if
!
        end do
      end do
!
!     +-------------------+
!     |  Mesh dimensions  |
!     +-------------------+
!
      Nnodes = Np_d
      Nbpnts = Np_b/2
      Nbm = maxval(bound_m)
!
      elgeom = 0
      elgeom(1) = count( ele_type_d == 2 )
      elgeom(2) = count( ele_type_d == 3 )
!
      write (*,*) '   - Nodes          ', 2*Nnodes
      if ( count( ele_type_d == 3 ) .gt. 0 ) write (*,*) '   - Quadrilaterals ', elgeom(2)
      if ( count( ele_type_d == 2 ) .gt. 0 ) write (*,*) '   - Trianlges      ', elgeom(1)
!
      open ( unit=11, file=trim(vtk_name) )
      write (11,'(a26)') '# vtk DataFile Version 3.0'
      write (11,'(a25)') 'Unstructured Grid Example'
      write (11,'(a5)') 'ASCII'
      write (11,'(a25)') 'DATASET UNSTRUCTURED_GRID'
!
!     +-----------------------------+
!     |  Writing nodes coordinates  |
!     +-----------------------------+
!
      write (11,*) 'POINTS', Nnodes, 'float'
      do i = 1, Nnodes
        write (11,*) rr(:,i), 0
      end do
!
!     +-----------------+
!     |  Writing cells  |
!     +-----------------+
!
      Ncel = elgeom(1) + elgeom(2)
      Ncsi = elgeom(1)*4 + elgeom(2)*5
      write (11,*) 'CELLS', Ncel, Ncsi
      do i = 1, Ne_d
        write (11,*) size(j_m_d(i)%vec), j_m_d(i)%vec - 1  
      end do
!
!     +---------------------+
!     |  Writing cell types |
!     +---------------------+
!
      write (11,*) 'CELL_TYPES', Ne_d
      do i = 1, Ne_d
        if ( ele_type_d(i) == 2 ) write (11,*) 5
        if ( ele_type_d(i) == 3 ) write (11,*) 9	
      end do
!
!
      return
      end subroutine save_vtk_mesh
!
!
!
!
!
      subroutine save_foam_mesh ( grid_name )
!     +-------------------------------------------------------------------------+
      use grid_1d
      use init_data
      use import_module
      use mesh_structure
      use nodes 
      use node_pair_structure
      use np_topology_gen
!
      implicit none
!
      character(len=64), intent(in)  :: grid_name
!
      type(d_i_v), dimension(:), pointer ::  m_c
      type(d_i_v), dimension(:), pointer ::  faces
      type(d_i_v), dimension(:), pointer ::  m_f
      type(d_i_v), dimension(:), pointer ::  f_m      
      integer, dimension(:,:), allocatable :: fbound
      integer, dimension(:,:), allocatable :: jb_jd
      integer, dimension(:,:), allocatable :: bfaces, edgesE
      integer, dimension(:), allocatable :: prd
      integer, dimension(:), allocatable :: etype, nds, bid, nodesE
      integer, dimension(:), allocatable :: bpoints2D, bound_p2D
      integer, dimension(:), allocatable :: bpoints2Dt, bound_p2Dt
      integer, dimension(:), allocatable :: bpoints3D, bound_p3D, locfac
      integer, dimension(:), allocatable :: periodic, NbfBid, srtFac
!
      integer, dimension(6) :: facemapPr
      integer, dimension(8) :: facemapHe
      integer, dimension(4,5) :: prism_faces
      integer, dimension(4,6) :: brick_faces
      integer, dimension(8) :: connect3D
      integer, dimension(4) :: connect
      integer, dimension(4) :: elgeom
      integer, dimension(2) :: nfac
!
      integer :: c, i, j, k, s, f, cpy, counter, Nbf, ndpel, bp, idx, b, cnt, m, cd, npf
      integer :: offset_prev, offset_next, NnodeP, index, cface, bcnt, nb1, nb2, fid
      integer :: n1, n2, n3, n4, np, prev, Npext, pBC, Nbpnts, cntD, Nbl, Ntf, fcnt
!
      real*8, dimension(2,2) :: rot_mat
      real*8, dimension(:,:), allocatable :: dx, rrR
      real*8, dimension(:), allocatable :: heigh
      real*8, dimension(:), pointer :: extrusion_path
      real*8, dimension(2) :: prv, nxt, avg, ed, nr
      real*8 :: off, extC, xcg, ycg, zcg, alpha
      real*8 :: pi = 3.14159265358979323
!
      character(len=64) :: foamfolder
      character(len=128) :: cmndfolder
      character(len=64) :: foamfakpat
      character(len=3) :: typel, facet
      character(len=1) :: quote
!
      logical, dimension(:), allocatable :: flag
      logical dir_e      
!     +-------------------------------------------------------------------------+
!
      quote = '"'
      foamfolder = trim(grid_name)//'-foam'
      foamfakpat = 'constant/polyMesh'
!
      inquire( file=foamfolder, exist=dir_e )
      if ( .not. dir_e ) then
        cmndfolder = ' mkdir '//trim(foamfolder)//';'// &
	             ' cd '//trim(foamfolder)//';'// &
		     ' mkdir constant;'// &
		     ' cd constant;'// &		     
		     ' mkdir polyMesh'
        call system(cmndfolder)
      end if
!
      call set_path_1D
      call build_grid_1D( grid_name, 1, extrusion_path )
      call read_npc_mesh ( grid_name )
!
!     +------------------+
!     | Inverting jd_jb  |
!     +------------------+
!
      allocate ( jb_jd(2,Np_d) )
      jb_jd = 0
!
      do i = 1, Np_d
        bcnt = 1
        do j = 1, Np_b
!
          if ( jd_jb(j) .eq. i ) then
            jb_jd(bcnt,i) = j
            bcnt = bcnt + 1
            if ( bcnt .gt. 2 ) exit
          end if
!
        end do
      end do
!
!     +-------------------+
!     |  Mesh dimensions  |
!     +-------------------+
!
      Npext = size(extrusion_path)
      Nnodes = Np_d * Npext
      NnodeP = Np_d * (Npext-1)
      Nbpnts = Np_b/2 * Npext
      Nbf = 2*Ne_d + Ne_b*(Npext - 1)
      Ntf = Nbf + (Nc_fv - Nc_b)
!
      elgeom = 0
      elgeom(1) = count( ele_type_d == 3 )
      elgeom(3) = count( ele_type_d == 2 )
!
      write (*,*) '   - Nodes     ', Nnodes
      if ( count( ele_type_d == 3 ) .gt. 0 ) write (*,*) '   - Hexahedra ', elgeom(1) * (Npext-1)
      if ( count( ele_type_d == 2 ) .gt. 0 ) write (*,*) '   - Prisms    ', elgeom(3) * (Npext-1)
      write (*,*) ''
      write (*,*) 'NOTE: Before using the openFoam mesh, please run the Foam utilities'
      write (*,*) 'checkMesh and renumberMesh to verify the mesh and ensure upper '
      write (*,*) 'triangular order for the faces. Also set writePrecision to 16.'
!
      open ( unit=11, file=trim(foamfolder)//'/'//trim(foamfakpat)//'/points' )
      call write_foam_header ( 11 )
      write (11,*) '    class       vectorField;'
      write (11,*) '    location    '//quote//trim(foamfakpat)//quote//';'
      write (11,*) '    object      points;'
      write (11,*) '}'
      call write_foam_separator ( 11 )
      write (11,*) ''
      write (11,*) Nnodes
      write (11,*) '('
      do i = 1, Np_d
        write (11,*) '(', rr(:,i), 1.0, ')'
      end do
      do i = 1, Np_d
        write (11,*) '(', rr(:,i), 0.0, ')'
      end do
!
      write (11,*) ')'
      write (11,*) ''
      call write_foam_separator ( 11 )
      close (11)
!
      allocate ( faces(0:Ntf-1) )
!
      open ( unit=11, file=trim(foamfolder)//'/'//trim(foamfakpat)//'/faces' )
      call write_foam_header ( 11 )
      write (11,*) '    class       faceList;'
      write (11,*) '    location    '//quote//trim(foamfakpat)//quote//';'
      write (11,*) '    object      faces;'
      write (11,*) '}'
      call write_foam_separator ( 11 )
      write (11,*) ''
      write (11,*) Ntf
      write (11,*) '('
!
!     +------------------+
!     |  Internal faces  |
!     +------------------+
!
      fcnt = -1
      do c = 1, Nc_fv
!
        if ( any( jb_jd(:,j_c_fv(1:2,c)) .eq. 0 ) ) then
	  n1 = j_c_fv(1,c) - 1
	  n2 = j_c_fv(2,c) - 1
	  n3 = n1 + Np_d
	  n4 = n2 + Np_d
	  fcnt = fcnt + 1
	  allocate ( faces(fcnt)%vec(4) )
	  faces(fcnt)%vec = (/ n1, n2, n4, n3 /)
	end if
!      
      end do
!
!     +--------------------------------------+
!     |  Boundary faces along 2D boundaries  |
!     +--------------------------------------+
!
      Nbl = maxval(bound_c)
      allocate ( NbfBid(Nbl+1) )
!
      do b = 1, Nbl
        cnt = 0
        do c = 1, Nc_b
	  if ( bound_c(c) .eq. b ) then
	    nb1 = j_c_b(1,c)
	    nb2 = j_c_b(2,c)
	    n1 = jd_jb(nb1) - 1
	    n2 = jd_jb(nb2) - 1
	    n3 = n1 + Np_d
	    n4 = n2 + Np_d
	    fcnt = fcnt + 1
	    allocate ( faces(fcnt)%vec(4) )
	    faces(fcnt)%vec = (/ n1, n2, n4, n3 /)
	    cnt = cnt + 1
	  else
	    cycle
	  end if
	end do
	NbfBid(b) = cnt
      end do
!
!     +----------------------------------------+
!     |  Boundary faces along symmetry planes  |
!     +----------------------------------------+
!
      do m = 1, Ne_d
!
        if ( ele_type_d(m) .eq. 2 ) then
!
	  n1 = j_m_d(m) % vec(1) - 1
	  n2 = j_m_d(m) % vec(2) - 1
	  n3 = j_m_d(m) % vec(3) - 1
!	  
	  fcnt = fcnt + 1
	  allocate ( faces(fcnt)%vec(3) )
	  faces(fcnt)%vec = (/ n1, n2, n3 /)
!
	else if ( ele_type_d(m) .eq. 3 ) then
!
	  n1 = j_m_d(m) % vec(1) - 1
	  n2 = j_m_d(m) % vec(2) - 1
	  n3 = j_m_d(m) % vec(3) - 1	  
	  n4 = j_m_d(m) % vec(4) - 1	
!	    
	  fcnt = fcnt + 1
	  allocate ( faces(fcnt)%vec(4) )
	  faces(fcnt)%vec = (/ n1, n2, n3, n4 /)
!
	end if
!
      end do
!
      do m = 1, Ne_d
!
        if ( ele_type_d(m) .eq. 2 ) then
!
	  n1 = j_m_d(m) % vec(1) - 1
	  n2 = j_m_d(m) % vec(2) - 1
	  n3 = j_m_d(m) % vec(3) - 1	  
!  
	  fcnt = fcnt + 1
	  allocate ( faces(fcnt)%vec(3) )	  
	  faces(fcnt)%vec = (/ n1+Np_d, n3+Np_d, n2+Np_d /)
!
	else if ( ele_type_d(m) .eq. 3 ) then
!
	  n1 = j_m_d(m) % vec(1) - 1
	  n2 = j_m_d(m) % vec(2) - 1
	  n3 = j_m_d(m) % vec(3) - 1	  
	  n4 = j_m_d(m) % vec(4) - 1
!	  
	  fcnt = fcnt + 1
	  allocate ( faces(fcnt)%vec(4) )	  
	  faces(fcnt)%vec = (/ n1+Np_d, n4+Np_d, n3+Np_d, n2+Np_d /)
!
	end if
!
      end do
      NbfBid(Nbl+1) = 2*Ne_d
!
      write (11,*) ')'
      write (11,*) ''
      call write_foam_separator ( 11 )
      close (11)      
!
      open ( unit=11, file=trim(foamfolder)//'/'//trim(foamfakpat)//'/boundary' )
      call write_foam_header ( 11 )
      write (11,*) '    class       polyBoundaryMesh;'
      write (11,*) '    location    '//quote//trim(foamfakpat)//quote//';'
      write (11,*) '    object      boundary;'
      write (11,*) '}'
      call write_foam_separator ( 11 )
      write (11,*) ''
      write (11,*) Nbl + 1
      write (11,*) '('
!
      allocate ( srtFac(Nbl+1) )
      srtFac = Ntf - Nbf
!     
      do b = 2, Nbl
        srtFac(b) = srtFac(b-1) + NbfBid(b-1)
      end do
!
      do b = 1, Nbl
        write (11,*) b_grid % edges(b) % edgename
	write (11,*) '{'
	write (11,*) '     type        patch;'
	write (11,*) '     nFaces    ', NbfBid(b), ';'
	write (11,*) '     startFace ', srtFac(b), ';'
	write (11,*) '}'
      end do
!
      write (11,*) 'defaultFaces'
      write (11,*) '{'
      write (11,*) '     type	     empty;'
      write (11,*) '     inGroups    1(empty);'
      write (11,*) '     nFaces    ', 2*Ne_d, ';'
      write (11,*) '     startFace ', srtFac(Nbl) + NbfBid(Nbl), ';'
      write (11,*) '}'
!
      write (11,*) ')'
      write (11,*) ''
      call write_foam_separator ( 11 )
      close (11)      
!
      open ( unit=11, file=trim(foamfolder)//'/'//trim(foamfakpat)//'/owner' )
      call write_foam_header ( 11 )
      write (11,*) '    class       labelList;'
      write (11,*) '    note        '//quote//'nPoints:', Nnodes, 'nCells:', Ne_d, &
        'nFaces:', Ntf, 'nInternalFaces:', Ntf - Nbf, quote//';'
      write (11,*) '    location    '//quote//trim(foamfakpat)//quote//';'
      write (11,*) '    object      owner;'
      write (11,*) '}'
      call write_foam_separator ( 11 )
      write (11,*) ''
      write (11,*) Ntf
      write (11,*) '('
!
      cd_cb = cd_cb_connectivity( j_c_fv, jd_jb,  j_c_b )
      allocate ( m_c(size_DIV(c_m_fv,3)) )
      m_c = invert_DIV ( c_m_fv )
      allocate ( m_f(0:Ntf-1) )
!
      fcnt = -1
      do c = 1, Nc_fv      
        if ( any( jb_jd(:,j_c_fv(1:2,c)) .eq. 0 ) ) then
	  write (11,*) minval(m_c(c)%vec) - 1
	  fcnt = fcnt + 1
	  allocate ( m_f(fcnt)%vec(1) )
	  m_f(fcnt)%vec = minval(m_c(c)%vec) - 1
	end if
      end do
!
      do b = 1, Nbl
         do c = 1, Nc_b
           if ( bound_c(c) .eq. b ) then
             cd = cd_cb(c)
             write (11,*) minval(m_c(cd)%vec) - 1
	     fcnt = fcnt + 1
	     allocate ( m_f(fcnt)%vec(1) )
	     m_f(fcnt)%vec = minval(m_c(cd)%vec) - 1
           end if
         end do
       end do
!
      do m = 1, Ne_d
        write (11,*) m - 1
	fcnt = fcnt + 1
	allocate ( m_f(fcnt)%vec(1) )
	m_f(fcnt)%vec = m - 1
      end do
      do m = 1, Ne_d
        write (11,*) m - 1
	fcnt = fcnt + 1
	allocate ( m_f(fcnt)%vec(1) )
	m_f(fcnt)%vec = m - 1
      end do
!
      write (11,*) ')'
      write (11,*) ''
      call write_foam_separator ( 11 )
      close (11)
!
      do i = 0, Ntf-1
        m_f(i)%vec = m_f(i)%vec + 1
      end do
      allocate ( f_m(size_DIV(m_f,3)) )
      f_m = invert_DIV ( m_f )
!
      open ( unit=11, file=trim(foamfolder)//'/'//trim(foamfakpat)//'/neighbour' )
      call write_foam_header ( 11 )
      write (11,*) '    class       labelList;'
      write (11,*) '    note        '//quote//'nPoints:', Nnodes, 'nCells:', Ne_d, &
        'nFaces:', Ntf, 'nInternalFaces:', Ntf - Nbf, quote//';'
      write (11,*) '    location    '//quote//trim(foamfakpat)//quote//';'
      write (11,*) '    object      neighbour;'
      write (11,*) '}'
      call write_foam_separator ( 11 )
      write (11,*) ''
      write (11,*) Ntf - Nbf
      write (11,*) '('
!
      do c = 1, Nc_fv      
        if ( any( jb_jd(:,j_c_fv(1:2,c)) .eq. 0 ) ) then
	  write (11,*) maxval(m_c(c)%vec) - 1
	end if
      end do
!
      write (11,*) ')'
      call write_foam_separator ( 11 )
      close (11)
!
      open ( unit=11, file=trim(foamfolder)//'/'//trim(foamfakpat)//'/faces' )
      call write_foam_header ( 11 )
      write (11,*) '    class       faceList;'
      write (11,*) '    location    '//quote//trim(foamfakpat)//quote//';'
      write (11,*) '    object      faces;'
      write (11,*) '}'
      call write_foam_separator ( 11 )
      write (11,*) ''
      write (11,*) Ntf
      write (11,*) '('
!
!     +-------------------------------------------------------------------+
!     |  Reordering face nodes to respect the normal leaving the element  |
!     +-------------------------------------------------------------------+
!
      do m = 1, size(f_m)      
!
        allocate ( locfac(size(f_m(m)%vec)) )
	locfac = f_m(m)%vec
        do f = 1, size(locfac)
	  f_m(m)%vec(f) = maxval(locfac)
	  locfac(maxloc(locfac)) = -1
	end do      
!
        fid = maxval(f_m(m)%vec) - 1
        npf = size(faces(fid)%vec)
!
	allocate ( nodesE(0:2*npf-1) )
	allocate ( edgesE(2,npf) )
	allocate ( rrR(2,3) )
!
        nodesE(npf:2*npf-1) = faces(fid)%vec
	nodesE(0:npf-1) = nodesE(npf:2*npf-1) - Np_d
!
	do i = 1, npf-1 	
	  edgesE(1,i) = nodesE(i-1)
          edgesE(2,i) = nodesE(i)	
	end do
	edgesE(1,npf) = nodesE(npf-1)
	edgesE(2,npf) = nodesE(0)
!
        xcg = sum( rr(1,nodesE(0:npf-1)+1) ) / npf
	ycg = sum( rr(2,nodesE(0:npf-1)+1) ) / npf 
!	  
        do f = 3, size(f_m(m)%vec)
	  fid = f_m(m)%vec(f) - 1
!
          cnt = 0
          floop: do i = 1, size(faces(fid)%vec)
	    do k = 0, npf-1
	      if ( nodesE(k) .eq. faces(fid)%vec(i) ) then
		cnt = cnt + 1
	        nfac(cnt) = nodesE(k)
		if (cnt .eq. 2) exit floop
	      end if
	    end do
	  end do floop
            
	  n1 = nfac(1)
	  n2 = nfac(2)
!
	  do k = 0, npf-1
	    if ( (nodesE(k) .ne. n1 .and. nodesE(k) .ne. n2) .or. &
	         (nodesE(k) .ne. n2 .and. nodesE(k) .ne. n1) ) then
	      n3 = nodesE(k)
	      exit
	    end if
	  end do	         
!
	  rrR(:,1) = rr(:,n1+1) - (/ xcg, ycg /)
	  rrR(:,2) = rr(:,n2+1) - (/ xcg, ycg /)
	  rrR(:,3) = rr(:,n3+1) - (/ xcg, ycg /)
!
	  ed = rrR(:,2) - rrR(:,1)
	  alpha = -atan2(ed(2),ed(1))
!	
          rot_mat(1,1) = cos(alpha)
          rot_mat(2,1) = sin(alpha)
          rot_mat(1,2) = -rot_mat(2,1)
          rot_mat(2,2) =  rot_mat(1,1)
!
          do j = 1, 3
            rrR(:,j) = matmul(rot_mat,rrR(:,j))
	  end do
!
          if ( rrR(2,3) .lt. 0 ) then
	    faces(fid)%vec = (/ n1, n2, n2+Np_d, n1+Np_d /)
	  else
  	    faces(fid)%vec = (/ n2, n1, n1+Np_d, n2+Np_d /)	
	  end if
!
	end do	
!
        deallocate (locfac)
	deallocate (nodesE)
	deallocate (edgesE)
	deallocate (rrR)
!	
      end do
!
      do f = 0, Ntf-1
        write (11,*) size(faces(f)%vec), '(', faces(f)%vec, ')'
      end do
!
      write (11,*) ')'
      write (11,*) ''
      call write_foam_separator ( 11 )
      close (11)
!
!
      contains
!
!
      subroutine write_foam_header ( unit )
!     +-------------------------------------------------------------------------+
      implicit none
      integer, intent(in) :: unit
!     +-------------------------------------------------------------------------+
!
      write (unit,*) '/*--------------------------------*- C++ -*----------------------------------*\'
      write (unit,*) '| =========                 |                                                 |'
      write (unit,*) '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'
      write (unit,*) '|  \\    /   O peration     | Version:  2.3.x                                 |'
      write (unit,*) '|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |'
      write (unit,*) '|    \\/     M anipulation  |                                                 |'
      write (unit,*) '\*---------------------------------------------------------------------------*/'
      write (unit,*) 'FoamFile'
      write (unit,*) '{'
      write (unit,*) '    version     2.0;'
      write (unit,*) '    format      ascii;'
!
      end subroutine write_foam_header
!
!
!
      subroutine write_foam_separator ( unit )
!     +-------------------------------------------------------------------------+
      implicit none
      integer, intent(in) :: unit
!     +-------------------------------------------------------------------------+
!
      write (unit,*) '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
!
      end subroutine write_foam_separator
!
      end subroutine save_foam_mesh
!
!
      end module export_module
