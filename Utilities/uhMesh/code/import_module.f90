!
      module import_module
!
!     +--------------------------------------------------------------------------------+
      use grid_2d, only: Ntrias, Nquads
      use grid_3d, only: Ntetra, Nhexas, Nprism, Npyram
!
      implicit none
      private
!      
      integer, dimension(14), parameter :: SU2_nodesE &
                                           = (/ 0,0,2,0,3,0,0,0,4,4,0,8,6,5 /)
      integer, dimension(14), parameter :: SU2_facesE &
                                           = (/ 0,0,2,0,3,0,0,0,4,4,0,6,5,5 /)
      integer, dimension(14), parameter :: SU2_edgesE &
                                           = (/ 0,0,1,0,3,0,0,0,4,6,0,12,9,8 /)
      integer, dimension(14), parameter :: SU2_NPM &
                                           = (/ 0,0,1,0,2,0,0,0,3,4,0,7,6,5 /)
                                           
      integer, parameter :: SU2_line          = 3,  &
                            SU2_triangle      = 5,  &
                            SU2_quadrilateral = 9,  &
                            SU2_tetrahedral   = 10, &
                            SU2_hexahedral    = 12, &
                            SU2_wedge         = 13, &
                            SU2_pyramid       = 14
!
      public :: read_npc_mesh, &
                read_su2_mesh, &
                read_fns_mesh, &
		export_su2_solution
!     +--------------------------------------------------------------------------------+
!
!
      contains
!
!
      subroutine read_npc_mesh ( grid_name )
!     +-----------------------------------------------+
      use nodes 
      use mesh_structure
      use node_pair_structure
      use np_topology_gen
!
      implicit none
      character (len=32), intent(in) :: grid_name
!     +-----------------------------------------------+
!
!     +-----------------+
!     |  Reading nodes  |
!     +-----------------+
!
      open ( unit=11, file='nodes.'//TRIM(grid_name) )  
      call read_nodes(11)                               
      close(11)       
!
!     +--------------------------------+
!     |  Reading element connectivity  |
!     +--------------------------------+
!
      open ( unit=11, file='grid.'//TRIM(grid_name) )
      call read_mesh(11)
      close (11)
!
!     +---------------------+
!     |  Edge connectivity  |
!     +---------------------+
!
      call node_pair_str_gen
!
      Nnodes = Np_d
      Ntrias = count( ele_type_d .eq. 2 )
      Nquads = count( ele_type_d .eq. 3 )
      Ntetra = count( ele_type_d .eq. 4 )
      Npyram = count( ele_type_d .eq. 5 )
      Nprism = count( ele_type_d .eq. 6 )
      Nhexas = count( ele_type_d .eq. 7 )
!
      return
      end subroutine read_npc_mesh 
!
!
!
!
!
      subroutine read_su2_mesh ( grid_name )
!     +---------------------------------------------------------+
      use nodes 
      use mesh_structure
      use node_pair_structure
      use np_topology_gen
!
      implicit none
      character (len=32), intent(in) :: grid_name
!
      character (len=32), dimension(:), allocatable :: bname
      character (len=1) :: kwrd
      integer :: i, j, et, dm, N_back, jcn
      integer :: k, Ne_l, N_bound, cnt
!     +---------------------------------------------------------+
!
      open (unit=11, file=trim(grid_name)//'.su2')
      read(11,*) kwrd, k_d
      read(11,*) kwrd, Ne_d
!
      allocate ( j_m_d(Ne_d) )
      allocate ( ele_type_d(Ne_d) )
      allocate ( ma_m_d(Ne_d) )
!
      do i = 1, Ne_d
        read(11,*) et          
        ele_type_d(i) = SU2_NPM(et)          
        allocate ( j_m_d(i)%vec(SU2_nodesE(et)) )         
        backspace(11)
        read(11,*) et, j_m_d(i)%vec, dm 
        j_m_d(i)%vec = j_m_d(i)%vec + 1
!
        if ( ele_type_d(i) .eq. 1 ) then
          allocate ( ma_m_d(i)%vec(2) )
        else if ( ele_type_d(i) .eq. 2 ) then
          allocate ( ma_m_d(i)%vec(3) )
        else if ( ele_type_d(i) .eq. 3 ) then
          allocate ( ma_m_d(i)%vec(4) )
        else if ( ele_type_d(i) .eq. 4 ) then
          allocate ( ma_m_d(i)%vec(4) )
        else if ( ele_type_d(i) .eq. 5 ) then
          allocate ( ma_m_d(i)%vec(5) )
        else if ( ele_type_d(i) .eq. 6 ) then
          allocate ( ma_m_d(i)%vec(5) )
        else if ( ele_type_d(i) .eq. 7 ) then
          allocate ( ma_m_d(i)%vec(6) )
        end if
        ma_m_d(i)%vec = 0
!
      end do
!
      read(11,*) kwrd, Np_d
      allocate ( rr(k_d,Np_d) )   
      do i = 1, Np_d
         read(11,*) rr(:,i)
      end do
!
      read(11,*) kwrd, N_bound
      allocate ( bname(N_bound) )
!
      Ne_b = 0
      Np_b = 0
!
      do i = 1, N_bound
        read(11,*) kwrd, bname(i)
        read(11,*) kwrd, Ne_l
        Ne_b = Ne_b + Ne_l
        do j = 1, Ne_l
          read(11,*) et
          Np_b = Np_b + SU2_nodesE(et)
        end do
      end do
!
      N_back = Ne_b + 2*N_bound
      do i = 1, N_back
        backspace(11)
      end do
!
      allocate ( j_m_b(Ne_b) )
      allocate ( ele_type_b(Ne_b) )
      allocate ( bound_m(Ne_b) )
      allocate ( bound_p(Np_b) )
      allocate ( jd_jb(Np_b) )
      allocate ( ma_m_b(Ne_b) )
!
      cnt = 1
      jcn = 1
      do i = 1, N_bound

        read(11,*) kwrd
        read(11,*) kwrd, Ne_l
        do j = 1, Ne_l
          read(11,*) et
          ele_type_b(cnt) = SU2_NPM(et)
          allocate ( j_m_b(cnt)%vec(SU2_nodesE(et)) )
          backspace(11)
          read(11,*) et, j_m_b(cnt)%vec
!
          j_m_b(cnt)%vec = j_m_b(cnt)%vec + 1
          bound_m(cnt) = i
!
          jd_jb(jcn:jcn-1+SU2_nodesE(et)) = j_m_b(cnt)%vec
          j_m_b(cnt)%vec = (/ (k, k=jcn,jcn-1+SU2_nodesE(et)) /)
          bound_p(jcn:jcn-1+SU2_nodesE(et)) = i
          jcn = jcn + SU2_nodesE(et)
!
          if ( ele_type_b(cnt) .eq. 1 ) then
            allocate ( ma_m_b(cnt)%vec(2) )
          else if ( ele_type_b(cnt) .eq. 2 ) then
            allocate ( ma_m_b(cnt)%vec(3) )
          else if ( ele_type_b(cnt) .eq. 3 ) then
            allocate ( ma_m_b(cnt)%vec(4) )
          end if
          ma_m_b(cnt)%vec = 0
!
          cnt = cnt + 1
        end do
!
      end do   
!
      close(11) 
!
      Nnodes = Np_d
      Ntrias = count( ele_type_d .eq. 2 )
      Nquads = count( ele_type_d .eq. 3 )
      Ntetra = count( ele_type_d .eq. 4 )
      Npyram = count( ele_type_d .eq. 5 )
      Nprism = count( ele_type_d .eq. 6 )
      Nhexas = count( ele_type_d .eq. 7 )
!
      return
      end subroutine read_su2_mesh
!
!
!
!
!
      subroutine read_fns_mesh ( grid_name )
!     +-------------------------------------------------------------------------+
      use init_data
      use mesh_structure
      use nodes 
      use node_pair_structure
      use np_topology_gen
!
      implicit none
!
      character(len=64), intent(in)  :: grid_name
!     +-------------------------------------------------------------------------+
!
      Nnodes = Np_d
      Ntrias = count( ele_type_d .eq. 2 )
      Nquads = count( ele_type_d .eq. 3 )
      Ntetra = count( ele_type_d .eq. 4 )
      Npyram = count( ele_type_d .eq. 5 )
      Nprism = count( ele_type_d .eq. 6 )
      Nhexas = count( ele_type_d .eq. 7 )
!
      return
      end subroutine read_fns_mesh 
!
!
!
!
!
      subroutine export_su2_solution ( sname )
      
      implicit none
      character (len=64), intent(in) :: sname
      character (len=1) :: aaa
      integer :: ef, i, Np, d
      real*8, dimension(:,:), allocatable :: sold
      real*8 :: dummy
      
      open ( unit=11, file=trim(sname) )
      Np = 0
      do
        read(11,*,iostat=ef)
	if ( ef .lt. 0 ) then
	  exit
	end if
        Np = Np + 1
      end do
      
      rewind(11)
      
      Np = Np-1
      allocate ( sold(4,Np) )
      
      read (11,*)
      do i = 1, Np
        read (11,*) d, dummy, dummy, sold(:,i)
      end do
       
      close (11)
       
      aaa = 'a'
      open ( unit=11, file='npc.solution' )
      do i = 1, 65+5
        write(11,*) aaa
      end do
      
      do i = 1, Np    
        write(11,*) i, sold(1,i)  
        write(11,*) sold(2:3,i)
        write(11,*) sold(4,i)
      enddo	
       
      close (11) 
       
      end subroutine export_su2_solution
!
      end module import_module
