!
      module init_data
!
!     +--------------------------------------------------------------------------------+
      use grid_types  !  b_grid types
      use list        !  list operations
      use curves      !  splines module
      use delaunay    !  delaunay basic module
      use smooth_bline
!
      implicit none
!
!     +----------------+
!     |  Type: VERTEX  |
!     +----------------+
!
      type vert
        type (point), pointer :: pnt
      end type vert
!
!     +--------------+
!     |  Type: EDGE  |
!     +--------------+
!
      type edge
!
        type (point), pointer :: pnt_head, axp_head, add_head
        type (point), pointer :: rem_head, gar_head
        type (edge),  pointer :: next, prev, ppp, old
        type (simplex), pointer :: spx_head, axs_head, int_head
        type (simplex), pointer :: ext_head, stc_head
!
        integer, dimension(:), pointer :: idata
        integer, dimension(:), pointer :: option
!
        integer :: index, ndata
        integer :: begin_vert, end_vert
        integer :: edge_kind, position, status
        integer :: struct, pnts_law, pnts_norm
        integer :: cl_open, adj_begin, adj_end
        integer :: leng_mod, leng_var
        integer :: cross_ovlp, auto_ovlp
        integer :: Npbms
!
        real*8, dimension(:,:), pointer :: hdata
        real*8, dimension(:), pointer :: height
        real*8, dimension(:), pointer :: sdata
!
        real*8 :: stretch, thick
        real*8 :: delta, leng_par
!
        character (len=3) :: curvetype
        character (len=32) :: edgename
!
!
      end type edge
!
!     +--------------+
!     |  Type: WAKE  |
!     +--------------+
!
      type wake
!
         type (point), pointer :: pnt_head, axp_head, add_head
         type (point), pointer :: rem_head, gar_head      
         type (point) :: origin
         type (simplex), pointer :: spx_head, axs_head, int_head
         type (simplex), pointer :: ext_head, stc_head
!
         real*8, dimension(:,:), pointer :: hdata
         real*8, dimension(:,:), pointer :: pdata
         real*8, dimension(:), pointer :: sdata
         real*8, dimension(:), pointer :: height         
!
         integer, dimension(:), pointer :: idata
         integer :: index, ndata, edge_kind
         integer :: begin_vert, end_vert
!
         character(len=1) :: wtype
!
      end type wake
!
!     +----------------+
!     |  Type: DOMAIN  |
!     +----------------+
!
      type domain
!
         type(vert), dimension(:), pointer :: verts
         type(edge), dimension(:), pointer :: edges
         type(wake), dimension(:), pointer :: wakes
!
         type(point),   pointer :: pnt_head, axp_head, add_head
         type(point),   pointer :: rem_head, gar_head
         type(simplex), pointer :: spx_head, axs_head, int_head
         type(simplex), pointer :: ext_head, stc_head
         type(point),   pointer :: points
!
         real*8 :: smooth
         real*8 :: gamma2
         real*8 :: gamma3
!
         integer :: number_of_verts
         integer :: number_of_edges
         integer :: number_of_layers
         integer :: number_of_wakes
         integer :: number_of_Uwakes
         integer :: number_of_Swakes
         integer :: done_layers
         integer :: insert_mtd
!
      end type domain
!
      type (domain), public :: grid
      type (domain), public :: back
      type (domain), public :: b_grid

      integer, parameter, public :: nd_b = 1
      integer, parameter, public :: nd_d = 2
      integer, parameter, public :: idf_cfg = 111
      integer, parameter, public :: idf_knt = 112
      integer, public :: number
!
      real (kind=8), dimension(:,:), allocatable :: backgrP
      real (kind=8), dimension(:,:), allocatable :: steinrP
      real*8, dimension(2) :: pbeg, pl
      real*8 :: stol, Re, uinf, rhoinf, visc, cf
      real*8 :: xfactor
      real*8 :: yfactor
      real*8 :: zfactor
      real*8 :: wlen
      real*8 :: wang
      real*8 :: h1, h2, wh

!
      integer :: m, w, npvl, Nstr_wakes, nwp
      integer :: max_lay
      integer :: steps
      integer :: lsmooth
      integer :: kmax
      integer :: nbackgr
      integer :: nsteinr
!
      character (len=32) :: floor_tag, roof_tag
      character (len=3) :: input_fmt
      character (len=5) :: ouput_fmt
      character (len=1) :: useypl      
!
      public :: load_topology, no_spaces, read_record
      private :: m, w, max_lay
!     +--------------------------------------------------------------------------------+
!
!
      contains
!
!
      subroutine read_record ( record_name, idx, act, gn )
!     +--------------------------------------------------------------------+
      use grid_1d
      implicit none
!
      character (*), intent(in) :: record_name
      integer, optional :: idx
      character (len=1), optional :: act
      character (len=64), optional :: gn
!
      character (len=26) :: keyword     
      character (len=1) :: recrd
      character (len=1) :: smooth
      character (len=20) :: command
      character (len=60) :: datafilename
!
      integer :: j, n, eof, eid
      integer :: Npassag, Gkernel, cnt, nbackg
!
      real*8, parameter :: pi=3.14159265358979d0
      real (kind=8), dimension(2,2) :: rot_mat
      real (kind=8), dimension(2)   :: xb, xe, x0
      real (kind=8) :: r1, r2, a1, a2
      real (kind=8) :: teta
!     +--------------------------------------------------------------------+
!
      n = 500
      do
        read (idf_cfg,'(a1)',iostat=eof) recrd
        if ( eof .lt. 0 ) then
          print*, ''
          print*, ' Error. End of file reached but'
          print*, trim(record_name)//' record was not found.'
          print*, ''
          stop
        else
!
          if ( recrd .eq. '#' ) then
            backspace (idf_cfg)
            read (idf_cfg,*) recrd, keyword
            if ( trim(keyword) .eq. trim(record_name) ) then
!
              select case ( trim(record_name) )
!              
                case ( 'UHM_DEFINITIONS' )
                read (idf_cfg,*) act
                read (idf_cfg,*) gn
                read (idf_cfg,*) meshtype
                read (idf_cfg,*) input_fmt
                read (idf_cfg,*) ouput_fmt
!
!
                case ( 'S_MESH_PARAMETERS' )
                read (idf_cfg,*) xfactor
                read (idf_cfg,*) yfactor
                if ( meshtype .eq. 3 ) read (idf_cfg,*) zfactor
!
!
                case ( 'G_MESH_PARAMETERS_BOUNDARY' )
                read (idf_cfg,*) b_grid % number_of_verts  
                read (idf_cfg,*) b_grid % number_of_edges  
                read (idf_cfg,*) b_grid % number_of_wakes
!
!     +----------------------------------------------------+
!     |  The following three parameters are:               |
!     |  - Height block. If high limits the height growth  |
!     |  - Condense angle                                  |
!     |  - Insert_mtd                                      |
!     +----------------------------------------------------+
!
                read (idf_cfg,*) b_grid % gamma2           
                read (idf_cfg,*) b_grid % gamma3           
                read (idf_cfg,*) b_grid % insert_mtd
!
!     +-------------------------------------------------------+
!     |  Automatic estimation for y+ at walls based on flat   |
!     |  plate correlations.                                  |
!     +-------------------------------------------------------+
!
		read (idf_cfg,*) useypl
		if ( useypl .eq. 'y' ) then
  		  read (idf_cfg,*) Re
		  read (idf_cfg,*) uinf
		  read (idf_cfg,*) rhoinf		
		  read (idf_cfg,*) visc
		  read (idf_cfg,*) npvl		
		end if
!
!
                case ( 'G_MESH_PARAMETERS_DOMAIN' )
                read (idf_cfg,*) steps
                read (idf_cfg,*) lsmooth
                read (idf_cfg,*) kmax
                read (idf_cfg,*) stol
!
!     +-------------------------------+
!     |  Backgrid and Steiner points  |
!     +-------------------------------+
!
                read (idf_cfg,*) nbackgr
                if ( nbackgr .gt. 0 .or. b_grid % number_of_Uwakes .gt. 0 ) then
                  if ( allocated(backgrP) ) deallocate (backgrP)
!
                  nbackg = nbackgr + sum( b_grid % wakes(:) % ndata )
!                 
                  allocate ( backgrP(3,nbackg) )
                  do j = 1, nbackgr
                    read (idf_cfg,*) backgrP(:,j)
                  end do
!
                  cnt = nbackgr
                  do j = 1, b_grid % number_of_wakes
                    if ( b_grid % wakes(j) % wtype .eq. 'u' ) then
                      backgrP(:,cnt+1:cnt+b_grid % wakes(j) % ndata) = &
                      b_grid % wakes(j) % pdata
                      cnt = cnt+ b_grid % wakes(j) % ndata
                    end if
                  end do
!
                  nbackgr = nbackg
!
                end if
                
                do j = 1, nbackgr
                write (112,*) backgrP(:,j)
                end do                
!
                read (idf_cfg,*) nsteinr
                if ( nsteinr .gt. 0 ) then
                  allocate ( steinrP(2,nsteinr) )
                  do j = 1, nsteinr
                    read (idf_cfg,*) steinrP(:,j)
                  end do
                end if
!
!
                case ( 'G_BOUNDARY_EDGE' )
                read (idf_cfg,*) eid                
                if ( eid .eq. idx ) then
!
                  b_grid % edges(eid) % index = eid
                  read (idf_cfg,'(a20)') command
                  call no_spaces( 20, command )
                  read (idf_cfg,*) b_grid % edges(eid) % edgename
                  read (idf_cfg,*) b_grid % edges(eid) % begin_vert, &
                                   b_grid % edges(eid) % end_vert
!
!     +------------------------------------------------------------+
!     |  Reading information on the geometry of the boundary edge  |
!     +------------------------------------------------------------+
!
                  select case ( command(1:3) )
!
!     +------------------------------------+
!     |  Drawing a boundary straight line  |
!     +------------------------------------+
!
                    case ('lin')
                    b_grid % edges(eid) % curvetype = 'lin'
                    read (idf_cfg,*) xb
                    read (idf_cfg,*) xe
		    !n = 1000000
                    n = 300
                    call gen_line( n, xb, xe )
		    n = 151
!
!     +-----------------------------------------------+
!     |  Drawing a boundary circle or arc or ellipse  |
!     +-----------------------------------------------+
!
                    case ('cir')
                    b_grid % edges(eid) % curvetype = 'cir'
                    read (idf_cfg,*) xb 
                    read (idf_cfg,*) r1, r2
                    read (idf_cfg,*) a1, a2
                    a1 = pi*a1/180.d0
                    a2 = pi*a2/180.d0
		    n = 300000
                    call gen_ell( n, xb, r1, r2, a1, a2 )
		    n = 500
!
!     +---------------------------------------+
!     |  Drawing a boundary from data points  |
!     +---------------------------------------+
!
                    case ('dat')
                    b_grid % edges(eid) % curvetype = 'dat'
                    read (idf_cfg,*) datafilename
                    call no_spaces( 60, datafilename )
                    read (idf_cfg,*) smooth, Npassag, Gkernel
                    read (idf_cfg,*) x0
                    read (idf_cfg,*) teta
!
                    teta = -teta*pi/180.d0
                    rot_mat(1,1) = cos(teta)
                    rot_mat(2,1) = sin(teta)
                    rot_mat(1,2) = -rot_mat(2,1)
                    rot_mat(2,2) =  rot_mat(1,1)
!
                    call gen_data_line( datafilename, x0, teta, rot_mat, &
                                        smooth, Npassag, Gkernel, eid )
!
!      +-----------------------------------------------------+
!      |  Drawing a boundary mesh spcified by user via file  |
!      +-----------------------------------------------------+
!
                    case ('bms')
                    b_grid % edges(eid) % curvetype = 'bms'
                    read (idf_cfg,*) datafilename
                    call no_spaces( 60, datafilename )
                    smooth = 'n'
                    Npassag = 0
                    Gkernel = 0
                    x0 = 0
                    teta = 0
                    rot_mat = 0
                    rot_mat(1,1) = 1
                    rot_mat(2,2) = 1
                    call gen_data_line( datafilename, x0, teta, rot_mat, &
                                        smooth, Npassag, Gkernel, eid )
!
                    case default
                    print*, ''; print*, 'Error. Unknown edge type'
                    print*, ''
!
                  end select
!
                else
                  cycle
                end if
!
!      +-----------------------------------------------+
!      |  Reading metrics associated to boundary edge  |
!      +-----------------------------------------------+
!
                read (idf_cfg,*) b_grid % edges(eid) % struct
                read (idf_cfg,*) b_grid % edges(eid) % delta
                read (idf_cfg,*) b_grid % edges(eid) % leng_mod
                read (idf_cfg,*) b_grid % edges(eid) % leng_var
                read (idf_cfg,*) b_grid % edges(eid) % leng_par
                read (idf_cfg,*) b_grid % edges(eid) % cross_ovlp
                read (idf_cfg,*) b_grid % edges(eid) % auto_ovlp
                read (idf_cfg,*) b_grid % edges(eid) % ndata
!
                n = b_grid % edges(eid) % ndata       
                b_grid % edges(eid) % edge_kind = 0
                b_grid % edges(eid) % position = 0
!
                if ( n .lt. 2 ) then
                  write (*,*) 'ERROR. Insufficent control points for edge '&
                    //trim(b_grid % edges(eid) % edgename)
                  write (*,*) ''; stop
                end if
!
                allocate ( b_grid % edges(eid) % idata(n) )
                allocate ( b_grid % edges(eid) % sdata(n) )
                allocate ( b_grid % edges(eid) % hdata(2*nd_d-1,n) )
                allocate ( b_grid % edges(eid) % height(n) )
                allocate ( b_grid % edges(eid) % option(n) )
                b_grid % edges(eid) % option = 0
!
                do j = 1, n
                  read(idf_cfg,*) b_grid % edges(eid) % sdata(j),   &  		  
                                  b_grid % edges(eid) % hdata(1,j), &
                                  b_grid % edges(eid) % idata(j),   &
                                  b_grid % edges(eid) % height(j),  &
                                  b_grid % edges(eid) % option(j)
!
                  b_grid % edges(eid) % hdata(2,j) = b_grid % edges(eid) % hdata(1,j)
                  b_grid % edges(eid) % hdata(3,j) = zero
!
                end do               
!
!     +-------------------------------------------------------+
!     |  Automatic estimation for y+ at walls based on flat   |
!     |  plate correlations.                                  |
!     +-------------------------------------------------------+
!
                if ( useypl .eq. 'y' ) then
		  cf = 0.026 / Re**(1.0/7.0)
		  b_grid % edges(eid) % height = (visc / ( sqrt(0.5*cf*uinf)*rhoinf )) / (npvl - 1)
		end if
!
!      +----------------------------------------------------+
!      |  Wakes are forced to be attached at the beginning  |
!      |  of the edge (not clear why...)                    |
!      +----------------------------------------------------+
!
                if ( b_grid % edges(eid) % option(n) .ne. 0 ) then
                   b_grid % edges(eid) % option(1) = b_grid % edges(eid) % option(n)
                   b_grid % edges(eid) % option(n) = 0
                end if
!
!
                case ( 'G_WAKE' )
                read (idf_cfg,*) eid                
                if ( eid .eq. idx ) then
!
                  read (idf_cfg,*) b_grid % wakes(eid) % wtype
                  b_grid % wakes(eid) % index = eid
!
                  if ( trim(b_grid % wakes(eid) % wtype) .eq. 's' ) then
!
                    read (idf_cfg,*) datafilename              
                    call no_spaces( 60, datafilename )      
                    read (idf_cfg,*) smooth, Npassag, Gkernel  
                    read (idf_cfg,*) teta                      
!
                    x0 = 0.d0
                    teta = -teta*pi/180.d0
                    rot_mat(1,1) = cos(teta)                
                    rot_mat(2,1) = sin(teta)                
                    rot_mat(1,2) = -rot_mat(2,1)            
                    rot_mat(2,2) =  rot_mat(1,1)
!
                    call gen_data_line( datafilename, x0, teta, rot_mat, &
                                        smooth, Npassag, Gkernel, eid )
!
                    read (idf_cfg,*) b_grid % wakes(eid) % ndata
                    b_grid % wakes(eid) % edge_kind = 0
                    n = b_grid % wakes(eid) % ndata       
!
                    if ( n .lt. 2 ) then
                      write (*,*) 'ERROR. Insufficent control points for wake ', w
                      write (*,*) ''; stop
                    end if
!
                    allocate ( b_grid % wakes(eid) % idata(n) )
                    allocate ( b_grid % wakes(eid) % sdata(n) )
                    allocate ( b_grid % wakes(eid) % hdata(2*nd_d-1,n) )
                    allocate ( b_grid % wakes(eid) % height(n) )
!
                    do j = 1, n
!
                       read (idf_cfg,*) b_grid % wakes(eid) % sdata(j),   &
                                        b_grid % wakes(eid) % hdata(1,j), &
                                        b_grid % wakes(eid) % idata(j),   &
                                        b_grid % wakes(eid) % height(j)
!
                       b_grid % wakes(eid) % hdata(2,j) = b_grid % wakes(eid) % hdata(1,j)
                       b_grid % wakes(eid) % hdata(3,j) = zero
!
                    end do
!
                  else
!
                    b_grid % number_of_Uwakes = b_grid % number_of_Uwakes + 1
                    read (idf_cfg,*) pbeg
                    read (idf_cfg,*) wlen
                    read (idf_cfg,*) wang
                    read (idf_cfg,*) h1, h2
!
                    nwp = int(wlen / ((h1+h2)/2))
                    b_grid % wakes(eid) % ndata = nwp
                    allocate ( b_grid % wakes(eid) % pdata(3,nwp) )
!                    
                    b_grid % wakes(eid) % pdata(:,1) = (/ pbeg, h1 /)
                    do j = 2, nwp
                      
                      wh = (h2 - h1)*j/nwp + h1
                      pl(1) = b_grid % wakes(eid) % pdata(1,j-1) + wh*cos(wang*3.14/180) 
                      pl(2) = b_grid % wakes(eid) % pdata(2,j-1) + wh*sin(wang*3.14/180)
!
                      b_grid % wakes(eid) % pdata(:,j) = (/ pl, wh /)
!                      
                    end do
!
                  end if
!
                else
                  cycle
                end if
!
                case ( 'G_EXTRUSION_PARAMETERS' )
                read (idf_cfg,*) floor_tag
                read (idf_cfg,*) roof_tag
                read (idf_cfg,*) Ns                                
!
	        allocate ( x_L(Ns), x_R(Ns) )
	        allocate ( metric_1(Ns), metric_2(Ns) )
	        allocate ( nm(Ns), str_type(Ns) )
! 
	        do j = 1, Ns
		  read (idf_cfg,*) x_L(j), x_R(j), nm(j), str_type(j), metric_1(j), metric_2(j)
	        enddo
!
                case DEFAULT
                print*, ''; print*, 'Error. Unknown record '
                print*, ''
!
              end select
!
              rewind (idf_cfg)
              exit
            else
              cycle
            end if
          else
            cycle
          end if
!
        end if
      end do
!
      return
      end subroutine read_record
!
!
!
      subroutine gen_data_line ( datafilename, x0, teta, rot_mat, smooth, &
                                 Npassag, Gkernel, eid )
!     +--------------------------------------------------------------------+
      implicit none
!
      character (len=60), intent(in) :: datafilename
      real*8, dimension(:), intent(in) :: x0
      real*8, intent(in) :: teta
      real*8, dimension(2,2), intent(in) :: rot_mat
      character (len=1), intent(in) :: smooth
      integer, intent(in) :: Npassag
      integer, intent(in) :: Gkernel
      integer, intent(in) :: eid
!
      real (kind=8), dimension(:,:), allocatable :: xdata
      integer :: j, n, idx, idf
      character (len=1) :: sort
!     +--------------------------------------------------------------------+
!
      idf = 113
      open ( idf, file=datafilename, form='formatted' )
!
      read (idf,*) n, sort
      call header( n, 0 )
      b_grid % edges(eid) % Npbms = n
!
      allocate ( xdata(2,n) )
      xdata = 0.d0
!
      do j = 1, n
        if ( sort .eq. 'y' ) then
          idx = n - (j-1)
        else
          idx = j
        end if
        read (idf,*) xdata(:,idx)  
      enddo
!
      close (idf)
!
      if ( smooth .eq. 'y' ) then
        call convolution ( Gkernel, Npassag, xdata )
      end if
!
      do j = 1, n
        xdata(:,j) = matmul(rot_mat,xdata(:,j)) + x0
        write (idf_knt,'(2(2x,e22.16))') xdata(:,j)
      end do
!
      deallocate ( xdata )      
!      
      return
      end subroutine gen_data_line
!
!
!
      subroutine init_boundary ( name )
!
!     +--------------------------------------------------------------------+
!     |  routine to generate geometry.anyname file starting from scratch.  |
!     |  generates a global file containing geometric data, as needed for  |
!     |  subroutine build_grid_2d                                          |
!     +--------------------------------------------------------------------+
      implicit none
!
      character(len=64), intent(in) :: name
!
      logical, parameter :: ascii = .true.
      character (len=60) :: outfilename
      integer :: i, j, k, idx
!     +---------------------------------------------------------------------+
!
      call read_record ( 'G_MESH_PARAMETERS_BOUNDARY' )
!
      if ( .not. associated(b_grid % verts) )  &
        allocate( b_grid % verts(b_grid % number_of_verts) )
      if ( .not. associated(b_grid % edges) )  &
        allocate( b_grid % edges(b_grid % number_of_edges) )
      if ( .not. associated(b_grid % wakes) )  &
        allocate( b_grid % wakes(b_grid % number_of_wakes) )
!
      m = b_grid % number_of_edges
      w = b_grid % number_of_wakes
!
      b_grid % number_of_Uwakes = 0
      b_grid % number_of_Swakes = 0
!
!     +-------------------------------------------------------------+
!     |  File containing geomtrical definitions for boundary edges  |
!     +-------------------------------------------------------------+
!
      k = len_trim(name)
      outfilename = 'knots.'//name(1:k)
      open ( idf_knt, file = outfilename, form = 'formatted' )
      write (idf_knt,'(6x,a6)') 'N_CURV'
      write (idf_knt,'(i12)') m
!
!     +------------------------------------------------------------+
!     |  Read EDGES metric and topology and create boundary lines  |
!     +------------------------------------------------------------+
!
      do i = 1, m
!
        call read_record ( 'G_BOUNDARY_EDGE', i )
!
        b_grid % edges(i) % pnt_head  =>  newhead_point ()
        b_grid % edges(i) % axp_head  =>  newhead_point ()
        b_grid % edges(i) % add_head  =>  newhead_point ()
        b_grid % edges(i) % rem_head  =>  newhead_point ()
        b_grid % edges(i) % gar_head  =>  newhead_point ()
!
        b_grid % edges(i) % spx_head  =>  newhead_simplex ()
        b_grid % edges(i) % axs_head  =>  newhead_simplex ()
        b_grid % edges(i) % int_head  =>  newhead_simplex ()
        b_grid % edges(i) % ext_head  =>  newhead_simplex ()
        b_grid % edges(i) % stc_head  =>  newhead_simplex ()
!
      end do   
      close (idf_knt)
!
      if ( w .ne. 0 ) then
!
!     +-----------------------------------------------------+
!     |  File containing geometrical definitions for wakes  |
!     +-----------------------------------------------------+
!
        outfilename = 'wknots.'//name(1:k)
        open ( unit=idf_knt, file=outfilename, form = 'formatted' )
        write (idf_knt,'(6x,a6)') 'N_CURV'
        write (idf_knt,'(i12)') w
!
!     +---------------------------------------------------------+
!     |  Read WAKES metric and topology and create wakes lines  |
!     +---------------------------------------------------------+
!
        do i = 1, w
!
          call read_record ( 'G_WAKE', i )
!
          if ( b_grid % wakes(i) % wtype .eq. 's' ) then
!
            b_grid % number_of_Swakes = b_grid % number_of_Swakes + 1
!
            b_grid % wakes(i) % pnt_head  =>  newhead_point ()
            b_grid % wakes(i) % axp_head  =>  newhead_point ()
            b_grid % wakes(i) % add_head  =>  newhead_point ()
            b_grid % wakes(i) % rem_head  =>  newhead_point ()
            b_grid % wakes(i) % gar_head  =>  newhead_point ()
!
            b_grid % wakes(i) % spx_head  =>  newhead_simplex ()
            b_grid % wakes(i) % axs_head  =>  newhead_simplex ()
            b_grid % wakes(i) % int_head  =>  newhead_simplex ()
            b_grid % wakes(i) % ext_head  =>  newhead_simplex ()
            b_grid % wakes(i) % stc_head  =>  newhead_simplex ()
          end if
!
        end do
        close (idf_knt)
!
      end if
!
      call make_geometry ( name )
!
      b_grid % pnt_head  =>  newhead_point ()
      b_grid % axp_head  =>  newhead_point ()
      b_grid % add_head  =>  newhead_point ()
      b_grid % rem_head  =>  newhead_point ()
      b_grid % gar_head  =>  newhead_point ()
!
      b_grid % spx_head  =>  newhead_simplex ()
      b_grid % axs_head  =>  newhead_simplex ()
      b_grid % int_head  =>  newhead_simplex ()
      b_grid % ext_head  =>  newhead_simplex ()
      b_grid % stc_head  =>  newhead_simplex ()
!
!     +-------------------------------------------------------+
!     | Computing the maximum number of quadrilateral layers  |
!     +-------------------------------------------------------+
!
      max_lay = 0   
      do i = 1, m
        max_lay = max( b_grid % edges(i) % struct, max_lay )   
      end do
      b_grid % number_of_layers = max_lay
!
      return
      end subroutine init_boundary
!
! 
!
!
!
      subroutine  connect_edges ( grid )
!
!     +-----------------------------------------------------+
      implicit none
!
      type(domain), intent(inout) :: grid
      integer :: i, j
      integer :: bv, bv1 
      integer :: ev, ev1
!     +-----------------------------------------------------+
!
!     +---------------------------------------+
!     |  Initialize the edge pointer to zero  |
!     +---------------------------------------+
!
      do i = 1, grid % number_of_edges
        nullify ( grid % edges(i) % next )
        nullify ( grid % edges(i) % prev )
        grid % edges(i) % status = 0
      enddo
!
      do i = 1, grid % number_of_edges
!
         bv = grid % edges(i) % begin_vert
         ev = grid % edges(i) % end_vert   
         grid % edges(i) % status = 1
!
         do j = 1, grid % number_of_edges
!
            if ( grid % edges(j) % status .eq. 1 )  cycle
!
            bv1 = grid % edges(j) % begin_vert
            ev1 = grid % edges(j) % end_vert
!
            if ( bv1 .eq. ev ) then
              grid % edges(j) % status = 1
              grid % edges(i) % next => grid % edges(j)
              grid % edges(j) % prev => grid % edges(i)    
            end if
!
!      +----------------------------------------------------------+
!      |  Edge i and egde j form a circle, i.e. the begin vertex  |
!      |  of one is also the end vertex of the other              |
!      +----------------------------------------------------------+
!
            if ( ev1 .eq. bv ) then
              grid % edges(i) %  prev => grid % edges(j)
            end if
!
         enddo
      enddo
!
!     +---------------------------------------+
!     |  The status of the edges is reset to  | 
!     |  zero (maybe for later use...)        |
!     +---------------------------------------+
!
      do i = 1, grid % number_of_edges
        grid % edges(i) % status = 0
      enddo
!
      return
      end subroutine connect_edges





  SUBROUTINE  no_spaces ( len_str, string )
  !----------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: len_str
  CHARACTER (LEN=len_str),INTENT(INOUT) :: string  
  INTEGER ::j,k
  !----------------------------------------------------------------------------------

  k=0
  DO j=1,len_str       
     IF ( string(j:j).EQ.' ')THEN
        k=k+1
     ELSE
        EXIT
     ENDIF
  ENDDO
  DO j=1,len_str-k
     string(j:j)=string(j+k:j+k)
  ENDDO
  DO j=len_str-k,len_str
     string(j:j)=' '
  ENDDO

  END SUBROUTINE  no_spaces


  
  

  SUBROUTINE gen_ell ( n, c, r1, r2, beg, end )
  !----------------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), DIMENSION(2), INTENT(IN) :: c
    REAL(KIND=8), INTENT(IN) :: r1, r2, beg, end

    INTEGER :: i
    REAL(KIND=8) :: delta
  !----------------------------------------------------------------------------------

    delta = (end-beg)/REAL(n-1,8)

    CALL header ( n, 2)
    WRITE(idf_knt,'(2f25.15)') (c(1) + r1*cos(beg+delta*i), c(2) +   &
                                r2*sin(beg+delta*i), i=0, n-1)

  END SUBROUTINE gen_ell
  




  SUBROUTINE gen_line ( n, beg, end )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL (KIND=8), DIMENSION(2), INTENT(IN) :: beg, end

    INTEGER :: i
    REAL (KIND=8), DIMENSION(2) :: delta

    delta = (end-beg)/REAL(n-1,8)

    CALL header ( n, 0)
    WRITE (idf_knt,'(2f25.15)') (beg+delta*i,i=0,n-1)

  END SUBROUTINE gen_line





  SUBROUTINE header ( n, m )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n, m

    WRITE (idf_knt,'(2a12)') '         DIM','      POINTS'
    WRITE (idf_knt,'(2i12)') 2,n
    WRITE (idf_knt,'( a12)') '         IDC'
    WRITE (idf_knt,'( i12)') m
    WRITE (idf_knt,'(2a12)') '           X','           Y'

  END SUBROUTINE header

  
  
  
  
      subroutine make_geometry( filename )
!
!     +--------------------------------------------------+
!     |  Reads a geometrical data point file, computes   |
!     |  the cubic spline passing through the nodes and  |
!     |  save the results (hermite form).                |
!     +--------------------------------------------------+
!
      implicit none
!
      character (len=64), intent(in) :: filename
!
      logical, parameter :: ascii = .true.
      integer :: n
      integer :: idf = 113
!     +--------------------------------------------------+
!
      n = len_trim(filename)
!
!     +------------------+
!     |  Boundary lines  |
!     +------------------+
!
      open ( unit=idf_knt, file='knots.'//filename(1:n), form='formatted' )
      call new_curves ( idf_knt, ascii, curv )
      close (idf_knt)
!
      open ( unit=idf, file='geometry.'//filename(1:n), form='formatted' )
      call save_curves ( idf, ascii, curv )
      close (idf)
!
!     +---------+
!     |  Wakes  |
!     +---------+
!
      if ( Nstr_wakes .ne. 0 ) then
!
        open ( unit=idf_knt, file='wknots.'//filename(1:n), form='formatted' )
        call new_curves ( idf_knt, ascii, wcurv )
        close (idf_knt)
!
        open ( unit=idf, file='wgeometry.'//filename(1:n), form='formatted' )
        call save_curves ( idf, ascii, wcurv )
        close (idf)
!
      end if
!
      return
      end subroutine make_geometry





  SUBROUTINE load_topology(idf) 

    !  READ THE 2D TOPOLOGY AND INTIALIZE THE VARIOUS LISTS

    IMPLICIT NONE

    REAL(KIND=8) , PARAMETER     :: pi = 3.1415927464102
    INTEGER      , INTENT(IN)    :: idf

    INTEGER :: i, j, n

    IF ( .NOT. ASSOCIATED ( b_grid % verts ) ) THEN
       ALLOCATE ( b_grid % verts(b_grid % number_of_verts) )
    ENDIF

    IF ( .NOT. ASSOCIATED ( b_grid % edges ) ) THEN
       ALLOCATE ( b_grid % edges(b_grid % number_of_edges) )
    ENDIF

    WRITE(*,*)' loading topology...'
    READ(idf,*) ! TOPOLOGY

    DO i=1,b_grid % number_of_edges

       READ (idf,*) ! #
       READ (idf,*) ! BEGIN
       READ (idf,*) ! INDEX, NDATA, BEGIN_VERT, END_VERT
       READ (idf,*) b_grid % edges(i) % index,         &
            b_grid % edges(i) % ndata,         &
            b_grid % edges(i) % begin_vert,    &
            b_grid % edges(i) % end_vert

       n = b_grid % edges(i) % ndata

       !by DD begin   
       b_grid % edges(i) % edge_kind = 0
       b_grid % edges(i) % position = 0
       !by DD end	 
       IF ( n < 2 ) THEN
          WRITE (*,*) 'load_topology error, valid ndata >= 2; I have ',n
          WRITE (*,*) '(stop)'
          STOP
       ENDIF

       ALLOCATE ( b_grid % edges(i) % idata(n) )
       ALLOCATE ( b_grid % edges(i) % sdata(n) )
       ALLOCATE ( b_grid % edges(i) % hdata(2*nd_d-1,n) )

       READ (idf,*) ! IDATA, SDATA, HDATA
       DO j=1,n
          READ (idf,*) b_grid % edges(i) % idata(j),   &
               b_grid % edges(i) % sdata(j),   &
                                !by DD begin
               b_grid % edges(i) % hdata(1,j)
          b_grid%edges(i)%hdata(2,j) = b_grid%edges(i)%hdata(1,j)
          b_grid%edges(i)%hdata(3,j) = zero
          !by DD end
       ENDDO
       READ (idf,*) ! END

       b_grid % edges(i) % pnt_head  =>  newhead_point ()
       b_grid % edges(i) % axp_head  =>  newhead_point ()
       b_grid % edges(i) % add_head  =>  newhead_point ()
       b_grid % edges(i) % rem_head  =>  newhead_point ()
       b_grid % edges(i) % gar_head  =>  newhead_point ()

       b_grid % edges(i) % spx_head  =>  newhead_simplex ()
       b_grid % edges(i) % axs_head  =>  newhead_simplex ()
       b_grid % edges(i) % int_head  =>  newhead_simplex ()
       b_grid % edges(i) % ext_head  =>  newhead_simplex ()
       b_grid % edges(i) % stc_head  =>  newhead_simplex ()

    ENDDO

    b_grid % pnt_head  =>  newhead_point ()
    b_grid % axp_head  =>  newhead_point ()
    b_grid % add_head  =>  newhead_point ()
    b_grid % rem_head  =>  newhead_point ()
    b_grid % gar_head  =>  newhead_point ()

    b_grid % spx_head  =>  newhead_simplex ()
    b_grid % axs_head  =>  newhead_simplex ()
    b_grid % int_head  =>  newhead_simplex ()
    b_grid % ext_head  =>  newhead_simplex ()
    b_grid % stc_head  =>  newhead_simplex ()


  END SUBROUTINE load_topology


END module init_data
