!
      program  uhMesh
!
!     +---------------------------------------------------------------+
      use export_module
      use grid_1d
      use grid_2d
      use grid_3d
      use grid_types, only: meshtype
      use import_module
      use init_data
      use mesh_transformation
      use nodes
!
      implicit none
!
      integer, parameter :: ONED = 1
      integer, parameter :: TWOD = 2
      integer, parameter :: THRD = 3
      integer :: narg
!
      character (len=1)  :: action
      character (len=64) :: directory = './'
      character (len=64) :: grid_name
      character (len=64) :: conf_name
      character (len=64) :: sol_name
!
      real(kind=8) :: st, ft
!     +---------------------------------------------------------------+
!
      call cpu_time(st)
      write (*,*) ''
      write (*,*) ' --------------------------------------------------'
      write (*,*) '      Unstructured-Hybrid Mesh generator '
      write (*,*) ' --------------------------------------------------'
      write (*,*) '   Department of Aerospace Engineering' 
      write (*,*) '   Politecnico di Milano. 1995-2008'
      write (*,*) ''
      write (*,*) '   CFD Laboratory, Mechanical Engineering'
      write (*,*) '   McGill University, 2009-2015'
      write (*,*) ''      
      write (*,*) '   ACE, Mechanical & Aerospace Engineering'
      write (*,*) '   University of Strathclyde, since 2016'
      write (*,*) ''      
      write (*,*) '   Authors: Stefano Rebay, Alberto Guardone,'
      write (*,*) '            Daniele Dussin, Marco Fossati'
      write (*,*) ' --------------------------------------------------'
      write (*,*) ''
!
      narg = iargc()
      if ( narg .lt. 1 ) then
        print*, ''; print*, '   Error. No configuration file supplied.'
        print*, '   Usage: uhMesh.exe configuration_file'
        print*, ''; stop
      else
        call getarg ( 1, conf_name )
      end if
!
!     +----------------------------------------------------------------+ 
!     |  Open configuration file and read general definitions          | 
!     +----------------------------------------------------------------+ 
!
      open ( unit=idf_cfg, file=trim(conf_name) )
      call read_record ( 'UHM_DEFINITIONS', narg, action, grid_name )
!
      if ( action .ne. 'g' ) then
        print*, ''
        print*, '   Reading '//trim(grid_name)//'.'//input_fmt
        select case ( input_fmt )
          case ( 'su2' ); call read_su2_mesh ( grid_name )
          case ( 'npc' ); call read_npc_mesh ( grid_name )
          case ( 'fns' ); call read_fns_mesh ( grid_name )
          case default 
          write(*,*) '   Error. Unsupported input format'; stop
        end select
      end if
!
      open (unit=911, file=trim(grid_name)//'.msg')
      select case ( action )
!
!     +-------------------------+
!     |  Mesh conversion        |
!     +-------------------------+
!
        case ( 'c' )
!       Nothing to be done in this case. Everything is 
!       managed via input and output formats
!
        select case ( ouput_fmt )
	case ( 'npc' ); call save_npc_mesh ( grid_name )
	case ( 'su2' ); call save_su2_mesh ( grid_name, meshtype )
        case default 
        write(*,*) '   Error. Unsupported input format'; stop
	end select
!
!     +-----------------------+
!     |  Mesh scaling         |
!     +-----------------------+
!
        case ( 's' )
        call read_record ( 'S_MESH_PARAMETERS' )
        call scale_mesh
        grid_name = trim(grid_name)//'_scaled'
!
!     +----------------------------------------+
!     |  Mesh generation                       |
!     +----------------------------------------+
!
        case ( 'g' )
        select case ( meshtype )
!
!     +------------------------+
!     |  One dimensional mesh  |
!     +------------------------+
!
          case ( ONED )
!          call read_path_1d
          print*, ''
          print*, 'WARNING: 1d meshes temporarily diasbled. Need to resolve the'
          print*, 'new reading parameters conflict with extrude mesh'
          print*, ''
          stop
          call build_grid_1d ( grid_name, 0 )
!
!     +------------------------+
!     |  Two dimensional mesh  |
!     +------------------------+
!
          case ( TWOD )
          call init_boundary ( grid_name )
          call build_grid_2d ( 1, directory, grid_name )
          call save_grid_2d ( 1, directory, grid_name )
          call plot_grid_2d ( 1, directory, grid_name, -1, 1, 0 )
!
!     +-----------------------------------+
!     |  Three dimensional extruded mesh  |
!     +-----------------------------------+
!
          case ( THRD )
          call init_boundary ( grid_name )
          call build_grid_2d ( 1, directory, grid_name )
          call save_grid_2d ( 1, directory, grid_name )
          call extrude_path ( grid_name )
          call save_grid_3d ( 1, grid_name )
          call plot_grid_3d ( 1, grid_name ) 
!
!
          case DEFAULT 
          write (*,*); write (*,*) 'ERROR. Unknown type of mesh'
          write (*,*); stop
!
        end select
      end select
      close (911)
!
      select case ( trim(ouput_fmt) ) 
!
        case ( 'npc' )                                                 
        write (*,*)                                                    
        write (*,*) '   Saving NODEPAIR mesh...'
        if ( action .ne. 'g' ) call save_npc_mesh ( grid_name )                       
        write (*,*) '   - Nodes        ', Nnodes                       
        if ( Nquads .gt. 0 ) write (*,*) '   - Quadrilaters ', Nquads  
        if ( Ntrias .gt. 0 ) write (*,*) '   - Triangles    ', Ntrias  
        if ( Nhexas .gt. 0 ) write (*,*) '   - Hexahedra    ', Nhexas  
        if ( Nprism .gt. 0 ) write (*,*) '   - Prisms       ', Nprism  
!
        case ( 'fns' )                                                 
        write (*,*) ''  
        write (*,*) ' --------------------------------------------------'   
        write (*,*) '   Saving FENSAP mesh...'                         
        call save_fns_mesh( grid_name )                                
!
        case ( 'su2' )                                                 
        write (*,*) ''                                                 
        write (*,*) ' --------------------------------------------------'   
        write (*,*) '   Saving SU2 mesh...'                            
        call save_su2_mesh( grid_name, meshtype )                                
!
        case ( 'vtk' )                                                 
        write (*,*) ''                                                 
        write (*,*) ' --------------------------------------------------'   
        write (*,*) '   Saving VTK mesh...'                            
        call save_vtk_mesh( grid_name )                                
!
        case ( 'foam' )                                                 
        write (*,*) ''                                                 
        write (*,*) ' --------------------------------------------------'   
        write (*,*) '   Saving FOAM mesh...'                            
        call save_foam_mesh( grid_name )                                
!
        case default                                                   
        write (*,*); write (*,*) 'ERROR. Unknown export format'        
        write (*,*); stop                                              
!
      end select                                                       
      close (idf_cfg)                                                  
!
      call cpu_time(ft)
!
      write (*,*) ''
      write (*,'(3x,a9,1x,e9.2,1x,a1)') 'CPU time:', ft - st,  's'
      write (*,*) ' --------------------------------------------------'         
      write (*,*) ''  
!
      end program  uhMesh
