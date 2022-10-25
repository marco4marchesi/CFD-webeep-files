!
      module grid_2d                                                                            
!
!     +---------------------------------------------------------------+  
      use grid_types     !  grid types                                                          
      use list           !  list operations                                                     
      use delaunay       !  delaunay basic module                                               
      use eucl_delaunay  !  delaunay euclid module                                              
      use metric_2d      !  metric module                                                       
      use back_2d        !  backgrid module                                                     
      use boundary_2d    !  boundary grid module                                                
      use front_2d       !  front grid module                                                   
      use arrays         !  arrays data structure                                               
      use refine_2d      !  refine module                                                       
      use check_grid     !  check module & quality properties                                   
      use init_data      !  utility for initial data generation                      
      use qgraph
!
      implicit none
!
      integer, public :: Ntrias, Nquads
!
      public :: build_grid_2d, &
                save_grid_2d,  &
                plot_grid_2d                                       
!     +---------------------------------------------------------------+  
!
!
      contains                                                                                  
!
!
      subroutine build_grid_2d ( idf, dir, fil )
!     +---------------------------------------------------------------+  
      implicit none                                                                             
!
      integer,            intent(in) :: idf                                            
      character (len=64), intent(in) :: dir, fil                                                
!
      integer :: choice
      integer :: d_end, f_end
      integer :: struct                              
!
      logical, parameter :: ascii = .true.                                                      
!     +---------------------------------------------------------------+
!
      write (*,*)'   Generating boundary grid...'
!
      d_end = last_c_leng (64, dir)                                                           
      f_end = last_c_leng (64, fil)                                                           
!
      choice = 0
      call connect_edges( b_grid )
!
      open ( unit=idf, file='geometry.'//fil(1:f_end) )
      call load_curves( idf, ascii, curv )
      close (idf)
!
      call delaunay_setup
      struct = 1                                                                              
      call boundary_grid ( idf, b_grid, dir, fil )
      call load_back                                                                          
!
      call read_record ( 'G_MESH_PARAMETERS_DOMAIN' )
      call build_back ( grid, back )                                   
!
      if ( steps .gt. 0 ) then                                                                   
!
        call front_grid ( grid, back, steps, choice )
        call build_arrays ( grid )
!
        if ( lsmooth .eq. 1 ) then
          write (*,*)'   Laplacian smoothing...'
          call smooth ( grid, kmax, stol )
          call update_geometry ( nd_d, grid % int_head )
        else
          write (*,*)'   no Laplacian smoothing'
        end if
!
      else
!
        write (*,*)'   Building triangulation arrays (backgrid only)...'
        call build_arrays ( grid )

      end if 

      end subroutine build_grid_2d                                                              
                                                         
  
  
  
      subroutine load_back
        
        implicit none                                                                           

        integer :: i                                                                            
        back%number_of_edges = grid %number_of_edges                                            
        back%number_of_verts = grid %number_of_verts                                            
        allocate(back%edges(back%number_of_edges))                                              
        allocate(back%verts(back%number_of_verts))                                              

        do i = 1, grid % number_of_edges
                                                                    
           back%edges(i)%index = grid%edges(i)%index                                            
           back%edges(i)%begin_vert = grid%edges(i)%begin_vert                                  
           back%edges(i)%end_vert = grid%edges(i)%end_vert                                      
           back%edges(i)%edge_kind = grid%edges(i)%edge_kind                                    
           back%edges(i)%position = grid%edges(i)%position                                      

           back % edges(i) % pnt_head  =>  newhead_point ()                                     
           back % edges(i) % axp_head  =>  newhead_point ()                                     
           back % edges(i) % add_head  =>  newhead_point ()                                     
           back % edges(i) % rem_head  =>  newhead_point ()                                     
           back % edges(i) % gar_head  =>  newhead_point ()                                     

           back % edges(i) % spx_head  =>  newhead_simplex ()                                   
           back % edges(i) % axs_head  =>  newhead_simplex ()                                   
           back % edges(i) % int_head  =>  newhead_simplex ()                                   
           back % edges(i) % ext_head  =>  newhead_simplex ()                                   
           back % edges(i) % stc_head  =>  newhead_simplex ()                                   

        enddo                                                                                   


        back % pnt_head  =>  newhead_point ()                                                   
        back % axp_head  =>  newhead_point ()                                                   
        back % add_head  =>  newhead_point ()                                                   
        back % rem_head  =>  newhead_point ()                                                   
        back % gar_head  =>  newhead_point ()                                                   

        back % spx_head  =>  newhead_simplex ()                                                 
        back % axs_head  =>  newhead_simplex ()                                                 
        back % int_head  =>  newhead_simplex ()                                                 
        back % ext_head  =>  newhead_simplex ()                                                 
        back % stc_head  =>  newhead_simplex ()                                                 

      end subroutine load_back                                                                  

  
      subroutine save_grid_2d ( idf, dir, fil )

        use nodes, only: Nnodes
        implicit none                                                                           

        integer, intent(in) :: idf
        character (len=64), intent(in) :: dir, fil                                              

        logical :: ascii = .false.                                                              
        integer :: d_end, f_end, p, m, n, j, offset                                             


        d_end = last_c_leng (64, dir)                                                           
        f_end = last_c_leng (64, fil)                                                           
                                  
        write(*,*) '   Building Delaunay arrays...'
        call build_arrays ( grid )

        if ( numbstrucedge .gt. 0 ) then                                                             
!
          write(*,*) '   Building quadrilateral arrays...'
          call structured_arrays
          write(*,*) '   Coupling arrays...'
          call couple_arrays
          call hyb_save_conv( fil, f_end, Nnodes, Ntrias, Nquads )
!
        else
!
          call uns_save_conv( fil, f_end, Nnodes, Ntrias )
!
        endif                                                                                   

        close (idf)                                                                            

        ascii=.true.
        open  (11,form='formatted',file='bnknots.'//fil(1:f_end)) 
        call new_curves( 11, ascii, curv )
        close (11)                                                                              

        open  (11,form='formatted',file='abscissa.'//fil(1:f_end))                              
          write(11,*) fil(1:f_end), ' - grid name'                                              
          write(11,*) b_grid % number_of_edges, ' - number of boundaries'                       
          write(11,*) size(b_bound_p), ' - number of points'                                    
          write(11,*)                                                                           
          write(11,*) '      node index          curve coordinate          boundary index'      
          write(11,*) ' (boundary notation)'                                                    

          offset = 0                                                                            

          do m = 1, b_grid % number_of_edges                                                    
                                                                                                
            p = size(curv(m)%s)                                                                 
            j = 2 + offset                                                                      
                                                                                                
            write(11,1112) 1 + offset, curv(m)%s(1)/curv(m)%s(p), m                             
                                                                                                
            do n = 2, p-1                                                                       
              write(11,1112) j,   curv(m)%s(n)/curv(m)%s(p), m                                  
              write(11,1112) j+1, curv(m)%s(n)/curv(m)%s(p), m                                  
              j = j+2                                                                           
            enddo                                                                               
                                                                                                
            write(11,1112) 2*(p-2)+2+offset, curv(m)%s(p)/curv(m)%s(p), m                       

            offset = offset + 2*(p-2)+2                                                         

            1112 format(4x,i7,15x,e15.7,12x,i7)                                                 

          enddo                                                                                 
                                                                                                
        close (11)                                                                              


      end subroutine save_grid_2d                                                               


  
      subroutine plot_grid_2d ( idf, dir, fil, metric_switch, &
        back_switch, ext_switch )

        implicit none                                                                           

        !  metric_switch:                                                                       
        !  < 0 no metric  plot                                                                  
        !  = 0    metric  plot                                                                  
        !  = 1    ellipse plot                                                                  
        !  = 2    length  plot                                                                  
        !  back_switch:                                                                         
        !  = 0 no back_grid  plot                                                               
        !  = 1    back_grid  plot                                                               
        !  ext_switch:                                                                          
        !  = 0 no external_grid  plot                                                           
        !  = 1    external_grid  plot                                                           

        integer, intent(in) :: idf, metric_switch, back_switch, ext_switch           
        character (len=64), intent(in) :: dir, fil                                              

        logical :: do_len = .false.                                                             
        integer :: d_end, f_end                                                                 
        character (len=64):: file_name                                                          

        d_end = last_c_leng (64, dir)                                                           
        f_end = last_c_leng (64, fil)                                                           



        if ( back_switch .eq. 1 ) then                                                               

           open  (idf,file='backgrid.'//fil(1:f_end)//'.mtv')                                   
           call g_plot (idf, back, back%int_head, back%ext_head)                                
           close (idf)                                                                          

           if (metric_switch >= 0)then                                                          
  
              if (metric_switch == 0) then                                                       
                 file_name = 'back_metric.'//fil(1:f_end)                                       
              else if(metric_switch == 1) then                                                  
                 file_name = 'back_ellipse.'//fil(1:f_end)                                      
              end if                                                                             
  
              open  (idf,file=dir(1:d_end)//file_name, status = 'unknown')                      
                call m_plot (idf,back%int_head, metric_switch)                                  
              close (idf)                                                                       
  
           end if                                                                                

        end if                                                                                   
!
!      +----------------------------------------------------+
!      |  The Delaunay part of the grid is now appended to  |
!      |  the quadrilateral part
!      +----------------------------------------------------+
!
        if ( numbstrucedge .gt. 0 ) then
!
          open ( unit=idf, file='mesh.'//fil(1:f_end)//'.mtv', &
            status='old', position='append' )
          call g_plot ( idf, grid, grid%int_head, grid%ext_head )                               
          close (idf)                                                                           
!
        else                                                                                    
!
!      +------------------------------------------------------------+
!      |  In this case no quadrilaters so the mesh.mtv file is new  |
!      +------------------------------------------------------------+
!
          open ( unit=idf, file='mesh.'//fil(1:f_end)//'.mtv', &
            status='unknown')
          write (idf,*) '$ data = curve2d'
          write (idf,*) '% equalscale = true'
          write (idf,*) '% boundary = true'
          write (idf,*) '% toplabel='' '' xlabel='' '' ylabel='' ''   '
          write (idf,*) '#%pointid = true'
          call g_plot ( idf, grid, grid%int_head, grid%ext_head ) 
          close (idf)
!
        end if                                                                                   

!         !by dd:begin-------------------------------------------                                 
!         if (metric_switch >= 0)then                                                             
!            if(metric_switch == 2) then                                                          
!               do_len = .true.                                                                   
!               file_name = 'len_grid.'//fil(1:f_end)                                             
!               !write (*,*) 'plotting  length grid ...', file_name                               
!               open  (idf,file=dir(1:d_end)//file_name, status = 'unknown')                      
!                 call g_plot (idf, grid, grid%int_head, grid%ext_head, do_len)                   
!               close (idf)                                                                       
!            else                                                                                 
!               if(metric_switch == 0) then                                                       
!                  file_name = 'metric_grid.'//fil(1:f_end)                                       
!                  !write (*,*) 'plotting  metric grid ...', file_name                            
!               else if(metric_switch == 1) then                                                  
!                  file_name = 'ellipse_grid.'//fil(1:f_end)                                      
!                  !write (*,*) 'plotting  ellipse grid ...', file_name                           
!               endif                                                                             
!               open  (idf,file=dir(1:d_end)//file_name, status = 'unknown')                      
!               call m_plot (idf, grid%int_head, metric_switch)                                   
!               close (idf)                                                                       
!            endif                                                                                
!         endif                                                                                   
!         !by dd:end (filename not defined unless metric_switch >=0!)                             
! 
!         if(ext_switch == 1) then                                                                
! !           write (*,*) 'plotting external grid...'                                             
!            open  (idf,file=dir(1:d_end)//'/extplot.'//fil(1:f_end))                             
!              call g_plot (idf, grid, grid%ext_head, grid%int_head)                              
!            close (idf)                                                                          
!         endif                                                                                   
! 
      return
      end subroutine plot_grid_2d



      subroutine g_plot (idf, grid, r_int_head, r_ext_head, do_len)                             

        ! plot the grid                                                                         
        ! dom_color -> accepted faces                                                           
        ! bad_color -> not accepted faces (l > tol_l)                                           
        ! not_color -> not accepted faces (l > max_l)                                           

        implicit none                                                                           

        logical, optional, intent(in) :: do_len                                                 
        integer          , intent(in) :: idf                                                    
        type (domain)    , intent(in) :: grid                                                   
        type (simplex)   , pointer    :: r_int_head, r_ext_head                                 

        !by dd:excluding riem_delaunay...                                                       
        real(kind=8), parameter       :: max_l  = 1.50d0    ! maximum distance among            
        ! accepted points 1.4                                                                   

        real(kind=8), parameter       :: tol_l  = 2.20d0    ! tolerance among                   
        ! accepted points 2.2                                                                   

        !by dd end                                                                              
        integer, parameter :: background = -1, foreground = 0,   &                              
             yellow = 1, blue = 2, green = 3, red = 4,   &                                      
             dark_blue = 5, orange = 6, pink = 7,   &                                           
             light_pink = 8, cyan = 9, brown = 9   !  plotmtv colors                            

        integer, parameter :: none = 0, dot = 1, cross = 2, capital_x = 3,   &                  
             white_square = 4, black_square = 5,   &                                            
             white_diamond = 6, black_diamond = 7,   &                                          
             white_triangle = 8, black_triangle = 9,   &                                        
             white_inv_triangle = 10, black_inv_triangle = 11,   &                              
             white_circle = 12, black_circle = 13   !   plotmtv markers                         

        integer, parameter :: dom_color = 0                                                     
        integer, parameter :: bad_color = red, not_color = orange                               
        integer, parameter :: bou_color = yellow                                                
        integer, parameter :: marker = none                                                     

        logical :: first                                                                        
        integer :: i, j, k, ipl, e, n_lines                                                     
        real (kind=8) :: length                                                                 
        type (point), pointer :: p1, p2                                                         
        type (simplex), pointer :: s_head, s                                                    
        real (kind=8), dimension(2) :: xpl                                                      


        !write (idf,*) '$ data = curve2d'                                                       
        !write (idf,*) '% equalscale = true'                                                    
        !write (idf,*) '% boundary = true'                                                      

        !**    goto 10 ! to print only the boundaries                                           

        !call q_g_plot(idf)                                                                     

        s_head  =>  r_ext_head; s  =>  s_head                                                   
        do; s  =>  s % next; if ( associated (s, s_head) ) exit                                 
           s % status = - s % status                                                            
        enddo                                                                                   

        if(.not. present(do_len)) then                                                          

           s_head  =>  r_int_head; s  =>  s_head; n_lines = 0                                   
           do; s  =>  s % next; if ( associated (s, s_head) ) exit                              
              if (s % status > 0) then                                                          
                 s % status = - s % status                                                      
                 do i=1,nd_d+1                                                                  
                    if ( s % opp(i) % spx % status > 0 ) then                                   
                       n_lines = n_lines+1; first = .true.                                      
                       write (idf,*)                                                            
                       write (idf,*) '% linecolor = ',dom_color,' markertype = ',marker         
                       do j=1,nd_d+1                                                            
                          if (j .ne. i) then                                                    
                             do k=1,nd_d                                                        
                                xpl(k) = s % opp(j) % pnt % x(k)                                
                             enddo                                                              
                             ipl = s % opp(j) % pnt % index                                     
                             write (idf,'(2e15.7,i8)') xpl, ipl                                 
                          endif                                                                 
                       enddo                                                                    
                    endif                                                                       
                 enddo                                                                          
              endif                                                                             
           enddo                                                                                

        else                                                                                    

           s_head  =>  r_int_head; s  =>  s_head; n_lines = 0                                   
           do; s  =>  s % next; if ( associated (s, s_head) ) exit                              
              if (s % status > 0) then                                                          
                 s % status = - s % status                                                      
                 do i=1,nd_d+1                                                                  
                    if ( s % opp(i) % spx % status > 0 ) then                                   
                       n_lines = n_lines+1; first = .true.                                      
                       write (idf,*)                                                            
                       write (idf,*) '% linecolor = ',dom_color,' markertype = ',marker         
                       j = mod(i,nd_d+1)+1                                                      
                       k = mod(j,nd_d+1)+1                                                      
                       p1 => s%opp(j)%pnt                                                       
                       p2 => s%opp(k)%pnt                                                       
                       length = segment_r_length(nd_d, p1, p2)                                  
                       if(length > tol_l) then                                                  
                          xpl = .5d0*(p1%x+p2%x)                                                
                          write (idf,*) '% linecolor = ', bad_color,' markertype = ', marker    

                          xpl = .5d0*(p1%x+p2%x)                                                
                          write (idf,*) '% linecolor = ', not_color,' markertype = ', marker    

                       else                                                                     
                          write (idf,*) '% linecolor = ', dom_color,' markertype = ', marker    
                       endif                                                                    
                       write (idf,'(2e15.7,i8)') p1%x, p1%index                                 
                       write (idf,'(2e15.7,i8)') p2%x, p2%index                                 
                    endif                                                                       
                 enddo                                                                          
              endif                                                                             
           enddo                                                                                

        endif                                                                                   

        100 format(a10,e15.7,1x,a3,e15.7,1x,a10,f6.2)                                               

        s  =>  s_head                                                                           
        do; s  =>  s % next; if ( associated (s, s_head) ) exit                                 
           s % status = abs(s % status)                                                         
        enddo                                                                                   

        s_head  =>  r_ext_head; s  =>  s_head                                                   
        do; s  =>  s % next; if ( associated (s, s_head) ) exit                                 
           s % status = abs(s % status)                                                         
        enddo

        10 do e = 1,grid % number_of_edges                                                           

           if ( grid % edges(e) % struct == 0 ) then
              s  =>  grid % edges(e) % int_head                                                 
              do; s  =>  s % next; if ( associated (s, grid % edges(e) % int_head) ) exit       

                 n_lines = n_lines+1                                                            

                 write (idf,*)                                                                  
                 write (idf,*) '% linecolor = ',bou_color,' markertype = ',marker               
                 do k=1,nd_d                                                                    
                    xpl(k) = s % opp(1) % pnt % ppp % x(k)                                      
                 enddo                                                                          
                 ipl = s % opp(1) % pnt % ppp % index                                           
                 write (idf,'(2e15.7,i8)') xpl, ipl                                             
                 do k=1,nd_d                                                                    
                    xpl(k) = s % opp(2) % pnt % ppp % x(k)                                      
                 enddo                                                                          
                 ipl = s % opp(2) % pnt % ppp % index                                           
                 write (idf,'(2e15.7,i8)') xpl, ipl                                             

              enddo                                                                             
           endif
        enddo                                                                                   

      end subroutine g_plot                                                                     


      end module grid_2d                                                                        



      function last_c_leng (len_str, string) result (leng)
!     +---------------------------------------------------------------+  
      implicit none

      integer, intent(in) :: len_str
      character (len=len_str), intent(in) :: string
      integer :: leng

      integer :: i
!     +---------------------------------------------------------------+  

      leng = len_str

      do i=1,len_str
        if ( string(i:i) .eq. ' ' ) then
          leng = i-1; exit
        end if
      end do

      end function last_c_leng

