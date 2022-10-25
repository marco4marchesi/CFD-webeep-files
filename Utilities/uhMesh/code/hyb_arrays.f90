  
  MODULE arrays

  !---------------------------------------------------------------------------------------
  USE convert
  USE schemi

  IMPLICIT NONE

  TYPE (domain), POINTER, PUBLIC :: array

  INTEGER, PUBLIC :: np_d_previous = 0, &
                     ns_d_previous = 0

  INTEGER, PUBLIC :: np_d,   &   ! Nj_d UNS
                     ns_d,   &   ! Nm_d UNS
 		     np_b,   &   ! Nj_b UNS
		     ns_b,   &   ! Nm_b UNS
		     g_np_d, &   ! Nj_d HYB
		     g_ns_d, &   ! Nm_d HYB
		     g_np_b      ! Nm_b HYB
		     
  ! NOTA 1 - Il numero di nodi di contorno Nj_b HYB e' pari a 2*g_np_b
  ! NOTA 2 - NP_B, dop arrays contiene per ogni lato anche i punti di box
  ! tre lati di 5 punti? => np_b = 21
  
  INTEGER, DIMENSION(:),             &
  ALLOCATABLE, PUBLIC :: jjdir,      &	 ! jd_jb UNS
			 b_jjdir,    &   ! jd_jb HYB
			 iflux,      &   ! bound_m UNS
			 b_bound_p,  &   ! bound_p/m HYB
			 g_m_t,      &   ! ele_type_d HYB
                         js_jns,     &
                         neighs,     &
			 neighs_low, &
			 m_t,	     &
                         iparent,    &
                         idir

  INTEGER, DIMENSION(:,:),           &
  ALLOCATABLE, PUBLIC :: jjs_b,      &   ! j_m_b UNS
			 b_jjs,      &   ! j_m_b HYB
			 neighs_b,   &   ! ma_m_b UNS
                         b_neighs_b, &   ! ma_m_b HYB
                         neigh,      &   ! ma_m_d UNS
			 g_neigh,    &   ! ma_m_d HYB
			 jloc,       &
			 mloc,       &
			 adj,	     &
			 jjs,	     &
			 b_jjs_b
						  
  INTEGER, DIMENSION(:,:), &
  POINTER, PUBLIC :: jj,   &   ! j_m_d UNS
                     g_jj, &   ! j_m_d HYB
		     s_jj, &
		     jj_previous
  
  REAL (KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: rr_b
  REAL (KIND=8), DIMENSION(:,:), POINTER,     PUBLIC :: rr,   &   ! rr UNS
	   	 				        g_rr, &   ! rr HYB
							s_rr, &
						        rr_previous

  PUBLIC  :: build_arrays, smooth
  PRIVATE :: p_vect_domain, p_vect_boundary, &
             s_vect_domain, s_vect_boundary
  !---------------------------------------------------------------------------------------
!
!
      contains
!
!
      subroutine structured_arrays
!
!     +---------------------------------------------------------------------------+
!     |  p % what indica il tipo di elemento che appoggia                         |
!     |           il proprio vertice di inizio su p;                              |
!     |           l'elemento si trova seguendo il verso convenzionale di          |
!     |           percorrenza della lista in cui si trova p                       |
!     |           sono previste 18 tipologie di elemento:                         |
!     |  0   nessun elemento o un elemento invisibile per interferenza            |
!     |  1   un rettangolo circondato da 4 elementi                               |
!     |  2   un rettangolo circondato da 3 elementi-esposto a destra              |
!     |  3   un rettangolo circondato da 3 elementi-esposto a sinistra            |
!     |  4   -                                                                    |
!     |  5   -                                                                    |
!     |  6   un triangolo di condensazione                                        |
!     |  7   due triangoli di condensazione(disabilitato=>6x2)                    |
!     |  8   un triangolo di spigolo + il quadrilatero a fianco (dx o sx)         |
!     |  9   due triangoli di spigolo + il quadrilatero seguente                  |
!     |  10  un quadrilatero suddiviso in tre parti (rarefazione)                 |
!     |         e circondato da 5 elementi                                        |
!     |  11  10 rovesciato (disabilitato)                                         |
!     |  12  un triangolo di rampa a sinistra                                     |
!     |  13  un triangolo di rampa a destra                                       |
!     |  14  una serie di elementi di scia                                        |
!     |  15  un rettangolo circondato da 1 elemento-esposto a dx, sopra, e a sx   |
!     |  16  un un quadrilatero diviso in tre parti                               |
!     |         e circondato esposto a dx ed eventualmente anche sopra            |
!     |  17  un quadrilatero diviso in tre parti ed esposto a sx/sopra            |
!     |  18  un quadrilatero diviso in tre parti esposto a dx, sx, sopra          |
!     |  19  Quattro triangoli di spigolo + il quadrilatero seguente              |
!     +---------------------------------------------------------------------------+
      implicit none
!
      type(point), pointer :: h, h1, h2
      type(point), pointer :: p, q, r, wkp
!
      integer :: m, i, j, k, l, w, z
      integer :: k_in, nlayer, wake_el
      integer :: sum, ind
      integer :: up_go, up_idx, up_typ, up_elms
!     +----------------------------------------------------------------------------+
!
!     +------------------------------------------------+
!     |  z is the number of lists of points points in  |
!     |  the first layer                               |
!     +------------------------------------------------+
!
      z = size( layer(1) % pnt_list )      
      nlayer = b_grid % number_of_layers
!
      allocate ( jloc(nlayer,z) )
      allocate ( mloc(nlayer,z) )
!
!     +--------------------------------------------------+
!     |  Loop over each list of nodes node (j) of        |
!     |  each layer (i) to count the number of elements  |
!     +--------------------------------------------------+
!
      k = 0
      sum = 0
!
      do i = 1, nlayer
        do j = 1, z
!
!     +------------------------------------------------+
!     |  mloc e' l'indice del primo elemento visibile  |
!     |  sulla lista j dello strato i                  |
!     +------------------------------------------------+
!
          mloc(i,j) = k + 1
          h => layer(i) % pnt_list(j)
          p => h          
!
!     +--------------------------------------------------+
!     |  The following loop counts the total number of   |
!     |  elements of the grid according to the variable  |
!     |  p % what that specifies for the first node of   |
!     |  the element the number of elements attached     |
!     +--------------------------------------------------+
!
          do; select case ( p % what )
!
            case ( 0)
            case ( 1);  k = k + 1    
            case ( 2);  k = k + 1    
            case ( 3);  k = k + 1    
            case ( 4)
            case ( 5)
            case ( 6);  k = k + 1    
            case ( 7)
            case ( 8);  k = k + 2
            case ( 9);  k = k + 3
            case (10);  k = k + 3
            case (11)
            case (12);  k = k + 1    
            case (13);  k = k + 1    
!
            case (14)
            if( (i .eq. 1) .and. (associated(p % wake)) ) then                           
               k_in = k                                                         
               wake_el = p % ppp1 % index - p % ppp % index
               k = k_in + wake_el + 1                          
            end if
!
            case (15);  k = k + 1     
            case (16);  k = k + 3     
            case (17);  k = k + 3     
            case (18);  k = k + 3
            case (19);  k = k + 5 
!
            end select
!
            p => p % next
            if ( associated(p,h) ) exit                                                
!
          end do
!
          layer(i) % num_el(j) = k - sum
          sum = sum + layer(i) % num_el(j)
!
          if ( layer(i) % num_el(j) .eq. 0 )  mloc(i,j) = 0
          if ( i .eq. nlayer )  layer(i+1) % num_el(j) = 0
!
        end do
      end do
!
!     +-----------------------------------------------+
!     |  k is now the total number of element in the  |
!     |  quadrilateral layers                         |
!     +-----------------------------------------------+
!
      m = k
!
!     +---------------------------------------------------------+
!     |  Tipo dell'elemento m-esimo (triangolo o quadrilatero)  |
!     |  quadrilatero = +3                                      |
!     |  triangolo con vertice all'insu' = -2                   |
!     |  triangolo con vertice all'ingiu' = +2                  |
!     +---------------------------------------------------------+
!
      allocate ( m_t(m) ) 
      m_t = 0
!
!     +-----------------------------------------+
!     |  Indici dei nodi dell'elemento m-esimo  |
!     +-----------------------------------------+
!
      allocate ( s_jj(4,m) ) 
      s_jj = 0
!
!     +-----------------------------------+
!     |  Coordinate dei punti con indice  |
!     +-----------------------------------+
!
      allocate ( s_rr(2,number) ) 
      s_rr = 0
!
!     +---------------------------------------------+
!     |  Loop over each list of nodes node (j) of   |
!     |  each layer (i) to define the connectivity  |
!     +---------------------------------------------+
!
      k = 0
      do i = 1, nlayer
        do j = 1, z
!
!     +------------------------------------------------+
!     |  h1 and h2 are the list of points along the    |
!     |  layer i and layer i + 1 that will define the  |
!     |  connectivity of the elements                  |
!     +------------------------------------------------+
!
          h1 => layer(i) % pnt_list(j)
          h2 => layer(i+1) % pnt_list(j)
!
          if ( layer(i) % num_el(j) .eq. 0 ) cycle
!
          p => h1
          q => h1 % next
          r => h2            
!
          do;  select case ( p % what )
!         
            case ( 0)
!
            if ( associated( r % old1 ) ) then
               if ( associated( r % old1, r % prev % old ) ) then
                  p => p % prev
                  q => q % prev
               else
                  if ( associated ( r % old, q % next ) ) then
                     p => p % next % next
                     q => p % next
                     cycle
                  else if ( associated ( r % old1 % next, q ) ) then
                     p => p % next
                     q => p % next
                     cycle
                  end if
               end if
            end if
!
            r => r % next
!
            case ( 1)
            k = k + 1
            m_t(k) = 3

            if ( r % what .eq. 0 ) r % mark = k

            s_jj(1,k) = p % index
            s_jj(2,k) = r % index
            r => r % next
            s_jj(3,k) = r % index
            s_jj(4,k) = q % index
!
            case ( 2)
            k = k + 1
            m_t(k) = 3

            p % mark = k
            if(r % what == 0) r % mark = k

            s_jj(1,k) = p % index
            s_jj(2,k) = r % index
            r => r % next
            s_jj(3,k) = r % index
            s_jj(4,k) = q % index
!
            case ( 3)
            k = k + 1
            m_t(k) = 3

            if(r % what == 0) r % mark = k

            s_jj(1,k) = p % index
            s_jj(2,k) = r % index
            r => r % next
            s_jj(3,k) = r % index
            s_jj(4,k) = q % index

            r % mark = k
!
            case ( 6)
            k = k + 1
            m_t(k) = -2
            s_jj(1,k) = p % index
            s_jj(2,k) = r % index
            s_jj(3,k) = q % index
!
            case ( 8)
            k = k + 2
            m_t(k-1) = 2
            m_t(k) = 3

            if ( r % what .eq. 0 ) r % mark = k-1

            s_jj(1,k-1) = p % index         
            s_jj(2,k-1) = r % index       
            r => r % next
            s_jj(3,k-1) = r % index
!
            if ( r % what .eq. 0 ) r % mark = k

            s_jj(1,k) = p % index  
            s_jj(2,k) = r % index
            r => r % next 
            s_jj(3,k) = r % index
            s_jj(4,k) = q % index
!
            case(10)                                                    
            k =  k + 3                                                  
            m_t(k-2) =  2                                               
            m_t(k-1) = -2                                               
            m_t(k)   =  2                                               

            if(r % what == 0) r % mark = k-2                            

            s_jj(1,k-2) = p % index                                     
            s_jj(2,k-2) = r % index                                     
            r => r % next                                               
            s_jj(3,k-2) = r % index                                     

            s_jj(1,k-1) = p % index                                     
            s_jj(2,k-1) = r % index                                     
            s_jj(3,k-1) = q % index                                     

            if(r % what == 0) r % mark = k                              

            s_jj(1,k) = r % index                                       
            r => r % next                                               
            s_jj(2,k) = r % index                                       
            s_jj(3,k) = q % index                                       
           
            case(12)                                                    
               k = k + 1                                                
               m_t(k) = -2                                              

               p % mark = k                                             
               s_jj(1,k) = p % index                                    
               r => r % next                                            
               s_jj(2,k) = r % index                                    
               s_jj(3,k) = q % index                                    

            case(13)                                                    
               k = k + 1                                                
               m_t(k) = -2                                              

               r % mark = k                                             
               s_jj(1,k) = p % index                                    
               s_jj(2,k) = r % index                                    
               s_jj(3,k) = q % index                                    
               r => r % next                                            



            case(9)                                                     
               k = k + 3                                                
               m_t(k-2) = 2                                             
               m_t(k-1) = 2                                             
               m_t(k)   = 3                                             

               if(r % what == 0) r % mark = k-2                         

               s_jj(1,k-2) = p % index                                  
               s_jj(2,k-2) = r % index                                  
               r => r % next                                            
               s_jj(3,k-2) = r % index                                  

               if(r % what == 0) r % mark = k-1                         

               s_jj(1,k-1) = p % index                                  
               s_jj(2,k-1) = r % index                                  
               r => r % next                                            
               s_jj(3,k-1) = r % index                                  

               if(r % what == 0) r % mark = k                           

               s_jj(1,k) = p % index                                    
               s_jj(2,k) = r % index                                    
               r => r % next                                            
               s_jj(3,k) = r % index                                    
               s_jj(4,k) = q % index                                    

            case(19)                                                     
               k = k + 5
               m_t(k-4) = 2                                             
               m_t(k-3) = 2                                             
               m_t(k-2) = 2
               m_t(k-1) = 2
               m_t(k)   = 3                                             

               if(r % what == 0) r % mark = k-4                         

               s_jj(1,k-4) = p % index                                  
               s_jj(2,k-4) = r % index                                  
               r => r % next                                            
               s_jj(3,k-4) = r % index                                  

               if(r % what == 0) r % mark = k-3                         

               s_jj(1,k-3) = p % index                                  
               s_jj(2,k-3) = r % index                                  
               r => r % next                                            
               s_jj(3,k-3) = r % index                                  

               if(r % what == 0) r % mark = k-2                         

               s_jj(1,k-2) = p % index                                  
               s_jj(2,k-2) = r % index                                  
               r => r % next                                            
               s_jj(3,k-2) = r % index                                  

               if(r % what == 0) r % mark = k-1                         

               s_jj(1,k-1) = p % index                                  
               s_jj(2,k-1) = r % index                                  
               r => r % next                                            
               s_jj(3,k-1) = r % index                                  

               if(r % what == 0) r % mark = k                           

               s_jj(1,k) = p % index                                    
               s_jj(2,k) = r % index                                    
               r => r % next                                            
               s_jj(3,k) = r % index                                    
               s_jj(4,k) = q % index                                    




            case(14)

               q => p % ppp                                             
               if(associated(p % wake))then                             
                  r => p % wake                                         
                  k_in = k                                              

                  do                                                    
                     k = k + 1                                          
                     m_t(k) = 3                                         

                     if ( associated(r % next) ) then                       
                        r => r % next                                   
                     else
!
!     +-----------------------------------------+
!     |  Four triangles at the end of the wake  |
!     +-----------------------------------------+
!
                        k = k + 4                                       
                        m_t(k-3) = 2
                        m_t(k-2) = 2
                        m_t(k-1) = 2
                        m_t(k) = 2                                      
!
!     +--------------------------------------+
!     |  Returning quadrilaters in the wake  |
!     +--------------------------------------+
!
                        do w = 1,(k - k_in - 4)                             
                          m_t(k + w) = 3
                        enddo                                           

                        m_t(k+w-1)= 3 ! per identificare la fine scia   
                        m_t(k+w) = 3                                    

                        exit                                            
                     endif                                              

                  enddo                                                 
                  k = k + w                                             

                  wake_el = p % ppp1 % index - p % ppp % index
                  q => p % ppp1                                         

                  !quadrilatero sull'altra faccia dello spigolo         
                  if ( q % what .eq. 0 ) then                                 
                     q % mark = k_in + wake_el + 1                      
                  endif                                                 
                  s_jj(1,k_in+wake_el+1) = p % index                    
                  s_jj(2,k_in+wake_el+1) = q % index                    
                  s_jj(3,k_in+wake_el+1) = q % next % index             
                  s_jj(4,k_in+wake_el+1) = p % next % index             

                  q => p % wake                                         
                  do w = 1,(wake_el - 4) / 2                                
                     
                     ! rettangoli 'sopra'                               
                     if(q % ppp1 % what == 0)then                       
                        q % ppp1 % mark = k_in + wake_el + 1 - w        
                     endif                                              
                     s_jj(1,k_in+wake_el+1-w) = q % index               
                     s_jj(2,k_in+wake_el+1-w) = q % ppp1 % index        
                     s_jj(3,k_in+wake_el+1-w) = q % ppp1 % next % index
                     s_jj(4,k_in+wake_el+1-w) = q % old % index         

                     ! rettangoli 'sotto'                               
                     if(q % old % ppp % what == 0)then                  
                        q % old % ppp % mark = k_in + w                 
                     endif                                              
                     s_jj(1,k_in+w) = q % old % index                   
                     s_jj(2,k_in+w) = q % old % ppp % index             
                     s_jj(3,k_in+w) = q % ppp % index                   
                     s_jj(4,k_in+w) = q % index
                                           
                     if(.not.(associated(q % next)))then                
                        exit                                            
                     else                                               
                        q => q % next                                   
                     endif
                                                                   
                  enddo                                                 

                  ! quattro triangoli finali                                
                  if(q % ppp % what == 0)then                           
                     q % ppp % mark = k_in + w + 1                      
                  endif                                                 
                  s_jj(1,k_in+w+1) = q % index                          
                  s_jj(2,k_in+w+1) = q % ppp % index                    
                  s_jj(3,k_in+w+1) = q % ppp % next % index
                  
                  wkp => q % ppp % next

                  if(q % ppp % what == 0)then                           
                     q % ppp % mark = k_in + w + 2
                  endif                                                 
                  s_jj(1,k_in+w+2) = q % index                          
                  s_jj(2,k_in+w+2) = wkp % index            
                  s_jj(3,k_in+w+2) = wkp % next % index               
               
                  wkp => wkp % next

                  if(q % ppp % what == 0)then                           
                     q % ppp % mark = k_in + w + 3                      
                  endif                                                 
                  s_jj(1,k_in+w+3) = q % index                          
                  s_jj(2,k_in+w+3) = wkp % index                    
                  s_jj(3,k_in+w+3) = wkp % next % index             
               
                  if(q % ppp1 % prev % what == 0)then                   
                     q % ppp1 % prev % mark = k_in + w + 4
                  endif                                                 
                  s_jj(1,k_in+w+4) = q % index                          
                  s_jj(2,k_in+w+4) = q % ppp1 % prev % index            
                  s_jj(3,k_in+w+4) = q % ppp1 % index
                                     
                  r => p % ppp1 % next                                  

              endif

            case(15)                                                    
               k = k + 1                                                
               m_t(k) = 3                                               
               p % mark = k                                             
               if(r % what == 0) r % mark = k                           

               r % next % mark = k                                      

               s_jj(1,k) = p % index                                    
               s_jj(2,k) = r % index                                    
               r => r % next                                            
               s_jj(3,k) = r % index                                    
               s_jj(4,k) = q % index                                    
            case(16)                                                    
               k =  k + 3                                               
               m_t(k-2) =  2                                            
               m_t(k-1) = -2                                            
               m_t(k)   =  2                                            

               p % mark = k-2                                           

               if(r % what == 0) r % mark = k-2                         

               s_jj(1,k-2) = p % index                                  
               s_jj(2,k-2) = r % index                                  
               r => r % next                                            
               s_jj(3,k-2) = r % index                                  

               s_jj(1,k-1) = p % index                                  
               s_jj(2,k-1) = r % index                                  
               s_jj(3,k-1) = q % index                                  

               if(r % what == 0) r % mark = k                           

               s_jj(1,k) = r % index                                    
               r => r % next                                            
               s_jj(2,k) = r % index                                    
               s_jj(3,k) = q % index                                    
            case(17)                                                    
               k =  k + 3                                               
               m_t(k-2) =  2                                            
               m_t(k-1) = -2                                            
               m_t(k)   =  2                                            

               if(r % what == 0) r % mark = k-2                         

               s_jj(1,k-2) = p % index                                  
               s_jj(2,k-2) = r % index                                  
               r => r % next                                            
               s_jj(3,k-2) = r % index                                  

               s_jj(1,k-1) = p % index                                  
               s_jj(2,k-1) = r % index                                  
               s_jj(3,k-1) = q % index                                  

               if(r % what == 0) r % mark = k                           

               r % next % mark = k                                      
               s_jj(1,k) = r % index                                    
               r => r % next                                            
               s_jj(2,k) = r % index                                    
               s_jj(3,k) = q % index                                    
            case(18)                                                    
               k =  k + 3                                               
               m_t(k-2) =  2                                            
               m_t(k-1) = -2                                            
               m_t(k)   =  2                                            

               p % mark = k-2                                           
               if(r % what == 0) r % mark = k-2                         

               s_jj(1,k-2) = p % index                                  
               s_jj(2,k-2) = r % index                                  
               r => r % next                                            
               s_jj(3,k-2) = r % index                                  

               s_jj(1,k-1) = p % index                                  
               s_jj(2,k-1) = r % index                                  
               s_jj(3,k-1) = q % index                                  

               if(r % what == 0) r % mark = k                           

               r % next % mark = k                                      
               s_jj(1,k) = r % index                                    
               r => r % next                                            
               s_jj(2,k) = r % index                                    
               s_jj(3,k) = q % index                                    
            end select                                                  
!
           p => p % next
           q => p % next
!
           if ( associated(p,h1) ) exit
!
            enddo
         enddo

      enddo

    !call wait
    ! -------- coordinate dei punti -------------
    do i = 1,nlayer + 1
       do j = 1,z
          h => layer(i)%pnt_list(j)
          p => h
          do
             if(associated(p % wake))then
                q => p % wake
                s_rr(:,q % index) = q % x
                do
                   if(associated(q % next))then
                      q => q % next
                      s_rr(:,q % index) = q % x
                   else
                      exit
                   endif
                enddo
             endif
             ind = p % index
             if(ind > 0)then 
                s_rr(:,ind) = p % x
             endif
             p => p % next
             if(associated(p,h))exit
          enddo
       enddo
    enddo
    !@#@#@#@#@ ripetizione delle coordinate del contorno fittizio?????????


    !---------- vettore degli elementi adiacenti -------
    allocate(adj(m,4))
    adj = 0 
    !------------ interno striscia -----------
    k = 1
    do i = 1, nlayer
       do j = 1, z

          h1 => layer(i) % pnt_list(j)
          h2 => layer(i+1) % pnt_list(j)

          if ( layer(i) % num_el(j) .eq. 0 ) cycle

          up_elms = layer(i+1) % num_el(j)

          p => h1
          q => h1 % next

          if(up_elms > 0)then
             up_go = mloc(i+1,j)-1
             r => h2
          else
             up_go = 0
          endif

          do

             select case(p % what)
             case(0)
                if(associated(r % old1)) then
                   if(associated(r % old1, r % prev % old))then
                      p => p % prev; q => q % prev
                   else
                      if(associated(r % old, p % next % next))then
                         p => p % next % next; cycle
                      elseif(associated(r % old1 % next, p % next))then
                         p => p % next; cycle
                      endif
                   endif
                endif
                call adj_up(up_idx, up_typ)

             case(1)
                adj(k,1) = k - 1 
                adj(k,3) = k + 1
                if(up_elms > 0)then

                   call adj_up(up_idx, up_typ)

                   adj(k,2) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                endif
                k = k + 1

             case(2)
                adj(k,1) = 0
                adj(k,3) = k + 1

                if(up_elms > 0)then

                   call adj_up(up_idx, up_typ)

                   adj(k,2) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k 
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                endif
                k = k + 1

             case(3)
                adj(k,1) = k - 1
                adj(k,3) = 0

                if(up_elms > 0)then

                   call adj_up(up_idx, up_typ)

                   adj(k,2) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k 
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                endif

                k = k + 1

             case(6)
                adj(k,1) = k + 1
                adj(k,3) = k - 1
                k = k + 1

             case(8)
                adj(k,2) = k + 1 
                adj(k,3) = k - 1
                adj(k+1,1) = k
                adj(k+1,3) = k + 2

                if(up_elms > 0)then

                   call adj_up(up_idx, up_typ)

                   adj(k,1) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                   k = k + 1

                   call adj_up(up_idx, up_typ)

                   adj(k,2) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                else
                   k = k + 1
                endif
                k = k + 1

             case(9)
                adj(k,2) = k + 1
                adj(k,3) = k - 1
                adj(k+1,2) = k + 2
                adj(k+1,3) = k 
                adj(k+2,1) = k + 1
                adj(k+2,3) = k + 3

                if(up_elms > 0)then

                   call adj_up(up_idx, up_typ)

                   adj(k,1) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                   k = k + 1

                   call adj_up(up_idx, up_typ)

                   adj(k,1) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                   k = k + 1

                   call adj_up(up_idx, up_typ)

                   adj(k,2) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                else
                   k = k + 2
                endif
                k = k + 1 

             case(10)
                adj(k,2) = k + 1
                adj(k,3) = k - 1
                adj(k+1,1) = k + 2 
                adj(k+1,3) = k 
                adj(k+2,1) = k + 3
                adj(k+2,2) = k + 1

                if(up_elms > 0)then

                   call adj_up(up_idx, up_typ)

                   adj(k,1) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k 
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                   k = k + 2 

                   call adj_up(up_idx, up_typ)

                   adj(k,3) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                else
                   k = k + 2
                endif
                k = k + 1

             case(12)
                adj(k,1) = k + 1
                adj(k,3) = 0

                call adj_up(up_idx, up_typ)

                k = k + 1

             case(13)
                adj(k,1) = 0
                adj(k,3) = k - 1
                call adj_up(up_idx, up_typ)
                k = k + 1

             case(14)
                k_in = k
                wake_el = p % ppp1 % index - p % ppp % index 
                ! ricalcolato x' potrebbero esserci piu' scie!
                do w = 1, (wake_el-2)/2 !gli elementi sotto
                   adj(k,1) = k_in + w - 2
                   adj(k,3) = k_in + w
                   adj(k,4) = k_in + wake_el - w

                   call adj_up(up_idx, up_typ)

                   adj(k,2) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k 
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k 
                   endif

                   k = k + 1
                enddo
                ! i due triangoli finali
                adj(k,2) = k + 1         ! il primo sotto
                adj(k,3) = k - 1

                call adj_up(up_idx, up_typ)

                adj(k,1) = up_idx
                if(up_typ == 3)then
                   adj(up_idx,4) = k 
                elseif(up_typ == -2)then
                   adj(up_idx,2) = k 
                endif

                k = k + 1
                ! il secondo sopra
                adj(k,2) = k + 1
                adj(k,3) = k - 1

                call adj_up(up_idx, up_typ)

                adj(k,1) = up_idx
                if(up_typ == 3)then
                   adj(up_idx,4) = k 
                elseif(up_typ == -2)then
                   adj(up_idx,2) = k 
                endif

                k = k + 1

                k_in = k
                do w = 1, (wake_el-2)/2 !gli elementi sopra
                   adj(k,1) = k_in + w - 2
                   adj(k,3) = k_in + w
                   adj(k,4) = k_in - 2 - w

                   call adj_up(up_idx, up_typ)
                   adj(k,2) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k 
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k 
                   endif

                   k = k + 1
                enddo
                !l'elemento che sale sul bordo di uscita di un profilo alare (es.)
                adj(k,1) = k - 1 
                adj(k,3) = k + 1

                call adj_up(up_idx, up_typ)

                adj(k,2) = up_idx
                if(up_typ == 3)then
                   adj(up_idx,4) = k
                elseif(up_typ == -2)then
                   adj(up_idx,2) = k
                endif


                k = k + 1
             case(15) ! un dentino!!
                adj(k,1) = 0
                adj(k,3) = 0
                if(up_elms > 0)then

                   call adj_up(up_idx, up_typ)

                   adj(k,2) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k 
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                endif

                k = k + 1

             case(16)
                adj(k,2) = k + 1
                adj(k,3) = 0
                adj(k+1,1) = k + 2 
                adj(k+1,3) = k 
                adj(k+2,1) = k + 3
                adj(k+2,2) = k + 1

                if(up_elms > 0)then

                   call adj_up(up_idx, up_typ)

                   adj(k,1) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k 
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                   k = k + 2 

                   call adj_up(up_idx, up_typ)

                   adj(k,3) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                else
                   k = k + 2
                endif

                k = k + 1
             case(17)
                adj(k,2) = k + 1
                adj(k,3) = k - 1
                adj(k+1,1) = k + 2 
                adj(k+1,3) = k 
                adj(k+2,1) = 0!k + 3
                adj(k+2,2) = k + 1
                !if(i<nlayer)then
                if(up_elms > 0)then

                   call adj_up(up_idx, up_typ)

                   adj(k,1) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k 
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                   k = k + 2 

                   call adj_up(up_idx, up_typ)

                   adj(k,3) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                else
                   k = k + 2
                endif

                k = k + 1
             case(18)
                adj(k,2) = k + 1
                adj(k,3) = 0
                adj(k+1,1) = k + 2 
                adj(k+1,3) = k 
                adj(k+2,1) = 0
                adj(k+2,2) = k + 1

                if(up_elms > 0)then

                   call adj_up(up_idx, up_typ)

                   adj(k,1) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k 
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                   k = k + 2 

                   call adj_up(up_idx, up_typ)

                   adj(k,3) = up_idx
                   if(up_typ == 3)then
                      adj(up_idx,4) = k
                   elseif(up_typ == -2)then
                      adj(up_idx,2) = k
                   endif
                else
                   k = k + 2
                endif

                k = k + 1
             end select

             p => p % next
             if(associated(p, h1))exit
          enddo

          q => p % prev ! correzione per il primo elemento dello strato
          k = k - 1

          if((q% what == 1).or.(q% what == 14).or.((q% what>5).and.(q% what<13)))then
             l = mloc(i,j)

             if(m_t(k) == 3)then
                adj(k,3) = l
             elseif(m_t(k) == 2)then

                adj(k,1) = l ! si trovera' solo il triangolo finale del caso 10!!!
                ! questa e' un'ipotesi forte!!!!!
             elseif(m_t(k) == -2)then
                adj(k,1) = l
             endif

             if(m_t(l) == 3)then
                adj(l,1) = k
             elseif(m_t(l) == 2)then
                adj(l,3) = k
             elseif(m_t(l) == -2)then
                adj(l,3) = k
             endif

          endif
          k = k + 1 

       enddo
    enddo
  contains
    subroutine adj_up(idx, typ)

      implicit none

      integer, intent(out)    :: idx, typ

      up_go = up_go + 1
      idx = up_go

      select case(r % what)
      case(1)
         typ = 3
      case(0)
         typ = 0 ; idx = 0      ; up_go = up_go - 1
      case(10)
         typ = -2; idx = idx + 1; up_go = idx + 1
      case(6)
         typ = -2
      case(2)
         typ = 3
      case(3)
         typ = 3
      case(12)
         typ = -2; 
      case(13)
         typ = -2
      case(8)
         typ = 3; idx = idx + 1; up_go = idx
      case(9)
         typ = 3; idx = idx + 2; up_go = idx
      case(14) ! si potrebbe anche togliere
         typ = 3
      case(15)
         typ = 3
      case(16:18)
         typ = -2; idx = idx + 1; up_go = idx + 1
      end select
      r => r % next
    end subroutine adj_up
      
      end subroutine structured_arrays







  SUBROUTINE couple_arrays()
  !-------------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER :: i, j, k, num, m, n, z, ns, check, l, v, &
             nlayer, s_ns_d, index, t, ns_fin

  TYPE(point),   POINTER :: pnt_head, add_head, p
  TYPE(simplex), POINTER :: int_head, s

  INTEGER,DIMENSION(3,3) :: mat
  !-------------------------------------------------------------------------------------

    mat(1,1) = 1;  mat(1,2) = 2;  mat(1,3) = 3
    mat(2,1) = 2;  mat(2,2) = 3;  mat(2,3) = 1
    mat(3,1) = 3;  mat(3,2) = 1;  mat(3,3) = 2
    
    ! Costruzione del vettore di interfaccia js_jns;
    ! indice strutturato corrispondente all'indice non strutturato
    ! dello stesso punto.
    ! la numerazione delle griglie strutturata e non strutturata
    ! e' gia' completata; i punti di frontiera tra le due griglie sono condivisi e
    ! quindi hanno una numerazione diversa

    ALLOCATE (js_jns(np_d))
    js_jns = 0
    
    ! Usare p % option come contenitore 1D dell'indice che deriva dalla 
    ! numerazione dell'ultimo layer

    DO i = 1, grid % number_of_edges
    
       ! Lato aperto e strutturato
       IF ((grid % edges(i) % begin_vert /= grid % edges(i) % end_vert) .AND. &
           (grid % edges(i) % struct == 1)) THEN

          index = grid % edges(i)% pnt_head % next % option
          grid % edges(i)% next % pnt_head % prev % option = index
          ! In pnt_head, il next e' la coda di add_head!!!!!

          index = grid % edges(i)% pnt_head % prev % option
          grid % edges(i)% prev % pnt_head % next % option = index

       ENDIF

    ENDDO
    
    
    DO i = 1, grid % number_of_edges

       pnt_head => grid % edges(i) % pnt_head
       p => pnt_head

       DO 
          p => p % next
          IF (ASSOCIATED(p, pnt_head)) EXIT

          js_jns(p % ppp % index) = p % option

       ENDDO

    ENDDO

    j = 0

    DO i = 1, grid % number_of_edges
    
       p => grid % edges(i) % pnt_head
       pnt_head => p % next % next
       DO
          j = j + 1
          p => p % prev
          IF (ASSOCIATED(p % old1))  neighs_low(j) = p % old1 % mark

          IF (ASSOCIATED(p, pnt_head)) EXIT
       ENDDO 

    ENDDO
    !-----------------------------------
    ! assegnare la nuova numerazione agli elementi di js_jns con indice nullo
    num = number
    !indice di partenza della numerazione sospesa all'ultimo layer
    DO i = 1,SIZE(js_jns)
       IF(js_jns(i) == 0)THEN
          num = num + 1
          js_jns(i) = num
       ENDIF
    ENDDO

    !-----------------------------------
    ! aggiornare gli indici dei nodi appartenenti ai simplessi
    ! il vettore si compila adessoo!!!
    DO i = 1,ns_d

       DO k = 1,nd_d+1 

          jj(k,i) = js_jns(jj(k,i))
       ENDDO

    ENDDO

    !-----------------------------------
    ! aggiornare la corrispondenza tra nodi 1D e 2D???????
    p => grid % pnt_head
    DO; p => p % next; IF(ASSOCIATED(p, grid % pnt_head))EXIT
       p % visited =.FALSE.
    ENDDO

    int_head => grid % int_head
    s => int_head
    DO;s => s % next;IF(ASSOCIATED(s,int_head))EXIT

       DO i = 1,nd_d + 1

          p => s % opp(i) % pnt

          IF(.NOT.(p % visited ))THEN

             p % index = js_jns(p % index)
             p % visited = .TRUE.

          ENDIF

       ENDDO

    ENDDO

    !----------------------------------
    ! unire i vettori delle coordinate(e degli indici!) rr,s_rr

    g_np_d = num

    ALLOCATE(g_rr(2,g_np_d))
    i = 0

    DO  
       i = i + 1
       ! coordinate dei punti strutturati
       g_rr(1,i) = s_rr(1,i)
       g_rr(2,i) = s_rr(2,i)

       IF( i == SIZE(s_rr,2))EXIT

    ENDDO

    i = 0
    DO 
       i = i + 1
       ! coordinate dei punti non strutturati
       g_rr(:,js_jns(i)) = rr(:,i)
       ! viene ripetuta l'assegnazione delle coordinate dei punti di contorno!

       IF(i == SIZE(js_jns))EXIT

    ENDDO

    ! WRITE(*,*)'unione vettori coordinate'
    !-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    ! generare il vettore jjdir
    ! indice di dominio corrispondente al punto di contorno

    g_np_b = 0

    DO i = 1,b_grid%number_of_edges
       g_np_b = g_np_b + numberof_point(b_grid % edges(i) % add_head)-1
    ENDDO

    ALLOCATE(b_jjdir(2*g_np_b))
    ALLOCATE(b_bound_p(2*g_np_b)) 

    j = 0  
    DO i = 1,b_grid % number_of_edges 
       v = 0

       add_head => b_grid % edges(i) % add_head
       p => add_head

       DO;p => p % next;IF(ASSOCIATED(p,add_head%prev))EXIT 
          v = v + 1
          j = j + 1                       ! non ha un ppp associato
          IF(((p % cnd == 0).AND.(b_grid % edges(i) % struct == 0)).OR. &
               (b_grid % edges(i) % struct >= 1))THEN

             b_jjdir(j) = p % ppp % index !punto di spalla
             b_bound_p(j) = i

          ELSEIF((p % cnd /= 0).AND.(b_grid % edges(i) % struct == 0))THEN

             b_jjdir(j) = p % ppp % ppp % index !punto intatto
             b_bound_p(j) = ABS(p % cnd)

             IF(ABS(p % cnd) /= i)THEN
                WRITE(*,*)'questa poi...cnd /= i,hybarr'
                STOP
             ENDIF

          ENDIF

          IF(v > 1) THEN ! si completa la coppietta 

             j = j+1
             b_jjdir(j) = b_jjdir(j-1)
             b_bound_p(j) = b_bound_p(j-1)

             ! due punti di contorno puntano allo stesso punto di dominio
          ENDIF
       ENDDO

       !write(*,*)'trattamento particolare'!---------------------------------------
       ! trattamento particolare per l'ultimo punto...
       j = j + 1
       b_bound_p(j) = i
       IF(ASSOCIATED(b_grid % edges(i) % next))THEN
          !write(*,*)'il lato continua'
          !il lato continua

          DO k = 1,b_grid %number_of_edges
             ! ricerca del lato contiguo
             IF( b_grid % edges(i) % end_vert == b_grid % edges(k) % begin_vert)EXIT
          ENDDO

          IF(b_grid % edges(k) % struct == 0)THEN

             !riemann normale e corretto
             IF(b_grid % edges(k) % add_head % next % cnd == 0)THEN
                ! riemann corretto
                b_jjdir(j) = b_grid % edges(k)% add_head % next % ppp % index
             ELSE 
                ! riemann normale
                b_jjdir(j) = b_grid % edges(k)% add_head % next % ppp % ppp % index
             ENDIF

          ELSE

             b_jjdir(j) = b_grid % edges(k)% add_head % next % ppp % index

          ENDIF

       ELSE

          DO k = 1,b_grid %number_of_edges
             ! ricerca del lato contiguo
             IF( b_grid % edges(i) % end_vert == b_grid % edges(k) % begin_vert)EXIT

          ENDDO

          b_jjdir(j) = b_grid % edges(k)% add_head % next % ppp % index

          IF((p % cnd /= 0).AND.(b_grid%edges(i)%struct == 0))THEN

             b_jjdir(j)=b_grid % edges(i) %  add_head % prev % ppp %ppp% index

          ENDIF

       ENDIF
       !------------------------------------------------------------------------
    ENDDO

    !-----------------------------------
    ! vettori jjs_b  e neighs_b 
    ! indice del punto appartenente al simplesso di contorno
    ! adiacenti 1D dei simplessi di contorno

    ALLOCATE(b_jjs_b(nd_b+1,g_np_b))
    ALLOCATE(b_neighs_b(nd_b +1,g_np_b))
    j = 1
    k = 0
    DO i = 1,b_grid % number_of_edges
       add_head => b_grid % edges(i) % add_head
       p => add_head % next

       k = k + 1
       b_neighs_b(1,k) = k+1
       b_neighs_b(2,k) = 0 ! ingresso
       b_jjs_b(1,k) = j
       b_jjs_b(2,k) = j+1

       j = j + 2
       DO 
          p => p % next
          k = k+1
          b_neighs_b(1,k) = k+1
          b_neighs_b(2,k) = k-1 ! include anche l'estremo iniziale(0)
          b_jjs_b(1,k) = j
          b_jjs_b(2,k) = j+1
          j = j + 2

          IF(ASSOCIATED(p,add_head % prev%prev))THEN! estremo finale(0)
             b_neighs_b(1,k) = 0
             EXIT
          ENDIF

       ENDDO
    ENDDO

    !-------------------------------------
    !interfaccia simplessi adiacenti
    ! lati chiusi, aperti , Delaunay 

    nlayer = b_grid % number_of_layers

    m = SIZE(m_t)
    s_ns_d = m

    CALL set_s2i(s_ns_d,grid % int_head)

    !rinumera i simplessi a partire dal numero di elementi strutturati

    n = 0
    DO i = 1,grid % number_of_edges
       CALL s_vect_boundary(grid % edges(i) % int_head,i,n)
    ENDDO
    ! per aggiornare neighs
    ! aggiornamento dei soli simplessi basati su lati strutturati



    DO i = 1,size(neigh,2)
       DO k = 1,3
          ! rinumera i simplessi 2D
          IF(neigh(k,i) /= 0)THEN

             neigh(k,i) = neigh(k,i) + m

          ENDIF

       ENDDO
    ENDDO

    z = 0
    ns = 0
    ns_fin = 0

    DO i = 1, grid % number_of_edges
       IF (grid % edges(i) % struct == 1) THEN
          
          z = z + 1
          ns = ns_fin

          DO
             ns = ns + 1
             IF (iflux(ns) == i) EXIT
             ns_fin = ns_fin + 1
          ENDDO
          ! ho trovato l'estremo iniziale
          DO
             ns_fin = ns_fin + 1
             IF (ns_fin == SIZE(iflux)) THEN     
                EXIT      
             ELSE
                IF (iflux(ns_fin + 1) > i) EXIT
             ENDIF
          ENDDO
          ! ho trovato l'estremo finale
          
          DO t = ns, ns_fin  ! sugli ns
          
             check = 0
             
             DO k = 1, 3 ! gira sugli opposti del simplesso
                IF (neigh(k, neighs(t)-m) == 0) THEN  
                
                   j = neighs_low(t)
                   
                   !print*, j
                   
                   neigh(k, neighs(t)-m) = j
                   
                   DO l = 1, SIZE(adj,2)
                      IF (adj(j,l) == 0) THEN
                         adj(j,l) = neighs(t)
                         EXIT
                      ENDIF
                   ENDDO
                   
                   EXIT
                
                ENDIF
             ENDDO
             
          ENDDO
          
          ! Lato aperto: e' necessario scambiare gli indici, perche' un quadrilatero 
          ! dell'ultimo strato di elementi E all'inizio di un lato aperto ha inizialmente 
          ! due indici nulli, per gli adiacenti a destra e sopra, ma solo quello sopra 
          ! e' interfacciato con un triangolo.           
          IF (grid % edges(i) % begin_vert /= grid % edges(i) % end_vert) THEN
             adj(neighs_low(ns),2) = adj(neighs_low(ns),1)
             adj(neighs_low(ns),1) = 0
          ENDIF

       ENDIF ! edges strutturati
    ENDDO ! sugli edges

    !@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
    ! UNIRE JJ indice dei punti appartenenti all' elemento   
    ALLOCATE( g_jj(4 , SIZE(s_jj,2)+SIZE(jj,2) ) )

    !quadrilateri
    DO i = 1,SIZE(s_jj,2)
       DO j = 1,ABS(m_t(i))+1
          g_jj(j,i) = s_jj(j,i)
       ENDDO
    ENDDO

    !triangoli
    DO k = i, i+SIZE(jj,2)-1
       DO j = 1,3
          g_jj(j,k) = jj(j,k-i+1)    
       ENDDO
    ENDDO

    !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    ! UNIRE  NEIGH ADJ  indici degli elementi adiacenti   
    !aggiornamento della numerazione contenuta nel vettore neigh degli adiacenti!!
    ALLOCATE( g_neigh(4,(SIZE(adj,1)+ SIZE(neigh,2)) ) )   
    !quadrilateri
    DO i = 1, SIZE(adj,1)
       DO j = 1,ABS(m_t(i))+1
          g_neigh(j,i) = adj(i,j)
       ENDDO
    ENDDO

    !triangoli
    DO k = i,i+SIZE(neigh,2)-1
       DO j = 1,3
          g_neigh(j,k) = neigh(j,k-i+1)
       ENDDO
    ENDDO

    !-#-#-#-##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#--#-#-#-#--#-
    ! ESTENDERE IL VETTORE M_T tipo di elemento
    ALLOCATE (g_m_t( SIZE(m_t) + ns_d))
    g_ns_d =  SIZE(m_t) + ns_d

    DO i = 1,SIZE(m_t)
       g_m_t(i) = ABS(m_t(i)) ! quadrilateri e triangoli
    ENDDO

    DO k = i,g_ns_d
       g_m_t(k) = 2  ! solo triangoli 
    ENDDO

  END SUBROUTINE couple_arrays
  
  


  
  SUBROUTINE build_arrays (grid)

    IMPLICIT NONE

    TYPE (domain), TARGET, INTENT(IN) :: grid
    INTEGER :: e, n

    array  =>  grid

    ! Initialization of  np_d_previous is needed to impose NULL status to
    ! previous grid pointers (otherwise error condition)

    IF (np_d_previous == 0) THEN
       NULLIFY(jj)
       NULLIFY(rr)
       NULLIFY(jj_previous)
       NULLIFY(rr_previous)
    ENDIF

    np_d_previous = np_d; ns_d_previous = ns_d


    np_d = 0; ns_d = 0 
    np_b = 0; ns_b = 0

    !progressive index assignement and elements count
    CALL set_p2i (np_d, grid % pnt_head) !np_d:beginning index


    CALL set_s2i (ns_d, grid % int_head)


    CALL clr_s2i (grid % ext_head) !set indexes to zero


    DO e=1,grid % number_of_edges
       CALL set_p2i (np_b, grid % edges(e) % pnt_head)
    ENDDO

    DO e=1,grid % number_of_edges
       CALL set_s2i (ns_b, grid % edges(e) % int_head)
    ENDDO


    IF ( ALLOCATED ( iparent )   ) DEALLOCATE ( iparent    )
    IF ( ALLOCATED ( idir )      ) DEALLOCATE ( idir       )
    IF ( ALLOCATED ( jjdir )     ) DEALLOCATE ( jjdir      )
    IF ( ALLOCATED ( iflux )     ) DEALLOCATE ( iflux      )
    IF ( ALLOCATED ( jjs )       ) DEALLOCATE ( jjs        )
    IF ( ALLOCATED ( neighs )    ) DEALLOCATE ( neighs     )
    IF ( ALLOCATED ( neighs_b )  ) DEALLOCATE ( neighs_b   )
    IF ( ALLOCATED ( neighs_low )) DEALLOCATE ( neighs_low )
    IF ( ALLOCATED ( neigh )     ) DEALLOCATE ( neigh      )

    IF ( ALLOCATED ( jjs_b ) )   DEALLOCATE ( jjs_b )
    IF ( ALLOCATED ( rr_b ) )    DEALLOCATE ( rr_b )

    IF ( ASSOCIATED ( jj_previous ) ) DEALLOCATE(jj_previous)
    IF ( ASSOCIATED ( rr_previous ) ) DEALLOCATE(rr_previous)

    IF ( ASSOCIATED ( jj ) ) THEN

       jj_previous  =>  jj
       rr_previous  =>  rr

       NULLIFY  ( jj )
       NULLIFY  ( rr )

    ENDIF

    ALLOCATE ( iparent(np_d))
    ALLOCATE ( idir(np_b)   )
    ALLOCATE ( jjdir(np_b)  )

    ALLOCATE ( iflux(ns_b)           )
    ALLOCATE ( jjs(nd_b+1, ns_b)     )
    ALLOCATE ( neighs(ns_b)          )
    ALLOCATE ( neighs_b(nd_b+1, ns_b))
    ALLOCATE ( neighs_low(ns_b)      ); neighs_low = 0

    ALLOCATE ( jj(nd_d+1, ns_d) )
    ALLOCATE ( neigh(nd_d+1, ns_d) )

    ALLOCATE ( rr(nd_d, np_d) )

    ALLOCATE ( jjs_b(nd_b+1, ns_b) )
    ALLOCATE ( rr_b(nd_b, np_b) )


    CALL p_vect_domain (grid % pnt_head)
    CALL s_vect_domain (grid % int_head)

    n = 0
    DO e=1,grid % number_of_edges
       CALL p_vect_boundary (grid % edges(e) % pnt_head, n)
       !completes a vector with all the boundary points
    ENDDO

    n = 0
    DO e=1,grid % number_of_edges
       CALL s_vect_boundary (grid % edges(e) % int_head, e, n)
    ENDDO


  END SUBROUTINE build_arrays


  !같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같
  SUBROUTINE p_vect_domain (head)

    IMPLICIT NONE

    TYPE (point), TARGET, INTENT(IN) :: head

    INTEGER :: k
    TYPE (point), POINTER :: p

    p  =>  head
    DO; p  =>  p % next; IF ( ASSOCIATED (p, head) ) EXIT
       iparent(p % index) = p % parent
       DO k=1,nd_d
          rr(k,p % index) = p % x(k)
       ENDDO
    END DO

  END SUBROUTINE p_vect_domain


  SUBROUTINE s_vect_domain (head)

    IMPLICIT NONE

    TYPE (simplex), TARGET, INTENT(IN) :: head
    TYPE (simplex), POINTER :: s
    INTEGER :: i


    s  =>  head
    DO; s  =>  s % next; IF ( ASSOCIATED (s, head) ) EXIT

       DO i=1,nd_d+1

          IF ( ASSOCIATED (s % opp(i) % pnt) ) THEN
             jj(i,s % index) = s % opp(i) % pnt % index
             !simplex points' global domain index
          ELSE
             WRITE(*,*)'WARNING: point not associated (s_vect_domain)'
             jj(i,s % index) = 0
          ENDIF

       ENDDO

       DO i=1,nd_d+1
          IF ( ASSOCIATED (s % opp(i) % spx) ) THEN
             neigh(i,s % index) = s % opp(i) % spx % index
             !adjacent simplices' global domain index
          ELSE
             neigh(i,s % index) = 0
          ENDIF

       ENDDO

    ENDDO

  END SUBROUTINE s_vect_domain

  !같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같�
  SUBROUTINE p_vect_boundary (head, n)

    IMPLICIT NONE

    TYPE(point), POINTER       :: head
    INTEGER    , INTENT(INOUT) :: n

    INTEGER :: i, k
    TYPE (point), POINTER :: p

    i = n; p  =>  head

    DO; p  =>  p % next; IF ( ASSOCIATED (p, head) ) EXIT
       i = i+1
       idir(i)   = p % cnd         
       jjdir(i)  = p % ppp % index !2D point index
       DO k=1,nd_b
          rr_b(k,i) = p % x(k)
       ENDDO
    END DO

    n = i

  END SUBROUTINE p_vect_boundary

  !같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같
  SUBROUTINE s_vect_boundary (head, e, n)

    IMPLICIT NONE

    TYPE (simplex), POINTER       :: head
    INTEGER       , INTENT(IN) :: e
    INTEGER       , INTENT(INOUT) :: n

    TYPE (simplex), POINTER :: s
    INTEGER :: i, j, k

    i = n; s  =>  head

    DO; s  =>  s % next; IF ( ASSOCIATED (s, head) ) EXIT

       i = i+1
       iflux(i) = e !1d simplex's belonging to edge

       IF ( ASSOCIATED(s % int) ) THEN
          neighs(i) = s % int % index                          
       ELSE
          neighs(i) = 0
          WRITE(*,*)'Boundary simplex without neighs() found'
       ENDIF

       DO j=1,nd_b+1
          IF ( ASSOCIATED (s % opp(j) % pnt) ) THEN
             IF ( idir(s % opp(j) % pnt % index) < 0 ) THEN ! AG: Boundary node
                DO k = 1, nd_b+1
                   IF ( k .NE. j) THEN  
                      s % opp(k) % spx % index = 0 ! AG: 
                   ENDIF
                ENDDO
             ENDIF !AG
          ENDIF
       ENDDO

       DO j=1,nd_b+1
          IF ( ASSOCIATED (s % opp(j) % pnt) ) THEN
             jjs_b(j,i) = s % opp(j) % pnt % index
             jjs(j,i)   = s % opp(j) % pnt % ppp % index

             neighs_b(j,i) = s % opp(j) % spx % index
          ELSE
             jjs_b(j,i) = 0
             jjs(j,i)   = 0
          ENDIF
       ENDDO

    ENDDO

    n = i


  END SUBROUTINE s_vect_boundary





  SUBROUTINE uns_save_conv ( name, name_length, Nnodes, Ntrias )
  !------------------------------------------------------------------------
  IMPLICIT NONE
  
  CHARACTER(LEN=64), INTENT(IN) :: name
  INTEGER,           INTENT(IN) :: name_length 
  INTEGER,           INTENT(INOUT) :: Nnodes
  INTEGER,           INTENT(INOUT) :: Ntrias
   
  INTEGER :: i, n, m
  INTEGER, DIMENSION(:), ALLOCATABLE :: jjdir_tmp
  !------------------------------------------------------------------------

  ! WARNING ADDED 10/01/2000
  ! Each boundary element has its proper and 
  ! independent nodes  -  CALL connect_edges

  Np_b = 2 * SIZE(jjs_b,2)
  ALLOCATE (jjdir_tmp(Np_b))

  n = 0
  DO m = 1, SIZE(jjs_b,2)

    i = jjs_b(1,m)    
    n = n+1
    jjdir_tmp(n) = jjdir(i)
    jjs_b(1,m) = n

    i = jjs_b(2,m)    
    n = n+1
    jjdir_tmp(n) = jjdir(i)
    jjs_b(2,m) = n

  ENDDO

  DEALLOCATE (jjdir)
  ALLOCATE (jjdir(Np_b))
  jjdir = jjdir_tmp


  OPEN (1,file='grid.'//name(1:name_length))
  CALL save_conv_mesh( 1, name, name_length, Ntrias )

  OPEN (1,file='nodes.'//name(1:name_length))
  CALL save_conv_nodes( 1, name, name_length, Nnodes )

  END SUBROUTINE uns_save_conv



  SUBROUTINE save_conv_mesh ( idf, name, name_length, Ntrias )
  !------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,            INTENT(IN) :: idf
  CHARACTER (LEN=64), INTENT(IN) :: name
  INTEGER,            INTENT(IN) :: name_length 
  INTEGER,            INTENT(INOUT) :: Ntrias

  INTEGER :: m, j_, f_
  !------------------------------------------------------------------------

  Ntrias = ns_d
  
  WRITE (idf,1000)
  WRITE (idf,1010)  name(1:name_length)
  WRITE (idf,1000); WRITE (idf,1020)
  WRITE (idf,1021)  ns_d, ns_b
  WRITE (idf,1000); WRITE (idf,1025); WRITE (idf,1000)
  WRITE (idf,1040); WRITE (idf,1041); WRITE (idf,1042)
  WRITE (idf,1043)

  DO m = 1, ns_d

    WRITE (idf,1046) m, 2
    DO j_ = 1, 3; WRITE (UNIT=idf,FMT=1047,ADVANCE='NO') jj(j_,m);    ENDDO
    WRITE(idf,*)
    DO f_ = 1, 3; WRITE (UNIT=idf,FMT=1047,ADVANCE='NO') neigh(f_,m); ENDDO
    WRITE(idf,*)

  ENDDO

  WRITE (idf,1000); WRITE (idf,1055); WRITE (idf,1000)
  WRITE (idf,1040); WRITE (idf,1065); WRITE (idf,1042) 
  WRITE (idf,1043)

  DO m = 1, ns_b

    WRITE (idf,1066) m, 1, iflux(m)
    DO j_ = 1, 2; WRITE (UNIT=idf,FMT=1047,ADVANCE='NO') jjs_b(j_,m);	 ENDDO
    WRITE(idf,*)
    DO f_ = 1, 2; WRITE (UNIT=idf,FMT=1047,ADVANCE='NO') neighs_b(f_,m); ENDDO
    WRITE(idf,*)

  ENDDO

  1000 FORMAT('###########################################################################')
  1010 FORMAT('#    NAME:      ',a15,'  					 #')
  1020 FORMAT('#  NE_D   NE_B								 #')
  1021 FORMAT(2i12)
  1025 FORMAT('#  **********  DOMAIN  **********					 #')
  1040 FORMAT('#  ++++  ELEMENTS ++++							 #')
  1041 FORMAT('#   IDX   TYPE								 #')
  1042 FORMAT('#   J_M  								 #')
  1043 FORMAT('#   MA_M 								 #')
  1046 FORMAT(2i12)
  1047 FORMAT(i12)
  1055 FORMAT('#  **********  BOUNDARY  **********					 #')
  1065 FORMAT('#   IDX   TYPE  BOUND							 #')
  1066 FORMAT(3i12)

  END SUBROUTINE save_conv_mesh



  SUBROUTINE save_conv_nodes ( idf, name, name_length, Nnodes )
  !------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,            INTENT(IN) :: idf
  CHARACTER (LEN=64), INTENT(IN) :: name
  INTEGER,            INTENT(IN) :: name_length 
  INTEGER,            INTENT(INOUT) :: Nnodes

  INTEGER :: i, m, n, p, NpB
  INTEGER, DIMENSION(:), ALLOCATABLE :: bound_p
  !------------------------------------------------------------------------

  ALLOCATE (bound_p(np_b))
  DO m = 1, ns_b 
    bound_p(jjs_b(:,m)) = iflux(m)		   
  ENDDO

  Nnodes = np_d

  WRITE (idf,1000); WRITE (idf,1010) name(1:name_length)
  WRITE (idf,1000); WRITE (idf,1020)
  WRITE (idf,1021)  nd_d, np_d, np_b 
  WRITE (idf,1000); WRITE (idf,1025)
  WRITE (idf,1000); WRITE (idf,1048)
  WRITE (idf,1050); WRITE (idf,1051)
  
  DO i = 1, np_d
    WRITE (idf,1052) i; WRITE (idf,1053) rr(:,i)
  ENDDO

  WRITE (idf,1000); WRITE (idf,1055)
  WRITE (idf,1000); WRITE (idf,1048)	 
  WRITE (idf,1075)
  
  DO i = 1, np_b
    WRITE (idf,1076) i, jjdir(i), bound_p(i)
  ENDDO    

  1000 FORMAT('###########################################################################')
  1010 FORMAT('#    NAME:      ',a15,'  					 #')
  1020 FORMAT('#    ND   NP_D	NP_B							 #')
  1021 FORMAT(5i12)
  1025 FORMAT('#  **********  DOMAIN  **********					 #')
  1048 FORMAT('#  ++++   NODES   ++++							 #')
  1050 FORMAT('#   IDX  								 #')
  1051 FORMAT('#    RR  								 #')
  1052 FORMAT(i12)
  1053 FORMAT(3(1x,e22.16))
  1055 FORMAT('#  **********  BOUNDARY  **********					 #')
  1075 FORMAT('#   IDX  JD_JB  BOUND							 #')
  1076 FORMAT(3i12)


  ! Marco Fossati. May 2008.
  ! Adding code to save curve coordinates for the
  ! nodes lying on the boundary lines. This information
  ! is necessary to grid adaption.
  OPEN  (idf, file = 'bnknots.'//name(1:name_length), form = 'formatted')
  WRITE (idf,'(6x,a6)') 'N_CURV'
  WRITE (idf,'(i12)') b_grid % number_of_edges

  DO m = 1, b_grid % number_of_edges
    
    NpB = COUNT(bound_p == m)
    p = size(curv(m)%s)
    
    WRITE (idf,'(2a12)') '	   DIM','      POINTS'
    WRITE (idf,'(2i12)') 2, NpB/2 + 1	     
    WRITE (idf,'( a12)') '	   IDC'

    IF (b_grid % edges(m) % curveType == 'lin')  WRITE (idf,'( i12)') 0
    IF (b_grid % edges(m) % curveType == 'ell')  WRITE (idf,'( i12)') 2
    IF (b_grid % edges(m) % curveType == 'cir')  WRITE (idf,'( i12)') 2
    IF (b_grid % edges(m) % curveType == 'dat')  WRITE (idf,'( i12)') 0
    IF (b_grid % edges(m) % curveType == 'bms')  WRITE (idf,'( i12)') 0

    WRITE (idf,'(2a12)') '	     X','	    Y'

    DO n = 1, SIZE(bound_p) - 1
      IF (bound_p(n) == m) THEN
  	IF (jjdir(n) /= jjdir(n+1)  .OR.  bound_p(n) /= bound_p(n+1)) THEN
  	  WRITE(idf,'(2(1x,e22.16))') rr(:,jjdir(n))
  	ENDIF
      ENDIF
    ENDDO
    
    IF (bound_p(SIZE(bound_p)) == m)  WRITE(idf,'(2(1x,e22.16))') rr(:,jjdir(SIZE(bound_p)))

  ENDDO
  
  ! Marco Fossati. May 2008
  ! Allocation only to count the number of boundary nodes
  ! to be used in writing file abscissa
  ALLOCATE (b_bound_p(SIZE(bound_p)))
    
  END SUBROUTINE save_conv_nodes



  SUBROUTINE hyb_save_conv ( name, name_length, Nnodes, Ntrias, Nquads )
  !------------------------------------------------------------------------
  IMPLICIT NONE

  CHARACTER (LEN=64), INTENT(IN) :: name
  INTEGER,            INTENT(IN) :: name_length
  INTEGER,            INTENT(INOUT) :: Nnodes, Ntrias, Nquads
  !------------------------------------------------------------------------

  OPEN  (1,file='grid.'//name(1:name_length))
  CALL hyb_save_conv_mesh ( 1, name, name_length, Ntrias, Nquads )
    
  OPEN  (1,file='nodes.'//name(1:name_length))
  CALL hyb_save_conv_nodes ( 1, name, name_length, Nnodes )

  END SUBROUTINE hyb_save_conv



  SUBROUTINE hyb_save_conv_mesh ( idf, name, name_length, Ntrias, Nquads )
  !------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,            INTENT(IN) :: idf
  CHARACTER (LEN=64), INTENT(IN) :: name
  INTEGER,            INTENT(IN) :: name_length 
  INTEGER,            INTENT(INOUT) :: Ntrias, Nquads

  INTEGER :: m, j_, f_, countT, countQ
  !------------------------------------------------------------------------

  WRITE (idf,1000); WRITE (idf,1010) name(1:name_length)
  WRITE (idf,1000); WRITE (idf,1020)
  WRITE (idf,1021)  g_ns_d, g_np_b
  WRITE (idf,1000); WRITE (idf,1025)
  WRITE (idf,1000); WRITE (idf,1040);
  WRITE (idf,1041); WRITE (idf,1042);
  WRITE (idf,1043)
  
  countT = 0
  countQ = 0

  DO m = 1, g_ns_d

    IF ( g_m_t(m) .eq. 2 ) THEN
      countT = countT + 1
    ELSE IF ( g_m_t(m) .eq. 3 ) THEN
      countQ = countQ + 1
    END IF

    WRITE (idf,1046) m, g_m_t(m)
    DO j_ = 1, g_m_t(m) + 1;  WRITE (UNIT=idf,FMT=1047,ADVANCE='NO') g_jj(j_,m);    ENDDO
    WRITE(idf,*)
    DO f_ = 1, g_m_t(m) + 1;  WRITE (UNIT=idf,FMT=1047,ADVANCE='NO') g_neigh(f_,m); ENDDO
    WRITE(idf,*)

  ENDDO

  Ntrias = countT
  Nquads = countQ

  WRITE (idf,1000); WRITE (idf,1055); WRITE (idf,1000);
  WRITE (idf,1040); WRITE (idf,1065); WRITE (idf,1042);
  WRITE (idf,1043)

  DO m = 1, g_np_b

    WRITE (idf,1066) m, 1, b_bound_p(2*m)
    DO j_ = 1, 2;  WRITE (UNIT=idf,FMT=1047,ADVANCE='NO') b_jjs_b(j_,m);    ENDDO
    WRITE(idf,*)
    DO f_ = 1, 2;  WRITE (UNIT=idf,FMT=1047,ADVANCE='NO') b_neighs_b(f_,m); ENDDO
    WRITE(idf,*)

  ENDDO

  1000 FORMAT('###########################################################################')
  1010 FORMAT('#    NAME:      ',a15,'  					 #')
  1020 FORMAT('#  NE_D   NE_B								 #')
  1021 FORMAT(2i12)
  1025 FORMAT('#  **********  DOMAIN  **********					 #')
  1040 FORMAT('#  ++++  ELEMENTS ++++							 #')
  1041 FORMAT('#   IDX   TYPE								 #')
  1042 FORMAT('#   J_M  								 #')
  1043 FORMAT('#   MA_M 								 #')
  1046 FORMAT(2i12)
  1047 FORMAT(i12)
  1055 FORMAT('#  **********  BOUNDARY  **********					 #')
  1065 FORMAT('#   IDX   TYPE  BOUND							 #')
  1066 FORMAT(3i12)

  END SUBROUTINE hyb_save_conv_mesh



  SUBROUTINE hyb_save_conv_nodes ( idf, name, name_length, Nnodes )
  !------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,            INTENT(IN) :: idf
  CHARACTER (LEN=64), INTENT(IN) :: name
  INTEGER,            INTENT(IN) :: name_length
  INTEGER,            INTENT(INOUT) :: Nnodes

  INTEGER :: i, n, m, NpB 
  REAL*8, DIMENSION(2) :: cc_prev, cc_curr
  REAL*8 :: scurr, length
  LOGICAL :: newc
  !------------------------------------------------------------------------

  Nnodes = g_np_d

  WRITE (idf,1000); WRITE (idf,1010) name(1:name_length)
  WRITE (idf,1000); WRITE (idf,1020)
  WRITE (idf,1021)  nd_d, g_np_d, 2*g_np_b 
  WRITE (idf,1000); WRITE (idf,1025)
  WRITE (idf,1000); WRITE (idf,1048)
  WRITE (idf,1050); WRITE (idf,1051)
  
  DO i = 1, g_np_d
    WRITE (idf,1052) i;  WRITE (idf,1053) g_rr(:,i)
  ENDDO

  WRITE (idf,1000); WRITE (idf,1055)
  WRITE (idf,1000); WRITE (idf,1048)	 
  WRITE (idf,1075)
  
  DO i = 1, 2*g_np_b
    WRITE (idf,1076) i, b_jjdir(i), b_bound_p(i)
  ENDDO

  1000 FORMAT('###########################################################################')
  1010 FORMAT('#    NAME:      ',a15,'  					 #')
  1020 FORMAT('#    ND   NP_D	NP_B							 #')
  1021 FORMAT(5i12)
  1025 FORMAT('#  **********  DOMAIN  **********					 #')
  1048 FORMAT('#  ++++   NODES   ++++							 #')
  1050 FORMAT('#   IDX  								 #')
  1051 FORMAT('#    RR  								 #')
  1052 FORMAT(i12)
  1053 FORMAT(3(1x,e22.16))
  1055 FORMAT('#  **********  BOUNDARY  **********					 #')
  1075 FORMAT('#   IDX  JD_JB  BOUND							 #')
  1076 FORMAT(3i12)


  ! Marco Fossati. May 2008.
  ! Adding code to save curve coordinates for the
  ! nodes lying on the boundary lines. This information
  ! is necessary to grid adaption.
  OPEN  (idf, file = 'bnknots.'//name(1:name_length), form = 'formatted')
  OPEN  (123, file = 'bnabsci.'//name(1:name_length), form = 'formatted')
  WRITE (idf,'(6x,a6)') 'N_CURV'
  WRITE (idf,'(i12)') b_grid % number_of_edges

  DO m = 1, b_grid % number_of_edges
    
    NpB = COUNT(b_bound_p == m)
    
    WRITE (idf,'(2a12)') '	   DIM','      POINTS'
    WRITE (idf,'(2i12)') 2, NpB/2 + 1
    WRITE (idf,'( a12)') '	   IDC'

    IF (b_grid % edges(m) % curveType == 'lin')  WRITE (idf,'( i12)') 0
    IF (b_grid % edges(m) % curveType == 'ell')  WRITE (idf,'( i12)') 2
    IF (b_grid % edges(m) % curveType == 'cir')  WRITE (idf,'( i12)') 2
    IF (b_grid % edges(m) % curveType == 'dat')  WRITE (idf,'( i12)') 0
    IF (b_grid % edges(m) % curveType == 'bms')  WRITE (idf,'( i12)') 0

    WRITE (idf,'(3a12)') '	     X','	    Y','	    S'

    newc = .true.
    length = 0.d0
    DO n = 1, SIZE(b_bound_p) - 1
      IF (b_bound_p(n) == m) THEN
  	IF (b_jjdir(n) /= b_jjdir(n+1)  .OR.  b_bound_p(n) /= b_bound_p(n+1)) THEN
          if ( newc ) then
            newc = .false.
            cc_prev = g_rr(:,b_jjdir(n))
          end if
          cc_curr = g_rr(:,b_jjdir(n))
          length = length + sqrt( sum(( cc_curr - cc_prev )**2) )
          cc_prev = g_rr(:,b_jjdir(n))
  	ENDIF  
      ENDIF
    ENDDO
    
    newc = .true.
    scurr = 0.d0
    DO n = 1, SIZE(b_bound_p) - 1
      IF (b_bound_p(n) == m) THEN
  	IF (b_jjdir(n) /= b_jjdir(n+1)  .OR.  b_bound_p(n) /= b_bound_p(n+1)) THEN
          if ( newc ) then
            newc = .false.
            cc_prev = g_rr(:,b_jjdir(n))
          end if
          cc_curr = g_rr(:,b_jjdir(n))
          scurr = scurr + sqrt( sum(( cc_curr - cc_prev )**2) )
          WRITE(123,*)g_rr(:,b_jjdir(n)), scurr / length
  	  WRITE(idf,'(3f20.10)') g_rr(:,b_jjdir(n)), scurr / length
          cc_prev = g_rr(:,b_jjdir(n))
  	ENDIF  
      ENDIF
    ENDDO
    
    IF (b_bound_p(SIZE(b_bound_p)) == m)  WRITE(idf,'(2f20.10)') g_rr(:,b_jjdir(SIZE(b_bound_p)))
    
  ENDDO

  END SUBROUTINE hyb_save_conv_nodes



  SUBROUTINE smooth(grid, kmax, tol)
  !------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER :: i, j, k
  INTEGER :: kmax
  REAL(KIND=8) :: tol

  REAL(KIND=8) :: e
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dr
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: rn, s

  TYPE(point), POINTER :: p

  TYPE(domain), TARGET, INTENT(IN) :: grid
  !------------------------------------------------------------------------

    NULLIFY(array)
    array => grid

    IF (.NOT. ASSOCIATED(array)) THEN
      WRITE (*,*); WRITE (*,*)  'ERROR. Array not associated.'
      WRITE (*,*); STOP
    ENDIF

    ALLOCATE (dr(nd_d, np_d))
    ALLOCATE (rn(np_d))
    ALLOCATE (s(nd_d))

    rn = 0.d0

    DO i = 1, ns_d
      DO j = 1, nd_d+1
        rn(jj(j,i)) = rn(jj(j,i))+2.
      ENDDO
    ENDDO

    DO i = 1, np_d;  rn(i) = 1./rn(i);  ENDDO

    DO k = 1, kmax
       dr = 0.

       DO i = 1, ns_d
         s = 0.d0
         DO j = 1, nd_d+1;  s(:) = s(:) + rr(:,jj(j,i));                           ENDDO
         DO j = 1, nd_d+1;  dr(:,jj(j,i)) = dr(:,jj(j,i)) + s(:) - rr(:,jj(j,i));  ENDDO
       ENDDO

       DO i = 1, np_d;  dr(:,i) = rn(i)*dr(:,i) - rr(:,i);  ENDDO
       DO i = 1, np_b;  dr(:,jjdir(i)) = 0.d0;              ENDDO
       DO i = 1, np_d;  rr(:,i) = rr(:,i) + dr(:,i);        ENDDO

       e = 0.

       DO i = 1, np_d
         DO j = 1, nd_d
            e = e + dr(j,i)**2
         ENDDO
       ENDDO

       !WRITE (*,'(3x,i5,e12.5)') k, e
       IF (e < tol) EXIT

    ENDDO

    DEALLOCATE (dr, rn, s)

    p  =>  array % pnt_head
    DO; p  =>  p % next; IF ( ASSOCIATED (p, array % pnt_head) ) EXIT
      p % x = rr(:,p % index)
    ENDDO

  END SUBROUTINE smooth

  END MODULE arrays

 
