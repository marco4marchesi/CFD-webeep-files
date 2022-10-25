    MODULE qgraph

    USE grid_types
    USE schemi

    contains

    SUBROUTINE q_g_plot ( idf )

    IMPLICIT NONE
  
    INTEGER, PARAMETER :: background = -1, foreground = 0,   &
         yellow = 1, blue = 2, green = 3, red = 4,   &
         dark_blue = 5, orange = 6, pink = 7,   &
         light_pink = 8, cyan = 9, brown = 9   !  PLOTMTV COLORS

    INTEGER, PARAMETER :: none = 0, dot = 1, cross = 2, capital_x = 3,   &
         white_square = 4, black_square = 5,   &
         white_diamond = 6, black_diamond = 7,   &
         white_triangle = 8, black_triangle = 9,   &
         white_inv_triangle = 10, black_inv_triangle = 11,   &
         white_circle = 12, black_circle = 13  !   PLOTMTV MARKERS

    INTEGER, PARAMETER :: dom_color = 0, bou_color = yellow
    
    INTEGER     :: marker, i, j, t, col, nlayer
    INTEGER      , INTENT(IN)   :: idf
    TYPE (point) , POINTER      :: p, head, q, r, s, h1, h2
    logical :: term

    WRITE (idf,*) '$ DATA = CURVE2D'
    WRITE (idf,*) '% equalscale = TRUE'
    WRITE (idf,*) '% boundary = TRUE'
    WRITE (idf,*) '% toplabel='' '' xlabel='' '' ylabel='' ''   '
    WRITE (idf,*) '#%pointID = TRUE'

    marker = none
    nlayer = b_grid % number_of_layers

    DO i = 1, nlayer
      DO j = 1, SIZE(layer(i) % pnt_list) 

      h1 => layer(i) % pnt_list(j)
      h2 => layer(i+1) % pnt_list(j)
       
      IF ( i .gt. layer(i) % max_lay(j) )  CYCLE
      
      p => h1 
      q => p % next

      r => h2
      s => r % next
      
      col = dom_color
     
     DO
      
      term = .true.
      SELECT CASE ( p % what )

      CASE (1) 
         CALL line_plot(r, p)
         CALL line_plot(r, s)
         r => s
         s => s % next
      
      CASE (0)   ! nessun elemento(ultimo fronte o interferenza
      !   r => s; s => s % next
         !IF(ASSOCIATED(r % old1))THEN
         ! write(*,*)'*-->', r % x
         ! r => s; s => s % next 
         !ENDIF
            IF(ASSOCIATED(r % old1)) THEN
              IF(ASSOCIATED(r % old1, r % prev % old))THEN
               p => p % prev; q => q % prev
               ! caso 10 invisibile          
              ELSE
               IF(ASSOCIATED(r % old, q % next))THEN
                ! caso 6 doppio 
                 p => p % next % next; q => p % next
                 CYCLE
               ELSEIF(ASSOCIATED(r % old1 % next, q))THEN
                ! caso 6 singolo 
                 p => p % next; q => p % next
                 CYCLE
               ENDIF
              ENDIF
            ENDIF
         r => s; s => s % next
  
      CASE (10) 
         CALL line_plot(r, p)
         CALL line_plot(r, s)
         CALL line_plot(s, p)
         CALL line_plot(s, q)
         CALL line_plot(s, s % next)
         r => s % next; s => r % next
      
      CASE (6) 
         CALL line_plot(r, p)
         
      CASE (2) 
         CALL line_plot(r, p)
         CALL line_plot(r, s)
         r => s; s => s % next
      
      CASE (3) 
         CALL line_plot(r, p)
         CALL line_plot(q, s)
         CALL line_plot(r, s)
         r => s; s => s % next
      
      CASE (12) 
         CALL line_plot(p, s)
         r => s; s => s % next
      
      CASE (13) 
         CALL line_plot(r, p)
         CALL line_plot(r, q)
         r => s ; s => s % next
      
      CASE (8) 
         CALL line_plot(r, p)
         CALL line_plot(r, s)
         r => s
         s => s % next
      
      CASE (9) 
         CALL line_plot(r, p)
         CALL line_plot(r, s)
         CALL line_plot(s, p)
         CALL line_plot(s, s % next)
         r => s % next; s => r % next
         CALL line_plot(r, p)
         CALL line_plot(r, s)
         r => s; s => s % next         
      
      CASE (14)                  ! estensione di scia
        IF((i==1).AND.(ASSOCIATED(p % wake)))THEN               
            q => p % wake
            r => p
            DO
              CALL line_plot(r, q)
              r => q
              IF(.NOT.ASSOCIATED(r % next)) EXIT
              q => r % next
            ENDDO
            s => p % ppp1
            r => p % ppp
            DO
             CALL line_plot(r, r % old)
             CALL line_plot(r, r % next)
             r => r % next
             IF(ASSOCIATED(r,s))EXIT
            ENDDO
            CALL line_plot(r, p)
            CALL line_plot(r, r % next)
            r => r % next; s => r % next
        ENDIF
      
      CASE (15)
         CALL line_plot(r, p)
         CALL line_plot(q, s)
         CALL line_plot(r, s)
         r => s; s => s % next
      
      CASE (16:18) 
         CALL line_plot(r, p)
         CALL line_plot(r, s)
         CALL line_plot(s, p)
         CALL line_plot(s, q)
         CALL line_plot(s, s % next)
         CALL line_plot(q, s % next)
         r => s % next; s => r % next
     
      CASE (19) 
         CALL line_plot(r, p)
         CALL line_plot(r, s)
         CALL line_plot(s, p)
         
         r => s; s => r % next
         CALL line_plot(r, s)
         CALL line_plot(s, p)
         
         r => s; s => r % next
         CALL line_plot(r, s)
         CALL line_plot(s, p)
         
         r => s; s => r % next
         CALL line_plot(r, s)
         CALL line_plot(s, p)

     END SELECT
     
       p => p % next
       q => p % next 
     
     IF ( ASSOCIATED(p,h1) .and. term ) EXIT

    ENDDO
   ENDDO
  ENDDO
  
!  BOUNDARIES....
  t = SIZE(layer(1) % pnt_list) 
  DO i = 1,t
    
    head => layer(1) % pnt_list(i)
    p => head
    DO 
    
      marker = none 
      col = 1
      q => p % next
      
      IF (q % index /= head % index)  CALL line_plot(p, p % next)
      p => p % next
      
      IF (ASSOCIATED(q,head)) EXIT
    
    ENDDO
  ENDDO


CONTAINS

  SUBROUTINE line_plot(p,q)
  
   IMPLICIT NONE 
  
   TYPE (point) ,POINTER :: p, q
    
   WRITE (idf,*)
   WRITE (idf,*) '% linecolor = ',col,' markertype = ',marker 
   WRITE (idf,'(2e15.7,i8)') p % x(1), p % x(2), p % index
   WRITE (idf,'(2e15.7,i8)') q % x(1), q % x(2), q % index 
   
  END SUBROUTINE line_plot
                   
END SUBROUTINE q_g_plot
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE q_g_plot_jump(idf, jdf)

    IMPLICIT NONE
  
    INTEGER, PARAMETER :: background = -1, foreground = 0,          &
                          yellow = 1, blue = 2, green = 3, red = 4, &
                          dark_blue = 5, orange = 6, pink = 7,      &
                          light_pink = 8, cyan = 9, brown = 9

    INTEGER, PARAMETER :: none = 0, dot = 1, cross = 2, capital_x = 3,       &
                          white_square = 4, black_square = 5,	             &
                          white_diamond = 6, black_diamond = 7,              &
                          white_triangle = 8, black_triangle = 9,            &
                          white_inv_triangle = 10, black_inv_triangle = 11,  &
                          white_circle = 12, black_circle = 13

    INTEGER, PARAMETER :: dom_color = 0
    INTEGER, PARAMETER :: bou_color = yellow
    INTEGER :: marker 
  
    INTEGER :: i,d,t,col, nP = 0
    INTEGER,INTENT(IN) :: idf, jdf
    TYPE(point), POINTER :: p, head, q

    
    WRITE (idf,*) '% equalscale = TRUE'
    WRITE (idf,*) '% toplabel='' '' xlabel='' '' ylabel='' ''	'
    WRITE (idf,*) '#%xmin=   xmax= '
    WRITE (idf,*) '#%xmin=   xmax= '
    WRITE (idf,*) '#%pointID = TRUE'

    marker = none

    !  DOMAIN....
    d = done_layers
    col = green
    
    t = SIZE(big_box)
    DO i = 1, t
   
      IF (big_box(i) % structured == 0) EXIT
        					     
      head => big_box(i) % virt_new_head	     
      p => head 				     

      DO					     
             					     
        !IF (ASSOCIATED(p % jump_f)) THEN
        !  col = green
        !  CALL line_plot(p, p % jump_f)
        !ENDIF					     
        					     
        IF (ASSOCIATED(p % jump_b)) THEN
          col = red
          CALL line_plot(p, p % jump_b)
	  nP = nP + 1			     
	  WRITE(jdf,*) p % index, p % jump_b % index
        ENDIF					     
        					     
        p => p % jump_f
        IF (ASSOCIATED(p, head)) EXIT 

      ENDDO					     

    ENDDO

    WRITE (jdf,*) 'Number of node-pairs at interface'
    WRITE (jdf,*) nP

    !  BOUNDARIES....
    col = yellow
    t = SIZE(layer(1) % pnt_list) 
    
    DO i = 1,t
    
      head => layer(1) % pnt_list(i)			   
      p => head 					   
      							   
      DO						   
      							   
        CALL line_plot(p, p % next)			   
        						   
        IF (ASSOCIATED(p % wake))THEN			   
        						   
          CALL line_plot(p, p % wake)			   
        						   
          IF (ASSOCIATED(p % wake % next)) THEN  	      
            q => p % wake % next 		             	   
          ELSE						      
            WRITE(*,*)'stop!  (wake)'
            WRITE(*,*)'scia costituita di un solo elemento?'  
          ENDIF  					      
        						   
          DO				     
         				    		   
            CALL line_plot(p, q) 	    		   
         				    		   
            IF (.NOT.ASSOCIATED(q%next)) THEN		   
              EXIT			    		   
            ELSE 			    		   
              q => q % next		    		   
            ENDIF			    		   
         				    		   
          ENDDO  			    		   

        ENDIF	 					   
        						   
        p => p % next
        IF (ASSOCIATED(p,head)) EXIT
       			   
      ENDDO						   

    ENDDO


    CONTAINS
  
  
    SUBROUTINE line_plot(p,q)
  
    IMPLICIT NONE  
    TYPE (point) ,POINTER :: p, q
     
    WRITE (idf,*)
    WRITE (idf,*) '% linecolor = ',col,' markertype = ',marker
    WRITE (idf,'(2e15.7,i8,1x,i2)') p % x(1), p % x(2), p % index, p % cnd
    WRITE (idf,'(2e15.7,i8,1x,i2)') q % x(1), q % x(2), q % index, q % cnd
    
    END SUBROUTINE line_plot


END SUBROUTINE q_g_plot_jump


END MODULE qgraph
