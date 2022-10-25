MODULE convert

   USE init_data

   IMPLICIT NONE

   TYPE p_pnter
      TYPE (point), POINTER :: pnter
   END TYPE p_pnter

   TYPE s_pnter
      TYPE (simplex), POINTER :: pnter
   END TYPE s_pnter


CONTAINS 


   !   ROUTINES TO SAVE, LOAD, AND COPY A GRID


   SUBROUTINE g_save (idf, ascii, grid)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: idf
      LOGICAL, INTENT(IN) :: ascii
      TYPE (domain), INTENT(IN) :: grid

      INTEGER :: e, m, n
      INTEGER :: np_b_dom, np_b_axp, ns_b_int, ns_b_ext, ns_b_axs
      INTEGER :: np_d_dom, np_d_axp, ns_d_int, ns_d_ext, ns_d_axs


      !  ###  POINT/SIMPLEX NUMBERS, LISTS, INDICES  ###
 
!write(*,*)'*'
      IF ( ascii ) THEN
         WRITE (idf,'(2(5x,a5))') 'VERTS','EDGES'
         WRITE (idf,'(2i10)') grid % number_of_verts, grid % number_of_edges
      ELSE
         WRITE (idf) grid % number_of_verts, grid % number_of_edges
      ENDIF

      !  Boundary point lists
!write(*,*)'*'
      n = 0

      IF ( ascii ) WRITE (idf,'(2(2x,a8))') 'NP_B_DOM','NP_B_AXP'
      DO e=1,grid % number_of_edges
	 m = n
	 CALL set_p2i (n, grid % edges(e) % pnt_head)
	 np_b_dom = n-m
!write(*,*)'*'
	 m = n
	 CALL set_p2i (n, grid % edges(e) % axp_head)
	 np_b_axp = n-m
!write(*,*)'*'
	 IF ( ascii ) THEN
	    WRITE (idf,'(2i10)') np_b_dom, np_b_axp
	 ELSE
	    WRITE (idf) np_b_dom, np_b_axp
	 ENDIF
      ENDDO

      !WRITE (*,*) ' b pnt number = ',n

      !  Boundary simplex lists
!write(*,*)'*'
      n = 0

      IF ( ascii ) WRITE (idf,'(3(2x,a8))') 'NS_B_INT','NS_B_EXT','NS_B_AXS'
      DO e=1,grid % number_of_edges
	 m = n
	 CALL set_s2i (n, grid % edges(e) % int_head)
	 ns_b_int = n-m

	 m = n
	 CALL set_s2i (n, grid % edges(e) % ext_head)
	 ns_b_ext = n-m

	 m = n
	 CALL set_s2i (n, grid % edges(e) % axs_head)
	 ns_b_axs = n-m

	 IF ( ascii ) THEN
	    WRITE (idf,'(3i10)') ns_b_int, ns_b_ext, ns_b_axs
	 ELSE
	    WRITE (idf) ns_b_int, ns_b_ext, ns_b_axs
	 ENDIF
      ENDDO

!      WRITE (*,*) ' b spx number = ',n

      !  Domain point lists

      n = 0

      m = n
      CALL set_p2i (n, grid % pnt_head)
      np_d_dom = n-m

      m = n
      CALL set_p2i (n, grid % axp_head)
      np_d_axp = n-m

      IF ( ascii ) THEN
         WRITE (idf,'(2(2x,a8))') 'NP_D_DOM','NP_D_AXP'
         WRITE (idf,'(2i10)') np_d_dom, np_d_axp
      ELSE
         WRITE (idf) np_d_dom, np_d_axp
      ENDIF

!      WRITE (*,*) ' d pnt number = ',n

      !  Domain simplex lists

      n = 0

      m = n
      CALL set_s2i (n, grid % int_head)
      ns_d_int = n-m

      m = n
      CALL set_s2i (n, grid % ext_head)
      ns_d_ext = n-m

      m = n
      CALL set_s2i (n, grid % axs_head)
      ns_d_axs = n-m

      IF ( ascii ) THEN
         WRITE (idf,'(3(2x,a8))') 'NS_D_INT','NS_D_EXT','NS_D_AXS'
         WRITE (idf,'(3i10)') ns_d_int, ns_d_ext, ns_d_axs
      ELSE
         WRITE (idf) ns_d_int, ns_d_ext, ns_d_axs
      ENDIF

!      WRITE (*,*) ' d spx number = ',n


      !  ###  SAVE DATA  ###


      !  Save boundary point data

      DO e=1,grid % number_of_edges
	     CALL p_save (idf, ascii, grid % edges(e) % pnt_head)
	     CALL p_save (idf, ascii, grid % edges(e) % axp_head)
      ENDDO

      !  Save boundary simplex data

      DO e=1,grid % number_of_edges
	     CALL s_save (idf, ascii, nd_b, grid % edges(e) % int_head)
	     CALL s_save (idf, ascii, nd_b, grid % edges(e) % ext_head)
	     CALL s_save (idf, ascii, nd_b, grid % edges(e) % axs_head)
      ENDDO

      !  Save domain point data

      CALL p_save (idf, ascii, grid % pnt_head)
      CALL p_save (idf, ascii, grid % axp_head)

      !  Save domain simplex data

      CALL s_save (idf, ascii, nd_d, grid % int_head)
      CALL s_save (idf, ascii, nd_d, grid % ext_head)
      CALL s_save (idf, ascii, nd_d, grid % axs_head)
write(*,*)'grid saved'

   END SUBROUTINE g_save
!같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같�

   SUBROUTINE g_load (idf, ascii, grid)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: idf
      LOGICAL, INTENT(IN) :: ascii
      TYPE (domain), INTENT(OUT) :: grid

      INTEGER :: e, m, n, np_b, np_d, ns_b, ns_d
      INTEGER :: np_d_dom, np_d_axp, ns_d_int, ns_d_ext, ns_d_axs

      INTEGER, DIMENSION(:), ALLOCATABLE :: np_b_dom, np_b_axp
      INTEGER, DIMENSION(:), ALLOCATABLE :: ns_b_int, ns_b_ext, ns_b_axs

      TYPE (p_pnter), DIMENSION(:), ALLOCATABLE :: i2p_b, i2p_d
      TYPE (s_pnter), DIMENSION(:), ALLOCATABLE :: i2s_b, i2s_d


      !  ###  POINT/SIMPLEX NUMBERS, LISTS, INDICES  ###


      IF ( ascii ) THEN
         READ  (idf,*)
         READ  (idf,*) grid % number_of_verts, grid % number_of_edges
      ELSE
         READ  (idf) grid % number_of_verts, grid % number_of_edges
      ENDIF
      !  modified by trive 24/2/97
      IF(ALLOCATED(np_b_dom)) THEN
         DEALLOCATE(np_b_dom, np_b_axp,ns_b_int, ns_b_ext, ns_b_axs) 
      ENDIF
      !  end

!by DD begin   
      ALLOCATE(grid%edges(grid%number_of_edges))      
!by DD end      
      
      ALLOCATE ( np_b_dom(grid % number_of_edges), np_b_axp(grid % number_of_edges) )
      ALLOCATE ( ns_b_int(grid % number_of_edges), ns_b_ext(grid % number_of_edges),   &
                 ns_b_axs(grid % number_of_edges) )

      !  Boundary point lists

      np_b = 0

      IF ( ascii ) READ  (idf,*)
      DO e=1,grid % number_of_edges
         IF ( ascii ) THEN
	         READ  (idf,*) np_b_dom(e), np_b_axp(e)
         ELSE
	         READ  (idf) np_b_dom(e), np_b_axp(e)
         ENDIF
	       np_b = np_b + np_b_dom(e) + np_b_axp(e)               !
        CALL p_list  (np_b_dom(e), grid % edges(e) % pnt_head)!
	       CALL p_list  (np_b_axp(e), grid % edges(e) % axp_head)!
      ENDDO
      !  Boundary simplex lists

      ns_b = 0

      IF ( ascii ) READ  (idf,*)
      DO e=1,grid % number_of_edges
         IF ( ascii ) THEN
	    READ  (idf,*) ns_b_int(e), ns_b_ext(e), ns_b_axs(e)
         ELSE
	    READ  (idf) ns_b_int(e), ns_b_ext(e), ns_b_axs(e)
         ENDIF
	 ns_b = ns_b + ns_b_int(e) + ns_b_ext(e) + ns_b_axs(e)

	 CALL s_list  (ns_b_int(e), grid % edges(e) % int_head)
	 CALL s_list  (ns_b_ext(e), grid % edges(e) % ext_head)
	 CALL s_list  (ns_b_axs(e), grid % edges(e) % axs_head)
      ENDDO
      !  Domain point lists

      np_d = 0

      IF ( ascii ) THEN
         READ  (idf,*)
         READ  (idf,*) np_d_dom, np_d_axp
      ELSE
         READ  (idf) np_d_dom, np_d_axp
      ENDIF

      np_d = np_d + np_d_dom + np_d_axp

      CALL p_list  (np_d_dom, grid % pnt_head)
      CALL p_list  (np_d_axp, grid % axp_head)
      !  Domain simplex lists

      ns_d = 0

      IF ( ascii ) THEN
         READ  (idf,*)
         READ  (idf,*) ns_d_int, ns_d_ext, ns_d_axs
      ELSE
         READ  (idf) ns_d_int, ns_d_ext, ns_d_axs
      ENDIF

      ns_d = ns_d + ns_d_int + ns_d_ext + ns_d_axs

      CALL s_list  (ns_d_int, grid % int_head)
      CALL s_list  (ns_d_ext, grid % ext_head)
      CALL s_list  (ns_d_axs, grid % axs_head)


      !  ###  SET I2P/I2S  ###


      ALLOCATE ( i2p_b(np_b), i2s_b(ns_b) )
      ALLOCATE ( i2p_d(np_d), i2s_d(ns_d) )


      !  Set i2p, boundary

      n = 0

      DO e=1,grid % number_of_edges
	 CALL set_i2p (n, grid % edges(e) % pnt_head, i2p_b)
	 CALL set_i2p (n, grid % edges(e) % axp_head, i2p_b)
      ENDDO

!      WRITE (*,*) 'b pnt number = ',n, np_b

      !  Set i2s, boundary

      n = 0

      DO e=1,grid % number_of_edges
	 CALL set_i2s (n, grid % edges(e) % int_head, i2s_b)
	 CALL set_i2s (n, grid % edges(e) % ext_head, i2s_b)
	 CALL set_i2s (n, grid % edges(e) % axs_head, i2s_b)
      ENDDO

      WRITE (*,*) 'b spx number = ',n, ns_b

      !  Set i2p, domain

      n = 0

      CALL set_i2p (n, grid % pnt_head, i2p_d)
      CALL set_i2p (n, grid % axp_head, i2p_d)

!      WRITE (*,*) 'd pnt number = ',n, np_d

      !  Set i2s, domain

      n = 0

      CALL set_i2s (n, grid % int_head, i2s_d)
      CALL set_i2s (n, grid % ext_head, i2s_d)
      CALL set_i2s (n, grid % axs_head, i2s_d)

!      WRITE (*,*) 'd spx number = ',n, ns_d


      !  ###  LOAD DATA  ###


      !  Load boundary point data

      n = 0

      DO e=1,grid % number_of_edges
	 m = n+1
	 n = n+np_b_dom(e)
	 CALL p_load (idf, ascii, m, n, i2p_b, i2s_b, i2p_d)

	 m = n+1
	 n = n+np_b_axp(e)
	 CALL p_load (idf, ascii, m, n, i2p_b, i2s_b, i2p_d)
      ENDDO

      !  Load boundary simplex data

      n = 0

      DO e=1,grid % number_of_edges
	 m = n+1
	 n = n+ns_b_int(e)
	 CALL s_load (idf, ascii, nd_b, m, n, i2p_b, i2s_b, i2s_d)

	 m = n+1
	 n = n+ns_b_ext(e)
	 CALL s_load (idf, ascii, nd_b, m, n, i2p_b, i2s_b, i2s_d)

	 m = n+1
	 n = n+ns_b_axs(e)
	 CALL s_load (idf, ascii, nd_b, m, n, i2p_b, i2s_b, i2s_d)
      ENDDO

      !  Load domain point data

      n = 0

      m = n+1
      n = n+np_d_dom
      CALL p_load (idf, ascii, m, n, i2p_d, i2s_d, i2p_d)

      m = n+1
      n = n+np_d_axp
      CALL p_load (idf, ascii, m, n, i2p_d, i2s_d, i2p_d)

      !  Load domain simplex data

      n = 0

      m = n+1
      n = n+ns_d_int
      CALL s_load (idf, ascii, nd_d, m, n, i2p_d, i2s_d, i2s_d)

      m = n+1
      n = n+ns_d_ext
      CALL s_load (idf, ascii, nd_d, m, n, i2p_d, i2s_d, i2s_d)

      m = n+1
      n = n+ns_d_axs
      CALL s_load (idf, ascii, nd_d, m, n, i2p_d, i2s_d, i2s_d)


      DEALLOCATE ( i2p_b, i2s_b, i2p_d, i2s_d )

   END SUBROUTINE g_load
!같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같

   SUBROUTINE g_copy (g_old, g_new)

      IMPLICIT NONE

      TYPE (domain), INTENT(IN) :: g_old
      TYPE (domain), INTENT(INOUT) :: g_new

      INTEGER :: ne, e, m, n, np_b, np_d, ns_b, ns_d
      TYPE (p_pnter), DIMENSION(:), ALLOCATABLE :: i2p_b, i2p_d
      TYPE (s_pnter), DIMENSION(:), ALLOCATABLE :: i2s_b, i2s_d


      !  ###  SET P2I AND ALLOCATE THE LISTS  ###


      ne = g_old % number_of_edges

      !  Boundary point lists

      np_b = 0

      DO e=1,ne
	 
   m = np_b
	 CALL set_p2i (np_b,   g_old % edges(e) % pnt_head)!set index to p (old)
	
   CALL p_list  (np_b-m, g_new % edges(e) % pnt_head)!allocate a new list
	
   m = np_b
	 CALL set_p2i (np_b,   g_old % edges(e) % axp_head)
	 CALL p_list  (np_b-m, g_new % edges(e) % axp_head)
      ENDDO

      !  Boundary simplex lists

      ns_b = 0

      DO e=1,ne
         m = ns_b
	 CALL set_s2i (ns_b,   g_old % edges(e) % int_head)
	 CALL s_list  (ns_b-m, g_new % edges(e) % int_head)
         m = ns_b
	 CALL set_s2i (ns_b,   g_old % edges(e) % ext_head)
	 CALL s_list  (ns_b-m, g_new % edges(e) % ext_head)
         m = ns_b
	 CALL set_s2i (ns_b,   g_old % edges(e) % axs_head)
	 CALL s_list  (ns_b-m, g_new % edges(e) % axs_head)
      ENDDO

      !  Domain point lists

      np_d = 0

      m = np_d
      CALL set_p2i (np_d,   g_old % pnt_head)
      CALL p_list  (np_d-m, g_new % pnt_head)
      m = np_d
      CALL set_p2i (np_d,   g_old % axp_head)
      CALL p_list  (np_d-m, g_new % axp_head)

      !  Domain simplex lists

      ns_d = 0

      m = ns_d
      CALL set_s2i (ns_d,   g_old % int_head)
      CALL s_list  (ns_d-m, g_new % int_head)
      m = ns_d
      CALL set_s2i (ns_d,   g_old % ext_head)
      CALL s_list  (ns_d-m, g_new % ext_head)
      m = ns_d
      CALL set_s2i (ns_d,   g_old % axs_head)
      CALL s_list  (ns_d-m, g_new % axs_head)


      !  ###  SET I2P AND I2S  ###


      ALLOCATE ( i2p_b(np_b), i2s_b(ns_b) )
      ALLOCATE ( i2p_d(np_d), i2s_d(ns_d) )


      !  Set i2p, boundary

      n = 0

      DO e=1,ne
	 CALL set_i2p (n, g_new % edges(e) % pnt_head, i2p_b)
	 CALL set_i2p (n, g_new % edges(e) % axp_head, i2p_b)
      ENDDO

!      WRITE (*,*) 'b pnt number = ',n, np_b

      !  Set i2s, boundary

      n = 0

      DO e=1,ne
	 CALL set_i2s (n, g_new % edges(e) % int_head, i2s_b)
	 CALL set_i2s (n, g_new % edges(e) % ext_head, i2s_b)
	 CALL set_i2s (n, g_new % edges(e) % axs_head, i2s_b)
      ENDDO

!      WRITE (*,*) 'b spx number = ',n, ns_b

      !  Set i2p, domain

      n = 0

      CALL set_i2p (n, g_new % pnt_head, i2p_d)
      CALL set_i2p (n, g_new % axp_head, i2p_d)

!      WRITE (*,*) 'd pnt number = ',n, np_d

      !  Set i2s, domain

      n = 0

      CALL set_i2s (n, g_new % int_head, i2s_d)
      CALL set_i2s (n, g_new % ext_head, i2s_d)
      CALL set_i2s (n, g_new % axs_head, i2s_d)

!      WRITE (*,*) 'd spx number = ',n, ns_d


      !  ###  COPY DATA  ###


      !  Copy boundary point data

      DO e=1,ne
	 CALL p_copy (g_old % edges(e) % pnt_head, i2p_b, i2s_b, i2p_d)
	 CALL p_copy (g_old % edges(e) % axp_head, i2p_b, i2s_b, i2p_d)
      ENDDO

      !  Copy boundary simplex data

      DO e=1,ne
	 CALL s_copy (g_old % edges(e) % int_head, nd_b, i2p_b, i2s_b, i2s_d)
	 CALL s_copy (g_old % edges(e) % ext_head, nd_b, i2p_b, i2s_b, i2s_d)
	 CALL s_copy (g_old % edges(e) % axs_head, nd_b, i2p_b, i2s_b, i2s_d)
      ENDDO

      !  Copy domain point data

      CALL p_copy (g_old % pnt_head, i2p_d, i2s_d, i2p_d)
      CALL p_copy (g_old % axp_head, i2p_d, i2s_d, i2p_d)

      !  Copy domain simplex data

      CALL s_copy (g_old % int_head, nd_d, i2p_d, i2s_d, i2s_d)
      CALL s_copy (g_old % ext_head, nd_d, i2p_d, i2s_d, i2s_d)
      CALL s_copy (g_old % axs_head, nd_d, i2p_d, i2s_d, i2s_d)


      DEALLOCATE ( i2p_b, i2s_b, i2p_d, i2s_d )


   END SUBROUTINE g_copy
!같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같

   !   ROUTINES TO SAVE, LOAD, AND COPY A SINGLE POINT/SIMPLEX


   SUBROUTINE p_save (idf, ascii, hp)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: idf
      LOGICAL, INTENT(IN) :: ascii
      TYPE (point), TARGET, INTENT (IN) :: hp

      INTEGER :: old, ppp, spx
      TYPE (point), POINTER :: p


      p  =>  hp

      IF ( ascii ) THEN

         DO; p  =>  p % next; IF ( ASSOCIATED (p, hp) ) EXIT

	    old = 0; ppp = 0; spx = 0

	    IF ( ASSOCIATED ( p % old ) ) old = p % old % index
	    IF ( ASSOCIATED ( p % ppp ) ) ppp = p % ppp % index
	    IF ( ASSOCIATED ( p % spx ) ) spx = p % spx % index

	    WRITE (idf,*)
	    WRITE (idf,*) p % visited
	    WRITE (idf,*) old, ppp, spx
	    WRITE (idf,*) p % cnd, p % ibk, p % parent
	    WRITE (idf,*) p % leng
	    WRITE (idf,*) p % x
	    WRITE (idf,*) p % metric

         ENDDO

      ELSE

         DO; p  =>  p % next; IF ( ASSOCIATED (p, hp) ) EXIT

	    old = 0; ppp = 0; spx = 0

	    IF ( ASSOCIATED ( p % old ) ) old = p % old % index
	    IF ( ASSOCIATED ( p % ppp ) ) ppp = p % ppp % index
	    IF ( ASSOCIATED ( p % spx ) ) spx = p % spx % index

	    WRITE (idf) p % visited
	    WRITE (idf) old, ppp, spx, p % cnd, p % ibk, p % parent
	    WRITE (idf) p % leng, p % x, p % metric

         ENDDO

      ENDIF

   END SUBROUTINE p_save

   !------------------------------------

   SUBROUTINE s_save (idf, ascii, nd, hs)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: idf, nd
      LOGICAL, INTENT(IN) :: ascii
      TYPE (simplex), TARGET, INTENT (IN) :: hs

      INTEGER :: i, int, ext
      TYPE (simplex), POINTER :: s
      INTEGER, DIMENSION(nd_d+1) :: s2p, s2s


      s  =>  hs

      IF ( ascii ) THEN

         DO; s  =>  s % next; IF ( ASSOCIATED (s, hs) ) EXIT

	        int = 0; ext = 0

	        IF ( ASSOCIATED ( s % int ) ) int = s % int % index
	        IF ( ASSOCIATED ( s % ext ) ) ext = s % ext % index

	        DO i=1,nd+1
	          IF ( ASSOCIATED ( s % opp(i) % pnt ) ) THEN
	            s2p(i) = s % opp(i) % pnt % index
	          ELSE
	            s2p(i) = 0
	          ENDIF

	          IF ( ASSOCIATED ( s % opp(i) % spx ) ) THEN
	            s2s(i) = s % opp(i) % spx % index
	          ELSE
	            s2s(i) = 0
	          ENDIF
	        ENDDO

	        WRITE (idf,*)
	        WRITE (idf,*) int, ext
	        WRITE (idf,*) (s2p(i),i=1,nd+1), (s2s(i),i=1,nd+1)
	        WRITE (idf,*) s % status
	        WRITE (idf,*) s % cmat
	        WRITE (idf,*) s % det
	        WRITE (idf,*) s % fmat

         ENDDO

      ELSE

         DO; s  =>  s % next; IF ( ASSOCIATED (s, hs) ) EXIT

	         int = 0; ext = 0

	         IF ( ASSOCIATED ( s % int ) ) int = s % int % index
	         IF ( ASSOCIATED ( s % ext ) ) ext = s % ext % index

	         DO i=1,nd+1
	           IF ( ASSOCIATED ( s % opp(i) % pnt ) ) THEN
	             s2p(i) = s % opp(i) % pnt % index
	           ELSE
	             s2p(i) = 0
	           ENDIF

	           IF ( ASSOCIATED ( s % opp(i) % spx ) ) THEN
	             s2s(i) = s % opp(i) % spx % index
	           ELSE
	             s2s(i) = 0
	           ENDIF
	         ENDDO

	         WRITE (idf) int, ext, (s2p(i),i=1,nd+1), (s2s(i),i=1,nd+1)
	         WRITE (idf) s % status, s % cmat
	         WRITE (idf) s % det, s % fmat

          ENDDO

      ENDIF

   END SUBROUTINE s_save

   !----------------------------------------------------

   SUBROUTINE p_load (idf, ascii, m, n, i2p, i2s, i2ppp)

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: idf, m, n
      LOGICAL, INTENT(IN) :: ascii
      TYPE (p_pnter), DIMENSION(:), INTENT(IN) :: i2p, i2ppp
      TYPE (s_pnter), DIMENSION(:), INTENT(IN) :: i2s

      INTEGER :: i, old, ppp, spx
      TYPE (point), POINTER :: p


      IF ( ascii ) THEN

         DO i=m,n

	    p  =>  i2p(i) % pnter

	    READ  (idf,*)
	    READ  (idf,*) p % visited
	    READ  (idf,*) old, ppp, spx
	    READ  (idf,*) p % cnd, p % ibk, p % parent
	    READ  (idf,*) p % leng
	    READ  (idf,*) p % x
	    READ  (idf,*) p % metric

	    IF ( old > 0 ) THEN
	       p % old  =>  i2p(old) % pnter
	    ELSE
	       NULLIFY ( p % old )
	    ENDIF

	    IF ( ppp > 0 ) THEN
	       p % ppp  =>  i2ppp(ppp) % pnter
	    ELSE
	       NULLIFY ( p % ppp )
	    ENDIF

	    IF ( spx > 0 ) THEN
	       p % spx  =>  i2s(spx) % pnter
	    ELSE
	       NULLIFY ( p % spx )
	    ENDIF

         ENDDO

      ELSE

         DO i=m,n

	    p  =>  i2p(i) % pnter

	    READ  (idf) p % visited
	    READ  (idf) old, ppp, spx, p % cnd, p % ibk, p % parent
	    READ  (idf) p % leng, p % x, p % metric

	    IF ( old > 0 ) THEN
	       p % old  =>  i2p(old) % pnter
	    ELSE
	       NULLIFY ( p % old )
	    ENDIF

	    IF ( ppp > 0 ) THEN
	       p % ppp  =>  i2ppp(ppp) % pnter
	    ELSE
	       NULLIFY ( p % ppp )
	    ENDIF

	    IF ( spx > 0 ) THEN
	       p % spx  =>  i2s(spx) % pnter
	    ELSE
	       NULLIFY ( p % spx )
	    ENDIF

         ENDDO

      ENDIF

   END SUBROUTINE p_load

   !-------------------------------------------------------

   SUBROUTINE s_load (idf, ascii, nd, m, n, i2p, i2s, i2sss)

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: idf, nd, m, n
      LOGICAL, INTENT(IN) :: ascii
      TYPE (p_pnter), DIMENSION(:), INTENT(IN) :: i2p
      TYPE (s_pnter), DIMENSION(:), INTENT(IN) :: i2s, i2sss

      INTEGER :: j, i, int, ext
      TYPE (simplex), POINTER :: s
      INTEGER, DIMENSION(nd_d+1) :: s2p, s2s


      IF ( ascii ) THEN

         DO j=m,n

	    s  =>  i2s(j) % pnter

	    READ  (idf,*)
	    READ  (idf,*) int, ext
	    READ  (idf,*) (s2p(i),i=1,nd+1), (s2s(i),i=1,nd+1)
	    READ  (idf,*) s % status
	    READ  (idf,*) s % cmat
	    READ  (idf,*) s % det
	    READ  (idf,*) s % fmat

	    IF ( int > 0 ) THEN
	       s % int  =>  i2sss(int) % pnter
	    ELSE
	       NULLIFY ( s % int )
	    ENDIF

	    IF ( ext > 0 ) THEN
	       s % ext  =>  i2sss(ext) % pnter
	    ELSE
	       NULLIFY ( s % ext )
	    ENDIF

	    DO i=1,nd+1
	       IF ( s2p(i) > 0 ) THEN
	          s % opp(i) % pnt  =>  i2p(s2p(i)) % pnter
	       ELSE
	          NULLIFY ( s % opp(i) % pnt )
	       ENDIF

	       IF ( s2s(i) > 0 ) THEN
	          s % opp(i) % spx  =>  i2s(s2s(i)) % pnter
	       ELSE
	          NULLIFY ( s % opp(i) % spx )
	       ENDIF
	    ENDDO

         ENDDO

      ELSE

         DO j=m,n

	    s  =>  i2s(j) % pnter

	    READ  (idf) int, ext, (s2p(i),i=1,nd+1), (s2s(i),i=1,nd+1)
	    READ  (idf) s % status, s % cmat
	    READ  (idf) s % det, s % fmat

	    IF ( int > 0 ) THEN
	       s % int  =>  i2sss(int) % pnter
	    ELSE
	       NULLIFY ( s % int )
	    ENDIF

	    IF ( ext > 0 ) THEN
	       s % ext  =>  i2sss(ext) % pnter
	    ELSE
	       NULLIFY ( s % ext )
	    ENDIF

	    DO i=1,nd+1
	       IF ( s2p(i) > 0 ) THEN
	          s % opp(i) % pnt  =>  i2p(s2p(i)) % pnter
	       ELSE
	          NULLIFY ( s % opp(i) % pnt )
	       ENDIF

	       IF ( s2s(i) > 0 ) THEN
	          s % opp(i) % spx  =>  i2s(s2s(i)) % pnter
	       ELSE
	          NULLIFY ( s % opp(i) % spx )
	       ENDIF
	    ENDDO

         ENDDO

      ENDIF

   END SUBROUTINE s_load

   !------------------------------------------

   SUBROUTINE p_copy (hp_old, i2p, i2s, i2ppp)

      IMPLICIT NONE

      TYPE (point), TARGET, INTENT (IN) :: hp_old
      TYPE (p_pnter), DIMENSION(:), INTENT(IN) :: i2p, i2ppp
      TYPE (s_pnter), DIMENSION(:), INTENT(IN) :: i2s

      INTEGER :: old, ppp, spx
      TYPE (point), POINTER :: pold, pnew


      pold  =>  hp_old

      DO; pold  =>  pold % next; IF ( ASSOCIATED (pold, hp_old) ) EXIT

	 pnew  =>  i2p(pold % index) % pnter

	 pnew % visited = pold % visited
	 pnew % cnd     = pold % cnd
	 pnew % ibk     = pold % ibk
	 pnew % parent  = pold % parent
	 pnew % leng    = pold % leng
	 pnew % x       = pold % x
	 pnew % metric  = pold % metric

	 old = 0; ppp = 0; spx = 0

	 IF ( ASSOCIATED ( pold % old ) ) old = pold % old % index
	 IF ( ASSOCIATED ( pold % ppp ) ) ppp = pold % ppp % index
	 IF ( ASSOCIATED ( pold % spx ) ) spx = pold % spx % index

	 IF ( old > 0 ) THEN
	    pnew % old  =>  i2p (old) % pnter
	 ELSE
	    NULLIFY ( pnew % old )
	 ENDIF

	 IF ( ppp > 0 ) THEN
	    pnew % ppp  =>  i2ppp (ppp) % pnter
	 ELSE
	    NULLIFY ( pnew % ppp )
	 ENDIF

	 IF ( spx > 0 ) THEN
	    pnew % spx  =>  i2s (spx) % pnter
	 ELSE
	    NULLIFY ( pnew % spx )
	 ENDIF

      ENDDO

   END SUBROUTINE p_copy

   !----------------------------------------------

   SUBROUTINE s_copy (hs_old, nd, i2p, i2s, i2sss)

      IMPLICIT NONE

      TYPE (simplex), TARGET, INTENT (IN) :: hs_old
      INTEGER, INTENT(IN) :: nd
      TYPE (p_pnter), DIMENSION(:), INTENT(IN) :: i2p
      TYPE (s_pnter), DIMENSION(:), INTENT(IN) :: i2s, i2sss

      INTEGER :: i, int, ext, p, s
      TYPE (simplex), POINTER :: sold, snew


      sold  =>  hs_old

      DO; sold  =>  sold % next; IF ( ASSOCIATED (sold, hs_old) ) EXIT

	 snew  =>  i2s(sold % index) % pnter

	 snew % status = sold % status
	 snew % cmat   = sold % cmat
	 snew % det    = sold % det
	 snew % fmat   = sold % fmat

	 int = 0; ext = 0

	 IF ( ASSOCIATED ( sold % int ) ) int = sold % int % index
	 IF ( ASSOCIATED ( sold % ext ) ) ext = sold % ext % index

	 IF ( int > 0 ) THEN
	    snew % int  =>  i2sss(int) % pnter
	 ELSE
	    NULLIFY ( snew % int )
	 ENDIF

	 IF ( ext > 0 ) THEN
	    snew % ext  =>  i2sss(ext) % pnter
	 ELSE
	    NULLIFY ( snew % ext )
	 ENDIF

	 DO i=1,nd+1

	    IF ( ASSOCIATED ( sold % opp(i) % pnt ) ) THEN
	       p = sold % opp(i) % pnt % index
	       snew % opp(i) % pnt  =>  i2p(p) % pnter
	    ELSE
	       NULLIFY ( snew % opp(i) % pnt )
	    ENDIF

	    IF ( ASSOCIATED ( sold % opp(i) % spx ) ) THEN
	       s = sold % opp(i) % spx % index
	       snew % opp(i) % spx  =>  i2s(s) % pnter
	    ELSE
	       NULLIFY ( snew % opp(i) % spx )
	    ENDIF
	 ENDDO

      ENDDO

   END SUBROUTINE s_copy


   !   ROUTINES TO ALLOCATE A NEW LIST


   SUBROUTINE p_list (n, p_head)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      TYPE (point), POINTER :: p_head

      INTEGER :: i
      TYPE (point), POINTER :: p

       p_head  =>  newhead_point ()

      DO i=1,n
	 p  =>  new_point ()
	 CALL pushtop_point (p_head, p)
      ENDDO

   END SUBROUTINE p_list

   !----------------------------
   SUBROUTINE s_list (n, s_head)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      TYPE (simplex), POINTER :: s_head

      INTEGER :: i
      TYPE (simplex), POINTER :: s

       s_head  =>  newhead_simplex ()

      DO i=1,n
	 s  =>  new_simplex ()
	 CALL pushtop_simplex (s_head, s)
      ENDDO

   END SUBROUTINE s_list


   !   ROUTINES TO RESET INDICES

   !------------------------
   SUBROUTINE clr_p2i (head)

      IMPLICIT NONE

      TYPE (point), TARGET, INTENT (IN) :: head

      TYPE (point), POINTER :: p

      p  =>  head

      DO; p  =>  p % next; IF ( ASSOCIATED (p, head) ) EXIT
	 p % index = 0
      END DO

   END SUBROUTINE clr_p2i

   !------------------------
   SUBROUTINE clr_s2i (head)

      IMPLICIT NONE

      TYPE (simplex), TARGET, INTENT (IN) :: head

      TYPE (simplex), POINTER :: s

      s  =>  head

      DO; s  =>  s % next; IF ( ASSOCIATED (s, head) ) EXIT
	 s % index = 0
      END DO

   END SUBROUTINE clr_s2i


   !   ROUTINES TO SET INDICES

   !---------------------------
   SUBROUTINE set_p2i (n, head)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: n
      TYPE (point), TARGET, INTENT (IN) :: head

      TYPE (point), POINTER :: p
      INTEGER :: i

      i = n; p  =>  head
     ! write(*,*)'P2I, head involved is',head % x
      DO; p  =>  p % next; IF ( ASSOCIATED (p, head) ) EXIT
         i = i+1; p % index = i
       ! WRITE(*,*)'P',p % index,p % x
      END DO
      n = i

   END SUBROUTINE set_p2i

   !---------------------------
   SUBROUTINE set_s2i (n, head)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: n
      TYPE (simplex), TARGET, INTENT (IN) :: head

      TYPE (simplex), POINTER :: s
      INTEGER :: i

      i = n;  s  =>  head
      DO; s  =>  s % next; IF ( ASSOCIATED (s, head) ) EXIT
         i = i+1; s % index = i
      END DO
      n = i

   END SUBROUTINE set_s2i


   !   ROUTINES TO SET POINTERS

   !--------------------------------
   SUBROUTINE set_i2p (n, head, i2p)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: n
      TYPE (point), POINTER :: head
      TYPE (p_pnter), DIMENSION(:), INTENT(INOUT) :: i2p

      INTEGER :: i
      TYPE (point), POINTER :: p

      i = n; p  =>  head
      DO; p  =>  p % next; IF ( ASSOCIATED (p, head) ) EXIT
	 i = i+1; i2p(i) % pnter  =>  p
      ENDDO
      n = i

   END SUBROUTINE set_i2p

   !--------------------------------
   SUBROUTINE set_i2s (n, head, i2s)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: n
      TYPE (simplex), TARGET, INTENT (IN) :: head
      TYPE (s_pnter), DIMENSION(:), INTENT (INOUT) :: i2s

      INTEGER :: i
      TYPE (simplex), POINTER :: s

      i = n; s  =>  head
      DO; s  =>  s % next; IF ( ASSOCIATED (s, head) ) EXIT
	 i = i+1; i2s(i) % pnter  =>  s
      ENDDO
      n = i

   END SUBROUTINE set_i2s


   !   ROUTINES NOT USED (ALTERNATE VERSIONS)

   
   SUBROUTINE set_i2p_a (head, i2p)

      IMPLICIT NONE

      TYPE (point), TARGET, INTENT(IN) :: head
      TYPE (p_pnter), DIMENSION(:), INTENT(INOUT) :: i2p

      INTEGER :: i
      TYPE (point), POINTER :: p

      p  =>  head
      DO i=1,SIZE(i2p)
	 p  =>  p % next; i2p(i) % pnter  =>  p
      ENDDO

   END SUBROUTINE set_i2p_a


   SUBROUTINE set_i2s_a (head, i2s)

      IMPLICIT NONE

      TYPE (simplex), TARGET, INTENT (IN) :: head
      TYPE (s_pnter), DIMENSION(:), INTENT (INOUT) :: i2s

      INTEGER :: i
      TYPE (simplex), POINTER :: s

      s  =>  head
      DO i=1,SIZE(i2s)
	 s  =>  s % next; i2s(i) % pnter  =>  s
      ENDDO

   END SUBROUTINE set_i2s_a


   SUBROUTINE set_i2p_b (n, head, i2p)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: n
      TYPE (point), TARGET, INTENT(IN) :: head
      TYPE (p_pnter), DIMENSION(:), INTENT(INOUT) :: i2p

      INTEGER :: i
      TYPE (point), POINTER :: p

      i = n; p  =>  head
      DO; p  =>  p % next; IF ( ASSOCIATED (p, head) ) EXIT
	 i = i+1; i2p(p % index) % pnter  =>  p
      ENDDO
      n = i

   END SUBROUTINE set_i2p_b


   SUBROUTINE set_i2s_b (n, head, i2s)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: n
      TYPE (simplex), TARGET, INTENT (IN) :: head
      TYPE (s_pnter), DIMENSION(:), INTENT (INOUT) :: i2s

      INTEGER :: i
      TYPE (simplex), POINTER :: s

      i = n; s  =>  head
      DO; s  =>  s % next; IF ( ASSOCIATED (s, head) ) EXIT
	 i = i+1; i2s(s % index) % pnter  =>  s
      ENDDO
      n = i

   END SUBROUTINE set_i2s_b


END MODULE convert
