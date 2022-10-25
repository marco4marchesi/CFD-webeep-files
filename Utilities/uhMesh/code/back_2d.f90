MODULE back_2d

  USE grid_types      !  Types module
  USE list            !  List module
  USE delaunay        !  Delaunay basic module
  USE eucl_delaunay   !  Delaunay Euclide module
  USE convert         !  Convert module
  USE metric_2d       !  Metric module
  use init_data, only: nbackgr, backgrP, read_record

  IMPLICIT NONE

  PUBLIC  :: build_back, eval_back

CONTAINS

  SUBROUTINE build_back ( grid, back )

    REAL(KIND=8), PARAMETER :: max_l = 1.5d0

    TYPE (domain), INTENT(IN)  :: grid
    TYPE (domain), INTENT(INOUT) :: back

    INTEGER                            :: i, j, n
    REAL(KIND=8)                       :: big, leng
    REAL(KIND=8) , DIMENSION(2*nd_d-1) :: a
    TYPE(point)  , POINTER             :: hp, p, q
    TYPE(simplex), POINTER             :: hs, s

    !  Boundary grid lengths
    !inserted points list
    big = HUGE(big); hp  =>  grid % pnt_head; p  =>  hp
    DO; p  =>  p % next; IF ( ASSOCIATED (p, hp) ) EXIT
       p % leng = big
    ENDDO

    DO i=1,grid % number_of_edges

       hs  =>  grid % edges(i) % int_head; s  =>  hs
       DO; s  =>  s % next; IF ( ASSOCIATED (s, hs) ) EXIT

          p  =>  s % opp(2) % pnt % ppp ! 2d points

          q  =>  s % opp(1) % pnt % ppp
          !  p        q
          !  @--------@
          !  2        1

          !  Euclidean metric => linear SQUARED length
          !  associated to p only

          leng = segment_length(nd_d, p, q)
          leng = SQRT(leng)

          p % leng = MIN( p % leng, leng )


          ! euclidean length is assigned to 2d point associated to 
          ! the 1d simplex's second point 


       ENDDO

    ENDDO

    !  Copy boundary grid into back grid

    CALL g_copy (grid, back)!shows boundary&domain,point&simplex number(grid_2d

    p  =>  hp; q  =>  back % pnt_head
    DO; p  =>  p % next; q  =>  q % next; IF ( ASSOCIATED (p, hp) ) EXIT
       p % ppp     =>  q
       q % leng    = p % leng
       q % ellipse = p % ellipse
       q % metric  = p % metric
    ENDDO


    !  Additional back grid points and lengths
    !   (domain.gridaname) 
    hs  =>  back % int_head

    !   Reads backgrid points
    call read_record ( 'G_MESH_PARAMETERS_DOMAIN' )
    WRITE(*,*)

    DO i = 1, nbackgr
       p  =>  new_point ()
       p % x(1:nd_d) = backgrP(1:nd_d,i)
       a(1) = backgrP(nd_d+1,i)
       a(2) = a(1)
       a(3) = 0 
!
       CALL gen_metr(a, p%ellipse, p%metric)
       p%leng = a(1)         

       CALL insert_point (nd_d, hs, p)
    ENDDO

  END SUBROUTINE build_back
  !같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같

  SUBROUTINE eval_back (p, testa, s_found, err_status)

    TYPE (point)  , TARGET      :: p
    TYPE (simplex), POINTER     :: testa
    TYPE (simplex), POINTER     :: s_found
    LOGICAL       , INTENT(OUT) :: err_status

    INTEGER                 :: i, j, k
    INTEGER                 :: status
    REAL (KIND=8)           :: a, b, c, leng2
    TYPE (simplex), POINTER :: s_back
    ! new point must be inserted inside the backgrid 

    err_status = .TRUE.
    a = 0.d0
    c = 0.d0
    s_back  =>  p % old % ppp % spx

    status = 0				!  0 to activate search on backgrid
    CALL walking_search (nd_d, p, s_back, status)

    IF(status == 1) THEN		! point located
       err_status = .FALSE.
       !WRITE(*,*)'Point located by Walking_search'
    ELSE IF(status == 0) THEN		! walking_search failed
       !write(*,*)'eval_back walking failed: use trivial search'
       CALL trivial_search (nd_d, testa, p, s_back, err_status)
       IF(err_status) THEN
          !WRITE(*,*)'Walking & Trivial search failed (eval_back)'
          RETURN
       ENDIF
       !WRITE(*,*)'Point located by Trivial_search - after walking failed'
    ELSE IF(status == -1) THEN		! if point not exist in grid
       !WRITE(*,*)'Walking says: point is not here! (eval_back)'
       CALL trivial_search (nd_d, testa, p, s_back, err_status)
       IF(err_status) THEN
          !WRITE(*,*)'Walking & Trivial search failed (eval_back)'
          RETURN
       ENDIF
    ENDIF

    IF (s_back%status == ext_status) THEN
       WRITE(*,*) 'ERROR.'
       WRITE(*,*) 'Element status external:', s_back%status
       WRITE(*,*) 'STOP'
       WRITE(*,*) ''
       STOP !research takes place in int_head list,where status=int

    ENDIF

    DO i=1,nd_d+1
       b = s_back % fmat(i,nd_d+1)
       leng2 = 0.d0
       DO j=1,nd_d
          b = b + s_back % fmat(i,j) * p % x(j)
       ENDDO
       b = MAX(0d0,b)  !  p must be internal to s
       IF ( b >= c ) THEN
          c = b
          k = i
       ENDIF
       a = a + b * s_back % opp(i) % pnt % leng
    ENDDO
    p % leng = a

    p % ppp  =>  s_back % opp(k) % pnt
    s_found => s_back

  END SUBROUTINE eval_back




END MODULE back_2d
