MODULE eucl_delaunay

  USE linsys      !  Direct Linear system solver

  USE grid_types  !  Grid types

  USE list        !  List operations

  USE delaunay    !  Delaunay Basic module

  IMPLICIT NONE
  
  PUBLIC  :: segment_length,segment_r_length, simplex_center, insert_point,&
            & insert_far_point

  PRIVATE :: delaunay_measure, delete_cavity, check_cavity

  REAL(KIND=8), PARAMETER, PUBLIC  :: e_tol = 1.50d0
  REAL(KIND=8), PARAMETER, PRIVATE :: m_tol = 0.50d0

  INTEGER,       DIMENSION(maxs)     , PRIVATE :: pivot
  REAL (KIND=8), DIMENSION(maxs)     , PRIVATE :: vector, center
  REAL (KIND=8), DIMENSION(maxs,maxs), PRIVATE :: matrix

CONTAINS
  !
  !********************************************************************
  !
  FUNCTION segment_length ( nd, a, b ) RESULT ( leng2 )

    ! +-------------------+
    ! |   EUCLID LENGTH   | SQUARED
    ! +-------------------+

    INTEGER, INTENT(IN) :: nd
    TYPE (point), INTENT(IN) :: a, b
    REAL(KIND=8) :: leng2

    INTEGER :: i
    REAL(KIND=8) :: r

    r = 0.d0
    DO i = 1, nd
      r = r + ( a % x(i) - b % x(i) )**2
    ENDDO
    leng2 = r

    IF ( leng2 .lt. 0.0 ) THEN
       WRITE (*,*) 'segment_length error, lenght**2 = ',leng2
       WRITE (*,*) '(stop)'
       STOP
    ENDIF

  END FUNCTION segment_length
  !
  !*******************************************************************
  !
!by DD:begin
  
  FUNCTION segment_r_length (nd, a, b) RESULT (leng2)

    ! Length is NOT squared

    INTEGER, INTENT(IN) :: nd
    TYPE (point), INTENT(IN) :: a, b
    REAL (KIND=8) :: leng2, l_1, l_2

    INTEGER :: j, k
    REAL (KIND=8) :: temp

    !  Calcul du vecteur distance (a-b)

    temp = zero
    DO j = 1, nd
       vector(j) = a % x(j) - b % x(j)
       temp = temp+vector(j)**2
    ENDDO

    temp = SQRT(temp)
    IF(temp < 1.d-12) THEN
       leng2 = zero
       RETURN
    ENDIF

    !leng2 = zero
    l_1   = zero
    l_2   = zero
    DO j = 1, nd
       DO k = 1, nd
          l_1   = l_1+vector(j)*a%metric(j,k)*vector(k)
          l_2 = l_2+vector(j)*b%metric(j,k)*vector(k)
          !leng2 = leng2+vector(j)*(a%metric(j,k)+b%metric(j,k))*vector(k)
       ENDDO
    ENDDO
    leng2 = .5d0*(SQRT(l_1)+SQRT(l_2))
    !leng2 = SQRT(.5d0*leng2)

  END FUNCTION segment_r_length

!by DD: end

  !
  !********************************************************************
  !
  SUBROUTINE simplex_center (nd, simp, centre, radius)

    ! +--------------------------------------+
    ! |   EUCLID SIMPLEX CENTER AND RADIUS   |
    ! +--------------------------------------+

    INTEGER, INTENT(IN) :: nd
    TYPE (simplex), TARGET, INTENT(INOUT) :: simp
    REAL (KIND=8), INTENT(OUT) :: radius
    REAL (KIND=8), DIMENSION(ND), INTENT(OUT) :: centre

    INTEGER :: i, j
    REAL (KIND=8) :: eps, sum, det
    REAL (KIND=8), DIMENSION(nd) :: x1, x2
    PARAMETER (eps = 1.0d-12)


    DO i=1,nd
       x1  =  simp % opp(i) % pnt % x
       x2  =  simp % opp(i+1) % pnt % x
       sum = 0.
       DO j=1,nd
          matrix(i,j) = x2(j)-x1(j)
          sum = sum + x2(j)**2-x1(j)**2
       ENDDO
       centre(i) = 0.5*sum
    ENDDO

    CALL lufact (nd, pivot, matrix, det)

    IF ( ABS(det) <= eps ) THEN
       WRITE (911,'(a30,e12.5)') 'simplex_center warning, det = ',det
       IF ( ABS(det) <= 0.0 ) THEN
          WRITE (*,*) ' ERROR. Negative area. Stop.'
          STOP
       ENDIF
    ENDIF

    CALL lusolv (nd, pivot, matrix, centre)

    sum = 0.
    DO i=1,nd+1
       x1  =  simp % opp(i) % pnt % x
       DO j=1,nd
          sum = sum + (centre(j)-x1(j))**2
       ENDDO
    ENDDO
    radius = sum/DBLE(nd+1)

  END SUBROUTINE simplex_center
  !
  !********************************************************************
  !
  FUNCTION delaunay_measure (nd, poin, simp) RESULT (measure)

    ! +--------------------------------------------+
    ! |   EUCLID MEASURE = DISTANCE**2/RADIUS**2   |
    ! +--------------------------------------------+

    INTEGER, INTENT(IN) :: nd
    TYPE (point), TARGET, INTENT(IN) :: poin
    TYPE (simplex), TARGET, INTENT(INOUT) :: simp
    REAL (KIND=8) :: measure

    INTEGER :: i, j
    REAL (KIND=8) :: eps, sum, det, radius, distance
    !REAL (KIND=8), DIMENSION(nd) :: x1, x2 ! by DD :casi 1D e 2D
    ! le coordinate dei punti sono sempre 2!!!
    REAL (KIND=8), DIMENSION(2) :: x1, x2
    PARAMETER (eps = 1.0d-12)


    DO i=1,nd
       x1  =  simp % opp(i) % pnt % x
       x2  =  simp % opp(i+1) % pnt % x
       sum = 0.
       DO j=1,nd
          matrix(i,j) = x2(j)-x1(j)
          sum = sum + x2(j)**2-x1(j)**2
       ENDDO
       center(i) = 0.5*sum
    ENDDO

    CALL lufact (nd, pivot, matrix, det)

    IF ( ABS(det) <= eps ) THEN
       WRITE (911,'(a30,e12.5)') 'simplex_center warning, det = ',det
       IF ( det <= 0.0 ) THEN
          WRITE (*,*) ' ERROR. Negative area. Stop.'
          STOP
       ENDIF
    ENDIF

    CALL lusolv (nd, pivot, matrix, center)

    sum = 0.
    DO i=1,nd+1
       x1  =  simp % opp(i) % pnt % x
       DO j=1,nd
          sum = sum + (center(j)-x1(j))**2
       ENDDO
    ENDDO
    radius = sum/DBLE(nd+1)

    distance = 0.
    DO j=1,nd
       distance = distance + (poin % x(j) - center(j))**2
    ENDDO

    measure = distance/radius

  END FUNCTION delaunay_measure
  !
  !********************************************************************
  !
  SUBROUTINE insert_point (nd, spx_head, pnew)

    ! +--------------------------+
    ! |    INSERT A NEW POINT    |
    ! +--------------------------+

    INTEGER       , INTENT(IN) :: nd
    TYPE (simplex), POINTER    :: spx_head
    TYPE (point)  , POINTER    :: pnew

    LOGICAL :: err_status

    IF ( ASSOCIATED (pnew % spx) ) THEN
       WRITE (*,*) 'New point is already connected (insert_point) -- STOP --'
       STOP
    END IF

    CALL delete_first    (nd, spx_head, pnew, err_status)
    IF(err_status) THEN
       WRITE(*,*)'Point not found (insert_point) -- STOP --'
       WRITE(*,*)'or new point outside the simplexes'' list (by DD)'
       write(*,*)'pnew', pnew % x, pnew % index
       STOP
    ENDIF

    CALL delete_cavity   (nd, pnew)

    CALL update_topology (nd, pnew)
    CALL update_geometry (nd, new_head)

    CALL jointop_simplex (gar_head, del_head)
    CALL jointop_simplex (spx_head, new_head)

  END SUBROUTINE insert_point
  !
  !********************************************************************
  !
  SUBROUTINE insert_far_point (nd, spx_head, pnew, far_point)

    ! +-----------------------------------------------+
    ! |    INSERT A NEW POINT, IF IT IS FAR ENOUGH    |
    ! +-----------------------------------------------+

    INTEGER       , INTENT(IN)  :: nd
    TYPE (simplex), POINTER     :: spx_head
    TYPE (point)  , POINTER     :: pnew
    LOGICAL       , INTENT(OUT) :: far_point

    LOGICAL :: err_status

    IF ( ASSOCIATED (pnew % spx) ) THEN
       WRITE (*,*) 'New point is already connected (insert_far_point) -- STOP --'
       STOP
    ENDIF

    CALL delete_first    (nd, spx_head, pnew, err_status)
    IF(err_status) THEN
       !WRITE(*,*)'Point not found (insert_far_point)'
       far_point = .FALSE.
       RETURN
    ENDIF

    CALL delete_cavity   (nd, pnew)
    CALL check_cavity    (nd, pnew, far_point)

    IF (far_point) THEN
       CALL update_topology (nd, pnew)
       CALL update_geometry (nd, new_head)

       CALL jointop_simplex (gar_head, del_head)
       CALL jointop_simplex (spx_head, new_head)
    ELSE
       CALL jointop_simplex (spx_head, del_head)
    ENDIF

  END SUBROUTINE insert_far_point
  !
  !********************************************************************
  !
  SUBROUTINE delete_cavity (nd, pnew)

    ! +-----------------------+
    ! |    CAVITY DELETION    |
    ! +-----------------------+

    INTEGER, INTENT(IN) :: nd
    TYPE (point), TARGET, INTENT(IN) :: pnew

    INTEGER :: i
    TYPE (simplex), POINTER :: simp, neig
    REAL (KIND=8), PARAMETER :: eps = 1.0d-10, limit = 2.d0 - eps
    REAL (KIND=8) :: pmeas, qmeas

    ! Status table:

    ! status  > 0 :: unvisited simplex
    ! status == 0 :: frozen simplex
    ! status  < 0 :: visited simplex

    ! Link the deleted simplices
    ! Set to "visited" the deleted and the adjacent simplices

    simp  =>  del_head
    DO; simp  =>  simp % next; IF ( ASSOCIATED (simp, del_head) ) EXIT
       DO i=1,nd+1
          neig  =>  simp % opp(i) % spx  !controlla gli adiacenti della base...
          IF ( neig % status > 0 ) THEN  !solo quelli non visitati,cioè quelli  
             neig % status = - neig % status !esterni alla base e li pone "visited"

             !  Euclidean metric
             pmeas = delaunay_measure ( nd, pnew, neig )
             qmeas = delaunay_measure ( nd, pnew, neig )

             IF ( ( pmeas+qmeas ) < limit ) THEN  !pto contenuto nel cerchio
                                                  !circoscritto del simpl
                CALL movebot_simplex (del_head, neig) 
               !mette sotto la testa i simplessi esterni alla base il cui ccc(cerchio
               !circoscritto)contiene il punto pnew
               
                                                  
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    ! pone "unvisited" i simplessi della base
    ! lasciando "visited" quelli esterni alla base

    simp  =>  del_head
    DO; simp  =>  simp % next; IF ( ASSOCIATED (simp, del_head) ) EXIT
       simp % status = ABS(simp % status)
    ENDDO

  END SUBROUTINE delete_cavity !poi update topology
  !
  !********************************************************************
  !
  SUBROUTINE check_cavity (nd, pnew, far_point)

    ! +------------------------------+
    ! |   CHECK THE DELETED CAVITY   |
    ! +------------------------------+

    INTEGER     , INTENT(IN)  :: nd
    TYPE (point), INTENT(IN)  :: pnew
    LOGICAL     , INTENT(OUT) :: far_point

    LOGICAL :: undelete
    INTEGER :: i
    TYPE (simplex), POINTER :: s
    REAL (KIND=8) :: pl2, leng2

    undelete = .FALSE.

    pl2 = (pnew%leng*m_tol)**2

    s  =>  del_head
    DO; s  =>  s % next; IF ( ASSOCIATED(s, del_head) ) EXIT

       DO i=1,nd+1

          !  Euclidean metric

          leng2 = segment_length(nd, s % opp(i) % pnt, pnew)   ! SQUARED

          IF ( leng2 < pl2 ) THEN
             undelete = .TRUE.
             EXIT
          ENDIF
       ENDDO
    ENDDO

    IF ( undelete ) THEN
       s  =>  del_head
       DO; s  =>  s % next; IF ( ASSOCIATED (s, del_head) ) EXIT
          DO i=1,nd+1
             s % opp(i) % spx % status = ABS(s % opp(i) % spx % status)
          ENDDO
       ENDDO
       far_point = .FALSE.
    ELSE
       far_point = .TRUE.
    ENDIF

  END SUBROUTINE check_cavity
  !
  !********************************************************************
  !
END MODULE eucl_delaunay
