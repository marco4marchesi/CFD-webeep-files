MODULE check_grid

  USE delaunay       ! Delaunay Basic Module

  USE eucl_delaunay  ! Delaunay Euclide Module

  USE back_2d        ! Background grid Module

  USE arrays         ! Data structure conversion Module

  IMPLICIT NONE

  INTEGER, PARAMETER :: max_loop = 1000000
  INTEGER, PARAMETER :: s_to_s = 1, s_to_p = 2, s_to_cmat_ii = 3
  INTEGER, PARAMETER :: s_to_cmat_ij = 4, b_neighs = 5

  REAL(kind=8), PARAMETER :: min_value = 1.d-12

  PUBLIC :: check_in_boundary, check_consistency, contability
  
CONTAINS
  
  !같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같�
 
  SUBROUTINE check_in_boundary(grid)

    IMPLICIT NONE

    TYPE(domain), INTENT(in) :: grid

    INTEGER :: i, j, i1, i2
    REAL(kind=8), DIMENSION(nd_d) :: delta


    edge_loop: DO i = 1, grid%number_of_edges

       i1 = grid%edges(i)%begin_vert
       i2 = grid%edges(i)%end_vert

       IF(i1 == i2) THEN

          !  Self closing boundary control

          delta = curv(i)%x(:, 1) - curv(i)%x(:, curv(i)%ns)
          IF(ABS(delta(1)) > min_value .OR. ABS(delta(2)) > min_value) THEN
             WRITE(*,*)
             WRITE(*,*)'*** Problems with self-closing boundary', i
             WRITE(*,*)'Different begin & end vertex'
             STOP
          ENDIF

       ELSE

          !  Open boundary pieces

          !  begin_vertex search in other boundaries end_vertex

          DO j = 1, grid%number_of_edges

             IF(grid%edges(j)%end_vert == i1) THEN

                delta = curv(i)%x(:, 1) - curv(j)%x(:, (curv(j)%ns))
                IF(ABS(delta(1)) > min_value .OR. ABS(delta(2)) > min_value) THEN
                   WRITE(*,*)
                   WRITE(*,*)'*** Problems with boundary', i,'and boundary', j
                   WRITE(*,*)'Begin & end vertex not consistent'
                   STOP
                ENDIF

             ENDIF

          ENDDO

          !  end_vertex search in other boundaries begin_vertex

          DO j = 1, grid%number_of_edges

             IF(grid%edges(j)%begin_vert == i2) THEN

                delta = curv(j)%x(:, 1) - curv(i)%x(:, curv(i)%ns)
                IF(ABS(delta(1)) > min_value .OR. ABS(delta(2)) > min_value) THEN
                   WRITE(*,*)
                   WRITE(*,*)'*** Problems with boundary', i,'and boundary', j
                   WRITE(*,*)'End vertex & begin not consistent'
                   STOP
                ENDIF

             ENDIF

          ENDDO

       ENDIF

    ENDDO edge_loop

    WRITE(*,*)'Boundary data checked.'

  END SUBROUTINE check_in_boundary
  
  !같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같
  
  SUBROUTINE check_consistency(grid)

    IMPLICIT NONE

    TYPE(domain), INTENT(IN) :: grid

    INTEGER :: i, j, k, ii
    INTEGER :: index
    REAL(kind=8) :: module_delta
    REAL(kind=8), DIMENSION(nd_d) :: delta

    TYPE(simplex), POINTER :: p_spx1, p_spx2

    WRITE(*,*)
    WRITE(*,*)'Grid data checking...'

    p_spx1 => grid%int_head

    spx_loop: DO i = 1, max_loop
       p_spx1 => p_spx1%next
       IF(ASSOCIATED(p_spx1, grid%int_head)) EXIT spx_loop

       !  check neigh 

       DO j = 1, nd_d+1

          p_spx2 => p_spx1%opp(j)%spx

          IF(ASSOCIATED(p_spx2)) THEN

             DO k = 1, nd_d+1
                IF(ASSOCIATED(p_spx2%opp(k)%spx, p_spx1)) EXIT
             ENDDO
             !by DD :!!!
             IF(k == nd_d+2) THEN
                CALL out_problem(s_to_s)
             ENDIF
             !by DD end
          ENDIF

       ENDDO

       !  check cmat

       DO j = 1, nd_d+1

          p_spx2 => p_spx1%opp(j)%spx

          IF(ASSOCIATED(p_spx2)) THEN

             index = p_spx1%cmat(j, j)

             IF(.NOT. ASSOCIATED(p_spx2%opp(index)%spx, p_spx1)) THEN
                CALL out_problem(s_to_cmat_ii)
             ENDIF

             DO k = 1, nd_d+1
                IF(k /= j) THEN

                   index = p_spx1%cmat(j, k)

                   delta = p_spx1%opp(k)%pnt%x-p_spx2%opp(index)%pnt%x

                   module_delta = zero

                   DO ii = 1, nd_d
                      module_delta = module_delta+delta(ii)**2
                   ENDDO
                   module_delta = SQRT(module_delta)

                   IF(module_delta > min_value) THEN
                      CALL out_problem(s_to_cmat_ij)
                   ENDIF

                ENDIF

             ENDDO

          ENDIF

       ENDDO

    ENDDO spx_loop
    IF(i >= max_loop) THEN
       WRITE(*,*)'spx_loop/check_consistency failure (max_loop =',max_loop,')'
       STOP
    ENDIF

    !  Check neighs

    bou_loop: DO i = 1, grid%number_of_edges
       p_spx1 => grid%edges(i)%int_head

       id_spx_loop: DO ii = 1, max_loop
          p_spx1 => p_spx1%next
          IF(ASSOCIATED(p_spx1, grid%edges(i)%int_head)) EXIT id_spx_loop

          !  Corresponding internal element

          p_spx2 => p_spx1%int

          k = 0
          DO j =1, nd_d+1
             IF(p_spx2%opp(j)%spx%index /= 0) k = k+1
          ENDDO
          IF(k == nd_d+1) THEN
             CALL out_problem(b_neighs)
          ENDIF

       ENDDO id_spx_loop
       IF(ii >= max_loop) THEN
          WRITE(*,*)'id_spx_loop/check_consistency loop failure (max_loop =',max_loop,')'
          STOP
       ENDIF

    ENDDO bou_loop

    !WRITE(*,*)'Grid data OK.'
    WRITE(*,*)
  CONTAINS
    !
    !-------------------------------------------------------------
    !
    SUBROUTINE out_problem(in_case)

      IMPLICIT NONE

      INTEGER, INTENT(in) :: in_case

      SELECT CASE(in_case)

      CASE(s_to_s)

         WRITE(*,*)
         WRITE(*,*)'blind neighbour'

      CASE(s_to_p)

         WRITE(*,*)
         WRITE(*,*)'point not seen'


      CASE(s_to_cmat_ii)

         WRITE(*,*)
         WRITE(*,*)'cmat diagonal fault'

      CASE(s_to_cmat_ij)

         WRITE(*,*)
         WRITE(*,*)'cmat extradiagonal fault'
         WRITE(*,*)'tolerance', min_value

      CASE(b_neighs)

         WRITE(*,*)
         WRITE(*,*)'neighs fault'
         WRITE(*,*)'Internal found, not pointing to 0'
         WRITE(*,*)'Boundary simplex:',p_spx1%index
         DO j = 1, nd_d
            WRITE(*,*)p_spx1%opp(j)%pnt%ppp%x
         ENDDO
         WRITE(*,*)'Corresponding domain simplex:',p_spx2%index
         WRITE(*,*)'neigh',(p_spx2%opp(j)%spx%index, j = 1, nd_d+1)
         DO j = 1, nd_d+1
            WRITE(*,*)p_spx2%opp(j)%pnt%x
         ENDDO

      END SELECT

    END SUBROUTINE out_problem
    !
    !-------------------------------------------------------------
    !
  END SUBROUTINE check_consistency

  
  !같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같
  
  SUBROUTINE contability(grid, verbose)

    !  Accepted values for "verbose" are:
    !  0 => no messages
    !  1 => only general information
    !  2 => detailed information


    IMPLICIT NONE

    INTEGER     , INTENT(IN) :: verbose
    TYPE(domain), INTENT(IN) :: grid

    INTEGER :: i1, i2, i3, i4, i5
    TYPE(simplex), POINTER :: s

    i1 = 0
    i2 = 0
    i3 = 0
    i4 = 0
    i5 = 0
    s => grid%int_head
    DO
       s => s%next
       IF(ASSOCIATED(s, grid%int_head)) EXIT

       IF(s%status == int_status) THEN
          i1 =i1 +1
       ELSE IF(s%status == und_status) THEN
          i2 =i2+1
          s%status = int_status
       ELSE IF(s%status == ext_status) THEN
          i3 =i3+1
       ELSE IF(s%status == stc_status) THEN
          i4 =i4+1
       ELSE
          i5 =i5+1
       ENDIF

    ENDDO
    
    IF(verbose == 2) THEN

       WRITE(*,*)'INTERNAL LIST flags'
       WRITE(*,*)'internal', i1
       WRITE(*,*)'undefined', i2
       WRITE(*,*)'external', i3
       WRITE(*,*)'stitch', i4
       WRITE(*,*)'unknown', i5
       WRITE(*,*)

    ENDIF

    i1 = 0
    i2 = 0
    i3 = 0
    i4 = 0
    i5 = 0
    s => grid%ext_head
    DO
       s => s%next
       IF(ASSOCIATED(s, grid%ext_head)) EXIT

       IF(s%status == int_status) THEN
          i1 =i1 +1
       ELSE IF(s%status == und_status) THEN
          i2 =i2+1
          s%status = ext_status
       ELSE IF(s%status == ext_status) THEN
          i3 =i3+1
       ELSE IF(s%status == stc_status) THEN
          i4 =i4+1
       ELSE
          i5 =i5+1
       ENDIF

    ENDDO

    IF(verbose == 2) THEN

       WRITE(*,*)'EXTERNAL LIST flags'
       WRITE(*,*)'internal', i1
       WRITE(*,*)'undefined', i2
       WRITE(*,*)'external', i3
       WRITE(*,*)'stitch', i4
       WRITE(*,*)'unknown', i5
       WRITE(*,*)

    ENDIF

    i1 = 0
    i2 = 0
    i3 = 0
    i4 = 0
    i5 = 0
    s => grid%axs_head
    DO
       s => s%next
       IF(ASSOCIATED(s, grid%axs_head)) EXIT

       IF(s%status == int_status) THEN
          i1 =i1 +1
       ELSE IF(s%status == und_status) THEN
          i2 =i2+1
          s%status = ext_status
       ELSE IF(s%status == ext_status) THEN
          i3 =i3+1
       ELSE IF(s%status == stc_status) THEN
          i4 =i4+1
       ELSE IF(s%status == frz_status) THEN
          i5 =i5+1
       ENDIF

    ENDDO

    IF(verbose == 2) THEN

       WRITE(*,*)'AUXILIARY LIST flags'
       WRITE(*,*)'internal', i1
       WRITE(*,*)'undefined', i2
       WRITE(*,*)'external', i3
       WRITE(*,*)'stitch', i4
       WRITE(*,*)'unknown', i5
       WRITE(*,*)
       
    ENDIF

    IF(verbose >= 1) THEN

       WRITE(*,*)'Grid details:'
       WRITE(*,*)'int_head', numberof_simplex(grid%int_head)
       WRITE(*,*)'ext_head', numberof_simplex(grid%ext_head)
       WRITE(*,*)'axs_head', numberof_simplex(grid%axs_head)
       WRITE(*,*)'stc_head', numberof_simplex(grid%stc_head)
       WRITE(*,*)'spx_head', numberof_simplex(grid%spx_head)
    ENDIF

  END SUBROUTINE contability
 
  
END MODULE check_grid



