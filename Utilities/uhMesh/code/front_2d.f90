MODULE front_2d

  USE grid_types     !  Grid types
  USE list           !  List operations
  USE delaunay       !  Delaunay basic module
  USE eucl_delaunay  !  Delaunay Euclid module
  USE back_2d        !  Backgrid module
  USE metric_2d      ! Metric module
  USE refine_2d      !  Refine module
  USE check_grid     !  Check module

  IMPLICIT NONE

  PUBLIC :: front_grid

CONTAINS
  !
  !********************************************************************
  !                                                      steiner
  SUBROUTINE front_grid (grid, back, steps, choice )!, file_name)

    TYPE (domain)     , INTENT(INOUT) :: grid
    TYPE (domain)     , INTENT(IN)    :: back
    INTEGER           , INTENT(IN)    :: steps, choice
    !CHARACTER(LEN=128), INTENT(IN)    :: file_name

    LOGICAL                :: display, update
    REAL(KIND=8)           :: max_leng, tmp_max_leng
    INTEGER                :: i, nfro, nadd
    TYPE(simplex), POINTER :: s

    CALL adjust_spx

    CALL set_status (grid % ext_head, frz_status)
    !pone lo status di tutti gli ext pari a frz

    CALL jointop_simplex (grid % spx_head, grid % int_head)
    !pushtop per int su spx,perchè altrimenti spx_h vuota!(e svuota int)

    CALL set_status (grid % spx_head, und_status)
    !pone lo status di tutti gli spx pari a und

    !  Insert Steiner points
    CALL steiner_points

    WRITE (*,*) '   Delaunay advancing front...'
    !WRITE (*,*)
    !WRITE (*,'(3(1x,a7), 5x, a10)') 'STEP','NFRO','NADD', 'MAX_LENGTH'

    max_leng = 0.d0

    DO i = 1, steps

       tmp_max_leng = max_leng
       max_leng     = 0.d0

       CALL uns_split()

       nfro = 0

       CALL find_points(choice)

       CALL add_frn_points()

       !WRITE (*,'(3i8, 5x, f10.2)') i, nfro, nadd, max_leng
       IF ( nadd .EQ. 0 ) EXIT

    ENDDO

    !  Only for information
    ! eseguito solo se non sono stati esauriti tutti i spx
    ! cioè se esiste ancora un fronte che separa accettati da non accettati
    IF(numberof_simplex(grid % spx_head) /= 0) THEN
       display = .FALSE.
       CALL uns_split( display )
    ENDIF

    update =.TRUE.
    !update =.FALSE.

    WRITE(*,*)

    IF(update)THEN
       CALL crn_points

       CALL add_crn_points

       CALL restore_lists

!       WRITE(*,*)'   ...corner points splitted'
    ELSE
!       WRITE(*,*)'   ...corner points NOT splitted'
    ENDIF

    IF(update)THEN

       !  Information on lists (0 no messages - 1 essential - 2 verbose)
       CALL contability(grid, 0)

       !  Actualisation of boundary-internal elements relationships
       CALL merge_ie(grid % spx_head, grid % int_head, grid % ext_head)

       DO i = 1,grid % number_of_edges
          CALL flag_shell (nd_d, grid%edges(i)%int_head, grid%spx_head, &
               & grid%edges(i)%stc_head)

          s => grid % edges(i) % stc_head
          IF(.NOT. ASSOCIATED(s%next, s)) THEN
             WRITE(*,*)
             WRITE(*,*) 'ERROR', i
             WRITE(*,*) 'elements', numberof_simplex(grid%edges(i)%stc_head)
             WRITE(*,*) 'STOP'
             WRITE(*,*)
             STOP
          ENDIF
       ENDDO

       CALL split_ie(nd_d, grid % spx_head, grid % int_head, grid % ext_head)

       !  Information on lists (0 no messages - 1 essential - 2 verbose)
       CALL contability(grid, 0)
    ELSE
       WRITE(*,*)'   ...boundary <-> inner domain NOT updated'
    ENDIF

  CONTAINS
    !
    !------------------------------------------------------------------
    !
    SUBROUTINE adjust_spx

      !  Update of referenced points: only internal elements are used 
      !  for reference
      INTEGER                 :: i
      TYPE (simplex), POINTER :: s1, b1, simp
      TYPE (point)  , POINTER :: p1

      s1 => grid%int_head

      b1 => back%int_head

      spx_loop: DO

         s1 => s1%next

         b1 => b1%next
         IF(ASSOCIATED(s1, grid%int_head)) EXIT spx_loop !int_h vuota!

         !  Adjusting pnt%spx using only internal elements

         DO i = 1, nd_d+1
            p1   => s1%opp(i)%pnt !i punti degli interni(grid)
            simp => p1%spx        !il loro simplesso di appartenenza (s1)!
            IF(ASSOCIATED(simp)) THEN  !se p1 appartiene ad un simplesso e 
               IF(.NOT. ASSOCIATED(s1, simp)) THEN !se quel simplesso è diverso 
                  p1%spx => s1    !da s1 corrente
               ENDIF              !assegno p1 al simplesso corrente
            ELSE
               WRITE(*,*)'   ...pnt%spx not associated? (adjust_spx)'
            ENDIF

            !  Adjusting pnt%ppp%spx using only internal elements (of backgrid)

            p1   => s1%opp(i)%pnt%ppp
            IF(ASSOCIATED(p1)) THEN
               simp => p1%spx !simpl di appartenenza del pto 2d 
               IF(ASSOCIATED(simp)) THEN   !corrispondente a 1d
                  IF(.NOT. ASSOCIATED(b1, simp)) THEN
                     p1%spx => b1
                  ENDIF
               ELSE
                  WRITE(*,*)'   ...pnt%ppp%spx not associated? (adjust_spx)'
               ENDIF
            ELSE
               write(*,*)'ciao!(by DD:verso di percorrenza dei lati?)'
               stop
            ENDIF
         ENDDO

      ENDDO spx_loop

    END SUBROUTINE adjust_spx
    !
    !------------------------------------------------------------------
    !
    SUBROUTINE set_status (head, value)

      TYPE (simplex), TARGET, INTENT(IN) :: head
      INTEGER, INTENT(IN) :: value
      TYPE (simplex), POINTER :: s

      s  =>  head
      DO; s  =>  s % next; IF (ASSOCIATED(s,head)) EXIT
         s % status = value
      ENDDO

    END SUBROUTINE set_status
    !
    !------------------------------------------------------------------
    !
    SUBROUTINE uns_split( display )

      !  Accepts or rejects elements accordingly to choice:

      !by DD: "choice" dependence removed

      !  choice = 0	Euclidian length
      !  choice = 1	Riemann length


      LOGICAL     , INTENT(IN), OPTIONAL :: display

      LOGICAL :: unsatisf
      INTEGER :: i, j, n_spx, n_int!, k
      REAL (KIND=8) :: dt2, dp2
      TYPE (simplex), POINTER :: s, t
      TYPE (point)  , POINTER :: pi, pj
      ! compila la nuova lista degli interni
      ! con dei confronti di lunghezze

      max_leng = 0.d0

      !  Euclidean metric

      s  =>  grid % spx_head % next
      DO; IF ( ASSOCIATED (s, grid % spx_head) ) EXIT
         t  =>  s % next
         unsatisf = .FALSE.
         !controllo le tre lunghezze del simpl
         !parto da 1 e controllo 2 e 3
         !quando parto da 2 devo controllare solo 3
         loop1: DO i=1,nd_d !1,2
            pi  =>  s % opp(i) % pnt 
            DO j=i+1,nd_d+1 !i=1->2,3  i=2->3,3
               pj  =>  s % opp(j) % pnt

               dp2 = segment_length (nd_d, pi, pj)!distanza^2 
               dt2 = (e_tol*0.5d0*(pi % leng + pj % leng))**2
               !0.75
               IF ( dp2 > dt2 ) THEN
                  max_leng = MAX(max_leng,dp2)
                  unsatisf = .TRUE.
                  EXIT loop1
               ENDIF

            ENDDO

         ENDDO loop1

         IF ( .NOT. unsatisf ) THEN !unsatisf =FALSE
            s % status = int_status
            CALL movetop_simplex (grid % int_head, s)
            !compila la lista dei simplessi accettati
            !relativamente al passo precedente
            !gli altri rimangono ancora in spx_h
         ENDIF

         s  =>  t
      ENDDO



!      IF(PRESENT(display)) THEN
!         n_spx = numberof_simplex(grid % spx_head)
!         n_int = numberof_simplex(grid % int_head)
!         IF(n_spx /= 0) THEN
!            WRITE(*,50)n_spx,'->',100.d0*DBLE(n_spx)/DBLE(n_spx+n_int), '%'
!            WRITE(*,*)
!         ENDIF                           !DBLE converte in  reale doppia precis.
!      ENDIF

      ! spx lista dei simplessi non accettati
      ! totale = non accettati + int

50    FORMAT('    *** not accepted', 1x, i8, 1x, a2, 1x, f5.2, 1x, a1)

    END SUBROUTINE uns_split
    !
    !------------------------------------------------------------------
    !
    SUBROUTINE find_points(choice)

      !  Rules for positioning internal points accordingly to choice:

      !  choice = 0	Euclidian length
      !  choice = 1	Riemann length

      INTEGER, INTENT(IN) :: choice

      SELECT CASE(choice)

      CASE(0)	! Euclidean metric
         CALL frn_points


      CASE DEFAULT 
         WRITE(*,*)'BAD OPTION: no construct with choice', choice
         STOP

      END SELECT

    END SUBROUTINE find_points
    !
    !------------------------------------------------------------------
    !
    SUBROUTINE frn_points

      !  Euclidean metric only

      LOGICAL :: err_status
      INTEGER :: m, i, j, k
      TYPE (simplex), POINTER :: s, s_back
      TYPE (point)  , POINTER :: p
      REAL (KIND=8) :: fn, srad2, leng2, num, den, frad2, heig2, dist2, t
      REAL (KIND=8) :: sum
      REAL (KIND=8), DIMENSION(nd_d) :: center

      m = 0
      fn = 1.d0/DBLE(nd_d)   !=0.5
      s  =>  grid % spx_head !lista dei simplessi non accettati

      spx_loop: DO

         s  =>  s % next
         IF ( ASSOCIATED (s, grid % spx_head) ) EXIT spx_loop
         CALL simplex_center (nd_d, s, center, srad2)! centro e raggio ccc

         sum = 0.0d0
         DO j= 1, nd_d+1
            sum = sum + s % opp(j) % pnt % leng
         ENDDO

         length_loop: DO i=1,nd_d+1

            IF ( s % opp(i) % spx % status /= und_status ) THEN

               !  Check if internal to the face(!?)
               num = s % fmat(i,nd_d+1)
               den = 0.0d0
               DO j=1,nd_d
                  num = num + s % fmat(i,j) * center(j)
                  den = den + s % fmat(i,j) ** 2
               ENDDO
               ! num  proporzionale ad una distanza
               IF ( num <= 0.d0 ) CYCLE !il centro del triangolo è esterno  
               !al triangolo stesso 
               !?distante o molto distante?  
               leng2 = (fn*(sum - s % opp(i) % pnt % leng))**2

               heig2 = num**2/den 
               frad2 = srad2-heig2
               dist2 = leng2-frad2

               t = SQRT(heig2/den)-SQRT(dist2/den)
               t = MAX (0.d0,t)

               IF (t >= 0.d0) THEN

                  m = m+1 !numero di front steps
                  k = 1+MOD(i,nd_d+1)
                  p => grid%gar_head%next

                  IF ( ASSOCIATED (p, grid % gar_head) ) THEN
                     p  =>  new_point()                  !lista vuota
                     CALL pushtop_point (grid % add_head, p)
                  ELSE                                   !lista esistente
                     CALL movetop_point (grid % add_head, p)
                  ENDIF

                  DO j=1,nd_d
                     p % x(j) = center(j) - s % fmat(i,j) * t
                     ! new point coordinates 
                     ! p % metric(j,j) = 0.00001d0/leng2
                  ENDDO

                  p % old  =>  s % opp(k) % pnt

                  CALL eval_back (p, back%int_head, s_back, err_status)

                  IF(err_status) THEN
                     !punto non rientrante nella backgrid
                     !WRITE(*,*)'   ...Research failed! Point not inserted'
                     !WRITE(*,*)'Point will not be inserted (frn_points)'
                     CALL movetop_point (grid % gar_head, p)
                     m = m-1
                     cycle spx_loop !passa al simplesso successivo
                  ENDIF

               ENDIF
            ENDIF
         ENDDO length_loop
      ENDDO spx_loop

      nfro = m

    END SUBROUTINE frn_points


    !
    !------------------------------------------------------------------
    !


    SUBROUTINE add_frn_points()

      !  Insert points accordingly to choice:

      !by DD: "choice" dependence removed 

      LOGICAL :: far
      INTEGER :: m, i
      TYPE (point), POINTER :: p, q

      i = 0
      m = 0
      p  =>  grid % add_head % next
      !


      !  Euclidean length

      DO
         IF ( ASSOCIATED (p, grid % add_head) ) EXIT
         q  =>  p % next

         CALL insert_far_point (nd_d, grid % spx_head, p, far)

         IF (far) THEN
            m = m+1
            CALL movetop_point (grid % pnt_head, p)

         ELSE

            CALL movetop_point (grid % gar_head, p)

         ENDIF

         p => q

      END DO

      nadd = m


    END SUBROUTINE add_frn_points
    !
    !------------------------------------------------------------------
    !
    SUBROUTINE crn_points

      INTEGER :: i, j, m
      TYPE (simplex), POINTER :: s, n


      m = 0; s  =>  grid % int_head
      !ricerca nella lista dei simplessi accettati
      DO; s  =>  s % next; IF ( ASSOCIATED (s, grid % int_head) ) EXIT

         j = 0

         DO i=1,nd_d+1
            n  =>  s % opp(i) % spx
            IF ( n % status == frz_status ) j = j+1
            ! conta i simplessi esterni della guaina dei contorni(backgrid)
            ! da ext_status erano stati posti a frz_status
         ENDDO


         IF ( j == nd_d ) THEN
            DO i=1,nd_d+1
               n  =>  s % opp(i) % spx
               IF ( ( n % status == int_status ) ) THEN
                  s % status = stc_status
                  m = m+1
                  CALL halve_dom (grid, i, s)
                  !mette i punti da inserire in grid%add_head
                  EXIT
               ENDIF
            ENDDO
         ELSEIF ( j == nd_d+1 ) THEN
            WRITE (*,*) 'refine_corners error, unpossible element found'
            WRITE (*,*) '(stop)'
            STOP
         ENDIF
      ENDDO

      s  =>  grid % int_head
      DO; s  =>  s % next; IF ( ASSOCIATED (s, grid % int_head) ) EXIT
         s % status = int_status
      ENDDO

      nadd = m

    END SUBROUTINE crn_points
    !
    !------------------------------------------------------------------
    !
    SUBROUTINE add_crn_points

      INTEGER :: m
      TYPE (point), POINTER :: p, q


      m = 0; p  =>  grid % add_head % next

      DO; IF ( ASSOCIATED (p, grid % add_head) ) EXIT
         q  =>  p % next
         m = m+1
         !write(*,*) 'trying to insert p',p % x
         CALL insert_point (nd_d, grid % int_head, p)
         CALL movetop_point (grid % pnt_head, p)
         write(*,*)'   ...corner point added', p % x
         p  =>  q
      END DO

      nadd = m

    END SUBROUTINE add_crn_points
    !
    !------------------------------------------------------------------
    !
    SUBROUTINE restore_lists

      TYPE (simplex), POINTER :: s

      s  =>  grid % spx_head
      DO; s  =>  s % next; IF ( ASSOCIATED (s, grid % spx_head) ) EXIT
         s % status = int_status
      ENDDO

      s  =>  grid % ext_head
      DO; s  =>  s % next; IF ( ASSOCIATED (s, grid % ext_head) ) EXIT
         s % status = ext_status
      ENDDO

      IF ( .NOT. ASSOCIATED (grid % spx_head % next, grid % spx_head) ) THEN
         CALL jointop_simplex (grid % int_head, grid % spx_head)
      ENDIF

    END SUBROUTINE restore_lists
    !
    !------------------------------------------------------------------
    !
    SUBROUTINE steiner_points
    use init_data, only: nsteinr, steinrP
      !  Additional grid points

      !CHARACTER(LEN=128), INTENT(IN) :: file_name

      LOGICAL  :: err_status
      INTEGER  :: i, j, n
      TYPE(point), POINTER :: p

      TYPE(point), TARGET    :: p_help 
      TYPE(simplex), POINTER :: s_back

      IF ( nsteinr > 0 )THEN
         WRITE(*,*)'Steiner points...'
      ENDIF
      DO i = 1, nsteinr

         p => grid%gar_head%next
         IF ( ASSOCIATED (p, grid % gar_head) ) THEN
            p  =>  new_point()
            CALL pushtop_point (grid % add_head, p)
         ELSE
            CALL movetop_point (grid % add_head, p)
         ENDIF

         p%x(1:nd_d) = steinrP(:,i)
         p%old => grid%pnt_head%next

         !by DD begin
         !  Find the backgrid simplex containing pnt
         !  Interpolation of metric of pnt over his backgrid 
         !  simplex

         err_status = .FALSE.

         CALL eval_back(p, back%int_head, s_back, err_status)

         IF( .not. err_status) THEN
            CALL int_metr_tria(s_back, p, p_help)
            ! Interpolate metric on back element
         ENDIF
         !by DD end

         IF(err_status) THEN

            WRITE(*,*)'    *** Steiner point NOT found!'
            CALL movetop_point (grid % gar_head, p)

         ELSE

            CALL insert_point (nd_d, grid % spx_head, p)
            CALL movetop_point (grid % pnt_head, p)

         ENDIF

      ENDDO

      RETURN

      !WRITE(*,*)'exit steiner_points'

1000  WRITE(*,*)' ** STEINER FILE DOES''NT EXIST **'
      RETURN
1100  WRITE(*,*)'UNEXPECTED END OF STEINER FILE'
      RETURN

    END SUBROUTINE steiner_points
    !
    !------------------------------------------------------------------
    !
  END SUBROUTINE front_grid
  !
  !********************************************************************
  !
END MODULE front_2d
