MODULE metric_2d
  
  USE curves          !  Curves module

  USE grid_types      !  Types module

  USE init_data       !  Boundary input module

  USE delaunay        !  Delaunay basic module

  USE eucl_delaunay   !  Delaunay Euclide module
    
  USE eigsys_2d       !  Eigenvalues & Eigenvectors module
  

  IMPLICIT NONE

  ! PARAMETERS

  REAL(kind=8), PARAMETER, PUBLIC :: pi = 3.1415927464102

  REAL(kind=8), PARAMETER :: tiny    = 1.d-8
  REAL(kind=8), PARAMETER :: epsilon = 1.d-12
  REAL(kind=8), PARAMETER :: lin_lim  = 2.d0
  REAL(kind=8), PARAMETER :: sin_lim  = 10.d0

  PUBLIC :: gen_metr, proj_metr, attr_metr, m_plot,  &
       & 	  int_metr_tria ,projection!, int_metr_seg

CONTAINS
  !
!********************************************************************
  !
  SUBROUTINE gen_versor(c, s0, vers)

    INTEGER                     , INTENT(IN)  :: c
    REAL(kind=8), DIMENSION(:)  , INTENT(IN)  :: s0
    REAL(kind=8), DIMENSION(:)  , INTENT(OUT) :: vers

    INTEGER      :: j
    REAL(KIND=8) :: mod_t

    CALL evalc_1(c, s0, vers,curv)

    mod_t = zero
    DO j = 1, nd_d
       mod_t = mod_t+vers(j)**2
    ENDDO
    mod_t = SQRT(mod_t)
    mod_t = one/mod_t
    vers  = vers*mod_t
    
  END SUBROUTINE gen_versor
  !
!********************************************************************
  !
  SUBROUTINE set_teta(a, t)

    REAL(kind=8)              , INTENT(OUT) :: a
    REAL(kind=8), DIMENSION(:), INTENT(IN)  :: t

    REAL(KIND=8) :: teta

    teta = ASIN(t(2))
    
    IF(t(1) >= zero) THEN

       IF(t(2) >= zero) THEN
          a    = teta
       ELSE
          a    = 2.d0*pi+teta
       ENDIF
    ELSE
       IF(t(2) >= zero) THEN
          a    = pi-teta
       ELSE
          a    = pi-teta
       ENDIF
    ENDIF

  END SUBROUTINE set_teta
  !
!********************************************************************
  !
  SUBROUTINE gen_metr(a, ellipse, mat)

    REAL(kind=8), DIMENSION(:)  , INTENT(IN)  :: a
    REAL(kind=8), DIMENSION(:)  , INTENT(OUT) :: ellipse 
    REAL(kind=8), DIMENSION(:,:), INTENT(OUT) :: mat 

    REAL(kind=8)                       :: teta
    REAL(kind=8), DIMENSION(nd_d,nd_d) :: p, l

    ellipse = a
    teta = ellipse(3)

    p(1,1) = COS(teta)
    p(2,1) = SIN(teta)
    p(1,2) = -p(2,1)
    p(2,2) =  p(1,1)

    l(1,1) = one/(ellipse(1)**2)
    l(2,2) = one/(ellipse(2)**2)
    l(1,2) = zero
    l(2,1) = zero

    mat = TRANSPOSE(p)
    mat = MATMUL(l,mat)
    mat = MATMUL(p,mat)

  END SUBROUTINE gen_metr
  !
!********************************************************************
  !
  SUBROUTINE proj_metr(c, edge_kind, s, a, h)

    INTEGER                   , INTENT(IN)  :: c, edge_kind
    REAL(kind=8)              , INTENT(IN)  :: s
    REAL(kind=8), DIMENSION(:), INTENT(IN)  :: a
    REAL(kind=8)              , INTENT(OUT) :: h

    INTEGER :: j, k
    REAL(kind=8), DIMENSION(nd_d-1)    :: s0
    REAL(kind=8), DIMENSION(nd_d)      :: t
    REAL(kind=8), DIMENSION(2*nd_d-1)  :: ell
    REAL(kind=8), DIMENSION(nd_d,nd_d) :: mat

    ! edge_kind = 0 isotropic metric
    ! edge_kind = 1 anisotropic metric with given angle
    ! edge_kind = 2 anisotropic metric aligned

    IF(edge_kind /= 1) THEN
       h = a(1)
    ELSE
       s0(1) = s

       CALL gen_versor(c, s0, t)

       CALL gen_metr(a, ell, mat)

       h = projection(mat, t)
       h = zero
       DO j = 1, nd_d
          DO k = 1, nd_d
             h = h+t(j)*mat(j,k)*t(k)
          ENDDO
       ENDDO
       h = one/SQRT(h)
    ENDIF

  END SUBROUTINE proj_metr
  !
!********************************************************************
  !
  FUNCTION projection(mat, t) RESULT(hh)

    REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: mat
    REAL(kind=8), DIMENSION(:)  , INTENT(IN) :: t
    REAL(kind=8) :: hh

    INTEGER :: j, k

    hh = zero

    DO j = 1, nd_d
       DO k = 1, nd_d
          hh = hh+t(j)*mat(j,k)*t(k)
          
       ENDDO
    ENDDO
 
    hh = one/SQRT(hh)

  END FUNCTION projection
  !
!********************************************************************
  !
  SUBROUTINE attr_metr(grid, i, pnt, i_l, i_r, pnt_2)

    INTEGER                   , INTENT(in) :: i, i_l, i_r
    TYPE(domain), TARGET      , INTENT(in) :: grid
    TYPE(point)               , POINTER    :: pnt,pnt_2

    INTEGER     :: edge_kind
    REAL(kind=8):: s1, s2, s_s
    REAL(kind=8), DIMENSION(nd_d)      :: t
    REAL(kind=8), DIMENSION(2*nd_d-1)  :: a1, a2, a

    ! edge_kind = 0 isotropic metric
    ! edge_kind = 1 anisotropic metric with given angle
    ! edge_kind = 2 anisotropic metric aligned

    ! position  = 0 external edge (clockwise)
    ! position  = 1 internal edge (counterclockwise)

    edge_kind = grid%edges(i)%edge_kind
    a1        = grid%edges(i)%hdata(:,i_l)! hdata transposed!hdatax,y,th in column
    s1        = grid%edges(i)%sdata(i_l)  !first point chord abscissa 

    IF(i_l == i_r) THEN

       IF(edge_kind == 2) THEN
          CALL gen_versor(i, pnt%x, t)
          IF(grid%edges(i)%position == 1) THEN
             t= -t
          ENDIF
          CALL set_teta(a1(3), t)
       ENDIF
       CALL gen_metr(a1, pnt_2%ellipse, pnt_2%metric) 

    ELSE

       a2 = grid%edges(i)%hdata(:,i_r)
       s2 = grid%edges(i)%sdata(i_r)

       s_s = (pnt%x(1)-s1)/(s2-s1)

       a = a1+s_s*(a2-a1)

       IF(edge_kind == 2) THEN
          call gen_versor(i, pnt%x, t)
          IF(grid%edges(i)%position == 1) THEN
             t= -t
          ENDIF
          call set_teta(a(3), t)
       ENDIF
       CALL gen_metr(a, pnt_2%ellipse, pnt_2%metric)

    ENDIF

  END SUBROUTINE attr_metr
  !
!********************************************************************
  !

  SUBROUTINE int_metr_tria(simp, pnt, p_star)

    TYPE (simplex), POINTER :: simp
    TYPE (point)  , TARGET  :: pnt, p_star

    INTEGER                         :: i
    REAL(KIND=8)                    :: csi
    REAL(KIND=8), DIMENSION(nd_d+1) :: p_dist
    REAL(KIND=8), DIMENSION(nd_d)   :: p2_p1, p3_p1, p_p3, p_p1
    TYPE(point) , POINTER           :: pnt1, pnt2, pnt3

    !  Test if point located on element face

    p_dist = zero
    CALL eval_facets(nd_d, pnt, simp, p_dist)

    DO i = 1, nd_d+1

       IF(p_dist(i) <= epsilon) THEN

         pnt1 => simp%opp(MOD(i, nd_d+1)+1)  %pnt
          pnt2 => simp%opp(MOD(i+1, nd_d+1)+1)%pnt

          p_p1  = pnt%x-pnt1%x
          p2_p1 = pnt2%x-pnt1%x

          IF( ABS(p2_p1(1)) > ABS(p2_p1(2)) ) THEN
             csi = p_p1(1)/p2_p1(1)
          ELSE
             csi = p_p1(2)/p2_p1(2)
          ENDIF

          ! Interpolate metric in p along pnt1-pnt2
          !by DD
          CALL int_metr_seg2(pnt1%metric, pnt2%metric, pnt%metric, csi)
           !pnt%metric=0.0
           
          RETURN

       ENDIF

    ENDDO

    pnt1 => simp%opp(1)%pnt
    pnt2 => simp%opp(2)%pnt
    pnt3 => simp%opp(3)%pnt

    !  Calculation of p_star, intersection of p1-p2 and p3-pnt

    p2_p1 = pnt2%x-pnt1%x
    p3_p1 = pnt3%x-pnt1%x
    p_p3  = pnt%x-pnt3%x

    csi = p2_p1(2)*p_p3(1)-p2_p1(1)*p_p3(2)
    IF(csi < epsilon) THEN
       WRITE(*,*)' csi1 denominator is singular - STOP- (int_metr_tria)'
       STOP
    ELSE
       csi = (p3_p1(2)*p_p3(1)-p3_p1(1)*p_p3(2))/csi
    ENDIF

    p_star%x = pnt1%x+csi*p2_p1

    ! Interpolate metric in p_star, along pnt1-pnt2
    !by DD
    CALL int_metr_seg2(pnt1%metric, pnt2%metric, p_star%metric, csi)
    !p_star%metric=0.0

    ! Interpolate metric in pnt, along p_star,pnt3

    p_p1  = pnt%x-p_star%x
    p2_p1 = pnt3%x-p_star%x
    
    IF( ABS(p2_p1(1)) > ABS(p2_p1(2)) ) THEN
       csi = p_p1(1)/p2_p1(1)
    ELSE
       csi = p_p1(2)/p2_p1(2)
    ENDIF
    !by DD
    CALL int_metr_seg2(p_star%metric, pnt3%metric, pnt%metric, csi)
    !pnt%metric=0.0
  
  END SUBROUTINE int_metr_tria
  !
!********************************************************************
  !
  
  SUBROUTINE int_metr_seg2(metric1, metric2, metric, csi)
    
    !  Interpolate metrics with simultaneous reduction; 
    !  linear // geometric // sinusoidal interpolation available
    !  with automatic switching for best results

    REAL (KIND=8)                , INTENT(IN)  :: csi
    REAL (KIND=8), DIMENSION(:,:), INTENT(IN)  :: metric1, metric2
    REAL (KIND=8), DIMENSION(:,:), INTENT(OUT) :: metric

    LOGICAL :: err, exit_status
    INTEGER :: i, j, k

    REAL(kind=8)                       :: c1, c
    REAL(KIND=8), DIMENSION(nd_d)      :: lambda1, lambda2, lambda
    REAL(KIND=8), DIMENSION(nd_d,nd_d) :: alambda, p1, p2, p, inv_p

    CALL invert_2d(metric1, p, err)

    IF (err) THEN
       WRITE(*,*)'metric0 is not invertible!'
       WRITE(*,*)'-- STOP -- (int_metr_seg)'
       STOP
    ELSE
       alambda = MATMUL(p, metric2)
    ENDIF

    p      = zero
    lambda = zero
    CALL eigensys_2d(alambda, lambda, p, exit_status)

    CALL ortogon(metric1, metric2, p)
    IF (exit_status) THEN  

 
       lambda1 = zero
       lambda2 = zero
       DO i = 1, nd_d
          DO j = 1, nd_d
             DO k = 1, nd_d
                lambda1(i) = lambda1(i)+p(j,i)*metric1(j,k)*p(k,i)
                lambda2(i) = lambda2(i)+p(j,i)*metric2(j,k)*p(k,i)
             ENDDO
          ENDDO
          lambda1(i) = one/SQRT(lambda1(i))
          lambda2(i) = one/SQRT(lambda2(i))
       ENDDO


       alambda = zero


       !  Linear interpolation

       alambda(1,1) = lambda1(1)+csi*(lambda2(1)-lambda1(1))
       alambda(2,2) = lambda1(2)+csi*(lambda2(2)-lambda1(2))



       alambda(1,1) = one/(alambda(1,1)**2)
       alambda(2,2) = one/(alambda(2,2)**2)

       CALL invert_2d(p, inv_p, err)

       IF (err) THEN
          WRITE(*,*)'eigenvectors is a Singular MATRIX'
          WRITE(*,*)p(1,1), p(1,2)
          WRITE(*,*)p(2,1), p(2,2)
          WRITE(*,*)'-- STOP -- (int_metr_seg)'
          STOP
       ENDIF

       p      = MATMUL(alambda, inv_p) 
       metric = MATMUL(TRANSPOSE(inv_p),p)


    ELSE

       !  metric a1 and a2 are proportional

       CALL eigensys_2d(metric1, lambda1, p1, exit_status)
       CALL eigensys_2d(metric2, lambda2, p2, exit_status)

       c1 = lambda1(1)/lambda2(1)
       c  = ABS(one-c1)

       IF(c <= 1.d-8) THEN
          metric = metric1
          RETURN
       ENDIF

       DO i = 1, nd_d
          lambda1(i) = one/SQRT(lambda1(i))
          lambda2(i) = one/SQRT(lambda2(i))
       ENDDO

       alambda = zero
       alambda(1,1) = lambda1(1)+csi*(lambda2(1)-lambda1(1))
       alambda(2,2) = lambda1(2)+csi*(lambda2(2)-lambda1(2))

       alambda(1,1) = one/(alambda(1,1)**2)
       alambda(2,2) = one/(alambda(2,2)**2)

       inv_p  = TRANSPOSE(p1)
       p      = MATMUL(alambda, inv_p) 
       metric = MATMUL(p1,p) 

    ENDIF

  END SUBROUTINE int_metr_seg2
  !
!********************************************************************
  !
  SUBROUTINE ortogon(metric1, metric2, p)
     
    !  Ortogonal eigenvectors to M1 an M2

    REAL (KIND=8), DIMENSION(:,:), INTENT(IN)    :: metric1, metric2
    REAL (KIND=8), DIMENSION(:,:), INTENT(INOUT) :: p
    INTEGER :: i, j

    REAL(KIND=8)                       ::  mod_p!min_diag
    REAL(KIND=8), DIMENSION(nd_d)      :: a, b
    REAL(KIND=8), DIMENSION(nd_d,nd_d) :: alambda


    DO i = 1, nd_d
       b(i) = zero
       DO j = 1, nd_d
          b(i) = b(i)+metric1(i,j)*p(j,1)
       ENDDO
    ENDDO

    DO i = 1, nd_d
       a(i) = zero
       DO j = 1, nd_d
          a(i) = a(i)+p(j,i)*b(j)
       ENDDO
    ENDDO

    IF(ABS(a(1)) > ABS(a(2))) THEN
       p(:,2) = p(:,2)-a(2)/a(1)*p(:,1)
       mod_p = SQRT(one+(a(2)/a(1))**2)
       p(:,2) = p(:,2)/mod_p
    ELSE!
       p(:,2) = p(:,1)-a(1)/a(2)*p(:,2)
       mod_p = SQRT(one+(a(1)/a(2))**2)
       p(:,2) = p(:,2)/mod_p
    ENDIF

    alambda = MATMUL(TRANSPOSE(p), metric1)
    alambda = MATMUL(alambda, p)
    
    alambda = MATMUL(TRANSPOSE(p), metric2)
    alambda = MATMUL(alambda, p)

  END SUBROUTINE ortogon
  !
!********************************************************************
  !
  SUBROUTINE m_plot (idf, r_int_head, iswitch)
    
    !by DD:grid and r_ext_head dependence removed  
    INTEGER, INTENT(IN) :: idf, iswitch
    
    TYPE (simplex), POINTER :: r_int_head 

    INTEGER :: i, i_x, i_y
    TYPE (point), POINTER :: p
    TYPE (simplex), POINTER :: s_head, s
    CHARACTER(len = 40) :: title
    CHARACTER(len = 1)  :: who1, who2
   
    !  For ellipse calculations
    LOGICAL                            :: exit_status
    REAL(KIND=8), DIMENSION(nd_d)      :: lambda
    REAL(KIND=8), DIMENSION(nd_d,nd_d) :: p_vec
    IF(iswitch == 0) THEN  !  Metric plot

       title = '% toplabel = Metric('   
       DO i_x = 1, nd_d
          who1(1:1) = CHAR(48+i_x)

          DO i_y = 1, nd_d

             !IF(MOD(i_x+i_y, nd_d) /= 0) CYCLE ! to plot diagonal terms only
             IF(i_x == nd_d .AND. i_y == 1) CYCLE ! to plot only one extradiag term

             who2(1:1)= CHAR(48+i_y)
             title = title(1:20)//who1//','//who2//')'

             WRITE (idf,*) '$ DATA = CONTCURVE'
             WRITE (idf,*) title
             WRITE (idf,*) '% contstyle = 2'
             WRITE (idf,*) '% cstep = 20'
             WRITE (idf,*) '% meshplot  = True'

             s_head  =>  r_int_head
             s  =>  s_head
             DO
                s  =>  s % next
                IF ( ASSOCIATED (s, s_head) ) EXIT
                WRITE (idf,*)
                DO i=1,nd_d+1
                   p => s%opp(i)%pnt
                   WRITE (idf,'(3e15.7,i8)') p%x, p%metric(i_x,i_y), p%index
                ENDDO
             ENDDO

          ENDDO
       ENDDO

    ELSE IF(iswitch == 1) THEN   ! Hdata plot

       title = '% toplabel = "Ellipse semiaxis('   
       DO i_x = 1, nd_d
          who1(1:1) = CHAR(48+i_x)

          title = title(1:31)//who1//')"'

          WRITE (idf,*) '$ DATA = CONTCURVE'
          WRITE (idf,*) title
          WRITE (idf,*) '% contstyle = 2'
          WRITE (idf,*) '% cstep = 20'
          WRITE (idf,*) '% meshplot  = True'

          s_head  =>  r_int_head
          s  =>  s_head
          DO
             s  =>  s % next
             IF ( ASSOCIATED (s, s_head) ) EXIT
             WRITE (idf,*)
             DO i=1,nd_d+1
                p => s%opp(i)%pnt

                !  Ellipse calculation

                IF(p%metric(1,1) /= zero) THEN

                   p_vec      = zero
                   lambda = zero
                   CALL eigensys_2d(p%metric, lambda, p_vec, exit_status)
                   p%ellipse(1) = one/SQRT(lambda(1))
                   p%ellipse(2) = one/SQRT(lambda(2))
                   CALL set_teta(p%ellipse(3), p_vec(:,1))

                ENDIF

                WRITE (idf,'(3e15.7,i8)') p%x, p%ellipse(i_x), p%index
             ENDDO
          ENDDO

       ENDDO

       title = '% toplabel = "Teta [degrees]"'   

       WRITE (idf,*) '$ DATA = CONTCURVE'
       WRITE (idf,*) title
       WRITE (idf,*) '% contstyle = 2'
       WRITE (idf,*) '% cstep = 20'
       WRITE (idf,*) '% meshplot  = True'

       s_head  =>  r_int_head
       s  =>  s_head
       DO
          s  =>  s % next
          IF ( ASSOCIATED (s, s_head) ) EXIT
          WRITE (idf,*)
          DO i=1,nd_d+1
             p => s%opp(i)%pnt
             WRITE (idf,'(3e15.7,i8)') p%x, p%ellipse(nd_d+1)*180.d0/pi, p%index
          ENDDO
       ENDDO

    ENDIF

  END SUBROUTINE m_plot
  !
!********************************************************************
  !

END MODULE metric_2d
