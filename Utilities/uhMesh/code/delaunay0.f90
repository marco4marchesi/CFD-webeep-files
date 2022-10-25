MODULE delaunay

  USE linsys      !  Direct Linear system solver

  USE grid_types  !  Grid types

  USE list        !  List operations

  IMPLICIT NONE

  !               ###  G L O B A L   V A R I A B L E S  ###

  REAL (KIND=8), PARAMETER, PUBLIC :: neg_tol = -1.d-12

  TYPE(simplex), POINTER  , PUBLIC :: bbl_head, del_head, new_head, gar_head

 
  INTEGER      , DIMENSION(maxs)     , PRIVATE :: pivot
  REAL (KIND=8), DIMENSION(maxs)     , PRIVATE :: vector
  REAL (KIND=8), DIMENSION(maxs,maxs), PRIVATE :: matrix

  !                  ###  R O U T I N E   I N D E X  ###

  !  1) GENERAL:

  PUBLIC  :: delaunay_setup, walking_search, trivial_search,  &
       & simplex_facets, simplex_normal, eval_facets

  !  2) BUBBLE AND SPLIT:

  PUBLIC  :: flag_shell, split_ie_1d, split_ie, merge_ie
  
  PRIVATE :: build_bubble, match_bubble_ndm1, swap_ie_ndm1, &
  &          match_bubble_nd, swap_ie_nd

  !  3) INSERTION:

  PUBLIC  :: set_delaunay, delete_first, update_geometry, update_topology
  
  PRIVATE :: delete_basis
!
!
      contains
!
!
      subroutine delaunay_setup
!
      bbl_head  =>  newhead_simplex ()
      del_head  =>  newhead_simplex ()
      new_head  =>  newhead_simplex ()
      gar_head  =>  newhead_simplex ()
!
      return
      end subroutine delaunay_setup
!
!
!
  SUBROUTINE walking_search (nd, pnew, simp, status)

    ! +------------------------------------------+
    ! |   SEARCH FOR A SIMPLEX CONTAINING PNEW   |
    ! +------------------------------------------+

    !  status in input determines where the search has to be performed: 
    !  0 -> backgrid
    !  1 -> maingrid

    !  status in output indicates the success or failure of the search:
    !  -1 -> point not contained in grid
    !   0 -> point not located (exceeded max_it)
    !   1 -> point located

    !USE predicates   !! Adaptive arithmetic? See later on

    INTEGER                , INTENT(IN)    :: nd
    TYPE (point)  , TARGET , INTENT(IN)    :: pnew
    TYPE (simplex), POINTER                :: simp
    INTEGER                , INTENT(INOUT) :: status

    INTEGER, PARAMETER :: max_it = 100000

    LOGICAL :: in_spx
    INTEGER :: i, n, m, avoid_status
    REAL (KIND=4):: h
    TYPE (simplex), POINTER :: s

    s  =>  simp
    NULLIFY (simp)

    IF(status == 0) THEN
       ! => back grid search
       avoid_status = ext_status
    ELSE IF(status == 1) THEN
       !=> main grid search
       avoid_status = frz_status
    ENDIF

    IF(s%status == avoid_status) THEN
       WRITE(*,*)'research is compromised. STOP (walking_search)'
       STOP
    ENDIF

    walk: DO m = 1, max_it

       CALL eval_facets (nd, pnew, s, vector)

       n = 0
       in_spx = .TRUE.
       DO i=1,nd+1
          IF (vector(i) <= neg_tol) THEN
             
             in_spx = .FALSE.

             IF(s%opp(i)%spx%status /= avoid_status) THEN
                n = n+1
                pivot(n) = i
             ENDIF

          ENDIF
       ENDDO

       IF (n == 0 ) THEN

          IF(in_spx) THEN

             !  Simplex found
             
             simp => s
             status = 1
             RETURN
             
          ELSE
             
             ! point lies on a forbidden simplex -> negative result

             status = -1
             RETURN

          ENDIF

       ELSE

          CALL RANDOM_NUMBER (h)
          s => s % opp(pivot(1+INT(n*h))) % spx

       ENDIF

    ENDDO walk

    !  Research failed after max_it steps

    status = 0
    WRITE(911,*)'Warning: walking search failed after', max_it, 'steps'



  END SUBROUTINE walking_search
  !
!********************************************************************
  !
  SUBROUTINE trivial_search (nd, spx_head, pnew, simp, err_status)

    ! +------------------------------------------+
    ! |   SEARCH FOR A SIMPLEX CONTAINING PNEW   |
    ! +------------------------------------------+

    INTEGER, INTENT(IN) :: nd
    TYPE (simplex), TARGET, INTENT(IN) :: spx_head
    TYPE (point), TARGET, INTENT(IN) :: pnew
    TYPE (simplex), POINTER :: simp
    LOGICAL, INTENT(OUT) :: err_status

    INTEGER :: i, n
    REAL (KIND=8) :: c_maxmin, c_min
    TYPE (simplex), POINTER :: s

    err_status = .TRUE.
    NULLIFY (simp); c_maxmin = -HUGE(c_maxmin)

    s  =>  spx_head
    DO; s  =>  s % next; IF ( ASSOCIATED(s, spx_head) ) EXIT

       IF ( s % status > 0 ) THEN
          CALL eval_facets (nd, pnew, s, vector)

          n = 0; c_min = +HUGE(c_min)
          DO i=1,nd+1
             IF ( vector(i) >= neg_tol ) n = n+1
             c_min = MIN(c_min, vector(i))
          ENDDO

          IF (n == nd+1) THEN
             simp  =>  s
             err_status = .FALSE.!point located by trivial search
             EXIT !DO
          ENDIF
         

          IF (c_maxmin < c_min) THEN
             c_maxmin  =  c_min
          ENDIF
       ENDIF

    ENDDO


  END SUBROUTINE trivial_search
  !
!********************************************************************
  !


  ! ### GENERAL ROUTINES (GEOMETRICAL)


  SUBROUTINE simplex_facets (nd, simp)

    ! +---------------------------------------------------------+
    ! |   FACET NORMAL VECTORS (BARYCENTRIC COORDINATE MATRIX)  |
    ! +---------------------------------------------------------+

    ! The barycentric coordinate matrix contains the nd+1
    ! coefficient of the facet nd-planes.

    ! 3D example: the equation of a 3-plane is ax + by + cz + d = 0,
    ! (a, b, c, d defined up to a constant value).
    
    ! The coefficients for face k are stored in simp % fmat's k-row
    ! face k opposite to point k(local number) 

    INTEGER, INTENT(IN) :: nd
    TYPE (simplex), TARGET, INTENT(INOUT) :: simp
    INTEGER :: i, j


    simp % fmat = 0.

    DO j=1,nd+1
       simp % fmat(j,j) = 1.
       DO i=1,nd
          matrix(i,j) = simp % opp(j) % pnt % x(i)
       ENDDO
       matrix(nd+1,j) = 1.
    ENDDO

    
    CALL lufact (nd+1, pivot, matrix, simp % det)
    
    
    
    DO j=1,nd+1
       CALL lusolv (nd+1, pivot, matrix, simp % fmat(:,j))
       !fmat:input by columns,output by columns
       !     but coefficients by rows!  
    ENDDO

  END SUBROUTINE simplex_facets
  !
!********************************************************************
  !
  SUBROUTINE simplex_normal (nd, simp_ndm1, normal)

    ! +---------------------------+
    ! |   SIMPLEX NORMAL VECTOR   |
    ! +---------------------------+

    ! WARNING: simp_ndm1 is a simplex in nd-1 dimensional space

    INTEGER, INTENT(IN) :: nd
    TYPE (simplex), TARGET, INTENT(IN) :: simp_ndm1
    REAL (KIND=8), DIMENSION(nd), INTENT(OUT) :: normal

    INTEGER :: i, j, k, m
    !REAL (KIND=8),DIMENSION(:),POINTER :: xk !by DD
    REAL (KIND=8),DIMENSION(nd) :: xk
    
    DO j=1,nd
       DO k=1,nd
          !xk  => simp_ndm1 % opp(k) % pnt % ppp % x
          xk  = simp_ndm1 % opp(k) % pnt % ppp % x
          
          DO i=0,nd-1
             m = MOD(i+j,nd+1)+1
             IF (m <= nd) THEN
                matrix(k,i+1) = xk(m)
             ELSE
                matrix(k,i+1) = 1.
             ENDIF
             
          ENDDO
          
       ENDDO
       CALL lufact (nd, pivot, matrix, normal(j))
    ENDDO
    
  END SUBROUTINE simplex_normal
  !
!********************************************************************
  !
  SUBROUTINE eval_facets (nd, poin, simp, x)

    ! +---------------------+
    ! |   AREA COORDINATE   |
    ! +---------------------+

    INTEGER, INTENT(IN) :: nd
    TYPE (point), INTENT(IN) :: poin
    TYPE (simplex), INTENT(IN) :: simp
    REAL (KIND=8), DIMENSION(nd+1), INTENT(OUT) :: x

    INTEGER :: i, j


    DO i=1,nd+1
       x(i) = simp % fmat(i,nd+1)
       DO j=1,nd
          x(i) = x(i) + simp % fmat(i,j) * poin % x(j)
       ENDDO
    ENDDO
    ! positive x(i) is given by inside points only(trivial&walking search)
  END SUBROUTINE eval_facets
  !
!********************************************************************
  !
  SUBROUTINE build_bubble (nd, bbl_head, poin)

    ! +------------------------+
    ! |   BUBBLE AROUND POIN   |
    ! +------------------------+
    
    ! point owners simplexes
    INTEGER, INTENT(IN) :: nd
    TYPE (simplex), POINTER :: bbl_head
    TYPE (point), TARGET, INTENT(IN) :: poin

    INTEGER :: i, j
    TYPE (simplex), POINTER :: simp, neig


    ! Initialize the surrounding list

    simp  =>  poin % spx

    IF ( .NOT. ASSOCIATED (simp) ) THEN
       WRITE (*,*) 'bubble error, unassociated pointer'
       WRITE (*,*) '(stop)'
       STOP
    ENDIF

    simp % status = - simp % status   !"visited"
   
    CALL movetop_simplex (bbl_head, simp)   
         
         
    ! Link the surrounding vertices to the list
 
    simp  =>  bbl_head
    DO; simp  =>  simp % next; IF ( ASSOCIATED (simp, bbl_head) ) EXIT
       DO i=1,nd+1
          neig  =>  simp % opp(i) % spx ! adjacent simpl...
          IF ( ASSOCIATED (neig) ) THEN
             IF ( neig % status > 0 ) THEN ! "unvisited"
                DO j=1,nd+1
                   IF ( ASSOCIATED (neig % opp(j) % pnt, poin) ) THEN
                         
                      neig % status = - neig % status  ! "visited"
                      CALL movebot_simplex (bbl_head, neig)
                      EXIT        
                   
                   ENDIF          
                ENDDO
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    ! Unflag the surrounding list

    simp  =>  bbl_head
    DO; simp  =>  simp % next; IF ( ASSOCIATED (simp, bbl_head) ) EXIT
       simp % status = ABS(simp % status)
      
    ENDDO

  END SUBROUTINE build_bubble
  !
!********************************************************************
  !
  SUBROUTINE match_bubble_ndm1 (nd, s_ndm1, bbl_head, s_nd_int, s_nd_ext)
!                                   * int_h element
    INTEGER, INTENT(IN) :: nd
    TYPE (simplex), INTENT(IN) :: s_ndm1
    TYPE (simplex), TARGET, INTENT(IN) :: bbl_head
    TYPE (simplex), POINTER :: s_nd_int, s_nd_ext

    LOGICAL :: first
    INTEGER :: i, j, k, kk
    TYPE (point), POINTER :: p
    TYPE (simplex), POINTER :: s
! searches in bbl_head the 2d simplices adjacent to 1d simplex

    first = .TRUE.; NULLIFY (s_nd_int, s_nd_ext)

    s  =>  bbl_head
    DO; s  =>  s % next; IF ( ASSOCIATED (s, bbl_head) ) EXIT
                !c'è solo il simpl di appartenenza
       k = 0
       DO i=1,nd
          kk = k
          p  =>  s_ndm1 % opp(i) % pnt % ppp !2d respective
          DO j=1,nd+1
             IF ( ASSOCIATED (p, s % opp(j) % pnt) ) THEN
                k = k+1 
                EXIT  ! if bubble is adjacent to boundary 
             ENDIF
          ENDDO
          IF (k == kk) EXIT
       ENDDO
       
       IF (k == nd) THEN  !found a base(not a vertex only)
          IF (first) THEN !
             s_nd_int  =>  s   
             first = .FALSE.
          ELSE
             s_nd_ext  =>  s   !exit cycle only when s_nd_ext defined
             EXIT              ! i.e.when for the second time k=2 
          ENDIF
       ENDIF

    ENDDO

  END SUBROUTINE match_bubble_ndm1
  !
!********************************************************************
  !
  SUBROUTINE swap_ie_ndm1 (nd, s_ndm1, s_nd_int, s_nd_ext)
 
    INTEGER, INTENT(IN) :: nd
    TYPE (simplex), POINTER :: s_ndm1, s_nd_int, s_nd_ext

    INTEGER :: i, j
    REAL (KIND=8) :: scal
    TYPE (simplex), POINTER :: s_nd_tmp


    CALL simplex_normal (nd, s_ndm1, vector)

    DO i=1,nd+1
       IF ( ASSOCIATED (s_nd_int % opp(i) % spx, s_nd_ext) ) THEN
          scal = 0.
          DO j=1,nd
             scal = scal + s_nd_int % fmat(i,j) * vector(j)
          ENDDO
          EXIT
       ENDIF
    ENDDO

    IF ( scal > 0.0 ) THEN
       s_nd_tmp  =>  s_nd_int
       s_nd_int  =>  s_nd_ext
       s_nd_ext  =>  s_nd_tmp
    ENDIF

  END SUBROUTINE swap_ie_ndm1
  !
!********************************************************************
  !
  SUBROUTINE match_bubble_nd (nd, k_nd, s_nd, bbl_head, s_nd_int, s_nd_ext)

    INTEGER, INTENT(IN) :: nd, k_nd
    TYPE (simplex), INTENT(IN) :: s_nd
    TYPE (simplex), TARGET, INTENT(IN) :: bbl_head
    TYPE (simplex), POINTER :: s_nd_int, s_nd_ext

    LOGICAL :: first
    INTEGER :: i, j, k, kk
    TYPE (point), POINTER :: p
    TYPE (simplex), POINTER :: s


    first = .TRUE.; NULLIFY (s_nd_int, s_nd_ext)

    s  =>  bbl_head
    DO; s  =>  s % next; IF ( ASSOCIATED (s, bbl_head) ) EXIT

       k = 0
       DO i=1,nd+1
          IF (i /= k_nd) THEN
             kk = k
             p  =>  s_nd % opp(i) % pnt
             DO j=1,nd+1
                IF ( ASSOCIATED (p, s % opp(j) % pnt) ) THEN
                   k = k+1
                   EXIT
                ENDIF
             ENDDO
             IF (k == kk) EXIT
          ENDIF
       ENDDO

       IF (k == nd) THEN
          IF (first) THEN
             s_nd_int  =>  s
             first = .FALSE.
          ELSE
             s_nd_ext  =>  s
             EXIT
          ENDIF
       ENDIF

    ENDDO

  END SUBROUTINE match_bubble_nd
  !
!********************************************************************
  !
  SUBROUTINE swap_ie_nd (nd, k_nd, s_nd, s_nd_int, s_nd_ext)

    INTEGER, INTENT(IN) :: nd, k_nd
    TYPE (simplex), INTENT(IN) :: s_nd
    TYPE (simplex), POINTER :: s_nd_int, s_nd_ext

    INTEGER :: i, j
    TYPE (simplex), POINTER :: s_nd_tmp
    REAL (KIND=8) :: scal


    DO i=1,nd+1
       IF ( ASSOCIATED (s_nd_int % opp(i) % spx, s_nd_ext) ) THEN
          scal = 0.
          DO j=1,nd
             scal = scal + s_nd % fmat(k_nd,j) * s_nd_int % fmat(i,j)
          ENDDO
          EXIT
       ENDIF
    ENDDO

    IF ( scal < 0.0 ) THEN
       s_nd_tmp  =>  s_nd_int
       s_nd_int  =>  s_nd_ext
       s_nd_ext  =>  s_nd_tmp
    ENDIF

  END SUBROUTINE swap_ie_nd
  !
!********************************************************************
  !
  SUBROUTINE flag_shell (nd, int_head_ndm1, spx_head, stc_head_ndm1)

    !     +--------------------------------+
    !     |   FLAG THE INT/EXT SIMPLICES   |
    !     +--------------------------------+

    INTEGER    , INTENT(IN) :: nd
    TYPE (simplex), POINTER :: int_head_ndm1, spx_head, stc_head_ndm1

    TYPE (simplex), POINTER :: s_ndm1, t_ndm1, s_nd_int, s_nd_ext

    s_ndm1  =>  int_head_ndm1 % next

    DO; IF ( ASSOCIATED (s_ndm1, int_head_ndm1) ) EXIT
       t_ndm1  =>  s_ndm1 % next

       CALL build_bubble      (nd, bbl_head, s_ndm1 % opp(1) % pnt % ppp)
       !                                     
       !bbl:first, belonging to simplex,"unvisited"
       !    last , adjacent simplices sharing point,"visited"
       
       CALL match_bubble_ndm1 (nd, s_ndm1, bbl_head, s_nd_int, s_nd_ext)
       !                                             puntatori in uscita
       CALL jointop_simplex   (spx_head, bbl_head)
       !                       pushtop  bbl in spx 
       IF ( ASSOCIATED (s_nd_int) .AND. ASSOCIATED (s_nd_ext) ) THEN
          ! both defined..
          CALL swap_ie_ndm1 (nd, s_ndm1, s_nd_int, s_nd_ext)
          ! swaps s_nd_int  with s_nd_ext
          ! normal vector?
          s_ndm1 % int  =>  s_nd_int 
          s_ndm1 % ext  =>  s_nd_ext

          s_nd_int % status = int_status
          s_nd_ext % status = ext_status
          ! guaina attorno al lato 1d:esterna ed interna,da usare in split_ie
       ELSE

 
          CALL movetop_simplex (stc_head_ndm1, s_ndm1)

       ENDIF

       s_ndm1  =>  t_ndm1
    ENDDO

  END SUBROUTINE flag_shell
  !
!********************************************************************
  !
  SUBROUTINE split_ie_1d (nd, spx_head, int_head, ext_head)

    !     +---------------------------------------------+
    !     |   SPLIT INT/EXT SIMPLEX LISTS --- 1D CASE   |
    !     +---------------------------------------------+
   ! int_head :simplex list,status> frz 
   ! ext_head :adjacents list,status=frz
   
    INTEGER       , INTENT(IN) :: nd
    TYPE (simplex), POINTER    :: spx_head, int_head, ext_head

    INTEGER :: i
    TYPE (simplex), POINTER :: s, t, n


    IF ( nd /= 1 ) THEN
       WRITE (*,*) 'split_1d error, unexpected nd value = ',nd
       WRITE (*,*) '(stop)'
       STOP
    ENDIF

    s  =>  spx_head % next
    DO; IF ( ASSOCIATED (s, spx_head) ) EXIT
       t  =>  s % next
       IF ( s % status > 0 ) THEN ! i.e./= frz_status
          s % status = int_status !=2
          CALL movetop_simplex (int_head, s)
       ENDIF
       s  =>  t
    ENDDO

    s  =>  int_head % next
    DO; IF ( ASSOCIATED (s, int_head) ) EXIT
       t  =>  s % next
       DO i=1,nd+1
          n  =>  s % opp(i) % spx 
          IF ( n % status == 0 ) THEN ! frz_status
             s % status = ext_status ! = 3
             CALL movetop_simplex (ext_head, s)
          ENDIF
       ENDDO
       s  =>  t
    ENDDO

  END SUBROUTINE split_ie_1d
  !
!********************************************************************
  !
  SUBROUTINE split_ie (nd, spx_head, int_head, ext_head)

    !     +---------------------------------+
    !     |   SPLIT INT/EXT SIMPLEX LISTS   |
    !     +---------------------------------+

    INTEGER       , INTENT(IN) :: nd
    TYPE (simplex), POINTER    :: spx_head, int_head, ext_head

    INTEGER :: i
    TYPE (simplex), POINTER :: s, t, n


    !  Link the shell simplices to the internal or external lists

!cerca nella lista spx_head gli int ed ext,come assegnati in split_ie_1d
    s  =>  spx_head % next
    DO; IF ( ASSOCIATED (s, spx_head) ) EXIT
       t  =>  s % next
       IF ( s % status == int_status ) THEN
          CALL movetop_simplex (int_head, s)
          ELSEIF ( s % status == ext_status ) THEN
          CALL movetop_simplex (ext_head, s)
       ENDIF
       s  =>  t
    ENDDO

    !  Link the non-shell internal simplices without crossing boundaries

!fra tutti gli int sposto gli adiacenti con und_status sotto la testa 
!e li pongo int_status

    s  =>  int_head % next
    DO; IF ( ASSOCIATED (s, int_head) ) EXIT
       t  =>  s % next
       DO i=1,nd+1
          n  =>  s % opp(i) % spx !int_head adjacents

          IF ( .NOT. ASSOCIATED (n) ) THEN
             WRITE (*,*) 'fill_shell error, I am lost !!!'
             WRITE (*,*) '(stop)'
             STOP
          ENDIF

          IF ( n % status == und_status ) THEN
             n % status = int_status
             CALL movebot_simplex  (int_head, n)
          ENDIF
       ENDDO
       s  =>  t
    ENDDO

    !  Link the remaining simplices (auxiliary type) as external

!mette ext_status i simplessi che non sono nè int nè ext nè frz(und o stc!)

    s  =>  spx_head
    DO; s  =>  s % next; IF ( ASSOCIATED (s, spx_head) ) EXIT
       IF ( s % status > 0 ) THEN
          s % status = ext_status
       ENDIF
    ENDDO

    IF ( .NOT. ASSOCIATED (spx_head % next, spx_head) ) THEN!cioè se la lista
                                                            !non è vuota
       CALL jointop_simplex (ext_head, spx_head)
       !                     pushtop di spx in ext
    ENDIF

  END SUBROUTINE split_ie
  !
!********************************************************************
  !
  SUBROUTINE merge_ie (spx_head, int_head, ext_head)

    TYPE (simplex), POINTER    :: spx_head, int_head, ext_head

    TYPE (simplex), POINTER :: s


    IF ( ASSOCIATED (spx_head % next, spx_head) ) THEN

       CALL jointop_simplex (spx_head, int_head)
       CALL jointop_simplex (spx_head, ext_head)

       s  =>  spx_head
       DO; s  =>  s % next; IF ( ASSOCIATED (s, spx_head) ) EXIT
          IF ( s % status > 0 ) s % status = und_status
       ENDDO

    ENDIF

  END SUBROUTINE merge_ie
  !
!********************************************************************
  !
  SUBROUTINE set_delaunay (nd, xlb, xrt, axp_head, axs_head, spx_head)

    ! +--------------------------------+
    ! |   DEFINE A NEW TRIANGULATION   |
    ! +--------------------------------+

    INTEGER                    , INTENT(IN) :: nd
    REAL (KIND=8), DIMENSION(:), INTENT(IN) :: xlb, xrt
    TYPE (point)               , POINTER    :: axp_head
    TYPE (simplex)             , POINTER    :: axs_head, spx_head

    INTEGER :: i, j, k
    REAL (KIND=8), PARAMETER :: parm1 = 0.25, parm2 = 2.*parm1+1.!=1.5
    TYPE (simplex), POINTER :: simpl, sfrst, slast
    TYPE (point)  , POINTER :: ip, jp
    
!write(*,*)' xlb,xrt',xlb,xrt
!write(*,*)'    ! ### SIMPLEX LISTS ###'

    DO i=1,nd+3 ! 4-5 
       simpl  =>  new_simplex ()
       CALL pushtop_simplex (axs_head, simpl) !axs testa in basso
    END DO


!write(*,*)'     ! ### BOX POINT LIST ### ' 

    DO i=1,nd+1 !2-3
       ip  =>  new_point ()
       CALL pushtop_point (axp_head, ip)!axp testa in basso,lista dei punti del 
    END DO                              !simplesso di supporto

    ip  =>  axp_head % next
    DO j=1,nd       
       ip % x(j) = xlb(j)-parm1*ABS(xrt(j)-xlb(j))
    END DO
    

    DO i=1,nd
       ip  =>  ip % next                !coordinate dei punti di supporto   
       DO j=1,nd
          ip % x(j) = xlb(j)-parm1*ABS(xrt(j)-xlb(j))     
       END DO
       ip % x(i) = xlb(i)+(parm2*nd-parm1)*(xrt(i)-xlb(i))
       
    END DO

!write(*,*)'     ! ### POINT TOPOLOGY ###'
!dati i punti,stabilire le connessioni simplesso-punti e viceversa
!
    sfrst  =>  axs_head % next
    slast  =>  axs_head % prev
 !write(*,*)' 1'  
    ip  =>  axp_head; simpl  =>  sfrst!=axs head%next
  ! write(*,*)' 2'
    DO i=1,nd+1                                    
       !appartenenza dei (tre)punti al simplesso(vista da entrambi) 
       ip  =>  ip % next                                        
       slast % opp(i) % pnt  =>  ip             
       ip % spx  =>  slast                              
                                                
       simpl  =>  simpl % next                   
       jp  =>  axp_head                        
      !  write(*,*)' +'
       DO j=1,nd+1                             
          jp  =>  jp % next                                  
          simpl % opp(j) % pnt  =>  jp
       END DO
      ! write(*,*)' ++'                                
!per slast corrispondenza simplesso<->punti
!per gli altri,a partire da frst solo corrispondenza simplesso->punti
       NULLIFY (simpl % opp(i) % pnt)!annulla il rif all'i-esimo pto
       NULLIFY (sfrst % opp(i) % pnt)
    END DO

!###sfrst e' un simplesso fantasma che serve per creare una vicinanza 
!###di simplessi che non esistono e rispettare cosi la struttura dei dati

!write(*,*)'     ! ### CONNECTIVITY MATRICES ###'

    simpl  =>  sfrst
    DO k=1,nd+1                !il ciclo inizia dopo sfrst!!
       simpl  =>  simpl % next
       DO i=1,nd+1                    !  simpl 1 2 3     1 2 0
          DO j=1,nd+1                 !        1 2 3  o  1 2 0
             simpl % cmat(i,j) = j    !        1 2 3     0 0 0 
          ENDDO
       ENDDO       
    END DO

    DO i=1,nd+1
       DO j=1,nd+1                    ! sfrst 0 0 0    slast 1 2 3  
          sfrst % cmat(i,j) = 0       !       0 0 0  e       1 2 3 
          slast % cmat(i,j) = j       !       0 0 0          1 2 3
       ENDDO
    ENDDO
    
!write(*,*)'     ! ### NEIGHBOUR TOPOLOGY ###'

    simpl  =>  sfrst

    DO i=1,nd+1
       simpl  =>  simpl % next

       DO j=1,nd+1                          !per ogni simpl dopo sfrst 
          simpl % opp(j) % spx  =>  sfrst   !i suoi opposti sono sempre sfrst
       END DO

       NULLIFY (sfrst % opp(i) % spx)

       simpl % opp(i) % spx  =>  slast      
       slast % opp(i) % spx  =>  simpl      
    END DO

!simpl4

!        opp(1)%spx slast    opp(1)%pnt p1
!(simpl3) opp(2)%spx sfrst    opp(2)%pnt p2    
!        opp(3)%spx sfrst    opp(3)%pnt 

!       opp(1)%spx sfrst    opp(1)%pnt p1
!simpl2 opp(2)%spx slast    opp(2)%pnt 
!       opp(3)%spx sfrst    opp(3)%pnt p3

!       opp(1)%spx sfrst    opp(1)%pnt 
!simpl1 opp(2)%spx sfrst    opp(2)%pnt p2
!       opp(3)%spx slast    opp(3)%pnt p3

!       opp(1)%spx          opp(1)%pnt 
!sfrst  opp(2)%spx          opp(2)%pnt 
!       opp(3)%spx          opp(3)%pnt 

! axs_head

!      opp(1)%spx simpl1    opp(1)%pnt p1
!slast opp(2)%spx simpl2    opp(2)%pnt p2
!      opp(3)%spx simpl3    opp(3)%pnt p3 

!i punti p1 p2 p3 appartengono tutti ad slast 

    ! ### SIMPLEX GEOMETRY ###

    simpl  =>  sfrst
    simpl % status = frz_status
    
    DO i=1,nd+1
       
       simpl  =>  simpl % next
       simpl % status = frz_status   !per simpl1,2,3
    END DO
   
    simpl  =>  simpl % next
    simpl % status = und_status !l'ultimo,in ogni caso non utilizzato
    
    CALL simplex_facets (nd, simpl) !sempre per l'ultimo
   
    CALL movetop_simplex (spx_head, simpl)
    
  END SUBROUTINE set_delaunay
  !
!********************************************************************
  !
  ! Cerca il simplesso nel quale si trova il nuovo punto per introdurlo
  ! nella lista del_head (delete_basis)
  SUBROUTINE delete_first (nd, spx_head, pnew, err_status)

    INTEGER               , INTENT(IN)  :: nd
    TYPE (simplex), TARGET, INTENT(IN)  :: spx_head
    TYPE (point)  , POINTER             :: pnew
    LOGICAL               , INTENT(OUT) :: err_status

    INTEGER :: status
    TYPE (simplex), POINTER :: simp

    err_status = .TRUE.
    IF ( ASSOCIATED (pnew % old) ) THEN

       simp  =>  pnew % old % spx

       status = 1			! 1 to activate search on maingrid
       CALL walking_search (nd, pnew, simp, status)
       IF(status == 1) THEN		! point has been located
          err_status = .FALSE.
          !WRITE(*,*)'Point located by Walking_search'
       ELSE IF(status == 0) THEN		! walking_search failed
          CALL trivial_search (nd, spx_head, pnew, simp, err_status)
          IF(err_status) THEN
             !WRITE(*,*)'Walking & Trivial search failed (delete_first)'
             RETURN
          ENDIF
          !WRITE(*,*)'Point located by Trivial_search - after walking failed'
       ELSE IF(status == -1) THEN
          !WRITE(*,*)'Walking says: point is not here! (delete_first)'
          RETURN
       ENDIF

    ELSE

       CALL trivial_search (nd, spx_head, pnew, simp, err_status)
       IF(err_status) THEN
          !WRITE(*,*)'Trivial search failed (delete_first)'
          RETURN
       ENDIF
       !WRITE(*,*)'Point located by Trivial_search'
    ENDIF
    
    CALL delete_basis (nd, pnew, simp)

  END SUBROUTINE delete_first
  !
!********************************************************************
  
  ! la base è costituita dal simplesso in cui si trova il nuovo punto 
  ! da inserire e dai suoi adiacenti,tutti posti con status<0(visited)
  
  SUBROUTINE delete_basis (nd, pnew, simp)

    INTEGER       , INTENT(IN) :: nd
    TYPE (point)  , POINTER    :: pnew
    TYPE (simplex), POINTER    :: simp

    INTEGER :: i, j, n
    TYPE (simplex), POINTER :: neig

    simp % status = - simp % status
    CALL movetop_simplex (del_head, simp) !il primo simplesso e' quello in cui
                                          !cade il punto 

    DO j=1,nd+1
       neig  =>  simp % opp(j) % spx
       IF ( neig % status > 0 ) THEN   !simpl non controllato
         CALL eval_facets (nd, pnew, neig, vector)

         n = 0
         DO i=1,nd+1
             IF ( vector(i) >= neg_tol ) n = n+1
         ENDDO

         IF ( n == nd+1 ) THEN        !il punto cade nel simplesso adiacente
           IF ( (neig%status /= und_status) .AND. (neig%status /= int_status) ) THEN
              WRITE(*,*)'delete_basis: element status:', neig%status
           ENDIF
             neig % status = - neig % status
             CALL movetop_simplex (del_head, neig)
             EXIT
          ENDIF
       ENDIF
    ENDDO
   
          
  END SUBROUTINE delete_basis !poi delete cavity in eucl_delaunay
  !
!********************************************************************
  !
  SUBROUTINE update_topology (nd, pnew)

    ! +------------------------------------------+
    ! |    UPDATE THE CAVITY TOPOLOGICAL DATA    |
    ! +------------------------------------------+

    INTEGER, INTENT(IN) :: nd
    TYPE (point), TARGET, INTENT(IN) :: pnew

    INTEGER :: mi, ni
    TYPE (simplex), POINTER :: del, adj, new

    ! Compute the new <---> adj data
    

    del  =>  del_head
    DO; del  =>  del % next; IF ( ASSOCIATED (del, del_head) ) EXIT
       DO ni=1,nd+1                      
          adj  =>  del % opp(ni) % spx   
          IF ( adj % status <= 0 ) THEN   
             mi = del % cmat(ni,ni)
             CALL adj_new  !connessione tra adj e il nuovo triangolo
          ENDIF            !opposto che nascerà dalla base
       ENDDO
    ENDDO

    ! Compute the new <---> new data

    del  =>  del_head
    DO; del  =>  del % next; IF ( ASSOCIATED (del, del_head) ) EXIT
       DO ni=1,nd+1
          adj  =>  del % opp(ni) % spx
          IF ( adj % status <= 0 ) THEN  
             new  =>  adj % opp(del % cmat(ni,ni)) % spx 
             DO mi=1,nd+1                                 
                IF (mi /= ni) CALL new_new 
             ENDDO                         
          ENDIF
       ENDDO
    ENDDO

    ! Compute the point ----> simplex 
    !    (appartenenza dei punti del simplesso al simplesso stesso!!!
    ! Reset to "unvisited" the adjacent simplices
    

    new  =>  new_head
    DO; new  =>  new % next; IF ( ASSOCIATED (new, new_head) ) EXIT
       DO mi=1,nd+1
          new % opp(mi) % pnt % spx  =>  new
          new % opp(mi) % spx % status = ABS(new % opp(mi) % spx % status)
       ENDDO
    ENDDO

  CONTAINS
!**********************************************
    SUBROUTINE adj_new

      INTEGER :: i, ni

      ! ###  NEW SIMPLEX  ###

      new  =>  gar_head % next  
      ! in insert_point ,gar è riempita con i del

      IF ( ASSOCIATED (new, gar_head) ) THEN   !uso per la prima volta new_head
                                               !ne definisco gli elementi,che in
         new  =>  new_simplex ()               !fine sono 6.   status?
         CALL pushtop_simplex (new_head, new)  !in grid_types è und

      ELSE  !dopo che il primo punto è stato inserito...
            !del_head ->in gar_head=>al secondo punto si arriva qui
            !new_head è come in new_head_simplex,cioè punta a se stesso
               
         new % status = und_status             !new ancora in gar_h;cambia stato
         CALL movetop_simplex (new_head, new)  !preleva da gar e mette in new_h
                                               !poi in gar_h la nuova del_h   
      ENDIF                                    !e new in spx(quindi new vuota)  

      ni = adj % cmat(mi,mi)

      ! ###  NEW POINT-POINTERS  ###

      DO i=1,nd+1
         new % opp(i) % pnt  =>  del % opp(i) % pnt
      ENDDO
      new % opp(ni) % pnt  =>  pnew  !cambia solo questo punto nel triangolo 
                                     !opposto ad adj    

      ! ###  NEW-ADJ SIMPLEX-POINTERS

      adj % opp(mi) % spx  =>  new
      new % opp(ni) % spx  =>  adj

      !   ###  NEW-ADJ CONNECTIVITY MATRIX  ###

      DO i=1,nd+1
         new % cmat(ni,i) = del % cmat(ni,i)
      ENDDO
      ! la connettività di adj rimane invariata,perchè adj vede new solo come 
      ! spostamento di un punto del simplesso appartenente alla base
      
      ! la connettività relativa agli altri due punti di new 
      ! è trattata in new_new 

    END SUBROUTINE adj_new

!***********************************************

! l'operazione di ridefinizione della triangolazione avviene su una base 
! di quattro triangoli che formano la cavità(vs. cavità effettiva...)

! i new sono visti dagli adj rispettivi come cambiamento di forma (spostamento 
! di un punto su pnew)dei del(rispettivi agli adj)
! 
! vedere la configurazione finale dei new 
! conservandone i rispettivi del di origine
!
! la corrispondenza new-new è 
! una corrispondenza dei rispettivi del-del
!
! a partire dall'angolo di new
! (del quale si vuole completare la riga della matrice di connettività new-new)
! si individua il del di interfaccia:
! - se è un petalo -> interfaccia diretta
! - oppure passaggio intermedio per il simplesso centrale della base
!   perchè è la SUA matrice di connettività ad interessare
!   (diventa un sostituto "vedente" del simplesso del di partenza)

    SUBROUTINE new_new

      TYPE (simplex), POINTER :: d, e
      INTEGER :: i, j, m, n, x, y

      j = 0; m = mi; n = ni

      DO i=1,nd+1
         IF ((i /= mi) .AND. (i /= ni)) THEN ! rileva l'indice che rimane
            j = j+1
            pivot(j) = i
         ENDIF
      ENDDO

      d  =>  del

      DO
         e  =>  d % opp(m) % spx
         IF ( e % status > 0 ) THEN  ! cioè finchè e non è un adj
            x = d % cmat(m,m)
            y = d % cmat(m,n)

            DO i=1,nd-1
               pivot(i) = d % cmat(m, pivot(i))
            ENDDO

            m = y  !=d % cmat(m,n)
            n = x  !=d % cmat(m,m)

            d  =>  e
         ELSE    ! quando raggiunge un petalo
            EXIT ! cioè quando per opp(m)trova un'interfaccia diretta 
         ENDIF   ! (già esistente)con un petalo 
      ENDDO

      new % opp(mi) % spx  =>  d % opp(m) % spx % opp(d % cmat(m,m)) % spx

      new % cmat(mi,mi) = n
      new % cmat(mi,ni) = m

      j = 0

      DO i=1,nd+1
         IF ((i /= mi) .AND. (i /= ni)) THEN
            j = j+1
            new % cmat(mi,i) = pivot(j)
         ENDIF
      ENDDO

    END SUBROUTINE new_new


  END SUBROUTINE update_topology
  !
!********************************************************************
  !
  SUBROUTINE update_geometry (nd, head)

    ! +----------------------------------------+
    ! |    UPDATE THE HEAD GEOMETRICAL DATA    |
    ! +----------------------------------------+

    INTEGER, INTENT(IN) :: nd
    TYPE (simplex), TARGET, INTENT(INOUT) :: head
    TYPE (simplex), POINTER :: s

    s  =>  head

    DO; s  =>  s % next; IF ( ASSOCIATED (s, head) ) EXIT
       CALL simplex_facets (nd, s) ! calcola fmat per i nuovi simplessi
    ENDDO

  END SUBROUTINE update_geometry


END MODULE delaunay
