MODULE grid_types

  IMPLICIT NONE

  integer :: meshtype
  INTEGER, PARAMETER :: maxp = 2, maxs = maxp+1

  INTEGER, PARAMETER :: nul_cnd = 0, nul_par = 0
  INTEGER, PARAMETER :: frz_status = 0, und_status = 1,   &
       &   int_status = 2, ext_status = 3,   &
       &   stc_status = 4

  REAL(KIND=8), PARAMETER :: zero = 0.0d0
  REAL(KIND=8), PARAMETER :: one  = 1.0d0
  
  !----------------------------------------------------------------------
  TYPE opposits
     TYPE (point),   POINTER :: pnt
     TYPE (simplex), POINTER :: spx
  END TYPE opposits
  
  !----------------------------------------------------------------------
  TYPE point

     LOGICAL :: visited
     INTEGER :: index, cnd, ibk, parent, what, opnd, option, mark
     
     REAL (KIND=8) :: leng, bou_leng1
     REAL (KIND=8) :: leng_max, leng_min, leng_min2o
     REAL (KIND=8) :: height, past_height
     REAL (KIND=8) :: theta, sum_theta

     TYPE (point),   POINTER :: next,  prev,  old,  ppp
     TYPE (point),   POINTER :: next1, prev1, old1, ppp1 
     TYPE (point),   POINTER :: wake 
     TYPE (point),   POINTER :: jump_f, jump_b 
     TYPE (simplex), POINTER :: spx

     REAL (KIND=8), DIMENSION (maxp)      :: x
     REAL (KIND=8), DIMENSION (2*maxp-1)  :: ellipse
     REAL (KIND=8), DIMENSION (maxp,maxp) :: metric
     REAL (KIND=8), DIMENSION (2,2) :: nv1

  END TYPE point
  
  !----------------------------------------------------------------------
  TYPE simplex
     LOGICAL       :: visited
     INTEGER       :: index, status
     REAL (KIND=8) :: det

     TYPE (simplex), POINTER :: prev, next, int, ext

     INTEGER,         DIMENSION (maxs,maxs) :: cmat
     REAL (KIND=8),   DIMENSION (maxs,maxs) :: fmat
     TYPE (opposits), DIMENSION (maxs)      :: opp
  END TYPE simplex
  !######    CMAT      #####
  !
  !i simplessi adiacenti(opposti,3) di un simplesso si identificano tramite 
  !un indice locale "i" da usare in opp(i)
  !per raggiungere velocemente i 9 adiacenti dei 3 opposti serve
  !la numerazione locale dei 3 opposti.
  !cmat(i,j)
  !colonna j ->informazioni sul simplesso opposto al nodo che localmente 
  !            ha indice j
  !cmat(j,j)=numero locale del nodo che sta di fronte a j
  !cmat(i,j)=numero locale del nodo i nella numerazione locale del simplesso 
  !          opposto a j
  !          /\
  !         / 3\       **numerazione antioraria
  !        /    \ 
  !       /1    2\   ->   cmat = .  .  .
  !       --------               .  .  . 
  !       \1    3/               1  3  2       
  !        \    /            
  !         \2 /
  !          \/             
  
  !----------------------------------------------------------------------
  TYPE stratum
     TYPE(point),DIMENSION(:),POINTER :: pnt_list
     INTEGER,DIMENSION(:),ALLOCATABLE :: num_el  ! numero di elementi per lato
     INTEGER,DIMENSION(:),ALLOCATABLE :: max_lay ! max numero di strati per lato
     !
  END TYPE stratum
  
  !----------------------------------------------------------------------
  TYPE box      ! definisce la struttura di supporto per il controllo
     ! di sovrapposizione dei fronti avanzanti
     !TYPE (point), POINTER :: corner1, corner2 
     REAL(KIND=8), DIMENSION(2) :: corner1, corner2  !coordinate allargate 
     TYPE(point), POINTER       :: head, c1, c2, virt_new_head
     ! punti iniziale e finale che
     ! definiscono gli estremi del box
     INTEGER :: level, structured, num_points, active, cross_ovlp, auto_ovlp
     ! indice di passaggio e identificativo del tipo
     ! di lato al quale il box fa riferimento
     TYPE (box), DIMENSION(:), POINTER :: small_box 
     ! riferimento ai sottoinsiemi di punti nei quali e' 
     ! suddiviso l'intero lato
  END TYPE box
  !----------------------------------------------------------------------

CONTAINS


  FUNCTION new_point () RESULT (p)

    IMPLICIT NONE

    TYPE (point), POINTER :: p

    INTEGER :: i


    ALLOCATE ( p )

    p % visited = .FALSE.
    p % index   = 0
    p % cnd     = nul_cnd
    p % ibk     = 0
    p % parent  = nul_par
    p % mark    = 0

    p % what   = 0  !######
    p % opnd   = 0  !######
    p % option = 0

    NULLIFY ( p % prev )
    NULLIFY ( p % next )
    NULLIFY ( p % old )
    NULLIFY ( p % ppp )
    NULLIFY ( p % spx )

    NULLIFY ( p % next1 ) !##########
    NULLIFY ( p % prev1 ) !##########
    NULLIFY ( p % old1 ) !##########
    NULLIFY ( p % ppp1 ) !##########
    NULLIFY ( p % wake ) !##########
    NULLIFY ( p % jump_f ) !##########
    NULLIFY ( p % jump_b) !##########

    p % leng    = zero
    p % leng_max   = zero
    p % leng_min   = zero      
    p % bou_leng1  = zero
    p % leng_min2o = zero
    p % x       = zero
    p % ellipse = zero
    p % metric  = zero
    p % nv1     = zero !############
    p % height = zero
    p % past_height = zero
    p % theta = zero
    p % sum_theta = zero
    DO i=1,maxp
       p % metric(i,i)  = one
    ENDDO

  END FUNCTION new_point
  !------------------------

  FUNCTION new_simplex () RESULT (s)

    IMPLICIT NONE

    TYPE (simplex), POINTER :: s

    INTEGER :: i


    ALLOCATE ( s )

    s % visited = .FALSE.
    s % index   = 0
    s % status  = und_status
    s % det     = zero

    NULLIFY ( s % prev )
    NULLIFY ( s % next )
    NULLIFY ( s % int )
    NULLIFY ( s % ext )

    DO i=1,maxs
       NULLIFY (s % opp(i) % pnt)
       NULLIFY (s % opp(i) % spx)
    ENDDO

    s % cmat    = 0
    s % fmat    = zero

  END FUNCTION new_simplex
  !°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  SUBROUTINE copy_point (p1,p2) 
    ! to copy p1 in p2
    IMPLICIT NONE

    TYPE (point), POINTER :: p1,p2

    INTEGER :: i

    p2 % visited = p1 % visited
    p2 % index   = p1 % index
    p2 % cnd     = p1 % cnd
    p2 % ibk     = p1 % ibk
    p2 % parent  = p1 % parent
    p2 % mark    = p1 % mark

    p2 % what   = p1 % what
    p2 % opnd   = p1 % opnd
    p2 % option = p1 % option

    NULLIFY ( p2 % prev )
    NULLIFY ( p2 % next )
    NULLIFY ( p2 % old )
    NULLIFY ( p2 % ppp )
    NULLIFY ( p2 % spx )

    NULLIFY ( p2 % next1 ) 
    NULLIFY ( p2 % prev1 ) 
    NULLIFY ( p2 % old1 ) 
    NULLIFY ( p2 % ppp1 ) 
    NULLIFY ( p2 % wake ) 

    IF(ASSOCIATED(p1 % prev))THEN 
       p2 % prev => p1 % prev    
    ENDIF
    IF(ASSOCIATED(p1 % next))THEN   
       p2 % next => p1 % next
    ENDIF
    IF(ASSOCIATED(p1 % old))THEN   
       p2 % old  => p1 % old
    ENDIF
    IF(ASSOCIATED(p1 % ppp))THEN  
       p2 % ppp  => p1 % ppp
       !escluso perchè non esiste ancora!
    ENDIF
    IF(ASSOCIATED(p1 % spx))THEN   
       p2 % spx  => p1 % spx
    ENDIF
    IF(ASSOCIATED(p1 % next1))THEN   
       p2 % next1  => p1 % next1
    ENDIF
    IF(ASSOCIATED(p1 % prev1))THEN   
       p2 % prev1  = p1 % prev1
    ENDIF
    IF(ASSOCIATED(p1 % old1))THEN   
       p2 % old1   => p1 % old1
    ENDIF
    IF(ASSOCIATED(p1 % ppp1))THEN   
       p2 % ppp1   => p1 % ppp1
    ENDIF
    IF(ASSOCIATED(p1 % wake))THEN   
       p2 % wake   => p1 % wake
    ENDIF
    IF(ASSOCIATED(p1 % jump_f))THEN   
       p2 % jump_f   => p1 % jump_f
    ENDIF
    IF(ASSOCIATED(p1 % jump_b))THEN   
       p2 % jump_b   => p1 % jump_b
    ENDIF

    p2 % leng    = p1 % leng
    p2 % leng_max   = p1 % leng_max
    p2 % leng_min   = p1 % leng_max
    p2 % bou_leng1  = p1 % bou_leng1
    p2 % leng_min2o = p1 % leng_min2o
    p2 % x       = p1 % x
    p2 % ellipse = p1 % ellipse
    p2 % metric  = p1 % metric

    p2 % nv1     = p1 % nv1
    p2 % height  = p1 % height
    p2 % past_height = p1 % past_height
    p2 % theta   = p1 % theta
    p2 % sum_theta   = p1 % sum_theta

    DO i=1,maxp
       p2 % metric(i,i)  = p1 % metric(i,i)
    ENDDO

  END SUBROUTINE copy_point
  !--------------------------------------------------------------
  FUNCTION new_box () RESULT (b)

    IMPLICIT NONE

    TYPE (box), POINTER :: b

    ALLOCATE ( b )

    b % corner1 = 0
    b % corner2 = 0

    NULLIFY(b % head)
    NULLIFY(b % c1) 
    NULLIFY(b % c2)
    NULLIFY(b % virt_new_head)

    !  b % active     = 0
    b % structured = 0
    b % num_points = 0
    b % level      = 0
    b % active     = 0
    b % cross_ovlp = 0
    b % auto_ovlp  = 0
  END FUNCTION new_box
  !------------------------


END MODULE grid_types
