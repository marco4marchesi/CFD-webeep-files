MODULE eigsys_2d

  ! PARAMETERS

  REAL(kind=8), PARAMETER, PRIVATE :: pi = 3.1415927464102
  REAL(kind=8), PARAMETER, PRIVATE :: zero = 0.d0
  REAL(kind=8), PARAMETER, PRIVATE :: one  = 1.d0
  REAL(kind=8), PARAMETER, PRIVATE :: tiny    = 1.d-8
  REAL(kind=8), PARAMETER, PRIVATE :: epsilon = 1.d-12

  PUBLIC :: regul_2d, invert_2d, eigensys_2d, eig_sys_2d

CONTAINS
  !
!********************************************************************
  !
  FUNCTION regul_2d(a) RESULT(ar)

    !  Transform a in ar symmetric matrix

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: a
    REAL(KIND=8), DIMENSION(2,2)              :: ar

    LOGICAL                      :: err
    REAL(KIND=8)                 :: mod_array 
    REAL(KIND=8), DIMENSION(2)   :: lambda 
    REAL(KIND=8), DIMENSION(2,2) :: p, alambda

    ar = a

    CALL eigensys_2d(ar, lambda, p, err)

    mod_array = ABS(ar(1,1)) + ABS(ar(2,2)) + ABS(ar(1,2)) + ABS(ar(2,1))
    alambda(1,1) = ABS(lambda(1)) + tiny*mod_array + epsilon
    alambda(2,2) = ABS(lambda(2)) + tiny*mod_array + epsilon
    alambda(1,2) = 0.d0
    alambda(2,1) = 0.d0

    CALL invert_2d(p, ar, err)

    IF(err) THEN
       WRITE(*,*)'p singular -> a not regularized (regul_2)'
       write(*,*)'-- STOP -- (eigsys_2d.f90)'
       STOP
    ELSE
       p  = MATMUL(p, alambda) 
       ar = MATMUL(p, ar) 
    ENDIF

  END FUNCTION regul_2d
  !
!********************************************************************
  !
  SUBROUTINE invert_2d(a, am1, err)

    ! Invert matrix a and store result in am1

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: a
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: am1
    LOGICAL                     , INTENT(OUT) :: err

    REAL(KIND=8):: det, mod_array

    mod_array = ABS(a(1,1)) + ABS(a(2,2)) + ABS(a(1,2)) + ABS(a(2,1))
    mod_array = mod_array**2
    det = a(1,1)*a(2,2) - a(1,2)*a(2,1)

    IF (ABS(det) .LE. tiny*mod_array ) THEN

       WRITE(*,*) ' INVERT_2d: determinant is zero'
       WRITE(*,*)a(1,1), a(1,2)
       WRITE(*,*)a(2,1), a(2,2)
       WRITE(*,*)'determinant:', det
       err = .TRUE.

    ELSE

       det      = one/det
       am1(1,1) = a(2,2)*det 
       am1(2,2) = a(1,1)*det
       am1(1,2) = -a(1,2)*det
       am1(2,1) = -a(2,1)*det
       err = .FALSE.

    ENDIF

  END SUBROUTINE invert_2d
  !
!********************************************************************
  !
  SUBROUTINE eigensys_2d(a, lambda, p, exit_status)

    !  Computation of Eigenvalues and EigenFactors for non singular matrices
    !  ** 2D ONLY **

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: a
    REAL(KIND=8), DIMENSION(:)  , INTENT(OUT) :: lambda 
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: p
    LOGICAL                     , INTENT(OUT) :: exit_status 

    !INTEGER      :: i, j
    REAL(KIND=8) :: disc, abs_mod, mod1, mod2, alpha, beta, qq, scale
    REAL(KIND=8), DIMENSION(2,2) :: aa !, lam , inv_p

    !  Scaling!!!

    aa = a
    scale = MIN(ABS(aa(1,1)), ABS(aa(1,2)), ABS(aa(2,2)))
    if(scale > one) then
       aa = aa/scale
    else
       scale = one
    endif

    mod2    = ABS(aa(1,2))+ABS(aa(2,1))
    mod1    = ABS(aa(1,1))+ABS(aa(2,2))+mod2

    alpha   = aa(1,1)*aa(2,2)-aa(1,2)*aa(2,1)
    beta    = aa(1,1)+aa(2,2)
    disc    = beta**2-4.d0*alpha
    abs_mod = tiny*mod1**2

    !  Check if NEGATIVE

    IF(disc <= -abs_mod) THEN 

       WRITE(*,*)'eigensys_2d: matrix has no real eigenvalue'
       WRITE(*,*)'discriminant:', disc
       STOP

    ELSE IF((mod2/mod1 < epsilon).OR. (disc < zero)) THEN

       !  Check if DIAGONAL

       IF (ABS(aa(1,1)-aa(2,2)) < epsilon) THEN

          ! ISOTROPE

          lambda = aa(1,1)
          exit_status = .FALSE.

       ELSE

          ! ANISOTROPE

          lambda(1) = aa(1,1)
          lambda(2) = aa(2,2)
          exit_status = .TRUE.

       ENDIF

       p(1,1) = one
       p(1,2) = zero
       p(2,1) = zero
       p(2,2) = one

    ELSE

       exit_status = .TRUE.

       ! Numerical Recipes way of lambda finding

       qq = .5d0*(beta+SIGN(1.d0, beta)*SQRT(disc))
       lambda(1) = qq
       lambda(2) = alpha/qq

       qq     = aa(1,1)-lambda(1)
       mod1   = SQRT(qq**2+aa(1,2)**2)
       p(1,1) = aa(1,2)/mod1 
       p(2,1) = -qq/mod1

       qq     = aa(2,2)-lambda(2)
       mod2   = SQRT(aa(2,1)**2+qq**2)
       p(1,2) = qq/mod2 
       p(2,2) = -aa(2,1)/mod2

    ENDIF

    !  Scaling back

    lambda = lambda*scale

    ! check the decomposition (developement test)


  END SUBROUTINE eigensys_2d
  !
!********************************************************************
  !
   SUBROUTINE eig_sys_2d(a, lambda, teta, exit_status)

    !  EIGENSYSTEM  *** 2D ONLY ***  NOT UPDATED*******

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: a
    REAL(KIND=8), DIMENSION(:)  , INTENT(OUT) :: lambda 
    REAL(KIND=8)                , INTENT(OUT) :: teta
    LOGICAL                     , INTENT(OUT) :: exit_status

    REAL(KIND=8):: disc, mod1, mod2, alpha, beta, e_diag, qq

    e_diag = ABS(a(1,2)) + ABS(a(2,1))
    mod1   = ABS(a(1,1)) + ABS(a(2,2)) + e_diag
    mod2   = tiny*mod1**2

    alpha  = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    beta   = a(1,1) + a(2,2)
    disc   = beta**2 - 4.d0*alpha

    IF(disc <= -mod2) THEN 

       WRITE(*,*) ' VECTP_VALP_2: no real eigen value'
       WRITE(*,*)'disc:', disc
       STOP

    ELSE IF(ABS(disc) < mod2) THEN

       lambda(1) = 0.5d0*beta
       lambda(2) = lambda(1)
       teta = zero
       exit_status = .FALSE.

    ELSE

       ! Numerical Recipes way of root finding

       qq = .5d0*(beta+SIGN(1.d0, beta)*SQRT(disc))
       lambda(1) = qq
       lambda(2) = alpha/qq

!**!  Standard, old, handy fashion...
!**
!**       lambda(1) = 0.5d0*(beta + SQRT(disc))
!**       lambda(2) = 0.5d0*(beta - SQRT(disc)) 

       alpha = 2.d0*a(1,2)/(lambda(1)-lambda(2))
       teta = .5d0*ASIN(alpha)

       IF(teta < zero) THEN
          alpha     = lambda(1)
          lambda(1) = lambda(2)
          lambda(2) = alpha
          alpha = lambda(1)*COS(teta)**2+lambda(2)*SIN(teta)**2-a(1,1)
          IF(ABS(alpha) > 1.d-5) THEN
             teta = teta+.5d0*pi
          ELSE
             teta = -teta
          ENDIF
       ENDIF

    ENDIF

  END SUBROUTINE eig_sys_2d
  !
!********************************************************************
  !
END MODULE eigsys_2d


