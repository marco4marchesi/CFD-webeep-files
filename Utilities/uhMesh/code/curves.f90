  MODULE curves

  !   1)  Building a Set of Curves
  !
  !   ...
  !
  !   2)  Loading and Saving
  !
  !   ...
  !
  !   3)  Evaluation
  !
  !   The curve parameter ranges from 0 to 1, and is proportional to
  !   the curvilinear abscissa of the polygonal line passing through
  !   the spline knots (i.e., if S denotes the curvilinear abscissa of
  !   the polygonal line through the knots, L its length, and P the
  !   parameter, P = S/L).
  
  USE banded

  !-----------------------------------------------------------------------
  IMPLICIT NONE

  TYPE curve
    INTEGER :: nd, ns
    REAL (KIND=8), DIMENSION(:),   POINTER :: s, h 
    REAL (KIND=8), DIMENSION(:,:), POINTER :: x, xs
  END TYPE curve

  TYPE(curve), DIMENSION(:), ALLOCATABLE :: curv, wcurv

  REAL(KIND=8), PARAMETER, PRIVATE :: zero = 0.d0, one = 1.d0

  INTEGER, PARAMETER :: chounp=1, natunp=2, claunp=3
  INTEGER :: n_curv


  PUBLIC :: new_curves, load_curves, save_curves, int_new_curve, &
            evalc_s, evalc_0, evalc_1, evalc_2, evalc_3, evalc_0123


  PRIVATE :: new_curve, load_curve, patch,save_curve, mkun,  & 
             choun, natun, claun, sysun, mkunp, sysunp,      &
             alcond, blcond, alpnts, blpnts, alun, blun, asun, bsun

  !PUBLIC::save_curve BY DD
  !-----------------------------------------------------------------------
  
  CONTAINS


  SUBROUTINE int_new_curve(data_in, g)

  ! Read the knot coordinates and compute the interpolant 
  ! univariate spline internal version -> data_in is an array 
  ! filled elsewhere. (ONLY NATURAL CONDITIONS IMPLEMENTED)
  !---------------------------------------------------------------------
  IMPLICIT NONE

  REAL(kind=8), DIMENSION(:, :), INTENT(IN)    :: data_in
  TYPE (curve),                  INTENT(INOUT) :: g

  REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: xs1, xsn, a, b, c
  INTEGER :: md, ms

  EXTERNAL choun, natun, claun
  !---------------------------------------------------------------------

  md = size(data_in, 2)
  ms = size(data_in, 1)
  
  !WRITE(*,*)'md,ms',md,ms

  ALLOCATE (g % s(ms))
  ALLOCATE (g % x(md,ms))
  ALLOCATE (g % xs(md,ms))

  ALLOCATE (xs1(md))
  ALLOCATE (xsn(md))

  g % nd = md
  g % ns = ms

  g % x  = data_in

  ALLOCATE (a(ms))
  ALLOCATE (b(ms))
  ALLOCATE (c(ms))

  CALL mkun (natunp, g % nd, g % ns, g % s, g % x, &
             xs1, xsn, a, b, c, g % xs)

  DEALLOCATE (c)
  DEALLOCATE (b)
  DEALLOCATE (a)

  DEALLOCATE (xsn)
  DEALLOCATE (xs1)

  END SUBROUTINE int_new_curve



  SUBROUTINE new_curves (idf, ascii, cur)
  !---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: idf
  LOGICAL, INTENT(IN) :: ascii
  TYPE(curve), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: cur

  INTEGER :: i
  !---------------------------------------------------------------------
  
  READ(idf,*)
  READ(idf,*) n_curv

  IF (ALLOCATED(cur)) DEALLOCATE(cur)
  ALLOCATE (cur(n_curv))

  DO i = 1, n_curv
    CALL new_curve (idf, ascii, cur(i))
  ENDDO

  END SUBROUTINE new_curves



  SUBROUTINE load_curves (idf, ascii, cur)
  !---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: idf
  LOGICAL, INTENT(IN) :: ascii
  TYPE(curve), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: cur

  INTEGER :: i
  !---------------------------------------------------------------------

  READ (idf,*)
  READ (idf,*) n_curv

  IF ( ALLOCATED(cur) )  DEALLOCATE(cur)

  ALLOCATE (cur(n_curv))

  DO i = 1, n_curv
     CALL load_curve (idf, ascii, cur(i))
  ENDDO

  END SUBROUTINE load_curves


  
  SUBROUTINE save_curves (idf, ascii,cur)
  !---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: idf
  LOGICAL, INTENT(IN) :: ascii
  TYPE(curve), DIMENSION(:), INTENT(INOUT) :: cur

  INTEGER :: i
  !---------------------------------------------------------------------

  WRITE (idf,'(4x,a6)') 'N_CURV'
  WRITE (idf,'(i10)') n_curv

  DO i = 1, n_curv
    CALL save_curve(idf, ascii, cur(i))
  ENDDO

  END SUBROUTINE save_curves












  SUBROUTINE evalc_s( i, n, s, cur )
  !---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: i
  INTEGER, INTENT(OUT) :: n
  
  REAL (KIND=8), DIMENSION(:), POINTER :: s
  TYPE(curve),   DIMENSION(:), INTENT(INOUT) :: cur
  !---------------------------------------------------------------------

  n  =   cur(i) % ns
  s  =>  cur(i) % s

  END SUBROUTINE evalc_s



  SUBROUTINE evalc_0 (i, se, xe, cur)

  ! Evaluate the 0th-derivative 
  ! of the geometry g
  !---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i
  REAL (KIND=8), DIMENSION(:), INTENT(IN)  :: se
  REAL (KIND=8), DIMENSION(:), INTENT(OUT) :: xe
  TYPE(curve),   DIMENSION(:), INTENT(INOUT) :: cur
  !---------------------------------------------------------------------

  CALL c_0 (cur(i) % ns, cur(i) % s, cur(i) % x, cur(i) % xs, se, xe)

  CONTAINS

        SUBROUTINE c_0 (n, s, x, xs, se, xe)

        ! Evaluate hermite univariate interpolation at point si
        ! si ranges from zero to one
        ! zeroth derivative xi = xi(si)

        ! Warning: to be updatated in order to eliminate "n"

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        REAL (KIND=8), DIMENSION(:),   INTENT(IN)  :: s, se
        REAL (KIND=8), DIMENSION(:,:), INTENT(IN)  :: x, xs
        REAL (KIND=8), DIMENSION(:),   INTENT(OUT) :: xe

        INTEGER       :: i1, i2, k
        REAL (KIND=8) :: uu, du, a1, b1, a2, b2, h01, h02, h03, h04


        uu = se(1)

        CALL patch (n, s, uu, du, i1, i2)

        a1  = uu
        b1  = 1.-uu
        a2  = a1**2
        b2  = b1**2

        h01 = (b1+3.*a1)*b2
        h02 = +a1*b2
        h03 = -a2*b1
        h04 = (a1+3.*b1)*a2

        DO k=1,SIZE(x,1)
           xe(k) = x(k,i1)*h01 + du*xs(k,i1)*h02 + du*xs(k,i2)*h03 + x(k,i2)*h04
        ENDDO

        END SUBROUTINE c_0

  END SUBROUTINE evalc_0


  
  SUBROUTINE evalc_1 (i, se, xe,cur)

  ! Evaluate the 1st-derivative 
  ! of the geometry g
  !---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i
  REAL (KIND=8), DIMENSION(:), INTENT(IN)  :: se
  REAL (KIND=8), DIMENSION(:), INTENT(OUT) :: xe
  TYPE(curve), DIMENSION(:),INTENT(INOUT):: cur !!!
  !---------------------------------------------------------------------

  CALL c_1 (cur(i) % ns, cur(i) % s, cur(i) % x, cur(i) % xs, se, xe)

  CONTAINS

        SUBROUTINE c_1 (n, s, x, xs, se, xe)

        ! Evaluate hermite univariate interpolation at point si

        ! si ranges from zero to one
        ! first derivative xi = (D/Ds)^1 (xi(si)) w.r.t. s

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        REAL (KIND=8), DIMENSION(:),   INTENT(IN)  :: s, se
        REAL (KIND=8), DIMENSION(:,:), INTENT(IN)  :: x, xs
        REAL (KIND=8), DIMENSION(:),   INTENT(OUT) :: xe

        INTEGER       :: i1, i2, k
        REAL (KIND=8) :: uu, du, d, a1, b1, ab, h11, h12, h13, h14


        uu = se(1)

        CALL patch (n, s, uu, du, i1, i2)

        d = s(n)/du   !   d   = one/du

        a1  = uu
        b1  = 1.-uu
        ab  = a1*b1

        h11 = -6.*ab
        h12 = +b1**2-2.*ab
        h13 = +a1**2-2.*ab
        h14 = +6.*ab

        DO k=1,SIZE(x,1)
           xe(k) = d*( x(k,i1)*h11 + du*xs(k,i1)*h12 + du*xs(k,i2)*h13 + x(k,i2)*h14 )
        ENDDO

        END SUBROUTINE c_1

  END SUBROUTINE evalc_1


  
  SUBROUTINE evalc_2 (i, se, xe,cur)

  ! Evaluate the 2nd-derivative 
  ! of the geometry g
  !---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i
  REAL (KIND=8), DIMENSION(:), INTENT(IN)  :: se
  REAL (KIND=8), DIMENSION(:), INTENT(OUT) :: xe
  TYPE(curve),   DIMENSION(:), INTENT(INOUT):: cur
  !---------------------------------------------------------------------

  CALL c_2 (cur(i) % ns, cur(i) % s, cur(i) % x, cur(i) % xs, se, xe)

  CONTAINS

        SUBROUTINE c_2 (n, s, x, xs, se, xe)

        ! Evaluate hermite univariate interpolation at point si
        ! si ranges from zero to one
        ! second derivative xi = (D/Ds)^2 (xi(si)) w.r.t. s

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        REAL (KIND=8), DIMENSION(:),   INTENT(IN)  :: s, se
        REAL (KIND=8), DIMENSION(:,:), INTENT(IN)  :: x, xs
        REAL (KIND=8), DIMENSION(:),   INTENT(OUT) :: xe

        INTEGER       :: i1, i2, k
        REAL (KIND=8) :: uu, du, d, a1, b1, h21, h22, h23, h24


        uu = se(1)

        CALL patch (n, s, uu, du, i1, i2)

        d = s(n)/du   !   d   = one/du
        d = d*d

        a1  = uu
        b1  = 1.-uu

        h21 = 6.*(a1-b1)
        h22 = 2.*a1-4.*b1
        h23 = 4.*a1-2.*b1
        h24 = 6.*(b1-a1)

        DO k=1,SIZE(x,1)
           xe(k) = d*( x(k,i1)*h21 + du*xs(k,i1)*h22 + du*xs(k,i2)*h23 + x(k,i2)*h24 )
        ENDDO

        END SUBROUTINE c_2

  END SUBROUTINE evalc_2


  
  SUBROUTINE evalc_3 (i, se, xe,cur)

  ! Evaluate the 3rd-derivative of the geometry g
  !---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i
  REAL (KIND=8), DIMENSION(:), INTENT(IN)  :: se
  REAL (KIND=8), DIMENSION(:), INTENT(OUT) :: xe
  TYPE(curve),   DIMENSION(:), INTENT(INOUT):: cur
  !---------------------------------------------------------------------
  
  CALL c_3 (cur(i) % ns, cur(i) % s, cur(i) % x, cur(i) % xs, se, xe)

  CONTAINS

        SUBROUTINE c_3 (n, s, x, xs, se, xe)

        ! Evaluate hermite univariate interpolation at point si
        ! si ranges from zero to one
        ! third derivative xi = (D/Ds)^3 (xi(si)) w.r.t. s

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        REAL (KIND=8), DIMENSION(:),   INTENT(IN)  :: s, se
        REAL (KIND=8), DIMENSION(:,:), INTENT(IN)  :: x, xs
        REAL (KIND=8), DIMENSION(:),   INTENT(OUT) :: xe

        INTEGER       :: i1, i2, k
        REAL (KIND=8) :: uu, du, d, h31, h32, h33, h34


        uu = se(1)

        CALL patch (n, s, uu, du, i1, i2)

        d = s(n)/du   !   d   = one/du
        d = d*d*d

        h31 = +12.
        h32 = +6.
        h33 = +6.
        h34 = -12.

        DO k=1,SIZE(x,1)
           xe(k) = d*( x(k,i1)*h31 + du*xs(k,i1)*h32 + du*xs(k,i2)*h33 + x(k,i2)*h34 )
        ENDDO

        END SUBROUTINE c_3

  END SUBROUTINE evalc_3



  SUBROUTINE evalc_0123 (i, se, x0e, x1e, x2e, x3e,cur)
  
  ! Evaluate the 0th,1st,2nd,3rd-derivatives 
  ! of the geometry g
  !---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i
  REAL (KIND=8), DIMENSION(:), INTENT(IN)  :: se
  REAL (KIND=8), DIMENSION(:), INTENT(OUT) :: x0e, x1e, x2e, x3e
  TYPE(curve),   DIMENSION(:), INTENT(INOUT):: cur
  !---------------------------------------------------------------------

  CALL c_0123 (cur(i) % ns, cur(i) % s, cur(i) % x, cur(i) % xs,   &
               se, x0e, x1e, x2e, x3e)

  CONTAINS

        SUBROUTINE c_0123 (n, s, x, xs, se, x0e, x1e, x2e, x3e)

        ! Evaluate hermite univariate interpolations at point si

        ! si ranges from zero to one
        ! zeroth, first, second and third derivatives w.r.t. s

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        REAL (KIND=8), DIMENSION(:),   INTENT(IN)  :: s, se
        REAL (KIND=8), DIMENSION(:,:), INTENT(IN)  :: x, xs
        REAL (KIND=8), DIMENSION(:),   INTENT(OUT) :: x0e, x1e, x2e, x3e

        INTEGER       :: i1, i2, k
        REAL (KIND=8) :: uu, du, d, e, f, a1, b1, ab, a2, b2,    &
             h01, h02, h03, h04, h11, h12, h13, h14,   &
             h21, h22, h23, h24, h31, h32, h33, h34


        uu = se(1)

        CALL patch (n, s, uu, du, i1, i2)

        d = s(n)/du   !   d   = one/du
        e = d*d
        f = e*d

        a1 = uu
        b1 = 1.-uu
        ab = a1*b1
        a2 = a1**2
        b2 = b1**2

        h01 = (b1+3.*a1)*b2
        h02 = +a1*b2
        h03 = -a2*b1
        h04 = (a1+3.*b1)*a2

        h11 = -6.*ab
        h12 = +b1**2-2.*ab
        h13 = +a1**2-2.*ab
        h14 = +6.*ab

        h21 = 6.*(a1-b1)
        h22 = 2.*a1-4.*b1
        h23 = 4.*a1-2.*b1
        h24 = 6.*(b1-a1)

        h31 = +12.
        h32 = +6.
        h33 = +6.
        h34 = -12.

        DO k=1,SIZE(x,1)
           x0e(k) =     x(k,i1)*h01 + du*xs(k,i1)*h02 + du*xs(k,i2)*h03 + x(k,i2)*h04
           x1e(k) = d*( x(k,i1)*h11 + du*xs(k,i1)*h12 + du*xs(k,i2)*h13 + x(k,i2)*h14 )
           x2e(k) = e*( x(k,i1)*h21 + du*xs(k,i1)*h22 + du*xs(k,i2)*h23 + x(k,i2)*h24 )
           x3e(k) = f*( x(k,i1)*h31 + du*xs(k,i1)*h32 + du*xs(k,i2)*h33 + x(k,i2)*h34 )
        ENDDO

        END SUBROUTINE c_0123

  END SUBROUTINE evalc_0123







  SUBROUTINE new_curve (idf, ascii, g)

  ! Read the knot coordinates and                                                 
  ! compute the interpolant univariate spline                                     
  ! input from file unit idf, to be opened and closed elsewhere                   
  !-------------------------------------------------------------------------      
  IMPLICIT NONE                                                                   

  INTEGER, INTENT(IN) :: idf                                                      
  LOGICAL, INTENT(IN) :: ascii                                                    
  TYPE (curve), INTENT(INOUT) :: g                                                

  INTEGER, PARAMETER :: nat=0, cla=1, per=2                                       
  INTEGER :: md, ms, idc                                                          
  REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: xs1, xsn, a, b, c, p, q             

  EXTERNAL choun, natun, claun                                                    
  !-------------------------------------------------------------------------      

  IF (ascii) THEN 
     READ (idf,*)                                                                 
     READ (idf,*) md, ms 
  ELSE                                                                            
     READ (idf) md, ms                                                             
  ENDIF                                                                           

  ALLOCATE (g % s(ms))                                                            
  ALLOCATE (g % x(md,ms))                                                         
  ALLOCATE (g % xs(md,ms))                                                        

  ALLOCATE (xs1(md))                                                              
  ALLOCATE (xsn(md))                                                              

  g % nd = md                                                                     
  g % ns = ms                                                                     

  IF (ascii) THEN                                                                 
     CALL alcond (idf, g % nd, idc, xs1, xsn)                                     
     CALL alpnts (idf, g % nd, g % ns, g % x)                                     
  ELSE                                                                            
     CALL blcond (idf, g % nd, idc, xs1, xsn)                                     
     CALL blpnts (idf, g % nd, g % ns, g % x)                                     
  ENDIF                                                                           

  ALLOCATE (a(ms))                                                                
  ALLOCATE (b(ms))                                                                
  ALLOCATE (c(ms))                                                                

  SELECT CASE (idc)                                                               

  CASE (nat)                                                                      

     CALL mkun (natunp, g % nd, g % ns, g % s, g % x, xs1, xsn, a, b, c, g % xs)  

  CASE (cla)                                                                      

     CALL mkun (claunp, g % nd, g % ns, g % s, g % x, xs1, xsn, a, b, c, g % xs)  

  CASE (per)                                                                      

     ALLOCATE (p(ms))                                                             
     ALLOCATE (q(ms))                                                             

     CALL mkunp (g % nd, g % ns, g % s, g % x, a, b, c, p, q, g % xs)             

     DEALLOCATE (q)                                                               
     DEALLOCATE (p)                                                               

  CASE DEFAULT                                                                    

     WRITE (*,*) ''                                                               
     WRITE (*,*) 'ERROR. Unknown univariate '                                     
     WRITE (*,*) 'spline end condition ',idc                                      
     WRITE (*,*) ''                                                               
     STOP                                                                         

  END SELECT                                                                      

  DEALLOCATE (c)                                                                  
  DEALLOCATE (b)                                                                  
  DEALLOCATE (a)                                                                  

  DEALLOCATE (xsn)                                                                
  DEALLOCATE (xs1)                                                                

  END SUBROUTINE new_curve



  SUBROUTINE load_curve (idf, ascii, g)

  ! Load univariate spline
  ! input from file unit idf, to be opened and closed elsewhere
  !-------------------------------------------------------------------------      
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: idf
  LOGICAL, INTENT(IN) :: ascii
  TYPE (curve), INTENT(INOUT) :: g

  INTEGER :: md, ms
  !-------------------------------------------------------------------------      

  IF (ascii) THEN                                            
    READ (idf,*)                                            
    READ (idf,*) md, ms                                      
  ELSE                                                       
    READ (idf) md, ms                                        
  ENDIF                                                      

  ALLOCATE (g % s(ms))                                     
  ALLOCATE (g % x(md,ms))                                  
  ALLOCATE (g % xs(md,ms))                                 

  g % nd = md                                                
  g % ns = ms                                                

  IF (ascii) THEN                                            
    CALL alun (idf, g % nd,  g % ns, g % s, g % x, g % xs)  
  ELSE                                                       
    CALL blun (idf, g % nd,  g % ns, g % s, g % x, g % xs)  
  ENDIF                                                      

  END SUBROUTINE load_curve



  SUBROUTINE save_curve (idf, ascii, g)

  ! Save univariate spline
  ! output on file unit idf, to be opened and closed elsewhere
  !-------------------------------------------------------------------------      
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: idf
  LOGICAL, INTENT(IN) :: ascii
  TYPE (curve), INTENT(INOUT) :: g
  !-------------------------------------------------------------------------      

  IF (ascii) THEN
  
     WRITE (idf,'(2(8x,a2))') 'ND','NS'
     WRITE (idf,'(2i10)') g % nd, g % ns

     CALL asun (idf, g % nd, g % ns, g % s, g % x, g % xs)
  
  ELSE
  
     WRITE (idf) g % nd, g % ns

     CALL bsun (idf, g % nd, g % ns, g % s, g % x, g % xs)
  
  ENDIF

  END SUBROUTINE save_curve





  SUBROUTINE patch (n, s, sx, ds, i1, i2)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL (KIND=8), DIMENSION(:), INTENT(IN)  :: s
    REAL (KIND=8), INTENT(INOUT) :: sx
    REAL (KIND=8), INTENT(OUT) :: ds
    INTEGER, INTENT(OUT) :: i1, i2

    INTEGER :: l, h, r


    IF ( ( sx .LT. zero ) .OR. ( sx .GT. one ) ) THEN

       WRITE (*,*) 'patch warning, parameter out of range'
       WRITE (*,*) 'parameter value = ',sx
       WRITE (*,*) 'parameter has been corrected'

       sx = MIN(one,MAX(zero,sx))

    ENDIF


    l = 1; r = n; sx = s(n)*sx

    DO; IF ( r-l == 1 ) EXIT
       h = (l+r)/2
       IF ( s(h) > sx ) THEN
          r = h
       ELSE
          l = h
       ENDIF
    ENDDO

    i1 = l; i2 = r; ds = s(r)-s(l); sx = (sx-s(l))/ds

  END SUBROUTINE patch


  SUBROUTINE mkun (conun,l,n,s,x,xs1,xsn,a,b,c,xs)

    ! Knot derivatives for cubic spline univariate interpolation

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: conun

    INTEGER :: l,n
    REAL (KIND=8) :: s,x,xs1,xsn,a,b,c,xs
    DIMENSION s(n),x(l,n),xs1(l),xsn(l),a(n),b(n),c(n),xs(l,n)

    CALL choun(l,n,x,s)

    IF ( conun == natunp ) then

       call natun(l,n,s,x,a,b,c,xs)

    else if ( conun == claunp ) then

       call claun(l,n,xs1,xsn,a,b,c,xs)

    endif

    CALL sysun  (l,n,s,x,a,b,c,xs)
    CALL luftr  (n,a,b,c)
    CALL soltr  (l,n,a,b,c,xs)

  END SUBROUTINE mkun


  SUBROUTINE choun (l,n,x,s)

    !     Chord parametrization for spline univariate interpolation

    !        input
    !     l       block DIMENSION
    !     n       knot number
    !     x       knot coordinates

    !        output
    !     s       chord abscissa

    IMPLICIT NONE

    INTEGER :: l,n
    REAL (KIND=8) :: x,s
    DIMENSION x(l,n),s(n)

    INTEGER :: i,k
    REAL (KIND=8) :: t


    s(1) = 0.
    DO i=2,n
       t = 0.
       DO k=1,l
          t = t + (x(k,i)-x(k,i-1))**2
       ENDDO
       s(i) = s(i-1) + sqrt(t)
    ENDDO

  END SUBROUTINE choun



  SUBROUTINE natun (l,n,s,f,a,b,c,r)

    !     System definitions for spline univariate interpolation
    !     natural end conditions

    !        input
    !     l       block DIMENSIONs
    !     n       knot number
    !     s       knot parametrization
    !     f       knot function values
    !by DD removed     dummy1  dummy variable
    ! "       "        dummyn  dummy variable

    !        output
    !     a,b,c   coefficients of the matrix
    !     r       right hand side

    IMPLICIT NONE

    INTEGER :: l,n
    REAL (KIND=8) :: s,f,a,b,c,r!,dummy1,dummyn
    DIMENSION s(n),f(l,n),a(n),b(n),c(n),r(l,n)!,dummy1(l),dummyn(l)

    INTEGER :: k


    b(1) = 2.
    c(1) = 1.
    a(n) = 1.
    b(n) = 2.

    DO k=1,l
       r(k,1) = 3.*(f(k,2)-f(k,1  ))/(s(2)-s(1  ))
       r(k,n) = 3.*(f(k,n)-f(k,n-1))/(s(n)-s(n-1))
    ENDDO

  END SUBROUTINE natun


  SUBROUTINE claun (l,n,rs1,rsn,a,b,c,r)

    !     System definitions for spline univariate interpolation
    !     clamped end conditions

    !        input
    !     l       block DIMENSIONs
    !     n       knot number
    !by DD removed    s       knot parametrization
    ! "       "       f       knot function values
    !     rs1     prescribed condition (s-slope) at node 1
    !     rsn     prescribed condition (s-slope) at node n

    !        output
    !     a,b,c   coefficients of the matrix
    !     r       right hand side

    IMPLICIT NONE

    INTEGER :: l,n
    REAL (KIND=8) :: rs1,rsn,a,b,c,r!s,f
    DIMENSION rs1(l),rsn(l),a(n),b(n),c(n),r(l,n)!s(n),f(l,n)

    INTEGER :: k


    b(1) = 1.
    c(1) = 0.
    a(n) = 0.
    b(n) = 1.

    DO k=1,l
       r(k,1) = rs1(k)
       r(k,n) = rsn(k)
    ENDDO

  END SUBROUTINE claun


  SUBROUTINE sysun (l,n,s,f,a,b,c,r)

    !     System definitions for spline univariate interpolation
    !     internal points

    !        input
    !     l       block DIMENSIONs
    !     n       knot number
    !     s       knot parametrization
    !     f       knot function values

    !        output
    !     a,b,c   coefficients of the matrix
    !     r       right hand side

    IMPLICIT NONE

    INTEGER :: l,n
    REAL (KIND=8) :: s,f,a,b,c,r
    DIMENSION s(n),f(l,n),a(n),b(n),c(n),r(l,n)

    INTEGER :: k,i
    REAL (KIND=8) :: fa,fc


    DO i=2,n-1
       a(i) = s(i+1)-s(i)
       c(i) = s(i)-s(i-1)
       b(i) = 2.*(a(i)+c(i))
       DO k=1,l
          fa     = f(k,i+1)-f(k,i)
          fc     = f(k,i)-f(k,i-1)
          r(k,i) = 3.*(fc*a(i)/c(i)+fa*c(i)/a(i))
       ENDDO
    ENDDO

  END SUBROUTINE sysun


  SUBROUTINE mkunp (l,n,s,x,a,b,c,p,q,xs)

    ! Knot derivatives for cubic spline univariate interpolation

    IMPLICIT NONE

    INTEGER :: l,n
    REAL (KIND=8) :: s,x,a,b,c,p,q,xs
    DIMENSION s(n),x(l,n),a(n),b(n),c(n),p(n),q(n),xs(l,n)

    CALL choun  (l,n,x,s)

    CALL sysunp (l,n,s,x,a,b,c,p,q,xs)
    CALL lufar  (n-1,a,b,c,p,q)
    CALL solar  (l,n-1,a,b,c,p,q,xs)
    CALL lastcp (l,xs(1,1),xs(1,n))

  END SUBROUTINE mkunp


  SUBROUTINE sysunp (l,n,s,f,a,b,c,p,q,r)

    !     System definitions for spline univariate interpolation
    !     internal points

    !        input
    !     l          block DIMENSIONs
    !     n          knot number
    !     s          knot parametrization
    !     f          knot function values

    !        output
    !     a,b,c,p,q  coefficients of the matrix
    !     r          right hand side

    IMPLICIT NONE

    INTEGER :: l,n
    REAL (KIND=8) :: s,f,a,b,c,p,q,r
    DIMENSION s(n),f(l,n),a(n),b(n),c(n),p(n),q(n),r(l,n)

    INTEGER :: k,i
    REAL (KIND=8) :: fa,fc


    a(1) = s(2)-s(1)
    c(1) = s(n)-s(n-1)
    b(1) = 2.*(a(1)+c(1))
    DO k=1,l
       fa     = f(k,2)-f(k,1)
       fc     = f(k,n)-f(k,n-1)
       r(k,1) = 3.*(fc*a(1)/c(1)+fa*c(1)/a(1))
    ENDDO

    DO i=2,n-1
       a(i) = s(i+1)-s(i)
       c(i) = s(i)-s(i-1)
       b(i) = 2.*(a(i)+c(i))
       DO k=1,l
          fa     = f(k,i+1)-f(k,i)
          fc     = f(k,i)-f(k,i-1)
          r(k,i) = 3.*(fc*a(i)/c(i)+fa*c(i)/a(i))
       ENDDO
    ENDDO

    q(1) = a(1)
    p(1) = c(n-1)
    DO i=2,n-1
       p(i) = 0.
       q(i) = 0.
    ENDDO

  END SUBROUTINE sysunp


  SUBROUTINE lastcp (n,a,b)

    IMPLICIT NONE

    INTEGER :: n,i
    REAL (KIND=8) :: a(n),b(n)

    DO i=1,n
       b(i) = a(i)
    ENDDO

  END SUBROUTINE lastcp





  SUBROUTINE alcond (unit,m,idc,ca,cb)

    ! Load end conditions; ascii file, free format
    ! input from file unit unit, to be opened and closed elsewhere

    IMPLICIT NONE

    INTEGER :: unit,m,idc
    REAL (KIND=8) :: ca,cb
    DIMENSION ca(m),cb(m)

    INTEGER :: clamp
    parameter (clamp=1)

    INTEGER :: i


    READ (unit,*)
    READ (unit,*) idc
    IF (idc .EQ. clamp) THEN
       READ (unit,*)
       READ (unit,*) (ca(i),i=1,m),(cb(i),i=1,m)
    ENDIF

  END SUBROUTINE alcond


  SUBROUTINE blcond (unit,m,idc,ca,cb)

    ! Load end conditions; binary file
    ! input from file unit unit, to be opened and closed elsewhere

    IMPLICIT NONE

    INTEGER :: unit,m,idc
    REAL (KIND=8) :: ca,cb
    DIMENSION ca(m),cb(m)

    INTEGER :: clamp
    parameter (clamp=1)

    INTEGER :: i


    READ (unit) idc
    IF (idc .EQ. clamp) then
       READ (unit) (ca(i),i=1,m),(cb(i),i=1,m)
    ENDIF

  END SUBROUTINE blcond





  SUBROUTINE alpnts (unit,m,n,x)

    ! Load knot coordinates; ascii file, free format
    ! input from file unit unit, to be opened and closed elsewhere

    IMPLICIT NONE

    INTEGER :: unit,m,n
    REAL (KIND=8) :: x
    DIMENSION x(m,n)

    INTEGER :: i,j


    READ(unit,*)
    DO i = 1, n
       READ(unit,*) (x(j,i),j=1,m)
    ENDDO

  END SUBROUTINE alpnts


  SUBROUTINE blpnts (unit,m,n,x)

    ! Load knot coordinates; binary file
    ! input from file unit unit, to be opened and closed elsewhere

    IMPLICIT NONE

    INTEGER :: unit,m,n
    REAL (KIND=8) :: x
    DIMENSION x(m,n)

    INTEGER :: i,j


    DO i=1,n
       READ (unit) (x(j,i),j=1,m)
    ENDDO

  END SUBROUTINE blpnts





  SUBROUTINE alun  (unit,m,n,s,x,xs)

    ! Load univariate spline; ascii file, free format
    ! input from file unit unit, to be opened and closed elsewhere

    IMPLICIT NONE

    INTEGER :: unit,m,n
    REAL (KIND=8) :: s,x,xs
    DIMENSION s(n),x(m,n),xs(m,n)

    INTEGER :: i,j


    READ  (unit,*)
    DO i=1,n
       READ  (unit,*) s(i),(x(j,i),j=1,m),(xs(j,i),j=1,m)
    ENDDO

  END SUBROUTINE alun


  SUBROUTINE blun  (unit,m,n,s,x,xs)

    ! Load univariate spline; binary file
    ! input from file unit unit, to be opened and closed elsewhere

    IMPLICIT NONE

    INTEGER :: unit,m,n
    REAL (KIND=8) :: s,x,xs
    DIMENSION s(n),x(m,n),xs(m,n)

    INTEGER :: i,j


    DO i=1,n
       READ  (unit) s(i),(x(j,i),j=1,m),(xs(j,i),j=1,m)
    ENDDO

  END SUBROUTINE blun





  SUBROUTINE asun  (unit,m,n,s,x,xs)

    ! Save univariate spline; ascii file (format for up to 3D problems)
    ! output on file unit unit, to be opened and closed elsewhere

    IMPLICIT NONE

    INTEGER :: unit,m,n
    REAL (KIND=8) :: s,x,xs
    DIMENSION s(n),x(m,n),xs(m,n)

    INTEGER :: i,j


    WRITE (unit,'(14x,a1,6(12x,a2,i1))') 'S',(' X',j,j=1,m),('XS',j,j=1,m)
    DO i=1,n
       WRITE (unit,'(7(1x,e15.9))') s(i),(x(j,i),j=1,m),(xs(j,i),j=1,m)
    ENDDO

  END SUBROUTINE asun


  SUBROUTINE bsun  (unit,m,n,s,x,xs)

    ! Save univariate spline; binary file
    ! output on file unit unit, to be opened and closed elsewhere

    IMPLICIT NONE

    INTEGER :: unit,m,n
    REAL (KIND=8) :: s,x,xs
    DIMENSION s(n),x(m,n),xs(m,n)

    INTEGER :: i,j


    DO i=1,n
       WRITE (unit) s(i),(x(j,i),j=1,m),(xs(j,i),j=1,m)
    ENDDO

  END SUBROUTINE bsun


END MODULE curves
