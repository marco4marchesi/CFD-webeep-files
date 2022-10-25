MODULE wakes_gen

  USE list
  USE grid_types

  IMPLICIT NONE
  PUBLIC :: insert_between


  CONTAINS


  SUBROUTINE insert_between ( p1, p2, q )
  
  IMPLICIT NONE
  TYPE(point), POINTER :: p1, p2, q

  IF ( .NOT. ASSOCIATED(p1 % next, p2) ) THEN
    WRITE(*,*)'error insert_between:points not consecutive'
    STOP
  END IF
  
  NULLIFY ( p1 % next )
  p1 % next => q
  q % prev => p1
  
  NULLIFY ( p2 % prev )
  p2 % prev => q
  q % next => p2

  END SUBROUTINE insert_between





  SUBROUTINE data_trail ( p1, p2, p3, h )

    IMPLICIT NONE

    TYPE(point) ,POINTER :: p1, p2, p3
    TYPE(point) ,POINTER :: q1, q2, q3
    TYPE(point) ,POINTER :: b1, b2, b3
    TYPE(point) ,POINTER :: h, this_point, next_point
    
    REAL(KIND=8),DIMENSION(2,2):: rot_mat 
    REAL(KIND=8),DIMENSION(2) :: q_med !,qin,qfin,vers
    REAL(KIND=8):: teta, leng, height  !,direction,distance,dist,t! dirin
    REAL(KIND=8):: leng_in, height_in  !,height_fin! par
    REAL(KIND=8),PARAMETER :: pi = 3.1415927464102d0
    real*8 :: mfangle

    INTEGER :: i

    ! the assumption will be that p2 has the same coordinates 
    ! as the first wake point
    p2 % what = 14

    leng_in = p2 % leng
    height_in = p2 % height 

    leng = leng_in
    height = height_in

   !--------------------------------- define a 90 rotation
    teta = 90
    teta = -teta*pi/180.d0 

    rot_mat(1,1) = COS(teta) ! rotazione oraria
    rot_mat(2,1) = SIN(teta) !(convenzione + antioraria
    rot_mat(1,2) = -rot_mat(2,1)
    rot_mat(2,2) =  rot_mat(1,1)
    !---------------------------------
    b1 => p1
    b2 => p2
    b3 => p3

    this_point => h
    ! the first wake point, coincident with the boundary point
    this_point => this_point % next                                     
    next_point => this_point % next 
                                   
    i = 0
    DO 
        i = i+1
       !--------------------------------- define the new points
       q1 => new_point()
       q2 => new_point()
       q3 => new_point()

       !--------------------------------- the pointers
       q1 % old => q2
       q2 % ppp => q1
       q3 % old => q2
       q2 % ppp1 => q3
       !--------------------------------- the carachteristic lengths
       leng = this_point % leng       
       height = this_point % height

       q1 % leng = leng ; q2 % leng = leng; q3 % leng = leng
       q1 % leng_max = leng ; q2 % leng_max = leng; q3 % leng_max = leng
       q1 % height = height ; q2 % height = height; q3 % height = height
       !--------------------------------- define the advancing direction 
       q_med = (next_point % x - this_point % x)

       q_med = q_med/SQRT(DOT_PRODUCT(q_med,q_med))
       !--------------------------------- point coordinates
       q2 % x = next_point % x 

       q1 % x = q2 % x + MATMUL(rot_mat,q_med)* (+height)
       q3 % x = q2 % x + MATMUL(rot_mat,q_med)* (-height)

       q1 % cnd = 0 ;q2 % cnd = 0 ;q3 % cnd = 0

       CALL insert_between ( b1, b3, q1 )
       CALL insert_between ( q1, b3, q3 )

       IF (i == 1) THEN
          b2 % wake => q2
          q2 % old => b2
       ELSEIF((i > 1))THEN
          b2 % next => q2
          q2 % old => b2
       ENDIF

       b1 => q1
       b2 => q2
       b3 => q3

       this_point => this_point % next
       next_point => next_point % next
       IF(ASSOCIATED(next_point,h % prev)) EXIT
!
    ENDDO
!
     NULLIFY ( b2 )
     q1 => new_point()
     q2 => new_point()
     q3 => new_point()
!
!     q1 % old => q2
!     q2 % ppp => q1
!     q3 % old => q2
!     q2 % ppp1 => q3
!
     q1 % old => b1 % old
     q2 % old => b1 % old
     q3 % old => b1 % old     
!
     q1 % cnd = 0; q2 % cnd = 0; q3 % cnd = 0
!    
     mfangle = -45.d0*PI/180.d0
     call rot( mfangle, this_point%x, q_med )
     q_med = q_med / SQRT(DOT_PRODUCT(q_med,q_med))
     q1 % x = this_point % x + q_med * height
    
     CALL insert_between( b1, b3, q1 )
     
!     NULLIFY(b2 % next)   
!     q1 % old => b1 % old
     
     leng = SQRT(DOT_PRODUCT(b1 % x - q1 % x,b1 % x - q1 % x))
     q1 % prev % leng = leng
     
     leng = SQRT(DOT_PRODUCT(b3 % x - q1 % x,b3 % x - q1 % x))
     q1 % leng = leng

     q1 % old % leng_max = ( q1 % prev % leng + q1 % leng )/2
     q1 % height = this_point % height

     mfangle = 0.d0*PI/180.d0
     call rot( mfangle, this_point%x, q_med )
     q_med = q_med / SQRT(DOT_PRODUCT(q_med,q_med))
     q2 % x = this_point % x + q_med * height
    
     CALL insert_between( q1, b3, q2 )

!     NULLIFY(b2 % next)
!     q2 % old => q1 % old
     
     leng = SQRT(DOT_PRODUCT(q1 % x - q2 % x,q1 % x - q2 % x))
     q2 % prev % leng = leng
     
     leng = SQRT(DOT_PRODUCT(b3 % x - q2 % x,b3 % x - q2 % x))
     q2 % leng = leng

     q2 % old % leng_max = ( q2 % prev % leng + q2 % leng )/2
     q2 % height = this_point % height

     mfangle = 45.d0*PI/180.d0
     call rot( mfangle, this_point%x, q_med )
     q_med = q_med / SQRT(DOT_PRODUCT(q_med,q_med))
     q3 % x = this_point % x + q_med * height
    
     CALL insert_between( q2, b3, q3 )

!     NULLIFY(b2 % next)    
!     q3 % old => q2 % old
     
     leng = SQRT(DOT_PRODUCT(q2 % x - q3 % x,q2 % x - q3 % x))
     q3 % prev % leng = leng
     
     leng = SQRT(DOT_PRODUCT(b3 % x - q3 % x,b3 % x - q3 % x))
     q3 % leng = leng

     q3 % old % leng_max = ( q3 % prev % leng + q3 % leng )/2
     q3 % height = this_point % height

     b1 % prev % leng = (2 * b1 % prev % leng / 3)
     b3 % leng = 2 * b3 % leng / 3

  END SUBROUTINE data_trail
  !같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같
  SUBROUTINE rot(theta,vecti,vecto)
    IMPLICIT NONE 
    REAL(KIND=8),INTENT(IN) ::theta
    REAL(KIND=8),DIMENSION(2),INTENT(IN) :: vecti
    REAL(KIND=8),DIMENSION(2),INTENT(OUT) :: vecto
    REAL(KIND=8),DIMENSION(2,2)::  rot_mat

    !convenzione + antioraria,angoli in radianti

    rot_mat(1,1) = COS(theta) 
    rot_mat(2,1) = SIN(theta) 
    rot_mat(1,2) = -rot_mat(2,1)
    rot_mat(2,2) = rot_mat(1,1)

    !vecto = MATMUL(rot_mat,vecti)
    vecto(1) = rot_mat(1,1)*vecti(1) + rot_mat(1,2)*vecti(2) 
    vecto(2) = rot_mat(2,1)*vecti(1) + rot_mat(2,2)*vecti(2) 

  END SUBROUTINE rot

  !같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같같
END MODULE wakes_gen
