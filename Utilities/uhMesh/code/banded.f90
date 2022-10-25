MODULE banded

   IMPLICIT NONE

CONTAINS


   SUBROUTINE luftr (n,a,b,c)

   !     lu factorization of a tridiagonal matrix

   !           input
   !     n              order of the matrix
   !     a (2:n)        lower diagonal coefficients
   !     b (1:n)        main diagonal coefficients
   !     c (1:n-1)      upper diagonal coefficients

   !           output
   !     a,b,c          lu factorized matrix to be solved by lustr or soltr

          IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL (KIND=8), DIMENSION(n), INTENT(INOUT) :: a,b,c

      INTEGER :: i


      DO i=2,n
         a(i) = a(i)/b(i-1)
         b(i) = b(i)-a(i)*c(i-1)
      ENDDO

   END SUBROUTINE luftr


   SUBROUTINE lustr (n,a,b,c,x)

   !     solution by forward and backward substitution of
   !     a tridiagonal matrix lu factorized with the SUBROUTINE luftr

   !           input
   !     n           order of the matrix
   !     a,b,c       lu factorized tridiagonal matrix produced by luftri
   !     x           right hand side

   !           output
   !     x           solution

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL (KIND=8), DIMENSION(n), INTENT(IN) :: a,b,c
      REAL (KIND=8), DIMENSION(n), INTENT(INOUT) :: x

      INTEGER :: i


      !     forward sweep

      DO i=2,n
         x(i) = x(i)-a(i)*x(i-1)
      ENDDO

      !     backward sweep

      x(n) = x(n)/b(n)

      DO i=n-1,1,-1
         x(i) = (x(i)-c(i)*x(i+1))/b(i)
      ENDDO

   END SUBROUTINE lustr


   SUBROUTINE soltr  (l,n,a,b,c,x)

   !     Block solution of a scalar lu factorized tridiagonal matrix

   !           input
   !     l           block DIMENSION
   !     n           order of the matrix
   !     a,b,c       lu factorized matrix produced by luftr
   !     x           right hand side

   !           output
   !     x           solution

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l,n
      REAL (KIND=8), DIMENSION(n), INTENT(IN) :: a,b,c
      REAL (KIND=8), DIMENSION(l,n), INTENT(INOUT) :: x

      INTEGER :: k,i


      DO i=2,n
         DO k=1,l
            x(k,i) = x(k,i)-a(i)*x(k,i-1)
         ENDDO
      ENDDO

      DO k=1,l
         x(k,n) = x(k,n)/b(n)
      ENDDO

      DO i=n-1,1,-1
         DO k=1,l
            x(k,i) = (x(k,i)-c(i)*x(k,i+1))/b(i)
         ENDDO
      ENDDO

   END SUBROUTINE soltr


   SUBROUTINE soltri (l,nu,nv,j,a,b,c,x)

   !     Block solution of a scalar lu factorized tridiagonal matrix
   !     for a j = constant row

   !        input
   !     l        block DIMENSION
   !     nu       knot row number
   !     nv       knot column number
   !     j        fixed row number
   !     a,b,c    lu factorized matrix produced by luftr
   !     x        right hand side

   !        output
   !     x        solution

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l,nu,nv,j
      REAL (KIND=8), DIMENSION(nu), INTENT(IN) :: a,b,c
      REAL (KIND=8), DIMENSION(l,nu,nv), INTENT(INOUT) :: x

      INTEGER :: i,k


      DO i=2,nu
         DO k=1,l
            x(k,i,j) = x(k,i,j)-a(i)*x(k,i-1,j)
         ENDDO
      ENDDO

      DO k=1,l
         x(k,nu,j) = x(k,nu,j)/b(nu)
      ENDDO

      DO i=nu-1,1,-1
         DO k=1,l
            x(k,i,j) = (x(k,i,j)-c(i)*x(k,i+1,j))/b(i)
         ENDDO
      ENDDO

   END SUBROUTINE soltri


   SUBROUTINE soltrj (l,nu,nv,i,a,b,c,x)

   !     Block solution of a scalar lu factorized tridiagonal matrix
   !     for a i = constant column

   !        input
   !     l        block DIMENSION
   !     nu       knot row number
   !     nv       knot column number
   !     i        fixed column number
   !     a,b,c    lu factorized matrix produced by luftr
   !     x        right hand side

   !        output
   !     x        solution

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l,nu,nv,i
      REAL (KIND=8), DIMENSION(nv), INTENT(IN) :: a,b,c
      REAL (KIND=8), DIMENSION(l,nu,nv), INTENT(INOUT) :: x

      INTEGER :: j,k


      DO j=2,nv
         DO k=1,l
            x(k,i,j) = x(k,i,j)-a(j)*x(k,i,j-1)
         ENDDO
      ENDDO

      DO k=1,l
         x(k,i,nv) = x(k,i,nv)/b(nv)
      ENDDO

      DO j=nv-1,1,-1
         DO k=1,l
            x(k,i,j) = (x(k,i,j)-c(j)*x(k,i,j+1))/b(j)
         ENDDO
      ENDDO

   END SUBROUTINE soltrj


   SUBROUTINE lufar (n,a,b,c,p,q)

   !     lu factorization of an arrow matrix

   !           input
   !     n              order of the matrix
   !     a (2:n)        lower diagonal coefficients
   !     b (1:n)        main diagonal coefficients
   !     c (1:n-1)      upper diagonal coefficients
   !     p (1:n-2)      bottom row coefficients
   !     q (1:n-2)      right column coefficients

   !           output
   !     a,b,c,p,q      lu factorized matrix to be solved by lusar or solar

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL (KIND=8), DIMENSION(n), INTENT(INOUT) :: a,b,c,p,q

      INTEGER :: i


      DO i=2,n-2
         a(i) = a(i)/b(i-1)
         b(i) = b(i)-a(i)*c(i-1)
         q(i) = q(i)-a(i)*q(i-1)
      ENDDO

      a(n-1) = a(n-1)/b(n-2)
      b(n-1) = b(n-1)-a(n-1)*c(n-2)
      c(n-1) = c(n-1)-a(n-1)*q(n-2)

      p(1) = p(1)/b(1)
      DO i=2,n-2
         p(i) = (p(i)-p(i-1)*c(i-1))/b(i)
      ENDDO

      a(n) = (a(n)-p(n-2)*c(n-2))/b(n-1)
      DO i=1,n-2
         b(n) = b(n)-p(i)*q(i)
      ENDDO
      b(n) = b(n)-a(n)*c(n-1)

   END SUBROUTINE lufar


   SUBROUTINE lusar (n,a,b,c,p,q,x)

   !     solution by forward and backward substitution of
   !     an arrow matrix lu factorized with the SUBROUTINE lufar

   !           input
   !     n              order of the matrix
   !     a,b,c,p,q      lu factorized arrow matrix produced by lufar
   !     x              right hand side

   !           output
   !     x              solution

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL (KIND=8), DIMENSION(n), INTENT(IN) :: a,b,c,p,q
      REAL (KIND=8), DIMENSION(n), INTENT(INOUT) :: x

      INTEGER :: i


      !     forward sweep

      DO i=2,n
         x(i) = x(i)-a(i)*x(i-1)
      ENDDO

      DO i=1,n-2
         x(n) = x(n)-p(i)*x(i)
      ENDDO

      !     backward sweep

      x(n) = x(n)/b(n)

      x(n-1) = (x(n-1)-c(n-1)*x(n))/b(n-1)

      DO i=n-2,1,-1
         x(i) = (x(i)-c(i)*x(i+1)-q(i)*x(n))/b(i)
      ENDDO

   END SUBROUTINE lusar


   SUBROUTINE solar  (l,n,a,b,c,p,q,x)

   !     Block solution of a scalar lu factorized arrow matrix

   !           input
   !     l              block DIMENSION
   !     n              order of the matrix
   !     a,b,c,p,q      lu factorized arrow matrix produced by lufar
   !     x              right hand side

   !           output
   !     x              solution

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l,n
      REAL (KIND=8), DIMENSION(n), INTENT(IN) :: a,b,c,p,q
      REAL (KIND=8), DIMENSION(l,n), INTENT(INOUT) :: x

      INTEGER :: k,i


      !     forward sweep

      DO i=2,n
         DO k=1,l
            x(k,i) = x(k,i)-a(i)*x(k,i-1)
         ENDDO
      ENDDO

      DO i=1,n-2
         DO k=1,l
            x(k,n) = x(k,n)-p(i)*x(k,i)
         ENDDO
      ENDDO

      !     backward sweep

      DO k=1,l
         x(k,n) = x(k,n)/b(n)
      ENDDO

      DO k=1,l
         x(k,n-1) = (x(k,n-1)-c(n-1)*x(k,n))/b(n-1)
      ENDDO

      DO i=n-2,1,-1
         DO k=1,l
            x(k,i) = (x(k,i)-c(i)*x(k,i+1)-q(i)*x(k,n))/b(i)
         ENDDO
      ENDDO

   END SUBROUTINE solar


END MODULE banded
