MODULE linsys
 

   IMPLICIT NONE


CONTAINS


   SUBROUTINE lufact (n, pivot, a, det)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      INTEGER, DIMENSION(:), INTENT(INOUT) :: pivot
      REAL (KIND=8), DIMENSION(:,:), INTENT(INOUT) :: a
      REAL (KIND=8), INTENT(OUT) :: det

      INTEGER :: nm1, i, j, k, kp1, m
      REAL (KIND=8) :: p, t


      det = 1.
      pivot(n) = 1

      IF (n.gt.1) THEN
         nm1 = n-1
         DO k=1,nm1
            kp1 = k+1

            ! find pivot p

            m = k
            DO i=kp1,n
               if (abs(a(i,k)).gt.abs(a(m,k))) m=i
            ENDDO
            pivot(k) = m
            IF (m.ne.k) pivot(n) = -pivot(n)
            p = a(m,k)
            a(m,k) = a(k,k)
            a(k,k) = p
            det = det*p

            IF (p.ne.0.0) THEN

               ! compute multipliers

               DO i=kp1,n
                  a(i,k) = -a(i,k)/p
               ENDDO

               ! interchange and eliminate by columns

               DO j=kp1,n
                  t = a(m,j)
                  a(m,j) = a(k,j)
                  a(k,j) = t
                  IF (t.ne.0.0) THEN
                     DO i=kp1,n
                        a(i,j) = a(i,j) + a(i,k)*t
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF

      det = det*a(n,n)*float(pivot(n))

   END SUBROUTINE lufact


   SUBROUTINE lusolv (n, pivot, a, b)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      INTEGER, DIMENSION(:), INTENT(IN) :: pivot
      REAL (KIND=8), DIMENSION(:,:), INTENT(IN) :: a
      REAL (KIND=8), DIMENSION(:), INTENT(INOUT) :: b

      INTEGER :: nm1, k, kb, kp1, km1, m, i
      REAL (KIND=8) :: s


      IF (n.gt.1) THEN

         ! forward elimination

         nm1 = n-1
         DO k=1,nm1
            kp1 = k+1
            m = pivot(k)
            s = b(m)
            b(m) = b(k)
            b(k) = s
            DO i=kp1,n
               b(i) = b(i) + a(i,k)*s
            ENDDO
         ENDDO

         ! back substitution

         DO kb=1,nm1
            km1 = n-kb
            k = km1+1
            b(k) = b(k)/a(k,k)
            s = -b(k)
            DO i=1,km1
               b(i) = b(i) + a(i,k)*s
            ENDDO
         ENDDO
      ENDIF

      b(1) = b(1)/a(1,1)

   END SUBROUTINE lusolv


END MODULE linsys
