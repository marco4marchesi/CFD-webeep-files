!
      module list
!
!     +-------------------------------------------------------------------+
!     |  The following procedures implement list operations               |
!     |  on items which are structures defined in module grid_types       |
!     |  The types presently considered are "point" and "simplex"         |
!     |  The only assumption made about the components of a generic item  |
!     |  is the existence of the following two components:                |
!     |  1) "item % prev" (pointer to the previous item of the list);     |
!     |  2) "item % next" (pointer to the next item of the list).         |
!     +-------------------------------------------------------------------+
!
      use grid_types
!
      implicit none
      public
!
!
      contains
!
!
!     +---------------------------------+
!     |  Section related to type POINT  |
!     +---------------------------------+
!
      function newhead_point () result ( newhead_point_r )
!
!     +--------------------------------------------------+
      implicit none
      type (point), pointer :: newhead_point_r
!     +--------------------------------------------------+

      allocate ( newhead_point_r )
!
      newhead_point_r % next  =>  newhead_point_r
      newhead_point_r % prev  =>  newhead_point_r
!
      end function newhead_point
!
!
!
!
!
   FUNCTION poptop_point (head) RESULT (linked_point)

      IMPLICIT NONE

      TYPE (point), POINTER :: head, linked_point

      linked_point  =>  head % next

   END FUNCTION poptop_point


   FUNCTION popbot_point (head) RESULT (linked_point)

      IMPLICIT NONE

      TYPE (point), POINTER :: head, linked_point

      linked_point  =>  head % prev

   END FUNCTION popbot_point


   SUBROUTINE extract_point (linked_point)

      IMPLICIT NONE

      TYPE (point), POINTER :: linked_point

      TYPE (point), POINTER :: p, n

! i collegamenti del punto estratto con i suoi adiacenti rimangono immutati
      p  =>  linked_point % prev
      n  =>  linked_point % next

      IF ( ASSOCIATED (p) .AND. ASSOCIATED (n) ) THEN
         IF ( ASSOCIATED (p,linked_point) .OR.  &
	      ASSOCIATED (n,linked_point) ) THEN
            WRITE (*,*) 'extract_point error, list is empty'
	    WRITE (*,*) '(stop)'
	    STOP
         ELSE
            p % next  =>  n
            n % prev  =>  p
         END IF
      ELSE
         WRITE (*,*) 'extract_point error, unlinked point'
	 WRITE (*,*) '(stop)'
	 STOP
      ENDIF

   END SUBROUTINE extract_point


   SUBROUTINE pushtop_point (head, unlinked_point)

      IMPLICIT NONE

      TYPE (point), POINTER :: head, unlinked_point
      TYPE (point), POINTER :: frst


      frst  =>  head % next

      head % next  =>  unlinked_point
      unlinked_point % prev  =>  head

      unlinked_point % next  =>  frst
      frst % prev  =>  unlinked_point

   END SUBROUTINE pushtop_point


   SUBROUTINE pushbot_point (head, unlinked_point)

      IMPLICIT NONE

      TYPE (point), POINTER :: head, unlinked_point
      TYPE (point), POINTER :: last


      last  =>  head % prev

      last % next  =>  unlinked_point
      unlinked_point % prev  =>  last

      unlinked_point % next  =>  head
      head % prev  =>  unlinked_point

   END SUBROUTINE pushbot_point


   SUBROUTINE movetop_point (head, linked_point)

      IMPLICIT NONE

      TYPE (point), POINTER :: head, linked_point


      CALL extract_point (linked_point)
      CALL pushtop_point (head, linked_point)

   END SUBROUTINE movetop_point


   SUBROUTINE movebot_point (head, linked_point)

      IMPLICIT NONE

      TYPE (point), POINTER :: head, linked_point


      CALL extract_point (linked_point)
      CALL pushbot_point (head, linked_point)

   END SUBROUTINE movebot_point


   SUBROUTINE jointop_point (fx_head, mv_head)

      IMPLICIT NONE

      TYPE (point), POINTER :: fx_head, mv_head
      TYPE (point), POINTER :: fx_frst, mv_frst, mv_last


      fx_frst  =>  fx_head % next
      mv_frst  =>  mv_head % next
      mv_last  =>  mv_head % prev

      fx_head % next  =>  mv_frst
      mv_frst % prev  =>  fx_head

      mv_last % next  =>  fx_frst
      fx_frst % prev  =>  mv_last

      mv_head % next  =>  mv_head
      mv_head % prev  =>  mv_head

   END SUBROUTINE jointop_point


   SUBROUTINE joinbot_point (fx_head, mv_head)

      IMPLICIT NONE

      TYPE (point), POINTER :: fx_head, mv_head
      TYPE (point), POINTER :: fx_last, mv_frst, mv_last


      fx_last  =>  fx_head % prev
      mv_frst  =>  mv_head % next
      mv_last  =>  mv_head % prev

      fx_last % next  =>  mv_frst
      mv_frst % prev  =>  fx_last

      mv_last % next  =>  fx_head
      fx_head % prev  =>  mv_last

      mv_head % next  =>  mv_head
      mv_head % prev  =>  mv_head

   END SUBROUTINE joinbot_point


   FUNCTION numberof_point (head) RESULT (numberof_point_r)

      IMPLICIT NONE

      TYPE (point), POINTER :: head
      INTEGER               :: numberof_point_r

      INTEGER n
      TYPE (point), POINTER :: p


      p  =>  head;  n = 0  

      DO
         p  =>  p % next;  n = n+1
         IF ( ASSOCIATED (p,head) ) EXIT
      END DO

      numberof_point_r = n-1

   END FUNCTION numberof_point


!            ######   TYPE (simplex)   ######


   FUNCTION newhead_simplex () RESULT (newhead_simplex_r)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: newhead_simplex_r


      ALLOCATE (newhead_simplex_r)

      newhead_simplex_r % next  =>  newhead_simplex_r
      newhead_simplex_r % prev  =>  newhead_simplex_r

   END FUNCTION newhead_simplex


   FUNCTION poptop_simplex (head) RESULT (linked_simplex)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: head, linked_simplex

      linked_simplex  =>  head % next

   END FUNCTION poptop_simplex


   FUNCTION popbot_simplex (head) RESULT (linked_simplex)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: head, linked_simplex

      linked_simplex  =>  head % prev

   END FUNCTION popbot_simplex


   SUBROUTINE extract_simplex (linked_simplex)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: linked_simplex
      TYPE (simplex), POINTER :: p, n


      p  =>  linked_simplex % prev
      n  =>  linked_simplex % next

      IF ( ASSOCIATED (p) .and. ASSOCIATED (n) ) THEN
         IF ( ASSOCIATED (p,linked_simplex) .OR.  &
	      ASSOCIATED (n,linked_simplex) ) THEN
            WRITE (*,*) 'extract_simplex error, list is empty'
	    WRITE (*,*) '(stop)'
	    STOP
         ELSE
            p % next  =>  n
            n % prev  =>  p
         END IF
      ELSE
         WRITE (*,*) 'extract_simplex error, unlinked simplex'
	 WRITE (*,*) '(stop)'
	 STOP
      ENDIF

   END SUBROUTINE extract_simplex


   SUBROUTINE pushtop_simplex (head, unlinked_simplex)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: head, unlinked_simplex
      TYPE (simplex), POINTER :: frst


      frst  =>  head % next

      head % next  =>  unlinked_simplex
      unlinked_simplex % prev  =>  head

      unlinked_simplex % next  =>  frst
      frst % prev  =>  unlinked_simplex

   END SUBROUTINE pushtop_simplex


   SUBROUTINE pushbot_simplex (head, unlinked_simplex)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: head, unlinked_simplex
      TYPE (simplex), POINTER :: last


      last  =>  head % prev

      last % next  =>  unlinked_simplex
      unlinked_simplex % prev  =>  last

      unlinked_simplex % next  =>  head
      head % prev  =>  unlinked_simplex

   END SUBROUTINE pushbot_simplex


   SUBROUTINE movetop_simplex (head, linked_simplex)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: head, linked_simplex

      CALL extract_simplex (linked_simplex)
      CALL pushtop_simplex (head, linked_simplex)

   END SUBROUTINE movetop_simplex


   SUBROUTINE movebot_simplex (head, linked_simplex)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: head, linked_simplex

      CALL extract_simplex (linked_simplex)
      CALL pushbot_simplex (head, linked_simplex)

   END SUBROUTINE movebot_simplex


   SUBROUTINE jointop_simplex (fx_head, mv_head)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: fx_head, mv_head
      TYPE (simplex), POINTER :: fx_frst, mv_frst, mv_last


      fx_frst  =>  fx_head % next
      mv_frst  =>  mv_head % next
      mv_last  =>  mv_head % prev

      fx_head % next  =>  mv_frst
      mv_frst % prev  =>  fx_head

      mv_last % next  =>  fx_frst
      fx_frst % prev  =>  mv_last

      mv_head % next  =>  mv_head
      mv_head % prev  =>  mv_head

   END SUBROUTINE jointop_simplex


   SUBROUTINE joinbot_simplex (fx_head, mv_head)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: fx_head, mv_head
      TYPE (simplex), POINTER :: fx_last, mv_frst, mv_last


      fx_last  =>  fx_head % prev
      mv_frst  =>  mv_head % next
      mv_last  =>  mv_head % prev

      fx_last % next  =>  mv_frst
      mv_frst % prev  =>  fx_last

      mv_last % next  =>  fx_head
      fx_head % prev  =>  mv_last

      mv_head % next  =>  mv_head
      mv_head % prev  =>  mv_head

   END SUBROUTINE joinbot_simplex


   FUNCTION numberof_simplex (head) RESULT (numberof_simplex_r)

      IMPLICIT NONE

      TYPE (simplex), POINTER :: head
      INTEGER :: numberof_simplex_r

      INTEGER n
      TYPE (simplex), POINTER :: p


      p  =>  head;  n = 0

      DO
         p  =>  p % next;  n = n+1
         IF ( ASSOCIATED (p,head) ) EXIT
      END DO

      numberof_simplex_r = n-1

   END FUNCTION numberof_simplex
!
      end module list
