MODULE refine_2d

  USE grid_types     !  Types module

!  USE domain_2d      !  Domain module

  USE boundary_2d    !  Boundary module

  USE delaunay       !  Delaunay basic module

  USE eucl_delaunay  !  Delaunay Euclide module

!  USE riem_delaunay  !  Delaunay Riemann module

  IMPLICIT NONE

  PUBLIC  ::  halve_dom ! refine_boundary, refine_domain

  REAL (KIND=8), PARAMETER, PRIVATE :: ref_length = 1.5d0

CONTAINS

  SUBROUTINE halve_dom (grid, idom, sdom)

    IMPLICIT NONE

    TYPE (domain), INTENT(INOUT) :: grid
    INTEGER, INTENT(IN) :: idom
    TYPE (simplex), TARGET, INTENT(IN) :: sdom

    INTEGER :: i
    REAL (KIND=8) :: factor
    TYPE (point), POINTER :: p


    factor = 1./DBLE(nd_d)

    p  =>  new_point()

    p % cnd = nul_cnd
    p % parent = sdom % index
    p % old  =>  sdom % opp(1) % pnt

    p % x = 0.; p % metric = 0.
    DO i=1,nd_d+1
       IF (i .NE. idom) THEN
          p % x = p % x + sdom % opp(i) % pnt % x
          p % metric = p % metric + sdom % opp(i) % pnt % metric
       ENDIF
    ENDDO
    p % x = factor * p % x; p % metric = factor * p % metric

    CALL pushtop_point (grid % add_head, p)

  END SUBROUTINE halve_dom
  !
  !********************************************************************
  !
END MODULE refine_2d
