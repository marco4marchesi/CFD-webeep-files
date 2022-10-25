!
      module grid_1d
!
!     +--------------------------------------------------------------------------------+
      use dynamic_vector
      implicit none

      type grid_type

        integer :: sd
        integer :: Nb
        
        ! nodes
        integer ::  nj_d, nj_b
        real(kind=8), dimension(:,:), pointer :: rr
        integer,      dimension(:),   pointer :: jd_jb
        integer,      dimension(:),   pointer :: bound_j
        
        ! elements
        integer ::  nm_d
        type (d_i_v), dimension(:),   pointer :: j_m_d
        type (d_i_v), dimension(:),   pointer :: m_j_d
        type (d_i_v), dimension(:),   pointer :: ma_m_d
        integer,      dimension(:),   pointer :: ele_type_d
        
        ! boundary elements
        integer ::  nm_b
        type (d_i_v), dimension(:),   pointer :: j_m_b
        type (d_i_v), dimension(:),   pointer :: m_j_b
        type (d_i_v), dimension(:),   pointer :: ma_m_b
        integer,      dimension(:),   pointer :: ele_type_b
        integer,      dimension(:),   pointer :: bound_m
        
        ! egdes
        integer                               :: nc_d, nc_b
        integer,      dimension(:,:), pointer :: j_c
        type (d_i_v), dimension(:),   pointer :: c_j
        type (d_i_v), dimension(:),   pointer :: m_c
        type (d_i_v), dimension(:),   pointer :: c_m
        integer,      dimension(:,:), pointer :: j_c_b
        integer,      dimension(:),   pointer :: cd_cb
        integer,      dimension(:),   pointer :: cb_cd
        integer,      dimension(:),   pointer :: bound_c

         ! metrics
         integer,      dimension(:,:), pointer :: jcd_fem
         integer,      dimension(:,:), pointer :: jcb_fem
         real(kind=8), dimension(:),   pointer :: cell
         real(kind=8), dimension(:,:), pointer :: eta
         real(kind=8), dimension(:,:), pointer :: chi_b
         real(kind=8), dimension(:,:), pointer :: xi_bp
        
      end type grid_type
!
      real(kind=8), dimension(:), allocatable :: x_L, x_R
      real(kind=8), dimension(:), allocatable :: metric_1
      real(kind=8), dimension(:), allocatable :: metric_2
!
      integer, dimension(:), allocatable :: Nm
      integer, dimension(:), allocatable :: str_type
      integer :: Ns, extr_dir
!     +--------------------------------------------------------------------------------+
!
!
      contains
!
!
!      subroutine read_path_1D(idf_cfg)
!     +--------------------------------------------------------------------------------+
!      implicit none
!      integer, intent(in) :: idf_cfg
!      integer :: i
!     +--------------------------------------------------------------------------------+
!
!      read (idf_cfg,*); read (idf_cfg,*); read (idf_cfg,*)
!      read (idf_cfg,*) extr_dir
!      read (idf_cfg,*) Ns
!
!      allocate ( x_L(Ns), x_R(Ns) )
!      allocate ( metric_1(Ns), metric_2(Ns) )
!      allocate ( nm(Ns), str_type(Ns) )
! 
!      do i = 1, Ns
!
!        read (idf_cfg,*); read(idf_cfg,*); read(idf_cfg,*)
!        read (idf_cfg,*) x_L(i), x_R(i)
!        read (idf_cfg,*) nm(i)
!        read (idf_cfg,*) str_type(i)
!
!        if ( str_type(i) == 2 ) then
!          read(idf_cfg,*) metric_1(i), metric_2(i)
!        else
!          read(idf_cfg,*) metric_1(i)
!        endif
!
!      enddo
! 
!      end subroutine read_path_1D
!
!
!
!
!
      subroutine set_path_1D
!     +--------------------------------------------------------------------------------+
      implicit none
      integer :: i
!     +--------------------------------------------------------------------------------+
!
      extr_dir = 3
      Ns = 1              
!
      allocate ( x_L(Ns), x_R(Ns) )
      allocate ( metric_1(Ns), metric_2(Ns) )
      allocate ( nm(Ns), str_type(Ns) )
! 
      x_L(1) = 0.d0
      x_R(1) = 0.1d0
      nm(1) = 1
      str_type(1) = 0
      metric_1(1) = 0.05
!
      end subroutine set_path_1D


 
  
      subroutine build_grid_1D ( name, option, path )
!     +--------------------------------------------------------------------------------+
      implicit none
!
      character(len=64), intent(in) :: name
      integer, intent(in) :: option
      real*8, dimension(:), pointer, intent(inout), optional :: path
!
      type (grid_type) :: grid
!
      integer, parameter :: MESH = 0
      integer, parameter :: EXPATH = 1
!     +--------------------------------------------------------------------------------+
!
      grid % sd = 1
      grid % Nb = 2
!
      call Create_Nodes( grid )
!
      if ( option .eq. EXPATH ) then
!
        allocate ( path(grid % Nj_d) )
        path = grid % rr(1,:)
        return
!
      else if ( option .eq. MESH ) then
!
        call Create_Elements( grid )
        call Write_Nodes( grid, name )
        call Write_Elements( grid, name )
!
        write(*,*) ' --------------------------------------------------------------------'  
        write(*,*) '   One dimensional grid generated'
        write(*,*) '   - Nodes', grid % Nj_d
        write(*,*) '   - Elements', grid % Nm_d
!
      end if
!
!
      contains

 
      subroutine  Create_Nodes ( grid )
!     +--------------------------------------------------------------------------------+
      implicit none
!
      type(grid_type), intent(inout) :: grid
     
      real*8, dimension(:), allocatable :: length
      real*8, dimension(:), allocatable :: dl
      real*8 :: csi, ss, ss1, ss2
      real*8 :: t, dt
!
      integer :: i, j, k, Nj
      integer, parameter :: UNIFORM      = 0
      integer, parameter :: STRETCH_BEG  = 1
      integer, parameter :: STRETCH_BOTH = 2
      integer, parameter :: STRETCH_END  = 3
!     +--------------------------------------------------------------------------------+
!
      grid % Nj_d = SUM(Nm) + 1
      allocate ( grid % rr(grid % sd, grid % Nj_d) )
      grid % rr(:,1) = x_L(1)
!
      allocate ( length(Ns) )
      length = 0.d0
      allocate ( dl(Ns) )
      dl = 0.d0

      ! Loop over the segments in which the domain 
      ! has been divided into.
      Nj = 1
      do i = 1, Ns
!
        length(i) = ABS(x_R(i) - x_L(i))
        dt = 1.d0 / FLOAT(Nm(i))
!
        select case ( str_type(i) )
!
        case ( UNIFORM )
!
          dl(i) = length(i) / Nm(i)
          k = 2
!
          do j = Nj + 1, Nj + Nm(i) - 1
            grid % rr(:,j) = x_L(i) + (k-1)*dl(i)
            k = k + 1
          end do
!
          grid % rr(:,Nj + Nm(i)) = x_R(i)
          Nj = Nj + Nm(i)
!          
        case ( STRETCH_BEG )
!
          ss = FLOAT(Nm(i)) * metric_1(i) / length(i)
          do j = 1, Nm(i) - 1
!
            k = Nj + j
            t  = FLOAT(j) * dt
            call Stretch1(ss, t, csi)
!
            grid % rr(:,k) = x_L(i) + csi * ABS(x_R(i) - x_L(i))
!
          end do
!
          grid % rr(:,Nj + Nm(i)) = x_R(i)
          Nj = Nj + Nm(i)
!
        case ( STRETCH_BOTH )
!        
          ss1 = FLOAT(Nm(i)) * metric_1(i) / length(i)
          ss2 = FLOAT(Nm(i)) * metric_2(i) / length(i)
!
          do j = 1, Nm(i) - 1
!
            k = Nj + j
            t  = FLOAT(j) * dt
            call Stretch2(ss1, ss2, t, csi)
!
            grid % rr(:,k) = x_L(i) + csi * ABS(x_R(i) - x_L(i))
!
          end do
!
          grid % rr(:,Nj + Nm(i)) = x_R(i)
          Nj = Nj + Nm(i)
!
        case ( STRETCH_END )
!
          ss = FLOAT(Nm(i)) * metric_1(i) / length(i)
          grid % rr(:, Nj + Nm(i)) = x_R(i)

          do j = 1, Nm(i) - 1
!
            k = Nj + Nm(i) - j         
            t  = FLOAT(j) * dt
            CALL Stretch1(ss, t, csi)
!
            grid % rr(:,k) = x_R(i) - csi * ABS(x_R(i) - x_L(i))
!
          end do 
!
          grid % rr(:,Nj) = x_L(i)
          Nj = Nj + Nm(i)
!
        end select
      end do
!
      ! Boundary nodes connectivity
      grid % Nj_b = 2
!
      allocate ( grid % jd_jb(grid % Nj_b) )   
      grid % jd_jb(1) = 1
      grid % jd_jb(2) = grid % Nj_d
!
      allocate ( grid % bound_j(grid % Nb) )
      grid % bound_j(1) = 1
      grid % bound_j(2) = 2
!
      end subroutine  Create_Nodes



        SUBROUTINE  Create_Elements(grid)
        !-----------------------------------------------------------------------
        IMPLICIT NONE

        TYPE(grid_type), INTENT(INOUT) :: grid

        INTEGER :: m
        !-----------------------------------------------------------------------   
        
        grid % Nm_d = SUM(Nm)

        ALLOCATE (grid % j_m_d(grid % Nm_d),       &
                  grid % ele_type_d(grid % Nm_d),  &
                  grid % ma_m_d(grid % Nm_d))
        
        DO m = 1, grid % Nm_d
          
          ALLOCATE (grid % j_m_d(m) % vec(2))
          
          grid % j_m_d(m) % vec(1) = m
          grid % j_m_d(m) % vec(2) = m + 1

          grid % ele_type_d(m) = 1

          ALLOCATE (grid % ma_m_d(m) % vec(2))

          IF (m == 1) THEN

            grid % ma_m_d(m) % vec(1) = 0     
            grid % ma_m_d(m) % vec(2) = m + 1
          
          ELSEIF (m == grid % Nm_d) THEN
          
            grid % ma_m_d(m) % vec(1) = m - 1
            grid % ma_m_d(m) % vec(2) = 0
          
          ELSE

            grid % ma_m_d(m) % vec(1) = m - 1
            grid % ma_m_d(m) % vec(2) = m + 1
                 
          ENDIF

        ENDDO

        ! Boundary Elements
        grid % Nm_b = 0
        
        END SUBROUTINE  Create_Elements
 
      END SUBROUTINE build_grid_1D



        SUBROUTINE Stretch1(S, T, CSI)
        !
        ! One sided stretching function. Vinokur (NASA Report 3313, 1980)
        !
        ! INPUT: s, t
        ! OUTPUT: csi
        ! 
        !-----------------------------------------------------------------------
        IMPLICIT NONE

        REAL(KIND=8):: S, T, CSI
        REAL(KIND=8):: PI, ONE, HLF, P, R, V, W
        !-----------------------------------------------------------------------
   
        PI  = 4.d00 * ATAN(1.0)
        HLF = 0.5d00
        ONE = 1.d00

        IF (S.LT.0.9999d00) THEN
   
          ! invert S=sin(2p)/2p to get p
          IF (S.LE.0.26938972d00) THEN
            P = HLF*PI*(ONE - S + S*S - (ONE-PI*PI/6.d00)*S*S*S  &
                + 6.794732d00*S*S*S*S - 13.205501d00*S*S*S*S*S   &
                + 11.726095d00*S*S*S*S*S*S)
          ELSE
            R = ONE - S
            P = HLF*SQRT(6.d00*R)*(ONE + 0.15d00*R + 0.057321429d00*R*R  &
                + 0.048774238d00*R*R*R - 0.053337753d00*R*R*R*R          &
                + 0.075845134d00*R*R*R*R*R)
          ENDIF
          
          ! Stretching function
          CSI = ONE + ATAN( (T-ONE)*TAN(P) ) / P               

        ELSE IF (S.GT.1.0001d00) THEN
   
          ! invert S=sinh(2p)/2p to get p
          IF (S.GE.2.7829628d00) THEN
            V = LOG(S)
            W = ONE/S - 0.028527431d00
            P = HLF*(V + (ONE+ONE/V)*LOG(2.d00*V) - 0.02041793d00  &
                + 0.24902722d00*W + 1.9496443d00*W*W               &
                - 2.6294547d00*W*W*W + 8.56795911d00*W*W*W*W)
          ELSE
            R = S - ONE
            P = HLF*SQRT(6.d00*R)*(ONE - 0.15d00*R + 0.057321429d00*R*R  &
                - 0.024907295d00*R*R*R + 0.0077424461d00*R*R*R*R          &
                - 0.0010794123d00*R*R*R*R*R)
          ENDIF
          
          ! Stretching function
          V = EXP(P)
          W = EXP(-P)
          R = (T-ONE)*(V - W)/(V + W)
          CSI = ONE + HLF * LOG( (ONE+R)/(ONE-R) ) / P               

        ELSE
   
          ! Stretching function near S=1.
          CSI = T*( ONE + HLF*(S-ONE)*(ONE-T)*(2.d00-T) )

        ENDIF


        END SUBROUTINE Stretch1





        SUBROUTINE Stretch2(SA,SB,T,CSI)
        !
        ! Two sided stretching function. Vinokur (NASA Report 3313, 1980)
        !
        ! INPUT: sa, sb, t
        ! OUTPUT: csi
        ! 
        !-----------------------------------------------------------------------
        IMPLICIT NONE

        REAL(KIND=8):: SA, SB, T, CSI
        REAL(KIND=8):: PI, ONE, HLF, AINV, B, P, R, U, V, W
        !-----------------------------------------------------------------------
   
        PI  = 4.d00 * ATAN(1.0)
        HLF = 0.5d00
        ONE = 1.d00

        AINV = SQRT(SB/SA)
        B    = SQRT(SA*SB)

        U = T/(AINV + (ONE-AINV)*T)

        IF (B.LT.0.9999d00) THEN
   
          ! invert B=sin(p)/p to get p
          IF (B.LE.0.26938972d00) THEN
            P = PI*(ONE - B + B*B - (ONE-PI*PI/6.d00)*B*B*B     &
                + 6.794732d00*B*B*B*B - 13.205501d00*B*B*B*B*B  &
                + 11.726095d00*B*B*B*B*B*B)
          ELSE
            R = ONE - B
            P = SQRT(6.d00*R)*(ONE + 0.15d00*R + 0.057321429d00*R*R  &
                + 0.048774238d00*R*R*R - 0.053337753d00*R*R*R*R      &
                + 0.075845134d00*R*R*R*R*R)
          ENDIF
          
          ! Stretching function
          CSI = HLF + ATAN( (2.d00*U-ONE)*TAN(HLF*P) ) / P               

        ELSE IF (B.GT.1.0001d00) THEN
   
          ! invert B=sinh(p)/p to get p
          IF (B.GE.2.7829628d00) THEN
            V = LOG(B)
            W = ONE/B - 0.028527431d00
            P = V + (ONE+ONE/V)*LOG(2.d00*V) - 0.02041793d00  &
                + 0.24902722d00*W + 1.9496443d00*W*W          &
                - 2.6294547d00*W*W*W + 8.56795911d00*W*W*W*W
          ELSE
            R = B - ONE
            P = SQRT(6.d00*R)*(ONE - 0.15d00*R + 0.057321429d00*R*R  &
                - 0.024907295d00*R*R*R + 0.0077424461d00*R*R*R*R      &
                - 0.0010794123d00*R*R*R*R*R)
          ENDIF
          
          ! Stretching function
          V = EXP(HLF*P)
          W = EXP(-HLF*P)
          R = (2.d00*U-ONE)*(V - W)/(V + W)
          CSI = HLF + LOG( (ONE+R)/(ONE-R) ) / P               

        ELSE
   
          ! Stretching function near B=1.
          CSI = U*( ONE - 2.d00*(B-ONE)*(U-HLF)*(ONE-U) )

        ENDIF

        END SUBROUTINE Stretch2



        SUBROUTINE  Write_Nodes(grid, grid_name)
        !-----------------------------------------------------------------------
        IMPLICIT NONE
   
        TYPE(grid_type),   INTENT(IN) :: grid
        CHARACTER(LEN=64), INTENT(IN) :: grid_name
   
        INTEGER :: j, funit = 11
        !----------------------------------------------------------------------- 

         
        OPEN (UNIT= funit, FILE= 'nodes.'//TRIM(grid_name), FORM= 'formatted', STATUS= 'unknown')
   
           WRITE(funit,1000)
           WRITE(funit,1010) grid_name

           WRITE(funit,1000)
           WRITE(funit,1020)
           WRITE(funit,1021) grid % sd, grid % Nj_d, grid % Nj_b 

           WRITE(funit,1000)
           WRITE(funit,1025)

           WRITE(funit,1000)
           WRITE(funit,1048)
           WRITE(funit,1050)
           WRITE(funit,1051)

           ! Volume nodes coordinates
           DO j = 1, grid % Nj_d
          
            WRITE(funit,1052) j
            WRITE(funit,1053) grid % rr(:,j)
          
           ENDDO

           WRITE(funit,1000)
           WRITE(funit,1055)

           WRITE(funit,1000)
           WRITE(funit,1048)     
           WRITE(funit,1075)

           ! Boundary nodes - domain nodes connectivity
           DO j = 1, grid % Nj_b
            WRITE(funit,1076) j, grid % jd_jb(j), grid % bound_j(j)
           ENDDO
   
        CLOSE(funit)    

1000    FORMAT('###########################################################################')
1010    FORMAT('#    NAME:      ',a15,'                                           #')
1020    FORMAT('#         ND        NP_D        NP_B                                      #')
1021    FORMAT(5i12)
1025    FORMAT('#  **********  DOMAIN  **********                                         #')
1048    FORMAT('#  ++++   NODES   ++++                                                    #')
1050    FORMAT('#        IDX                                                              #')
1051    FORMAT('#         RR                                                              #')
1052    FORMAT(i12)
1053    FORMAT(3e18.9)
1055    FORMAT('#  **********  BOUNDARY  **********                                       #')
1075    FORMAT('#        IDX       JD_JB       BOUND                                      #')
1076    FORMAT(3i12) 

      END SUBROUTINE  Write_Nodes





      SUBROUTINE  Write_Elements(grid, grid_name)
        !-----------------------------------------------------------------------
        IMPLICIT NONE
        
        TYPE(grid_type),   INTENT(IN) :: grid
        CHARACTER(LEN=64), INTENT(IN) :: grid_name
        
        INTEGER :: j, m, f, funit = 11
        !----------------------------------------------------------------------- 

        OPEN (UNIT=funit, FILE='grid.'//TRIM(grid_name), FORM='formatted', STATUS='unknown')  

          WRITE(funit,2000)
          WRITE(funit,2010) grid_name(1:LEN_TRIM(grid_name))

          WRITE(funit,2000)
          WRITE(funit,2020)
          WRITE(funit,2021) grid % Nm_d, grid % Nm_b

          WRITE(funit,2000)
          WRITE(funit,2025)

          WRITE(funit,2000)
          WRITE(funit,2040)
          WRITE(funit,2041) 
          WRITE(funit,2042) 
          WRITE(funit,2043) 

          ! Element of the domain
          ! ---------------------
          DO m = 1, grid % Nm_d
          
            WRITE(funit,2046) m, grid % ele_type_d(m)
            
            DO j = 1, SIZE( grid % j_m_d(m) % vec ) 
               WRITE(UNIT= funit, FMT= 2047, ADVANCE= 'NO') grid % j_m_d(m) % vec(j)
            ENDDO
            
            WRITE(funit,*)
            
            DO f = 1, SIZE(grid % ma_m_d(m) % vec) 
               WRITE(UNIT= funit, FMT= 2047, ADVANCE= 'NO') grid % ma_m_d(m) % vec(f)
            ENDDO
            
            WRITE(funit,*)
             
          ENDDO


          WRITE(funit,2000)
          WRITE(funit,2055)

          WRITE(funit,2000)
          WRITE(funit,2040)     
          WRITE(funit,2065)
          WRITE(funit,2042) 
          WRITE(funit,2043) 

          ! Element of the boundary
          ! -----------------------
          DO m = 1, grid % Nm_b
          
             WRITE(funit,2066) m, grid % ele_type_b(m), grid % bound_m(m)
          
             DO j = 1, SIZE(grid % j_m_b(m) % vec)
                WRITE(UNIT= funit, FMT= 2047, ADVANCE= 'NO')  grid % j_m_b(m) % vec(j)
             ENDDO
          
             WRITE(funit,*)
          
             DO f = 1, SIZE(grid % ma_m_b(m) % vec)
                WRITE(UNIT= funit, FMT= 2047, ADVANCE= 'NO') grid % ma_m_b(m) % vec(f)
             ENDDO
          
             WRITE(funit,*)
          
          ENDDO
          
        CLOSE (funit)

2000    FORMAT('###########################################################################')
2010    FORMAT('#    NAME:      ',a15,'                                           #')
2020    FORMAT('#       NE_D        NE_B                                                  #')
2021    FORMAT(2i12)
2025    FORMAT('#  **********  DOMAIN  **********                                         #')
2040    FORMAT('#  ++++  ELEMENTS ++++                                                    #')
2041    FORMAT('#        IDX        TYPE                                                  #')
2042    FORMAT('#        J_M                                                              #')
2043    FORMAT('#       MA_M                                                              #')
2046    FORMAT(2i12)
2047    FORMAT(i12)
2055    FORMAT('#  **********  BOUNDARY  **********                                       #')
2065    FORMAT('#        IDX        TYPE       BOUND                                      #')
2066    FORMAT(3i12)
 
      END SUBROUTINE  Write_Elements
!
!
      END MODULE grid_1d
