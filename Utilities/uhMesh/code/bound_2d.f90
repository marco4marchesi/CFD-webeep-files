!     
      module boundary_2d
!
!     +---------------------------------------------------------------------+
      use curves          !  curves module
      use grid_types      !  types module
      use list            !  list module
      use delaunay        !  delaunay basic module
      use eucl_delaunay   !  delaunay euclid module
      use metric_2d       !  metric module
      use schemi
      use qgraph
!
      implicit none
!
      integer :: numbstrucedge
!
      public :: boundary_grid, &
                add_points_1d, &
                add_points_12d, &
                stitch_boundary, &
                halve_bou, &
                numbstrucedge
!     +---------------------------------------------------------------------+
!
!
      contains
!
!
      subroutine boundary_grid ( idf, b_grid, dir, fil )
!     +---------------------------------------------------------------------+
      implicit none
!
      integer,           intent(in)    :: idf
      type (domain),     intent(inout) :: b_grid
      character(len=64), intent(in)    :: dir, fil  
!
      integer :: z, i
      integer :: d_end, f_end, nse
!
      logical, parameter :: ascii = .true.
!     +---------------------------------------------------------------------+
!
      d_end = last_c_leng ( 64, dir )
      f_end = last_c_leng ( 64, fil )
!
      open  ( unit=idf, file='blines.'//fil(1:f_end)//'.mtv' )
      call bou_r_points ( b_grid, idf )
      close (idf)
!
!      if ( b_grid % number_of_wakes .gt. 0 ) then
      if ( b_grid % number_of_Uwakes .gt. 0 ) then
        open  ( unit=idf, file='wakes.'//fil(1:f_end)//'.mtv' )
        call wake_r_points ( b_grid, idf )
        close (idf)
      end if
!
!     +------------------------------------------------+
!     |  Loop to count the number of structured edges  |
!     +------------------------------------------------+
!
      numbstrucedge = 0
      do nse = 1, size( b_grid % edges )
        numbstrucedge = numbstrucedge + &
          b_grid % edges(nse) % struct
      end do
!
!     +---------------------------------------+
!     |  Different operations between HYBRID  |
!     |  and UNSTRUCTURED mesh                |
!     +---------------------------------------+
!  
      if ( numbstrucedge .gt. 0 ) then
        write (*,*)
        write (*,*) '   Generating quadrilaterals layers...'
        call structured_mesh ( b_grid, z )
      else
        allocate ( h_many_edg(b_grid % number_of_edges) )
        h_many_edg = 0
      end if
!
!      open ( unit=idf, file=dir(1:d_end)//'/mesh.'// &
!        fil(1:f_end)//'.mtv',status='REPLACE' )
!      call q_g_plot(idf)
!      close (idf)
!
      write (*,*)
      write (*,*) '   Interface with Delaunay solver...'
      call define_interface
!
!     +----------------------------------------+
!     |  Different operations between HYBRID   |
!     |  and UNSTRUCTURED mesh. Namely in the  |
!     |  case of structured edges the quadrila |
!     |  teral layers are written in the file  |
!     +----------------------------------------+
!  
      if ( numbstrucedge .gt. 0 ) then
!
        open ( unit=idf, file=dir(1:d_end)//'/mesh.'// &
          fil(1:f_end)//'.mtv',status='REPLACE' )
        call q_g_plot(idf)
        close (idf)
!
        open ( unit=idf, file='qt-interface.'//fil(1:f_end)//'.mtv' )
        open ( unit=idf+1, file='np-interface.'//fil(1:f_end) )
        call q_g_plot_jump ( idf, idf+1 )
        close (idf)
        close (idf+1) 
!
      end if
!
      call set_grid ( grid )
      call add_points_12d ( grid )
!
!
      contains
!
!
      subroutine bou_r_points ( b_grid, idf1 )
!
!     +-------------------------------------------------------------------------+
      type (domain), intent(inout) :: b_grid
      integer, intent(in) :: idf1
!
      type (point), pointer :: h1, p1
      real*8, dimension(:), allocatable :: new_s
      real*8, dimension(:), pointer :: s, h
      real*8 :: d1, d2, sdata
!
      integer :: n, ns, j, k
      integer :: z, z1, z2
      integer :: option
!     +-------------------------------------------------------------------------+
!
      nullify ( s, h )
!
      write (idf1,*) '% toplabel='' '' xlabel='' '' ylabel='' ''   '
      write (idf1,*) '#%xmin=   xmax= '
      write (idf1,*) '#%xmin=   xmax= '
      write (idf1,*) '%pointid = true' 
      write (idf1,*) '%equalscale = true'
!
      do i = 1, b_grid % number_of_edges
!
         h1  =>  b_grid % edges(i) % add_head
!
         if ( b_grid % edges(i) % curvetype .eq. 'bms' ) then
!
           ns = b_grid % edges(i) % Npbms ! numero di punti finale
           allocate ( new_s(ns) )
           call points_1d ( h1, s, h, ns )
!
           n = 0
           p1 => h1
           do; p1  =>  p1 % next; if ( associated (p1, h1) ) exit
              n = n+1
              new_s(n) = curv(i) % s(n) / curv(i) % s(ns) ! da 0 a 1, Ascissa curvilinea normalizzata
           enddo
!
         else
!
           call metric_r_1d ( b_grid, i, s, h )
           call points_1d ( h1, s, h )
!
           ns = numberof_point(h1) ! numero di punti finale
           allocate ( new_s(ns) )
!
           n = 0
           p1 => h1
           do; p1  =>  p1 % next; if ( associated (p1, h1) ) exit
              n = n+1
              new_s(n) = p1 % x(1) ! da 0 a 1, Ascissa curvilinea normalizzata
           enddo
!
         end if
!
!
         do j = 1, b_grid % edges(i) % ndata
!
            option = b_grid % edges(i) % option(j)
            sdata = b_grid % edges(i) % sdata(j)
!
            if ( option .ne. 0 ) then
               do z = 1,n
                  z1 = z
                  z2 = z + 1
                  if (sdata == new_s(z1))then
                     !beccato al primo tentativo
                     h1 % next% option = option
                     exit 
                  elseif((sdata < new_s(z2)))then
                     ! a quale z è più vicino?
                     p1 => h1
                     k = 1
                     do
                        p1 => p1 % next
                        !write(*,*) 'abscissa = ',p1%x(1), ' vs ', new_s(z1) 
                        if ( p1 % x(1) > new_s(z1) )  exit
                        k = k + 1
                     enddo
                     d2 = z2 - sdata 
                     d1 =  sdata - z1
                     if (d2 >= d1) then
                        p1 % option = option
                     else
                        p1 % next % option = option
                     endif
                     exit
                  endif
               enddo
            endif
!
         enddo
!
!
         write (idf1,*) '% linecolor = 1    markertype = 12'        
         k = 0
         p1  =>  h1
!
         do; p1  =>  p1 % next; if ( associated (p1, h1) ) exit
!
           k = k + 1
           if ( b_grid % edges(i) % curvetype .eq. 'bms' ) then 
             p1 % x = curv(i) % x(:,k)
           else
             call evalc_0 ( i, p1%x, p1%x, curv )
           end if
!
           p1 % cnd = i
           p1 % old  =>  p1 % prev
           
           write (idf1,'(1x,2(2x,e15.9),2x,i10)') p1 % x, int(new_s(k)*10**3)
!
         enddo
         write (idf1,*)
!
         h1 % next % cnd = -i
         h1 % prev % cnd = -i
!
         nullify ( h1 % next % old )
         write (*,'(4x,a6,i2,1x,a16,1x,i5.5)') '- edge', i, 'Riemann points =', n
!
         if ( b_grid % edges(i) % struct .ge. 1 ) then 

            call p_height( b_grid, i, new_s, h )       
            p1 => h1
            n = 0 

            do; p1  =>  p1 % next; if ( associated (p1, h1) ) exit
               n = n+1
               p1 % height = h(n)         
            enddo

         endif
         deallocate (new_s)
!
      enddo
!
      end subroutine bou_r_points
!
!
!
!
!
    SUBROUTINE wake_r_points ( b_grid, idf1 )

      INTEGER, INTENT(IN) :: idf1
      INTEGER :: n, ns, k
      TYPE(domain), INTENT(INOUT) :: b_grid
      TYPE (point) , POINTER :: h1, p1
      REAL*8, DIMENSION(:), POINTER :: s, h
      REAL*8, DIMENSION(:), ALLOCATABLE :: new_s
      !REAL(KIND=8) :: d1, d2!, sdata

      NULLIFY ( s, h)

      !  1d points
      !write(*,*)'Warning: n majored to obtain l<1'
      ! write(*,*)'options assigned to the most near point!'

      WRITE (idf1,*) '% toplabel='' '' xlabel='' '' ylabel='' ''   '
      WRITE (idf1,*) '#%xmin=   xmax= '
      WRITE (idf1,*) '#%xmin=   xmax= '
      WRITE (idf1,*) '%pointID = TRUE' 

      DO i=1,b_grid % number_of_Swakes
         h1  =>  b_grid % wakes(i) % add_head

         !write(*,*) 'before wake_metric_r_1d'
         !final number of points determining(metric imposal)
         CALL wake_metric_r_1d (b_grid,i, s, h)
         !write(*,*) 'before points_1d'
         CALL points_1d (h1, s, h)

         ns = numberof_point(h1) ! numero di punti finale

         ALLOCATE (new_s(ns))

         n = 0
         p1 => h1
         DO; p1  =>  p1 % next; IF ( ASSOCIATED (p1, h1) ) EXIT
            n = n+1
            new_s(n) = p1 % x(1) ! da 0 a 1 ,ascissa curvilinea normalizzata
         ENDDO
         !--------------------------------------------------------        
         !         !write(*,*)'options assigned to the most near point!'
         !         DO j =1,b_grid % wakes(i) % ndata
         !          option = b_grid % wakes(i) % option(j)
         !          sdata = b_grid % wakes(i) % sdata(j)
         !          IF (option /= 0) THEN
         !           
         !           DO z = 1,n
         !             z1 = z
         !             z2 = z + 1
         !             IF (sdata == new_s(z1))THEN
         !               !beccato al primo tentativo
         !                h1 % next% option = option
         !		EXIT 
         !             ELSEIF((sdata < new_s(z2)))THEN
         !               ! a quale z è più vicino?
         !                p1 => h1
         !                k = 1
         !                DO
         !                 p1 => p1 % next
         !                 
         !                 IF(p1 % x(1) > new_s(z1))THEN
         !                  EXIT
         !                 ENDIF
         !                 k = k + 1
         !                ENDDO
         !                d2 = z2 - sdata 
         !                d1 =  sdata - z1
         !                IF (d2 >= d1) THEN
         !                  p1 % option = option
         !                ELSE
         !                  p1 % next % option = option
         !                ENDIF
         !
         !                EXIT
         !             ENDIF
         !           ENDDO
         !          ENDIF
         !         ENDDO
         !         !WRITE(*,*)'options done'
         !---------------------------------
         WRITE (idf1,*) '% linecolor = 1    markertype = 12'  
         !WRITE (idf,*) '% pointID = TRUE'  
         k = 0
         p1  =>  h1

         DO; p1  =>  p1 % next; IF ( ASSOCIATED (p1, h1) ) EXIT

            CALL evalc_0(i,p1%x,p1%x,wcurv)
            k = k + 1
            p1 % cnd = i
            p1 % old  =>  p1 % prev
            WRITE(idf1,*) p1%x, INT(new_s(k)*10**3)
         ENDDO
         WRITE(idf1,*) !spazio vuoto


         h1 % next % cnd = -i
         h1 % prev % cnd = -i
         NULLIFY (h1 % next % old)

         WRITE (*,'(4x,a6,i2,1x,a16,1x,i5.5)') '- wake', i, 'Riemann points =', n
         !IF(b_grid % wakes(i) % struct >= 1)THEN 
         CALL wake_p_height(b_grid,i,new_s,h)       
         p1 => h1
         n = 0 

         DO; p1  =>  p1 % next; IF ( ASSOCIATED (p1, h1) ) EXIT
            n = n+1
            p1 % height = h(n)         
         ENDDO
         !ENDIF 
         DEALLOCATE (new_s)
      ENDDO


    END SUBROUTINE wake_r_points
    !-------------------------------------------------------------------
    SUBROUTINE p_height (b_grid,c, new_s, h)

      TYPE(domain),INTENT(INOUT)::b_grid
      INTEGER, INTENT(IN) :: c
      REAL*8, DIMENSION(:), POINTER ::  h

      INTEGER :: ns, ms, i, j, j0, jm
      REAL*8 :: ss, s0, h0, si, hi, fs, fh, sh, dh, tj
      INTEGER,       DIMENSION(:)    , POINTER :: idata
      REAL*8, DIMENSION(:)    , POINTER :: sdata
      REAL*8, DIMENSION(:)    , POINTER :: height
      REAL*8, DIMENSION(:)     :: new_s

      ns = numberof_point(b_grid % edges(c)% add_head)

      ! h vettore delle altezze dei punti finale di contorno
      IF ( ASSOCIATED ( h ) ) DEALLOCATE ( h )
      ALLOCATE ( h(ns) )

      idata  =>  b_grid % edges(c) % idata
      sdata  =>  b_grid % edges(c) % sdata
      height =>  b_grid % edges(c) % height

      ms = SIZE(idata)        ! punti di imposizione metrica

      j0 = 1
      ss = new_s(ns)          ! ascissa ultimo punto(normalizzata),1
      s0 = ss*sdata(1)        ! 0

      h0 = height(1)
      h(j0) = h0

      DO i=2,ms

         IF(height(i) <= 0) CYCLE
         ! per by-passare i punti nei quali non voglio imporre l'altezza

         jm = j0
         si = ss*sdata(i)     ! ascissa normalizzata dei pti imposizione metrica
         hi = height(i)

         DO j=j0+1,ns         ! ns numero di punti del lato
            IF ( new_s(j) >= si ) THEN
               jm = j-1       ! mi fermo prima del punto di imposizione
               EXIT           ! il nuovo jm
            ENDIF
         ENDDO

         IF ( jm == j0 ) THEN ! lato costituito di due soli punti????
            WRITE (*,*) 'build_h error for si = ',si
            WRITE (*,*) '--STOP-- (p_height)'
            STOP
         ENDIF

         SELECT CASE (idata(i))

         CASE (1) ! Linear interpolation

            fs = 1.0/(si-s0)
            dh = hi-h0

            DO j=j0+1,jm      ! mi fermo prima del punto di imposizione
               tj   = fs*(new_s(j)-s0)
               h(j) = h0+tj*dh
            ENDDO

         CASE (2) ! Geometrical interpolation

            fs = 1.0/(si-s0)
            fh = hi/h0
            DO j=j0+1,jm
               tj   = fs*(new_s(j)-s0)
               h(j) = h0*(fh**tj)
            ENDDO

         CASE (3) ! Sinusoidal interpolation

            fs = pi/(si-s0)
            sh = 0.5*(hi+h0); dh = 0.5*(hi-h0)
            DO j=j0+1,jm
               tj   = fs*(new_s(j)-s0)
               h(j) = sh-dh*cos(tj)
            ENDDO

         END SELECT

         j0 = jm; s0 = si; h0 = hi 

      ENDDO
      h(ns) = height(ms)  ! semplice trasferimento dell'informazione disponibile

    END SUBROUTINE p_height
    !---------------------------------------------------------------------------

    SUBROUTINE metric_r_1d ( b_grid, c, s, h )

      TYPE(domain),INTENT(INOUT)::b_grid
      INTEGER, INTENT(IN) :: c
      REAL*8, DIMENSION(:), POINTER :: s, h

      INTEGER :: ns, ms, i, j, j0, jm, edge_kind
      REAL*8 :: ss, s0, h0, si, hi, fs, fh, sh, dh, tj
      INTEGER,       DIMENSION(:)    , POINTER :: idata
      REAL*8, DIMENSION(:)    , POINTER :: sdata
      REAL*8, DIMENSION(:,:)  , POINTER :: hdata


      CALL evalc_s ( c, ns, s, curv )

      IF ( ASSOCIATED ( h ) ) DEALLOCATE ( h )
      ALLOCATE ( h(ns) )

      edge_kind = b_grid % edges(c) % edge_kind

      idata  =>  b_grid % edges(c) % idata
      sdata  =>  b_grid % edges(c) % sdata
      hdata  =>  b_grid % edges(c) % hdata

      ms = SIZE(idata)

      !  Compute h(s) --- interpolated local spacing at curve knots

      j0 = 1
      ss = s(ns)
      s0 = ss*sdata(1)

      !  First point projected metric

      CALL proj_metr(c, edge_kind, sdata(1), hdata(:,1), h0)

      h(j0) = h0

      DO i=2,ms

         jm = j0
         si = ss*sdata(i)
         CALL proj_metr(c, edge_kind, sdata(i), hdata(:,i), hi)

         DO j=j0+1,ns
            IF ( s(j) >= si ) THEN
               jm = j-1
               EXIT
            ENDIF
         ENDDO

         IF ( jm == j0 ) THEN
            WRITE (*,*) 'build_h error for si = ',si
            WRITE (*,*) '--STOP-- (metric_r_1d)'
            STOP
         ENDIF

         SELECT CASE (idata(i))

         CASE (1) ! Linear interpolation

            fs = 1.0/(si-s0)
            dh = hi-h0

            DO j=j0+1,jm
               tj   = fs*(s(j)-s0)
               h(j) = h0+tj*dh
            ENDDO

         CASE (2) ! Geometrical interpolation

            fs = 1.0/(si-s0)
            fh = hi/h0

            DO j=j0+1,jm
               tj   = fs*(s(j)-s0)
               h(j) = h0*(fh**tj)
            ENDDO

         CASE (3) ! Sinusoidal interpolation

            fs = pi/(si-s0)
            sh = 0.5*(hi+h0); dh = 0.5*(hi-h0)

            DO j=j0+1,jm
               tj   = fs*(s(j)-s0)
               h(j) = sh-dh*cos(tj)
            ENDDO

         END SELECT

         j0 = jm; s0 = si; h0 = hi 

      ENDDO

      CALL proj_metr(c, edge_kind, sdata(ms), hdata(:,ms), h(ns))

      !  Integral of 1/h(s) --- riemannian curvilinear abscissa at knots

      hi   = 1.0/h(1)
      h(1) = 0.0
      DO i = 2, ns
         h0   = hi
         hi   = 1.0/h(i)
         dh   = 0.5*(s(i)-s(i-1))*(hi+h0)
         !**         if(dh > .5d0) then
         !**            write(*,*)'need a point! dh > 1/2 (metric_r_1d)'
         !**         endif
         h(i) = h(i-1) + dh
      ENDDO

    END SUBROUTINE metric_r_1d
    !
    !********************************************************************
    !
    SUBROUTINE wake_p_height (b_grid,c, new_s, h)

      TYPE(domain),INTENT(INOUT)::b_grid
      INTEGER, INTENT(IN) :: c
      REAL*8, DIMENSION(:), POINTER ::  h

      INTEGER :: ns, ms, i, j, j0, jm
      REAL*8 :: ss, s0, h0, si, hi, fs, fh, sh, dh, tj
      INTEGER,       DIMENSION(:)    , POINTER :: idata
      REAL*8, DIMENSION(:)    , POINTER :: sdata
      REAL*8, DIMENSION(:)    , POINTER :: height
      REAL*8, DIMENSION(:)     :: new_s

      ns = numberof_point(b_grid % wakes(c)% add_head)

      ! h vettore delle altezze dei punti finale di contorno
      IF ( ASSOCIATED ( h ) ) DEALLOCATE ( h )
      ALLOCATE ( h(ns) )

      idata  =>  b_grid % wakes(c) % idata
      sdata  =>  b_grid % wakes(c) % sdata
      height =>  b_grid % wakes(c) % height

      ms = SIZE(idata)        ! punti di imposizione metrica

      j0 = 1
      ss = new_s(ns)          ! ascissa ultimo punto(normalizzata),1
      s0 = ss*sdata(1)        ! 0

      h0 = height(1)
      h(j0) = h0

      DO i=2,ms

         IF(height(i) <= 0) CYCLE
         ! per by-passare i punti nei quali non voglio imporre l'altezza

         jm = j0
         si = ss*sdata(i)     ! ascissa normalizzata dei pti imposizione metrica
         hi = height(i)

         DO j=j0+1,ns         ! ns numero di punti del lato
            IF ( new_s(j) >= si ) THEN
               jm = j-1       ! mi fermo prima del punto di imposizione
               EXIT           ! il nuovo jm
            ENDIF
         ENDDO

         IF ( jm == j0 ) THEN ! lato costituito di due soli punti????
            WRITE (*,*) 'build_h error for si = ',si
            WRITE (*,*) '--STOP-- (p_height)'
            STOP
         ENDIF

         SELECT CASE (idata(i))

         CASE (1) ! Linear interpolation

            fs = 1.0/(si-s0)
            dh = hi-h0

            DO j=j0+1,jm      ! mi fermo prima del punto di imposizione
               tj   = fs*(new_s(j)-s0)
               h(j) = h0+tj*dh
            ENDDO

         CASE (2) ! Geometrical interpolation

            fs = 1.0/(si-s0)
            fh = hi/h0
            DO j=j0+1,jm
               tj   = fs*(new_s(j)-s0)
               h(j) = h0*(fh**tj)
            ENDDO

         CASE (3) ! Sinusoidal interpolation

            fs = pi/(si-s0)
            sh = 0.5*(hi+h0); dh = 0.5*(hi-h0)
            DO j=j0+1,jm
               tj   = fs*(new_s(j)-s0)
               h(j) = sh-dh*cos(tj)
            ENDDO

         END SELECT

         j0 = jm; s0 = si; h0 = hi 

      ENDDO
      h(ns) = height(ms)  ! semplice trasferimento dell'informazione disponibile

    END SUBROUTINE wake_p_height
    !
    !********************************************************************
    !
    SUBROUTINE wake_metric_r_1d (b_grid,c, s, h)

      TYPE(domain),INTENT(INOUT)::b_grid
      INTEGER, INTENT(IN) :: c
      REAL*8, DIMENSION(:), POINTER :: s, h

      INTEGER :: ns, ms, i, j, j0, jm, edge_kind
      REAL*8 :: ss, s0, h0, si, hi, fs, fh, sh, dh, tj
      INTEGER,       DIMENSION(:)    , POINTER :: idata
      REAL*8, DIMENSION(:)    , POINTER :: sdata
      REAL*8, DIMENSION(:,:)  , POINTER :: hdata

      CALL evalc_s (c, ns, s,wcurv)

      IF ( ASSOCIATED ( h ) ) DEALLOCATE ( h )
      ALLOCATE ( h(ns) )

      edge_kind = b_grid % wakes(c) % edge_kind

      sdata  =>  b_grid % wakes(c) % sdata
      idata  =>  b_grid % wakes(c) % idata
      hdata  =>  b_grid % wakes(c) % hdata

      ms = SIZE(idata)

      !  Compute h(s) --- interpolated local spacing at curve knots

      j0 = 1
      ss = s(ns)
      s0 = ss*sdata(1)

      !  First point projected metric

      CALL proj_metr(c, edge_kind, sdata(1), hdata(:,1), h0)

      h(j0) = h0

      DO i=2,ms

         jm = j0
         si = ss*sdata(i)
         CALL proj_metr(c, edge_kind, sdata(i), hdata(:,i), hi)

         DO j=j0+1,ns
            IF ( s(j) >= si ) THEN
               jm = j-1
               EXIT
            ENDIF
         ENDDO

         IF ( jm == j0 ) THEN
            WRITE (*,*) 'build_h error for si = ',si
            WRITE (*,*) '--STOP-- (metric_r_1d)'
            STOP
         ENDIF

         SELECT CASE (idata(i))

         CASE (1) ! Linear interpolation

            fs = 1.0/(si-s0)
            dh = hi-h0

            DO j=j0+1,jm
               tj   = fs*(s(j)-s0)
               h(j) = h0+tj*dh
            ENDDO

         CASE (2) ! Geometrical interpolation

            fs = 1.0/(si-s0)
            fh = hi/h0

            DO j=j0+1,jm
               tj   = fs*(s(j)-s0)
               h(j) = h0*(fh**tj)
            ENDDO

         CASE (3) ! Sinusoidal interpolation

            fs = pi/(si-s0)
            sh = 0.5*(hi+h0); dh = 0.5*(hi-h0)

            DO j=j0+1,jm
               tj   = fs*(s(j)-s0)
               h(j) = sh-dh*cos(tj)
            ENDDO

         END SELECT

         j0 = jm; s0 = si; h0 = hi 

      ENDDO

      CALL proj_metr(c, edge_kind, sdata(ms), hdata(:,ms), h(ns))

      !  Integral of 1/h(s) --- riemannian curvilinear abscissa at knots

      hi   = 1.0/h(1)
      h(1) = 0.0
      DO i=2,ns
         h0   = hi
         hi   = 1.0/h(i)
         dh   = 0.5*(s(i)-s(i-1))*(hi+h0)
         !**         if(dh > .5d0) then
         !**            write(*,*)'need a point! dh > 1/2 (metric_r_1d)'
         !**         endif
         h(i) = h(i-1) + dh
      ENDDO

    END SUBROUTINE wake_metric_r_1d
    !
    !********************************************************************
    !
    SUBROUTINE points_1d ( head, s, h, Nbms )

      TYPE(point)               , POINTER :: head
      REAL(KIND=8), DIMENSION(:), POINTER :: s, h
      INTEGER, OPTIONAL :: nBMS


      INTEGER :: ns, n, j0, i, j
      REAL*8 :: hh, ss, dh, hi, si
      TYPE (point), POINTER :: p


      !  Subdivide the curve into L segments of unit riemannian
      !  length, L being the riemannian length of the curve.

!      ns = SIZE(s)      
!      hh = h(ns)

      if ( present(Nbms) ) then

        n = Nbms - 1
        ns = Nbms
        hh = 1.0        

        if ( associated(s) ) deallocate (s)
        allocate( s(ns) )
        s = 1.d0
        
        if ( associated(h) ) deallocate (h)
        allocate( h(ns) )
        h = 1.d0

      else

        ns = SIZE(s)      
        hh = h(ns)
        n  = MAX(INT(hh),1)

      end if

      dh = hh/float(n)
      if ( dh > 1.d0 ) then
         n = n+1
         dh = hh / float(n)
      endif
      ss = 1.0/s(ns)
      
      j0 = 1
      p  =>  new_point ()
      p % x(1) = 0.0

      CALL pushbot_point (head, p)

      DO i = 2, n
         hi = dh*(i-1)
         DO j = j0+1, ns
            IF ( h(j) >= hi ) THEN
               j0 = j-1
               si = s(j-1) + (hi-h(j-1))*(s(j)-s(j-1))/(h(j)-h(j-1))

               p  =>  new_point ()
               p % x(1) = si*ss   !'Riemann abscissa'

               CALL pushbot_point (head, p)
               EXIT
            ENDIF
         ENDDO
      ENDDO

      p  =>  new_point ()
      p % x(1) = 1.0

      CALL pushbot_point (head, p)

    END SUBROUTINE points_1d
    !
    !********************************************************************
    !
    SUBROUTINE set_grid (grid)

      INTEGER :: i
      REAL*8, DIMENSION(:), ALLOCATABLE :: xlb, xrt
      TYPE (point), POINTER :: hap, p
      TYPE (simplex), POINTER :: has, hs
      TYPE(domain),INTENT(INOUT)::grid
      !  1d grids

      ALLOCATE ( xlb(nd_b) )
      ALLOCATE ( xrt(nd_b) )

      xlb = 0.
      xrt = 1.

      DO i=1,grid % number_of_edges

         ! 1d box definition
         hap  =>  grid % edges(i) % axp_head  ! auxiliary point head
         has  =>  grid % edges(i) % axs_head  ! auxiliary simplex head
         hs   =>  grid % edges(i) % spx_head  ! simplex head

         CALL set_delaunay (nd_b, xlb, xrt, hap, has, hs)

      ENDDO

      DEALLOCATE ( xlb )
      DEALLOCATE ( xrt )

      !  2d grid

      ALLOCATE ( xlb(nd_d) )
      ALLOCATE ( xrt(nd_d) )

      xlb = +HUGE(xlb)
      xrt = -HUGE(xrt)

      !  find max and min  coordinates values
      !  (box must include every point)
      p  =>  grid % add_head      
      DO; p  =>  p % next; IF ( ASSOCIATED (p, grid % add_head) ) EXIT
         xlb = MIN(xlb, p % x); xrt = MAX(xrt, p % x)
      END DO
      ! 2d box definition
      hap  =>  grid % axp_head
      has  =>  grid % axs_head
      hs   =>  grid % spx_head

      CALL set_delaunay (nd_d, xlb, xrt, hap, has, hs)

      DEALLOCATE ( xlb )
      DEALLOCATE ( xrt )

    END SUBROUTINE set_grid
!
      end subroutine boundary_grid
!
!
!
!
!
      subroutine add_points_1d ( grid )
!
!     +----------------------------------------------------+
      type (domain), intent(inout) :: grid
!
      type (simplex), pointer :: hs, hi, he
      type (point), pointer :: p, q
!
      integer :: i
!     +----------------------------------------------------+
!
      do i = 1, grid % number_of_edges
!
         hs => grid % edges(i) % spx_head
         hi => grid % edges(i) % int_head
         he => grid % edges(i) % ext_head
!
         call merge_ie ( hs, hi, he )
         ! 
         ! int_head and ext_head pushtop in spx_head
         ! put status=und to every simplex
!
         p => grid % edges(i) % add_head % next
         do
!
            if ( associated (p, grid % edges(i) % add_head) ) exit
            q  =>  p % next
!
            call insert_point ( nd_b, hs, p )
            call movetop_point ( grid % edges(i) % pnt_head, p )
            !inserted points list
            p => q
!
         end do
!
!     +---------------------------------------------------+
!     |  Moves in int the simplices having status /= frz  |
!     |  ext the simplices having status = frz            |
!     +---------------------------------------------------+
!
         call split_ie_1d ( nd_b, hs, hi, he )
!
      end do

      end subroutine add_points_1d
!
!
!
!
!
      subroutine add_points_12d ( grid )
!
!     +----------------------------------------------------+
      implicit none
!
      type(domain), intent(inout) :: grid
!
      type(simplex), pointer :: hs, hi, he
      type(simplex), pointer :: hbi, hbs
      type(point), pointer :: p, q
!
      integer :: i, nstc, k
!     +----------------------------------------------------+
!
      call add_points_1d ( grid )
!
      hs => grid % spx_head
      hi => grid % int_head
      he => grid % ext_head

      k = 0

      add_2d_points: do

         k = k + 1

         call merge_ie ( hs, hi, he )

         p => grid % add_head % next

         do; if ( associated (p, grid % add_head) ) exit
            q  =>  p % next

            call insert_point  (nd_d, grid % spx_head, p)
            call movetop_point (grid % pnt_head, p)
            !inserted points list  
            p  =>  q
         end do


         do i = 1, grid % number_of_edges
           hbi  =>  grid % edges(i) % int_head
           hbs  =>  grid % edges(i) % stc_head
           call flag_shell ( nd_d, hbi, hs, hbs )
         end do


         call stitch_boundary ( grid, nstc )

         if ( nstc .gt. 0 ) then

            write (*,*) 'stitch points = ',nstc          
            call add_points_1d ( grid )

         else
!
            ! gets from spx_h internal and external boundary adjacent simplices
            ! then spreads the internal and external status
            ! except for frozen status simplices(box adjacents)
            call split_ie (nd_d, hs, hi, he)
            exit

         end if

      end do add_2d_points

      end subroutine add_points_12d
!
!
!
!
!
      subroutine stitch_boundary ( grid, nadd )
!
!     +---------------------------------------------+
      type (domain), intent(inout) :: grid
      integer, intent(out) :: nadd
!
      type (simplex), pointer :: head, s
      integer :: e, g, n
!     +---------------------------------------------+
!
      n = 0
      do e = 1, grid % number_of_edges
!
         g = grid % edges(e) % index
!
         head => grid % edges(e) % stc_head
         s => head
!
         do
!
           s  =>  s % next
           if ( associated (s, head) ) exit
!
           n = n + 1
           call halve_bou ( grid, e, g, s )
!
         end do
!
      end do
      nadd = n
!
      end subroutine stitch_boundary
!
!
!
!
!
  SUBROUTINE halve_bou (grid, ibou, igeo, sbou)

    TYPE (domain), INTENT(INOUT) :: grid
    TYPE (simplex), TARGET, INTENT(IN) :: sbou
    INTEGER :: ibou, igeo

    INTEGER :: i
    REAL*8 :: factor
    TYPE (point), POINTER :: p_dom, p_bou

    factor = 1./DBLE(nd_b+1)

    p_bou  =>  new_point()
    p_dom  =>  new_point()

    p_bou % cnd = ibou

    IF(ASSOCIATED(sbou % int)) THEN		! by trive
       p_bou % parent = sbou % int % index
    ELSE
       p_bou % parent = 0
    ENDIF

    p_bou % old  =>  sbou % opp(1) % pnt
    p_bou % ppp  =>  p_dom
    p_bou % x = 0.
    DO i=1,nd_b+1
       p_bou % x = p_bou % x + sbou % opp(i) % pnt % x
    ENDDO

    p_bou % x = factor * p_bou % x

    p_dom % cnd = ibou

    IF(ASSOCIATED(sbou % int)) THEN		! by trive
       p_dom % parent = sbou % int % index
    ELSE
       p_dom % parent = 0
    ENDIF

    p_dom % old  =>  sbou % opp(1) % pnt % ppp

    p_dom % metric = 0.


    DO i=1,nd_b+1
       !write(*,*)'sbou % opp(',i,') % pnt % ppp % x',sbou % opp(i) % pnt % ppp % x

       p_dom % metric = p_dom % metric + sbou % opp(i) % pnt % ppp % metric
    ENDDO

    p_dom % metric = factor * p_dom % metric

    CALL evalc_0       (igeo, p_bou % x, p_dom % x,curv)

    CALL pushtop_point (grid % edges(ibou) % add_head, p_bou)

    CALL pushtop_point (grid % add_head, p_dom)

  END SUBROUTINE halve_bou

  !-------------------------------------------------------------------------

  SUBROUTINE define_interface

    IMPLICIT NONE

    INTEGER :: i, j, k, n, nedg, z, bv, t
    TYPE(point), POINTER :: b_head, p, q, head, p1, h1, p2 
    REAL(KIND=8) :: s, ds
    TYPE(edge), POINTER :: edg

    !modificare le liste di punti relative ai lati da triangolare 
    !adiacenti ai lati quadrangolati(restrizioni)
    !###=> usare solo le informazioni dei punti,PER TUTTI!
    !### interfaccia solo per ADD_POINTS_12D

    !++++++++++++++++++++++++++++++++++DETERMINAZIONE DEL NUMERO DI LATI FINALE
    nedg = SIZE( h_many_edg )
    n = 0

    DO i = 1, b_grid % number_of_edges
      IF ( b_grid % edges(i) % status .eq. 0 ) n = n + 1
    END DO

    DO i = 1, nedg
      IF ( h_many_edg(i) .ne. 0 )  n = n + 1
    END DO

    grid % number_of_edges = n
    grid % number_of_verts = n
    
    ALLOCATE ( grid%edges(n) )
    ALLOCATE ( grid%verts(n) )

    DO i = 1, n

       grid % edges(i) % index = i

       grid % edges(i) % pnt_head  =>  newhead_point ()
       grid % edges(i) % axp_head  =>  newhead_point ()
       grid % edges(i) % add_head  =>  newhead_point ()
       grid % edges(i) % rem_head  =>  newhead_point ()
       grid % edges(i) % gar_head  =>  newhead_point ()

       grid % edges(i) % spx_head  =>  newhead_simplex ()
       grid % edges(i) % axs_head  =>  newhead_simplex ()
       grid % edges(i) % int_head  =>  newhead_simplex ()
       grid % edges(i) % ext_head  =>  newhead_simplex ()
       grid % edges(i) % stc_head  =>  newhead_simplex ()

    ENDDO

    grid % pnt_head  =>  newhead_point ()
    grid % axp_head  =>  newhead_point ()
    grid % add_head  =>  newhead_point ()
    grid % rem_head  =>  newhead_point ()
    grid % gar_head  =>  newhead_point ()

    grid % spx_head  =>  newhead_simplex ()
    grid % axs_head  =>  newhead_simplex ()
    grid % int_head  =>  newhead_simplex ()
    grid % ext_head  =>  newhead_simplex ()
    grid % stc_head  =>  newhead_simplex ()

    k = 0
    z = 0
    i = 1
    n = 1

    DO WHILE ( i .le. b_grid % number_of_edges )

       edg => b_grid % edges(i)
       k = k+1

       IF ( n .eq. 1 ) THEN
         grid % edges(k) % begin_vert = k
         bv = k
       END IF
       
       IF ( associated(edg % next) ) THEN
       
          IF ( edg % struct .ge. 1 ) THEN
             n = 1
             z = z+1
             j = h_many_edg(z)-1

             DO t = 1, j
                n = n+1 
                edg => edg % next
             ENDDO
             
             IF ( associated(edg % next) )THEN

                ! si continua il lato ,ma in modo non strutturato
                grid % edges(k) % end_vert = k+1
                grid % edges(k+1) % begin_vert = k+1
                i = i + n 
                n = n+1      


             ELSE

                grid % edges(k) % end_vert = k
                !il lato è self-closed
                i = i + n 
                !write(*,*) 'i next',i
                n = 1
             ENDIF
             
          ELSE

             grid % edges(k) % end_vert = k+1
             grid % edges(k+1) % begin_vert = k+1
             n = n+1
             i = i+1
             
          ENDIF
       
       ELSE

          IF ( edg % struct .ge. 1 ) z = z+1
          grid % edges(k) % end_vert = bv
          i = i+1 
          n = 1
       
       ENDIF
    ENDDO

    DO i = 1, grid % number_of_edges
       !WRITE(*,*)'  edge',i,' beg vert',grid %edges(i)%begin_vert,' end vert',grid %edges(i)%end_vert
       grid % edges(i) % struct = 0
    ENDDO

    !DO i=1,b_grid%number_of_edges
    !  WRITE(*,*)'edg',i,'status',b_grid%edges(i)%status
    !ENDDO 

    z = 0
    j = 0 
    i = 1

    if ( (b_grid % edges(1) % status .ne. 0  .or.  &
         b_grid % edges(b_grid % number_of_edges) % status .ne. 0 ) .and. &
	 any(b_grid % edges % status .eq. 0) ) then
      write (*,*) ''
      write (*,*) 'WARNING: Your mesh might have an open strucutred edge'
      write (*,*) 'Please ensure that such edge is not the first or last'
      write (*,*) 'edge in the series of edges since this might cause problems'
      write (*,*) 'A structured open edge should not appear as first or last'
      write (*,*) 'in the list of edges.'
      write (*,*) ''
    end if

    DO 
       IF ( b_grid % edges(i) % status .eq. 0 ) THEN
       
          j=j+1
          WRITE(*,'(4x,a6,1x,i2,1x,a14)') '- edge',j,'type: Delaunay'
          grid % edges(j) % struct = 0

          b_grid % edges(i) % ppp => grid % edges(j)  
          grid % edges(j) % old   => b_grid % edges(i)

          head => grid%edges(j)%add_head
          b_head => b_grid%edges(i)%add_head
          p => b_head
          DO;p => p% next;IF(ASSOCIATED(p,b_head))EXIT
             q => new_point()
             q % x = p % x
             p % ppp => q 
             q % old1 => p !NOVITA' 10-02-2004
             q % option = p % index

             q % cnd = j
             CALL pushbot_point(head,q)
          ENDDO
          head % next % cnd = -j
          head % prev % cnd = -j
          NULLIFY(head % next % old)
          i = i+1
          
       ELSE
          
          j = j+1
          WRITE (*,'(4x,a6,1x,i2,1x,a11)') '- edge',j,'type: Layer'
          z = z+1

          grid % edges(j) % struct = 1
          head => grid % edges(j) % add_head
          
          !    b_head => layer(done_layers)%pnt_list(z)
          IF ( SIZE(big_box) > 0 ) THEN
             b_head => big_box(z) % virt_new_head
             IF(b_head % opnd == 1)THEN
                b_head => layer(1) % pnt_list(z)      
             ENDIF
          ELSE
             b_head => layer(done_layers)%pnt_list(z)
          ENDIF

          p => b_head

          DO 
             q => new_point()
             q % x = p % x
             p % ppp => q 
             q % old1 => p 
             q % option = p % index
             q % cnd = j
             CALL pushbot_point(head,q)

             p => p % jump_f
             IF(ASSOCIATED(p,b_head))THEN

                IF(b_head% opnd == 0) THEN !ripetizione del primo elemento

                   q => new_point()
                   q % x = b_head % x
                   q % old1 => head % next % old1
                   q % option = head % index
                   q % cnd = j
                   CALL pushbot_point(head,q)
                ENDIF

                EXIT

             ENDIF
          ENDDO

          head % next % cnd = -j
          head % prev % cnd = -j
          NULLIFY(head % next % old)
          i = i + h_many_edg(z)
          
       ENDIF
       
       IF ( i .gt. b_grid % number_of_edges ) EXIT

    ENDDO

    CALL connect_edges2dir( grid )
    CALL adjust_bound

    !write(*,*)'points 2d'!---------------------- variante originale

    DO i =1,grid %number_of_verts !vertex points
       p2=>new_point()
       grid % verts(i) % pnt => p2
       CALL pushtop_point(grid % add_head,p2)
    ENDDO

    DO i = 1, grid % number_of_edges !edge points

       h1 => grid %edges(i)%add_head
       p1 => h1 % next
       p2 => grid % verts(grid %edges(i)%begin_vert)%pnt
       p1 % ppp =>p2

       p2 % cnd = p1 % cnd
       p2 % x = p1 % x

       h1 => h1 % prev

       DO;p1 => p1 % next;IF(ASSOCIATED(p1,h1))EXIT
          p2 => new_point()
          p1 % ppp => p2 
          p2 % cnd = p1 % cnd
          p2 % x = p1 % x

          CALL pushtop_point(grid%add_head,p2)
       ENDDO


       p2 => grid % verts(grid % edges(i)%end_vert)% pnt
       p1 % ppp => p2 

       p2 % cnd = p1 % cnd

       IF ( grid % edges(i) % begin_vert .ne. grid % edges(i) % end_vert ) THEN
          p2 % x = p1 % x
       END IF

    ENDDO

    !----------------------------------------------------------------------
    !normalizzazione delle coordinate dei punti di grid%edges(i)% add_head
    !write(*,*)'normalizzare' 

    DO i = 1, grid % number_of_edges
    
       n = numberof_point(grid % edges(i) % add_head)

       ds = 1 / float(n-1)

       s = 0
       head => grid %edges(i)%add_head
       p1 => head % next

       p1 % x(1) = s

       DO; p1 => p1 % next; IF(ASSOCIATED(p1,head))EXIT
          s = s + ds

          p1 % x(1) = s

       ENDDO
       p1 % x(1) = 1

    ENDDO

  END SUBROUTINE define_interface
  !----------------------------------------------------------------------




  SUBROUTINE n_circ (head,s)
    IMPLICIT NONE
    TYPE(point),POINTER::p,head
    INTEGER ,INTENT(OUT):: s
    p => head
    s = 1
    DO;p => p % next; IF (ASSOCIATED(p,head))EXIT
       s = s + 1
    ENDDO

  END SUBROUTINE n_circ





  SUBROUTINE connect_edges2dir(grid)

    IMPLICIT NONE

    TYPE( domain ), INTENT(INOUT) :: grid
    INTEGER                       :: i, j, bv, ev, bv1, ev1
    
    !#### SUPPONENDO UN VALORE EFFETTIVO DI ORDINE PER index
    DO i = 1, grid % number_of_edges
      NULLIFY ( grid % edges(i) % next )
      NULLIFY ( grid % edges(i) % prev )
      grid % edges(i) % status = 0
    ENDDO
    
    DO i = 1, grid % number_of_edges

       bv = grid % edges(i) % begin_vert
       ev = grid % edges(i) % end_vert
       grid % edges(i) % status = 1

       DO j = 1, grid % number_of_edges       
       
          IF ( grid % edges(j) % status .eq. 1 ) CYCLE

          bv1 = grid % edges(j) % begin_vert
          ev1 = grid % edges(j) % end_vert

          IF ( bv1 .eq. ev ) THEN ! il lato continua
             grid % edges(j) % status = 1
             grid % edges(i) % next => grid % edges(j)
             grid % edges(j) % prev => grid % edges(i)    
          END IF
          
          IF ( ev1 .eq. bv ) THEN ! il lato chiude
             grid % edges(i) % prev => grid % edges(j)
             grid % edges(j) % next => grid % edges(i)! ** diff con init_data
             !il primo elemento del blocco punta anch all'ultimo(circolarità minima)
          END IF
          
       ENDDO

    ENDDO
    
    DO i = 1, grid % number_of_edges
       grid % edges(i) % status = 0
    ENDDO
    !un lato self-closing non ha riferimenti
    !un lato collegato ha il riferimento 
  END SUBROUTINE connect_edges2dir





  FUNCTION last_c_leng (len_str, string) RESULT (leng)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: len_str
    CHARACTER (LEN=len_str), INTENT(IN) :: string
    INTEGER :: leng

    INTEGER :: i

    leng = len_str

    DO i=1,len_str
       IF ( string(i:i) .EQ. ' ' ) THEN
          leng = i-1; EXIT
       ENDIF
    ENDDO

  END FUNCTION last_c_leng


END MODULE boundary_2d
