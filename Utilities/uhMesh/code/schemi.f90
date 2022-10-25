!
      module schemi
!
!     +-------------------------------------------------------------+
      use list
      use eucl_delaunay
      use init_data
      use grid_types
      use wakes_gen
      use fork
      use interfer
!
      implicit none
!
      type(stratum), dimension(:), pointer, public:: layer
      type(box), dimension(:), pointer, public :: big_box
!
      integer, dimension(:), allocatable, public :: which_e
      integer, dimension(:), allocatable, public :: h_many_edg
      integer, dimension(:), allocatable, public :: delaun_edges
!
      integer :: done_layers
      public :: wait
!     +-------------------------------------------------------------+
!
!
      contains
!
!
      subroutine structured_mesh ( grid, z )
!
!     +--------------------------------------------------------------------+
!     |   subroutine dedicata a:                                           |         
!     |   organizzazione dei lati contigui nelle liste                     |       
!     |   richiesta della legge di variazione delle altezze                |       
!     |   inserimento ricorsivo degli strati di punti                      |       
!     |   numerazione dei punti inseriti                                   |       
!     +--------------------------------------------------------------------+
      implicit none
!
      type(domain), intent(inout) :: grid
      integer, intent(out) :: z
!
      type(point), pointer :: h, h_new, ptmp
      type(point), pointer :: p, head, tail
      type(point), pointer :: virt_h_new
!
      type(box), dimension(:), pointer :: pnt2bigbox                      
      type(box), pointer :: continent, pnt2box              
!
      integer, dimension(:), allocatable :: leng_mod, leng_var
      integer, dimension(:), allocatable :: cross_ovlp, auto_ovlp           
!
      integer :: i, j, n
      integer :: r, s, t
      integer :: nlayer, insert_mtd, num_small_boxes
      integer ::del, size_big_box, min_pt, enlarge, res                      
!
      real(kind=8), dimension(:), allocatable  :: delta, leng_par
      real(kind=8) :: gamma2, gamma3                  
      real(kind=8), parameter :: pi = 3.1415927464102d0          
!     +----------------------------------------------------------------------+
!
      number = 0
      nlayer = grid % number_of_layers
      n = grid % number_of_edges
!
      allocate ( delta(n) )       
      allocate ( leng_mod(n) )    
      allocate ( leng_var(n) )    
      allocate ( leng_par(n) )    
      allocate ( cross_ovlp(n) )
      allocate ( auto_ovlp(n) )
!
      delta = 0
      leng_mod = 0
      leng_var = 0
      leng_par = 0
      cross_ovlp = 0
      auto_ovlp  = 0
!
!     +-----------------------------------------------------------------------------+
!     |  Organizzazione dei contorni in liste:                                      |
!     |  z = numero di contorni indipendenti, numero di liste per strato            |
!     |  Initial_box_setup:                                                         |  
!     |  il lato continente e' il contorno esterno del dominio; tutti i lati        |
!     |  interni sono affetti da interferenza di box nei sui confronti, se chiuso   |
!     |  su se stesso, per cui bisogna frammentarlo in piccoli box da includere     |
!     |  nella lista di big_box iniziale; se si chiude su se' stesso ma e' composto |   
!     |  nel file boundary da piu lati, non e' necessaria la suddivisione in        |
!     |  small box                                                                  |
!     +-----------------------------------------------------------------------------+
!
      call organize_bou ( nlayer+1, grid, z, delta, leng_mod, &
        leng_var, leng_par, cross_ovlp, auto_ovlp )
!
      r = 0
      nullify ( continent )
!
!     +---------------------+
!     |  Self-closing edge  |
!     +---------------------+
!
      if ( grid % edges(1) % begin_vert .eq. &
        grid % edges(1) % end_vert ) then
!
        continent => new_box()
!
!     +-----------------------+
!     |  Non-structured edge  |
!     +-----------------------+
!
        if ( grid % edges(1) % struct .eq. 0 ) then
!
          call jumpy ( grid % edges(1) % add_head )
!
          continent % head => grid % edges(1) % add_head
          continent % c1 => grid % edges(1) % add_head % next
          continent % c2 => grid % edges(1) % add_head % prev   
          continent % cross_ovlp = grid % edges(1) % cross_ovlp
          continent % auto_ovlp  = grid % edges(1) % auto_ovlp
!
!     +-------------------+
!     |  Structured edge  |
!     +-------------------+
!
        else
!
          head => layer(1) % pnt_list(1) 
          call jumpy ( head ) 
          nullify ( head )
!
          continent % head => layer(1) % pnt_list(1)
          continent % c1 => layer(1) % pnt_list(1)
          continent % c2 => layer(1) % pnt_list(1) % prev 
          continent % cross_ovlp = cross_ovlp(1)
          continent % auto_ovlp  = auto_ovlp(1)
!
        end if
!
        if ( grid % edges(1) % struct .eq. 0 ) then
          continent % structured = 0
        else
          ! Vale per un solo lato strutturato che in piu' si chiude su se' stesso
          ! non posso attribuire il valore del numero di strati, ma solo 1/0
          continent % structured = 1 
        end if

        if ( continent % structured .eq. 0 ) then
           r = 1
           min_pt = 5
           call big2small_box ( continent, min_pt, num_small_boxes )
        else
           num_small_boxes = 1
        end if
!
!     +-------------+
!     |  Open edge  |
!     +-------------+
!
      else
!
         ! annulla il -1 che si trova calcolando size_big_box
         num_small_boxes = 1 
!
      end if

    del = 0
    DO t=1,SIZE(delaun_edges)
       IF(delaun_edges(t) /= 0) del = del + 1
    ENDDO
    ! essendo noto il numero complessivo di box, i lati originali con l'eventuale
    ! frazionamento del continente, si alloca l'array di box big_box

    size_big_box = z+del+num_small_boxes-1 

    ALLOCATE(big_box(size_big_box))
    DO t=1,size_big_box
       big_box(t) = new_box()
    ENDDO

    ! big_box e' organizzato come una lista ordinata di box; prima vengono
    ! quelli strutturati, poi quelli delaunay. Gli small box del continente
    ! vanno messi tra i big e per primi nei rispettivi gruppi(strutt o del)
    IF(num_small_boxes >= 1)THEN
       !IF(num_small_boxes > 1)THEN
       IF(ASSOCIATED(continent))THEN
          IF(continent % structured == 1)THEN
             !write(*,*)'continente strutturato'
             big_box(1) % structured = 0
             IF(num_small_boxes > 1)THEN! condizione solo per delaunay
                !write(*,*)'$$$$'
                DO s=2, num_small_boxes 
                   big_box(s) = continent % small_box(s)
                   big_box(s) % structured = 1
                ENDDO
             ENDIF
          ELSE
             !write(*,*) 'continente da',z+1,'a',z+num_small_boxes
             DO s=1,num_small_boxes
                big_box(s+z) = continent % small_box(s)
                big_box(s+z) % structured = 2 
             ENDDO
          ENDIF
       ENDIF
    ENDIF
    ! ora vengono definiti i big_box dei lati strutturati e delaunay, facendo
    ! attenzione a saltare gli eventuali small_box del continente appena
    ! inserito; quelli delaunay sono riconoscibili dal tag structured = 2
    t=1
    !write(*,*)'zeta', z
    DO s=1,size_big_box

       IF( big_box(s) % structured == 1) CYCLE

       head => layer(1) % pnt_list(t) 
       CALL jumpy(head)

       head => layer(1) % pnt_list(t) 
       tail => layer(1) % pnt_list(t) % prev
       enlarge = 1

       CALL define_struct_box(head, tail, big_box(s), enlarge)
       big_box(s) % c1 => head
       big_box(s) % c2 => tail
       big_box(s) % cross_ovlp = cross_ovlp(t)
       big_box(s) % auto_ovlp = auto_ovlp(t)

       t=t+1
       IF (t>z) EXIT
    ENDDO

    t=s+1 !r=1!r=0
    r = 0 
    IF(ASSOCIATED(continent))THEN
       IF(continent % structured == 0) r = 1
    ENDIF
    IF(size_big_box > t)THEN
       DO s=t, size_big_box
          IF(big_box(s) % structured == 2) THEN
             big_box(s) % structured = 0
             CYCLE
          ENDIF
          r = r + 1
          !write(*,*)'delaun_edges(r)',r,delaun_edges(r)
          CALL jumpy(grid % edges(delaun_edges(r)) % add_head)

          head => grid % edges(delaun_edges(r)) % add_head % next
          tail => grid % edges(delaun_edges(r)) % add_head % prev

          CALL define_del_box(head, tail, big_box(s))
          NULLIFY(big_box(s) % head)
          big_box(s) % head => head % prev ! e' una correzione necessaria 
          ! per avere un fit di define_del_box alle operazioni sul level = 1
          big_box(s) % c1 => head
          big_box(s) % c2 => tail
          big_box(s) % cross_ovlp = grid % edges(delaun_edges(r)) % cross_ovlp
          big_box(s) % auto_ovlp = grid % edges(delaun_edges(r)) % auto_ovlp

          IF(r>SIZE(delaun_edges)) EXIT
       ENDDO
    ENDIF
!
!     +--------------------------------------+
!     |  Loop over the boundary nodes with   |
!     |  boundary layer option               |
!     +--------------------------------------+
!
    done_layers = 0
    DO j = 1, z

       h => layer(1) % pnt_list(j)
!
!     +--------------------------------------------+
!     |  In sequence do:                           |
!     |  1. calcolo distanze fra punti consecutivi |
!     |     del contorno, serve da base per la     |
!     |     lunghezza di riferimento               |
!     |  2. lungh. max di riferim. per ogni punto  |
!     +--------------------------------------------+
!
       CALL base_leng ( h )
       CALL calc_leng_max ( leng_mod(j), h )
!
!     +---------------------------------------------+
!     |  La distanza tra i nodi viene trasferita    |
!     |  immutata ai punti degli strati successivi  |
!     +---------------------------------------------+
!
       p => h
       DO 
          p % cnd = 0
          p % bou_leng1 = p % leng
          p => p % next
          IF (ASSOCIATED(p,h)) EXIT
       END DO
!
!     +-----------------------------------------------------------------------+
!     |  Calcola il numero di strati necessario per avere chiusure e          |
!     |  assegna il tag cond  per effettuare eventuali interventi correttivi  |
!     +-----------------------------------------------------------------------+
!
       CALL find_fork ( h, 1, nlayer, delta(j) )
!
    END DO
    done_layers = 1
!
    res = 0
    pnt2bigbox => big_box
!
!     +------------------------------------------+
!     |  controllo di interferenza di base       |
!     |  utile per eventuali errori in boundary  |
!     +------------------------------------------+
!
    CALL check_interference ( pnt2bigbox, res ) 
                                                
    IF ( res .gt. 0 ) THEN
      WRITE (*,*)
      WRITE (*,*) 'ERROR. Boundaries overlap.'
      WRITE (*,*) 'Stop after first stage check_interference.'
      WRITE (*,*); STOP
    ENDIF
!
!     +----------------------------------------------------+
!     |  Computing the advancing front for boundary layer  |
!     +----------------------------------------------------+
!
    gamma2 = grid % gamma2
    gamma3 = grid % gamma3
    insert_mtd = grid % insert_mtd
!
    WRITE (*,*) '   Quadrilaterals advancing front...'

    DO i = 1, nlayer   
       DO j = 1, z

          h => layer(i) % pnt_list(j)         ! testa di lista dello strato attuale
          h_new => layer(i+1) % pnt_list(j)   ! testa di lista del nuovo strato

          layer(i+1) % max_lay(j) = layer(i) % max_lay(j)
          IF ( h % opnd .eq. 1 )  h_new % opnd = h % opnd  ! da eliminare?     

          IF ( i .eq. 1 ) THEN
            CALL first_layer ( grid, h, h_new, insert_mtd )
          ELSE
            CALL other_layer ( h, h_new, insert_mtd )
          END IF
!
!      +------------------------------------------------------------+
!      |  In order do:                                              |
!      |  1. Numerazione (anche per scie)                           |
!      |  2. Calcolo delle lunghezze dei segmenti del nuovo strato  |
!      |  3. Calcolo altezze di inserimento                         |
!      |  4. Aggiornamento della lunghezza massima di riferimento   |
!      +------------------------------------------------------------+
!
          CALL numbers ( h )
          CALL base_leng ( h_new )
          CALL update_height ( h_new, delta(j) )
          CALL update_leng_max ( h_new, leng_var(j), leng_par(j), i, nlayer )
!
          if ( i .gt. 1 )  then
!
!      +-----------------------+
!      |  CONCAVITY treatment  |
!      +-----------------------+
!          
            CALL condense ( h_new, gamma3, delta(j) )
!
            if ( i .le. nlayer ) then 
             
               p => h_new
               do 
                  call angle( p % prev, p, p % next, p % theta )
                  p => p % next
                  if ( associated(p, h_new) ) exit
               end do
!
!      +-----------------------+
!      |  CONVEXITY treatment  |
!      +-----------------------+
!
               call rarefactions( h_new )
                       
            end if
          end if
!
!      +-------------------------------+
!      |  CONCAVITY treatment forecast |
!      +-------------------------------+
!
          ! Previsione delle future intersezioni sul fronte concavo
          CALL find_fork( h_new, i, nlayer, delta(j) )
          
          ! Blocco eventuale delle altezze di inserimento          
          CALL check_stretch( h_new, delta(j), gamma2 )

          ! Connettivita' impiegata per il controllo di interferenza
          CALL jump_and_mark( h_new )
          
          IF ( i .eq. 1 ) THEN
          
             virt_h_new => h_new
             big_box(j) % virt_new_head => h_new
             
             IF ( h_new % opnd == 1 )  CALL update_jump_ar(h_new)

          ELSE          
	  
             IF (h % mark == 0)  big_box(j) % virt_new_head => h_new

             ! E' una testa che puo' non essere aggiornata con l'avanzamento del fronte
             ! se interessata da inteferenza
             virt_h_new => big_box(j) % virt_new_head

             ! h_new eredita le informazioni relative all'interferenza dello strato 
             ! precedente, i jump devono essere modificati di conseguenza             
             IF (i < layer(i) % max_lay(j)+1) then
	       CALL update_jump_ar ( h_new )
	     end if
             
          ENDIF

          pnt2box => big_box(j)
          IF ( i > layer(i) % max_lay(j)+1 ) THEN
             ! Do nothing...
          ELSE
             ! I box sono basati sulle liste parallele (rispetto alle liste
             ! next/prev) dei jump          
             CALL update_big_box(pnt2box, virt_h_new)
          ENDIF

       ENDDO ! sui diversi lati

       res = 0
       pnt2bigbox => big_box

       ! Controllo di interferenza dei fronti 0=self_check disabilitato
       CALL check_interference(pnt2bigbox, res)
       
       DO j = 1, z
          IF ( i .gt. layer(i) % max_lay(j) ) THEN

             ! Superato il numero di strati prescritto per il lato, si ipotizza per i
             ! successivi un'interferenza che ne coinvolge tutti i punti
             h_new => layer(i+1) % pnt_list(j)
             p => h_new

             DO
                p % mark = 1
                p => p % next 
                IF (ASSOCIATED(p, h_new)) EXIT
             ENDDO
             
             res = 1
             
          ENDIF
       ENDDO
       
       
       ! Interrompe forzatamente e indirettamente il proseguimento del fronte,
       ! simulando una interferenza che permette di attribuire a ciascun lato un
       ! numero diverso di strati da inserire
       IF ( res .ge. 1 ) THEN
          
          IF ( i .eq. 1 ) THEN
             WRITE(*,*)''          
             WRITE(*,*)'ERROR. Overlapping detected after FIRST layer insertion:'
             WRITE(*,*)'Try increasing the gaps between boundaries or the elements'' height.'
             WRITE(*,*)''
             STOP
          ENDIF
          
          DO j = 1, z
          
             h_new => layer(i+1) % pnt_list(j)
             h     => layer(i)   % pnt_list(j)
             
             IF ((h_new % mark == 1) .AND. &
                 (ASSOCIATED(h_new, big_box(j)%virt_new_head))) &
                
                big_box(j) % virt_new_head => h

             CALL update_jump_ar(h_new)
          
          ENDDO
          
          done_layers = i + 1 ! togliere e ripristinare ! !
          
       ELSE
          done_layers = i + 1
       ENDIF



       IF ( i .eq. nlayer ) THEN
          
          DO j = 1, z
          
             h_new => layer(i+1) % pnt_list(j)
             
             ! Assegnazione degli indici sull'ultimo
             ! strato di punti immessi
             CALL numbers(h_new)                 

             head  => big_box(j) % virt_new_head
             p => head

             DO
                p % mark = 0
                p => p % jump_f 
                IF (ASSOCIATED(p, head)) EXIT
             ENDDO

          ENDDO

       ENDIF

    ENDDO   ! Loop over layers

  END SUBROUTINE structured_mesh



  !-------------------------------------------------------------------------
  SUBROUTINE jump_and_mark(head)

    IMPLICIT NONE
    TYPE(point), POINTER :: head, p

    CALL jumpy(head)        
    ! duplicazione della lista next/prev in jump_f/jump_b
    ! mark = 1 => pto/segmento interessato da interferenza

    p => head
    DO                       
       IF(ASSOCIATED(p % old1))THEN
          p % mark = MAX( p % old1 % mark, p % old % mark)
       ELSE
          p % mark = p % old % mark
       ENDIF
       p => p % next; IF(ASSOCIATED(p, head))EXIT
    ENDDO
    ! i mark dello strato prec. vengono
    ! trasferiti alla lista h_new

  END SUBROUTINE jump_and_mark
  !-------------------------------------------------------------------------
  SUBROUTINE update_jump_ar(head)
    ! rileva i mark (punti interessati da interferenza del fronte)
    ! viene utilizzata due volte:
    ! - prima di effettuare il controllo di interferenza, per avere una lista-jump 
    !   con i vecchi mark ereditati aderente alla griglia 'erosa' (per il calcolo 
    !   corretto dei box)
    ! - dopo il controllo di interferenza, per aggiornare la lista-jump con
    !   l'introduzione dei nuovi mark

    ! viene implementata percorrendo la lista nei due versi, per trattare la zona a
    ! mark=1 sempre dall'esterno

    ! previsti due comportamenti quando si verifica interferenza: alcuni elementi
    ! dello strato scompaiono e lasciano un gradino scoperto oppure spianato da un
    ! triangolo
    IMPLICIT NONE

    TYPE(point),POINTER :: head,p,q
    INTEGER :: mark, nxtmark, what


    p => head                             ! ciclo di andata, segue i next

    IF((p % opnd == 1).AND.(p % mark == 0))THEN
       p % old % what = 2
       p % jump_b => p % old
       p % jump_b % jump_f => p
       p => p % next
    ENDIF
    mark = p % mark
    DO
       nxtmark = p % next % mark

       IF(mark < nxtmark)THEN              ! si va incontro ad una zona con mark=1
          IF(p % prev % mark > mark)THEN
             ! c'e' un buchetto che conviene riempire...
             p % mark = 1
             p => p % next; IF(ASSOCIATED(p,head))EXIT
             mark = nxtmark
             CYCLE
          ENDIF
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF(.NOT.ASSOCIATED(p % old1))THEN
             IF(p % old % what == 10)THEN
                ! il mark e' in mezzo ad un quaTRIlatero
                ! ora sono posizionato all'inizio del quaTRIlatero
                p % old % what = 0
                IF(p % old % prev % what == 1)THEN
                   ! quaTRIlatero preceduto da un rettangolo
                   p % old % prev % what = 3 ! qudrilatero exp a SX
                ELSEIF(p % old % prev % what == 10)THEN
                   ! quaTRIlatero preceduto da un quaTRIlatero
                   p % old % prev % what = 17! quaTRIlatero esposto a SX
                ELSE
                   !    WRITE(*,*)'uncovered upd_jump_ar+ at p',p%x
                   !write(*,*)'old prev what=', p%old%prev%what,p%old%prev%index
                   !    write(*,*)'quaTRIlater0 preceduto da qcsa /=da quad/quaTRIlat.'
                   !   STOP
                ENDIF
                ! bisogna eliminare il punto introdotto per rarefazione
                q => p % next
                CALL delete_point(q)
             ENDIF
          ELSE
             IF((p % old1 % what == 10).OR.(p%old1%what==17))THEN
                ! il mark alla fine del quaTRIlatero
                ! ora sono posizionato in mezzo al quaTRIlatero
                q => p 
                p => p % prev
                p % old % what = 0
                IF(p % old % prev % what == 1)THEN
                   ! quaTRIlatero preceduto da un rettangolo
                   p % old % prev % what = 3 ! quadrilatero esposto a SX
                ELSEIF(p % old % prev % what == 10)THEN
                   ! quaTRIlatero preceduto da un quaTRIlatero
                   p % old % prev % what = 17! quaTRIlatero esposto a SX
                ELSE
                   !WRITE(*,*)'uncovered update_jump_ar+ at p',p%x
                   !write(*,*)'old prev what=', p%old%prev%what,p%old%prev%index
                   !STOP
                ENDIF
                ! bisogna eliminare il punto introdotto per rarefazione
                CALL delete_point(q)
             ENDIF
          ENDIF
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          p % jump_f => p % old 
          p % jump_f % jump_b => p 
          p % what = 0
          p % jump_f % what = 0
          what = p % jump_f % prev % what
          IF(what == 10)THEN
             p % jump_f % prev % what = 17 
          ELSEIF(what == 1)THEN
             p % jump_f % prev % what = 3 
          ELSEIF(what == 2)THEN
             p % jump_f % prev % what = 15
          ELSE
             ! WRITE(*,*)'uncovered update_jump_ar++ at p',p % x
             !       write(*,*)'jump_f prev what=',what,p%jump_f%prev%index
          ENDIF
       ELSEIF((mark > 0).AND.ASSOCIATED(p % jump_f, p % old))THEN 
          IF(ASSOCIATED(p % old % jump_b, p))THEN
             p % old % jump_b => p % old % prev
             p % jump_f => p % next
             p % what = 0
             p % old % jump_b % what = 0
          ENDIF
       ENDIF
       p => p % next
       IF(ASSOCIATED(p,head))EXIT
       mark = nxtmark
    ENDDO
    !-------------------------------------
    p => head                             ! ciclo di ritorno, segue i prev
    IF((p % opnd == 1).AND.(p % prev % mark == 0))THEN
       p => p % prev
       p % old % what = 0
       p % jump_f => p % old
       p % jump_f % jump_b => p
       what = p % old % prev % what
       IF(what == 2)THEN
          p % prev % old % what = 15
       ELSEIF(what == 1)THEN
          p % prev % old % what = 3
       ELSEIF(what == 3)THEN

       ELSE
          write(*,*) ''
	  write(*,*) ' ERROR: Exception found in rarefaction for quads.'
	  write(*,*) ' Try increasing the ratio of up/low to trigger triangles insertion'
	  write(*,*) ''
	  stop
       ENDIF
       p => p % prev
    ENDIF
    mark = p % mark
    DO
       nxtmark = p % prev % mark

       IF(mark < nxtmark)THEN          ! si va incontro ad una zona con mark=1

          !----------------------------------------------------------------------
          IF(.NOT.ASSOCIATED(p % old1))THEN
             !IF(p % old % prev % what == 10)THEN
             IF((p % old % prev % what == 10).OR.(p % old % prev % what == 16))THEN
                ! il mark e' in mezzo al quaTRIlatero
                ! ora sono posizionato alla fine del quaTRIlatero 
                p % old % prev % what = 0
                what = p % old % what
                IF(what == 1)THEN
                   p % old % what = 2
                ELSEIF(what == 3)THEN
                   p % old % what = 15
                ELSEIF(what == 10)THEN
                   p % old % what = 16
                ELSEIF(what == 17)THEN
                   p % old % what = 18
                ELSE
                   !  write(*,*)'uncovered update_jump_ar- at p',p % x
                   !  write(*,*)'old what=', p%old%what,p%old %index
                ENDIF
                q => p % prev
                CALL delete_point(q)
             ENDIF
          ELSE
             IF((p % old1 % what==10).OR.(p%old1 % what==16))THEN
                ! il mark e' all'inizio del quaTRIlatero
                ! ora sono posizionato in mezzo al quaTRIlatero
                q => p
                p => p % next ! torno 'indietro'
                p % old % prev % what = 0
                what = p % old % what
                IF(what == 1)THEN
                   p % old % what = 2
                ELSEIF(what == 3)THEN
                   p % old % what = 15
                ELSEIF(what == 10)THEN
                   p % old % what = 16
                ELSEIF(what == 17)THEN
                   p % old % what = 18
                ELSE
                   !   write(*,*)'uncovered update_jump_ar- at p',p % x
                   !   write(*,*)'old what=', p % old % what, p % old % index
                ENDIF
                CALL delete_point(q)
             ENDIF
          ENDIF


          !----------------------------------------------------------------------

          p % jump_b => p % old 
          p % jump_b % jump_f => p 
          p % prev % what = 0
          what = p % jump_b % what
          IF(what == 1)THEN
             p % jump_b % what = 2
          ELSEIF(what == 2)THEN
          ELSEIF(what == 3)THEN
             p % jump_b % what = 15 ! e' un dentino!!!!!!!!!!!!!!!!!!!!
          ELSEIF(what == 10)THEN
             p % jump_b % what = 16
          ELSEIF(what == 16)THEN

          ELSEIF(what == 17)THEN
             p % jump_b % what = 18
          ELSEIF(what == 15)THEN

          ELSE
             ! WRITE(*,*)'uncovered update_jump_ar-- at p', p % x
             !  write(*,*)'jump_b what=',p%jump_b%what, p%jump_b%index
          ENDIF
          p % jump_b % prev % what = 0
       ELSEIF((mark > 0).AND.ASSOCIATED(p % jump_b, p % old))THEN 

          IF(ASSOCIATED(p % jump_b, p % old))THEN
             !       write(*,*)'  correzione-', p % x
             p % jump_b => p % prev
             p % old % jump_f => p % old % next
             p % what = 0
             p % old % what = 0
          ENDIF
       ENDIF
       IF((mark == nxtmark).AND.(mark == 1))THEN
          p % old % prev % what = 0
          IF(ASSOCIATED(p % old1))THEN 
             p % old1 % prev % what = 0 
             IF(ASSOCIATED(p % old % prev % prev, p % old1))THEN
                p % old1 % what = 0 
                ! e' per completare il caso dei due triangoli di condensazione affiancati
             ENDIF
          ENDIF
       ENDIF
       p => p % prev
       IF(ASSOCIATED(p,head))EXIT
       mark = nxtmark
    ENDDO
  END SUBROUTINE update_jump_ar
  !-------------------------------------------------------------------------
  SUBROUTINE delete_point(q)

    IMPLICIT NONE
    TYPE(point), POINTER :: p, q, r

    p => q % prev
    r => q % next

    p % next => r
    r % prev => p

    p % jump_f => r
    r % jump_b => p
    NULLIFY(p % old % ppp1)
    NULLIFY(r % old % ppp1)

    DEALLOCATE(q)

  END SUBROUTINE delete_point





  SUBROUTINE jumpy(head)
    ! copia la connettivita' di lista in jump
    IMPLICIT NONE

    TYPE(point),POINTER :: head,p

    p => head
    DO
       p % jump_f => p % next
       p % jump_b => p % prev
       p => p % next
       IF (ASSOCIATED(p, head))EXIT
    ENDDO

  END SUBROUTINE jumpy





  SUBROUTINE check_stretch ( head, delta, gamma )
!
! subroutine di blocco della crescita delle altezze di inserimento
! quando lo stretch degli elementi è poco adeguato all'esecuzione di
! condensazioni; sono comunque possibili condensazioni 'difficili' 
! a causa di chiusure ampie della griglia e vengono gestite da un ulteriore
! parametro di stretch 'locale'
!
    IMPLICIT NONE

    TYPE(point),POINTER :: head,p
    REAL(KIND=8),PARAMETER :: pi = 3.1415927464102d0
    REAL(KIND=8) :: stretch, max_stretch
    REAL(KIND=8),INTENT(IN) :: gamma
    REAL(KIND=8),INTENT(INOUT) :: delta

    p => head
    max_stretch = 1 / (2 * TAN(gamma*pi / 360))

    DO
       stretch = p % height / p % leng

       IF ( (p % theta .lt. 0) .AND. (stretch .gt. max_stretch) ) THEN

          delta = 0
          EXIT

       ENDIF
       p => p % next
       IF(ASSOCIATED(p,head))EXIT 

    ENDDO

  END SUBROUTINE check_stretch





  !! Assign an id to the points contained into the head list
  SUBROUTINE numbers ( head )

    IMPLICIT NONE
    TYPE(point),POINTER :: head, p, q

    ! assegnazione dell'indice anche alle spalle dei lati aperti!!!!!
    ! quindi Delaunay ne assegnerà un altro!!!!!! 

    p => head
    DO 
!
       IF ( p % mark .eq. 0 ) THEN
         number = number + 1
         p % index = number
       END IF
!
       IF ( ASSOCIATED(p % wake) ) THEN
          ! numerazione anche per il filamento di scia
          q => p % wake
          DO
             number= number + 1
             q % index = number
             IF(.NOT.ASSOCIATED(q % next))EXIT
             q => q % next
          END DO
       END IF
!
       p => p % next
       IF ( ASSOCIATED(p, head) ) EXIT

    END DO
!
  END SUBROUTINE numbers




   SUBROUTINE first_layer ( grid, first, head, insert_mtd )
!
!  This method inserts the first layer of structured cells on a given boundary;
!  the method differs from the other_layer method for specific boundary handling.
!
!  Input parameters are:
!  - grid: to retrieve the wake geometry data
!  - first: the pointer to the current layer point list
!  - head: the pointer to the next layer (empty) point list 
!  - insert_mtd: the insertion method specifier (along the normal or with adaptation)
!
!  Specific boundary handling consists of wake data insertion (depending on the p%option parameter
!  as it is given in the boundary.name file as 'wake type') and of bi/trifurcation (depending on 
!  the p%theta parameter)
!
   IMPLICIT NONE

   TYPE(domain), INTENT(INOUT):: grid
   TYPE(point), POINTER :: first, head
   INTEGER, INTENT(IN) :: insert_mtd
     
   TYPE(point), POINTER :: p, q, r, pp
   TYPE(point), POINTER :: s, t, u
   TYPE(point), POINTER :: wake_head

   REAL(KIND=8), DIMENSION(2,2) :: RM
   REAL(KIND=8), DIMENSION(2) :: q_med, qmed
   REAL(KIND=8), DIMENSION(2) :: q1, q_med1
   REAL(KIND=8), DIMENSION(2) :: q2, q_med2
   REAL(KIND=8), DIMENSION(2) :: rot_vec1
   REAL(KIND=8), DIMENSION(2) :: rot_vec2

   REAL(KIND=8) :: theta, theta1, theta2, te_alpha
   REAL(KIND=8) :: l1, l2, alpha, beta1, beta2
   REAL(KIND=8), PARAMETER    :: pi = 3.1415927464102d0
   real*8 :: mfangle

   INTEGER :: cnd, option
   CHARACTER(LEN=4) :: optn

   p => first ! primo punto dello strato esistente
   q => head  ! primo punto del nuovo strato

   DO
   
      theta = p % theta
      cnd = p % cnd 
      option = p % option 

      IF ( p % option .gt. 0 ) THEN

         ! Scia. Punto di inizio della scia sul primo strato di punti
         p % next % cnd = 0
         p % prev % cnd = 0
         p % what = 14

         IF (.NOT. ASSOCIATED(p, first)) THEN     

            ! Punti successivi al primo della lista
            q => new_point()
            optn = 'frst'
            CALL insert_next_evo ( optn, p, q, r, head )

         ELSE
         
            ! Primo punto della lista 
            q1 = p % nv1(1,:)  
            q1 = q1 / SQRT(DOT_PRODUCT(q1, q1))

            q % x = p % x + q1 * p % height 

         ENDIF

         p % ppp => q
         q % old => p
         q % cnd = 0

         ! Punto di completamento della scia sul primo strato di punti
         r => new_point()
         optn ='scnd'
         
         CALL insert_next_evo ( optn, p, r, r, head )

         p % ppp1 => r
         r % old => p
         r % cnd = 0

         q1 = p % nv1(1,:)
         q2 = p % nv1(2,:)
         
         q_med = q1

         ! Rotazione di q1 fino ad essere allineato con la bisettrice          
         mfangle = p % theta / 2.d0
         CALL rot(mfangle, q_med, q_med)

         ! Wake coordinates are relative to a reference point!
         ! assign the right wake coordinates adding the attachment 
         ! reference point (p) coordinates
         wake_head => grid % wakes(p % option) % add_head

         qmed = wake_head % next % next %x
         qmed = qmed / SQRT(DOT_PRODUCT(qmed, qmed));

         pp => wake_head % next

         DO
            pp % x = pp % x + p % x
            pp => pp%next
            IF (ASSOCIATED(pp, wake_head)) EXIT
         ENDDO

         CALL angle(p%prev, p, wake_head%next%next, beta1) 
         mfangle = -(PI+beta1)/2.d0
         CALL rot(mfangle, qmed, rot_vec1)          

         CALL angle(wake_head%next%next, p, p%next, beta2)
         mfangle = (PI+beta2)/2.d0 
         CALL rot( mfangle, qmed, rot_vec2)          

         q % x = p % x + rot_vec1*p%height/ ABS(COS( beta1/2))
         r % x = p % x + rot_vec2*p%height/ ABS(COS( beta2/2))

         CALL data_trail ( q, p, r, wake_head )


      ELSE IF ( (p % theta > (pi - 130*pi/180)) .AND. &
                (p % theta < (pi -  60*pi/180)) ) THEN  ! Triforcazione

         p % next % cnd = 0
         p % prev % cnd = 0
         p % what = 9

         IF ( .NOT. ASSOCIATED(p,first) ) THEN  ! p is NOT the first point of the list add two 
         ! points along the two local normal directions

            q => new_point()      
            r => new_point()

            optn = 'both'
            CALL insert_next_evo( optn, p, q, r, head )

         ELSE  ! p is the first point of the list then still add two points 
         ! (the difference with the above is not very clear)

            q1 = p % nv1(1,:)
            q1 = q1/SQRT(DOT_PRODUCT(q1,q1))
            q % x = p % x + q1 * p % height

            q2 = p % nv1(2,:)
            q2 = q2/SQRT(DOT_PRODUCT(q2,q2))
            r => new_point()

            optn = 'scnd'
            CALL insert_next_evo( optn, p, r, r, head )

         END IF
         
         p % ppp => q
         p % ppp1 => r

         ! The attribute 'old' refers to the corresponding point on the previous layer
         ! the correspondence is one-to-one if no bi-trifurcation is present

         q % old => p
         r % old => p

         q % cnd = 0
         r % cnd = 0

         ! The attribute 'nv1' gives the normal versors departing from point p. One
         ! is along the normal of the left adjacent element, the other is along the normal
         ! of the right adjacent element

         q1 = p % nv1(1,:)
         q2 = p % nv1(2,:)

         ! Insert a point between the new points q and r that are obtained  from
         ! the source point p laying on the previous layer

         q_med = (q1 + q2) / 2
         q_med = q_med / SQRT(DOT_PRODUCT(q_med,q_med)) 

         s => new_point()
         s % x = p % x + q_med * p % height

         CALL insert_between( q, r, s )
         s % old => p

         s % cnd = 0 
         q % cnd = 0
         r % cnd = 0

         p % bou_leng1 = SQRT( segment_length(2, q, s) )


      ELSE IF ( p % theta > (pi - 60*pi/180) ) THEN  ! Pentaforcazione


         p % next % cnd = 0
         p % prev % cnd = 0
         p % what = 19

         IF ( .NOT. ASSOCIATED(p,first) )THEN     

            q => new_point()      
            r => new_point()

            optn = 'both'
            CALL insert_next_evo(optn, p, q, r, head)

         ELSE

            q1 = p % nv1(1,:)
            q1 = q1/SQRT(DOT_PRODUCT(q1,q1))
            q % x = p % x + q1 * p % height

            q2 = p % nv1(2,:)
            q2 = q2/SQRT(DOT_PRODUCT(q2,q2))
            r => new_point()

            optn = 'scnd'
            CALL insert_next_evo(optn, p, r, r, head)

         ENDIF
         
         p % ppp => q
         p % ppp1 => r

         q % old => p
         r % old => p

         q % cnd = 0
         r % cnd = 0

         q1 = p % nv1(1,:)          
         q2 = p % nv1(2,:)
         
         alpha = acos( DOT_PRODUCT(q1,q2) )
         
         mfangle = 0.25*alpha
         call rot( mfangle, q1, q_med )
         q_med = q_med / SQRT(DOT_PRODUCT(q_med,q_med))
         s => new_point()
         s % x = p % x + q_med * p % height
         call insert_between(q,r,s)
         s % old => p
         s % cnd = 0
         p % bou_leng1 = SQRT(segment_length(2,q,s))

         mfangle = 0.5*alpha
         call rot( mfangle, q1, q_med )
         q_med = q_med / SQRT(DOT_PRODUCT(q_med,q_med))
         t => new_point()
         t % x = p % x + q_med * p % height
         call insert_between(s,r,t)
         t % old => p
         t % cnd = 0
         p % bou_leng1 = p % bou_leng1 + SQRT(segment_length(2,s,t))

         mfangle = 0.75*alpha
         call rot( mfangle, q1, q_med )
         q_med = q_med / SQRT(DOT_PRODUCT(q_med,q_med))
         u => new_point()
         u % x = p % x + q_med * p % height
         call insert_between(t,r,u)
         u % old => p
         u % cnd = 0
         p % bou_leng1 = p % bou_leng1 + SQRT(segment_length(2,t,u))


      ELSE IF ( p % opnd > 0 ) THEN
                
         ! Inizio di un lato aperto
         p % what = 2    
         q => head
         q2 = p % nv1(2,:)
         q2 = q2/SQRT(DOT_PRODUCT(q2,q2))
         q % x = p % x + q2 * p % height
         p % ppp => q
         q % old => p
         
         
         
      ELSE IF ( p % next % opnd > 0 ) THEN
         
         ! Fine di un lato aperto
         p % what = 3
         q => new_point()
         optn = 'frst'
         CALL insert_next_evo(optn,p,q,r,head)

         p % ppp => q
         q % old => p

      ELSE

         ! Punto qualsiasi
         p % what = 1
         
         IF (.NOT. ASSOCIATED(p, first)) THEN     
         
            q => new_point()
            optn = 'half'

            CALL insert_next_evo(optn,p,q,r,head)

         ELSE
         
            q1 = p % nv1(1,:)
            q2 = p % nv1(2,:)

             
             IF ( ( p % x(1) .le. 0.725  .and. p % x(1) .ge. 0.69  .and. p % x(2) .le.  0.04   .and. p%x(2) .ge.  0.013 ) ) THEN !.or. & ! cove
!                  ( p % x(1) .le. 1.137  .and. p % x(1) .ge. 1.127 .and. p % x(2) .le. -0.141  .and. p%x(2) .ge. -0.148 ) .or. & ! blunt
!                  ( p % x(1) .le. 0.058  .and. p % x(1) .ge. 0.043 .and. p % x(2) .le. -0.0281 .and. p%x(2) .ge. -0.038 ) .or. & ! mainSlat
!                  ( p % x(1) .le. 0.887  .and. p % x(1) .ge. 0.870 .and. p % x(2) .le.  0.032  .and. p%x(2) .ge.  0.02) ) THEN   ! flapMain
!             
                q_med = (q1 + q2)/2
                q_med = q_med/(SQRT(DOT_PRODUCT(q_med, q_med)))
! 
                IF ( p % theta /= 0 ) THEN
                   q % x = p % x + q_med * p % height/ABS(COS(p % theta/2))   
                ELSE
                   q % x = p % x + q_med * p % height
                ENDIF            
!            
            ELSE IF ( insert_mtd .eq. 1 ) THEN !optn = give
!            IF ( p % theta .gt. pi - 130*pi/180.0 ) THEN !optn = give
            
               !theta1 = (p % theta - p % prev % theta)/3 
               !theta2 = (p % theta - p % next % theta)/3 
               theta1 = 20*pi/180d0
               theta2 = theta1
               CALL rot(theta1,q1,q_med1) 
               CALL rot(-theta2,q2,q_med2)

               l1 = p % prev % leng
               l2 = p % leng
               !alpha = 1+(2*(l1-l2)/(l1+l2))
               !q_med = ((alpha/2)**2)*q_med1+(1-((alpha/2))**2)*q_med2
               alpha =(2*(l1-l2)/(l1+l2))

               q_med = (1-alpha)*q_med1 + (1 + alpha)*q_med2

               q_med = q_med/(SQRT(DOT_PRODUCT(q_med,q_med)))
               q % x = p % x + q_med * p % height / ABS(COS((theta1+theta2)/2))
            ELSE
               q_med = (q1 + q2)/2
               q_med = q_med/(SQRT(DOT_PRODUCT(q_med, q_med)))

               IF(p % theta /= 0)THEN
                  q % x = p % x + q_med * p % height/ABS(COS(p % theta/2))   
               ELSE
                  q % x = p % x + q_med * p % height
               ENDIF
            ENDIF
         
         ENDIF
         
         p % ppp => q
         q % old => p

         IF ( p % cnd > 0 ) THEN
            ! Diminuzione per aperture (non piu' attivata)  
            q % cnd = p % cnd - 1
         ELSE IF ( p % cnd < 0 ) THEN
            ! Incremento per chiusure (-1 => condensazione)
            q % cnd = p % cnd + 1
         END IF

      END IF

      p => p % next
      IF ( ASSOCIATED(p, first) ) EXIT

   END DO

   END SUBROUTINE first_layer






    SUBROUTINE other_layer (first,head,insert_mtd)
!
!   Input parameters are:
!   - first: the pointer to the current layer point list
!   - head: the pointer to the next layer (empty) point list 
!   - insert_mtd: the insertion method specifier (along the normal or with adaptation)
!
    IMPLICIT NONE

    TYPE(point), POINTER       :: p, q, r
    TYPE(point), POINTER       :: first, head    

    REAL(KIND=8)               :: theta, te_alpha, theta1,theta2,l1,l2,alpha
    REAL(KIND=8), PARAMETER    :: pi = 3.1415927464102d0
    REAL(KIND=8), DIMENSION(2) :: q1, q2, q_med,q_med1,q_med2

    INTEGER                    :: cnd, option
    INTEGER,INTENT(IN)         :: insert_mtd

    CHARACTER(LEN=4) :: optn
!
    p => first  ! primo punto dello strato esistente
    q => head   ! primo punto del nuovo strato

    DO

       theta = p % theta
       cnd = p % cnd 
       option = p % option  ! indica un intervento particolare da svolgere sul contorno

       IF ( p % opnd > 0 ) THEN
       
          ! Inizio di un lato aperto
          p % what = 2    
          q => head
          q2 = p % nv1(2,:)
          q2 = q2/SQRT(DOT_PRODUCT(q2,q2))
          q % x = p % x + q2 * p % height
          p % ppp => q
          q % old => p
       
       ELSE IF ( p % next % opnd > 0 ) THEN

          ! Fine di un lato aperto
          p % what = 3
          q => new_point()
          optn = 'frst'
          CALL insert_next_evo(optn,p,q,r,head)

          p % ppp => q
          q % old => p

       ELSE

          ! Punto qualsiasi
          p % what = 1
          
          IF (.NOT. ASSOCIATED(p,first)) THEN

             IF ( ( p % x(1) .le. 0.725  .and. p % x(1) .ge. 0.69  .and. p % x(2) .le.  0.04   .and. p%x(2) .ge.  0.013 ) ) THEN !.or. & ! cove
!                  ( p % x(1) .le. 1.137  .and. p % x(1) .ge. 1.127 .and. p % x(2) .le. -0.141  .and. p%x(2) .ge. -0.148 ) .or. & ! blunt
!                  ( p % x(1) .le. 0.046  .and. p % x(1) .ge. 0.037 .and. p % x(2) .le. -0.0229 .and. p%x(2) .ge. -0.0276 ) .or. & ! mainSlat
!                  ( p % x(1) .le. 0.887  .and. p % x(1) .ge. 0.870 .and. p % x(2) .le.  0.032  .and. p%x(2) .ge.  0.02) ) THEN   ! flapMain
!             
                 optn = 'half'
!             
             ELSE IF ( insert_mtd .eq. 1 ) THEN
!             IF ( insert_mtd .eq. 1 ) THEN
                optn = 'give' ! opzione inserimento fuori dalla normale
             ELSE
                optn = 'half'
             ENDIF
             
             q => new_point()     
             CALL insert_next_evo ( optn, p, q, r, head )

          ELSE
          
             q1 = p % nv1(1,:)
             q2 = p % nv1(2,:)

             IF ( ( p % x(1) .le. 0.725  .and. p % x(1) .ge. 0.69  .and. p % x(2) .le.  0.04   .and. p%x(2) .ge.  0.013 ) ) THEN !.or. & ! cove
!                  ( p % x(1) .le. 1.137  .and. p % x(1) .ge. 1.127 .and. p % x(2) .le. -0.141  .and. p%x(2) .ge. -0.148 ) .or. & ! blunt
!                  ( p % x(1) .le. 0.046  .and. p % x(1) .ge. 0.037 .and. p % x(2) .le. -0.0229 .and. p%x(2) .ge. -0.0276 ) .or. & ! mainSlat
!                  ( p % x(1) .le. 0.887  .and. p % x(1) .ge. 0.870 .and. p % x(2) .le.  0.032  .and. p%x(2) .ge.  0.02) ) THEN   ! flapMain
!             
                 q_med = (q1 + q2)/2
                 q_med = q_med/(SQRT(DOT_PRODUCT(q_med, q_med)))
! 
                 IF ( p % theta /= 0 ) THEN
                    q % x = p % x + q_med * p % height/ABS(COS(p % theta/2))   
                 ELSE
                    q % x = p % x + q_med * p % height
                 ENDIF
!            
             ELSE IF ( insert_mtd .eq. 1 ) THEN
!             IF ( insert_mtd .eq. 1 ) THEN
             
                theta1 = 20*pi/180d0
                theta2 = theta1
                CALL rot(theta1,q1,q_med1) 
                CALL rot(-theta2,q2,q_med2)

                l1 = p % prev % leng
                l2 = p % leng

                alpha =(2*(l1-l2)/(l1+l2))

                q_med = (1-alpha)*q_med1 + (1 + alpha)*q_med2

                q_med = q_med/(SQRT(DOT_PRODUCT(q_med,q_med)))
                q % x = p % x + q_med * p % height / ABS(COS((theta1+theta2)/2))
             
             ELSE
             
                q_med = (q1 + q2)/2
                q_med = q_med/(SQRT(DOT_PRODUCT(q_med, q_med)))

                IF ( p % theta /= 0 ) THEN
                   q % x = p % x + q_med * p % height/ABS(COS(p % theta/2))   
                ELSE
                   q % x = p % x + q_med * p % height
                ENDIF
             
             ENDIF
          ENDIF
          
          p % ppp => q
          q % old => p

          IF ( p % cnd > 0 ) THEN
             q % cnd = p % cnd - 1  ! diminuzione per aperture( non piu' attivata)  
          ELSE IF ( p % cnd < 0 ) THEN
             q % cnd = p % cnd + 1  ! incremento per chiusure (-1 => condensazione)
          ENDIF

       ENDIF

       p => p % next
       IF ( ASSOCIATED(p, first) ) EXIT

    ENDDO

  END SUBROUTINE other_layer


  !°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

  SUBROUTINE organize_bou (nlayer,grid,z,delta,mod,var,par,cross,auto)

    IMPLICIT NONE

    INTEGER                     :: i, j, k, m, t, e, bv, bv1, ev, ev1, s, sum
    INTEGER , INTENT(IN)        :: nlayer
    INTEGER , INTENT(OUT)       :: z
    TYPE(domain), INTENT(INOUT) :: grid
    TYPE(point), POINTER        :: add_head, head, p, q
    TYPE(point), POINTER        :: pnt_list
    LOGICAL                     :: found
    REAL(KIND=8)                :: last_height, ref_delta, ref_mod, ref_var, &
         ref_par
    REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: save_delta

    INTEGER, INTENT(INOUT),DIMENSION(:)      :: mod, var, cross, auto 
    REAL(KIND=8),INTENT(INOUT), DIMENSION(:) :: par, delta
    INTEGER :: ref_cross, ref_auto 
    !............................................
    ! procedura di collegamento dei lati consecutivi che si vogliono grigliare in
    ! modo strutturato
    ! si effettua trovando la coincidenza dei vertici di inizio e di fine lato

!    WRITE(*,*)'   ...initial front description '
    ALLOCATE(layer(nlayer))
    ! array di strati

    m = grid % number_of_edges

    IF (ALLOCATED(h_many_edg)) DEALLOCATE(h_many_edg)
    ALLOCATE(h_many_edg(m))
    h_many_edg = 0
    ! individuati i lati indipendenti, h_many_edg informa di QUANTI lati originali
    ! sono costituiti

    IF (ALLOCATED(which_e)) DEALLOCATE(which_e)
    ALLOCATE(which_e(m))
    which_e = 0
    ! informa di QUALI lati originali sono costituiti  

    DO i=1,m

       grid % edges(i) % status = 0
       ! individua(1) i lati gia' presi in considerazione

       grid % edges(i) % cl_open = 0
       ! individua(1) il lato che chiude un lato indipendente finale

       WRITE (*,'(4x,a6,1x,i2,1x,a9,i2,2x,a1,1x,i2)') '- edge', i, 'vertices ', &
         grid % edges(i) % begin_vert, '&', grid % edges(i) % end_vert

    ENDDO
!    WRITE(*,*)
    !k ordine progressivo dei lati
    !which_e corrispondente a livello di grid%edges() per k
    ! risponde alla domanda:quali lati strutturati?(in ordine di percorrenza)
    !h_many_edg(m) quanti lati costituiscono il contorno chiuso 
    !(m dimensione massima)
    ! risponde alla domanda: quanti lati consecutivi strutturati per l'i-esimo 
    !                        contorno indipendente?

    k = 0

    z = 0
    !individua il lato finale indipendente

    DO i=1,m
       IF ((grid%edges(i)%struct >= 1).AND.(grid%edges(i)%status == 0))THEN

          k = k+1
          which_e(k) = i !grid%edges(e)....

          bv = grid % edges(i) % begin_vert
          ev = grid % edges(i) % end_vert
          found = .FALSE.      
          DO j=1,m

             IF((grid%edges(j)%struct >= 1).AND.(grid%edges(j)%status == 0)) THEN

                bv1 = grid % edges(j) % begin_vert
                ev1 = grid % edges(j) % end_vert

                IF ((i == j))THEN         
                   grid % edges(i)% status = k
                   IF (bv1 == ev1)THEN !self-closing edge
                      z=z+1
                      h_many_edg(z)=1
                      grid % edges(i) % status = k 
                      !WRITE(*,*)'   ...closed & structured boundary',z,'is made of', &
                      !          h_many_edg(z),'edge(s)'
                      found = .TRUE.
                      EXIT
                   ENDIF
                ENDIF

                IF (bv1 == ev)THEN  !next edge


                   k=k+1
                   grid % edges(j) % status = k
                   which_e(k)=j          
                   ev = grid % edges(j) % end_vert 
                   IF (ev1 == bv)THEN !next edge is also the final edge


                      z = z+1
                      sum=0

                      DO s=1,z-1
                         sum = sum + h_many_edg(s)
                      ENDDO
                      found=.TRUE.
                      h_many_edg(z)= k-sum
                      !WRITE(*,*)'   ...closed & structured boundary',z,'is made of', &
                       !          h_many_edg(z),'edge(s)'

                   ENDIF
                ENDIF

             ENDIF

          ENDDO
          IF (.NOT.found)THEN        
             z = z + 1
             sum=0

             DO s=1,z-1
                sum = sum + h_many_edg(s)
             ENDDO
             found=.TRUE.
             h_many_edg(z)= k-sum
             !WRITE(*,*)'   ...open & structured boundary',z,'is made of',h_many_edg(z),'edge(s)'
             grid % edges(which_e(k))%cl_open = 1    
             !indica la fine di un lato che non si chiude su se stesso

          ENDIF


       ENDIF
    ENDDO

!    WRITE(*,*)'which_e',which_e
!    WRITE(*,*)'h_many_edg',h_many_edg
    !-------------------------------------------
    ALLOCATE(delaun_edges(m))
    delaun_edges = 0 
    j = 0 
    DO i=1,m
       IF (grid % edges(i) % struct == 0)THEN
          j = j + 1
          delaun_edges(j) = i
          ! identifico in un vettore i  lati non strutturati
       ENDIF
    ENDDO
    !WRITE(*,*) 'delaun_edges', delaun_edges


    !Sezione di controllo dell'uniformità di dati 
    !   dei lati costituenti un lato indipendente  
    ! delta,leng_mod,leng_var,leng_par
    !write(*,*)'test begin'

    ALLOCATE(save_delta(m))

    t = 1
    DO j=1,z

       e = which_e(t)
       ref_delta = grid % edges(e) % delta
       ref_mod = grid % edges(e) % leng_mod
       ref_var = grid % edges(e) % leng_var
       ref_par = grid % edges(e) % leng_par
       ref_cross = grid % edges(e) % cross_ovlp       
       ref_auto = grid % edges(e) % auto_ovlp       

       DO k=t,t+h_many_edg(j)-1
          e = which_e(k)
          IF ((grid % edges(e) % delta    /= ref_delta).OR. &
               (grid % edges(e) % leng_mod /= ref_mod)  .OR. &
               (grid % edges(e) % leng_var /= ref_var)  .OR. &
               (grid % edges(e) % leng_par /= ref_par)) THEN

             WRITE(*,*)
             WRITE(*,*)'ERROR. Configuration file:' 
             WRITE(*,*)'adjacent edges ',which_e(t),'and',which_e(k)
             WRITE(*,*)'must have the same delta/leng_parameters. STOP' 
             WRITE(*,*)
             STOP 
          ENDIF
          IF((grid % edges(e) % cross_ovlp /= ref_cross).OR. &
               (grid % edges(e) % auto_ovlp  /= ref_auto ))THEN

             WRITE(*,*)
             WRITE(*,*)'ERROR. Configuration file:'
             WRITE(*,*)'adjacent edges ',which_e(t),'and',which_e(k)
             WRITE(*,*)'must have the same cross/auto overlapping check parameters. STOP'
             WRITE(*,*)
             STOP
          ENDIF
          IF(k == (t+h_many_edg(j)-1))THEN
             ! completato il controllo sul lato indip j senza pb
             delta(j) = ref_delta
             mod(j)   = ref_mod
             var(j)   = ref_var
             par(j)   = ref_par

             cross(j) = ref_cross
             auto(j)  = ref_auto
          ENDIF
       ENDDO
       t = t + h_many_edg(j)

    ENDDO
    DEALLOCATE(save_delta) 
    !-------------------------------------------  
    DO i=1,nlayer

       ALLOCATE( layer(i) % pnt_list(z))   
       ALLOCATE( layer(i) % num_el(z))
       ALLOCATE( layer(i) % max_lay(z))

       t = 1
       DO j=1,z   !per ogni lato indipendente...

          pnt_list => layer(i)%pnt_list(j)     
          !testa di lista per accumulare i punti dei vari lati in
          !un solo self-closing lato      

          pnt_list = new_point()

          !  layer(i) % max_lay(j) = grid % edges(which_e(t))%struct
          !  write(*,*)'copied lay ij',i,j, layer(i)%max_lay(j)

          IF (i > 1) THEN 
             layer(i) % max_lay(j) = layer(i-1) % max_lay(j)
             !write(*,*)'copied lay ij',i,j, layer(i)%max_lay(j)
             CYCLE
          ELSE
             layer(i) % max_lay(j) = grid % edges(which_e(t))%struct
             !write(*,*)'copied lay ij',i,j, layer(i)%max_lay(j)
          ENDIF

          head => pnt_list
          head % opnd = 0
          e = which_e(t)
          p => grid % edges(e) % add_head % next 

          DO k=t,t+h_many_edg(j)-1

             e = which_e(k)
             !last_height evita che ad un vertice non venga applicata la media delle
             ! altezze ciascuna relativa ad un lato diverso;altrimenti si attribuisce
             ! l'altezza del lato che sta per iniziare        
             !PUNTI INIZIALI         
             add_head => grid % edges(e)%add_head
             IF (k == t)THEN  !lato di inizio
                !*****************
                p => add_head % next  !decide se includere il primo elemento 
                head % x = p % x
                p % ppp => head!################
                head % height = p % height

                IF (p % option /= 0)THEN
                   !trasferimento delle informazioni riguardanti la scia

                   head % option = p % option
                ENDIF

             ELSEIF (k > t) THEN !lato continuativo
                p => add_head
             ENDIF

             s=0

             DO
                !per non ricopiare l'ultimo
                p => p % next
                IF(grid % edges(e)%cl_open == 0)THEN
                   IF(ASSOCIATED(p%next,add_head))THEN
                      last_height = p % height
                      IF(k == t + h_many_edg(j)-1)THEN
                         head % height = (head % height + last_height)/2
                      ENDIF
                      EXIT
                   ENDIF
                ELSE
                   head % opnd = 1 !serve a base_leng 
                   IF(ASSOCIATED(p,add_head)) THEN
                      last_height = 0
                      EXIT
                   ENDIF
                ENDIF
                s=s+1
                q => new_point()
                q % x = p % x
                p % ppp => q!#############
                IF((s == 1).AND.(k > t)) THEN 
                   ! lato continuativo
                   q % height = (p % height + last_height)/2
                ELSE
                   q % height = p % height
                ENDIF
                IF (p % option /= 0)THEN
                   !trasferimento delle informazioni riguardanti la scia
                   !WRITE(*,*)'option',p % option,'transferred',p % x 
                   q % option = p % option
                ENDIF

                CALL insert_next(head,q)

             ENDDO


          ENDDO !sui lati del contorno indipendente
          t=t+h_many_edg(j)

       ENDDO!su ogni contorno indipendente
    ENDDO!su ogni strato

!    WRITE(*,*) '   done.'

  END SUBROUTINE organize_bou





   SUBROUTINE base_leng ( head )
!
! Compute the distance between consecutive points of the same layer/list.
! Point p stores its distance from point p%next.
!
    IMPLICIT NONE
    TYPE(point), POINTER :: head, p, p2, p3

    !calcola le lunghezze tra i punti dello strato
    !attribuisce al punto la MINIMA fra le due adiacenti

    p => head
    DO

       !p1 => p % prev
       p2 => p
       p3 => p % next

       p % leng = SQRT(segment_length(2,p2,p3))
       p => p % next
       
       IF ( ASSOCIATED(p,head) ) THEN
          IF ( head % opnd .eq. 1 )THEN ! procedura per lati aperti
             head % leng = SQRT(segment_length(2,head%next,head))
             head % prev % leng = SQRT(segment_length(2,head%prev,head%prev%prev))
          END IF
          EXIT
       END IF

    ENDDO

  END SUBROUTINE base_leng





  SUBROUTINE calc_leng_max ( mod, head )
!
! Compute the maximum allowed length for the quadrilateral cells; this value, together
! with the height, allows to fix the local maximum cell stretch.
!
    IMPLICIT NONE
!
    INTEGER, INTENT(IN) :: mod
    TYPE(point), POINTER :: head
!
    TYPE(point), POINTER :: p
    REAL(KIND=8) :: lmin, lmax
    REAL(KIND=8) :: sum
    INTEGER :: n
!
    n = 0
    sum = 0
!
    p => head
    lmin = p % leng 
    lmax = p % leng

    DO
       n = n+1                   ! conta gli elementi
       sum = sum + p % leng      ! somma gli elementi 
       lmin = MIN(lmin,p % leng) ! minimo
       lmax = MAX(lmax,p % leng) ! massimo
       p => p % next; IF(ASSOCIATED(p,head))EXIT
    ENDDO
    sum = sum / n

!   LUNGHEZZA LOCALE
    IF ( mod .eq. 1 ) THEN
!
       p => head
       DO
          p % leng_max = p % leng
          p => p % next; IF (ASSOCIATED(p,head)) EXIT
       END DO

!   LUNGHEZZA MEDIATA ( su tre punti )
    ELSE IF ( mod .eq. 2 ) THEN
!
       p => head
       DO
          p % leng_max = ( p % prev % leng + p % leng + p % next % leng )/ 3
          p => p % next; IF (ASSOCIATED(p,head)) EXIT
       END DO
!       
!   LUNGHEZZA GLOBALE ( media sul lato )
    ELSE IF ( mod .eq. 3 ) THEN
!
       p => head
       DO
          p % leng_max = sum
          p => p % next; IF (ASSOCIATED(p,head)) EXIT
       END DO
!
    END IF

  END SUBROUTINE calc_leng_max



  SUBROUTINE update_leng_max ( head, law, par, lay, nlayer )

    IMPLICIT NONE
    TYPE(point),POINTER :: head,p,q
    INTEGER ::  law,lay,nlayer, cnt
    REAL(KIND=8),INTENT(IN) :: par
    REAL(KIND=8),PARAMETER :: pi = 3.1415927464102d0
!
!    +-------------------+
!    |  Update bou_leng  |
!    +-------------------+
!
    p => head
    DO 
       p % bou_leng1 = p % old % bou_leng1
       IF ( (lay .eq. 1) .AND. (p % old % theta .gt. (pi - 60*pi/180)) ) THEN
       
          IF ( ASSOCIATED( p%old%ppp, p%prev%prev ) ) THEN
            p % bou_leng1 = p % old % leng
          END IF
       
       END IF
       p => p % next
       IF (ASSOCIATED(p,head)) EXIT
    END DO

    p => head
!
!    +-------------------------------+
!    |  Constant leng_max = K * L_0  |
!    +-------------------------------+
!
    IF ( law .eq. 1 ) THEN

       DO

          IF ( lay .eq. 1 ) THEN
             p % leng_max = p % old % leng_max * par 
          ELSE
             p % leng_max = p % old % leng_max
          ENDIF

          p => p % next
          IF(ASSOCIATED(p,head))EXIT

       END DO
!
!    +--------------------------------------+
!    |  Progressive leng_max = K * L_(i-1)  |
!    +--------------------------------------+
!
    ELSE IF ( law .eq. 2 ) THEN
!

       cnt = 0
       DO
!
          IF ( (lay .eq. 1) .AND. (ASSOCIATED(p % old % wake)) ) THEN
             
             q => p
             DO
                p % leng_max = p % leng
                IF ( ASSOCIATED(p % next % old, q % old) ) EXIT
                p => p % next
                cnt = cnt + 1
                !write (123,'(1x,i5,2(1x,e15.6),2x,2(1x,e15.6),1x,a1)') cnt, head%x, p%x, 'w'
             ENDDO

          ELSE          
             p % leng_max = p % old % leng_max * par
             cnt = cnt + 1
          END IF
          
          p => p % next
          !write (123,'(1x,i5,2(1x,e15.6),2x,2(1x,e15.6))') cnt, head%x, p%x
          IF ( ASSOCIATED(p, head) ) EXIT

       ENDDO
!
!    +-------------------------------------------------------------+
!    |  Interpolation between layer_0 and assigned leng_max (par)  |
!    +-------------------------------------------------------------+
!
    ELSE IF ( law .eq. 3 ) THEN 

       DO
          p % leng_max = p % bou_leng1 + (par - 1)*p % bou_leng1 * float(lay)/float(nlayer)
          p => p % next
          IF(ASSOCIATED(p,head))EXIT
       ENDDO
!
    ENDIF

  END SUBROUTINE update_leng_max


  
  
  
  SUBROUTINE update_height(head, delta)
  !! Update the thickness of cells to be inserted at each of the points contained into the list.
  !! delta is usually a growth parameter for boundary layer meshes.

    IMPLICIT NONE

    TYPE(point),POINTER :: head,p
    INTEGER :: optn
    REAL(KIND=8)   :: height
    REAL(KIND=8),INTENT(IN)   :: delta

    ! Forced to be 1
    optn = 1

    p => head 
    DO

       IF (ASSOCIATED(p % old)) THEN

          IF (optn == 1) THEN

             height = (1 + delta) * p % old % height
             p % height = height
             p % past_height = p % old % past_height + p % old % height

          ELSEIF (optn == 2)THEN
          
             p % height = ((1+ delta)**2) * p % old % height
             p % leng_max = ((1+ delta)**2) * p % old % leng_max
             WRITE(*,*)'aggiornare past height '
             STOP
          
          ENDIF
          
       ELSE
          
          WRITE(*,*)'warning:old not associated to p:',p % x
          WRITE(*,*)'head',head% x
          p % height = (p % prev % height + p % next % height )/2 

          WRITE(*,*)'height =',p % height 
          WRITE(*,*)'aggiornare past height** '
          STOP
           
       ENDIF

       p => p % next
       
       IF (ASSOCIATED(p, head))  EXIT

    ENDDO

  END SUBROUTINE update_height
  !----------------------------------------------------------------------
  !! Compute the concavity/convexity angle for each point contained into the list
  SUBROUTINE update_angle(head)
    IMPLICIT NONE
    TYPE(point),POINTER :: head,p1,p2,p3
    REAL(KIND=8) ::theta
    p1 => head % prev
    p2 => head 
    p3 => head % next
    DO
       CALL angle(p1,p2,p3,theta)
       p2 % theta = theta

       p1 => p1 % next
       p2 => p2 % next
       p3 => p3 % next

       IF(ASSOCIATED(p2,head))THEN
          EXIT
       ENDIF

    ENDDO

  END SUBROUTINE update_angle
  !-----------------------------------------------------------------------
  !! debug function to suspend the program execution  
  SUBROUTINE wait
    IMPLICIT NONE

    WRITE (*,*)'continue?'
    READ(*,*)

  END SUBROUTINE wait





  SUBROUTINE adjust_bound

    IMPLICIT NONE
    TYPE(point),POINTER :: h1,h2,h3 !teste di lista dei punti
    !h2 centrale strutturata gr%edg(i)%
    !h1 laterale b_gr%edg(i)%
    TYPE(point),POINTER :: sub_h1, sub_h3!, sub_h2
    TYPE(point),POINTER :: c,p1,p3,p_list,p,q,r,head,qh2

    REAL(KIND=8) :: r1,r2,theta,btheta,proj,nproj!,projmed,dn
    REAL(KIND=8),DIMENSION(2) :: l_c1,l_c3,vh1,vh3,displ,vect
    !  REAL(KIND=8),DIMENSION(2,2) :: rot_mat
    REAL(KIND=8),PARAMETER :: pi = 3.1415927464102d0
    real*8 :: mfangle

    INTEGER :: i,k,z,extracted,kmax,store_opt!,sign

    z = 0
    DO i = 1, grid % number_of_edges

       IF ( grid % edges(i) % struct .eq. 1 )THEN
          z = z+1
          p_list => layer( b_grid % number_of_layers + 1 ) % pnt_list(z)

          IF ( p_list % opnd .eq. 1 ) THEN

             h1 => grid % edges(i) % prev % add_head
             h2 => grid % edges(i) % add_head
             h3 => grid % edges(i) % next % add_head 

             p => h2
             DO
                p => p % next
                !       write(*,*)'h2 ', p % x, p % index, p % cnd, p % option
                IF(ASSOCIATED(p,h2))EXIT
             ENDDO
             ! liste dei contorni 1D  
             sub_h1 => grid % edges(i) % prev % old % add_head 
             sub_h3 => grid % edges(i) % next % old % add_head
             !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             ! determino lo strato di appartenenza della (eventuale) virt_head
             kmax = b_grid % number_of_layers
             p => p_list
             DO
                !        write(*,*)'p_list', p % x,'with mark',p % mark
                IF((p % mark == 0).OR.(ASSOCIATED(p, h2))) EXIT
                p => p % old; kmax = kmax - 1
             ENDDO
             store_opt = p % index
             p_list => p
             ! qh2 e' la testa di lista dell'ultimo strato di punti
             ! puo' anche essere una virt_head, qindi devo portarmi su di essa
             qh2 => p_list

             !      WRITE(*,*)' edge',i,':'
             !punto di frontiera tra i due tipi di lato
             c  => h1 % prev
             !penultimo punto del lato precedente non strutturato 
             p1 => h1 % prev % prev
             !secondo punto del lato strutturato
             p3 => h2 % next % next 

             !       write(*,*)'c' , c%x
             !       write(*,*)'p1', p1%x
             !       write(*,*)'p3', p3%x

             r1   = SQRT(segment_length(2,c,p1))
             r2   = SQRT(segment_length(2,c,p3))
             l_c1 = (p1 % x - c % x)/r1 !normalized vectors
             l_c3 = (p3 % x - c % x)/r2

             CALL angle ( p1, c, p3, theta )
             CALL angle ( p1, c, layer(1) % pnt_list(z) % next, btheta )
             !      write(*,*)'theta, btheta',theta*180/pi, btheta * 180/pi

             IF((btheta < 0).AND.(theta > 0)) theta = -theta

             !------------------------------------------------------- inizio lato    
             IF((btheta*180/pi > -170).AND.(btheta*180/pi < -65))THEN

                !WRITE(*,*)'       beginning side:  points moved on adjacent Delaunay edge'
                grid % edges(i) % adj_begin = 1

                head => qh2 ! lista strutturata
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! determino l'inclinazione della  retta sulla quale proietto i punti 
                IF((btheta*180/pi > -170).AND.(btheta*180/pi < -90))THEN
                   ! lato sovrapposto alla spalla,rotazione antioraria di l_c3
                   mfangle = pi/2.d0 
                   CALL rot( mfangle,l_c3,displ)
                   vect = displ * head % past_height * TAN(theta)
                   IF(SQRT(DOT_PRODUCT(vect,vect)) > head % leng)THEN
                      WRITE(*,*)'****la sovrapposizione sull''inizio del lato aperto'
                      WRITE(*,*)'    supera la seconda colonna di punti!STOP'  
                      STOP
                   ENDIF
                ELSEIF((btheta*180/pi >= -90).AND.(btheta*180/pi < -65))THEN
                   ! lato non sovrapposto,rotazione oraria di l_c3
                   mfangle = -pi/2.d0 
                   CALL rot( mfangle,l_c3,displ)

                ENDIF
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                ! cerco il punto sul lato delaunay che sta appena sopra la proiezione
                ! della testa aperta&strutturata sul lato delaunay stesso
                p => h1
                q => sub_h1
                nproj = 0
                proj = r2
                k = 0
                DO; p => p % prev; q => q % prev;IF(ASSOCIATED(p,h1))EXIT !all'indietro..
                   !proseguo sul contorno adiacente
                   k = k + 1
                   vh1 = p % x - c % x
                   proj = MAX(nproj,proj)
                   nproj = SQRT(vh1(1)**2 + vh1(2)**2)*ABS(COS(theta))
                   !         write(*,*)'  p', p % x
                   !         write(*,*)'  proj,nproj',proj, nproj
                   !proiezione sulla spalla del contorno adiacente

                   IF ((nproj > proj).AND.(k > 2))THEN

                      r => new_point()
                      r %  x = qh2 % x + displ * qh2 % past_height * TAN(theta)
                      !write(*,*)'*'
                      IF(SQRT(segment_length(2,r,p)) < qh2 % height/3 )THEN
                         !            WRITE(*,*)'  ROSICCHIO ANCORA UN PO'''
                         p => p % prev
                      ENDIF
                      DEALLOCATE(r)
                      EXIT
                   ENDIF
                ENDDO  !proiezione sulla spalla
                
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! tutti i punti che stanno 'sotto' la proiezione appena trovata vengono
                ! eliminati dalla lista dei punti del lato delaunay, sia 2D(p) che 1D(q)
                ! attraverso la funzione extract_point
                ! i nuovi punti proiettati dal blocco strutturato vengono inseriti solo
                ! nella lista delaunay 1D; quella 2D rimane mutilata e viene passata
                ! come contorno al triangolatore.

                ! inverto il senso di percorrenza
                p => p % next; q => q % next
                extracted = 0
                DO 
                   IF(ASSOCIATED(p,h1 % prev))THEN
                      EXIT !WRITE(*,*)'no points extracted from Delaunay edge'
                   ENDIF
                   ! per estrarre i punti che rientrano nella proiezione della spalla sul
                   ! contorno adiacente
                   !         write(*,*)'deleting p', p % x
                   CALL extract_point(p); CALL extract_point(q)

                   extracted = extracted + 1

                   p => p % next; q => q % next

                   IF(ASSOCIATED(p,h1 % prev)) EXIT
                ENDDO
                !       WRITE(*,*)' ',extracted,'points deleted from Delaunay edge'
                !estraggo fino al penultimo;cambiero' coordinate all'ultimo.
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! elimino i punti da proiettare dalla lista delaunay che viene passata 
                ! al triangolatore -> senza i fianchi del blocco strutturato
                !p => h2 % next
                !DO k= 2, kmax+1
                !  p => p % next
                !  write(*,*)'...extracting',k,' p', p % x,p % cnd,p % option
                !  CALL extract_point(p)
                !ENDDO
                p => h2 
                DO k= 1, kmax
                   p => p % next
                   !         write(*,*)'...extracting',k,' p', p % x,p % cnd,p % option
                   CALL extract_point(p)
                ENDDO
                h2 % next % cnd = -i
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! proietto i punti del blocco strutturato sul contorno delaunay
                ! a partire dallo strato piu' avanzato
                DO k =1,kmax !b_grid%number_of_layers !cambiamento coordinate
                   !         write(*,*)'8  ', head % x, displ, head % past_height, theta 
                   head %  x = head % x + displ * head % past_height * TAN(theta)
                   !         write(*,*)'8' , head % x, head % past_height

                   IF(ASSOCIATED(head,qh2))THEN
                      h1 % prev % x = head % x
                      h2 % next % x = head % x
                      !         h2 % next % option = store_opt
                      !write(*,*)'store_opt',store_opt
                      !write(*,*)' last not-extracted is moved to',head%x,h2%next%cnd,h2%next%option
                      !cambiamento di coordinate dell'ultimo non estratto
                   ENDIF
                   !introduzione dei punti proiettati solo nella lista delaunay 1D 
                   q => new_point() 
                   CALL pushbot_point(sub_h1 % prev,q)
                   q % x = head % x
                   q % ppp => head!????
                   !      write(*,*)'9' 

                   head => head % old ! passaggio al layer successivo
                   sub_h1 % prev % cnd = 0
                   !      write(*,*)'10' 
                ENDDO

                NULLIFY(head % next % old)
                !      write(*,*)'11' 
             ELSE
                ! no  edge modifying
                ! -------------------------------------ADD SIDE POINTS TO THE LAYER
                WRITE(*,*)'       beginning side:  points added to b.layer front', layer(1)%max_lay(z)
                grid % edges(i) % adj_begin = 0
                h2 % next % cnd = i

                head => qh2
                p => head !riemann type

                !DO k=1,b_grid%number_of_layers
                !   p => p % old
                !   q => new_point()
                !   q % x = p % x

                !  WRITE(*,*)'p',p % x,p % index, k
                !  q % cnd = i
                !  q % option = p % index
                !in progressione dall'alto verso il contorno solido
                !  CALL pushtop_point(h2,q) !fifo
                !ENDDO
                h2 % next % cnd = -i

                !NULLIFY(head % next % old)

             ENDIF
             !#############################################################################    
             !------------------------------------------------------- fine lato    
             p_list => layer(b_grid%number_of_layers+1)% pnt_list(z)
             !      write(*,*)'12' 
             c => h3 % next                           
             p1 => h2 % prev % prev 
             p3 => h3 % next % next  

             !       write(*,*)'c' , c%x
             !       write(*,*)'p1', p1%x
             !       write(*,*)'p3', p3%x


             r1   = SQRT(segment_length(2,c,p1))
             r2   = SQRT(segment_length(2,c,p3))
             l_c1 = (p1 % x - c % x)/r1 !normalized vectors
             l_c3 = (p3 % x - c % x)/r2

             CALL angle(p1, c, p3, theta)
             CALL angle(layer(1) % pnt_list(z) % prev % prev, c, p3, btheta)
             !      write(*,*)'theta, btheta',theta*180/pi, btheta * 180/pi

             IF((btheta < 0).AND.(theta > 0)) theta = -theta
             !per eliminare le incertezze di segno quando l'angolo è prossimo a 180

             IF((btheta*180/pi > -170).AND.(btheta*180/pi < -65))THEN    
                !WRITE(*,*)'  sign<0,SIDE POINTS DISPLACEMENT'
                !       WRITE(*,*)'       final side:      points moved on adjacent Delaunay edge'
                grid % edges(i) % adj_end = 1

                !        write(*,*)'qh2',qh2 % x, qh2 % mark
                !        write(*,*)'qh2 next',qh2 % next % x, qh2 %next% mark
                !        write(*,*)'qh2 prev',qh2 % prev%x, qh2 %prev% mark
                !        write(*,*)'qh2 jump_f',qh2 % jump_f % x, qh2 %jump_f% mark
                !        write(*,*)'qh2 jump_b',qh2 % jump_b%x, qh2 %jump_b% mark
                !        write(*,*)'plist',p_list % x, p_list % mark

                kmax = b_grid % number_of_layers
                p => p_list % prev
                DO
                   !          write(*,*)'p_list end', p % x,'with mark',p % mark
                   IF((p % mark == 0).OR.(ASSOCIATED(p, h2 % prev))) EXIT
                   p => p % old
                   kmax = kmax - 1
                ENDDO
                store_opt = p % index
                qh2 => p
                head => qh2
                !        write(*,*)'end head', head % x, head % mark

                IF((btheta*180/pi > -170).AND.(btheta*180/pi < -90))THEN
                   ! lato sovrapposto alla spalla,rotazione oraria di l_c1
                   mfangle = -pi/2.d0 
                   CALL rot( mfangle, l_c1, displ)
                   vect = displ * head % past_height * TAN(theta)
                   IF(SQRT(DOT_PRODUCT(vect,vect)) > head % leng)THEN
                      WRITE(*,*)'****la sovrapposizione sulla fine del lato aperto'
                      WRITE(*,*)'    supera la penultima colonna di punti!STOP'  
                      STOP
                   ENDIF
                ELSEIF((btheta*180/pi >= -90).AND.(btheta*180/pi < -65))THEN
                   ! lato non sovrapposto,rotazione antioraria di l_c1
                   mfangle = pi/2.d0 
                   CALL rot( mfangle, l_c1, displ)
                ENDIF

                p => h3
                q => sub_h3
                nproj = 0
                proj = r1
                k = 0
                !WRITE(*,*)'be',p % x,r1
                DO;p => p % next;q => q % next;IF(ASSOCIATED(p,h3))EXIT !in avanti...
                   k = k + 1
                   vh3 = p % x - c % x
                   proj = MAX(proj,nproj)
                   ! proj   altezza della spalla(pti inseriti in direzione normale!!!!!)
                   nproj = SQRT(vh3(1)**2 + vh3(2)**2)*ABS(COS(theta))
                   ! nproj  proiezione sulla spalla del contorno adiacente

                   IF ((nproj > proj).AND.(k > 3))THEN
                      !write(*,*)'qh2 h' ,qh2 %prev% height /3
                      r => new_point()
                      r %  x = qh2 %prev% x + displ*qh2%prev%past_height*TAN(theta)
                      !write(*,*)'segm',SQRT(segment_length(2,r,p % prev))
                      IF(SQRT(segment_length(2,r,p)) < qh2 %prev% height/3 )THEN
                         !            WRITE(*,*)'  ROSICCHIO ANCORA UN PO'''
                         p => p % next
                      ENDIF
                      DEALLOCATE(r)
                      EXIT
                      EXIT

                   ENDIF
                ENDDO  !proiezione sulla spalla

                p => p % prev 
                q => q % prev
                extracted = 0
                DO 
                   IF(ASSOCIATED(p,h3 % next))THEN
                      !WRITE(*,*)'no points extracted from Delaunay edge'
                      EXIT ! all'indietro...
                   ENDIF

                   !         write(*,*)'deleting p', p % x
                   CALL extract_point(p)
                   CALL extract_point(q)
                   extracted = extracted + 1 
                   p => p % prev
                   q => q % prev
                   IF(ASSOCIATED(p,h3 % next))EXIT
                ENDDO
                !        WRITE(*,*)' ',extracted,'points deleted from Delaunay edge'

                !p => h2 % prev
                !DO k= 2, kmax+1
                !  p => p % prev
                !  write(*,*)'...extracting',k,' p', p % x, p % cnd,p % option
                !  CALL extract_point(p)
                !ENDDO
                p => h2 
                DO k= 1, kmax
                   p => p % prev
                   !         write(*,*)'...extracting',k,' p', p % x,p % cnd,p % option
                   CALL extract_point(p)
                ENDDO
                h2 % prev % cnd = -i


                DO k =1,kmax!b_grid%number_of_layers !cambiamento coordinate
                   !       write(*,*)'1   ', head % x, displ, head % past_height, theta 
                   head % x = head % x + displ*head%past_height*TAN(theta)
                   !       write(*,*)'2', head % x

                   IF(ASSOCIATED(head,qh2))THEN!%prev))THEN
                      h3 % next % x = head % x
                      h2 % prev % x = head % x
                      !          h2 % prev % option = store_opt
                      !write(*,*)'store_opt',store_opt
                      !write(*,*)' last not-extracted is moved to',head%x,h2%prev%cnd,h2%prev%option
                   ENDIF
                   !       write(*,*)'3' 
                   q => new_point()
                   CALL pushtop_point(sub_h3 % next,q)
                   q % x = head % x
                   q % ppp => head

                   head => head % old ! passaggio al layer successivo
                   !       write(*,*)'4' 
                ENDDO
                sub_h3 % next % ppp => head
                sub_h3 % next % cnd = 0
                !       write(*,*)'5' 

                NULLIFY(head % next % old)

             ELSE
                WRITE(*,*)'       final side:      points added to b.layer front', layer(1)%max_lay(z)
                ! no  edge modifying
                ! add side points to the layer
                grid % edges(i) % adj_end = 0
                h2 % prev % cnd = i

                head => qh2
                p => head % prev !riemann type

                !DO k=1,b_grid%number_of_layers
                !   p => p % old
                !   q => new_point()
                !   q % x = p % x
                !   WRITE(*,*)'p',p % x,p % index,  k
                !   q % cnd = i
                !   q % option = p % index
                !  !in progressione dall'alto verso il contorno solido
                !   CALL pushbot_point(h2,q) !fifo
                !ENDDO
                h2 % next % cnd = -i

                !NULLIFY(head % next % old)
             ENDIF

          ENDIF

       ENDIF

    ENDDO

  END SUBROUTINE adjust_bound





END MODULE schemi
