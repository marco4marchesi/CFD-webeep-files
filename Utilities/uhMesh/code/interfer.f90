MODULE interfer

  USE grid_types
  USE eucl_delaunay

  IMPLICIT NONE

  CONTAINS


  SUBROUTINE check_interference(bi_box, gath_final_res)
  !-----------------------------------------------------------------------------------
  IMPLICIT NONE

  TYPE(box), POINTER, DIMENSION(:) :: bi_box
  INTEGER,            INTENT(OUT)  :: gath_final_res

  TYPE(box), POINTER :: box1, box2
  INTEGER :: size_bi_box, k, kmax, i, j, i_min, &
             res, num1, num2, min_points, final_res
  !-----------------------------------------------------------------------------------

    size_bi_box = SIZE(bi_box) 
    kmax = (size_bi_box * (size_bi_box-1) / 2) 

    k = 0
    gath_final_res = 0
    ! self_check indica l'attivazione del controllo di auto-interferenza
    ! che viene effettuato sugli small_block e quindi cio' deve essere noto a
    ! livello dei big_block, che devono in questo caso essere x forza divisi 

    DO i = 1, size_bi_box

       box1 => bi_box(i)
       IF (box1 % cross_ovlp == 0) CYCLE

       IF (box1 % auto_ovlp == 1) THEN
          i_min = i 
       ELSE
          i_min = i + 1 
       ENDIF

       DO j = i_min, size_bi_box

          box2 => bi_box(j)
          
          IF (box2 % cross_ovlp == 0) CYCLE
          IF ((box1 % structured == 0) .AND. (box2 % structured == 0)) CYCLE

          CALL box_compare2(box1, box2, res)

          IF (res == 1) THEN

             ! dopo aver controllato la sovrapposizione, i big_box interessati 
             ! sono divisi in parti piu' piccole, se questo non e' gia' stato fatto per
             ! interferenze con altri blocchi(active = 1)

             min_points = 3
             ! numero di punti minimo utilizzato per la suddivisione dei box

             IF (box1 % active == 1) THEN
                num1 = SIZE(box1 % small_box)
             ELSE
                CALL big2small_box(box1, min_points, num1)
                box1 % active = 1
             ENDIF

             IF (box2 % active == 1) THEN
                num2 = SIZE(box2 % small_box)
             ELSE
                CALL big2small_box(box2, min_points, num2)
                box2 % active = 1
             ENDIF

             CALL small_box_interference(box1, box2, num1, num2, final_res)

             IF (final_res > 0) THEN
                gath_final_res = gath_final_res + 1
                WRITE(*,*) '   WARNING: Interference detected. Grid generation may fail...'
             ENDIF

          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE check_interference
  
  
  
  
  
  !-----------------------------------------------------------------------------
  SUBROUTINE small_box_interference(box1, box2, num1, num2, gath_outcome)

    IMPLICIT NONE

    TYPE(box), POINTER    :: box1, box2, sm_box1, sm_box2
    INTEGER, INTENT(IN)   :: num1, num2
    INTEGER, INTENT(OUT)  :: gath_outcome
    INTEGER               :: res, imin, num_box, i, j, imax, tnum1, tnum2, &
         outcome, check

    ! subroutine dedicata al controllo dell'interferenza tra small_boxes, viene
    ! eseguito, attivato da check, il controllo di self_interference

    gath_outcome = 0
    num_box = num1 + num2

    imax = num_box

    DO i=1, num1
       !------------------ 
       IF(ASSOCIATED(box1,box2))THEN ! auto interference
          imin = num1 + i + 1
       ELSE                          ! cross interference
          imin = num1 + 1
       ENDIF
       !-----------------
       sm_box1 => box1 % small_box(i)
       DO j = imin, imax

          IF((ASSOCIATED(box1,box2)).AND.(j-imin+1 <= i+1)) CYCLE

          sm_box2 => box2 % small_box(j-imin+1)

          IF(check == 1)THEN
             res = 1
          ELSE     
             CALL box_compare2(sm_box1, sm_box2, res)
          ENDIF

          IF(res == 1) THEN

             IF(sm_box1 % active == 1) THEN
                tnum1 = SIZE(sm_box1 % small_box)
             ELSE
                CALL small2tiny_box(sm_box1, tnum1)
                sm_box1 % active = 1
             ENDIF

             IF(sm_box2 % active == 1) THEN
                tnum2 = SIZE(sm_box2 % small_box)
             ELSE
                CALL small2tiny_box(sm_box2, tnum2)
                sm_box2 % active = 1
             ENDIF

             CALL tiny_box_interference(sm_box1, sm_box2, tnum1, tnum2, check, outcome)

             IF(outcome > 0)THEN
                gath_outcome = gath_outcome + 1
             ENDIF
          ENDIF
       ENDDO

    ENDDO

  END SUBROUTINE small_box_interference
  !-----------------------------------------------------------------------------
  SUBROUTINE tiny_box_interference(sm_box1,sm_box2,num1,num2,check,gath_outcome)

    IMPLICIT NONE

    TYPE(box), POINTER    :: sm_box1, sm_box2, tiny_box1, tiny_box2
    INTEGER, INTENT(IN)   :: check, num1, num2
    INTEGER, INTENT(OUT)  :: gath_outcome
    INTEGER               :: res, imin, num_box, i, j, imax, outcome

    ! subroutine dedicata al controllo dell'interferenza tra small_boxes, viene
    ! eseguito, attivato da check, il controllo di self_interference

    gath_outcome = 0
    num_box = num1 + num2
    IF(check == 1)THEN
       imax = num1
    ELSE
       imax = num1
    ENDIF

    DO i=1, imax

       IF(check == 1)THEN
          imin = num1 + 1
       ELSE
          imin = num1 + 1
       ENDIF

       tiny_box1 => sm_box1 % small_box(i)
       DO j=imin, num_box
          IF(ASSOCIATED(sm_box1 % c2, sm_box2 % c1).AND.(i==num1).AND.(j==imin))CYCLE

          tiny_box2 => sm_box2 % small_box(j-imin+1)
          CALL box_compare2(tiny_box1, tiny_box2, res)

          IF(res == 1) THEN
             CALL solve_tinybox_interference(tiny_box1, tiny_box2, outcome)
             IF(outcome == 1)THEN
                gath_outcome = gath_outcome + 1
             ENDIF
          ENDIF
       ENDDO

    ENDDO

  END SUBROUTINE tiny_box_interference
  !-----------------------------------------------------------------------------
  SUBROUTINE solve_tinybox_interference(box1, box2, out)

    IMPLICIT NONE
    TYPE(box), POINTER   :: box1, box2, segm_box, point_box
    TYPE(point), POINTER :: p, q
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: trig_qnt
    REAL(KIND=8)         ::  x1,y1, y2, xp, yp, delta,&
         h_lim, h_real, ymin, ymax, leng
    INTEGER, INTENT(OUT) :: out
    INTEGER              :: num_points,i, mark

    ! subroutine per il controllo definitivo della possibile interferenza
    ! individuata a livello di small_box; 

    ! e' necessario conoscere gli angoli di inclinazione assoluta dei segmenti di
    ! ciascun box; per i lati strutturati basterebbe calcolare il primo ed
    ! aggiungere per gli altri i relativi theta gia' calcolati; per i lati delaunay
    ! bisogna invece calcolare tutto da zero.
    ! SBAGLIATO! se ci sono due lati strutturati ho gia' tutto, altrimenti uso il
    ! lato delaunay come lato dei punti, cosi' non devo calcolare niente di piu'.

    out = 0
    segm_box => box1
    point_box => box2
    IF((box1 % structured+box2 % structured)== 1)THEN
       IF(point_box % structured == 1)THEN
          !    write(*,*)'exchange, box1 = point_box'
          point_box => box1
          segm_box => box2
       ENDIF
    ENDIF

    ! ipotizzo che iniziare da un box o da un altro sia indifferente, quanto ad
    ! onere di calcolo (sara' vero?) 

    ! ciclo sui punti del point_box
    !calcolo il vettore dei seni e coseni dei segmenti
    mark = 0
    q => segm_box % c1

    num_points = segm_box % num_points
    ALLOCATE(trig_qnt((num_points-1),2))

    DO i=1,num_points-1

       leng = SQRT(segment_length(2,q,q%jump_f))
       trig_qnt(i,1) = (q % jump_f % x(2) - q % x(2))/ leng ! seno
       trig_qnt(i,2) = (q % jump_f % x(1) - q % x(1))/ leng ! coseno

       q => q % jump_f 
    ENDDO

    p => point_box % c1
    DO
       ! ciclo sui segmenti del segment_box
       q => segm_box % c1   
       xp = p % x(1)
       yp = p % x(2)

       DO i=1,num_points-1
          ! verifico la posizione di P rispetto al segmento rappresentato come retta
          ! ax + by + c = 0 e di estremi P1 e P2
          ! (x-x1)sinA-(y-y1)cosA-h_lim=delta
          ! ovvero   h_real-h_lim=delta
          ! con delta<=0 il punto e' a distanza normale dal segmento <= h_lim, 
          ! l'altezza di inserimento del nuovo elemento(si prende una media pesata 
          ! con quella del punto in questione, perche' deve rimanere lo spazio per 
          ! due elementi piu' qualcos'altro)
          x1 = q % x(1)
          y1 = q % x(2)

          h_real = (xp - x1)*trig_qnt(i,1) - (yp - y1)*trig_qnt(i,2)
          !                        sinA                            cosA

          h_lim = (MAX(q % height, q % jump_f % height) + p % height)*3/2

          delta = h_real - h_lim

          ! verifico la 'prospicienza' del punto rispetto al segmento
          IF(delta <= 0)THEN
             ! il punto oltrepassa la retta identificata dal segmento, ma e' necessario
             ! verificare che il punto non sia comunque troppo distante (lateralmente)
             ! dal segmento stesso

             y2 = q % jump_f % x(2)

             ymax = y2 - h_real*trig_qnt(i,2)+h_lim*trig_qnt(i,1)
             ymin = y1 - h_real*trig_qnt(i,2)-h_lim*trig_qnt(i,1)
             !                  cosA                sinA

             delta = MAX(ymin, ymax)
             ymin = MIN(ymin, ymax)
             ymax = delta

             ! sono necessarie delle correzioni (estensioni)
             IF((yp <= ymax).AND.(yp >= ymin))THEN
                ! il punto P sta proprio di fronte al segmento.
                ! e' indifferente confrontare le y o le x, siamo su di una retta.
                out = 1 
                !        write(*,*)'              PERDINCIBACCO!'
                p % mark = 1
                q % mark = 1

             ENDIF

          ENDIF
          q => q % jump_f
       ENDDO
       p => p % jump_f
       IF(ASSOCIATED(p, point_box % c2))EXIT
    ENDDO

  END SUBROUTINE solve_tinybox_interference
  !-----------------------------------------------------------------------------
  SUBROUTINE box_compare2(box1, box2, res)

    IMPLICIT NONE

    TYPE(box),  POINTER      :: box1, box2
    INTEGER,    INTENT(OUT)  :: res
    REAL(KIND=8)  :: c1x, c1y, c2x, c2y, cxmin, cxmax, cymin, cymax, &
         p1x, p1y, p2x, p2y, pxmin, pxmax, pymin, pymax, &
         dcx, dcy, dpx, dpy, d1, d2

    ! controllo se le coordinate del secondo box sono comprese negli intervalli di
    ! x e y definiti dal primo box
    res = 0
    IF((box1 % c2 % x(1) == box2 % c1 % x(1)).AND. &
         (box1 % c2 % x(2) == box2 % c1 % x(2))) THEN
       !    write(*,*)'blocchi CONSECUTIVI '
       !RETURN
    ELSE
       ! controlla che i due box non siano consecutivi sullo stesso contorno,
       ! altrimenti ammettendo la self-interference si avrebbero sempre 
       ! controlli positivi...
       c1x = box1 % corner1(1)
       c1y = box1 % corner1(2)
       c2x = box1 % corner2(1)
       c2y = box1 % corner2(2)

       cxmin = MIN(c1x, c2x)  ;  cymin = MIN(c1y, c2y)
       cxmax = MAX(c1x, c2x)  ;  cymax = MAX(c1y, c2y)

       p1x = box2 % corner1(1)  
       p1y = box2 % corner1(2)  
       p2x = box2 % corner2(1) 
       p2y = box2 % corner2(2)

       pxmin = MIN(p1x, p2x)  ;  pymin = MIN(p1y, p2y)
       pxmax = MAX(p1x, p2x)  ;  pymax = MAX(p1y, p2y)

       dcx = cxmax - cxmin ! positiva per forza
       dpx = pxmax - pxmin
       d1 = pxmin - cxmin

       res = 0
       IF(d1 > 0)THEN
          IF(d1 < dcx)THEN
             res = 1
             ! pxmin si trova tra cxmin e cxmax 
          ENDIF
       ELSE
          IF(dpx > ABS(d1))THEN
             res = 1
             ! pxmax supera cxmin
          ENDIF
       ENDIF

       IF(res == 1)THEN ! si e' verificata un'alternanza delle xmin, xmax
          dcy = cymax - cymin
          dpy = pymax - pymin
          d2 = pymin - cymin
          IF(d2 > 0)THEN
             IF(d2 < dcy)THEN
                res = 2
                ! pxmin si trova tra cxmin e cxmax 
             ENDIF
          ELSE
             IF(dpy > ABS(d2))THEN
                res = 2
                ! pxmax supera cxmin
             ENDIF
          ENDIF
       ENDIF

       IF(res == 2)THEN
          res = 1
          !write(*,*)'              overlap!!!!'
       ELSE
          res = 0
       ENDIF
    ENDIF
  END SUBROUTINE box_compare2
  !------------------------------------------------------------------------------
  SUBROUTINE box_compare(box1, box2, res)
    IMPLICIT NONE
    TYPE(box), POINTER :: box1, box2
    INTEGER, INTENT(OUT) :: res
    ! INTEGER         :: xcross, ycross
    REAL(KIND=8)  :: c1x, c1y, c2x, c2y, xmin, xmax, ymin, ymax, &
         p1x, p1y, p2x, p2y

    ! xcross = 0
    ! ycross = 0

    ! controllo se le coordinate del secondo box sono comprese negli intervalli di
    ! x e y definiti dal primo box

    c1x = box1 % corner1(1)
    c1y = box1 % corner1(2)
    c2x = box1 % corner2(1)
    c2y = box1 % corner2(2)

    xmin = MIN(c1x, c2x)  ;  ymin = MIN(c1y, c2y)
    xmax = MAX(c1x, c2x)  ;  ymax = MAX(c1y, c2y)

    p1x = box2 % corner1(1)  
    p1y = box2 % corner1(2)  
    p2x = box2 % corner2(1) 
    p2y = box2 % corner2(2)
    ! IF((box2%corner1(1) > c1x).AND.((box2%corner1(1) < c2x)).AND.&
    !    (box2%corner1(2) > c1y).AND.((box2%corner1(2) < c2y)))THEN
    ! IF(((box2%corner(1) > MIN(c1x, c2x)).AND.(box2%corner(1) < MAX(c1x, c2x)))&
    !.OR.((box2%corner(2) > MIN(c1y, c2y)).AND.(box2%corner(2) < MAX(c1y, c2y))))&
    IF(((( p1x < xmax).AND.(p1x > xmin)).AND.(( p1y < ymax).AND.(p1y > ymin)))&
         .OR.&
         ((( p2x < xmax).AND.(p2x > xmin)).AND.(( p2y < ymax).AND.(p2y > ymin))))&
         THEN
       res = 1 
       write(*,*)'     yes1'
       write(*,*)'xmin', xmin,'xmax',xmax
       write(*,*)'ymin', ymin,'ymax',ymax
       write(*,*)'p1x', p1x,'p1y',p1y
       write(*,*)'p2x', p2x,'p2y',p2y
    ELSE 
       write(*,*)'      no1'
       write(*,*)'xmin', xmin,'xmax',xmax
       write(*,*)'ymin', ymin,'ymax',ymax
       write(*,*)'p1x', p1x,'p1y',p1y
       write(*,*)'p2x', p2x,'p2y',p2y  
       res = 0 
    ENDIF
    IF(res == 0)THEN 
       xmin = MIN(p1x, p2x)  ;  ymin = MIN(p1y, p2y)
       xmax = MAX(p1x, p2x)  ;  ymax = MAX(p1y, p2y)

       IF(((( c1x < xmax).AND.(c1x > xmin)).AND.(( c1y < ymax).AND.(c1y > ymin)))&
            .OR.&
            ((( c2x < xmax).AND.(c2x > xmin)).AND.(( c2y < ymax).AND.(c2y > ymin))))&
            THEN
          res = 1 
          write(*,*)'     yes2'
          write(*,*)'xmin', xmin,'xmax',xmax
          write(*,*)'ymin', ymin,'ymax',ymax
          write(*,*)'c1x', c1x,'c1y',c1y
          write(*,*)'c2x', c2x,'c2y',c2y
       ELSE 
          write(*,*)'      no2'
          write(*,*)'xmin', xmin,'xmax',xmax
          write(*,*)'ymin', ymin,'ymax',ymax
          write(*,*)'c1x', c1x,'c1y',c1y
          write(*,*)'c2x', c2x,'c2y',c2y
          res = 0 
       ENDIF
    ENDIF
  END SUBROUTINE box_compare
  !-----------------------------------------------------------------------------
  SUBROUTINE big2small_box(big, min_points, num_small)
    IMPLICIT NONE
    TYPE(box), POINTER :: big, small
    TYPE(point), POINTER :: b_head, tail, p
    INTEGER, INTENT(IN)  ::  min_points
    INTEGER, INTENT(OUT) :: num_small
    INTEGER              :: mark, pts, elems, i, j, level, enlarge, opnd

    ! suddivisione di un box contenente un lato intero in box piu' piccoli
    ! min_points e' un parametro in ingresso, eventualmente ottimizzabile 

    b_head=> big % head 
    ! calcolo il numero k di punti del box, tenendo conto della forma diversa
    ! delle liste tra lati strutturati e delaunay

    level = big % level
    ! il level indica il grado di suddivisione di un box 
    ! 0 => box iniziale di un intero lato; 1 => box iniziale con 1 suddivisione
    ! viene aggiornato nella subroutine big2small
    ! unica eccezione il lato continente che puo' arrivare fino a 2

    opnd = 0
    IF ( big % structured == 1) THEN
       IF(level == 1)THEN                  ! chiamata su un continent
          tail => big % c2
       ELSEIF(level == 0)THEN
          tail => b_head
          IF(b_head % opnd == 1)THEN        ! lato strutturato aperto
             opnd = 1
             tail => tail % jump_b
          ENDIF
       ENDIF
       mark = 1
       !----------------------------------
    ELSEIF(big % structured == 0)THEN
       IF(level == 1)THEN ! bigtosmall chiamata per suddividere il continent
          tail => big % c2    
       ELSEIF(level == 0)THEN
          tail => b_head % jump_b ! la testa dei lati delaunay e' vuota
          b_head => b_head % jump_f ! la testa dei lati delaunay e' vuota!
          IF((b_head % x(1) == tail % x(1)).AND.(b_head % x(2) == tail % x(2)))THEN 
             opnd = 0
          ENDIF
       ENDIF
       mark = 0
    ENDIF

    CALL enumerate(b_head, tail, pts)
    IF ( big % structured == 1) THEN
       IF(level == 1)THEN ! chiamata su un continent
          elems = pts - 1

       ELSEIF(level == 0)THEN
          IF(opnd == 0) pts = pts - 1
          elems = pts
          IF(opnd == 1) elems = elems - 1

       ENDIF
    ELSEIF(big % structured == 0)THEN
       IF(level == 1)THEN ! chiamata su un continent
          elems = pts - 1

       ELSEIF(level == 0)THEN
          IF(opnd == 0) pts = pts - 1
          elems = pts
          IF(opnd == 1) elems = elems - 1   

       ENDIF
    ENDIF

    num_small =INT(elems/(min_points-1)) 
    !----------------------------------

    IF(num_small < 1) THEN
       write(*,*) 'low number of points 4 box splitting. stop(big2small_box)'
       STOP
    ENDIF


    ALLOCATE(big % small_box(num_small))

    p => b_head

    DO i=1, num_small
       small => big % small_box(i)
       small % structured = mark
       small % head => p
       small % c1 => p  ! e' la "testa di lista" dello small_box, il primo corner
       small % level = small % level + 1
       small % active = 0
       small % cross_ovlp = big % cross_ovlp
       small % auto_ovlp  = big % auto_ovlp
       j = 0

       DO
          j = j + 1 
          IF(j == min_points) THEN
             IF(i == num_small)THEN
                small % c2 => tail

                IF(elems == pts)THEN
                   small % num_points = pts - (num_small-1)*(min_points-1) + 1
                ELSE
                   small % num_points = pts - (num_small-1)*(min_points-1)
                ENDIF

             ELSE
                small % c2 => p
                small % num_points = min_points

             ENDIF

             IF(mark == 1)THEN
                enlarge = 1
                CALL define_struct_box(small % c1, small % c2, small, enlarge)
             ELSEIF(mark == 0)THEN
                CALL define_del_box(small % c1, small % c2, small)
             ENDIF

             EXIT

          ENDIF
          IF(j == pts)EXIT
          p => p % jump_f
       ENDDO
    ENDDO
  END SUBROUTINE big2small_box
  !-----------------------------------------------------------------------------
  SUBROUTINE enumerate (head, tail, points)

    IMPLICIT NONE

    TYPE(point), POINTER :: head, tail, p

    INTEGER, INTENT(OUT) :: points

    p => head
    points = 1

    DO
       p => p % jump_f
       points = points + 1
       IF(ASSOCIATED(p,tail))EXIT
    ENDDO

  END SUBROUTINE enumerate
  !-----------------------------------------------------------------------------
  SUBROUTINE small2tiny_box(small, num_tiny)

    IMPLICIT NONE

    TYPE(box), POINTER   :: small, tin
    TYPE(point), POINTER :: b_head, tail, p, q, p_xmin, p_xmax, p_ymin, p_ymax

    INTEGER, INTENT(OUT) :: num_tiny
    INTEGER              :: mark, pts, elems, i

    REAL(KIND=8)         :: xmin, xmax, ymin, ymax, dx, dy

    ! suddivisione di un box contenente un lato intero in box piu' piccoli
    ! min_points e' un parametro in ingresso, eventualmente ottimizzabile 
    ! write(*,*) '    entering small2tiny' 
    b_head=> small % head 
    tail => small % c2
    IF(small % structured == 1)THEN
       mark = 1
    ELSE 
       mark = 0
    ENDIF
    CALL enumerate(b_head, tail, pts)

    elems = pts - 1
    num_tiny = elems 

    ALLOCATE(small % small_box(num_tiny))

    p => b_head

    DO i=1, num_tiny
       tin => small % small_box(i)
       tin % structured = mark
       tin % head => p
       tin % c1 => p  ! e' la "testa di lista" dello small_box, il primo corner
       tin % active = 0
       IF(i == num_tiny) EXIT
       tin % c2 => p % jump_f
       tin % num_points = 2
       p => p % jump_f
    ENDDO
    tin % c2 => tail 
    tin % num_points = pts - (num_tiny - 1)

    DO i=1, num_tiny
       tin => small % small_box(i)
       p => tin % c1
       q => tin % c2
       xmin= p % x(1); xmax = p % x(1)
       ymin= p % x(2); ymax = p % x(2)
       p_xmin => p; p_xmax => p
       p_ymin => p; p_ymax => p

       IF(q % x(1) > xmax) THEN
          xmax = q % x(1)
          p_xmax => q
       ENDIF
       IF(q % x(1) < xmin) THEN
          xmin = q % x(1)
          p_xmin => q
       ENDIF

       IF(q % x(2) > ymax) THEN
          ymax = q % x(2)
          p_ymax => q
       ENDIF
       IF(q % x(2) < ymin) THEN
          ymin = q % x(2)
          p_ymin => q
       ENDIF
       IF(tin % structured == 1)THEN
          xmin = xmin - p_xmin % height
          xmax = xmax + p_xmax % height
          ymin = ymin - p_ymin % height
          ymax = ymax + p_ymax % height
       ELSE
          ! alcuni segmenti delaunay (approx orizzontali o verticali) definiscono box
          ! molto sottili che ritardano la rilevazione dell'interferenza
          ! l'idea e' quella di avere un box circa quadrato
          dx = xmax - xmin
          dy = ymax - ymin
          IF(dx > 2*dy)THEN ! orizzontale
             ymin = (ymin + ymax - dx)/2
             ymax = ymin + dx
          ELSEIF(2*dx < dy)THEN ! verticale
             xmin = (xmin + xmax - dy)/2
             xmax = xmin + dy
          ENDIF
       ENDIF

       tin % corner1(1) = xmin
       tin % corner1(2) = ymin

       tin % corner2(1) = xmax
       tin % corner2(2) = ymax

    ENDDO


  END SUBROUTINE small2tiny_box
  !-----------------------------------------------------------------------------
  SUBROUTINE define_struct_box( head, tail, str_bx, enlarge)
    IMPLICIT NONE
    TYPE(point), POINTER     :: head, tail, p, p_xmin, p_xmax, p_ymin, p_ymax
    TYPE(box), INTENT(INOUT) :: str_bx
    REAL(KIND=8)             :: xmin, xmax, ymin, ymax
    INTEGER, INTENT(IN)      :: enlarge

    ! attraverso la testa di lista di un lato di partenza del fronte avanzante si
    ! accede a tutti i suoi punti; il box e' definito dai due corner = i punti di
    ! coordinate cartesiane estreme, che si individuano intersecando tra loro le
    ! rette a xmin,ymin e xmax,ymax. Le liste (prelevate dal layer) hanno la testa
    ! di lista che rappresenta un punto effettivo, il primo punto del lato
    p => head
    xmin= p % x(1); xmax = p % x(1)
    ymin= p % x(2); ymax = p % x(2)
    p_xmin => p; p_xmax => p
    p_ymin => p; p_ymax => p

    p => p % jump_f

    DO
       IF(p % x(1) > xmax) THEN
          xmax = p % x(1)
          p_xmax => p
       ENDIF
       IF(p % x(1) < xmin) THEN
          xmin = p % x(1)
          p_xmin => p
       ENDIF

       IF(p % x(2) > ymax) THEN
          ymax = p % x(2)
          p_ymax => p
       ENDIF
       IF(p % x(2) < ymin) THEN
          ymin = p % x(2)
          p_ymin => p
       ENDIF

       IF(ASSOCIATED(p,tail))EXIT
       p => p % jump_f
    ENDDO

    ! le coordinate minime e massime vengon allargate con il valore locale di
    ! altezza dei quadrilateri che verranno introdotti
    IF(enlarge == 1)THEN
       xmin = xmin - p_xmin % height
       xmax = xmax + p_xmax % height
       ymin = ymin - p_ymin % height
       ymax = ymax + p_ymax % height
    ENDIF

    str_bx % corner1(1) = xmin
    str_bx % corner1(2) = ymin

    str_bx % corner2(1) = xmax
    str_bx % corner2(2) = ymax

    str_bx % head => head
    str_bx % structured = 1

  END SUBROUTINE define_struct_box
  !-----------------------------------------------------------------------------
  SUBROUTINE define_del_box(head, tail, del_bx )
    IMPLICIT NONE
    TYPE(point), POINTER     :: head, tail, p
    TYPE(box), INTENT(INOUT) :: del_bx
    REAL(KIND=8)             :: xmin, xmax, ymin, ymax

    ! attraverso la testa di lista di un lato si
    ! accede a tutti i suoi punti; il box e' definito dai due corner = i punti di
    ! coordinate cartesiane estreme, che si individuano intersecando tra loro le
    ! rette a xmin,ymin e xmax,ymax. Le liste circolari (prelevate da grid) 
    p => head 
    xmin= p % x(1); xmax = p % x(1)
    ymin= p % x(2); ymax = p % x(2)

    p => p % jump_f

    DO
       IF(p % x(1) > xmax) THEN
          xmax = p % x(1)
       ENDIF
       IF(p % x(1) < xmin) THEN
          xmin = p % x(1)
       ENDIF

       IF(p % x(2) > ymax) THEN
          ymax = p % x(2)
       ENDIF
       IF(p % x(2) < ymin) THEN
          ymin = p % x(2)
       ENDIF

       IF(ASSOCIATED(p,tail))EXIT
       p => p % jump_f
    ENDDO

    del_bx % corner1(1) = xmin
    del_bx % corner1(2) = ymin

    del_bx % corner2(1) = xmax
    del_bx % corner2(2) = ymax

    del_bx % head => head

  END SUBROUTINE define_del_box
  !-----------------------------------------------------------------------------
  SUBROUTINE update_big_box(big, head)

    IMPLICIT NONE

    TYPE(box),   POINTER :: big, small
    TYPE(point), POINTER :: head, tail
    INTEGER              :: i, enlarge

    ! nei controlli successivi al primo, bisogna:
    !   - ridefinire i box strutturati, sia big che small
    !   - mantenere i box delaunay, sia big che small

    ! cancellare gli small box 

    IF(big % active == 1)THEN

       DO i=1, SIZE(big % small_box)
          small => big % small_box(i)

          IF(small % active /= 0)THEN
             DEALLOCATE(small % small_box)
             small % active = 0
          ENDIF

       ENDDO

       DEALLOCATE(big % small_box)
       big % active = 0

    ENDIF
    ! azzerare i big box

    ! ridefinire i big box attraverso le teste di lista del nuovo layer
    tail => head % jump_b

    enlarge = 1
    CALL define_struct_box(head, tail, big, enlarge) 

  END SUBROUTINE update_big_box
  !-----------------------------------------------------------------------------

END MODULE interfer

