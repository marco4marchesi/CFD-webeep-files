MODULE fork

  USE eucl_delaunay
  USE grid_types
  USE wakes_gen
  USE init_data
  
  IMPLICIT NONE


  CONTAINS

  SUBROUTINE find_fork ( h, iter, nlayer, deltaj )

    IMPLICIT NONE

    TYPE (point), POINTER :: h
    INTEGER, INTENT(IN) :: iter
    INTEGER, INTENT(IN) :: nlayer
    REAL(KIND=8),INTENT(IN)  :: deltaj
    
    TYPE (point), POINTER :: p
    REAL(KIND=8) :: n, sum, sum1, h_need2c, h_n ! delta
    REAL(KIND=8) :: theta, height
    REAL(KIND=8), PARAMETER :: pi = 3.1415927464102d0
    INTEGER :: i

    !delta = grid % delta   ! (1+delta) = costante moltiplicativa delle altezze

    p => h
    DO ! calcolo degli angoli

       CALL angle( p % prev, p, p % next, theta )

       p % theta = theta
       p => p % next
       IF (ASSOCIATED(p,h))EXIT

    ENDDO

    p => h
    DO     
    
       IF (p % opnd == 1)THEN ! niente biforcazioni agli estremi
          ! punto iniziale del lato aperto
          p % theta = 0
          p % prev %theta = 0
          p => p % next
          CYCLE                ! passa ai punti successivi    
       ENDIF
       
       IF((p % cnd == 0).OR.(p % cnd < 0)) THEN 
          ! punti non impegnati o interessati da chiusura   
          IF ((p%theta == 0).AND.(p % next % theta == 0))THEN!||||||||||||||||

             IF(p% cnd /= 0 )THEN
                WRITE(*,*)'aia!'
                CALL waitme
             ENDIF
             n = 0  
          ELSE!|||||||||||||||||||||||||||||||||||||||||||||||||

             IF (p % theta /= 0) THEN

                h_need2c = (2 * MIN(p % leng, p%prev%leng) *COS(p % theta/2))/ &
                     SQRT(2 * (1 - COS(p % theta)) )
                ! altezza necessaria perche' la minima delle due lunghezze adiacenti 
                ! al punto si riduca a zero nel passare ad uno strato successivo
             ENDIF
             !----------------------------------------------------    
             ! calcolo quanti strati (i) sono necessari per superare l'altezza
             ! determinata ( si arriva fino a sum1), considerando anche la legge 
             ! di variazione delle altezze

             height = p % height
             sum = p % past_height
             i = 0        
             DO   
                i = i+1
                sum1 = sum +  height
                height =  height *(1 + deltaj) ! nuova altezza,passo seguente

                IF (p % theta < 0)THEN  ! chiusure--------

                   h_n = p % past_height + h_need2c
                   ! 'altezza' rispetto al contorno

                   IF ( sum1 >= (p%past_height + h_need2c) )EXIT

                   !ELSE                    ! aperture---------
                   !h_need2o = 0 ! da togliere
                   !h_n = p % past_height + h_need2o
                   !IF ( sum1 >= (p % past_height+h_need2o) )EXIT

                ENDIF

                IF (i > nlayer-iter+1)EXIT

                sum = sum1 
             ENDDO
             !----------------------------------------------------
             IF ( i > nlayer - iter + 1 )THEN!""""""""""""""""""""""

                p % cnd = 0 

             ELSE!"""""""""""""""""""""""

                IF(p % theta < 0)THEN !p % theta <0  chiusura del contorno

                   !WRITE(*,*)'chiusura'
                   IF (sum == sum1 )THEN
                      p % cnd = -1
                      ! per evitare la singolarita' di sum -sum1
                      CYCLE
                   ENDIF

                   IF ((h_n - sum)/(sum1-sum) >= 0.65)THEN
                      p % cnd = -i              !ritardo della condensazione
                   ELSE
                      p % cnd = - MAX(1,i - 1)  !anticipo della condensazione
                   ENDIF
                   ! il punto di chiusura si trova all'altezza complessiva h_n
                   ! che e' compresa tra sum e sum1
                   ! si distribuiscono i casi intermedi immaginando che la chiusura
                   ! avvenga allo strato precedente(sum) o successivo(sum1)

                ELSE

                   p % cnd = 0

                ENDIF
             ENDIF!""""""""""""""""""""""""

          ENDIF!||||||||||||||||||||||||||||||||||||||||||||||||||||


       ENDIF
       p => p % next

       IF (p % next % opnd == 1)THEN !niente biforcazioni agli estremi
          ! punto finale del lato aperto
          p % theta = 0
          EXIT

       ELSEIF(ASSOCIATED(p,h))THEN

          EXIT

       ENDIF


    ENDDO


  END SUBROUTINE find_fork

  !°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

      subroutine angle ( p1, p2, p3, theta )
!
      implicit none
!
      type (point), pointer :: p1, p2, p3
      real*8, intent(out) :: theta 
!
      real*8, dimension(2) :: q1, q2
      real*8, dimension(2) :: l_23, l_12
!
      real*8 :: r1, r2, leng_q12
      real*8, parameter :: pi = 3.1415927464102d0
!
      integer :: sign
!
      r1 = dsqrt( segment_length ( 2, p1, p2 ) ) 
      r2 = dsqrt( segment_length ( 2, p2, p3 ) )
!
      if ( (r1 .lt. 1.d-15).or.(r2 .lt. 1.d-15) ) then
         write (*,*) 'ERROR. Computing angles:' 
         write (*,*) 'Nodes duplication at', p2 % x
         write (*,*) ''
         stop
      end if

      l_23 = p3 % x - p2 % x
      l_12 = p2 % x - p1 % x

      l_23  = l_23 / r2
      l_12  = l_12 / r1
      ! vettori normalizzati a uno e allineati con i due segmenti

      theta = 90
      theta = -theta * pi / 180.d0 
      ! rotazione oraria(positiva antioraria)

      call rot( theta, l_12, q1 )
      call rot( theta, l_23, q2 )
      ! q1 e q2, versori normali ai due segmenti convergenti nel punto

      p2 % nv1(1,:) = q1
      p2 % nv1(2,:) = q2

      leng_q12 = ( q2(2) - q1(2) )**2 + ( q2(1) - q1(1) )**2
      ! distanza tra le due 'punte'

      if ( abs(leng_q12 - 4) .lt. 1.d-05 ) then

         !write(*,*)'warning(angle forced to ''pi'')!',p2 % x
         !write(*,*)'...hope it is an edge or an open straight side...'
         theta = pi
         ! i due versori sono allineati con verso opposto
      else 

         theta = acos(1 - ((leng_q12) / 2 )) !rad !teorema di carnot

         call theta_sign ( q1, q2, theta, sign )
         theta = sign * theta
         ! determina la concavita' o convessita' del contorno
      endif
!
      return
      end subroutine angle

  !°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

  SUBROUTINE theta_sign ( q1, q2, theta, sign )

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(2), INTENT(IN) :: q1
    REAL(KIND=8), DIMENSION(2), INTENT(IN) :: q2
    REAL(KIND=8), INTENT(IN) :: theta
    INTEGER, INTENT(OUT) :: sign

    REAL(KIND=8), DIMENSION(2) :: qr
    REAL(KIND=8) :: l1, l2

    CALL rot( theta, q1, qr)
    l1 = SQRT(( q2(1) - qr(1))**2 + ( q2(2) - qr(2) )**2)
    ! rotazione antioraria del versore normale q1; si ottiene qr
    ! calcolo della distanza l1 fra le punte di qr e q2

    CALL rot( -theta, q1, qr)
    l2 = SQRT(( q2(1) - qr(1) )**2 + (q2(2) - qr(2) )**2)
    ! rotazione oraria del versore normale q1; si ottiene qr
    ! calcolo della distanza l2 fra le punte di qr e q2

    IF ( l1 .le. l2 ) THEN
       sign = 1
       !apertura del contorno,rotazione ulteriore negativa
    ELSE IF ( l1 .gt. l2 ) THEN
       sign = -1
       !chiusura del contorno,rotazione ulteriore positiva,theta+ antiorario
    ENDIF

  END SUBROUTINE theta_sign

  !°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°




  SUBROUTINE insert_next(head,q)

    IMPLICIT NONE
    TYPE(point),POINTER ::head,p,q

    ! inserisce un nuovo punto in avanti (next) nella lista circolare 
    ! condizione necessaria: esistenza del punto head 

    IF (.NOT. ASSOCIATED(head % next)) THEN

       head % next => q
       q % prev => head

       head % prev => q
       q % next => head
       RETURN

    ENDIF

    p => head % prev
    q % prev => p
    p % next => q

    q % next => head
    head % prev => q 

  END SUBROUTINE insert_next





  SUBROUTINE insert_next_evo ( option, p, q, r, head )
  
  ! Inserts points in the direction normal to the wall
  !-------------------------------------------------------------------
  IMPLICIT NONE
  
  CHARACTER (LEN=4) :: option
  TYPE(point), POINTER :: head, p, q
  TYPE(point), POINTER, OPTIONAL :: r
  
  REAL(KIND=8),DIMENSION(2):: q_med
  REAL(KIND=8),DIMENSION(2):: q1, q_med1
  REAL(KIND=8),DIMENSION(2):: q2, q_med2

  REAL(KIND=8) :: theta1, theta2
  REAL(KIND=8) :: alpha, l1, l2

  REAL(KIND=8),PARAMETER :: pi = 3.1415927464102d0
  !-------------------------------------------------------------------
   
    SELECT CASE ( option )
    
    CASE ('frst')  ! Insert only along the first normal (no stretch)

       q1 = p % nv1(1,:)
       q % x = p % x + q1 * p % height
       
       CALL insert_next( head, q )


    CASE ('scnd')  ! Insert only along the second normal (no stretch)

       q2 = p % nv1(2,:)
       q % x = p % x + q2 * p % height
       
       CALL insert_next( head, q )


    CASE ('both')  ! Insert both (no stretch)

       IF (present(r)) THEN
       
          q1 = p % nv1(1,:)
          q2 = p % nv1(2,:)

          q % x = p % x + q1 * p % height
          r % x = p % x + q2 * p % height

          CALL insert_next( head, q )
          CALL insert_next( head, r )
       
       ELSE
          WRITE(*,*)'too few arguments for insert_next_evo both'
          STOP
       ENDIF


    CASE ('half')  ! Insert intermediate (stretch)

       q1 = p % nv1(1,:)  
       q2 = p % nv1(2,:)  
       
       ! Media delle due direzioni normali       
       q_med = (q1 + q2)/2 
       q_med = q_med/(SQRT(DOT_PRODUCT(q_med,q_med)))


       IF (p % theta /= 0)THEN
          q % x = p % x + q_med * p % height / ABS(COS(p % theta / 2))
       ELSE
          q % x = p % x + q_med * p % height
       ENDIF

       CALL insert_next(head,q) 


    CASE ('give')  ! Analogo di 'half' ma con scostamento rispetto alla normale

       q1 = p % nv1(1,:)  
       q2 = p % nv1(2,:)
        
       !q_med = (q1 + q2)/2 
       !q_med = q_med/(SQRT(DOT_PRODUCT(q_med,q_med))) 

       !theta1 = (p % theta - p % prev % theta)/3 !+ 2*pi/180.d0 
       !theta2 = (p % theta - p % next % theta)/3 !+ 2*pi/180.d0 
       
       theta1 = 20*pi/180d0
       theta2 = theta1
       
       CALL rot( theta1, q1, q_med1) ! + antiorario
       CALL rot(-theta2, q2, q_med2)

       ! q_med = (q_med1 + q_med2)/2

       l1 = p % prev % leng
       l2 = p % leng

       !alpha = 1+(2*(l1-l2)/(l1+l2))
       alpha =(2*(l1-l2)/(l1+l2))

       !q_med = ((alpha/2)**2)*q_med1+(1-((alpha/2))**2)*q_med2

       q_med = (1-alpha)*q_med1 + (1 + alpha)*q_med2
       q_med = q_med/(SQRT(DOT_PRODUCT(q_med,q_med)))

       !if((p%x(1)<1.005).and.(p % x(1) >0.99))then
       !write(*,*)'***theta ',p%theta
       !write(*,*)'theta1,theta2',theta1,theta2
       !write(*,*)'q1,q1',q1,q2
       !write(*,*)'qmed1,qmed2',q_med1,q_med2
       !write(*,*)'l1,l2,alpha',l1,l2,alpha
       !write(*,*)'QMED-height',q_med,p % height
       !endif
       !q % x = p % x + q_med * p % height / ABS(COS(p % theta / 2))

       q % x = p % x + q_med * p % height / ABS(COS((theta1+theta2)/2))
       
       CALL insert_next(head, q)
       
    END SELECT

  END SUBROUTINE insert_next_evo
  
  
  
  
  
  !°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  SUBROUTINE waitme(string)
    IMPLICIT NONE
    CHARACTER(LEN=10),OPTIONAL::string

    IF(present(string))then

       WRITE(*,*)string,';continue?'

    ELSE

       WRITE (*,*)'continue?'

    ENDIF
    READ(*,*)

  END SUBROUTINE waitme
  !-------------------------------------------------------------------------
  SUBROUTINE condense(new_head,gamma_lim3, deltaj)!,j)

    IMPLICIT NONE

    TYPE(point),POINTER :: p1, p2, p3, o1, o2, o3, q, new_head! m

    REAL(KIND=8):: theta_lim, thetad_now, theta3, height_aft,&
         height_now, t, theta2, height_lim, height, &
         theta5, theta4, upper_min, below_min, theta_cpl1, theta_cpl2, leng
    REAL(KIND=8),PARAMETER    :: pi = 3.1415927464102d0

    REAL(KIND=8),INTENT(IN) :: gamma_lim3, deltaj
    INTEGER                   :: testa, n, nh
    INTEGER,DIMENSION(2)      :: prior        

    LOGICAL,PARAMETER         :: ascii =.TRUE.
    !..........................................
    ! scorre la lista dei punti inseriti 'condensando' 2 o 3 punti ogni volta che
    ! l'indicatore cnd (stabilito da find_fork) passa da -1 a 0 inserendo un nuovo
    ! strato di punti
    ! due criteri di condensazione: per concavita' ampie o ridotte 
    !..........................................


    nh = 0 
    p1 => new_head
    DO
       nh = nh + 1
       p1 => p1 % next
       IF(ASSOCIATED(p1,new_head))EXIT
    ENDDO

    p1 => new_head % prev
    p2 => new_head
    p3 => new_head % next

    o1 => p1 % old
    o2 => p2 % old
    o3 => p3 % old
    !write(*,*)'condensating...'

    testa = 0


    n = 0

    DO 
       n = n+1
       IF(n > nh)THEN
          !WRITE(*,*)'WARNING condense: forced exit!'
          EXIT

       ENDIF

       IF(p2 % opnd == 1)THEN
          ! per evitare le condensazioni su tre punti che coinvolgano gli estremi dei
          ! lati aperti
          n = n + 1    ! estremo iniziale
          nh = nh - 2  ! estremo finale

          p1 => p1 % next % next 
          p2 => p2 % next % next
          p3 => p3 % next % next

          o1 => p1 % old
          o2 => p2 % old
          o3 => p3 % old

          CYCLE

       ENDIF

       IF(o2 % theta <= 0) THEN!*****************************

          height = SQRT(segment_length(2,o2,p2))

          IF (o2%cnd <=0)THEN

             height = p2 % height

             prior(1) = 0 
             prior(2) = 0

             ! priorita' di intervento di condensazione
             ! = 1 si ha intersezione di due altezze gia' avvenuta 
             !     con questo inserimento di punti
             ! = 2 l'intersezione avverra' prima che sia raggiunta la meta'(vd.t) della
             !     successiva altezza
             ! = 0 l'intersezione avviene dopo l'inserimento del successivo strato 

             !**********************************************************
             ! 1° Calcolo 

             theta_lim = (pi - ABS(o1 % theta))/2
             ! angolo o2 - o1 - p1 , corrispondente all'angolo o2 - o1 - p2
             !   se p2 coincidesse con p1, cioe' se fosse avvenuta l'intersezione

             theta2 = (pi - ABS(o2 % theta))/2
             ! angolo p2 - o2 - o1

             CALL angle(o1,p2,o2,theta3)

             theta3 = (pi - ABS(theta3))
             ! angolo o1 - p2 - o2

             thetad_now = pi - theta2 - theta3 
             ! angolo o2 - o1 - p2
             ! <  theta_lim  intersezione non avvenuta
             ! >= theta_lim        "      avvenuta

             t = 1 + 0.5 *(1 + deltaj)
             IF( thetad_now >= theta_lim)THEN 
                !write(*,*)'intersezione SUBITO'
                prior(1) = 1

             ELSE  ! l'intersezione avverra' a breve(prossimi/o strati/o)

                height_lim = o1%leng*SIN(theta_lim)*SIN(theta2)/SIN(theta3)
                ! distanza del punto di intersezione dal segmento o1 - o2

                height_now = SQRT(segment_length(2,p2,o2)) * SIN(theta2 )
                ! distanza del punto p2 dal segmento o1 - o2

                theta4 = (thetad_now+theta_lim)/2

                theta_cpl1 = pi - theta2 - theta4

                below_min = MIN(theta2, theta4, theta_cpl1)

                leng=SQRT(segment_length(2,p1,p2))

                height=o1%height*(1+deltaj)

                theta_cpl2=ASIN(leng*SIN(theta2)/SQRT(leng**2+(height/SIN(theta2))**2 &
                     - 2*(leng*height*COS(theta2)/SIN(theta2))))

                theta5 = pi -theta2 - theta_cpl2

                upper_min= theta_cpl2       

                IF(upper_min <= (gamma_lim3*pi/180d0))THEN
                   !       write(*,*)'intersezione VICINA1 '
                   prior(1) = 1                  
                ELSE
                   !write(*,*)'intersezione lontana1'
                   prior(1) = 0
                ENDIF
             ENDIF

             ! 2° Calcolo
             theta_lim = (pi - ABS(o3 % theta) )/2
             theta2 = (pi - ABS(o2 % theta))/2
             CALL angle(o3,p2,o2,theta3)
             theta3 = (pi - ABS(theta3))
             thetad_now = pi - theta2 - theta3 


             IF( thetad_now >= theta_lim)THEN ! intersezione avvenuta
                !write(*,*)'intersezione subito'
                prior(2) = 1
             ELSE  ! l'intersezione avverra' a breve
                height_lim = o2%leng*SIN(theta_lim)*SIN(theta2)/SIN(theta3)

                height_now = SQRT(segment_length(2,p2,o2)) * SIN(theta2)

                height_aft = t * height_now


                theta4 = (thetad_now+theta_lim)/2
                theta_cpl1 = pi - theta2 - theta4
                below_min = MIN(theta2, theta4, theta_cpl1)


                leng=SQRT(segment_length(2,p2,p3))
                height=o3%height*(1+deltaj)
                theta_cpl2=ASIN(leng*SIN(theta2)/SQRT(leng**2+(height/SIN(theta2))**2  &
                     - 2*(leng*height*COS(theta2)/SIN(theta2))))
                theta5 = pi -theta2 - theta_cpl2


                upper_min= theta_cpl2

                IF(upper_min <= (gamma_lim3*pi/180d0))THEN
                   !        write(*,*)'intersezione VICINA'
                   prior(2) = 1                 
                ELSE
                   !write(*,*)'intersezione lontana2'
                   prior(2) = 0
                ENDIF
             ENDIF

             !**********************************************************
             IF((prior(1) /= 0).OR.(prior(2)/= 0))THEN!°°°°°°°°°°°°°°    
                IF((prior(1) /= 0).AND.(prior(2)== 0))THEN

                   !        WRITE(*,*)'substitute1'
                   CALL substitute_point2(p1, p2, q, new_head, testa)
                   nh = nh - 1

                ELSEIF((prior(1) == 0).AND.(prior(2)/= 0))THEN

                   !        WRITE(*,*)'substitute2'
                   CALL substitute_point2(p2, p3, q, new_head, testa)
                   nh = nh - 1

                   !ELSEIF((prior(1) /= 0).AND.(prior(2)/= 0))THEN
                ELSE!((prior(1) /= 0).AND.(prior(2)/= 0))THEN  
                   !        WRITE(*,*)'substitute3'
                   CALL substitute_point (p1, p2, p3, q, new_head, testa)
                   nh = nh - 2

                ENDIF

                q % cnd = 0 
                q % leng = MIN(SQRT(segment_length(2,q,q % next)) ,&
                     SQRT(segment_length(2,q,q % prev)))
                q % leng_max    = (q % prev % leng_max  + q % next % leng_max )/2 
                q % bou_leng1   = (q % prev % bou_leng1 + q % next % bou_leng1)/2
                q % height      =  height
                q % past_height = (q % prev % past_height + q % next % past_height)/2

                p1 => q
                p2 => q % next 
                p3 => q % next % next  

             ENDIF

          ENDIF
       ENDIF!**************************************************

       IF( n == nh) EXIT

       p1 => p1 % next 
       p2 => p2 % next
       p3 => p3 % next

       o1 => p1 % old
       o2 => p2 % old
       o3 => p3 % old

    ENDDO

  END SUBROUTINE condense

  !-------------------------------------------------------------------------

  SUBROUTINE substitute_point(p1,p2,p3,q,new_head,testa)

    IMPLICIT NONE

    TYPE(point),POINTER   :: p1, p2, p3, b1, b2, o1, o2, o3, q, head, new_head
    INTEGER,INTENT(OUT)         :: testa
    !EAL(KIND=8)                :: c, theta
    !EAL(KIND=8),DIMENSION(2,2) :: rot_mat
    REAL(KIND=8),DIMENSION(2)   ::  coo !displ,
    REAL(KIND=8),PARAMETER      :: pi = 3.1415927464102d0

    !.......................................
    ! sostituzione di tre punti con uno 
    ! in presenza di una testa di lista si mantiene sempre la testa
    ! altrimenti si eliminano i tre punti e se ne introduce uno nuovo
    !.......................................
    o1 => p1 % old
    o2 => p2 % old
    o3 => p3 % old

    head => new_head % old

    b1 => p1 % prev
    b2 => p3 % next

    NULLIFY (o1 % ppp) ;NULLIFY (o2 % ppp) ;NULLIFY (o3 % ppp) 
    testa = 0

    o1 % what = 6
    o2 % what = 6

    IF(ASSOCIATED(o1,head).OR.ASSOCIATED(o2,head).OR.ASSOCIATED(o3,head))THEN
       q => new_head 
       ! il punto e ' gia' inserito
    ELSE
       q => new_point()
    ENDIF

    o1 % ppp => q;   o2 % ppp => q;   o3 % ppp => q

    !write(*,*)'o1 x',o1 % x
    !write(*,*)'o2 x',o2 % x
    !write(*,*)'o3 x',o3 % x

    coo= (p1 % x + p2 %x + p3 % x)/3
    ! coordinata media

    IF(ASSOCIATED(o1,head))THEN!-------------------------------

       testa = 1
       DEALLOCATE(p2); DEALLOCATE(p3)
       new_head % next => b2
       b2 % prev => new_head

       new_head % old1 => o1
       new_head % old => o3

    ELSEIF(ASSOCIATED(o2,head))THEN!---------------------------

       testa = 2
       DEALLOCATE(p1); DEALLOCATE(p3)

       new_head % next => b2 ; new_head % prev => b1
       b2 % prev => new_head ; b1 % next => new_head

       new_head % old1 => o1
       new_head % old => o3

    ELSEIF(ASSOCIATED(o3,head))THEN!---------------------------

       testa = 3
       DEALLOCATE(p2); DEALLOCATE(p1)

       new_head % prev => b1
       b1 % next => new_head 

       new_head % old1 => o1
       new_head % old => o3

    ELSE
       testa = 0
       DEALLOCATE(p1); DEALLOCATE(p2); DEALLOCATE(p3)

       b1 % next => b2 ; b2 % prev => b1 !serve a insert_between
       CALL insert_between(b1,b2,q)

       q % old1 => o1
       q % old => o3

    ENDIF

    q % x = coo

  END SUBROUTINE substitute_point

  !------------------------------------------------------------

  SUBROUTINE substitute_point2(p1,p2,q,new_head,testa)

    IMPLICIT NONE

    TYPE(point),POINTER         :: p1, p2, b1, b2, o1, o2, q, head, new_head

    INTEGER, INTENT(OUT)        :: testa

    REAL(KIND=8),DIMENSION(2)   :: coo 
    REAL(KIND=8),PARAMETER      :: pi = 3.1415927464102d0

    ! sostituzione di due punti con uno solo

    o1 => p1 % old
    o2 => p2 % old

    head => new_head % old

    b1 => p1 % prev
    b2 => p2 % next

    NULLIFY (o1 % ppp) ; NULLIFY (o2 % ppp) 
    testa = 0 
    o1 % what = 6

    IF(ASSOCIATED(o1, head).OR.ASSOCIATED(o2, head))THEN
       q => new_head

    ELSE
       q => new_point()
    ENDIF

    o1 % ppp => q;   o2 % ppp => q

    !write(*,*)'o1 x',o1 % x
    !write(*,*)'o2 x',o2 % x

    coo= (p1 % x + p2 %x )/2

    IF(ASSOCIATED(o1, head))THEN!-------------------------------

       testa = 1
       DEALLOCATE(p2)

       new_head % next => b2
       b2 % prev => new_head

       new_head % old1 => o1
       new_head % old => o2

    ELSEIF(ASSOCIATED(o2, head))THEN!---------------------------

       testa = 2
       DEALLOCATE(p1)

       new_head % next => b2 ; new_head % prev => b1
       b2 % prev => new_head ; b1 % next => new_head

       new_head % old1 => o1
       new_head % old => o2

    ELSE
       testa = 0
       DEALLOCATE(p1); DEALLOCATE(p2)

       b1 % next => b2 ; b2 % prev => b1 !serve a insert_between
       CALL insert_between(b1, b2, q)

       q % old1 => o1
       q % old => o2

    ENDIF

    q % x = coo

  END SUBROUTINE substitute_point2

  !--------------------------------------------------------------------------

  SUBROUTINE point_more1 ( p1, q, p2 )
!
!   inserimento di un punto tra due che hanno distanza maggiore di quella di
!   riferimento. collegamenti nello strato e al di fuori( verso il precedente ) 
!
    IMPLICIT NONE

    TYPE(point), POINTER :: p1, p2, q
    TYPE(point), POINTER :: o1, o2
    REAL(KIND=8) :: theta

    o1 => p1 % old
    o2 => p2 % old 

    o1 % what = 10

    o1 % ppp  => p1
    o1 % ppp1 => q
    o2 % ppp => p2
    o2 % ppp1 => q

    CALL eq_ang ( p1, q, p2 ) ! percorrendo i tre punti in ordine, si ha la stessa
    ! variazione di angolo fra i segmenti 
    !CALL cirq(p1,q,p2)
    ! dispone i punti lungo un arco di cerchio

    !  Re-assign pointers, next and prev of p1, p2 ad q such that
    !  q is between p1 and p2
    CALL insert_between( p1, p2, q )

    q % old1 => o1
    q % old => o2

    q % leng  = SQRT(segment_length(2, q, p2))
    p1 % leng = SQRT(segment_length(2, q, p1))

    q % leng_max  = (p1 % leng_max  + p2 % leng_max ) /2
    q % bou_leng1 = (p1 % bou_leng1 + p2 % bou_leng1)/2
    q % height = (p1 % height + p2 % height)/2
    q % past_height = (p1 % past_height + p2 % past_height)/2

    CALL angle ( p1, q, p2, theta )
    q % theta = theta

    CALL angle ( p1 % prev, p1, q, theta )
    p1 % theta = theta

    CALL angle ( q, p2, p2 % next, theta )
    p2 % theta = theta


  END SUBROUTINE point_more1
  
  
  
  
  
  SUBROUTINE cirq(p,q,r)
    IMPLICIT NONE

    TYPE(point),POINTER       :: p,q,r
    REAL(KIND=8)              :: radius1,radius2,theta1,theta2,gamma1,gamma2
    REAL(KIND=8),DIMENSION(2) :: vect
    REAL(KIND=8),PARAMETER    :: pi = 3.1415927464102d0
    ! subroutine di calcolo della coordinata di un punto di condensazione
    ! disposizione del punto secondo un arco di cerchio 

    !angolo parziale di passaggio diretto
    CALL angle(p % prev,p,r,theta1)
    CALL angle(p,r,r%next,theta2)

    gamma1 = 2 * ABS(theta1)
    gamma2 = 2 * ABS(theta2)
    !raggio del cerchio 'osculatore'
    radius1 = SQRT(segment_length(2,p,r)) / (SQRT(2*(1- COS(gamma1))))
    radius2 = SQRT(segment_length(2,p,r)) / (SQRT(2*(1- COS(gamma2))))

    gamma1 = MIN(gamma1,gamma2)
    IF (gamma1 == gamma2)THEN
       radius1 = radius2
    ENDIF
    vect = (r%x - p%x)/SQRT(segment_length(2,p,r))


    CALL rot((gamma1-pi)/2,vect,vect)

    vect = vect * radius1 * gamma1 /2

    q % x = p % x + vect

  END SUBROUTINE cirq





  SUBROUTINE eq_ang(p,q,r)

    IMPLICIT NONE

    TYPE(point),POINTER       :: p, q, r

    REAL(KIND=8)              :: leng, theta1, theta2, gamma
    REAL(KIND=8),DIMENSION(2) :: vect 
    REAL(KIND=8),PARAMETER    :: pi = 3.1415927464102d0
    ! subroutine calcolo coordinata del nuovo punto
    ! lo dispone in modo tale che il suo angolo e' identico ai due adiacenti

    ! p, r  i due punti esistenti
    ! q     punto da inserire

    leng = dsqrt( segment_length(2,p,r) )

    CALL angle(p % prev, p, r, theta1)
    CALL angle(p, r, r%next, theta2)

    gamma = MIN(ABS(theta1),ABS(theta2))/2

    vect = (r % x - p % x) / SQRT(segment_length(2, p, r))
    !versore allineato con i due punti di partenza
    CALL rot(-gamma, vect, vect)
    ! rotazione oraria(il nuovo punto va all'esterno della congiungente p - r
    vect = vect * leng /(2* COS(gamma)) 
    q % x = p % x + vect

  END SUBROUTINE eq_ang





  SUBROUTINE rarefactions ( head )

    IMPLICIT NONE

    TYPE (point), POINTER :: head
    TYPE (point), POINTER :: p, q, r
    INTEGER :: head_involved

    head_involved = 0
    p => head
    DO

       IF ( p % theta .gt. 1.d-15 ) THEN 
      !IF((p % theta > 0).AND.(p % mark == 0))THEN 
!
!    +------------------------+
!    |  caso del lato aperto  |
!    +------------------------+
!
          IF ( (p % opnd .gt. 0) .AND. (ASSOCIATED(p % next, head)) ) EXIT
!
!    +----------------------------------------------------------------+
!    |  Confronto di lunghezza massima, aggiunge un punto se il caso  |
!    +----------------------------------------------------------------+
!
          IF ( (p % leng - p % leng_max) .gt. 1.d-08 ) THEN
!
             IF ( ASSOCIATED(p % old1) .OR. ASSOCIATED(p % next % old1) ) THEN
             !IF((p % old % mark ==1).OR.ASSOCIATED(p % old1).OR.ASSOCIATED(p % next % old1))THEN
             !write(*,*)'esizialissimo caso!!!'
                p => p % next
                IF (ASSOCIATED(p, head)) EXIT
                CYCLE
             END IF
!
             r => p % next
             q => new_point()
             CALL point_more1 ( p, q, r ) 
!
!    +---------------------------------------------------------------------------+
!    |  se un intervento avviene sulla testa bisogna trattare opportunamente un  |
!    |  intervento analogo alla fine( per evitare spostamenti di ppp1)           |
!    |  infatti i puntatori sono ppp e ppp1 ma qui i collegamenti sarebbero 3!!  |
!    +---------------------------------------------------------------------------+
!
             IF ( ASSOCIATED(p, head) )  head_involved = 1
!
!    +---------------------------------------------------------------------------+
!    |  con questo si ovvia temporaneamente al pasticcio di una testa compresa   |
!    |  tra due aggiunte di punti!!!!                                            |
!    +---------------------------------------------------------------------------+
!
             IF ( ASSOCIATED(p % next % next, head) .AND. (head_involved .eq. 1) ) THEN
                NULLIFY(head % ppp1)
                head % old % ppp1 => head % next
                ! ristabilisce una successione naturale dei ppp1
             END IF
!
             p => q
!
!    +---------------------------------------------------------------------------+
!            poi si controlla il successivo di q!
!
!           ELSEIF((p%leng < p%leng_min2o).AND.&
!                  (p % next % leng < p % next % leng_min2o))THEN
!              
!              ! confronto di lunghezza minima( non attivato, leng_min2o = 0!)
!              
!              IF((p%next%opnd > 0))EXIT
!              ! caso del lato aperto
!              
!              CALL point_less1(p,p % next,p % next % next)
!              ! elimina p % next       
!              
!              ! mi trovo naturalmente a controllare il successivo
!    +---------------------------------------------------------------------------+
!
          END IF
       END IF

       p => p % next
       IF ( ASSOCIATED(p, head) ) EXIT

    END DO

  END SUBROUTINE rarefactions





  SUBROUTINE point_less1(p1,p2,p3)

    IMPLICIT NONE

    TYPE(point),POINTER :: p1,p2,p3,o2
    REAL(KIND=8) :: theta
    !..................................
    ! riduce di uno il numero di punti
    ! dopo un confronto di lunghezza minima 

    o2 => p2 % old 

    o2 % what = 11

    o2 % ppp => p3
    o2 % ppp1 => p1

    CALL extract_point(p2)
    DEALLOCATE(p2)

    p1 % old1 => o2
    p3 % old1 => o2

    p1 % leng = SQRT(segment_length(2,p1,p3))

    CALL angle(p1 % prev,p1,p3,theta)

    IF(p1 % opnd == 0)THEN

       p1 % theta = theta

    ENDIF

    CALL angle(p1,p3,p3 % next,theta)

    IF(p3 % opnd == 0)THEN

       p3 % theta = theta

    ENDIF


  END SUBROUTINE point_less1

  !############################################################################
  !     UNUSED PROCEDURES   !
  !############################################################################
  SUBROUTINE fuse_top(head)
    IMPLICIT NONE
    TYPE(point),POINTER :: head,p1,p2,q,endh



    REAL(KIND=8),PARAMETER :: pi = 3.1415927464102d0

    INTEGER :: c1,c2,n,sis
    REAL(KIND=8) :: k,leng1,theta

    write(*,*)'fuse_top '
    write(*,*)'aggiorna cnd!!'
    stop
    p1 => head
    p2 => head % next

    endh => p1 

    n = 0
    do 

       c1 = p1 % cnd ; c2 = p2 % cnd
       !if ((c1 /= 0).AND.(c2 /= 0))THEN  
       if ((c1 == 1).AND.(c2 == -1))then
          q => new_point()
          CALL insert_between(p1,p2,q)
          !inserzione di punti mantenendo le tangenti
          sis = 0

          call angle(p1%prev,p1,p2,theta)
          p1 % theta = theta
          call angle(p1,p2,p2%next,theta)
          p2 % theta = theta
          leng1= SQRT(segment_length(2,p1,p2))/2

          k = leng1*TAN((p1%theta+p2%theta)/3) 

          q % x = (p1%x+p2%x)/2 + (p1%nv1(2,:)*k*0.75)

          p1 % cnd = 0
          p2 % cnd = 0 
          q % cnd = 0
          q % height = (p1 % height + p2 % height)/2
          q % past_height =(p1 %past_height + p2 % past_height)/2 
          q % leng_max = (p1 % leng_max  + p2 % leng_max )/2

          write(*,*)'   p1',p1 % x,p1 % leng

          write(*,*)'   p2',p2 % x,p2 % leng
          write(*,*)'   q',q % x,q%leng
          write(*,*)'-----'
          p1 % next1 => p2
          p2 % prev1 => p1
          call angle(p1,q,p2,theta)
          q % theta = theta
          n = n+1
          p1 => p2  
          p2 => p2 % next 
       endif

       p1 => p1 % next 
       p2 => p2 % next
       if (associated(p1,endh))exit

    enddo
    write(*,*)'         ',n,'points'
  END SUBROUTINE fuse_top
  !°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  SUBROUTINE fuse_bot(new_head)
    IMPLICIT NONE

    TYPE(point),POINTER :: new_head,p,q,m!,r,s
    REAL(KIND=8),PARAMETER :: pi = 3.1415927464102d0
    INTEGER :: n
    REAL(KIND=8) :: leng
    !triangolo cuscinetto tra due strati ,riduce di uno il numero di punti

    p => new_head

    write(*,*)'fuse bot'

    n = 0

    DO 
       IF(p % theta < 0) THEN!*****************************

          IF ((p % cnd == 1))THEN!°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

             write(*,*)'p % prev',p %prev % x
             write(*,*)'p',p % x
             write(*,*)'p % next',p %next % x
             n = n+1
             write(*,*)'n',n
             NULLIFY(p % next % prev)
             NULLIFY(p % next % prev1)
             p % next % prev1 => p
             p % next % prev => p % prev

             NULLIFY(p % prev % next)
             NULLIFY(p % prev % next1)
             p % prev % next1 => p
             p % prev % next => p % next

             !call angle(p % prev % prev,p % prev,p % next)
             !call angle(p % prev,p % next,p % next % next)
             p % prev % cnd = 0
             p % next % cnd = 0
             p % cnd = 0

             leng = SQRT(segment_length(2,p % prev,p % next))
             p % prev% leng =MIN(leng,p % prev % leng)
             p % next% leng =MIN(leng,p % next % leng)

             write(*,*)'p pre next',p % prev % next%x
             write(*,*)'p nex prev',p % next % prev %x
             if (associated(p,new_head))then
                !cambiare l'identità della head
                write(*,*)'head involved!'
                m => p % prev
                q => new_point()

                q % x = m % x
                q % old => m % old
                q % prev  => m % prev 
                q % prev1 => m % prev1
                q % next => m % next
                !q % next1 => m % next1
                write(*,*)'*'
                m % x = new_head % x
                m % old => new_head % old
                m % prev => new_head
                NULLIFY (m % prev1)
                m % next => q % next 
                NULLIFY (m % next1)
                write(*,*)'**'
                new_head % x = q % x
                new_head % old => q % old
                new_head % prev => q % prev
                NULLIFY (new_head % prev % next)
                new_head % prev % next => new_head
                new_head % prev1 => q % prev1
                new_head % next => q % next
                new_head % next1 => m
                write(*,*)'***'
                new_head % prev % cnd = 0
                !s => new_point()
                !r => new_point()

                ! call copy_point(p % prev,r)
                !write(*,*)'copy1'
                !call copy_point(p,s)
                !call copy_point(r,p)
                !call copy_point(s,p%prev)
                write(*,*)'done,head',new_head % x,new_head % next % x
             endif
             p => new_head !% next
             !if (associated(p,new_head))then

             !write(*,*)'exit previsto,n',n
             !exit
             !endif 
          ENDIF!°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
       ENDIF!**************************************************


       p => p % next


       IF(ASSOCIATED(p,new_head))then
          write(*,*)'fuse_bot done:',n,'bridges'

          EXIT
       ENDIF
    ENDDO
    write(*,*)'fine'
  END SUBROUTINE fuse_bot

  SUBROUTINE condense2(new_head)

    IMPLICIT NONE
    TYPE(point),POINTER :: p1,p2,p3,p4,p5,p6,new_head,q,p,endp

    REAL(KIND=8) :: lmax,lmin,th1,th2,lmax1,lmin1,seg,gamma,radius,lcirc,theta
    REAL(KIND=8),PARAMETER :: pi = 3.1415927464102d0
    INTEGER :: testa,i
    !gestione della condensazione su 5 punti invece di 3

    p1 => new_head % prev % prev
    p2 => new_head % prev
    p3 => new_head
    p4 => new_head % next
    p5 => new_head % next % next

    write(*,*)'condensating..2.'
    testa = 0
    i = 0
    endp => p4
    DO 
       i = i + 1

       IF(p3 % theta < 0) THEN!*****************************

          IF ((p3%cnd /= 0))then
             write(*,*)'p3 % cnd ',p3 % cnd ,p3 % x

             !lunghezze : percorso diretto lmin1
             !           percorso lungo i tratti lmax1
             !   lmin   minima lunghezza  di segmento incontrata
             !   lmax   massima lunghezza  di segmento incontrata  
             lmin1 = SQRT(segment_length(2,p1,p5))
             lmax1 = 0
             p => p1
             q => p1 % next
             th2 = ABS(p % theta)
             lmax = SQRT(segment_length(2,p1,p1%prev))
             lmin = lmax
             do 
                seg = SQRT(segment_length(2,p,q))
                lmax1 = lmax1 + seg 
                th2 = th2 + ABS(q % theta)! angolo di passaggio lungo la curva 
                lmax =MAX(lmax,seg)
                lmin=MIN(lmin,seg)
                p => p % next
                q => q % next
                if (associated(p,p5)) exit
             enddo
             !-------------------------
             call angle(p1 % prev,p1,p5,theta)
             p1 % theta = theta
             !angolo parziale di passaggio diretto
             th1 = ABS(p1 % theta)
             gamma = 2 * th1
             !raggio del cerchio 'osculatore'
             radius = lmin1 / (SQRT(2*(1- COS(gamma))))
             lcirc = radius * gamma ! lunghezza arco di circonferenza
             call angle(p1,p5,p5 % next,theta)
             p5 % theta = theta
             th1 = th1 + ABS(p5 % theta)
             !angolo complessivo di passaggio diretto

             write(*,*)'   tlmax1',lmax1,'tlmin1',lmin1
             write(*,*)'   th2',th2*180/pi,'th1',th1*180/pi
             write(*,*)'   lmax',lmax,'lmin',lmin
             write(*,*)'   radius',radius,'lcirc',lcirc
             write(*,*)'   lmax1-lcirc',lmax1-lcirc
             if (new_head % index /= 0)then

                write(*,*)'p1',p1 % cnd ,p1 % x  
                write(*,*)'p2',p2 % cnd ,p2 % x
                write(*,*)'p3',p3 % cnd ,p3 % x
                write(*,*)'p4',p4 % cnd ,p4 % x
                write(*,*)'p5',p5 % cnd ,p5 % x  
                call waitme
                if(associated(p1,new_head))then
                   testa = 1
                endif
                if(associated(p2,new_head))then
                   CALL change_head(new_head,p1)
                   testa = 2

                   p1 => new_head

                endif

                if(associated(p3,new_head))then
                   CALL change_head(new_head,p1)
                   testa = 3

                   p1 => new_head

                endif

                if(associated(p4,new_head))then
                   CALL change_head(new_head,p1)
                   testa = 4

                   p1 => new_head
                endif
                if (associated(p5,new_head))then
                   testa = 5
                endif
                p2 => p1 % next ;p3 => p2 % next ;p4 => p3 % next 
                p5 => p4 % next

                p6 => new_point()

                if (p3 % theta > 0)then
                   p2 % next => p3
                   p3 % prev => p2

                   p2 % ppp => p6
                   p3 % ppp => p6
                   p4 % ppp => p6
                   p6 % old => p4
                   p6 % old1 => p2
                else  

                   p2 % next => p4
                   p2 % next1 => p3

                   p4 % prev => p2
                   p4 % prev1 => p3


                   p1 % next1 => p2

                   p5 % prev1 => p4

                   p2 % ppp => p6
                   p4 % ppp => p6
                   p6 % old => p4
                   p6 % old1 => p2
                endif

                p6 % next => p5
                p5 % prev => p6
                p6 % prev => p1
                p1 % next => p6



                p6 % x = (p1 % x + p5 % x + p3 % x)/ 3

                CALL angle(p1,p6,p5,theta)
				
                p6 % theta = theta

                p6 % leng = MIN(SQRT(segment_length(2,p6,p1)),&
                     SQRT(segment_length(2,p6,p5)))
                p6 % leng_max = (p1 % leng_max+ p5 % leng_max)/2
                p6 % height =(p1 % height+ p5 % height)/2                
                p6 % past_height =(p1 %past_height+ p5 % past_height)/2

                p1 % cnd = 0
                p5 % cnd = 0
                p6 % cnd = 0
                if((i == 1).AND.(testa > 0))then
                   endp => p1
                elseif((i == 2).AND.(testa == 5))then 
                   endp => p1
                elseif((i == 3).AND.(testa == 5))then
                   endp => p1
                else 
                   endp => new_head   
                endif



                p1 => p1 % next
                p2 => p1 % next
                p3 => p2 % next ; p4 => p3 % next
                p5 => p4 % next 

             endif

          ENDIF!°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
       ENDIF!**************************************************

       p1 => p1 % next 
       p2 => p2 % next
       p3 => p3 % next
       p4 => p4 % next 
       p5 => p5 % next

       if(associated(p4,endp))exit


    ENDDO

    call waitme

  END SUBROUTINE condense2
  !°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  SUBROUTINE change_head(head,poin)
    IMPLICIT NONE
    TYPE(point),POINTER :: head,poin,q,r,a,b,c,d
    INTEGER :: order
    ! si assume che la testa non abbia
    ! prev1, next1 ,  

    if (associated(head%prev,poin))then
       ! a ,1-poin , 2-head , d
       a => poin % prev
       d => head % next
       order = 1
       !sub_order = 0
    elseif(associated(head % next,poin))then
       ! a ,1-head , 2-poin , d  
       a => head % prev
       d => poin % next
       order = 2
       !sub_order = 0
    else
       ! a , head , b ,*********, c , poin , d
       a => head % prev
       b => head % next
       c => poin % prev
       d => poin % next
       order = 3
    endif
    !cambiare l'identità della head
    !q,r punti di appoggio,raccolgono le informazioni di point,head
    q => new_point()
    r => new_point()
    ! al punto che riceve i dati viene annullato ogni puntatore
    !    copy    (poin --> q )

    CALL copy_point(poin,q)
    CALL copy_point(head,r)
    CALL copy_point(q,head) 
    CALL copy_point(r,poin)
    ! ma gli adiacenti dei due punti estirpati mantengono i puntatori 
    ! ai vecchi punti
    if (order == 1) then
       a % next => head
       d % prev => poin

       head % next => poin
       poin % prev => head

    elseif(order == 2) then
       a % next => poin
       d % prev => head

       head % prev => poin
       poin % next => head

    elseif(order == 3) then
       d % prev => head
       c % next => head
       a % next => poin
       b % prev => poin

    endif
    !esclusi i ppp1!!!!!!!!!!!!!!!!!!!!!!!

    NULLIFY(poin % old % ppp) ! q vecchio poin
    poin % old % ppp => poin ! nuovo poin
    if (associated(poin % old1))then
       if (.NOT.associated(poin % old1 % next,poin % old))then
          ! q risulta da una biforcazione
          NULLIFY(poin % old1 % next % ppp)
          poin % old1 % next % ppp => poin 

       endif
       NULLIFY(poin % old1 % ppp)
       poin % old1 % ppp => poin 

    endif


    NULLIFY(head % old % ppp) ! r vecchia head
    head % old % ppp => head ! nuova head
    if (associated(head % old1))then
       if (.NOT.associated(head % old1 % next,head % old))then
          ! r risulta da una biforcazione
          NULLIFY(head % old1 % next % ppp)
          head % old1 % next % ppp => head

       endif
       NULLIFY(head % old1 % ppp)
       head % old1 % ppp => head
    endif
    ! write(*,*)'head old',head % old % x
    ! write(*,*)'head ',head % x
    ! write(*,*)'poin old',poin % old % x
    ! write(*,*)'poin prev',poin% prev % x
    ! write(*,*)'poin',poin % x
    ! write(*,*)'poin next',poin% next% x
    ! write(*,*)'poin prevnext',poin%prev% next% x
    DEALLOCATE(q)
    DEALLOCATE(r)

  END SUBROUTINE change_head
  !°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  SUBROUTINE check_layer(head,test)
    IMPLICIT NONE

    !INTEGER,INTENT(IN) :: iter
    INTEGER,INTENT(OUT) :: test
    TYPE(point),POINTER :: head,p1,p2,p3,o2
    REAL(KIND=8),PARAMETER :: pi = 3.1415927464102d0

    write(*,*)'checking layer'
    p1 => head % prev
    p2 => head
    p3 => head % next
    !OPND!!!!!!!!!!!!!!!!!!!!!!!
    test = 0 
    DO 
       ! no variazioni di segno
       ! no variazioni eccessive di angolo: theta < 90°  , theta > -90°
       IF(ASSOCIATED(p2 % old))THEN 
          o2 => p2 % old
          IF ((o2 % theta < 0 ).AND.(p2 % theta > 0))THEN
             write(*,*)'o2',o2 % x,o2 % theta * 180 / pi
             write(*,*)'p2 ',p2 % x, p2 % theta * 180 / pi
             test = 1
             EXIT
          ENDIF

       ENDIF
       p1 => p1 % next
       p2 => p2 % next
       p3 => p3 % next
       IF (ASSOCIATED(p2,head))THEN
          EXIT
       ENDIF
    ENDDO

  END SUBROUTINE check_layer
  !°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  SUBROUTINE addi(head)
    IMPLICIT NONE
    TYPE(point),POINTER :: head,p,o,q,r

    p => head
    o => p % old
    DO

       IF ((o % cnd == 1).AND. (.NOT.(o % visited )))THEN  
          r => p % next
          q => new_point()
          CALL point_more1(p,q,r)
          p => q
          o % visited = .TRUE.
       ENDIF
       p => p % next
       IF(ASSOCIATED(p, head))THEN
          EXIT
       ENDIF
       o => p % old
    ENDDO
  END SUBROUTINE addi
END MODULE fork
