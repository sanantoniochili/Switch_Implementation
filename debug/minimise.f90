  
  ! Common Conjugate Gradient

  subroutine minimise(xc,fc,gc,pvect,glast)
    implicit none

    integer                                          :: i
    real,  dimension(2),   intent(inout)             :: xc
    real,  dimension(2),   intent(in)                :: gc
    real,                  intent(in)                :: fc
    real,  dimension(2),   intent(inout)             :: pvect
    real,  dimension(2)                              :: plast
    real,  dimension(2)                              :: rk
    real,                  intent(inout)             :: glast
    real                                             :: gn
    real                                             :: beta

  call gnorm(gc, gn)
  do i = 1,2
    rk(i) = - gc(i)/gn ! residual
  end do

  do i = 1,2
    xc(i) = xc(i) + alpha*pvect(i)
  end do

  beta = gn/glast ! ||cur_grad|| / ||prev_grad||

  ! Calculate new direction
  do i = 1,2
    pvect(i) = rk(i) + beta*pvect(i)
  end do

  print*, glast
  glast = gn ! save last gradient norm 
  end subroutine minimise


! Calculate gn = ||grad||

  subroutine gnorm(gc,gn)
    implicit none

    integer                                             :: i
    real, dimension(2), intent(in)                      :: gc
    real                                                :: grad
    real,               intent(out)                     :: gn

    grad = 0
    do i = 1,2
      grad = gc(i)**2
    end do

    gn = sqrt(grad)
  end subroutine gnorm