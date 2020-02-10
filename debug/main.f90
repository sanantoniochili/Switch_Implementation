  program main
    implicit none

    real,  dimension(2)            :: xc
    real,  dimension(2)            :: gc
    real,  dimension(2)            :: pvect

    real                           :: fc
    real                           :: glast
    integer                        :: i
    
    do i=1,2 ! initial point x0 = (0,0)
      xc(i) = 0
      pvect(i) = 1
    end do

    glast = 1.0

    call energy(xc, fc, gc)
    call minimise(xc, gc, pvect, glast)
  
    ! do i = 1,2
    !   print *,gc(i)
    ! end do    

  end program main