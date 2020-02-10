  subroutine energy(xc,fc,gc)
    implicit none

    integer                                          :: i
    real,  dimension(2),   intent(in)                :: xc
    real,  dimension(2),   intent(inout)             :: gc
    real,                  intent(inout)             :: fc  

    ! f(x1,x2) = (x1-2)^2/4 + (x2-1)^2

    fc = (xc(1)-2)**2/4.0 + (xc(2)-1)*2

    ! grad = (x1/2-1, 2x2-2)

    gc(1) = xc(1)/2 - 1
    gc(2) = 2*xc(2) - 2
    
  end subroutine energy