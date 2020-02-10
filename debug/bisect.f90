subroutine bisect(xc,fc,nfc)
    implicit none
    real :: args

    ! mid = left+right/2.0
    mult = fc*nfc

    if (mult < 0) then
      right = mid
    else
      left = mid      
    endif
    
end subroutine bisect