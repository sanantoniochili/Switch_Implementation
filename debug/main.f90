  program main
    real, dimension(3) :: xc
    real, dimension(3) :: gc

    call minimise(1,2,3)
    write(*,*) "ok"
  end