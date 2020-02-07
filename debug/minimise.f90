  program minimise(xc,fc,gc)
    integer, parameter :: i4  = selected_int_kind(5)
    integer, parameter :: dp  = kind(1.0d0)

  integer(i4)                                :: nvar=3
!   real(dp),             intent(inout)         :: fc
!   real(dp),             intent(inout)         :: gc(nvar)
!   real(dp),             intent(inout)         :: xc(nvar)
! !
! !  Local variables
! !

!   real(dp), dimension(:),   allocatable       :: gg

!   !
!   !  Set local number of variables
!   !
!     maxnvar = max(nvar,1_i4)

  !
  !  Allocate local memory
  !
  ! allocate(gg(maxnvar),stat=status)

! !************************
! !  Conjugate gradients  *
! !************************
!   if (lfrst.or.mod(jcyc,nupdate).eq.0) then
!     do i = 1,nvar
!       xlast(i) = xc(i)
!       glast(i) = - gsca*gc(i)
!       pvect(i) = glast(i)
!     enddo
!     lfrst = .false.
!   else
!     ggg = 0.0_dp
!     dggg = 0.0_dp
!     do i = 1,nvar
!       ggg = ggg + glast(i)*glast(i)
!       dggg = dggg + (glast(i) + gsca*gc(i))*gc(i)*gsca
!     enddo
!     gam = dggg/ggg
!     do i = 1,nvar
!       xlast(i) = xc(i)
!       glast(i) = - gsca*gc(i)
!       pvect(i) = glast(i) + gam*pvect(i)
!     enddo
!   endif
!   pnorm = sqrt(ddot(nvar,pvect,1_i4,pvect,1_i4))
!   call linmin(xc,alp,pvect,nvar,fc,okf,gg,imode)
  
!
!  Free local memory
!
  ! deallocate(xvar,stat=status)


  return
  end

