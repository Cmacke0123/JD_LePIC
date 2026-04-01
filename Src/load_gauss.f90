subroutine load_gauss(vx,vy,vt,rnd)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:       AUG/07
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Gaussian temperature loading
!     NOTE:          
!     --------------------------------------------------------------
  implicit none
  real(kind=8):: theta,vp,vx,vy,vt,rnd(2)
  include 'constants.h'

  ! Gaussian loading
  vp = vt*dsqrt( -dlog(1-rnd(1)) )
  
  ! Update transverse normalized momentum
  theta = 2.d0*pi*rnd(2)
  
  ! Calculate new normalized velocity
  vx = vp*dcos(theta)
  vy = vp*dsin(theta)
  
  return
end subroutine load_gauss
