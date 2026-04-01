subroutine tridiag(n,h,u,b)
!     ==============================================================
!     VERSION:         
!     LAST MOD:   Sep/21   
!     MOD AUTHOR: Numerical Recipes   
!     COMMENTS:   Modified by G. Fubiani   
!     NOTE:       Published in Numerical Recipes (Fortran 77), 
!                 2nd edition, Sec. 2.4      
!     --------------------------------------------------------------
  implicit none
  integer:: ix,is,n(2)
  real(kind=8):: h(3),aw,ae,ac
  real(kind=8):: bet,gam(n(1)+1),b(n(1)+1,n(2)+1),u(0:n(1)+2,0:n(2)+2)
  include 'constants.h'
  include 'particle_info.h'

  ! aw*u(i-1) + ac*u(i) + ae*u(i+1)= b(i)
  aw= eps0/(h(1)*h(1))
  ae= aw
  ac= - ( aw + ae )

  is= 2
  if(flag_nmn.eq.1) is=1
  
  ! Initialisation
  if(is.eq.2) b(2,1)= b(2,1) - aw*u(1,1)
  b(n(1),1)= b(n(1),1) - ae*u(n(1)+1,1)

  bet= ac
  u(is,1)= b(is,1)/bet

  do ix=is+1,n(1)
     gam(ix)= ae/bet
     if(flag_nmn.eq.1 .and. ix.eq.2) then
        gam(ix)= 2.d0*ae/bet
     endif
     bet= ac - aw*gam(ix)
     if(bet.eq.0.) then 
        print*, 'Tridiag failed!'
        call stop_calculation
     endif
     u(ix,1)= ( b(ix,1) - aw*u(ix-1,1) )/bet
  enddo

  do ix=n(1)-1,is,-1
     u(ix,1)=u(ix,1)-gam(ix+1)*u(ix+1,1)
  enddo

  return

end subroutine tridiag
