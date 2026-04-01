subroutine iterative_solver(u,b,bcnd,h,n,ncycl,eps,k,ktot,res,ng,rank,nproc)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      July/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Solve diffusion equation using a V-shaped
!                    multi-grid method. SOR algorithm is used for 
!                    relaxation.
!     NOTE:          u(x,y,t) is defined as u(0:nx+2,0:ny+2) where
!                    1 and n+1 are for the boundary conditions.
!     --------------------------------------------------------------
  implicit none
  integer:: k,n(3),ng,ncycl,bcnd(0:n(1)+2,0:n(2)+2),rank,nproc
  real(kind=8):: h(3),u(0:n(1)+2,0:n(2)+2),b(n(1)+1,n(2)+1)
  real(kind=8):: omega,eps,ktot,u0,alpha
  real(kind=8):: res,res_s ! residuals
  logical:: converged
  include 'constants.h'
  include 'mg.h'

  ! Parameters for convergence test
  u0= 1.d0 ! 1V
  alpha= ABS((n(1)+1)*(n(2)+1)*ac)*u0

  ! Initialize variables and arrays
  ktot=0
  res_s=0
  omega= 1.4d0

  ! Calculate Solution
  k=1
  converged=.FALSE.
  do while (.NOT.converged)
   
     call mg(u,b,bcnd,h,res,n,omega,ng,alpha*eps,k,ktot,rank,nproc)
     res= res/alpha

     ! Some warnings 
     if(k.eq.1) then 
        res_s=res
     else
        if(ABS(res/res_s).gt.10.d0) then
           if(rank.eq.0) then 
              print*, ' '
              print*, 'Warning: PDE solver is diverging'
              print*, 'Abort calculation ...'
              print*, 'Final sol: k=',k,'ratio=',res/res_s
              call stop_calculation
           endif
        endif
     endif

     ! Test for convergence
     if( res.le.eps .or. k.eq.ncycl ) converged=.TRUE.

     ! Print info. on screen
     if(MOD(k,1000).eq.0) then
        if(rank.eq.0) then
           print*, 'k=',k,'res=',res
           print*, 'Equivalent SOR iterations:',nint(ktot)
        endif
     endif

     k=k+1

  enddo ! end dowhile

  return
end subroutine iterative_solver
