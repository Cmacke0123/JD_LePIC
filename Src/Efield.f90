subroutine calc_Efield(n,h,phi,Ei,bcnd)
!     ===================================================================
!     VERSION:         0.3
!     LAST MOD:      Sep/20
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     NOTE:    
!     -------------------------------------------------------------------
  implicit none
  integer:: ix,iy,n(3),bcnd(0:n(1)+2,0:n(2)+2)
  real(kind=8):: h(3),phi(0:n(1)+2,0:n(2)+2),Ei(2,0:n(1)+2,0:n(2)+2)
  include 'particle_info.h'

  !
  ! Interior points
  !
  !$OMP PARALLEL
  !$OMP DO
  do iy=1,n(2)+1
     do ix=1,n(1)+1

        ! Interior points
        if( bcnd(ix,iy).le.0 ) then

           ! second order correct in dx
           Ei(1,ix,iy)= -( phi(ix+1,iy)-phi(ix-1,iy) )/(2.d0*h(1))

           ! second order correct in dy
           Ei(2,ix,iy)= -( phi(ix,iy+1)-phi(ix,iy-1) )/(2.d0*h(2))

        endif

     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  !
  ! Boundary conditions
  !
  !$OMP PARALLEL
  !$OMP DO
  do iy=1,n(2)+1
     do ix=1,n(1)+1

        ! Skip everything 
        if( bcnd(ix,iy).eq.-1 ) goto 10

        !
        ! Neumann BCs
        !
        if( bcnd(ix,iy).eq.-2 ) then ! Along (Oy), LHS only
           
           ! Calculate Ex
           Ei(1,1,iy)= 0.d0

           if( iy.gt.1 .and. iy.lt.(n(2)+1) ) then
              ! Calculate Ey
              Ei(2,1,iy)= -( phi(1,iy+1)-phi(1,iy-1) )/(2.d0*h(2))
              goto 10
           endif

        endif

        !
        ! Periodic BCs
        !
        if( bcnd(ix,iy).eq.0 .or. bcnd(ix,iy).eq.-2 ) then

           if( iy.eq.1 ) then ! Along (Ox)
              ! Calculate Ey
              Ei(2,ix,1)= -( phi(ix,2)-phi(ix,n(2)) )/(2.d0*h(2))
              Ei(2,ix,0)= Ei(2,ix,n(2))
              Ei(2,ix,n(2)+1)= Ei(2,ix,1)
              Ei(2,ix,n(2)+2)= Ei(2,ix,2)
              ! Calculate Ex
              Ei(1,ix,0)= Ei(1,ix,n(2))
              Ei(1,ix,n(2)+2)= Ei(1,ix,2)
              goto 10
           endif

        endif

        !
        ! Walls
        !
        if( bcnd(ix,iy).ge.1 ) then

           ! Ex
           Ei(1,ix,iy)= 0.d0
           ! Ey
           Ei(2,ix,iy)= 0.d0
           
           ! Wall on the LHS
           if( bcnd(ix+1,iy).le.0 ) then          
              ! Ex
              Ei(1,ix,iy)= 2.d0*Ei(1,ix+1,iy) - Ei(1,ix+2,iy)   
           endif

           ! Wall on the RHS
           if( bcnd(ix-1,iy).le.0 ) then          
              ! Ex
              Ei(1,ix,iy)= 2.d0*Ei(1,ix-1,iy) - Ei(1,ix-2,iy)  
           endif

           ! Wall at the bottom
           if( bcnd(ix,iy+1).le.0 ) then          
              ! Ey
              Ei(2,ix,iy)= 2.d0*Ei(2,ix,iy+1) - Ei(2,ix,iy+2)   
           endif

           ! Wall at the top
           if( bcnd(ix,iy-1).le.0 ) then          
              ! Ey
              Ei(2,ix,iy)= 2.d0*Ei(2,ix,iy-1) - Ei(2,ix,iy-2) 
           endif

        endif

10      continue       
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  return
  
end subroutine calc_Efield
