subroutine directsolver(nx,ny,h,jp1,jp2,NxNy,phi,rho,bcnd,flag_pbc,flag,mpi_rank)
  !     ==============================================================
  !      This routine calls Pardiso
  !     --------------------------------------------------------------
  use mod_directsolver
  implicit none
  integer:: flag,flag_pbc,mpi_rank
  integer:: nx,ny,jp1,jp2,NxNy
  integer:: bcnd(0:nx+2,0:ny+2)
  real(kind=8):: h(3),phi(0:nx+2,0:ny+2),rho(NxNy)
  include 'constants.h'

  ! Flag periodic conditions
  switch_perio=0
  if(flag_pbc.eq.1) switch_perio=1

  dx= h(1)
  dy= h(2)
  dx2i=eps0/(dx*dx)
  dy2i=eps0/(dy*dy)

  !
  ! Call PDE solver
  !
  if(flag.eq.1) call initializepardiso(nx,ny,NxNy,jp1,jp2,bcnd,mpi_rank)  
  if(flag.eq.2) call runpardiso(nx,ny,NxNy,jp1,jp2,phi,rho)

end subroutine directsolver
