! switch_perio = 0 : no periodic boundary conditions for y direction - NxNy=n1x*n1y,jp1=0,jp2=ny
! switch_perio = 1 : periodic boundary conditions for y direction NxNy=n1x*ny,jp1=1,jp2=ny
module mod_system
  use sparse_matrices
  implicit none
  
  real(kind=8)   ::    dx,dy,dx2i,dy2i
  integer(kind=4) :: nprt
  integer switch_perio
 
  !----------------------------  PARDISO  ---------------------------
  integer			:: indx_init=1
  integer                   :: msglvl,mtype,phase,maxfct,mnum,nrhs,idum,solver,iparm(64),error    	
  integer (kind=8)          :: pt(64)	
  real (kind=8)             :: ddum,dparm(64)
  type(sparseMatrix)	:: M_final
  
end module mod_system
