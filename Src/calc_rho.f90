subroutine np_periodic(n,np,bcnd,ntype,nproc)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Oct/20
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  use omp_lib
  implicit none
  integer:: ix,n(3),ntype,iproc,nproc
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2)
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,ntype,nproc)

  ! Periodic boundary conditions (performed in parallel)
  !$OMP PARALLEL PRIVATE(ix,iproc)
  ! Get processor id (from 0 to nproc-1)
  iproc= omp_get_thread_num() + 1

  do ix=1,n(1)+1  ! along (Ox)
     if( bcnd(ix,1).eq.0 ) then
        np(ix,1,:,iproc)= 0.5d0*( np(ix,1,:,iproc) + &
             np(ix,n(2)+1,:,iproc) )
        np(ix,0,:,iproc)= np(ix,n(2),:,iproc)
     endif
     
     if( bcnd(ix,n(2)+1).eq.0 ) then
        np(ix,n(2)+2,:,iproc)= np(ix,2,:,iproc)
        np(ix,n(2)+1,:,iproc)= np(ix,1,:,iproc)   
     endif
  enddo
  !$OMP END PARALLEL
  
  return
end subroutine np_periodic

subroutine calc_rhs(n,np,rhs,ntype,nproc,nproc_mpi)
!     ==============================================================
!     VERSION:         0.6
!     LAST MOD:      Oct/20
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  use omp_lib
  implicit none
  include 'mpif.h'
  integer ierr
  integer:: ix,iy,ptype,n(3),ntype,iproc,nproc,nproc_mpi
  ! Particle arrays
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,ntype,nproc)
  real(kind=8):: rhs(n(1)+1,n(2)+1)
  include 'particle_info.h'

  rhs=0.d0 ! Initialization

  !$OMP PARALLEL
  !$OMP DO
  do iy=1,n(2)+1
     do ix=1,n(1)+1
        do iproc=1,nproc 
           do ptype=1,ntype
              rhs(ix,iy)= rhs(ix,iy) - &
                   charge(ptype)*np(ix,iy,ptype,iproc)
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL  
     
  if(nproc_mpi.gt.1) then
     call MPI_ALLREDUCE(MPI_IN_PLACE, rhs, (n(1)+1)*(n(2)+1), &
          MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  endif
     
 return
end subroutine calc_rhs

subroutine calc_rhs_par(n,np,rhs_par,bcnd,phi,jpl,jpr,NxNy,ntype,nproc,nproc_mpi)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Oct/20
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  use omp_lib
  implicit none
  include 'mpif.h'
  integer ierr
  integer:: ix,iy,ptype,n(3),ntype,iproc,nproc,nproc_mpi,&
       jpr,jpl,NxNy,itot
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2)
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,ntype,nproc), &
       phi(0:n(1)+2,0:n(2)+2)
  real(kind=8):: rhs_par(NxNy)
  include 'particle_info.h'

     rhs_par= 0.d0

     !$OMP PARALLEL PRIVATE(itot)
     !$OMP DO
     do iy=jpl,jpr
        do ix=0,n(1)
           if(flag_pbc.eq.0) then
              itot=ix+iy*(n(1)+1)+1
           else
              itot=ix+(iy-1)*(n(1)+1)+1
           endif
           ! Must shift by one cell index by convention
           if(bcnd(ix+1,iy+1).le.0) then 
              do iproc=1,nproc 
                 do ptype=1,ntype
                    rhs_par(itot)= rhs_par(itot) - &
                         charge(ptype)*np(ix+1,iy+1,ptype,iproc)
                 enddo
              enddo
           else ! update Dirichlet BC's
              rhs_par(itot)= phi(ix+1,iy+1)/real(nproc_mpi)
           endif
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL  
     
     if(nproc_mpi.gt.1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE, rhs_par, NxNy, &
             MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     endif

  return
end subroutine calc_rhs_par

subroutine dens_red(n,np,np_red,bcnd,ntype,nproc,nproc_mpi)
!     ==============================================================
!     VERSION:         0.4
!     LAST MOD:      Oct/16
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  use omp_lib
  implicit none
  include 'mpif.h'
  integer ierr
  integer:: ix,iy,ptype,n(3),ntype,iproc,nproc,nproc_mpi
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2)
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,ntype,nproc), &
       np_red(0:n(1)+2,0:n(2)+2,ntype)
  real(kind=8), allocatable:: np_red_tmp(:,:,:) 
  include 'particle_info.h'

  ! Initialize arrays
  np_mx(1:ntype)= 0.d0
  np_red=0.d0

  !
  ! Reduction
  !
  !$OMP PARALLEL
  !$OMP DO
  do iy=1,n(2)+1
     do ix=1,n(1)+1
        do iproc=1,nproc 
           do ptype=1,ntype
              np_red(ix,iy,ptype)= np_red(ix,iy,ptype) + &
                   np(ix,iy,ptype,iproc)
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL  

  if(nproc_mpi.gt.1) then
     allocate( np_red_tmp(0:n(1)+2,0:n(2)+2,ntype) )
     np_red_tmp= 0.d0
     call MPI_ALLREDUCE(np_red(:,:,:), np_red_tmp, (n(1)+3)*(n(2)+3)*ntype, &
          MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     np_red= np_red_tmp
     deallocate( np_red_tmp )
  endif

  do ptype=1,ntype
     !$OMP PARALLEL REDUCTION(MAX:np_mx)
     !$OMP DO
     do iy=1,n(2)+1
        do ix=1,n(1)+1
           ! Calculate maximum density
           if( bcnd(ix,iy).eq.-1 ) & ! inside simulation domain
                np_mx(ptype)= MAX(np_mx(ptype), np_red(ix,iy,ptype))
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL  
  enddo

  return
end subroutine dens_red
