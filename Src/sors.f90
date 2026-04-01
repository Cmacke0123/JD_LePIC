subroutine sor_rb(u,b,h,bcnd,res,n,omega,eps,ksor,kmg,dig,rank,nproc)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      July/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Solve poisson equation using a weighted Gauss-
!                    Seidel method (Successive Over-Relaxation).
!     NOTE:          u(x,y) is defined as u(0:nx+2,0:ny+2) where
!                    1 and n+1 are for the boundary conditions.
!                    Red-black ordering is used
!                    Residual is not exact but combine new and old
!                    values. This avoid another n^2 loop.
!     --------------------------------------------------------------
  implicit none
  include 'mpif.h'
  integer:: i,j,is,jl,jr,k,ksor,kmg,n(3),nr,dig,m,u_0,u_n2p1, &
       flag
  integer ierr,rank,nproc,status(MPI_STATUS_SIZE)
  integer:: bcnd(0:n(1)+2,0:n(2)+2),shift(0:nproc-1), &
       length(0:nproc-1)
  real(kind=8):: h(3),u(0:n(1)+2,0:n(2)+2),b(n(1)+1,n(2)+1), &
       u_pbc(0:n(1)+2,2),mr
  real(kind=8), allocatable:: u_tmp(:,:),b_tmp(:,:)
  integer, allocatable:: bcnd_tmp(:,:)
  parameter (u_0=1, u_n2p1=2)
  ! Solver constants
  real(kind=8):: omega,rij,res,res_tmp,res_s,eta,eps ! residual
  include 'particle_info.h'
  include 'constants.h'
  include 'mg.h'

  !
  ! Compute size of local block
  !
  m= n(2)/nproc
  mr= real(n(2))/real(nproc)

  flag=0
  if( (real(m)-mr).eq.0.d0 .and. (real(m/2)-mr/2.d0).lt.0 .and. &
       dig.lt.0 ) flag=1

  if( (m-mr).eq.0.d0 ) then 
     ! At least 1 node per proc
     jl= rank*m+1
     if( (rank+1).eq.nproc ) then    
        jr= (rank+1)*m+1
     else
        jr= (rank+1)*m
     endif
  else ! Less than 1 node per proc
     jl= 1
     jr= n(2)+1
  endif

  if(flag.eq.1) allocate ( u_tmp(0:n(1)+2,0:n(2)+2), &
       b_tmp(n(1)+1,n(2)+1), bcnd_tmp(0:n(1)+2,0:n(2)+2) )

  !
  ! Initialization & multi-grid parameters
  !
  res_s=0.d0
  ksor=0
  nr=10
  if( kmg.eq.1 .and. dig.eq.0 ) nr=10**6
  eta=0.6d0

  !
  ! SOR iteration
  !
  do k=1,nr

     res=0.d0

     !
     ! Red ordering based on old back values
     ! i=2 => j=1,3,5, etc.
     ! i=3 => j=2,4,6, etc.
     !

     !$OMP PARALLEL PRIVATE(is)
     !$OMP DO
     do j=jl,jr
     
        if(mod(j,2).eq.0) then
           is=1
        else
           is=2
        endif

        do i=is,n(1)+1,2
    
           if(bcnd(i,j).le.0) then
              
              ! Neumann BC's (LHS only)
              if( bcnd(i,j).eq.-2 ) then 
                 u(i,j)= (1.d0-omega)*u(i,j) + &
                      omega*(  b(i,j) - an*u(i,j+1) - & 
                      as*u(i,j-1) - 2.d0*ae*u(i+1,j)  )/ac
                 goto 100
              endif

              ! Periodic BCs (calculate from 1->n)
              if( bcnd(i,j).eq.0 ) then
                 if( j.eq.n(2)+1 ) goto 100
              endif

              ! Update u(i,j) [iter. k+1] inside domain     
              u(i,j)= (1.d0-omega)*u(i,j) + &
                   omega*(  b(i,j) - aw*u(i-1,j) - an*u(i,j+1) - & 
                   as*u(i,j-1) - ae*u(i+1,j)  )/ac

100           continue
           endif
           
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL

     if(m.lt.1) goto 110 ! Skip inter-node communications for m<1

     !
     ! Send/receive ghost nodes to/from neighbour processes
     !
     if(MOD(rank,2).eq.1) then ! rank is odd
        ! Send ghost value u(:,j=rank*m+1) to rank-1
        call MPI_SEND(u(0,rank*m+1),n(1)+3, MPI_REAL8, rank-1, &
             0, MPI_COMM_WORLD, ierr)
        ! Received ghost value from rank-1 stored in u(:,rank*m-1:rank*m)
        call MPI_RECV(u(0,rank*m-1),(n(1)+3)*2, MPI_REAL8, rank-1, &
             0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
     
        if(rank.lt.nproc-1) then
           ! Send ghost value u(:,(rank+1)*m-1:(rank+1)*m) to rank+1
           call MPI_SEND(u(0,(rank+1)*m-1),(n(1)+3)*2, MPI_REAL8, rank+1, &
                0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
           ! Received ghost value from rank+1 stored in u(:,j=(rank+1)*m+1)
           call MPI_RECV(u(0,(rank+1)*m+1),n(1)+3, MPI_REAL8, rank+1, &
                0, MPI_COMM_WORLD, status, ierr)
        endif
     else ! Rank is even 
        if(rank.gt.0) then
           call MPI_RECV(u(0,rank*m-1),(n(1)+3)*2, MPI_REAL8, rank-1, & 
                0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
           call MPI_SEND(u(0,rank*m+1),n(1)+3, MPI_REAL8, rank-1, &
                0, MPI_COMM_WORLD, ierr)
        endif

        if(rank.lt.nproc-1) then
           call MPI_RECV(u(0,(rank+1)*m+1),n(1)+3, MPI_REAL8, rank+1, &
                0, MPI_COMM_WORLD, status, ierr)
           call MPI_SEND(u(0,(rank+1)*m-1),(n(1)+3)*2, MPI_REAL8, rank+1, &
                0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
        endif
     endif
     
110  continue

     !
     ! Send/receive periodic boundary values (update j=0 and n+1)
     !
     if(flag_pbc.eq.0) goto 120

     if( nproc.gt.1 .and. (real(m)-mr).eq.0.d0 ) then
        if(rank.eq.0) then
           ! Send periodic BC u(:,1) to nprocs-1
           call MPI_SEND(u(0,1),n(1)+3, MPI_REAL8, nproc-1, &
                0, MPI_COMM_WORLD, ierr)
        endif
        if(rank.eq.(nproc-1)) then
           ! Received periodic value from rank #0 : u(:,n(2)+1)=u(:,1)
           call MPI_RECV(u_pbc(0,u_n2p1),n(1)+3, MPI_REAL8, 0, &
                0, MPI_COMM_WORLD, status, ierr)
           ! Send periodic BC u(:,1) to rank #0
           call MPI_SEND(u(0,n(2)),n(1)+3, MPI_REAL8, 0, &
                1, MPI_COMM_WORLD, ierr)
        endif
        if(rank.eq.0) then
           ! Received periodic value from nproc-1 : u(:,0)=u(:,n(2))
           call MPI_RECV(u_pbc(0,u_0),n(1)+3, MPI_REAL8, nproc-1, &
                1, MPI_COMM_WORLD, status, ierr)
        endif
     endif
       
     do i=0,n(1)+2
        if( bcnd(i,1).le.0 ) then
           if( nproc.gt.1 .and. (real(m)-mr).eq.0.d0 ) then
              if(rank.eq.0) u(i,0)= u_pbc(i,u_0)
           else
              u(i,0)= u(i,n(2))
           endif
        endif
        if( bcnd(i,n(2)+1).le.0 ) then
           if( nproc.gt.1 .and. (real(m)-mr).eq.0.d0 ) then
              if(rank.eq.nproc-1) u(i,n(2)+1)= u_pbc(i,u_n2p1)
           else
              u(i,n(2)+1)= u(i,1)
           endif
        endif
     enddo

120 continue

     !
     ! Black ordering based on new red values
     ! i=2 => j=2,4,6, etc.
     ! i=3 => j=1,3,5, etc.
     !     

     !$OMP PARALLEL PRIVATE(is)
     !$OMP DO
     do j=jl,jr
        
        if(mod(j,2).eq.0) then
           is=2
        else
           is=1
        endif
        
        do i=is,n(1)+1,2
        
           if(bcnd(i,j).le.0) then
              
              ! Neumann BC's (LHS only)
              if( bcnd(i,j).eq.-2 ) then 
                 u(i,j)= (1.d0-omega)*u(i,j) + &
                      omega*(  b(i,j) - an*u(i,j+1) - & 
                      as*u(i,j-1) - 2.d0*ae*u(i+1,j)  )/ac
                 goto 130
              endif

              ! Periodic BCs (calculate from 1->n)
              if( bcnd(i,j).eq.0 ) then
                 if( j.eq.n(2)+1 ) goto 130
              endif
              
              ! Update u(i,j) [iter. k+1] inside domain     
              u(i,j)= (1.d0-omega)*u(i,j) + &
                   omega*(  b(i,j) - aw*u(i-1,j) - an*u(i,j+1) - & 
                   as*u(i,j-1) - ae*u(i+1,j)  )/ac
              
130           continue
           endif
        
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL

     if(m.lt.1) goto 140 ! Skip inter-node communications for m<1

     !
     ! Send/receive ghost nodes to/from neighbour processes 
     !
     if(MOD(rank,2).eq.1) then ! rank is odd
        ! Send ghost value u(:,j=rank*m+1) to rank-1
        call MPI_SEND(u(0,rank*m+1),n(1)+3, MPI_REAL8, rank-1, &
             0, MPI_COMM_WORLD, ierr)
        ! Received ghost value from rank-1 stored in u(:,rank*m-1:rank*m)
        call MPI_RECV(u(0,rank*m-1),(n(1)+3)*2, MPI_REAL8, rank-1, &
             0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
        
        if(rank.lt.nproc-1) then
           ! Send ghost value u(:,(rank+1)*m-1:(rank+1)*m) to rank+1
           call MPI_SEND(u(0,(rank+1)*m-1),(n(1)+3)*2, MPI_REAL8, rank+1, &
                0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
           ! Received ghost value from rank+1 stored in u(:,j=(rank+1)*m+1)
           call MPI_RECV(u(0,(rank+1)*m+1),n(1)+3, MPI_REAL8, rank+1, &
                0, MPI_COMM_WORLD, status, ierr)
        endif
     else ! Rank is even 
        if(rank.gt.0) then
           call MPI_RECV(u(0,rank*m-1),(n(1)+3)*2, MPI_REAL8, rank-1, &
                0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
           call MPI_SEND(u(0,rank*m+1),n(1)+3, MPI_REAL8, rank-1, &
                0, MPI_COMM_WORLD, ierr)
        endif

        if(rank.lt.nproc-1) then
           call MPI_RECV(u(0,(rank+1)*m+1),n(1)+3, MPI_REAL8, rank+1, &
                0, MPI_COMM_WORLD, status, ierr)
           call MPI_SEND(u(0,(rank+1)*m-1),(n(1)+3)*2, MPI_REAL8, rank+1, &
                0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
        endif
     endif

140  continue

     !
     ! Send/receive periodic boundary values (update j=0 and n+1)
     !
     if(flag_pbc.eq.0) goto 150

     if( nproc.gt.1 .and. (real(m)-mr).eq.0.d0 ) then
        if(rank.eq.0) then
           ! Send periodic BC u(:,1) to nprocs-1
           call MPI_SEND(u(0,1),n(1)+3, MPI_REAL8, nproc-1, &
                0, MPI_COMM_WORLD, ierr)
        endif
        if(rank.eq.(nproc-1)) then
           ! Received periodic value from rank #0 : u(:,n(2)+1)=u(:,1)
           call MPI_RECV(u_pbc(0,u_n2p1),n(1)+3, MPI_REAL8, 0, &
                0, MPI_COMM_WORLD, status, ierr)
           ! Send periodic BC u(:,n(2)) to rank #0
           call MPI_SEND(u(0,n(2)),n(1)+3, MPI_REAL8, 0, &
                1, MPI_COMM_WORLD, ierr)
        endif
        if(rank.eq.0) then
           ! Received periodic value from nproc-1 : u(:,0)=u(:,n(2))
           call MPI_RECV(u_pbc(0,u_0),n(1)+3, MPI_REAL8, nproc-1, &
                1, MPI_COMM_WORLD, status, ierr)
        endif
     endif
       
     do i=0,n(1)+2
        if( bcnd(i,1).le.0 ) then
           if( nproc.gt.1 .and. (real(m)-mr).eq.0.d0 ) then
              if(rank.eq.0) u(i,0)= u_pbc(i,u_0)
           else
              u(i,0)= u(i,n(2))
           endif
        endif
        if( bcnd(i,n(2)+1).le.0 ) then
           if( nproc.gt.1 .and. (real(m)-mr).eq.0.d0 ) then
              if(rank.eq.nproc-1) u(i,n(2)+1)= u_pbc(i,u_n2p1)
           else
              u(i,n(2)+1)= u(i,1)
           endif
        endif
     enddo

150  continue

     !
     ! Calculate residual
     !
     
     res=0.d0
     !$OMP PARALLEL PRIVATE(rij) REDUCTION(+:res)
     !$OMP DO
     do j=jl,jr
        do i=1,n(1)+1
           
           rij=0.d0
           
              ! Interior points
           if(bcnd(i,j).eq.-1) then
              rij= b(i,j) - ( aw*u(i-1,j) + an*u(i,j+1) +  &
                   as*u(i,j-1) + ae*u(i+1,j) + ac*u(i,j) )
           endif

           res= res + ABS(rij)
           
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL

     !
     ! Calculate max residual between all processes
     !
     if( nproc.gt.1 .and. (real(m)-mr).eq.0.d0 ) then
        call MPI_ALLREDUCE(res, res_tmp, 1, MPI_REAL8, MPI_SUM, &
             MPI_COMM_WORLD, ierr)
        res= res_tmp
     endif

     !
     ! Convergence test
     !
     ksor= ksor + 1 

     if( dig.eq.0 ) then
        if( res.le.eps ) exit
     else
        if( res_s.gt.0.d0 .and. (res/res_s).le.eta ) exit
     endif

     res_s=res

  enddo ! End iterative loop 

  !
  ! Concatenate arrays u, bcnd, b when m=1 and move toward coarser grids
  !
  if(flag.eq.1) then  

     ! u() & bcnd()
     shift(0)=0
     do i=0,nproc-1
        length(i)=m
        if(i.ge.1) shift(i)=i*m+1
     enddo
     length(0)=length(0)+1
     length(nproc-1)=length(nproc-1)+2
     
     length= length*(n(1)+3)
     shift= shift*(n(1)+3)

     jl=rank*m+1
     if(rank.eq.0) jl=jl-1
     jr=(rank+1)*m
     if(rank.eq.nproc-1) jr=jr+1
     
     call MPI_ALLGATHERV(u(0:n(1)+2,jl:jr),length(rank),MPI_REAL8,u_tmp, &
          length,shift,MPI_REAL8,MPI_COMM_WORLD,ierr)
     u= u_tmp

     call MPI_ALLGATHERV(bcnd(0:n(1)+2,jl:jr),length(rank),MPI_INTEGER,bcnd_tmp, &
          length,shift,MPI_INTEGER,MPI_COMM_WORLD,ierr)
     bcnd= bcnd_tmp

     ! b()
     shift(0)=0
     do i=0,nproc-1
        length(i)=m
        if(i.ge.1) shift(i)=i*m
     enddo
     length(nproc-1)=length(nproc-1)+1
     
     length= length*(n(1)+1)
     shift= shift*(n(1)+1)

     jl=rank*m+1
     jr=(rank+1)*m
     if(rank.eq.nproc-1) jr=jr+1

     call MPI_ALLGATHERV(b(1:n(1)+1,jl:jr),length(rank),MPI_REAL8,b_tmp, &
          length,shift,MPI_REAL8,MPI_COMM_WORLD,ierr)
     b= b_tmp

  endif

  return
end subroutine sor_rb
