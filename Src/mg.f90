subroutine mg(u,b,bcnd,h,res,n,omega,ng,eps,k,ktot,rank,nproc)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      July/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Solve poisson equation using a V-shaped
!                    multi-grid method. SOR algorithm is used for 
!                    relaxation.
!     NOTE:          1) u(x,y,t) is defined as u(0:nx+2,0:ny+2) where
!                    1 and n+1 are for the boundary conditions.
!                    2) Graphical representation of MG scheme:
!                     
!                             Fine grid (n=8)
!                    |---|---|---|---|---|---|---|---|
!                    1   2   3   4   5   6   7   8   9
!                    1(BC)                          n+1(BC)  
!                             
!                          First level coarsed grid
!                    |-------|-------|-------|-------|
!                    1       2       3       4       5
!                    1(BC)                         n/2+1(BC)  
!
!     --------------------------------------------------------------
  implicit none
  integer:: i,k,ksor,ig,n(2),ng,iter(ng),il,rank,nproc
  real(kind=8):: h(3),u(0:n(1)+2,0:n(2)+2),b(n(1)+1,n(2)+1)
  ! Define multigrid arrays
  real(kind=8):: r2(n(1)/2+1,n(2)/2+1),e2(0:n(1)/2+2,0:n(2)/2+2), &
       r4(n(1)/4+1,n(2)/4+1),e4(0:n(1)/4+2,0:n(2)/4+2), &
       r8(n(1)/8+1,n(2)/8+1),e8(0:n(1)/8+2,0:n(2)/8+2), &
       r16(n(1)/16+1,n(2)/16+1),e16(0:n(1)/16+2,0:n(2)/16+2), &
       r32(n(1)/32+1,n(2)/32+1),e32(0:n(1)/32+2,0:n(2)/32+2), &
       r64(n(1)/64+1,n(2)/64+1),e64(0:n(1)/64+2,0:n(2)/64+2), &
       r128(n(1)/128+1,n(2)/128+1),e128(0:n(1)/128+2,0:n(2)/128+2), &
       r256(n(1)/256+1,n(2)/256+1),e256(0:n(1)/256+2,0:n(2)/256+2), &
       r512(n(1)/512+1,n(2)/512+1),e512(0:n(1)/512+2,0:n(2)/512+2), &
       r1024(n(1)/1024+1,n(2)/1024+1),e1024(0:n(1)/1024+2,0:n(2)/1024+2), &
       r2048(n(1)/2048+1,n(2)/2048+1),e2048(0:n(1)/2048+2,0:n(2)/2048+2)
  integer:: bcnd(0:n(1)+2,0:n(2)+2),bcnd2(0:n(1)/2+2,0:n(2)/2+2), &
       bcnd4(0:n(1)/4+2,0:n(2)/4+2),bcnd8(0:n(1)/8+2,0:n(2)/8+2), &
       bcnd16(0:n(1)/16+2,0:n(2)/16+2),bcnd32(0:n(1)/32+2,0:n(2)/32+2), &
       bcnd64(0:n(1)/64+2,0:n(2)/64+2),bcnd128(0:n(1)/128+2,0:n(2)/128+2), &
       bcnd256(0:n(1)/256+2,0:n(2)/256+2),bcnd512(0:n(1)/512+2,0:n(2)/512+2), &
       bcnd1024(0:n(1)/1024+2,0:n(2)/1024+2),bcnd2048(0:n(1)/2048+2,0:n(2)/2048+2)
  integer:: n2(2),n4(2),n8(2),n16(2),n32(2),n64(2),n128(2),n256(2), &
       n512(2),n1024(2),n2048(2)
  real(kind=8):: h2(3),h4(3),h8(3),h16(3),h32(3),h64(3),h128(3),h256(3), &
       h512(3),h1024(3),h2048(3)
  real(kind=8):: ktot,eps,omega
  real(kind=8):: res ! residual


  ! Initialize some arrays
  do i=1,ng
     iter(i)=0
  enddo

  ! Define multi-grid arrays
  n2=n/2
  h2=2.d0*h
  n4=n/4
  h4=4.d0*h
  n8=n/8
  h8=8.d0*h
  n16=n/16
  h16=16.d0*h
  n32=n/32
  h32=32.d0*h
  n64=n/64
  h64=64.d0*h
  n128=n/128
  h128=128.d0*h
  n256=n/256
  h256=256.d0*h
  n512=n/512
  h512=512.d0*h
  n1024=n/1024
  h1024=1024.d0*h
  n2048=n/2048
  h2048=2048.d0*h

  !
  ! Start iteration 
  !
  do ig=1,ng*2-1 ! 2*ng-1 for V shape iteration
     
     ! Level #1
     il=1
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(u,b,h,bcnd,res,n,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(u,b,h,bcnd,bcnd2,r2,e2,n,rank,nproc)
     endif

     ! Level #2
     il=2
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then
        call sor_rb(e2,r2,h2,bcnd2,res,n2,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e2,r2,h2,bcnd2,bcnd4,r4,e4,n2,rank,nproc)
        if( ig.gt.ng ) call prolongation(u,e2,bcnd,n,rank,nproc)
     endif

     ! Level #3
     il=3
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then
        call sor_rb(e4,r4,h4,bcnd4,res,n4,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e4,r4,h4,bcnd4,bcnd8,r8,e8,n4,rank,nproc)
        if( ig.gt.ng ) call prolongation(e2,e4,bcnd2,n2,rank,nproc)
     endif

     ! Level #4
     il=4
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e8,r8,h8,bcnd8,res,n8,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e8,r8,h8,bcnd8,bcnd16,r16,e16,n8,rank,nproc)
        if( ig.gt.ng ) call prolongation(e4,e8,bcnd4,n4,rank,nproc)
     endif

     ! Level #5
     il=5
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e16,r16,h16,bcnd16,res,n16,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e16,r16,h16,bcnd16,bcnd32,r32,e32,n16,rank,nproc)
        if( ig.gt.ng ) call prolongation(e8,e16,bcnd8,n8,rank,nproc)
     endif

     ! Level #6
     il=6
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e32,r32,h32,bcnd32,res,n32,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e32,r32,h32,bcnd32,bcnd64,r64,e64,n32,rank,nproc)
        if( ig.gt.ng ) call prolongation(e16,e32,bcnd16,n16,rank,nproc)
     endif

     ! Level #7
     il=7
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e64,r64,h64,bcnd64,res,n64,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e64,r64,h64,bcnd64,bcnd128,r128,e128,n64,rank,nproc)
        if( ig.gt.ng ) call prolongation(e32,e64,bcnd32,n32,rank,nproc)
     endif

     ! Level #8
     il=8
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e128,r128,h128,bcnd128,res,n128,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e128,r128,h128,bcnd128,bcnd256,r256,e256,n128,rank,nproc)
        if( ig.gt.ng ) call prolongation(e64,e128,bcnd64,n64,rank,nproc)
     endif

     ! Level #9
     il=9
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e256,r256,h256,bcnd256,res,n256,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e256,r256,h256,bcnd256,bcnd512,r512,e512,n256,rank,nproc)
        if( ig.gt.ng ) call prolongation(e128,e256,bcnd128,n128,rank,nproc)
     endif

     ! Level #10
     il=10
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e512,r512,h512,bcnd512,res,n512,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e512,r512,h512,bcnd512,bcnd1024,r1024,e1024,n512,rank,nproc)
        if( ig.gt.ng ) call prolongation(e256,e512,bcnd256,n256,rank,nproc)
     endif

     ! Level #11
     il=11
     if( (ig.lt.ng.and.ig.eq.il) .or. &
          (ig.gt.ng.and.ig.eq.(2*ng-il)) ) then 
        call sor_rb(e1024,r1024,h1024,bcnd1024,res,n1024,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(il)= iter(il) + ksor
        if( ig.lt.ng ) call restriction(e1024,r1024,h1024,bcnd1024,bcnd2048,r2048,e2048,n1024,rank,nproc)
        if( ig.gt.ng ) call prolongation(e512,e1024,bcnd512,n512,rank,nproc)
     endif

     ! Level max.
     if( ig.eq.ng ) then
        if(ng.eq.1) call sor_rb(u,b,h,bcnd,res,n,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.2) call sor_rb(e2,r2,h2,bcnd2,res,n2,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.3) call sor_rb(e4,r4,h4,bcnd4,res,n4,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.4) call sor_rb(e8,r8,h8,bcnd8,res,n8,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.5) call sor_rb(e16,r16,h16,bcnd16,res,n16,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.6) call sor_rb(e32,r32,h32,bcnd32,res,n32,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.7) call sor_rb(e64,r64,h64,bcnd64,res,n64,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.8) call sor_rb(e128,r128,h128,bcnd128,res,n128,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.9) call sor_rb(e256,r256,h256,bcnd256,res,n256,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.10) call sor_rb(e512,r512,h512,bcnd512,res,n512,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.11) call sor_rb(e1024,r1024,h1024,bcnd1024,res,n1024,omega,eps,ksor,k,ig-ng,rank,nproc)
        if(ng.eq.12) call sor_rb(e2048,r2048,h2048,bcnd2048,res,n2048,omega,eps,ksor,k,ig-ng,rank,nproc)
        iter(ig)= iter(ig) + ksor

        if(ng.eq.2) call prolongation(u,e2,bcnd,n,rank,nproc)
        if(ng.eq.3) call prolongation(e2,e4,bcnd2,n2,rank,nproc)
        if(ng.eq.4) call prolongation(e4,e8,bcnd4,n4,rank,nproc)
        if(ng.eq.5) call prolongation(e8,e16,bcnd8,n8,rank,nproc)
        if(ng.eq.6) call prolongation(e16,e32,bcnd16,n16,rank,nproc)
        if(ng.eq.7) call prolongation(e32,e64,bcnd32,n32,rank,nproc)
        if(ng.eq.8) call prolongation(e64,e128,bcnd64,n64,rank,nproc)
        if(ng.eq.9) call prolongation(e128,e256,bcnd128,n128,rank,nproc)
        if(ng.eq.10) call prolongation(e256,e512,bcnd256,n256,rank,nproc)
        if(ng.eq.11) call prolongation(e512,e1024,bcnd512,n512,rank,nproc)
        if(ng.eq.12) call prolongation(e1024,e2048,bcnd1024,n1024,rank,nproc)
     endif
     
  enddo
  
  do i=1,ng
     ktot= ktot + real(iter(i))/2.d0**(2*(i-1))
  enddo
  
  return
end subroutine mg


subroutine restriction(u,b,h,bcnd,bcnd2h,r2h,e2h,n,rank,nproc)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      July/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      h -> 2h grid mapping using a 9 point bilinear 
!                    interpolation
!     NOTE:          1) We do not need to compute 1 and n/2+1
!                    2) For the stencyl [a() array] restriction, see. 
!                       A Multigrid Tutorial, 2nd Ed., SIAM, p. 130-131
!                       Note that their domain range from 0 to n 
!                       (we use 1 to n+1) consequently we must at 1 to 
!                       each indexes in the book in order to properly 
!                       convert to our conventions. 
!
!                          W   E   W   E   W   E   W 
!                      W   E   W   E   W   E   W   E
!                    |-*-|-*-|-*-|-*-|-*-|-*-|-*-|-*-|
!                    1   2   3   4   5   6   7   8   9
!                    1(BC)      Fine grid           n+1(BC)  
!                             
!                        W       E       W       E
!                    |---*---|---*---|---*---|---*---|
!                    1       2       3       4       5
!                    1(BC)     Coarsed grid       n/2+1(BC)  
! 
!                    Legend= W==aw, E==ae, ac= -( aw + ae )
!                            aw==a^i-1/2, ae==a^i+1/2
!
!     ---------------------------------------------------------------
  implicit none
  include 'mpif.h'
  integer:: i,j,ir,jr,jm,jp,n(2),rank,nproc,m
  integer ierr,status(MPI_STATUS_SIZE)
  real(kind=8):: h(3),u(0:n(1)+2,0:n(2)+2),b(n(1)+1,n(2)+1), &
       e2h(0:n(1)/2+2,0:n(2)/2+2),r2h(n(1)/2+1,n(2)/2+1), &
       r_pbc(n(1)+1),mr
  integer:: bcnd(0:n(1)+2,0:n(2)+2),bcnd2h(0:n(1)/2+2,0:n(2)/2+2)
  real(kind=8):: re,r1,r2,r3,r4,r5,r6,r7,r8,r9
  include 'particle_info.h'

  !
  ! Compute size of local block
  !
  m= (n(2)/2)/nproc
  mr= (real(n(2))/2.d0)/real(nproc)

  if( (real(m)-mr).eq.0.d0 ) then 
     ! At least 1 node per proc  
     jm= rank*m+1
     if( (rank+1).eq.nproc ) then    
        jp= (rank+1)*m+1
     else
        jp= (rank+1)*m
     endif
  else ! Less than 1 node per proc
     jm= 1
     jp= n(2)/2+1
  endif

  !
  ! Calculate redidual @ j=n(2) for periodic BCs
  !
  if(flag_pbc.eq.0) goto 90

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ir)
  !$OMP DO
  do i=1,n(1)/2+1
     ir=2*i-1
     if( bcnd(ir,n(2)+1).le.0 ) then
        if(ir.ne.n(1)+1) call getres(u,b,h,ir+1,n(2),n,r_pbc(ir+1),flag_nmn)
        call getres(u,b,h,ir,n(2),n,r_pbc(ir),flag_nmn)
        if(ir.ne.1) call getres(u,b,h,ir-1,n(2),n,r_pbc(ir-1),flag_nmn)
     endif
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  if( nproc.gt.1 .and. (real(m)-mr).eq.0.d0 ) then
     if(rank.eq.(nproc-1)) then
        ! Send periodic BC r(:,n(2)) to rank #0
        call MPI_SEND(r_pbc,n(1)+1, MPI_REAL8, 0, &
             0, MPI_COMM_WORLD, ierr)
     endif
     if(rank.eq.0) then
        ! Received periodic value from nproc-1 : r(:,0)=r(:,n(2))
        call MPI_RECV(r_pbc,n(1)+1, MPI_REAL8, nproc-1, &
             0, MPI_COMM_WORLD, status, ierr)
     endif
  endif

90 continue

  !
  ! Sweep the coarse grid (bilinear 9 point interpolation)
  !

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ir,jr,r1,r2,r3,r4,&
  !$OMP  r5,r6,r7,r8,r9,re)
  !$OMP DO
  do j=jm,jp
     do i=1,n(1)/2+1
        
        ! Address of the fine grid point
        ir=2*i-1
        jr=2*j-1

        ! Update bcnd() array for higher grid level
        bcnd2h(i,j)=bcnd(ir,jr)

        if(ir.eq.1) bcnd2h(0,j)=bcnd(0,jr)
        if(ir.eq.n(1)+1) bcnd2h(n(1)/2+2,j)=bcnd(n(1)+2,jr)
        if(jr.eq.1) bcnd2h(i,0)=bcnd(ir,0)
        if(jr.eq.n(2)+1) bcnd2h(i,n(2)/2+2)=bcnd(ir,n(2)+2)

        ! Computes the residuals around fine grid point
        if(bcnd2h(i,j).le.0) then

           ! Initialization
           r1=0.d0
           r2=0.d0
           r3=0.d0
           r4=0.d0
           r5=0.d0
           r6=0.d0
           r7=0.d0
           r8=0.d0
           r9=0.d0

           call getres(u,b,h,ir,jr,n,r9,flag_nmn)

           if( ir.ne.(n(1)+1) .and. jr.ne.(n(2)+1) ) then
              if(bcnd(ir+1,jr+1).le.0) &
                   call getres(u,b,h,ir+1,jr+1,n,r1,flag_nmn)
           endif

           if( ir.ne.(n(1)+1) ) then
              if(bcnd(ir+1,jr).le.0) &
                   call getres(u,b,h,ir+1,jr,n,r2,flag_nmn)
           endif

           if(bcnd(ir+1,jr-1).le.0) then
              if( ir.ne.(n(1)+1) .and. jr.ne.1 ) then
                 call getres(u,b,h,ir+1,jr-1,n,r3,flag_nmn)
              else
                 if( jr.eq.1 ) r3= r_pbc(ir+1)
              endif
           endif

           if(bcnd(ir,jr-1).le.0) then
              if( jr.ne.1 ) then
                 call getres(u,b,h,ir,jr-1,n,r4,flag_nmn)
              else
                 r4= r_pbc(ir)
              endif
           endif

           if(bcnd(ir-1,jr-1).le.0 .and. ir.ne.1 ) then
              if(jr.ne.1) then
                 call getres(u,b,h,ir-1,jr-1,n,r5,flag_nmn)
              else
                 r5= r_pbc(ir-1)
              endif
           endif

           if(bcnd(ir-1,jr).le.0) then
              if( ir.ne.1 ) call getres(u,b,h,ir-1,jr,n,r6,flag_nmn)
           endif

           if(bcnd(ir-1,jr+1).le.0) then
              if( ir.ne.1 .and. jr.ne.(n(2)+1) ) &
                   call getres(u,b,h,ir-1,jr+1,n,r7,flag_nmn)
           endif

           if(bcnd(ir,jr+1).le.0) then
              if( jr.ne.(n(2)+1) ) call getres(u,b,h,ir,jr+1,n,r8,flag_nmn)
           endif

           ! Compute residual 
           re=(4.0d0*r9+2.0d0*(r2+r4+r6+r8)+r1+r3+r5+r7)/16.0d0
           r2h(i,j)=re

        else

           r2h(i,j)=0.d0
           
        endif
   
        ! Zero initial solution on coarse grid
        e2h(i,j)=0.0d0
        
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  ! Zero the remaining unknowns on the coarse grid
  if( (real(m)-mr).eq.0.d0 ) then
     e2h(:,rank*m)=0.d0
     e2h(:,(rank+1)*m+1:(rank+1)*m+2)=0.d0
     
     e2h(0,rank*m:(rank+1)*m+2)=0.d0
     e2h(n(1)/2+2,rank*m:(rank+1)*m+2)=0.d0
  else
     e2h(:,0)=0.d0
     e2h(:,n(2)/2+2)=0.d0
     
     e2h(0,:)=0.d0
     e2h(n(1)/2+2,:)=0.d0
  endif

  if(m.lt.1) goto 100 ! Skip inter-node communications for m<1

  !
  ! Send ghost nodes for r2h to processes of higher rank 
  !
  if(MOD(rank,2).eq.1) then ! rank is odd
     ! Received ghost value from rank-1 stored in r2h(:,j=rank*m)
     call MPI_RECV(r2h(1,rank*m),n(1)/2+1, MPI_REAL8, rank-1, &
          0, MPI_COMM_WORLD, status, ierr)
     if(rank.lt.nproc-1) then
        ! Send ghost value r2h(:,j=(rank+1)*m) to rank+1
        call MPI_SEND(r2h(1,(rank+1)*m),n(1)/2+1, MPI_REAL8, rank+1, &
             0, MPI_COMM_WORLD, ierr)
     endif
  else ! Rank is even 
     if(rank.gt.0) then
        call MPI_RECV(r2h(1,rank*m),n(1)/2+1, MPI_REAL8, rank-1, &
             0, MPI_COMM_WORLD, status, ierr)
     endif
     if(rank.lt.nproc-1) then
        call MPI_SEND(r2h(1,(rank+1)*m),n(1)/2+1, MPI_REAL8, rank+1, &
                0, MPI_COMM_WORLD, ierr)
     endif
  endif

  !
  ! Send/receive ghost nodes for bcnd2h to/from neighbour processes
  !
  if(MOD(rank,2).eq.1) then ! rank is odd
     ! Send ghost value bcnd2h(:,j=rank*m+1) to rank-1
     call MPI_SEND(bcnd2h(0,rank*m+1),n(1)/2+3, MPI_INTEGER, rank-1, &
          0, MPI_COMM_WORLD, ierr)
     ! Received ghost value from rank-1 stored in bcnd2h(:,rank*m-1:rank*m)
     call MPI_RECV(bcnd2h(0,rank*m-1),(n(1)/2+3)*2, MPI_INTEGER, rank-1, &
          0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
     
     if(rank.lt.nproc-1) then
        ! Send ghost value bcnd2h(:,(rank+1)*m-1:(rank+1)*m) to rank+1
        call MPI_SEND(bcnd2h(0,(rank+1)*m-1),(n(1)/2+3)*2, MPI_INTEGER, rank+1, &
             0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
        ! Received ghost value from rank+1 stored in bcnd2h(:,j=(rank+1)*m+1)
        call MPI_RECV(bcnd2h(0,(rank+1)*m+1),n(1)/2+3, MPI_INTEGER, rank+1, &
             0, MPI_COMM_WORLD, status, ierr)
     endif
  else ! Rank is even 
     if(rank.gt.0) then
        call MPI_RECV(bcnd2h(0,rank*m-1),(n(1)/2+3)*2, MPI_INTEGER, rank-1, & 
                0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
        call MPI_SEND(bcnd2h(0,rank*m+1),n(1)/2+3, MPI_INTEGER, rank-1, &
             0, MPI_COMM_WORLD, ierr)
     endif
     
        if(rank.lt.nproc-1) then
           call MPI_RECV(bcnd2h(0,(rank+1)*m+1),n(1)/2+3, MPI_INTEGER, rank+1, &
                0, MPI_COMM_WORLD, status, ierr)
           call MPI_SEND(bcnd2h(0,(rank+1)*m-1),(n(1)/2+3)*2, MPI_INTEGER, rank+1, &
                0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
        endif
     endif

100 continue

  return
end subroutine restriction


subroutine getres(uorg,rhsorg,h,ir,jr,n,r,flag)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      July/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Calculate residual
!     NOTE:          
!     --------------------------------------------------------------
  implicit none
  integer:: ir,jr,n(2),flag
  real(kind=8):: h(3),uorg(0:n(1)+2,0:n(2)+2),rhsorg(n(1)+1,n(2)+1)
  real(kind=8):: uij,ue,uw,un,us,urhs,r
  include 'constants.h'
  include 'mg.h'

  
  ! data points used in the calculation of the residual
  uij=uorg(ir,jr)
  un=uorg(ir,jr+1)
  us=uorg(ir,jr-1)
  ue=uorg(ir+1,jr)
  if( ir.eq.1 .and. flag.eq.1 ) then
     ! Neumann BC's
     uw=0.d0
     aw=0
     ae=2.d0*ae
  else
     uw=uorg(ir-1,jr)
  endif
  urhs=rhsorg(ir,jr)
  
  ! ...and the residual
  r= urhs - ( aw*uw + ae*ue + an*un + as*us + ac*uij )

  return
end subroutine getres

subroutine prolongation(u,e2h,bcnd,n,rank,nproc)
!     ==============================================================
!     VERSION:         0.3
!     LAST MOD:      Dec/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      2h -> h grid mapping using a 9 point bilinear 
!                    (Full Weighted) interpolation      
!     NOTE:          We do not need to compute 0 and n/2+1  
!     --------------------------------------------------------------
  implicit none
  include 'mpif.h'
  integer:: i,j,ir,jr,jm,jp,n(2),rank,nproc,m
  real(kind=8):: u(0:n(1)+2,0:n(2)+2),e2h(0:n(1)/2+2,0:n(2)/2+2), &
       r1,r2,r3,r4,mr
  integer:: bcnd(0:n(1)+2,0:n(2)+2)

  !
  ! Compute size of local block
  !
  m=(n(2)/2)/nproc
  mr= (real(n(2))/2.d0)/real(nproc)

  if( (real(m)-mr).eq.0.d0 ) then     
     jm= rank*m+1    
     jp= (rank+1)*m+1
  else ! Less than 1 node per proc
     jm= 1
     jp= n(2)/2+1
  endif

  !
  ! Sweep the coarse grid
  !

  !$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(ir,jr,r1,r2,r3,r4)
  !$OMP DO
  do j=jm,jp
     do i=1,n(1)/2+1
        
        ! Load all the required data points
        r1=e2h(i,j)
        r2=e2h(i,j-1)
        r3=e2h(i-1,j-1)
        r4=e2h(i-1,j)

        ! Compute fine grid address
        ir=2*i-1
        jr=2*j-1

        ! Prolong
        if(bcnd(ir,jr).le.0) u(ir,jr)= u(ir,jr) + r1
           
        if(bcnd(ir,jr-1).le.0) &
             u(ir,jr-1)= u(ir,jr-1) + (r1+r2)/2.0d0        

        if(bcnd(ir-1,jr-1).le.0) &
             u(ir-1,jr-1)= u(ir-1,jr-1) + (r1+r2+r3+r4)/4.0d0
           
        if(bcnd(ir-1,jr).le.0 ) &
             u(ir-1,jr)= u(ir-1,jr) + (r1+r4)/2.0d0 
           
     enddo ! repeat over the whole coarse grid
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  return
end subroutine prolongation

