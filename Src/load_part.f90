subroutine load_part(n,h,bcnd,np,vxp,ntype,nmax,kq,ni0,np_tot,nproc,&
     iseed,sum_dEk,Nh,n_icp,mpi_rank,nproc_mpi)
!     ==============================================================
!     VERSION:         0.3
!     LAST MOD:      FEB/15
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:   
!     NOTE:          
!     --------------------------------------------------------------
  use omp_lib 
  implicit none
  include 'mpif.h'
  include 'particle_info.h'
  integer ierr,mpi_rank,nproc_mpi
  integer:: ntype,ptype,nmax,n(3),iproc,nproc,Nh(n_icp,nproc),n_icp
  real(kind=8):: h(3)
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2),np_tot(ntype,nproc), &
       sum_np_tot_OMP(ntype),sum_np_tot(ntype)
  real(kind=8):: vxp(6,nmax,ntype,nproc)
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,ntype,nproc),ni0(npart),&
       sum_dEk(n_icp,nproc),kq(0:n(1)+2,0:n(2)+2)
  ! Macroscopic parameters 
  integer:: iseed(nproc),iseed_OMP
  
  if(mpi_rank.eq.0) print*, 'Loading particles ...' 
  
  !
  ! Loop over OpenMP processes
  !
  
  !$OMP PARALLEL PRIVATE(iproc,ptype,iseed_OMP)
  ! Get processor id (from 0 to nproc-1)
  iproc= omp_get_thread_num() + 1
  iseed_OMP= iseed(iproc)
  call load_part_OMP(n,h,bcnd,np,vxp,ntype,nmax,kq,ni0,np_tot,&
       nproc,iseed_OMP,sum_dEk,Nh,n_icp,nproc_mpi,iproc)
  iseed(iproc)=iseed_OMP
  !$OMP END PARALLEL

  sum_np_tot= 0
  sum_np_tot_OMP= SUM(np_tot,DIM=2)
  call MPI_ALLREDUCE(sum_np_tot_OMP, sum_np_tot, ntype, MPI_INTEGER, MPI_SUM, &
       MPI_COMM_WORLD, ierr)
     
  !
  ! Write info on screen
  !
  if(mpi_rank.eq.0) then
     write(*,100) pname(1:ntype)
100  format(' Particles: ',10(a5,1x))
     write(*,101) sum_np_tot ! Sum over the second index (nproc)
101  format(' np=',10(1x,i9))
  endif
  
  return
end subroutine load_part

subroutine load_part_OMP(n,h,bcnd,np,vxp,ntype,nmax,kq,ni0,&
     np_tot,nproc,iseed,sum_dEk,Nh,n_icp,nproc_mpi,iproc)
!     ==============================================================
!     VERSION:         0.3
!     LAST MOD:      OCT/16
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:   Loop order not optimized. This does not matter 
!                 much in terms of calculation time because this
!                 subroutine is only called once at the beginning.   
!     NOTE:          
!     --------------------------------------------------------------
  implicit none
  include 'mpif.h'
  include 'particle_info.h'
  integer nproc_mpi
  integer:: ix,iy,j,jmax,k,n_icp,i_icp
  integer:: ntype,ptype,nmax,n(3),iproc,nproc,Nh(n_icp,nproc)
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2),np_tot(ntype,nproc)
  real(kind=8):: h(3),vxp(6,nmax,ntype,nproc),x,y,z,vx,vy,vt,dEk,rnd(2),&
       ran2,ki(4),kp,px,py,ni0(npart),sum_dEk(n_icp,nproc),vz_sav(ntype)
  real(kind=8):: np(0:n(1)+2,0:n(2)+2,ntype,nproc), &
       kq(0:n(1)+2,0:n(2)+2)
  ! Macroscopic parameters 
  integer:: iseed

  !
  ! Initialize particle counter, variables & arrays
  !
  k=0
  np_tot(:,iproc)=0
  vz_sav=0.d0
  np(:,:,:,iproc)=0.d0
  sum_dEk(:,iproc)=0.d0
  Nh(:,iproc)=0

  ! Total number of particles
  jmax= NINT(real(np_cell*n_cell)/real(nproc_mpi)/real(nproc))
  do j=1,jmax
                 
     ! random loading
70   continue

     ! Flattop profile
     rnd(1)=ran2(iseed)
     x= rnd(1)*x_load
     rnd(1)=ran2(iseed)
     y= rnd(1)*ymax
     rnd(1)=ran2(iseed)
     z= rnd(1)*zmax
        
     ! Get particle left grid index
     ix= INT( x/h(1) ) + 1
     iy= INT( y/h(2) ) + 1

     ! Load uniquely inside simulation domain
     if( bcnd(ix,iy).ge.1 .and. bcnd(ix+1,iy).ge.1 .and. &
          bcnd(ix+1,iy+1).ge.1 .and. bcnd(ix,iy+1).ge.1 ) goto 70
     
     do ptype=1,ntype
        
        ! Probability to inject ptype particle
        rnd(1)= ran2(iseed)
        if( rnd(1).gt.ni0(ptype) ) goto 80
        
        ! Add particle to counter
        np_tot(ptype,iproc)= np_tot(ptype,iproc) + 1
        
        ! Get thermal velocity
        vt= vt0(ptype)
           
        ! Particle index
        k= np_tot(ptype,iproc)
        
        ! Warning
        if(k.gt.nmax) then
           print*, 'k > nmax in load_part'
           print*, 'Abort calculation ...'
           call stop_calculation
        endif
        
        ! Get particle location
        vxp(1,k,ptype,iproc)= x ! same location for all particles
        vxp(2,k,ptype,iproc)= y
        vxp(3,k,ptype,iproc)= z
        
        rnd(1)= ran2(iseed)
        rnd(2)= ran2(iseed)
        call load_gauss(vx,vy,vt,rnd)
        vxp(4,k,ptype,iproc)= vx
        vxp(5,k,ptype,iproc)= vy
        if(vz_sav(ptype).eq.0.d0) then
           rnd(1)= ran2(iseed)
           rnd(2)= ran2(iseed)
           call load_gauss(vx,vy,vt,rnd)
           vxp(6,k,ptype,iproc)= vx
           vz_sav(ptype)= vy
        else
           vxp(6,k,ptype,iproc)= vz_sav(ptype)
           vz_sav(ptype)= 0.d0
        endif
           
        !
        ! Calculate initial density
        !
        px=( ix*h(1) - x )/h(1)
        
        if(flag_1D.eq.0) then
           py=( iy*h(2) - y )/h(2)
        
           kp=Nm(ptype)/(h(1)*h(2)*zmax)
           ki(1)= kp*px*py
           ki(2)= kp*(1.d0-px)*py
           ki(3)= kp*(1.d0-px)*(1.d0-py)
           ki(4)= kp*px*(1.d0-py)       
        else
           ! 1D case
           kp=Nm(ptype)/h(1) ! per m^2
           ki(1)= kp*px
           ki(2)= kp*(1.d0-px)
        endif
        
        ! Charge assigned to the bottom left grid point
        np(ix,iy,ptype,iproc)= np(ix,iy,ptype,iproc) + kq(ix,iy)*ki(1)

        ! Charge assigned to the bottom right grid point
        np(ix+1,iy,ptype,iproc)= np(ix+1,iy,ptype,iproc) + kq(ix+1,iy)*ki(2)

        if(flag_1D.eq.0) then
           ! Charge assigned to the top right grid point
           np(ix+1,iy+1,ptype,iproc)= np(ix+1,iy+1,ptype,iproc) + kq(ix+1,iy+1)*ki(3)

           ! Charge assigned to the top left grid point
           np(ix,iy+1,ptype,iproc)= np(ix,iy+1,ptype,iproc) + kq(ix,iy+1)*ki(4)
        endif

        ! Electrons Maxwellian heating
        if( ptype.eq.1 .and. Pabs(1).gt.0.d0 ) then 

           if(flag_c.eq.0) then ! Slit
              if( x.lt.xl_pow .or. x.gt.xr_pow .or. &
                   y.lt.yl_pow .or. y.gt.yr_pow ) goto 90
           else ! Disk
              if( ((x-xa)**2 + (y-ymax/2.d0)**2).gt.dr**2 ) goto 90
           endif
           
           ! Calculate kinetic energy of macroparticle
           dEk= 0.5d0*Nm(ptype)*mass(ptype)*( &
                vxp(4,k,ptype,iproc)*vxp(4,k,ptype,iproc) + &
                vxp(5,k,ptype,iproc)*vxp(5,k,ptype,iproc) + &
                vxp(6,k,ptype,iproc)*vxp(6,k,ptype,iproc) )
           
           i_icp=INT(n_icp*y/ymax) + 1
           sum_dEk(i_icp,iproc)= sum_dEk(i_icp,iproc) + dEk
           Nh(i_icp,iproc)= Nh(i_icp,iproc) + 1
           
90         continue
        endif
                
80      continue
     enddo
  enddo

  return

end subroutine load_part_OMP

