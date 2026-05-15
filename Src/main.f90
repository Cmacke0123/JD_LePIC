!     ================================================================
!     VERSION:         4.7.4
!     LAST MOD:       Sep/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      A 2.5D explicit parallel hybrid OpenMP/MPI PIC code
!                   
!     NOTE:         Domain extend from ix(iy)=1 to ix(iy)=nx(ny)+1
!                   ix(iy)=1 & ix(iy)=nx(ny)+1 are boundary conditions, 
!                   unknown interior points are located between ix(iy)=2 
!                   and ix(iy)=nx(ny)
!                   Example: nx=5 & along (Ox)
!                   BC                  BC
!                    |---|---|---|---|---|
!                    1   2   3   4   5   6
!                   (0)                (xmax)  dx=xmax/nx
!
!
! Copyright or © or Copr. Gwenael Fubiani (2022/11/22)
!
! gwenael.fubiani@cnrs.fr
!
! This software is a computer program whose purpose is to model low
! temperature plasmas with a Particle-In-Cell algorihm in 2.5D-3V
! dimensions.
!
! This software is governed by the CeCILL-C license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL-C
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 

! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
!
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-C license and that you accept its terms.
!     ----------------------------------------------------------------
program main
  use omp_lib 
  implicit none
  include 'mpif.h'
  include 'particle_info.h'
  include 'constants.h'
  ! MPI
  integer ierr,mpi_rank,nproc_mpi
  ! Loop indexes
  integer:: it,ix,iy,i,j,i_rg,ptype,niter,jl,jr,m,n_icp,i_icp,&
       NxNy,jpl,jpr,itot,flag_pardiso,nmax_ps,flag_nopart,iyl,iyr,flag_diag
  ! Simulation parameters (integers)
  integer:: n(3),na(3),n_B(2),ncycle,nmax,nmax_tmp,ngrid,ntype,ntmax,nsav, &
       ng,iproc,nproc,ctime(0:10),MSTIMER,B_dir,namlen,&
       sum_Nh_tmp,sum_Nh,n_neu,flag_sav,flag_wrt,ns_convP,nseq,igrid,&
       flag_updatephi
  real(kind=8):: h(3),Vc,h_B(2)
  parameter (nmax_tmp=2*10**8)
  ! Physical scales
  real(kind=8):: lbd_d,wp,kt,wc,ni0(npart),Te,sum_dEk_tmp,&
       sum_dEk_tot,I_inj
  ! Arrays
  integer, allocatable:: bcnd(:,:),iseed(:),shift(:),length(:),cnt_b(:,:)
  integer(kind=1), allocatable:: flag_dead(:,:,:)
  real(kind=8), allocatable:: phi(:,:),rhs(:,:),sum_dEk(:,:), &
       vxp(:,:,:,:),np(:,:,:,:),Ei(:,:,:),p_mac(:,:,:,:),&
       P_loss(:,:,:),kq(:,:),Bi(:,:,:),rhs_par(:),&
       np_red(:,:,:),phi_tmp(:,:),uB(:,:,:),vt(:)
  ! Sorting
  integer:: pl_max,nsort
  integer, allocatable:: Plist(:,:,:),N_inj(:,:)
  real(kind=8), allocatable:: sorting(:,:)
  ! Macroscopic parameters
  integer, allocatable:: np_tot(:,:),&
       Nh(:,:),sum_np_tot_tmp(:),sum_np_tot(:),ip_ps(:)
  real(kind=8), allocatable:: cnt_col(:,:,:,:),p_mts(:,:,:,:,:),&
       phi_avg(:,:,:),data_pavg(:,:,:,:),cnt_dead(:),sum_q(:,:),Vgrd(:),&
       sum_q_y(:,:,:,:),sum_q_red_y(:,:,:),sum_I_sav(:,:,:),phase_space(:,:,:),&
       ss2D(:,:,:,:,:)
  ! Collisions
  integer:: ncol_mx,npt_mx
  parameter (ncol_mx=1200, npt_mx=1500)
  integer:: sig_list(npart,ncol_mx),col_info(ncol_mx,10), &
       sig_type(ncol_mx),scol_rank(npart,ncol_mx)
  real(kind=8):: sig(npt_mx,ncol_mx),sig_Er(npt_mx),sig_Eex(ncol_mx,2), &
       scol_info(ncol_mx,4),sigv_mx(npart,ncol_mx),kd
  ! Simulation parameters (reals)
  integer:: cnt_i,cnt_f,cnt_rate,navg,n_mts,cnt_avg(2),lgh,sav_np,iseed_OMP,&
       np_pos,np_pos0,np_dup
  real(kind=8):: res,eps,ksor,time,tmax,sum_time,a1,a2,xl_rg(nm_rg), &
       xr_rg(nm_rg),yl_rg,yr_rg,dtime,cnt_dead_tmp,&
       Pext,Pabs_cor,sum_q_tmp,Ca,phi0_sav,phi0_RF,f0_RF,phi1_RF,f1_RF
  character:: rname*20,name*20,pnum*1,pnum_bck*3

  CALL MPI_Init(ierr)                             ! starts MPI
  CALL MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)  ! get current process id
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc_mpi, ierr) ! get # of procs

  call system_clock(cnt_i,cnt_rate)

  !
  ! Write code info on screen
  !
  if(mpi_rank.eq.0) call introduction

  !
  ! Constants
  !
  qe=1.60217646e-19 ! Coulombs
  c=2.99792458d8 ! m/s
  eps0=8.854187817d-12 ! S.I.
  pi=4.d0*datan(1.d0)
  amu=1.66053886d-27


  !
  ! Get number of OMP threads
  !
  nproc= omp_get_max_threads()

  if(mpi_rank.eq.0) then
     print'(1x,"Number of OpenMP threads= ",i3)',nproc
     print'(1x,"Number of MPI threads= ",i3)',nproc_mpi
  endif

  !
  ! Read input parameters
  !
  call read_input(n,na,tmax,nsav,eps,kt,rname,ngrid,ng,n_icp, &
       xl_rg,xr_rg,yl_rg,yr_rg,I_inj,np_dup,mpi_rank,Ca,n_B,&
       phi0_RF,f0_RF,phi1_RF,f1_RF)

  flag_sav= 0
  if(nsav.lt.0) then
     flag_sav= 1
     nsav= ABS(nsav)
  endif

  flag_bak= 0
  if(nbak.lt.0) then
     flag_bak= 1
     nbak= ABS(nbak)
  endif
  
  ! Correct for nmax
  nmax= NINT(nmax_tmp/real(nproc)/real(nproc_mpi))

  if( mpi_rank.eq.0 .and. k_eps0.ne.1 ) &
       print*, 'eps0 HAS BEEN RE-SCALED: k_eps0=',k_eps0
  eps0=k_eps0*eps0
  
  ! Warning
  if(k_eps0.eq.0) then
     print*, 'Please correct k_eps0'
     call stop_calculation
  endif

  !
  ! Read cross sections from files
  !
  call read_reactions(sig,sig_Er,sig_list,sig_Eex,ncol_mx,sig_type, &
       npt_mx,rname,col_info,scol_rank,scol_info,sigv_mx,ni0,ntype,n_neu,mpi_rank)

  ntype= ntype - n_neu
  ! Neutral density
  np_mx(ntype+1:ntype+n_neu)= ni0(ntype+1:ntype+n_neu)

  ! Electron fraction
  if(tag_neg.gt.0) ni0(1)= ni0(1) - ni0(tag_neg)

  flag_diffsrc=0
  if(( I_inj.gt.0.d0 .and. (SUM(charge(1:ntype)*ni0(1:ntype))).ne.0.d0 .or. ABS(opt_inj).ge.3) ) then
     flag_diffsrc=1
     if(mpi_rank.eq.0) &
          print*, 'Potential energy from particles injected in the volume is accounted for'
  endif

  !
  ! Allocate arrays
  !
  allocate (  phi(0:n(1)+2,0:n(2)+2), &
       vxp(6,nmax,ntype,nproc),np(0:n(1)+2,0:n(2)+2,ntype,nproc), &
       Ei(2,0:n(1)+2,0:n(2)+2),p_mac(ntype,2,0:ngrid,nproc),P_loss(4,ntype,nproc), &
       np_red(0:n(1)+2,0:n(2)+2,ntype),kq(0:n(1)+2,0:n(2)+2), &
       Bi(4,0:n_B(1)+2,0:n_B(2)+2),sum_dEk(n_icp,nproc),&
       uB(na(1)+1,na(2)+1,ntype),vt(n_icp) )
  allocate (  phi_tmp(0:n(1)+2,0:n(2)+2), &
       shift(0:nproc_mpi-1), length(0:nproc_mpi-1) )
  allocate ( np_tot(ntype,nproc),Nh(n_icp,nproc),sum_np_tot_tmp(ntype),sum_np_tot(ntype), &
       N_inj(ntype,nproc),cnt_dead(nproc),sum_q(ntype,nproc),&
       sum_q_y(ngrid,0:n(2)+2,ntype,nproc),sum_q_red_y(ngrid,0:n(2)+2,ntype),&
       sum_I_sav(ngrid,0:n(2)+2,ntype) )
  allocate ( cnt_col(3,ncol_mx,nproc,nm_rg),iseed(nproc) )
  allocate ( bcnd(0:n(1)+2,0:n(2)+2) )
  allocate ( flag_dead(nmax,ntype,nproc),cnt_b(2,nproc) )
  allocate ( phi_avg(3,0:n(1)+2,0:n(2)+2), &
       data_pavg(7,na(1)+1,na(2)+1,ntype), Vgrd(ngrid),&
       ss2D(2,0:n(1)+2,0:n(2)+2,ntype,nproc) )
  nmax_ps= NINT(80000./real(nproc))
  allocate ( phase_space(nmax_ps,4,nproc),ip_ps(nproc) )

  pl_max= n(1) + n(1)*n(2) + 1 ! (n(1)-1) + n(1)*(n(2)-1) + 1
  allocate ( Plist(0:pl_max,ntype,nproc) )

  !
  ! Initialize variables
  !
  cnt_col=0.d0
  time=0.d0
  Ei=0.d0
  Bi=0.d0
  phi=0.d0
  ctime=0
  flag_dead=0
  Bmax=0.d0
  cnt_avg=0
  p_mac=0.d0
  P_loss=0.d0
  data_pavg=0.d0
  phi_avg=0.d0
  N_inj=0
  cnt_dead=0.d0
  flag_wrt=0
  Pabs_cor=1.d0
  sum_q= 0.d0
  sum_q_y= 0.d0
  cnt_b= 0
  cnt_plt= 0
  sum_I_sav= 0.d0
  ip_ps=0
  ss2D= 0.d0
  th_B=0.d0

  !
  ! Generate boundary counditions
  !
  call generate_boundary(phi,n,h,bcnd,Vgrd,ngrid,n_icp,mpi_rank)

  ! Domain decomposition for Poisson solver
  m= n(2)/nproc_mpi

  !
  ! PDE solvers
  !
  if(opt_solver.eq.1 .or. opt_solver.eq.3) then ! MG & Thomas alg.

     allocate ( rhs(n(1)+1,n(2)+1) )

     if(opt_solver.eq.1) then ! MG
        ncycle= 20000 ! PDE max cycle
        
        ! Calculate maximum number multigrid levels
        a1= alog10(real(n(1)))/alog10(2.)
        a2= alog10(real(n(2)))/alog10(2.)
        if( ABS(NINT(a1)-a1).le.1d-3 .and. ABS(NINT(a2)-a2).le.1d-3 ) then
           ng= NINT(MIN(a1,a2) )
        else
           a1= real(n(1))/2**(ng-1)
           a2= real(n(2))/2**(ng-1)
           ! Warning
           if( ABS(NINT(a1)-a1).gt.0 .or. ABS(NINT(a2)-a2).gt.0 ) then
              print*, 'Incorrect number of grid levels in MG'
              print*, 'Please correct ...'
              call stop_calculation
           endif
        endif
     endif
     
     ! Direct solver
  else ! Parameters for Pardiso

     if(flag_pbc.eq.0) then
        NxNy=(n(1)+1)*(n(2)+1)
        jpl=0
        jpr=n(2)
     else
        NxNy=(n(1)+1)*n(2)
        jpl=1
        jpr=n(2)
     endif

     allocate(rhs_par(NxNy))

     ! Setup Dirichlet BC's values 
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
           if(bcnd(ix+1,iy+1).ge.1) then 
              rhs_par(itot)= phi(ix+1,iy+1)
           endif
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL 

  endif

  !
  ! Coordinates of sub-regions
  !
  do i=1,n_rg
     if(xl_rg(i).lt.0) xl_rg(i)= 0
     if(xr_rg(i).gt.xmax) xr_rg(i)= xmax
     ixl_rg(i)= NINT(xl_rg(i)/h(1)) + 1
     ixr_rg(i)= NINT(xr_rg(i)/h(1)) + 1
  enddo
  if(yl_rg.lt.0) yl_rg= 0
  if(yr_rg.gt.ymax) yr_rg= ymax
  iyl_rg= NINT(yl_rg/h(2)) + 1
  iyr_rg= NINT(yr_rg/h(2)) + 1 


  !
  ! Define charge per cell coefficient
  !
  kq= 1.d0
  !$OMP PARALLEL
  !$OMP DO
  do ix=0,n(1)+2
     do iy=0,n(2)+2
        if(bcnd(ix,iy).ne.-1) kq(ix,iy)= 2.d0
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  if(flag_1D.eq.1)  kq(2:n(1),:)= 1.d0

  ! Write to a file
  open(30,file='DATA/kq.dat',form='UNFORMATTED')
  write(30) n(1),n(2)
  write(30) kq(1:n(1)+1,1:n(2)+1)
  close(30)

  !
  ! Open (some) output files
  !
  open(12,file='DATA/time.dat')
  open(40,file='DATA/residual.dat')

  !
  ! Magnetic field profile
  !
  if(flag_B.eq.1) then
     h_B(1)= xmax/n_B(1)
     h_B(2)= ymax/n_B(2)
  
     flag_gridB=0
     if( n_B(1).ne.n(1) .or. n_B(2).ne.n(2) ) flag_gridB=1
  
     do i=1,nB

        ! Skip
        if(B_file(i).eq.'s'.or.B_file(i).eq.'S') goto 10

        ! Extract file name
        name= B_name(i)
        do namlen=20,1,-1
           if( name(namlen:namlen).ne.' ' ) goto 20
        enddo
20      continue

        ! B-field map from an external file
        if( B_file(i).eq.'y'.or. B_file(i).eq.'Y' ) then      
           do B_dir=1,3
              if(B_dir.eq.1) name= name(1:namlen)//'Bx.dat'
              if(B_dir.eq.2) name= name(1:namlen)//'By.dat'
              if(B_dir.eq.3) name= name(1:namlen)//'Bz.dat'
              call read_Bfield_map(Bi,n_B,name,namlen+6,B_dir,B_scale(i),mpi_rank)
           enddo
        endif

        ! B-field map from gaussian profile/PG current
        if( B_file(i).eq.'n'.or. B_file(i).eq.'N' ) then 
           call find_B_dir(B_info(i),B_dir)
           if( name(1:1).eq. 'g' .or. name(1:1).eq. 'G' ) then
              if(mpi_rank.eq.0) print*, 'Gaussian magnetic filter field profile'
              call gaussian_Bfield(Bi,n_B,h_B,B0(i),x0(i),dL(i),B_dir)
           endif
           if( name(1:1).eq. 'p' .or. name(1:1).eq. 'P' ) then
              if(mpi_rank.eq.0) print*, 'Magnetic filter field generated by a PG current' 
              call PG_current(Bi,n_B,h_B,B0(i),B_dir)
           endif
           if( name(1:1).eq. 'e' .or. name(1:1).eq. 'E' ) then
              if(mpi_rank.eq.0) print*, 'Cusp magnetic field profile is generated along: ',B_info(i) 
              call EE_magnets(Bi,n_B,h_B,B0(i),x0(i),y0(i),dL(i),B_dir)
           endif
        endif

10      continue
     enddo
 
     ! Calculate ||B||
     !$OMP PARALLEL REDUCTION(MAX:Bmax)
     !$OMP DO
     do iy=0,n_B(2)+2
        do ix=0,n_B(1)+2
           Bi(4,ix,iy)= dsqrt( Bi(1,ix,iy)*Bi(1,ix,iy) + &
                Bi(2,ix,iy)*Bi(2,ix,iy) + &
                Bi(3,ix,iy)*Bi(3,ix,iy) )
           Bmax= MAX(Bmax,Bi(4,ix,iy))
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL 

     if(flag_icp.eq.1) then
        !$OMP PARALLEL PRIVATE(i)
        do i=2,n_icp
           !$OMP DO
           do iy=1,n_B(2)/n_icp+1
              do ix=1,n_B(1)+1
                 Bi(:,ix,(i-1)*n_B(2)/n_icp+iy)= Bi(:,ix,iy)
              enddo
           enddo
           !$OMP END DO NOWAIT
        enddo
        !$OMP END PARALLEL 
     endif

     ! Write to files
     if( n_B(1).gt.1 .and. mpi_rank.eq.0) then      
        open(30,file='DATA/Bx.dat',form='UNFORMATTED')
        open(31,file='DATA/By.dat',form='UNFORMATTED')
        open(32,file='DATA/Bz.dat',form='UNFORMATTED')
        do i=1,3
           write(30+(i-1)) n_B(1),n_B(2)
           write(30+(i-1)) Bi(i,1:n_B(1)+1,1:n_B(2)+1)*1.d4
           close(30+(i-1))
        enddo
     endif

     if(opt_inj.eq.4 .and. flag_nmn.eq.1 .and. ABS(Bi(1,1,1)).gt.0.d0) then
        th_B= datan(Bi(2,1,1)/Bi(1,1,1)) ! arctg(By/Bx)
     endif
     
  endif

  !
  ! Simulation parameters
  !
  lbd_d=dsqrt( (eps0*Ti(1))/(n0*ABS(charge(1))) ) ! Debye length
  wp=dsqrt( n0*charge(1)**2./(eps0*mass(1)) ) ! Plasma frequency

  ! Cell volume
  Vc= h(1)*h(2)*zmax
  if(flag_1D.eq.1) Vc= h(1) ! per m^2
  
  do ptype=1,ntype+n_neu
     ! Thermal velocities [sqrt(2*kT/m)]
     vt0(ptype)=dsqrt(2.d0*qe*Ti(ptype)/ABS(mass(ptype)))
     ! Charged macroparticle weight
     if(ABS(charge(ptype)).gt.0) then
        Nm(ptype)=n0*Vc/ABS(np_cell)
     else ! Gas macroparticle weight
        Nm(ptype)=ngas*Vc/ABS(np_cell)
     endif
  enddo
  dt= kt*MIN(h(1),h(2))/vt0(1) ! Time step k*dr/vte, k is arbitrary

  ntmax= 40000001
  nsort=10 ! Frequency of calls to sorting subroutine
  ns_coll= 1*nsort ! Frequency of calls to collision subroutine
  ns_heat=4 ! Frequency of calls to heating subroutine
  ns_inj=1 ! Frequency for particle injection
  if(flag_bak.eq.0) then
     navg=5*nsort ! Frequency for averaging
  else
     navg= nsav
  endif
  ! Save sequence of profiles every tseq seconds
  if(tseq.gt.0) then 
     nseq= MAX(NINT((tseq/dt)/nsav),1)*nsav 
     tseq= nseq*dt ! correct tseq
  endif
  nudt=nu_h*(ns_heat*dt)  ! % of electrons which should be heated
  nu_uplim(1)=4.d7 ! Maximum collision frequency (e-)
  nu_uplim(2:ntype)=5.d6 ! (ions)
  ns_Hm= 1 ! Negative ion injection frequency on the PE
  eheat_type=1 ! Electron velocity sampled from Maxwellian: 1= full velocity, 2= only increment dv
  if(Pabs(1).le.0.d0) eheat_type=1 ! when fixed Te is set
  flag_gpp=0
  if( opt_inj.lt.0 .and. Pabs(1).gt.0.d0) flag_gpp=1 ! Gaussian power profile

  ! Floating end-plates (3rd dimension) or dielectric BC's
  if(flag_float.eq.1) phi0_sav= phi0

  if(mpi_rank.eq.0) then
     print'(1x,"Frequency of calls to collision subroutine= ",i3)',ns_coll
     print'(1x,"Frequency of calls to electron heating subroutine= ",i3)',ns_heat
     if(flag_heat.eq.1) then
        if(eheat_type.eq.1) print'(1x,"Electron velocity sampled from a Maxwellian: velocity is replaced")'
        if(eheat_type.eq.2) print'(1x,"Electron velocity sampled from a Maxwellian: an incremental dv is added")'
     endif
     print'(1x,"Frequency of calls to sort subroutine= ",i3)',nsort
     if(flag_float.eq.1) print'(1x,"End-plates (3rd dim) are floating")'
     if(flag_nmn.eq.1) print'(1x,"Left-hand-side boundary condition set to Neumann type")'
     if(flag_pbc.eq.1) print'(1x,"Top/bottom boundary conditions set as periodic")'
     if( flag_pbc.eq.1 .and. flag_spec.eq.1 ) print'(1x,"Specular reflexion for particles at the top and bottom boundaries")'

     ! Warning
     if(nudt.gt.1.d0) then
        print*, 'Warning, nu*dt>1, please correct ...'
        print*, 'nu*dt=',nudt
        call stop_calculation
     endif
  endif

  !
  ! Set random seed number
  !
  do iproc=1,nproc
     iseed(iproc)= 123456*iproc*(10*mpi_rank+1)
  enddo

  !
  ! Load particles
  !

  ! RF heating & injection box
  if(xr_pow.gt.xmax) xr_pow= xmax
  if(xl_pow.lt.0) xl_pow= 0.d0
  if(yr_pow.gt.ymax) yr_pow= ymax
  if(yl_pow.lt.0) yl_pow= 0.d0

  if( xl_pow.eq.0.d0 ) then
     ! Half-gaussian power profile option  
     xa=0.d0
     dr= xr_pow
  endif

  ! Location of the plasma electrode
  if(x_load.gt.xmax) x_load= xmax
  nx_PE= INT( x_load/h(1) ) + 1

  if(mpi_rank.eq.0) then
     print'(1x,"Heating region (cm): x<=",f6.1,", x>=",f6.1", y<=",f6.1", y>=",f6.1)',&
          xl_pow*1.d2,xr_pow*1.d2,yl_pow*1.d2,yr_pow*1.d2
     print'(1x,"Particle loading region: nx<= ",i5)',nx_PE
  endif
     
  ! Estimate # of cells inside simulation domain
  n_cell=0
  !$OMP PARALLEL REDUCTION(+:n_cell)
  !$OMP DO
  do iy=1,n(2)
     do ix=1,nx_PE
        if( bcnd(ix,iy).eq.-1 .or. bcnd(ix+1,iy).eq.-1 .or. &
             bcnd(ix+1,iy+1).eq.-1 .or. bcnd(ix,iy+1).eq.-1 ) &
             n_cell= n_cell + 1        
     enddo
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  if(flag_restart.eq.0) then
     if(np_cell.gt.0) then
        call load_part(n,h,bcnd,np,vxp,ntype,nmax,kq,ni0,np_tot,&
             nproc,iseed,sum_dEk,Nh,n_icp,mpi_rank,nproc_mpi)
     else
        np= 0.d0
        np_tot(1:ntype,:)= 0
     endif
  else
     call restart(n,h,np,vxp,nmax,ntype,kq,time,nproc,&
          np_tot,sum_dEk,Nh,n_icp,np_dup,mpi_rank,nproc_mpi,Pabs_cor)
     ! Load potential on dielectric surfaces
     if(flag_dielec.eq.1) then
        if(flag_restart.eq.2) then
           open(41,file='SAV_DATA/phi.bak',form='UNFORMATTED')
        else
           open(41,file='DATA.BAK/phi.bak',form='UNFORMATTED')
        endif
        read(41) phi(0:n(1)+2,0:n(2)+2)
        close(41)
     endif
     if( tseq.gt.0.d0 .and. time.ge.tseq_init .and. time.le.tseq_final ) &
          cnt_plt= NINT((time-tseq_init)/tseq)
  endif

  !
  ! Print info. on screen
  ! 
  if(mpi_rank.eq.0) then
     print*, 'lde (mm)=',lbd_d*1.d3
     print*, 'lpe (mm)=',2*pi*c/wp*1.d3
     print*, 'vt (m/s)=',vt0(1:ntype) ! Electron thermal velocity
     print*, '<|v|> (m/s)=',2.d0/dsqrt(pi)*vt0(1) ! Electron electron speed
     print*, 'dx (mm)=',h(1)*1.d3,'dx/lde=',h(1)/lbd_d
     if(flag_1D.eq.0) print*, 'dy (mm)=',h(2)*1.d3,'dy/lde=',h(2)/lbd_d
     print*, 'Nm=',Nm(1)
     print*, 'dt (s)=',dt,', wpe*dt=',wp*dt
     print*, 'vte*dt/dx=',vt0(1)*dt/MIN(h(1),h(2))

     if(Bmax.gt.0.d0) then
        if(n_B(1).eq.1) then
           print'(1x,"B-field is constant")'
           print'(1x,"Bx(T)= ",es10.2,", By(T)= ",es10.2,", Bz(T)= ",es10.2)', Bi(1,1,1), Bi(2,1,1), Bi(3,1,1)
        endif
     
        wc=qe*Bmax/mass(1)
        a1= 2.*datan(wc*dt/2.)/(wc*dt)
        a2= dsqrt(1+(wc*dt/2.)**2)
        print'(1x,"Bm(T)= ",es10.2,", wc*dt= ",f5.2,", wc (/s)= ",es12.4,", & 
             re (mm)= ",f5.2)', Bmax,wc*dt,wc,vt0(1)/wc*1.d3
        print'(1x,"err: w*/wc= ",f5.2,", r*/re= ",f5.2)', a1,a2
     endif
     
     if(opt_solver.eq.1) then
        print'(1x,"Multigrid PDE Solver is used, eps= ",es10.2)',eps
        print'(1x,"ng= ",i2)',ng
     endif
     if(opt_solver.eq.2) print'(1x,"Direct Solver Pardiso is used")'
     print'(1x,"nx=",i5,", ny=",i5," # of cells: ",i8)',n(1),n(2),n_cell 
     if(navg.eq.nsav) then
        print'(1x,"Data will not be averaged!")'
     else
        print'(1x,"Frequency for averaging= ",i3)',navg
     endif
     if(flag_sav.eq.0 .and. tseq.gt.0 .and. time.le.tseq_final) &
          print'(1x,"Frequency for saving sequence of profiles= ",i6,", dt(us)=",f5.2,", cnt_plt=",i4)', &
          nseq,tseq*1.d6,cnt_plt

     if(Pabs(1).gt.0) then
        if(flag_gpp.eq.0) then
           print*, 'Absorbed power profile is flattop'
        else
           print*, 'Absorbed power profile is Gaussian'
        endif
     endif

     if(opt_inj.lt.0 .and. I_inj.gt.0) print'(1x,"Particle injection profile is Gaussian, dr(mm)=",f6.2)',dr*1.d3
     if(flag_1D.eq.1) print*, '1D calculation option is set.'
     if(flag_sec.ne.0) then
        print'(1x,"Secondary particle emission on grid xg=",f6.2,1x,"cm, ig=",i4,", dir=",i2)',xg_sec*1.d2,INT(xg_sec/h(1))+1,dir_sec
     else
        if(jne.gt.0.d0) &
             print'(1x,"Particle flux emmited on grid located at xg=",f6.2,1x,"cm, ig=",i4)',xg1*1.d2,ind_g
     endif
     if(plt_src.eq.1) print'(1x,"Source/sink calculated from cross section option on")'
     if(opt_inj.eq.4 .and. flag_nmn.eq.1) then
        print'(1x,"Particle refluxing on the LHS boundary")'
        if(flag_B.eq.1) print'(1x,"Angle of magnetic field with respect to the (OX) axis (in degrees)=",f6.2)',th_B*180.d0/pi
     endif
  endif
  
  ! Iterative process to find Pabs : initialization
  if(flag_convP.eq.0) then
     ns_convP=10000
     ! Save initial number of positive ions
     np_pos0= SUM(np_tot(2:ntype,1:nproc))
     if(tag_neg.gt.0) np_pos0= np_pos0 - SUM(np_tot(tag_neg,1:nproc))
     if(mpi_rank.eq.0) print*, 'Pabs, I_inj or Vgrd will be calculated through an iterative process'
  endif

  ! Reset timer        
  ctime(0)= MSTIMER()
  ctime=0

  call system_clock(cnt_f)
  dtime=REAL(cnt_f-cnt_i)/REAL(cnt_rate)
  if(mpi_rank.eq.0) print'(" *** Startup time (s): ",f5.1," ***")', dtime

  if(mpi_rank.eq.0) then
     print*, ' '
     print*, 'Running simulation ...'
  endif

  !
  ! Start iteration 
  !
  do it=1,ntmax

     time= time + dt
     if(time.ge.tmax) exit ! Stop calculation

     ctime(0)= MSTIMER()

     ! Stat without any background plasma
     flag_nopart= 0
     if(SUM(np_tot(1:ntype,1:nproc)).eq.0) flag_nopart= 1

     flag_updatephi= 0
     ! Iterate to find Pabs, I_inj or V
     if( flag_convP.eq.0 .and. MOD(it,ns_convP).eq.0 ) then
        np_pos= SUM(np_tot(2:ntype,1:nproc))
        if(tag_neg.gt.0) np_pos= np_pos - SUM(np_tot(tag_neg,1:nproc))
        if( ABS(opt_inj).ne.2 .and. I_inj.gt.0.d0 ) then
           ! Inject particle beam
           I_inj=I_inj*real(np_pos0)/real(np_pos)
           if(nproc_mpi.gt.1) call MPI_Bcast(I_inj,1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        endif
        ! Inject external power
        if(Pabs(1).gt.0.d0) then
           Pabs=Pabs*real(np_pos0)/real(np_pos)
           if(nproc_mpi.gt.1) call MPI_Bcast(Pabs,n_icp, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)     
        endif
        ! Secondary electron emission
        if( Pabs(1).le.0.d0 .and. I_inj.le.0.d0 .and. (gam_sec.gt.0.d0 .or. jne.gt.0.d0) ) then
           ! Update flag
           flag_updatephi= 1
           ! Iterative Vgrd
           Vgrd(igrid_sec)=Vgrd(igrid_sec)*real(np_pos0)/real(np_pos)
        endif
     endif

     if(flag_RFpot.eq.1) then
        ! Warning
        if(flag_updatephi.eq.1) then
           print*, 'Cathode potential update used twice : for iterative proceduce and RF. Please correct...'
           call stop_calculation
        endif
        ! Update flag
        flag_updatephi=1
        ! Calculate RF potential
        Vgrd(igrid_sec)= phi0_RF*dsin(2.d0*pi*f0_RF*time) + &
             phi1_RF*dsin(2.d0*pi*f1_RF*time)
     endif

     if(flag_updatephi.eq.1) then
        ! Update cathode potential
        if(nproc_mpi.gt.1) call MPI_Bcast(Vgrd(igrid_sec),1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        !$OMP PARALLEL PRIVATE(igrid)
        !$OMP DO
        do iy=0,n(2)+2
           do ix=0,n(1)+2
              igrid= bcnd(ix,iy)
              if(igrid.eq.igrid_sec) &
                   phi(ix,iy)= Vgrd(igrid_sec)                
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
     endif
     
     !
     ! Dielectric boundary conditions
     !
     if(flag_dielec.eq.1) then
     
        ! Reduction
        sum_q_red_y= 0.d0
        !$OMP PARALLEL PRIVATE(kd)
        !$OMP DO       
        do iy=0,n(2)+2
           kd= 1.d0
           if( iy.eq.1 .or. iy.eq.n(2)+1 ) kd= 2.d0 ! half-cell
           do iproc=1,nproc
              sum_q_red_y(:,iy,:)= sum_q_red_y(:,iy,:) + kd*sum_q_y(:,iy,:,iproc)
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 
        
        call MPI_ALLREDUCE(MPI_IN_PLACE, sum_q_red_y, (ntype*ngrid*(n(2)+3)), &
             MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        
        ! Re-initialize surface charge counters on dielectrics
        sum_q_y=0.d0
        
        ! Save profiles
        sum_I_sav= sum_I_sav + sum_q_red_y/dt

        iyl= 2
        iyr= n(2)
        if(flag_pbc.eq.1) then
           iyl= 0
           iyr= n(2)+2
        endif
        if(flag_1D.eq.1) then
           iyl= 1
           iyr= 2
        endif

        ! Update potential at boundary
        !$OMP PARALLEL PRIVATE(igrid)
        !$OMP DO
        do iy=iyl,iyr
           do ix=0,n(1)+2
              ! dV/dt= C*dQ/dt : V^k+1= V^k + C*(Q^k+1 - Q^k) (1)
              ! V^1= V_0 + C*Q^1
              ! V^2 = V^1 + C*(Q^2-Q^1)
              !    = V_0 + C*Q^2
              ! ...
              ! V^n = V_0 + C*Q^n (2) : Q^n is the total charge at time k= n
              ! We chose to implement (1) with no loss of generality.
              igrid= bcnd(ix,iy) ! LHS
              if(igrid.ge.1) then
                 ! Ca == F/m2
                 if(dtype(igrid).eq.1) phi(ix,iy)= phi(ix,iy) +  &
                      SUM(sum_q_red_y(igrid,iy,1:ntype))/(Ca*h(2)) 
              endif
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL  
        
     endif

     !
     ! Calculate rho
     !
     if(flag_pbc.eq.1) call np_periodic(n,np,bcnd,ntype,nproc)
     if(opt_solver.eq.1 .or. opt_solver.eq.3) call calc_rhs(n,np,rhs,ntype,nproc,nproc_mpi)
     if(opt_solver.eq.2) call calc_rhs_par(n,np,rhs_par,bcnd,phi,jpl,jpr,NxNy,ntype,nproc,nproc_mpi)
     ctime(1)= ctime(1) + MSTIMER()
     
     !
     ! Calculate Potential
     !  
     if(opt_solver.eq.1) then ! MG
        call iterative_solver(phi,rhs,bcnd,h,n,ncycle,eps,niter,ksor,res,ng,mpi_rank,nproc_mpi)
        
        ! Concatenate phi
        if(it.eq.1) then
           shift(0)=0
           do i=0,nproc_mpi-1
              length(i)=m
              if(i.ge.1) shift(i)=i*m+1
           enddo
           length(0)=length(0)+1
           length(nproc_mpi-1)=length(nproc_mpi-1)+2
           
           length= length*(n(1)+3)
           shift= shift*(n(1)+3)
           
           jl=mpi_rank*m+1
           if(mpi_rank.eq.0) jl=jl-1
           jr=(mpi_rank+1)*m
           if(mpi_rank.eq.nproc_mpi-1) jr=jr+1
        endif
        
        call MPI_ALLGATHERV(phi(0:n(1)+2,jl:jr),length(mpi_rank),MPI_REAL8,phi_tmp, &
             length,shift,MPI_REAL8,MPI_COMM_WORLD,ierr)
        phi= phi_tmp
     endif
     
     if(opt_solver.eq.2) then ! Pardiso
        if(it.eq.1) then
           flag_pardiso=1 ! Initialize
           call directsolver(n(1),n(2),h,jpl,jpr,NxNy,phi,rhs_par,bcnd,flag_pbc,flag_pardiso,mpi_rank)
           flag_pardiso=2 ! Run pardiso
        endif
        call directsolver(n(1),n(2),h,jpl,jpr,NxNy,phi,rhs_par,bcnd,flag_pbc,flag_pardiso,mpi_rank)
     endif
     
     if(opt_solver.eq.3) then
        call tridiag(n,h,phi,rhs)
        phi(:,2)= phi(:,1)
     endif
     
     ctime(2)= ctime(2) + MSTIMER()
     
     !
     ! Calculate electric field
     !
     if(opt_solver.ne.3) then
        call calc_Efield(n,h,phi,Ei,bcnd)
     else
        call calc_Efield1D(n,h,phi,Ei)
        Ei(1,:,2)= Ei(1,:,1) ! Ex
     endif
     
     ctime(1)= ctime(1) + MSTIMER()
     
     !
     ! Sort particles inside vxp()
     !
     if( flag_nopart.eq.0 .and. MOD(it,nsort).eq.1) then
        
        ! Allocate array
        allocate ( sorting(6,nmax*nproc) )
        
        do ptype=1,ntype
           if( SUM(np_tot(ptype,1:nproc)).gt.0 ) then
              call part_sorting(vxp,sorting,nmax,ntype,n,h,Plist,pl_max,nproc,np_tot,ptype,mpi_rank)
           endif
        enddo
        
        ! deallocate
        deallocate ( sorting )
        
        ctime(3)= ctime(3) + MSTIMER()
     endif
     
     !
     ! Calculate averages
     !
     if( flag_nopart.eq.0 .and. MOD(it,navg).eq.1 ) then 

        n_mts=4
        allocate ( p_mts(n_mts,na(1)+1,na(2)+1,ntype,nproc) )
        
        !$OMP PARALLEL PRIVATE(iproc,ptype)
        iproc= omp_get_thread_num() + 1
        
        ! Particle moments
        do ptype=1,ntype
           call part_moments(n,na,h,vxp,nmax,ntype,kq,nproc,np_tot, &
                iproc,ptype,p_mts,n_mts)           
        enddo ! enddo over ptype particle
        !$OMP END PARALLEL

        if(plt_src.eq.1) then
           flag_diag=1
           ss2D= 0.d0
           call collisions(it,vxp,n,h,ntype,nmax,sig,sig_Er,sig_list,sig_Eex,&
                ncol_mx,npt_mx,cnt_col,P_loss,Plist,pl_max,sigv_mx,col_info,&
                np_red,flag_dead,nproc,np_tot,iseed,nproc_mpi,mpi_rank,ss2D,flag_diag)
        endif
        
        call calc_avg(n,h,na,np,p_mts,data_pavg,ntype,n_mts,ss2D,nproc,nproc_mpi)

        if(plt_src.eq.0) ss2D= 0.d0
        deallocate( p_mts )
        
        ! Calculate average value of the potential
        !$OMP PARALLEL
        !$OMP DO
        do iy=0,n(2)+2
           do ix=0,n(1)+2
              phi_avg(1,ix,iy)= phi_avg(1,ix,iy) + phi(ix,iy)
              phi_avg(2,ix,iy)= phi_avg(2,ix,iy) + Ei(1,ix,iy)
              phi_avg(3,ix,iy)= phi_avg(3,ix,iy) + Ei(2,ix,iy)
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        ! Update counter
        cnt_avg(1)= cnt_avg(1) + 1

        ! Calculate average value of uB
        do ptype=2,ntype
           !$OMP PARALLEL PRIVATE(Te)
           !$OMP DO
           do iy=1,na(2)+1
              do ix=1,na(1)+1
                 Te= data_pavg(2,ix,iy,1)/real(cnt_avg(1))
                 uB(ix,iy,ptype)= dsqrt( qe*Te/mass(ptype) )
              enddo
           enddo
           !$OMP END DO NOWAIT
           !$OMP END PARALLEL
        enddo
        
        ctime(4)= ctime(4) + MSTIMER()
     endif

     !
     ! Collisions
     !
     if( it.gt.1 .and. ncol.gt.0 .and. MOD(it,ns_coll).eq.1 ) then
        call dens_red(n,np,np_red,bcnd,ntype,nproc,nproc_mpi)
        ! WARNING: must be right before part_mover() due to flag_dead(), i.e.,
        ! do not insert sorting() in between as flag_dead() is a particle list
        flag_diag=0
        call collisions(it,vxp,n,h,ntype,nmax,sig,sig_Er,sig_list,sig_Eex,&
             ncol_mx,npt_mx,cnt_col,P_loss,Plist,pl_max,sigv_mx,col_info,&
             np_red,flag_dead,nproc,np_tot,iseed,nproc_mpi,mpi_rank,ss2D,flag_diag)
        ctime(5)= ctime(5) + MSTIMER() 
     endif

     !
     ! Move particles
     !
     if( MOD(it,ns_heat).eq.0 .and. Pabs(1).gt.0.d0 ) then
        do i_icp=1,n_icp
           ! Reduction
           sum_dEk_tot= SUM(sum_dEk(i_icp,1:nproc)) ! sum over OMP proc
           if(nproc_mpi.gt.1) then ! sum over MPI proc
              sum_dEk_tmp=0.d0
              call MPI_ALLREDUCE(sum_dEk_tot, sum_dEk_tmp, 1, MPI_REAL8, MPI_SUM, &
                   MPI_COMM_WORLD, ierr)
              sum_dEk_tot= sum_dEk_tmp
           endif

           sum_Nh= SUM(Nh(i_icp,1:nproc))
           if(nproc_mpi.gt.1) then
              sum_Nh_tmp=0
              call MPI_ALLREDUCE(sum_Nh, sum_Nh_tmp, 1, MPI_INTEGER, MPI_SUM, &
                   MPI_COMM_WORLD, ierr)
              sum_Nh= sum_Nh_tmp
           endif

           if(eheat_type.eq.1) then ! update electron velocity (sampled from Maxwellian)
              Te= (2.d0/3.d0)*( sum_dEk_tot + Pabs(i_icp)/nu_h )/(Nm(1)*qe*real(sum_Nh))
           else ! add a small dv incrementally 
              Te= (2.d0/3.d0)*(Pabs(i_icp)/nu_h )/(Nm(1)*qe*real(sum_Nh))
           endif

           ! Calculated thermal velocity
           vt(i_icp)=dsqrt(2.d0*qe*Te/ABS(mass(1)))
        enddo
     endif

     ! Floating end-plates (third dimension)
        if( it.gt.1 .and. flag_float.eq.1 ) then
              sum_q_tmp=0.d0
              call MPI_ALLREDUCE(SUM(sum_q(1:ntype,1:nproc)), sum_q_tmp, 1, MPI_REAL8, MPI_SUM, &
                   MPI_COMM_WORLD, ierr)
              ! Ca == F/m2
              phi0= phi0_sav + sum_q_tmp/(Ca*xmax*ymax*2.d0)
        endif


     !$OMP PARALLEL PRIVATE(iproc,ptype,iseed_OMP,sav_np)
     iproc= omp_get_thread_num() + 1
     iseed_OMP= iseed(iproc)

     ! Initialize counters and arrays
     Nh(:,iproc)=0
     sum_dEk(:,iproc)=0.d0
        
     ! Electron heating
     if( MOD(it,ns_heat).eq.0 .and. flag_heat.eq.1 ) then
        if(flag_inj.eq.1) vt= vt0(1) ! Constant electron temperature
        call eheating(vxp,nmax,ntype,nproc,iseed_OMP,&
             np_tot,vt,n_icp,iproc,P_loss,Pabs_cor)
     endif
     
     ! Push particles
     if(flag_nopart.eq.0) then
        do ptype=1,ntype
           ! Do not consider neutrals in this subroutine
           if( ABS(charge(ptype)).gt.0.d0 ) then
              if(ABS(opt_inj).eq.2) sav_np= SUM(p_mac(ptype,np_loss,0:ngrid,iproc))
              if(flag_1D.eq.0) then
                 call part_mover(it,n,na,h,Ei,Bi,p_mac,P_loss,vxp,bcnd,nmax,ntype,&
                      ngrid,flag_dead,nproc,np_tot,phi,uB,iproc,ptype,iseed_OMP,&
                      cnt_dead,sum_q,sum_q_y,phase_space,ip_ps,nmax_ps,&
                      n_B,h_B,ss2D)
              else
                 call part_mover1D(it,n,h,Ei,Bi,p_mac,P_loss,vxp,&
                      bcnd,nmax,ntype,ngrid,flag_dead,nproc,np_tot,&
                      iproc,ptype,iseed,cnt_dead,sum_q_y,n_B,h_B,ss2D)
                 np(:,2,ptype,iproc)= np(:,1,ptype,iproc)
              endif
              if(ABS(opt_inj).eq.2) N_inj(ptype,iproc)= N_inj(ptype,iproc) + &
                   ( SUM(p_mac(ptype,np_loss,0:ngrid,iproc)) - sav_np )
           endif
        enddo
     endif
     
     ! Particle injection
     if( MOD(it,ns_inj).eq.0 .and. flag_inj.eq.1 ) then
        do ptype=1,ntype
           if( ABS(charge(ptype)).gt.0.d0 ) &          
                call part_injection(n,h,bcnd,vxp,ss2D,ntype,nmax,ni0,I_inj,&
                np_tot,nproc,iseed_OMP,nproc_mpi,N_inj,iproc,ptype,P_loss,phi)
        enddo
        ! Reset counter 
        N_inj(:,iproc)= 0
     endif

     ! Calculate particle density
     do ptype=1,ntype
        ! Initialize particle density array
        np(:,:,ptype,iproc)=0.d0
        if(np_tot(ptype,iproc).gt.0) then
           if(flag_1D.eq.0) then
              call charge_deposition(n,h,vxp,nmax,ntype,kq,np,nproc,np_tot,sum_dEk,&
                   Nh,iproc,ptype,n_icp)
           else
              call charge_deposition1D(n,h,vxp,nmax,ntype,kq,np,nproc,np_tot,sum_dEk,&
                   Nh,iproc,ptype,n_icp)
           endif
        endif
     enddo
          
     iseed(iproc)=iseed_OMP
     !$OMP END PARALLEL
     ctime(6)= ctime(6) + MSTIMER()

     ! Negative ion emmiter on the surface of the extraction electrode 
     if( jne.gt.0.d0 .and. (ns_Hm.eq.1 .or. MOD(it,ns_Hm).eq.1) ) then
        call part_flux_emitter(it,vxp,n,h,ntype,nmax,bcnd,ss2D,nproc,&
             iseed,np_tot,nproc_mpi,mpi_rank,cnt_b,n_icp,P_loss)
        ctime(7)= ctime(7) + MSTIMER()
     endif

     !
     ! Backup simulation data
     !
     if( it.gt.1 .and. MOD(it,nbak).eq.1 ) then 

        if(mpi_rank.eq.0) print*, 'Backing up simulation data ...'
           
        write (pnum_bck,'(i3)'),mpi_rank
        if(mpi_rank.le.9) lgh=3
        if(mpi_rank.ge.10 .and. mpi_rank.le.99 ) lgh=2
        if(mpi_rank.ge.100 .and. mpi_rank.le.999 ) lgh=1

        ! Create backup files
        if(flag_wrt.eq.0) then
           open(41+mpi_rank,file='DATA.BAK/particles'//pnum_bck(lgh:3)//'.bak',form='UNFORMATTED')
           if( flag_dielec.eq.1 .and. mpi_rank.eq.0 ) &
                open(43,file='DATA.BAK/phi.bak',form='UNFORMATTED')
           flag_wrt= 1
        else
           open(41+mpi_rank,file='DATA.BAK2/particles'//pnum_bck(lgh:3)//'.bak',form='UNFORMATTED')
           if( flag_dielec.eq.1 .and. mpi_rank.eq.0 ) &
                open(43,file='DATA.BAK2/phi.bak',form='UNFORMATTED')
           flag_wrt= 0
        endif

        ! Save particles
        write(41+mpi_rank) time,gtype,Pabs_cor
        write(41+mpi_rank) np_tot(1:ntype,1:nproc)
        do iproc=1,nproc
           do ptype=1,ntype
              write(41+mpi_rank) vxp(:,1:np_tot(ptype,iproc),ptype,iproc)
           enddo
        enddo
        close(41+mpi_rank)

        ! Save dielectric potential
        if( flag_dielec.eq.1 .and. mpi_rank.eq.0 ) then
           write(43) phi(0:n(1)+2,0:n(2)+2)
           close(43)
        endif

        if(mpi_rank.eq.0) print*, 'Done!'
        
        !
        ! Save ion properties on num_grd
        !
        if(flag_ps.eq.2) then 
           if(mpi_rank.eq.0) then
              print*, 'Saving phase space distribution ...'
              open(41,file='DATA/beam_phase_space.dat')
              write(41,*) '# y(m), yp(rad)'
              do iproc=1,nproc
                 do i=1,ip_ps(iproc)
                    write(41,'(3(es14.6,1x),(i2,1x))') ( phase_space(i,j,iproc), j=1,3 ), &
                         NINT(phase_space(i,4,iproc))
                 enddo
              enddo
              close(41)
           endif

           flag_ps=3
           if(mpi_rank.eq.0) print*, 'Done!'
        endif

        !
        ! Save particle distribution function
        !
        if(flag_pdf.eq.1) then
           if(mpi_rank.eq.0) then 

              print*, 'Saving particle distribution function ...'
           
              do i_rg=1,n_rg
                 write (pnum,'(i1)'),i_rg
                 open(41+i_rg,file='DATA/xv'//pnum//'.dat')
              enddo
              
              do iproc=1,nproc
                 do i=1,np_tot(ptype_pdf,iproc)
                    do i_rg=1,n_rg
                       
                       ! Save velocity for each predefined sub-regions
                       if ( vxp(1,i,ptype_pdf,iproc).ge.xl_rg(i_rg) .and. &
                            vxp(1,i,ptype_pdf,iproc).le.xr_rg(i_rg) .and. &
                            vxp(2,i,ptype_pdf,iproc).ge.yl_rg .and. &
                            vxp(2,i,ptype_pdf,iproc).le.yr_rg ) then  
                          write(41+i_rg,'(3(es14.6,1x))') vxp(4,i,ptype_pdf,iproc),&
                               vxp(5,i,ptype_pdf,iproc),vxp(6,i,ptype_pdf,iproc)
                          goto 30
                       endif
                       
                    enddo
30                  continue
                 enddo
              enddo
           
              do i_rg=1,n_rg
                 write (pnum,'(i1)'),i_rg
                 close(41+i_rg)
              enddo
              
              print*, 'Done!'
           endif
        endif
        
        ctime(8)= ctime(8) + MSTIMER()
        
     endif
     
     !
     ! Save data
     !
     if( MOD(it,nsav).eq.1 ) then

        ! Sum np_tot over all processes
        sum_np_tot= SUM(np_tot,DIM=2) ! sum over OMP proc
        if(nproc_mpi.gt.1) then ! sum over MPI proc
           sum_np_tot_tmp=0
           call MPI_ALLREDUCE(sum_np_tot, sum_np_tot_tmp, ntype, MPI_INTEGER, MPI_SUM, &
                MPI_COMM_WORLD, ierr)
           sum_np_tot= sum_np_tot_tmp
        endif

        if(mpi_rank.eq.0) then

           ! Save energy conservation, particle errors & calculation time
           sum_time= real(SUM(ctime(1:10)))
           if(it.eq.1) write(12,'(a105)') &
                '# time(s), nsteps, ctime{E/rho, poisson, sorting, avg/write, MC, mover, &
                neg. ions, bck}(%), total(ms)'
           write(12,100) time,it,( real(ctime(i))/sum_time*100., i=1,8 ), &
                real(SUM(ctime(1:10)))/real(it)
100        format(es15.8,2x,i8,8(2x,es10.2),(2x,f7.1))
           
           ! Print on screen
           write(*,101) it,time*1.d6, (sum_np_tot(ptype), ptype=1,ntype)
101        format(' it= ',i8,', t (us)= ',f6.1,', np= ',10(i9,2x))
           
           ! Chamfered apertures
           if( flag_chmfrd.eq.1 .and. tag_b.gt.0 ) &
                print'(1x,"Chamfered over front grid surface current ratio=",f6.2)', &
                real(SUM(cnt_b(2,1:nproc)))/real(SUM(cnt_b(1,1:nproc)))
           
           ! Save PDE solver residual & # of iterations
           if(opt_solver.eq.1) write(40,*) res,niter,nint(ksor)          
           
        endif
       
        cnt_avg(2)= nsav
        if( it.eq.1 .or. flag_sav.eq.1 ) cnt_avg(2)= it
        flag_bak=0
        if( flag_sav.eq.0 .and. tseq.gt.0.d0 .and. &
             time.ge.tseq_init .and. time.le.tseq_final) then
           ! Save sequence of profiles in Macho format.
           if(MOD(it,nseq).eq.1 ) flag_bak=2 
        endif
        call write_data(it,time,n,na,h,p_mac,P_loss,phi_avg,data_pavg,ntype,ngrid,&
             cnt_col,ncol_mx,sig_list,nproc,cnt_avg,mpi_rank,nproc_mpi,Pext,Vgrd,I_inj,n_B)

        ! Floating end-plates
        if( flag_float.eq.1 .and. mpi_rank.eq.0 ) print'(1x,"End-plates potential (V):",f8.2)',phi0

        ! Negative ion current lost through collisions
        if( it.gt.1 .and. tag_neg.gt.0 ) then 
              cnt_dead_tmp=0.d0
              call MPI_REDUCE(SUM(cnt_dead,DIM=1), cnt_dead_tmp, 1, MPI_REAL8, MPI_SUM, &
                   0,MPI_COMM_WORLD, ierr)
              if( mpi_rank.eq.0 .and. ABS(cnt_dead_tmp).gt.0.d0 ) &
                   print'(1x,"Neg. ion current lost through collisions (A)=",es10.2)', &
                   cnt_dead_tmp/(real(cnt_avg(2))*dt)
        endif

        ! Size of phase space macroparticle sample 
        if( it.gt.1 .and. flag_ps.eq.1 ) then 
           if( mpi_rank.eq.0 ) &
                   print'(1x,"Size of phase space macroparticle sample:",i6)', &
                   SUM(ip_ps(1:nproc))
        endif

        ! Plot current profile on dielectric
        if( flag_dielec.eq.1 .and. mpi_rank.eq.0 ) then 
           do igrid= 1,ngrid
              if(dtype(igrid).eq.1 ) then
                 write (pnum,'(i1)'),igrid
                 open(15,file='DATA/jg'//pnum//'_diel.dat')
                 do iy=1,n(2)+1
                    write(15,*) (iy-1)*h(2),( sum_I_sav(igrid,iy,ptype)/(h(2)*real(cnt_avg(2))), ptype=1,ntype )
                 enddo
                 close(15)
              endif
           enddo
        endif

        ! Reset arrays
        if(flag_sav.eq.0) then
           phi_avg=0.d0
           data_pavg=0.d0
           cnt_avg(1)=0
           p_mac=0.d0
           P_loss=0.d0
           cnt_dead=0.d0
           cnt_col(1:2,:,:,:)=0.d0
           if(flag_dielec.eq.1) sum_I_sav= 0.d0
        endif
        
        ! Adjust the heating frequency - Gaussian power profile
        if ( flag_gpp.eq.1 .and. Pext.gt.0.d0 ) then 
           Pabs_cor= Pabs_cor*(Pabs(1)/Pext)
           ! Warning
           if((nudt*9.d0*Pabs_cor).gt.1) then
              print*, 'nu0*dt>1 for Gaussian power profile, please correct...'
              call stop_calculation
           endif
        endif

        ! Time lag for writing data and calculating averages
        ctime(4)= ctime(4) + MSTIMER()

        if( it.gt.1 .and. mpi_rank.eq.0 ) &
             write(*,103) ( real(ctime(i))/real(it), i=1,6 ),&
             real(ctime(8))/real(it), &
             real(SUM(ctime(1:10)))/real(it)
103     format(' <t> (ms): E/rho= ',f7.1,', poisson= ',f7.1,', sorting= ',f7.1, &
             ', avg/write= ',f7.1,', MC= ',f7.1,', mover= ',f7.1,', bck= ',f7.1,', total= ',f7.1)

     endif

  enddo ! enddo over ntmax iterations

  ! 
  ! Close external files
  !
  close(12)
  close(40)

  call MPI_Finalize(ierr)

end program main

subroutine introduction
!     ==============================================================
!     VERSION:         4.7.4
!     LAST MOD:       Sep/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Display code info
!     NOTE:              /
!     --------------------------------------------------------------
  implicit none

  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) '+                                                                          +'
  write(*,*) '+ 2.5D-3V explicit parallel particle in cell code                          +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ LePIC 2.5D algorithm version: 4.7.4                                      +'
  write(*,*) '+ Last modifications: Sep/2024                                             +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ Copyright or © or Copr. Gwenael Fubiani (2022/11/22)                     +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ gwenael.fubiani@cnrs.fr                                                  +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ This software is a computer program whose purpose is to model low        +'
  write(*,*) '+ temperature plasmas with a Particle-In-Cell algorihm in 2.5D-3V          +'
  write(*,*) '+ dimensions.                                                              +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ This software is governed by the CeCILL-C license under French law and   +'
  write(*,*) '+ abiding by the rules of distribution of free software.  You can  use,    +'
  write(*,*) '+ modify and/ or redistribute the software under the terms of the CeCILL-C +'
  write(*,*) '+ license as circulated by CEA, CNRS and INRIA at the following URL        +'
  write(*,*) '+ "http://www.cecill.info".                                                +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ As a counterpart to the access to the source code and  rights to copy,   +'
  write(*,*) '+ modify and redistribute granted by the license, users are provided only  +'
  write(*,*) "+ with a limited warranty and the software's author, the holder of the     +"
  write(*,*) '+ economic rights,  and the successive licensors  have only  limited       +'
  write(*,*) '+ liability.                                                               +'
  write(*,*) '+                                                                          +'
  write(*,*) "+ In this respect, the user's attention is drawn to the risks associated   +"
  write(*,*) '+ with loading,  using,  modifying and/or developing or reproducing the    +'
  write(*,*) '+ software by the user in light of its specific status of free software,   +'
  write(*,*) '+ that may mean  that it is complicated to manipulate,  and  that  also    +'
  write(*,*) '+ therefore means  that it is reserved for developers  and  experienced    +'
  write(*,*) '+ professionals having in-depth computer knowledge. Users are therefore    +'
  write(*,*) "+ encouraged to load and test the software's suitability as regards their  +"
  write(*,*) '+ requirements in conditions enabling the security of their systems and/or +'
  write(*,*) '+ data to be ensured and,  more generally, to use and operate it in the    +'
  write(*,*) '+ same conditions as regards security.                                     +'
  write(*,*) '+                                                                          +'
  write(*,*) '+ The fact that you are presently reading this means that you have had     +'
  write(*,*) '+ knowledge of the CeCILL-C license and that you accept its terms.         +'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) '                                                 '
  
  return
end subroutine introduction
