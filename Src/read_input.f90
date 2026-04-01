subroutine read_input(n,na,tmax,nsav,eps,kt,rname,ngrid,ng, &
     n_icp,xl_rg,xr_rg,yl_rg,yr_rg,I_inj,np_dup,mpi_rank,Ca,n_B,&
     phi0_RF,f0_RF,phi1_RF,f1_RF)
!     ==============================================================
!     VERSION:         0.6
!     LAST MOD:      Sep/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  integer:: iB,n(3),na(3),nsav,flag_read,ngrid,i_rg,i_icp,n_icp,ng,&
       mpi_rank,opt_endplts,np_dup,n_B(2)
  real(kind=8):: eps,tmax,kt,xl_rg(nm_rg),xr_rg(nm_rg), &
       yl_rg,yr_rg,I_inj,Ca,phi0_RF,f0_RF,phi1_RF,f1_RF
  character:: end_file*3,ans1*1,ans2*1,rname*20

  ! Initialize variables & arrays
  xl_rg=0.d0
  xr_rg=0.d0
  flag_grd=0

  ! Open input file
  open(10,file='input_dir/conditions.inp')

  read(10,*,end=999) rname,flag_ionrep,fEk ! name of file storing reactions, ionization scheme, Ek ratio between elec.

  read(10,*,end=999) Ti(1) ! read electron temperature (eV)

  !
  ! Simulation parameters
  !
  read(10,*,err=999) ng  ! # of levels in MG
  read(10,*,err=999) eps ! Min. precision in residual
  opt_solver=1
  read(10,*,err=999) n(1),n(2),opt_solver ! n: nx=2**n, ny=2**n. Number of grid points, opt_solver: 1 (MG), 2 (Pardiso)
  flag_1D= 0
  if(n(2).eq.1) then 
     flag_1D= 1
     opt_solver=3 ! Thomas algorithm
  endif
  
  read(10,*,err=999) xl_pow,xr_pow,yl_pow,yr_pow ! particle heating region [cm]
  xl_pow= xl_pow*1.d-2 ! convert to meters
  xr_pow= xr_pow*1.d-2
  yl_pow= yl_pow*1.d-2
  yr_pow= yr_pow*1.d-2
  flag_c=0 ! Slit
  if( (xl_pow.lt.0.d0 .or. xr_pow.lt.0.d0) .and. flag_1D.eq.0 ) then
     flag_c=1 ! Disk
  endif
  xl_pow= ABS(xl_pow)
  xr_pow= ABS(xr_pow)
  xa= (xl_pow+xr_pow)/2.d0
  dr= (xr_pow-xl_pow)/2.d0

  read(10,*,err=999) x_load ! particle loading region [cm]
  x_load= x_load*1.d-2

  read(10,*,err=999) nB ! # of magnetic field files
  flag_B_pos=0
  if(nB.lt.0) flag_B_pos=1
  nB= ABS(nB)
  read(10,*,err=999) n_B(1),n_B(2) ! Grid for the magnetic field map
  if( n(2).eq.1 ) n_B(2)= 1 ! 1D option
  flag_B=0
  do iB=1,nB
     ! Ext. file (y/n/s), name, scaling, B(G), direction, L/x0/y0 (cm)
     read(10,*,err=999) B_file(iB),B_name(iB),B_scale(iB),B0(iB),B_info(iB),&
          dL(iB),x0(iB),y0(iB)
     B0(iB)= B0(iB)*1.d-4 ! convert into Tesla
     if(B_file(iB).ne.'s' .or. B_file(iB).ne.'S') flag_B=1
  enddo
  if( n_B(1).le.1 .or. flag_B.eq.0 .or. nB.eq.0 ) n_B= 1  ! Option for B=cste
  if(flag_B.eq.0) nB=0
  dL= dL*1.d-2
  x0= x0*1.d-2
  y0= y0*1.d-2
  
  ! # of drivers (<0 = iterative P), Pabs(W) (n times, <0 for a fixed Te), heat. freq. (<0 == off)
  read(10,*,err=999) n_icp,(Pabs(i_icp),i_icp=1,ABS(n_icp)),nu_h
  flag_convP=1
  if(n_icp.eq.-1) flag_convP=0
  flag_icp=0
  if(n_icp.ge.2) flag_icp=1
  n_icp= ABS(n_icp)
  if(flag_1D.eq.1) n_icp=1
  flag_heat=1
  if(nu_h.lt.0.d0) flag_heat=0
  ! Inject current I(A) (<0 == off), options (1=fixed, 2= equal to ion losses, 3= beam along Z, <0= gaussian)
  read(10,*,err=999) I_inj, opt_inj  
  flag_inj=1
  if(I_inj.lt.0.d0) then 
     flag_inj=0
     I_inj= 0.d0
  else
     if(ABS(opt_inj).eq.2) I_inj= 0.d0
  endif
 
  read(10,*,err=999) tmax ! s
  read(10,*,err=999) kt ! time coeff. kt ; dt= kt*dx/vt
  read(10,*,err=999) n0 ! m-3
  read(10,*,err=999) ngas ! m-3
  read(10,*,err=999) np_cell ! number of particles per cell
  read(10,*,err=999) ngrid ! number of wall labels
  read(10,*,err=999) nsav ! Frequency for data saving
  read(10,*,err=999) nbak,cnt_plt,tseq_init,tseq_final ! Frequency for simulation backup
  tseq= ABS(tseq_final-tseq_init)/real(cnt_plt)
  read(10,*,err=999) size_na ! Size ratio of arrays used for saving data
  if(flag_1D.eq.1) size_na=1
  na=n/size_na
  read(10,*,err=999) k_eps0
  plt_src=0
  if(k_eps0.lt.0) then
     ! Flag for plotting source/sink terms (0: count macroparticle production, 1: n*nu is evaluated)
     plt_src=1
     k_eps0= ABS(k_eps0)
  endif
     
  ! Save pdf (y/n), ptype, yl & yr, # of regions n_rg, xl & xr for each n_rg
  read(10,*,err=999) ans1,ptype_pdf,yl_rg,yr_rg,n_rg,(xl_rg(i_rg),xr_rg(i_rg),i_rg=1,n_rg) 
  yl_rg= yl_rg*1.d-2
  yr_rg= yr_rg*1.d-2
  xl_rg= xl_rg*1.d-2
  xr_rg= xr_rg*1.d-2

  if( n_rg.gt.nm_rg .and. mpi_rank.eq.0 ) then ! Warning
     print*, 'Warning nm_rg is too small, please correct ...'
     call stop_calculation
  endif

  flag_pdf=0
  if(ans1.eq.'y'.or.ans1.eq.'Y') flag_pdf=1
   
  read(10,*,err=999) ans1,omp_rank_max,mpi_rank_max,np_dup ! restart (yes/no/external/gauss),  # of OMP threads, # of MPI, # for particle duplication
  flag_restart=0
  if(ans1.eq.'y'.or.ans1.eq.'Y') flag_restart=1
  if(ans1.eq.'e'.or.ans1.eq.'E') flag_restart=2
  if(flag_restart.eq.0) np_dup=1

  
  read(10,*,err=999) ans1,xg1,xg2,Lg,Lh,n_holes,ans2,ind_g
  if(ans1.eq.'y'.or.ans1.eq.'Y') flag_grd(1)=1
  if(ans1.eq.'o'.or.ans1.eq.'O') flag_grd(1)=2
  if(ans2.eq.'y'.or.ans1.eq.'Y') flag_grd(2)=1

  read(10,*,err=999) jne,ptype_fx,THm,num_grd
  flag_Tp=0
  if(THm.lt.0) then
     flag_Tp=1
     THm= ABS(THm)
  endif
  flag_ps=0
  if(num_grd.lt.0) then
     flag_ps=1
     num_grd= ABS(num_grd)
  endif
  read(10,*,err=999) hz,zmax,phi0,opt_endplts,Ca ! h=ns/<n> [<0 length=1m and no losses], length in 3rd dim (cm), bias(V), opts (1=metal, 2=float), Capacitance (F/m2)
  flag_float=0
  if(opt_endplts.eq.2) flag_float=1
  zmax= ABS(zmax)*1.d-2
  if(hz.lt.0 .or. flag_1D.eq.1) zmax=1.d0 ! per meter or m^2

  read(10,*,err=999) gam_sec,igrid_sec,ans1,phi0_RF,f0_RF,phi1_RF,f1_RF ! secondary emission coeff, grid index (<0, fixed Is instead), phi0(V), f0(Hz), phi1(V), f1(Hz)
  
  flag_RFpot=0
  if(ans1.eq.'y'.or.ans1.eq.'Y') flag_RFpot=1
  
  ! Warning
  if( gam_sec.gt.0.d0 .and. jne.gt.0.d0 .and. ind_g.eq.igrid_sec ) then
     if(mpi_rank.eq.0) print*, 'gam_sec>0 and jne>0 options incompatible on the same wall index, please correct...'
     call stop_calculation
  endif
  
  flag_sec=0
  if(gam_sec.gt.0.d0) then
     ! Fixed secondary emission coefficient gamma
     flag_sec=1
     ! Fixed I_sec= |I_inj| with adaptative gamma
     if(igrid_sec.lt.0) flag_sec=2
  else
     gam_sec= 0.d0
  endif
  igrid_sec= ABS(igrid_sec)
  
  flag_read=0
  do while( flag_read.eq.0 )
     read(10,*,err=999) end_file
     if( end_file.eq.'END' .or. end_file.eq.'end'  ) flag_read=1
  enddo
 
  if(flag_read.eq.1) then
     if(mpi_rank.eq.0) print*, 'Input file was read correctly ...'
     close(10)
     return
  endif

999 if(mpi_rank.eq.0) print*, 'Input file was not read correctly!'  
  call stop_calculation

end subroutine read_input
