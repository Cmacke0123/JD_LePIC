!
! Counters and seed number
!
  integer:: npart,np_cell,flag_pdf, &
       flag_restart,n_cell,nbak,n_rg, &
       flag_pbc,cnt_seg,tag_neg,num_grd,&
       tag_neu,flag_heat,flag_inj,flag_nmn,&
       opt_inj,flag_B_pos,flag_convP,eheat_type,&
       flag_gpp,flag_dielec,mpi_rank_max,omp_rank_max,&
       flag_chmfrd,tag_b,flag_Tp,cnt_plt,flag_bak,flag_spec,&
       flag_float,dtype(10),opt_solver,flag_ps,igrid_sec,dir_sec,&
       flag_sec,flag_ionrep,flag_1D,flag_diffsrc,flag_icp,flag_gridB,&
       plt_src,flag_RFpot
  parameter (npart=32)

!
! Collisions
!
  integer:: ncol,p_ncol(npart),sig_npt_mx,nscol,p_nscol(npart)
  real(kind=8):: np_mx(npart),nu_max(npart),nu_uplim(npart),&
            gam_sec,xg_sec

!
! Particle arrays & flags
!
  integer:: np_loss,P_w,ind_nre,ind_nby,ind_Eth, &
            ind_dE,cnt_neg
  parameter ( np_loss=1, P_w=2, ind_nre=1, &
       ind_nby=2, ind_Eth=1, ind_dE=2 )
  real(kind=8):: charge(npart),mass(npart),vt0(npart),Ti(npart),phi0,chmfrd_array(4)

!
! Magnetic field
!
  integer:: nB
  real(kind=8):: B_scale(5),B0(5),dL(5),x0(5),y0(5),th_B
  character:: B_file(5)*1,B_name(5)*20,B_info(5)*2

!
! Particle and simulation info
!
  integer:: ns_heat,ns_coll,gtype,nm_rg,n_holes,flag_grd(2),&
          ind_g,ptype_pdf,ns_Hm,flag_B,nm_icp,nx_PE,ns_inj,flag_c,&
          ptype_fx
  parameter ( nm_rg=5, nm_icp=4 )
  integer:: ixl_rg(nm_rg),ixr_rg(nm_rg),iyl_rg,iyr_rg,size_na
  real(kind=8):: Pabs(nm_icp),n0,ngas,Nm(npart),dt,xmax,ymax,zmax,x_load, &
            nu_h,Lg,Lh,xg1,xg2,Bmax,nudt,hz,xa,dr,tseq,tseq_init,&
            tseq_final,fEk,Lgrd
  real(kind=8):: k_eps0,jne,RNeta,THm,xl_pow,xr_pow,yl_pow,yr_pow
  character:: pname(npart)*6

!
! Common blocks
!
  common /part_info_int/ np_cell,ncol,sig_npt_mx,flag_pdf,n_cell,nbak,&
       flag_restart,flag_B,ns_heat,ns_coll,p_ncol,nscol,&
       p_nscol,gtype,n_rg,flag_pbc,n_holes,flag_grd,ind_g,nB,cnt_seg,tag_neg,&
       cnt_neg,tag_neu,num_grd,ptype_pdf,ns_Hm,ixl_rg,ixr_rg,iyl_rg,iyr_rg,&
       nx_PE,size_na,flag_heat,flag_inj,ns_inj,flag_nmn,opt_inj,&
       flag_B_pos,flag_c,flag_convP,eheat_type,flag_gpp,flag_dielec,&
       mpi_rank_max,omp_rank_max,flag_chmfrd,tag_b,flag_Tp,cnt_plt,&
       flag_bak,flag_spec,flag_float,dtype,opt_solver,flag_ps,igrid_sec,&
       dir_sec,flag_sec,flag_ionrep,ptype_fx,flag_1D,flag_diffsrc,flag_icp,&
       flag_gridB,plt_src,flag_RFpot
  common /part_info_dp_1/ charge,mass,Ti,vt0,k_eps0,B_scale, &
       B0,dL,x0,y0,jne,RNeta,THm,phi0,xl_pow,xr_pow,chmfrd_array,yl_pow,yr_pow,&
       th_B
  common /part_info_dp_2/ Pabs,n0,ngas,Nm,dt,xmax,ymax,zmax,x_load, &
       nu_h,Lg,Lh,xg1,xg2,Bmax,nudt,hz,xa,dr,tseq,tseq_init,&
       tseq_final,gam_sec,xg_sec,fEk,Lgrd
  common /part_info_dp_3/ np_mx,nu_max,nu_uplim
  common /part_info_ch/ pname,B_file,B_name,B_info
