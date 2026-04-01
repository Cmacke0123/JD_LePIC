subroutine collisions(it,vxp,n,h,ntype,nmax,sig,sig_Er,sig_list,sig_Eex, &
     ncol_mx,npt_mx,cnt_col,P_loss,Plist,pl_max,sigv_mx,col_info,&
     np_red,flag_dead,nproc,np_tot,iseed,nproc_mpi,mpi_rank,ss2D,flag_diag)
!     ===================================================================
!     VERSION:         0.8
!     LAST MOD:      Nov/23
!     MOD AUTHOR:    G. Fubiani
!
!     COMMENTS:     
!     NOTES:     (1) When sorting collision frequencies, indx(1) gives the 
!                    index of the smallest frequency for the list of p_ncol 
!                    reactions associated with incident particle ptype. 
!                    The index of the associated cross section as found in 
!                    the input file is sig_list(ptype,indx(#))
!
!                (2) For any type of inelastic collisions we assume: 
!                    |p1|=|p2|=...=|pn|,
!                    where |p1|=m1*v1. We then have
!                    Ekr'= 0.5 mu v'r**2= 0.5 mu vr**2 - Eth 
!                    and, 
!                    0.5*m1*v1**2 + 0.5*m2*v2**2 + ... + 0.5*mn*vn**2 = Ekr'
!                    We obtain 
!                    0.5*m1*v1**2= (1/m1*Ekr')/(1/m1+1/m2+...1/mn),
!                    0.5*m2*v2**2= (1/m2*Ekr')/(1/m1+1/m2+...1/mn), 
!                    etc ...
!                    We now know |v1|,|v2|, ...,|vn|
!                    We then assume equal angle differences between momentum
!                    vectors p: taking th1=rnd*2*pi and 
!                    th2=th1+2*pi/n, 
!                    th3=th2+2*pi/n, 
!                    together with an angle phi identical of all particles
!                    etc ...
!                    For instance with n=3 and th1=10° then
!                    th2=130° and th3=250° with phi=20°.
!                    we deduce the coordinates of v1 as
!                    e1x= cos(th1)
!                    e1y= sin(th1)*dsin(phi)
!                    e1z= sin(th1)*dcos(phi)
!                    and,
!                    v1x= vx_cm + v1*e1x,
!                    v1y= vy_cm + v1*e1y,
!                    v1z= vz_cm + v1*e1z,
!                    in the laboratory frame.
!                    Note than |m1v1|=|m2v2|= ...  so
!                    Sum ex=0, Sum ey=0 and Sum ez=0
!                    
!                (3) For dissociation: provide in CM equal angle separately for
!                    the electrons and for the neutral/ion momentums. Share CM
!                    energy minus threshold between all of them and then add
!                    the extra energy associated with the dissociative process
!                    to the heavy particles in the center of mass
!
!     -------------------------------------------------------------------
  use omp_lib
  implicit none
  include 'mpif.h'
  include 'particle_info.h'
  include 'constants.h'
  ! MPI
  integer ierr,nproc_mpi,mpi_rank,flag_diag
  ! Collision parameters, arrays & variables
  integer nmax,ptype,n(3),ntype,iproc,nproc,Nc(ntype,nproc), &
       np_add(ntype,nproc),Nc_tmp,iseed(nproc), &
       np_tot(ntype,nproc),sum_np_tot,sum_np_tot_tmp,&
       np_tot_rg(ntype,nm_rg),i_rg,err_coll(ntype,nproc),iseed_OMP,it
  integer:: npt_mx,icol,ncol_mx,sig_list(npart,ncol_mx),col_info(ncol_mx,10),n_re,&
       ttype,ind_col,sav_ttype(ncol_mx)
  integer(kind=1):: flag_dead(nmax,ntype,nproc)
  real(kind=8):: h(3),sig(npt_mx,ncol_mx),sig_Er(npt_mx),sig_Eex(ncol_mx,2),&
       sigv_mx(npart,ncol_mx),ran2,rnd,k,Pmax,err_coll_pc,nu_max_OMP(nproc),dNc
  ! Particle arrays & variables
  integer:: pl_max,Plist(0:pl_max,ntype,nproc)
  real(kind=8):: vxp(6,nmax,ntype,nproc), &
       np_red(0:n(1)+2,0:n(2)+2,ntype)
  ! Macroscopic parameters
  real(kind=8):: cnt_col(3,ncol_mx,nproc,nm_rg),P_loss(4,ntype,nproc),&
       ss2D(2,0:n(1)+2,0:n(2)+2,ntype,nproc)

  !
  ! Set collision parameters
  !

  ! Initialize variables & arrays
  np_add= 0
  nu_max= 0.d0
  Nc= 0
  sum_np_tot= 0
  sum_np_tot_tmp= 0
  err_coll=0

  ! Loop over particle species
  do ptype=1,ntype

     ! Jump to end of the loop for collisionless particles or negative masses
     if( p_ncol(ptype).eq.0 .or. mass(ptype).le.0 ) goto 100

     ! Calculate nu_max & number of collisions
     do icol=1,p_ncol(ptype)
        ind_col= sig_list(ptype,icol)
        n_re= col_info(ind_col,1)
        ttype= col_info(ind_col,ind_nby+n_re)
        sav_ttype(ind_col)= ttype
        nu_max(ptype)= nu_max(ptype) + np_mx(ttype)*sigv_mx(ptype,ind_col)
     enddo

     nu_max(ptype)= MIN(nu_max(ptype),nu_uplim(ptype))
     sum_np_tot= SUM(np_tot(ptype,1:nproc))

     if(nproc_mpi.gt.1) then
        call MPI_ALLREDUCE(sum_np_tot, sum_np_tot_tmp, 1, MPI_INTEGER, MPI_SUM, &
             MPI_COMM_WORLD, ierr)
        sum_np_tot= sum_np_tot_tmp
     endif

     ! Total number of collisions (summed over all threads for better statistics)
     Pmax= nu_max(ptype)*real(ns_coll)*dt

     ! Warning
     if(Pmax.gt.1.d0) then 
        print*, 'Collision probability greater than 1, please correct ...'
        print*, 'Origin: "collisions" subroutine'
        call stop_calculation
     endif

     ! Number of collisions per MPI and OMP thread
     dNc= sum_np_tot*Pmax/real(nproc_mpi*nproc)
     Nc_tmp= INT(dNc)

     ! Correct for round-off errors
     rnd= ran2(iseed(1))
     if( rnd.le.(dNc-Nc_tmp) ) Nc_tmp= Nc_tmp + 1
     Nc(ptype,:)= Nc_tmp

     ! Warning
     if( SUM(np_tot(1,1:nproc))/n_cell.ge.5 .and. ptype.eq.1 .and. Nc_tmp.le.1 ) then
        if(mpi_rank.eq.0) then
           print*, 'Collision time step too small for electrons, dNc=',dNc
           print*, 'Please correct; abort simulation ...'
        endif
        call stop_calculation
     endif

     ! Calculate total # of ptype particles per sub-regions
     k=h(1)*h(2)*zmax/Nm(ptype)
     if(flag_1D.eq.1) k=h(1)/Nm(ptype)
     do i_rg=1,n_rg
        if(flag_1D.eq.0) then 
           np_tot_rg(ptype,i_rg)= k*SUM(np_red(ixl_rg(i_rg):ixr_rg(i_rg), &
                iyl_rg:iyr_rg,ptype))
        else
           ! 1D case
           np_tot_rg(ptype,i_rg)= k*SUM(np_red(ixl_rg(i_rg):ixr_rg(i_rg),1,ptype))
        endif
     enddo

100  continue
  enddo ! end-loop over ptype particles

  !
  ! Perform loop in parallel
  !
  !$OMP PARALLEL PRIVATE(iproc,ptype,iseed_OMP)
  ! Get processor id (from 0 to nproc-1)
  iproc= omp_get_thread_num() + 1
  iseed_OMP= iseed(iproc)

  do ptype=1,ntype
     if( p_ncol(ptype).eq.0 .or. mass(ptype).le.0 .or. &
          Nc(ptype,iproc).eq.0 .or. np_tot(ptype,iproc).eq.0 ) goto 110

     ! Correct nu_max for rounding errors
     nu_max_OMP(iproc)= real(Nc(ptype,iproc))/real(np_tot(ptype,iproc))/(real(ns_coll)*dt)
    
     call collision_OMP(vxp,n,h,ntype,nmax,sig,sig_Er,sig_list,sig_Eex,ncol_mx,&
          npt_mx,cnt_col,P_loss,col_info,np_red,flag_dead,nproc,np_tot,iseed_OMP,&
          Plist,pl_max,np_add,Nc,ptype,sav_ttype,iproc,np_tot_rg,nu_max_OMP,err_coll,ss2D,&
          flag_diag)

110 enddo ! end-loop over ptype particles
  iseed(iproc)= iseed_OMP
  !$OMP END PARALLEL

  ! Add new particles to particle counter
  np_tot= np_tot + np_add

  ! Warning
  if( MAXVAL(np_tot).gt.nmax ) then
     print*, 'nmax is too small, please correct'
     print*, 'Abort calculation ...'
     call stop_calculation
  endif

  ! Warning (DSMC)
  if(MOD(it,ns_coll*250).eq.1) then
     do ptype=1,ntype
        if(SUM(Nc(ptype,1:nproc)).gt.0) then
           err_coll_pc= real(SUM(err_coll(ptype,1:nproc)))/real(SUM(Nc(ptype,1:nproc)))
           if( err_coll_pc.gt.0.1 .and. mpi_rank.eq.0 ) then
              print*, '>10% of collisions are missed due to an insufficient number of ppc'
              print'(1x,"Err(%)= ",f8.2,", ptype=",i3)',err_coll_pc*100.d0,ptype
           endif
        endif
     enddo
  endif

  return
end subroutine collisions

subroutine collision_OMP(vxp,n,h,ntype,nmax,sig,sig_Er, &
     sig_list,sig_Eex,ncol_mx,npt_mx,cnt_col,P_loss,col_info,np_red,&
     flag_dead,nproc,np_tot,iseed,Plist,pl_max,np_add,Nc,ptype,&
     sav_ttype,iproc,np_tot_rg,nu_max_OMP,err_coll,ss2D,flag_diag)
!     ===================================================================
!     VERSION:         0.7
!     LAST MOD:      Nov/23
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:     
!     NOTES:   
!     -------------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  ! Collision parameters, arrays & variables
  integer ic,ip,it,ix,iy,ib,nmax,ptype,n(3),ntype,iproc,nproc,&
       Nc(ntype,nproc),flag_coll,c_ind,np_add(ntype,nproc),flag_ttype(npart),&
       sav_it(npart),ind_ce,iseed,np_tot(ntype,nproc),&
       np_tot_rg(ntype,nm_rg),i_rg,err_coll(ntype,nproc),ici,ict, &
       flag_err(npart),Nc_tmp,flag_diag
  integer:: ipt,npt,npt_mx,icol,ncol_mx,sig_list(npart,ncol_mx),&
       col_info(ncol_mx,10),n_re,n_by,ttype,ind_col,sav_ttype(ncol_mx),&
       indx(ncol_mx),i_by,i_re,btype,ctype,flag_add(10),opt_add, &
       tproc,cnt_proc,bproc,sav_tproc(npart),cnt_ppc,ipt_L,ipt_R,ipt_M,&
       opt_scat,cnt_ele,flag_ion
  integer(kind=1):: flag_dead(nmax,ntype,nproc)
  real(kind=8):: h(3),sig(npt_mx,ncol_mx),sig_Er(npt_mx),Ek_L,Ek_R,sig_L, &
       sig_R,sig_p,nu(ncol_mx),sig_Eex(ncol_mx,2),sort_arr(ncol_mx),Ek_M, &
       vzb_sav,vbx,vby,vbz,Ek_b(npart),ki(4),k1,k2
  real(kind=8):: ran2,rnd(2),sum_nu,cos_th,sin_th,th,phi,ex,ey,ez, &
       Eth,Ee,rt,th_add,nu_max_OMP(nproc),cos_ph,sin_ph,ex1,ey1,ez1,cos_th_s,phi_s
  ! Particle arrays & variables
  integer:: pl_max,Plist(0:pl_max,ntype,nproc)
  real(kind=8):: vxp(6,nmax,ntype,nproc),vx(2),vy(2),vz(2),vp,sum_mass,&
       vr(npart),Ekr(npart),mu(npart),vx_cm(npart),vy_cm(npart),&
       vz_cm(npart),np_t,np_red(0:n(1)+2,0:n(2)+2,ntype),vz_sav,v2old,&
       xp,yp,px,py,Vc
  ! Macroscopic parameters
  real(kind=8):: cnt_col(3,ncol_mx,nproc,nm_rg),P_loss(4,ntype,nproc),&
       ss2D(2,0:n(1)+2,0:n(2)+2,ntype,nproc)
  ! Extra
  character:: answer*1

  real(kind=8) :: Ek_before, Ek_after, Ek_diff


  ! Initialization
  vz_sav=0.d0
  vzb_sav=0.d0
  Vc= h(1)*h(2)*zmax

  !
  ! Loop over colliding particles
  !
  Nc_tmp= Nc(ptype,iproc)
  ! Use subroutine for diagnostic only
  if(flag_diag.eq.1) Nc_tmp= np_tot(ptype,iproc)
  do ic=1,Nc_tmp

     ! Get incident particle index
     if(flag_diag.eq.0) then
        rnd(1)= ran2(iseed)
        ip= INT(np_tot(ptype,iproc)*rnd(1)) + 1
        ! Warning
        if( ip.eq.(np_tot(ptype,iproc)+1) ) then
           print*, 'ip= np_tot + 1, please correct...'
           call stop_calculation
        endif
     else
        ip= ic
     endif
     
     ! Find grid location
     xp= vxp(1,ip,ptype,iproc)
     yp= vxp(2,ip,ptype,iproc)
     ix= INT( xp/h(1) ) + 1
     iy= INT( yp/h(2) ) + 1

     if(flag_diag.eq.1) then
        px=( ix*h(1) - xp )/h(1)
        py=( iy*h(2) - yp )/h(2)
                
        ! Bilinear interpolation
        ki(1)= px*py
        ki(2)= (1.d0-px)*py
        ki(3)= (1.d0-px)*(1.d0-py)
        ki(4)= px*(1.d0-py)                          
     endif
     
     ! Get cell index
     ici= (ix-1) + n(1)*(iy-1) + 1
        
     ! Get velocity of selected particle
     vx(1)= vxp(4,ip,ptype,iproc)
     vy(1)= vxp(5,ip,ptype,iproc)
     vz(1)= vxp(6,ip,ptype,iproc)

     ! Initialize variables & arrays
     flag_add=0
     flag_ttype= 0
     sav_it= 0
     Ekr= 0.d0
     vr= 0.d0
     c_ind= 0
     sum_nu= 0.d0
     nu= 0.d0
     sav_tproc= 0
     flag_err= 0
     sort_arr=0.d0
        
     !
     ! Extract sigma for each collision processes
     !
     do icol=1,p_ncol(ptype)

        ! Index of collision and target type
        ind_col= sig_list(ptype,icol)
        ttype= sav_ttype(ind_col)

        ! If no target particle was found in a previous iteration 
        if(flag_err(ttype).eq.1) goto 105 ! jump to next reaction

        ! If target particle has not yet been selected
        if(flag_ttype(ttype).eq.0) then

           ! Monte-Carlo collisions (MC)
           if(charge(ttype).eq.0) goto 107
        
           ! Target particle OMP proc same as incident particle
           tproc= iproc

           ! Target particle cell index same as incident particle
           ict= ici

           ! Look for a target particle
           cnt_ppc= Plist(ict,ttype,tproc) - Plist(ict-1,ttype,tproc)
           if(cnt_ppc.eq.0) then ! no target particles in iproc
              
              ! # of ppc summed over OMP procs (same cell index as incident particle)
              cnt_ppc= SUM(Plist(ict,ttype,1:nproc) - Plist(ict-1,ttype,1:nproc))
              if(cnt_ppc.gt.0) goto 104 ! bypass other tests
              
              ! Nothing in ic, find target particle in neighboring cells (4 in 2D)
              if(cnt_ppc.eq.0) then ! ix + 1
                 ict= ici + 1
                 if( ict.gt.0 .and. ict.le.pl_max ) &
                      cnt_ppc= SUM(Plist(ict,ttype,1:nproc) - Plist(ict-1,ttype,1:nproc))
              else
                 goto 104 ! bypass other tests
              endif
              
              if(cnt_ppc.eq.0) then ! ix - 1
                 ict= ici - 1
                 if( ict.gt.0 .and. ict.le.pl_max ) &
                      cnt_ppc= SUM(Plist(ict,ttype,1:nproc) - Plist(ict-1,ttype,1:nproc))
              else
                 goto 104
              endif

              if(cnt_ppc.eq.0) then ! iy + 1
                 ict= ici + n(1)
                 if( ict.gt.0 .and. ict.le.pl_max ) &
                      cnt_ppc= SUM(Plist(ict,ttype,1:nproc) - Plist(ict-1,ttype,1:nproc))
              else
                 goto 104
              endif

              if(cnt_ppc.eq.0) then ! iy - 1
                 ict= ici - n(1)
                 if( ict.gt.0 .and. ict.le.pl_max ) &
                      cnt_ppc= SUM(Plist(ict,ttype,1:nproc) - Plist(ict-1,ttype,1:nproc))
              else
                 goto 104
              endif
                            
              ! No target particle was found 
              if(cnt_ppc.eq.0) then                 
                 flag_err(ttype)= 1
                 ! Count errors
                 if( charge(ttype).ge.0.d0 .and. ix.le.nx_PE ) &
                      err_coll(ptype,iproc)= err_coll(ptype,iproc) + 1
                 ! Jump to next reaction
                 goto 105 
              endif
              
104           continue
       
              ! If a target particle was indeed found
              cnt_proc= 0
106           cnt_ppc= Plist(ict,ttype,tproc) - Plist(ict-1,ttype,tproc)
              ! Find in which OMP proc is the particle
              if( cnt_ppc.eq.0 .and. cnt_proc.lt.nproc ) then
                 tproc= tproc + 1
                 cnt_proc= cnt_proc + 1
                 if(tproc.gt.nproc) tproc=1 
                 goto 106
              endif
    
           endif

           ! Extract target particle index
           rnd(1)= ran2(iseed)
           it= Plist(ict,ttype,tproc)*(1.d0-rnd(1)) + &
                Plist(ict-1,ttype,tproc)*rnd(1) + 1
           sav_it(ttype)= it
           sav_tproc(ttype)= tproc

           !
           ! Properties of target particle
           !
107        continue

           if(charge(ttype).eq.0.d0) then
              ! Extract velocity of neutral from a Maxwellian
              rnd(1)= ran2(iseed)
              rnd(2)= ran2(iseed)
              call load_gauss(vx(2),vy(2),vt0(ttype),rnd)
              if(vz_sav.eq.0.d0) then
                 rnd(1)= ran2(iseed)
                 rnd(2)= ran2(iseed)
                 call load_gauss(vz(2),vz_sav,vt0(ttype),rnd)
              else
                 vz(2)= vz_sav
                 vz_sav= 0.d0
              endif
           else
              ! Extract velocity of target particle (DSMC)
              vx(2)= vxp(4,it,ttype,tproc)
              vy(2)= vxp(5,it,ttype,tproc)
              vz(2)= vxp(6,it,ttype,tproc)
           endif

           ! Relative velocity, energy in center of mass (CM) frame & CM velocity
           mu(ttype)= ABS(mass(ptype))*ABS(mass(ttype))/ &
                ( ABS(mass(ptype)) + ABS(mass(ttype)) ) ! ptype vs. ttype
              
           vr(ttype)= dsqrt( (vx(1) - vx(2))*(vx(1) - vx(2)) + &
                (vy(1) - vy(2))*(vy(1) - vy(2)) + &
                (vz(1) - vz(2))*(vz(1) - vz(2)) )

           Ekr(ttype)= ( 0.5d0*mu(ttype)*vr(ttype)**2. )/qe ! eV   

           Ek_before = 0.5d0*abs(mass(ptype)*(vx(1)**2.+vy(1)**2.+vz(1)**2.)) &
                  + 0.5d0*abs(mass(ttype)*(vx(2)**2.+vy(2)**2.+vz(2)**2.))
           
           sum_mass= ABS(mass(ptype)) + ABS(mass(ttype)) 
           vx_cm(ttype)= ( ABS(mass(ptype))*vx(1) + ABS(mass(ttype))*vx(2) )/sum_mass
           vy_cm(ttype)= ( ABS(mass(ptype))*vy(1) + ABS(mass(ttype))*vy(2) )/sum_mass
           vz_cm(ttype)= ( ABS(mass(ptype))*vz(1) + ABS(mass(ttype))*vz(2) )/sum_mass

           if( Ekr(ttype).gt.sig_Er(sig_npt_mx) ) then
              print*, 'Warning: Ek_r > max. allowed in cross-section arrays ...'
              print*, 'Ek_r (eV)=',Ekr(ttype)
              print*, 'Do you want to abort calculation? (y/n)'
              read(*,*) answer
              if( answer.eq.'y' .or. answer.eq.'Y' ) call stop_calculation
           endif

           flag_ttype(ttype)=1

        endif ! end select target particle 

        ! Dichotomy
        npt= sig_npt_mx
        
        ipt_L= 1
        ipt_R= npt
        
        do ipt=1,5
           ipt_M= INT(real(ipt_R-ipt_L)/2.d0) + ipt_L
           Ek_L= sig_Er(ipt_L)
           Ek_R= sig_Er(ipt_R)
           Ek_M= sig_Er(ipt_M)
           if( Ekr(ttype).ge.Ek_L .and. Ekr(ttype).le.Ek_M ) ipt_R= ipt_M
           if( Ekr(ttype).gt.Ek_M .and. Ekr(ttype).le.Ek_R ) ipt_L= ipt_M
        enddo
        
        ! Warning
        if(ipt_L.eq.ipt_R) print*, 'Warning: ipt_L= ipt_R in collision_OMP()'
        
        ! Scan cross-section
        ipt_loop: do ipt=ipt_L+1,ipt_R,1
           
           ! Same energy range for all cross-sections
           Ek_L= sig_Er(ipt-1)
           Ek_R= sig_Er(ipt)
    
           ! Linear interpolation
           if( Ekr(ttype).gt.Ek_L .and. Ekr(ttype).le.Ek_R ) then

              sig_L= sig(ipt-1,ind_col)
              sig_R= sig(ipt,ind_col)
                 
              Eth= sig_Eex(ind_col,ind_Eth)
              
              ! Calculate sigma
              sig_p = sig_L + ( Ekr(ttype) - Ek_L )* &
                   ( sig_R - sig_L )/( Ek_R - Ek_L )
                 
              ! Deduce normalized collision frequency (ix_t=ix, iy_t=iy)
              if(charge(ttype).eq.0.d0) then
                 ! MC
                 np_t= np_mx(ttype)
              else
                 ! DSMC
                 np_t= 0.25d0*( np_red(ix,iy,ttype) + np_red(ix+1,iy,ttype) + &
                      np_red(ix,iy+1,ttype) + np_red(ix+1,iy+1,ttype) )                   
              endif

              nu(icol)= np_t*sig_p*vr(ttype)/nu_max_OMP(iproc)
              sum_nu= sum_nu + nu(icol)

              ! Calculate averaged frequencies               
              i_rg_loop: do i_rg=1,n_rg
                 if ( ix.ge.ixl_rg(i_rg) .and. ix.le.ixr_rg(i_rg) .and. &
                      iy.ge.iyl_rg .and. iy.le.iyr_rg ) then
                    cnt_col(1,ind_col,iproc,i_rg)= cnt_col(1,ind_col,iproc,i_rg) + &
                         nu(icol)*nu_max_OMP(iproc)
                    cnt_col(2,ind_col,iproc,i_rg)= cnt_col(2,ind_col,iproc,i_rg) + 1.d0
                    exit i_rg_loop ! exit loop
                 endif
              enddo i_rg_loop

              ! Calculate source term's 2D map 
              if(flag_diag.eq.1) then
                 ! # of reactants & byproducts
                 n_re= col_info(ind_col,1)
                 n_by= col_info(ind_col,2)
                 ! Reaction type (ionization, etc.)
                 rt= col_info(ind_col,2+n_re+n_by+1)

                 ! Store total source and sink terms from ionization, dissociation and charge exchange. 
                 if(rt.ne.1 .and. rt.ne.3) then
                    flag_add(1:2)=-1
                    k2= nu(icol)*nu_max_OMP(iproc)

                    ! Loop over byproduct particles
                    do i_by=1,n_by
                       ! Extract byproduct type                    
                       btype= col_info(ind_col,2+n_re+i_by)
                       if(charge(btype).ne.0) then ! Neutrals not accounted for
                          if(btype.ne.ptype .and. btype.ne.ttype) then
                             ! Density element
                             k1= Nm(btype)/Vc
                             k1= k1*k2

                             ! Source term 
                             ss2D(1,ix,iy,btype,iproc)= ss2D(1,ix,iy,btype,iproc) + ki(1)*k1
                             ss2D(1,ix+1,iy,btype,iproc)= ss2D(1,ix+1,iy,btype,iproc) + ki(2)*k1
                             ss2D(1,ix+1,iy+1,btype,iproc)= ss2D(1,ix+1,iy+1,btype,iproc) + ki(3)*k1
                             ss2D(1,ix,iy+1,btype,iproc)= ss2D(1,ix,iy+1,btype,iproc) + ki(4)*k1
                          endif
                       endif

                       if(btype.eq.ptype) flag_add(1)=0
                       if(btype.eq.ttype) flag_add(2)=0
                                                 
                    enddo

                    ! Loop over reactants
                    do i_re=1,2
                       if(flag_add(i_re).eq.-1) then
                          if(i_re.eq.1) btype=ptype
                          if(i_re.eq.2) btype=ttype
                          
                          if(charge(btype).ne.0) then ! Neutrals not accounted for
                             ! Density element
                             k1= Nm(btype)/Vc
                             k1= k1*k2
                             
                             ! Sink term
                             ss2D(2,ix,iy,btype,iproc)= ss2D(2,ix,iy,btype,iproc) + ki(1)*k1
                             ss2D(2,ix+1,iy,btype,iproc)= ss2D(2,ix+1,iy,btype,iproc) + ki(2)*k1
                             ss2D(2,ix+1,iy+1,btype,iproc)= ss2D(2,ix+1,iy+1,btype,iproc) + ki(3)*k1
                             ss2D(2,ix,iy+1,btype,iproc)= ss2D(2,ix,iy+1,btype,iproc) + ki(4)*k1
                          endif
                       endif
                    enddo
                    
                 endif
              endif
              
              ! Temporary store collision frequencies
              sort_arr(icol)= nu(icol)

              exit ipt_loop ! exit loop

           endif
           
        enddo ipt_loop

105     continue
     enddo
        
     !
     ! Sort collision frequencies
     !
     call indexx(p_ncol(ptype),sort_arr,indx)                          

     !
     ! Check for collision type
     !
     rnd(1)=ran2(iseed)

     flag_coll= 0
     ! Check for null-collisions first
     if( rnd(1).le.sum_nu ) then
        
        sum_nu= 0.d0 ! re-initialize
        do icol=1,p_ncol(ptype)              
           
           ! Sum over collision frequencies starting with smallest
           sum_nu= sum_nu + nu(indx(icol))
           
           ! Check if collision occured
           if( rnd(1).le.sum_nu ) then
              flag_coll= 1
              ! Save reaction index
              c_ind= sig_list(ptype,indx(icol)) 

              exit ! exit loop
           endif
           
        enddo
        
     endif

     if(flag_diag.eq.1) flag_coll= 0
     
     !
     ! Perform collision process
     !
     if(flag_coll.eq.1) then ! flag_coll=0 is null collision

        ! # of reactants & byproducts
        n_re= col_info(c_ind,1)
        n_by= col_info(c_ind,2)
        ! Reaction type (ionization, etc.)
        rt= col_info(c_ind,2+n_re+n_by+1)
        ! Target particle type
        ttype= sav_ttype(c_ind)
        ! Target particle index
        it= sav_it(ttype)
        ! Target particle processor
        tproc= sav_tproc(ttype)

        !
        ! All reactions but charge exchange
        !
        Ek_after = 0.0d0
        if( rt.ne.4 ) then                 

           ! Get threshold energy 
           Eth= sig_Eex(c_ind,ind_Eth)
                            
           ! Remove Eth from impacting particle Ek in CM
           Ee= 0.5d0*mu(ttype)*vr(ttype)**2. - Eth*qe
              
           ! Scattering option
           opt_scat= 1 

           ! Check ionization scheme
           if( rt.eq.2 .and. flag_ionrep.ge.1 .and. n_by.eq.3 ) then
              ! Simplified repartition of energy after an ionization event 
              ! (1) Draw velocity for the ion from a Maxwellian distribution 
              ! in the laboratory frame. (2) Given percentage for the incident 
              ! and byproduct electron. 
              ! This method neither conserve energy nor momentum
              ! We further assume in this over simplified model that the neutral
              ! mass in infinite, i.e., v_cm~ 0 and Ee~ 1/2*me*ve^2 - Eth
              ! Note that the ionization frequency was correctly calculated in the
              ! CM frame.
              cnt_ele=0
              flag_ion=0
              Ek_b= 0.d0
              
              do i_by=1,n_by
                 ! Extract byproduct type                    
                 btype= col_info(c_ind,2+n_re+i_by)
                 
                 ! Electrons
                 if(btype.eq.1) then 
                    cnt_ele= cnt_ele+1
                    if(cnt_ele.eq.1) Ek_b(i_by)= fEk
                    if(cnt_ele.eq.2) Ek_b(i_by)= 1.d0-fEk
                 endif

                 ! Ion
                 if( btype.ne.1 .and. charge(btype).ne.0 ) flag_ion=1                 
              enddo ! end-loop over byproducts 
              
              ! 3-body event with 2 byproduct electrons and 1 ion.
              if( cnt_ele.eq.2 .and. flag_ion.eq.1 ) then 
                 opt_scat= 2
                 ! Neglect contribution from neutral in Ee
                 Ee= 0.5d0*mass(1)*( vx(1)*vx(1) + vy(1)*vy(1) + &
                      vz(1)*vz(1) ) - Eth*qe
                 ! May happen is rare occasions due to the inclusion of the neutral energy above in the selection of collision events.  
                 if(Ee.lt.0.d0) Ee=0.d0
              endif
              
           else ! Scheme which conserves energy and momentum
              
              ! Sum mass inverses
              sum_mass= SUM(1.d0/ABS(mass(col_info(c_ind,2+n_re+1:2+n_re+n_by))))
           
              ! Calculate scattering angles for momentum vectors in frame (Ox',Oy',Oz')
              rnd(1)= ran2(iseed)
              rnd(2)= ran2(iseed)
              
              th_add= 2.d0*pi/real(n_by)           
              cos_th= 1.d0 - 2.d0*rnd(1)           
              th= dacos(cos_th)
              phi= 2.d0*pi*rnd(2)
              cos_ph= dcos(phi)
              sin_ph= dsin(phi)

              ! Calculate scattering angles for unit vector (Ox',Oy',Oz') -> (Ox_CM,Oy_CM,Oz_CM) (CM frame)          
              rnd(1)= ran2(iseed)
              rnd(2)= ran2(iseed)
           
              cos_th_s= 1.d0 - 2.d0*rnd(1)            
              phi_s= 2.d0*pi*rnd(2)
           endif ! end-if ionization scheme

           ! Warning
           if(Ee.lt.0.d0) then 
              print*, 'Warning: negative collision energy was found ...'
              print*, 'Ee (eV)=',Ee/qe,', Eth=',Eth,', c_ind=',c_ind,', flag_ionrep=',flag_ionrep
              print*, 'Abort calculation ...' 
              call stop_calculation
           endif
          
           ! Initialize flag_add
           flag_add(1:2)=-1
           flag_add(n_re+1:n_re+n_by)=0
           
           !
           ! Loop over byproducts
           !
           do i_by=1,n_by

              ! Extract byproduct type                    
              btype= col_info(c_ind,2+n_re+i_by)

              ! Do not add particles with a negative mass
              if(mass(btype).le.0) goto 110 

              !
              ! Add particle to list
              !
              opt_add=1

              ! Set bproc
              bproc= iproc

              ! Update incident particles characteristics
              if( i_by.eq.1 .and. btype.eq.ptype ) then 
                 ib=ip
                 opt_add=0
                 ! Keep track of particle removal
                 flag_add(1)=0
              endif

              ! Elastic collisions & excitation events
              if( rt.eq.1 .or. rt.eq.3 ) then
                 if(btype.eq.ptype) ib=ip
                 if(btype.eq.ttype) then 
                    ib=it
                    bproc= tproc
                 endif
                 opt_add=0
                 ! Keep track of particle removal
                 flag_add(1:2)=0
              endif

              ! Create new particle
              if(opt_add.eq.1) then
                 np_add(btype,bproc)= np_add(btype,bproc) + 1
                 ! New particle array index
                 ib= np_tot(btype,bproc) + np_add(btype,bproc)
                 ! We assume that (ix_t,iy_t)= (ix,iy)
                 if(plt_src.eq.0) ss2D(1,ix,iy,btype,bproc)= ss2D(1,ix,iy,btype,bproc) + 1.
                 ! Same location as incident particle
                 vxp(1,ib,btype,bproc)= vxp(1,ip,ptype,bproc)
                 vxp(2,ib,btype,bproc)= vxp(2,ip,ptype,bproc)
                 vxp(3,ib,btype,bproc)= vxp(3,ip,ptype,bproc)
                 ! Keep track of particle creation
                 flag_add(n_re+i_by)=1
              endif

              !
              ! Calculate byproduct velocity
              !
              if(opt_add.eq.0) then
                 ! Save old velocity
                 v2old= vxp(4,ib,btype,bproc)*vxp(4,ib,btype,bproc) + &
                      vxp(5,ib,btype,bproc)*vxp(5,ib,btype,bproc) + &
                      vxp(6,ib,btype,bproc)*vxp(6,ib,btype,bproc)
              else
                 v2old= 0.d0
              endif

              ! Method which both conserves energy and momentum
              if(opt_scat.eq.1) then
                 ! Theta angle between momentum vectors (equidistant)
                 cos_th= dcos(th)
                 sin_th= dsin(th)
                 th= th + th_add
                 
                 ! Velocity in laboratory frame
                 ex= cos_th
                 ey= sin_th*sin_ph 
                 ez= sin_th*cos_ph              
                 call SCATTER(ex1,ey1,ez1,ex,ey,ez,cos_th_s,phi_s) 
                 
                 vp= dsqrt( 2.d0*(Ee/sum_mass)/ABS(mass(btype))**2. )
                 vbx= vx_cm(ttype) + vp*ex1
                 vby= vy_cm(ttype) + vp*ey1
                 vbz= vz_cm(ttype) + vp*ez1
              endif
              
              ! Method which fixes energy repartition (does not conserve momentum)
              ! me >> mi, we hence assume (1) v_cm~ 0 and (2) Ee~ 1/2*me*ve^2 - Eth
              if(opt_scat.eq.2) then
                 
                 ! Electrons  
                 if(btype.eq.1) then 
                    ! Velocity in laboratory frame
                    rnd(1)= ran2(iseed)
                    rnd(2)= ran2(iseed)           
                    cos_th= 1.d0 - 2.d0*rnd(1)    
                    sin_th= 2.d0*dsqrt(rnd(1)*(1.d0-rnd(1)))       
                    phi= 2.d0*pi*rnd(2)
                    ex1= cos_th
                    ey1= sin_th*dsin(phi) 
                    ez1= sin_th*dcos(phi)              
                    
                    ! Ek_b is a given percentage of Ee for electrons
                    vp= dsqrt( 2.d0*(Ee*Ek_b(i_by)/ABS(mass(btype))) )
                    vbx= vp*ex1 
                    vby= vp*ey1
                    vbz= vp*ez1                   
                 else
                    ! Ion            
                    rnd(1)= ran2(iseed)
                    rnd(2)= ran2(iseed)
                    
                    ! Draw the ion kinetic energy from a Maxwellian 
                    call load_gauss(vbx,vby,vt0(btype),rnd)
                    if(vzb_sav.eq.0.d0) then
                       rnd(1)= ran2(iseed)
                       rnd(2)= ran2(iseed)
                       call load_gauss(vbz,vzb_sav,vt0(btype),rnd)
                    else
                       vbz= vzb_sav
                       vzb_sav= 0.d0
                    endif
                 endif
              endif

              ! Update velocity
              vxp(4,ib,btype,bproc)= vbx
              vxp(5,ib,btype,bproc)= vby
              vxp(6,ib,btype,bproc)= vbz

              Ek_after = Ek_after + 0.5d0*abs(mass(btype)*(vbx**2.+vby**2.+vbz**2.))

              ! Total power lost/gained per time step
              P_loss(3,btype,bproc)= P_loss(3,btype,bproc) + 0.5d0*Nm(btype)*mass(btype)*( &
                   vxp(4,ib,btype,bproc)*vxp(4,ib,btype,bproc) + &
                   vxp(5,ib,btype,bproc)*vxp(5,ib,btype,bproc) + &
                   vxp(6,ib,btype,bproc)*vxp(6,ib,btype,bproc) - v2old ) 
              
110           continue
           enddo ! end-Loop over byproducts

         !   if (rt.eq.1 .and. bproc.eq.1) then
         !       Ek_diff= Ek_before - Ek_after
         !       print*, ' '
         !       print*, '------------   Collision diagnostic ------------'
         !       print*, pname(ptype), " + ", pname(ttype), " -> ", &
         !            (pname(col_info(c_ind,2+n_re+i_by)), i_by=1,n_by)
         !       print*, 'Collision number: ', c_ind
         !       print*, 'Threshold energy = ', Eth, "eV"
               
         !       print*, 'Energy before collision = ', Ek_before/qe, "eV"
         !       print*, 'Energy after collision = ', Ek_after/qe, "eV"
         !       if(Ek_diff.gt.0.d0) then
         !          print*, 'Energy loss = ', Ek_diff/qe, "eV"
         !       else
         !          print*, 'Energy gain = ', abs(Ek_diff)/qe, "eV"
         !       endif
         !       print*, '------------------------------------------------'
         !       print*, ' '
         !       call stop_calculation
         !   endif
        
           !
           ! Reactant particles
           !
           do i_re=1,n_re

              ! Reactant type
              ctype= (col_info(c_ind,2+i_re))

              ! Do not kill particles with a negative mass
              if(mass(ctype).le.0) then 
                 flag_add(i_re)=0
                 goto 120
              endif

              if( flag_add(i_re).eq.-1 )  then                       
                 ! Kill target particle (done in mover)
                 if(ctype.eq.ptype) then 
                    ib=ip
                    bproc=iproc
                 endif
                 if(ctype.eq.ttype) then 
                    ib=it
                    bproc=tproc
                 endif
                 flag_dead(ib,ctype,bproc)=1
                 ! Sink term
                 if(plt_src.eq.0) ss2D(2,ix,iy,ctype,bproc)= ss2D(2,ix,iy,ctype,bproc) + 1.
              endif

120           continue
           enddo
              
        endif ! end-if all reaction types except charge exchange events

        !
        ! Charge exchange collisions
        !
        if( rt.eq.4 ) then
              
           ! Initialization
           Eth= 0.d0
           flag_add(1:n_re+n_by)= 0

           do i_by=1,n_by
              
              ! Extract byproduct type                    
              btype= col_info(c_ind,2+n_re+i_by)
              
              ! Do not add particles with a negative mass
              if(mass(btype).le.0) goto 130
              
              if(btype.eq.ptype) then
                 ib= ip
                 ind_ce= 1
                 bproc=iproc
              else
                 ib= sav_it(ttype)
                 ind_ce= 2
                 bproc= sav_tproc(tproc)
              endif
 
              ind_ce= ind_ce + 1
              if(ind_ce.eq.3) ind_ce= 1
              
              ! Save old velocity
              v2old=  vxp(4,ib,btype,bproc)*vxp(4,ib,btype,bproc) + &
                   vxp(5,ib,btype,bproc)*vxp(5,ib,btype,bproc) + &
                   vxp(6,ib,btype,bproc)*vxp(6,ib,btype,bproc)

              ! Update velocity of target particle
              vxp(4,ib,btype,bproc)= vx(ind_ce)
              vxp(5,ib,btype,bproc)= vy(ind_ce)
              vxp(6,ib,btype,bproc)= vz(ind_ce)

              ! Total power lost/gained per time step
              P_loss(3,btype,bproc)= P_loss(3,btype,bproc) + 0.5d0*Nm(btype)*mass(btype)*( &
                   vxp(4,ib,btype,bproc)*vxp(4,ib,btype,bproc) + &
                   vxp(5,ib,btype,bproc)*vxp(5,ib,btype,bproc) + &
                   vxp(6,ib,btype,bproc)*vxp(6,ib,btype,bproc) - v2old ) 
               
               Ek_after = Ek_after + 0.5d0*abs(mass(btype)*(vx(ind_ce)**2.+vy(ind_ce)**2.+vz(ind_ce)**2.))
                            
130           continue
           enddo
              
        endif ! End-if charge exchange collisions
        
        if (rt.eq.1 .and. bproc.eq.1) then
               Ek_diff= Ek_before - Ek_after
               print*, ' '
               print*, '------------   Collision diagnostic ------------'
               print*, pname(ptype), " + ", pname(ttype), " -> ", &
                    (pname(col_info(c_ind,2+n_re+i_by)), i_by=1,n_by)
               print*, 'Collision number: ', c_ind
               print*, 'Threshold energy = ', Eth, "eV"
               
               print*, 'Energy before collision = ', Ek_before/qe, "eV"
               print*, 'Energy after collision = ', Ek_after/qe, "eV"
               if(Ek_diff.gt.0.d0) then
                  print*, 'Energy loss = ', Ek_diff/qe, "eV"
               else
                  print*, 'Energy gain = ', abs(Ek_diff)/qe, "eV"
               endif
               print*, '------------------------------------------------'
               print*, ' '
               call stop_calculation
         endif

        ! Count total number of actual collision events
        i_rg_loop2: do i_rg=1,n_rg
           if ( ix.ge.ixl_rg(i_rg) .and. ix.le.ixr_rg(i_rg) .and. &
                iy.ge.iyl_rg .and. iy.le.iyr_rg ) then
              cnt_col(3,c_ind,iproc,i_rg)= cnt_col(3,c_ind,iproc,i_rg) + &
                   1.d0/np_tot_rg(ptype,i_rg)
              exit i_rg_loop2 ! exit loop
           endif
        enddo i_rg_loop2

     endif ! end-if a collision occured
     
  enddo ! end-loop over Nc collisions
  
  return
end subroutine collision_OMP

SUBROUTINE SCATTER(vx1,vy1,vz1,vx,vy,vz,costheta,phi)
!     ===================================================================
!     VERSION:         0.1
!     LAST MOD:      April/19
!     MOD AUTHOR:    G. Hagelaar
!     COMMENTS:     
!     NOTES:   Turn vector vx,vy,vz over a scattering angle theta and azimuthal angle phi
!     -------------------------------------------------------------------
  implicit none
  real(kind=8):: costheta,sintheta,cosphi,sinphi,phi,vx,vy,vz,v,vv,vx1,vy1,vz1
  
  sintheta=dsqrt(dmax1(1d0-costheta**2,0d0))
  sinphi=dsin(phi)
  cosphi=dcos(phi)
  v=dsqrt(vx**2+vy**2+vz**2)
  IF (v==0) RETURN
  IF (dabs(vy)>dabs(vz)) THEN
     vv=dsqrt(vx**2+vy**2)
     vx1=vx*costheta+(vy*v*sinphi+vx*vz*cosphi)/vv*sintheta
     vy1=vy*costheta+(-vx*v*sinphi+vy*vz*cosphi)/vv*sintheta
     vz1=vz*costheta-vv*cosphi*sintheta
  ELSE
     vv=dsqrt(vx**2+vz**2)
     vx1=vx*costheta+(vz*v*sinphi-vy*vx*cosphi)/vv*sintheta
     vy1=vy*costheta+vv*cosphi*sintheta
     vz1=vz*costheta-(vx*v*sinphi+vy*vz*cosphi)/vv*sintheta
  ENDIF
return
END SUBROUTINE SCATTER
