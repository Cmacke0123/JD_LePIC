subroutine part_injection(n,h,bcnd,vxp,ss2D,ntype,nmax,ni0,&
     I_inj,np_tot,nproc,iseed,nproc_mpi,N_inj,iproc,ptype,&
     P_loss,phi)
!     ==============================================================
!     VERSION:         0.5
!     LAST MOD:      Mar/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:    
!     NOTE:          
!     --------------------------------------------------------------
  implicit none
  include 'mpif.h'
  include 'particle_info.h'
  include 'constants.h'
  integer nproc_mpi
  integer:: ix,iy,i,k
  integer:: ntype,ptype,nmax,n(3),iproc,nproc
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2),np_tot(ntype,nproc),&
       N_inj(ntype,nproc),N_inj_tmp,flag
  real(kind=8):: h(3),vxp(6,nmax,ntype,nproc),x,y,z,vx,vy,vz,vt,vb(ntype),rnd(2),&
       ran2,ni0(npart),I_inj,dN_inj,vz_sav,x_sav,P_loss(4,ntype,nproc),&
       ss2D(2,0:n(1)+2,0:n(2)+2,ntype,nproc),vb_tmp
  real(kind=8):: phi(0:n(1)+2,0:n(2)+2),phip,ki(4),px,py,dt_tmp,&
       vmax(ntype),fmax(ntype)
  ! Macroscopic parameters 
  integer:: iseed

  ! Initialize particle counter, variables & arrays
  vz_sav=0.d0
  x_sav= 0.d0
  flag= 0

  ! Parameters for the shifted Maxwellian flux MC subroutine
  ! opt=3 : inject electron beam along (OZ)
  ! opt=4 : electron/ion emission off a cathode along +X
  vb(ptype)=0.d0
  vmax(ptype)=0.d0
  fmax(ptype)=0.d0
  if(ABS(opt_inj).eq.3 .or. ABS(opt_inj).eq.4) then
     vb(ptype)= dsqrt(2.d0*qe*THm/ABS(mass(ptype)))
     if(ABS(opt_inj).eq.4) then
        ! Option specific for divertor modeling
        if(ptype.eq.1 .and. flag_nmn.eq.1) vb(ptype)= 0.d0  ! Half-Maxwellian flux
        call init_shifted_maxwellian_flux(vb,vmax,fmax,ptype,ntype)
     endif
  endif

  ! # of particles per specie per time step to inject 
  if(ABS(opt_inj).eq.2) then
     ! Re-inject an electron-ion pair for each positive ion lost
     if( ptype.eq.1 .or. (ptype.eq.tag_neg .and. tag_neg.gt.0) ) then 
        flag= 1
        ! # of negative charges to re-inject
        N_inj_tmp= SUM(N_inj(2:ntype,iproc))
        if(tag_neg.gt.0) N_inj_tmp= N_inj_tmp - &
             N_inj(tag_neg,iproc)

        ! # percentage ratio for selected negatively charged specie
        dN_inj=ni0(ptype)*N_inj_tmp  
     endif
  else
     ! Inject a fixed particle current 
     flag= 1
     ! Number of particles to inject per time step per OMP and MPI threads
     dN_inj= I_inj*ns_inj*dt*ni0(ptype)/(qe*Nm(1))/real(nproc_mpi)/real(nproc)
  endif
  
  if(flag.eq.1) then
     N_inj_tmp= INT(dN_inj)
     ! Correct for round-off errors
     rnd(1)= ran2(iseed)
     if( rnd(1).le.(dN_inj-N_inj_tmp) ) N_inj_tmp= N_inj_tmp + 1
     N_inj(ptype,iproc)= N_inj_tmp
  endif

  ! Inject electron beam
  if( ABS(opt_inj).eq.3 .or. &
       (ABS(opt_inj).eq.4 .and. flag_nmn.eq.0) ) N_inj(2:ntype,iproc)=0 

  ! Inject particles
  do i=1,N_inj(ptype,iproc)
                 
     ! Random location
70   if(opt_inj.gt.0) then 
60      rnd(1)=ran2(iseed)
        x= rnd(1)*(xl_pow-xr_pow) + xr_pow ! inside [x<,x>]
        if(flag_c.eq.1) then ! Disk
           rnd(1)=ran2(iseed)
           y= rnd(1)*(2.d0*dr) + ymax/2.d0-dr
           if( ((x-xa)**2 + (y-ymax/2.d0)**2).gt.dr**2 ) goto 60
        endif
     endif

     ! Gaussian loading
     if(opt_inj.lt.0) then
        if(flag_c.eq.0) then ! Slit
50         if(x_sav.eq.0) then
              rnd(1)= ran2(iseed)
              rnd(2)= ran2(iseed)
              ! Only along (Ox)
              call load_gauss(x,x_sav,dr,rnd)
              x= x + xa 
              x_sav= x_sav + xa
           else
              x = x_sav
              x_sav= 0.d0
           endif
           if( flag_nmn.eq.1 .and. xa.eq.0.d0) x=ABS(x)
           if( x.lt.(xa-4.d0*dr) .or. x.gt.(xa+4.d0*dr) .or. &
                x.lt.0.d0 .or. x.gt.xmax ) goto 50
        else ! Disk           
40         rnd(1)= ran2(iseed)
           rnd(2)= ran2(iseed)
           call load_gauss(x,y,dr,rnd)
           x= x + xa 
           y= y + ymax/2.d0           
           if( ((x-xa)**2 + (y-ymax/2.d0)**2).gt.4.d0*dr**2 .or.  &
                x.lt.0.d0 .or. x.gt.xmax .or. &
                y.lt.0.d0 .or. y.gt.ymax ) goto 40           
        endif
     endif

     if(flag_c.eq.0) then ! Slit
        rnd(1)=ran2(iseed)
        y= rnd(1)*(yl_pow-yr_pow) + yr_pow ! inside [y<,y>]
     endif
     rnd(1)=ran2(iseed)
     z= rnd(1)*zmax
     
     ! Get particle left grid index
     ix= INT( x/h(1) ) + 1
     iy= INT( y/h(2) ) + 1
     
     ! Load uniquely inside simulation domain
     if( bcnd(ix,iy).ge.1 .and. bcnd(ix+1,iy).ge.1 .and. &
          bcnd(ix+1,iy+1).ge.1 .and. bcnd(ix,iy+1).ge.1 ) goto 70
     
     ! Add particle to counter
     np_tot(ptype,iproc)= np_tot(ptype,iproc) + 1
          
     ! Particle index
     k= np_tot(ptype,iproc)
        
     ! Warning
     if(k.gt.nmax) then
        print*, 'k > nmax in part_injection'
        print*, 'Abort calculation ...'
        call stop_calculation
     endif
        
     ! Thermal velocity
     vt= vt0(ptype)
     ! Option to use Tb as Te for the plasma electron source term
     if(ABS(opt_inj).le.2) then
        if( flag_Tp.eq.1 .and. ptype.eq.1 ) &
             vt= dsqrt(2.d0*qe*THm/ABS(mass(ptype)))
     endif

     ! Shifted Maxwellian flux distribution with u=vB and T=Te
     if(ABS(opt_inj).eq.4) then                
        ! Emission from the LHS boundary along (OX)
        vb_tmp= vb(ptype)*dcos(th_B)
        ! Shifted Maxwellian flux distribution along (OX)
75      call shifted_maxwellian_flux(vx,vb_tmp,vt,fmax(ptype),iseed)     
        rnd(1)= ran2(iseed)
        rnd(2)= ran2(iseed)
        call load_gauss(vy,vz,vt,rnd)
        vy= vy + vb(ptype)*dsin(th_B) ! Shifted Maxwellian along (OY)
        
        ! Spread the position over the distance traveled during one time step ns_inj*dt
        rnd(1)= ran2(iseed)
        dt_tmp= rnd(1)*(real(ns_inj)*dt)
        x= vx*dt_tmp ! Inject on the LHS at x=0
        
        ! Check bounds
        if(x.lt.0.d0 .or. x.gt.xmax) goto 75
        ix= FLOOR( x/h(1) ) + 1           
     else
        ! Gaussian distribution
        rnd(1)= ran2(iseed)
        rnd(2)= ran2(iseed)
        call load_gauss(vx,vy,vt,rnd)
        if(vz_sav.eq.0.d0) then
           rnd(1)= ran2(iseed)
           rnd(2)= ran2(iseed)
           call load_gauss(vx,vy,vt,rnd)
           vz= vx 
           vz_sav= vy
        else
           vz= vz_sav 
           vz_sav= 0.d0
        endif
     endif
     
     ! Save particle 6d-coordinates
     vxp(1,k,ptype,iproc)= x
     vxp(2,k,ptype,iproc)= y
     vxp(3,k,ptype,iproc)= z        
     vxp(4,k,ptype,iproc)= vx
     vxp(5,k,ptype,iproc)= vy
     if(ABS(opt_inj).eq.3) vz= vz + vb(ptype) ! Shifted Maxwellian along (OZ)
     vxp(6,k,ptype,iproc)= vz
     
     if(plt_src.eq.0) ss2D(1,ix,iy,ptype,iproc)= ss2D(1,ix,iy,ptype,iproc) + 1.
     
     P_loss(4,ptype,iproc)= P_loss(4,ptype,iproc) + 0.5d0*Nm(ptype)*mass(ptype)*( &
          vxp(4,k,ptype,iproc)*vxp(4,k,ptype,iproc) + &
          vxp(5,k,ptype,iproc)*vxp(5,k,ptype,iproc) + &
          vxp(6,k,ptype,iproc)*vxp(6,k,ptype,iproc) ) 

     ! if div j_tot is not null one must add the term
     ! int_V [e*phi*S_tot*dV] which comes from the volume integration
     ! of the plasma energy equation; S_i= n_i*nu; n_i= sum_k delta(x-x_i)
     ! See Fundamentals of Plasma Physics, Golant et al. p. 173
     if(flag_diffsrc.eq.1) then
        px=( ix*h(1) - x )/h(1)
        py=( iy*h(2) - y )/h(2)

        ki(1)= px*py
        ki(2)= (1.d0-px)*py
        ki(3)= (1.d0-px)*(1.d0-py)
        ki(4)= px*(1.d0-py)

        phip= ki(1)*phi(ix,iy) + &
             ki(2)*phi(ix+1,iy) + &
             ki(3)*phi(ix+1,iy+1) + &
             ki(4)*phi(ix,iy+1)

        P_loss(1,ptype,iproc)= P_loss(1,ptype,iproc) - Nm(ptype)*charge(ptype)*phip
     endif
     
  enddo

  return

end subroutine part_injection

subroutine shifted_maxwellian_flux(v,vb,vt,fmax,iseed)  
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Apr/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:    This subroutine is designed to generate a shifted
!                  Maxwellian flux distribution
!     NOTE:          
!     --------------------------------------------------------------
  use mpi
  implicit none
  integer iseed
  real(kind=8):: v,vb,vt,fmax,f,ran2,vm,vp,rnd(2)

  vm= MAX(0.d0,vb-4.d0*vt) ! vmin
  vp= vb+4.d0*vt ! vmax
  
90 rnd(1)= ran2(iseed)
  rnd(2)= ran2(iseed)
  v= vm + rnd(1)*(vp-vm)
  f= v*dexp(-(v-vb)**2/vt**2)
  
  ! acceptance/rejection
  if(rnd(2).gt.f/fmax) goto 90

  return
end subroutine shifted_maxwellian_flux

subroutine init_shifted_maxwellian_flux(vb,vmax,fmax,ptype,ntype)  
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Apr/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:    This subroutine is designed to generate the
!                  parameters for the shifted Maxwellian flux MC
!                  calculation.
!
!     NOTE:          
!     --------------------------------------------------------------
  use mpi
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer ptype,ntype
  real(kind=8):: vb(ntype),vmax(ntype),fmax(ntype)

  ! Parameters for injection/rejection MC calculation.
  vmax(ptype)= (vb(ptype) + dsqrt(vb(ptype)**2 + 2.d0*vt0(ptype)**2))/2.d0
  fmax(ptype)= vmax(ptype)*dexp(-(vmax(ptype)-vb(ptype))**2/vt0(ptype)**2)
  
  return
end subroutine init_shifted_maxwellian_flux
