subroutine part_mover1D(istep,n,h,Ei,Bi,p_mac,P_loss,vxp,&
     bcnd,nmax,ntype,ngrid,flag_dead,nproc,np_tot,&
     iproc,ptype,iseed,cnt_dead,sum_q_y,n_B,h_B)
!     ==============================================================
!     VERSION:         0.4
!     LAST MOD:      April/24
!     MOD AUTHOR
!     NOTE:  1D version of the particle mover
!                Interpolation |------|-------------|
!                              ip      xp           ip+1 
!                              <------><------------>
!                                 1-p         p  
!
!                 Indexes in vxp(): 1== x
!                                   2== y    
!                                   3== z
!                                   4== vx
!                                   5== vy
!                                   6== vz                
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: n(3),n_B(2),nmax,ntype,ngrid,igrid
  integer:: ix,iy,i,ptype,istep,flag_BS,flag_lost,iproc,&
       nproc,ptype_sec,ip_sec,n_sec
  ! Field arrays
  real(kind=8):: h(3),h_B(2),Ei(2,0:n(1)+2,0:n(2)+2),Exp,ki(2), &
       Bi(4,0:n_B(1)+2,0:n_B(2)+2)
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2)
  integer(kind=1):: flag_dead(nmax,ntype,nproc)
  real(kind=8):: vxp(6,nmax,ntype,nproc),xp_new,vpx_new,&
       vpy_new,vpz_new,Eki,px,sum_q_y(ngrid,0:n(2)+2,ntype,nproc),&
       vx_sec,vy_sec,vz_sec,vt
  real(kind=8):: v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,Bpx,Bpy,Bpz,B_tot,&
       k1,k2,vb(ntype),vmax(ntype),fmax(ntype),vbx,dt_tmp
  ! Macroscopic parameters
  integer:: np_tot(ntype,nproc),np_lost(ntype,nproc),i_shift,iseed
  real(kind=8):: p_mac(ntype,2,0:ngrid,nproc),P_loss(4,ntype,nproc),&
       ran2,rnd(2),cnt_dead(nproc)

  ! Initialize particle arrays
  np_lost(ptype,iproc)=0

  ! Parameters for refluxing on the Neumann boundary condition
  vb(ptype)=0.d0
  vmax(ptype)=0.d0
  fmax(ptype)=0.d0
  if(ABS(opt_inj).eq.4 .and. flag_nmn.eq.1) then
     vb(ptype)= dsqrt(2.d0*qe*THm/ABS(mass(ptype)))
     ! Option specific for divertor modeling
     if(ptype.eq.1 .and. flag_nmn.eq.1) vb(ptype)= 0.d0  ! Half-Maxwellian flux
     call init_shifted_maxwellian_flux(vb,vmax,fmax,ptype,ntype)
  endif

  ! Set coefficients for Boris solver
  k1= dt*charge(ptype)/(2.d0*mass(ptype))

  !
  ! Move particles
  !
  do i=1,np_tot(ptype,iproc),1

     ! Initialize
     Eki=0.d0

     !
     ! Save actual velocity & position (avoid calling arrays many times)
     !
     xp_new= vxp(1,i,ptype,iproc)
     vpx_new= vxp(4,i,ptype,iproc)
     vpy_new= vxp(5,i,ptype,iproc)
     vpz_new= vxp(6,i,ptype,iproc)

     !
     ! Remove particles killed during a collision in the previous step 
     !
     if(flag_dead(i,ptype,iproc).eq.1) then
        ! Add another lost particle
        np_lost(ptype,iproc)= np_lost(ptype,iproc) + 1
        ! Re-initialize 
        flag_dead(i,ptype,iproc)= 0       
        ! Counter 
        if(ptype.eq.tag_neg) cnt_dead(iproc)= cnt_dead(iproc) + Nm(ptype)*charge(ptype)
        ! Calculate kinetic energy of macroparticle
        Eki= 0.5d0*Nm(ptype)*mass(ptype)*( vpx_new*vpx_new + &
             vpy_new*vpy_new + vpz_new*vpz_new ) 
        ! Total power lost per time step
        P_loss(3,ptype,iproc)= P_loss(3,ptype,iproc) - Eki
        ! Jump to the end of the loop
        goto 110
     endif

     !
     ! Get particle left grid index
     !
     ix= INT( xp_new/h(1) ) + 1
     iy= 1

     !
     ! Calculate fields at macroparticle location
     !
     px=( ix*h(1) - xp_new )/h(1)

     ki(1)= px
     ki(2)= (1.d0-px)

     Exp= ki(1)*Ei(1,ix,iy) + &
          ki(2)*Ei(1,ix+1,iy)
          
     !
     ! Calculate new velocity
     !
     flag_BS=0
     if(flag_B.eq.1) then
        
        if(n_B(1).eq.1) then ! Constant B-field

           Bpx= Bi(1,1,1)
           Bpy= Bi(2,1,1)
           Bpz= Bi(3,1,1)
           
           flag_BS=1
        else ! B-field map
           
           if(flag_gridB.eq.1) then
              ! Get particle left grid index on B-field map
              ix= INT( xp_new/h_B(1) ) + 1
           
              ! Calculate B-fields at position x and y on the map
              px=( ix*h_B(1) - xp_new )/h_B(1)

              ki(1)= px
              ki(2)= (1.d0-px)
           endif
        
           ! Bx
           Bpx= ki(1)*Bi(1,ix,iy) + &
                ki(2)*Bi(1,ix+1,iy)
           
           ! By
           Bpy= ki(1)*Bi(2,ix,iy) + &
                ki(2)*Bi(2,ix+1,iy)
           
           ! Bz
           Bpz= ki(1)*Bi(3,ix,iy) + &
                ki(2)*Bi(3,ix+1,iy)
         
           flag_BS=1           
        endif
     endif

     ! Positive ions not magnetized
     if( flag_B_pos.eq.1 .and. charge(ptype).gt.0.d0 ) flag_BS=0

     if(flag_BS.eq.1) then 
        ! Boris scheme
        B_tot= dsqrt( Bpx*Bpx + Bpy*Bpy + Bpz*Bpz )  
        k2= k1*( 2.d0/(1.d0 + (k1*B_tot)*(k1*B_tot)) )

        v1x= vpx_new + k1*Exp
        v1y= vpy_new
        v1z= vpz_new

        v3x= v1x + k1*( v1y*Bpz - v1z*Bpy )
        v3y= v1y + k1*( v1z*Bpx - v1x*Bpz )
        v3z= v1z + k1*( v1x*Bpy - v1y*Bpx )

        v2x= v1x + k2*( v3y*Bpz - v3z*Bpy )
        v2y= v1y + k2*( v3z*Bpx - v3x*Bpz )
        v2z= v1z + k2*( v3x*Bpy - v3y*Bpx ) 

        vpx_new= v2x + k1*Exp
        vpy_new= v2y
        vpz_new= v2z
     else 
        ! Electrostatic Leap-Frog solver
        vpx_new= vpx_new + 2.d0*k1*Exp           
     endif

     ! Calculate new position  
     xp_new= xp_new + dt*vpx_new

     ! Get particle new grid location
     ix= FLOOR( xp_new/h(1) ) + 1

     ! Correct for particles out of bounds
     call check_outofbounds(ix,iy,n)

     !
     ! Check if the particle did not leave the simulation domain
     !
     flag_lost= 0
     if( (xp_new.le.0.d0 .and. flag_nmn.eq.0) .or. xp_new.ge.xmax ) then 
        flag_lost= 1
     endif

     !
     ! Lost particles
     !
     if(flag_lost.ge.1) then
        ! Add another lost particle
        np_lost(ptype,iproc)= np_lost(ptype,iproc) + 1
        
        ! Calculate kinetic energy of macroparticle
        Eki= 0.5d0*Nm(ptype)*mass(ptype)*( vpx_new*vpx_new + &
             vpy_new*vpy_new + vpz_new*vpz_new ) 
        
        ! Total lost power per time step
        P_loss(1,ptype,iproc)= P_loss(1,ptype,iproc) + Eki
        
        ! Get grid index
        igrid= bcnd(ix,iy)

        ! Dielectrics        
        if(flag_dielec.eq.1) then
           if( dtype(igrid).eq.1 ) &
                sum_q_y(igrid,iy,ptype,iproc)= sum_q_y(igrid,iy,ptype,iproc) + &
                Nm(ptype)*charge(ptype) 
        endif

        ! Secondary particle emission
        if(charge(ptype).gt.0 .and. flag_sec.ge.1) then
           if(igrid.eq.igrid_sec) then
              rnd(1)= ran2(iseed)
              ptype_sec = 1 ! Electrons
              n_sec= INT(gam_sec)
              if(rnd(1) .le. (gam_sec-n_sec)) n_sec= n_sec+1
              if(n_sec.eq.0) goto 130
              ! Generate secondary particle
              do ip_sec=1,n_sec

                 ! Create a new electron
                 np_tot(ptype_sec,iproc)= np_tot(ptype_sec,iproc) + 1
                 i_shift= np_tot(ptype_sec,iproc)

                 ! Use THm for electron temperature
                 vt= dsqrt(2.d0*qe*THm/ABS(mass(ptype_sec))) 
                 rnd(1)=ran2(iseed)
                 vx_sec = -sign(1.d0,vpx_new)*vt*dsqrt( -dlog(1-rnd(1)) )
                 rnd(1)= ran2(iseed)
                 rnd(2)= ran2(iseed)
                 ! Gaussian loading (flux normal to the grid surface)
                 call load_gauss(vy_sec,vz_sec,vt,rnd)
                 
                 ! y-location of impacting ion.
                 vxp(1,i_shift,ptype_sec,iproc)= xg_sec
                 vxp(4,i_shift,ptype_sec,iproc)= vx_sec 
                 vxp(5,i_shift,ptype_sec,iproc)= vy_sec 
                 vxp(6,i_shift,ptype_sec,iproc)= vz_sec
                
                 ! Count power injected into the plasma electrons
                 P_loss(4,ptype_sec,iproc)= P_loss(4,ptype_sec,iproc) + 0.5d0*Nm(ptype_sec)*mass(ptype_sec)*( &
                      vx_sec*vx_sec + vz_sec*vz_sec + vz_sec*vz_sec )
                 
              enddo
130           continue
           endif
        endif

        ! Total number of particle lost at the wall per time step
        p_mac(ptype,np_loss,igrid,iproc)= p_mac(ptype,np_loss,igrid,iproc) + 1
        
        ! Total lost power at the wall per time step
        p_mac(ptype,P_w,igrid,iproc)= p_mac(ptype,P_w,igrid,iproc) + Eki
        
     endif ! end-if flag_lost=1
     
     !
     ! Particles inside the simulation domain
     !
     if(flag_lost.eq.0) then

        ! Neumann BCs (LHS only)
        if(flag_nmn.eq.1) then
           if( xp_new.le.0.d0 ) then
              if(ABS(opt_inj).ne.4) then
                 ! Specular reflection
                 xp_new= -xp_new
                 vpx_new= -vpx_new
              else                 
                 ! Refluxing
                 vbx= vb(ptype)*dcos(th_B)
                 vt= vt0(ptype)

                 ! Calculate kinetic energy of macroparticle
                 Eki= 0.5d0*Nm(ptype)*mass(ptype)*( vpx_new*vpx_new + &
                      vpy_new*vpy_new + vpz_new*vpz_new ) 
                 ! Total power lost
                 P_loss(1,ptype,iproc)= P_loss(1,ptype,iproc) + Eki
                 
                 ! Shifted Maxwellian flux distribution
75               call shifted_maxwellian_flux(vpx_new,vbx,vt,fmax(ptype),iseed)     
                 rnd(1)= ran2(iseed)
                 rnd(2)= ran2(iseed)
                 call load_gauss(vpy_new,vpz_new,vt,rnd)
                 vpy_new= vpy_new + vb(ptype)*dsin(th_B) ! Shifted Maxwellian along (OY)
                 
                 ! Spread the position over the distance traveled during one time step ns_inj*dt
                 rnd(1)= ran2(iseed)
                 dt_tmp= rnd(1)*dt
                 xp_new= vpx_new*dt_tmp ! Inject on the LHS at x=0
                 
                 ! Sanity check
                 if(xp_new.lt.0.d0 .or. xp_new.gt.xmax) goto 75
                 
                 ! Calculate kinetic energy of macroparticle
                 Eki= 0.5d0*Nm(ptype)*mass(ptype)*( vpx_new*vpx_new + &
                      vpy_new*vpy_new + vpz_new*vpz_new ) 
                 ! Total power injected
                 P_loss(4,ptype,iproc)= P_loss(4,ptype,iproc) + Eki                
              endif
           endif
        endif
        
        ! Shift location of particle inside array
        i_shift= i-np_lost(ptype,iproc)
        vxp(1,i_shift,ptype,iproc)= xp_new
        vxp(4,i_shift,ptype,iproc)= vpx_new
        vxp(5,i_shift,ptype,iproc)= vpy_new
        vxp(6,i_shift,ptype,iproc)= vpz_new

     endif

110  continue
  enddo ! end loop over np_tot(ptype) particles

  !
  ! Update particle counter
  !
  np_tot(ptype,iproc)= np_tot(ptype,iproc) - np_lost(ptype,iproc)

  return
 
end subroutine part_mover1D

subroutine calc_Efield1D(n,h,phi,Ei)
!     ===================================================================
!     VERSION:         0.1
!     LAST MOD:      Sep/21
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     NOTE:    
!     -------------------------------------------------------------------
  implicit none
  integer:: ix,n(3)
  real(kind=8):: h(3),phi(0:n(1)+2,0:n(2)+2),Ei(2,0:n(1)+2,0:n(2)+2)
  include 'particle_info.h'

  !
  ! Interior points
  !
  !$OMP PARALLEL
  !$OMP DO
  do ix=1,n(1)+1
     ! second order correct in dx
     Ei(1,ix,1)= -( phi(ix+1,1)-phi(ix-1,1) )/(2.d0*h(1))
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
  
  !
  ! Boundary conditions
  !
  
  ! Walls  
  Ei(1,1,1)= 2.d0*Ei(1,2,1) - Ei(1,3,1)  ! Ex, LHS
  Ei(1,n(1)+1,1)= 2.d0*Ei(1,n(1),1) - Ei(1,n(1)-1,1) ! Ex, RHS  

  ! Neumann BCs
  if( flag_nmn.eq.1 ) Ei(1,1,1)= 0.d0 ! Ex= 0
      
  return
  
end subroutine calc_Efield1D


subroutine charge_deposition1D(n,h,vxp,nmax,ntype,kq,np,nproc,np_tot,sum_dEk,&
     Nh,iproc,ptype,n_icp)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Apr/24
!     MOD AUTHOR     G. Fubiani
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: n(3),nmax,ntype,Nh(n_icp,nproc)
  integer:: ix,iy,i,ptype,iproc,nproc,n_icp
  real(kind=8):: h(3),ki(2),vxp(6,nmax,ntype,nproc),xp_new,vpx_new,&
       vpy_new,vpz_new,Eki,px,sum_dEk(n_icp,nproc)
  ! Particle densities
  real(kind=8):: k1,kq(0:n(1)+2,0:n(2)+2),&
       np(0:n(1)+2,0:n(2)+2,ntype,nproc)
  ! Macroscopic parameters
  integer:: np_tot(ntype,nproc)
  
  do i=1,np_tot(ptype,iproc),1
     
     !
     ! Save actual velocity & position (avoid calling arrays many times)
     !
     xp_new= vxp(1,i,ptype,iproc)
     vpx_new= vxp(4,i,ptype,iproc)
     vpy_new= vxp(5,i,ptype,iproc)
     vpz_new= vxp(6,i,ptype,iproc)
     
     !
     ! Get particle left grid index
     !
     ix= INT( xp_new/h(1) ) + 1
     iy= 1
     
     !
     ! Calculate density
     !
     px=( ix*h(1) - xp_new )/h(1)
     k1=Nm(ptype)/h(1)
     ki(1)= k1*px
     ki(2)= k1*(1.d0-px)
     
     ! Charge assigned to the grid node i
     np(ix,iy,ptype,iproc)= np(ix,iy,ptype,iproc) + &
          kq(ix,iy)*ki(1)
     
     ! Charge assigned to the grid node i+1
     np(ix+1,iy,ptype,iproc)= np(ix+1,iy,ptype,iproc) + &
          kq(ix+1,iy)*ki(2)
     
     ! Electrons Maxwellian heating
     if( ptype.eq.1 .and. Pabs(1).gt.0.d0 ) then
        if( xp_new.ge.xl_pow .or. xp_new.le.xr_pow  ) then
           ! Store number of particles
           Nh(1,iproc)= Nh(1,iproc) + 1
           if(eheat_type.eq.1) then ! Replace particle velocity
              ! Calculate kinetic energy of macroparticle
              Eki= 0.5d0*Nm(ptype)*mass(ptype)*( vpx_new*vpx_new + &
                   vpy_new*vpy_new + vpz_new*vpz_new )
              ! Store energy
              sum_dEk(1,iproc)= sum_dEk(1,iproc) + Eki
           endif
        endif
     endif
      
  enddo
  
end subroutine charge_deposition1D

