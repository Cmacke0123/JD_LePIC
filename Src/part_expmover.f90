subroutine part_mover(istep,n,na,h,Ei,Bi,p_mac,P_loss,vxp,bcnd,&
     nmax,ntype,ngrid,flag_dead,nproc,np_tot,phi,uB,iproc,ptype,&
     iseed,cnt_dead,sum_q,sum_q_y,phase_space,ip_ps,nmax_ps,n_B,h_B,ss2D)
!     ==============================================================
!     VERSION:         0.12
!     LAST MOD:      Mar/24
!     MOD AUTHOR
!     NOTE:       Interpolation |------|-------------|
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
  integer:: n(3),na(3),n_B(2),nmax,ntype,ngrid,igrid
  integer:: ix,iy,ix_c,iy_c,ix_w,iy_w,i,ptype,istep,flag_B_tmp,flag_lost,iproc,&
       nproc,nmax_ps,ip_ps(nproc),ip_ps_tmp,igrid_tmp,dtype_tmp,ptype_sec,&
       ip_sec,n_sec
  ! Field arrays
  real(kind=8):: h(3),ha(3),h_B(2),Ei(2,0:n(1)+2,0:n(2)+2),Exp,Eyp,ki(4), &
       Bi(4,0:n_B(1)+2,0:n_B(2)+2)
  ! Particle arrays
  integer:: bcnd(0:n(1)+2,0:n(2)+2)
  integer(kind=1):: flag_dead(nmax,ntype,nproc)
  real(kind=8):: vxp(6,nmax,ntype,nproc),xp_new,yp_new,zp_new,vpx_new,&
       vpy_new,vpz_new,Eki,px,py,&
       sum_q_y(ngrid,0:n(2)+2,ntype,nproc),vx_sec,vy_sec,vz_sec,vt
  real(kind=8):: k1,k2,k3,k4,phi(0:n(1)+2,0:n(2)+2),&
       uB(na(1)+1,na(2)+1,ntype),phase_space(nmax_ps,4,nproc)
  real(kind=8):: v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,Bpx,Bpy,Bpz,B_tot,phip,&
       vb(ntype),vmax(ntype),fmax(ntype),vbx,dt_tmp
  ! Macroscopic parameters
  integer:: np_tot(ntype,nproc),np_lost(ntype,nproc),i_shift,iseed
  real(kind=8):: p_mac(ntype,2,0:ngrid,nproc),P_loss(4,ntype,nproc),&
       ran2,rnd(2),nu_3D,cnt_dead(nproc),sum_q(ntype,nproc),&
       ss2D(2,0:n(1)+2,0:n(2)+2,ntype,nproc)

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

  ! Grid size for uB()
  ha=size_na*h

  !
  ! Move particles
  !
  do i=1,np_tot(ptype,iproc),1

     !
     ! Save actual velocity & position (avoid calling arrays many times)
     !
     xp_new= vxp(1,i,ptype,iproc)
     yp_new= vxp(2,i,ptype,iproc)
     zp_new= vxp(3,i,ptype,iproc)
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
     iy= INT( yp_new/h(2) ) + 1

     !
     ! Calculate fields at position x,y
     !
     px=( ix*h(1) - xp_new )/h(1)
     py=( iy*h(2) - yp_new )/h(2)

     ki(1)= px*py
     ki(2)= (1.d0-px)*py
     ki(3)= (1.d0-px)*(1.d0-py)
     ki(4)= px*(1.d0-py)

     Exp= ki(1)*Ei(1,ix,iy) + &
          ki(2)*Ei(1,ix+1,iy) + &
          ki(3)*Ei(1,ix+1,iy+1) + &
          ki(4)*Ei(1,ix,iy+1)
     
     Eyp= ki(1)*Ei(2,ix,iy) + &
          ki(2)*Ei(2,ix+1,iy) + &
          ki(3)*Ei(2,ix+1,iy+1) + &
          ki(4)*Ei(2,ix,iy+1)
     
     !
     ! Calculate new velocity
     !
     flag_B_tmp=0
     if(flag_B.eq.1)then

        if(n_B(1).eq.1) then ! Constant B-field

           Bpx= Bi(1,1,1)
           Bpy= Bi(2,1,1)
           Bpz= Bi(3,1,1)
           
           flag_B_tmp=1
        else ! B-field map
           
           if(flag_gridB.eq.1) then
              ! Get particle left grid index on B-field map
              ix= INT( xp_new/h_B(1) ) + 1
              iy= INT( yp_new/h_B(2) ) + 1
           endif
        
           if( Bi(4,ix,iy).gt.1.d-5 ) then

              if(flag_gridB.eq.1) then
                 ! Calculate B-fields at position x and y on the map
                 px=( ix*h_B(1) - xp_new )/h_B(1)
                 py=( iy*h_B(2) - yp_new )/h_B(2)

                 ki(1)= px*py
                 ki(2)= (1.d0-px)*py
                 ki(3)= (1.d0-px)*(1.d0-py)
                 ki(4)= px*(1.d0-py)
              endif

              ! Bx
              Bpx= ki(1)*Bi(1,ix,iy) + &
                   ki(2)*Bi(1,ix+1,iy) + &
                   ki(3)*Bi(1,ix+1,iy+1) + &
                   ki(4)*Bi(1,ix,iy+1)
              
              ! By
              Bpy= ki(1)*Bi(2,ix,iy) + &
                   ki(2)*Bi(2,ix+1,iy) + &
                   ki(3)*Bi(2,ix+1,iy+1) + &
                   ki(4)*Bi(2,ix,iy+1)
           
              ! Bz
              Bpz= ki(1)*Bi(3,ix,iy) + &
                   ki(2)*Bi(3,ix+1,iy) + &
                   ki(3)*Bi(3,ix+1,iy+1) + &
                   ki(4)*Bi(3,ix,iy+1)
  
              flag_B_tmp=1           
           endif
        endif
     endif

     ! Positive ions not magnetized
     if( flag_B_pos.eq.1 .and. charge(ptype).gt.0.d0 ) flag_B_tmp=0

     if(flag_B_tmp.eq.1) then ! Boris scheme

        B_tot= dsqrt( Bpx*Bpx + Bpy*Bpy + Bpz*Bpz )         
        k2= k1*( 2.d0/(1.d0 + (k1*B_tot)*(k1*B_tot)) )

        v1x= vpx_new + k1*Exp
        v1y= vpy_new + k1*Eyp
        v1z= vpz_new

        v3x= v1x + k1*( v1y*Bpz - v1z*Bpy )
        v3y= v1y + k1*( v1z*Bpx - v1x*Bpz )
        v3z= v1z + k1*( v1x*Bpy - v1y*Bpx )

        v2x= v1x + k2*( v3y*Bpz - v3z*Bpy )
        v2y= v1y + k2*( v3z*Bpx - v3x*Bpz )
        v2z= v1z + k2*( v3x*Bpy - v3y*Bpx ) 

        vpx_new= v2x + k1*Exp
        vpy_new= v2y + k1*Eyp
        vpz_new= v2z

     else ! Electrostatic Leap-Frog solver

        vpx_new= vpx_new + 2.d0*k1*Exp
        vpy_new= vpy_new + 2.d0*k1*Eyp
           
     endif

     ! Calculate new position  
     xp_new= xp_new + dt*vpx_new
     yp_new= yp_new + dt*vpy_new
     ! Track z-location in 2.5D
     zp_new= zp_new + dt*vpz_new

     ! Periodic BCs
     k4= INT(yp_new/ymax) ! 0 (y<0) or 1 usually
     if(flag_pbc.eq.1) then
        if( yp_new.ge.ymax ) then
           ! Re-inject particle at bottom of simulation box
           if(flag_spec.eq.0) then
              yp_new= yp_new - k4*ymax
           else ! Specular reflection 
              yp_new= (k4+1)*ymax - yp_new
              vpy_new= -vpy_new
           endif
        endif
        if( yp_new.le.0.d0 ) then 
           if(flag_spec.eq.0) then
              yp_new= (1-k4)*ymax + yp_new
           else
              yp_new= k4*ymax - yp_new
              vpy_new= -vpy_new
           endif
        endif
     endif

     ! Get particle new grid location
     ix= FLOOR( xp_new/h(1) ) + 1
     iy= FLOOR( yp_new/h(2) ) + 1

     ! Correct for particles out of bounds
     call check_outofbounds(ix,iy,n)

     !
     ! Check if the particle did not leave the simulation domain
     !
     flag_lost= 0
     if( bcnd(ix,iy).ge.1 .and. bcnd(ix+1,iy).ge.1 .and. &
          bcnd(ix+1,iy+1).ge.1 .and. bcnd(ix,iy+1).ge.1 ) then 
        flag_lost= 1
        goto 120
     endif

     ! Negative ions and Neumann BC's (no specular reflexion)
     if( charge(ptype).lt.0 .and. mass(ptype).ge.amu ) then
        if( xp_new.lt.0.d0 .and. flag_nmn.eq.1 ) then 
           flag_lost=3 
           goto 120
        endif
     endif

     !
     ! 3rd virtual dimension
     !
     if(hz.lt.0) goto 120 ! No particle losses

     ! Electrons & negative ions
     if(charge(ptype).lt.0) then
              
        if( zp_new.le.0.d0 .or. zp_new.ge.zmax ) then               
           ! Calculate kinetic energy of macroparticle
           Eki= 0.5d0*mass(ptype)*vpz_new*vpz_new/qe

           ! Check if negative charge hits the wall
           if( Eki.ge.(phi(ix,iy)-phi0) ) then 
              flag_lost= 2
           else ! Particle is reflected back
              vpz_new= -vpz_new
              if(zp_new.le.0.d0) then ! Modify z-location
                 zp_new=0.02d0*zmax
              else
                 zp_new=0.98d0*zmax
              endif
           endif
        endif
        
        ! Positive ions (Bohm frequency)
     else
        
        if( na(1).eq.n(1) ) then
           ix_c= ix
           iy_c= iy
        else ! coarse grid coordinates
           ix_c= INT( xp_new/ha(1) ) + 1
           iy_c= INT( yp_new/ha(2) ) + 1
        endif

        ! hz*uB*2/L
        nu_3D=hz*2.d0*uB(ix_c,iy_c,ptype)/zmax
        rnd(1)= ran2(iseed)
        if( rnd(1).le.nu_3D*dt ) flag_lost= 2
        
     endif

     !
     ! Lost particles
     !
120  if(flag_lost.ge.1) then

        ! Add another lost particle
        np_lost(ptype,iproc)= np_lost(ptype,iproc) + 1
        
        ! Calculate kinetic energy of macroparticle
        Eki= 0.5d0*Nm(ptype)*mass(ptype)*( vpx_new*vpx_new + &
             vpy_new*vpy_new + vpz_new*vpz_new ) 
        
        ! Total lost power per time step
        P_loss(1,ptype,iproc)= P_loss(1,ptype,iproc) + Eki

        ! Sink term
        if(plt_src.eq.0 .and. flag_lost.eq.1) then
           ix_w= NINT( xp_new/h(1) ) + 1
           iy_w= NINT( yp_new/h(2) ) + 1
           call check_outofbounds(ix_w,iy_w,n)
           ss2D(2,ix_w,iy_w,ptype,iproc)= ss2D(2,ix_w,iy_w,ptype,iproc) + 1.
        endif

        ! Particle lost on Dirichlet boundaries
        if(flag_lost.eq.1) then
        
           ! Get grid index
           igrid= bcnd(ix,iy)
        
           ! Dielectrics        
           if(flag_dielec.eq.1) then
              dtype_tmp= dtype(igrid)
              igrid_tmp= bcnd(ix,iy+1)
              dtype_tmp= MAX(dtype_tmp,dtype(igrid_tmp))
              
              if( dtype_tmp.eq.1 ) then
                 k3= Nm(ptype)*charge(ptype)
                 py=( iy*h(2) - yp_new )/h(2)
                 sum_q_y(igrid,iy,ptype,iproc)= sum_q_y(igrid,iy,ptype,iproc) + k3*py
                 sum_q_y(igrid,iy+1,ptype,iproc)= sum_q_y(igrid,iy+1,ptype,iproc) + k3*(1.d0-py)           
              endif
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
                    vt= dsqrt(2.d0*qe*ABS(THm)/ABS(mass(ptype_sec))) 
                    rnd(1)=ran2(iseed)
                    !  dir_sec or -sign(1.d0,vpx_new)
                    vx_sec = -sign(1.d0,vpx_new)*vt*dsqrt( -dlog(1-rnd(1)) )
                    rnd(1)= ran2(iseed)
                    rnd(2)= ran2(iseed)
                    ! Gaussian loading (flux normal to the grid surface)
                    call load_gauss(vy_sec,vz_sec,vt,rnd)
                    
                    ! y-location of impacting ion.
                    vxp(1,i_shift,ptype_sec,iproc)= xg_sec
                    vxp(2,i_shift,ptype_sec,iproc)= yp_new 
                    vxp(3,i_shift,ptype_sec,iproc)= zp_new 
                    vxp(4,i_shift,ptype_sec,iproc)= vx_sec 
                    vxp(5,i_shift,ptype_sec,iproc)= vy_sec 
                    vxp(6,i_shift,ptype_sec,iproc)= vz_sec
                    
                    ! Count power injected into the plasma electrons
                    P_loss(4,ptype_sec,iproc)= P_loss(4,ptype_sec,iproc) + 0.5d0*Nm(ptype_sec)*mass(ptype_sec)*( &
                         vx_sec*vx_sec + vz_sec*vz_sec + vz_sec*vz_sec )
                    
                 enddo
130              continue
              endif
           endif
        endif

        ! 3rd dimension
        if(flag_lost.ge.2) then  
           igrid=0
           ! Floating end-plates
           if(flag_float.eq.1) sum_q(ptype,iproc)= sum_q(ptype,iproc) + Nm(ptype)*charge(ptype)

           if(flag_lost.eq.2) then ! 2.5D
              px=( ix*h(1) - xp_new )/h(1)
              py=( iy*h(2) - yp_new )/h(2)
              
              ki(1)= px*py
              ki(2)= (1.d0-px)*py
              ki(3)= (1.d0-px)*(1.d0-py)
              ki(4)= px*(1.d0-py)
              
              phip= ki(1)*phi(ix,iy) + &
                   ki(2)*phi(ix+1,iy) + &
                   ki(3)*phi(ix+1,iy+1) + &
                   ki(4)*phi(ix,iy+1)

              ! Power lost when removing particles inside a potentiel well (U*I)
              P_loss(1,ptype,iproc)= P_loss(1,ptype,iproc) + Nm(ptype)*charge(ptype)*phip
           endif
   
        endif

        ! Total number of particle lost at the wall per time step
        p_mac(ptype,np_loss,igrid,iproc)= p_mac(ptype,np_loss,igrid,iproc) + 1
        
        ! Total lost power at the wall per time step
        p_mac(ptype,P_w,igrid,iproc)= p_mac(ptype,P_w,igrid,iproc) + Eki
        
        ! Ion properties impacting num_grd
        if(flag_ps.eq.1) then
           if(igrid.eq.num_grd .and. mass(ptype).ge.amu) then              
              ! Save phase-space distribution  
              if(ip_ps(iproc).lt.nmax_ps) then                        
                 ip_ps(iproc)= ip_ps(iproc) + 1
                 ip_ps_tmp= ip_ps(iproc)
                 phase_space(ip_ps_tmp,1,iproc)= yp_new
                 phase_space(ip_ps_tmp,2,iproc)= vpx_new
                 phase_space(ip_ps_tmp,3,iproc)= vpy_new
                 phase_space(ip_ps_tmp,4,iproc)= ptype
              else
                 flag_ps=2
              endif
           endif
        endif

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
        vxp(2,i_shift,ptype,iproc)= yp_new
        vxp(3,i_shift,ptype,iproc)= zp_new
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
 
end subroutine part_mover


subroutine charge_deposition(n,h,vxp,nmax,ntype,kq,np,nproc,np_tot,sum_dEk,&
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
  integer:: ix,iy,i,ptype,iproc,nproc,n_icp,i_icp
  real(kind=8):: h(3),ki(4),vxp(6,nmax,ntype,nproc),xp_new,yp_new,zp_new,vpx_new,&
       vpy_new,vpz_new,Eki,px,py,sum_dEk(n_icp,nproc)
  ! Particle densities
  real(kind=8):: k1,kq(0:n(1)+2,0:n(2)+2),&
       np(0:n(1)+2,0:n(2)+2,ntype,nproc)
  ! Macroscopic parameters
  integer:: np_tot(ntype,nproc)

  do i=1,np_tot(ptype,iproc),1
     
     ! 
     ! Save particle 6D-coordinates
     !
     xp_new= vxp(1,i,ptype,iproc)
     yp_new= vxp(2,i,ptype,iproc)
     zp_new= vxp(3,i,ptype,iproc)
     vpx_new= vxp(4,i,ptype,iproc)
     vpy_new= vxp(5,i,ptype,iproc)
     vpz_new= vxp(6,i,ptype,iproc)
     
     !
     ! Get particle left grid index
     !
     ix= INT( xp_new/h(1) ) + 1
     iy= INT( yp_new/h(2) ) + 1
     
     
     !
     ! Calculate density
     !
     px=( ix*h(1) - xp_new )/h(1)
     py=( iy*h(2) - yp_new )/h(2)
     k1=Nm(ptype)/(h(1)*h(2)*zmax)
     ki(1)= k1*px*py
     ki(2)= k1*(1.d0-px)*py
     ki(3)= k1*(1.d0-px)*(1.d0-py)
     ki(4)= k1*px*(1.d0-py)
     
     ! Charge assigned to the grid node 00
     np(ix,iy,ptype,iproc)= np(ix,iy,ptype,iproc) + &
          kq(ix,iy)*ki(1)
     
     ! Charge assigned to the grid node +0
     np(ix+1,iy,ptype,iproc)= np(ix+1,iy,ptype,iproc) + &
          kq(ix+1,iy)*ki(2)
     
     ! Charge assigned to the grid node ++
     np(ix+1,iy+1,ptype,iproc)= np(ix+1,iy+1,ptype,iproc) + &
          kq(ix+1,iy+1)*ki(3)
     
     ! Charge assigned to the grid node 0+
     np(ix,iy+1,ptype,iproc)= np(ix,iy+1,ptype,iproc) + &
          kq(ix,iy+1)*ki(4)
     
     ! Electrons Maxwellian heating
     if( ptype.eq.1 .and. Pabs(1).gt.0.d0 ) then
        if(flag_c.eq.0) then ! Slit 
           if( xp_new.lt.xl_pow .or. xp_new.gt.xr_pow .or. &
                yp_new.lt.yl_pow .or. yp_new.gt.yr_pow ) goto 140
        else ! Disk
           if( ((xp_new-xa)**2 + (yp_new-ymax/2.d0)**2).gt.dr**2 ) goto 140
        endif
        ! Index of driver 
        i_icp=INT(n_icp*yp_new/ymax) + 1
        ! Store number of particles
        Nh(i_icp,iproc)= Nh(i_icp,iproc) + 1
        if(eheat_type.eq.1) then ! Replace particle velocity
           ! Calculate kinetic energy of macroparticle
           Eki= 0.5d0*Nm(ptype)*mass(ptype)*( vpx_new*vpx_new + &
                vpy_new*vpy_new + vpz_new*vpz_new )
           ! Store energy
           sum_dEk(i_icp,iproc)= sum_dEk(i_icp,iproc) + Eki
        endif
140     continue
     endif
     
  enddo
  
end subroutine charge_deposition

subroutine check_outofbounds(ix,iy,n)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Apr/24
!     MOD AUTHOR     G. Fubiani
!     --------------------------------------------------------------
  implicit none
  integer ix,iy,n(3)

  if(ix.lt.0) ix=0
  if(ix.gt.n(1)+1) ix=n(1)+1
  if(iy.lt.0) iy=0
  if(iy.gt.n(2)+1) iy=n(2)+1

  return
end subroutine check_outofbounds
