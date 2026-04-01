subroutine part_flux_emitter(istep,vxp,n,h,ntype,nmax,bcnd,ss2D,nproc,&
     iseed,np_tot,nproc_mpi,mpi_rank,cnt_b,n_icp,P_loss)
!     ===================================================================
!     VERSION:         0.7
!     LAST MOD:      Nov/23
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:     
!     NOTES:    
!     -------------------------------------------------------------------
  use omp_lib 
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer mpi_rank,nproc_mpi ! MPI
  integer:: istep,n(3),nmax,ntype,Nh,iproc,nproc,iseed(nproc),iseed_OMP,n_icp
  integer:: bcnd(0:n(1)+2,0:n(2)+2),np_tot(ntype,nproc),cnt_b(2,nproc)
  real(kind=8):: h(3),vxp(6,nmax,ntype,nproc),js,ran2,rnd(2),dNh,&
       th_chrd,d_chrd,x1,x2,y1,y2,P_loss(4,ntype,nproc),I_sec,&
       ss2D(2,0:n(1)+2,0:n(2)+2,ntype,nproc)
  
  ! Initialization
  js= jne*10.d0 ! A/m-2
  Lgrd=n_icp*(Lg-Lh*n_holes) ! effective grid length

  ! Chamfered grid
  if(ptype_fx.eq.tag_neg .and. flag_chmfrd.eq.1) then
     x1=  chmfrd_array(1)
     x2=  chmfrd_array(2)
     y1=  chmfrd_array(3)
     y2=  chmfrd_array(4)
     d_chrd= dsqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) )
     th_chrd= dacos( (x2-x1)/d_chrd )
     Lgrd=2.d0*(y1+d_chrd) !  effective grid length

     if( istep.eq.1 .and. mpi_rank.eq.0) print'(1x,"Aperture is chamfered: L_ch(mm)=",f6.2,&
          ", th(degrees)=",f6.1,", Lgrd(mm)=",f6.2)',2.d0*d_chrd*1.d3,th_chrd*180./pi,Lgrd*1.d3

     if(Lgrd.gt.ymax) then
        if(mpi_rank.eq.0) print*, 'Lgrd>ymax, please correct...'
        call stop_calculation
     endif
  endif

  ! Number of ions per MPI proc per OMP proc
  I_sec= js*zmax*Lgrd
  dNh=I_sec*(real(ns_Hm)*dt)/(Nm(ptype_fx)*qe)/real(nproc_mpi)/real(nproc)
  Nh= INT(dNh)
  
  ! Correct for round-off errors
  rnd(1)= ran2(iseed(1))
  if( rnd(1).le.(dNh-Nh) ) Nh= Nh + 1

  !$OMP PARALLEL PRIVATE(iproc,iseed_OMP)
  ! Get processor id (from 0 to nproc-1)
  iproc= omp_get_thread_num() + 1
  iseed_OMP= iseed(iproc)
  call part_flux_OMP(vxp,n,h,ntype,nmax,bcnd,ss2D,nproc,iseed_OMP,np_tot,Nh,&
       iproc,th_chrd,cnt_b,n_icp,P_loss)
  iseed(iproc)= iseed_OMP
  !$OMP END PARALLEL
  
  ! Write on screen
  if( istep.eq.1 .and. mpi_rank.eq.0) then
     print'(1x,"Particle current generated on the electrode: Isurf(A)= ",&
          es10.2,", jsurf(A/m2)= ",f6.1)',dNh*real(nproc*nproc_mpi)*qe*Nm(ptype_fx)/(real(ns_Hm)*dt),js
  endif
  
  return
end subroutine part_flux_emitter

subroutine part_flux_OMP(vxp,n,h,ntype,nmax,bcnd,ss2D,nproc,&
     iseed,np_tot,Nh,iproc,th_chrd,cnt_b,n_icp,P_loss)
!     ===================================================================
!     VERSION:         0.8
!     LAST MOD:      Nov/23
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:     
!     NOTES:    
!     -------------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: n(3),nmax,ntype,ix,iy,i,Nh,iproc,nproc,iseed,ptype,n_icp
  integer:: bcnd(0:n(1)+2,0:n(2)+2),np_tot(ntype,nproc),np_tot_tmp,&
       cnt_b(2,nproc),cnt_lhs,cnt_rhs
  real(kind=8):: h(3),vxp(6,nmax,ntype,nproc),ran2,rnd(2),x_inj,y_inj,z_inj,&
       vt,vx,vy,vz,th_chrd,th_chrdtmp,x_injtmp,y1,d0,d1,P_loss(4,ntype,nproc),dir,&
       vx_tmp,vy_tmp,vz_tmp,dt_tmp,ss2D(2,0:n(1)+2,0:n(2)+2,ntype,nproc)
  
  ! Inject @ xg1
  x_injtmp= xg1

  ! Initialize
  cnt_lhs=0
  cnt_rhs=0

  ! Chamfered grid height
  if(ptype_fx.eq.tag_neg .and. flag_chmfrd.eq.1) y1=  chmfrd_array(3)

  ! Loop over negative ions
  do i=1,Nh
   
90   x_inj= x_injtmp   
     rnd(1)= ran2(iseed)
     if(n_icp.eq.1) then
        y_inj= (ymax-Lg)/2.d0 + Lg*rnd(1)
     else
        y_inj= ymax*(1.d0-rnd(1))
     endif

     ! Define ptype
     ptype= ptype_fx

     ! Chamfered surface model only for negative ions
     if( ptype_fx.eq.tag_neg .and. flag_chmfrd.eq.1 .and. ((y_inj.le.ymax/2.d0 .and. &
          y_inj.gt.y1) .or. (y_inj.gt.ymax/2.d0 .and. ABS(ymax-y_inj).gt.y1)) ) then
        ! Define ptype
        if(tag_b.gt.0) ptype= tag_b
        ! Angle with respect to the (Ox) axis: bottom side
        th_chrdtmp= th_chrd 
        ! Lengths
        d0= y_inj-y1
        d1= y1
        if(y_inj.gt.ymax/2.d0) then 
           ! Angle with respect to the (Ox) axis: top side (-theta)
           th_chrdtmp= -th_chrd 
           ! Lengths
           d0= ABS(y_inj-ymax) - y1
           d1= ymax - y1
        endif
        x_inj= x_injtmp + d0*dcos(th_chrdtmp)
        y_inj= d1 + d0*dsin(th_chrdtmp)
    endif

    ix= INT( x_inj/h(1) ) + 1
    iy= INT( y_inj/h(2) ) + 1
     
    rnd(1)= ran2(iseed)
    z_inj= zmax*(1.d0-rnd(1))
    
    if( (bcnd(ix,iy).ge.1 .and. &
         bcnd(ix,iy+1).ge.1) .or. &
         (bcnd(ix+1,iy).ge.1 .and. &         
         bcnd(ix+1,iy+1).ge.1)  ) then
       ! Find direction of flux
       if(bcnd(ix-1,iy).eq.-1) cnt_lhs= cnt_lhs+1
       if(bcnd(ix+1,iy).eq.-1) cnt_rhs= cnt_rhs+1
       if(cnt_lhs.ne.cnt_rhs) then
          if(cnt_lhs.gt.cnt_rhs) then
             dir=1
          else
             dir= -1
          endif
       else
          ! Failed to find direction, look for another macroparticle
          goto 90
       endif
    else
       ! Do not load inside apertures
       goto 90
    endif
    
     ! Count current ratio chamfered vs. flat surface
     if(ptype.eq.tag_b) then
        cnt_b(2,iproc)= cnt_b(2,iproc) + 1 
     else
        cnt_b(1,iproc)= cnt_b(1,iproc) + 1 
     endif

     ! Create ion
     np_tot(ptype,iproc)= np_tot(ptype,iproc) + 1
     np_tot_tmp= np_tot(ptype,iproc)

     ! Do not reset sour 
     ss2D(1,ix,iy,ptype,iproc)= ss2D(1,ix,iy,ptype,iproc) + 1.
          
     ! Gaussian loading (flux normal to the PE)
     vt= dsqrt(2.d0*qe*THm/ABS(mass(ptype))) 
     
     rnd(1)=ran2(iseed)
     vx = dir*vt*dsqrt( -dlog(1-rnd(1)) )
     
     rnd(1)= ran2(iseed)
     rnd(2)= ran2(iseed)
     call load_gauss(vy,vz,vt,rnd)
     
     ! Rotation with respect to the (Ox) axis
     th_chrdtmp= pi/2.d0 ! Normal to the left-hand-side
     if( ptype_fx.eq.tag_neg .and. flag_chmfrd.eq.1 .and. ((y_inj.le.ymax/2.d0 .and. &
          y_inj.gt.y1) .or. (y_inj.gt.ymax/2.d0 .and. ABS(ymax-y_inj).gt.y1)) ) then
        th_chrdtmp= th_chrd  ! Normal to the bottom chamfered grid surface
        if(y_inj.gt.ymax/2.d0) th_chrdtmp= pi/2.d0+th_chrd ! top chamfered grid side
     endif

     vx_tmp= -vx*dsin(th_chrdtmp) + vy*dcos(th_chrdtmp) 
     vy_tmp= vx*dcos(th_chrdtmp) + vy*dsin(th_chrdtmp)
     vz_tmp= vz
     
     vxp(4,np_tot_tmp,ptype,iproc)= vx_tmp
     vxp(5,np_tot_tmp,ptype,iproc)= vy_tmp
     vxp(6,np_tot_tmp,ptype,iproc)= vz_tmp

     rnd(1)=ran2(iseed)
     dt_tmp= rnd(1)*(real(ns_Hm)*dt)
     vxp(1,np_tot_tmp,ptype,iproc)= x_inj + vx_tmp*dt_tmp
     vxp(2,np_tot_tmp,ptype,iproc)= y_inj + vy_tmp*dt_tmp
     vxp(3,np_tot_tmp,ptype,iproc)= z_inj + vz_tmp*dt_tmp

     P_loss(4,ptype,iproc)= P_loss(4,ptype,iproc) + 0.5d0*Nm(ptype)*mass(ptype)*( &
          vxp(4,np_tot_tmp,ptype,iproc)*vxp(4,np_tot_tmp,ptype,iproc) + &
          vxp(5,np_tot_tmp,ptype,iproc)*vxp(5,np_tot_tmp,ptype,iproc) + &
          vxp(6,np_tot_tmp,ptype,iproc)*vxp(6,np_tot_tmp,ptype,iproc) ) 
     
  enddo ! end-loop over Nh

  return
end subroutine part_flux_OMP
