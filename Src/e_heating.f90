subroutine eheating(vxp,nmax,ntype,nproc,iseed,np_tot,vt,&
     n_icp,iproc,P_loss,Pabs_cor)
!     ==============================================================
!     VERSION:         0.5
!     LAST MOD:      Jan/24
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:   Electrons Maxwellian heating
!                 1== xp
!                 4== vpx
!                 5== vpy
!                 6== vpz
!
!                Estimates for nu - Gaussian power profile
!                <nu>= Int_D n(x,y)*nu*exp[-(x^2+y^2)/L^2]*dx*dy/ 
!                      Int_D n(x,y)*dx*dy 
!                n(x,y)=Sum_i dirac_delta(x-xi,y-yi), D is a disk of
!                radius dr. if n(x,y)= cste, Int_D dx*dy~ pi*dr**2
!                dr, L, <nu> are input parameters, L= dr/3 (3 sigma).
!                We get nu= 9*<nu>. For arbitrary n(x,y) profiles, nu is  
!                calculated iteratively using the previous value of the 
!                absorbed power (correcting factor Pabs_cor).
!                
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  include 'constants.h'
  integer:: nmax,ntype,iseed,np_tot(ntype,nproc),np_tot_tmp,&
       i,ptype,iproc,nproc,n_icp,i_icp
  real(kind=8):: vxp(6,nmax,ntype,nproc),xp_new,yp_new,vt(n_icp),&
       vz_sav,ran2,rnd(2),rp2,v2old,P_loss(4,ntype,nproc),dvx,dvy,dvz,&
       nu_tmp,Pabs_cor,flag_eh

  ! Only electrons
  ptype= 1

  ! Initialization
  vz_sav= 0.d0
  nu_tmp= nudt/dt ! nu_h*(ns_heat*dt)

  !
  ! Maxwellian heating
  !
  np_tot_tmp= np_tot(ptype,iproc)
  do i=1,np_tot_tmp,1

     ! Get particle left grid index
     xp_new= vxp(1,i,ptype,iproc)
     yp_new= vxp(2,i,ptype,iproc)

     ! Add heating by external source
     flag_eh=0
     if(flag_c.eq.0) then ! Slit
        rp2= (xp_new-xa)**2
        if( rp2.le.dr**2 .and. yp_new.ge.yl_pow .and. yp_new.le.yr_pow ) flag_eh=1
     else ! Disk
        rp2= (xp_new-xa)**2 + (yp_new-ymax/2.d0)**2
        if(rp2.le.dr**2) flag_eh=1
     endif
     
     if(flag_eh.eq.1) then
        ! Option for gaussian power profile
        if(flag_gpp.eq.1) nu_tmp=(nudt/dt)*9.d0*Pabs_cor*dexp(-9.d0*rp2/dr**2)
     else
        goto 100
     endif

     ! Use Maxwellian with self-consistent temperature 
     rnd(1)= ran2(iseed)
     if( rnd(1).le.nu_tmp*dt ) then
             
        v2old= vxp(4,i,ptype,iproc)*vxp(4,i,ptype,iproc) + &
             vxp(5,i,ptype,iproc)*vxp(5,i,ptype,iproc) + &
             vxp(6,i,ptype,iproc)*vxp(6,i,ptype,iproc)

        i_icp=INT(n_icp*yp_new/ymax) + 1

        ! Calculate new velocity 
        rnd(1)= ran2(iseed)
        rnd(2)= ran2(iseed)
        call load_gauss(dvx,dvy,vt(i_icp),rnd)
        if(vz_sav.eq.0.d0) then
           rnd(1)= ran2(iseed)
           rnd(2)= ran2(iseed)
           call load_gauss(dvz,vz_sav,vt(i_icp),rnd)
        else
           dvz= vz_sav
           vz_sav= 0.d0
        endif

        if(eheat_type.eq.1) then
           vxp(4,i,ptype,iproc)=  dvx
           vxp(5,i,ptype,iproc)=  dvy
           vxp(6,i,ptype,iproc)=  dvz
        else
           vxp(4,i,ptype,iproc)=  vxp(4,i,ptype,iproc) + dvx
           vxp(5,i,ptype,iproc)=  vxp(5,i,ptype,iproc) + dvy
           vxp(6,i,ptype,iproc)=  vxp(6,i,ptype,iproc) + dvz
        endif
        
        P_loss(2,ptype,iproc)= P_loss(2,ptype,iproc) + 0.5d0*Nm(ptype)*mass(ptype)*( &
             vxp(4,i,ptype,iproc)*vxp(4,i,ptype,iproc) + &
             vxp(5,i,ptype,iproc)*vxp(5,i,ptype,iproc) + &
             vxp(6,i,ptype,iproc)*vxp(6,i,ptype,iproc) - v2old ) 

     endif
     
100  continue
  enddo ! end loop over np_tot(ptype) particles

  return
 
end subroutine eheating
