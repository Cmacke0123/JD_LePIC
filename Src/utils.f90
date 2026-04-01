subroutine part_moments(n,na,h,vxp,nmax,ntype,kq,nproc,np_tot, &
     iproc,ptype,p_mts,n_mts)
!     ==============================================================
!     VERSION:         0.3
!     LAST MOD:      APR/20
!     MOD AUTHOR
!     NOTE:        Indexes in vxp(): 1== x
!                                    2== y    
!                                    3== z
!                                    4== vx
!                                    5== vy
!                                    6== vz                
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  integer:: ix,iy,i,ptype,n(3),na(3),nmax,ntype,iproc,nproc,n_mts, &
       Ekg,j1,j2,j3,ix_f,iy_f
  parameter ( Ekg=1, j1=2, j2=3, j3=4 ) 
  real(kind=8):: h(3),ha(3),ki(8),k1,k2,k_tmp
  ! Particle arrays
  real(kind=8):: vxp(6,nmax,ntype,nproc),xp,yp,vx,vy,vz,Eki,px,py,&
       kq(0:n(1)+2,0:n(2)+2),p_mts(n_mts,na(1)+1,na(2)+1,ntype,nproc)
  ! Macroscopic parameters
  integer:: np_tot(ntype,nproc)

  ! Initialization
  p_mts(:,:,:,ptype,iproc)=0.d0 
  ha= size_na*h
  if(flag_1D.eq.0) then 
     k1=Nm(ptype)/(ha(1)*ha(2)*zmax)
  else
     k1=Nm(ptype)/ha(1) ! per m^2
  endif

  ! Loop over ptype particles
  do i=1,np_tot(ptype,iproc)

     vx= vxp(4,i,ptype,iproc)
     vy= vxp(5,i,ptype,iproc)
     vz= vxp(6,i,ptype,iproc)
     xp= vxp(1,i,ptype,iproc) - vx*dt/2.d0
     yp= vxp(2,i,ptype,iproc)
     if(flag_1D.eq.0) yp= yp - vy*dt/2.d0
             
     ! Neumann BCs (LHS only)
     if(flag_nmn.eq.1) then
        if( xp.le.0.d0 ) then
           ! Specular reflection
           xp= -xp
           vx= -vx
        endif
     endif

     ! Periodic BCs (off in 1D)
     if(flag_pbc.eq.1) then
        k2= INT(yp/ymax) ! 0 (y<0) or 1 usually
        if( yp.ge.ymax ) then
           ! Re-inject particle at bottom of simulation box
           if(flag_spec.eq.0) then
              yp= yp - k2*ymax ! k2= 1 : -1.
           else ! Specular reflection 
              yp= (k2+1)*ymax - yp ! 2.
              vy= -vy
           endif             
        endif
        if( yp.le.0.d0 ) then 
           if(flag_spec.eq.0) then
              yp= (1-k2)*ymax + yp ! k2= 0 : 1.
           else
              yp= k2*ymax - yp ! 0.
              vy= -vy
           endif
        endif
     endif

     ! Calculate particle fluxes and kinetic energies
     if( xp.gt.xmax .or. xp.lt.0.d0 .or. &
          yp.gt.ymax .or. yp.lt.0.d0 ) goto 100

     Eki= vx*vx + vy*vy + vz*vz
     
     ix= INT( xp/ha(1) ) + 1
     px=( ix*ha(1) - xp )/ha(1)
             
     if(flag_1D.eq.0) then
        iy= INT( yp/ha(2) ) + 1
        py=( iy*ha(2) - yp )/ha(2)
     
        ki(1)= k1*px*py
        ki(2)= k1*(1.d0-px)*py
        ki(3)= k1*(1.d0-px)*(1.d0-py)
        ki(4)= k1*px*(1.d0-py)
     else
        ! 1D case
        iy= 1
        ki(1)= k1*px
        ki(2)= k1*(1.d0-px)
     endif

     ix_f= size_na*ix - (size_na-1)
     iy_f= size_na*iy - (size_na-1)

     ! Grid node 00
     k_tmp= kq(ix_f,iy_f)*ki(1)
     p_mts(j1,ix,iy,ptype,iproc)= p_mts(j1,ix,iy,ptype,iproc) + k_tmp*vx
     p_mts(j2,ix,iy,ptype,iproc)= p_mts(j2,ix,iy,ptype,iproc) + k_tmp*vy
     p_mts(j3,ix,iy,ptype,iproc)= p_mts(j3,ix,iy,ptype,iproc) + k_tmp*vz
     p_mts(Ekg,ix,iy,ptype,iproc)= p_mts(Ekg,ix,iy,ptype,iproc) + k_tmp*Eki

     ! Grid node +0
     k_tmp= kq(ix_f+1,iy_f)*ki(2)
     p_mts(j1,ix+1,iy,ptype,iproc)= p_mts(j1,ix+1,iy,ptype,iproc) + k_tmp*vx
     p_mts(j2,ix+1,iy,ptype,iproc)= p_mts(j2,ix+1,iy,ptype,iproc) + k_tmp*vy
     p_mts(j3,ix+1,iy,ptype,iproc)= p_mts(j3,ix+1,iy,ptype,iproc) + k_tmp*vz
     p_mts(Ekg,ix+1,iy,ptype,iproc)= p_mts(Ekg,ix+1,iy,ptype,iproc) + k_tmp*Eki

     if(flag_1D.eq.0) then
        ! Grid node ++
        k_tmp= kq(ix_f+1,iy_f+1)*ki(3)
        p_mts(j1,ix+1,iy+1,ptype,iproc)= p_mts(j1,ix+1,iy+1,ptype,iproc) + k_tmp*vx
        p_mts(j2,ix+1,iy+1,ptype,iproc)= p_mts(j2,ix+1,iy+1,ptype,iproc) + k_tmp*vy
        p_mts(j3,ix+1,iy+1,ptype,iproc)= p_mts(j3,ix+1,iy+1,ptype,iproc) + k_tmp*vz
        p_mts(Ekg,ix+1,iy+1,ptype,iproc)= p_mts(Ekg,ix+1,iy+1,ptype,iproc) + k_tmp*Eki
        
        ! Grid node 0+
        k_tmp= kq(ix_f,iy_f+1)*ki(4)
        p_mts(j1,ix,iy+1,ptype,iproc)= p_mts(j1,ix,iy+1,ptype,iproc) + k_tmp*vx
        p_mts(j2,ix,iy+1,ptype,iproc)= p_mts(j2,ix,iy+1,ptype,iproc) + k_tmp*vy
        p_mts(j3,ix,iy+1,ptype,iproc)= p_mts(j3,ix,iy+1,ptype,iproc) + k_tmp*vz
        p_mts(Ekg,ix,iy+1,ptype,iproc)= p_mts(Ekg,ix,iy+1,ptype,iproc) + k_tmp*Eki
     endif

100  enddo

  return
end subroutine part_moments

subroutine calc_avg(n,h,na,np,p_mts,data_pavg,ntype,n_mts,ss2D,&
     nproc,nproc_mpi)
!     ==============================================================
!     VERSION:         0.4
!     LAST MOD:      DEC/23
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:  Calculate average value of particle moments
!     NOTE: 
!     --------------------------------------------------------------
  implicit none
  include 'mpif.h'
  include 'particle_info.h'
  include 'constants.h'
  integer:: ix,iy,ix_f,iy_f,ptype,n(3),na(3),iproc,nproc,ntype,&
       n_mts,nproc_mpi,ierr
  integer:: np_avg,Tp_avg,j1_avg,j2_avg,j3_avg,sour_avg,sink_avg,&
       alpha,beta
  parameter ( np_avg=1, Tp_avg=2, j1_avg=3, j2_avg=4, j3_avg=5, &
       sour_avg=6, sink_avg=7 ) 
  real(kind=8):: h(3),np(0:n(1)+2,0:n(2)+2,ntype,nproc), &
       p_mts(n_mts,na(1)+1,na(2)+1,ntype,nproc),&
       data_pavg(7,na(1)+1,na(2)+1,ntype),&
       Vc,k,np_tmp,jx_tmp,jy_tmp,jz_tmp,Tp_tmp,&
       ss2D(2,0:n(1)+2,0:n(2)+2,ntype,nproc)
  real(kind=8),allocatable:: data_pavg_tmp(:,:,:,:),data_pavg_MPI(:,:,:,:)

  allocate( data_pavg_tmp(7,0:na(1)+2,0:na(2)+2,ntype) )

  ! Initialize
  data_pavg_tmp=0.d0
  
  alpha= size_na
  beta= size_na-1

  Vc= h(1)*h(2)*zmax
  if(flag_1D.eq.1) Vc= h(1) ! per m^2 in 1D

  !
  ! Reduction
  !
  do iproc=1,nproc
     do ptype=1,ntype

        ! Charge density element
        k=Nm(ptype)/Vc

        !$OMP PARALLEL
        !$OMP DO
        do iy=1,na(2)+1
           do ix=1,na(1)+1
              ix_f= alpha*ix - beta
              iy_f= alpha*iy - beta
              data_pavg_tmp(np_avg,ix,iy,ptype)= data_pavg_tmp(np_avg,ix,iy,ptype) + &
                   np(ix_f,iy_f,ptype,iproc) 
              data_pavg_tmp(j1_avg,ix,iy,ptype)= data_pavg_tmp(j1_avg,ix,iy,ptype) + &
                   p_mts(2,ix,iy,ptype,iproc) ! jx
              data_pavg_tmp(j2_avg,ix,iy,ptype)= data_pavg_tmp(j2_avg,ix,iy,ptype) + &
                   p_mts(3,ix,iy,ptype,iproc) ! jy
              data_pavg_tmp(j3_avg,ix,iy,ptype)= data_pavg_tmp(j3_avg,ix,iy,ptype) + &
                   p_mts(4,ix,iy,ptype,iproc) ! jz
              data_pavg_tmp(Tp_avg,ix,iy,ptype)= data_pavg_tmp(Tp_avg,ix,iy,ptype) + &
                   p_mts(1,ix,iy,ptype,iproc) ! n*<v²>
              if(plt_src.eq.0) then
                 data_pavg_tmp(sour_avg,ix,iy,ptype)= data_pavg_tmp(sour_avg,ix,iy,ptype) + &
                      ss2D(1,ix_f,iy_f,ptype,iproc)*k
                 data_pavg_tmp(sink_avg,ix,iy,ptype)= data_pavg_tmp(sink_avg,ix,iy,ptype) + &
                      ss2D(2,ix_f,iy_f,ptype,iproc)*k
              else
                 data_pavg_tmp(sour_avg,ix,iy,ptype)= data_pavg_tmp(sour_avg,ix,iy,ptype) + &
                      ss2D(1,ix_f,iy_f,ptype,iproc)
                 data_pavg_tmp(sink_avg,ix,iy,ptype)= data_pavg_tmp(sink_avg,ix,iy,ptype) + &
                      ss2D(2,ix_f,iy_f,ptype,iproc)
              endif
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

     enddo
  enddo

  if(nproc_mpi.gt.1) then
     allocate( data_pavg_MPI(7,0:na(1)+2,0:na(2)+2,ntype) )
     data_pavg_MPI=0.d0
     call MPI_ALLREDUCE(data_pavg_tmp(:,:,:,:), data_pavg_MPI, 7*(na(1)+3)*(na(2)+3)*ntype, &
          MPI_REAL8, MPI_SUM,MPI_COMM_WORLD, ierr)
     data_pavg_tmp= data_pavg_MPI
     deallocate( data_pavg_MPI )     
  endif

  !
  ! Calculate Tp
  !
  do ptype=1,ntype
        
     !$OMP PARALLEL PRIVATE(np_tmp,jx_tmp,jy_tmp,jz_tmp,Tp_tmp)
     !$OMP DO
     do iy=1,na(2)+1
        do ix=1,na(1)+1
           np_tmp= data_pavg_tmp(np_avg,ix,iy,ptype)
           data_pavg(np_avg,ix,iy,ptype)= data_pavg(np_avg,ix,iy,ptype) + np_tmp
           jx_tmp= data_pavg_tmp(j1_avg,ix,iy,ptype)
           data_pavg(j1_avg,ix,iy,ptype)= data_pavg(j1_avg,ix,iy,ptype) + jx_tmp
           jy_tmp= data_pavg_tmp(j2_avg,ix,iy,ptype)
           data_pavg(j2_avg,ix,iy,ptype)= data_pavg(j2_avg,ix,iy,ptype) + jy_tmp
           jz_tmp= data_pavg_tmp(j3_avg,ix,iy,ptype)
           data_pavg(j3_avg,ix,iy,ptype)= data_pavg(j3_avg,ix,iy,ptype) + jz_tmp
           if(np_tmp.gt.0) then
              Tp_tmp= mass(ptype)/qe*(1.d0/3.d0)*( data_pavg_tmp(Tp_avg,ix,iy,ptype)/np_tmp - &
                   (jx_tmp*jx_tmp + jy_tmp*jy_tmp + jz_tmp*jz_tmp)/(np_tmp*np_tmp) )
           else
              Tp_tmp= 0.d0
           endif
           ! Issue which may occur in corners or lack of statistics
           if(Tp_tmp.gt.0) data_pavg(Tp_avg,ix,iy,ptype)= data_pavg(Tp_avg,ix,iy,ptype) + Tp_tmp
           data_pavg(sour_avg,ix,iy,ptype)= data_pavg(sour_avg,ix,iy,ptype) + &
                data_pavg_tmp(sour_avg,ix,iy,ptype)
           data_pavg(sink_avg,ix,iy,ptype)= data_pavg(sink_avg,ix,iy,ptype) + &
                data_pavg_tmp(sink_avg,ix,iy,ptype)
        enddo
     enddo
     !$OMP END DO NOWAIT
     !$OMP END PARALLEL

  enddo
  
  return
end subroutine calc_avg

FUNCTION MSTIMER()
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Mar/10
!     MOD AUTHOR:    G. Hagelaar, G. Fubiani
!     COMMENTS:
!     NOTE: 
!     --------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(10):: cdat,ctim,czon
  INTEGER:: time(8),t_new,t_old,MSTIMER
  SAVE t_old

  CALL date_and_time(cdat,ctim,czon,time)
  t_new=((time(5)*60+time(6))*60+time(7))*1000+time(8)
  IF (t_old.EQ.0.d0) t_old=t_new
  MSTIMER=t_new-t_old
  t_old=t_new

  ! This happens at midnight ...
  if(MSTIMER.lt.0) MSTIMER=0

  RETURN

END FUNCTION MSTIMER
