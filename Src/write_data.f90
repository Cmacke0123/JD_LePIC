subroutine write_data(it,time,n,na,h,p_mac,P_loss,phi_avg,data_pavg,ntype,&
     ngrid,cnt_col,ncol_mx,sig_list,nproc,cnt_avg,mpi_rank,nproc_mpi,Pext,Vgrd,&
     I_inj,n_B)
!     ==============================================================
!     VERSION:         0.5
!     LAST MOD:      Dec/23
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:
!     NOTE: 
!     --------------------------------------------------------------
  use omp_lib
  implicit none
  include 'mpif.h'
  include 'particle_info.h'
  include 'constants.h'
  integer ierr,flag_stop
  integer:: ix,iy,ptype,n(3),na(3),nproc,ntype,i_rg,cnt_avg(2),n_B(2)
  integer:: np_avg,Tp_avg,sour_avg,iproc,j1_avg,j2_avg,j3_avg,ngrid,igrid,&
       icol,ncol_mx,it,sig_list(npart,ncol_mx),mpi_rank,nproc_mpi,every,i_pl,&
       sink_avg
  parameter ( np_avg=1, Tp_avg=2, j1_avg=3, j2_avg=4, j3_avg=5, sour_avg=6,&
       sink_avg=7 ) 
  real(kind=8):: h(3),p_mac(ntype,2,0:ngrid,nproc), &
       phi_avg(3,0:n(1)+2,0:n(2)+2), &
       data_pavg(7,na(1)+1,na(2)+1,ntype),&
       P_loss(4,ntype,nproc),time,cnt_col(3,ncol_mx,nproc,nm_rg), &
       sum_nu(ntype,nm_rg),P_loss_tmp(4,ntype),cnt_col_tmp(3,ncol_mx,nm_rg),&
       p_mac_tmp(ntype,2,0:ngrid),dh,Ip,Im,Pwall,Pext,Pinj,Pcoll,I_mw,&
       Ek_ew,Ek_iw,Vgrd(ngrid),Ptmp(ntype),I_inj,Isec(ngrid),Ip_cath
  real(kind=8),allocatable:: ld(:,:),cnt_col_red(:,:,:),p_mac_red(:,:,:),&
       P_loss_red(:,:)       
  character, save:: name(70)*15,pnum*1,name0D(16)*15,strg*3,plnum*3,corrnum*3

  allocate( ld(na(1)+1,na(2)+1) )
  
  ! Initialize
  flag_stop=0
  ld=0.d0

  ! 1D case
  if(flag_1D.eq.1) data_pavg(:,:,2,:)= data_pavg(:,:,1,:)
  
  !
  ! Calculate average collision frequency  
  !
  cnt_col_tmp= 0

  ! Reduction 
  do i_rg=1,n_rg
     do iproc=1,nproc
        do icol=1,ncol
           ! Average frequency calculated over all collision candidates
           cnt_col_tmp(1:2,icol,i_rg)= cnt_col_tmp(1:2,icol,i_rg) + &
                cnt_col(1:2,icol,iproc,i_rg)
           ! Actual collision events
           cnt_col_tmp(3,icol,i_rg)= cnt_col_tmp(3,icol,i_rg) + &
                cnt_col(3,icol,iproc,i_rg)/(real(it)*dt)
        enddo
     enddo
  enddo

  if(nproc_mpi.gt.1) then
     allocate( cnt_col_red(3,ncol_mx,nm_rg) )
     cnt_col_red=0.d0
     call MPI_REDUCE(cnt_col_tmp(:,:,:), cnt_col_red, 3*ncol_mx*nm_rg, &
          MPI_REAL8, MPI_SUM, 0,MPI_COMM_WORLD, ierr)
     cnt_col_tmp= cnt_col_red
     deallocate( cnt_col_red )
  endif

  !
  ! Calculate average power & current loss 
  !

  P_loss_tmp= 0.d0
  p_mac_tmp= 0.d0

  ! Reduction
  do iproc=1,nproc
     do ptype=1,ntype
        P_loss_tmp(:,ptype)= P_loss_tmp(:,ptype) + &
             P_loss(:,ptype,iproc)/real(cnt_avg(2))/dt
        p_mac_tmp(ptype,np_loss,:)= p_mac_tmp(ptype,np_loss,:) + &
             p_mac(ptype,np_loss,:,iproc)*charge(ptype)*Nm(ptype)/real(cnt_avg(2))/dt
        p_mac_tmp(ptype,P_w,:)= p_mac_tmp(ptype,P_w,:) + &
             p_mac(ptype,P_w,:,iproc)/real(cnt_avg(2))/dt
     enddo
  enddo

  if(nproc_mpi.gt.1) then
     allocate( p_mac_red(ntype,2,0:ngrid) )
     p_mac_red=0.d0
     call MPI_REDUCE(p_mac_tmp(:,:,:), p_mac_red, ntype*2*(ngrid+1), &
          MPI_REAL8, MPI_SUM, 0,MPI_COMM_WORLD, ierr) ! mpi_rank=0 only
     p_mac_tmp= p_mac_red
     deallocate( p_mac_red )

     allocate( P_loss_red(4,ntype) )
     P_loss_red=0.d0
     call MPI_REDUCE(P_loss_tmp(:,:), P_loss_red, 4*ntype, &
          MPI_REAL8, MPI_SUM, 0,MPI_COMM_WORLD, ierr)
     P_loss_tmp= P_loss_red
     deallocate( P_loss_red )
  endif

  Ip= SUM(p_mac_tmp(2:ntype,np_loss,0:ngrid))
  Ip_cath= SUM(p_mac_tmp(2:ntype,np_loss,igrid_sec))
  if(tag_neg.gt.0) then
     Ip= Ip - SUM(p_mac_tmp(tag_neg,np_loss,0:ngrid))
     Ip_cath= Ip_cath - p_mac_tmp(tag_neg,np_loss,igrid_sec)
  endif
  ! Isec >0 as it corresponds to electrons drawn from the power supply
  Isec=0.d0
  ! Fixed secondary electron emission coefficient gamma
  if(flag_sec.eq.1) &
       Isec(igrid_sec)= gam_sec*Ip_cath

  ! Inject a fixed current of secondary electrons with an adaptative gamma
  if(flag_sec.eq.2) then
     ! Send to all MPI threads
     if(nproc_mpi.gt.1) call MPI_Bcast(Ip_cath,1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     ! Warning
     if(Ip_cath.le.0) then
        if(mpi_rank.eq.0) then
           print*, 'Warning: Ip_cath <=0, please correct...'
           flag_stop=1
        endif
     endif
     Isec(igrid_sec)= ABS(I_inj)
     ! Update secondary emission yield
     if(it.gt.1) gam_sec= Isec(igrid_sec)/Ip_cath
  endif

  igrid= ind_g
  ! Inject a uniform particle current on the electrode surface
  if(flag_grd(1).eq.0) igrid= igrid_sec
  Isec(igrid)= Isec(igrid) + jne*10.d0*zmax*Lgrd

  ! Add secondary electron current emitted on cathode
  p_mac_tmp(1,np_loss,igrid_sec)= p_mac_tmp(1,np_loss,igrid_sec) + Isec(igrid_sec)
  if(flag_grd(1).gt.0) p_mac_tmp(1,np_loss,ind_g)= p_mac_tmp(1,np_loss,ind_g) + Isec(ind_g)

  !
  ! Only on rank 0  
  !
  if(mpi_rank.eq.0) then

     dh= dsqrt(h(1)*h(2))
     if(flag_1D.eq.1) dh= h(1) ! 1D

     ! Calculate Temperature & Debye length
     do ptype=1,ntype
        !$OMP PARALLEL
        !$OMP DO
        do iy=1,na(2)+1
           do ix=1,na(1)+1
              if(data_pavg(np_avg,ix,iy,ptype).gt.0.d0) then
                 ! Electron Debye length
                 if(ptype.eq.1) then 
                    ld(ix,iy)= dsqrt( (eps0*data_pavg(Tp_avg,ix,iy,1))/ &
                         (data_pavg(np_avg,ix,iy,1)*qe) )                                        
                    ! Calculate ratio with dr
                    if(ld(ix,iy).gt.0) ld(ix,iy)=dh/ld(ix,iy)
                 endif
              endif
           enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL       
     enddo

     if(flag_pbc.eq.1) then 
        data_pavg(sour_avg,:,na(2)+1,:)= data_pavg(sour_avg,:,1,:)
        data_pavg(sink_avg,:,na(2)+1,:)= data_pavg(sink_avg,:,1,:)
     endif

     !
     ! Write data
     !
     Im= 0.d0
     if(tag_neg.gt.0) Im= SUM(p_mac_tmp(tag_neg,np_loss,0:ngrid))
     do ptype=1,ntype ! Include contribution of the bias voltage on the total power
        Ptmp(ptype)= SUM(p_mac_tmp(ptype,np_loss,1:ngrid)*Vgrd(1:ngrid))
     enddo
     Pwall= -SUM(P_loss_tmp(1,1:ntype)) -SUM(Ptmp(1:ntype)) ! -sum(I)*Vbias
     Pext= SUM(P_loss_tmp(2,1:ntype))
     Pcoll= SUM(P_loss_tmp(3,1:ntype))
     Pinj= SUM(P_loss_tmp(4,1:ntype))
     I_mw= SUM(p_mac_tmp(1,np_loss,0:ngrid))+Im
     if(SUM(p_mac_tmp(1,np_loss,0:ngrid)).ne.0.d0) &
          Ek_ew= SUM(p_mac_tmp(1,P_w,0:ngrid))/ABS(SUM(p_mac_tmp(1,np_loss,0:ngrid)))
     Ek_iw= 0.d0
     if(SUM(ABS(p_mac_tmp(2:ntype,np_loss,0:ngrid))).gt.0.d0) &
          Ek_iw= SUM(p_mac_tmp(2:ntype,P_w,0:ngrid))/SUM(ABS(p_mac_tmp(2:ntype,np_loss,0:ngrid))) 
     write(*,102) Pwall,Pext,Pcoll,Pinj,I_mw,Ip,Ek_ew,Ek_iw  ! print on screen
     
     open(22,file='DATA/Ptot.dat',access='APPEND')
     if( it.eq.1 .and. flag_restart.eq.0 ) write(22,'("# Time (s), Pwall (W), Pabs, Pcoll, Pinj, I_mw (A), I_pw, Ek_ew (eV), Ek_iw, I_inj(A), Vgrd")')
     write(22,101) time,Pwall,Pext,Pcoll,Pinj,I_mw,Ip,Ek_ew,Ek_iw,I_inj,Vgrd(igrid_sec) 
     close(22)
     
     ! Collision frequencies per specie per sub-regions
     if(it.gt.1) then
        do i_rg=1,n_rg
           
           write (pnum,'(i1)'),i_rg
           do icol=1,ncol
              if(cnt_col_tmp(2,icol,i_rg).gt.0) then
                 cnt_col_tmp(1,icol,i_rg)= cnt_col_tmp(1,icol,i_rg)/cnt_col_tmp(2,icol,i_rg)
              else
                 cnt_col_tmp(1,icol,i_rg)= 0.d0
              endif
           enddo
           
           do ptype=1,ntype
              sum_nu(ptype,i_rg)= SUM(cnt_col_tmp(3,sig_list(ptype,1:p_ncol(ptype)),i_rg))
              
              ! Warning
              if(sum_nu(ptype,i_rg).gt.0.25d0*nu_uplim(ptype)) then
                 print*, 'Warning: sum_nu/nu_uplim=',sum_nu(ptype,i_rg)/nu_uplim(ptype),'>25% for ptype=',ptype
                 print*, 'nu_max=',nu_max(ptype),', nu_uplim=',nu_uplim(ptype)
                 flag_stop=1
              endif
              
              ! Warning
              if( (sum_nu(ptype,i_rg)*real(ns_coll)*dt).gt.0.2d0 ) then
                 print*, 'Warning: sum_nu*dt > 20% for ptype=',ptype
              endif
              
           enddo
           
           open(14,file='DATA/frequency'//pnum//'.out')
           
           write(14,*) '# null collision frequency for each species'
           write(14,'(10(1x,es15.8))') ( nu_max(ptype),ptype=1,ntype )
           write(14,*) '# sum of collision frequencies for each species'
           write(14,'(10(1x,es15.8))') ( sum_nu(ptype,i_rg),ptype=1,ntype )
           write(14,*) '# reaction index/frequencies for each species'
           
           do icol=1,ncol
              do ptype=1,ntype
                 strg='no'
                 if(ptype.eq.ntype) strg='yes'
                 if(icol.le.p_ncol(ptype)) then
                    write(14,'((i3,1x,es15.8))',advance=strg) sig_list(ptype,icol), &
                         cnt_col_tmp(1,sig_list(ptype,icol),i_rg)
                 else
                    write(14,'((i3,1x,es15.8))',advance=strg) 0,0.
                 endif
              enddo
           enddo
           
           close(14)           
        enddo
     endif

  endif

  ! Broadcast Pext to other MPI threads
  if(nproc_mpi.gt.1) call MPI_Bcast(Pext, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

  ! Save potential and Debye length
  if(mpi_rank.eq.0) then

     ! Get data file names
     if(it.eq.1) then
        open(14,file='Macho/namelist.inp')
        if(Bmax.gt.0.d0 .and. n_B(1).gt.1) then 
           write(14,*) 'B.dat'
           write(14,*) 'Bz.dat'
        endif
        write(14,*) 'bcnd.dat'
        write(14,*) 'kq.dat'
        write(14,*) 'phi.dat'
        write(14,*) 'Ex.dat'
        write(14,*) 'Ey.dat'
        write(14,*) 'dr.dat'
     endif
        
     do ptype=1,ntype
        write (pnum,'(i1)'),ptype
           
        name(7*ptype-6)='n'//pnum//'.dat'
        name(7*ptype-5)='T'//pnum//'.dat'
        name(7*ptype-4)='sour'//pnum//'.dat'
        name(7*ptype-3)='j'//pnum//'x.dat'
        name(7*ptype-2)='j'//pnum//'y.dat'
        name(7*ptype-1)='j'//pnum//'z.dat'
        name(7*ptype)='sink'//pnum//'.dat'

        if(it.eq.1) then
           write(14,*) 'n'//pnum//'.dat'
           write(14,*) 'T'//pnum//'.dat'
           write(14,*) 'sour'//pnum//'.dat'
           write(14,*) 'sink'//pnum//'.dat'
           if(flag_1D.eq.0) then
              write(14,*) 'j'//pnum//'.dat'
           else
              write(14,*) 'j'//pnum//'x.dat'
           endif

           name0D(2*ptype-1)='Iw'//pnum//'.dat' 
           name0D(2*ptype)='Pw'//pnum//'.dat'  
        endif
     enddo        

     if(it.eq.1) close(14)
        
     open(14,file='DATA/dr.dat',form='UNFORMATTED')
     open(15,file='DATA/phi.dat',form='UNFORMATTED')
     open(16,file='DATA/Ex.dat',form='UNFORMATTED')
     open(17,file='DATA/Ey.dat',form='UNFORMATTED')
    
     write(14) na(1),na(2)
     write(15) n(1),n(2)
     write(16) n(1),n(2)
     write(17) n(1),n(2)
     write(14) ld(1:na(1)+1,1:na(2)+1)
     write(15) phi_avg(1,1:n(1)+1,1:n(2)+1)/real(cnt_avg(1))
     write(16) phi_avg(2,1:n(1)+1,1:n(2)+1)/real(cnt_avg(1))
     write(17) phi_avg(3,1:n(1)+1,1:n(2)+1)/real(cnt_avg(1))
     close(14)
     close(15)
     close(16)
     close(17)
     
     ! Write density, temperature, source term and currents for each species
     do ptype=1,ntype

        open(14,file='DATA/'//name(7*ptype-6),form='UNFORMATTED')
        open(15,file='DATA/'//name(7*ptype-5),form='UNFORMATTED')
        open(16,file='DATA/'//name(7*ptype-4),form='UNFORMATTED')
        open(17,file='DATA/'//name(7*ptype-3),form='UNFORMATTED')
        open(18,file='DATA/'//name(7*ptype-2),form='UNFORMATTED')
        open(19,file='DATA/'//name(7*ptype-1),form='UNFORMATTED')
        open(20,file='DATA/'//name(7*ptype),form='UNFORMATTED')

        write(14) na(1),na(2)
        write(15) na(1),na(2)
        write(16) na(1),na(2)
        write(17) na(1),na(2)
        write(18) na(1),na(2)
        write(19) na(1),na(2)
        write(20) na(1),na(2)
        write(14) data_pavg(np_avg,1:na(1)+1,1:na(2)+1,ptype)/real(cnt_avg(1))
        write(15) data_pavg(Tp_avg,1:na(1)+1,1:na(2)+1,ptype)/real(cnt_avg(1))
        if(plt_src.eq.0) then
           write(16) data_pavg(sour_avg,1:na(1)+1,1:na(2)+1,ptype)/real(cnt_avg(2)*dt)
        else
           write(16) data_pavg(sour_avg,1:na(1)+1,1:na(2)+1,ptype)/real(cnt_avg(1))
        endif
        write(17) data_pavg(j1_avg,1:na(1)+1,1:na(2)+1,ptype)/real(cnt_avg(1))
        write(18) data_pavg(j2_avg,1:na(1)+1,1:na(2)+1,ptype)/real(cnt_avg(1))
        write(19) data_pavg(j3_avg,1:na(1)+1,1:na(2)+1,ptype)/real(cnt_avg(1))
        if(plt_src.eq.0) then
           write(20) data_pavg(sink_avg,1:na(1)+1,1:na(2)+1,ptype)/real(cnt_avg(2)*dt)
        else
           write(20) data_pavg(sink_avg,1:na(1)+1,1:na(2)+1,ptype)/real(cnt_avg(1))
        endif    
        close(14)
        close(15)
        close(16)
        close(17)
        close(18)
        close(19)
        close(20)
        
     enddo ! end loop over ptype
     
     ! Save profiles for a given time window without averaging
     if(flag_bak.eq.2) then 
        every= INT(na(1)/256)
        if(every.eq.0) every=1

        cnt_plt= cnt_plt+1
        if(cnt_plt.lt.10) then 
           write (plnum,'(i1)'),cnt_plt
           i_pl=1
        endif
        if( cnt_plt.ge.10 .and. cnt_plt.lt.100 ) then 
           write (plnum,'(i2)'),cnt_plt
           i_pl=2
        endif
        if( cnt_plt.ge.100 .and. cnt_plt.lt.1000 ) then 
           write (plnum,'(i3)'),cnt_plt
           i_pl=3
        endif
        
        do ptype=1,ntype
           
           if(i_pl.eq.1) corrnum= '_00'
           if(i_pl.eq.2) corrnum= '_0'
           if(i_pl.eq.3) corrnum= '_'
           
           ! Filenames
           write (pnum,'(i1)'),ptype
           name(4*ptype-3)='n'//pnum//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco'
           name(4*ptype-2)='T'//pnum//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco'                      
           name(4*ptype-1)='j'//pnum//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco'      
           name(4*ptype)='phi'//corrnum(1:3-i_pl+1)//plnum(1:i_pl)//'.mco'           
           
           open(15,file='Macho/DATA_seq/'//name(4*ptype-3))
           open(16,file='Macho/DATA_seq/'//name(4*ptype-2))
           open(17,file='Macho/DATA_seq/'//name(4*ptype-1))
           open(18,file='Macho/DATA_seq/'//name(4*ptype))
           
           write(15,*) na(1)/every,na(2)/every
           write(16,*) na(1)/every,na(2)/every
           write(17,*) na(1)/every,na(2)/every
           write(18,*) na(1)/every,na(2)/every
           do iy=na(2)+1,1,-1*every
              write(15,103) ( data_pavg(np_avg,ix,iy,ptype)/real(cnt_avg(1)), ix=1,na(1)+1,every )
              write(16,103) ( data_pavg(Tp_avg,ix,iy,ptype)/real(cnt_avg(1)), ix=1,na(1)+1,every )
              write(17,103) ( dsqrt( data_pavg(j1_avg,ix,iy,ptype)**2 + &
                   data_pavg(j2_avg,ix,iy,ptype)**2 ), ix=1,na(1)+1,every )
              write(18,103) ( phi_avg(1,ix,iy)/real(cnt_avg(1)), ix=1,na(1)+1,every )
           enddo
           
           write(17,*) 'vector'
           do iy=na(2)+1,1,-1*every ! Theta= arctg(jy/jx)
              write(17,103) ( datan2(data_pavg(j2_avg,ix,iy,ptype),&
                   data_pavg(j1_avg,ix,iy,ptype) ), ix=1,na(1)+1,every )  
           enddo
           
           close(15)
           close(16)
           close(17)
           close(18)
           
        enddo ! end loop over ptype
     endif

     !
     ! Write space integrated data
     !
     do ptype=1,ntype        
        open(22,file='DATA/'//name0D(2*ptype-1),access='APPEND')
        open(23,file='DATA/'//name0D(2*ptype),access='APPEND')
        
        write(22,101) time,( p_mac_tmp(ptype,np_loss,igrid), igrid=0,ngrid )
        write(23,101) time,( p_mac_tmp(ptype,P_w,igrid), igrid=0,ngrid )

        close(22)
        close(23)             
     enddo

     ! Save some local values of phi, ne and Te
     open(24,file='DATA/phi_Te_ne_cntr.dat',access='APPEND')
     if( it.eq.1 .and. flag_restart.eq.0 ) write(24,'("# Time(s), phi(xm/2,ym/2), Te(xm/2,ym/2), ne at ym/2 and &
          xm/2, xm/10, 2xm/10, 3xm/10, 4xm/10, 6xm/10, 7xm/10, 8xm/10, 9xm/10 ")')
     write(24,101) time, phi_avg(1,n(1)/2+1,n(2)/2+1)/real(cnt_avg(1)),&
          data_pavg(Tp_avg,na(1)/2+1,na(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg(np_avg,na(1)/2+1,na(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg(np_avg,NINT(na(1)/10.)+1,na(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg(np_avg,NINT(2*na(1)/10.)+1,na(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg(np_avg,NINT(3*na(1)/10.)+1,na(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg(np_avg,NINT(4*na(1)/10.)+1,na(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg(np_avg,NINT(6*na(1)/10.)+1,na(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg(np_avg,NINT(7*na(1)/10.)+1,na(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg(np_avg,NINT(8*na(1)/10.)+1,na(2)/2+1,1)/real(cnt_avg(1)),&
          data_pavg(np_avg,NINT(9*na(1)/10.)+1,na(2)/2+1,1)/real(cnt_avg(1))
     close(24)

  endif

  deallocate(ld)

101 format(20(1x,es16.8))
102 format(' Pwall (W)= ',es10.2,', Pabs (W)= ',es10.2,', Pcoll (W)= ',es10.2,&
         ', Pinj (W)= ',es10.2,', I_w (A)= ',2(es10.2,2x),', Ek_w (eV)= ',2(es10.2,2x))
103 format(800(e18.6,1x))

  if(nproc_mpi.gt.1) call MPI_Bcast(flag_stop,1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
  if(flag_stop.eq.1) call stop_calculation
     
  return

end subroutine write_data
