program press
  implicit none
  integer ix,iy,n,nx,ny,ne,nH2p,nHm,nH3p, &
       nHp,every,iy0,nx_new,opt,cnt_grd
  parameter (n=2100,ne=1,nHm=2,nH2p=3,nH3p=4,nHp=5)
  real(kind=8):: dx,dy,ni(n,n,5),np_nm(n,n),dP(n,n,2),Ei(n,n,2), &
       sum_np(n,n),sum_nm(n),Te(n,n),phi(n,n),nm_ne(n,n),Bz(n,n),&
       vD(n,n,2),vE(n,n,2),jx(n,n,5),jy(n,n,5),ue(n,n,2),ui(n,n,2),&
       amu,qe,mi,ui_RMS,vE_RMS,ue_RMS

  open(8,file='input.dat',status='OLD')
  read(8,*) opt
  read(8,*) dx,dy
  read(8,*) nx_new
  read(8,*) iy0
  close(8)

  dx=dx*1.d-3
  dy=dy*1.d-3

  ni=0.d0
  np_nm=0.d0
  sum_np=0.d0
  sum_nm=0.d0
  dP=0.d0
  Ei=0.d0
  Bz=0.d0
  nm_ne=0.d0
  ue=0.d0
  ui=0.d0
  amu= 1.66d-27
  qe= 1.6d-19
  mi= amu
  ui_RMS=0.d0
  vE_RMS=0.d0
  ue_RMS=0.d0
  cnt_grd= 0

  if(opt.ge.2) open(7,file='../../DATA/Bz.dat',status='OLD',form='UNFORMATTED') ! Bz
  open(8,file='../../DATA/phi.dat',status='OLD',form='UNFORMATTED') ! phi
  open(9,file='../../DATA/T1.dat',status='OLD',form='UNFORMATTED') ! Te
  open(10,file='../../DATA/n1.dat',status='OLD',form='UNFORMATTED') ! ne
  open(11,file='../../DATA/n2.dat',status='OLD',form='UNFORMATTED') ! nH2+
  open(30,file='../../DATA/j1x.dat',status='OLD',form='UNFORMATTED') ! je_x
  open(31,file='../../DATA/j1y.dat',status='OLD',form='UNFORMATTED') ! je_y
  open(32,file='../../DATA/j2x.dat',status='OLD',form='UNFORMATTED') ! jH2+_x
  open(33,file='../../DATA/j2y.dat',status='OLD',form='UNFORMATTED') ! jH2+_y
  if(opt.eq.1 .or. opt.eq.2) then
     open(12,file='../../DATA/n3.dat',status='OLD',form='UNFORMATTED') ! nH-
     open(13,file='../../DATA/n4.dat',status='OLD',form='UNFORMATTED') ! nH3+
     open(14,file='../../DATA/n5.dat',status='OLD',form='UNFORMATTED') ! nH+
     open(34,file='../../DATA/j3x.dat',status='OLD',form='UNFORMATTED') ! jH-_x
     open(35,file='../../DATA/j3y.dat',status='OLD',form='UNFORMATTED') ! jH-_y
     open(36,file='../../DATA/j4x.dat',status='OLD',form='UNFORMATTED') ! jH3+_x
     open(37,file='../../DATA/j4y.dat',status='OLD',form='UNFORMATTED') ! jH3+_y
     open(38,file='../../DATA/j5x.dat',status='OLD',form='UNFORMATTED') ! jH+_x
     open(39,file='../../DATA/j5y.dat',status='OLD',form='UNFORMATTED') ! jH+_y
  endif
  open(15,file='n2_n1.mco')
  open(16,file='np.mco')
  if(opt.ge.2) then
     open(17,file='vDe.mco')
     open(18,file='vE.mco')
  endif
  open(19,file='ue.mco')
  open(20,file='ui.mco')
  if(iy0.gt.0) open(21,file='rho_1D.dat')
  if(opt.eq.1 .or. opt.eq.2) open(22,file='nHm_ne.mco')
  open(23,file='vDex.mco')
  open(24,file='vDey.mco')
  open(25,file='vEx.mco')
  open(26,file='vEy.mco')
  open(27,file='uex.mco')
  open(28,file='uey.mco')
  open(29,file='uix.mco')
  open(50,file='uiy.mco')
  open(51,file='uBp.mco')

  !
  ! Read file
  !
  if(opt.gt.1) then
     read(7) nx,ny
     read(7) Bz(1:nx+1,1:ny+1)
  endif
  read(8) nx,ny
  read(8) phi(1:nx+1,1:ny+1)
  read(9) nx,ny
  read(9) Te(1:nx+1,1:ny+1)
  read(10) nx,ny
  read(10) ni(1:nx+1,1:ny+1,ne)
  read(11) nx,ny
  read(11) ni(1:nx+1,1:ny+1,nH2p)
  read(30) nx,ny
  read(30) jx(1:nx+1,1:ny+1,ne)
  read(31) nx,ny
  read(31) jy(1:nx+1,1:ny+1,ne)
  read(32) nx,ny
  read(32) jx(1:nx+1,1:ny+1,nH2p)
  read(33) nx,ny
  read(33) jy(1:nx+1,1:ny+1,nH2p)
  if(opt.eq.1 .or. opt.eq.2) then
     read(12) nx,ny
     read(12) ni(1:nx+1,1:ny+1,nHm)
     read(13) nx,ny
     read(13) ni(1:nx+1,1:ny+1,nH3p)
     read(14) nx,ny
     read(14) ni(1:nx+1,1:ny+1,nHp)
     read(34) nx,ny
     read(34) jx(1:nx+1,1:ny+1,nHm)
     read(35) nx,ny
     read(35) jy(1:nx+1,1:ny+1,nHm)
     read(36) nx,ny
     read(36) jx(1:nx+1,1:ny+1,nH3p)
     read(37) nx,ny
     read(37) jy(1:nx+1,1:ny+1,nH3p)
     read(38) nx,ny
     read(38) jx(1:nx+1,1:ny+1,nHp)
     read(39) nx,ny
     read(39) jy(1:nx+1,1:ny+1,nHp)
  endif

  every= NINT(real(nx)/real(nx_new))
  print'(1x,"nx=",i5,", ny=",i5,", every=",i5)', nx,ny,every

  Bz= Bz*1.d-4 ! In tesla

  do iy=2,ny
     do ix=2,nx
        
        ! Electric field
        Ei(ix,iy,1)= -( phi(ix+1,iy) - phi(ix-1,iy) )/(2.d0*dx)
        Ei(ix,iy,2)= -( phi(ix,iy+1) - phi(ix,iy-1) )/(2.d0*dy)

        ! Bz
        if(opt.ge.2) then 
           
           if(ni(ix,iy,ne).gt.0.d0) then 
              ! Pressure gradient (Flux_i= -( diP_n + Ei )
              dP(ix,iy,1)= ( ni(ix+1,iy,ne)*Te(ix+1,iy) - &
                   ni(ix-1,iy,ne)*Te(ix-1,iy) )/(2.d0*dx)
              dP(ix,iy,1)= dP(ix,iy,1)/ni(ix,iy,ne)
              dP(ix,iy,2)= ( ni(ix,iy+1,ne)*Te(ix,iy+1) - &
                   ni(ix,iy-1,ne)*Te(ix,iy-1) )/(2.d0*dy)
              dP(ix,iy,2)= dP(ix,iy,2)/ni(ix,iy,ne)
           endif
           
           if( Bz(ix,iy).gt.0.d0) then
              ! vDe= -(nabla PxB)/(q*n*B^2)
              vD(ix,iy,1)= dP(ix,iy,2)/Bz(ix,iy)
              vD(ix,iy,2)= -dP(ix,iy,1)/Bz(ix,iy)
              ! vE= ExB/B^2
              vE(ix,iy,1)= Ei(ix,iy,2)/Bz(ix,iy)
              vE(ix,iy,2)= -Ei(ix,iy,1)/Bz(ix,iy)
           endif
           
        endif

        ! Density ratio
        if(opt.ge.1 .and. SUM(ni(ix,iy,ne:nHm)).gt.0.d0) then 
           np_nm(ix,iy)= SUM(ni(ix,iy,nH2p:nHp))/SUM(ni(ix,iy,ne:nHm)) - 1.d0
           nm_ne(ix,iy)= ni(ix,iy,nHm)/ni(ix,iy,ne)
        endif
        if(opt.eq.0 .and. ni(ix,iy,ne).gt.0.d0) &
             np_nm(ix,iy)= ni(ix,iy,nH2p)/ni(ix,iy,ne) - 1.d0
           
        ! Sum n+
        sum_np(ix,iy)= SUM(ni(ix,iy,nH2p:nHp))           

        if(ni(ix,iy,ne).gt.0.d0) then
           ue(ix,iy,1)= jx(ix,iy,ne)/ni(ix,iy,ne)
           ue(ix,iy,2)= jy(ix,iy,ne)/ni(ix,iy,ne)
        endif

        if(sum_np(ix,iy).gt.0.d0) then
           ui(ix,iy,1)= SUM(jx(ix,iy,nH2p:nHp))/sum_np(ix,iy)
           ui(ix,iy,2)= SUM(jy(ix,iy,nH2p:nHp))/sum_np(ix,iy)
        endif

     enddo
     if(iy0.gt.0) sum_nm(ix)= SUM(ni(ix,iy0,ne:nHm))
  enddo
  
     !
     ! Print in Macho format
     !
     write(15,*) nx/every,ny/every
     write(16,*) nx/every,ny/every
     if(opt.ge.2) then
        write(17,*) nx/every,ny/every
        write(18,*) nx/every,ny/every
     endif
     write(19,*) nx/every,ny/every
     write(20,*) nx/every,ny/every
     if(opt.eq.1 .or. opt.eq.2) &
          write(22,*) nx/every,ny/every
     do iy=ny+1,1,-1*every
        write(15,*) ( np_nm(ix,iy), ix=1,nx+1,every )  
        write(16,*) ( sum_np(ix,iy), ix=1,nx+1,every ) 
        if(opt.ge.2) then
           write(17,*) ( dsqrt( vD(ix,iy,1)**2. + &
                vD(ix,iy,2)**2. ), ix=1,nx+1,every ) 
           write(23,*) ( vD(ix,iy,1), ix=1,nx+1,every ) 
           write(24,*) ( vD(ix,iy,2), ix=1,nx+1,every )   
           write(18,*) ( dsqrt( vE(ix,iy,1)**2. + &
                vE(ix,iy,2)**2. ), ix=1,nx+1,every ) 
           write(25,*) ( vE(ix,iy,1), ix=1,nx+1,every ) 
           write(26,*) ( vE(ix,iy,2), ix=1,nx+1,every )  
        endif
        write(19,*) ( dsqrt( ue(ix,iy,1)**2. + &
             ue(ix,iy,2)**2. ), ix=1,nx+1,every ) 
        write(27,*) ( ue(ix,iy,1), ix=1,nx+1,every )  
        write(28,*) ( ue(ix,iy,2), ix=1,nx+1,every )  
        write(20,*) ( dsqrt( ui(ix,iy,1)**2. + &
             ui(ix,iy,2)**2. ), ix=1,nx+1,every )
        write(29,*) ( ui(ix,iy,1), ix=1,nx+1,every )  
        write(50,*) ( ui(ix,iy,2), ix=1,nx+1,every )  
        if(opt.eq.1 .or. opt.eq.2) &
             write(22,*) ( nm_ne(ix,iy), ix=1,nx+1,every ) 
        write(51,*) ( dsqrt(qe*Te(ix,iy)/mi) , ix=1,nx+1,every )  
     enddo

     if(opt.ge.2) then
        write(17,*) 'vector'
        do iy=ny+1,1,-1*every
           write(17,*) (datan2(vD(ix,iy,2),vD(ix,iy,1)), ix=1,nx+1,every)  
        enddo
        write(18,*) 'vector'
        do iy=ny+1,1,-1*every
           write(18,*) (datan2(vE(ix,iy,2),vE(ix,iy,1)), ix=1,nx+1,every)  
        enddo
     endif
     write(19,*) 'vector'
     do iy=ny+1,1,-1*every
        write(19,*) (datan2(ue(ix,iy,2),ue(ix,iy,1)), ix=1,nx+1,every)  
     enddo
     write(20,*) 'vector'
     do iy=ny+1,1,-1*every
        write(20,*) (datan2(ui(ix,iy,2),ui(ix,iy,1)), ix=1,nx+1,every)  
     enddo
     
     if(iy0.gt.0) then
        do ix=1,nx+1
           write(21,*) ix,sum_np(ix,iy0),sum_nm(ix)
        enddo
     endif

     do ix=1,nx+1
        do iy=1,ny+1
           if(opt.ge.2) then
              vE_RMS= vE_RMS + vE(ix,iy,1)**2 + vE(ix,iy,2)**2
           endif

           if(sum_np(ix,iy).gt.0.d0) cnt_grd= cnt_grd + 1
           ui_RMS= ui_RMS + ui(ix,iy,1)**2 + ui(ix,iy,2)**2
           ue_RMS= ue_RMS + ue(ix,iy,1)**2 + ue(ix,iy,2)**2
        enddo
     enddo

     if(opt.ge.2) then 
        vE_RMS= dsqrt(vE_RMS/cnt_grd)
        print*, 'vE_RMS=',vE_RMS
     endif
     ui_RMS= dsqrt(ui_RMS/cnt_grd)
     ue_RMS= dsqrt(ue_RMS/cnt_grd)
     print*, 'ui_RMS=',ui_RMS,', ue_RMS=',ue_RMS

     close(8)
     close(9)
     close(10)
     close(11)
     close(12)
     close(13)
     close(14)
     close(15)
     close(16)
     if(opt.ge.2) then
        close(17)
        close(18)
     endif
     close(19)
     close(20)
  if(iy0.gt.0) close(21)
  if(opt.eq.1 .or. opt.eq.2) close(22)
  close(23)
  close(24)
  close(25)
  close(26)
  close(27)
  close(28)
  close(29)
  close(50)
  close(51)

end program press
