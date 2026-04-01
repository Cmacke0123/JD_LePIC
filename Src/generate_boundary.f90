subroutine generate_boundary(u,n,h,bcnd,V,ngrid,n_icp,mpi_rank)
!     ==============================================================
!     VERSION:         0.7
!     LAST MOD:       Mar/22
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Generates boundary condition matrix
!     NOTE:              
!     --------------------------------------------------------------
  implicit none
  integer:: i,j,k,flag,istart,iend,jstart,ind,ind_sav,ind_max, &
       ngrid,j0,j1,i0,i1,ih,jm,jp,cor,n_icp,mpi_rank,check_grd
  integer:: n(3),jmax(0:n(1)+2),bcnd(0:n(1)+2,0:n(2)+2), &
       i_tmp,j_tmp,flag_dc
  real(kind=8):: h(3),u(0:n(1)+2,0:n(2)+2),x1,x2,y1,y2,ytmp,V(ngrid), &
       xl,xr,yl,yr,Rc,yg_sec
  include 'particle_info.h'

  ind_max= 0
  flag_chmfrd=0
  check_grd=0

  ! Read potential values
  open(10,file='input_dir/boundary.inp')
  do i=1,ngrid
     read(10,*,end=1000) V(i)
  enddo
  
1000 continue
  close(10)


  ! Warning
  if( (i-1).lt.ngrid ) then
     if(mpi_rank.eq.0) then
        print*, 'Insufficient number of wall labels found in file boundary.inp '
        print*, 'please correct ...'
     endif
     call stop_calculation
  endif

  
  ! Open input files
  open(10,file='input_dir/geometry.inp')

  ! Define box size
  read(10,*) xl,xr,yr
  flag_dc=0
  if(yr.lt.0.d0) then 
     flag_dc=1
     yr= ABS(yr)
  endif
  yl=-yr

  xmax= (xr-xl)
  ymax= (yr-yl)

  ! Define dx,dy
  h(1)=(xr-xl)/n(1)
  h(2)=(yr-yl)/n(2)   

  ! Initialize simulation box
  bcnd=-1


  ! Draw a circle
  if(flag_dc.eq.1) then
     Rc= MIN(xmax/2.d0,ymax/2.d0)
     do j=1,n(2)+1
        do i=1,n(1)+1
           if( (((i-1)*h(1)-xmax/2.d0)**2 +  &
                ((j-1)*h(2)-ymax/2.d0)**2).ge.Rc**2 ) then 
              bcnd(i,j)= 1
              u(i,j)= V(1)
           endif
        enddo
     enddo
     bcnd(0,:)=bcnd(1,:)
     bcnd(n(1)+2,:)=bcnd(n(1)+1,:)
     bcnd(:,0)=bcnd(:,1)
     bcnd(:,n(2)+2)=bcnd(:,n(2)+1)
     goto 999
  endif

  dtype=0
  flag_pbc=0
  flag_spec=0
  flag_nmn=0
  flag_dielec=0
  flag=0
  k=1
  do while(flag.eq.0)

     if(k.eq.1) then
        read(10,*,end=999) x1,y1,ind
        ind_max= MAX(ind_max,ind)
        read(10,*,end=999) x2,y2,ind
        ind_max= MAX(ind_max,ind)
     else
        read(10,*,end=999) x2,y2,ind
        ind_max= MAX(ind_max,ind)
        if(ind.ne.ind_sav) then
           x1= x2
           y1= y2
           read(10,*,end=999) x2,y2,ind
           ind_max= MAX(ind_max,ind)
        endif
     endif

     ! Chamfered grid (front side only)
     if( x2.ne.x1 .and. y2.lt.y1 .and. flag_chmfrd.eq.0 ) then
        flag_chmfrd=1
        chmfrd_array(1)= x1*1.d-2 ! in meters
        chmfrd_array(2)= x2*1.d-2
        chmfrd_array(3)= (ymax/2.d0-y1)*1.d-2
        chmfrd_array(4)= (ymax/2.d0-y2)*1.d-2
     endif
     
     ! Secondary particle emission
     xg_sec= 0.d0
     if(flag_sec.ge.1 .and. igrid_sec.eq.ind) then
        i= nint(x1/h(1))+1
        xg_sec= (i-1)*h(1)*1.d-2 ! in meters
        yg_sec= (y1+y2)/2.d0*1.d-2
        Lg= 2.d0*abs((y2-y1))
        Lh= 0.d0
        ! Warning 
        if((x2-x1).ne.0) then
           if(mpi_rank.eq.0) then
              print*, 'Secondary particle production can only be implemented in x=cste planes'
              print*, 'please correct...'
           endif
           call stop_calculation
        endif
        check_grd=1
     endif

     ! Periodic BCs
     if( ind.eq.0 .or. ind.eq.-3 ) then 
        if(ind.eq.-3) flag_spec=1
        flag_pbc=1
        ind=0
     endif

     ! Neumann BCs
     if(ind.eq.-2) flag_nmn=1

     ! Dielectric BCs
     if(ind.le.-4) then 
        flag_dielec=1
        ind= ABS(ind)
        ind_max= MAX(ind_max,ind)
        if(x1.eq.x2) dtype(ind)=1
        if(y1.eq.y2) then           
           if(mpi_rank.eq.0) &
                print*, 'Dielectrics can only be implemented in x=cste planes, please correct...'
           call stop_calculation
        endif
     endif

     ! Warning
     if( ind_max.gt.ngrid ) then
        if(mpi_rank.eq.0) then
           print*, '# of grid labels > ngrid '
           print*, 'please correct ...'
        endif
        call stop_calculation
     endif

     ! Set istart, iend
     istart= nint((x1-xl)/h(1))+1
     iend= nint((x2-xl)/h(1))+1

     ! Correct for lines parallel to coordinate axis
     if((x2-x1).eq.0) then
        istart= nint((x1-xl)/h(1))+1
        iend= istart
        jstart= MIN(nint((y1-yl)/h(2))+1,nint((y2-yl)/h(2))+1) 
     endif

     if((y2-y1).eq.0) then
        istart= nint((x1-xl)/h(1))+1
        iend= nint((x2-xl)/h(1))+1
        jstart= nint((y1-yl)/h(2))+1
     endif

     ! Correct istart & iend if necessary
     if( ind.eq.0 ) istart= istart+1

     if( istart.eq.1 ) istart= 0
     if( iend.eq.(n(1)+1) ) iend= n(1)+2

     ! Start drawing
     do i=istart,iend
        
        ! Oblique
        if(x2.ne.x1.and.y2.ne.y1) then
           ytmp= (y2-y1)/(x2-x1)*((i-1)*h(1)-x1) + y1
           jstart= nint((ytmp-yl)/h(2))+1
        endif
        
        jmax(i)=jstart
        
        do j=jstart,n(2)+2

           ! Fill bcnd
           bcnd(i,j)= ind
           bcnd(i,2+n(2)-j)= ind

           ! Initialize u-vector
           if(ind.ge.1) then 
              u(i,j)= V(ind)
              u(i,2+n(2)-j)= V(ind)
           else
              u(i,j)= 0.d0
              u(i,2+n(2)-j)= 0.d0              
           endif

           ! Store wall coordinates
           if( (iend-istart).le.1 ) then
              i_tmp=istart
              j_tmp=j
           else
              i_tmp=i
              j_tmp=jstart
           endif
           if(i_tmp.eq.0) i_tmp=1
           if(j_tmp.eq.n(2)+2) j_tmp=n(2)+1
        enddo

     enddo

     ! Save boundary infos
     x1=x2
     y1=y2
     ind_sav=ind
     k=k+1

  enddo

999 continue
  close(10)

  !
  ! Draw grid and create holes
  !
  if(xg1.gt.xmax) xg1=xmax
  if(xg2.gt.xmax) xg2=xmax
  if(Lg.gt.ymax) Lg= ymax
  if(Lh.ge.Lg) then
     if(mpi_rank.eq.0) then
        print*, 'Warning Lh>Lg, please correct ...'
     endif
     call stop_calculation
  endif

  ! No grid case
  if(flag_grd(1).eq.0) goto 20 

  ! Lower edge of grid
  if( ((ymax-Lg)/2.).gt.0 ) then
     j0= INT( ((ymax-Lg)/2.)/h(2) ) + 1
  else
     j0= 1
  endif

  ! Half-size index of a single hole
  ih= INT(Lh/h(2))

  ! Update value of Lh 
  Lh= ih*h(2)
  if(mpi_rank.eq.0) print'(1x,"Updated hole radius (cm):",f8.2," cm")',Lh

  ! Lower edge of holes region
  y1= ymax/2. - (real(n_holes)+0.5)*Lh
  j1= INT(y1/h(2)) + 1
  
  if(flag_grd(2).eq.0) ind= ind_g

  ! Warning
  if( flag_grd(2).eq.1 .and. (ind+n_holes+1).gt.ngrid ) then 
     if(mpi_rank.eq.0) then
        print'(" Warning: ngrid must be set to :",1x,i2)',ind+n_holes+1
        print*, 'Please correct ...'
     endif
     call stop_calculation
  endif

  if(j0.gt.j1) then
     if(mpi_rank.eq.0) then
        print'(" Warning: grid size too small for # of holes:")'
        print*, 'Please correct ...'
     endif
     call stop_calculation
  endif

  i0= NINT(xg1/h(1))+1
  i0=MAX(i0,2)
  i1= NINT(xg2/h(1))+1
  i1=MIN(i1,n(1))

  if(xg1.eq.xmax) then
     i0= n(1)+1
     i1= i0+1
  endif

  ! Loop over hole number
  do k=0,n_holes,1

     if(flag_grd(2).eq.1) ind=ind+1

     ! Loop over grid width
     do i=i0,i1

        i_tmp=i0
        if(i.gt.NINT(real(i0+i1)/2.)) i_tmp=i1
  
        ! Draw grid
        if(k.eq.0) then
           if(n_holes.eq.0) then
              ! No holes then draw everything
              jm= j0
              jp= n(2)+2-j0
              do j= jm,jp
                 bcnd(i,j)= ind_g
                 u(i,j)= V(ind_g)
              enddo
              goto 10 ! jump to end of i-loop
           else ! Draw lower edge
              jm= j0
              jp= j1
              do j= jm,jp
                 bcnd(i,j)= ind
                 u(i,j)= V(ind)
              enddo
           endif
        endif
        
        if(k.eq.n_holes) then ! Draw upper edge
              jm= j1+(2*k+1)*ih
              jp= n(2)+2-j0
           do j= jm,jp
              bcnd(i,j)= ind
              u(i,j)= V(ind)
           enddo
        endif

        ! Interstice between holes
        jm= j1+(2*k*ih)
        jp= j1+(2*k+1)*ih
        cor=1
        if(ABS(jp-jm).ge.4) cor=0
        do j= jm,jp 
           bcnd(i,j)= ind
           u(i,j)= V(ind)
           ! Interior surface of a hole (2 grid depth)
           j_tmp= j
           if( k.gt.0 .and. i.gt.i0 .and. i.lt.i1 .and. &
                j.le.(jm+1-cor) ) j_tmp= jm
           if( k.lt.n_holes .and. i.gt.i0 .and. i.lt.i1 .and. &
                j.ge.(jp-2+cor) ) j_tmp= jp-1
        enddo

        ! Hole interiors
        if(k.lt.n_holes) then
           jm= j1+(2*k+1)*ih
           jp= j1+(2*k+2)*ih-1
           do j= jm,jp
              bcnd(i,j)=-1 
              ! Impose a hole potential
              if(flag_grd(1).eq.2) then
                 bcnd(i,j)= ind_g-1
                 u(i,j)= V(ind_g-1)
              endif     
           enddo
        endif
        
10      continue
     enddo
  enddo
  
  !
  ! Add multi-drivers to simulation domain
  !
20 continue
  if(flag_icp.eq.1) then
     do i=0,n(1)+2 ! Shrink
        do j=1,n(2)/n_icp+1
           bcnd(i,j)= bcnd(i,j*n_icp)
           u(i,j)= u(i,j*n_icp)
        enddo
        
        do k=2,n_icp
           do j=1,n(2)/n_icp+1 ! Add drivers
              bcnd(i,k*n(2)/n_icp-j)= bcnd(i,j)
              u(i,k*n(2)/n_icp-j)= u(i,j)
           enddo
        enddo
     enddo
  
     bcnd(:,1)= bcnd(:,n(2)+1)
     u(:,1)= u(:,n(2)+1)
  
     ymax=ymax*real(n_icp)
     h(2)=h(2)*real(n_icp)
  endif

  !
  ! 1D simulations
  !
  if(flag_1D.eq.1) then
     flag_pbc=0
     bcnd(2:n(1),:)=-1
  endif
  
  !
  ! Write boundaries in files
  !
30 continue
  open(12,file='DATA/bcnd.dat',form='UNFORMATTED')
  write(12) n(1),n(2)
  write(12) dble(bcnd(1:n(1)+1,1:n(2)+1))
  close(12)

  ! Convert units in meters
  h= h*1.d-2
  xmax= xmax*1.d-2
  ymax= ymax*1.d-2
  Lg= Lg*1.d-2
  Lh= Lh*1.d-2
  xg1= xg1*1.d-2
  xg2= xg2*1.d-2

  ! Direction of secondary emission
  if(flag_sec.ge.1) then
     i_tmp= INT(xg_sec/h(1)) + 1
     j_tmp= INT(yg_sec/h(2)) + 1
     if(bcnd(i_tmp-1,j_tmp).ge.1) then
        dir_sec= 1 ! Flux toward pos. x
     else
        dir_sec= -1
     endif
  endif

  ! Warning
  if(flag_sec.ge.1 .and. check_grd.eq.0) then
     if(mpi_rank.eq.0) then
        print*, 'Emmissive cathode wall index has not been defined. please correct...'
     endif
     call stop_calculation           
  endif
  
  return
end subroutine generate_boundary
