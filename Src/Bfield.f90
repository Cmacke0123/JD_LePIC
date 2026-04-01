subroutine read_Bfield_map(Bi,n,name,namlen,B_dir,scaling,mpi_rank)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      Oct/11
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  implicit none
  include 'particle_info.h'
  integer:: i,ix,iy,nx,ny,n(3),B_dir,namlen,nmax,mpi_rank
  parameter (nmax=3*10**7)
  real(kind=8):: Bi(4,0:n(1)+2,0:n(2)+2),B,scaling
  character:: name*20

  if(mpi_rank.eq.0) print*, 'Opening: ',name(1:namlen)
  open(30,file='mag_field_files/'//name(1:namlen),status='OLD')

  nx= 0
  ny= 0

  ! Generate B-field map 
  do i=1,nmax
     read(30,*,end=999) ix,iy,B
     B= B*scaling
     Bi(B_dir,ix,iy)= Bi(B_dir,ix,iy) + B   
     nx= MAX(nx,ix)
     ny= MAX(ny,iy)
  enddo

999 continue

  if(mpi_rank.eq.0) &
       print'(1x,"nx= ",i4,", ny= ",i4,", scaling= ",f5.2)', nx,ny,scaling

  close(30)

  return

end subroutine read_Bfield_map

subroutine gaussian_Bfield(Bi,n,h,B0,x0,Lx,B_dir)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      Oct/11
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  implicit none
  integer:: ix,iy,n(3),B_dir
  real(kind=8):: h(3),Bi(4,0:n(1)+2,0:n(2)+2),x,B0,x0,Lx

  ! Generate B-field map 
  do ix=0,n(1)+2
     x= (ix-1)*h(1)
     do iy=0,n(2)+2
        Bi(B_dir,ix,iy)= Bi(B_dir,ix,iy) + & ! update B-field array
             B0*dexp( -(x-x0)**2/(2*Lx**2) )
     enddo
  enddo

  return

end subroutine gaussian_Bfield

subroutine PG_current(Bi,n,h,B0,B_dir)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      Jun/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:    1D fit for the magnetic filter field profile 
!                  in ELISE device, IPP Garching.
!                  Franzen et al., Plasma Phys. Control. Fusion 56 
!                  (2014) p. 025007, Fig 5.
!     --------------------------------------------------------------
  implicit none
  integer:: ix,iy,n(3),B_dir
  real(kind=8):: h(3),Bi(4,0:n(1)+2,0:n(2)+2), &
       a(5),b(2),x,B0,B_tmp
     
  a(1) = 0.0175596       
  a(2) = 0.0227362       
  a(3) = 0.00193113      
  a(4) = -8.51235e-05    
  a(5) = 9.58673e-07     

  b(1) = 15.7416
  b(2) = -0.377132

  ! Generate B-field map 
  do ix=0,n(1)+2

     x= (ix-1)*h(1)*1.d2 ! in cm
     if(x.le.39.d0) then 
        B_tmp=  a(1) + a(2)*x + a(3)*x**2 + a(4)*x**3 + &
             a(5)*x**4
     else
        B_tmp=  b(1) + b(2)*x
     endif

     do iy=0,n(2)+2
        Bi(B_dir,ix,iy)= Bi(B_dir,ix,iy) + B0*B_tmp ! update B-field array
     enddo

  enddo

  return

end subroutine PG_current

subroutine EE_magnets(Bi,n,h,B0,x0,y0,d,B_dir)
!     ==============================================================
!     VERSION:         0.3
!     LAST MOD:      Oct/19
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      d is the distance between the extraction magnets
!                    placed parallel to the (Ox) axis for B_dir=='Bx'
!                    and along (Oy) for B_dir=='By'. x0 is the position 
!                    of the magnets along (Ox), y0 along (Oy). 
!     --------------------------------------------------------------
  implicit none
  integer:: ix,iy,n(3),flag_simB,B_dir
  real(kind=8):: h(3),Bi(4,0:n(1)+2,0:n(2)+2),x,y,B0,x0,y0,d,ix0,iy0,&
       Bx,By
  include 'constants.h'

  ! Options
  ix0= 1
  iy0= 1
  flag_simB= 0
  if( B0.lt.0 .and. d.gt.0) then 
     flag_simB= 1
     if(B_dir.eq.1) iy0= n(2)/2+1
     if(B_dir.eq.2) ix0= n(1)/2+1
  endif
  if( B0.lt.0 .and. d.lt.0 ) flag_simB= 2
  d= ABS(d)
  B0= ABS(B0)

  ! Generate B-field map 
  do ix=n(1)+1,ix0,-1
     x= (ix-1)*h(1)
     do iy=n(2)+1,iy0,-1
        y= (iy-1)*h(2)

        if(B_dir.eq.1) then
           Bx= B0*dsin(pi*(x-x0)/d)*dexp(-pi*(y0-y)/d)
           By= -B0*dcos(pi*(x-x0)/d)*dexp(-pi*(y0-y)/d)
        endif
        if(B_dir.eq.2) then
           Bx= B0*dsin(pi*(y-y0)/d)*dexp(-pi*(x0-x)/d)
           By= B0*dcos(pi*(y-y0)/d)*dexp(-pi*(x0-x)/d)
        endif

        Bi(1,ix,iy)= Bi(1,ix,iy) + Bx ! Bx
        Bi(2,ix,iy)= Bi(2,ix,iy) + By ! By

        ! Mirror
        if(flag_simB.eq.1) then 
           if(B_dir.eq.1) then
              Bi(1,ix,n(2)+2-iy)= Bi(1,ix,iy) 
              Bi(2,ix,n(2)+2-iy)= -Bi(2,ix,iy) 
           endif
           if(B_dir.eq.2) then
              Bi(1,n(1)+2-ix,iy)= -Bi(1,ix,iy) 
              Bi(2,n(1)+2-ix,iy)= Bi(2,ix,iy) 
           endif
        endif

        ! Magnetic field lines in opposite direction
        if(flag_simB.eq.2) then 
           if(B_dir.eq.1) then 
              Bx= -B0*dsin(pi*(x-x0)/d)*dexp(-pi*y/d)
              By= -B0*dcos(pi*(x-x0)/d)*dexp(-pi*y/d)              
           endif
           if(B_dir.eq.2) then 
              Bx= B0*dsin(pi*(y-y0)/d)*dexp(-pi*x/d)
              By= -B0*dcos(pi*(y-y0)/d)*dexp(-pi*x/d)
           endif
           Bi(1,ix,iy)= Bi(1,ix,iy) + Bx ! Bx
           Bi(2,ix,iy)= Bi(2,ix,iy) + By ! By
        endif

     enddo
  enddo

  return

end subroutine EE_magnets

subroutine find_B_dir(B_info,B_dir)
!     ==============================================================
!     VERSION:         0.2
!     LAST MOD:      Oct/11
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS: 
!     --------------------------------------------------------------
  implicit none
  integer:: B_dir
  character:: B_info*2

  ! Find direction of magnetic field
  if( B_info.eq.'Bx' .or. B_info.eq.'BX' .or. &
       B_info.eq.'bx' .or. B_info.eq.'bX' ) B_dir=1
  if( B_info.eq.'By' .or. B_info.eq.'BY' .or. &
       B_info.eq.'by' .or. B_info.eq.'bY' ) B_dir=2
  if( B_info.eq.'Bz' .or. B_info.eq.'BZ' .or. &
       B_info.eq.'bz' .or. B_info.eq.'bZ' ) B_dir=3
  
  return

end subroutine find_B_dir
