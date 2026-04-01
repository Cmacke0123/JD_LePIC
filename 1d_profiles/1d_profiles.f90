program prof
  implicit none
  integer:: i,ic,opt,ix,iy,ix0,iy0,n1,n2,nx,ny
  parameter (n1=1024,n2=1536)
  real(kind=8):: dx1,dx2,dx3,k,Itot,ni(1:n1+1,1:n2+1)
  character:: name*20,ans*1

  Itot= 0.d0

  open(9,file='namelist.inp')
  read(9,*,end=1000) name
  print*, 'Reading file: ',name

  do ic=1,20
     if( name(ic:ic).eq.'.' .or. name(ic:ic).eq.' ' ) exit
  enddo

  ! Open file
  open(10,file='../DATA/'//name(1:ic-1)//'.dat',status='OLD',form='UNFORMATTED')

  ! Read file
  read(10,end=1000) nx,ny
  read(10,end=1000) ni(1:nx+1,1:ny+1)

  ! Plot options
  print*, 'along Ox=1, Oy=2?'
  read(*,*) opt

  if(opt.eq.1) then
     print*, 'ix=?'
     read(*,*) ix0
     open(11,file='1d_prof_ix.dat')
  else
     print*, 'iy=?'
     read(*,*) iy0
     open(11,file='1d_prof_iy.dat')
  endif
  
  k=1.d0
  print*, 'Do you want to calculate currents (y/n)?'
  read(*,*) ans
  if(ans.eq.'y'.or.ans.eq.'Y') then
     if(opt.eq.1) print*, 'dx, dy, z (mm)=?'
     if(opt.eq.2) print*, 'dy, dx, z (mm)=?'
     read(*,*) dx1,dx2,dx3
     dx1=dx1*1.d-3
     dx2=dx2*1.d-3
     dx3=dx3*1.d-3

     k=1.6e-19*dx1
  endif

  ! Draw plot
  do iy=1,ny+1
     do ix=1,nx+1
        
        if( opt.eq.1 .and. ix.eq.ix0 ) then
           write(11,*) iy,ni(ix,iy)*k
           Itot= Itot + ni(ix,iy)*k*dx2*dx3
        endif
        if( opt.eq.2 .and. iy.eq.iy0 ) then
           write(11,*) ix,ni(ix,iy)*k
           Itot= Itot + ni(ix,iy)*k*dx2*dx3
        endif
        
     enddo
  enddo

  goto 1001

  ! Messages
1000 continue
  print*, 'Error opening file'

1001 continue
  print*, 'End of file reached'

  if(ans.eq.'y'.or.ans.eq.'Y') print*, 'Itot (A)=',Itot

  if(opt.eq.1) then
     print*, '<.>=',ix0,sum(ni(ix0,1:ny))/ny
  else
     print*, '<.>=',sum(ni(1:nx,iy0))/nx
  endif

end program prof
