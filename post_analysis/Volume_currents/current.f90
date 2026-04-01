program current
  implicit none
  integer:: i,ix,iy,iz,ix0,n,num,nx,ny
  parameter (n=3200)
  real(kind=8):: time,ni(0:n,0:n),dx,dy,dz,k,Itot
  character:: dummy,pnum*1

  Itot= 0.d0
  ni= 0.d0

  open(11,file='input.dat',status='OLD')

  read(11,*) dx,dy,dz

  print*, 'ptype=?'
  read(*,*) num

  write(pnum,'(i1)'),num

  open(12,file='../../DATA/sour'//pnum//'.dat',status='OLD',form='UNFORMATTED')

  dx=dx*1.d-3
  dy=dy*1.d-3
  dz=dz*1.d-3
  k=1.6d-19*dx*dy*dz

  read(12) nx,ny
  read(12) ni(1:nx+1,1:ny+1)
  print'(1x,"nx=",i5,", ny=",i5,", every=",i5)', nx,ny

  print*, 'ix < ix0, ix0?'
  read(*,*) ix0

  do iy=1,ny+1
     do ix=1,nx+1
        if(ix.le.ix0) Itot= Itot + ni(ix,iy)*k
     enddo
  enddo

  print*, 'Itot (A)=',Itot

  close(11)
  close(12)

end program current
