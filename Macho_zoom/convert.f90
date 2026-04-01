program convert
  implicit none
  integer:: i,ix,iy,nx,ny,n,np,ic,nfiles,flag_read,flag,every,lghx,lghy
  parameter (n=3200)
  real(kind=8):: ni(0:n,0:n,2),time,th(0:n,0:n)
  character:: name*20,dummy

  open(9,file='namelist.inp')
  
  print*, 'Plot every n points?'
  read(*,*) every

  flag_read=0
  do while(flag_read.eq.0)

     read(9,*,end=1002) name
     print*, 'Reading file: ',name

     flag=0

     do ic=1,20
        if( name(ic:ic).eq.'.' .or. name(ic:ic).eq.' ' ) exit
     enddo

     !
     ! Open and create file names
     !
     open(12,file=name(1:ic-1)//'.mco')
     open(13,file=name(1:ic-1)//'_1D.dat')
     
     ! Vectors
     if( name(ic-2:ic-1).eq.'j1' .or. name(ic-2:ic-1).eq.'j2' .or. &
          name(ic-2:ic-1).eq.'j3' .or. name(ic-2:ic-1).eq.'j4' .or. &
          name(ic-2:ic-1).eq.'j5' .or. name(ic-2:ic-1).eq.'j6' .or. &
          name(ic-2:ic-1).eq.'j7' .or. name(1:ic-1).eq.'B' ) then
        open(10,file='../DATA/'//name(1:ic-1)//'x.dat',form='UNFORMATTED')
        open(11,file='../DATA/'//name(1:ic-1)//'y.dat',form='UNFORMATTED')
        flag=1
     else
        open(10,file='../DATA/'//name(1:ic-1)//'.dat',form='UNFORMATTED')
     endif
  
     !
     ! Read file
     !

     read(10) nx,ny
     read(10) ni(1:nx+1,1:ny+1,1)
     if(flag.eq.1) then 
        read(11) nx,ny
        read(11) ni(1:nx+1,1:ny+1,2)
     endif

     print'(1x,"nx=",i5,", ny=",i5,", every=",i5)', nx,ny,every

     !
     ! Print in Macho format
     !
     lghx= 160 !256
     lghy= 160 !512
     write(12,*) lghx/every,lghy/every
     do iy=1536+lghy/2,1536-lghy/2,-1*every
        if(flag.eq.0) write(12,*) ( ni(ix,iy,1), ix=2048-lghx,2048,every  ) 
        if(flag.eq.1) write(12,*) ( dsqrt( ni(ix,iy,1)**2. + &
             ni(ix,iy,2)**2. ), ix=2048-lghx,2048,every  ) 
     enddo
     if(flag.eq.1) then
        write(12,*) 'vector'
        do iy=1536+lghy/2,1536-lghy/2,-1*every
           write(12,*) (datan2(ni(ix,iy,2),ni(ix,iy,1)), ix=2048-lghx,2048,every )  
        enddo
     endif

     !
     ! Save line plot
     !
     do iy=1536-768,1536+768
        if(flag.eq.0) write(13,*) iy,ni(2048,iy,1) 
        if(flag.eq.1) write(13,*) iy,dsqrt( ni(2048,iy,1)**2. + &
             ni(2048,iy,2)**2. )
     enddo

     close(10)
     close(11)
     close(12)
     close(13)

  enddo

1002 continue
  print*, 'End of file reached'
  
  close(9)

end program convert
