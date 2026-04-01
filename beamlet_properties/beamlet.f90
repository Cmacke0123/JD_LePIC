program neg_ions
  implicit none
  integer:: i,nmax,iy,iy_min,iy_max,cnt_h,dih,Nh
  parameter (nmax=10**6)
  real(kind=8):: j_ion(7000),j_tmp,j_tot,j_avg,j_beam
  character ans*1
  
  open(10,file='../1d_profiles/1d_prof_ix.dat')
  open(12,file='impacts_avg.dat')

  iy_min=10**10
  iy_max=0
  j_ion=0.d0
  j_beam=0.d0
  cnt_h=0
  j_tot=0.d0

  do i=1,nmax
     read(10,*,end=999) iy,j_tmp
     j_ion(iy)= j_tmp
     j_tot= j_tot + j_tmp

     if(ABS(j_tmp).gt.0) then
        iy_min=MIN(iy_min,iy)
        iy_max=MAX(iy_max,iy)
     endif
     
  enddo

999 close(10)

print*, 'iy<=',iy_min,', iy>=',iy_max

  print*, 'Number of apertures along (Oy)?'
  read(*,*) Nh

  dih=NINT(real(iy_max-iy_min)/real(Nh))
  iy_min= iy_min - NINT(real(dih)/4.d0)
  iy_max= iy_max + NINT(real(dih)/4.d0)
  dih=NINT(real(iy_max-iy_min)/real(Nh))

  print*, 'Size of aperture:',dih

  j_avg= 1.d0/real(Nh) ! Normalized

  do iy=iy_min,iy_max
     if( cnt_h.eq.dih .or. iy.eq.iy_max ) then
           print*, '<j>(%)=',((j_beam-j_avg)/j_avg)*100.,', iy<, iy>=',iy-dih,iy
           write(12,*) iy,((j_beam-j_avg)/j_avg)*100.
           j_beam=0.d0
           cnt_h=0
     else
        ! Calculate average/hole
        j_beam= j_beam + j_ion(iy)/j_tot
        cnt_h= cnt_h + 1
     endif
  enddo

  close(11)
  close(12)

end program neg_ions
