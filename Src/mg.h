!
! Matrix coefficients
!
  real(kind=8):: aw,ae,an,as,ac

! 
! Coefficients used in Poisson's equation and residual
! h(1)==dx and h(2)==dy
! ac= - an - as - aw - ae
!
  aw= eps0/(h(1)*h(1))
  ae= aw
  an= eps0/(h(2)*h(2))
  as= an
  ac= - ( an + as + aw + ae )
