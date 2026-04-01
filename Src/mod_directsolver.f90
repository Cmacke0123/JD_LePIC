module mod_directsolver
  use mod_system
  
  implicit none
  
contains
  !-----------------------------------------------------
  subroutine initializepardiso(nx,ny,NxNy,jp1,jp2,bcnd,mpi_rank)
    !     ==============================================================
    !     Initialize Pardiso solver
    !     --------------------------------------------------------------
    use sparse_matrices    
    implicit none    
    integer:: i,j,itot,mpi_rank
    type(sparseMatrix):: M, M_CL
    type(sparseMatrix_LL):: M_LL, M_CL_LL    
    type(mySM_Row_LL), pointer:: row_ptr
    integer:: NxNy,jp1,jp2,nx,ny
    integer:: bcnd(0:nx+2,0:ny+2)
    
    ! Initialization - only one time
    IF (indx_init.EQ.1) THEN
       
       ! Matrix coefficients
       do j=jp1,jp2
          do i=0,nx
             IF (switch_perio.EQ.0) THEN
                itot=i+j*(nx+1)+1
             ELSE
                itot=i+(j-1)*(nx+1)+1
             ENDIF
             CALL mySM_LL_getRowPtr(M_LL, itot, row_ptr)
                         
             IF (i.GE.1.AND.i.LE.(nx-1).AND.j.GE.(jp1+1).AND.j.LE.(jp2-1)) then
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,-2d0*dx2i-2.0d0*dy2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot-1,1.0d0*dx2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot+1,1.0d0*dx2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot-(nx+1),1.0d0*dy2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot+(nx+1),1.0d0*dy2i))
             ENDIF
             
             ! i= 0
             IF (i.EQ.0) THEN
                 IF (bcnd(i+1,j+1).ge.1) THEN
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,1.0d0))
                ENDIF
                 IF (bcnd(i+1,j+1).eq.-2) THEN
                   IF (j.GE.(jp1+1).AND.j.LE.(jp2-1)) THEN
                      CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,-2.0d0*dx2i-2.0d0*dy2i)) 
                      CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot+1,2.0d0*dx2i)) 
                      CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot-(nx+1),1.0d0*dy2i)) 
                      CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot+(nx+1),1.0d0*dy2i))
                   ELSEIF (j.EQ.jp1) THEN
                      CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,-2.0d0*dx2i-1.0d0*dy2i)) 
                      CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot+1,2.0d0*dx2i)) 
                      CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot+(nx+1),1.0d0*dy2i)) 
                   ELSEIF (j.EQ.jp2) THEN
                      CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,-2.0d0*dx2i-1.0d0*dy2i)) 
                      CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot+1,2.0d0*dx2i)) 
                      CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot-(nx+1),1.0d0*dy2i))
                   ENDIF
                ENDIF
             ENDIF

             ! i= nx
             IF (i.EQ.nx) CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,1.0d0))
             
             ! j=jp1, i>0 & <nx
             IF (j.EQ.jp1.and.i.NE.0.AND.i.NE.nx) THEN
                 IF (bcnd(i+1,j+1).ge.1) THEN
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,1.0d0))
                ENDIF
                IF (switch_perio.EQ.1) THEN
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,-2.0d0*dx2i-2.0d0*dy2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot-1,1.0d0*dx2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot+1,1.0d0*dx2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(NxNy-(nx+1)+i+1,1.0d0*dy2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot+(nx+1),1.0d0*dy2i))
                ENDIF
             ENDIF
             
             ! j=jp2, i>0 & <nx
             IF (j.EQ.jp2.and.i.NE.0.AND.i.NE.nx) THEN
                IF (bcnd(i+1,j+1).ge.1) THEN 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,1.0d0))
                ENDIF
                IF (switch_perio.EQ.1) THEN
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,-2.0d0*dx2i-2.0d0*dy2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot-1,1.0d0*dx2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot+1,1.0d0*dx2i)) 
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot-(nx+1),1.0d0*dy2i))
                   CALL mySM_rll_inc(row_ptr, mySM_mk_elt(i+1,1.0d0*dy2i))
                ENDIF
             ENDIF
             
          ENDDO
       ENDDO
       
       ! Dirichle BC's inside simulation domain
       CALL mySM_spm_LL_to_CSR(M_LL,M)
       DO j=jp1,jp2
          DO i=0,nx
             IF (switch_perio.EQ.0) THEN
                itot=i+j*(nx+1)+1
             ELSE
                itot=i+(j-1)*(nx+1)+1
             ENDIF
             CALL mySM_LL_getRowPtr(M_CL_LL, itot, row_ptr)
             IF (bcnd(i+1,j+1).ge.1) then
                CALL mySM_rll_inc(row_ptr, mySM_mk_elt(itot,1.0d0))
             ENDIF
          ENDDO
       ENDDO
       CALL mySM_spm_LL_to_CSR(M_CL_LL,M_CL)
       CALL mySM_mergematrices_by_rows(M, M_CL, M_final)


       !-------pardiso init
       mtype = 11
       solver = 0
       pt=0 ! important !
       
       CALL pardisoinit(pt,mtype,iparm)
       
       maxfct= 1      !Number of numerical factorizations in memory
       mnum = 1       !Actual matrix to factorize	
       nrhs = 1       !Number of right-hand sides
       msglvl = 0     !no information printed to the screen (default=1)

       !.. Reordering and Symbolic Factorization, This step also allocates all memory that is necessary for the factorization

       phase = 11 ! only reordering and symbolic factorization
!        iparm(33) = 0
       
       CALL pardiso (pt, maxfct, mnum, mtype, phase, NxNy, M_final%a, M_final%ia, M_final%ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)

       if(mpi_rank.eq.0) then
          write(*,*) '-------------------------------------------------------------'
          write(*,*) ' Informations - PARDISO  '
          write(*,*) 'Reordering completed ... '
          if (error .NE. 0) then
             write(*,*) 'The following ERROR was detected: ', error
             STOP
          endif
          write(*,*) 'Number of nonzeros in factors = ',iparm(18)
          write(*,*) 'Number of factorization MFLOPS = ',iparm(19)
          write(*,*) '-------------------------------------------------------------'
       endif
       
       !.. Factorization.
       phase = 22 ! only factorization

       CALL pardiso (pt, maxfct, mnum, mtype, phase, NxNy, M_final%a, M_final%ia, M_final%ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)

       IF(mpi_rank.eq.0) THEN
          WRITE(*,*) 'Factorization completed ... '
          IF (error .NE. 0) THEN
             WRITE(*,*) 'The following ERROR was detected: ', error
             STOP
          ENDIF
       ENDIF
       
       phase = 33 ! only solve
       indx_init=0
       
    ENDIF ! end-if flag indx_init=1

    CALL mySM_dealloc(M)   

  end subroutine initializepardiso

  subroutine runpardiso(nx,ny,NxNy,jp1,jp2,phi,rhs)
    !     ==============================================================
    !     Run Pardiso. Call initializepardiso first.
    !     --------------------------------------------------------------
    implicit none
    integer :: i,j,itot
    double precision, dimension(NxNy) :: rhs,phi_temp
    real(kind=8):: phi(0:nx+2,0:ny+2)
    integer:: NxNy,jp1,jp2,nx,ny
    
    CALL pardiso (pt, maxfct, mnum, mtype, phase, NxNy, M_final%a, M_final%ia, M_final%ja,idum, nrhs, iparm, msglvl, rhs,phi_temp, error)
    
    do j=jp1,jp2
       do i=0,nx						
          IF (switch_perio.EQ.0) THEN
             itot=i+j*(nx+1)+1
          ELSE
             itot=i+(j-1)*(nx+1)+1
          ENDIF
          phi(i+1,j+1)=phi_temp(itot)
       enddo
    enddo
    
    if (switch_perio.NE.0) then
       do i=0,nx
          phi(i+1,1)=phi(i+1,ny+1)		
          phi(i+1,ny+2)=phi(i+1,2)		
       enddo
    endif
    
  end subroutine runpardiso
  
end module mod_directsolver
