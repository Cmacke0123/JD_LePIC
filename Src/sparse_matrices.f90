

!include 'mkl_pardiso.f90'
!include 'mkl_service.f90'

!pour toute la gestion des matrices creuses
module sparse_matrices
    !use mkl_pardiso
    !use global_data
    !use file_handling
    !use MKL_SERVICE
    
    implicit none


    
    type sparseMatrix !stockage en format CSR, one-based indexing, 3 array-variation
        double precision, allocatable :: a(:)  ! les éléments non nuls de a
        integer, allocatable          :: ia(:) ! index premier élement de chaque ligne dans a. Le dernier élément(nb lignes+1) contient le nombre d'élements de a+1
        integer, allocatable          :: ja(:) ! les colonnes des élements de a
    endtype sparseMatrix
    
    type mySM_Elt
        double precision                :: val
        integer                         :: key
    endtype mySM_Elt
    
    ! Type de définition d'une nouvelle ligne dans une matrice sparse
    ! Ne pas tenter de remplir le type 'à la main', mais utiliser une approcher du type :
    ! definir une variable mySM_row :: r
    ! definir un mySM_elt :: elt et lui donner l'élément de matrice qu'on veut ajouter
    ! faire successivement r=r+elt autant de fois que nécessaire
    ! Pour vider, simplement faire deallocate r
    type mySM_Row 
        type(mySM_Elt), allocatable     :: r(:)    
    endtype mySM_Row

    
    type sparseMatrix_LL !structure intermédiaire pour la création incrémentielle (ligne par ligne) d'une matrice creuse. Appeler ensuite la fonction adéquate pour la transformer dans le format CSR
        type(mySM_Row_LL), pointer      :: first_row=>null()
        type(mySM_Row_LL), pointer      :: last_get=>null() ! le dernier row renvoyé, pour acceler la création incrémentielle 
        integer                         :: last_get_i=0 ! la ligne du dernier renvoyé
    endtype sparseMatrix_LL
    
    type mySM_Row_LL ! linked list sur les lignes de matrices sparse pour une création efficace
        type(mySM_Row)                  :: row
        type(mySM_Row_LL), pointer      :: next_ptr => null()
    end type mySM_Row_LL
    
    !interface operator(+)
    !    module procedure mySM_Row_add_elt
    !end interface
    
    ! solver mkl pardiso
    integer, parameter                                      :: p_nsys = 1 ! le nombre de systèmes d'équations de structure différente qu'on va chercher à résoudre avec pardiso. Dans le but de ne faire qu'une analyse de structure
    integer, parameter                                      :: p_sys_pardiso=1! les id des différents systèmes d'équations
    logical,dimension(p_nsys)                               :: p_sys_cgs_precond_ok=.false. ! est ce qu'on a déjà calculé la "forme" de la matrice une première fois ?
    
    !TYPE(MKL_PARDISO_HANDLE), private,dimension(p_nsys,64)  :: p_pt ! utiliser p_nsys*64 pointeurs permet de conserver p_nsys structures de matrices
    !integer, private,dimension(p_nsys,64)                   :: p_iparm
    
    
    !TYPE(MKL_PARDISO_HANDLE), private,dimension(64)         :: p_pt            ! OLD
    !integer, private,dimension(64)                          :: p_iparm
    integer, private                                        :: p_error, p_msglvl
    
contains
    ! renvoie le pointeur vers une structure mySM_Row_LL. Remplir simplement le membre row par la suite
    subroutine mySM_LL_getRowPtr(spm_LL, row_n, row_LL_ptr)
        type(sparseMatrix_LL), INTENT(INOUT)    :: spm_LL
        integer, INTENT(IN)                     :: row_n
        type(mySM_Row_LL), pointer, INTENT(OUT) :: row_LL_ptr
        type(mySM_Row_LL), pointer              :: cur
        integer                                 :: i,start
    
        if(.not. associated(spm_LL%first_row)) then
            allocate(spm_LL%first_row)
            cur=>spm_LL%first_row
            start=1
        else
            cur=>spm_LL%first_row
            
            ! accelération si on profite d'une création par lignes consécutives :
            if(spm_LL%last_get_i.NE.0) then
                if(spm_LL%last_get_i.EQ.row_n) then
                    ! on renvoie le dernier récupéré :
                    row_LL_ptr=>spm_LL%last_get
                    return
                else if(row_n.GT. spm_LL%last_get_i) then
                    ! la ligne qu'on veut récupérer vient après, pas besoin de tout retraverser :
                    start=spm_LL%last_get_i
                    cur=>spm_LL%last_get
                else
                    !pas de chance, faut tout reprendre :
                    start=1
                end if
            end if
        end if
        

        
        do i=start,row_n-1
            if(.not.(associated(cur%next_ptr))) then
                allocate(cur%next_ptr)
                cur=>cur%next_ptr
            else
                cur=>cur%next_ptr
            end if
            
        end do
        row_LL_ptr=>cur
        spm_LL%last_get=>cur
        spm_LL%last_get_i=row_n
            
    end subroutine

    ! convertit la matrice sous forme Linked List en format CSR
    ! désalloue également  (si dealloc = FALSE, optionnel , dealloc = TRUE (default))
    subroutine mySM_spm_LL_to_CSR(spm_LL, spm_CSR,dealloc)
        type(sparseMatrix_LL), INTENT(INOUT)   :: spm_LL
        type(sparseMatrix),    INTENT(OUT)     :: spm_CSR
        logical, optional,INTENT(IN)           :: dealloc
        type(mySM_Row_LL), pointer             :: cur,cur2
        
        integer                                :: nrows,i, nelts, row_i, elt_i
        logical                                :: next
        logical                                :: do_cleanup
        
        if(present(dealloc)) then
            do_cleanup=dealloc
        else
            do_cleanup=.TRUE.
        end if
        
        
        ! recherche du nombre de lignes et du nombre d'éléments de la matrice :
        nrows=1
        nelts=0
        cur=>spm_LL%first_row
        if(allocated(cur%row%r)) then
            nelts=nelts+ubound(cur%row%r,1) ! nombre d'éléments dans cette ligne
        else
            nelts=nelts+1 ! on mettra un zéro en début de ligne pour la matrice csr
        end if
        
        
        do while(associated(cur%next_ptr))
            cur=>cur%next_ptr
            nrows=nrows+1
            if(allocated(cur%row%r)) then
                nelts=nelts+ubound(cur%row%r,1) ! nombre d'éléments dans cette ligne
            else
                nelts=nelts+1 ! on mettra un zéro en début de ligne pour la matrice csr
            end if
        end do
        
        ! allocation de la matrice csr: 
        allocate(spm_CSR%ia(1:(nrows+1)))
        allocate(spm_CSR%a(1:nelts))
        allocate(spm_CSR%ja(1:nelts))
        spm_CSR%ia=0
        spm_CSR%a=0.0d0
        spm_CSR%ja=0

        ! on retraverse une seconde fois la matrice en copiant cette fois les éléments :
        row_i=1
        elt_i=1
        
        cur=>spm_LL%first_row
        
        next=.TRUE.
        
        do while(next)
            spm_CSR%ia(row_i)=elt_i
            if(allocated(cur%row%r)) then
                ! on copie la ligne :
                do i=1, ubound(cur%row%r,1) !par construction, les éléments sont triés par ordre de colonne croissant
                    spm_CSR%ja(elt_i)=cur%row%r(i)%key
                    spm_CSR%a(elt_i)=cur%row%r(i)%val
                    elt_i=elt_i+1
                end do
                ! et on désalloue la ligne pour nettoyage :
                if(do_cleanup) then
                    deallocate(cur%row%r)
                end if
            else
                spm_CSR%ja(elt_i)=1
                spm_CSR%a(elt_i)=0.0d0
                elt_i=elt_i+1
            end if
            
            if(associated(cur%next_ptr)) then
                next=.TRUE.
                cur2=>cur%next_ptr
                if(do_cleanup) then
                    deallocate(cur)
                end if
                cur=>cur2
            else
                next=.FALSE. ! si il n'y a plus d'éléments, on arrête l'itération
                if(do_cleanup) then ! mais attention à quand même nettoyer la ligne actuelle
                    deallocate(cur)
                end if
            end if
            row_i=row_i+1
        end do
        
        if(do_cleanup) then
            nullify(spm_LL%first_row)
            nullify(spm_LL%last_get)
        end if
        
        if((elt_i.NE.(nelts+1)).OR. (row_i.NE.(nrows+1))) then
            print *, 'convert sp_LL to CSR : problème !'
        else
            spm_CSR%ia(row_i)=nelts+1
        end if
        
        
    end subroutine
    

    ! Affecte la matrice sparse sp_m à la diagonale vecteur spécifiée
    ! suppose le 1 based indexing 
    subroutine mySM_vectorDiag(sp_m,diag)
        type(sparseMatrix), intent(OUT)                 :: sp_m
        double precision, dimension(:),INTENT(IN)       :: diag
        
        integer                                         ::u,i
        
        u=ubound(diag,1)
        
        allocate(sp_m%a(1:u))
        allocate(sp_m%ja(1:u))
        allocate(sp_m%ia(1:u+1))
        
        sp_m%a=diag
        do i=1,u
            sp_m%ja(i)=i
            sp_m%ia(i)=i
        end do
        
        sp_m%ia(u+1)=u+1
    end subroutine
    
    subroutine mySM_getRow(sp_m, nrow,SM_row)
        type(sparseMatrix), intent(IN)                  :: sp_m
        integer,INTENT(IN)                              :: nrow
        type(mySM_Row), INTENT(OUT)                     :: SM_row
        integer                                         :: i,i1,i2
        
        if(nrow.LE.(ubound(sp_m%ia,1)-1)) then
            i1=sp_m%ia(nrow)
            i2=sp_m%ia(nrow+1)-1
            do i=i1,i2
                !SM_row=SM_row+mySM_mk_elt(sp_m%ja(i),sp_m%a(i))                
                CALL mySM_r_inc(SM_row,mySM_mk_elt(sp_m%ja(i),sp_m%a(i)))
            end do
        else
            print *, ' mySM_getRow : erreur : nrow .GT. nrows(sp_m)'
            print *, '',nrow,' > ',ubound(sp_m%ia,1)-1 
            STOP
        end if
        
    end subroutine
    subroutine mySM_setRow(sp_m, nrow,SM_row)
        type(sparseMatrix), intent(INOUT)               :: sp_m
        integer, intent(IN)                             :: nrow ! ligne
        type(mySM_Row), INTENT(IN)                      :: SM_row
        integer, allocatable                            :: cols(:)
        double precision,allocatable                    :: vals(:)
        integer                                         :: nvals
        
        integer                                         :: m,s,ds,nzeros,i,ua,uia,uib
        double precision,allocatable                    :: b(:)
        integer, allocatable                            :: ib(:), jb(:)
        
   
        
        if(.NOT.allocated(SM_row%r)) then
            return
        end if
        
        ! extraction de la table colonne et valeurs
        nvals=ubound(SM_row%r,1)
        allocate(cols(1:nvals))
        allocate(vals(1:nvals))
        do i=1,nvals
            cols(i)=SM_row%r(i)%key
            vals(i)=SM_row%r(i)%val
        end do
        
        
    
        ! nombre actuel de lignes
        if(.not. allocated(sp_m%ia)) then
            m=1
            ALLOCATE(sp_m%a(1))
            ALLOCATE(sp_m%ia(1:2))
            ALLOCATE(sp_m%ja(1))
            sp_m%a(1)=0.0d0
            sp_m%ia(1)=1
            sp_m%ia(2)=2
            sp_m%ja(1)=1
        else
            m=ubound(sp_m%ia,1)-1 ! on suppose le 1-indexing
        end if
        if(nrow.LE.0) then
            print *, 'sparse matrix, index négatif !!!!'
            STOP
        end if
        
        if(nrow.LE.m) then
            ! ok on est bien "dans" la matrice, reste plus qu'à modifier la ligne correspondante
            ! (en particulier, virer les 0 éventuels si la matrice n'a pas été construite séquentiellement ligne par ligne
            ! et qu'il nous en fallait un pour remplir en attendant)
            ! on va faire simple : c'est un set donc de toute façon on supprime la ligne
            
            ! on veut calculer le nombre d'élements que contiendra la matrice après l'insertion de la ligne
            ! mais si la ligne contient d'autres éléments insérés préalablement, il faut les supprimer
            ds=nvals-(sp_m%ia(nrow+1)-sp_m%ia(nrow))   ! fonctionne également pour la dernière ligne grace à l'élément supplémentaire à m+1
            ua=ubound(sp_m%a,1)
            s=ua+ds
            
            uia=ubound(sp_m%ia,1)
            
              ! on peut maintenant effectuer les allocations 
            allocate(b(1:s))
            allocate(ib(1:uia))
            allocate(jb(1:s))
            
            uib=uia
            
            
            ! on copie la partie intacte de la matrice originale
            if(nrow.GE.2) then ! le début 
                ib(1:(nrow-1))=sp_m%ia(1:(nrow-1))
                b(1:(sp_m%ia(nrow)-1))=sp_m%a(1:(sp_m%ia(nrow)-1))
                jb(1:(sp_m%ia(nrow)-1))=sp_m%ja(1:(sp_m%ia(nrow)-1))
            end if
            if(nrow.LT.m) then ! la fin
                ib((nrow+1):(uib-1))=sp_m%ia((nrow+1):(uia-1))+spread(ds,1,(uia-1)-(nrow+1)+1)
                ib(uib)=s+1
                jb((sp_m%ia(nrow+1)+ds):s)=sp_m%ja(sp_m%ia(nrow+1):ua)
                b(ib(nrow+1):s)=sp_m%a(sp_m%ia(nrow+1):ua)
            end if
            
            ! il reste à compléter la nouvelle ligne avec les valeurs à insérer :
            ib(nrow)=sp_m%ia(nrow)
            ib(uia)=s+1
            jb(ib(nrow):(ib(nrow+1)-1))=cols
            b(ib(nrow):(ib(nrow+1)-1))=vals

        else
            ! ici on ajoute plus loin que la 'fin' actuelle de la matrice
            ! on détermine si il faut ajouter des 0 avant la ligne qu'on veut insérer :
            nzeros=nrow-m-1
            ua=ubound(sp_m%a,1)
            s=ua+nzeros+nvals
            
            uia=ubound(sp_m%ia,1)
            
            allocate(b(1:s))
            allocate(ib(1:(uia+nzeros+1)))
            allocate(jb(1:s))
            b(1:ua)=sp_m%a(1:ua)
            jb(1:ua)=sp_m%ja(1:ua)
            ib(1:(uia-1))=sp_m%ia(1:(uia-1))
            if(nzeros.GT.0) then
                b((ua+1):(ua+nzeros))=0.0d0
                jb((ua+1):(ua+nzeros))=1 ! on met les zeros en début de ligne
                do i=1,nzeros
                    ib(uia+i-1)=ua+i
                end do
            end if
            ! on ajoute maintenant la ligne de valeurs
            b((ua+nzeros+1):s)=vals
            jb((ua+nzeros+1):s)=cols
            ib(uia+nzeros)=s-(nvals-1)
            ib(ubound(ib,1))=s+1
        end if
        
        ! transfert vers la variable de stockage : (et desallocation de b)
        CALL MOVE_ALLOC(b,sp_m%a)
        CALL MOVE_ALLOC(ib,sp_m%ia)
        CALL MOVE_ALLOC(jb,sp_m%ja)
        
        deallocate(cols)
        deallocate(vals)
    end subroutine mySM_setRow
    
    subroutine empty_mySM_Row(row)
        type(mySM_Row), intent(inout)   :: row
        integer                         :: err
        if(allocated(row%r)) then
            deallocate(row%r, STAT=err)
            if(err.NE.0) then
                print *, 'erreur dealloc row, arret'
                STOP
            end if
        end if
    end subroutine empty_mySM_Row
    
    subroutine mySM_dealloc(sp_m)
    ! à appeler en sortie de programme pour chaque matrice
        type(sparseMatrix)                              :: sp_m
        integer                                         :: st1,st2,st3
        
        st1=0
        st2=0
        st3=0
        
        if(allocated(sp_m%a)) then
            deallocate(sp_m%a,STAT=st1)
        end if
        if(allocated(sp_m%ja)) then
            deallocate(sp_m%ja,STAT=st2)
        end if
        if(allocated(sp_m%ia)) then
            deallocate(sp_m%ia,STAT=st3)
        end if
        
        if((st1.NE.0).OR.(st2.NE.0).OR.(st3.NE.0)) then
            print *, "erreur lors du dealloc d'une matrice sparse, arrêt"
            stop
        end if
        
        
    end subroutine mySM_dealloc
    
    ! equivalent de l'opération row+= elt que j'aurais aimé pouvoir surcharger... malheureusement, fortran ne veut pas me laisser faire
    ! à noter que faire row=row+elt produit probablement des memory leaks important (sur le row de RHS, qui se retrouve dans la nature sans deallocation)
    ! c'est pour ça que je pense que l'overload de + vers mySM_Row_add_elt(row,elt) est une mauvaise idée
    subroutine mySM_r_inc(row, elt)
        type(mySM_row), INTENT(INOUT)           :: row
        type(mySM_elt), INTENT(IN)              :: elt
        type(mySM_row)                          :: tmp_row
        
        tmp_row=mySM_row_add_elt(row,elt)
        CALL move_alloc(tmp_row%r,row%r)
        
    end subroutine
    
    subroutine mySM_rLL_inc(row_LL_ptr, elt)
        type(mySM_row_LL), pointer, INTENT(INOUT)  :: row_LL_ptr
        type(mySM_elt), INTENT(IN)              :: elt
        type(mySM_row)                          :: tmp_row
        
        tmp_row=mySM_row_add_elt(row_LL_ptr%row,elt)
        CALL move_alloc(tmp_row%r,row_LL_ptr%row%r)
        
        
    end subroutine
    
    function mySM_Row_add_elt(row, elt)
    ! ajoute un element à la ligne (simplement un tableau d'éléments, mais ordonné)
    ! utiliser UNIQUEMENT ceci pour créer une nouvelle ligne à insérer dans une matrice sparse
    ! l'intérêt étant que les valeurs sont effectivement classées par ordre croissant de colonne
    ! ce qui est nécessaire dans le format CSR
    
        type(mySM_Row), INTENT(IN)             :: row
        type(mySM_Row)                         :: mySM_Row_add_elt
        type(mySM_Elt), INTENT(IN)             :: elt
        
        integer                                :: i,ub
        
        
        if(elt%val.EQ.0.0d0) then
            if(allocated(row%r)) then
                allocate(mySM_Row_add_elt%r(1:ubound(row%r,1)))
                mySM_Row_add_elt%r=row%r
            end if
            return
        end if
        
        
        if(.NOT.allocated(row%r)) then
            allocate(mySM_Row_add_elt%r(1))
            mySM_Row_add_elt%r(1)=elt
            return
        else
            ub=ubound(row%r,1)
            do i=1,ub
                if(row%r(i)%key.EQ.elt%key) then
                    if(row%r(i)%val.NE.-elt%val) then 
                        allocate(mySM_Row_add_elt%r(1:ub))
                        mySM_Row_add_elt%r(1:ub)=row%r
                        mySM_Row_add_elt%r(i)%val=mySM_Row_add_elt%r(i)%val+elt%val
                    else ! on va obtenir un zéro
                        if(ub.GT.1) then 
                            allocate(mySM_Row_add_elt%r(1:(ub-1)))
                            mySM_Row_add_elt%r(1:(i-1))=row%r(1:(i-1))
                            mySM_Row_add_elt%r(i:(ub-1))=row%r((i+1):ub)
                       end if 
                    end if
                    return
                elseif(row%r(i)%key.GT.elt%key) then
                    ! on sait que notre vecteur est trié par ordre croissant, il ne peut rien y avoir plus loin d'interessant
                    ! on va donc ajouter ici même le nouvel élement, ce qui va garantir l'ordre voulu
                    allocate(mySM_Row_add_elt%r(1:(ub+1)))
                    if(i.GT.1) then
                       mySM_Row_add_elt%r(1:(i-1))=row%r(1:(i-1))
                    end if
                    mySM_Row_add_elt%r(i)=elt
                    
                    if(i.LE.ub) then
                       mySM_Row_add_elt%r((i+1):(ub+1))=row%r(i:ub) 
                    end if
                    return
                    
                end if
            
            end do
            ! si on arrive ici, il suffit d'ajouter l'élément à la fin
             allocate(mySM_Row_add_elt%r(ub+1))
             mySM_Row_add_elt%r(1:ub)=row%r
             mySM_Row_add_elt%r(ub+1)=elt

        end if
        
    end function mySM_Row_add_elt

    function mySM_mk_elt(col, val)
        integer,INTENT(IN)                         :: col
        double precision,INTENT(IN)                :: val
        type(mySM_elt)                             :: mySM_mk_elt
        
        mySM_mk_elt%key=col
        mySM_mk_elt%val=val
    end function mySM_mk_elt
    

    
    ! à utiliser à des fins de Debug, probablement assez lourd(contient un max)
    ! possible d'améliorer en prenant le max que sur les éléments sp_m%ja(sp_m%ia-1), et en supposant les colonnes ordonnées
    ! par ordre croissant, ce qui est normalement le cas, mais est ce vraiment utile ? Ce genre de fonction ne devrait pas servir
    ! dans le code où a priori , on sait ce qu'on fait
    subroutine mySM_size(sp_m, m,n)
        type(sparseMatrix),INTENT(IN)               :: sp_m  
        integer, INTENT(OUT)                        :: m,n
        
        if(allocated(sp_m%a)) then
            m=ubound(sp_m%ia,1)-1
            n=maxval(sp_m%ja)
        else
            m=0
            n=0
        end if
    
    end subroutine
    
    ! écit la matrice sparse spécifiée dans un fichier
    subroutine mySM_mattofile(sp_a, fname)
        type(sparseMatrix), INTENT(IN)          :: sp_a
        character(len=*), INTENT(IN)                  :: fname
        double precision, allocatable           :: v(:)
        integer                                 :: col,i,j
        type(mySM_Row)                          :: r
        
        
        
        col=maxval(sp_a%ja)

        allocate(v(1:col))
        
        OPEN(UNIT=999, FILE=TRIM(ADJUSTL(fname)), STATUS='REPLACE', recl=4*sizeof(v))
        
        do i=1, ubound(sp_a%ia,1)-1
            v=0.0d0
            CALL mySM_getRow(sp_a,i,r)
            do j=1,ubound(r%r,1)
                v(r%r(j)%key)=r%r(j)%val
            end do
            
            write(999,*) v(1:col)
            
            CALL empty_mySM_Row(r)
        end do
        
        deallocate(v)
        CLOSE(999)
    
    end subroutine
    ! copie la matrice a dans b, s'occupe de l'allocation de b
    subroutine mySM_copy(sp_a,sp_b)
        type(sparseMatrix), INTENT(IN)          :: sp_a
        type(sparseMatrix), INTENT(OUT)         :: sp_b
        
        CALL mySM_dealloc(sp_b)
        allocate(sp_b%ia(1:ubound(sp_a%ia,1)))
        allocate(sp_b%ja(1:ubound(sp_a%ja,1)))
        allocate(sp_b%a(1:ubound(sp_a%a,1)))
        
        sp_b%ia=sp_a%ia
        sp_b%ja=sp_a%ja
        sp_b%a=sp_a%a
    end subroutine
    
    ! fusionne les matrices a et b, en prenant par défaut les lignes de A, mais celles de B si non nulles
    ! le résultat va dans c
    subroutine mySM_mergeMatrices_by_Rows(sp_a, sp_b,sp_c)
        type(sparseMatrix), INTENT(IN)                  :: sp_a
        type(sparseMatrix), INTENT(IN)                  :: sp_b
        type(sparseMatrix), INTENT(OUT)                 :: sp_c
        
        integer                                         :: i,nrows,maxelts,ma,mb, n
        type(sparseMatrix)                              :: temp_sp_c
        
        logical                                         :: copyA
        
        ma=ubound(sp_a%ia,1)-1
        mb=ubound(sp_b%ia,1)-1
        nrows=max(ma,mb)
        !nb maximum d'éléments qu'on pourrait trouver :
        maxelts=ubound(sp_a%a,1)+ubound(sp_b%a,1)
        
        ! on alloue une matrice temporaire avec ces deux tailles :
        allocate(temp_sp_c%ia(1:(nrows+1)))
        allocate(temp_sp_c%ja(1:maxelts))
        allocate(temp_sp_c%a(1:maxelts))
            
        i=1
        n=1
        do while((i.LE.ma) .OR. (i.LE.mb))
            ! on parcourt les matrices ligne par ligne :
            
            ! test pour savoir d'où on doit copier la ligne
            if(i.LE.ma) then
                if(i.LE.mb) then
                    if(ANY(sp_b%a(sp_b%ia(i):(sp_b%ia(i+1)-1)).NE.0.0d0)) then
                        copyA=.FALSE. ! il y a un élément non nul dans B, on substitue la ligne
                    else
                        copyA=.TRUE. ! pas d'élément non nul dans B, on copie A
                    end if
                else
                    copyA=.TRUE. ! on dépasse la dernière ligne de B, on copie A nécessairement
                end if
            else
                copyA=.FALSE. ! copie de B si on a dépassé la dernière ligne de A
            end if
            
            
            
            if(.not.copyA) then
                !il y a un élément non nul dans cette ligne de B : on la copie donc dans C
                temp_sp_c%ia(i)=n
                temp_sp_c%a(n:(n+sp_b%ia(i+1)-sp_b%ia(i)-1))=sp_b%a(sp_b%ia(i):(sp_b%ia(i+1)-1))
                temp_sp_c%ja(n:(n+sp_b%ia(i+1)-sp_b%ia(i)-1))=sp_b%ja(sp_b%ia(i):(sp_b%ia(i+1)-1))
                n=n+sp_b%ia(i+1)-sp_b%ia(i)
            else
                ! on prend ici la ligne de A
                temp_sp_c%ia(i)=n
                temp_sp_c%a(n:(n+sp_a%ia(i+1)-sp_a%ia(i)-1))=sp_a%a(sp_a%ia(i):(sp_a%ia(i+1)-1))
                temp_sp_c%ja(n:(n+sp_a%ia(i+1)-sp_a%ia(i)-1))=sp_a%ja(sp_a%ia(i):(sp_a%ia(i+1)-1))
                n=n+sp_a%ia(i+1)-sp_a%ia(i)
            end if

            i=i+1
        end do
        ! dernier élément de ic (à la ligne nrows +1) :
        temp_sp_c%ia(i)=n !(nb elements +1)
        
        ! puis redimensionnement à la vraie taille :
        allocate(sp_c%ia(1:(nrows+1)))
        allocate(sp_c%a(1:(n-1)))
        allocate(sp_c%ja(1:(n-1)))
        
        sp_c%ia=temp_sp_c%ia
        sp_c%ja=temp_sp_c%ja(1:(n-1))
        sp_c%a=temp_sp_c%a(1:(n-1))
        
        deallocate(temp_sp_c%a)
        deallocate(temp_sp_c%ja)
        deallocate(temp_sp_c%ia)
    
    end subroutine
      
  end module sparse_matrices
