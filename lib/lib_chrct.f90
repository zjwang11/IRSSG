module chrct 

    use lib_comms
    implicit none 

contains

subroutine irssg_chrct(num_sym, SU2, &
                       num_bands, bot_band, top_band, ene_bands, &
                       dim_basis, num_basis, & 
                       coeffup_basis, coeffdn_basis, & 
                       chrct_sum, deg_set, numdeg_set, &
                       map_mat, map_vec, G_phase_pw, tolE_)

    !!! This subroutine calculates the characters
    !!! from bot_band to top_band on one kpoints for all operations.

    integer,     intent(in)  :: num_sym
    complex(dp), intent(in)  :: SU2(2,2,num_sym)

    

    integer,     intent(in)  :: bot_band, top_band 
    integer,     intent(in)  :: num_bands 
    real(dp),    intent(in)  :: ene_bands(num_bands)

    integer,     intent(in)  :: dim_basis, num_basis 
    complex(dp), intent(in)  :: coeffup_basis(dim_basis, num_bands)
    complex(dp), intent(in)  :: coeffdn_basis(dim_basis, num_bands)
    
    complex(dp), intent(out) :: chrct_sum(num_bands, MAXSYM)
    integer,     intent(out) :: deg_set(num_bands)
    integer,     intent(out) :: numdeg_set(num_bands)

    real(dp), optional,    intent(in)  :: tolE_


    complex(dp), intent(in), optional :: map_mat(dim_basis,dim_basis,num_sym)
    integer,     intent(in), optional :: map_vec(dim_basis,num_sym)
    complex(dp), intent(in), optional :: G_phase_pw(dim_basis,num_sym)
    
    integer     :: i, j, j2, iband, iope, kope 
    integer     :: num_deg, cnt_deg_set 
    complex(dp) :: csum
    integer  :: m
    complex(dp) :: gphase

    complex(dp) :: coeffup_norm(dim_basis, num_bands)
    complex(dp) :: coeffdn_norm(dim_basis, num_bands)

    complex(dp) :: coeffup_tmp(dim_basis, MAXDG)
    complex(dp) :: coeffdn_tmp(dim_basis, MAXDG)

    real(dp) :: tolE
    if (present(tolE_)) then 
        tolE = tolE_
    else
        tolE = 1.E-3_dp
    endif

    chrct_sum = cmplx(7.777777_dp, 7.777777_dp, dp)

    ! normalization
    coeffup_norm = 0.d0
    coeffdn_norm = 0.d0 
    do j = 1, num_bands
        csum = cmplx_0
        csum = dot_product(coeffup_basis(:,j), coeffup_basis(:,j))
        csum = csum + dot_product(coeffdn_basis(:,j),coeffdn_basis(:,j))
        csum = zsqrt(csum)
        coeffup_norm(:,j) = coeffup_basis(:,j)/csum
        coeffdn_norm(:,j) = coeffdn_basis(:,j)/csum 
    enddo 
    !coeffup_norm = coeffup_basis
    !coeffdn_norm = coeffdn_basis 

    ! calculate characters
    deg_set     = 0
    cnt_deg_set = 0
    iband = bot_band 

    DO WHILE (iband <= top_band)

        ! find the degenerate sets
        num_deg = 1
        cnt_deg_set = cnt_deg_set + 1
        if (iband < top_band) then
            do while (ene_bands(iband+num_deg)-ene_bands(iband) <= tolE)
                num_deg = num_deg + 1
                if ((iband+num_deg) >  top_band) exit 
            enddo 
        endif 

        do i = iband, iband+num_deg-1 
            deg_set(i)    = cnt_deg_set
            numdeg_set(i) = num_deg

            if (save_kcount <= MAXKSAVE) save_numdeg(i,save_kcount) = num_deg 
        enddo 
        
        if (num_deg > MAXDG) then 
            write(6, '(A,I3)') "WARNING: num_deg is larger than MAXDG for the band:", iband
            iband = iband + num_deg 
            cycle 
        endif 

        ! for each operations in little group
        do iope = 1, num_sym
            kope = iope
           
            ! rotation on spinor part

            j2 = 0
            coeffup_tmp = cmplx_0
            coeffdn_tmp = cmplx_0 
            do j = iband, iband + num_deg - 1
                j2 = j2 + 1
                do i = 1, num_basis 
                    coeffup_tmp(i,j2) = (SU2(1,1,kope)*coeffup_norm(i,j)) + &
                                        (SU2(1,2,kope)*coeffdn_norm(i,j))
                    coeffdn_tmp(i,j2) = (SU2(2,1,kope)*coeffup_norm(i,j)) + &
                                        (SU2(2,2,kope)*coeffdn_norm(i,j))
                enddo 
            enddo 


    
            ! rotation effect on the other part encoded in map_basis
            ! which should be provided in the main code
            chrct_sum(iband, kope) = cmplx_0
            do j = 0, num_deg - 1
                csum = cmplx_0
                if (present(map_mat)) then 

                    csum = dot_product(coeffup_norm(:,iband+j), &
                            matmul(map_mat(:,:,kope), coeffup_tmp(:,1+j))) + &
                           dot_product(coeffdn_norm(:,iband+j), &
                            matmul(map_mat(:,:,kope), coeffdn_tmp(:,1+j))) 
                else if (present(map_vec)) then 
                    do i = 1, num_basis  
                        m = map_vec(i,kope)
                        gphase = G_phase_pw(i,kope)
                        csum = csum + dconjg(coeffup_norm(m,iband+j))*coeffup_tmp(i,1+j)*gphase + &
                                      dconjg(coeffdn_norm(m,iband+j))*coeffdn_tmp(i,1+j)*gphase 
                    enddo 
                else
                    stop "The basis transformation under symmetry operations are not given"
                endif 
            
                chrct_sum(iband, kope) = chrct_sum(iband, kope) + csum 
            enddo 

            ! set all degenerate states the same character
            do j = iband+1, iband+num_deg-1
                chrct_sum(j, kope) = chrct_sum(iband, kope)
            enddo 

        enddo ! iope

        iband = iband + num_deg 

    ENDDO ! while  

end subroutine irssg_chrct


subroutine irssg_reps(num_litt_group_unitary,litt_group_unitary,num_bands,bot_band, top_band,EE,&
                            irrep_num, ch_table,irrep_name_list,&
                             chrct_set, deg_set, numdeg_set, kphase, numreps_set)
    integer, intent(in)            :: num_litt_group_unitary
    integer, intent(in)            :: litt_group_unitary(num_litt_group_unitary)
    integer, intent(in)            :: num_bands
    integer, intent(in)            :: bot_band, top_band
    integer, intent(in)             :: irrep_num
    complex(dp), intent(in)        :: chrct_set(num_bands, MAXSYM)
    complex(dp), intent(in)        :: kphase(num_litt_group_unitary)
    integer,     intent(in)        :: deg_set(num_bands)
    integer,     intent(in)        :: numdeg_set(num_bands)
    integer,     intent(out)       :: numreps_set(num_bands)
    complex(dp), intent(in)        :: ch_table(irrep_num, num_litt_group_unitary)
    character(len=15)               :: irrep_name_list(irrep_num)
    real(dp), intent(in)           ::  EE(num_bands)

    integer :: kope, iope
    complex(dp) :: chrct_phase(num_bands, num_litt_group_unitary)
    integer :: iband

    integer :: i1,j1,i2,j2,jl
    logical     :: match 
    real(dp) :: diff
    integer :: nir
    character(len=40) :: fmt

    integer :: record

    record = 0
    
    write(fmt,'(I0)')num_litt_group_unitary
#ifdef TEST
    write(10450,'(A,'//trim(adjustl(fmt))//'I20)')'     ',litt_group_unitary(1:num_litt_group_unitary)
#endif   

    write(*,*)
    write(*,'(A)')'Band Characters and Co-representations:'
    write(*,'(A)',advance='no')' bnd ndg    eigval'
    write(180,'(A)',advance='no')' bnd ndg    eigval'
    write(180,'(1I10,'//trim(adjustl(fmt))//'I13)')litt_group_unitary(1:num_litt_group_unitary)
    if (num_litt_group_unitary<=6) then
        write(*,'(1I10,'//trim(adjustl(fmt))//'I13)')litt_group_unitary(1:num_litt_group_unitary)
    else
        write(*,'(1I10,2I13,A,1I11,2I13)')litt_group_unitary(1),litt_group_unitary(2:3),'        ... ...',litt_group_unitary(num_litt_group_unitary-2:num_litt_group_unitary)
    endif

    write(fmt,'(I0)')num_litt_group_unitary*2
    
    nir = irrep_num
    ! chrct_phase(:,:) = conjg(chrct_set(:,:))
    chrct_phase(:,:) = chrct_set(:,1:num_litt_group_unitary)
    do iope = 1, num_litt_group_unitary
        kope = iope
        ! chrct_phase(:,kope) = chrct_phase(:,kope)*conjg(kphase(kope))*conjg(kphase_difftau(kope))
        chrct_phase(:,kope) = chrct_phase(:,kope)*kphase(kope)!*kphase_difftau(kope)
    enddo 

    iband = bot_band 
    do while (iband <= top_band) 
        do i1 = 1,irrep_num
            match = .true. 
            do iope = 1, num_litt_group_unitary
                kope = iope
                diff = abs(chrct_phase(iband,kope)-ch_table(i1,kope))
                ! write(*,*)i1,diff - (dble(nir)*TOLREP)
                if (diff > (dble(nir)*TOLREP)) then 
                    match = .false.
                    exit
                endif 
            enddo 
            if (match) exit 
        enddo 
        if (match) then 
            numreps_set(iband) = 1
            if (save_kcount <= MAXKSAVE) then 
                save_chrct(iband, save_kcount) = i1 
                save_numrep(iband, save_kcount) = 1
            endif 

#ifdef TEST
            write(10450,'(1I5,'//trim(adjustl(fmt))//'F10.5,1I5)')iband,chrct_phase(iband,1:num_litt_group_unitary),i1
#endif

            write(*,'(1I4,1I3,1E12.4E2,A)',advance='no')iband,nint(real(chrct_phase(iband,1))),EE(iband),' '
            write(180,'(1I4,1I3,1E12.4E2,A)',advance='no')iband,nint(real(chrct_phase(iband,1))),EE(iband),' '
            record = iband

            if(num_litt_group_unitary<=6) then
                do jl = 1, num_litt_group_unitary
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo
            else
                do jl = 1, 3
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo

                write(*,'(A)',advance='no')'   ... ...   '
                do jl = num_litt_group_unitary-2, num_litt_group_unitary
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo
            endif

            write(*,'(A)')'   = '//trim(adjustl(irrep_name_list(i1)))

            ! write(180,'(1I4,A,1E12.4E2,A)',advance='no')iband,' ',EE(iband),' '

            do jl = 1, num_litt_group_unitary
                if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                    write(180,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                else
                    write(180,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                endif
            enddo

            write(180,'(A)')'   = '//trim(adjustl(irrep_name_list(i1)))
        endif 

        if (.not.match) then 
    two:do i1 = 1, irrep_num
        do j1 = i1, irrep_num
            match = .true.
            do iope = 1, num_litt_group_unitary
                kope = iope
                diff = abs(chrct_phase(iband,kope)-ch_table(i1,kope) &
                                                    -ch_table(j1,kope))

                if (diff > (dble(nir)*TOLREP)) then 
                    match = .false.
                    exit
                endif 
            enddo 
            if (match) exit two 
        enddo 
        enddo two 
        if (match) then 
            numreps_set(iband) = 2
            if (save_kcount <= MAXKSAVE) then 
                save_chrct(iband,  save_kcount)  = i1
                save_chrct(iband+1,save_kcount)  = j1
                save_numrep(iband, save_kcount)  = 2
                save_numrep(iband+1,save_kcount) = 2
            endif 
#ifdef TEST
            write(10450,'(1I5,'//fmt//'F10.5,2I5)')iband,chrct_phase(iband,1:num_litt_group_unitary),i1,j1
#endif
            
            write(*,'(1I4,1I3,1E12.4E2,A)',advance='no')iband,nint(real(chrct_phase(iband,1))),EE(iband),' '
            write(180,'(1I4,1I3,1E12.4E2,A)',advance='no')iband,nint(real(chrct_phase(iband,1))),EE(iband),' '
            record = iband

            if(num_litt_group_unitary<=6) then
                do jl = 1, num_litt_group_unitary
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo
            else
                do jl = 1, 3
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo

                write(*,'(A)',advance='no')'   ... ...   '
                do jl = num_litt_group_unitary-2, num_litt_group_unitary
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo
            endif

            write(*,'(A)')'   = '//trim(adjustl(irrep_name_list(i1)))//' + '//trim(adjustl(irrep_name_list(j1)))

            ! write(180,'(1I4,A,1E12.4E2,A)',advance='no')iband,' ',EE(iband),' '
            do jl = 1, num_litt_group_unitary
                if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                    write(180,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                else
                    write(180,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                endif
            enddo

            write(180,'(A)')'   = '//trim(adjustl(irrep_name_list(i1)))//' + '//trim(adjustl(irrep_name_list(j1)))

        endif 
        endif 

    if (.not.match) then 
    thr:do i1 = 1, irrep_num
        do j1 = i1, irrep_num
        do i2 = j1, irrep_num
            match = .true.
            do iope = 1, num_litt_group_unitary
                kope = iope
                diff = abs(chrct_phase(iband,kope)-ch_table(i1,kope) &
                                                -ch_table(j1,kope) &
                                                -ch_table(i2,kope))

                if (diff > (dble(nir)*TOLREP)) then 
                    match = .false.
                    exit
                endif 
            enddo 
            if (match) exit thr
        enddo 
        enddo 
        enddo thr  
        if (match) then 
            numreps_set(iband) = 3
            if (save_kcount <= MAXKSAVE) then 
                save_chrct(iband, save_kcount) = i1
                save_chrct(iband+1,save_kcount) = j1
                save_chrct(iband+2,save_kcount) = i2
                save_numrep(iband, save_kcount) = 3
                save_numrep(iband+1, save_kcount) = 3
                save_numrep(iband+2, save_kcount) = 3
            endif 
#ifdef TEST
            write(10450,'(1I5,'//fmt//'F10.5,3I5)')iband,chrct_phase(iband,1:num_litt_group_unitary),i1,j1,i2
#endif

            write(*,'(1I4,1I3,1E12.4E2,A)',advance='no')iband,nint(real(chrct_phase(iband,1))),EE(iband),' '
            write(180,'(1I4,1I3,1E12.4E2,A)',advance='no')iband,nint(real(chrct_phase(iband,1))),EE(iband),' '
            record = iband

            if(num_litt_group_unitary<=6) then
                do jl = 1, num_litt_group_unitary
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo
            else
                do jl = 1, 3
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo

                write(*,'(A)',advance='no')'   ... ...   '
                do jl = num_litt_group_unitary-2, num_litt_group_unitary
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo
            endif

            write(*,'(A)')'   = '//trim(adjustl(irrep_name_list(i1)))//' + '//trim(adjustl(irrep_name_list(j1)))//' + '//trim(adjustl(irrep_name_list(i2)))


            ! write(180,'(1I4,A,1E12.4E2,A)',advance='no')iband,' ',EE(iband),' '
            do jl = 1, num_litt_group_unitary
                if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                    write(180,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                else
                    write(180,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                endif
            enddo

            write(180,'(A)')'   = '//trim(adjustl(irrep_name_list(i1)))//' + '//trim(adjustl(irrep_name_list(j1)))//' + '//trim(adjustl(irrep_name_list(i2)))

        endif 
        endif 
    if (.not.match) then 
    fou:do i1 = 1, irrep_num
        do j1 = i1, irrep_num
        do i2 = j1, irrep_num
        do j2 = i2, irrep_num
            match = .true.
            do iope = 1, num_litt_group_unitary
                kope = iope
                diff = abs(chrct_phase(iband,kope)-ch_table(i1,kope) &
                                                    -ch_table(j1,kope) &
                                                    -ch_table(i2,kope) &
                                                    -ch_table(j2,kope))
                if (diff > (dble(nir)*TOLREP)) then 
                    match = .false.
                    exit
                endif 
            enddo 
            if (match) exit fou 
        enddo 
        enddo 
        enddo 
        enddo fou
        if (match) then 
            numreps_set(iband) = 4
            if (save_kcount <= MAXKSAVE) then 
                save_chrct(iband, save_kcount)   = i1
                save_chrct(iband+1, save_kcount) = j1
                save_chrct(iband+2, save_kcount) = i2 
                save_chrct(iband+3, save_kcount) = j2 
                save_numrep(iband, save_kcount) = 4
                save_numrep(iband+1, save_kcount) = 4
                save_numrep(iband+2, save_kcount) = 4
                save_numrep(iband+3, save_kcount) = 4
            endif 
#ifdef TEST
            write(10450,'(1I5,'//fmt//'F10.5,4I5)')iband,chrct_phase(iband,1:num_litt_group_unitary),i1,j1,i2,j2
#endif

            write(*,'(1I4,1I3,1E12.4E2,A)',advance='no')iband,nint(real(chrct_phase(iband,1))),EE(iband),' '
            write(180,'(1I4,1I3,1E12.4E2,A)',advance='no')iband,nint(real(chrct_phase(iband,1))),EE(iband),' '
            record = iband

            if(num_litt_group_unitary<=6) then
                do jl = 1, num_litt_group_unitary
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo
            else
                do jl = 1, 3
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo

                write(*,'(A)',advance='no')'   ... ...   '
                do jl = num_litt_group_unitary-2, num_litt_group_unitary
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo
            endif

            write(*,'(A)')'   = '//trim(adjustl(irrep_name_list(i1)))//' + '//trim(adjustl(irrep_name_list(j1)))//' + '//trim(adjustl(irrep_name_list(i2)))//' + '//trim(adjustl(irrep_name_list(j2)))

            ! write(180,'(1I4,A,1E12.4E2,A)',advance='no')iband,' ',EE(iband),' '

            do jl = 1, num_litt_group_unitary
                if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                    write(180,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                else
                    write(180,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                endif
            enddo

            write(180,'(A)')'   = '//trim(adjustl(irrep_name_list(i1)))//' + '//trim(adjustl(irrep_name_list(j1)))//' + '//trim(adjustl(irrep_name_list(i2)))//' + '//trim(adjustl(irrep_name_list(j2)))



        endif 
        endif 

        if (.not. match) then
#ifdef TEST
            write(10450,'(1I5,'//fmt//'F10.5,A)')iband,chrct_phase(iband,1:num_litt_group_unitary),'  ??'
#endif

            write(*,'(1I4,1I3,1E12.4E2,A)',advance='no')iband,nint(real(chrct_phase(iband,1))),EE(iband),' '
            write(180,'(1I4,1I3,1E12.4E2,A)',advance='no')iband,nint(real(chrct_phase(iband,1))),EE(iband),' '
            record = iband

            if(num_litt_group_unitary<=6) then
                do jl = 1, num_litt_group_unitary
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo
            else
                do jl = 1, 3
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo

                write(*,'(A)',advance='no')'   ... ...   '
                do jl = num_litt_group_unitary-2, num_litt_group_unitary
                    if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                    else
                        write(*,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                    endif
                enddo
            endif

            write(*,'(A)')'   = ??'


            ! write(180,'(1I4,A,1E12.4E2,A)',advance='no')iband,' ',EE(iband),' '
            do jl = 1, num_litt_group_unitary
                if (aimag(chrct_phase(iband,jl))>=0.0_dp) then
                    write(180,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'+',aimag(chrct_phase(iband,jl))+1e-6,'i'
                else
                    write(180,'(A,1F5.2,A,1F4.2,A)',advance='no')'  ',real(chrct_phase(iband,jl))+1e-6,'-',-aimag(chrct_phase(iband,jl))+1e-6,'i'
                endif
            enddo

            write(180,'(A)')'   = ??'

        endif
        iband = iband + numdeg_set(iband)
    enddo

endsubroutine


end module chrct
