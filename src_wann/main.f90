module driver_wann
contains
subroutine run_wann()
    use get_ssg_mod
    use comms_wann 
    use init_wann
    use wave_data_wann
    use linear_rep
    use comprel
    use kgroup_mod, only: kgroup
    implicit none 

    integer            :: kkk

    ! command argument 
    integer            :: narg, iarg, lens, stat 
    character(len=100) :: arg, cmd 

    integer,     allocatable :: rot_lg(:,:,:)
    real(dp),     allocatable :: spin_rot_lg(:,:,:)
    real(dp),    allocatable :: SO3_lg(:,:,:)
    real(dp),    allocatable :: tau_lg(:,:)
    complex(dp), allocatable :: SU2_lg(:,:,:)
    complex(dp), allocatable :: kphase(:)
    complex(dp), allocatable :: tilte_mat(:,:,:)
    complex(dp), allocatable :: Gphase(:,:)
    integer,     allocatable :: tilte_vec(:,:)
    integer,     allocatable :: litt_group(:)
    integer,     allocatable :: time_reversal_lg(:)
    integer     ::           numLG
    integer     ::          numULG
    integer, allocatable :: litt_group_unitary(:)
    integer, allocatable :: rot_unitary_lg(:,:,:)
    real(dp), allocatable :: SO3_unitary_lg(:,:,:)
    complex(dp), allocatable :: SU2_unitary_lg(:,:,:)
    real(dp), allocatable :: tau_unitary_lg(:,:)
    integer :: Ncoirrep
    integer :: Nirrep
    integer, allocatable :: op_order(:)
    complex(dp), allocatable :: ch_table(:,:)
    complex(dp), allocatable :: Ch_table2(:,:)

    complex(dp), allocatable :: ch_unitary_table(:,:)
    complex(dp), allocatable :: Ch_table1(:,:)

    integer, allocatable :: irrep_coirrep_relation(:,:)

    character(len=15),allocatable :: irrep_name_list(:)

    integer:: kchoose_start, kchoose_end
    logical :: kchoose_flag
    logical :: spin_no_reversal
    integer, allocatable :: torsion(:)
    integer, allocatable :: discriminant_value(:)

    character(len=4) :: k_name
    character(len=15) :: k_frac_symbol

    integer :: i, j
    complex(dp), allocatable :: phase(:)
    character(len=100) :: ssg_label
    real(dp) :: tolE
    integer :: spg

    

interface
subroutine irrep_ssg(numULG,litt_group_unitary, rot_lg, tau_lg,&
                      SO3_lg, SU2_lg, &
                      KKK, WK, kphase, &
                      num_bands, m, n, ene_bands, &
                      dim_basis, num_basis, &
                      coeffa, coeffb, Ncoirrep, ch_table, &
                      irrep_name_list, &
                      G_phase_pw, rot_vec_pw, rot_mat_tb, tolE_)

    integer, parameter :: dp = 8
    integer, intent(in) :: numULG
    integer,     intent(in) :: rot_lg(3,3,numULG)
    real(dp),    intent(in) :: tau_lg(3,numULG)
    real(dp),    intent(in) :: SO3_lg(3,3,numULG)
    complex(dp), intent(in) :: SU2_lg(2,2,numULG)
    integer,     intent(in) :: litt_group_unitary(numULG)
    integer,     intent(in) :: KKK 
    real(dp),    intent(in) :: WK(3)
    complex(dp), intent(in) :: kphase(numULG)

    integer,     intent(in) :: num_bands, m, n 
    real(dp),    intent(in) :: ene_bands(num_bands) 

    integer,     intent(in) :: dim_basis, num_basis  
    complex(dp), intent(in) :: coeffa(dim_basis, num_bands)
    complex(dp), intent(in) :: coeffb(dim_basis, num_bands)

    integer, intent(in) :: Ncoirrep
    character(len=15), intent(in) :: irrep_name_list(Ncoirrep)
    complex(dp), intent(in) :: ch_table(Ncoirrep,numULG)
    complex(dp), intent(in), optional :: G_phase_pw(dim_basis, numULG)
    integer,     intent(in), optional :: rot_vec_pw(dim_basis, numULG)
    complex(dp), intent(in), optional :: rot_mat_tb(dim_basis, dim_basis, numULG)
    real(dp),intent(in) :: tolE_


end subroutine irrep_ssg
end interface 
    kchoose_flag=.false.
    kchoose_start = 0
    kchoose_end = 0
    bot_band = 0
    top_band = 0
    tolE = 1.0e-3
    call get_command(cmd)
    write(*,*) 'Current command : ', trim(cmd) 
    narg = command_argument_count()
    write(*,*) 'Argument count : ', narg 

    iarg = 1
    do while (.true.)
        call get_command_argument(iarg, arg, lens, stat)
        if (len_trim(arg) == 0) exit 
        if (trim(arg) == '-nk') then
            kchoose_flag = .true.
            iarg = iarg + 1
            call get_command_argument(iarg, arg, lens, stat)
            read(arg, *) kchoose_start
            iarg = iarg + 1
            call get_command_argument(iarg, arg, lens, stat)
            read(arg, *) kchoose_end
        endif
        if (trim(arg) == '-nb') then 
            iarg = iarg + 1
            call get_command_argument(iarg, arg, lens, stat)
            read(arg, *) bot_band
            iarg = iarg + 1
            call get_command_argument(iarg, arg, lens, stat)
            read(arg, *) top_band
        endif 
        if (trim(arg) == '-tolE') then 
            iarg = iarg + 1
            call get_command_argument(iarg, arg, lens, stat)
            read(arg, *) tolE
        endif

        iarg = iarg + 1
    enddo 
    !
    allocate(rot(3,3,max_nsym))  ;  rot=0
    allocate(tau(3,max_nsym))     ;  tau=0.0_dp
    allocate(SU2(2,2,max_nsym)) ; SU2=cmplx(0.0_dp,0.0_dp,kind=dp)
    allocate(time_reversal(max_nsym)) ; time_reversal=0
    allocate(spin_rot(3,3,max_nsym)) ; spin_rot=0.0_dp
    !
    ssg_label = ''
    !
    call get_ssg_op(spg, ssg_label, isSpinor, num_sym, rot, tau, SU2, spin_rot, time_reversal)
    !
    call load_bilbao(spg)
    !
    call read_tbbox()
    !
    if(isSpinPola) then
        call read_HmnR_spin_polarized()
    else
        call read_HmnR()
    endif
    !
    call output_ssg_operator(ssg_label, isSpinor, num_sym, rot, tau, SU2, spin_rot, time_reversal)
    !
    call setarray()
    !
    allocate(rot_lg(3,3,num_sym))                      ; rot_lg=0
    allocate(SO3_lg(3,3,num_sym))                      ; SO3_lg=0.d0
    allocate(SU2_lg(2,2,num_sym))                      ; SU2_lg=0.d0 
    allocate(tau_lg(3,num_sym))                      ; tau_lg=0.d0 
    allocate(kphase(num_sym))                       ; kphase=0.d0
    allocate(litt_group(num_sym))                   ; litt_group=0        
    allocate(time_reversal_lg(num_sym))                ; time_reversal_lg=0  
    allocate(SO3_unitary_lg(3,3,num_sym))              ; SO3_unitary_lg=0.d0
    allocate(SU2_unitary_lg(2,2,num_sym))              ; SU2_unitary_lg=0.d0
    allocate(rot_unitary_lg(3,3,num_sym))         ; rot_unitary_lg=0
    allocate(tau_unitary_lg(3,num_sym))                ; tau_unitary_lg=0.d0
    allocate(litt_group_unitary(num_sym))            ; litt_group_unitary=0   
    allocate(spin_rot_lg(3,3,num_sym))                 ; spin_rot_lg=0.0_dp
    allocate(tilte_mat(num_sumorb,num_sumorb,num_sym)) ; tilte_mat = 0.d0

    !
    if(.not. kchoose_flag) then
        kchoose_start = 1
        kchoose_end = num_k
    endif
    !

    
    if(kchoose_end > num_k) then
        stop  'kpoints chosen is larger than the number of kpoints!!!'
    endif


    do kkk = kchoose_start, kchoose_end
        
        WK(:) = k(:,kkk)
        do i = 1, 3
            if (abs(3.d0*WK(i)-1.d0).lt.0.0002d0) WK(i) = 1.d0/3.d0
        enddo 
        
        call kgroup(num_sym, num_sym, rot, time_reversal, WK, litt_group, numLG)

        numULG = 0
        do i=1, numLG
            rot_lg(:,:,i) = rot(:,:,litt_group(i))
            SO3_lg(:,:,i) = matmul(matmul(br2, rot_lg(:,:,i)), br4)
            SU2_lg(:,:,i) = SU2(:,:,litt_group(i))
            spin_rot_lg(:,:,i) = spin_rot(:,:,litt_group(i))
            time_reversal_lg(i) = time_reversal(litt_group(i))
            tau_lg(:,i) = tau(:,litt_group(i))

            if (time_reversal_lg(i)==1) then
                numULG = numULG + 1
                litt_group_unitary(numULG) = litt_group(i)
                SO3_unitary_lg(:,:,numULG) = SO3_lg(:,:,i)
                SU2_unitary_lg(:,:,numULG) = SU2_lg(:,:,i)
                rot_unitary_lg(:,:,numULG) = rot_lg(:,:,i)
                tau_unitary_lg(:,numULG) = tau_lg(:,i)
            endif
            
        enddo

        
        isSpinor = .true.
        if (nspin==2) then
            call judge_up_down_relation(numLG,SU2_lg,time_reversal_lg,spin_no_reversal)
            if (.not. spin_no_reversal) then
                isSpinor = .true.
                call get_WF_spin_polarized(3)
            else
                isSpinor = .false.
                call get_WF_spin_polarized(1)
            endif
        else
            call get_WF()
        endif

        
        write(*,'(A)')'********************************************************************************'
        write(*,*)
        write(*,*)

        write(154,'(A)')'********************************************************************************'


        write(180,'(A)')'********************************************************************************'
        write(180,*)
        write(180,*)

        call get_k_name( WK ,k_name, k_frac_symbol)


        if (nspin==2) then
            if (spin_no_reversal) then
                write(*,'(A,1I3,A)')'knum = ',kkk,' spin = up'
                write(180,'(A,1I3,A)')'knum = ',kkk,' spin = up'
            else
                write(*,'(A)')'There is a SSG operator flipping spin in the k-little group.'
                write(180,'(A)')'There is a SSG operator flipping spin in the k-little group.'
                write(*,'(A,1I3,A)')'knum = ',kkk,' spin = up & down'
                write(180,'(A,1I3,A)')'knum = ',kkk,' spin = up & down'
            endif
        else
            write(*,'(A,1I3)')'knum = ',kkk
            write(180,'(A,1I3)')'knum = ',kkk
        endif

        write(*,'(A,3F9.6)')'K = ',WK
        write(*,'(A)')'kname = '//trim(adjustl(k_name))//'    Fractional coordinate: '//trim(adjustl(k_frac_symbol))//' (given in the conventional basis)'
        ! write(154,'(A,3F9.6)')'K = ',WK
        write(154,'(A)')'kname = '//trim(adjustl(k_name))//'    Fractional coordinate: '//trim(adjustl(k_frac_symbol))//' (given in the conventional basis)'
        write(180,'(A,3F9.6)')'K = ',WK
        write(180,'(A)')'kname = '//trim(adjustl(k_name))//'    Fractional coordinate: '//trim(adjustl(k_frac_symbol))//' (given in the conventional basis)'

        call tb_setup_ssg(WK, &
                numULG, &
                ncenter, pos_center, &
                num_sumorb, num_sumorb, norb_center, orbt, &
                rot_unitary_lg, SO3_unitary_lg, tau_unitary_lg, &
                kphase, tilte_mat)

        allocate(op_order(numLG))
        allocate(ch_table(numLG,numLG))
        allocate(ch_unitary_table(numLG,numLG))
        allocate(discriminant_value(numLG))
        allocate(torsion(numLG))
        allocate(phase(numLG))
 
        call get_ch_from_op(numLG, spin_rot_lg(:,:,1:numLG), rot_lg(:,:,1:numLG), &
                            tau_lg(:,1:numLG), SU2_lg(:,:,1:numLG), time_reversal_lg(1:numLG), &
                            WK, op_order, Ncoirrep, ch_table, phase, &
                            Nirrep, ch_unitary_table, torsion)

        allocate(Ch_table2(Ncoirrep,numULG))
        allocate(Ch_table1(Nirrep,numULG))
        allocate(irrep_coirrep_relation(2,Ncoirrep))
        allocate(irrep_name_list(Ncoirrep))


        Ch_table2(1:Ncoirrep,1:numULG) = ch_table(1:Ncoirrep,1:numULG)
        Ch_table1(1:Nirrep,1:numULG) = ch_unitary_table(1:Nirrep,1:numULG)

        deallocate(ch_table)
        deallocate(ch_unitary_table)
        
        

        call find_relation_between_unitary_irrep_coirrep(numULG,Ncoirrep,Ch_table2,Nirrep,Ch_table1,irrep_coirrep_relation,&
                                                        discriminant_value,torsion)

        call output_character_table(k_name,numLG,numULG,litt_group,op_order,Ncoirrep,Ch_table2,phase, &
                                Nirrep,Ch_table1,irrep_coirrep_relation,irrep_name_list,discriminant_value)


        call get_comprel(numLG,numULG,litt_group,op_order,Ncoirrep,irrep_name_list,Ch_table2)

        do i=1,Ncoirrep
            do j=1,numULG
                Ch_table2(i,j) = Ch_table2(i,j) * phase(j)
            enddo
        enddo


        if (nspin==1 .or. spin_no_reversal) then
            call irrep_ssg( numULG,litt_group_unitary,&
                            rot_unitary_lg, tau_unitary_lg, SO3_unitary_lg, SU2_unitary_lg, &
                            kkk, WK, kphase, &
                            num_bands, bot_band, top_band, EE, &
                            num_sumorb, num_sumorb, &
                            coeffa, coeffb, &
                            Ncoirrep, Ch_table2, &
                            irrep_name_list, &
                            rot_mat_tb=tilte_mat, tolE_=tolE)
        else
            call irrep_ssg( numULG,litt_group_unitary,&
                    rot_unitary_lg, tau_unitary_lg, SO3_unitary_lg, SU2_unitary_lg, &
                    kkk, WK, kphase, &
                    num_bands*2, bot_band*2-1, top_band*2, EE, &
                    num_wann, num_wann, &
                    coeffa, coeffb, &
                    Ncoirrep, Ch_table2, &
                    irrep_name_list, &
                    rot_mat_tb=tilte_mat, tolE_=tolE)
        endif

        deallocate(Ch_table2)
        deallocate(Ch_table1)
        deallocate(op_order)
        deallocate(irrep_coirrep_relation)
        deallocate(irrep_name_list)
        deallocate(discriminant_value)
        deallocate(torsion)
        deallocate(phase)

    enddo 

    call end_get_comprel()

    if (isSpinPola) then
    do kkk = kchoose_start, kchoose_end
        WK(:) = k(:,kkk)
        do i = 1, 3
            if (abs(3.d0*WK(i)-1.d0).lt.0.0002d0) WK(i) = 1.d0/3.d0
        enddo 
        
        call kgroup(num_sym, num_sym, rot, time_reversal, WK, litt_group, numLG)

        numULG = 0

        do i=1, numLG
            rot_lg(:,:,i) = rot(:,:,litt_group(i))
            SO3_lg(:,:,i) = matmul(matmul(br2, rot_lg(:,:,i)), br4)
            SU2_lg(:,:,i) = SU2(:,:,litt_group(i))
            spin_rot_lg(:,:,i) = spin_rot(:,:,litt_group(i))
            time_reversal_lg(i) = time_reversal(litt_group(i))
            tau_lg(:,i) = tau(:,litt_group(i))

            if (time_reversal_lg(i)==1) then
                numULG = numULG + 1
                litt_group_unitary(numULG) = litt_group(i)
                SO3_unitary_lg(:,:,numULG) = SO3_lg(:,:,i)
                SU2_unitary_lg(:,:,numULG) = SU2_lg(:,:,i)
                rot_unitary_lg(:,:,numULG) = rot_lg(:,:,i)
                tau_unitary_lg(:,numULG) = tau_lg(:,i)
            endif
            
        enddo

        call judge_up_down_relation(numLG,SU2_lg,time_reversal_lg,spin_no_reversal)

        if (.not. spin_no_reversal) then
            isSpinor = .true.
            call get_WF_spin_polarized(3)
        else
            isSpinor = .false.
            call get_WF_spin_polarized(2)
        endif


        write(*,'(A)')'********************************************************************************'
        write(*,*)
        write(*,*)

        write(180,'(A)')'********************************************************************************'
        write(180,*)
        write(180,*)

        call get_k_name( WK ,k_name, k_frac_symbol)
        
        if (spin_no_reversal) then
            write(*,'(A,1I3,A)')'knum = ',kkk,' spin = down'
            write(180,'(A,1I3,A)')'knum = ',kkk,' spin = down'
        else
            write(*,'(A)')'There is a SSG operator flipping spin in the k-little group.'
            write(*,'(A,1I3,A)')'knum = ',kkk,' spin = up & down'

            write(180,'(A)')'There is a SSG operator flipping spin in the k-little group.'
            write(180,'(A,1I3,A)')'knum = ',kkk,' spin = up & down'
        endif

        write(*,'(A,3F9.6)')'K = ',WK
        write(*,'(A)')'kname = '//trim(adjustl(k_name))//'    Fractional coordinate: '//trim(adjustl(k_frac_symbol))//' (given in the conventional basis)'
        write(180,'(A,3F9.6)')'K = ',WK
        write(180,'(A)')'kname = '//trim(adjustl(k_name))//'    Fractional coordinate: '//trim(adjustl(k_frac_symbol))//' (given in the conventional basis)'
        
        call tb_setup_ssg(WK, &
                numULG, &
                ncenter, pos_center, &
                num_sumorb, num_sumorb, norb_center, orbt, &
                rot_unitary_lg, SO3_unitary_lg, tau_unitary_lg, &
                kphase, tilte_mat)
                    
        allocate(op_order(numLG))
        allocate(ch_table(numLG,numLG))
        allocate(ch_unitary_table(numLG,numLG))
        allocate(discriminant_value(numLG))
        allocate(torsion(numLG))
        allocate(phase(numLG))

        call get_ch_from_op(numLG, spin_rot_lg(:,:,1:numLG), rot_lg(:,:,1:numLG), &
                            tau_lg(:,1:numLG), SU2_lg(:,:,1:numLG), time_reversal_lg(1:numLG), &
                            WK, op_order, Ncoirrep, ch_table, phase, Nirrep, ch_unitary_table, torsion)

        allocate(Ch_table2(Ncoirrep,numULG))
        allocate(Ch_table1(Nirrep,numULG))
        allocate(irrep_coirrep_relation(2,Ncoirrep))
        allocate(irrep_name_list(Ncoirrep))

        Ch_table2(1:Ncoirrep,1:numULG) = ch_table(1:Ncoirrep,1:numULG)
        Ch_table1(1:Nirrep,1:numULG) = ch_unitary_table(1:Nirrep,1:numULG)

        deallocate(ch_table)
        deallocate(ch_unitary_table)



        call find_relation_between_unitary_irrep_coirrep(numULG,Ncoirrep,Ch_table2,Nirrep,Ch_table1,irrep_coirrep_relation,&
                                                        discriminant_value,torsion)

        call output_character_table(k_name,numLG,numULG,litt_group,op_order,Ncoirrep,Ch_table2,phase,&
                                Nirrep,Ch_table1,irrep_coirrep_relation,irrep_name_list,discriminant_value)

        do i=1,Ncoirrep
            do j=1,numULG
                Ch_table2(i,j) = Ch_table2(i,j) * phase(j)
            enddo
        enddo

        if (spin_no_reversal) then
            call irrep_ssg( numULG,litt_group_unitary,&
                            rot_unitary_lg, tau_unitary_lg, SO3_unitary_lg, SU2_unitary_lg, &
                            kkk, WK, kphase, &
                            num_bands, bot_band, top_band, EE, &
                            num_wann, num_wann, &
                            coeffa, coeffb, &
                            Ncoirrep, Ch_table2, &
                            irrep_name_list, &
                            rot_mat_tb=tilte_mat, tolE_=tolE)
        else
            call irrep_ssg( numULG,litt_group_unitary,&
                    rot_unitary_lg, tau_unitary_lg, SO3_unitary_lg, SU2_unitary_lg, &
                    kkk, WK, kphase, &
                    num_bands*2, bot_band*2-1, top_band*2, EE, &
                    num_sumorb, num_sumorb, &
                    coeffa, coeffb, &
                    Ncoirrep, Ch_table2, &
                    irrep_name_list, &
                    rot_mat_tb=tilte_mat, tolE_=tolE)
        endif


        deallocate(Ch_table2)
        deallocate(Ch_table1)
        deallocate(op_order)
        deallocate(irrep_coirrep_relation)
        deallocate(irrep_name_list)
        deallocate(discriminant_value)
        deallocate(torsion)
        deallocate(phase)
    enddo 
    endif

    call downarray()

    deallocate(rot_lg)
    deallocate(SO3_lg)
    deallocate(SU2_lg)
    deallocate(tau_lg)
    deallocate(kphase)
    deallocate(tilte_mat)
    deallocate(litt_group)
    deallocate(time_reversal_lg)
    deallocate(rot_unitary_lg)
    deallocate(SO3_unitary_lg)
    deallocate(SU2_unitary_lg)
    deallocate(litt_group_unitary)
    deallocate(tau_unitary_lg)
    deallocate(spin_rot_lg)
    WRITE(6,*) 
    WRITE(6,*) "*****************************************************"
    WRITE(6,*) 
    WRITE(6,*) "TOTAL END"


end subroutine run_wann
end module driver_wann
