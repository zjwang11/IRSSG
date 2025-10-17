module comprel
use lib_params, only:dp
integer, save :: num_litt_group_old
integer, allocatable, save :: op_order_old(:)
integer, save :: num_litt_group_unitary_old
character(len=15), allocatable, save :: irrep_name_old(:)
integer, save :: irrep_num_old
complex(dp), allocatable, save :: ch_table_old(:,:)
integer, allocatable, save :: litt_group_old(:)
public :: get_comprel
public :: end_get_comprel
contains



subroutine get_comprel(num_litt_group,num_litt_group_unitary,litt_group,op_order,irrep_num,irrep_name,ch_table)
    implicit none
    integer , intent(in) :: num_litt_group
    integer, intent(in) :: num_litt_group_unitary
    integer, intent(in) :: op_order(num_litt_group)
    integer :: relation
    integer, intent(in) :: irrep_num
    character(len=15), intent(in) :: irrep_name(irrep_num)
    complex(dp), intent(in) :: ch_table(irrep_num,num_litt_group_unitary)
    integer, intent(in) :: litt_group(num_litt_group)
    integer :: op_order_(num_litt_group)
    integer :: op_order_old_(num_litt_group_old)
    integer :: i
    if (.not. allocated(op_order_old)) then
        open(unit=202, file="comprel.log", status="replace", action="write")
        close(202)
        allocate(op_order_old(num_litt_group))
        op_order_old = op_order
        num_litt_group_old = num_litt_group
        allocate(irrep_name_old(irrep_num))
        irrep_name_old = irrep_name
        irrep_num_old = irrep_num
        num_litt_group_unitary_old = num_litt_group
        allocate(ch_table_old(num_litt_group,num_litt_group_unitary))
        ch_table_old = ch_table

        allocate(litt_group_old(num_litt_group))
        litt_group_old = litt_group

        return
    else
        open(unit=202, file="comprel.log", status="unknown", action="write",access="sequential", position="append")
    end if

    do i = 1, num_litt_group
        op_order_(i) = litt_group(op_order(i))
    end do

    do i = 1, num_litt_group_old
        op_order_old_(i) = litt_group_old(op_order_old(i))
    end do



    call judge_little_group_relation(num_litt_group_old,num_litt_group,op_order_old_,op_order_,relation)
    if (relation == 2) then
        return
        ! write(202,*) 'Compatibility relations: '
        do i = 1, irrep_num
            write(202,'(A)') trim(adjustl(irrep_name(i)))//'->'//trim(adjustl(irrep_name_old(i)))
        end do
    elseif (relation == 1) then
        ! write(202,*) 'Compatibility relations: '
        write(202,'(A)')'========================'
        call decomposition_rep(num_litt_group_old,num_litt_group,num_litt_group_unitary_old,num_litt_group_unitary,op_order_old_,op_order_,irrep_num_old,irrep_num,irrep_name_old,irrep_name,ch_table_old,ch_table)
    elseif (relation == -1) then
        ! write(202,*) 'Compatibility relations: '
        write(202,'(A)')'========================'
        call decomposition_rep(num_litt_group,num_litt_group_old,num_litt_group_unitary,num_litt_group_unitary_old,op_order_,op_order_old_,irrep_num,irrep_num_old,irrep_name,irrep_name_old,ch_table,ch_table_old)
    endif


    deallocate(irrep_name_old)
    allocate(irrep_name_old(irrep_num))
    irrep_num_old = irrep_num
    irrep_name_old = irrep_name

    deallocate(op_order_old)
    allocate(op_order_old(num_litt_group))
    op_order_old = op_order
    num_litt_group_old = num_litt_group
    num_litt_group_unitary_old = num_litt_group_unitary

    deallocate(ch_table_old)
    allocate(ch_table_old(num_litt_group,num_litt_group_unitary))
    ch_table_old = ch_table

    deallocate(litt_group_old)
    allocate(litt_group_old(num_litt_group))
    litt_group_old = litt_group
    close(202)
endsubroutine


subroutine end_get_comprel()
    implicit none
    if (allocated(op_order_old)) then
        deallocate(op_order_old)
    end if
    if (allocated(irrep_name_old)) then
        deallocate(irrep_name_old)
    end if
    if (allocated(ch_table_old)) then
        deallocate(ch_table_old)
    end if
    if (allocated(litt_group_old)) then
        deallocate(litt_group_old)
    end if
end subroutine end_get_comprel

subroutine decomposition_rep(num_litt_group1,num_litt_group2,num_litt_group_unitary_1,num_litt_group_unitary_2,op_order1,op_order2,irrep_num1,irrep_num2,irrep_name1,irrep_name2,ch_table1,ch_table2)
    implicit none
    integer, intent(in) :: num_litt_group1, num_litt_group2
    integer, intent(in) :: num_litt_group_unitary_1, num_litt_group_unitary_2
    integer, intent(in) :: op_order1(num_litt_group1), op_order2(num_litt_group2)
    integer, intent(in) :: irrep_num1, irrep_num2
    character(len=15), intent(in) :: irrep_name1(irrep_num1), irrep_name2(irrep_num2)
    complex(dp), intent(in) :: ch_table1(irrep_num1,num_litt_group_unitary_1), ch_table2(irrep_num2,num_litt_group_unitary_2)
    integer :: i, j, k, j1, j2, j3, j4, j5, j6, l
    logical :: flag

    do i=1, irrep_num2
        do j=1, irrep_num1
            flag = .true.
            do k=1, num_litt_group_unitary_1
                if (abs(ch_table1(j,k) - ch_table2(i,find_index(op_order1(k),num_litt_group2,op_order2)))>1e-2) then
                    flag = .false.
                    exit
                end if
            end do
            if (flag) then
                exit
            endif
        end do
        if (flag) then
            write(202,'(A)') trim(adjustl(irrep_name2(i)))//'->'//trim(adjustl(irrep_name1(j)))
        else
            do j1=1, irrep_num1
                do j2 = 1, irrep_num1
                    flag = .true.
                    do k=1, num_litt_group_unitary_1
                        l = find_index(op_order1(k),num_litt_group2,op_order2)
                        if (abs(ch_table1(j1,k) + ch_table1(j2,k) - ch_table2(i,l))>1e-2) then
                            flag = .false.
                            exit
                        end if
                    end do
                    if (flag) then
                        exit
                    endif
                enddo
                if (flag) then
                    exit
                endif
            end do

            if (flag) then
                write(202,'(A)') trim(adjustl(irrep_name2(i)))//'->'//trim(adjustl(irrep_name1(j1)))//'+'//trim(adjustl(irrep_name1(j2)))
            else
                do j1=1, irrep_num1
                    do j2 = 1, irrep_num1
                        do j3 = 1, irrep_num1
                            flag = .true.
                            do k=1, num_litt_group_unitary_1
                                l = find_index(op_order1(k),num_litt_group2,op_order2)
                                if (abs(ch_table1(j1,k) + ch_table1(j2,k) + ch_table1(j3,k) - ch_table2(i,l))>1e-2) then
                                    flag = .false.
                                    exit
                                end if
                            end do
                            if (flag) then
                                exit
                            endif
                        enddo
                        if (flag) then
                            exit
                        endif
                    enddo
                    if (flag) then
                        exit
                    endif
                end do
                if (flag) then
                    write(202,'(A)') trim(adjustl(irrep_name2(i)))//'->'//trim(adjustl(irrep_name1(j1)))//'+'//trim(adjustl(irrep_name1(j2)))//&
                                    '+'//trim(adjustl(irrep_name1(j3)))
                else
                    do j1=1, irrep_num1
                        do j2 = 1, irrep_num1
                            do j3 = 1, irrep_num1
                                do j4 = 1, irrep_num1
                                    flag = .true.
                                    do k=1, num_litt_group_unitary_1
                                        l = find_index(op_order1(k),num_litt_group2,op_order2)
                                        if (abs(ch_table1(j1,k) + ch_table1(j2,k) + ch_table1(j3,k) + ch_table1(j4,k) - ch_table2(i,l))>1e-2) then
                                            flag = .false.
                                            exit
                                        end if
                                    end do
                                    if (flag) then
                                        exit    
                                    endif
                                enddo
                                if (flag) then
                                    exit
                                endif
                            enddo
                            if (flag) then
                                exit
                            endif
                        enddo
                        if (flag) then
                            exit
                        endif
                    end do
                    if (flag) then
                        write(202,'(A)') trim(adjustl(irrep_name2(i)))//'->'//trim(adjustl(irrep_name1(j1)))//'+'//trim(adjustl(irrep_name1(j2)))//&
                                            '+'//trim(adjustl(irrep_name1(j3)))//'+'//trim(adjustl(irrep_name1(j4)))
                    else
                        do j1=1, irrep_num1
                            do j2 = 1, irrep_num1
                                do j3 = 1, irrep_num1
                                    do j4 = 1, irrep_num1
                                        do j5 = 1, irrep_num1
                                            flag = .true.
                                            do k=1, num_litt_group_unitary_1
                                                l = find_index(op_order1(k),num_litt_group2,op_order2)
                                                if (abs(ch_table1(j1,k) + ch_table1(j2,k) + ch_table1(j3,k) + ch_table1(j4,k) + ch_table1(j5,k) - ch_table2(i,l))>1e-2) then
                                                    flag = .false.
                                                    exit
                                                end if
                                            end do
                                            if (flag) then
                                                exit
                                            endif
                                        enddo
                                        if (flag) then
                                            exit
                                        endif
                                    enddo
                                    if (flag) then
                                        exit
                                    endif
                                enddo
                                if (flag) then
                                    exit
                                endif
                            enddo
                            if (flag) then
                                exit
                            endif
                        end do
                        if (flag) then
                            write(202,'(A)') trim(adjustl(irrep_name2(i)))//'->'//trim(adjustl(irrep_name1(j1)))//'+'//trim(adjustl(irrep_name1(j2)))//&
                                                '+'//trim(adjustl(irrep_name1(j3)))//'+'//trim(adjustl(irrep_name1(j4)))//'+'//trim(adjustl(irrep_name1(j5)))
                        else
                            do j1=1, irrep_num1
                                do j2 = 1, irrep_num1
                                    do j3 = 1, irrep_num1
                                        do j4 = 1, irrep_num1
                                            do j5 = 1, irrep_num1
                                                do j6 = 1, irrep_num1
                                                    flag = .true.
                                                    do k=1, num_litt_group_unitary_1
                                                        l = find_index(op_order1(k),num_litt_group2,op_order2)
                                                        if (abs(ch_table1(j1,k) + ch_table1(j2,k) + ch_table1(j3,k) + ch_table1(j4,k) + ch_table1(j5,k) + ch_table1(j6,k) - ch_table2(i,l))>1e-2) then
                                                            flag = .false.
                                                            exit
                                                        end if
                                                    end do  
                                                    if (flag) then
                                                        exit
                                                    endif
                                                end do
                                                if (flag) then
                                                    exit
                                                endif
                                            end do
                                            if (flag) then
                                                exit
                                            endif
                                        end do  
                                        if (flag) then
                                            exit
                                        endif
                                    end do
                                    if (flag) then
                                        exit
                                    endif
                                end do
                                if (flag) then
                                    exit
                                endif
                            end do
                            if (flag) then
                                write(202,'(A)') trim(adjustl(irrep_name2(i)))//'->'//trim(adjustl(irrep_name1(j1)))//'+'//trim(adjustl(irrep_name1(j2)))//&
                                                    '+'//trim(adjustl(irrep_name1(j3)))//'+'//trim(adjustl(irrep_name1(j4)))//'+'//trim(adjustl(irrep_name1(j5)))//'+'//trim(adjustl(irrep_name1(j6)))
                            else
                                continue
                            endif
                        endif
                    endif
                endif
            endif
            
        endif
    end do


endsubroutine

subroutine judge_little_group_relation(num_litt_group1,num_litt_group2,op_order1,op_order2,relation)
    implicit none
    integer, intent(in) :: num_litt_group1, num_litt_group2
    integer, intent(in) :: op_order1(num_litt_group1), op_order2(num_litt_group2)
    integer, intent(out) :: relation
    integer :: i, j, k

    logical :: flag

    relation = 0

    if (mod(num_litt_group1,num_litt_group2) /= 0 .and. (mod(num_litt_group2,num_litt_group1) /= 0)) then
        relation = 0
        return
    end if

    if (num_litt_group1 == num_litt_group2) then
        if (all(op_order1 == op_order2)) then
            relation = 2
            return
        end if
    end if
    
    if (num_litt_group1 < num_litt_group2) then
        do i = 1, num_litt_group1
            flag = .false.
            do j = 1, num_litt_group2
                if (op_order1(i) == op_order2(j)) then
                    flag = .true.
                    exit
                end if
            end do
            if (.not. flag) then
                relation = 0
                return
            end if
        end do
        relation = 1
        return
    else
        do i = 1, num_litt_group2
            flag = .false.
            do j = 1, num_litt_group1
                if (op_order2(i) == op_order1(j)) then
                    flag = .true.
                    exit
                end if
            end do
            if (.not. flag) then
                relation = 0
                return
            end if
        end do
        relation = -1
        return
    end if


endsubroutine judge_little_group_relation



    integer function find_index(val, dim, arr) result(pos)
        implicit none
        integer,           intent(in) :: val
        integer,           intent(in) :: dim
        integer, intent(in) :: arr(dim)
        integer :: i
        pos = 0                              ! 未找到默认 0
        do i = 1, dim
            if (arr(i) == val) then
                pos = i
                return
            end if
        end do
    end function find_index
endmodule comprel