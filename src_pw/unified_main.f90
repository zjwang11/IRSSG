program irssg
    use driver_pw,   only: run_pw
    use driver_wann, only: run_wann
    implicit none

    integer :: narg, iarg, stat, lens
    character(len=256) :: arg
    logical :: use_wann

    use_wann = .false.
    narg = command_argument_count()
    iarg = 1
    do while (iarg <= narg)
        call get_command_argument(iarg, arg, lens, stat)
        if (len_trim(arg) == 0) exit
        if (trim(arg) == '--wann') then
            use_wann = .true.
        else if (trim(arg) == '--mode') then
            if (iarg + 1 <= narg) then
                call get_command_argument(iarg+1, arg, lens, stat)
                if (trim(arg) == 'wann') use_wann = .true.
            end if
        else if (index(arg, '--mode=') == 1) then
            if (trim(arg(8:)) == 'wann') use_wann = .true.
        end if
        iarg = iarg + 1
    end do

    if (use_wann) then
        call run_wann()
    else
        call run_pw()
    end if

end program irssg

