module init_wann

    use comms_wann 
    implicit none 
    public :: setarray
    public :: downarray
contains

subroutine setarray()

    integer :: irot

    if (.not. isSpinPola) then
        allocate (EE(num_bands))                         ; EE=0.d0
        allocate (coeffa(num_sumorb, num_wann))           ; coeffa=0.d0
        allocate (coeffb(num_sumorb, num_wann))           ; coeffb=0.d0
    else
        allocate (EE(num_bands*2))                         ; EE=0.d0
        allocate (coeffa(num_sumorb, num_wann*2))           ; coeffa=0.d0
        allocate (coeffb(num_sumorb, num_wann*2))           ; coeffb=0.d0
        allocate (order(num_bands*2))                      ; order = 0
    endif



end subroutine

subroutine downarray()
    if(allocated(EE    ))    deallocate (EE    )
    if(allocated(coeffa))    deallocate (coeffa)
    if(allocated(coeffb))    deallocate (coeffb)

end subroutine


subroutine read_HmnR_spin_polarized()

    integer :: i, j, ir 
    integer :: itmp, jtmp 

    real(dp) :: r1tmp, r2tmp

    integer, parameter :: hr = 1000


    open(unit=hr, file=hrfile_up, status='old')
    read(hr, *)
    read(hr, *) num_wann
    read(hr, *) num_block
    num_sumorb = num_wann
    if (bot_band == 0) bot_band = 1
    if (top_band == 0) top_band = num_wann

    allocate(deg_block(num_block));                deg_block = 0
    allocate(vec_block(3, num_block));             vec_block = 0
    allocate(HmnR_up(num_wann, num_wann, num_block)); HmnR_up = czero

    read(hr, "(15I5)") deg_block(:)
    do ir = 1, num_block
        do j = 1, num_wann
        do i = 1, num_wann
            read(hr, *) vec_block(:,ir), itmp, jtmp, r1tmp, r2tmp
            HmnR_up(i,j,ir) = cmplx(r1tmp, r2tmp, dp)
        enddo 
        enddo

        if (deg_block(ir) == 1) cycle
        HmnR_up(:,:,ir) = HmnR_up(:,:,ir)/deg_block(ir)
        deg_block(ir) = 1
    enddo 

    close(hr)

    deallocate(deg_block)
    deallocate(vec_block)

    open(unit=hr, file=hrfile_dn, status='old')
    read(hr, *)
    read(hr, *) num_wann
    num_bands = num_wann
    read(hr, *) num_block
    num_sumorb = num_wann
    if (bot_band == 0) bot_band = 1
    if (top_band == 0) top_band = num_wann
    allocate(deg_block(num_block));                deg_block = 0
    allocate(vec_block(3, num_block));             vec_block = 0
    allocate(HmnR_dn(num_wann, num_wann, num_block)); HmnR_dn = czero

    read(hr, "(15I5)") deg_block(:)
    do ir = 1, num_block
        do j = 1, num_wann
        do i = 1, num_wann
            read(hr, *) vec_block(:,ir), itmp, jtmp, r1tmp, r2tmp
            HmnR_dn(i,j,ir) = cmplx(r1tmp, r2tmp, dp)
        enddo 
        enddo

        if (deg_block(ir) == 1) cycle
        HmnR_dn(:,:,ir) = HmnR_dn(:,:,ir)/deg_block(ir)
        deg_block(ir) = 1
    enddo 

    close(hr)
end subroutine read_HmnR_spin_polarized


subroutine read_HmnR()

    integer :: i, j, ir 
    integer :: itmp, jtmp 

    real(dp) :: r1tmp, r2tmp

    integer, parameter :: hr = 1000

    open(unit=hr, file=hrfile, status='old')
    read(hr, *)
    read(hr, *) num_wann
    read(hr, *) num_block
    num_bands = num_wann
    num_sumorb = num_wann/2
    if (bot_band == 0) bot_band = 1
    if (top_band == 0) top_band = num_wann

    allocate(deg_block(num_block));                deg_block = 0
    allocate(vec_block(3, num_block));             vec_block = 0
    allocate(HmnR(num_wann, num_wann, num_block)); HmnR = czero
    
    read(hr, "(15I5)") deg_block(:)
    do ir = 1, num_block
        do j = 1, num_wann
        do i = 1, num_wann
            read(hr, *) vec_block(:,ir), itmp, jtmp, r1tmp, r2tmp
            HmnR(i,j,ir) = cmplx(r1tmp, r2tmp, dp)
        enddo 
        enddo

        if (spincov == 2) then
            call reorder_oddeven(HmnR(1:num_wann,1:num_wann,ir),num_wann)
        endif
        
        if (deg_block(ir) == 1) cycle
        HmnR(:,:,ir) = HmnR(:,:,ir)/deg_block(ir)
        deg_block(ir) = 1
        
    enddo 

    close(hr)
    
end subroutine read_HmnR


subroutine read_tbbox()

    character(len=22)  :: chtp22

    integer            :: i, j, i1, i2, irot 
    integer            :: ierr 
    
    integer            :: kmesh, Nk 
    real(dp), allocatable, save :: node_kpoints(:,:)
    
    integer            :: itmp
    real(dp)           :: rtmp(8), rotmt(3,3)
    real(dp)           :: tmp(3), tmpr(3)
    complex(dp)        :: crotmt(3,3), protmt(3,3), drotmt(5,5), srotmt(2,2), frotmt(7,7)
    complex(dp)        :: cmat3(3,3), cmat5(5,5), cmat7(7,7)
    real(dp) :: dummy1, dummy2, dummy3

    integer, parameter :: tbbox = 1001


    open(unit=tbbox, file='tbbox.in',form='formatted',status='old')

    isSpinPola = .false.
    chtp22 = 'spinpol'
    call get_key_para_cht(chtp22, tbbox, casename)
    casename = trim(adjustl(casename))
    if (casename(1:1) == 't' .or. casename(1:1) == 'T') then 
        isSpinPola = .true.
        nspin = 2
    else
        isSpinPola = .false.
        nspin = 1
    endif

    if (isSpinPola) then 
        hrfile_up = 'wannier90.up_hr.dat'
        chtp22 = 'hr_up_name'
        call get_key_para_cht(chtp22, tbbox, casename)
        hrfile_up = trim(casename)


        hrfile_dn = 'wannier90.dn_hr.dat'
        chtp22 = 'hr_dn_name'
        call get_key_para_cht(chtp22, tbbox, casename)
        hrfile_dn = trim(casename)
    else
        hrfile = 'wannier90_hr.dat'
        chtp22 = 'hr_name'
        call get_key_para_cht(chtp22, tbbox, casename)
        hrfile = trim(adjustl(casename))
    endif

    chtp22 = 'proj'
    call get_key_para_loc(chtp22, tbbox)

    chtp22 = 'orbt'
    orbt = -120
    call get_key_para_intct(chtp22, tbbox, orbt)
    
    if (orbt == 1 .or. orbt == -120) then 
        write(6, *) "Orbital convention 1: &
                     s, px, py, pz, xy, yz, zx, x2-y2, 3z2-r2 !!!"
    elseif (orbt == 2) then 
        write(6, *) "Orbital convention 2: &
                     s, pz, px, py, 3z2-r2, xz, yz, x2-y2, xy(Wannier90) !!!"
    else
        write(*,*) "Error: orbt !!!!"
        stop 
    endif 

    chtp22 = 'spincov'
    spincov = -120
    if (orbt == -120) then
        call get_key_para_int(chtp22, tbbox, spincov)
    else
        call get_key_para_intct(chtp22, tbbox, spincov)
    endif
    if (spincov == 1 .or. spincov == -120) then 
        write(6, *) "Spin convention 1: &
                     up up .../ dn dn ... (Wannier90.1.x) !!!"
    elseif (spincov == 2) then 
        write(6, *) "Spin convention 2: &
                     up dn/ up dn/ ... (Wannier90.3.x) !!!"
    else
        write(*,*) "Error: spincov !!!!"
        stop 
    endif 

    chtp22 = 'ntau'

    if (spincov == -120) then
        call get_key_para_int(chtp22, tbbox, ncenter) 
    else
        call get_key_para_intct(chtp22, tbbox, ncenter) 
    endif

    allocate(type_center(ncenter), norb_center(ncenter), startorb_center(ncenter))
    allocate(pos_center(3, ncenter))
    startorb_center = 0
    do i = 1, ncenter 
        read(tbbox, *) pos_center(:, i), dummy1, dummy2, dummy3, type_center(i), norb_center(i)
        if (i /= 1) startorb_center(i) = startorb_center(i-1) + norb_center(i-1)
    enddo 

    ! isInv         = .false.
    ! isSymmorphic  = .false. 
    chtp22 = 'unit_cell'
    call get_key_para_loc(chtp22, tbbox)
    do i = 1, 3
        read(tbbox, *) br2(:,i)
    enddo 
    lattice = br2
    call invreal33(br2,br4)


    ! kpath 
    chtp22 = 'kpoint'
    call get_key_para_loc(chtp22, tbbox)
    chtp22 = 'kmesh'
    call get_key_para_int(chtp22, tbbox, kmesh)
    chtp22 = 'Nk'
    call get_key_para_int(chtp22, tbbox, Nk)
    allocate(node_kpoints(3,Nk)); node_kpoints = 0.d0
    call get_key_para_nvec(tbbox, Nk, node_kpoints)

    allocate(len_k((Nk-1)*kmesh+1)); len_k = 0.d0 
    allocate(k(3,(Nk-1)*kmesh+1)); k = 0.d0 
    do i = 1, Nk - 1
        tmp = (node_kpoints(:,i+1)-node_kpoints(:,i))/kmesh
        do j = 1, kmesh
            k(:,(i-1)*kmesh+j) = node_kpoints(:,i)+tmp*(j-1)
            tmpr = matmul(tmp(:), br4(:,:))
            if (i/=1 .or. j/=1) then 
                len_k((i-1)*kmesh+j) = len_k((i-1)*kmesh+j-1) + &
                                  dsqrt(tmpr(1)**2+tmpr(2)**2+tmpr(3)**2)
            endif 
        enddo 
    enddo 
    k(:,(Nk-1)*kmesh+1) = node_kpoints(:,Nk)
    len_k = 2.d0*PI*len_k 
    num_k = (Nk-1)*kmesh + 1

    close(tbbox)    


    ! isComplexWF = .true. 
    ! if (.not.isSpinor .and. isInv) isComplexWF = .false. 
    ! if (     isSymmorphic) write(6,'(A19)',advance='NO') ' Symmorphic crystal'
    ! if (.not.isSymmorphic) write(6,'(A23)',advance='NO') ' Non-symmorphic crystal'

    ! if (     isInv) write(6,'(A24)') ' with inversion symmetry'
    ! if (.not.isInv) write(6,'(A27)') ' without inversion symmetry'

    ! if (     isComplexWF) write(6, '(A23)') ' Complex eigenfunctions'
    ! if (.not.isComplexWF) write(6, '(A20)') ' Real eigenfunctions'

    ! if (     isSpinor) write(6, '(A45)') ' Spin-orbit eigenfunctions (->time inversion)'
    ! if (.not.isSpinor) write(6, '(A29)') ' No spin-orbit eigenfunctions'

      WRITE(6,590)
      WRITE(6,592)
      WRITE(6,595) ((br2(I1,I2),I2=1,3),I1=1,3)
      WRITE(6,591)
      WRITE(6,596)  (br4( 1,I2),I2=1,3),' : g1/2pi'
      WRITE(6,596)  (br4( 2,I2),I2=1,3),' : g2/2pi'
      WRITE(6,596)  (br4( 3,I2),I2=1,3),' : g3/2pi'
   
      RETURN
!529  FORMAT(1X,A10,/,1X,A4,' lattice')
 529  FORMAT(1X,A10,/)
 590  FORMAT(//,' Transformations:',/, &
             ' Direct lattice vectors in Cartesian coord. system (BR2)')
 591  FORMAT(' Reciprocal lattice vectors in Cartesian coord. system (BR4)')
 592  FORMAT('        t1 :            t2 :            t3 : ')
 595  FORMAT(3(3F16.8,/))
 596  FORMAT((3F16.8),A9)
end subroutine read_tbbox

subroutine reorder_oddeven(H,dim)
  implicit none
  integer,intent(in) :: dim
  complex(dp), intent(inout) :: H(dim,dim)
  integer :: n, i, p
  integer, allocatable :: perm(:)
  complex(dp), allocatable :: tmp(:,:)

  n = dim

  allocate(perm(n))
  p = 0
  do i = 1, n, 2     
    p = p + 1
    perm(p) = i
  end do
  do i = 2, n, 2    
    p = p + 1
    perm(p) = i
  end do

  allocate(tmp(n,n))
  tmp = H(perm, perm)  
  H = tmp

  deallocate(tmp, perm)
end subroutine reorder_oddeven
end module init_wann
