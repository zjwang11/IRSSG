! irssg
! shengzhang221@mails.ucas.ac.cn

subroutine irrep_ssg(num_litt_group_unitary, litt_group_unitary, rot, tau,&
                      SO3, SU2, &
                      KKK, WK, kphase, &
                      num_bands, m, n, ene_bands, &
                      dim_basis, num_basis, &
                      coeffa, coeffb, irrep_num, ch_table, &
                      irrep_name_list, &
                      G_phase_pw, rot_vec_pw, rot_mat_tb, tolE_)

    use lib_comms
    use chrct
    implicit none 
    intrinsic :: nint

    ! number of space-group operations (module the integer lattice translations) 
    ! should not be changed for different k
    integer,     intent(in) :: num_litt_group_unitary
    integer,     intent(in) :: litt_group_unitary(num_litt_group_unitary)

    ! the rotation part of space-group operations with respect to primitive lattice vectors
    ! should not be changed for different k
    integer,     intent(in) :: rot(3,3,num_litt_group_unitary)

    ! the translation part of space-group operations with respect to primitive lattice vectors
    ! should not be changed for different k
    real(dp),    intent(in) :: tau(3,num_litt_group_unitary)

    ! the rotation part of space-group operations given in Cartesian coordinates
    ! should not be changed for different k
    real(dp),    intent(in) :: SO3(3,3,num_litt_group_unitary)

    ! the rotation part of space-group operations given in spin-1/2 space
    ! should not be changed for different k
    complex(dp), intent(in) :: SU2(2,2,num_litt_group_unitary)

    ! the sequential number of the given k-point
    integer,     intent(in) :: KKK

    ! the coordinate of the k-point with respect to primitive reciprocal lattice vectors
    real(dp),    intent(in) :: WK(3)

    ! the k-dependent phase factors due to the translation part of space-group operations
    complex(dp), intent(in) :: kphase(num_litt_group_unitary)

    ! the total number of bands
    integer,     intent(in) :: num_bands

    ! the IRs of the set of bands [m,n] are computed
    integer,     intent(in) :: m, n 

    ! the energy of bands at the k-point
    real(dp),    intent(in) :: ene_bands(num_bands) 

    ! the reserved number of the PW/TB basis
    ! if rot_vec_pw is given, dim_basis should be larger than the PW number of any k-points
    ! if rot_mat_tb is given, one should set dim_basis = num_basis
    integer,     intent(in) :: dim_basis

    ! the number of PW or orthogonal TB basis for the k-point 
    ! (note: the PW basis numbers for different k-points are usually different)
    integer,     intent(in) :: num_basis  

    ! the coefficient of spinor up part of wave functions at the given k-point 
    ! (note: coeffup_basis(1:num_basis,1:num_bands) is nonzero)
    complex(dp), intent(in) :: coeffa(dim_basis, num_bands)

    ! the coefficent of spinor down part of wave functions at the given k-point if spinor is .true.
    ! (note: coeffdn_basis(1:num_basis,1:num_bands) is nonzero)
    complex(dp), intent(in) :: coeffb(dim_basis, num_bands)
 
    ! the phase factor dependent by the PW vectors
    complex(dp), intent(in), optional :: G_phase_pw(dim_basis, num_litt_group_unitary)

    ! the transformation vectors of space-group operations, which send the jth PW to the j'th PW
    integer,     intent(in), optional :: rot_vec_pw(dim_basis, num_litt_group_unitary)

    ! the transformation matrices of space-group operations in orthogonal TB basis
    complex(dp), intent(in), optional :: rot_mat_tb(dim_basis, dim_basis, num_litt_group_unitary)

    integer, intent(in) :: irrep_num

    complex(dp), intent(in) :: ch_table(irrep_num,num_litt_group_unitary)

    character(len=15), intent(in) :: irrep_name_list(irrep_num)

    real(dp), intent(in), optional :: tolE_
    real(dp) :: tolE
    
    ! calculated character tables, without outer k phase
    complex(dp)              :: chrct_set(num_bands, MAXSYM)

    ! degenerate informations of input bands
    integer                  :: deg_set(num_bands)
    integer                  :: numdeg_set(num_bands)

    ! representation names
    character(len=20)        :: reps_set(num_bands)
    integer                  :: numreps_set(num_bands)

    integer  :: i, j

    if (present(tolE_)) then 
        tolE = tolE_
    else
        tolE = 1.E-3_dp
    endif

    if (.not. allocated(save_chrct)) then
        save_kcount = 0
        allocate(save_chrct(num_bands, MAXKSAVE))
        allocate(save_numdeg(num_bands, MAXKSAVE))
        allocate(save_numrep(num_bands, MAXKSAVE))
        allocate(save_ktype(MAXKSAVE))
        save_chrct = -7
        save_numdeg = 0
        save_numrep = 0
        save_ktype = 0
    endif

    save_kcount = save_kcount + 1
    ! write(*,*) SU2,coeffa,coeffb,ene_bands,G_phase_pw,rot_vec_pw,chrct_set

    if (present(rot_mat_tb)) then
        call irssg_chrct(num_litt_group_unitary, SU2, &
                            num_bands, m, n, ene_bands, &
                            dim_basis, num_basis, &
                            coeffa, coeffb, &
                            chrct_set, deg_set, numdeg_set, &
                            map_mat=rot_mat_tb, tolE_=tolE) 
    elseif (present(rot_vec_pw)) then
        call irssg_chrct(num_litt_group_unitary, SU2, &
                            num_bands, m, n, ene_bands, &
                            dim_basis, num_basis, &
                            coeffa, coeffb, &
                            chrct_set, deg_set, numdeg_set, &
                            map_vec=rot_vec_pw, G_phase_pw=G_phase_pw, tolE_=tolE) 
    endif
    
    call irssg_reps(num_litt_group_unitary,litt_group_unitary,num_bands,m, n, ene_bands,&
                            irrep_num, ch_table, irrep_name_list, &
                            chrct_set, deg_set, numdeg_set, kphase, numreps_set)

end subroutine irrep_ssg


subroutine pw_setup_ssg(WK, num_litt_group_unitary,&
                    dim_basis, num_basis, Gvec, &
                    rot_unitary, tau_unitary, &
                    kphase, Gphase_pw, rot_vec_pw)

    use lib_comms
    implicit none 

    ! the coordinates of the k-point with respect to primitive reciprocal lattice vectors
    real(dp),    intent(in)  :: WK(3)

    ! the number of space-group operations 
    integer,     intent(in)  :: num_litt_group_unitary
    
    ! the translation part of space-group operations with respect to primitive lattice vectors
    real(dp),    intent(in) :: tau_unitary(3,num_litt_group_unitary)

    ! the reserved number of the PW_basis (dim_basis >= num_basis)
    integer,     intent(in)  :: dim_basis 

    ! the number of PW for the k-point (note: num_basis for different k-points are usually different)
    integer,     intent(in)  :: num_basis 

    ! the plane-wave G-vector with respected to reciprocal lattice vectors
    integer,     intent(in)  :: Gvec(3, dim_basis)

    ! the rotation part of space-group operations with respect to primitive lattice vectors
    integer,     intent(in) :: rot_unitary(3,3,num_litt_group_unitary)

    ! the k-dependent phase factors due to the translation part of space-group operations
    complex(dp), intent(out) :: kphase(num_litt_group_unitary)

    ! the phase factor dependent by the PW vectors
    complex(dp), intent(out) :: Gphase_pw(dim_basis, num_litt_group_unitary)

    ! the transformation vectors of Rs, which send the jth PW to the j'th PW
    integer,     intent(out) :: rot_vec_pw(dim_basis, num_litt_group_unitary)




    ! parameter used inside the library 
    integer :: i, j, irot

    
    complex(dp) :: so3tmp(3,3), su2tmp(2,2)


    integer  :: ind, ikv 
    real(dp) :: RKV(3)
    real(dp) :: diff, ang
    real(dp) :: KV(3,dim_basis)

    integer  :: invrot(3,3,num_litt_group_unitary)

    real(dp) :: A(3,3)

    do irot = 1, num_litt_group_unitary
        call invmati(rot_unitary(:,:,irot), invrot(:,:,irot)) 
    enddo 
    
    ! get kphase
    kphase = 0.d0 
    do irot = 1, num_litt_group_unitary
        ang = -2.d0*PI*( WK(1)*tau_unitary(1,irot) &
                          +WK(2)*tau_unitary(2,irot) &
                          +WK(3)*tau_unitary(3,irot))
        kphase(irot) = cmplx(dcos(ang), dsin(ang))
    enddo 
    
    ! get Gphase and rot_vec_pw
    KV = 0.d0 
    do ikv = 1, num_basis 
        KV(:,ikv) = dble(Gvec(:,ikv)) + WK(:)
    enddo 


    Gphase_pw = 0.d0
    rot_vec_pw = 0
    do irot = 1, num_litt_group_unitary
        do ikv = 1, num_basis 
            
            do j = 1, 3
                RKV(j) = KV(1,ikv)*dble(invrot(1,j,irot)) &
                        +KV(2,ikv)*dble(invrot(2,j,irot)) &
                        +KV(3,ikv)*dble(invrot(3,j,irot))

            enddo 

            ind = 1
            do while ( (abs(RKV(1)-KV(1,ind)) &
                       +abs(RKV(2)-KV(2,ind)) &
                       +abs(RKV(3)-KV(3,ind))) > 1.d-3 .and. (ind <= num_basis) )
                ind = ind + 1
            enddo 

            if (ind == num_basis + 1) then 
                write(6,*) "cannot find (k+G)inv(Ri)"
                stop
            endif 

            ang = 0.d0 
            do j = 1, 3
                ang = ang - 2.d0*PI*tau_unitary(j, irot)*(KV(j,ind)-WK(j))!!!!
            enddo 

            Gphase_pw(ikv, irot)  = cmplx(dcos(ang), dsin(ang))
            rot_vec_pw(ikv, irot) = ind 

        enddo 
    enddo     


end subroutine pw_setup_ssg



! setup for Tight-binding part of irrep_bcs
! prepare kphase, rot_mat_tb
subroutine tb_setup_ssg(WK, &
                    num_litt_group_unitary,&
                    num_atom, atom_position, &
                    dim_basis, num_basis, angularmom, orbt, &
                    rot_unitary, SO3_unitary, tau_unitary, &
                    kphase, rot_mat_tb)

    use lib_comms
    use mathlib, only:O3_to_axis_angle, Dmatrix
    implicit none 

    ! the coordinates of the k-point with respected to primitive reciprocal lattice vectors
    real(dp),    intent(in)  :: WK(3)

    ! the number of space-group operations 
    integer,     intent(in)  :: num_litt_group_unitary

    ! the determination of rotation part of space-group operations
    real(dp)  :: det(num_litt_group_unitary)
    
    ! the rotation angles of space-group operations
    real(dp)  :: angle(num_litt_group_unitary)

    ! the rotation axis of space-group operations in Cartesian coordinates
    real(dp)  :: axis(3, num_litt_group_unitary)

    ! the translation part of space-group operations with respect to primitive lattice vectors
    real(dp),    intent(in)  :: tau_unitary(3,num_litt_group_unitary)

    ! the number of atoms in the TB Hamiltonian
    integer,     intent(in)  :: num_atom

    ! the coordinates of atoms with respect to primitive lattice vectors
    real(dp),    intent(in)  :: atom_position(3, num_atom)

    ! the reserved number of the TB basis (dim_basis = num_basis)
    integer,     intent(in)  :: dim_basis 

    ! the number of orthogonal local orbitals for the k-point
    integer,     intent(in)  :: num_basis 

    ! the local orbital information on each atom. Detailed explainations can be found in Table S3
    ! could be 1,3,4,5,6,7,8,9 for s,p,s+p,d,s+d,f,p+d,s+p+d 
    integer,     intent(in)  :: angularmom(num_atom)

    ! the convention of the local orbitals on each atom
    ! if orbt=1, local orbitals are in the order of Table S3
    ! if orbt=2, local orbitals are in the order as implemented in Wannier90
    integer,     intent(in)  :: orbt

    ! the rotation part of space-group operations with respect to primitive lattice vectors
    integer,     intent(in) :: rot_unitary(3,3,num_litt_group_unitary)

    ! the rotation part of space-group operations in Cartesian coordinates
    real(dp),    intent(in) :: SO3_unitary(3,3,num_litt_group_unitary)

    ! the k-dependent phase factors due to the translation part
    complex(dp), intent(out) :: kphase(num_litt_group_unitary)

    ! the transformation matrices of Rs in the orthogonal TB basis
    complex(dp), intent(out) :: rot_mat_tb(dim_basis, dim_basis, num_litt_group_unitary)


    ! parameter used inside the library 
    integer :: i, j, irot, iatom, jatom 
    
    complex(dp) :: srotmt(2,2), protmt(3,3), drotmt(5,5), frotmt(7,7)
    complex(dp) :: cmat3(3,3), cmat5(5,5), cmat7(7,7)

    real(dp) :: RPOS(3), diffr(3)
    real(dp) :: diff, ang

    integer  :: invrot(3,3,num_litt_group_unitary)

    integer :: norb_i, norb_j, startorb_i, startorb_j
    complex(dp) :: phk 

    allocate(rot_orb(mix2l, mix2l, num_atom, num_litt_group_unitary))
    allocate(startorb(num_atom))
    allocate(rot_atom(num_atom, num_litt_group_unitary))
    allocate(rot_phivec(3,num_atom,num_litt_group_unitary))
    do irot = 1, num_litt_group_unitary
        call invmati(rot_unitary(:,:,irot), invrot(:,:,irot)) 
        call O3_to_axis_angle(SO3_unitary(:,:,irot), axis(:,irot), angle(irot), det(irot))
        angle(irot) = angle(irot)*180.d0/PI
        ! write(150,*)irot
        ! write(150,*)SO3_unitary(:,:,irot)
        ! write(150,*)axis(:,irot),angle(irot), det(irot)

        call Dmatrix(rnx=axis(1,irot),rny=axis(2,irot),rnz=axis(3,irot),degree=angle(irot),twoja1=3,Dmat=protmt)
        if (orbt == 2) then
            protmt = matmul(protmt, pwann)
            cmat3  = transpose(pwann)
            protmt = matmul(cmat3, protmt) 
        endif 

        call Dmatrix(rnx=axis(1,irot),rny=axis(2,irot),rnz=axis(3,irot),degree=angle(irot),twoja1=5,Dmat=drotmt)
        if (orbt == 2) then 
            drotmt = matmul(drotmt, dwann)
            cmat5 = transpose(dwann)
            drotmt = matmul(cmat5, drotmt)
        endif 

        call Dmatrix(rnx=axis(1,irot),rny=axis(2,irot),rnz=axis(3,irot),degree=angle(irot),twoja1=7,Dmat=frotmt)  !! modified by rhzhang on 23/12/2023
        if (orbt == 2) then 
            frotmt = matmul(frotmt, fwann)
            cmat7 = transpose(fwann)
            frotmt = matmul(cmat7, frotmt)
        endif

        call Dmatrix(rnx=axis(1,irot),rny=axis(2,irot),rnz=axis(3,irot),degree=angle(irot),twoja1=2,Dmat=srotmt)

        do iatom = 1, num_atom
            if     (angularmom(iatom) == 1) then 
                rot_orb(1,1,iatom,irot) = 1.d0 
            elseif (angularmom(iatom) == 3) then 
                rot_orb(1:3,1:3,iatom,irot) = real(protmt, dp)
                if (det(irot) < 1.d-2) rot_orb(1:3,1:3,iatom,irot) = -real(protmt, dp)
            elseif (angularmom(iatom) == 4) then 
                rot_orb(1,1,iatom,irot) = 1.d0 
                rot_orb(2:4,2:4,iatom,irot) = real(protmt, dp)
                if (det(irot) < 1.d-2) rot_orb(2:4,2:4,iatom,irot) = -real(protmt, dp)
            elseif (angularmom(iatom) == 5) then 
                rot_orb(1:5,1:5,iatom,irot) = real(drotmt, dp)
            elseif (angularmom(iatom) == 6) then 
                rot_orb(1,1,iatom,irot) = 1.d0
                rot_orb(2:6,2:6,iatom,irot) = real(drotmt, dp)
            elseif (angularmom(iatom) == 7) then 
                rot_orb(1:7,1:7,iatom,irot) = real(frotmt, dp)
                if (det(irot) < 1.d-2) rot_orb(1:7,1:7,iatom,irot) = -real(frotmt, dp)
            elseif (angularmom(iatom) == 8) then 
                rot_orb(1:3,1:3,iatom,irot) = real(protmt, dp)
                rot_orb(4:8,4:8,iatom,irot) = real(drotmt, dp)
                if (det(irot) < 1.d-2) rot_orb(1:3,1:3,iatom,irot) = -real(protmt, dp)
            elseif (angularmom(iatom) == 9) then 
                rot_orb(1,1,iatom,irot) = 1.d0 
                rot_orb(2:4,2:4,iatom,irot) = real(protmt, dp)
                rot_orb(5:9,5:9,iatom,irot) = real(drotmt, dp)
                if (det(irot) < 1.d-2) rot_orb(2:4,2:4,iatom,irot) = -real(protmt, dp)
            elseif (angularmom(iatom) ==12) then 
                rot_orb(1:5,1:5,iatom,irot) = real(drotmt, dp)
                rot_orb(6:12,6:12,iatom,irot) = real(frotmt, dp)
                if (det(irot) < 1.d-2) rot_orb(6:12,6:12,iatom,irot) = -real(frotmt, dp)
            else
                write(*,*) "Please set the correct angularmom for each atom"
                stop
            endif 
        enddo 
    enddo

    startorb = 0
    do iatom = 1, num_atom
        if (iatom /= 1) startorb(iatom) = startorb(iatom-1) + angularmom(iatom-1)
    enddo 

    rot_atom = 0
    do irot = 1, num_litt_group_unitary
        do iatom = 1, num_atom 
            RPOS(:) = matmul(rot_unitary(:,:,irot), atom_position(:,iatom)) + tau_unitary(:,irot)
            do jatom = 1, num_atom
                diffr = RPOS(:) - atom_position(:,jatom)
                if (abs(diffr(1)-nint(diffr(1))) .lt. 0.1d-2 .and. &
                    abs(diffr(2)-nint(diffr(2))) .lt. 0.1d-2 .and. &
                    abs(diffr(3)-nint(diffr(3))) .lt. 0.1d-2 ) then 
                    exit
                endif 
            enddo
            if (jatom == num_atom + 1) then 
                write(*,*) "rotation on atoms error", iatom
                stop 
            endif 
            rot_atom(iatom, irot) = jatom
            rot_phivec(:,iatom,irot) = -matmul(invrot(:,:,irot),tau_unitary(:,irot)) &
                                        +matmul(invrot(:,:,irot),nint(diffr(:))) 
        enddo 
    enddo 


    ! get kphase
    kphase = 0.d0 
    do irot = 1, num_litt_group_unitary
        ang = -2.d0*PI*( WK(1)*tau_unitary(1,irot) &
                          +WK(2)*tau_unitary(2,irot) &
                          +WK(3)*tau_unitary(3,irot))
        kphase(irot) = cmplx(dcos(ang), dsin(ang))
    enddo 

    ! get rot_mat_tb
    rot_mat_tb = 0.d0
    do irot = 1, num_litt_group_unitary
        do iatom = 1, num_atom 
            jatom = rot_atom(iatom, irot)
            norb_i = angularmom(iatom)
            norb_j = angularmom(jatom)
            startorb_i = startorb(iatom)
            startorb_j = startorb(jatom)
            if (norb_i /= norb_j) then 
                write(*,*) "error symmetry"
                stop
            endif 
            phk = -dot_product(WK(:),rot_phivec(:,iatom,irot))
            rot_mat_tb((/(startorb_j + i, i=1, norb_j)/), (/(startorb_i + j, j=1, norb_i)/), irot) = &
                rot_orb(1:norb_j, 1:norb_i, iatom, irot)*exp(2.d0*PI*cmplx_i*phk)
        enddo 
    enddo 
    deallocate(rot_orb)
    deallocate(startorb)
    deallocate(rot_atom)
    deallocate(rot_phivec)

end subroutine tb_setup_ssg

subroutine get_k_name(kpoints, kname, k_frac_symbol)

    use bilbao
    implicit none
    real(dp), intent(in) :: kpoints(3)
    character(len=4), intent(out) :: kname
    ! character(len=15), intent(out) :: kpath_name
    character(len=15), intent(out) :: k_frac_symbol
    real(dp) :: k_tmp(3)
    k_tmp = kpoints
    kname = ''
    k_frac_symbol = ''
    call bilbao_getkid(k_tmp, kname, k_frac_symbol)

end subroutine get_k_name

subroutine load_bilbao(sgn)

    use bilbao
    implicit none
    integer, intent(in) :: sgn
    call bilbao_read(sgn)
end subroutine load_bilbao
