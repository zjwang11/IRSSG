module linear_rep
    use lib_params, only:dp
    implicit none
    public :: get_ch_from_op
    ! integer, parameter :: dp = kind(1d0)
contains
    ! Helper function to compare the tau values
    logical function identity_tau(tau1, tau2, tol)
        real(dp), intent(in) :: tau1(3), tau2(3)
        real(dp), intent(in) :: tol
        real(dp) :: diff(3)
        ! Use nearest-integer wrapping to be robust against values near 1.0 ± eps
        diff = tau1 - tau2
        diff = diff - nint(diff)
        identity_tau = sqrt(dot_product(diff,diff)) < tol
    end function identity_tau

    ! Rearranging function based on time-reversal order
    subroutine rearrange_ops(num_litt_group, time_reversal, order)
        integer, intent(in) :: num_litt_group
        integer, intent(in) :: time_reversal(num_litt_group)
        
        integer, intent(out) :: order(num_litt_group)

        integer :: i, pos

        pos = 0
        do i = 1, num_litt_group
            if (time_reversal(i) == 1) then
                pos = pos + 1
                order(pos) = i
            end if
        end do

        do i = 1, num_litt_group
            if (time_reversal(i) == -1) then
                pos = pos + 1
                order(pos) = i
            end if
        end do

        ! write(*,*) order

    end subroutine rearrange_ops

    ! Main function for finding the index
    subroutine find_index(num_litt_group,rot_new, tau_new, su2_new, t_new, space_rot, space_tau, su2_list, time_reversal, index)
        integer, intent(in) :: num_litt_group
        integer, intent(in) :: rot_new(3,3)
        integer, intent(in) :: space_rot(3,3,num_litt_group)
        real(dp), intent(in) :: tau_new(3)
        integer, intent(in) :: t_new
        real(dp), intent(in) :: space_tau(3,num_litt_group)
        complex(dp), intent(in) :: su2_new(2,2)
        complex(dp), intent(in) :: su2_list(2,2,num_litt_group)
        integer, intent(in) :: time_reversal(num_litt_group)
        integer, intent(out) :: index
        integer :: i
        logical :: id
        real(dp) :: tol
        ! Find matching index
        tol = 1e-3
        index = 0
        do i = 1, num_litt_group
            if (all(abs(rot_new - space_rot(:,:,i)) < tol) .and. abs(t_new - time_reversal(i)) < tol) then
                id = identity_tau(tau_new, space_tau(:,i),tol)
                if (id) then
                    if (maxval(abs(su2_new - su2_list(:,:,i))) < tol) then
                        index = i
                        return
                    else if (maxval(abs(su2_new + su2_list(:,:,i))) < tol) then
                        index = -i
                        return
                    end if
                end if
            end if
        end do
        print*, 'Error: Cannot find index of multiple'
        stop
    end subroutine find_index



subroutine get_table(num_litt_group, space_rot, space_tau, spin, su2, time_reversal, &
                     multiple_table, factor_su2)
   implicit none
   integer,               intent(in) :: num_litt_group
   integer,               intent(in) :: space_rot(3,3,num_litt_group)
   real(dp),              intent(in) :: space_tau(3,num_litt_group)
   real(dp),               intent(in) :: spin(3,3,num_litt_group)
   complex(dp),           intent(in) :: su2(2,2,num_litt_group)
   integer,               intent(in) :: time_reversal(num_litt_group)
   integer,               intent(out):: multiple_table(num_litt_group,num_litt_group)
   complex(dp),           intent(out):: factor_su2   (num_litt_group,num_litt_group)
   integer          :: i, j, t_new, index
   integer          :: rot_new(3,3)
   real(dp)         :: tau_new(3)
   complex(dp)      :: su2_new(2,2)

   do i = 1, num_litt_group
      do j = 1, num_litt_group
         rot_new = matmul( space_rot(:,:,i), space_rot(:,:,j) )
         tau_new = matmul( space_rot(:,:,i), space_tau(:,j) ) + space_tau(:,i)

         if (time_reversal(i) < 0) then
            su2_new = matmul( su2(:,:,i), conjg( su2(:,:,j) ) )
         else
            su2_new = matmul( su2(:,:,i), su2(:,:,j) )
         end if

         t_new = time_reversal(i) * time_reversal(j)
         call find_index(num_litt_group, rot_new, tau_new, su2_new, t_new, &
                         space_rot, space_tau, su2, time_reversal, index)

         if (index < 0) then
            multiple_table(i,j) = -index                
            factor_su2  (i,j) = cmplx(-1.0_dp, 0.0_dp)  
         else
            multiple_table(i,j) =  index
            factor_su2  (i,j) = cmplx( 1.0_dp, 0.0_dp)
         end if
      end do
   end do
#ifdef TEST
    do i=1,num_litt_group
        write(1340,'(10000I4)')multiple_table(i,1:num_litt_group)*nint(real(factor_su2(i,1:num_litt_group)))
    enddo
#endif
end subroutine get_table

subroutine get_trans_factor(num_litt_group, time_reversal, space_rot, space_tau, kpoint, factor)
      integer,              intent(in)  :: num_litt_group
      integer,              intent(in)  :: time_reversal(num_litt_group)
      integer,              intent(in)  :: space_rot(3,3,num_litt_group)
      real(dp),             intent(in)  :: space_tau(3,num_litt_group)
      real(dp),             intent(in)  :: kpoint(3)
      complex(dp),          intent(out) :: factor(num_litt_group,num_litt_group)

      integer  :: i, j
      real(dp) :: vec(3), dotk
      real(dp) :: PI=3.1415926535

      factor = cmplx(0.0_dp,0.0_dp,dp)

      do i = 1, num_litt_group
         do j = 1, num_litt_group
            ! vec = R_i * tau_j  -  t_i * tau_j
            vec  = matmul( real(space_rot(:,:,i),dp), space_tau(:,j) )                          &
                    - real(time_reversal(i),dp) * space_tau(:,j)

            dotk = sum( kpoint * vec )*2*PI        ! k · vec
            factor(i,j) = exp( cmplx( 0.0_dp, -dotk, kind=dp ) )   ! e^{-i·dotk}
         end do
      end do
   end subroutine get_trans_factor

    subroutine regular_rep(num_litt_group, mul_table, su2_table, trans_table, rep_mat)
    implicit none
    integer,               intent(in)  :: num_litt_group
    integer,               intent(in)  :: mul_table(num_litt_group,num_litt_group)
    complex(dp),           intent(in)  :: su2_table (num_litt_group,num_litt_group)
    complex(dp),           intent(in)  :: trans_table(num_litt_group,num_litt_group)
    complex(dp),           intent(out) :: rep_mat(num_litt_group,num_litt_group,num_litt_group)
    complex(dp) :: factor(num_litt_group,num_litt_group)
    integer      :: g, i, m

    factor = su2_table * trans_table  

    rep_mat = cmplx(0.0_dp, 0.0_dp)

    do g = 1, num_litt_group
        do i = 1, num_litt_group  
            m = mul_table(g,i)       
        
            rep_mat(m, i, g) = factor(g, i)
        end do
    end do
    end subroutine regular_rep

    subroutine generate_random_matrix(g, mat)
    integer, intent(in)              :: g
    complex(dp), intent(out)         :: mat(g,g)
    real(dp)                         :: r(g,g), s(g,g)

    if (g == 2) then
      ! Fixed 2x2 matrix
      mat(1,1) = cmplx(1.23456_dp, 2.14793_dp, kind=dp)
      mat(1,2) = cmplx(1.57993_dp, 4.22339_dp, kind=dp)
      mat(2,1) = cmplx(3.4563277_dp, 1.778431_dp, kind=dp)
      mat(2,2) = cmplx(2.3399532_dp, 5.114779_dp, kind=dp)
    else
      call random_number(r)   
      call random_number(s)   

      mat = cmplx(10.0_dp * r, 10.0_dp * s, kind=dp)
    end if

  end subroutine generate_random_matrix


subroutine diag_hermitian(A_in, n, w, Z)
   implicit none

   !==== Interface ====
   integer,                  intent(in)  :: n              ! Matrix order
   complex(dp),              intent(in)  :: A_in(n,n)      ! Input Hermitian matrix
   real(dp),                 intent(out) :: w(n)           ! Real eigenvalues (ascending)
   complex(dp),              intent(out) :: Z(n,n)         ! Columns = normalized eigenvectors

   !==== Locals ====
   complex(dp), allocatable :: work(:)
   real(dp),    allocatable :: rwork(:)
   complex(dp)              :: work_query
   integer                  :: lwork, info

   Z = A_in

   lwork = -1
   allocate(rwork(max(1,3*n-2)))   
   call zheev('V','U', n, Z, n, w, work_query, lwork, rwork, info)
   if (info /= 0) then
      print *, 'zheev workspace query failed, info=', info
      stop
   end if
   lwork = int(real(work_query))     
   allocate(work(lwork))

   call zheev('V','U', n, Z, n, w, work, lwork, rwork, info)
   if (info /= 0) then
      print *, 'zheev failed, info=', info
      stop
    else
#ifdef TEST
      write(380,*)w(1:n)
#endif
   end if

   deallocate(work, rwork)
end subroutine diag_hermitian


subroutine get_irreducible_rep(num_litt_group, time_reversal, k_point, factor, mul_table, tau_list, rep_mat, irrep_num, ch_table, phase, torsion_out)
    implicit none
    ! Input parameters
    integer, intent(in)    :: num_litt_group              ! Total number of symmetry operations (group elements)
    integer, intent(in)    :: time_reversal(num_litt_group) ! Array indicating time-reversal (1 for unitary, -1 for antiunitary) for each operation (already rearranged)
    real(dp), intent(in)   :: k_point(3)                  ! k-point vector (Bloch wavevector)
    complex(dp), intent(in) :: factor(num_litt_group, num_litt_group)
    integer, intent(in)    :: mul_table(num_litt_group, num_litt_group)
    real(dp), intent(in)   :: tau_list(3, num_litt_group)
    complex(dp), intent(in) :: rep_mat(num_litt_group, num_litt_group, num_litt_group)
    ! Output parameters
    integer, intent(out)   :: irrep_num                   ! Number of irreducible representations found
    complex(dp), intent(out) :: ch_table(num_litt_group, num_litt_group)
    complex(dp), intent(out) :: phase(num_litt_group)     ! Phase factors for each unitary operation
    ! Local variables
    integer :: nu               ! nu: number of unitary operations (time_reversal = 1)
    integer :: i, j, m, block_idx, ir   ! loop indices
    complex(dp), allocatable :: H(:,:), H_sum(:,:), H_hermitian(:,:)    ! matrices used for random Hamiltonian and its Hermitian form
    complex(dp), allocatable :: eigenvectors(:,:)        ! matrix of eigenvectors
    real(dp),    allocatable :: eigenvalues(:)           ! array of eigenvalues
    complex(dp), allocatable :: E_adj(:,:), E_conj(:,:)  ! E_adj = conjugate transpose of eigenvector matrix, E_conj = elementwise conjugate of eigenvector matrix
    complex(dp), allocatable :: block_rep_all(:,:,:)     ! block_rep_all(:,:,i) will store the block-diagonal form of rep_mat(:,:,i) in the eigenbasis
    integer :: rep_size_list(num_litt_group)             ! Array to hold sizes of degenerate blocks (repetitions)
    integer :: nb                                       ! Number of degenerate blocks (i.e., number of candidate irreps before filtering duplicates)
    complex(dp) :: sum_ch, kappa1, kappa2_omega, kappa2_mul  ! variables for character orthonormality checks
    complex(dp) :: val
    logical :: check_ir, fail_flag
    integer :: loop_num
    integer :: unique_count                             ! Number of unique irreducible characters found (rep_num before sorting)
    integer :: sorted_indices(num_litt_group)           ! Array to hold sorted indices for irreps
    real(dp) :: tor                                     ! Used for computing torsion (Kramers degeneracy indicator)
    integer,intent(out),optional :: torsion_out(num_litt_group)                  ! Torsion values for each irrep (2 for Kramers degenerate, 0 for unitary, etc.)
    integer :: torsion(num_litt_group)
    ! Character tables (temporary storage)
    complex(dp), allocatable :: char_mat(:,:), char_now_mat(:,:) 
    ! Determine group order and number of unitary operations

    complex(dp) :: tmp_col(num_litt_group)
    logical :: complete_flag
    real(dp) :: PI=3.1415926535


    nu = 0
    do i = 1, num_litt_group
        if (time_reversal(i) == 1) nu = nu + 1
    end do

    ! Allocate dynamic arrays based on group size
    allocate(H(num_litt_group, num_litt_group), H_sum(num_litt_group, num_litt_group), H_hermitian(num_litt_group, num_litt_group))
    allocate(eigenvalues(num_litt_group), eigenvectors(num_litt_group, num_litt_group))
    allocate(E_adj(num_litt_group, num_litt_group), E_conj(num_litt_group, num_litt_group))
    allocate(block_rep_all(num_litt_group, num_litt_group, num_litt_group))
    allocate(char_mat(num_litt_group, num_litt_group), char_now_mat(num_litt_group, num_litt_group))
    ch_table = cmplx(0.0_dp, 0.0_dp, dp)  ! Initialize output table to zero

    check_ir = .false.
    loop_num = 0

    ! Attempt up to 5 random Hermitian matrices to achieve a full irreducible decomposition
    do while (.not. check_ir .and. loop_num < 5)
        loop_num = loop_num + 1

        ! 1. Generate a random complex matrix H of size (num_litt_group x num_litt_group)
        call generate_random_matrix(num_litt_group, H)

        ! 2. Compute H_sum = Σ_{i}( R_i * H * R_i^{-1} ) for unitary ops,
        !    and H_sum = Σ_{i}( R_i * conj(H) * R_i^{-1} ) for antiunitary ops.
        H_sum = cmplx(0.0_dp, 0.0_dp, dp)         ! initialize sum matrix to zero
        do i = 1, num_litt_group
            if (time_reversal(i) == 1) then
                ! Unitary operation: conjugation by rep_mat(:,:,i)
                H_sum = H_sum + matmul( rep_mat(:,:,i), matmul(H, conjg(transpose(rep_mat(:,:,i)))) )
            else
                ! Antiunitary operation: conjugation by rep_mat and complex conjugation of H
                H_sum = H_sum + matmul( rep_mat(:,:,i), matmul(conjg(H), conjg(transpose(rep_mat(:,:,i)))) )
            end if
        end do

        ! 3. Make H_sum Hermitian: H_hermitian = H_sum + (H_sum)†
        H_hermitian = H_sum + conjg(transpose(H_sum))
        ! Now H_hermitian is Hermitian (self-adjoint), ensuring real eigenvalues.

        ! 4. Compute eigenvalues and eigenvectors of H_hermitian.
        !    (Using an external routine or LAPACK is necessary; not implemented here.)
        !    The eigenvectors will be used as the transformation to block-diagonalize the representation.
        ! TODO: Use LAPACK (e.g., ZHEEV) or similar to solve for eigenvalues/vectors of H_hermitian
        call diag_hermitian(H_hermitian, num_litt_group, eigenvalues, eigenvectors)  ! Placeholder call

        ! 5. Sort eigenvalues in descending order (and reorder eigenvectors accordingly)
        do i = 1, num_litt_group-1
            m = i
            do j = i+1, num_litt_group
                if (eigenvalues(j) > eigenvalues(m)) m = j
            end do
            if (m /= i) then
                ! Swap eigenvalues(i) and eigenvalues(m)
                val = cmplx(eigenvalues(i), 0.0_dp, dp)
                eigenvalues(i) = eigenvalues(m)
                eigenvalues(m) = real(val)  ! real(val) gets original eigenvalues(i)
                ! Swap corresponding eigenvector columns i and m
                tmp_col(1:num_litt_group) = eigenvectors(1:num_litt_group, i)  ! use E_conj as temp storage for column
                eigenvectors(1:num_litt_group, i) = eigenvectors(1:num_litt_group, m)
                eigenvectors(1:num_litt_group, m) = tmp_col(1:num_litt_group)
            end if
        end do

        ! 6. Determine degeneracy structure of eigenvalues to identify block sizes.
        nb = 0
        rep_size_list(:) = 0      ! reset or initialize
        m = 1
        do i = 2, num_litt_group
            if (abs(eigenvalues(i) - eigenvalues(i-1)) <= 1e-3_dp) then
                m = m + 1   ! continue counting equal eigenvalues (degeneracy)
            else
                nb = nb + 1
                rep_size_list(nb) = m
                m = 1       ! reset count for next block
            end if
        end do
        nb = nb + 1
        rep_size_list(nb) = m
        ! (At this point, rep_size_list(1:nb) contains the sizes of each degenerate eigenvalue block, and sum(rep_size_list(1:nb)) = num_litt_group.)

        ! 7. Form the transformation matrices:
        !    Let E = eigenvectors (columns are eigenvecs). Compute E_adj = E† and E_conj = element-wise conjugate of E.
        E_adj = conjg(transpose(eigenvectors))
        E_conj = conjg(eigenvectors)

        ! 8. Compute block_rep_all(:,:, i) = E_adj * rep_mat(:,:,i) * E    for unitary i
        !                              and = E_adj * rep_mat(:,:,i) * E_conj for antiunitary i.
        do i = 1, num_litt_group
            if (time_reversal(i) == 1) then
                block_rep_all(:,:, i) = matmul( E_adj, matmul(rep_mat(:,:,i), eigenvectors) )
            else
                block_rep_all(:,:, i) = matmul( E_adj, matmul(rep_mat(:,:,i), E_conj) )
            end if
        end do

        ! 9. Calculate character contributions for each unitary operation across each degenerate block.
        !    We sum the diagonal elements of each block of block_rep_all for unitary operations.
        char_mat = cmplx(0.0_dp, 0.0_dp, dp)
        unique_count = 0
        ! Iterate over each degenerate block defined by rep_size_list
        m = 1   ! m will track the starting index of the current block in the matrix
        
        do block_idx = 1, nb
            ! Size of this block
            j = rep_size_list(block_idx)
            ! Sum diagonal elements of this block for each unitary operation
            do i = 1, nu   ! loop over unitary operations (columns 1..nu correspond to time_reversal = 1)
                sum_ch = cmplx(0.0_dp, 0.0_dp, dp)
                ! Sum diagonal entries from m to m+j-1 for matrix block_rep_all of operation i
                do ir = 0, j-1
                    sum_ch = sum_ch + block_rep_all(m+ir, m+ir, i)
                end do
                char_mat(block_idx, i) = sum_ch
            end do
            m = m + j   ! advance to the starting index of the next block
        end do

        ! 10. Filter out duplicate irreducible characters:
        !     Some blocks may correspond to identical irreducible representations (degenerate eigenvalues can repeat the same irrep).
        unique_count = 0
        do block_idx = 1, nb
            ! Check if the character vector of this block (char_mat(block_idx, 1:nu)) is already in char_now_mat
            if (unique_count == 0) then
                unique_count = unique_count + 1
                char_now_mat(unique_count, 1:nu) = char_mat(block_idx, 1:nu)
            else
                fail_flag = .false.
                do ir = 1, unique_count
                    ! Compute norm difference between existing char_now_mat(ir,:) and current char_mat(block_idx,:)
                    sum_ch = cmplx(0.0_dp, 0.0_dp, dp)
                    do j = 1, nu
                        ! Sum of squared differences (use Euclidean norm for comparison)
                        sum_ch = sum_ch + (char_now_mat(ir, j) - char_mat(block_idx, j)) * conjg(char_now_mat(ir, j) - char_mat(block_idx, j))
                    end do
                    if (sqrt(real(sum_ch)) < 1e-4_dp) then
                        fail_flag = .true.    ! found a match (character already listed)
                        exit
                    end if
                end do
                if (.not. fail_flag) then
                    unique_count = unique_count + 1
                    char_now_mat(unique_count, 1:nu) = char_mat(block_idx, 1:nu)
                end if
            end if
        end do
        irrep_num = unique_count   ! provisional count of irreps found

        ! 11. Compute the **linear character** including the Bloch phase factor for each irreducible representation.
        !     We multiply each character by exp(-i * k_point · tau) for each unitary operation.
        phase = cmplx(0.0_dp, 0.0_dp, dp)  ! Initialize phase factors
        do j = 1, nu
            phase(j) = exp( cmplx(0.0_dp, - 2*PI*dot_product(k_point, tau_list(:, j)), dp) )
        end do

        do ir = 1, irrep_num
            do j = 1, nu
                ch_table(ir, j) = char_now_mat(ir, j) * phase(j)
            end do
        end do

        ! 12. Sort the irreducible representations (rows of char_now_mat/ch_table) by the custom lexicographical order:
        !     Sort primarily by real parts of the character values (for each unitary op in sequence, tolerance 1e-2),
        !     and secondarily by imaginary parts (tolerance 1e-3).
        !     This ensures a consistent ordering of the irreps.
        sorted_indices = 0
        do i = 1, irrep_num
            sorted_indices(i) = 0
        end do
        
        call sort_perm(irrep_num, nu, char_now_mat(1:irrep_num,1:nu), sorted_indices(1:irrep_num))

        ! do i = 1, irrep_num
        !     if (sorted_indices(i) == 0) sorted_indices(i) = i
        ! end do

        ! Reorder char_now_mat and ch_table according to sorted_indices
        if (irrep_num > 1) then
            if (allocated(E_conj)) then
                deallocate(E_conj)
            endif

            if (allocated(E_adj)) then
                deallocate(E_adj)
            endif

            allocate(E_conj(irrep_num, nu))         ! reuse E_conj as a temporary matrix for swapping
            allocate(E_adj(irrep_num, nu))
            ! Copy sorted data to temp arrays
            do i = 1, irrep_num
                E_conj(i, 1:nu) = char_now_mat(sorted_indices(i), 1:nu)
                E_adj(i, 1:nu)  = ch_table(sorted_indices(i), 1:nu)
            end do
            ! Move sorted data back to char_now_mat and ch_table
            do i = 1, irrep_num
                char_now_mat(i, 1:nu) = E_conj(i, 1:nu)
                ch_table(i, 1:nu)     = E_adj(i, 1:nu)
            end do
            deallocate(E_conj, E_adj)
        end if

        ch_table(1:irrep_num,1:nu) = char_now_mat(1:irrep_num,1:nu)

        ! 13. Compute torsion for each irreducible representation:
        !     (Torsion indicates the degeneracy under time-reversal: e.g., 2 for Kramers doublet, 0 for no degeneracy)
        if (any(time_reversal(:) == -1)) then
            do ir = 1, irrep_num
                tor = 0.0_dp
                do i = 1, num_litt_group
                    if (time_reversal(i) > 0) then
                        ! Only include unitary operations in the sum
                        if (i <= nu) then
                            tor = tor + (abs(char_now_mat(ir, i))**2)
                        endif
                    end if
                end do
                tor = tor / real(num_litt_group, dp) * 2.0_dp
                if (abs(tor - nint(tor)) > 1e-3_dp) then
                    print *, 'Warning: torsion not an integer (possible error) for irrep', ir
                end if
                torsion(ir) = nint(tor)   ! round to nearest integer (should be 0, 1, or 2 typically)
            end do
        else
            do ir = 1, irrep_num
                torsion(ir) = 0
            end do
        end if
#ifdef TEST
        write(1223,*)torsion(1:irrep_num)
#endif
        ! 14. Verify the orthonormality of characters (irreducibility check):
        fail_flag = .false.
        if (any(time_reversal(:) == -1)) then
            ! For magnetic (antiunitary present) groups
            do ir = 1, irrep_num
                sum_ch = cmplx(0.0_dp, 0.0_dp, dp)
                ! Sum_{i=1 to nu} [ |χ_ir(i)|^2 + χ_ir( mul_table(nu+i, nu+i) ) * ω_i ] / (2*nu)
                do i = 1, nu
                    kappa1      = char_now_mat(ir, i) * conjg(char_now_mat(ir, i))
                    kappa2_omega = factor(nu+i, nu+i)
                    m = mul_table(nu+i, nu+i)    ! result of (antiunitary_i * antiunitary_i)
                    if (m >= 1 .and. m <= nu) then
                        kappa2_mul = char_now_mat(ir, m)
                    else
                        kappa2_mul = cmplx(0.0_dp, 0.0_dp, dp)  ! if the product yields an antiunitary, treat character as 0
                    end if
                    sum_ch = sum_ch + (kappa1 + kappa2_mul * kappa2_omega) / (2.0_dp * real(nu, dp))
                end do
#ifdef TEST
                write(190,*)sum_ch
#endif
                if (abs(sum_ch - 1.0_dp) > 1e-2_dp) then
                    fail_flag = .true.
                    exit
                end if
            end do
        else
            ! For purely unitary groups (no time-reversal elements)
            do ir = 1, irrep_num
                sum_ch = cmplx(0.0_dp, 0.0_dp, dp)
                ! Sum_{i=1 to nu} |χ_ir(i)|^2 / nu
                do i = 1, nu
                    kappa1 = char_now_mat(ir, i) * conjg(char_now_mat(ir, i))
                    sum_ch = sum_ch + kappa1 / real(nu, dp)
                end do
                if (abs(sum_ch - 1.0_dp) > 1e-2_dp) then
                    fail_flag = .true.
                    exit
                end if
            end do
        end if

        if (.not. fail_flag) then
            ! All orthonormality checks passed: irreducible decomposition found
            check_ir = .true.
        else
            ! Failed orthonormality check, try again with a new random matrix if attempts remain
            print *, "retry!!"
            if (loop_num == 5) print *, "something wrong"
        end if
    end do  ! end of while loop for attempts

    ! If check_ir is false after 5 attempts, we proceed with the last obtained result (with a warning)
    if (.not. check_ir) then
        print *, "Warning: Irreducible representation check failed after 5 attempts."
    end if

    ! At this point, ch_table(1:irrep_num, 1:nu) contains the linear characters for each irreducible representation.
    ! (The remaining entries of ch_table are zero.)

    ! Deallocate local allocatable arrays (they will also auto-deallocate when subroutine exits)
    if (allocated(E_adj)) then
        deallocate(E_adj)
    endif

    if (allocated(E_conj)) then
        deallocate(E_conj)
    endif

    if (present(torsion_out)) then
        ! If torsion_out is provided, copy torsion values to it
            torsion_out(1:num_litt_group) = torsion(1:num_litt_group)
    end if

    deallocate(H, H_sum, H_hermitian, eigenvalues, eigenvectors, block_rep_all, char_mat, char_now_mat)
end subroutine get_irreducible_rep


    ! Main entry point for get_ch_from_op
    subroutine get_ch_from_op(num_litt_group,spin_list, rot_list, tau_list, su2_list, time_reversal, kpoint, order, &
                            irrep_num, ch_table, phase, irrep_unitary_num, ch_unitary_table, torsion)
        integer, intent(in) :: num_litt_group
        integer, intent(in) :: time_reversal(num_litt_group)
        integer, intent(in) :: rot_list(3,3,num_litt_group)
        real(dp), intent(in) :: spin_list(3,3,num_litt_group)
        real(dp), intent(in) :: tau_list(3,num_litt_group)
        complex(dp), intent(in) :: su2_list(2,2,num_litt_group)
        real(dp), intent(in) :: kpoint(3)
        
        integer,intent(out) :: irrep_num
        integer,intent(out) :: irrep_unitary_num
        complex(dp), intent(out) :: ch_table(num_litt_group,num_litt_group)
        integer, intent(out) :: order(num_litt_group)
        complex(dp), intent(out) :: ch_unitary_table(num_litt_group,num_litt_group)


        integer :: time_reversal_(num_litt_group)
        integer :: rot_list_(3,3,num_litt_group)
        real(dp) :: spin_list_(3,3,num_litt_group)
        real(dp) :: tau_list_(3,num_litt_group)
        complex(dp) :: su2_list_(2,2,num_litt_group)
        integer :: i,j
        complex(dp) :: time_rev_mat(2,2)
        complex(dp),allocatable :: rep_mat(:,:,:)

        integer :: num_litt_group_unitary
        complex(dp),allocatable :: factor(:,:)
        complex(dp),allocatable :: trans_table(:,:)
        integer,allocatable  :: mul_table(:,:)
        integer,intent(out) :: torsion(num_litt_group)
        complex(dp), intent(out) :: phase(num_litt_group)

        allocate(mul_table(num_litt_group,num_litt_group))
        allocate(rep_mat(num_litt_group,num_litt_group,num_litt_group))
        allocate(factor(num_litt_group,num_litt_group))
        allocate(trans_table(num_litt_group,num_litt_group))


        time_rev_mat = 0.0_dp
        time_rev_mat(1,2) = cmplx(-1.0_dp,0.0_dp,dp)
        time_rev_mat(2,1) = cmplx(1.0_dp,0.0_dp,dp)

        ! Call functions to process and get representations
        call rearrange_ops(num_litt_group, time_reversal, order)

        num_litt_group_unitary = 0
        do i = 1, num_litt_group
            rot_list_(:,:,i)   = rot_list(:,:, order(i))
            spin_list_(:,:,i)  = spin_list(:,:,     order(i))
            tau_list_(:,i)     = tau_list(:,  order(i))
            time_reversal_(i)      = time_reversal(order(i))
            if(time_reversal_(i) == -1) then
                su2_list_(:,:,i)   = matmul(su2_list(:,:,     order(i)),time_rev_mat)
            else
                su2_list_(:,:,i)   = su2_list(:,:,     order(i))
                num_litt_group_unitary = num_litt_group_unitary + 1
            endif
        end do

        call get_table(num_litt_group, rot_list_, tau_list_, spin_list_, su2_list_, time_reversal_, mul_table, factor)
        call get_trans_factor(num_litt_group, time_reversal_, rot_list_, tau_list_, kpoint, trans_table)
        call regular_rep(num_litt_group, mul_table, factor, trans_table, rep_mat)
        
        factor = factor * trans_table

        call get_irreducible_rep(num_litt_group, time_reversal_, kpoint, factor, mul_table, tau_list_, rep_mat, irrep_num, ch_table, phase, torsion)
        deallocate(mul_table,factor,trans_table,rep_mat)


        allocate(mul_table(num_litt_group_unitary,num_litt_group_unitary))
        allocate(rep_mat(num_litt_group_unitary,num_litt_group_unitary,num_litt_group_unitary))
        allocate(factor(num_litt_group_unitary,num_litt_group_unitary))
        allocate(trans_table(num_litt_group_unitary,num_litt_group_unitary))


        call get_table(num_litt_group_unitary, rot_list_, tau_list_, spin_list_, su2_list_, time_reversal_, mul_table, factor)
        call get_trans_factor(num_litt_group_unitary, time_reversal_, rot_list_, tau_list_, kpoint, trans_table)
        call regular_rep(num_litt_group_unitary, mul_table, factor, trans_table, rep_mat)
        factor = factor * trans_table
        call get_irreducible_rep(num_litt_group_unitary, time_reversal_, kpoint, factor, mul_table, tau_list_, &
                        rep_mat, irrep_unitary_num, ch_unitary_table(1:num_litt_group_unitary,1:num_litt_group_unitary),phase)
        

        ! call get_discriminant_value(num_litt_group,num_litt_group_unitary,ch_unitary_table(1:irrep_unitary_num,1:num_litt_group_unitary),irrep_unitary_num,&
        !                                 time_reversal_,rot_list_,tau_list_,su2_list_, &
        !                                 rot_list_,tau_list_,su2_list_,discriminant_value)

        deallocate(mul_table,factor,trans_table,rep_mat)
    end subroutine get_ch_from_op

    subroutine output_character_table(kname,num_litt_group,num_litt_group_unitary,litt_group,order_op,&
                                irrep_num,character_table,phase,irrep_unitary_num,character_unitary_table,&
                                irrep_coirrep_relation,irrep_name_list_reduce,discriminant_value,is_msg)
        character(len=4) :: kname
        integer, intent(in) :: num_litt_group_unitary
        integer, intent(in) :: num_litt_group
        integer, intent(in) :: order_op(num_litt_group)
        integer, intent(in) :: irrep_num
        integer, intent(in) :: irrep_unitary_num
        complex(dp), intent(in) :: character_table(irrep_num,num_litt_group_unitary)
        complex(dp), intent(in) :: character_unitary_table(irrep_unitary_num,num_litt_group_unitary)
        integer, intent(in) :: irrep_coirrep_relation(2,irrep_num)
        integer, intent(in) :: litt_group(num_litt_group)
        complex(dp), intent(in) :: phase(num_litt_group_unitary)
        
        integer :: i,j
        logical :: all_phase_one
        real(dp) :: tol_phase
        integer :: jj
        character(len=100) :: irrep_name
        character(len=15) :: irrep_name_list(irrep_num)
        character(len=5) :: irrep_unitary_name_list(irrep_unitary_num)
        character(len=15),intent(out) :: irrep_name_list_reduce(irrep_num)
        integer, intent(in) :: discriminant_value(num_litt_group)
        logical, intent(in) :: is_msg
        character(len=4) :: tmp
        character(len=1) :: msg_prefix


        if (is_msg) then
            msg_prefix = 'M'
        else
            msg_prefix = ''
        endif

        write(*,*)
        write(154,*)
        write(tmp,'(I0)') num_litt_group_unitary
        write(*,'(A)')trim(adjustl(tmp))//' unitary elements of the little group:'
        write(154,'(A)')trim(adjustl(tmp))//' unitary elements of the little group:'
        do i=1,num_litt_group_unitary
            write(*,'(1I5)',advance='no')litt_group(order_op(i))
            write(154,'(1I5)',advance='no')litt_group(order_op(i))
            if (mod(i,15)==0 .and. i/=num_litt_group_unitary) then
                write(*,*)
                write(154,*)
            endif
        enddo

        write(*,*)
        write(154,*)
        if (num_litt_group_unitary == num_litt_group) then
            write(*,'(A)',advance='no')'No anti-unitary elements of the little group.'
            write(154,'(A)',advance='no')'No anti-unitary elements of the little group.'
        else
            write(tmp,'(I0)') num_litt_group - num_litt_group_unitary
            write(*,'(A)')trim(adjustl(tmp))//' anti-unitary elements of the little group:'
            write(154,'(A)')trim(adjustl(tmp))//' anti-unitary elements of the little group:'
            do i=num_litt_group_unitary+1,num_litt_group
                write(*,'(1I5)',advance='no')litt_group(order_op(i))
                write(154,'(1I5)',advance='no')litt_group(order_op(i))
                if (mod(i-num_litt_group_unitary,15)==0 .and. i/=num_litt_group) then
                    write(*,*)
                    write(154,*)
                endif
            enddo
        endif
        write(*,*)
        write(*,*)
        write(154,*)
        write(154,*)

        irrep_name_list(:) = ''
        irrep_unitary_name_list(:) = ''
        irrep_name_list_reduce(:) = ''

#ifdef TEST
        write(178,'(A)')'Character table for unitary group'
        write(178,'(A,1000I20)')'   ',litt_group(order_op(1:num_litt_group_unitary))
        do i=1,irrep_unitary_num
            write(178,'(1I3,1000F10.5)')i,character_unitary_table(i,1:num_litt_group_unitary)+cmplx(1e-6_dp,1e-6_dp,dp)
        enddo
        write(178,*)
        write(178,'(A)')'Character table for complete group'
        write(178,'(A,1000I20)')'   ',litt_group(order_op(1:num_litt_group_unitary))
        do i=1,irrep_num
            write(178,'(1I3,1000F10.5)',advance='no')i,character_table(i,1:num_litt_group_unitary)+cmplx(1e-6_dp,1e-6_dp,dp)
            if (irrep_coirrep_relation(2,i) == 0) then
                write(178,'(A,1I10)')' ',irrep_coirrep_relation(1,i)
            else
                write(178,'(A,1I10,A,1I10)')' ',irrep_coirrep_relation(1,i),'+',irrep_coirrep_relation(2,i)
            endif
        enddo
#endif


        write(*,'(A)')'Character table for unitary group'
        write(*,'(A)',advance='no')'      tor'
        write(*,'(1I8,1000I13)')litt_group(order_op(1:num_litt_group_unitary))
        do i=1,irrep_unitary_num
            irrep_name = ''
            write(irrep_name,*)i
            write(*,'(A4,A,1I3,A)',advance='no')trim(msg_prefix)//trim(adjustl(kname))//trim(adjustl(irrep_name)),' ',discriminant_value(i),'   '
            irrep_unitary_name_list(i) = trim(msg_prefix)//trim(adjustl(kname))//trim(adjustl(irrep_name))
            do j=1,num_litt_group_unitary
                if (aimag(character_unitary_table(i,j)) >= 0.0_dp) then
                    write(*,'(1F5.2,A,1F4.2,A)',advance='no')real(character_unitary_table(i,j))+1e-6,'+',aimag(character_unitary_table(i,j))+1e-6,'i'
                else
                    write(*,'(1F5.2,A,1F4.2,A)',advance='no')real(character_unitary_table(i,j))+1e-6,'-',-aimag(character_unitary_table(i,j))+1e-6,'i'
                endif
                write(*,'(A)',advance='no')'  '
            enddo
            write(*,*)
        enddo
        
        ! Only print phase line if not all phases are 1
        tol_phase = 1.0e-4_dp
        all_phase_one = .true.
        do jj=1,num_litt_group_unitary
            if (abs(real(phase(jj)) - 1.0_dp) > tol_phase .or. abs(aimag(phase(jj))) > tol_phase) then
                all_phase_one = .false.
                exit
            endif
        enddo

        if (.not. all_phase_one) then
            write(*,'(A)',advance='no')' phase     '
            do j=1,num_litt_group_unitary
                if (aimag(phase(j)) >= 0.0_dp) then
                    write(*,'(1F5.2,A,1F4.2,A)',advance='no')real(phase(j))+1e-6,'+',aimag(phase(j))+1e-6,'i'
                else
                    write(*,'(1F5.2,A,1F4.2,A)',advance='no')real(phase(j))+1e-6,'-',-aimag(phase(j))+1e-6,'i'
                endif
                write(*,'(A)',advance='no')'  '
            enddo
            write(*,*)
        endif

        ! write(154,'(A)')'Character table for unitary group'
        ! write(154,'(1000I13)')litt_group(order_op(1:num_litt_group_unitary))
        ! do i=1,irrep_unitary_num
        !     irrep_name = ''
        !     write(irrep_name,*)i
        !     write(154,'(A4,A)',advance='no')trim(adjustl(kname))//trim(adjustl(irrep_name)),'   '
        !     do j=1,num_litt_group_unitary
        !         if (aimag(character_unitary_table(i,j)) >= 0.0_dp) then
        !             write(154,'(1F5.2,A,1F4.2,A)',advance='no')real(character_unitary_table(i,j)),'+',aimag(character_unitary_table(i,j)),'i'
        !         else
        !             write(154,'(1F5.2,A,1F4.2,A)',advance='no')real(character_unitary_table(i,j)),'-',-aimag(character_unitary_table(i,j)),'i'
        !         endif
        !         write(154,'(A)',advance='no')'  '
        !     enddo
        !     write(154,*)
        ! enddo

        write(154,'(A)')'Character table for unitary group'
        write(154,'(A)',advance='no')'      tor'
        write(154,'(1I8,1000I13)')litt_group(order_op(1:num_litt_group_unitary))
        do i=1,irrep_unitary_num
            irrep_name = ''
            write(irrep_name,*)i
            write(154,'(A4,A,1I3,A)',advance='no')trim(msg_prefix)//trim(adjustl(kname))//trim(adjustl(irrep_name)),' ',discriminant_value(i),'   '
            irrep_unitary_name_list(i) = trim(msg_prefix)//trim(adjustl(kname))//trim(adjustl(irrep_name))
            do j=1,num_litt_group_unitary
                if (aimag(character_unitary_table(i,j)) >= 0.0_dp) then
                    write(154,'(1F5.2,A,1F4.2,A)',advance='no')real(character_unitary_table(i,j))+1e-6,'+',aimag(character_unitary_table(i,j))+1e-6,'i'
                else
                    write(154,'(1F5.2,A,1F4.2,A)',advance='no')real(character_unitary_table(i,j))+1e-6,'-',-aimag(character_unitary_table(i,j))+1e-6,'i'
                endif
                write(154,'(A)',advance='no')'  '
            enddo
            write(154,*)
        enddo

        irrep_name = ''
        do i=1,irrep_num
            irrep_name=trim(adjustl(irrep_name))//','
            if (irrep_coirrep_relation(2,i) == 0) then
                irrep_name_list_reduce(i) = trim(adjustl(irrep_unitary_name_list(irrep_coirrep_relation(1,i))))
                irrep_name=trim(adjustl(irrep_name))//trim(adjustl(irrep_unitary_name_list(irrep_coirrep_relation(1,i))))
            else
                irrep_name_list_reduce(i) = trim(adjustl(irrep_unitary_name_list(irrep_coirrep_relation(1,i))))//&
                                        trim(adjustl(irrep_unitary_name_list(irrep_coirrep_relation(2,i))))
                irrep_name=trim(adjustl(irrep_name))//trim(adjustl(irrep_name_list_reduce(i)))
            endif
            
            
        enddo

        ! if (.not. all_phase_one) then
        !     write(154,'(A)',advance='no')' phase     '
        !     do j=1,num_litt_group_unitary
        !         if (aimag(phase(j)) >= 0.0_dp) then
        !             write(154,'(1F5.2,A,1F4.2,A)',advance='no')real(phase(j))+1e-6,'+',aimag(phase(j))+1e-6,'i'
        !         else
        !             write(154,'(1F5.2,A,1F4.2,A)',advance='no')real(phase(j))+1e-6,'-',-aimag(phase(j))+1e-6,'i'
        !         endif
        !         write(154,'(A)',advance='no')'  '
        !     enddo
        !     write(154,*)
        ! endif

        if (num_litt_group/=num_litt_group_unitary) then 
            write(*,'(A)')'Coirreps for complete group: '//irrep_name(2:)
            write(154,'(A)')'Coirreps for complete group: '//irrep_name(2:)
            
        else
            write(*,'(A)')'Irreps for complete group: '//irrep_name(2:)
            ! write(154,'(A)')'There are no anti-unitary operations in the little group.'
            write(154,'(A)')'Irreps for complete group: '//irrep_name(2:)
        endif

        ! write(154,*)
        ! write(154,'(A)')'Character table for complete group'
        !  write(154,'(1000I13)')litt_group(order_op(1:num_litt_group_unitary))
        ! do i=1,irrep_num
        !     irrep_name = ''
        !     write(irrep_name,*)i
        !     write(154,'(A4,A)',advance='no')'A'//trim(adjustl(irrep_name)),'   '
        !     irrep_name_list(i) = 'A'//trim(adjustl(irrep_name))
        !     do j=1,num_litt_group_unitary

        !         if (aimag(character_table(i,j)) >= 0.0_dp) then
        !             write(154,'(1F5.2,A,1F4.2,A)',advance='no')real(character_table(i,j)),'+',aimag(character_table(i,j)),'i'
        !         else
        !             write(154,'(1F5.2,A,1F4.2,A)',advance='no')real(character_table(i,j)),'-',-aimag(character_table(i,j)),'i'
        !         endif
        !         write(154,'(A)',advance='no')'  '
        !     enddo

        !     if (irrep_coirrep_relation(2,i) == 0) then
        !         write(154,'(A,A)')' = ',trim(adjustl(irrep_unitary_name_list(irrep_coirrep_relation(1,i))))
        !     else
        !         write(154,'(A,A,A)')' = ',trim(adjustl(irrep_unitary_name_list(irrep_coirrep_relation(1,i)))),&
        !                                 ' + ',trim(adjustl(irrep_unitary_name_list(irrep_coirrep_relation(2,i))))
        !     endif

        ! enddo


    endsubroutine

    subroutine find_relation_between_unitary_irrep_coirrep(num_litt_group_unitary,irrep_num,ch_table,irrep_unitary_num,ch_unitary_table,irrep_coirrep_relation,&
                                                            discriminant_value,torsion)
        integer, intent(in) :: num_litt_group_unitary
        integer, intent(in) :: irrep_num
        integer, intent(in) :: irrep_unitary_num
        complex(dp), intent(in) :: ch_table(irrep_num,num_litt_group_unitary)
        complex(dp), intent(in) :: ch_unitary_table(irrep_unitary_num,num_litt_group_unitary)
        integer, intent(out) :: irrep_coirrep_relation(2,irrep_num)
        integer, intent(out) :: discriminant_value(irrep_unitary_num)
        integer, intent(in) :: torsion(irrep_num)
        
        integer :: i,j1,j2,k

        logical :: flag1,flag2
        
        irrep_coirrep_relation = 0
        do i=1,irrep_num
            flag1 = .false.
            do j1=1,irrep_unitary_num
                flag2 = .true.
                
                do k=1,num_litt_group_unitary
                    if (abs(ch_unitary_table(j1,k) - ch_table(i,k))>1e-4) then
                        flag2 = .false.
                        exit
                    endif

                enddo

                if (flag2) then
                    flag1 = .true.
                    exit
                endif
            enddo
            
            if (flag1) then
                irrep_coirrep_relation(1,i) = j1
                discriminant_value(j1) = torsion(i)
            else
            two:do j1=1,irrep_unitary_num
                    do j2=1,irrep_unitary_num
                        flag2 = .true.
                        do k=1,num_litt_group_unitary
                            if (abs(ch_unitary_table(j1,k) + ch_unitary_table(j2,k) - ch_table(i,k))>1e-4) then
                                flag2 = .false.
                                exit
                            endif
                        enddo

                        if (flag2) then
                            flag1 = .true.
                            exit two
                        endif
                    enddo
                enddo two
                if (flag1) then
                    irrep_coirrep_relation(1,i) = j1
                    irrep_coirrep_relation(2,i) = j2
                    discriminant_value(j1) = torsion(i)
                    discriminant_value(j2) = torsion(i)
                else
                    write(*,*)'ERROR!!! Reduce failed!!!'
                endif
            endif
        enddo
    endsubroutine


    subroutine get_discriminant_value(num_litt_group,num_litt_group_unitary,ch_unitary_table,irrep_unitary_num,&
                                        time_reversal,rot,tau,SU2,&
                                        rot_unitary,tau_unitary,SU2_unitary,discriminant_value)
                
        integer, intent(in) :: num_litt_group
        integer, intent(in) :: num_litt_group_unitary
        integer, intent(in) :: irrep_unitary_num
        complex(dp), intent(in) :: ch_unitary_table(irrep_unitary_num,num_litt_group_unitary)
        integer, intent(in) :: time_reversal(num_litt_group)
        integer, intent(in) :: rot(3,3,num_litt_group)
        integer, intent(in) :: rot_unitary(3,3,num_litt_group_unitary)
        real(dp), intent(in) :: tau(3,num_litt_group)
        real(dp), intent(in) :: tau_unitary(3,num_litt_group_unitary)
        complex(dp), intent(in) :: SU2(2,2,num_litt_group)
        complex(dp), intent(in) :: SU2_unitary(2,2,num_litt_group_unitary)
        integer, intent(out) :: discriminant_value(irrep_unitary_num)

        integer :: i,j
        complex(dp) :: sum_ch_2
        integer :: rot_tmp(3,3)
        complex(dp) :: SU2_tmp(2,2)
        real(dp) :: tau_tmp(3)
        integer :: index

        integer :: time_reversal_tmp(num_litt_group_unitary)
        intrinsic :: nint

        time_reversal_tmp = 1

        do i=1,irrep_unitary_num
            sum_ch_2 = cmplx(0.0_dp,0.0_dp,dp)
            do j=1,num_litt_group
                if (time_reversal(j) == -1) then
                    rot_tmp = matmul(rot(:,:,j),rot(:,:,j))
                    SU2_tmp = matmul(SU2(:,:,j),conjg(SU2(:,:,j)))
                    tau_tmp = tau(:,j) + matmul(rot(:,:,j),tau(:,j))
                    call find_index(num_litt_group_unitary,rot_tmp, tau_tmp, SU2_tmp, 1, rot_unitary, tau_unitary, SU2_unitary, time_reversal_tmp, index)
                    sum_ch_2 = sum_ch_2 + sign(1,index)*ch_unitary_table(i,abs(index))
                    
                endif
                
            enddo

            sum_ch_2 = sum_ch_2/num_litt_group_unitary

            discriminant_value(i) = nint(real(sum_ch_2))

            if (abs(sum_ch_2-discriminant_value(i))>1e-2) then
                write(*,*)'WARNING!!! Reduce may fail!!!'
            endif
        enddo

    endsubroutine

    logical function is_less(ii, jj, irrep_num, num_litt_group_unitary, char_now_mat)  result(less)
    implicit none
    integer, intent(in) :: ii, jj
    integer, intent(in) :: irrep_num, num_litt_group_unitary
    complex(dp), intent(in) :: char_now_mat(irrep_num,num_litt_group_unitary)
    integer :: b
    real(dp), parameter :: tol_r = 1e-2_dp, tol_i = 1e-3_dp

    less = .false.

    if     ( real(char_now_mat(ii,1)) < real(char_now_mat(jj,1)) - tol_r ) then
        less = .true.; return
    elseif ( real(char_now_mat(ii,1)) > real(char_now_mat(jj,1)) + tol_r ) then
        return                   
    end if


    do b = 2, num_litt_group_unitary

        if     ( real(char_now_mat(ii,b)) > real(char_now_mat(jj,b)) + tol_r ) then
            less = .true.; return
        elseif ( real(char_now_mat(ii,b)) < real(char_now_mat(jj,b)) - tol_r ) then
            return
        end if

        if     ( aimag(char_now_mat(ii,b)) > aimag(char_now_mat(jj,b)) + tol_i ) then
            less = .true.; return
        elseif ( aimag(char_now_mat(ii,b)) < aimag(char_now_mat(jj,b)) - tol_i ) then
            return
        end if
    end do

    stop 'error character table sorting: should not reach here!'

    end function is_less

    subroutine sort_perm(irrep_num, num_litt_group, char_now_mat, perm)
        implicit none
        integer :: i, j, key
        integer, intent(in) :: irrep_num, num_litt_group
        complex(dp), intent(in) :: char_now_mat(irrep_num,num_litt_group)
        integer, intent(out) :: perm(irrep_num)
        perm = [(i, i = 1, irrep_num)]
        
        do i = 2, irrep_num
            key = perm(i)
            j = i - 1
            do while (j >= 1)
                if (is_less(key, perm(j), irrep_num, num_litt_group, char_now_mat)) then
                    perm(j+1) = perm(j)
                    j = j - 1
                else
                    exit
                end if
            end do
            perm(j+1) = key
        end do
    end subroutine sort_perm
end module linear_rep
