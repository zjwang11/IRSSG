module fixed_subspace_mod
    use lib_params, only: dp
    implicit none
    public :: fixed_subspace
contains

    subroutine fixed_subspace(R, n, tol, dim, basis)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: R(3,3,n)
        real(dp), intent(in) :: tol
        integer, intent(out) :: dim
        real(dp), intent(out) :: basis(3,3)

        real(dp), allocatable :: A(:,:), U(:,:), VT(:,:), work(:)
        real(dp) :: S(3)
        integer :: i, j, k, info, lwork, m
        external :: dgesvd

        if (n <= 0) then
            dim = 3
            basis = 0.0_dp
            basis(1,1) = 1.0_dp
            basis(2,2) = 1.0_dp
            basis(3,3) = 1.0_dp
            return
        endif

        m = 3 * n
        allocate(A(m,3))
        A = 0.0_dp

        do i = 1, n
            do j = 1, 3
                do k = 1, 3
                    A(3*(i-1)+j, k) = R(j,k,i)
                end do
                A(3*(i-1)+j, j) = A(3*(i-1)+j, j) - 1.0_dp
            end do
        end do

        allocate(U(m,m), VT(3,3))
        allocate(work(1))
        call dgesvd('A','A',m,3,A,m,S,U,m,VT,3,work,-1,info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(max(1,lwork)))

        call dgesvd('A','A',m,3,A,m,S,U,m,VT,3,work,lwork,info)
        if (info /= 0) then
            dim = 0
            basis = 0.0_dp
            deallocate(A, U, VT, work)
            return
        endif

        dim = 0
        basis = 0.0_dp
        do i = 1, 3
            if (S(i) < tol) then
                dim = dim + 1
                basis(:,dim) = VT(i,:)
            endif
        end do

        deallocate(A, U, VT, work)
    end subroutine fixed_subspace

end module fixed_subspace_mod
