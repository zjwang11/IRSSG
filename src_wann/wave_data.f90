module wave_data_wann

    use comms_wann
    implicit none

    complex(dp), allocatable, save :: Hk(:,:)
    complex(dp), allocatable, save :: evec(:,:)

contains

subroutine get_WF()

    integer :: i, info 

    allocate(  Hk(num_wann,num_wann))
    allocate(evec(num_wann,num_wann))

    info = 0
    call get_Hk()
    evec = czero 
    call myzheev(num_wann, num_wann, Hk, EE, evec ,info)
    coeffa = czero; coeffb = czero
    coeffa(1:num_sumorb, 1:num_wann) = evec(1:num_sumorb, 1:num_wann)
    coeffb(1:num_sumorb, 1:num_wann) = evec(num_sumorb+1:2*num_sumorb, 1:num_wann)

    deallocate(Hk)
    deallocate(evec)
end subroutine get_WF 


subroutine get_WF_spin_polarized(spin_case)

    integer :: i, info 
    integer,intent(in) :: spin_case ! 1 up; 2 down; 3 all

    if (spin_case <= 2) then
        allocate(  Hk(num_wann,num_wann))
        allocate(evec(num_wann,num_wann))
    else 
        allocate(  Hk(num_wann*2,num_wann*2))
        allocate(evec(num_wann*2,num_wann*2))    
    endif

    info = 0

    call get_Hk_spin_polarized(spin_case)
    
    evec = czero 

    if (spin_case <= 2) then
        call myzheev(num_wann, num_wann, Hk, EE, evec ,info)
        coeffa = czero; coeffb = czero
        write(*,*)num_sumorb,num_wann
        coeffa(1:num_sumorb, 1:num_wann) = evec(1:num_sumorb, 1:num_wann)

    else
        call myzheev(num_wann*2, num_wann*2, Hk, EE, evec ,info)

        coeffa = czero; coeffb = czero
        coeffa(1:num_sumorb, 1:num_wann*2) = evec(1:num_sumorb, 1:num_wann*2)
        coeffb(1:num_sumorb, 1:num_wann*2) = evec(num_sumorb+1:2*num_sumorb, 1:num_wann*2)
        ! call bubble_sort_with_index(num_bands*2,EE,order)
    endif
    
    

    deallocate(Hk)
    deallocate(evec)
end subroutine get_WF_spin_polarized

subroutine get_Hk()

    integer  :: i
    real(dp) :: kdotr 

    Hk = czero 
    do i = 1, num_block
        kdotr = dot_product(WK(:), vec_block(:,i))
        Hk(1:num_wann,1:num_wann) = Hk(1:num_wann,1:num_wann) + &
                HmnR(1:num_wann, 1:num_wann, i)*exp(CI*kdotr*2.d0*PI)
    enddo 
    ! call reorder_oddeven(Hk,num_wann)
end subroutine get_Hk 

subroutine get_Hk_spin_polarized(spin_case)

    integer  :: i
    real(dp) :: kdotr 
    integer,intent(in) :: spin_case ! 1 up; 2 down; 3 all


    Hk = czero 
    do i = 1, num_block
        kdotr = dot_product(WK(:), vec_block(:,i))
        if (spin_case == 1) then
            Hk(1:num_wann,1:num_wann) = Hk(1:num_wann,1:num_wann) + &
                    HmnR_up(1:num_wann, 1:num_wann, i)*exp(CI*kdotr*2.d0*PI)
        elseif (spin_case == 2) then
            Hk(1:num_wann,1:num_wann) = Hk(1:num_wann,1:num_wann) + &
                        HmnR_dn(1:num_wann, 1:num_wann, i)*exp(CI*kdotr*2.d0*PI)
        else
            Hk(1:num_wann,1:num_wann) = Hk(1:num_wann,1:num_wann) + &
                        HmnR_up(1:num_wann, 1:num_wann, i)*exp(CI*kdotr*2.d0*PI)
            Hk(num_wann+1:2*num_wann,num_wann+1:2*num_wann) = Hk(num_wann+1:2*num_wann,num_wann+1:2*num_wann) + &
                        HmnR_dn(1:num_wann, 1:num_wann, i)*exp(CI*kdotr*2.d0*PI)
        endif
    enddo 


end subroutine get_Hk_spin_polarized 

subroutine myzheev(lda, N, Hmn, eval, evec, info)

    integer,     intent(in)  :: lda
    integer,     intent(in)  :: N
    complex(dp), intent(in)  :: Hmn(lda, N)
    real(dp),    intent(out) :: eval(N)
    complex(dp), intent(out) :: evec(lda, N)
    integer,     intent(out) :: info 

    ! dummy integer variables for lapack call
    integer     :: LWORK
    real(dp)    :: RWORK(3*N)
    COMPLEX(dp) :: WORK(2*N)
    complex(dp) :: A(N,N)

    LWORK = 2*N

    A      = czero
    A(:,:) = Hmn(1:N,1:N)
    call ZHEEV('V','U',N,A,N,eval,WORK,LWORK,RWORK,INFO)
    evec(1:N,1:N) = A(:,:)

    return 

end subroutine myzheev 


subroutine bubble_sort_with_index(nband, A, idx)
    integer, intent(in) :: nband
      real(dp)   , intent(inout) :: A(nband)
      integer, intent(out) :: idx(nband)

      integer :: n, i, j, tmp_i
      real    :: tmp_r

      n = size(A)
      idx = [(i, i=1, n)]  

      do i = 1, n-1
         do j = 1, n-i
            if (A(idx(j)) > A(idx(j+1))) then
               tmp_i      = idx(j)
               idx(j)     = idx(j+1)
               idx(j+1)   = tmp_i
            end if
         end do
      end do
      A = A(idx)
   end subroutine bubble_sort_with_index


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

end module wave_data_wann 
