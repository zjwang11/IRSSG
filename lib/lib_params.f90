module lib_params
    ! Basic parameters used by library modules
    implicit none 

    integer,     public, parameter :: dp = 8
    
    real(dp),    public, parameter :: pi = 3.141592653589793238462643383279d0

    real(dp),    public            :: lattice(3,3)

end module lib_params

