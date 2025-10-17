module comms  
    use lib_params, only: dp, pi, lattice
    implicit none

    integer,     public, parameter :: MAXSYM = 96

    integer,     public, parameter :: MAXIR = 48

    integer,     public, parameter :: MAXIRDG = 4

    integer,     public, parameter :: MAXKTYPE = 40


    logical,     public            :: isComplexWF
    logical,     public            :: isSpinor
    logical,     public            :: isInv
    logical,     public            :: isSpinPola
    logical,     public            :: isNoncollinear

    real(dp),    public            :: br2(3,3), br4(3,3)

    integer,     public            :: num_sym
    
    integer,     public            :: num_k

    integer,     public            :: num_bands

    integer,     public            :: bot_band

    integer,     public            :: top_band
    
    integer,     public            :: nspin 

    logical,     public            :: spinor

    integer,     public            :: max_plane
    
    integer,    allocatable, save :: rot_input(:,:,:)
    real(dp),    allocatable, save :: tau_input(:,:)
    complex(dp),    allocatable, save :: SU2_input(:,:,:)
    integer,    allocatable, save :: time_reversal_input(:)
    real(dp),    allocatable, save :: spin_rot_input(:,:,:)

    ! from WAVECAR
    integer    , allocatable, save, public ::  igall(:,:) 
    real(dp)   , allocatable, save, public ::  KV(:,:) 

    complex(dp), allocatable, save, public ::  coeffa(:,:),coeffb(:,:)
    real(dp)   , allocatable, save, public ::  EE(:)

    real(dp),                 save, public ::  WK(3)

    integer,                  save, public :: ncnt, nplane
    integer, allocatable, save, public :: order(:)
    real(dp), allocatable, save, public :: kpoints_from_outcar(:,:) 

end module comms

