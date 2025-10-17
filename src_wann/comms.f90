module comms_wann  
    use lib_params, only: dp, pi, lattice
    implicit none
    
    integer,     public, parameter :: MAXSYM = 96

    complex(dp),       public, parameter :: czero = (0._dp, 0._dp)
    complex(dp),       public, parameter :: cone  = (1._dp, 0._dp)
    complex(dp),       public, parameter :: CI    = (0._dp, 1._dp)

    logical,     public            :: isComplexWF
    logical,     public            :: isSpinor
    logical,     public            :: isInv
    logical,     public            :: isSpinPola
    logical,     public            :: isNoncollinear

    real(dp),    public            :: br2(3,3), br4(3,3)

    integer,           public, save      :: sgn
    integer,           public, save      :: bot_band
    integer,           public, save      :: top_band

    character(len=22), public, save      :: casename 
    character(len=22), public, save      :: hrfile_dn
    character(len=22), public, save      :: hrfile_up
    character(len=22), public, save      :: hrfile 
    
    integer,     public            :: num_k

    integer,     public            :: num_bands
    
    integer,     public            :: nspin 

    logical,     public            :: spinor
    
    integer,    allocatable, save :: rot_input(:,:,:)
    real(dp),    allocatable, save :: tau_input(:,:)
    complex(dp),    allocatable, save :: SU2_input(:,:,:)
    integer,    allocatable, save :: time_reversal_input(:)
    real(dp),    allocatable, save :: spin_rot_input(:,:,:)

    integer,           public, save      :: nrot
    integer,           public ,save      :: ncenter
    integer,           public, save      :: orbt 
    integer,           public, save      :: spincov

    integer,           public, save      :: num_sym
    integer,           public, save      :: num_wann
    integer,           public, save      :: num_sumorb
    integer,           public, save      :: num_block

    integer,  allocatable, public, save  :: type_center(:)
    integer,  allocatable, public, save  :: norb_center(:)
    integer,  allocatable, public, save  :: startorb_center(:)
    real(dp), allocatable, public, save  :: pos_center(:,:)
    real(dp), allocatable, public, save  :: rot_orb(:,:,:,:)

    real(dp), allocatable, public, save  :: len_k(:)
    real(dp), allocatable, public, save  :: k(:,:)

    integer,  allocatable, public, save  :: deg_block(:)
    integer,  allocatable, public, save  :: vec_block(:,:)
    complex(dp), allocatable, public, save :: HmnR(:,:,:)
    complex(dp), allocatable, public, save :: HmnR_up(:,:,:)
    complex(dp), allocatable, public, save :: HmnR_dn(:,:,:)

    ! from WAVECAR
    integer    , allocatable, save, public ::  igall(:,:) 
    real(dp)   , allocatable, save, public ::  KV(:,:) 

    complex(dp), allocatable, save, public ::  coeffa(:,:),coeffb(:,:)
    real(dp)   , allocatable, save, public ::  EE(:)

    real(dp),                 save, public ::  WK(3)

    integer,                  save, public :: ncnt, nplane
    integer, allocatable, save, public :: order(:)

end module comms_wann
