!************************************
!  Module for GULP-KIM interaction  *
!************************************
!
!  10/12 Created from modules.f90
!   4/14 Modified for OpenKIM v1.4
!   7/14 pkim_model now kept as integer pointer
!   7/16 pkim_model changed to an array of type c_ptr
!   8/18 Re-written for OpenKIM v2.0
!   8/18 Datatypes standardised so that for OpenKIM builds
!        the default types are C compatible
!   9/18 Particle energies added
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, September 2018
!

!
!  OpenKIM models
!
  module kim_models
    use datatypes
#ifdef KIM
    use kim_simulator_headers_module
#endif

    logical                                      :: lkim_model = .false.         ! If true then OpenKIM model to be used
    integer(i4),      parameter                  :: kim_len = 10000              ! Length of KIM descriptor string
    integer(i4),      parameter                  :: maxkimmodel = 10             ! Maximum number of KIM models
    integer(i4)                                  :: nkimcutoffs = 0              ! Number of KIM cutoffs
    integer(i4)                                  :: nkimmodel = 0                ! Number of KIM models
    integer(i4)                                  :: nlibnkimmodel = 0            ! Number of KIM models prior to libraries
    character(len=256,kind=c_char)               :: kim_model_name(maxkimmodel)  ! Names of KIM models
!
#ifdef KIM
    type(kim_model_handle_type)                  :: kim_model(maxkimmodel)       ! Model handle
    type(kim_compute_arguments_handle_type)      :: kim_compute_arguments(maxkimmodel) ! Compute arguments handle
#endif
!
    integer(i4)                                  :: kim_energy_units(maxkimmodel) ! Integer indicating energy units
    integer(i4)                                  :: kim_length_units(maxkimmodel) ! Integer indicating length units
    integer(i4),                     parameter   :: kim_ndim = 3                 ! Integer containing number of dimensions
    integer(i4)                                  :: kim_numat                    ! Integer containing number of atoms for KIM
    integer(i4)                                  :: kim_maxneighbour             ! Maximum number of neighbours
    integer(i4)                                  :: kim_num_nlist(maxkimmodel)   ! Number of neighbour lists
    integer(i4), dimension(:,:),   pointer       :: kim_ncontributing_only       ! Flag to whether neighbours of non-contributing are needed
    integer(i4), dimension(:,:,:), pointer       :: kim_neighbourlist            ! Neighbour list in KIM form
    integer(i4), dimension(:,:),   pointer       :: kim_nspec                    ! Array to hold species types for KIM per species
    integer(i4), dimension(:),     pointer       :: kim_nspecies                 ! Array to hold species types for KIM
    integer(i4), dimension(:,:),   pointer       :: kim_ncontributing            ! Array for flag for contributing particles
!
    real(dp),    dimension(:,:), pointer         :: kim_coord  => null()         ! Array to hold coordinates for KIM
    real(dp),                             target :: kim_energy                   ! Variable to hold KIM energy
    real(dp),    dimension(:,:), pointer         :: kim_forces => null()         ! Array to hold coordinates for KIM
    real(dp),    dimension(:),   pointer         :: kim_particle => null()       ! Array to hold particle energies for KIM
    real(dp),                             target :: kim_virial(6)                ! Variable to hold KIM virial
!
    real(dp)                                     :: kim_cutoff                   ! Variable to hold cutoff for KIM models
    real(dp),    dimension(:,:), pointer         :: kim_cutoffs => null()        ! Array to hold cutoffs for KIM models
    real(dp)                                     :: kim_influence                ! Variable to hold overall influence distance
    real(dp)                                     :: kim_influences(maxkimmodel)  ! Variable to hold influence distance
    real(dp)                                     :: kim_energy_conversion(maxkimmodel) ! Conversion factor for energy
    real(dp)                                     :: kim_length_conversion(maxkimmodel) ! Conversion factor for length
!
    logical                                      :: kim_any_ncontributing_only(maxkimmodel) ! If true then at least one list needs neighbours of non-contributing particles
    logical                                      :: kim_forcesOK                 ! If true then model supports forces
    logical                                      :: kim_forcesRequired           ! If true then forces must be requested from the model
    logical                                      :: kim_particleOK               ! If true then model supports particle energies
    logical                                      :: kim_particleNeeded           ! If true then GULP needs the particle energies
    logical                                      :: kim_particleRequired         ! If true then particle energies must be requested from the model
    logical                                      :: kim_virialOK                 ! If true then model supports virial
    logical                                      :: kim_virialNeeded             ! If true then GULP needs the virial
    logical                                      :: kim_virialRequired           ! If true then virial must be requested from the model
    logical,     dimension(:,:), pointer         :: lkim_model_cfg_OK => null()  ! Stores logical that configuration is OK for KIM model
!
  end module kim_models
!
!  OpenKIM module for cell images within influence distance
!
  module kim_cellimages
    use datatypes
    integer(i4)                                  :: kim_maxncell = 0             ! Maximum number of cell images - array dimension
    integer(i4)                                  :: kim_midcell                  ! Number of central cell for cutoff
    integer(i4)                                  :: kim_midicell                 ! Number of central cell for influence
    integer(i4)                                  :: kim_ncell                    ! Total number of cell images for cutoff
    integer(i4)                                  :: kim_ncells(3)                ! Number of cells per dimension for cutoff
    integer(i4)                                  :: kim_nicell                   ! Total number of cell images for influence
    integer(i4)                                  :: kim_nicells(3)               ! Number of cells per dimension for influence
    integer(i4), dimension(:,:),        pointer  :: kim_icell => null()          ! Cell indices for vector
    real(dp),    dimension(:,:),        pointer  :: kim_cellvector => null()     ! Vectors for cell images
  end module kim_cellimages
