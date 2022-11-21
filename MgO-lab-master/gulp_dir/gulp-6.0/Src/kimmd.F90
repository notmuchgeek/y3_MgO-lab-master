  subroutine kimmd(ekim,esregion12,esregion2,eattach,lgrad1)
!
!  Calculates the energy and up to first derivatives for OpenKIM potentials.
!
!  On entry : 
!
!  nkim            = integer reference to number of KIM model for this call
!  lgrad1          = if .true. calculate the first derivatives
!
!  On exit :
!
!  ekim            = the value of the energy contribution from the nkim'th model
!  eattach         = sum of particle energies for the growth slice (if applicable)
!  eregion12       = currently unchanged since we can't separate region 1-2 contributions
!  eregion2        = sum of particle energies for region 2 (if applicable)
!
!  10/12 Created from bondordermd.f90
!   4/14 Modifications for OpenKIM version 1.4 made
!   7/14 Modifications for OpenKIM version 1.7.3 made
!   8/16 Forces now subtracted to change sign to be derivative
!   2/18 Trace added
!   8/18 Modified for version 2 of OpenKIM
!   9/18 Particle energies added
!   9/18 Region and attachment energies added as arguments
!   9/19 nrel2 / nrelat changed to nrela2f / nrelf2a
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
!  Copyright Curtin University 2019
!
!  Julian Gale, CIC, Curtin University, September 2019
!
#ifdef KIM
  use kim_compute_arguments_module, only : kim_get_argument_support_status, &
                                           kim_set_argument_pointer, &
                                           kim_set_callback_pointer
  use kim_support_status_module,    only : kim_support_status_type
  use kim_model_module,             only : kim_compute, &
                                           kim_compute_arguments_create, &
                                           kim_compute_arguments_destroy
#endif
  use datatypes
  use current
  use iochannels
#ifdef KIM
  use configurations, only : nregions, nregionno, lsliceatom
  use control,        only : lsiteenergy
  use derivatives,    only : xdrv, ydrv, zdrv, rstrd
  use energies,       only : siteenergy, eregion2region
  use kim_cellimages, only : kim_nicell
  use kim_functions,  only : get_neigh
  use kim_functions,  only : set_kim_neighbours
  use kim_functions,  only : set_kim_neighbour_list
  use kim_models
  use symmetry,       only : lstr
#endif
  use parallel
  use reallocate
  use times
#ifdef TRACE
  use trace,          only : trace_in, trace_out
#endif
  implicit none
!
!  Passed variables
!
  real(dp),      intent(out)                       :: ekim               ! KIM model energy on return
  real(dp),      intent(inout)                     :: esregion12         ! Region 1-2 interaction energy
  real(dp),      intent(inout)                     :: esregion2          ! Region 2-2 interaction energy
  real(dp),      intent(inout)                     :: eattach            ! Energy of growth slice
  logical,       intent(in)                        :: lgrad1             ! If true, then compute the forces
#ifdef KIM
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ii
  integer(i4)                                      :: imid
  integer(i4)                                      :: KIMerror           ! KIM error flag
  integer(i4)                                      :: nm
  integer(i4)                                      :: nregioni
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: conversion
  real(dp)                                         :: e2local
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp),                                pointer :: noforces => null()
  real(dp),                                pointer :: noparticle => null()
  type(kim_support_status_type)                    :: support_status
#endif
#ifdef TRACE
  call trace_in('kimmd')
#endif
!
#ifdef KIM
  t1 = g_cpu_time()
!
!  Initialise energy
!
  ekim = 0.0_dp
!
!  Call neighbour list set up for KIM potentials - common to all models
!
  call set_kim_neighbours
!*********************
!  Loop over models  *
!*********************
  do nm = 1,nkimmodel
    if (lkim_model_cfg_OK(ncf,nm)) then
!
!  Create compute arguments
!
      call kim_compute_arguments_create(kim_model(nm),kim_compute_arguments(nm),KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error creating KIM compute argument',0_i4)
        call stopnow('kimmd')
      endif
!
!  Check whether model supports forces/virial and whether they are required
!
      call kim_get_argument_support_status(kim_compute_arguments(nm), &
        kim_compute_argument_name_partial_forces,support_status,KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error checking status of KIM argument - forces',0_i4)
        call stopnow('kimmd')
      endif
      if (support_status.eq.kim_support_status_not_supported) then
        kim_forcesOK = .false.
        kim_forcesRequired = .false.
      else
        kim_forcesOK = .true.
        if (support_status.eq.kim_support_status_required) then
          kim_forcesRequired = .true.
        else
          kim_forcesRequired = .false.
        endif
      endif
!
      call kim_get_argument_support_status(kim_compute_arguments(nm), &
        kim_compute_argument_name_partial_virial,support_status,KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error checking status of KIM argument - virial',0_i4)
        call stopnow('kimmd')
      endif
      if (support_status.eq.kim_support_status_not_supported) then
        kim_virialOK = .false.
        kim_virialRequired = .false.
      else
        kim_virialOK = .true.
        if (support_status.eq.kim_support_status_required) then
          kim_virialRequired = .true.
        else
          kim_virialRequired = .false.
        endif
      endif
!
!  Check whether model is compatible with force/virial requests from GULP
!
      if (lgrad1.and..not.kim_forcesOK) then
        call outerror('KIM model does not support forces',0_i4)
        call stopnow('kimmd')
      endif
      if (kim_virialNeeded.and..not.kim_virialOK) then
        call outerror('KIM model does not support virial for strains',0_i4)
        call stopnow('kimmd')
      endif
!
!  Check whether model supports particle energies and whether they are required
!
      call kim_get_argument_support_status(kim_compute_arguments(nm), &
        kim_compute_argument_name_partial_particle_energy,support_status,KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error checking status of KIM argument - particle energy',0_i4)
        call stopnow('kimmd')
      endif
      if (support_status.eq.kim_support_status_not_supported) then
        kim_particleOK = .false.
        kim_particleRequired = .false.
      else
        kim_particleOK = .true.
        if (support_status.eq.kim_support_status_required) then
          kim_particleRequired = .true.
        else
          kim_particleRequired = .false.
        endif
      endif
!
!  Check whether model is compatible with particle energy requests from GULP
!
      if (kim_particleNeeded.and..not.kim_particleOK) then
        call outerror('KIM model does not support particle energies',0_i4)
        call stopnow('kimmd')
      endif
!
!  Set total number of particles needed for OpenKIM
!
      kim_numat = numat*kim_nicell
!
!  Reallocate arrays for data
!
      call realloc(kim_nspecies,kim_numat,KIMerror)
      if (KIMerror.ne.0) call outofmemory('setkim','kim_nspecies')
      call realloc(kim_ncontributing,kim_numat,kim_num_nlist(nm),KIMerror)
      if (KIMerror.ne.0) call outofmemory('setkim','kim_ncontributing')
      call realloc(kim_coord,kim_ndim,kim_numat,KIMerror)
      if (KIMerror.ne.0) call outofmemory('setkim','kim_coord')
      call realloc(kim_forces,kim_ndim,kim_numat,KIMerror)
      if (KIMerror.ne.0) call outofmemory('setkim','kim_forces')
      call realloc(kim_particle,kim_numat,KIMerror)
      if (KIMerror.ne.0) call outofmemory('setkim','kim_particle')
!
      call realloc(kim_neighbourlist,kim_maxneighbour,kim_numat,kim_num_nlist(nm),KIMerror)
      if (KIMerror.ne.0) call outofmemory('setkim','kim_neighbourlist')
!
!  Associate memory for KIM
!
      call kim_set_argument_pointer(kim_compute_arguments(nm), &
        kim_compute_argument_name_number_of_particles,kim_numat,KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error associating KIM pointer - kim_numat',0_i4)
        call stopnow('kimmd')
      endif
!
      call kim_set_argument_pointer(kim_compute_arguments(nm), &
        kim_compute_argument_name_particle_species_codes,kim_nspecies, &
        KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error associating KIM pointer - kim_nspecies',0_i4)
        call stopnow('kimmd')
      endif
!
      call kim_set_argument_pointer(kim_compute_arguments(nm), &
        kim_compute_argument_name_particle_contributing,kim_ncontributing, &
        KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error associating KIM pointer - kim_ncontributing',0_i4)
        call stopnow('kimmd')
      endif
!
      call kim_set_argument_pointer(kim_compute_arguments(nm), &
        kim_compute_argument_name_coordinates,kim_coord,KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error associating KIM pointer - kim_coord',0_i4)
        call stopnow('kimmd')
      endif
!
      call kim_set_argument_pointer(kim_compute_arguments(nm), &
        kim_compute_argument_name_partial_energy,kim_energy,KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error associating KIM pointer - kim_energy',0_i4)
        call stopnow('kimmd')
      endif
!
      if (lstr.and.kim_virialNeeded) then
        call kim_set_argument_pointer(kim_compute_arguments(nm), &
          kim_compute_argument_name_partial_virial,kim_virial,KIMerror)
        if (KIMerror.ne.0) then
          call outerror('error associating KIM pointer - kim_virial',0_i4)
          call stopnow('kimmd')
        endif
      else
!
!  If virial is not needed then pass a null pointer unless this is a required argument
!
        if (kim_virialRequired) then
          call kim_set_argument_pointer(kim_compute_arguments(nm), &
            kim_compute_argument_name_partial_virial,kim_virial,KIMerror)
          if (KIMerror.ne.0) then
            call outerror('error associating KIM pointer - kim_virial',0_i4)
            call stopnow('kimmd')
          endif
        elseif (kim_virialOK) then
          call kim_set_argument_pointer(kim_compute_arguments(nm), &
            kim_compute_argument_name_partial_virial,noforces,KIMerror)
          if (KIMerror.ne.0) then
            call outerror('error associating KIM pointer - kim_virial',0_i4)
            call stopnow('kimmd')
          endif
        endif
      endif
!
      if (lgrad1) then
        call kim_set_argument_pointer(kim_compute_arguments(nm), &
          kim_compute_argument_name_partial_forces,kim_forces,KIMerror)
        if (KIMerror.ne.0) then
          call outerror('error associating KIM pointer - kim_forces',0_i4)
          call stopnow('kimmd')
        endif
      else
!
!  If forces are not needed then pass a null pointer unless this is a required argument
!
        if (kim_forcesRequired) then
          call kim_set_argument_pointer(kim_compute_arguments(nm), &
            kim_compute_argument_name_partial_forces,kim_forces,KIMerror)
          if (KIMerror.ne.0) then
            call outerror('error associating KIM pointer - kim_forces',0_i4)
            call stopnow('kimmd')
          endif
        elseif (kim_forcesOK) then
          call kim_set_argument_pointer(kim_compute_arguments(nm), &
            kim_compute_argument_name_partial_forces,noforces,KIMerror)
          if (KIMerror.ne.0) then
            call outerror('error associating KIM pointer - kim_forces',0_i4)
            call stopnow('kimmd')
          endif
        endif
      endif
!
!  Particle energies
!
      if (lsiteenergy.or.nregions(ncf).gt.1) then
        call kim_set_argument_pointer(kim_compute_arguments(nm), &
          kim_compute_argument_name_partial_particle_energy,kim_particle,KIMerror)
        if (KIMerror.ne.0) then
          call outerror('error associating KIM pointer - kim_particle',0_i4)
          call stopnow('kimmd')
        endif
      else
!
!  If particle energies are not needed then pass a null pointer unless this is a required argument
!
        if (kim_particleRequired) then
          call kim_set_argument_pointer(kim_compute_arguments(nm), &
            kim_compute_argument_name_partial_particle_energy,kim_particle,KIMerror)
          if (KIMerror.ne.0) then
            call outerror('error associating KIM pointer - kim_particle',0_i4)
            call stopnow('kimmd')
          endif
        elseif (kim_particleOK) then
          call kim_set_argument_pointer(kim_compute_arguments(nm), &
            kim_compute_argument_name_partial_particle_energy,noparticle,KIMerror)
          if (KIMerror.ne.0) then
            call outerror('error associating KIM pointer - kim_particle',0_i4)
            call stopnow('kimmd')
          endif
        endif
      endif
!
!  Populate KIM arrays with coordinates and atom properties
!
      call set_kim_neighbour_list(nm)
!
!  Set pointer in KIM object to neighbour list routine and object
!
      call kim_set_callback_pointer(kim_compute_arguments(nm), &
        kim_compute_callback_name_get_neighbor_list,kim_language_name_fortran, &
        c_funloc(get_neigh),c_loc(kim_neighbourlist),KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error associating KIM pointer - get_neigh',0_i4)
        call stopnow('kimmd')
      endif
!*********************************************
!  Compute KIM energy and optionally forces  *
!*********************************************
!
!  Call KIM model compute
!
      call kim_compute(kim_model(nm),kim_compute_arguments(nm),KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error in KIM model compute',0_i4)
        call stopnow('kimmd')
      endif
!****************************************
!  Copy energy, forces and virial back  *
!****************************************
      ekim = ekim + kim_energy*kim_energy_conversion(nm)
!
!  Particle energies (optional)
!
      if (lsiteenergy) then
        imid = 0
        do ii = 1,kim_nicell
          do i = 1,numat
            siteenergy(i) = siteenergy(i) + kim_particle(imid+i)*kim_energy_conversion(nm)
          enddo
          imid = imid + numat
        enddo
      endif
!
!  Region contributions
!
      if (nregions(ncf).gt.1) then
        imid = 0
        e2local = 0.0_dp
        do ii = 1,kim_nicell
          do i = 1,numat
            if (lsliceatom(nsft + nrelf2a(i))) then
              eattach = eattach + kim_particle(imid+i)*kim_energy_conversion(nm)
            endif
            nregioni = nregionno(nsft+nrelf2a(i))
            if (nregioni.ne.1) then
              e2local = e2local + kim_particle(imid+i)*kim_energy_conversion(nm)
            endif
            eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + &
              kim_particle(imid+i)*kim_energy_conversion(nm)
          enddo
          imid = imid + numat
        enddo
        esregion2 = esregion2 + e2local
!
!  Correction region 1 energy by subtracting region 2 contribution
!
        ekim = ekim - e2local
      endif
!
      if (lgrad1) then
        conversion = kim_energy_conversion(nm)/kim_length_conversion(nm)
        imid = 0
        do ii = 1,kim_nicell
          do i = 1,numat
            xdrv(i) = xdrv(i) - kim_forces(1,imid+i)*conversion
            ydrv(i) = ydrv(i) - kim_forces(2,imid+i)*conversion
            zdrv(i) = zdrv(i) - kim_forces(3,imid+i)*conversion
          enddo
          imid = imid + numat
        enddo
        if (lstr) then
          select case(ndim)
            case(1)
              rstrd(1) = rstrd(1) + kim_virial(1)*kim_energy_conversion(nm)
            case(2)
              rstrd(1) = rstrd(1) + kim_virial(1)*kim_energy_conversion(nm)
              rstrd(2) = rstrd(2) + kim_virial(2)*kim_energy_conversion(nm)
              rstrd(3) = rstrd(3) + kim_virial(3)*kim_energy_conversion(nm)
            case(3)
              rstrd(1) = rstrd(1) + kim_virial(1)*kim_energy_conversion(nm)
              rstrd(2) = rstrd(2) + kim_virial(2)*kim_energy_conversion(nm)
              rstrd(3) = rstrd(3) + kim_virial(3)*kim_energy_conversion(nm)
              rstrd(4) = rstrd(4) + kim_virial(4)*kim_energy_conversion(nm)
              rstrd(5) = rstrd(5) + kim_virial(5)*kim_energy_conversion(nm)
              rstrd(6) = rstrd(6) + kim_virial(6)*kim_energy_conversion(nm)
          end select
        endif
      endif
!
!  Destroy compute arguments
!
      call kim_compute_arguments_destroy(kim_model(nm),kim_compute_arguments(nm),KIMerror)
      if (KIMerror.ne.0) then
        call outerror('error destroying KIM compute argument',0_i4)
        call stopnow('kimmd')
      endif
    endif
!
!  End of loop over models
!
  enddo
!
!  NB: Need to handle regions & eattach!
!
  t2 = g_cpu_time()
  tkim = tkim + t2 - t1
#else
  ekim = 0.0_dp
#endif
#ifdef TRACE
  call trace_out('kimmd')
#endif
!
  return
  end
