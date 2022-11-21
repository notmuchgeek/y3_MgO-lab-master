  subroutine setkim
!
!  Initialises OpenKIM models for each configuration 
!
!  10/12 Created
!   7/14 Modifications for OpenKIM v1.4 started
!   7/16 Modifications for OpenKIM v1.7.3 added
!   8/18 Modified for version 2 of OpenKIM
!   9/18 Multiple KIM models added
!  11/18 KIM subroutine names updated due to changes in beta release
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
!  Julian Gale, Curtin University, November 2018
!
  use, intrinsic :: iso_c_binding
#ifdef KIM
  use kim_model_module,        only : kim_get_units, &
                                      kim_get_species_support_and_code, &
                                      kim_get_number_of_neighbor_lists, &
                                      kim_get_neighbor_list_values, &
                                      kim_get_influence_distance
  use kim_species_name_module, only : kim_from_string
  use kim_unit_system_module,  only : kim_length_unit_type, &
                                      kim_energy_unit_type, &
                                      kim_charge_unit_type, &
                                      kim_temperature_unit_type, &
                                      kim_time_unit_type
  use configurations, only : ncfg, ncellmaxcfg, ncellmincfg, nregions
#endif
  use control
  use current
#ifdef KIM
  use g_constants,    only : autoangs, autoev, kcaltoev, evtoj
#endif
  use iochannels
  use reallocate
  use kim_models
  use parallel
#ifdef KIM
  use species,        only : nspec, natspec
#endif
  implicit none
#ifdef KIM
!
!  Local variables
!
  character(len=2,kind=c_char)                   :: kim_string       ! Element symbol for KIM atom typing
  character(len=2)                               :: type_string      ! Element symbol for atom typing
  integer(i4)                                    :: ii
  integer(i4)                                    :: nc
  integer(i4)                                    :: nco
  integer(i4)                                    :: nm
  integer(i4)                                    :: ns
  integer(i4)                                    :: speciesOK        ! Used to check species
  integer(i4)                                    :: unitsOK          ! Used to check units
  integer(i4)                                    :: KIMerror         ! KIMerror flag
  integer(i4),                      allocatable  :: itmp(:)
  real(dp)                                       :: cut
  real(dp),                         allocatable  :: tmp(:)
  type(kim_species_name_type)                    :: kim_species_name 
  type(kim_length_unit_type)                     :: length_unit
  type(kim_energy_unit_type)                     :: energy_unit
  type(kim_charge_unit_type)                     :: charge_unit
  type(kim_temperature_unit_type)                :: temperature_unit
  type(kim_time_unit_type)                       :: time_unit
!
!  Loop over models
!
  do nm = 1,nkimmodel
!
!  Set default unit flags for KIM - only energy and length are needed for now
!
!  kim_length_units : 1 => Angstroms; 2 => Bohr; 3 => nm
!  kim_energy_units : 1 => eV; 2 => Hartree; 3 => J; 4 => kcal/mol
!
    kim_length_units = 1
    kim_energy_units = 1
    kim_length_conversion = 1.0_dp
    kim_energy_conversion = 1.0_dp
!
!  Create KIM model to test units
!
    call kim_model_create(kim_numbering_one_based, &
                          kim_length_unit_a, &
                          kim_energy_unit_ev, &
                          kim_charge_unit_e, &
                          kim_temperature_unit_k, &
                          kim_time_unit_ps, &
                          trim(kim_model_name(nm)), &
                          unitsOK, kim_model(nm), KIMerror)
    if (KIMerror.ne.0) then
      call outerror('error calling KIM model create',0_i4)
      call stopnow('setkim')
    endif
!
!  Check whether units were compatible
!
    if (unitsOK.eq.0) then
!
!  Find out what units can be used
!
      call kim_get_units(kim_model(nm),length_unit,energy_unit, &
         charge_unit,temperature_unit,time_unit)
!
      if (length_unit.eq.kim_length_unit_a) then
        kim_length_units(nm) = 1
        kim_length_conversion(nm) = 1.0_dp
      elseif (length_unit.eq.kim_length_unit_bohr) then
        kim_length_units(nm) = 2
        kim_length_conversion(nm) = autoangs
      elseif (length_unit.eq.kim_length_unit_nm) then
        kim_length_units(nm) = 3
        kim_length_conversion(nm) = 10.0_dp
      else
        call outerror('length units are incompatible with KIM model',0_i4)
        call stopnow('setkim')
      endif
!
      if (energy_unit.eq.kim_energy_unit_ev) then
        kim_energy_units(nm) = 1
        kim_energy_conversion(nm) = 1.0_dp
      elseif (energy_unit.eq.kim_energy_unit_hartree) then
        kim_energy_units(nm) = 2
        kim_energy_conversion(nm) = autoev
      elseif (energy_unit.eq.kim_energy_unit_j) then
        kim_energy_units(nm) = 3
        kim_energy_conversion(nm) = 1.0_dp/evtoj
      elseif (energy_unit.eq.kim_energy_unit_kcal_mol) then
        kim_energy_units(nm) = 4
        kim_energy_conversion(nm) = kcaltoev
      else
        call outerror('energy units are incompatible with KIM model',0_i4)
        call stopnow('setkim')
      endif
    endif
!
!  Allocate array to hold species type mapping to OpenKIM model
!
    call realloc(kim_nspec,nspec,nkimmodel,KIMerror)
    if (KIMerror.ne.0) call outofmemory('setkim','kim_nspec')
!
!  Initialise species numbers to zero to trap for unsupported species
!
    kim_nspec(1:nspec,nm) = 0
!
!  Loop over species looking for valid types in the model
!
    do ns = 1,nspec
!
!  Get element symbol
!
      call label(natspec(ns),0_i4,type_string)
      kim_string = type_string
!
!  Find KIM species
!
      call kim_from_string(kim_string, kim_species_name)
!
!  Translate to code for species in KIM
!
      call kim_get_species_support_and_code(kim_model(nm), &
        kim_species_name, speciesOK, kim_nspec(ns,nm), KIMerror)
!
!  Error check on call
!
      if (KIMerror.ne.0) then
        call outerror('error calling KIM model get species',0_i4)
        call stopnow('setkim')
      endif
!
!  If species is not supported then stop
!
      if (speciesOK.ne.1) then
        call outerror('species not supported in KIM model',0_i4)
        call stopnow('setkim')
      endif
    enddo
!
!  Get influence distance information
!
    call kim_get_influence_distance(kim_model(nm),kim_influences(nm))
!
!  Get number of neighbour lists
!
    call kim_get_number_of_neighbor_lists(kim_model(nm),kim_num_nlist(nm))
!
!  Allocate array to hold cutoffs for neighbour lists
!
    call realloc(kim_cutoffs,kim_num_nlist(nm),nkimmodel,KIMerror)
    if (KIMerror.ne.0) call outofmemory('setkim','kim_cutoffs')
    call realloc(kim_ncontributing_only,kim_num_nlist(nm),nkimmodel,KIMerror)
    if (KIMerror.ne.0) call outofmemory('setkim','kim_ncontributing_only')
!
!  Get cutoff information
!
    allocate(itmp(kim_num_nlist(nm)))
    allocate(tmp(kim_num_nlist(nm)))
!
    call kim_get_neighbor_list_values(kim_model(nm),tmp, &
      itmp,KIMerror)
    if (KIMerror.ne.0) then
      call outerror('error in KIM get neighbor list call',0_i4)
      call stopnow('setkim')
    endif
!
    kim_cutoffs(1:kim_num_nlist(nm),nm) = tmp(1:kim_num_nlist(nm))
    kim_ncontributing_only(1:kim_num_nlist(nm),nm) = itmp(1:kim_num_nlist(nm))
!
    deallocate(tmp)
    deallocate(itmp)
!
!  Set flag as to whether any neighbour list needs neighbours of non-contributing particles
!
    kim_any_ncontributing_only(nm) = .false.
    do nc = 1,kim_num_nlist(nm)
      if (kim_ncontributing_only(nc,nm).eq.0) kim_any_ncontributing_only(nm) = .true.
    enddo
  enddo
!
!  Find maximum cutoff
!
  kim_cutoff = 0.0_dp
  do nm = 1,nkimmodel
    cut = maxval(kim_cutoffs(1:kim_num_nlist(nm),nm))
    kim_cutoff = max(kim_cutoff,cut)
  enddo
  kim_influence = maxval(kim_influences(1:nkimmodel))
!
!  Set a flag as to whether the virial will be needed by GULP at any point
!
  kim_virialNeeded = (ndim.gt.0.and.(lprop.or.lstressout))
  if (.not.kim_virialNeeded) then
!
!  Check for configurations with cell parameters being optimised
!
    do ii = 1,ncfg
      nco = ncellmaxcfg(ii) - ncellmincfg(ii)
      if (nco.gt.0) kim_virialNeeded = .true.
    enddo
  endif
!
!  Set a flag as to whether particle energies will be needed by GULP at any point
!
  kim_particleNeeded = lsiteenergy
  if (.not.kim_particleNeeded) then
!
!  Check for configurations with multiple regions
!
    do ii = 1,ncfg
      if (nregions(ii).gt.1) kim_particleNeeded = .true.
    enddo
  endif
#else
  kim_cutoff    = 0.0_dp
  kim_influence = 0.0_dp
  kim_num_nlist = 0
#endif
!
  return
  end
