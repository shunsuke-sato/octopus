
!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

module ion_interaction_m
  use comm_m
  use geometry_m
  use global_m
  use distributed_m
  use loct_math_m
  use messages_m
  use mpi_m
  use parser_m
  use periodic_copy_m
  use profiling_m
  use simul_box_m
  use species_m
  use ssys_ionic_m
  use unit_system_m

  implicit none

  private
  public ::                           &
    ion_interaction_t,                &
    ion_interaction_init,             &
    ion_interaction_end,              &
    ion_interaction_add_subsys_ionic, &
    ion_interaction_calculate,        &
    ion_interaction_test

  type ion_interaction_t
    type(ssys_ionic_t), pointer :: subsys_ionic !< Subsystems ionic term.
    FLOAT                       :: alpha
  end type ion_interaction_t

  integer, parameter ::            &
    ION_COMPONENT_REAL     = 1,    &
    ION_COMPONENT_SELF     = 2,    &
    ION_COMPONENT_FOURIER  = 3,    &
    ION_NUM_COMPONENTS     = 3
  
contains

  subroutine ion_interaction_init(this)
    type(ion_interaction_t), intent(out)   :: this

    PUSH_SUB(ion_interaction_init)

    !%Variable EwaldAlpha
    !%Type float
    !%Default 1.1313708
    !%Section Hamiltonian
    !%Description
    !% The value 'Alpha' that controls the splitting of the Coulomb
    !% interaction in the Ewald sum used to calculation the ion-ion
    !% interaction for periodic systems. This value affects the speed
    !% of the calculation, normally users do not need to modify it.
    !%End
    call parse_variable('EwaldAlpha', CNST(0.21), this%alpha)
    
    nullify(this%subsys_ionic)

    POP_SUB(ion_interaction_init)
  end subroutine ion_interaction_init
  
  ! ---------------------------------------------------------
  
  subroutine ion_interaction_end(this)
    type(ion_interaction_t), intent(inout) :: this

    PUSH_SUB(ion_interaction_end)

    this%alpha = -CNST(1.0)
    
    nullify(this%subsys_ionic)

    POP_SUB(ion_interaction_end)
  end subroutine ion_interaction_end
  
  ! ---------------------------------------------------------
  
  subroutine ion_interaction_add_subsys_ionic(this, subsys_ionic)
    type(ion_interaction_t),    intent(inout) :: this
    type(ssys_ionic_t), target, intent(in)    :: subsys_ionic
    
    PUSH_SUB(ion_interaction_add_subsys_ionic)
    
    this%subsys_ionic => subsys_ionic
    
    POP_SUB(ion_interaction_add_subsys_ionic)
  end subroutine ion_interaction_add_subsys_ionic
  
  ! ---------------------------------------------------------
  !> For details about this routine, see
  !! http://www.tddft.org/programs/octopus/wiki/index.php/Developers:Ion-Ion_interaction
  subroutine ion_interaction_calculate(this, geo, sb, ignore_external_ions, energy, force, energy_components, force_components)
    type(ion_interaction_t),  intent(inout) :: this
    type(geometry_t), target, intent(in)    :: geo
    type(simul_box_t),        intent(in)    :: sb
    logical,                  intent(in)    :: ignore_external_ions
    FLOAT,                    intent(out)   :: energy
    FLOAT,                    intent(out)   :: force(:, :)
    FLOAT, optional,          intent(out)   :: energy_components(:)
    FLOAT, optional,          intent(out)   :: force_components(:, :, :)
    
    FLOAT, allocatable:: r(:), f(:)
    FLOAT :: rr, dd, zi, zj, epsilon, sigma
    integer :: iatom, jatom, natom, iindex, jindex
    type(species_t), pointer :: spci, spcj
    type(profile_t), save :: ion_ion_prof
    logical,  allocatable :: in_box(:)

    PUSH_SUB(ion_interaction_calculate)
    call profiling_in(ion_ion_prof, "ION_ION_INTERACTION")

    SAFE_ALLOCATE(r(1:sb%dim))
    SAFE_ALLOCATE(f(1:sb%dim))

    if(present(energy_components)) then
      ASSERT(ubound(energy_components, dim = 1) == ION_NUM_COMPONENTS)
      energy_components = CNST(0.0)
    end if      

    if(present(force_components)) then
      ASSERT(all(ubound(force_components) == (/sb%dim, geo%natoms, ION_NUM_COMPONENTS/)))
      force_components = CNST(0.0)
    end if 

    energy = M_ZERO
    if(associated(this%subsys_ionic))then
      ! Get the subsystems interaction energy.
      call ssys_ionic_calc(this%subsys_ionic)
      call ssys_ionic_get(this%subsys_ionic, energy, except=(/"live"/))
    end if
    force(1:sb%dim, 1:geo%natoms) = M_ZERO

    if(simul_box_is_periodic(sb)) then

      ASSERT(geo%ncatoms==0)

      spci => geo%atom(1)%species
      ! This depends on the area, but we should check if it is fully consistent.        
      if( species_type(spci) == SPECIES_JELLIUM_SLAB ) then
        energy = energy + &
          M_PI*species_zval(spci)**2/(M_FOUR*sb%lsize(1)*sb%lsize(2))*(sb%lsize(3) - species_jthick(spci)/M_THREE)
      else
        call ion_interaction_periodic(this, geo, sb, energy, force, energy_components, force_components)
      end if

    else
      
      natom = geo%natoms + geo%ncatoms

      if(ignore_external_ions) then
        SAFE_ALLOCATE(in_box(1:natom))
        do iatom = 1, geo%natoms
          in_box(iatom) = simul_box_in_box(sb, geo, geo%atom(iatom)%x)
        end do
        do iatom = 1, geo%ncatoms
          in_box(geo%natoms + iatom) = simul_box_in_box(sb, geo, geo%catom(iatom)%x)
        end do
      end if
      
      ! only interaction inside the cell
      do iatom = 1, geo%natoms
        if(ignore_external_ions) then
          if(.not. in_box(iatom)) cycle
        end if
        
        if(associated(this%subsys_ionic))then
          ! Calculate the interaction forces.
          call ssys_ionic_interaction(this%subsys_ionic, geo%atom(iatom), f, except=(/"live"/) )
          force(1:sb%dim,iatom) = force(1:sb%dim,iatom) + f(1:sb%dim)
        end if

        spci => geo%atom(iatom)%species
        zi = species_zval(spci)

        select case(species_type(spci))
        case(SPECIES_JELLIUM)
          energy = energy + (M_THREE/M_FIVE)*zi**2/species_jradius(spci)
          ! The part depending on the simulation sphere is neglected
          
        case(SPECIES_JELLIUM_SLAB)
          energy = energy - M_PI*zi**2/(M_FOUR*sb%lsize(1)*sb%lsize(2))*species_jthick(spci)/M_THREE
          ! The part depending on the simulation box transverse dimension is neglected
        end select
        
        do jatom = iatom + 1, geo%natoms

          if(ignore_external_ions) then
            if(.not. in_box(jatom)) cycle
          end if
          
          spcj => geo%atom(jatom)%species

          r(1:sb%dim) = geo%atom(iatom)%x(1:sb%dim) - geo%atom(jatom)%x(1:sb%dim)

          rr = sqrt(sum(r**2))
          iindex = species_index(spci)
          jindex = species_index(spcj)
          
          select case(geo%ionic_interaction_type(iindex, jindex))
          case(INTERACTION_COULOMB)
            zj = species_zval(spcj)
            !the force
            dd = zi*zj/rr
            f(1:sb%dim) = (dd/rr**2)*r(1:sb%dim)
            force(1:sb%dim,iatom) = force(1:sb%dim,iatom) + f(1:sb%dim)
            force(1:sb%dim,jatom) = force(1:sb%dim,jatom) - f(1:sb%dim)
            !energy
            energy=energy + dd

          case(INTERACTION_LJ)
            epsilon= geo%ionic_interaction_parameter(LJ_EPSILON, iindex, jindex)
            sigma  = geo%ionic_interaction_parameter(LJ_SIGMA,   iindex, jindex)
            dd = (sigma/rr)**6

            !the force
            f(1:sb%dim) = (CNST(24.0)*epsilon*(dd/rr**2)*(CNST(2.0)*dd - M_ONE))*r
            force(1:sb%dim, iatom) = force(1:sb%dim,iatom) + f(1:sb%dim)
            force(1:sb%dim, jatom) = force(1:sb%dim,jatom) - f(1:sb%dim)

            !energy
            energy = energy + CNST(4.0)*epsilon*dd*(dd - M_ONE)
          end select
          
        end do !jatom
      end do !iatom
      
      do iatom = 1, geo%natoms

        if(ignore_external_ions) then
          if(.not. in_box(iatom)) cycle
        end if
        
        do jatom = 1, geo%ncatoms
          if(ignore_external_ions) then
            if(.not. in_box(geo%natoms+jatom)) cycle
          end if
          
          r(1:sb%dim) = geo%atom(iatom)%x(1:sb%dim) - geo%catom(jatom)%x(1:sb%dim)
          rr = sqrt(sum(r**2))
          !INTERACTION_COULOMB
          zi = species_zval(geo%atom(iatom)%species)
          zj = geo%catom(jatom)%charge
          !the force
          dd = zi*zj/rr
          force(1:sb%dim,iatom) = force(1:sb%dim,iatom) + (dd/rr*2)*r(1:sb%dim)
          !energy
          energy = energy + dd
          
        end do !jatom
      end do !iatom
      
    end if

    SAFE_DEALLOCATE_A(in_box)
    SAFE_DEALLOCATE_A(r)
    SAFE_DEALLOCATE_A(f)
    
    call profiling_out(ion_ion_prof)
    
    POP_SUB(ion_interaction_calculate)
  end subroutine ion_interaction_calculate

  ! ---------------------------------------------------------
  
  subroutine ion_interaction_periodic(this, geo, sb, energy, force, energy_components, force_components)
    type(ion_interaction_t),   intent(in)    :: this
    type(geometry_t),          intent(in)    :: geo
    type(simul_box_t),         intent(in)    :: sb
    FLOAT,                     intent(out)   :: energy
    FLOAT,                     intent(out)   :: force(:, :) !< sb%dim, geo%natoms
    FLOAT, optional,           intent(out)   :: energy_components(:)
    FLOAT, optional,           intent(out)   :: force_components(:, :, :)

    FLOAT :: rr, xi(1:MAX_DIM), zi, zj, ereal, efourier, eself, erfc, rcut
    integer :: iatom, jatom, icopy, atom_start, atom_end
    type(periodic_copy_t) :: pc
    integer :: ix, iy, iz, isph, ss, idim
    FLOAT   :: gg(1:MAX_DIM), gg2, gx
    FLOAT   :: factor, charge
    CMPLX   :: sumatoms, tmp(1:MAX_DIM), aa
    CMPLX, allocatable :: phase(:)
    type(profile_t), save :: prof_short, prof_long

    PUSH_SUB(ion_interaction_periodic)

    ! Check http://www.duke.edu/~kz10/file/ewald.pdf for the equations
    ! implemented here.

    if(any(geo%ionic_interaction_type /= INTERACTION_COULOMB)) then
      message(1) = "Cannot calculate non-Coulombic interaction for periodic systems."
      call messages_fatal(1)
    end if

    ereal = M_ZERO

    force(1:sb%dim, 1:geo%natoms) = M_ZERO

    ! if the system is periodic we have to add the energy of the
    ! interaction with the copies
    
    rcut = CNST(6.0)/this%alpha

    call profiling_in(prof_short, "EWALD_SHORT")
    
    ! the short-range part is calculated directly
    do iatom = geo%atoms_dist%start, geo%atoms_dist%end
      if (.not. species_represents_real_atom(geo%atom(iatom)%species)) cycle
      zi = species_zval(geo%atom(iatom)%species)

      call periodic_copy_init(pc, sb, geo%atom(iatom)%x, rcut)
      
      do icopy = 1, periodic_copy_num(pc)
        xi(1:sb%dim) = periodic_copy_position(pc, sb, icopy)
        
        do jatom = 1,  geo%natoms
          zj = species_zval(geo%atom(jatom)%species)
          rr = sqrt( sum( (xi(1:sb%dim) - geo%atom(jatom)%x(1:sb%dim))**2 ) )
          
          if(rr < CNST(1e-5)) cycle
          
          erfc = M_ONE - loct_erf(this%alpha*rr)

          ! energy
          ereal = ereal + M_HALF*zj*zi*erfc/rr
          
          ! force
          force(1:sb%dim, jatom) = force(1:sb%dim, jatom) - &
            zj*zi*(xi(1:sb%dim) - geo%atom(jatom)%x(1:sb%dim))*&
            (erfc/rr + M_TWO*this%alpha/sqrt(M_PI)*exp(-(this%alpha*rr)**2))/rr**2
        end do
        
      end do
      
      call periodic_copy_end(pc)
    end do

    if(geo%atoms_dist%parallel) then
      call comm_allreduce(geo%atoms_dist%mpi_grp%comm, ereal)
      call comm_allreduce(geo%atoms_dist%mpi_grp%comm, force)
    end if

    if(present(force_components)) then
      force_components(1:sb%dim, 1:geo%natoms, ION_COMPONENT_REAL) = force(1:sb%dim, 1:geo%natoms)
    end if

    call profiling_out(prof_short)

    call profiling_in(prof_long, "EWALD_LONG")
    
    ! self-interaction
    eself = M_ZERO
    charge = M_ZERO
    do iatom = 1, geo%natoms
      zi = species_zval(geo%atom(iatom)%species)
      charge = charge + zi
      eself = eself - this%alpha/sqrt(M_PI)*zi**2
    end do

    ! And the long-range part, using an Ewald sum
    SAFE_ALLOCATE(phase(1:geo%natoms))

    ! get a converged value for the cutoff in g
    rcut = huge(rcut)
    do idim = 1, sb%dim
      rcut = min(rcut, sum(sb%klattice(1:sb%dim, idim)**2))
    end do

    rcut = sqrt(rcut)
    
    isph = ceiling(CNST(9.5)*this%alpha/rcut)

    ! First the G = 0 term (charge was calculated previously)
    efourier = -M_PI*charge**2/(M_TWO*this%alpha**2*sb%rcell_volume)

    do ix = -isph, isph
      do iy = -isph, isph
        do iz = -isph, isph
          
          ss = ix**2 + iy**2 + iz**2
          
          if(ss == 0 .or. ss > isph**2) cycle

          gg(1:sb%dim) = ix*sb%klattice(1:sb%dim, 1) + iy*sb%klattice(1:sb%dim, 2) + iz*sb%klattice(1:sb%dim, 3)
          gg2 = sum(gg(1:sb%dim)**2)

          ! g=0 must be removed from the sum
          if(gg2 < M_EPSILON) cycle
          
          gx = -CNST(0.25)*gg2/this%alpha**2

          if(gx < CNST(-36.0)) cycle

          factor = M_TWO*M_PI/sb%rcell_volume*exp(gx)/gg2

          if(factor < epsilon(factor)) cycle

          sumatoms = M_Z0
          !$omp parallel do private(iatom, gx, aa) reduction(+:sumatoms)
          do iatom = 1, geo%natoms
            gx = sum(gg(1:sb%dim)*geo%atom(iatom)%x(1:sb%dim))
            aa = species_zval(geo%atom(iatom)%species)*TOCMPLX(cos(gx), sin(gx))
            phase(iatom) = aa 
            sumatoms = sumatoms + aa
          end do
          
          efourier = efourier + TOFLOAT(factor*sumatoms*conjg(sumatoms))
          
          do iatom = 1, geo%natoms
            tmp(1:sb%dim) = M_ZI*gg(1:sb%dim)*phase(iatom)
            force(1:sb%dim, iatom) = force(1:sb%dim, iatom) &
              - factor*TOFLOAT(conjg(tmp(1:sb%dim))*sumatoms + tmp(1:sb%dim)*conjg(sumatoms))
          end do
          
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(phase)

    if(present(energy_components)) then
      energy_components(ION_COMPONENT_REAL) = ereal
      energy_components(ION_COMPONENT_SELF) = eself
      energy_components(ION_COMPONENT_FOURIER) = efourier      
    end if

    if(present(force_components)) then
      force_components(1:sb%dim, 1:geo%natoms, ION_COMPONENT_FOURIER) = &
        force(1:sb%dim, 1:geo%natoms) - force_components(1:sb%dim, 1:geo%natoms, ION_COMPONENT_REAL)
    end if

    energy = ereal + efourier + eself
    
    call profiling_out(prof_long)
    
    POP_SUB(ion_interaction_periodic)
  end subroutine ion_interaction_periodic

  ! --------------------------------------------------------------
  
  subroutine ion_interaction_test(geo, sb)
    type(geometry_t),         intent(in)    :: geo
    type(simul_box_t),        intent(in)    :: sb

    type(ion_interaction_t) :: ion_interaction
    FLOAT :: energy
    FLOAT, allocatable :: force(:, :), force_components(:, :, :)
    FLOAT :: energy_components(1:ION_NUM_COMPONENTS)
    integer :: iatom, idir
    
    PUSH_SUB(ion_interaction_test)

    call ion_interaction_init(ion_interaction)

    SAFE_ALLOCATE(force(1:sb%dim, 1:geo%natoms))
    SAFE_ALLOCATE(force_components(1:sb%dim, 1:geo%natoms, ION_NUM_COMPONENTS))
    
    call ion_interaction_calculate(ion_interaction, geo, sb, .false., energy, force, &
      energy_components = energy_components, force_components = force_components)

    call messages_write('Ionic energy        =')
    call messages_write(energy, fmt = '(f20.10)')
    call messages_info()

    call messages_write('Real space energy   =')
    call messages_write(energy_components(ION_COMPONENT_REAL), fmt = '(f20.10)')
    call messages_info()

    call messages_write('Self energy         =')
    call messages_write(energy_components(ION_COMPONENT_SELF), fmt = '(f20.10)')
    call messages_info()

    call messages_write('Fourier energy      =')
    call messages_write(energy_components(ION_COMPONENT_FOURIER), fmt = '(f20.10)')
    call messages_info()
    
    call messages_info()

    do iatom = 1, geo%natoms
      call messages_write('Ionic force         atom')
      call messages_write(iatom)
      call messages_write(' =')
      do idir = 1, sb%dim
        call messages_write(force(idir, iatom), fmt = '(f20.10)')
      end do
      call messages_info()

      call messages_write('Real space force    atom')
      call messages_write(iatom)
      call messages_write(' =')
      do idir = 1, sb%dim
        call messages_write(force_components(idir, iatom, ION_COMPONENT_REAL), fmt = '(f20.10)')
      end do
      call messages_info()

      call messages_write('Fourier space force atom')
      call messages_write(iatom)
      call messages_write(' =')
      do idir = 1, sb%dim
        call messages_write(force_components(idir, iatom, ION_COMPONENT_FOURIER), fmt = '(f20.10)')
      end do
      call messages_info()

      call messages_info()
    end do

    SAFE_DEALLOCATE_A(force)
    SAFE_DEALLOCATE_A(force_components)
    
    call ion_interaction_end(ion_interaction)
    
    POP_SUB(ion_interaction_test)
  end subroutine ion_interaction_test
    
end module ion_interaction_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End: