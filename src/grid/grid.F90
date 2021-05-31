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

#include "global.h"

module grid_oct_m
  use affine_coordinates_oct_m
  use box_image_oct_m
  use cartesian_oct_m
  use coordinate_system_oct_m
  use cube_oct_m
  use curv_briggs_oct_m
  use curv_gygi_oct_m
  use curv_modine_oct_m
  use derivatives_oct_m
  use global_oct_m
  use ions_oct_m
  use mesh_oct_m
  use mesh_init_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use namespace_oct_m
  use nl_operator_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use simul_box_oct_m
  use stencil_oct_m
  use stencil_cube_oct_m
  use symmetries_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                &
    grid_t,                &
    grid_init_stage_1,     &
    grid_init_stage_2,     &
    grid_end,              &
    grid_write_info

  type grid_t
    ! Components are public by default
    type(simul_box_t)           :: sb
    type(mesh_t)                :: mesh
    type(derivatives_t)         :: der
    class(coordinate_system_t), pointer :: coord_system
    type(stencil_t)             :: stencil

    type(symmetries_t)          :: symm
  end type grid_t

  integer, parameter :: &
    CURV_AFFINE  = 1,   &
    CURV_GYGI    = 2,   &
    CURV_BRIGGS  = 3,   &
    CURV_MODINE  = 4

contains

  !-------------------------------------------------------------------
  subroutine grid_init_stage_1(gr, namespace, ions, space)
    type(grid_t),      intent(inout) :: gr
    type(namespace_t), intent(in)    :: namespace
    type(ions_t),      intent(inout) :: ions
    type(space_t),     intent(in)    :: space

    type(stencil_t) :: cube
    integer :: enlarge(1:MAX_DIM)
    type(block_t) :: blk
    integer :: idir, cv_method
    FLOAT :: def_h, def_rsize
    FLOAT :: grid_spacing(1:MAX_DIM)

    PUSH_SUB(grid_init_stage_1)

    call simul_box_init(gr%sb, namespace, ions, space)

    call symmetries_init(gr%symm, namespace, ions, space)

    call ions%grid_defaults(def_h, def_rsize)
    
    ! initialize to -1
    grid_spacing = -M_ONE

    !%Variable Spacing
    !%Type float
    !%Section Mesh
    !%Description
    !% The spacing between the points in the mesh. This controls the
    !% quality of the discretization: smaller spacing gives more
    !% precise results but increased computational cost.
    !%
    !% When using curvilinear coordinates, this is a canonical spacing
    !% that will be changed locally by the transformation. In periodic
    !% directions, your spacing may be slightly different than what
    !% you request here, since the box size must be an integer
    !% multiple of the spacing.
    !%
    !% The default value is defined by the species if only default pseudopotentials are used
    !% or by the image resolution if <tt>BoxShape = box_image</tt>. Otherwise, there is
    !% no default.
    !%
    !% It is possible to have a different spacing in each one of the Cartesian directions
    !% if we define <tt>Spacing</tt> as block of the form
    !%
    !% <tt>%Spacing
    !% <br>&nbsp;&nbsp;spacing_x | spacing_y | spacing_z
    !% <br>%</tt>
    !%End

    if(parse_block(namespace, 'Spacing', blk) == 0) then
      if(parse_block_cols(blk,0) < space%dim) call messages_input_error(namespace, 'Spacing')
      do idir = 1, space%dim
        call parse_block_float(blk, 0, idir - 1, grid_spacing(idir), units_inp%length)
        if(def_h > M_ZERO) call messages_check_def(grid_spacing(idir), .true., def_h, 'Spacing', units_out%length)
      end do
      call parse_block_end(blk)
    else
      call parse_variable(namespace, 'Spacing', -M_ONE, grid_spacing(1), units_inp%length)
      grid_spacing(1:space%dim) = grid_spacing(1)
      if(def_h > M_ZERO) call messages_check_def(grid_spacing(1), .true., def_h, 'Spacing', units_out%length)
    end if

    if (associated(gr%sb%box)) then
      select type (box => gr%sb%box)
      type is (box_image_t)
        do idir = 1, space%dim
          ! default grid_spacing is determined from the pixel size such that one grid point = one pixel.
          if(grid_spacing(idir) < M_ZERO) then
            grid_spacing(idir) = box%pixel_size(idir)
          end if
        end do
      end select
    end if

    if (any(grid_spacing(1:space%dim) < M_EPSILON)) then
      if (def_h > M_ZERO .and. def_h < huge(def_h)) then
        call ions%grid_defaults_info()
        do idir = 1, space%dim
          grid_spacing(idir) = def_h
          write(message(1), '(a,i1,3a,f6.3)') "Info: Using default spacing(", idir, &
            ") [", trim(units_abbrev(units_out%length)), "] = ",                        &
            units_from_atomic(units_out%length, grid_spacing(idir))
          call messages_info(1)
        end do
        ! Note: the default automatically matches the 'recommended' value compared by messages_check_def above.
      else
        message(1) = 'Either:'
        message(2) = "   *) variable 'Spacing' is not defined and"
        message(3) = "      I can't find a suitable default"
        message(4) = "   *) your input for 'Spacing' is negative or zero"
        call messages_fatal(4)
      end if
    end if

    !%Variable PeriodicBoundaryMask
    !%Type block
    !%Section Mesh
    !%Description
    !% (Experimental) Defines a mask for which periodic boundaries are replaced by zero boundary conditions.
    !%End
    if(parse_block(namespace, 'PeriodicBoundaryMask', blk) < 0) then
      gr%mesh%masked_periodic_boundaries = .false.
    else
      gr%mesh%masked_periodic_boundaries = .true.
      call parse_block_string(blk, 0, 0, gr%mesh%periodic_boundary_mask)
      call messages_experimental('PeriodicBoundaryMask')
    end if

    ! Initialize coordinate system

    !%Variable CurvMethod
    !%Type integer
    !%Default curv_uniform
    !%Section Mesh::Curvilinear
    !%Description
    !% The relevant functions in octopus are represented on a mesh in real space.
    !% This mesh may be an evenly spaced regular rectangular grid (standard mode),
    !% or else an adaptive or curvilinear grid. We have implemented
    !% three kinds of adaptive meshes, although only one is currently working,
    !% the one invented by F. Gygi (<tt>curv_gygi</tt>). The code will stop if any of
    !% the other two is invoked. All are experimental with domain parallelization.
    !%Option curv_affine 1
    !% Regular, uniform rectangular grid.
    !%Option curv_gygi 2
    !% The deformation of the grid is done according to the scheme described by
    !% F. Gygi [F. Gygi and G. Galli, <i>Phys. Rev. B</i> <b>52</b>, R2229 (1995)].
    !%Option curv_briggs 3
    !% The deformation of the grid is done according to the scheme described by
    !% Briggs [E.L. Briggs, D.J. Sullivan, and J. Bernholc, <i>Phys. Rev. B</i> <b>54</b> 14362 (1996)]
    !% (NOT WORKING).
    !%Option curv_modine 4
    !% The deformation of the grid is done according to the scheme described by
    !% Modine [N.A. Modine, G. Zumbach and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289 (1997)]
    !% (NOT WORKING).
    !%End
    call parse_variable(namespace, 'CurvMethod', CURV_AFFINE, cv_method)
    if (.not. varinfo_valid_option('CurvMethod', cv_method)) call messages_input_error(namespace, 'CurvMethod')
    call messages_print_var_option(stdout, "CurvMethod", cv_method)

    ! FIXME: The other two methods are apparently not working
    if (cv_method > CURV_GYGI) then
      call messages_experimental('Selected curvilinear coordinates method')
    end if

    select case (cv_method)
    case (CURV_BRIGGS)
      gr%coord_system => curv_briggs_t(namespace, space%dim, gr%sb%lsize(1:space%dim), grid_spacing(1:space%dim))
    case (CURV_GYGI)
      gr%coord_system => curv_gygi_t(namespace, space%dim, ions%natoms, ions%pos)
    case (CURV_MODINE)
      gr%coord_system => curv_modine_t(namespace, space%dim, ions%natoms, ions%pos, gr%sb%lsize(1:space%dim), &
        grid_spacing(1:space%dim))
    case (CURV_AFFINE)
      if (ions%latt%nonorthogonal) then
        gr%coord_system => affine_coordinates_t(namespace, space%dim, ions%latt%rlattice_primitive)
      else
        gr%coord_system => cartesian_t(namespace, space%dim)
      end if
    end select

    ! initialize derivatives
    call derivatives_init(gr%der, namespace, space, gr%sb%latt, gr%coord_system%local_basis)
    ! the stencil used to generate the grid is a union of a cube (for
    ! multigrid) and the Laplacian.
    call stencil_cube_get_lapl(cube, space%dim, order = 2)
    call stencil_union(space%dim, cube, gr%der%lapl%stencil, gr%stencil)
    call stencil_end(cube)

    enlarge = 0
    enlarge(1:space%dim) = 2
    enlarge = max(enlarge, gr%der%n_ghost)

    call mesh_init_stage_1(gr%mesh, namespace, space, gr%sb, gr%coord_system, grid_spacing, enlarge)
    call mesh_init_stage_2(gr%mesh, space, gr%sb, gr%stencil)

    POP_SUB(grid_init_stage_1)
  end subroutine grid_init_stage_1


  !-------------------------------------------------------------------
  subroutine grid_init_stage_2(gr, namespace, space, mc)
    type(grid_t), target, intent(inout) :: gr
    type(namespace_t),    intent(in)    :: namespace
    type(space_t),        intent(in)    :: space
    type(multicomm_t),    intent(in)    :: mc

    PUSH_SUB(grid_init_stage_2)

    call mesh_init_stage_3(gr%mesh, namespace, space, gr%stencil, mc)

    call nl_operator_global_init(namespace)
    call derivatives_build(gr%der, namespace, space, gr%mesh)

    ! print info concerning the grid
    call grid_write_info(gr, stdout)

    POP_SUB(grid_init_stage_2)
  end subroutine grid_init_stage_2


  !-------------------------------------------------------------------
  subroutine grid_end(gr)
    type(grid_t), intent(inout) :: gr

    class(coordinate_system_t), pointer :: coord_system

    PUSH_SUB(grid_end)

    call nl_operator_global_end()

    call derivatives_end(gr%der)

    ! We need to take a pointer here, otherwise we run into a gfortran bug.
    coord_system => gr%coord_system
    SAFE_DEALLOCATE_P(coord_system)

    call mesh_end(gr%mesh)

    call symmetries_end(gr%symm)

    call stencil_end(gr%stencil)

    POP_SUB(grid_end)
  end subroutine grid_end


  !-------------------------------------------------------------------
  subroutine grid_write_info(gr, iunit)
    type(grid_t),     intent(in) :: gr
    integer,          intent(in) :: iunit

    PUSH_SUB(grid_write_info)

    if(.not.mpi_grp_is_root(mpi_world)) then
      if(debug%info) call messages_debug_newlines(6)
      POP_SUB(grid_write_info)
      return
    end if

    call messages_print_stress(iunit, "Grid")

    message(1) = 'Simulation Box:'
    call messages_info(1, iunit)
    call gr%sb%write_info(iunit)

    message(1) = "Main mesh:"
    call messages_info(1, iunit)
    call mesh_write_info(gr%mesh, iunit)

    if (gr%mesh%use_curvilinear) then
      call gr%coord_system%write_info(iunit)
    end if
    
    call messages_print_stress(iunit)

    POP_SUB(grid_write_info)
  end subroutine grid_write_info

end module grid_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
