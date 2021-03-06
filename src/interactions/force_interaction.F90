!! Copyright (C) 2020 M. Oliveira
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

module force_interaction_oct_m
  use interaction_with_partner_oct_m

  implicit none

  private
  public ::                &
    force_interaction_t

  type, extends(interaction_with_partner_t), abstract :: force_interaction_t
    integer :: dim = 0       !< spatial dimensions
    integer :: system_np = 0 !< number of particles in the system that the forces are acting on

    FLOAT, allocatable, public :: force(:,:)
  end type force_interaction_t

end module force_interaction_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
