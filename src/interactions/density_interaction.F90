!! Copyright (C) 2021 N. Tancogne-Dejean
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

module density_interaction_oct_m
  use interaction_with_partner_oct_m

  implicit none

  private
  public ::                &
    density_interaction_t

  type, extends(interaction_with_partner_t), abstract :: density_interaction_t
    FLOAT, allocatable, public :: density(:,:)
  end type density_interaction_t

end module density_interaction_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
