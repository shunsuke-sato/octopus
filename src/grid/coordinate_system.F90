!! Copyright (C) 2021 M. Oliveira
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

module coordinate_system_oct_m
  implicit none

  private
  public :: coordinate_system_t

  type, abstract :: coordinate_system_t
    logical :: local_basis  !< Do the basis vectors depend on the position, i.e., is the basis local? (false for Cartesian and affine, true for curvilinear coordinates in general)
    logical :: orthogonal   !< Are the basis vectors orthogonal?
    integer :: dim          !< Dimension of the space
    FLOAT :: min_mesh_scaling_product !< product of the smallest scaling :: min(distance between the grid points / spacing)
  contains
    procedure(coordinate_system_chi2x),      deferred :: chi2x
    procedure(coordinate_system_x2chi),      deferred :: x2chi
    procedure(coordinate_system_det_jac),    deferred :: det_jac
    procedure(coordinate_system_write_info), deferred :: write_info
  end type coordinate_system_t

  abstract interface
    ! ---------------------------------------------------------
    subroutine coordinate_system_chi2x(this, chi, xx)
      import coordinate_system_t
      class(coordinate_system_t), target, intent(in)  :: this
      FLOAT,                              intent(in)  :: chi(:)
      FLOAT,                              intent(out) :: xx(:)
    end subroutine coordinate_system_chi2x

    ! ---------------------------------------------------------
    subroutine coordinate_system_x2chi(this, xx, chi)
      import coordinate_system_t
      class(coordinate_system_t), target, intent(in)  :: this
      FLOAT,                              intent(in)  :: xx(:)
      FLOAT,                              intent(out) :: chi(:)
    end subroutine coordinate_system_x2chi

    ! ---------------------------------------------------------
    FLOAT function coordinate_system_det_jac(this, xx, chi) result(jdet)
      import coordinate_system_t
      class(coordinate_system_t), intent(in)  :: this
      FLOAT,                      intent(in)  :: xx(:)
      FLOAT,                      intent(in)  :: chi(:)
    end function coordinate_system_det_jac

    ! ---------------------------------------------------------
    subroutine coordinate_system_write_info(this, unit)
      import coordinate_system_t
      class(coordinate_system_t), intent(in) :: this
      integer,                    intent(in) :: unit
    end subroutine coordinate_system_write_info
  end interface

end module coordinate_system_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
