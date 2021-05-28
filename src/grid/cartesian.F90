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

module cartesian_oct_m
  use affine_coordinates_oct_m
  use basis_vectors_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::         &
    cartesian_t,    &
    cartesian_copy

  type, extends(affine_coordinates_t) :: cartesian_t
    private
  contains
    procedure :: chi2x => cartesian_chi2x
    procedure :: x2chi => cartesian_x2chi
    procedure :: det_jac => cartesian_det_jac
    procedure :: write_info => cartesian_write_info
  end type cartesian_t

  interface cartesian_t
    procedure cartesian_constructor
  end interface cartesian_t

contains

  ! ---------------------------------------------------------
  function cartesian_constructor(namespace, dim) result(cart)
    type(namespace_t), intent(in)  :: namespace
    integer,           intent(in)  :: dim
    class(cartesian_t), pointer :: cart

    integer :: idir
    FLOAT :: basis_vectors(dim, dim)

    PUSH_SUB(cartesian_constructor)

    SAFE_ALLOCATE(cart)

    cart%dim = dim
    cart%local_basis = .false.
    cart%orthogonal = .true. 
    basis_vectors = M_ZERO
    do idir = 1, dim
      basis_vectors(idir, idir) = M_ONE
    end do
    cart%basis = basis_vectors_t(namespace, dim, basis_vectors)

    POP_SUB(cartesian_constructor)
  end function cartesian_constructor

  ! -------------------------------------------------------------- 
  subroutine cartesian_copy(this_out, this_in)
    type(cartesian_t), intent(inout) :: this_out
    type(cartesian_t), intent(in)    :: this_in

    PUSH_SUB(cartesian_copy)

    this_out%dim = this_in%dim
    this_out%local_basis = this_in%local_basis

    POP_SUB(cartesian_copy)
  end subroutine cartesian_copy

  ! ---------------------------------------------------------
  subroutine cartesian_chi2x(this, chi, xx)
    class(cartesian_t), target, intent(in)  :: this
    FLOAT,                      intent(in)  :: chi(:)
    FLOAT,                      intent(out) :: xx(:)

    ! no PUSH_SUB, called too often

    xx = chi

  end subroutine cartesian_chi2x

  ! ---------------------------------------------------------
  subroutine cartesian_x2chi(this, xx, chi)
    class(cartesian_t), target, intent(in)  :: this
    FLOAT,                      intent(in)  :: xx(:)
    FLOAT,                      intent(out) :: chi(:)

    ! no PUSH_SUB, called too often

    chi = xx

  end subroutine cartesian_x2chi

  ! ---------------------------------------------------------
  FLOAT function cartesian_det_jac(this, xx, chi) result(jdet)
    class(cartesian_t),    intent(in)  :: this
    FLOAT,                 intent(in)  :: xx(:)
    FLOAT,                 intent(in)  :: chi(:)

    ! No PUSH_SUB, called too often

    jdet = M_ONE

  end function cartesian_det_jac

  ! ---------------------------------------------------------
  subroutine cartesian_write_info(this, unit)
    class(cartesian_t), intent(in) :: this
    integer,            intent(in) :: unit

    PUSH_SUB(cartesian_write_info)

    write(message(1), '(a)')  '  Using Cartesian coordinates'
    call messages_info(1, unit)

    POP_SUB(cartesian_write_info)
  end subroutine cartesian_write_info

end module cartesian_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
