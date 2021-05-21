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

!> This module implements the curvilinear coordinates given in
!! E.L. Briggs, D.J. Sullivan, and J. Bernholc, PRB 54 14362 (1996)
!!
!! It assumes that the Oxygen atom is located at x0=0 (see Eq. (12))

module curv_briggs_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                     &
    curv_briggs_t,              &
    curv_briggs_init,           &
    curv_briggs_end,            &
    curv_briggs_copy,           &
    curv_briggs_chi2x,          &
    curv_briggs_jacobian_inv

  type curv_briggs_t
    private
    FLOAT, allocatable :: lsize(:) !< size of the box
    FLOAT :: beta                  !< adjustable parameter between 0 and 1 that controls the degree of scaling
  end type curv_briggs_t

contains

  ! ---------------------------------------------------------
  subroutine curv_briggs_init(cv, namespace, dim, lsize, spacing, min_scaling_product)
    type(curv_briggs_t), intent(out) :: cv
    type(namespace_t),   intent(in)  :: namespace
    integer,             intent(in)  :: dim
    FLOAT,               intent(in)  :: lsize(1:dim)
    FLOAT,               intent(in)  :: spacing(1:dim)
    FLOAT,               intent(out) :: min_scaling_product
  
    SAFE_ALLOCATE(cv%lsize(1:dim))
    cv%lsize(1:dim) = lsize(1:dim)

    call parse_variable(namespace, 'CurvBriggsBeta', M_HALF, cv%beta)

    if (cv%beta < M_ZERO .or. cv%beta > M_ONE) then
      message(1) = 'The parameter "CurvBriggsBeta" must lie between 0 and 1.'
      call messages_fatal(1)
    end if

    call curv_briggs_min_scaling(cv, dim, spacing, min_scaling_product)

  end subroutine curv_briggs_init

  ! ---------------------------------------------------------
  subroutine curv_briggs_end(cv)
    type(curv_briggs_t), intent(inout) :: cv

    PUSH_SUB(curv_briggs_end)

    SAFE_DEALLOCATE_A(cv%lsize)

    POP_SUB(curv_briggs_end)
  end subroutine curv_briggs_end

  ! ---------------------------------------------------------
  subroutine curv_briggs_copy(this_out, this_in)
    type(curv_briggs_t), intent(inout) :: this_out
    type(curv_briggs_t), intent(in)    :: this_in

    PUSH_SUB(curv_briggs_copy)

    SAFE_ALLOCATE_SOURCE_A(this_out%lsize, this_in%lsize)
    this_out%beta = this_in%beta

    POP_SUB(curv_briggs_copy)
  end subroutine curv_briggs_copy

  ! ---------------------------------------------------------
  subroutine curv_briggs_chi2x(cv, dim, chi, x)
    type(curv_briggs_t), intent(in)  :: cv
    integer,             intent(in)  :: dim
    FLOAT,               intent(in)  :: chi(:)  !< chi(dim)
    FLOAT,               intent(out) :: x(:)    !< x(dim)

    x = chi - cv%lsize*cv%beta/(M_TWO*M_PI)*sin(M_TWO*M_PI*chi/cv%lsize)

  end subroutine curv_briggs_chi2x

  ! ---------------------------------------------------------
  subroutine curv_briggs_jacobian_inv(cv, dim, chi, jac)
    type(curv_briggs_t), intent(in)  :: cv
    integer,             intent(in)  :: dim
    FLOAT,               intent(in)  :: chi(:)    !< x(dim)
    FLOAT,               intent(out) :: jac(:,:)  !< jac(dim,dim), the Jacobian

    integer :: idim

    jac = M_ZERO
    do idim = 1, dim
      jac(idim, idim) = M_ONE - cv%beta*cos(M_TWO*M_PI*chi(idim)/cv%lsize(idim))
    end do

  end subroutine curv_briggs_jacobian_inv

  ! ---------------------------------------------------------
  subroutine curv_briggs_min_scaling(cv, dim, spacing, min_scaling_product)
    type(curv_briggs_t), intent(in)  :: cv
    integer,             intent(in)  :: dim
    FLOAT,               intent(in)  :: spacing(1:dim)
    FLOAT,               intent(out) :: min_scaling_product

    integer :: idim

    min_scaling_product = M_ONE
    do idim = 1, dim
      ! corresponds to the distance of grid points at [+spacing/2,-spacing/2]
      min_scaling_product = min_scaling_product * (M_ONE / &
        (M_ONE - cv%lsize(idim) * cv%beta / (M_PI * spacing(idim)) * sin(M_PI * spacing(idim) / cv%lsize(idim))))
    end do
  end subroutine curv_briggs_min_scaling

end module curv_briggs_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
