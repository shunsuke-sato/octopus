!! Copyright (C) 2017 Shunsuke A. Sato
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

module lattice_dynamics_oct_m
  implicit none

  private


  type lattice_dynamics_t
    private
    logical          :: move_lattice

  end type ion_dynamics_t

contains

  ! ---------------------------------------------------------
  subroutine lattice_dynamics_init(this)
    type(lattice_dynamics_t), intent(out)   :: this

    PUSH_SUB(lattice_dynamics_init)


    !%Variable MoveIons
    !%Type logical
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable controls whether atoms are moved during a time
    !% propagation run. The default is yes when the ion velocity is
    !% set explicitly or implicitly, otherwise is no.
    !%End
    call parse_variable('MoveLattice', .false., this%move_lattice)
    call messages_print_var_value(stdout, 'MoveLattice', this%move_lattice)



    POP_SUB(lattice_dynamics_init)
  end subroutine lattice_dynamics_init


  ! ---------------------------------------------------------
  subroutine lattice_dynamics_end(this)
    type(ion_dynamics_t), intent(inout) :: this

    PUSH_SUB(lattice_dynamics_end)


    POP_SUB(lattice_dynamics_end)
  end subroutine lattice_dynamics_end


  ! ---------------------------------------------------------


end module lattice_dynamics_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
