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

program octopus
  use global_oct_m
  use calc_mode_par_oct_m
  use command_line_oct_m
  use io_oct_m
  use loct_oct_m
  use mpi_oct_m
  use parser_oct_m
  use profiling_oct_m
  use run_oct_m
  use string_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use messages_oct_m

  use io_binary_oct_m
  use io_function_oct_m

  implicit none

  character(len=256) :: config_str
  integer :: inp_calc_mode, ierr
  type(block_t) :: blk

  integer :: np,np_part,nst,nkp,nst_gs
  real(8) :: aLx,aLy,aLz,Omega
  character(999) :: dir_name,filename
  real(8),allocatable :: occ(:,:),occ_gs(:,:),kvec(:,:),kweight(:)
  complex(8),allocatable :: wfn(:,:),wfn_gs(:,:),read_ff(:)
  integer :: ik,ist,ist_gs,offset=0
  integer :: inum, inum_gs
  real(8),allocatable :: occ_ex(:,:), num_ex_elec_each_k(:)
  real(8) :: num_elec, num_ex_elec
  real(8) :: ss


  call global_init()
  call messages_init()


!  dir_name="/scratch/ssato/shift_current/BaTiO3/nonlinear_response/ps_b_NK8NL32_p1d10_w4.0ev_dt_0.03_p/"
  read(*,*)dir_name

  call read_mesh
  call read_state_info
  allocate(occ(nst,nkp),occ_gs(nst_gs,nkp),kvec(3,nkp),kweight(nkp))
  allocate(occ_ex(nst_gs,nkp),num_ex_elec_each_k(nkp))
  allocate(wfn(np,nst),wfn_gs(np,nst_gs))
  allocate(read_ff(np))
  call read_occ_kvec

  inum=0; inum_gs=0
  do ik=1,nkp
     write(*,"(A,2x,I7,A,I7)")"ik=",ik,"/",nkp
! td
     do ist=1,nst
        inum = inum+1
        write(filename,'(i10.10)') inum
        filename=trim(dir_name)//"restart/td/"//trim(filename)//".obf"
        call io_binary_read(trim(filename), np, read_ff, ierr, offset = offset)
        wfn(:,ist)=read_ff(:)
     end do

! gs
     do ist=1,nst_gs
        inum_gs = inum_gs+1
        write(filename,'(i10.10)') inum_gs
        filename=trim(dir_name)//"restart/gs/"//trim(filename)//".obf"
        call io_binary_read(trim(filename), np, read_ff, ierr, offset = offset)
        wfn_gs(:,ist)=read_ff(:)
     end do

     ss = 0d0
     do ist_gs=1,nst_gs
        ss = 0d0
        do ist=1,nst
           ss = ss + abs( sum( conjg(wfn_gs(:,ist_gs))*wfn(:,ist) )/np*Omega )**2
        end do
        occ_ex(ist_gs,ik) = ss
     end do
  end do

  num_elec = dble(nst)
  ss = 0d0
  do ik = 1,nkp
     ss = ss + sum(occ_ex(1:nst,ik))*kweight(ik)
     num_ex_elec_each_k(ik) = num_elec - sum(occ_ex(1:nst,ik)) 
  end do
  num_ex_elec = num_elec - ss

  open(20,file="occ_ex.dat")
  write(20,"(A,2x,2es26.16e3)")"num_elec, num_ex_elec",num_elec,num_ex_elec
  do ik = 1,nkp
     do ist = 1,nst_gs
        write(20,"(I7,2x,I7,e26.16e3)")ik,ist,occ_ex(ist,ik)
     end do
  end do
  close(20)


  open(20,file="num_ex_at_k.dat")
  do ik = 1,nkp,2
        write(20,"(I7,2x,999e26.16e3)")ik ,num_ex_elec_each_k(ik),num_ex_elec_each_k(ik+1) &
             ,min(num_ex_elec_each_k(ik)/num_ex_elec_each_k(ik+1),num_ex_elec_each_k(ik+1)/num_ex_elec_each_k(ik))
  end do
  close(20)
!  open(20,file="/scratch/ssato/shift_current/BaTiO3/nonlinear_response/ps_b_NK8NL32_p1d10_w4.0ev_dt_0.03_p/restart/td/0000000001.obf",form="unformatted")
!  read(20)wfn(1:np+1,1)
!  close(20)
!  write(*,*)sum( abs(wfn(1:np,1))**2)/np*Omega
!  write(*,*)wfn(1,1)
!  write(*,*)wfn(2,1)
!  write(*,*)wfn(3,1)
!  write(*,*)wfn(4,1)
!  write(*,*)wfn(5,1)
!  write(*,*)wfn(6,1)

  call messages_end()
  call global_end()

contains
!===================================================================================================
  subroutine read_mesh
    implicit none
    character(999) :: file_name,ctmp

    file_name=trim(dir_name)//"restart/td/mesh"
    open(20,file=file_name)
    read(20,*); read(20,*); read(20,*); read(20,*)
    read(20,*)ctmp, np
    read(20,*)ctmp, np_part
    read(20,*); read(20,*); read(20,*); read(20,*); read(20,*); read(20,*)
    read(20,*); read(20,*); read(20,*); read(20,*); read(20,*); read(20,*)
    read(20,*); read(20,*); read(20,*); read(20,*)
    read(20,*)ctmp, aLx, aLy,aLz
    aLx = 2d0*aLx; aLy = 2d0*aLy; aLz = 2d0*aLz
    close(20)
    Omega = aLx*aLy*aLz

    write(*,"(A,2x,I7)")"# np=",np
    write(*,"(A,2x,I7)")"# np_part=",np_part
    write(*,"(A,2x,3es16.6e3)")"# aLx,aLy,aLz=",aLx,aLy,aLz
    write(*,"(A,2x,es16.6e3)")"# Omega=",Omega
  end subroutine read_mesh
!===================================================================================================
  subroutine read_state_info
    implicit none
    character(999) :: file_name,ctmp
    
    file_name=trim(dir_name)//"restart/td/states"
    open(20,file=file_name)
    read(20,*)ctmp,nst
    read(20,*)
    read(20,*)ctmp,nkp
    
    file_name=trim(dir_name)//"restart/gs/states"
    open(20,file=file_name)
    read(20,*)ctmp,nst_gs
    write(*,"(A,2x,3I5)")"# nst, nst_gs, nkp=",nst,nst_gs, nkp
    
  end subroutine read_state_info
!===================================================================================================
  subroutine read_occ_kvec
    implicit none
    integer :: ik,ist
    real(8) :: rtmp(999)
    character(999) :: file_name,ctmp(999)

! time-dependent
    file_name=trim(dir_name)//"restart/td/occs"
    open(20,file=file_name)
    read(20,*); read(20,*)
    do ik=1,nkp
       do ist=1,nst
          read(20,*)occ(ist,ik),ctmp(1),rtmp(1),ctmp(2),rtmp(2),ctmp(3) &
               ,kvec(1,ik),ctmp(4),kvec(2,ik),ctmp(5),kvec(3,ik),ctmp(6),kweight(ik)
       end do
    end do

    do ik=1,nkp
       do ist=1,nst
          write(*,"(2I5,2x,999e16.6e3)")ik,ist,occ(ist,ik),kvec(:,ik),kweight(ik)
       end do
    end do

! ground state
    file_name=trim(dir_name)//"restart/gs/occs"
    open(20,file=file_name)
    read(20,*); read(20,*)
    do ik=1,nkp
       do ist=1,nst_gs
          read(20,*)occ_gs(ist,ik),ctmp(1),rtmp(1),ctmp(2),rtmp(2),ctmp(3) &
               ,kvec(1,ik),ctmp(4),kvec(2,ik),ctmp(5),kvec(3,ik),ctmp(6),kweight(ik)
       end do
    end do

!    do ik=1,nkp
!       do ist=1,nst
!          write(*,"(2I5,2x,999e16.6e3)")ik,ist,occ(ist,ik),kvec(:,ik),kweight(ik)
!       end do
!    end do
!
!    do ik=1,nkp
!       do ist=1,nst_gs
!          write(*,"(2I5,2x,999e16.6e3)")ik,ist,occ_gs(ist,ik),kvec(:,ik),kweight(ik)
!       end do
!    end do
    
  end subroutine read_occ_kvec






!  call getopt_init(ierr)
!  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
!  if(ierr  ==  0) call getopt_octopus(trim(config_str))
!  call getopt_end()
!
!  call global_init()
!  call messages_init()
!
!  !%Variable ReportMemory
!  !%Type logical
!  !%Default no
!  !%Section Execution::Debug
!  !%Description
!  !% If true, after each SCF iteration <tt>Octopus</tt> will print
!  !% information about the memory the code is using. The quantity
!  !% reported is an approximation to the size of the heap and
!  !% generally it is a lower bound to the actual memory <tt>Octopus</tt> is
!  !% using.
!  !%End
!  call parse_variable('ReportMemory', .false., conf%report_memory)
!
!  ! need to find out calc_mode already here since some of the variables here (e.g.
!  ! periodic dimensions) can be different for the subsystems
!
!  !%Variable CalculationMode
!  !%Type integer
!  !%Default gs
!  !%Section Calculation Modes
!  !%Description
!  !% Decides what kind of calculation is to be performed.
!  !%Option gs 01
!  !% Calculation of the ground state.
!  !%Option unocc 02
!  !% Calculation of unoccupied/virtual KS states. Can also be used for a non-self-consistent
!  !% calculation of states at arbitrary k-points, if <tt>density.obf</tt> from <tt>gs</tt>
!  !% is provided in the <tt>restart/gs</tt> directory.
!  !%Option td 03
!  !% Time-dependent calculation (experimental for periodic systems).
!  !%Option go 05
!  !% Optimization of the geometry.
!  !%Option opt_control 07
!  !% Optimal control.
!  !%Option em_resp 08
!  !% Calculation of the electromagnetic response: electric
!  !% polarizabilities and hyperpolarizabilities and magnetic
!  !% susceptibilities (experimental for periodic systems).
!  !%Option casida 09
!  !% Excitations via Casida linear-response TDDFT; for finite systems only.
!  !%Option vdw 11
!  !% Calculate van der Waals coefficients.
!  !%Option vib_modes 12
!  !% Calculation of the vibrational modes.
!  !%Option one_shot 14
!  !% Obsolete. Use <tt>gs</tt> with <tt>MaximumIter = 0</tt> instead.
!  !%Option kdotp 15
!  !% Calculation of effective masses by <math>\vec{k} \cdot \vec{p}</math> perturbation theory (experimental).
!  !%Option dummy 17
!  !% This calculation mode does nothing. Useful for debugging, testing and benchmarking.  
!  !%Option invert_ks 18
!  !% Invert the Kohn-Sham equations (experimental).
!  !%Option recipe 99
!  !% Prints out a tasty recipe.
!  !%End
!  if(parse_block('CalculationMode', blk) == 0) then
!    call messages_write('The datasets mode has been deprecated,', new_line = .true.)
!    call messages_write('please use several Octopus runs.')
!    call messages_fatal()
!  end if
!
!  call parse_variable('CalculationMode', CM_GS, inp_calc_mode)
!  if(.not.varinfo_valid_option('CalculationMode', inp_calc_mode)) call messages_input_error('CalculationMode')
!
!  ! Now we can initialize the I/O
!  call io_init()
!
!  call calc_mode_par_init()
!  
!  ! now we declare octopus as running
!  call io_switch_status('running')
!  
!  call profiling_init()
!  
!  call print_header()
!  
!  ! now we really start
!  call run(inp_calc_mode)
!  
!#if defined(HAVE_MPI)
!  ! wait for all processors to finish
!  call MPI_Barrier(mpi_world%comm, mpi_err)
!#endif
!  
!  ! run finished successfully
!  call io_switch_status('finished')
!  call io_end()
!  
!  call profiling_end()
!  
!  call calc_mode_par_end()
!  
!  call print_date("Calculation ended on ")
!  call print_walltime()
!
!  call messages_end()
!  call global_end()
!
!contains
!
!  subroutine print_walltime()
!    integer :: days, hours, min, sec, usec
!
!    call loct_gettimeofday(sec, usec)
!    call epoch_time_diff(sec, usec)
!    
!    days  = sec / 86400
!    hours = (sec / 3600) - (days * 24)
!    min   = (sec / 60) - (days * 1440) - (hours * 60)
!    sec   = modulo(sec, 60)
!
!    message(2) = ''
!    if(days  > 0) write(message(2), '(i3,a)') days, ' days,'
!    if(hours > 0.or.message(2) /= '') &
!      write(message(2), '(a,1x,i2.2,a)') trim(message(2)), hours, 'h'
!    if(min   > 0.or.message(1) /= '') &
!      write(message(2), '(a,1x,i2.2,a)') trim(message(2)), min, 'm'
!    write(message(2), '(a,1x,i2.2,a,i3,a)') trim(message(2)), sec, '.', usec/1000, 's'
!    message(1) = str_center('Walltime: ' // trim(message(2)), 70)
!    call messages_info(1)
!
!  end subroutine print_walltime

end program octopus

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
