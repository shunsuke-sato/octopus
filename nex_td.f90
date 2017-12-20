program main
  implicit none

  real(8),parameter :: ev=27.21138505d0

  character(len=256) :: config_str
  integer :: inp_calc_mode, ierr

  integer :: np,np_part,nst,nkp,nst_gs
  real(8) :: aLx,aLy,aLz,Omega
  real(8)  :: rlattice(3,3)
  character(999) :: dir_name,filename
  real(8),allocatable :: occ(:,:),occ_gs(:,:),kvec(:,:),kweight(:)

  integer :: ik,ist,ist_gs,offset=0
  integer :: inum, inum_gs
  real(8),allocatable :: occ_ex(:,:), num_ex_elec_each_k(:)
  real(8),allocatable :: eps_ex(:,:),occ_zero(:,:)
  real(8) :: num_elec, num_ex_elec
  real(8) :: ss
  integer :: ik_t, ist_t
  
  real(8),parameter  :: beta = 1d0/(0.1d0/ev)
  real(8) :: eps_Fermi
  character(999) :: ctmp


  dir_name = "./"

  call read_mesh
  call read_state_info
  allocate(occ(nst,nkp),occ_gs(nst_gs,nkp),kvec(3,nkp),kweight(nkp))
  allocate(occ_ex(nst_gs,nkp),num_ex_elec_each_k(nkp),eps_ex(nst_gs,nkp))
  allocate(occ_zero(nst_gs,nkp))
  call read_occ_kvec

  open(20,file="occ_ex.dat")
  read(20,*)
  do ik = 1,nkp
     do ist = 1,nst_gs
        read(20,*)ik_t,ist_t,occ_ex(ist,ik)
     end do
  end do
  close(20)


  open(20,file="static/eigenvalues")
  read(20,*);   read(20,*)
  read(20,*);   read(20,*)
  read(20,*)

  do ik = 1,nkp
  
     read(20,*)
     do ist = 1,nst_gs
        read(20,*)ist_t,ctmp,eps_ex(ist,ik)
     end do
  end do
  
  close(20)


  open(20,file='nex_eps_dist.out')
  do ik = 1,nkp
     do ist = 1,nst_gs
        write(20,"(I7,2x,I7,2e26.16e3)")ik,ist,occ_ex(ist,ik),eps_ex(ist,ik)
     end do
  end do
  close(20)

  call Fermi_Dirac_distribution

contains
  subroutine read_mesh
    implicit none
    character(999) :: file_name,ctmp
    real(8) :: rvec(3)

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
    read(20,*)
    read(20,*)
    read(20,*)ctmp, rlattice(1,1),rlattice(2,1),rlattice(3,1)
    read(20,*)ctmp, rlattice(1,2),rlattice(2,2),rlattice(3,2)
    read(20,*)ctmp, rlattice(1,3),rlattice(2,3),rlattice(3,3)

    close(20)

    rlattice(:,1) = rlattice(:,1)*aLx
    rlattice(:,2) = rlattice(:,2)*aLy
    rlattice(:,3) = rlattice(:,3)*aLz

    rvec(1) = rlattice(2,1)*rlattice(3,2) - rlattice(3,1)*rlattice(2,2)
    rvec(2) = rlattice(3,1)*rlattice(1,2) - rlattice(1,1)*rlattice(3,2)
    rvec(3) = rlattice(1,1)*rlattice(2,2) - rlattice(2,1)*rlattice(1,2)

    Omega = abs(sum(rvec(:)*rlattice(:,3)))

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

!    do ik=1,nkp
!       do ist=1,nst
!          write(*,"(2I5,2x,999e16.6e3)")ik,ist,occ(ist,ik),kvec(:,ik),kweight(ik)
!       end do
!    end do

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

  subroutine Fermi_Dirac_distribution
    implicit none
    integer,parameter :: niter_max = 500
    real(8) :: mu_min, mu_max, mu
    real(8) :: tot_elec_num, elec_num,ph_num
    integer :: iter

    tot_elec_num = 0d0
    do ik = 1,nkp
       tot_elec_num = tot_elec_num + sum(occ_ex(:,ik))*kweight(ik)
    end do

    write(*,*)"Total # of electrons =",tot_elec_num

    mu_min = minval(occ_ex)
    mu_max = maxval(occ_ex)

    iter = 0
    do 
       iter = iter + 1
       mu = 0.5d0*(mu_min+mu_max)

       do ik = 1,nkp
          do ist = 1,nst_gs

             occ_zero(ist,ik) = 2d0/(exp((eps_ex(ist,ik)-mu)*beta)+1d0)

          end do
       end do

       elec_num = 0d0
       do ik = 1,nkp
          elec_num = elec_num + sum(occ_zero(:,ik))*kweight(ik)
       end do

       if(elec_num < tot_elec_num)then
          mu_min = mu
       else
          mu_max = mu
       end if
       if(mu_max-mu_min < 1d-12)exit
       
    end do
    write(*,*)"Fermi-Dirac distirbution"
    write(*,*)"# of iteration",iter
    write(*,*)"mu_max - mu_min (ev)",(mu_max-mu_min)*ev
    write(*,*)"elec_num",elec_num,tot_elec_num


    ph_num = 0d0
    do ik = 1,nkp
       ph_num = ph_num + sum(abs(occ_ex(:,ik)-occ_zero(:,ik)))*kweight(ik)
    end do
    ph_num = 0.5d0*ph_num

    write(*,*)"ph_num",ph_num

  end subroutine Fermi_Dirac_distribution


end program main
