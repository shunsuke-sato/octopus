program main
  implicit none

  real(8),parameter :: pi = 4d0*atan(1d0)
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
  
  real(8),parameter  :: beta = 1d0/(0.01d0/ev)
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
        write(20,"(I7,2x,I7,999e26.16e3)")ik,ist,eps_ex(ist,ik)&
             ,occ_ex(ist,ik),occ_zero(ist,ik),kweight(ik)
     end do
  end do
  close(20)

  call Fermi_Dirac_distribution


  call write_DOS

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

    close(20)
    
    file_name=trim(dir_name)//"restart/gs/states"
    open(20,file=file_name)
    read(20,*)ctmp,nst_gs
    write(*,"(A,2x,3I5)")"# nst, nst_gs, nkp=",nst,nst_gs, nkp

    close(20)
    
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

    close(20)

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

    close(20)

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

    eps_Fermi = mu

    ph_num = 0d0
    do ik = 1,nkp
       ph_num = ph_num + sum(abs(occ_ex(:,ik)-occ_zero(:,ik)))*kweight(ik)
    end do
    ph_num = 0.5d0*ph_num

    write(*,*)"ph_num",ph_num

  end subroutine Fermi_Dirac_distribution

  subroutine write_dos
    implicit none
    integer,parameter :: Nw = 4096
    real(8),parameter :: gamma = 0.1d0/ev
    real(8) :: eps_max, eps_min
    real(8) :: ww_max, ww_min, dw, ww
    real(8) :: dos_zero(0:Nw),dos_ex(0:Nw)
    real(8) :: dos_zero_s(0:Nw),dos_ex_s(0:Nw)
    integer :: ik, ist, iw, iw2
    real(8) :: ff

    eps_max = maxval(eps_ex)
    eps_min = minval(eps_ex)

    ww_min = eps_min - 0.5d0*(eps_max-eps_min)
    ww_max = eps_max + 0.5d0*(eps_max-eps_min)
    dw = (ww_max - ww_min)/Nw

    dos_zero = 0d0
    do ik = 1,nkp
       do ist = 1,nst_gs
          iw = aint( (eps_ex(ist,ik)-ww_min)/dw )
          dos_zero(iw) = dos_zero(iw) + occ_zero(ist,ik)*kweight(ik)
          dos_ex(iw)   = dos_ex(iw)   + occ_ex(ist,ik)  *kweight(ik)

       end do
    end do

    dos_zero = dos_zero/dw
    dos_ex   = dos_ex/dw

    dos_zero_s = 0d0
    dos_ex_s   = 0d0

    do iw = 0,Nw

       do iw2 =0,Nw
          ww = (iw-iw2)*dw
          ff = (gamma/pi)/(ww**2+gamma**2)*dw
          dos_zero_s(iw2) = dos_zero_s(iw2) + dos_zero(iw)*ff 
          dos_ex_s(iw2)   = dos_ex_s(iw2)   + dos_ex(iw)*ff 

       end do

    end do

    write(*,*)"dos-raw integral",sum(dos_zero)*dw,sum(dos_ex)*dw
    write(*,*)"dos-sme integral",sum(dos_zero_s)*dw,sum(dos_ex_s)*dw
    
    open(20,file='nex_w_dist.out')
    do iw = 0,Nw
       ww = ww_min + dw*iw - eps_Fermi
       write(20,"(999e26.16e3)")ww,dos_ex_s(iw),dos_ex_s(iw)-dos_zero_s(iw)

    end do
    close(20)


  end subroutine write_dos


end program main
