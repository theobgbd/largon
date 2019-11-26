!!! Computional CHemistry - PCCP Program
!!! MD-Metropolis code - NVT ensemble
!!! Lennard-Jones Argon
!!! Theo Beigbeder 2019
!
program largon
  !
  use routines
  use omp_lib
  !
  implicit none
  !
  integer              :: natoms          ! number of atoms
  integer              :: nsteps          ! number of configurations to sample

  real,    allocatable :: positions(:,:)      ! atomic positions
  real,    allocatable :: positions_tmp(:,:)  ! tmp atomic positions


  real,    allocatable :: charges(:)      ! charges
  real,    allocatable :: masses(:)       ! masses
  real                 :: cell(3)         ! cell size

  ! input parameters
  ! all of them have a reasonable default value, set in read_input()
  real           :: tstep          ! simulation timestep
  real           :: temperature    ! temperature
  real           :: friction       ! friction for Langevin dynamics (for NVE, use 0)
  real           :: listcutoff     ! cutoff for neighbour list
  integer        :: nstep          ! number of steps
  integer        :: nconfig        ! stride for output of configurations
  integer        :: nstat          ! stride for output of statistics
  integer        :: maxneighbour   ! maximum average number of neighbours per atom
  integer        :: idum           ! seed
  logical        :: wrapatoms      ! if true, atomic coordinates are written wrapped in minimal cell
  character(256) :: outputfile     ! name of file with final configuration (xyz)
  character(256) :: trajfile       ! name of the trajectory file (xyz)
  character(256) :: statfile       ! name of the file with statistics
  !
  real    :: eref        ! energy of the reference  configuration
  real    :: enew        ! energy of the new sampled configuration
  integer :: iatom
  integer :: i,j,k,n, id
  real :: internal_energy,internal_energy_tmp, sigma, eps, bfact, deltaE
  real, parameter :: kb = 1.38d-23
  !
  call read_input(temperature,tstep,friction, &
                  listcutoff,nstep, natoms,nconfig,nstat, &
                  wrapatoms, &
                  outputfile,trajfile,statfile, &
                  maxneighbour,idum)
  !
  write(*,*) "Number of configurations :", nstep
  write(*,*) "Number of atoms          :", natoms
  !
  open(10,file=outputfile)
  open(20,file=trajfile)
  !
  ! Set parallel env
  !$omp parallel
  id = omp_get_thread_num()
  if ( id == 0 ) then
    write ( *, '(a,i4)' ) &
      '  The number of threads available is ', omp_get_num_threads()
  end if
  !$omp end parallel
  !
  ! allocation of dynamical arrays
  allocate(positions(3,natoms))
  allocate(positions_tmp(3,natoms))
  allocate(charges(natoms))
  !
  ! Atom setting override
  sigma  = 1
  eps    = 0.5
  charges(:) = 0.0d0
  cell(:) = (/10,10,10/)
  !
  ! Initialize with a random config
  !
  do i=1, natoms
    do k=1,3
        !
        positions(k,i) = rand() * cell(k)
        !
      end do
      !
  end do 
  !
  !MC-Metropolis sampling algorithm
  !
  config_loop : do n=1, nstep
    !
    internal_energy  = total_energy(natoms,positions,charges, sigma, eps)
    !
    call random_conf(natoms, cell,positions_tmp) ! making a new random config
    !
    ! Compute the total energy for the system
    internal_energy_tmp = total_energy(natoms,positions_tmp,charges, sigma, eps)
    !
    deltaE = internal_energy_tmp - internal_energy
    !
    !write(*,*) " ------------------"
    !write(*,*) "config # :", n 
    !write(*,*) "Internal energy ref.: ", internal_energy
    !write(*,*) "Internal energy new : ", internal_energy_tmp
    !write(*,*) "Energy diff     : ", deltaE
    !
    ! Metropolis test
    !
    if (deltaE .lt. 0) then
      !
      positions(:,:) = positions_tmp(:,:)
      write(*,*) "Selected !"
      !
      write(20,*) natoms
      write(20,*) "frame", n
      do, i=1,natoms
        write(20,*) "Ar", ( positions_tmp(j,i), j=1,3 )
      end do
      !
      write(10,*) n, internal_energy, internal_energy_tmp, deltaE  
      !
    else 
      !
      call random_number(bfact)
      !
      if(bfact .lt. exp(-deltaE/kb*temperature)) then
        positions(:,:) = positions_tmp(i,j) 
        write(*,*) "Selected !"
        !
        write(20,*) natoms
        write(20,*) "frame", n
        !
        do, i=1,natoms
          write(20,*) "Ar", ( positions_tmp(j,i), j=1,3 )
        end do
        !
        write(10,*) n, internal_energy, internal_energy_tmp, deltaE  
      else
        cycle
        !write(*,*) "Rejected !"
      end if
      !
    end if
    !    !
  end do config_loop  
  !
  return
  !
end program largon


