!!! Computional Chemistry - PCCP Program
!!! MD-Metropolis code - NVT sampling
!!! Lennard-Jones Argon
!!! Theo Beigbeder 2019
!!!
program largon
  !
  use routines
  use omp_lib
  !
  implicit none
  !
  integer              :: natoms          ! number of atoms
  integer              :: nsteps          ! number of configurations to sample

  double precision,    allocatable :: positions(:,:)      ! atomic positions
  double precision,    allocatable :: positions_tmp(:,:)  ! tmp atomic positions


  double precision,    allocatable :: charges(:)      ! charges
  double precision,    allocatable :: masses(:)       ! masses
  double precision                 :: cell(3)         ! cell size

  ! input parameters
  ! all of them have a reasonable default value, set in read_input()
  double precision           :: tstep          ! simulation timestep
  double precision           :: temperature    ! temperature
  double precision           :: friction       ! friction for Langevin dynamics (for NVE, use 0)
  double precision           :: listcutoff     ! cutoff for neighbour list
  double precision           :: init_energy_lim! threshold engery for the inital configuration energy
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
  double precision    :: eref        ! energy of the reference  configuration
  double precision    :: enew        ! energy of the new sampled configuration
  integer :: iatom
  integer :: i,j,k,n,m, id
  double precision :: internal_energy,internal_energy_tmp, sigma, eps, bfact, deltaE
  double precision, parameter :: kb = 1.38d-23
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
  !
  !$omp parallel
  !$ id = omp_get_thread_num()
  !$ if ( id == 0 ) then
  !$  write ( *, '(a,i4)' ) &
  !$    '  The number of threads available is ', omp_get_num_threads()
  !$ end if
  !$omp end parallel
  !
  ! allocation of dynamical arrays
  allocate(positions(3,natoms))
  allocate(positions_tmp(3,natoms))
  allocate(charges(natoms))
  !
  ! Atom setting override
  sigma  = 1
  eps    = 2
  charges(:) = 0.0d0
  cell(:) = (/10,10,10/)
  !
  ! Initialize with a random config
  !
  write(*,*) 'Selecting initial configuration ...'
  !
  init_energy_lim = 1.0D20
  !
  init_conformity : do ! infinite loop until a correct config is found
    !
    do i=1, natoms
      do k=1,3
        !
        call random_number(positions)
        positions(1,:) = positions(1,:) * cell(1)
        positions(2,:) = positions(2,:) * cell(2)
        positions(3,:) = positions(3,:) * cell(3)
        !
      end do
      !
    end do
    !
    internal_energy  = total_energy(natoms,positions,charges, cell, sigma, eps)
    !
    if (internal_energy .gt. init_energy_lim ) then
      !write(*,*) 'Initial configuration rejected !'
      cycle
    else 
      write(*,*) 'Initial configuration selected !'
      write(*,*) 'Initial energy : ', internal_energy 
      exit
    end if
    !
  end do init_conformity
  !
  !MC-Metropolis sampling algorithm
  !
  n = 0
  m = 0
  !
  config_loop : do 
    !
    internal_energy  = total_energy(natoms,positions,charges, cell, sigma, eps)
    !
    call random_conf(natoms, cell,positions,positions_tmp) ! making a new random config
    !
    ! Compute the total energy for the system
    internal_energy_tmp = total_energy(natoms,positions_tmp,charges, cell, sigma, eps)
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
      !write(*,*) "Selected !"
      positions(:,:) = positions_tmp(:,:)
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
        positions(:,:) = positions_tmp(:,:)
        !write(*,*) "Selected !"
        !
        write(20,*) natoms
        write(20,*) "frame", n
        !
        do, i=1,natoms
          write(20,*) "Ar", ( positions_tmp(j,i), j=1,3 )
        end do
        !
        write(10,*) n, internal_energy, internal_energy_tmp, deltaE
        !
      else
        !
        cycle
        !
      end if
      !
    end if
    !
    if (n .gt. nstep) then
      exit
    end if
    !
    n = n + 1
    m = m + 1
    !
  end do config_loop
  !
  write(*,*) "Done !"
  !
  ! deallocation of dynamical arrays
  deallocate(positions)
  deallocate(positions_tmp)
  deallocate(charges)
  !
  return
  !
end program largon
