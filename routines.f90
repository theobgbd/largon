!!! routines LJ MD
module routines
!
use omp_lib
!
implicit none
!
contains
!
subroutine read_input(temperature,tstep,friction, &
                      listcutoff,nstep, natoms,&
                      nconfig,nstat, &
                      wrapatoms, &
                      outputfile,trajfile,statfile, &
                      maxneighbours,idum)
  implicit none
  real,           intent(out) :: temperature
  real,           intent(out) :: tstep
  real,           intent(out) :: friction
  real,           intent(out) :: listcutoff
  integer,        intent(out) :: nstep
  integer,        intent(out) :: natoms
  integer,        intent(out) :: nconfig
  integer,        intent(out) :: nstat
  logical,        intent(out) :: wrapatoms
  integer,        intent(out) :: maxneighbours
  character(256), intent(out) :: outputfile
  character(256), intent(out) :: trajfile
  character(256), intent(out) :: statfile
  integer, intent(out) :: idum
  integer :: iostat
  character(256) :: line,keyword,keyword1
  integer :: i
  logical :: foundsharp
! default values
  temperature=1.0
  tstep=0.005
  friction=0.0
  listcutoff=3.0
  nstep=1
  nconfig=10
  nstat=1
  maxneighbours=1000
  idum=0
  wrapatoms=.false.
  statfile=""
  trajfile=""
  outputfile=""

  do
    read(*,"(a)",iostat=iostat) line
! when the file finishes, exit the loop
    if(iostat/=0) exit
! delete everything past an eventual "#" comment
    foundsharp=.false.
    do i=1,len(line)
      if(line(i:i)=="#") foundsharp=.true.
      if(foundsharp) line(i:i)=" "
    end do
! if the remaining line is empty, skip it
    if(len_trim(line)==0) cycle
! read the first word from line
    read(line,*) keyword
! the second word is then read to the proper variable
    select case(keyword)
    case("temperature")
      read(line,*) keyword1,temperature
    case("tstep")
      read(line,*) keyword1,tstep
    case("friction")
      read(line,*) keyword1,friction
    case("listcutoff")
      read(line,*) keyword1,listcutoff
    case("nstep")
      read(line,*) keyword1,nstep
    case("natoms")
      read(line,*) keyword1,natoms
    case("nconfig")
      read(line,*) keyword1,nconfig,trajfile
    case("nstat")
      read(line,*) keyword1,nstat,statfile
    case("maxneighbours")
      read(line,*) keyword1,maxneighbours
    case("wrapatoms")
      read(line,*) keyword1,wrapatoms
    case("outputfile")
      read(line,*) keyword1,outputfile
    case("trajfile")
      read(line,*) keyword1,trajfile
    case("seed")
      read(line,*) keyword1,idum
      idum=-idum ! idum for ran1() needs to be negative
    case default
! an unknown word will stop the execution
      write(0,*) "Unknown keyword :",trim(keyword)
      stop
    end select
  end do
  if(outputfile=="") then
    write(0,*) "Specify output file"
    stop
  end if
  if(trajfile=="") then
    write(0,*) "Specify traj file"
    stop
  end if
end subroutine read_input
!
! Reading the initial configuration
subroutine read_positions(inputfile,natoms,positions,cell)
! read positions and cell from a file called inputfile
! natoms (input variable) and number of atoms in the file should be consistent
  implicit none
  character(*), intent(in)  :: inputfile
  integer,      intent(in)  :: natoms
  real,         intent(out) :: positions(3,natoms)
  real,         intent(out) :: cell(3)
  integer :: iatom
  character(100) :: atomname
  open(10,file=inputfile)
  read(10,*)
  read(10,*) cell
  do iatom=1,natoms
    read(10,*) atomname,positions(:,iatom)
  end do
! note: atomname is read but not used
  close(10)
end subroutine read_positions
!
!
! Lennard-Jones potential function
function lj_pot(r,sigma,eps) result(u)
    !
    real , intent(in)    :: r, sigma, eps
    real                          :: u 
    !
    u = 4*eps*((sigma/r)**12 - (sigma/r)**6 )
    !
end function lj_pot
!
! Coulombic potential function
function el_pot(r,c1,c2) result(u)
    !
    real , intent(in)   :: r,c1,c2
    real,  parameter 	:: const = 9.0d9   
    real                :: u 
    !
    u = const * (c1*c2) / r
    !
end function el_pot
!
! Function for total energy
function total_energy(natoms,positions,charges, sigma, eps) result(u)
  ! External
  integer,        intent(in)    :: natoms
  real,        	  intent(in)    :: eps, sigma
  real ,          intent(in)    :: positions(3,natoms)      ! atomic positions
  real ,          intent(in)    :: charges  (natoms)      ! atomic positions
  real 			                :: u
  ! Internal
  integer 						:: i,j,k 
  real 		:: r 
  !
  U = 0.0D0
  !
  !$omp parallel do private(i,j,r) reduction(+:u)
  !
  do i=1,natoms
    !
    do j=1, natoms
      !
      if (i .ne. j) then 
        !
        r = sqrt(    (positions(1,i) - positions(1,j))**2 &
                 +   (positions(2,i) - positions(2,j))**2 &
                 +   (positions(3,i) - positions(3,j))**2)
        !
        U = U + el_pot(r,charges(i), charges(j)) + lj_pot(r,sigma,eps)
        ! 
      end if
      !
    end do
    !
  end do
  !
  !$omp end parallel do 
  !  
  return 
  !
end function total_energy
!
! Random position subroutine
subroutine random_conf(natoms, cell,positions)
	integer,        intent(in)    	 :: natoms
  	real ,          intent(inout)    :: positions(3,natoms)
 	real ,          intent(inout)    :: cell(3)

  	integer :: i,k
	!
	do i=1, natoms
    do k=1,3
        !
        call random_number(positions(k,i))
        positions(k,i) = positions(k,i) * cell(k) 
        !
      end do
      !
  end do 
end subroutine random_conf
!
end module routines

