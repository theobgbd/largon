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
  double precision,           intent(out) :: temperature
  double precision,           intent(out) :: tstep
  double precision,           intent(out) :: friction
  double precision,           intent(out) :: listcutoff
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
  double precision,         intent(out) :: positions(3,natoms)
  double precision,         intent(out) :: cell(3)
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
    double precision , intent(in)    :: r, sigma, eps
    double precision                          :: u
    !
    u = 4*eps*((sigma/r)**12 - (sigma/r)**6 )
    !
end function lj_pot
!
! Coulombic potential function
function el_pot(r,c1,c2) result(u)
    !
    double precision , intent(in)   :: r,c1,c2
    double precision,  parameter 	:: const = 9.0d9
    double precision                :: u
    !
    u = const * (c1*c2) / r
    !
end function el_pot
!
! Function for total energy
function total_energy(natoms,positions,charges,cell, sigma, eps) result(u)
  ! External
  integer,        intent(in)    :: natoms
  double precision,        	  intent(in)    :: eps, sigma
  double precision ,          intent(in)    :: positions(3,natoms)      ! atomic positions
  double precision ,          intent(in)    :: charges  (natoms)      ! atomic positions
  double precision     			                :: u, dx, dy, dz, rx, ry, rz
  double precision                 :: cell(3)         ! cell size
  ! Internal
  integer 						:: i,j,k
  double precision 		:: r
  !
  U = 0.0D0
  !
  !$omp parallel do private(i,j,r,rx,ry,rz,dx,dy,dz) reduction(+:U)
  !
  do i=1,natoms
    !
    do j=1, natoms
      !
      if (i .ne. j) then
        !
        dx = positions(1,i) - positions(1,j)
        dy = positions(2,i) - positions(2,j)
        dz = positions(3,i) - positions(3,j)
        rx = dx - cell(1)*dfloat(idint(dx/(cell(1)/2)))
        ry = dy - cell(2)*dfloat(idint(dy/(cell(2)/2)))
        rz = dz - cell(3)*dfloat(idint(dz/(cell(3)/2)))
        !
        r = sqrt( rx**2 + ry**2 + rz**2 )
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
subroutine random_conf(natoms, cell,positions_old, positions)
	integer,        intent(in)    	 :: natoms
  double precision ,          intent(inout)    :: positions(3,natoms), positions_old(3,natoms)
 	double precision ,          intent(inout)    :: cell(3)
  double precision :: newatom(3)
  double precision :: atomid
	integer :: i,k
	!
  call random_number(atomid)
  k = int(atomid*natoms)
  !
	call random_number(newatom)
  newatom(:) = newatom(:) * cell(:)
  !write(*,*) k, newatom
  positions(:,:) = positions_old(:,:)
  positions(:,k) = newatom(:)
  !
end subroutine random_conf
!
end module routines
