MODULE xyz_types

TYPE atoms
   CHARACTER(LEN=100)            :: filename  ! name of the xyz file
   INTEGER                       :: label     ! read and write unit specifier
   INTEGER                       :: na, ntyp  ! number of atoms, number of types of atoms
   INTEGER                       :: ibrav     ! quantum espresso unit cell specifier
   CHARACTER(LEN=4)              :: coord_typ ! cartesian or crystal
   CHARACTER(LEN=4), ALLOCATABLE :: atl(:,:)  ! The atomic label and symbol
   REAL*8, ALLOCATABLE           :: ats(:,:)  ! The array of atoms 
   REAL*8, DIMENSION(3,3)        :: unit_cell ! Unit cell (cartesian)
END TYPE atoms

TYPE atom
   CHARACTER(LEN=4)       :: symbol         ! Atomic symbol
   REAL*8                 :: magnitude      ! Length of vector coordinate
   REAL*8, DIMENSION(1:3) :: coordinate     ! Cartesian coordinate
   REAL*8, DIMENSION(1:3) :: tmp_coordinate ! Temporary coordinate variable
END TYPE atom

CONTAINS

SUBROUTINE bond_length(bond)

   IMPLICIT NONE

   TYPE(atom) :: bond

   bond%magnitude = SQRT(bond%coordinate(1)**2 + bond%coordinate(2)**2 + bond%coordinate(3)**2)

   RETURN
END SUBROUTINE

END MODULE
