MODULE xyz_io

CONTAINS

SUBROUTINE read_xyz(xyz_in)

   USE xyz_types
   
   IMPLICIT NONE
   
   TYPE(atoms) xyz_in
   INTEGER     :: i
   
   ! Assign label and open file
   xyz_in%label = 20
   OPEN(UNIT=xyz_in%label,FILE=TRIM(xyz_in%filename))
   
   ! Read the number of atoms and allocate arrays
   READ(xyz_in%label,*) xyz_in%na
   ALLOCATE(xyz_in%atl(xyz_in%na,2))
   ALLOCATE(xyz_in%ats(xyz_in%na,3))
   
   ! Read the atoms from file
   DO i = 1, xyz_in%na
      READ(xyz_in%label, *) xyz_in%atl(i,1), xyz_in%ats(i,1:3)
   END DO
   
   CLOSE(xyz_in%label)
   
   RETURN
END SUBROUTINE

SUBROUTINE write_xyz(xyz_out, print_atomic_labels)

   USE xyz_types
   
   IMPLICIT NONE
   
   TYPE(atoms) xyz_out
   INTEGER :: i, j
   LOGICAL :: print_atomic_labels
   
   
   ! Assign write label and open file
   xyz_out%label = 30
   OPEN(UNIT=xyz_out%label,FILE=TRIM(xyz_out%filename))
   
   ! Write the data to file
   WRITE(xyz_out%label,*) xyz_out%na
   WRITE(xyz_out%label,*)
   
   ! whether to print atomic labels or user labels
   j = 1
   IF (print_atomic_labels .EQV. .TRUE.) THEN
      j = 2
   END IF
   
   DO i = 1, xyz_out%na
      WRITE(xyz_out%label, *) xyz_out%atl(i,j), xyz_out%ats(i,1:3)
   END DO
   
   CLOSE(xyz_out%label)

   RETURN
END SUBROUTINE

SUBROUTINE composition(bulk)

   USE xyz_types, ONLY : atoms

   TYPE(atoms) :: bulk
   INTEGER     :: ia

   bulk%ntyp = 0
   DO ia = 1, bulk%na
      
   END DO

   RETURN   
END SUBROUTINE

SUBROUTINE write_lattice(system)

   USE xyz_types, ONLY : atoms
   
   IMPLICIT NONE

   TYPE(atoms) :: system
   INTEGER     :: i, j

   OPEN(unit=45,file=TRIM(system%filename)//'.lat')

   DO i = 1,3
      WRITE(45,*) (system%unit_cell(i,j), j=1,3)
   END DO

   CLOSE(45)
   RETURN
END SUBROUTINE

END MODULE
