PROGRAM protonate

USE cell_data
USE xyz_types
USE xyz_io
USE xyz_ops

   IMPLICIT NONE

   TYPE(atoms)       :: flake
   TYPE(atom)        :: r060, r120
   TYPE(atom)        :: h, c, bond
   REAL*8            :: x000, x060, x120, width, vacuum(3)
   REAL*8, PARAMETER :: deg_060=60., deg_120=120., cutoff=1.6, &
                        ch_bond=1.0892830999999994
   INTEGER           :: i,j 
   LOGICAL           :: exst

   ! Read and check filename
   CALL getarg(1, flake%filename)
   INQUIRE(FILE=flake%filename, EXIST=exst)
   IF (.NOT. exst) THEN
      WRITE(*,101) 1
      STOP
   END IF

   READ(*,*) width

   ! Read coordinates from xyz
   ! and centre them
   CALL read_xyz(flake)

   ! Change edge atoms to H
   DO i = 1, flake%na

      ! rotated atomic coordinates
      r060%coordinate(:) = flake%ats(i,:)
      r120%coordinate(:) = flake%ats(i,:)
      CALL rotate_atom_z(r060, deg_060)
      CALL rotate_atom_z(r120, deg_120)

      x000 = flake%ats(i,1)
      x060 = r060%coordinate(1)       
      x120 = r120%coordinate(1)       

      IF (ABS(x000) .GT. width) THEN
               flake%atl(i,1) = 'H'
      ELSE IF (ABS(x060) .GT. width) THEN
               flake%atl(i,1) = 'H'
      ELSE IF (ABS(x120) .GT. width) THEN
               flake%atl(i,1) = 'H'
      END IF

   END DO

   ! Loop over H atoms and correct C-H bond length
   DO i = 1, flake%na

      IF (TRIM(flake%atl(i,1)) .EQ. 'H') THEN
         ! Find the nearest C atom
         h%coordinate(:) = flake%ats(i,:)
         DO j = 1, flake%na

            IF (TRIM(flake%atl(j,1)) .EQ. 'C') THEN
               c%coordinate(:) = flake%ats(j,:)
               bond%coordinate(:) = c%coordinate(:) - h%coordinate
               CALL bond_length(bond)
               IF (bond%magnitude .LT. cutoff) THEN
                  flake%ats(i,:) = c%coordinate(:) -&
                                   ch_bond*(bond%coordinate(:)/bond%magnitude)
   
               END IF
            END IF

         END DO
      END IF

   END DO

   vacuum(1) = 10.
   vacuum(2) = 10.
   vacuum(3) = 30.
   ! Centre in a box
   CALL centre_in_box(flake, vacuum)
   
   flake%filename = 'test.xyz'
   CALL write_xyz(flake, .FALSE.) 
   CALL write_lattice(flake)

101 FORMAT('Error ',I1,': Input file does not exist')
END PROGRAM
