PROGRAM flake 

USE cell_data
USE xyz_types
USE xyz_io
USE xyz_ops

   IMPLICIT NONE

   TYPE(atoms)       :: big 
   TYPE(atom)        :: r060, r120
   REAL*8            :: x000, x060, x120, width
   REAL*8, PARAMETER :: deg_060=60., deg_120=120.
   !REAL*8, DIMENSION(2), PARAMETER :: shift=(0.0, 0.0)
   REAL*8, DIMENSION(2), PARAMETER :: shift=(3.0855065590976949E-002,0.06092965570996256) 
   INTEGER           :: i
   LOGICAL           :: exst

   ! Read and check filename
   CALL getarg(1, big%filename)
   INQUIRE(FILE=big%filename, EXIST=exst)
   IF (.NOT. exst) THEN
      WRITE(*,101) 1
      STOP
   END IF

   READ(*,*) width

   ! Read coordinates from xyz
   ! and centre them
   CALL read_xyz(big)
   CALL centre(big, shift)


   ! Open new file for flake coordinates
   OPEN(UNIT=83,file='check_flake.xyz')

   DO i = 1, big%na

      ! Rotate the atomic vector
      r060%coordinate(:) = big%ats(i,:)
      r120%coordinate(:) = big%ats(i,:)
      CALL rotate_atom_z(r060, deg_060)
      CALL rotate_atom_z(r120, deg_120)

      x000 = big%ats(i,2)
      x060 = r060%coordinate(2)       
      x120 = r120%coordinate(2)       
      IF (ABS(x000) .LT. width) THEN 
         IF (ABS(x060) .LT. width) THEN
            IF (ABS(x120) .LT. width) THEN
               WRITE(83,*) big%atl(i,1), big%ats(i,1), big%ats(i,2), big%ats(i,3)
            END IF
         END IF
      END IF

   END DO

   CLOSE(83)

101 FORMAT('Error ',I1,': Input file does not exist')
END PROGRAM
