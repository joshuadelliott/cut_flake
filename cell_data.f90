MODULE cell_data

REAL*8, PARAMETER      :: bohrtoang = 0.529177
REAL*8, DIMENSION(1:6) :: celldim

CONTAINS

SUBROUTINE create_vectors(bulk)

   USE xyz_types

   IMPLICIT NONE

   TYPE(atoms)            :: bulk
   REAL*8                 :: term1, term2, sr2, sr3
!   REAL*8, DIMENSION(1:6) :: celldim

   IF (bulk%ibrav == 5) THEN
   
      WRITE(*,*) '   ibrav = 5 detected: provide'
      WRITE(*,*) '   celldim(1) celldim(4)'
   
      READ(*,*) celldim(1), celldim(4)
   
      ! From espresso latgen.f90
      IF (celldim(4) <= -0.5 .OR. celldim(4) >= 1.0) THEN
         WRITE(*,102) 2
         STOP
      END IF
   
      term1 = SQRT(1.0 + 2.0*celldim(4))
      term2 = SQRT(1.0 - celldim(4))
   
      sr2 = SQRT(2.)
      sr3 = SQRT(3.)
   
      bulk%unit_cell(1,1) = celldim(1)*term2/sr2
      bulk%unit_cell(2,2) = sr2*celldim(1)*term2/sr3
      bulk%unit_cell(2,3) = celldim(1)*term1/sr3
      bulk%unit_cell(1,2) = -1.*bulk%unit_cell(1,1)/sr3
      bulk%unit_cell(1,3) =  1.*bulk%unit_cell(2,3)
      bulk%unit_cell(2,1) =  0.0 
      bulk%unit_cell(3,1) = -1.*bulk%unit_cell(1,1)
      bulk%unit_cell(3,2) =  1.*bulk%unit_cell(1,2)
      bulk%unit_cell(3,3) =  1.*bulk%unit_cell(2,3)
   
   END IF

   bulk%unit_cell(:,:) = bulk%unit_cell(:,:)
   
   ! Print Unit cell vectors in Angstroms
   WRITE(*,*) 'Unit cell vectors (Ang):'
   WRITE(*,103) bulk%unit_cell(1,1:3)*bohrtoang
   WRITE(*,103) bulk%unit_cell(2,1:3)*bohrtoang
   WRITE(*,103) bulk%unit_cell(3,1:3)*bohrtoang
   
   102 FORMAT('Error ',I1,': celldim(4) is incorrect')
   103 FORMAT(3(F10.7 '  '))
   RETURN
END SUBROUTINE

END MODULE
