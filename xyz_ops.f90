MODULE xyz_ops

CONTAINS

SUBROUTINE rotate_atoms_x(atom_block, angle_deg, geometric_centre)

   USE xyz_types
   
   IMPLICIT NONE
   
   TYPE(atom) vector
   TYPE(atoms) atom_block
   
   INTEGER :: i
   
   REAL*8, DIMENSION(3,3) :: rotation=0.0
   REAL*8, DIMENSION(3)   :: geometric_centre
   REAL*8, PARAMETER      :: pi=4.0*ATAN(1.0)
   REAL*8                 :: angle_deg, angle_rad
   
   ! Convert to radians
   angle_rad = (pi/180.0)*angle_deg
   
   ! Initialize rotation matrix
   rotation(1,1) = 1.0
   rotation(2,2) = COS(angle_rad)
   rotation(2,3) = -1.*SIN(angle_rad)
   rotation(3,2) = SIN(angle_rad)
   rotation(3,3) = COS(angle_rad)
   
   ! Rotate geometric centre
   geometric_centre = MATMUL(rotation,geometric_centre)
   
   DO i = 1, atom_block%na
      ! Ax = b
      vector%coordinate(:) = atom_block%ats(i,:)
      vector%tmp_coordinate(:) = MATMUL(rotation,vector%coordinate)
   
      atom_block%ats(i,:) = vector%tmp_coordinate
   END DO
   
   RETURN
END SUBROUTINE

SUBROUTINE rotate_atoms_y(atom_block, angle_deg, geometric_centre)

   USE xyz_types
   
   IMPLICIT NONE
   
   TYPE(atom) vector
   TYPE(atoms) atom_block
   
   INTEGER :: i
   
   REAL*8, DIMENSION(3,3) :: rotation=0.0
   REAL*8, DIMENSION(3)   :: geometric_centre
   REAL*8, PARAMETER      :: pi=4.0*ATAN(1.0)
   REAL*8                 :: angle_deg, angle_rad
   
   ! Convert to radians
   angle_rad = (pi/180.0)*angle_deg
   
   ! Initialize rotation matrix
   rotation(1,1) = COS(angle_rad)
   rotation(1,3) = SIN(angle_rad)
   rotation(2,2) = 1.0
   rotation(3,1) = -1.*SIN(angle_rad)
   rotation(3,3) = COS(angle_rad)
   
   ! Rotate geometric centre
   geometric_centre = MATMUL(rotation,geometric_centre)
   
   DO i = 1, atom_block%na
      ! Ax = b
      vector%coordinate(:) = atom_block%ats(i,:)
      vector%tmp_coordinate(:) = MATMUL(rotation,vector%coordinate)
   
      atom_block%ats(i,:) = vector%tmp_coordinate
   END DO
   
   RETURN
END SUBROUTINE
   
SUBROUTINE rotate_atoms_z(atom_block, angle_deg, geometric_centre)
   
   USE xyz_types
   
   IMPLICIT NONE
   
   TYPE(atom) vector
   TYPE(atoms) atom_block
   
   INTEGER :: i
   
   REAL*8, DIMENSION(3,3) :: rotation=0.0
   REAL*8, DIMENSION(3)   :: geometric_centre
   REAL*8, PARAMETER      :: pi=4.0*ATAN(1.0)
   REAL*8                 :: angle_deg, angle_rad
   
   ! Convert to radians
   angle_rad = (pi/180.0)*angle_deg
   
   ! Initialize rotation matrix
   rotation(1,1) = COS(angle_rad)
   rotation(1,2) = -1.*SIN(angle_rad)
   rotation(2,1) = SIN(angle_rad)
   rotation(2,2) = COS(angle_rad)
   rotation(3,3) = 1.0
   
   ! Rotate geometric centre
   geometric_centre = MATMUL(rotation,geometric_centre)
   
   DO i = 1, atom_block%na
      ! Ax = b
      vector%coordinate(:) = atom_block%ats(i,:)
      vector%tmp_coordinate(:) = MATMUL(rotation,vector%coordinate)
   
      atom_block%ats(i,:) = vector%tmp_coordinate
   END DO
   
   RETURN
END SUBROUTINE

SUBROUTINE rotate_atom_z(vector, angle_deg)
! To use the routine the system must already be at the geometric centre
! HINT: use subroutine centre
   USE xyz_types
   
   IMPLICIT NONE
   
   TYPE(atom) vector
   
   INTEGER :: i
   
   REAL*8, DIMENSION(3,3) :: rotation=0.0
   REAL*8, PARAMETER      :: pi=4.0*ATAN(1.0)
   REAL*8                 :: angle_deg, angle_rad
   
   ! Convert to radians
   angle_rad = (pi/180.0)*angle_deg
   
   ! Initialize rotation matrix
   rotation(1,1) = COS(angle_rad)
   rotation(1,2) = -1.*SIN(angle_rad)
   rotation(2,1) = SIN(angle_rad)
   rotation(2,2) = COS(angle_rad)
   rotation(3,3) = 1.0
   
   vector%tmp_coordinate(:) = MATMUL(rotation,vector%coordinate)
   vector%coordinate(:) = vector%tmp_coordinate
   
   RETURN
END SUBROUTINE

SUBROUTINE crys_to_cart(bulk)
! Routine converts crystal coorindates to cartesian coordinates
! in a general way. Converts the coordinate unit from alat
! to Angstrom, stores new coordinates in the bulk object and changes
! the flag to cart
! Josh Elliott - 22/01/2018
   USE xyz_types
   USE cell_data, ONLY : bohrtoang, celldim

   IMPLICIT NONE

   TYPE(atoms)         :: bulk
   TYPE(atom)          :: at
   INTEGER             :: ia
   
   DO ia = 1, bulk%na

      ! Tempoary storage of individual atom
      at%symbol = bulk%atl(ia, 1)
      at%tmp_coordinate(1:3) = bulk%ats(ia,1:3)

      ! conversion of crystal (direct) coordinate to cartesian
      at%coordinate(1) = at%tmp_coordinate(1) * bulk%unit_cell(1,1) + &
                         at%tmp_coordinate(2) * bulk%unit_cell(2,1) + &
                         at%tmp_coordinate(3) * bulk%unit_cell(3,1)  

      at%coordinate(2) = at%tmp_coordinate(1) * bulk%unit_cell(1,2) + &
                         at%tmp_coordinate(2) * bulk%unit_cell(2,2) + &
                         at%tmp_coordinate(3) * bulk%unit_cell(3,2)  

      at%coordinate(3) = at%tmp_coordinate(1) * bulk%unit_cell(1,3) + &
                         at%tmp_coordinate(2) * bulk%unit_cell(2,3) + &
                         at%tmp_coordinate(3) * bulk%unit_cell(3,3)  

      bulk%ats(ia,1:3) = at%coordinate(1:3)*bohrtoang

      WRITE(*,105) ia, bulk%ats(ia,1:3)
   END DO

   bulk%coord_typ = 'cart'

   RETURN
105 FORMAT(I3, 3(F7.3 "   "))
END SUBROUTINE

SUBROUTINE minimum_distance(ap1, ap2, bulk, length, adj_cell)
! Subroutine calculates this distance between two atoms in a structure. If
! adj_cell set to true, this distance is minimised by considering all possible
! periodic replicas.
   USE xyz_types
   USE cell_data, ONLY : bohrtoang
   IMPLICIT NONE

   TYPE(atoms)            :: bulk
   TYPE(atom)             :: ap1, ap2
   LOGICAL                :: adj_cell

   REAL*8, INTENT(OUT)    :: length
   REAL*8                 :: tmp_length

   INTEGER                :: v1, v2, v3

   length = (ap1%coordinate(1)-ap2%coordinate(1))**2 + &
            (ap1%coordinate(2)-ap2%coordinate(2))**2 + &
            (ap1%coordinate(3)-ap2%coordinate(3))**2

   length = SQRT(length)

   IF (adj_cell) THEN ! we should check periodic replicas of ap2
   ! very bad programming we look in all 26 adjacent cells

l1:   DO v1 = -1, 1
l2:      DO v2 = -1, 1
l3:         DO v3 = -1, 1

               ap2%tmp_coordinate(1:3) = ap2%coordinate(1:3) + &
                                         REAL(v1)*bulk%unit_cell(1,1:3)*bohrtoang + & 
                                         REAL(v2)*bulk%unit_cell(2,1:3)*bohrtoang + &
                                         REAL(v3)*bulk%unit_cell(3,1:3)*bohrtoang


               tmp_length = (ap1%coordinate(1)-ap2%tmp_coordinate(1))**2 + &
                        (ap1%coordinate(2)-ap2%tmp_coordinate(2))**2 + &
                        (ap1%coordinate(3)-ap2%tmp_coordinate(3))**2
            
               tmp_length = SQRT(tmp_length)

               IF (tmp_length < length) length = tmp_length

            END DO l3
         END DO l2
      END DO l1

   END IF

   RETURN
END SUBROUTINE

SUBROUTINE centre(cart, shift)

   USE xyz_types

   IMPLICIT NONE

   TYPE(atoms) :: cart
   REAL*8, DIMENSION(2) :: shift
   REAL*8      :: cm(3)
   INTEGER     :: i

   cm(:) = 0.

   DO i = 1, cart%na
      cm(1) = cm(1) + cart%ats(i,1)
      cm(2) = cm(2) + cart%ats(i,2)
      cm(3) = cm(3) + cart%ats(i,3)
   END DO

   cm(:) = cm(:)/REAL(cart%na)

   DO i = 1, cart%na
      cart%ats(i,1) = cart%ats(i,1) - cm(1) - shift(1)
      cart%ats(i,2) = cart%ats(i,2) - cm(2) - shift(2)
      cart%ats(i,3) = cart%ats(i,3) - cm(3)
   END DO

   RETURN
END SUBROUTINE

SUBROUTINE centre_in_box(system, vac)
! The routine works for systems located with centre of mass
! at the origin
   USE xyz_types

   IMPLICIT NONE

   TYPE(atoms) :: system
   REAL*8      :: vac(3), box(3,3)
   INTEGER     :: i

   ! Set up box cell vectors
   system%unit_cell = 0.
   
   system%unit_cell(1,1) = 2.*(vac(1) + MAXVAL(system%ats(:,1)))
   system%unit_cell(2,2) = 2.*(vac(2) + MAXVAL(system%ats(:,2)))
   system%unit_cell(3,3) = 2.*(vac(3) + MAXVAL(system%ats(:,3)))

   DO i = 1, system%na
      system%ats(i,1) = system%ats(i,1) + 0.5*system%unit_cell(1,1)
      system%ats(i,2) = system%ats(i,2) + 0.5*system%unit_cell(2,2)
      system%ats(i,3) = system%ats(i,3) + 0.5*system%unit_cell(3,3)
   END DO

   RETURN
END SUBROUTINE

END MODULE
