Two executables to (1) cut hexagonal flakes from infinite graphene sheets
and (2) protonate the boarder.

# Makefile is setup for gfortran, build executables with

	make flake
	make protonate
	make all

clean directory with

	make trash

# Usage:

	flake.x "my_infinite_sheet.xyz"
	then provide a width for the flake

	You must edit the header of output "check_flake.xyz" by hand (sorry for that)

	protonate.x "check_flake.xyz"
	then provide the "width" of the outermost C atoms, which will be converted to H

# Notes:

	There are variables in the codes (flake.f90, protonate.f90) that can be modified
	on the fly. If your output flake looks weird check the following

	# flake.f90
	shift(:)  (repositions the centre of mass)
	# protonate.f90
	cutoff    (sets the upper bound for C-C bond lengths)
	ch_bond   (sets the length of the new C-H bond)
	vacuum(:) (sets the vacuum size in the x, y and z directions)

