# many-body-lib
A series of .hpp files which allows the generation of different Hamiltonians used in many-body phyiscs, using either a fermionic, bosonic or tensorproduct basis
# Name defenitions:
State is object containg one lattice and the identification
In a state:
LatticeIt iterates throug a lattice= a vector with particles
In a Basis:
Lattice is an object containg lattice array with elctrons + dim + id
BasisType is the tuple containg lattice + id + position. The type of the basis
Latticeit is same
Basis is set od basistypes
Basisit iterator poiting to basis element

 In tensor product:
 BasisType is the tuple containing left and right identiefiers and position= id
 LTLattice and R are the types of the left and right basis'
 Basis set containing basistype
 basisIt iterates through basis
 