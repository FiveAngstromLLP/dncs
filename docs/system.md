## System

The `System` struct represents a molecular system composed of particles (atoms) and encapsulates information about the system's topology, force field parameters, and interactions. It is designed to handle the initialization and management of molecular systems for simulations or analyses.

### Fields

- `seq: String`: Represents the sequence of the molecular system, such as a protein or nucleic acid sequence.

- `forcefield: ForceField`: Contains the force field parameters used for the system, including atom types, residues, and non-bonded interactions.

- `particles: Particles`: A vector of `Atom` instances representing all the atoms in the system.

- `dihedral: Vec<(Atom, Atom, Atom, Atom)>`: A list of tuples representing dihedral angles in the system. Each tuple contains four atoms that define a dihedral.

- `dihedral_angle: Vec<(Atom, Atom, Atom, Atom)>`: Similar to `dihedral`, but specifically used for angle calculations.

- `firstbonded: Vec<Particles>`: A list where each element is a vector of atoms that are directly bonded (first neighbors) to a specific atom.

- `secondbonded: Vec<Particles>`: A list where each element is a vector of atoms that are second-bonded (bonded to first-bonded atoms) to a specific atom.

- `nonbonded: Vec<Particles>`: A list where each element is a vector of atoms that are non-bonded to a specific atom.

- `bonded1_4: Vec<Particles>`: A list where each element is a vector of atoms that are 1-4 bonded (bonded through three bonds) to a specific atom.

- `hydrogen: Vec<(Atom, Atom)>`: A list of atom pairs representing potential hydrogen bonds in the system.

### Methods

- `new(seq: &str, forcefield: ForceField) -> Self`:
  Creates a new `System` instance from a given sequence string and force field. It initializes the particles based on the sequence and sets up the data structures for bonded and non-bonded interactions.

- `from_pdb(file: &str, forcefield: ForceField) -> Self`:
  Creates a new `System` instance by parsing a PDB file. It extracts atom information from the file and initializes the system accordingly.

- `get_atomtype(&mut self)`:
  Assigns atom types to each atom in the system based on the force field definitions. Updates atom properties like `typeid`, `atomtype`, `epsilon`, `sigma`, and `charge`.

- `init_parameters(&mut self)`:
  Initializes various parameters for the system by calling methods to assign atom types, determine neighbors, set energy parameters, and identify hydrogen bonds.

- `get_energyparameters(&mut self)`:
  Sets the energy parameters (`epsilon` and `sigma`) for each atom according to its type and position in the sequence (e.g., special handling for terminal atoms).

- `get_dihedral(&mut self)`:
  Calculates the dihedral angles in the system, including both backbone and side-chain dihedrals, and populates the `dihedral` and `dihedral_angle` lists.

- `get_hydrogen_bonded(&mut self)`:
  Identifies pairs of atoms that can form hydrogen bonds and populates the `hydrogen` list. Searches for hydrogen donors and acceptors within the system.

- `get_neighbours(&mut self)`:
  Determines the bonded and non-bonded neighbors for each atom by creating `Neighbour` instances and updating the neighbor lists.

- `to_pdb(&self, filename: &str)`:
  Exports the system's atom information to a PDB file, allowing visualization or further analysis using standard molecular modeling tools.

## Neighbor

The `Neighbour` struct is responsible for calculating and storing neighbor atoms for a specific atom within a molecular system. It categorizes neighbors into different levels based on bonding and spatial relationships.

### Fields

- `atom: Atom`: The reference atom for which neighbors are being calculated.

- `polymer: Particles`: A vector of all atoms in the system, representing the molecular chain.

- `firstbonded: Particles`: A vector of atoms that are directly bonded (first neighbors) to the reference atom.

- `secondbonded: Particles`: A vector of atoms that are bonded to the first-bonded atoms but not directly to the reference atom (second neighbors).

- `nonbonded: Particles`: A vector of atoms that are not directly bonded to the reference atom or its immediate neighbors (non-bonded interactions).

- `bonded1_4: Particles`: A vector of atoms that are connected to the reference atom via three bonds (1-4 interactions).

### Methods

- `new(polymer: Particles, atom: Atom) -> Self`:
  Constructs a new `Neighbor` instance for a given atom and the entire polymer chain. Initializes the neighbor lists.

- `get_neighbours(&mut self)`:
  Initiates the neighbor calculation by determining non-bonded and 1-4 bonded atoms using helper methods.

- `first_bonded(&self) -> Vec<Atom>`:
  Identifies atoms that are directly bonded to the reference atom. Uses connection data (`ALLCONN`) to find bonds based on residue type and atom names.

- `second_bonded(&self) -> (Vec<Atom>, Vec<Atom>)`:
  Identifies second-bonded atoms by finding atoms bonded to the first-bonded atoms. Returns a tuple containing both first and second bonded vectors.

- `third_bonded(&self) -> (Vec<Atom>, Vec<Atom>, Vec<Atom>)`:
  Identifies third-bonded atoms by finding atoms bonded to the second-bonded atoms. Returns a tuple with first, second, and third bonded vectors.

- `fourth_bonded(&self) -> Vec<Atom>`:
  Determines atoms that are beyond third neighbors, effectively considering them as fourth-bonded or non-bonded atoms.

- `formated_fourth(&self) -> Vec<Atom>`:
  Formats the list of fourth-bonded atoms by selecting unique atoms and managing the ordering, avoiding duplicates and ensuring proper indexing.

- `get_bonded1_4(&mut self)`:
  Calculates atoms involved in 1-4 bonded interactions with the reference atom. Considers specific residue and atom combinations to exclude certain interactions based on molecular structure.

- `get_nonbonded(&mut self)`:
  Determines non-bonded atoms relative to the reference atom, which are atoms separated by more than 4 bonds from the reference atom.

## Example

```rust
use libdncs::*;

const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AMBERFB15;

fn main() {
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    for atom in sys.particles.iter() {
        println!("{:?}", atom); // Atom
    }
    for (a, b, c, d) in sys.dihedral.iter() {
        println!("{} {} {} {}", a.name, b.name, c.name, d.name);  // Dihedrals
    }
}
```
