

## Sobol

A struct that generates Sobol sequences for sampling high-dimensional spaces. These sequences are useful in Monte Carlo simulations and numerical integration due to their uniformity properties.

### Fields

- `count: usize` - The current count of Sobol points generated.
- `total: usize` - The total number of Sobol points that can be generated (based on the bit size `SIZE`).
- `current: Vec<usize>` - The current vector of Sobol indices.
- `direction: Vec<Vec<usize>>` - A vector of direction numbers used in Sobol sequence generation.

### Methods

- `new(dimension: usize) -> Self` - Creates a new Sobol instance with the specified `dimension` (should be between 1 and 21201). Initializes the `current` and `direction` fields.

- `get_direction(d: usize) -> Vec<Vec<usize>>` - Generates direction numbers for the Sobol sequence based on the dimension. This method reads from a precomputed direction number file.

### Trait Implementations

- `Iterator for Sobol` - Implements an iterator that produces Sobol sequences as vectors of `f64` numbers. The iterator returns a vector of Sobol points scaled to the range [0, 1], with added uniform noise.

### Example

```rust
use libdncs::*;

fn main() {
    let s = Sobol::new(10);
    for point in s.take(10) {
        println!("{:?}", point);
    }
}
```



## RotateAtDihedral

A struct used to rotate molecular systems around a specified dihedral angle.

### Fields

- `system: System` - The molecular system that is being rotated.
- `rotated: Particles` - The particle system after the rotation has been applied.

### Methods

- `new(system: System) -> Self` - Creates a new `RotateAtDihedral` instance, initializing it with a cloned version of the input `system`.

- `from_pdb(file: &str) -> Vec<f64>` - Reads a molecular system from a PDB file and calculates its dihedral angles. Returns a vector of the calculated dihedral angles.

- `rotate(angle: Vec<f64>)` - Rotates the atoms in the system according to the provided vector of dihedral angles.

- `rotated_energy(angle: Vec<f64>) -> f64` - Rotates the system and calculates the energy of the rotated system using the Amber force field.

- `dihedral_angle(a: &Atom, b: &Atom, c: &Atom, d: &Atom) -> f64` - Calculates the dihedral angle between four atoms in 3D space.

- `energy() -> f64` - Returns the potential energy of the current configuration of the system.

- `to_pdbstring(model: usize, energy: f64) -> String` - Converts the rotated molecular system into a PDB format string, including the energy of the system.


### Example
```rust
use libdncs::*;

const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AMBERFB15;

fn main() {
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    let mut rotate = RotateAtDihedral::new(sys.clone());
    let sobol = Sobol::new(sys.dihedral.len());
    for angle in sobol.skip(50).take(1) {
        println!("{:?}", angle);
        rotate.rotate(angle);
        for atom in rotate.rotated.iter() {
            println!("{:?}", atom);
        }
    }
}
```


## Sampler

A struct that generates samples of molecular conformations by rotating dihedral angles and calculating the resulting energies.

### Fields

- `system: System` - The molecular system being sampled.
- `rotate: RotateAtDihedral` - The `RotateAtDihedral` instance used for rotating the system.
- `angles: Vec<Vec<f64>>` - A vector containing the sampled dihedral angles.
- `sample: Vec<System>` - A vector containing the sampled molecular systems.
- `energy: Vec<f64>` - A vector containing the energies corresponding to each sampled system.

### Methods

- `new(system: System) -> Self` - Creates a new `Sampler` instance initialized with a molecular system.

- `sample(maxsample: usize)` - Generates samples by rotating the dihedral angles using the Sobol sequence. The method calculates the energy for each sampled configuration and stores the angles, energy, and sampled systems.

- `rotatesample(angle: Vec<f64>)` - Reinitializes the system and rotates it according to the provided dihedral angles.

- `conformational_sort()` - Sorts the sampled conformations based on their Boltzmann weights at 300K. The more thermodynamically favorable conformations are ordered first.

- `write_sampled_angles(filename: &str)` - Writes the sampled dihedral angles and corresponding energies to a specified file.

- `to_pdb(filename: &str)` - Writes the sampled molecular systems to a PDB file, including the energy of each conformation.

- `to_pdbfiles(foldername: &str)` - Writes each sampled conformation to a separate PDB file in the specified folder.


## Example
```rust
use libdncs::*;

// Configuration
const NAME: &str = "DNCS";
const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AMBERFB15;
const NO_OF_SAMPLE: usize = 10;

fn main() {
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    let mut sample = Sampler::new(sys);
    sample.sample(NO_OF_SAMPLE);
    sample.write_angles(&format!("{}.out", NAME));
    sample.to_pdb(&format!("{}.pdb", NAME));
}
```
