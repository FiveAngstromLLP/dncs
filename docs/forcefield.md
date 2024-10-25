## Amber

The `Amber` struct is used to calculate molecular system energies using Amber force field parameters. This Structure performs calculations for molecular energy including non-bonded interactions (Lennard-Jones and electrostatic forces), bond stretching energies using harmonic potentials, angle bending energies between three bonded atoms, torsional rotations around bonds, and hydrogen bonding interactions. Below are its properties and methods, with mathematical representations included for clarity.

### Fields
- `system`: An instance of the `System` struct that encapsulates information about the particles in the system and their interaction parameters.

### Methods
- `new(system: System) -> Self` creates a new `Amber` instance from a `System` containing particles and their interaction properties, returning a new `Amber` instance.

- `energy(&self) -> f64` computes and returns the total system energy in kcal/mol: \\[ E_{\text{total}} = E_{\text{non-bonded}} + E_{\text{bonds}} + E_{\text{angles}} + E_{\text{torsions}} + E_{\text{hydrogen bonds}} \\]

- `nonbonded_energy(&self, iatom: &Atom) -> f64` calculates the non-bonded energy (in kg·Å²/s²) for an atom including Lennard-Jones and electrostatic interactions: \\[ E_{\text{non-bonded}} = \sum_{j \neq i} \left( E_{\text{LJ}}(i, j) + E_{\text{electrostatic}}(i, j) \right) \\]

- `harmonic_bond_force(&self, iatom: &Atom) -> f64` calculates the harmonic bond energy (in kcal/mol) for an atom: \\[ E_{\text{bond}} = \frac{1}{2} k (r - r_0)^2 \\] where \\( k \\) is the bond force constant, \\( r \\) is current bond length, \\( r_0 \\) is equilibrium length.

- `harmonic_angle_force(&self, iatom: &Atom) -> f64` computes the harmonic angle energy (in kcal/mol) for an atom: \\[ E_{\text{angle}} = \frac{1}{2} k_{\theta} (\theta - \theta_0)^2 \\] where \\( k_{i\theta} \\) is angle force constant, \\( \theta \\) is current angle, \\( \theta_0 \\) is equilibrium angle.

- `periodic_torsional_force(&self) -> f64` returns the total periodic torsional energy in kcal/mol: \\[ E_{\text{torsion}} = \sum_{n} \frac{1}{2} k_n (1 + \cos(n \phi - \gamma)) \\] where \\( V_n \\) is torsional force constant, \\( n \\) is periodicity, \\( \phi \\) is torsion angle, \\( \gamma \\) is phase.

- `hydrogen_bond_energy(&self) -> f64` returns the total hydrogen bond energy in kcal/mol.

- `lennard_jones_energy(i: &Atom, j: &Atom) -> f64` calculates Lennard-Jones energy (in kg·Å²/s²) between atoms: \\[ E_{\text{LJ}}(i, j) = 4\epsilon \left( \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^{6} \right) \\] where \\( r \\) is atom distance, \\( \epsilon \\) is potential well depth, \\( \sigma \\) is zero-potential distance.

- `electrostatic_energy(i: &Atom, j: &Atom) -> f64` calculates electrostatic energy (in kg·Å²/s²) between atoms: \\[ E_{\text{electrostatic}}(i, j) = \frac{k_e \cdot q_i \cdot q_j}{r} \\] where \\( k_e \\) is Coulomb's constant, \\( q_i \\), \\( q_j \\) are atom charges, \\( r \\) is distance.

- `angle(i: &Atom, j: &Atom, k: &Atom) -> f64` calculates angle in degrees between three atoms: \\[ \theta = \arccos \left( \frac{\vec{b} \cdot \vec{c}}{|\vec{b}| |\vec{c}|} \right) \\] where \\( \vec{b} \\), \\( \vec{c} \\) are vectors from atom positions.

- `distance(i: &Atom, j: &Atom) -> f64` returns Euclidean distance in angstroms between atoms: \\[ d = \sqrt{(x_j - x_i)^2 + (y_j - y_i)^2 + (z_j - z_i)^2} \\]

## Example

```rust
use libdncs::*;

// Configuration
const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AMBERFB15;

fn main() {
    // System
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    let amber = Amber::new(sys);
    let eng = amber.energy();
    println!("Energy: {} KCal/Mol", eng);
}
```
