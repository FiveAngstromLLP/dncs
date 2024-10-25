## Minimizer

A structure for performing energy minimization and conformational analysis of molecular systems.

### Fields
- `sample: Sampler` - The sampler object containing the initial molecular system and sampling parameters.

- `minimized: Vec<System>` - A vector containing the minimized molecular systems for each sampled conformation.

- `energy: Vec<f64>` - The potential energies (in `KCal/mol`) corresponding to each minimized conformation.

- `angles: Vec<Vec<f64>>` - The optimized dihedral angles (in `degrees`) for each minimized conformation.

### Methods

- `new` - Creates a new Minimizer instance initialized with the provided Sampler (sample: Sampler) and returns a new Minimizer instance (Self). The vectors for minimized systems, energies, and angles are pre-allocated based on the number of samples in the Sampler.

- `diff` - Calculates the energy of a given dihedral angle configuration using a provided function (`Arc<dyn Fn(Vec<f64>) -> f64>`) that calculates energy from dihedral angles (wrapped in an Arc for thread safety) and a vector of dihedral angles to evaluate (`Vec<f64>`). Returns a tuple containing the energy value at the given angles (`f64`) and a vector of energy gradients with respect to each angle (`Vec<f64>`). The function performs numerical differentiation of the energy function with respect to the angles and is executed asynchronously using tokio tasks.

- `minimize` - Performs energy minimization on all sampled conformations using the L-BFGS algorithm. Optimizes the dihedral angles to find local energy minima. Updates the minimized systems, energies, and angles vectors with the optimization results. The minimization is performed asynchronously using tokio tasks. Returns no value (async).

- `conformational_sort` - Sorts the conformations based on their Boltzmann weights at `300K`. Reorders the minimized systems, energies, and angles vectors according to the conformational weights in descending order. Higher weights indicate more thermodynamically favorable conformations. Returns no value (async).

- `write_angles` - Writes the optimized dihedral angles and corresponding energies to a specified file (`&str`). Each line contains the conformation index, energy, and comma-separated angles. Returns an async `Result<()>` indicating success or failure of the file operation.

- `to_pdb` - Writes all conformations to a PDB format file specified by filename (`&str`). Includes both the initial system and all minimized conformations. Each MODEL entry contains atomic coordinates and the corresponding energy. Returns an async `Result<()>` indicating success or failure of the file operation.


### Example
```rust
use libdncs::*;

// Configuration
const NAME: &str = "DNCS";
const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AMBERFB15;
const NO_OF_SAMPLE: usize = 10;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // System
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    // Sample
    let mut sample = Sampler::new(sys);
    println!("Generating Samples..");
    sample.sample(NO_OF_SAMPLE);
    // Minimizer
    let mut minimizer = Minimizer::new(sample);
    println!("Executing Minimizer..");
    minimizer.minimize().await;
    minimizer.conformational_sort().await;
    minimizer.write_angles(&format!("{}.out", NAME)).await?;
    minimizer.to_pdb("minimized.pdb").await?;
    println!("Completed!");
    Ok(())
}
```
