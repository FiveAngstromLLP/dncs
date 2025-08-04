// use std::sync::Arc;

use libdncs::*;

// const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AmberFB15;

fn main() {
    let start_time = std::time::Instant::now();

    let mut sys = System::new("AA", FORCE_FIELD.init());
    sys.to_pdb("output.pdb");

    let init_start_time = std::time::Instant::now();
    sys.init_parameters();
    let init_elapsed = init_start_time.elapsed();

    println!("Parameter initialization took: {:?}", init_elapsed);

    // let mut rotate = RotateAtDihedral::new(Arc::new(sys.clone()));
    // let sobol = Sobol::new(sys.dihedral.len());
    // for angle in sobol.skip(50).take(1) {
    //     println!("{:?}", angle);
    //     rotate.rotate(angle);
    //     for atom in rotate.rotated.iter() {
    //         println!("{:?}", atom);
    //     }
    // }

    let total_elapsed = start_time.elapsed();
    println!("Total execution time: {:?}", total_elapsed);
}
