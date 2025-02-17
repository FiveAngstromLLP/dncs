use std::sync::Arc;

use libdncs::*;

const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AmberFB15;

fn main() {
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    let mut rotate = RotateAtDihedral::new(Arc::new(sys.clone()));
    let sobol = Sobol::new(sys.dihedral.len());
    for angle in sobol.skip(50).take(1) {
        println!("{:?}", angle);
        rotate.rotate(angle);
        for atom in rotate.rotated.iter() {
            println!("{:?}", atom);
        }
    }
}
