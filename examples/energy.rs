use std::sync::Arc;

use libdncs::*;

// Configuration
const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AMBERFB15;

fn main() {
    // System
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    let amber = Amber::new(Arc::new(sys));
    let eng = amber.energy();
    println!("Energy: {} KCal/Mol", eng);
}
