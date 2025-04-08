use libdncs::*;
use std::sync::Arc;
const FORCE_FIELD: FF = FF::AmberFB15;

fn main() {
    let mut sys = System::new("AAAAAAAAAA", FORCE_FIELD.init());
    sys.init_parameters();

    let amber = Amber::new(Arc::new(sys.clone()));
    assert!(&sys.particles.len() == &amber.system.particles.len());
    let eng = amber.energy();
    println!("Energy: {} kJ/mol", eng);
}
