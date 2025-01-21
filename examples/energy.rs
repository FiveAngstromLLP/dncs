use std::sync::Arc;

use libdncs::*;

// Configuration
const SEQUENCE: &str = "YGGFM";
const FORCE_FIELD: FF = FF::AMBERFB15;

fn main() {
    let sys = System::new("BB", FORCE_FIELD.init());
    for i in sys.particles.iter() {
        println!("{:?}", i);
    }

    // System
    let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.init_parameters();
    let amber = Amber::new(Arc::new(sys));
    let eng = amber.energy();
    println!("Energy: {} KCal/Mol", eng);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_name() {
        let sys = System::new("BB", FORCE_FIELD.init());
        println!("{}", sys.particles);
    }
}
