use std::sync::Arc;

use libdncs::*;

// Configuration

// const SEQUENCE: &str = "AA";
const FORCE_FIELD: FF = FF::AmberFB15;

fn main() {
    // let mut sys = System::from_pdb("Result/sample/sample_0000.pdb", FORCE_FIELD.init());
    let mut sys = System::new("AA", FORCE_FIELD.init());
    // for i in sys.particles.iter() {
    //     println!("{:?}", i);
    // }

    // System
    // let mut sys = System::new(SEQUENCE, FORCE_FIELD.init());
    sys.to_pdb("experiment.pdb");
    sys.init_parameters();
    let amber = Amber::new(Arc::new(sys));
    let eng = amber.energy();
    println!("Energy: {} kJ/mol", eng);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_name() {
        let sys = System::new("BA", FORCE_FIELD.init());
        println!("{}", sys.particles);
    }
}
