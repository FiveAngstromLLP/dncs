use libdncs::*;
use std::sync::Arc;
const FORCE_FIELD: FF = FF::AmberFB15;

fn main() {
    let mut sys = System::from_pdb("experiment.pdb", FORCE_FIELD.init());
    // let sys = System::new("YGGFM", FORCE_FIELD.init());
    // sys.to_pdb("experiment.pdb");

    // let pdb: Vec<Atom> = std::fs::read_to_string("experiment.pdb")
    //     .unwrap()
    //     .lines()
    //     .filter(|x| x.starts_with("ATOM"))
    //     .map(String::from)
    //     .map(Atom::new)
    //     .collect();

    // let
    sys.init_parameters();

    // let amber = Amber::new(Arc::new(sys));
    // let eng = amber.energy();
    // println!("Energy: {} kJ/mol", eng);
    //
    // let mut sample =
}
