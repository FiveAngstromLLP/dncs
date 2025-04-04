extern crate libdncs;

use libdncs::{System, FF};

fn main() {
    let mut system = System::new("AA", FF::Amber99SB.init());
    system.get_neighbours();

    use std::fs;
    use std::io::Write;
    use std::path::Path;

    // Delete the neighbours directory if it exists, then create a new one
    let neighbours_dir = Path::new("examples/neighbours");
    if neighbours_dir.exists() {
        fs::remove_dir_all(neighbours_dir).expect("Failed to delete neighbours directory");
    }
    fs::create_dir_all(neighbours_dir).expect("Failed to create neighbours directory");

    // For each bond distance (0 to 5)
    for bond_index in 1..6 {
        // Create or overwrite the file
        let file_path = neighbours_dir.join(format!("bond_{}.log", bond_index));
        let mut file = fs::File::create(file_path).expect("Failed to create file");

        // Write all atom pairs with this bond distance
        for (i, particle_i) in system.particles.iter().enumerate() {
            write!(file, "{}:{} \t>> ", i + 1, particle_i.name).expect("Failed to write to file");

            for (j, particle_j) in system.particles.iter().enumerate() {
                // Check if particles are bonded at the specified bond distance
                let is_bonded = match bond_index {
                    1 => system.firstbonded.get(i).map_or(false, |neighbors| {
                        neighbors.iter().any(|n| n.serial == particle_j.serial)
                    }),
                    2 => system.secondbonded.get(i).map_or(false, |neighbors| {
                        neighbors.iter().any(|n| n.serial == particle_j.serial)
                    }),
                    3 => system.thirdbonded.get(i).map_or(false, |neighbors| {
                        neighbors.iter().any(|n| n.serial == particle_j.serial)
                    }),
                    4 => system.bonded1_4.get(i).map_or(false, |neighbors| {
                        neighbors.iter().any(|n| n.serial == particle_j.serial)
                    }),
                    5 => system.nonbonded.get(i).map_or(false, |neighbors| {
                        neighbors.iter().any(|n| n.serial == particle_j.serial)
                    }),
                    _ => false,
                };

                if is_bonded {
                    write!(file, "{}:{}; ", j + 1, particle_j.name)
                        .expect("Failed to write to file");
                }
            }

            writeln!(file).expect("Failed to write to file");
        }
    }
}
