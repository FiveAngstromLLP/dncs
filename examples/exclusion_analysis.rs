use std::sync::Arc;
use libdncs::*;

fn main() {
    println!("Analyzing non-bonded exclusions...");
    
    let mut sys = System::from_pdb("Result/6RRO/sample_0854.pdb", FF::AmberFB15.init());
    sys.init_parameters();
    
    println!("Total atoms: {}", sys.particles.len());
    
    // Check the first few atoms' nonbonded lists
    println!("\nNon-bonded interaction ranges for first 10 atoms:");
    for i in 0..10.min(sys.particles.len()) {
        let atom = &sys.particles[i];
        let nb_list = &sys.nonbonded[atom.serial - 1];
        print!("Atom {}: ", atom.serial);
        
        if nb_list.is_empty() {
            println!("no non-bonded interactions");
        } else {
            for chunk in nb_list.chunks(2) {
                if chunk.len() == 2 {
                    print!("{}..{} ", chunk[0].serial, chunk[1].serial);
                }
            }
            println!();
        }
    }
    
    // Let's see what pairs are actually calculated vs. which clashes exist
    println!("\nLet's see which close pairs are actually included in non-bonded calculations:");
    
    // First collect all the short distance pairs again
    let mut short_pairs = Vec::new();
    for i in 0..sys.particles.len() {
        for j in (i+1)..sys.particles.len() {
            let atom_i = &sys.particles[i];
            let atom_j = &sys.particles[j];
            
            let dx = atom_j.position[0] - atom_i.position[0];
            let dy = atom_j.position[1] - atom_i.position[1];
            let dz = atom_j.position[2] - atom_i.position[2];
            let distance = (dx*dx + dy*dy + dz*dz).sqrt() * 0.1;
            
            if distance < 0.15 {
                short_pairs.push((i, j, distance));
            }
        }
    }
    
    println!("Found {} short distance pairs", short_pairs.len());
    
    // Check which of these are actually calculated in non-bonded
    let mut calculated_pairs = 0;
    let mut excluded_pairs = 0;
    
    for (i, j, dist) in &short_pairs {
        let atom_i = &sys.particles[*i];
        let atom_j = &sys.particles[*j];
        
        // Check if this pair would be calculated
        let mut is_calculated = false;
        
        // Check in atom_i's nonbonded list
        for chunk in sys.nonbonded[atom_i.serial - 1].chunks(2) {
            if chunk.len() == 2 {
                let start = chunk[0].serial;
                let end = chunk[1].serial;
                if atom_j.serial >= start && atom_j.serial <= end {
                    is_calculated = true;
                    break;
                }
            }
        }
        
        if is_calculated {
            calculated_pairs += 1;
            println!("CALCULATED: Atoms {}({}) - {}({}) dist={:.4}nm", 
                     atom_i.serial, atom_i.name, atom_j.serial, atom_j.name, dist);
        } else {
            excluded_pairs += 1;
        }
    }
    
    println!("\nSummary:");
    println!("Total short pairs: {}", short_pairs.len());
    println!("Excluded from non-bonded: {}", excluded_pairs);
    println!("Calculated in non-bonded: {}", calculated_pairs);
}