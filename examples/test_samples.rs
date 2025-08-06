use std::sync::Arc;
use libdncs::*;

fn main() {
    println!("Testing energy calculation for specific samples...");
    
    // Test Sample 854 (the one with clashes)
    println!("\n=== Sample 854 (with clashes) ===");
    let mut sys854 = System::from_pdb("Result/6RRO/sample_0854.pdb", FF::AmberFB15.init());
    sys854.init_parameters();
    let energy854 = Amber::new(Arc::new(sys854)).energy();
    println!("Total energy: {} kJ/mol", energy854);
    
    // Test Sample 706 (better geometry)
    println!("\n=== Sample 706 (better geometry) ===");
    let mut sys706 = System::from_pdb("Result/6RRO/sample_0706.pdb", FF::AmberFB15.init());
    sys706.init_parameters();
    let energy706 = Amber::new(Arc::new(sys706)).energy();
    println!("Total energy: {} kJ/mol", energy706);
    
    println!("\n=== Comparison ===");
    println!("Energy difference: {} kJ/mol", energy854 - energy706);
    if energy854 < energy706 {
        println!("Sample 854 has LOWER energy (incorrect - it has clashes!)");
    } else {
        println!("Sample 706 has lower energy (correct)");
    }
}