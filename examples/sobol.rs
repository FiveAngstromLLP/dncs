use libdncs::Sobol;
use rayon::iter::{ParallelBridge, ParallelIterator};

fn main() {
    let s = Sobol::new(1);
    s.take(100)
        .into_iter()
        .par_bridge()
        .for_each(|i| println!("{:?}", i));
}
