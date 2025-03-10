use libdncs::Sobol;
use rayon::iter::{ParallelBridge, ParallelIterator};

fn main() {
    let s = Sobol::new(10);
    s.take(10)
        .into_iter()
        .par_bridge()
        .for_each(|i| println!("{:?}", i));
}
