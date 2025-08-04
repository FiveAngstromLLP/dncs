use libdncs::Sobol;

fn main() {
    let s = Sobol::new(10, libdncs::Method::Fold);
    for i in 0..100 {
        println!("{:?}", s.get_index(i));
    }
}
