#![allow(dead_code)]
use rayon::prelude::*;
use std::sync::LazyLock;

const DIRECTION: LazyLock<Vec<String>> = LazyLock::new(|| {
    include_str!("../data/new-joe-kuo-6.21201")
        .lines()
        .map(String::from)
        .collect()
});

const SIZE: usize = 32;

pub struct Sobol {
    count: usize,
    total: usize,
    current: Vec<usize>,
    direction: Vec<Vec<usize>>,
}

impl Sobol {
    pub fn new(dimension: usize) -> Self {
        assert!(
            dimension >= 1 && dimension <= 21201,
            "DIMESNSION must in range (1-21201)"
        );
        let mut sobol = Self {
            count: 0,
            total: 1 << SIZE,
            current: vec![0; dimension],
            direction: vec![vec![0; SIZE]; dimension],
        };
        sobol.direction = sobol.get_direction(dimension);
        sobol
    }
    fn get_direction(&self, d: usize) -> Vec<Vec<usize>> {
        (1..=d)
            .into_par_iter()
            .map(|d| match d {
                1 => (1..=SIZE).map(|i| 1 << (SIZE - i)).collect(),
                _ => {
                    let mut val = vec![0; SIZE];
                    let direction: Vec<usize> = DIRECTION[d]
                        .split_whitespace()
                        .skip(2)
                        .map(|f| f.parse().unwrap())
                        .collect();
                    for i in 0..SIZE {
                        let s = direction.len() - 1;
                        if i < s {
                            for j in 0..s {
                                val[j] = direction[j + 1] << (SIZE - j - 1)
                            }
                        } else {
                            for k in s..SIZE {
                                val[k] = val[k - s] ^ (val[k - s] >> s);
                                for l in 1..s {
                                    let a = (direction[0] >> (s - l)) & 1;
                                    val[k] ^= a * val[k - l];
                                }
                            }
                        }
                    }
                    val
                }
            })
            .collect()
    }
}

impl Iterator for Sobol {
    type Item = Vec<f64>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.count < self.total {
            let rmz = (self.count ^ (self.total - 1)).trailing_zeros() as usize;
            let val: Vec<usize> = self.current.clone();
            self.count += 1;
            self.current = val
                .par_iter()
                .enumerate()
                .map(|(d, m)| self.direction[d][rmz] ^ m)
                .collect();
            let num: Vec<f64> = val
                .par_iter()
                .map(|i| *i as f64 / f64::powi(2.0, SIZE as i32))
                .collect();
            // let num = bakers_transform(num);
            Some(num)
        } else {
            None
        }
    }
}

#[inline]
fn bakers_transform(input: Vec<f64>) -> Vec<f64> {
    let m = input.len();
    if m < 2 {
        input
    } else {
        let k = ((m as f64 * input[0]).floor() as usize + 1).min(m);
        let first = m as f64 * input[0] - (k - 1) as f64;
        let transformed: Vec<f64> = std::iter::once(first)
            .chain(
                input[1..]
                    .iter()
                    .map(|&xi| (xi + (k - 1) as f64) / m as f64),
            )
            .collect();
        transformed
    }
}

mod test {
    #[allow(unused)]
    use super::*;
    #[test]
    fn sobol() {
        let s = Sobol::new(10);
        for i in s.take(10) {
            println!("{:?}", i)
        }
    }
}
