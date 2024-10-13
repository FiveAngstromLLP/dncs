#![allow(unused_variables)]
use liblbfgs_sys::*;
use std::os::raw::{c_int, c_void};

// Define our function to minimize: f(x) = x^2 + 1
extern "C" fn evaluate(
    instance: *mut c_void,
    x: *const f64,
    g: *mut f64,
    n: c_int,
    step: f64,
) -> f64 {
    unsafe {
        let x = std::slice::from_raw_parts(x, n as usize);
        let g = std::slice::from_raw_parts_mut(g, n as usize);

        // Calculate f(x) = x^2 + 1
        let fx = 1.0 / x[0] - 100.0;

        // Calculate gradient: f'(x) = 2x
        g[0] = 2.0 * x[0];

        fx
    }
}

fn main() {
    let n = 1; // Dimension of the problem
    let mut x = vec![2.0]; // Initial guess
    let mut f = 0.0; // Will store the minimum function value

    let mut param = lbfgs_parameter_t {
        m: 0,
        epsilon: 0.0,
        past: 0,
        delta: 0.0,
        max_iterations: 0,
        linesearch: 0,
        max_linesearch: 0,
        min_step: 0.0,
        max_step: 0.0,
        ftol: 0.0,
        wolfe: 0.0,
        gtol: 0.0,
        xtol: 0.0,
        orthantwise_c: 0.0,
        orthantwise_start: 0,
        orthantwise_end: 0,
    };
    unsafe {
        lbfgs_parameter_init(&mut param);
    }

    let result = unsafe {
        lbfgs(
            n,
            x.as_mut_ptr(),
            &mut f,
            Some(evaluate),
            None,
            std::ptr::null_mut(),
            &mut param as *mut lbfgs_parameter_t,
        )
    };

    println!("LBFGS optimization result: {:?}", result);
    println!("Optimal x: {:?}", x);
    println!("Minimum function value: {}", f);
}
