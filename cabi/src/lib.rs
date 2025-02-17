use libdncs::*;
use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::str::FromStr;
use std::sync::Arc;

#[no_mangle]
pub extern "C" fn sample(
    seq: *const c_char,
    folder: *const c_char,
    ff: *const c_char,
    sample: usize,
    method: *const c_char,
    temp: f64,
) -> i32 {
    let seq_str = unsafe {
        if seq.is_null() {
            eprintln!("Error: Sequence pointer is null");
            return 1;
        }
        CStr::from_ptr(seq).to_str().unwrap()
    };

    let folder_str = unsafe {
        if folder.is_null() {
            eprintln!("Error: Folder pointer is null");
            return 1;
        }
        CStr::from_ptr(folder).to_str().unwrap()
    };

    let ff_str = unsafe {
        if ff.is_null() {
            eprintln!("Error: FF pointer is null");
            return 1;
        }
        CStr::from_ptr(ff).to_str().unwrap()
    };

    let method_str = unsafe {
        if method.is_null() {
            eprintln!("Error: Method pointer is null");
            return 1;
        }
        CStr::from_ptr(method).to_str().unwrap()
    };

    let method = match Method::from_str(method_str) {
        Ok(m) => m,
        Err(_) => {
            eprintln!("Error: Invalid method name: {}", method_str);
            return 2;
        }
    };

    let mut system = System::new(seq_str, FF::from_str(ff_str).unwrap().init());
    system.init_parameters();
    let mut sampler = Sampler::new(Arc::new(system), method, folder_str.to_string());
    sampler.sample(sample, temp);

    0 // Success
}

#[no_mangle]
pub extern "C" fn pdb_to_angle(filename: *const c_char) -> *mut c_char {
    let filename_str = unsafe {
        if filename.is_null() {
            eprintln!("Error: Filename pointer is null");
            return std::ptr::null_mut(); // Return null pointer on error
        }
        CStr::from_ptr(filename).to_str().unwrap()
    };

    let angle = RotateAtDihedral::from_pdb(filename_str);
    let angle_csv = angle
        .iter()
        .map(|&a| a.to_string())
        .collect::<Vec<String>>()
        .join(", ");

    match CString::new(angle_csv) {
        Ok(c_str) => c_str.into_raw(),
        Err(e) => {
            eprintln!("Error creating CString: {}", e);
            std::ptr::null_mut()
        }
    }
}

#[no_mangle]
pub extern "C" fn free_string(s: *mut c_char) {
    unsafe {
        if s.is_null() {
            return;
        }
        let _ = CString::from_raw(s); // Drop CString to free memory
    }
}
