use libdncs::*;
use std::str::FromStr;
use std::sync::Arc;

#[no_mangle]
pub extern "C" fn sample(
    seq: *const std::os::raw::c_char,
    folder: *const std::os::raw::c_char,
    ff: *const std::os::raw::c_char,
    sample: usize,
    method: Method,
) {
    let seq = unsafe { std::ffi::CStr::from_ptr(seq).to_str().unwrap() };
    let folder = unsafe { std::ffi::CStr::from_ptr(folder).to_str().unwrap() };
    let ff = unsafe { std::ffi::CStr::from_ptr(ff).to_str().unwrap() };

    let system = System::new(seq, FF::from_str(ff).unwrap().init());
    let mut sampler = Sampler::new(Arc::new(system), method, folder.to_string());
    sampler.sample(sample);
}

#[no_mangle]
pub extern "C" fn pdb_to_angle(
    filename: *const std::os::raw::c_char,
) -> *const std::os::raw::c_char {
    let filename = unsafe { std::ffi::CStr::from_ptr(filename).to_str().unwrap() };
    let angle = RotateAtDihedral::from_pdb(filename);
    let angle_csv = angle
        .iter()
        .map(|&a| a.to_string())
        .collect::<Vec<String>>()
        .join(", ");
    let c_str = std::ffi::CString::new(angle_csv).unwrap();
    c_str.into_raw()
}
