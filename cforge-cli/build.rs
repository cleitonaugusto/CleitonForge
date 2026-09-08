//! Makes the optional q1tsim build runnable.
//!
//! q1tsim ships as a Rust dylib, so linking it leaves the binary loading
//! libq1tsim from `target/<profile>/deps` and libstd from the toolchain, at run
//! time, from directories the loader does not search. Without an rpath the
//! binary dies before `main` with a shared-object error.
//!
//! Only the feature build needs this. The default build links everything
//! statically and gets no extra link arguments.

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    if std::env::var_os("CARGO_FEATURE_Q1TSIM").is_none() {
        return;
    }

    // Sibling `deps/`, relative to wherever the binary ends up.
    println!("cargo:rustc-link-arg=-Wl,-rpath,$ORIGIN/deps");

    // libstd comes from the toolchain, whose location is only known here.
    if let Ok(out) =
        std::process::Command::new(std::env::var("RUSTC").unwrap_or_else(|_| "rustc".into()))
            .arg("--print")
            .arg("target-libdir")
            .output()
    {
        if out.status.success() {
            let dir = String::from_utf8_lossy(&out.stdout).trim().to_string();
            if !dir.is_empty() {
                println!("cargo:rustc-link-arg=-Wl,-rpath,{dir}");
            }
        }
    }
}
