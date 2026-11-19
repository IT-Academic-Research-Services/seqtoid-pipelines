use std::{
    env,
    error::Error,
    ffi::OsString,
    fs::{self, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
};

fn sidecar(prefix: &Path, suffix: &str) -> PathBuf {
    let mut name: OsString = prefix.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        return Err(
            "usage: write_prefilter OUTPUT_PREFIX QUERY_KEY TARGET_KEY PREF_SCORE".into()
        );
    }

    let prefix = Path::new(&args[1]);
    let query_key: u64 = args[2].parse()?;
    let target_key: u64 = args[3].parse()?;
    let score: i32 = args[4].parse()?;

    let index_path = sidecar(prefix, ".index");
    let type_path = sidecar(prefix, ".dbtype");

    // Refuse to overwrite an existing experiment.
    for path in [prefix, index_path.as_path(), type_path.as_path()] {
        if path.exists() {
            return Err(format!("output already exists: {}", path.display()).into());
        }
    }

    // One query entry, containing one candidate, terminated by one NUL.
    let mut data = format!("{target_key}\t{score}\t0\n").into_bytes();
    data.push(0);

    let mut output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(prefix)?;
    output.write_all(&data)?;
    output.flush()?;

    // Indexed length includes the terminating NUL.
    fs::write(&index_path, format!("{query_key}\t0\t{}\n", data.len()))?;

    // DBTYPE_PREFILTER_RES = 7; uncompressed, little-endian.
    fs::write(&type_path, 7_u32.to_le_bytes())?;

    println!(
        "Wrote {}: query={}, target={}, pref_score={}, diagonal=0, entry_bytes={}",
        prefix.display(), query_key, target_key, score, data.len()
    );
    Ok(())
}
