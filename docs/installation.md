## Quick Setup
1. **Install Rust**: If `cargo -h` fails, install via [rustup.rs](https://rustup.rs).
2. **Clone Repository**: `git clone https://github.com/IT-Academic-Research-Services/seqtoid-pipelines.git`
3. **Navigate**: `cd seqtoid-pipelines`
4. **Build**:
- Standard: `cargo build --release`
- Optimized for your CPU: `RUSTFLAGS="-C target-cpu=native" cargo build --release`
- Optimized for the UCSF cluster nodes: `RUSTFLAGS="-C target-cpu=znver3" cargo build --release`
5. **Make Executable**: `chmod u+x target/release/seqtoid-pipelines`
6. **Install Globally** (optional): `sudo mv target/release/seqtoid-pipelines /usr/local/bin`
7. **Test**: From another directory, run `seqtoid-pipelines --help` to verify.
