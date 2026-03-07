# Rust Superkmers Project Guide

## Build & Test Commands
- Build: `cargo build`
- Test all: `cargo test`
- Run single test: `cargo test main`
- Benchmark: `./bench.sh` or `RUSTFLAGS='-C target-cpu=native' cargo +nightly bench`
- Format code: `cargo fmt`
- Check code: `cargo clippy`

## Code Style Guidelines
- **Naming**: Use snake_case for variables/functions, CamelCase for types/structs
- **Imports**: Group by external crates first, then internal modules
- **Error handling**: Use `expect()` with descriptive messages for file ops, `unwrap()` for trusted operations
- **Types**: Include explicit type annotations, especially for public interfaces
- **Comments**: Use doc comments (`///`) for public API, regular comments (`//`) for implementation details
- **Module structure**: Separate modules by implementation approach, main library in lib.rs
- **Function design**: Prefer small, focused functions with clear inputs/outputs
- **Performance**: For this bioinformatics project, performance is critical - benchmark different approaches

Follow standard Rust idioms and patterns from the official Rust style guide.