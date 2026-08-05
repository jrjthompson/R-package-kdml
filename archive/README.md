# KDML MCMC prototype archive

`KDML-MCMC-Sampler-reference.zip` preserves the small source and specification
set from the standalone prototype that was used while porting the sampler into
the KDML R package.

Included:

- Rust and CUDA source, Cargo metadata, and the standalone CLI;
- the mathematical model specification in TeX and PDF form;
- the original kernel-parity test;
- the standalone R analysis entry points; and
- the original compressed mixed-type datasets.

Excluded:

- generated experiment results under `local/`;
- Cargo build artifacts under `rust/target/`; and
- report caches and rendered experimental outputs.

The package implementation is maintained independently in `R/`, `src/`, and
`tests/`; it does not read this archive at build time or runtime.

Archive SHA-256:
`1B78282F125E4308F217CE1C6B4538C5A44DC874A4E2003042E4ADC3C6D4AC97`.
