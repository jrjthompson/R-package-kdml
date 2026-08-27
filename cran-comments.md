
## kdml -- Update

This update adds a compiled (CUDA and CPU) Markov chain Monte Carlo (MCMC) algorithm for
selecting feature-specific kernels and bandwidths for the product- and
summation-similarity distances. It also adds supporting posterior
diagnostics and distance reconstruction tools. In addition, it corrects bugs
in the `dkps()`, `kss()`, and `dkss()` functions affecting mixed-type data,
single-column inputs, nominal category counts, ordinal variables, and maximum
similarity cross-validation.

## Test environments

- Local: macOS Tahoe 26.5.1, R 4.5.2 (aarch64-apple-darwin20)
- win-builder R-release: Windows Server 2022 x64 (build 20348), R 4.6.1
  (2026-06-24 ucrt, x86_64-w64-mingw32)
- win-builder R-oldrelease: Windows Server 2022 x64 (build 20348), R 4.5.3
  (2026-03-11 ucrt, x86_64-w64-mingw32)
- win-builder R-devel: Windows Server 2022 x64 (build 20348), R Under development
  (unstable) (2026-08-17 r90424 ucrt, x86_64-w64-mingw32)

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
