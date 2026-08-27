# kdml 1.2.0

- Added CPU and optional CUDA MCMC backends for kernel
  and bandwidth selection with DKPS and DKSS.
- Added posterior summaries, convergence diagnostics, plots, and distance construction from retained MCMC states.
- Fixed some minor typos and bugs such as mixed-type column indexing, one-column categorical handling, nominal
  kernel naming, and invalid-bandwidth penalties.
- Added test coverage for mixed-type distances, MSCV bandwidth selection, diagnostics, and CUDA and CPU sampling.
- Reworked the vignette to show an example analysis and new package features.
