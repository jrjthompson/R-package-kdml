\name{cuda_info}
\alias{cuda_info}
\title{Inspect CUDA Support for the KDML Sampler}
\description{
Return build, runtime, and device information for the optional CUDA backend.
}
\usage{
cuda_info()
}
\details{
CUDA is optional. If it is unavailable, the package still works normally using
the CPU backend. The \code{available} value is \code{TRUE} only when CUDA
support was included during installation and the current R session can access
a compatible CUDA device and runtime.
}
\value{
A named list with components:
\item{compiled}{logical; whether this installation contains CUDA support.}
\item{available}{logical; whether the CUDA backend is currently usable.}
\item{compile_version}{integer CUDA toolkit version used to build the
sidecar, or zero when CUDA was not compiled.}
\item{compile_version_string}{character representation of the compile-time
CUDA version or its availability status.}
\item{runtime_version}{integer CUDA runtime version, or \code{NA_integer_}
when it cannot be queried.}
\item{device_count}{integer number of devices that could be queried.}
\item{devices}{a list with one element per device. Each available element
contains the device index, name, global memory in bytes, compute-capability
major and minor versions, multiprocessor count, driver version, and runtime
version. A device element can be \code{NULL} if its properties could not be
queried.}
}
\seealso{
\code{\link{mcmc.dkps}}, \code{\link{mcmc.dkss}}
}
\examples{
info <- cuda_info()
info[c("compiled", "available", "device_count")]
}
\author{
Daniel Krasnov \email{daniel.krasnov@mail.mcgill.ca}
}
