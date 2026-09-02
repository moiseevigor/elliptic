# Strict device stub for the GPU code paths

`gpuArray.m` is a stand-in for an `ocl` (OpenCL) device array that behaves
like the real one where it matters for correctness of the kernels:

- elementwise operators work between two device arrays or with a plain
  scalar; mixing a device array with a host **matrix** errors exactly like
  `ocl` does (`binary operator not implemented for 'ocl matrix' by 'matrix'`),
- logical indexing and indexed assignment with logical/device masks error
  (`ocl` has no logical indexing),
- `isreal` is false, `gather` rejects host arrays.

`testGpuStrict.m` puts this directory in front of the path, switches
`elliptic_config('gpu', true)` and compares every GPU path with the CPU
path.  The identity stubs used before (plain `gpuArray = @(x) x`) let three
host/device mixing defects through to the L4 hardware runs; this stub
catches them on a laptop.
