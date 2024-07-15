# Tester for Polynomial Smoothers

Test for the Polynomial Smoothers in AMG4PSBLAS

To clone the repository with the needed submodules do:
```bash
git clone --recurse-submodules git@github.com:Cirdans-Home/polynomialsmoothers.git
```
this will populate also the `psblas` and  `amg4psblas` folders with the relevant code.

## Contributors
- Pasqua D'Ambra,
- Fabio Durastante,
- Salvatore Filippone,
- Stefano Massei,
- Stephen Thomas

## Examples

The folder `3dlaplacian` contains the Finite Difference 3D Poisson Benchmark, while the folder `anisotropy` contains the FEM 2D diffusion with rotated anisotropy. To use the code you have to compile and install PSBLAS and AMG4PSBLAS, then the code in these two folders can be compiled and linked against it. Example of scripts for running the test in HPC environment is contained in the logfiles folder.
