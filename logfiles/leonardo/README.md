# README

The `Makefile` makes hard-link to the executables in the folder `3dlaplacian` and `anisotropy`. You have
to adjust it to your configuration.

The SLURM batch script to launch the experiments can be obtained by doing:
```bash
./gen_script.sh <script to launch> <number of polynomial iters>
```
and analogously for the `gen_script_3dlap.sh` and other variants.

## Logfiles

The folders contain the complete logfiles of the runs used to generate the results included in the paper.
The MATLAB files do the parsing of the results and the generation of the figures.
