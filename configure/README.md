# Machines

Each of the subfolders here contains the configure script for different machines.

| Machine  | Specifics | GPU and Accelerators |
|----------|-----------|----------------------|
| Alex | 2 x AMD EPYC 7713 ("Milan", "Zen3"), 2 x 64 cores @2.0 GHz <br> Main memory per node 1 TB | 8 x Nvidia A100 (40 GB HBM2) |
| Karolina | 720x 2x AMD 7H12, 64 cores, 2,6 GHz, 92,160 cores in total <br> 72x 2x AMD 7763, 64 cores, 2,45 GHz, 9,216 cores in total <br> 32x Intel Xeon-SC 8628, 24 cores, 2,9 GHz, 768 cores in total <br> 2x 2x AMD 7452, 32 cores, 2,35 GHz, 128 cores in total <br>  256 GB / 1 TB (GPU) / 24 TB fat node, 320 GB HBM2 (8 x 40 GB) GPU, Infiniband HDR 200 Gb/s | 72x 8x NVIDIA A100 GPU, 576 GPU in total <br> 36x 2x AMD 7H12, 64 cores, 2,6 GHz, 4,608 cores in total |
| Piz Daint | Intel® Xeon® E5-2690 v3 @ 2.60GHz (12 cores, 64GB RAM) - 5704 Nodes <br> Two Intel® Xeon® E5-2695 v4 @ 2.10GHz (2 x 18 cores, 64/128 GB RAM) - 1813 Nodes | NVIDIA Tesla P100 16GB |
| Toeplitz | 4 Nodes Intel(R) Xeon(R) CPU E5-2650 v4 @ 2.20GHz with 2 threads per core, 12 cores per socket and 2 socket with 250 Gb and 1 Node with Intel(R) Xeon(R) CPU E5-2643 v4 @ 3.40GHz with 2 threads per core, 6 cores per socket and 2 socket, 125 Gb | NO |
| Euler | Workstation with 12th Gen Intel(R) Core(TM) i7-12700 with 2 threads per core, 12 cores per socket, and 1 socket with 16 Gb  |  NVIDIA T1000 |
| x580gd | Laptop with Intel(R) Core(TM) i7-8750H CPU @ 2.20GHz with 2 threads per core, 6 cores per socket, and 1 socket with 16 Gb  |  NVIDIA GeForce GTX 1050 |
