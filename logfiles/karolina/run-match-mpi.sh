#!/usr/bin/bash
#SBATCH --job-name anisopb
#SBATCH --account DD-23-167
#SBATCH --partition qcpu
#SBATCH --time 01:00:00
#SBATCH --nodes @NNODES@

# Load environment
module load gcccuda/2022b gompi/2022b CMake/3.24.3-GCCcore-12.2.0 OpenBLAS/0.3.21-GCC-12.2.0
module list

idim=@SIZE@
psize=@NTASKS@
theta=30
epsilon=100

srun -n ${psize} ./anisopsblas >> matchcpu/log_cheby4_match_l1jac_${idim}_task_${psize}_thr_1.txt 2>&1  <<EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
CSR                     ! matrix storage format
${idim}                 ! Discretization grid size
${theta}                ! Theta Coefficient (Degree)
${epsilon}              ! Epsilon of the anisotropy
FCG                     ! Krylov Solver
2                       ! Stopping criterion
500                     ! Maximum number of iterations
1                       ! Trace of FCG
500                     ! Restart (RGMRES e BICGSTAB)
1.d-7                   ! Tolerance
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-VMATCH-4CHEB4-30L1JAC ! verbose description of the prec
ML                      ! Preconditioner type
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                    ! smoother type
1                       ! (pre-)smoother / 1-lev prec sweeps
4                       ! degree for polynomial smoother
POLY_LOTTES             ! polynomial variant
0                       ! number of overlap layers
HALO                    ! restriction  over application of AS
NONE                    ! prolongation over application of AS
L1-JACOBI               ! local subsolver
1                       ! inner solver sweeps
LLK                     ! AINV variant
0                       ! Fill level P for ILU(P) and ILU(T,P)
1                       ! Inverse Fill level P for INVK
1.d-4                   ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                    ! Second (post) smoother, ignored if NONE
1                       ! Number of sweeps for (post) smoother
4                       ! degree for polynomial smoother
POLY_LOTTES_BETA        ! Polynomial variant
0                       ! Number of overlap layers for AS preconditioner
HALO                    ! AS restriction operator: NONE HALO
NONE                    ! AS prolongation operator: NONE SUM AVG
L1-JACOBI               ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                       ! Inner solver sweeps (GS and JACOBI)
LLK                     ! AINV variant
0                       ! Fill level P for ILU(P) and ILU(T,P)
8                       ! Inverse Fill level P for INVK
1.d-4                   ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Multilevel parameters %%%%%%%%%%%%%%%%
VCYCLE                  ! AMG cycle type
1                       ! number of 1lev/outer sweeps
-3                      ! Max Number of levels in a multilevel preconditioner; if <0, lib default
-3                      ! Target coarse matrix size per process; if <0, lib default
SMOOTHED                ! Type of aggregation: SMOOTHED UNSMOOTHED
COUPLED                     ! Parallel aggregation: DEC, SYMDEC, COUPLED
MATCHBOXP                    ! aggregation measure SOC1, MATCHBOXPcall read_data(prec%aggr_size,inp_unit) ! Requested size of the aggregates for MATCHBOXP
4                       ! Requested size of the aggregates for MATCHBOXP
NATURAL                 ! Ordering of aggregation NATURAL DEGREE
-1.5                    ! Coarsening ratio, if < 0 use library default
FILTER                  ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0               ! Smoothed aggregation threshold, ignored if < 0
-2                      ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025              ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI               ! coarsest-lev solver
L1-JACOBI               ! coarsest-lev subsolver
FCG                     ! type of Krylov method
DIST                    ! coarsest mat layout
0                       ! fill-in for incompl LU
1.d-2                   ! Threshold for ILUT
30                      ! sweeps for GS/JAC subsolver
%%%%%%%%%%%  Dump parms %%%%%%%%%%%%%%%%%%%%%%%%%%
F                       ! Dump on file?
1                       ! Minimum level to dump
10                      ! Maximum level to dump
T                       ! Dump coarse matrices
T                       ! Dump restrictors
T                       ! Dump prolongators
T                       ! Dump smoothers
T                       ! Dump coarse solver
T                       ! Dump using global numbering?
EOF

srun -n ${psize} ./anisopsblas >> matchcpu/log_optcheby4_match_l1jac_${idim}_task_${psize}_thr_1.txt 2>&1  <<EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
CSR                     ! matrix storage format
${idim}                 ! Discretization grid size
${theta}                ! Theta Coefficient (Degree)
${epsilon}              ! Epsilon of the anisotropy
FCG                     ! Krylov Solver
2                       ! Stopping criterion
500                     ! Maximum number of iterations
1                       ! Trace of FCG
500                     ! Restart (RGMRES e BICGSTAB)
1.d-7                   ! Tolerance
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-VMATCH-4OPTCHEB4-30L1JAC ! verbose description of the prec
ML                      ! Preconditioner type
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                    ! smoother type
1                       ! (pre-)smoother / 1-lev prec sweeps
4                       ! degree for polynomial smoother
POLY_LOTTES_BETA        ! polynomial variant
0                       ! number of overlap layers
HALO                    ! restriction  over application of AS
NONE                    ! prolongation over application of AS
L1-JACOBI               ! local subsolver
1                       ! inner solver sweeps
LLK                     ! AINV variant
0                       ! Fill level P for ILU(P) and ILU(T,P)
1                       ! Inverse Fill level P for INVK
1.d-4                   ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                    ! Second (post) smoother, ignored if NONE
1                       ! Number of sweeps for (post) smoother
4                       ! degree for polynomial smoother
POLY_LOTTES_BETA        ! Polynomial variant
0                       ! Number of overlap layers for AS preconditioner
HALO                    ! AS restriction operator: NONE HALO
NONE                    ! AS prolongation operator: NONE SUM AVG
L1-JACOBI               ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                       ! Inner solver sweeps (GS and JACOBI)
LLK                     ! AINV variant
0                       ! Fill level P for ILU(P) and ILU(T,P)
8                       ! Inverse Fill level P for INVK
1.d-4                   ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Multilevel parameters %%%%%%%%%%%%%%%%
VCYCLE                  ! AMG cycle type
1                       ! number of 1lev/outer sweeps
-3                      ! Max Number of levels in a multilevel preconditioner; if <0, lib default
-3                      ! Target coarse matrix size per process; if <0, lib default
SMOOTHED                ! Type of aggregation: SMOOTHED UNSMOOTHED
COUPLED                     ! Parallel aggregation: DEC, SYMDEC, COUPLED
MATCHBOXP                    ! aggregation measure SOC1, MATCHBOXPcall read_data(prec%aggr_size,inp_unit) ! Requested size of the aggregates for MATCHBOXP
4                       ! Requested size of the aggregates for MATCHBOXP
NATURAL                 ! Ordering of aggregation NATURAL DEGREE
-1.5                    ! Coarsening ratio, if < 0 use library default
FILTER                  ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0               ! Smoothed aggregation threshold, ignored if < 0
-2                      ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025              ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI               ! coarsest-lev solver
L1-JACOBI               ! coarsest-lev subsolver
FCG                     ! type of Krylov method
DIST                    ! coarsest mat layout
0                       ! fill-in for incompl LU
1.d-2                   ! Threshold for ILUT
30                      ! sweeps for GS/JAC subsolver
%%%%%%%%%%%  Dump parms %%%%%%%%%%%%%%%%%%%%%%%%%%
F                       ! Dump on file?
1                       ! Minimum level to dump
10                      ! Maximum level to dump
T                       ! Dump coarse matrices
T                       ! Dump restrictors
T                       ! Dump prolongators
T                       ! Dump smoothers
T                       ! Dump coarse solver
T                       ! Dump using global numbering?
EOF

srun -n ${psize} ./anisopsblas >> matchcpu/log_optcheby1_match_l1jac_${idim}_task_${psize}_thr_1.txt 2>&1  <<EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
CSR                     ! matrix storage format
${idim}                 ! Discretization grid size
${theta}                ! Theta Coefficient (Degree)
${epsilon}              ! Epsilon of the anisotropy
FCG                     ! Krylov Solver
2                       ! Stopping criterion
500                     ! Maximum number of iterations
1                       ! Trace of FCG
500                     ! Restart (RGMRES e BICGSTAB)
1.d-7                   ! Tolerance
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-VMATCH-4OPTCHEB1-30L1JAC ! verbose description of the prec
ML                      ! Preconditioner type
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                    ! smoother type
1                       ! (pre-)smoother / 1-lev prec sweeps
4                       ! degree for polynomial smoother
POLY_LOTTES_NEW         ! polynomial variant
0                       ! number of overlap layers
HALO                    ! restriction  over application of AS
NONE                    ! prolongation over application of AS
L1-JACOBI               ! local subsolver
1                       ! inner solver sweeps
LLK                     ! AINV variant
0                       ! Fill level P for ILU(P) and ILU(T,P)
1                       ! Inverse Fill level P for INVK
1.d-4                   ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                    ! Second (post) smoother, ignored if NONE
1                       ! Number of sweeps for (post) smoother
4                       ! degree for polynomial smoother
POLY_LOTTES_BETA        ! Polynomial variant
0                       ! Number of overlap layers for AS preconditioner
HALO                    ! AS restriction operator: NONE HALO
NONE                    ! AS prolongation operator: NONE SUM AVG
L1-JACOBI               ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                       ! Inner solver sweeps (GS and JACOBI)
LLK                     ! AINV variant
0                       ! Fill level P for ILU(P) and ILU(T,P)
8                       ! Inverse Fill level P for INVK
1.d-4                   ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Multilevel parameters %%%%%%%%%%%%%%%%
VCYCLE                  ! AMG cycle type
1                       ! number of 1lev/outer sweeps
-3                      ! Max Number of levels in a multilevel preconditioner; if <0, lib default
-3                      ! Target coarse matrix size per process; if <0, lib default
SMOOTHED                ! Type of aggregation: SMOOTHED UNSMOOTHED
COUPLED                     ! Parallel aggregation: DEC, SYMDEC, COUPLED
MATCHBOXP                    ! aggregation measure SOC1, MATCHBOXPcall read_data(prec%aggr_size,inp_unit) ! Requested size of the aggregates for MATCHBOXP
4                       ! Requested size of the aggregates for MATCHBOXP
NATURAL                 ! Ordering of aggregation NATURAL DEGREE
-1.5                    ! Coarsening ratio, if < 0 use library default
FILTER                  ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0               ! Smoothed aggregation threshold, ignored if < 0
-2                      ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025              ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI               ! coarsest-lev solver
L1-JACOBI               ! coarsest-lev subsolver
FCG                     ! type of Krylov method
DIST                    ! coarsest mat layout
0                       ! fill-in for incompl LU
1.d-2                   ! Threshold for ILUT
30                      ! sweeps for GS/JAC subsolver
%%%%%%%%%%%  Dump parms %%%%%%%%%%%%%%%%%%%%%%%%%%
F                       ! Dump on file?
1                       ! Minimum level to dump
10                      ! Maximum level to dump
T                       ! Dump coarse matrices
T                       ! Dump restrictors
T                       ! Dump prolongators
T                       ! Dump smoothers
T                       ! Dump coarse solver
T                       ! Dump using global numbering?
EOF
