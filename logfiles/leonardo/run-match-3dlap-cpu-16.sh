#!/usr/bin/bash -l
#SBATCH --job-name 3dlap
#SBATCH --partition boost_usr_prod
#SBATCH --time 01:00:00
#SBATCH --nodes 4
#SBATCH --ntasks=16
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=4
#SBATCH --export=NONE
#SBATCH -A CNHPC_1465132


unset SLURM_EXPORT_ENV

# Load environment
module load openmpi/4.1.6--gcc--12.2.0 openblas/0.3.24--gcc--12.2.0 gcc/12.2.0 cuda/12.1
module list

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/leonardo/home/userexternal/pdambra0/polynomialsmoothers/install/lib

cd ${HOME}/polynomialsmoothers/logfiles/leonardo

degiter=4
idim=8034
psize=16

mpirun -np ${psize} ./3dlaplacian >> 3dlap/match/cpu/log_hgs_match_l1jac_${idim}_task_${psize}_thr_1.txt 2>&1  <<EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
CSR                     ! matrix storage format
${idim}                 ! Discretization grid size
CONST                   ! PDECOEFF: CONST, EXP, BOX, GAUSS Coefficients of the PDE
FCG                     ! Krylov Solver
2                       ! Stopping criterion
00500                     ! Maximum number of iterations
1                       ! Trace of FCG
30                      ! Restart (RGMRES e BICGSTAB)
1.d-7                   ! Tolerance
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-VSMATCH-HGS-30L1JAC ! verbose description of the prec
ML                      ! Preconditioner type
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
FBGS                    ! smoother type
1                       ! (pre-)smoother / 1-lev prec sweeps
${degiter}              ! degree for polynomial smoother
POLY_NEW                ! polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
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
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
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
COUPLED                 ! Parallel aggregation: DEC, SYMDEC, COUPLED
MATCHBOXP               ! aggregation measure SOC1, MATCHBOXP 
8                       ! Requested size of the aggregates for MATCHBOXP
NATURAL                 ! Ordering of aggregation NATURAL DEGREE
-1.5                    ! Coarsening ratio, if < 0 use library default
FILTER                  ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0               ! Smoothed aggregation threshold, ignored if < 0
-2                      ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025              ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI               ! coarsest-lev solver
L1-JACOBI               ! coarsest-lev subsolver
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

mpirun -np ${psize} ./3dlaplacian >> 3dlap/match/cpu/log_l1jacobi_match_l1jac_${idim}_task_${psize}_thr_1.txt 2>&1  <<EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
CSR                     ! matrix storage format
${idim}                 ! Discretization grid size
CONST                   ! PDECOEFF: CONST, EXP, BOX, GAUSS Coefficients of the PDE
FCG                     ! Krylov Solver
2                       ! Stopping criterion
00500                     ! Maximum number of iterations
1                       ! Trace of FCG
30                      ! Restart (RGMRES e BICGSTAB)
1.d-7                   ! Tolerance
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-VSMATCH-L1JAC-30L1JAC ! verbose description of the prec
ML                      ! Preconditioner type
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
L1-JACOBI               ! smoother type
4                       ! (pre-)smoother / 1-lev prec sweeps
1                       ! degree for polynomial smoother
POLY_NEW                ! polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
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
4                       ! Number of sweeps for (post) smoother
4                       ! degree for polynomial smoother
POLY_LOTTES_BETA        ! Polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
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
COUPLED                 ! Parallel aggregation: DEC, SYMDEC, COUPLED
MATCHBOXP               ! aggregation measure SOC1, MATCHBOXP 
8                       ! Requested size of the aggregates for MATCHBOXP
NATURAL                 ! Ordering of aggregation NATURAL DEGREE
-1.5                    ! Coarsening ratio, if < 0 use library default
FILTER                  ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0               ! Smoothed aggregation threshold, ignored if < 0
-2                      ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025              ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI               ! coarsest-lev solver
L1-JACOBI               ! coarsest-lev subsolver
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
