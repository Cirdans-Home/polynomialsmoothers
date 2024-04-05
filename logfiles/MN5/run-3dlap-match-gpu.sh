#!/usr/bin/bash -l
#SBATCH --job-name 3dlap@SIZE@
#SBATCH --partition a100
#SBATCH --time 01:00:00
#SBATCH --nodes @NNODES@
#SBATCH --gres=gpu:a100:@NGPUS@
#SBATCH --ntasks=@NTASKS@
#SBATCH --export=NONE
@MORETHANONENODE@

unset SLURM_EXPORT_ENV

# Load environment
module load openmpi/4.1.6-gcc11.2.0-cuda gcc/11.2.0 cuda/12.1.1 
module list

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hpc/ihpc/ihpc100h/polynomialsmoothers/install/lib

cd ${HOME}/polynomialsmoothers/logfiles/alex

degiter=@DEGITER@
idim=@SIZE@
psize=@NTASKS@

mpirun -np ${psize} ./3dlaplacian >> 3dlap/match/${degiter}/log_cheby4_match_l1jac_${idim}_task_${psize}_thr_1.txt 2>&1  <<EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
HLG                         ! Storage format CSR COO JAD
${idim}                     ! IDIM; domain size. Linear system size is IDIM**3
CONST                       ! PDECOEFF: CONST, EXP, BOX, GAUSS Coefficients of the PDE
FCG                         ! Iterative method: BiCGSTAB BiCGSTABL BiCG CG CGS FCG GCR RGMRES
2                           ! ISTOPC
00500                       ! ITMAX
1                           ! ITRACE
30                          ! IRST (restart for RGMRES and BiCGSTABL)
1.d-7                       ! EPS
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-VSMATCH-${degiter}CHEBY4-30L1JAC ! Longer descriptive name for preconditioner (up to 20 chars)
ML                        ! Preconditioner type: NONE JACOBI GS FBGS BJAC AS ML POLY
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
${degiter}                  ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                     ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
1                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                        ! Second (post) smoother, ignored if NONE
6                           ! Number of sweeps for (post) smoother
1                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Multilevel parameters %%%%%%%%%%%%%%%%
VCYCLE                      ! Type of multilevel CYCLE: VCYCLE WCYCLE KCYCLE MULT ADD
1                           ! Number of outer sweeps for ML
-3                          ! Max Number of levels in a multilevel preconditioner; if <0, lib default
-3                          ! Target coarse matrix size per process; if <0, lib default
SMOOTHED                    ! Type of aggregation: SMOOTHED UNSMOOTHED
COUPLED                     ! Parallel aggregation: DEC, SYMDEC, COUPLED
MATCHBOXP                   ! aggregation measure SOC1, MATCHBOXP
8                           ! Requested size of the aggregates for MATCHBOXP
NATURAL                     ! Ordering of aggregation NATURAL DEGREE
-1.5                        ! Coarsening ratio, if < 0 use library default
FILTER                      ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0                    ! Smoothed aggregation threshold, ignored if < 0
-2                           ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025                  ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
DIST                        ! Coarsest-level matrix distribution: DIST  REPL
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%  Dump parms %%%%%%%%%%%%%%%%%%%%%%%%%%
F                           ! Dump preconditioner on file
1                           ! Min level 
20                          ! Max level 
T                           ! Dump AC
T                           ! Dump RP
F                           ! Dump TPROL
F                           ! Dump SMOOTHER
F                           ! Dump SOLVER
F                           ! Global numering ? 
EOF

mpirun -np ${psize} ./3dlaplacian >> 3dlap/match/${degiter}/log_optcheby4_match_l1jac_${idim}_task_${psize}_thr_1.txt 2>&1  <<EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
HLG                         ! Storage format CSR COO JAD
${idim}                     ! IDIM; domain size. Linear system size is IDIM**3
CONST                       ! PDECOEFF: CONST, EXP, BOX, GAUSS Coefficients of the PDE
FCG                         ! Iterative method: BiCGSTAB BiCGSTABL BiCG CG CGS FCG GCR RGMRES
2                           ! ISTOPC
00500                       ! ITMAX
1                           ! ITRACE
30                          ! IRST (restart for RGMRES and BiCGSTABL)
1.d-7                       ! EPS
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-VSMATCH-${degiter}OCHEBY4-30L1JAC ! Longer descriptive name for preconditioner (up to 20 chars)
ML                        ! Preconditioner type: NONE JACOBI GS FBGS BJAC AS ML POLY
%
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
${degiter}                  ! degree for polynomial smoother
POLY_LOTTES_BETA            ! Polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                     ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
1                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                        ! Second (post) smoother, ignored if NONE
6                           ! Number of sweeps for (post) smoother
1                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Multilevel parameters %%%%%%%%%%%%%%%%
VCYCLE                      ! Type of multilevel CYCLE: VCYCLE WCYCLE KCYCLE MULT ADD
1                           ! Number of outer sweeps for ML
-3                          ! Max Number of levels in a multilevel preconditioner; if <0, lib default
-3                          ! Target coarse matrix size per process; if <0, lib default
SMOOTHED                    ! Type of aggregation: SMOOTHED UNSMOOTHED
COUPLED                     ! Parallel aggregation: DEC, SYMDEC, COUPLED
MATCHBOXP                   ! aggregation measure SOC1, MATCHBOXP
8                           ! Requested size of the aggregates for MATCHBOXP
NATURAL                     ! Ordering of aggregation NATURAL DEGREE
-1.5                        ! Coarsening ratio, if < 0 use library default
FILTER                      ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0                    ! Smoothed aggregation threshold, ignored if < 0
-2                           ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025                  ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
DIST                        ! Coarsest-level matrix distribution: DIST  REPL
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%  Dump parms %%%%%%%%%%%%%%%%%%%%%%%%%%
F                           ! Dump preconditioner on file
1                           ! Min level 
20                          ! Max level 
T                           ! Dump AC
T                           ! Dump RP
F                           ! Dump TPROL
F                           ! Dump SMOOTHER
F                           ! Dump SOLVER
F                           ! Global numering ? 
EOF

mpirun -np ${psize} ./3dlaplacian >> 3dlap/match/${degiter}/log_optcheby1_match_l1jac_${idim}_task_${psize}_thr_1.txt 2>&1  <<EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
HLG                         ! Storage format CSR COO JAD
${idim}                     ! IDIM; domain size. Linear system size is IDIM**3
CONST                       ! PDECOEFF: CONST, EXP, BOX, GAUSS Coefficients of the PDE
FCG                         ! Iterative method: BiCGSTAB BiCGSTABL BiCG CG CGS FCG GCR RGMRES
2                           ! ISTOPC
00500                       ! ITMAX
1                           ! ITRACE
30                          ! IRST (restart for RGMRES and BiCGSTABL)
1.d-7                       ! EPS
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-VSMATCH-${degiter}OCHEBY1-30L1JAC ! Longer descriptive name for preconditioner (up to 20 chars)
ML                        ! Preconditioner type: NONE JACOBI GS FBGS BJAC AS ML POLY
%
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
${degiter}                  ! degree for polynomial smoother
POLY_NEW                    ! Polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                     ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
1                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                        ! Second (post) smoother, ignored if NONE
6                           ! Number of sweeps for (post) smoother
1                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Multilevel parameters %%%%%%%%%%%%%%%%
VCYCLE                      ! Type of multilevel CYCLE: VCYCLE WCYCLE KCYCLE MULT ADD
1                           ! Number of outer sweeps for ML
-3                          ! Max Number of levels in a multilevel preconditioner; if <0, lib default
-3                          ! Target coarse matrix size per process; if <0, lib default
SMOOTHED                    ! Type of aggregation: SMOOTHED UNSMOOTHED
COUPLED                     ! Parallel aggregation: DEC, SYMDEC, COUPLED
MATCHBOXP                   ! aggregation measure SOC1, MATCHBOXP
8                           ! Requested size of the aggregates for MATCHBOXP
NATURAL                     ! Ordering of aggregation NATURAL DEGREE
-1.5                        ! Coarsening ratio, if < 0 use library default
FILTER                      ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0                    ! Smoothed aggregation threshold, ignored if < 0
-2                           ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025                  ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
DIST                        ! Coarsest-level matrix distribution: DIST  REPL
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%  Dump parms %%%%%%%%%%%%%%%%%%%%%%%%%%
F                           ! Dump preconditioner on file
1                           ! Min level 
20                          ! Max level 
T                           ! Dump AC
T                           ! Dump RP
F                           ! Dump TPROL
F                           ! Dump SMOOTHER
F                           ! Dump SOLVER
F                           ! Global numering ? 
EOF

mpirun -np ${psize} ./3dlaplacian >> 3dlap/match/${degiter}/log_l1jacobi_match_l1jac_${idim}_task_${psize}_thr_1.txt 2>&1  <<EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
HLG                         ! Storage format CSR COO JAD
${idim}                     ! IDIM; domain size. Linear system size is IDIM**3
CONST                       ! PDECOEFF: CONST, EXP, BOX, GAUSS Coefficients of the PDE
FCG                         ! Iterative method: BiCGSTAB BiCGSTABL BiCG CG CGS FCG GCR RGMRES
2                           ! ISTOPC
00500                       ! ITMAX
1                           ! ITRACE
30                          ! IRST (restart for RGMRES and BiCGSTABL)
1.d-7                       ! EPS
%%%%%%%%%%%  Main preconditioner choices %%%%%%%%%%%%%%%%
ML-VSMATCH-${degiter}L1JAC-30L1JAC ! Longer descriptive name for preconditioner (up to 20 chars)
ML                        ! Preconditioner type: NONE JACOBI GS FBGS BJAC AS ML POLY
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
${degiter}                  ! Number of sweeps for smoother
1                           ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                     ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
1                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                        ! Second (post) smoother, ignored if NONE
6                           ! Number of sweeps for (post) smoother
1                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_EST_POWER      ! Algorithm to estimate spectral radius (ignored if next larger than 0)
1.0                     ! Spectral radius estimate
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Multilevel parameters %%%%%%%%%%%%%%%%
VCYCLE                      ! Type of multilevel CYCLE: VCYCLE WCYCLE KCYCLE MULT ADD
1                           ! Number of outer sweeps for ML
-3                          ! Max Number of levels in a multilevel preconditioner; if <0, lib default
-3                          ! Target coarse matrix size per process; if <0, lib default
SMOOTHED                    ! Type of aggregation: SMOOTHED UNSMOOTHED
COUPLED                     ! Parallel aggregation: DEC, SYMDEC, COUPLED
MATCHBOXP                   ! aggregation measure SOC1, MATCHBOXP
8                           ! Requested size of the aggregates for MATCHBOXP
NATURAL                     ! Ordering of aggregation NATURAL DEGREE
-1.5                        ! Coarsening ratio, if < 0 use library default
FILTER                      ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0                    ! Smoothed aggregation threshold, ignored if < 0
-2                           ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025                  ! Thresholds
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
DIST                        ! Coarsest-level matrix distribution: DIST  REPL
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%  Dump parms %%%%%%%%%%%%%%%%%%%%%%%%%%
F                           ! Dump preconditioner on file
1                           ! Min level 
20                          ! Max level 
T                           ! Dump AC
T                           ! Dump RP
F                           ! Dump TPROL
F                           ! Dump SMOOTHER
F                           ! Dump SOLVER
F                           ! Global numering ? 
EOF
