#!/usr/bin/bash -l
#SBATCH --job-name 3dlap
#SBATCH --partition boost_usr_prod
#SBATCH --time 01:00:00
#SBATCH --nodes @NNODES@
#SBATCH --gres=gpu:@NGPUS@
#SBATCH --ntasks=@NTASKS@
#SBATCH --ntasks-per-node=4
#SBATCH --export=NONE
#SBATCH -A CNHPC_1465132
@MORETHANONENODE@

unset SLURM_EXPORT_ENV

# Load environment
module load openmpi/4.1.6--gcc--12.2.0 openblas/0.3.24--gcc--12.2.0 gcc/12.2.0 cuda/12.1
module list

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/leonardo/home/userexternal/pdambra0/polynomialsmoothers/install/lib

cd ${HOME}/polynomialsmoothers/logfiles/leonardo

idim=@SIZE@
psize=@NGPUS@

srun ./3dlaplacian_multi >> 3dlap/match/multi/log_poly_match_l1jac_${idim}_task_${psize}_thr_1.txt 2>&1  <<EOF
%%%%%%%%%%%  General  arguments % Lines starting with % are ignored.
HLG                         ! Storage format CSR COO JAD
${idim}                     ! IDIM; domain size. Linear system size is IDIM**3
CONST                       ! PDECOEFF: CONST, EXP, BOX, GAUSS Coefficients of the PDE
FCG                         ! Iterative method: BiCGSTAB BiCGSTABL BiCG CG CGS FCG GCR RGMRES
2                           ! ISTOPC
00100                       ! ITMAX
1                           ! ITRACE
30                          ! IRST (restart for RGMRES and BiCGSTABL)
1.d-6                       ! EPS
%%%%%%%%%%%   %%%%%%%%%%%%%%%%
ML                        ! Preconditioner type: NONE JACOBI GS FBGS BJAC AS ML POLY
%%%%%%%%%%%  Multilevel parameters %%%%%%%%%%%%%%%%
VCYCLE                      ! Type of multilevel CYCLE: VCYCLE WCYCLE KCYCLE MULT ADD
1                           ! Number of outer sweeps for ML
-3                          ! Max Number of levels in a multilevel preconditioner; if <0, lib default
-3                          ! Target coarse matrix size per process; if <0, lib default
DIST                        ! Coarsest-level matrix distribution: DIST  REPL
SMOOTHED                    ! Type of aggregation: SMOOTHED UNSMOOTHED
COUPLED                     ! Parallel aggregation: DEC, SYMDEC, COUPLED
MATCHBOXP                   ! aggregation measure SOC1, MATCHBOXP
8                           ! Requested size of the aggregates for MATCHBOXP
NATURAL                     ! Ordering of aggregation NATURAL DEGREE
-1.5                        ! Coarsening ratio, if < 0 use library default
FILTER                      ! Filtering of matrix:  FILTER NOFILTER
-0.0100d0                   ! Smoothed aggregation threshold, ignored if < 0
-2                          ! Number of thresholds in vector, next line ignored if <= 0
0.05 0.025                  ! Thresholds
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-4CHEBY4-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
4                           ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
4                           ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-6CHEBY4-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
6                           ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-8CHEBY4-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
8                           ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
8                           ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-10CHEBY4-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
10                          ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_BA
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                   ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
1                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                        ! Second (post) smoother, ignored if NONE
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                   ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-12CHEBY4-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
12                          ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_LOTTES                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-4CHEBY4O-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
4                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
4                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-6CHEBY4O-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
6                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-8CHEBY4O-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
8                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
8                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-10CHEBY4O-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
10                          ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_BA
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                   ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
1                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                        ! Second (post) smoother, ignored if NONE
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                   ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-12CHEBY4O-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
12                          ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_LOTTES_BETA                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-4CHEBY1O-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
4                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
4                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-6CHEBY1O-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
6                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-8CHEBY1O-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
8                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
8                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-10CHEBY1O-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
10                          ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_BA
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                   ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
1                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Second smoother, always ignored for non-ML  %%%%%%%%%%%%%%%%
NONE                        ! Second (post) smoother, ignored if NONE
1                           ! Number of sweeps for (post) smoother
10                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                   ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-12CHEBY1O-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
POLY                        ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
1                           ! Number of sweeps for smoother
12                          ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
12                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-4L1JAC-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
4                           ! Number of sweeps for smoother
1                          ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-6L1JAC-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
6                           ! Number of sweeps for smoother
1                          ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-8L1JAC-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
8                           ! Number of sweeps for smoother
1                          ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-10L1JAC-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
10                          ! Number of sweeps for smoother
1                          ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
ML-VMATCH-12L1JAC-L1JAC30   ! Longer descriptive name for preconditioner (up to 20 chars)
%%%%%%%%%%%  First smoother (for all levels but coarsest) %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Smoother type JACOBI FBGS GS BWGS BJAC AS POLY  r 1-level, repeats previous.
12                           ! Number of sweeps for smoother
1                          ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_BA
1.0
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
1                           ! Number of sweeps for (post) smoother
6                           ! degree for polynomial smoother
POLY_NEW                 ! Polynomial variant
POLY_RHO_EST_POWER
1.0
0                           ! Number of overlap layers for AS preconditioner
HALO                        ! AS restriction operator: NONE HALO
NONE                        ! AS prolongation operator: NONE SUM AVG
L1-JACOBI                  ! Subdomain solver for BJAC/AS: JACOBI GS BGS ILU ILUT MILU MUMPS SLU UMF
1                           ! Inner solver sweeps (GS and JACOBI)
LLK                         ! AINV variant
0                           ! Fill level P for ILU(P) and ILU(T,P)
8                           ! Inverse Fill level P for INVK
1.d-4                       ! Threshold T for ILU(T,P)
%%%%%%%%%%%  Coarse level solver  %%%%%%%%%%%%%%%%
L1-JACOBI                   ! Coarsest-level solver: MUMPS UMF SLU SLUDIST JACOBI GS BJAC
ILU                         ! Coarsest-level subsolver for BJAC: ILU ILUT MILU UMF MUMPS SLU
1                           ! Coarsest-level fillin P for ILU(P) and ILU(T,P)
1.d-4                       ! Coarsest-level threshold T for ILU(T,P)
30                          ! Number of sweeps for JACOBI/GS/BJAC coarsest-level solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Next preconditioner name and smoothers  %%%%%%%%%%%%%%%%
END-OF-TESTS       ! Longer descriptive name for preconditioner (up to 20 chars)
EOF
