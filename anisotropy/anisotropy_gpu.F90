program anisopsblas
  use psb_base_mod
  use psb_ext_mod
  use psb_krylov_mod
  use psb_util_mod
#if defined(OPENMP)
  use omp_lib
#endif
#ifdef CUDA_MODE
  use psb_cuda_mod
#endif
  use amg_prec_mod
  use data_input
  implicit none

  ! do i want debug prints?
  logical :: doidebug = .false.

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt, pdecoeff
  integer(psb_ipk_) :: ne
  real(psb_dpk_)    :: theta,eps
  common /ex54_theta/ theta
  integer(psb_epk_) :: system_size

  ! miscellaneous
  real(psb_dpk_) :: t1, t2, tprec, thier, tslv

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a
  type(psb_ldspmat_type) :: aglob
  type(amg_dprec_type)  :: prec
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense vectors
  type(psb_d_vect_type) :: xvec,b,r
  ! parallel environment
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: iam, np, nl

  ! solver parameters
  integer(psb_ipk_)        :: iter, itmax,itrace, istopc, irst, nlv
  integer(psb_epk_) :: amatsize, precsize, descsize
  real(psb_dpk_)   :: err, resmx, resmxp

  ! Krylov solver data
  type solverdata
    character(len=40)  :: kmethd      ! Krylov solver
    integer(psb_ipk_)  :: istopc      ! stopping criterion
    integer(psb_ipk_)  :: itmax       ! maximum number of iterations
    integer(psb_ipk_)  :: itrace      ! tracing
    integer(psb_ipk_)  :: irst        ! restart
    real(psb_dpk_)     :: eps         ! stopping tolerance
  end type solverdata
  type(solverdata)       :: s_choice

  ! preconditioner data
  type precdata

    ! preconditioner type
    character(len=40)  :: descr       ! verbose description of the prec
    character(len=10)  :: ptype       ! preconditioner type

    integer(psb_ipk_)  :: outer_sweeps ! number of outer sweeps: sweeps for 1-level,
    ! AMG cycles for ML
    ! general AMG data
    character(len=32)  :: mlcycle      ! AMG cycle type
    integer(psb_ipk_)  :: maxlevs     ! maximum number of levels in AMG preconditioner

    ! AMG aggregation
    character(len=32)  :: aggr_prol    ! aggregation type: SMOOTHED, NONSMOOTHED
    character(len=32)  :: par_aggr_alg ! parallel aggregation algorithm: DEC, SYMDEC
    character(len=32)  :: aggr_type   ! Type of aggregation SOC1, SOC2, MATCHBOXP
    integer(psb_ipk_)  :: aggr_size   ! Requested size of the aggregates for MATCHBOXP
    character(len=32)  :: aggr_ord    ! ordering for aggregation: NATURAL, DEGREE
    character(len=32)  :: aggr_filter ! filtering: FILTER, NO_FILTER
    real(psb_dpk_)     :: mncrratio  ! minimum aggregation ratio
    real(psb_dpk_), allocatable :: athresv(:) ! smoothed aggregation threshold vector
    integer(psb_ipk_)  :: thrvsz      ! size of threshold vector
    real(psb_dpk_)     :: athres      ! smoothed aggregation threshold
    integer(psb_ipk_)  :: csizepp     ! minimum size of coarsest matrix per process

    ! AMG smoother or pre-smoother; also 1-lev preconditioner
    character(len=32)  :: smther      ! (pre-)smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps     ! (pre-)smoother / 1-lev prec. sweeps
    integer(psb_ipk_)  :: degree       ! degree for polynomial smoother
    character(len=32)  :: pvariant     ! polynomial  variant
    character(len=32)  :: prhovariant ! how to estimate rho(M^{-1}A)
    real(psb_dpk_)     :: prhovalue   ! if previous is set value, we set it from this one
    integer(psb_ipk_)  :: novr        ! number of overlap layers
    character(len=32)  :: restr       ! restriction over application of AS
    character(len=32)  :: prol        ! prolongation over application of AS
    character(len=32)  :: solve       ! local subsolver type: ILU, MILU, ILUT,
                                      ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    integer(psb_ipk_)  :: ssweeps      ! inner solver sweeps
    character(len=32)  :: variant     ! AINV variant: LLK, etc
    integer(psb_ipk_)  :: fill        ! fill-in for incomplete LU factorization
    integer(psb_ipk_)  :: invfill     ! Inverse fill-in for INVK
    real(psb_dpk_)     :: thr         ! threshold for ILUT factorization

    ! AMG post-smoother; ignored by 1-lev preconditioner
    character(len=32)  :: smther2      ! post-smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps2     ! post-smoother sweeps
    integer(psb_ipk_)  :: degree2      ! degree for polynomial smoother
    character(len=32)  :: pvariant2    ! polynomial  variant
    character(len=32)  :: prhovariant2 ! how to estimate rho(M^{-1}A)
    real(psb_dpk_)     :: prhovalue2   ! if previous is set value, we set it from this one
    integer(psb_ipk_)  :: novr2        ! number of overlap layers
    character(len=32)  :: restr2       ! restriction  over application of AS
    character(len=32)  :: prol2        ! prolongation over application of AS
    character(len=32)  :: solve2       ! local subsolver type: ILU, MILU, ILUT,
                                       ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    integer(psb_ipk_)  :: ssweeps2     ! inner solver sweeps
    character(len=32)  :: variant2     ! AINV variant: LLK, etc
    integer(psb_ipk_)  :: fill2        ! fill-in for incomplete LU factorization
    integer(psb_ipk_)  :: invfill2     ! Inverse fill-in for INVK
    real(psb_dpk_)      :: thr2         ! threshold for ILUT factorization

    ! coarsest-level solver
    character(len=32)  :: cmat        ! coarsest matrix layout: REPL, DIST
    character(len=32)  :: csolve      ! coarsest-lev solver: BJAC, SLUDIST (distr.
                                      ! mat.); UMF, MUMPS, SLU, ILU, ILUT, MILU
                                      ! (repl. mat.)
    character(len=32)  :: csbsolve    ! coarsest-lev local subsolver: ILU, ILUT,
                                      ! MILU, UMF, MUMPS, SLU
    character(len=32)  :: krmeth      ! Krylov method for subsolve
    integer(psb_ipk_)  :: cfill       ! fill-in for incomplete LU factorization
    real(psb_dpk_)     :: cthres      ! threshold for ILUT factorization
    integer(psb_ipk_)  :: cjswp       ! sweeps for GS or JAC coarsest-lev subsolver

    ! Dump data
    logical            :: dump = .false.
    integer(psb_ipk_)  :: dlmin       ! Minimum level to dump
    integer(psb_ipk_)  :: dlmax       ! Maximum level to dump
    logical            :: dump_ac         = .false.
    logical            :: dump_rp         = .false.
    logical            :: dump_tprol      = .false.
    logical            :: dump_smoother   = .false.
    logical            :: dump_solver     = .false.
    logical            :: dump_global_num = .false.

  end type precdata
  type(precdata)       :: p_choice

  ! Discretization variables
  real(psb_dpk_) :: h,blb(2),ev(2),thk,sg(3,9),x,y,coord(2,4)
  real(psb_dpk_) :: ss(4,4),shp(3,9),xsj,dd(2,2),a1,a2,val(1),ssvec(16)
  integer(psb_ipk_) :: nr,nlr,nnz,nel,ndf,f2,ll,f4,j1,ki,kj,i1
  integer(psb_lpk_) :: M,nt,geq,Istart,Iend,qj,qi,idx(4),irw(1)
  integer(psb_lpk_), allocatable :: vl(:)
  integer(psb_lpk_), allocatable     :: irow(:),icol(:),myidx(:)


  ! other variables
  integer(psb_ipk_)  :: info, i, k, nth, idxgpu, tngpus, xgpu
  character(len=20)  :: name,ch_err
  character(len=40)  :: matrixname

  ! To work with different matrix formats we use mold objects
  type(psb_d_ell_sparse_mat), target    :: aell
  type(psb_d_csr_sparse_mat), target   :: acsr
  type(psb_d_coo_sparse_mat), target   :: acoo
  type(psb_d_hll_sparse_mat), target    :: ahll
  type(psb_d_hdia_sparse_mat), target  :: ahdia
  type(psb_d_dns_sparse_mat), target   :: adns
  class(psb_d_base_sparse_mat), pointer :: amold => acsr
  class(psb_d_base_sparse_mat), pointer :: amold2 => acsr

  type(psb_d_base_vect_type), target   :: dvect
  class(psb_d_base_vect_type), pointer :: vmold => dvect
  type(psb_i_base_vect_type), target   :: ivect
  class(psb_i_base_vect_type), pointer :: imold => ivect
#ifdef CUDA_MODE
  type(psb_d_cuda_elg_sparse_mat), target   :: aelg
  type(psb_d_cuda_csrg_sparse_mat), target  :: acsrg
  type(psb_d_cuda_hlg_sparse_mat), target   :: ahlg
  type(psb_d_cuda_hdiag_sparse_mat), target :: ahdiag
  type(psb_d_vect_cuda), target         :: dvgpu
  type(psb_i_vect_cuda), target         :: ivgpu
#endif

  info=psb_success_

  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)
#ifdef CUDA_MODE
  call psb_cuda_init(ctxt,iam)
#endif

! openmp stuff
#if defined(OPENMP)
  !$OMP parallel shared(nth)
  !$OMP master
  nth = omp_get_num_threads()
  !$OMP end master
  !$OMP end parallel
#else
  nth = 1
#endif

  if (iam < 0) then
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif

  if(psb_get_errstatus() /= 0) goto 9999
  name='Rotated Anisotropy'
  call psb_set_errverbosity(itwo)
  !
  ! Hello world
  !
  if (iam == psb_root_) then
    write(psb_out_unit,*) 'Welcome to AMG4PSBLAS version: ',amg_version_string_
    write(psb_out_unit,*) 'This is the ',trim(name),' test program'
  end if

  !
  ! Reading parameters from input file
  !
  call get_parms(ctxt,afmt,ne,theta,eps,s_choice,p_choice)

#ifdef CUDA_MODE
  select case(psb_toupper(afmt))
  case('ELG')
    amold => aelg
  case('HLG')
    call psi_set_hksz(32)
    amold => ahlg
  case('HDIAG')
    amold => ahdiag
  case('CSRG')
    amold => acsrg
  case('ELL')
    amold => aell
  case('HLL')
    call psi_set_hksz(32)
    amold => ahll
  case('HDIA')
    amold => ahdia
  case('CSR')
    amold => acsr
  case('DNS')
    amold => adns
  case default
    write(*,*) 'Unknown format defaulting to HLG'
    amold => ahlg
  end select
  vmold  => dvgpu
  imold  => ivgpu
#else
  select case(psb_toupper(afmt))
  case('ELL')
    amold => aell
  case('HLL')
    call psi_set_hksz(32)
    amold => ahll
  case('HDIA')
    amold => ahdia
  case('CSR')
    amold => acsr
  case('DNS')
    amold => adns
  case default
    write(*,*) 'Unknown format defaulting to CSR'
    amold => acsr
  end select
#endif

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Solving ",I7," x ",I7," grid with θ = ",F16.4," ϵ = ",F16.4)') &
      & ne,ne,theta,eps
  end if


  !
  ! Building the discretization
  !
  h = 2.0/ne
  M = (ne+1)
  M = M*(ne+1)
  ! Block-row distribution *n   = bs*(Nbs/size + ((Nbs % size) > rank)); Nbs = *N/bs;
  ! nt = (M+np-1)/np !-1
  nt = M/np
  if ( iam < modulo(M,np) ) nt = nt + 1
  nr = max(0,min(nt,M-(iam*nt)))
  nnz = 7*((M+np-1)/np) !-1
  nt = nr
  call psb_sum(ctxt,nt)
  if (nt /= M) then
    write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,m
    info = -1
    call psb_barrier(ctxt)
    call psb_abort(ctxt)
    return
  end if

  call psb_cdall(ctxt,desc_a,info,nl=nr)
  myidx = desc_a%get_global_indices()
  nlr = size(myidx)

  !
  ! Parameters
  !
  f2 = 2
  f4 = 4
  theta = theta / 57.2957795
  ki = 2
  blb(1) = 0.0
  blb(2) = 0.0
  ev(1) = 1.0
  ev(2) = eps*ev(1)
  !
  ! Compute the matrix and right-hand-side vector that define the linear system
  !                                   Ax = b.
  !
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call psb_spall(a,desc_a,info,nnz,bldmode=psb_matbld_remote_,dupl=psb_dupl_add_)
  call psb_geall(b,desc_a,info)
  call psb_geall(xvec,desc_a,info)

  thk = 1.0              ! thickness
  nel = 4                ! nodes per element (quad)
  ndf = 1
  call int2d(f2,sg)
  Istart = myidx(1) -1
  Iend = myidx(nlr)
  if (doidebug) write(*,*)'Process ',iam,' Istart ',Istart,' Iend ',Iend
  do geq=Istart,Iend-1,1
     qj = geq/(ne+1); qi = mod(geq,(ne+1))
     x = h*qi - 1.0; y = h*qj - 1.0 ! lower left corner (-1,-1)
     if (qi < ne .and. qj < ne) then
        coord(1,1) = x;   coord(2,1) = y
        coord(1,2) = x+h; coord(2,2) = y
        coord(1,3) = x+h; coord(2,3) = y+h
        coord(1,4) = x;   coord(2,4) = y+h
        ! form stiff
        ss = 0.0
        if (doidebug) then
          write(0,*) 'spins : coord1 ',coord(1,:)
          write(0,*) 'spins : coord2 ',coord(2,:)
        end if
        do ll = 1,4
           call shp2dquad(sg(1,ll),sg(2,ll),coord,shp,xsj,f2)
           xsj = xsj*sg(3,ll)*thk
           if (doidebug) then
             write(0,*) 'spins : xsj ',xsj
           end if
           call thfx2d(ev,coord,shp,dd,f2,f2,f4,ex54_psi ,geq)
           if (doidebug) then
             write(0,*) 'spins : dd ',dd(1,:)
             write(0,*) 'spins : dd ',dd(2,:)
           end if

           j1 = 1
           do kj = 1,nel
              a1 = (dd(1,1)*shp(1,kj) + dd(1,2)*shp(2,kj))*xsj
              a2 = (dd(2,1)*shp(1,kj) + dd(2,2)*shp(2,kj))*xsj
              ! Compute residual
              ! p(j1) = p(j1) - a1*gradt(1) - a2*gradt(2)
              ! Compute tangent
              if (doidebug) then
                write(0,*) 'spins : a1a2 ',a1,a2
              end if

              i1 = 1
              do ki = 1,nel
                 ss(i1,j1) = ss(i1,j1) + a1*shp(1,ki) + a2*shp(2,ki)
                 i1 = i1 + ndf
              end do
              j1 = j1 + ndf
           end do
        enddo

        idx(1) = geq;
        idx(2) = geq+1;
        idx(3) = geq+(ne+1)+1
        idx(4) = geq+(ne+1)
!!$        write(0,*) iam,' IDX ',idx(1:4)+1
        call idx2psblas(idx,irow,icol)
        if (qj > 0) then
!            call MatSetValues(Amat,f4,idx,f4,idx,ss,ADD_VALUES,ierr)
          ssvec = reshape(ss,(/16/))

          call psb_spins(16, irow, icol, ssvec, a, desc_a, info)
!!$          do ki=1,16
!!$            if ((irow(ki) == 62)) then
!!$              write(0,*) iam,'Touching : ',irow(ki), icol(ki),ssvec(ki)
!!$            end if
!!$          end do
        else                !     a BC
          do ki=1,4,1
            do kj=1,4,1
              if (ki<3 .or. kj<3) then
                if (ki==kj) then
                  ss(ki,kj) = .1*ss(ki,kj)
                else
                  ss(ki,kj) = 0.0
                endif
              endif
            enddo
          enddo
          !            call MatSetValues(Amat,f4,idx,f4,idx,ss,ADD_VALUES,ierr)
          ssvec = reshape(ss,(/16/))
          call psb_spins(16, irow, icol, ssvec, a, desc_a, info)
!!$          do ki=1,16
!!$            if ((irow(ki) == 62)) then
!!$              write(0,*) iam,'Touching : ',irow(ki), icol(ki),ssvec(ki)
!!$            end if
!!$          end do
        endif               ! BC
        if (doidebug) then
          write(0,*) 'Spins : idx(1) ',idx(1),' ss ',ss(1,1:4)
          write(0,*) 'Spins : idx(2) ',idx(2),' ss ',ss(2,1:4)
          write(0,*) 'Spins : idx(3) ',idx(3),' ss ',ss(3,1:4)
          write(0,*) 'Spins : idx(4) ',idx(4),' ss ',ss(4,1:4)
        end if

     endif                  ! add element
     if (qj > 0) then      ! set rhs
        val(1) = h*h*exp(-100*((x+h/2)-blb(1))**2)*exp(-100*((y+h/2)-blb(2))**2)
!         call VecSetValues(bvec,one,geq,val,INSERT_VALUES,ierr)
        irw(1) = geq+1
        call psb_geins(1,irw,val,b,desc_a,info)
        val(1) = dzero
        call psb_geins(1,irw,val,xvec,desc_a,info)
     endif
  enddo
  flush(0)
  call psb_cdasb(desc_a, info, mold=imold)
  call psb_spasb(a, desc_a, info, mold=amold)
  call psb_geasb(xvec, desc_a, info, mold=vmold)
  call psb_geasb(b, desc_a, info, mold=vmold)

  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='anisopsblas'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (p_choice%dump) then

    call psb_gather(aglob,a,desc_a,info)
    if (iam == 0) then
      write(matrixname,'("A-np-",I1,".mtx")') np
      call mm_mat_write(aglob,"Anisotropy",info,filename=trim(matrixname))
    end if
  end if

  if (iam == psb_root_) &
       & write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if ( (iam == psb_root_).and.(a%is_dev()) ) then
    write(psb_out_unit,'("Matrix is on device memory")')
  elseif ( (iam == psb_root_).and.(a%is_host()) ) then
    write(psb_out_unit,'("Matrix is on host memory")')
  end if
  if ( (iam == psb_root_).and.(a%is_sync()) ) &
      & write(psb_out_unit,'("Matrix is in sync state")')
  if (iam == psb_root_) &
       & write(psb_out_unit,'(" ")')

  !
  ! Preconditioner
  !
  call prec%init(ctxt,p_choice%ptype,info)
  select case(trim(psb_toupper(p_choice%ptype)))
  case ('NONE','NOPREC')
    ! Do nothing, keep defaults

  case ('JACOBI','L1-JACOBI','GS','FWGS','FBGS')
    ! 1-level sweeps from "outer_sweeps"
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)

  case ('BJAC','POLY')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    call prec%set('solver_sweeps',   p_choice%ssweeps,   info)
    call prec%set('poly_degree',     p_choice%degree,    info)
    call prec%set('poly_variant',    p_choice%pvariant,  info)
    if (psb_toupper(p_choice%solve)=='MUMPS') &
         & call prec%set('mumps_loc_glob','local_solver',info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)

  case('AS')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_ovr',         p_choice%novr,    info)
    call prec%set('sub_restr',       p_choice%restr,   info)
    call prec%set('sub_prol',        p_choice%prol,    info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    call prec%set('solver_sweeps',   p_choice%ssweeps,   info)
    if (psb_toupper(p_choice%solve)=='MUMPS') &
         & call prec%set('mumps_loc_glob','local_solver',info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)

  case ('ML')
    ! multilevel preconditioner

    call prec%set('ml_cycle',        p_choice%mlcycle,    info)
    call prec%set('outer_sweeps',    p_choice%outer_sweeps,info)
    if (p_choice%csizepp>0)&
         & call prec%set('min_coarse_size_per_process', p_choice%csizepp,    info)
    if (p_choice%mncrratio>1)&
         & call prec%set('min_cr_ratio',   p_choice%mncrratio, info)
    if (p_choice%maxlevs>0)&
         & call prec%set('max_levs',    p_choice%maxlevs,    info)
    if (p_choice%athres >= dzero) &
         & call prec%set('aggr_thresh',     p_choice%athres,  info)
    if (p_choice%thrvsz>0) then
      do k=1,min(p_choice%thrvsz,size(prec%precv)-1)
        call prec%set('aggr_thresh',     p_choice%athresv(k),  info,ilev=(k+1))
      end do
    end if

    call prec%set('aggr_prol',       p_choice%aggr_prol,   info)
    call prec%set('par_aggr_alg',    p_choice%par_aggr_alg,   info)
    call prec%set('aggr_type',       p_choice%aggr_type, info)
    call prec%set('aggr_size',       p_choice%aggr_size, info)

    call prec%set('aggr_ord',        p_choice%aggr_ord,   info)
    call prec%set('aggr_filter',     p_choice%aggr_filter,info)


    call prec%set('smoother_type',   p_choice%smther,     info)
    call prec%set('smoother_sweeps', p_choice%jsweeps,    info)
    call prec%set('poly_degree',     p_choice%degree,    info)
    call prec%set('poly_variant',    p_choice%pvariant,  info)
    if (p_choice%prhovalue > dzero) then
        call prec%set('poly_rho_ba', p_choice%prhovalue, info)
    else
        call prec%set('poly_rho_estimate', p_choice%prhovariant,info)
    end if

    select case (psb_toupper(p_choice%smther))
    case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI','L1-FBGS')
      ! do nothing
    case default
      call prec%set('sub_ovr',         p_choice%novr,       info)
      call prec%set('sub_restr',       p_choice%restr,      info)
      call prec%set('sub_prol',        p_choice%prol,       info)
      select case(trim(psb_toupper(p_choice%solve)))
      case('INVK')
        call prec%set('sub_solve',       p_choice%solve,   info)
      case('INVT')
        call prec%set('sub_solve',       p_choice%solve,   info)
      case('AINV')
        call prec%set('sub_solve',       p_choice%solve,   info)
        call prec%set('ainv_alg', p_choice%variant,   info)
      case default
        call prec%set('sub_solve',       p_choice%solve,   info)
        if (psb_toupper(p_choice%solve)=='MUMPS') &
          & call prec%set('mumps_loc_glob','local_solver',info)
      end select
      call prec%set('solver_sweeps',   p_choice%ssweeps,    info)
      call prec%set('sub_fillin',      p_choice%fill,       info)
      call prec%set('inv_fillin',      p_choice%invfill,    info)
      call prec%set('sub_iluthrs',     p_choice%thr,        info)
    end select

    if (psb_toupper(p_choice%smther2) /= 'NONE') then
      call prec%set('smoother_type',   p_choice%smther2,   info,pos='post')
      call prec%set('smoother_sweeps', p_choice%jsweeps2,  info,pos='post')
      call prec%set('poly_degree',     p_choice%degree2,   info,pos='post')
      call prec%set('poly_variant',    p_choice%pvariant2, info,pos='post')
      if (p_choice%prhovalue > dzero) then
        call prec%set('poly_rho_ba', p_choice%prhovalue, info, pos='post')
      else
        call prec%set('poly_rho_estimate', p_choice%prhovariant,info, pos='post')
      end if

      select case (psb_toupper(p_choice%smther2))
      case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI','L1-FBGS')
        ! do nothing
      case default
        call prec%set('sub_ovr',         p_choice%novr2,     info,pos='post')
        call prec%set('sub_restr',       p_choice%restr2,    info,pos='post')
        call prec%set('sub_prol',        p_choice%prol2,     info,pos='post')
        select case(trim(psb_toupper(p_choice%solve2)))
        case('INVK')
          call prec%set('sub_solve',       p_choice%solve2,   info)
        case('INVT')
          call prec%set('sub_solve',       p_choice%solve2,   info)
        case('AINV')
          call prec%set('sub_solve',       p_choice%solve2,   info)
          call prec%set('ainv_alg', p_choice%variant2,   info)
        case default
          call prec%set('sub_solve',       p_choice%solve2,   info, pos='post')
          if (psb_toupper(p_choice%solve2)=='MUMPS') &
               & call prec%set('mumps_loc_glob','local_solver',info)
        end select
        call prec%set('solver_sweeps',   p_choice%ssweeps2,  info,pos='post')
        call prec%set('sub_fillin',      p_choice%fill2,     info,pos='post')
        call prec%set('inv_fillin',      p_choice%invfill2,  info,pos='post')
        call prec%set('sub_iluthrs',     p_choice%thr2,      info,pos='post')
      end select
    end if

    call prec%set('coarse_solve',    p_choice%csolve,    info)
    if (psb_toupper(p_choice%csolve) == 'BJAC') &
         &  call prec%set('coarse_subsolve', p_choice%csbsolve,  info)
    if (psb_toupper(p_choice%csolve) == 'KRM') then
	call prec%set('krm_method',p_choice%krmeth,info)
	call prec%set('krm_sub_solve',p_choice%csbsolve,  info)
    end if
    call prec%set('coarse_mat',      p_choice%cmat,      info)
    call prec%set('coarse_fillin',   p_choice%cfill,     info)
    call prec%set('coarse_iluthrs',  p_choice%cthres,    info)
    call prec%set('coarse_sweeps',   p_choice%cjswp,     info)

  end select

  ! build the preconditioner
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%hierarchy_build(a,desc_a,info)
  thier = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_hierarchy_bld')
    goto 9999
  end if
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%smoothers_build(a,desc_a,info,amold=amold,vmold=vmold,imold=imold)
  tprec = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_smoothers_bld')
    goto 9999
  end if

  call psb_amx(ctxt, thier)
  call psb_amx(ctxt, tprec)

  if(iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Preconditioner: ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Preconditioner time: ",es12.5)')thier+tprec
    write(psb_out_unit,'(" ")')
  end if
  call prec%memory_use(info,iout=psb_out_unit)
  flush(psb_out_unit)

  if (p_choice%dump) then
    call prec%dump(info,istart=p_choice%dlmin,iend=p_choice%dlmax,&
         & ac=p_choice%dump_ac,rp=p_choice%dump_rp,tprol=p_choice%dump_tprol,&
         & smoother=p_choice%dump_smoother, solver=p_choice%dump_solver, &
         & global_num=p_choice%dump_global_num)
  end if

  !
  ! iterative method parameters
  !
  call prec%allocate_wrk(info,vmold=vmold)
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call psb_krylov(s_choice%kmethd,a,prec,b,xvec,s_choice%eps,&
       & desc_a,info,itmax=s_choice%itmax,iter=iter,err=err,itrace=s_choice%itrace,&
       & istop=s_choice%istopc,irst=s_choice%irst)
  call psb_barrier(ctxt)
  tslv = psb_wtime() - t1

  call prec%free_wrk(info)
  call psb_amx(ctxt,tslv)

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='solver routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ctxt)
  tslv = psb_wtime() - t1
  call psb_amx(ctxt,tslv)
  !  write(0,*) iam,' after krylov',psb_errstatus_fatal()
  ! compute residual norms
  call psb_geall(r,desc_a,info)
  call r%zero()
  call psb_geasb(r,desc_a,info)
!  write(0,*) iam,' before geaxpby',psb_errstatus_fatal()
  call psb_geaxpby(done,b,dzero,r,desc_a,info)
  call psb_spmm(-done,a,xvec,done,r,desc_a,info)
  resmx  = psb_genrm2(r,desc_a,info)
  resmxp = psb_geamax(r,desc_a,info)
  !write(0,*) iam,' before descr',psb_errstatus_fatal()
  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  system_size = desc_a%get_global_rows()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  call psb_sum(ctxt,precsize)
  call prec%descr(info,iout=psb_out_unit)
  call prec%memory_use(info,iout=psb_out_unit)
#ifdef CUDA_MODE
!idxgpu = psb_cuda_getdevice()
!tngpus = psb_cuda_getdevicecount()
!if (idxgpu==0) then 
!  call psb_sum(ctxt,tngpus)
!else
!  xgpu = 0
!  call psb_sum(ctxt,xgpu)
!tngpus = xgpu
!end if
#endif
!$OMP PARALLEL
!$OMP SINGLE
 if ((iam == psb_root_)) then
    write(psb_out_unit,'("Computed solution on ",i8," processors")')  np
    write(psb_out_unit,'("Linear system size                 : ",i12)') system_size
    write(psb_out_unit,'("Theta                              : ",F16.5)') theta
    write(psb_out_unit,'("Anisotropy eps                     : ",F16.5)') eps
    write(psb_out_unit,'("Number of threads                  : ",i12)') nth
!   write(psb_out_unit,'("Number of gpus                     : ",i12)') tngpus
    write(psb_out_unit,'("Krylov method                      : ",a)') trim(s_choice%kmethd)
    write(psb_out_unit,'("Preconditioner                     : ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Iterations to convergence          : ",i12)')    iter
    write(psb_out_unit,'("Relative error estimate on exit    : ",es12.5)') err
    write(psb_out_unit,'("Number of levels in hierarchy      : ",i12)')    prec%get_nlevs()
    write(psb_out_unit,'("Time to build hierarchy            : ",es12.5)') thier
    write(psb_out_unit,'("Time to build smoothers            : ",es12.5)') tprec
    write(psb_out_unit,'("Total time for preconditioner      : ",es12.5)') tprec+thier
    write(psb_out_unit,'("Time to solve system               : ",es12.5)') tslv
    write(psb_out_unit,'("Time per iteration                 : ",es12.5)') tslv/iter
    write(psb_out_unit,'("Total time                         : ",es12.5)') tslv+tprec+thier
    write(psb_out_unit,'("Residual 2-norm                    : ",es12.5)') resmx
    write(psb_out_unit,'("Residual inf-norm                  : ",es12.5)') resmxp
    write(psb_out_unit,'("Total memory occupation for A      : ",i12)') amatsize
    write(psb_out_unit,'("Total memory occupation for DESC_A : ",i12)') descsize
    write(psb_out_unit,'("Total memory occupation for PREC   : ",i12)') precsize
    write(psb_out_unit,'("Storage format for A               : ",a  )') a%get_fmt()
    write(psb_out_unit,'("Storage format for DESC_A          : ",a  )') desc_a%get_fmt()

  end if
!$OMP END SINGLE
!$OMP END PARALLEL
  !
  !  cleanup storage and exit
  !
!  write(0,*) iam,' before spfree',psb_errstatus_fatal()
  call psb_spfree(a,desc_a,info)
!  write(0,*) iam,' before gefree 1',psb_errstatus_fatal()
  call psb_gefree(b,desc_a,info)
!  write(0,*) iam,' before gefree 2',psb_errstatus_fatal()
  call psb_gefree(xvec,desc_a,info)
!  write(0,*) iam,' before cdfree ',psb_errstatus_fatal()
  call psb_cdfree(desc_a,info)
!  write(0,*) iam,' before precfree ',psb_errstatus_fatal()
  call prec%free(info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
#ifdef CUDA_MODE
  call psb_cuda_exit()
#endif
  call psb_exit(ctxt)
  stop

9999 continue
  call psb_error(ctxt)

contains

  subroutine get_parms(ctxt,afmt,idim,theta,eps,solve,prec)

    implicit none

    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: idim
    real(psb_dpk_)      :: theta
    real(psb_dpk_)      :: eps
    character(len=*)    :: afmt
    type(solverdata)    :: solve
    type(precdata)      :: prec
    integer(psb_ipk_)   :: iam, nm, np, inp_unit
    character(len=1024)   :: filename

    call psb_info(ctxt,iam,np)

    if (iam == psb_root_) then
      if (command_argument_count()>0) then
        call get_command_argument(1,filename)
        inp_unit = 30
        open(inp_unit,file=filename,action='read',iostat=info)
        if (info /= 0) then
          write(psb_err_unit,*) 'Could not open file ',filename,' for input'
          call psb_abort(ctxt)
          stop
        else
          write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
        end if
      else
        inp_unit=psb_inp_unit
      end if
      ! read input data
      !
      call read_data(afmt,inp_unit)            ! matrix storage format
      call read_data(idim,inp_unit)            ! Discretization grid size
      call read_data(theta,inp_unit)           ! Theta Coefficient
      call read_data(eps,inp_unit)             ! Eps Coefficient
      ! Krylov solver data
      call read_data(solve%kmethd,inp_unit)    ! Krylov solver
      call read_data(solve%istopc,inp_unit)    ! stopping criterion
      call read_data(solve%itmax,inp_unit)     ! max num iterations
      call read_data(solve%itrace,inp_unit)    ! tracing
      call read_data(solve%irst,inp_unit)      ! restart
      call read_data(solve%eps,inp_unit)       ! tolerance
      ! preconditioner type
      call read_data(prec%descr,inp_unit)      ! verbose description of the prec
      call read_data(prec%ptype,inp_unit)      ! preconditioner type
      ! First smoother / 1-lev preconditioner
      call read_data(prec%smther,inp_unit)     ! smoother type
      call read_data(prec%jsweeps,inp_unit)    ! (pre-)smoother / 1-lev prec sweeps
      call read_data(prec%degree,inp_unit)     ! (pre-)smoother / 1-lev prec sweeps
      call read_data(prec%pvariant,inp_unit)   ! Variant of polynomial smoother
      call read_data(prec%prhovariant,inp_unit)! Variant of spectral radius estimate
      call read_data(prec%prhovalue,inp_unit)  ! User defined spectral radius estimate
      call read_data(prec%novr,inp_unit)       ! number of overlap layers
      call read_data(prec%restr,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve,inp_unit)      ! local subsolver
      call read_data(prec%ssweeps,inp_unit)    ! inner solver sweeps
      call read_data(prec%variant,inp_unit)    ! AINV variant
      call read_data(prec%fill,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%invfill,inp_unit)    !Inverse fill-in for INVK
      call read_data(prec%thr,inp_unit)        ! threshold for ILUT
      ! Second smoother/ AMG post-smoother (if NONE ignored in main)
      call read_data(prec%smther2,inp_unit)     ! smoother type
      call read_data(prec%jsweeps2,inp_unit)    ! (post-)smoother sweeps
      call read_data(prec%degree2,inp_unit)    ! (post-)smoother sweeps
      call read_data(prec%pvariant2,inp_unit)   ! Variant of polynomial smoother
      call read_data(prec%prhovariant2,inp_unit)! Variant of spectral radius estimate
      call read_data(prec%prhovalue2,inp_unit)  ! User defined spectral radius estimate   
      call read_data(prec%novr2,inp_unit)       ! number of overlap layers
      call read_data(prec%restr2,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol2,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve2,inp_unit)      ! local subsolver
      call read_data(prec%ssweeps2,inp_unit)    ! inner solver sweeps
      call read_data(prec%variant2,inp_unit)    ! AINV variant
      call read_data(prec%fill2,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%invfill2,inp_unit)    !Inverse fill-in for INVK
      call read_data(prec%thr2,inp_unit)        ! threshold for ILUT
      ! general AMG data
      call read_data(prec%mlcycle,inp_unit)     ! AMG cycle type
      call read_data(prec%outer_sweeps,inp_unit) ! number of 1lev/outer sweeps
      call read_data(prec%maxlevs,inp_unit)    ! max number of levels in AMG prec
      call read_data(prec%csizepp,inp_unit)       ! min size coarsest mat
      ! aggregation
      call read_data(prec%aggr_prol,inp_unit)    ! aggregation type
      call read_data(prec%par_aggr_alg,inp_unit)    ! parallel aggregation alg
      call read_data(prec%aggr_type,inp_unit)   ! type of aggregation
      call read_data(prec%aggr_size,inp_unit) ! Requested size of the aggregates for MATCHBOXP
      call read_data(prec%aggr_ord,inp_unit)    ! ordering for aggregation
      call read_data(prec%mncrratio,inp_unit)   ! minimum aggregation ratio
      call read_data(prec%aggr_filter,inp_unit) ! filtering
      call read_data(prec%athres,inp_unit)      ! smoothed aggr thresh
      call read_data(prec%thrvsz,inp_unit)      ! size of aggr thresh vector
      if (prec%thrvsz > 0) then
        call psb_realloc(prec%thrvsz,prec%athresv,info)
        call read_data(prec%athresv,inp_unit)   ! aggr thresh vector
      else
        read(inp_unit,*)                        ! dummy read to skip a record
      end if
      ! coasest-level solver
      call read_data(prec%csolve,inp_unit)      ! coarsest-lev solver
      call read_data(prec%csbsolve,inp_unit)    ! coarsest-lev subsolver
      call read_data(prec%krmeth,inp_unit)      ! type of Krylov method
      call read_data(prec%cmat,inp_unit)        ! coarsest mat layout
      call read_data(prec%cfill,inp_unit)       ! fill-in for incompl LU
      call read_data(prec%cthres,inp_unit)      ! Threshold for ILUT
      call read_data(prec%cjswp,inp_unit)       ! sweeps for GS/JAC subsolver
      ! dump
      call read_data(prec%dump,inp_unit)       ! Dump on file?
      call read_data(prec%dlmin,inp_unit)      ! Minimum level to dump
      call read_data(prec%dlmax,inp_unit)      ! Maximum level to dump
      call read_data(prec%dump_ac,inp_unit)
      call read_data(prec%dump_rp,inp_unit)
      call read_data(prec%dump_tprol,inp_unit)
      call read_data(prec%dump_smoother,inp_unit)
      call read_data(prec%dump_solver,inp_unit)
      call read_data(prec%dump_global_num,inp_unit)

      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if
    end if

    call psb_bcast(ctxt,afmt)
    call psb_bcast(ctxt,idim)
    call psb_bcast(ctxt,theta)
    call psb_bcast(ctxt,eps)

    call psb_bcast(ctxt,solve%kmethd)
    call psb_bcast(ctxt,solve%istopc)
    call psb_bcast(ctxt,solve%itmax)
    call psb_bcast(ctxt,solve%itrace)
    call psb_bcast(ctxt,solve%irst)
    call psb_bcast(ctxt,solve%eps)

    call psb_bcast(ctxt,prec%descr)
    call psb_bcast(ctxt,prec%ptype)

    ! broadcast first (pre-)smoother / 1-lev prec data
    call psb_bcast(ctxt,prec%smther)
    call psb_bcast(ctxt,prec%jsweeps)
    call psb_bcast(ctxt,prec%degree)
    call psb_bcast(ctxt,prec%pvariant)
    call psb_bcast(ctxt,prec%prhovariant)
    call psb_bcast(ctxt,prec%prhovalue)
    call psb_bcast(ctxt,prec%novr)
    call psb_bcast(ctxt,prec%restr)
    call psb_bcast(ctxt,prec%prol)
    call psb_bcast(ctxt,prec%solve)
    call psb_bcast(ctxt,prec%ssweeps)
    call psb_bcast(ctxt,prec%variant)
    call psb_bcast(ctxt,prec%fill)
    call psb_bcast(ctxt,prec%invfill)
    call psb_bcast(ctxt,prec%thr)
    ! broadcast second (post-)smoother
    call psb_bcast(ctxt,prec%smther2)
    call psb_bcast(ctxt,prec%jsweeps2)
    call psb_bcast(ctxt,prec%degree2)
    call psb_bcast(ctxt,prec%pvariant2)
    call psb_bcast(ctxt,prec%prhovariant2)
    call psb_bcast(ctxt,prec%prhovalue2)
    call psb_bcast(ctxt,prec%novr2)
    call psb_bcast(ctxt,prec%restr2)
    call psb_bcast(ctxt,prec%prol2)
    call psb_bcast(ctxt,prec%solve2)
    call psb_bcast(ctxt,prec%ssweeps2)
    call psb_bcast(ctxt,prec%variant2)
    call psb_bcast(ctxt,prec%fill2)
    call psb_bcast(ctxt,prec%invfill2)
    call psb_bcast(ctxt,prec%thr2)

    ! broadcast AMG parameters
    call psb_bcast(ctxt,prec%mlcycle)
    call psb_bcast(ctxt,prec%outer_sweeps)
    call psb_bcast(ctxt,prec%maxlevs)
    call psb_bcast(ctxt,prec%csizepp)

    call psb_bcast(ctxt,prec%aggr_prol)
    call psb_bcast(ctxt,prec%par_aggr_alg)
    call psb_bcast(ctxt,prec%aggr_type)
    call psb_bcast(ctxt,prec%aggr_size)
    call psb_bcast(ctxt,prec%aggr_ord)
    call psb_bcast(ctxt,prec%aggr_filter)
    call psb_bcast(ctxt,prec%mncrratio)
    call psb_bcast(ctxt,prec%thrvsz)
    if (prec%thrvsz > 0) then
      if (iam /= psb_root_) call psb_realloc(prec%thrvsz,prec%athresv,info)
      call psb_bcast(ctxt,prec%athresv)
    end if
    call psb_bcast(ctxt,prec%athres)

    call psb_bcast(ctxt,prec%cmat)
    call psb_bcast(ctxt,prec%csolve)
    call psb_bcast(ctxt,prec%csbsolve)
    call psb_bcast(ctxt,prec%krmeth)
    call psb_bcast(ctxt,prec%cfill)
    call psb_bcast(ctxt,prec%cthres)
    call psb_bcast(ctxt,prec%cjswp)
    ! dump
    call psb_bcast(ctxt,prec%dump)
    call psb_bcast(ctxt,prec%dlmin)
    call psb_bcast(ctxt,prec%dlmax)

    call psb_bcast(ctxt,prec%dump_ac)
    call psb_bcast(ctxt,prec%dump_rp)
    call psb_bcast(ctxt,prec%dump_tprol)
    call psb_bcast(ctxt,prec%dump_smoother)
    call psb_bcast(ctxt,prec%dump_solver)
    call psb_bcast(ctxt,prec%dump_global_num)


  end subroutine get_parms

  subroutine idx2psblas(ind,irowind,icolind)

    implicit none

    integer(psb_lpk_), intent(in) :: ind(4)
    integer(psb_lpk_), allocatable, intent(inout) :: irowind(:)
    integer(psb_lpk_), allocatable, intent(inout) :: icolind(:)

    if (.not.allocated(irowind)) then
      allocate(irowind(16))
    end if
    if (.not.allocated(icolind)) then
      allocate(icolind(16))
    end if

    irowind(1:4) = ind(1:4)+1
    irowind(5:8) = ind(1:4)+1
    irowind(9:12) = ind(1:4)+1
    irowind(13:16) = ind(1:4)+1

    icolind(1:4) = ind(1)+1
    icolind(5:8) = ind(2)+1
    icolind(9:12) = ind(3)+1
    icolind(13:16) = ind(4)+1

  end subroutine idx2psblas

  subroutine thfx2d(ev,xl,shp,dd,ndm,ndf,nel,dir,geq)
  implicit  none

  integer(psb_ipk_) ::   ndm,ndf,nel,i
  integer(psb_lpk_) ::   geq
  real(psb_dpk_) :: ev(2),xl(ndm,nel),shp(3,*),dir
  real(psb_dpk_) :: xx,yy,psi,cs,sn,c2,s2,dd(2,2)

  xx       = 0.0
  yy       = 0.0
  do i = 1,nel
    xx       = xx       + shp(3,i)*xl(1,i)
    yy       = yy       + shp(3,i)*xl(2,i)
  end do
  psi = dir(xx,yy)
!     Compute thermal flux
  cs  = cos(psi)
  sn  = sin(psi)
  c2  = cs*cs
  s2  = sn*sn
  cs  = cs*sn

  dd(1,1) = c2*ev(1) + s2*ev(2)
  dd(2,2) = s2*ev(1) + c2*ev(2)
  dd(1,2) = cs*(ev(1) - ev(2))
  dd(2,1) = dd(1,2)
! if (geq==0)  write(0,*) 'thfx2d p:',psi,ev(1),ev(2),xx,yy
! flux(1) = -dd(1,1)*gradt(1) - dd(1,2)*gradt(2)
! flux(2) = -dd(2,1)*gradt(1) - dd(2,2)*gradt(2)

  end subroutine thfx2d

! !     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! !     shp2dquad - shape functions - compute derivatives w/r natural coords.
! !     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine shp2dquad(s,t,xl,shp,xsj,ndm)
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Shape function routine for 4-node isoparametric quads
!
!      Inputs:
!         s,t       - Natural coordinates of point
!         xl(ndm,*) - Nodal coordinates for element
!         ndm       - Spatial dimension of mesh

!      Outputs:
!         shp(3,*)  - Shape functions and derivatives at point
!                     shp(1,i) = dN_i/dx  or dN_i/dxi_1
!                     shp(2,i) = dN_i/dy  or dN_i/dxi_2
!                     shp(3,i) = N_i
!         xsj       - Jacobian determinant at point
!-----[--.----+----.----+----.-----------------------------------------]
  implicit  none
  integer(psb_ipk_)  :: ndm
  real(psb_dpk_) :: xo,xs,xt, yo,ys,yt, xsm,xsp,xtm
  real(psb_dpk_) :: xtp, ysm,ysp,ytm,ytp
  real(psb_dpk_) :: s,t, xsj,xsj1, sh,th,sp,tp,sm
  real(psb_dpk_) :: tm, xl(ndm,4),shp(3,4)

!     Set up interpolations

  sh = 0.5*s
  th = 0.5*t
  sp = 0.5 + sh
  tp = 0.5 + th
  sm = 0.5 - sh
  tm = 0.5 - th
  shp(3,1) =   sm*tm
  shp(3,2) =   sp*tm
  shp(3,3) =   sp*tp
  shp(3,4) =   sm*tp

!     Set up natural coordinate functions (times 4)

  xo =  xl(1,1)-xl(1,2)+xl(1,3)-xl(1,4)
  xs = -xl(1,1)+xl(1,2)+xl(1,3)-xl(1,4) + xo*t
  xt = -xl(1,1)-xl(1,2)+xl(1,3)+xl(1,4) + xo*s
  yo =  xl(2,1)-xl(2,2)+xl(2,3)-xl(2,4)
  ys = -xl(2,1)+xl(2,2)+xl(2,3)-xl(2,4) + yo*t
  yt = -xl(2,1)-xl(2,2)+xl(2,3)+xl(2,4) + yo*s

!     Compute jacobian (times 16)

  xsj1 = xs*yt - xt*ys

!     Divide jacobian by 16 (multiply by .0625)

  xsj = 0.0625*xsj1
  if (xsj1.eq.0.0) then
     xsj1 = 1.0
  else
     xsj1 = 1.0/xsj1
  endif

!     Divide functions by jacobian

  xs  = (xs+xs)*xsj1
  xt  = (xt+xt)*xsj1
  ys  = (ys+ys)*xsj1
  yt  = (yt+yt)*xsj1

!     Multiply by interpolations

  ytm =  yt*tm
  ysm =  ys*sm
  ytp =  yt*tp
  ysp =  ys*sp
  xtm =  xt*tm
  xsm =  xs*sm
  xtp =  xt*tp
  xsp =  xs*sp

!     Compute shape functions

  shp(1,1) = - ytm+ysm
  shp(1,2) =   ytm+ysp
  shp(1,3) =   ytp-ysp
  shp(1,4) = - ytp-ysm
  shp(2,1) =   xtm-xsm
  shp(2,2) = - xtm-xsp
  shp(2,3) = - xtp+xsp
  shp(2,4) =   xtp+xsm

  end subroutine shp2dquad
!
! !     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! !     int2d
! !     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine int2d(l,sg)
!-----[--.----+----.----+----.-----------------------------------------]
!     Purpose: Form Gauss points and weights for two dimensions

!     Inputs:
!     l       - Number of points/direction

!     Outputs:
!     sg(3,*) - Array of points and weights
!-----[--.----+----.----+----.-----------------------------------------]
  implicit  none
  integer(psb_ipk_) :: l,i,lr(9),lz(9)
  real(psb_dpk_) :: g,third,sg(3,*)
  data      lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
  data      third / 0.3333333333333333 /

!     2x2 integration
  g = sqrt(third)
  do i = 1,4
     sg(1,i) = g*lr(i)
     sg(2,i) = g*lz(i)
     sg(3,i) = 1.0
  end do

  end subroutine int2d

!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     ex54_psi - anusotropic material direction
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  function ex54_psi(x,y)
  implicit  none
  real(psb_dpk_) :: ex54_psi
  real(psb_dpk_) :: x,y,theta
  common /ex54_theta/ theta
  ex54_psi = theta
  if (theta < 0.) then     ! circular
     if (y==0) then
        ex54_psi = 2.0*atan(1.0)
     else
        ex54_psi = atan(-x/y)
     endif
  endif
  end function ex54_psi

end program
