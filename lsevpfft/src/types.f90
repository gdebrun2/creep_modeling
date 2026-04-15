module types
    use global
    implicit none

    type voxel_type

        double precision :: wgtc = 0 ! current-config voxel weight (vol frac)
        double precision :: sg(3, 3) = 0 ! Cauchy stress tensor
        double precision :: velgrad(3, 3) = 0 ! velocity gradient L = grad(v)
        integer :: grain = 0 ! grain id
        integer :: phase = 0 ! phase id
        logical :: gas = .false. ! dummy gas voxel flag
        double precision :: cg66(6, 6) = 0 ! elastic stiffness (Voigt 6x6)
        double precision :: sg66(6, 6) = 0 ! elastic compliance (Voigt 6x6)
        double precision :: ag(3, 3) = 0 ! crystal->sample rotation matrix
        double precision :: gacumgr = 0 ! accumulated slip (hardening scalar)
        double precision :: crss(nsysmx, 2) = 0 ! slip resistance per system (fwd/bwd)
        double precision :: trialtau(nsysmx, 2) = 0 ! trial slip resistance (fwd/bwd)
        double precision :: kin(nsysmx) = 0 ! kinematic hardening/backstress
        double precision :: eth(3, 3) = 0 ! eigenstrain/thermal strain tensor
        double precision :: defgrad(3, 3) = 0 ! deformation gradient F
        double precision :: defgrade(3, 3) = 0 ! elastic deformation gradient Fe
        double precision :: defgradp(3, 3) = 0 ! plastic deformation gradient Fp
        double precision :: defgradinv(3, 3) = 0 ! inverse deformation gradient F^{-1}
        double precision :: edotp(3, 3) = 0 ! plastic strain-rate tensor
        double precision :: epvm = 0 ! accumulated von Mises plastic strain
        double precision :: sgt(3, 3) = 0 ! previous-step stress (elastic update)
        double precision :: disgradsym(3, 3) = 0 ! symmetric displacement gradient
        double precision :: detF = 0 ! det(F) Jacobian (local volume change)
        double precision :: defgradini(3, 3) = 0 ! initial deformation gradient F0
        double precision :: defgradiniinv(3, 3) = 0 ! inverse initial deformation gradient
        double precision :: xgr(3) = 0 ! current voxel center coordinates
        double precision :: gamdot(nsysmx) = 0 ! slip rates per system
        double precision :: sch(5, nsysmx) = 0 ! Schmid tensors (5-comp basis)
        double precision :: c066mod(6, 6) = 0 ! stiffness mapped by current F
        double precision :: c066modGoperr066mod(6, 6) = 0 ! AL operator term for stabilization
        double precision :: cgas(6, 6) = 0 ! gas-phase compliance operator (Voigt)
        double precision :: sgPK1(3, 3) = 0 ! PK1-like stress fluctuation for FFT
        double precision :: dvelgradref(3, 3) = 0 ! velgrad correction (reference config)
        double precision :: velgradold(3, 3) = 0 ! previous velgrad (mixing/relaxation)
        double precision :: ept(3, 3) = 0 ! accumulated plastic strain tensor
        double precision :: Lpinter(3, 3) = 0 ! plastic velocity gradient (intermediate)
        double precision :: defgradinc(3, 3) = 0 ! incremental F = exp(L*dt)
        double precision :: defgradeinc(3, 3) = 0 ! incremental Fe update
        integer :: igamdotmx(3) = 0 ! indices of top-3 |gamdot| systems
        double precision :: velgrad_snapshot(3, 3) = 0
        double precision :: velgradold_snapshot(3, 3) = 0
        double precision :: sg_snapshot(3, 3) = 0
        double precision :: edotp_snapshot(3, 3) = 0
        double precision :: gamdot_snapshot(nsysmx) = 0
        double precision :: trialtau_snapshot(nsysmx, 2) = 0

    end type voxel_type

    type micro_type

        double precision :: wgt = 0 ! reference voxel weight
        integer :: npts1 = 0, npts2 = 0, npts3 = 0 ! grid points in x,y,z
        double precision :: c066(6, 6) = 0 ! macroscopic tangent stiffness (Voigt)
        double precision :: s066(6, 6) = 0 ! macroscopic tangent compliance (Voigt)
        double precision :: c0(3, 3, 3, 3) = 0 ! macroscopic stiffness tensor (4th order)
        double precision :: s0(3, 3, 3, 3) = 0 ! macroscopic compliance tensor (4th order)
        double precision :: Goperr0(3, 3, 3, 3) = 0 ! Green operator kernel at origin (real-space)
        double precision :: wphcsol = 0 ! solid volume fraction (current config)
        double precision :: wphcgas = 0 ! gas volume fraction (current config)
        double precision :: defgradavg(3, 3) ! macroscopic average deformation gradient
        double precision :: defgradinvavgc(3, 3) ! avg F^{-1} (all voxels)
        double precision :: defgradinvavgcs(3, 3) ! avg F^{-1} (solid voxels)
        double precision :: defgradinvavgcg(3, 3) ! avg F^{-1} (gas voxels)
        double precision :: dvelgradavg(3, 3) = 0 ! BC correction to avg velgrad
        double precision :: velgradavg(3, 3) = 0 ! average velocity gradient
        double precision :: velgradavgsol(3, 3) = 0 ! avg velgrad in solid
        double precision :: velgradavggas(3, 3) = 0 ! avg velgrad in gas
        double precision :: sgavg(3, 3) = 0 ! average Cauchy stress (all voxels)
        double precision :: sgavgsol(3, 3) = 0 ! average Cauchy stress (solid only)
        double precision :: sgsolavg(3, 3) = 0 ! average solid stress tensor (output)
        double precision :: de_ds(3, 3, 3, 3) = 0 ! macro compliance tensor for BCs
        double precision :: svmavg = 0 ! von Mises stress of sgavg
        double precision :: dvmavg = 0 ! von Mises strain-rate of velgradavg(sym)
        double precision :: edotvmmx = 0 ! max voxel von Mises strain-rate
        double precision :: evmavg = 0 ! accumulated von Mises strain
        double precision :: epavg(3, 3) = 0 ! average plastic strain tensor
        double precision :: eavg(3, 3) = 0 ! average total strain tensor
        double precision :: eelavg(3, 3) = 0 ! average elastic strain tensor
        double precision :: edotpavg(3, 3) = 0 ! average plastic strain-rate tensor
        double precision :: evmpavg = 0 ! average von Mises plastic strain
        double precision :: dvmpavg = 0 ! average von Mises plastic strain-rate
        double precision :: lvmavg = 0 ! von Mises norm of velgradavg
        double precision, allocatable :: Goper(:, :, :, :, :, :, :) ! Fourier-space Green operator
        double precision, allocatable :: xnode(:, :, :, :) ! nodal coordinates (current)

        type(voxel_type), allocatable :: voxel(:, :, :) ! per-voxel state array

    end type micro_type

    type phase_type

        logical :: igas = .false. ! phase is dummy gas flag

        double precision :: dnca(3, nsysmx) = 0 ! slip plane normals (crystal basis)
        double precision :: dbca(3, nsysmx) = 0 ! slip directions (crystal basis)
        double precision :: schca(5, nsysmx) = 0 ! Schmid tensors in crystal basis (5)
        double precision :: tau(nsysmx, 3) = 0 ! tau params (tau0f, tau0b, tau1)
        double precision :: hard(nsysmx, nsysmx) = 0 ! latent hardening matrix
        double precision :: thet(nsysmx, 0:1) = 0 ! hardening rates (thet0/thet1)
        double precision :: nrs(nsysmx) = 0 ! rate sensitivity exponent per system
        double precision :: gamd0(nsysmx) = 0 ! reference shear rate per system
        double precision :: twsh(nsysmx) = 0 ! twinning shear per system
        double precision :: hpfac(nsysmx) = 0 ! Hall-Petch factor per system
        integer :: nsm(nmodmx) = 0 ! number of systems per mode
        integer :: nmodes = 0 ! number of deformation modes
        integer :: nsyst = 0 ! total number of systems
        integer :: ntwmod = 0 ! number of twinning modes
        integer :: ntwsys = 0 ! number of twinning systems
        integer :: isectw(nsysmx) = 0 ! twinning section id per system
        double precision :: C(3, 3, 3, 3) = 0 ! crystal elastic stiffness tensor

    end type phase_type

    type process_type

        integer :: ithermo = 0 ! eigenstrain/thermal strain flag
        integer :: iudot(3, 3) = 0 ! imposed udot components mask
        integer :: idsim(6) = 0 ! imposed strain-rate components mask (Voigt)
        integer :: iscau(6) = 0 ! imposed stress components mask (Voigt)
        integer :: ictrl = 0 ! timestep/strain control mode
        integer :: nsteps = 0 ! number of steps in this process
        integer :: itmax = 0 ! max FFT equilibrium iterations per step
        integer :: irecover = 0 ! recover initial stress from file flag
        integer :: isave = 0 ! save state flag
        integer :: iupdate = 0 ! update kinematics/stiffness flag
        integer :: iuphard = 0 ! update hardening within step flag
        integer :: iwtex = 0 ! write texture output flag
        integer :: iwfields = 0 ! write field outputs flag
        integer :: iwstep = 0 ! write field outputs cadence
        integer :: igamma = 0 ! derivative operator selection (spectral)
        double precision :: udot(3, 3) = 0 ! imposed macroscopic velocity gradient
        double precision :: dsim(3, 3) = 0 ! symmetric part of udot
        double precision :: tomtot(3, 3) = 0 ! antisymmetric part of udot
        double precision :: dvm = 0 ! von Mises magnitude of dsim
        double precision :: scauchy(3, 3) = 0 ! imposed macroscopic Cauchy stress
        double precision :: tdot = 0 ! current time increment dt
        double precision :: tdotref = 0 ! reference dt
        double precision :: tdotmin = 0 ! minimum allowed dt
        double precision :: devmmx = 0 ! max allowed von Mises increment
        double precision :: error = 0, erroral = 0 ! equilibrium / AL tolerances
        double precision :: xc0 = 0 ! reference stiffness scaling factor
        double precision :: devmapp = 0 ! target vm increment (strain control)
        double precision :: dsim_snapshot(3, 3) = 0
        double precision :: scauchy_snapshot(3, 3) = 0

    end type process_type

    type props_type

        integer :: npts1 = 0, npts2 = 0, npts3 = 0 ! grid points in x,y,z
        integer :: nph = 0 ! number of phases
        integer :: imicrosave = 0 ! restart microstructure step id
        integer :: nproc = 0 ! number of loading processes
        integer :: idispini = 0 ! initial displacement/defgrad option
        double precision :: delt(3) = 4 ! voxel spacing (dx,dy,dz)
        double precision :: wgt = 0 ! reference voxel weight
        double precision :: complgas = 0 ! gas compliance scaling
        double precision :: wphsol = 0 ! solid volume fraction (reference)
        double precision :: wphgas = 0 ! gas volume fraction (reference)
        logical :: res = .false. ! restart flag
        logical :: modAL = .true. ! modified-AL toggle

        type(phase_type), allocatable :: phase(:) ! per-phase material data
        type(process_type), allocatable :: process(:) ! per-process loading/BC data

    end type props_type

    type diag_type
        ! NOTE: diagnostics default to 0.0 unless explicitly updated.
        ! For min(detF), we use HUGE() as a sentinel during accumulation and
        ! convert it to 0.0 if a phase (solid/gas) is absent.

        ! min detF
        double precision :: detF_min = 0.d0
        double precision :: detF_min_sol = 0.d0
        double precision :: detF_min_gas = 0.d0
        double precision :: detF_min_prev = 0.0d0
        double precision :: detF_min_ratio = 1.0d0
        double precision :: detF_min_relstep = 0.0d0
        double precision :: detF_min_sol_prev = 0.0d0
        double precision :: detF_min_sol_ratio = 1.0d0
        double precision :: detF_min_sol_relstep = 0.0d0
        double precision :: detF_min_gas_prev = 0.0d0
        double precision :: detF_min_gas_ratio = 1.0d0
        double precision :: detF_min_gas_relstep = 0.0d0

        ! kinematic increment maxima from defgradinc (per step)
        double precision :: volinc_max_sol = 0.0d0        ! |log(det(Finc))|
        double precision :: volinc_max_gas = 0.0d0
        double precision :: volinc_max_sol_prev = 0.0d0
        double precision :: volinc_max_gas_prev = 0.0d0
        double precision :: volinc_max_sol_relstep = 0.0d0
        double precision :: volinc_max_gas_relstep = 0.0d0
        double precision :: stretchinc_max_sol = 0.0d0    ! ||Uinc - I||_F (polar decomposition)
        double precision :: stretchinc_max_gas = 0.0d0
        double precision :: stretchinc_max_sol_prev = 0.0d0
        double precision :: stretchinc_max_gas_prev = 0.0d0
        double precision :: stretchinc_max_sol_relstep = 0.0d0
        double precision :: stretchinc_max_gas_relstep = 0.0d0
        double precision :: rotang_max_sol = 0.0d0        ! rotation angle (radians) from Rinc
        double precision :: rotang_max_gas = 0.0d0
        double precision :: rotang_max_sol_prev = 0.0d0
        double precision :: rotang_max_gas_prev = 0.0d0

        ! per-step (macrostep) velgrad change diagnostics comparing the converged field
        ! to velgrad_snapshot (snapshot at macrostep start), so localized collapse shows up.
        double precision :: dL_step_rms_sol = 0.0d0
        double precision :: dL_step_rms_gas = 0.0d0
        double precision :: dL_step_rms = 0.0d0
        double precision :: dL_step_max_sol = 0.0d0
        double precision :: dL_step_max_gas = 0.0d0
        double precision :: dL_step_max = 0.0d0

        ! magnitude of the macroscopic BC correction to avg velgrad (Frobenius norm)
        double precision :: dvelgradavg_mag = 0.0d0

        ! final normalized field errors for this macrostep (as reported by the solver loop)
        double precision :: erre = 0.0d0
        double precision :: errs = 0.0d0
        double precision :: errsbc = 0.0d0

        ! reference timestep for ramping (as stored/updated in props%process(ipr)%tdotref)
        double precision :: tdotref = 0.0d0

        ! max von Mises strain increment over voxels (currently based on micro%edotvmmx, i.e. solid max)
        double precision :: devm_inc_max = 0.0d0

        ! macroscopic stress diagnostics (z-axis component)
        double precision :: s33_target = 0.0d0
        double precision :: s33 = 0.0d0
        double precision :: s33_sol = 0.0d0

        ! operator conditioning (set on finalize; "max" means worst-case)
        double precision :: cond_c066mod_max_sol = 0.0d0
        double precision :: cond_c066mod_max_gas = 0.0d0
        double precision :: cond_cgas_max = 0.0d0

        integer :: imin_sol(3) = 0, imin_gas(3) = 0

    contains
        procedure, pass(self) :: write_diagnostics

    end type diag_type

    public :: init_micro
    public :: update_diagnostics

contains

    subroutine init_micro(micro, props)

        implicit none
        type(micro_type), intent(inout) :: micro ! microstructure state to initialize
        type(props_type), intent(in) :: props ! input properties/configuration
        integer :: ip1, ip2, ip3, npts3node ! grid indices and z-node count

        micro%npts1 = props%npts1
        micro%npts2 = props%npts2
        micro%npts3 = props%npts3

        allocate (micro%Goper(3, 3, 3, 3, props%npts1, props%npts2, props%npts3))
        if (props%npts3 == 1) then
            npts3node = 1
        else
            npts3node = props%npts3 + 1
        end if
        allocate (micro%xnode(3, props%npts1 + 1, props%npts2 + 1, npts3node))
        micro%Goper = 0
        do ip3 = 1, npts3node
        do ip2 = 1, props%npts2 + 1
        do ip1 = 1, props%npts1 + 1
            micro%xnode(1, ip1, ip2, ip3) = float(ip1) - 0.5
            micro%xnode(2, ip1, ip2, ip3) = float(ip2) - 0.5
            if (npts3node == 1) then
                micro%xnode(3, ip1, ip2, ip3) = 0.0
            else
                micro%xnode(3, ip1, ip2, ip3) = float(ip3) - 0.5
            end if
        end do
        end do
        end do

        allocate (micro%voxel(props%npts1, props%npts2, props%npts3))
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1
            micro%voxel(ip1, ip2, ip3)%wgtc = props%wgt
            micro%voxel(ip1, ip2, ip3)%sg = 0
            micro%voxel(ip1, ip2, ip3)%velgrad = 0
            micro%voxel(ip1, ip2, ip3)%grain = 0
            micro%voxel(ip1, ip2, ip3)%phase = 0
            micro%voxel(ip1, ip2, ip3)%gas = .false.
            micro%voxel(ip1, ip2, ip3)%cg66 = 0
            micro%voxel(ip1, ip2, ip3)%sg66 = 0
            micro%voxel(ip1, ip2, ip3)%ag = 0
            micro%voxel(ip1, ip2, ip3)%gacumgr = 0
            micro%voxel(ip1, ip2, ip3)%crss = 0
            micro%voxel(ip1, ip2, ip3)%trialtau = 0
            micro%voxel(ip1, ip2, ip3)%kin = 0
            micro%voxel(ip1, ip2, ip3)%eth = 0
            micro%voxel(ip1, ip2, ip3)%edotp = 0
            micro%voxel(ip1, ip2, ip3)%epvm = 0
            micro%voxel(ip1, ip2, ip3)%sgt = 0
            micro%voxel(ip1, ip2, ip3)%disgradsym = 0
            micro%voxel(ip1, ip2, ip3)%gamdot = 0
            micro%voxel(ip1, ip2, ip3)%sch = 0
            micro%voxel(ip1, ip2, ip3)%c066mod = 0
            micro%voxel(ip1, ip2, ip3)%c066modGoperr066mod = 0
            micro%voxel(ip1, ip2, ip3)%cgas = 0
            micro%voxel(ip1, ip2, ip3)%sgPK1 = 0
            micro%voxel(ip1, ip2, ip3)%dvelgradref = 0
            micro%voxel(ip1, ip2, ip3)%velgradold = 0
            micro%voxel(ip1, ip2, ip3)%ept = 0
            micro%voxel(ip1, ip2, ip3)%Lpinter = 0
            micro%voxel(ip1, ip2, ip3)%igamdotmx = 0
            micro%voxel(ip1, ip2, ip3)%defgrad = id3
            micro%voxel(ip1, ip2, ip3)%defgrade = id3
            micro%voxel(ip1, ip2, ip3)%defgradp = id3
            micro%voxel(ip1, ip2, ip3)%defgradinv = id3
            micro%voxel(ip1, ip2, ip3)%detF = 1.0
            micro%voxel(ip1, ip2, ip3)%defgradini = id3
            micro%voxel(ip1, ip2, ip3)%defgradiniinv = id3
            micro%voxel(ip1, ip2, ip3)%defgradinc = id3
            micro%voxel(ip1, ip2, ip3)%defgradeinc = id3
            micro%voxel(ip1, ip2, ip3)%xgr(1) = float(ip1)
            micro%voxel(ip1, ip2, ip3)%xgr(2) = float(ip2)
            micro%voxel(ip1, ip2, ip3)%xgr(3) = float(ip3)
        end do
        end do
        end do

        micro%defgradavg = id3
        micro%defgradinvavgc = id3
        micro%defgradinvavgcs = id3
        micro%defgradinvavgcg = id3

    end subroutine init_micro

    subroutine update_diagnostics(micro, diag, reset, finalize, scan_cond)

        use tensor_functions, only: polar_dcmp
        implicit none

        type(micro_type), intent(in) :: micro
        type(diag_type), intent(inout) :: diag
        logical, intent(in), optional :: reset, finalize, scan_cond

        logical :: do_reset, do_finalize, do_scan_cond
        integer :: i, j, k, lwork, info, ii
        integer :: imin_sol(3), imin_gas(3)

        ! kinematics scratch
        double precision :: Finc(3, 3), R(3, 3), U(3, 3)
        double precision ::  Q(3, 3), w(3), work(64), tmp(3, 3), logU(3, 3)
        double precision :: Jinc, volinc, stretchinc, rotang, arg

        ! conditioning scratch
        double precision :: condA

        ! step-to-step velgrad change scratch
        double precision :: dL33(3, 3), dL, L0, denom, rel, wgtc
        double precision :: sum_rel2_sol, sum_w_sol, max_rel_sol
        double precision :: sum_rel2_gas, sum_w_gas, max_rel_gas
        double precision :: sum_rel2_all, sum_w_all, max_rel_all
        double precision :: Lref, Lfloor
        double precision :: detF_min_all, detF_min_sol_all, detF_min_gas_all
        double precision :: volinc_max_sol_all, volinc_max_gas_all
        double precision :: stretchinc_max_sol_all, stretchinc_max_gas_all
        double precision :: rotang_max_sol_all, rotang_max_gas_all
        double precision :: cond_c066mod_max_sol_all, cond_c066mod_max_gas_all, cond_cgas_max_all
        double precision :: detF_min_loc, detF_min_sol_loc, detF_min_gas_loc
        integer :: imin_sol_loc(3), imin_gas_loc(3)
        double precision :: sum_rel2_sol_loc, sum_w_sol_loc, max_rel_sol_loc
        double precision :: sum_rel2_gas_loc, sum_w_gas_loc, max_rel_gas_loc
        double precision :: sum_rel2_all_loc, sum_w_all_loc, max_rel_all_loc
        double precision :: volinc_max_sol_loc, volinc_max_gas_loc
        double precision :: stretchinc_max_sol_loc, stretchinc_max_gas_loc
        double precision :: rotang_max_sol_loc, rotang_max_gas_loc
        double precision :: cond_c066mod_max_sol_loc, cond_c066mod_max_gas_loc, cond_cgas_max_loc

        do_reset = .false.; if (present(reset)) do_reset = reset
        do_finalize = .false.; if (present(finalize)) do_finalize = finalize
        do_scan_cond = .false.; if (present(scan_cond)) do_scan_cond = scan_cond

        if (do_reset) then
            diag%detF_min = huge(1.0d0)
            diag%detF_min_sol = huge(1.0d0)
            diag%detF_min_gas = huge(1.0d0)
            diag%detF_min_ratio = 1.0d0
            diag%detF_min_relstep = 0.0d0
            diag%detF_min_sol_ratio = 1.0d0
            diag%detF_min_sol_relstep = 0.0d0
            diag%detF_min_gas_ratio = 1.0d0
            diag%detF_min_gas_relstep = 0.0d0

            diag%volinc_max_sol = 0.0d0
            diag%volinc_max_gas = 0.0d0
            diag%volinc_max_sol_relstep = 0.0d0
            diag%volinc_max_gas_relstep = 0.0d0
            diag%stretchinc_max_sol = 0.0d0
            diag%stretchinc_max_gas = 0.0d0
            diag%stretchinc_max_sol_relstep = 0.0d0
            diag%stretchinc_max_gas_relstep = 0.0d0
            diag%rotang_max_sol = 0.0d0
            diag%rotang_max_gas = 0.0d0

            diag%dL_step_rms_sol = 0.0d0
            diag%dL_step_rms_gas = 0.0d0
            diag%dL_step_rms = 0.0d0
            diag%dL_step_max_sol = 0.0d0
            diag%dL_step_max_gas = 0.0d0
            diag%dL_step_max = 0.0d0

            diag%dvelgradavg_mag = 0.0d0
            diag%erre = 0.0d0
            diag%errs = 0.0d0
            diag%errsbc = 0.0d0
            diag%tdotref = 0.0d0
            diag%devm_inc_max = 0.0d0
            diag%s33_target = 0.0d0
            diag%s33 = 0.0d0
            diag%s33_sol = 0.0d0

            diag%cond_c066mod_max_sol = 0.0d0
            diag%cond_c066mod_max_gas = 0.0d0
            diag%cond_cgas_max = 0.0d0
        end if

        if (do_finalize) then
            ! Recompute per-step diagnostics from scratch (independent of whether reset was called).
            detF_min_all = huge(1.0d0)
            detF_min_sol_all = huge(1.0d0)
            detF_min_gas_all = huge(1.0d0)
            volinc_max_sol_all = 0.0d0
            volinc_max_gas_all = 0.0d0
            stretchinc_max_sol_all = 0.0d0
            stretchinc_max_gas_all = 0.0d0
            rotang_max_sol_all = 0.0d0
            rotang_max_gas_all = 0.0d0

            imin_sol = 0
            imin_gas = 0

            cond_c066mod_max_sol_all = 0.0d0
            cond_c066mod_max_gas_all = 0.0d0
            cond_cgas_max_all = 0.0d0

            sum_rel2_sol = 0.0d0; sum_w_sol = 0.0d0; max_rel_sol = 0.0d0
            sum_rel2_gas = 0.0d0; sum_w_gas = 0.0d0; max_rel_gas = 0.0d0
            sum_rel2_all = 0.0d0; sum_w_all = 0.0d0; max_rel_all = 0.0d0
            Lref = max(micro%lvmavg, tiny(1.0d0))
            Lfloor = max(1.0d-12 * Lref, tiny(1.0d0))

            !$OMP PARALLEL DEFAULT(NONE) &
            !$OMP SHARED(micro, do_scan_cond, Lfloor, detF_min_all, detF_min_sol_all, detF_min_gas_all, imin_sol, imin_gas, &
            !$OMP        volinc_max_sol_all, volinc_max_gas_all, stretchinc_max_sol_all, stretchinc_max_gas_all, &
            !$OMP        rotang_max_sol_all, rotang_max_gas_all, cond_c066mod_max_sol_all, cond_c066mod_max_gas_all, cond_cgas_max_all, &
            !$OMP        sum_rel2_sol, sum_w_sol, max_rel_sol, sum_rel2_gas, sum_w_gas, max_rel_gas, sum_rel2_all, sum_w_all, max_rel_all) &
            !$OMP PRIVATE(i, j, k, lwork, info, ii, Finc, R, U, Q, w, work, tmp, logU, Jinc, volinc, stretchinc, rotang, arg, &
            !$OMP         condA, dL33, dL, L0, denom, rel, wgtc, detF_min_loc, detF_min_sol_loc, detF_min_gas_loc, &
            !$OMP         imin_sol_loc, imin_gas_loc, sum_rel2_sol_loc, sum_w_sol_loc, max_rel_sol_loc, &
            !$OMP         sum_rel2_gas_loc, sum_w_gas_loc, max_rel_gas_loc, sum_rel2_all_loc, sum_w_all_loc, max_rel_all_loc, &
            !$OMP         volinc_max_sol_loc, volinc_max_gas_loc, stretchinc_max_sol_loc, stretchinc_max_gas_loc, &
            !$OMP         rotang_max_sol_loc, rotang_max_gas_loc, cond_c066mod_max_sol_loc, cond_c066mod_max_gas_loc, cond_cgas_max_loc)

            lwork = size(work)
            detF_min_loc = huge(1.0d0)
            detF_min_sol_loc = huge(1.0d0)
            detF_min_gas_loc = huge(1.0d0)
            imin_sol_loc = 0
            imin_gas_loc = 0

            sum_rel2_sol_loc = 0.0d0; sum_w_sol_loc = 0.0d0; max_rel_sol_loc = 0.0d0
            sum_rel2_gas_loc = 0.0d0; sum_w_gas_loc = 0.0d0; max_rel_gas_loc = 0.0d0
            sum_rel2_all_loc = 0.0d0; sum_w_all_loc = 0.0d0; max_rel_all_loc = 0.0d0

            volinc_max_sol_loc = 0.0d0
            volinc_max_gas_loc = 0.0d0
            stretchinc_max_sol_loc = 0.0d0
            stretchinc_max_gas_loc = 0.0d0
            rotang_max_sol_loc = 0.0d0
            rotang_max_gas_loc = 0.0d0

            cond_c066mod_max_sol_loc = 0.0d0
            cond_c066mod_max_gas_loc = 0.0d0
            cond_cgas_max_loc = 0.0d0

            !$OMP DO COLLAPSE(3) SCHEDULE(static)
            do k = 1, micro%npts3
            do j = 1, micro%npts2
            do i = 1, micro%npts1

                dL33 = micro%voxel(i, j, k)%velgrad - micro%voxel(i, j, k)%velgrad_snapshot
                dL = sqrt(sum(dL33**2))
                L0 = sqrt(sum(micro%voxel(i, j, k)%velgrad_snapshot**2))
                denom = max(L0, Lfloor)
                rel = dL / denom
                wgtc = micro%voxel(i, j, k)%wgtc

                sum_rel2_all_loc = sum_rel2_all_loc + (rel * rel) * wgtc
                sum_w_all_loc = sum_w_all_loc + wgtc
                max_rel_all_loc = max(max_rel_all_loc, rel)

                detF_min_loc = min(detF_min_loc, micro%voxel(i, j, k)%detF)

                if (micro%voxel(i, j, k)%gas) then
                    if (micro%voxel(i, j, k)%detF <= detF_min_gas_loc) then
                        detF_min_gas_loc = micro%voxel(i, j, k)%detF
                        imin_gas_loc = (/i, j, k/)
                    end if
                    sum_rel2_gas_loc = sum_rel2_gas_loc + (rel * rel) * wgtc
                    sum_w_gas_loc = sum_w_gas_loc + wgtc
                    max_rel_gas_loc = max(max_rel_gas_loc, rel)
                else
                    if (micro%voxel(i, j, k)%detF <= detF_min_sol_loc) then
                        detF_min_sol_loc = micro%voxel(i, j, k)%detF
                        imin_sol_loc = (/i, j, k/)
                    end if
                    sum_rel2_sol_loc = sum_rel2_sol_loc + (rel * rel) * wgtc
                    sum_w_sol_loc = sum_w_sol_loc + wgtc
                    max_rel_sol_loc = max(max_rel_sol_loc, rel)
                end if

                Finc = micro%voxel(i, j, k)%defgradinc
                Jinc = determinant33(Finc)
                volinc = abs(log(max(Jinc, tiny(1.0d0))))

                call polar_dcmp(Finc, R, U)
                arg = 0.5d0 * ((R(1, 1) + R(2, 2) + R(3, 3)) - 1.0d0)
                arg = max(-1.0d0, min(1.0d0, arg))
                rotang = acos(arg) * 180 / pi
                Q = U
                call dsyev('V', 'U', 3, Q, 3, w, work, lwork, info)
                if (info /= 0) then
                    stretchinc = huge(1.0d0)
                else
                    tmp = 0.0d0
                    do ii = 1, 3
                        tmp(:, ii) = Q(:, ii) * log(max(w(ii), tiny(1.0d0)))
                    end do
                    logU = matmul(tmp, transpose(Q))
                    stretchinc = sqrt(sum(logU**2))
                end if

                if (micro%voxel(i, j, k)%gas) then
                    volinc_max_gas_loc = max(volinc_max_gas_loc, volinc)
                    rotang_max_gas_loc = max(rotang_max_gas_loc, rotang)
                    stretchinc_max_gas_loc = max(stretchinc_max_gas_loc, stretchinc)
                else
                    volinc_max_sol_loc = max(volinc_max_sol_loc, volinc)
                    rotang_max_sol_loc = max(rotang_max_sol_loc, rotang)
                    stretchinc_max_sol_loc = max(stretchinc_max_sol_loc, stretchinc)
                end if

                if (do_scan_cond) then
                    condA = cond1_est_6x6(micro%voxel(i, j, k)%c066mod)
                    if (micro%voxel(i, j, k)%gas) then
                        cond_c066mod_max_gas_loc = max(cond_c066mod_max_gas_loc, condA)
                        condA = cond1_est_6x6(micro%voxel(i, j, k)%cgas)
                        cond_cgas_max_loc = max(cond_cgas_max_loc, condA)
                    else
                        cond_c066mod_max_sol_loc = max(cond_c066mod_max_sol_loc, condA)
                    end if
                end if

            end do
            end do
            end do
            !$OMP END DO

            !$OMP CRITICAL(update_diagnostics_combine)
            detF_min_all = min(detF_min_all, detF_min_loc)

            if (detF_min_sol_loc < detF_min_sol_all) then
                detF_min_sol_all = detF_min_sol_loc
                imin_sol = imin_sol_loc
            end if
            if (detF_min_gas_loc < detF_min_gas_all) then
                detF_min_gas_all = detF_min_gas_loc
                imin_gas = imin_gas_loc
            end if

            sum_rel2_all = sum_rel2_all + sum_rel2_all_loc
            sum_w_all = sum_w_all + sum_w_all_loc
            max_rel_all = max(max_rel_all, max_rel_all_loc)

            sum_rel2_sol = sum_rel2_sol + sum_rel2_sol_loc
            sum_w_sol = sum_w_sol + sum_w_sol_loc
            max_rel_sol = max(max_rel_sol, max_rel_sol_loc)

            sum_rel2_gas = sum_rel2_gas + sum_rel2_gas_loc
            sum_w_gas = sum_w_gas + sum_w_gas_loc
            max_rel_gas = max(max_rel_gas, max_rel_gas_loc)

            volinc_max_sol_all = max(volinc_max_sol_all, volinc_max_sol_loc)
            volinc_max_gas_all = max(volinc_max_gas_all, volinc_max_gas_loc)
            rotang_max_sol_all = max(rotang_max_sol_all, rotang_max_sol_loc)
            rotang_max_gas_all = max(rotang_max_gas_all, rotang_max_gas_loc)
            stretchinc_max_sol_all = max(stretchinc_max_sol_all, stretchinc_max_sol_loc)
            stretchinc_max_gas_all = max(stretchinc_max_gas_all, stretchinc_max_gas_loc)

            cond_c066mod_max_sol_all = max(cond_c066mod_max_sol_all, cond_c066mod_max_sol_loc)
            cond_c066mod_max_gas_all = max(cond_c066mod_max_gas_all, cond_c066mod_max_gas_loc)
            cond_cgas_max_all = max(cond_cgas_max_all, cond_cgas_max_loc)
            !$OMP END CRITICAL(update_diagnostics_combine)

            !$OMP END PARALLEL

            diag%detF_min = detF_min_all
            diag%detF_min_sol = detF_min_sol_all
            diag%detF_min_gas = detF_min_gas_all
            diag%volinc_max_sol = volinc_max_sol_all
            diag%volinc_max_gas = volinc_max_gas_all
            diag%stretchinc_max_sol = stretchinc_max_sol_all
            diag%stretchinc_max_gas = stretchinc_max_gas_all
            diag%rotang_max_sol = rotang_max_sol_all
            diag%rotang_max_gas = rotang_max_gas_all
            diag%imin_sol = imin_sol
            diag%imin_gas = imin_gas
            diag%cond_c066mod_max_sol = cond_c066mod_max_sol_all
            diag%cond_c066mod_max_gas = cond_c066mod_max_gas_all
            diag%cond_cgas_max = cond_cgas_max_all

            ! Convert sentinel minima to 0.0 if a phase is absent.
            if (diag%detF_min >= huge(1.0d0) * 0.5d0) diag%detF_min = 0.0d0
            if (diag%detF_min_sol >= huge(1.0d0) * 0.5d0) diag%detF_min_sol = 0.0d0
            if (diag%detF_min_gas >= huge(1.0d0) * 0.5d0) diag%detF_min_gas = 0.0d0

            if (sum_w_all > 0.0d0) then
                diag%dL_step_rms = sqrt(sum_rel2_all / sum_w_all)
            else
                diag%dL_step_rms = 0.0d0
            end if
            diag%dL_step_max = max_rel_all

            if (sum_w_sol > 0.0d0) then
                diag%dL_step_rms_sol = sqrt(sum_rel2_sol / sum_w_sol)
            else
                diag%dL_step_rms_sol = 0.0d0
            end if
            diag%dL_step_max_sol = max_rel_sol

            if (sum_w_gas > 0.0d0) then
                diag%dL_step_rms_gas = sqrt(sum_rel2_gas / sum_w_gas)
            else
                diag%dL_step_rms_gas = 0.0d0
            end if
            diag%dL_step_max_gas = max_rel_gas

            diag%dvelgradavg_mag = sqrt(sum(micro%dvelgradavg**2))

            if (.not. do_scan_cond) then
                if (diag%imin_sol(1) > 0) then
                    diag%cond_c066mod_max_sol = cond1_est_6x6(micro%voxel(diag%imin_sol(1), diag%imin_sol(2), diag%imin_sol(3))%c066mod)
                end if
                if (diag%imin_gas(1) > 0) then
                    diag%cond_c066mod_max_gas = cond1_est_6x6(micro%voxel(diag%imin_gas(1), diag%imin_gas(2), diag%imin_gas(3))%c066mod)
                    diag%cond_cgas_max = cond1_est_6x6(micro%voxel(diag%imin_gas(1), diag%imin_gas(2), diag%imin_gas(3))%cgas)
                end if
            end if

            ! Step-to-step relative changes (and update prev).
            diag%detF_min_ratio = safe_ratio(diag%detF_min, diag%detF_min_prev)
            diag%detF_min_relstep = safe_rel_change(diag%detF_min, diag%detF_min_prev)
            diag%detF_min_prev = diag%detF_min

            diag%detF_min_sol_ratio = safe_ratio(diag%detF_min_sol, diag%detF_min_sol_prev)
            diag%detF_min_sol_relstep = safe_rel_change(diag%detF_min_sol, diag%detF_min_sol_prev)
            diag%detF_min_sol_prev = diag%detF_min_sol

            diag%detF_min_gas_ratio = safe_ratio(diag%detF_min_gas, diag%detF_min_gas_prev)
            diag%detF_min_gas_relstep = safe_rel_change(diag%detF_min_gas, diag%detF_min_gas_prev)
            diag%detF_min_gas_prev = diag%detF_min_gas

            diag%volinc_max_sol_relstep = safe_rel_change(diag%volinc_max_sol, diag%volinc_max_sol_prev)
            diag%volinc_max_sol_prev = diag%volinc_max_sol
            diag%volinc_max_gas_relstep = safe_rel_change(diag%volinc_max_gas, diag%volinc_max_gas_prev)
            diag%volinc_max_gas_prev = diag%volinc_max_gas

            diag%stretchinc_max_sol_relstep = safe_rel_change(diag%stretchinc_max_sol, diag%stretchinc_max_sol_prev)
            diag%stretchinc_max_sol_prev = diag%stretchinc_max_sol
            diag%stretchinc_max_gas_relstep = safe_rel_change(diag%stretchinc_max_gas, diag%stretchinc_max_gas_prev)
            diag%stretchinc_max_gas_prev = diag%stretchinc_max_gas
            diag%rotang_max_sol_prev = diag%rotang_max_sol
            diag%rotang_max_gas_prev = diag%rotang_max_gas
        end if

    contains

        pure double precision function one_norm_6x6(A) result(anorm)
            double precision, intent(in) :: A(6, 6)
            double precision :: colsum
            integer :: i6, j6
            anorm = 0.0d0
            do j6 = 1, 6
                colsum = 0.0d0
                do i6 = 1, 6
                    colsum = colsum + abs(A(i6, j6))
                end do
                anorm = max(anorm, colsum)
            end do
        end function one_norm_6x6

        double precision function cond1_est_6x6(Ain) result(cond)
            ! 1-norm condition number estimate via LAPACK dgetrf + dgecon.
            double precision, intent(in) :: Ain(6, 6)
            double precision :: A(6, 6), work_local(24), anorm, rcond
            integer :: ipiv(6), iwork(6), info_local

            A = Ain
            anorm = one_norm_6x6(Ain)

            call dgetrf(6, 6, A, 6, ipiv, info_local)
            if (info_local /= 0) then
                cond = huge(1.0d0)
                return
            end if

            call dgecon('1', 6, A, 6, anorm, rcond, work_local, iwork, info_local)
            if (info_local /= 0 .or. rcond <= 0.0d0) then
                cond = huge(1.0d0)
            else
                cond = 1.0d0 / rcond
            end if
        end function cond1_est_6x6

        pure double precision function safe_rel_change(x, x_prev) result(rel_change)
            double precision, intent(in) :: x, x_prev
            if (x_prev <= tiny(1.0d0) .or. x_prev >= huge(1.0d0) * 0.5d0) then
                rel_change = 0.0d0
            else
                rel_change = abs(x - x_prev) / x_prev
            end if
        end function safe_rel_change

        pure double precision function safe_ratio(x, x_prev) result(ratio)
            double precision, intent(in) :: x, x_prev
            if (x_prev <= tiny(1.0d0) .or. x_prev >= huge(1.0d0) * 0.5d0) then
                ratio = 1.0d0
            else
                ratio = x / x_prev
            end if
        end function safe_ratio

        pure double precision function determinant33(A) result(detA)
            double precision, intent(in) :: A(3, 3)
            detA = A(1, 1) * (A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2)) &
                   - A(1, 2) * (A(2, 1) * A(3, 3) - A(2, 3) * A(3, 1)) &
                   + A(1, 3) * (A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1))
        end function determinant33

    end subroutine update_diagnostics

    subroutine write_diagnostics(self, output_dir, iter_local, imacroloop_local, ipr_local, sim_time, wall_time, tdot, init)
        implicit none

        class(diag_type), intent(in) :: self
        character(len=*), intent(in) :: output_dir
        integer, intent(in) :: iter_local, imacroloop_local, ipr_local
        double precision, intent(in) :: sim_time, wall_time, tdot
        logical, intent(in), optional :: init

        character(len=512) :: csv_path
        character(len=2048) :: header
        integer :: u, ios
        logical :: do_init, exists
        integer :: nout

        nout = len_trim(output_dir)
        if (nout == 0) stop 'write_diagnostics: output_dir is empty'

        if (output_dir(nout:nout) == '/') then
            csv_path = trim(output_dir)//'diagnostics.csv'
        else
            csv_path = trim(output_dir)//'/diagnostics.csv'
        end if

        do_init = .false.
        if (present(init)) do_init = init

        header = 'iter,macrostep,ipr,sim_time,wall_time,tdot,'// &
                 'vm_max_inc,'// &
                 'tdotref,'// &
                 'erre,errs,errsbc,'// &
                 's33_target,s33,s33_sol,'// &
                 'detF_min_sol,detF_min_gas,detF_min_sol_rel,detF_min_gas_rel,'// &
                 'volinc_max_sol,volinc_max_gas,volinc_max_sol_rel,volinc_max_gas_rel,'// &
                 'stretchinc_max_sol,stretchinc_max_gas,stretchinc_max_sol_rel,stretchinc_max_gas_rel,'// &
                 'rotang_max_sol,rotang_max_gas,'// &
                 'dL_step_rms_sol,dL_step_rms_gas,dL_step_rms,dL_step_max_sol,dL_step_max_gas,dL_step_max,'// &
                 'dvelgradavg_mag,cond_c066mod_max_sol,cond_c066mod_max_gas,cond_cgas_max,'// &
                 'imin_sol_i,imin_sol_j,imin_sol_k,imin_gas_i,imin_gas_j,imin_gas_k'

        if (do_init) then
            open (newunit=u, file=trim(csv_path), status='replace', action='write', iostat=ios)
            if (ios /= 0) stop 'write_diagnostics: failed to open diagnostics.csv for writing'
            write (u, '(A)') trim(header)
            close (u)
        end if

        inquire (file=trim(csv_path), exist=exists)
        if (.not. exists) then
            open (newunit=u, file=trim(csv_path), status='replace', action='write', iostat=ios)
            if (ios /= 0) stop 'write_diagnostics: failed to create diagnostics.csv'
            write (u, '(A)') trim(header)
            close (u)
        end if

        open (newunit=u, file=trim(csv_path), status='old', position='append', action='write', iostat=ios)
        if (ios /= 0) stop 'write_diagnostics: failed to open diagnostics.csv for append'

        write (u, '(*(g0,:,","))') iter_local, imacroloop_local, ipr_local, sim_time, wall_time, tdot, &
            self%devm_inc_max, &
            self%tdotref, &
            self%erre, self%errs, self%errsbc, &
            self%s33_target, self%s33, self%s33_sol, &
            self%detF_min_sol, self%detF_min_gas, self%detF_min_sol_relstep, self%detF_min_gas_relstep, &
            self%volinc_max_sol, self%volinc_max_gas, self%volinc_max_sol_relstep, self%volinc_max_gas_relstep, &
            self%stretchinc_max_sol, self%stretchinc_max_gas, self%stretchinc_max_sol_relstep, self%stretchinc_max_gas_relstep, &
            self%rotang_max_sol, self%rotang_max_gas, &
            self%dL_step_rms_sol, self%dL_step_rms_gas, self%dL_step_rms, &
            self%dL_step_max_sol, self%dL_step_max_gas, self%dL_step_max, &
            self%dvelgradavg_mag, self%cond_c066mod_max_sol, self%cond_c066mod_max_gas, self%cond_cgas_max, &
            self%imin_sol(1), self%imin_sol(2), self%imin_sol(3), self%imin_gas(1), self%imin_gas(2), self%imin_gas(3)

        flush (u)
        close (u)
    end subroutine write_diagnostics

end
