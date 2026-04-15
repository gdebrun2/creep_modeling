# Mimosa Codebase Reference

## Table of Contents

- [1. Scope, Sources, and Reconciliation Rules](#1-scope-sources-and-reconciliation-rules)
- [2. Repository Architecture](#2-repository-architecture)
- [3. Python Orchestration and the Input Contract](#3-python-orchestration-and-the-input-contract)
  - [3.1 Python Driver Entry Point](#31-python-driver-entry-point)
  - [3.2 FFT Input Writer and Simulation Templates](#32-fft-input-writer-and-simulation-templates)
  - [3.3 Material-file writers](#33-material-file-writers)
  - [3.4 Microstructure metadata and load conversion](#34-microstructure-metadata-and-load-conversion)
  - [3.5 Run registration and output layout](#35-run-registration-and-output-layout)
- [4. Fortran Data Model](#4-fortran-data-model)
  - [4.1 Voxel Type](#41-voxel-type)
  - [4.2 Micro Type](#42-micro-type)
  - [4.3 Phase Type](#43-phase-type)
  - [4.4 Process Type](#44-process-type)
  - [4.5 Diagnostics Type](#45-diagnostics-type)
- [5. Input Parsing and Phase Initialization](#5-input-parsing-and-phase-initialization)
  - [5.1 Load Input](#51-load-input)
  - [5.2 Current igas Semantics in the Parser](#52-current-igas-semantics-in-the-parser)
  - [5.3 Consequence for gas material files](#53-consequence-for-gas-material-files)
  - [5.4 Crystal plasticity file reader](#54-crystal-plasticity-file-reader)
  - [5.5 Hall-Petch directional offset during initialization](#55-hall-petch-directional-offset-during-initialization)
- [6. Kinematic and Constitutive Notation](#6-kinematic-and-constitutive-notation)
  - [6.1 Current-configuration averaging](#61-current-configuration-averaging)
  - [6.2 PK1-like fluctuation field](#62-pk1-like-fluctuation-field)
- [7. Top-Level Solver Algorithm](#7-top-level-solver-algorithm)
  - [7.1 Command-line controls](#71-command-line-controls)
  - [7.2 One macrostep in pseudocode](#72-one-macrostep-in-pseudocode)
- [8. Spectral Equilibrium Solve](#8-spectral-equilibrium-solve)
  - [8.1 Green operator construction](#81-green-operator-construction)
  - [8.2 Push-forwarded stiffness and modified AL operator](#82-push-forwarded-stiffness-and-modified-al-operator)
  - [8.3 PK1 convolution and correction back to current configuration](#83-pk1-convolution-and-correction-back-to-current-configuration)
  - [8.4 Adaptive Mixing of Velgrad](#84-adaptive-mixing-of-velgrad)
- [9. Mixed Macroscopic Boundary Conditions](#9-mixed-macroscopic-boundary-conditions)
  - [9.1 What stresses are averaged](#91-what-stresses-are-averaged)
  - [9.2 Stress boundary-condition mismatch](#92-stress-boundary-condition-mismatch)
  - [9.3 Macro compliance relation](#93-macro-compliance-relation)
  - [9.4 State 6x6 EVPSC](#94-state-6x6-evpsc)
  - [9.5 Average-velocity-gradient correction](#95-average-velocity-gradient-correction)
- [10. Solid Crystal Plasticity Model](#10-solid-crystal-plasticity-model)
  - [10.1 Slip geometry and Schmid tensors](#101-slip-geometry-and-schmid-tensors)
  - [10.2 Resolved shear stress and sign-dependent resistance](#102-resolved-shear-stress-and-sign-dependent-resistance)
  - [10.3 Exact slip law used by the code](#103-exact-slip-law-used-by-the-code)
  - [10.4 Plastic strain-rate tensor](#104-plastic-strain-rate-tensor)
  - [10.5 Consistent derivative of the slip law](#105-consistent-derivative-of-the-slip-law)
  - [10.6 Elastic strain-rate contribution](#106-elastic-strain-rate-contribution)
  - [10.7 Local nonlinear constitutive solve](#107-local-nonlinear-constitutive-solve)
  - [10.8 Hardening law actually used by the source](#108-hardening-law-actually-used-by-the-source)
  - [10.9 When hardening is applied](#109-when-hardening-is-applied)
  - [10.10 Macro tangent update](#1010-macro-tangent-update)
- [11. Current Gas Model](#11-current-gas-model)
  - [11.1 What the current source actually implements](#111-what-the-current-source-actually-implements)
  - [11.2 Gas compliance operator](#112-gas-compliance-operator)
  - [11.3 Gas Stress Update in Solve Res](#113-gas-stress-update-in-solve-res)
  - [11.4 Why gas material files do not matter today](#114-why-gas-material-files-do-not-matter-today)
  - [11.5 Historical and Planned igas 2 Behavior](#115-historical-and-planned-igas-2-behavior)
- [12. Accepted-Step Finite-Strain Updates](#12-accepted-step-finite-strain-updates)
  - [12.1 Immediate state accumulation](#121-immediate-state-accumulation)
  - [12.2 Plastic deformation gradient](#122-plastic-deformation-gradient)
  - [12.3 Total deformation gradient](#123-total-deformation-gradient)
  - [12.4 Elastic deformation gradient](#124-elastic-deformation-gradient)
  - [12.5 Current elastic stiffness](#125-current-elastic-stiffness)
  - [12.6 Objective rotation and stress push-forward](#126-objective-rotation-and-stress-push-forward)
- [13. Timestep Control, Substepping, and Current Iteration Safeguards](#13-timestep-control-substepping-and-current-iteration-safeguards)
  - [13.1 Strain-controlled timestep rescaling](#131-strain-controlled-timestep-rescaling)
  - [13.2 Kinematic substepping and rollback](#132-kinematic-substepping-and-rollback)
  - [13.3 Step rejection path](#133-step-rejection-path)
- [14. Diagnostics and Output Files](#14-diagnostics-and-output-files)
  - [14.1 Results CSV](#141-results-csv)
  - [14.2 Diagnostics CSV](#142-diagnostics-csv)
  - [14.3 What the diagnostics mean mathematically](#143-what-the-diagnostics-mean-mathematically)
  - [14.4 Python-side loaders](#144-python-side-loaders)
- [15. Calibration, Loss Functions, and Postprocessing](#15-calibration-loss-functions-and-postprocessing)
  - [15.1 What the docs recommend conceptually](#151-what-the-docs-recommend-conceptually)
  - [15.2 What the current calibration code actually computes](#152-what-the-current-calibration-code-actually-computes)
  - [15.3 Creep objective](#153-creep-objective)
  - [15.4 SRJ objective](#154-srj-objective)
  - [15.5 Combined loss](#155-combined-loss)
  - [15.6 Postprocessing](#156-postprocessing)
- [16. Material Parameters: Meaning, Priors, and Current Defaults](#16-material-parameters-meaning-priors-and-current-defaults)
  - [16.1 Parameter meanings](#161-parameter-meanings)
  - [16.2 Documentation priors versus active code defaults](#162-documentation-priors-versus-active-code-defaults)
  - [16.3 Gas parameter defaults](#163-gas-parameter-defaults)
- [17. The 0D Bulk Surrogate Note](#17-the-0d-bulk-surrogate-note)
- [18. Code Map by Algorithmic Responsibility](#18-code-map-by-algorithmic-responsibility)
- [19. Baseline Summary](#19-baseline-summary)


## 2. Repository Architecture

At a high level, `mimosa` is an end-to-end workflow around a large-strain elastic-viscoplastic FFT crystal plasticity solver.

The repository splits into two cooperating layers:

- `lsevpfft/src`: the Fortran solver, which owns the constitutive update, spectral equilibrium solve, finite-strain kinematics, timestep adaptation, diagnostics, and CSV/VTK field output.
- `utils`: the Python layer, which writes input files, allocates run directories, launches the solver, parses outputs, performs calibration, computes roughness/slip postprocessing from VTK fields, and generates standardized plots.

The effective execution chain is:

1. Python chooses a simulation case and parameter set.
2. Python writes `cuel.sx`, `cupl2.sx`, optional gas material files, and `fft.in`.
3. The Fortran solver reads `fft.in`, the phase files, and the microstructure.
4. The Fortran solver advances the microstructure in time and writes `results.csv`, `diagnostics.csv`, and optional VTK fields.
5. Python reloads those outputs into typed arrays and generates plots.

## 3. Python Orchestration and the Input Contract

### 3.1 Python Driver Entry Point

The main orchestration entry point is `utils/run.py`.

Its responsibilities are:

- allocate a per-run directory under `results/<sim_case>/<run_id>/`;
- load an optional calibration row from `results/calibrations/calibration_index.csv`;
- write solid material files using `write_cuel` and `write_cupl2`;
- optionally write gas material files when `igas == 2`;
- generate `fft.in` through `utils/write_fft.py`;
- invoke `LS-EVPFFT` with command-line flags such as `--nthreads`, `--in`, `--out`, `--fftw_save`, `--mix`, `--substep`, and `--modal`;
- always attempt postprocessing in a `finally` block.

The core configuration object is `RunConfig`, which inherits from `FftConfig` and adds run-directory, solver, and postprocessing metadata.

### 3.2 FFT Input Writer and Simulation Templates

`utils/write_fft.py` dispatches to `utils/sim_cases/*.py` based on `cfg.sim_case`.

This file is the Python-side contract writer for `fft.in`:

- `sim_case` selects which template is used;
- `igas`, `complgas`, `load`, `microstructure`, `cuel`, `cupl2`, `cuel_gas`, and `cupl2_gas` are inserted into the template text;
- `write_text` writes the finished file with a trailing newline.

The templates encode the structure expected by `IO_functions.f90::load_input`, so the Python and Fortran layers are tightly coupled through file layout, not through an API object.

### 3.3 Material-file writers

`utils/io_utils.py` writes the two material-file types used by the solver:

1. `cuel.sx`

   - `write_cuel(..., iso=0)` writes a full cubic 6x6 stiffness matrix with entries `C11`, `C12`, `C44`.
   - `write_cuel(..., iso=1)` writes isotropic `Young, nu`.

2. `cupl2.sx`

   `write_cupl2` writes a single FCC `{111}<110>` slip mode with:

   - `nrsx`
   - `gamd0x`
   - `tau0xf`
   - `tau0xb`
   - `tau1x`
   - `thet0`
   - `thet1`
   - `hselfx`
   - `hlatex`

The Python writer is consistent with the current Fortran crystal-file reader in `IO_functions.f90::data_crystal`.

### 3.4 Microstructure metadata and load conversion

`utils/data_utils.py::MicrostructureInfo.load` parses the microstructure file as:

- columns 3,4,5: voxel coordinates `x,y,z`
- column 6: grain id
- column 7: phase id

and infers:

- `nx, ny, nz`
- `nphase`
- `solid_frac`
- `gas_frac`

`utils/calc_utils.py::convert_load` then applies the repository convention:

$$
\sigma_{33}^{\text{written to fft.in}} \approx f_s \, \sigma_{33}^{\text{desired solid average}}
$$

where `f_s` is the reference solid volume fraction.

The intent of this helper is to compensate for highly compliant gas padding. It assumes the imposed macroscopic stress in `fft.in` is interpreted as a domain average while gas carries negligible traction.

That assumption is only approximately consistent with the current Fortran implementation:

- during the PK1 fluctuation build, `micro%sgavg` is first computed over all voxels;
- during boundary-condition correction, `calc_velgradavg_corr` rebuilds `micro%sgavg` from solid contributions only, without dividing by `micro%wphcsol`.

So the code is effectively using a domain-average stress with gas omitted, not a strict all-voxel average and not a normalized solid-only average either.

### 3.5 Run registration and output layout

Every run is registered in `results/run_index.csv` through `utils/io_utils.py::register_run`, which records:

- `run_id`
- `sim_case`
- `calib_id`
- `igas`
- `load`
- `complgas`
- `substep`
- `mix`
- `modal`
- `fft`
- `nthreads`
- `run_dir`
- `microstructure`
- `notes`

The expected run layout is:

- `<run_dir>/fft.in`
- `<run_dir>/cuel.sx`
- `<run_dir>/cupl2.sx`
- `<run_dir>/cuel_gas.sx` when written
- `<run_dir>/cupl2_gas.sx` when written
- `<run_dir>/output.txt`
- `<run_dir>/evpfft_outputs/results.csv`
- `<run_dir>/evpfft_outputs/diagnostics.csv`
- `<run_dir>/evpfft_outputs/microstr_cell_fields_stp_*.vtk`
- `<run_dir>/plots/`

## 4. Fortran Data Model

The central derived types are declared in `lsevpfft/src/types.f90`.

### 4.1 Voxel Type

Each voxel stores:

- kinematics: `velgrad`, `defgrad`, `defgradinv`, `defgrade`, `defgradp`, `detF`, `defgradinc`, `defgradeinc`
- stress and strain measures: `sg`, `sgt`, `edotp`, `ept`, `disgradsym`, `epvm`
- crystal plasticity state: `gamdot`, `crss`, `trialtau`, `kin`, `sch`, `ag`, `gacumgr`
- local constitutive operators: `cg66`, `sg66`, `c066mod`, `c066modGoperr066mod`, `cgas`
- spectral quantities: `sgPK1`, `dvelgradref`
- current volume weight: `wgtc`
- gas flag: `gas`
- snapshot fields used for rollback and diagnostics

### 4.2 Micro Type

The full microstructure object stores:

- macro tangents: `c066`, `s066`, `c0`, `s0`, `de_ds`
- Green operator data: `Goper`, `Goperr0`
- current-configuration averages: `velgradavg`, `velgradavgsol`, `velgradavggas`, `sgavg`, `sgavgsol`
- finite-strain averages: `defgradavg`, `defgradinvavgc`, `defgradinvavgcs`, `defgradinvavgcg`
- current solid/gas volume fractions: `wphcsol`, `wphcgas`
- output averages: `epavg`, `eavg`, `eelavg`, `edotpavg`, `sgsolavg`, `evmavg`, `evmpavg`, `dvmavg`, `dvmpavg`, `svmavg`, `lvmavg`
- nodal coordinates `xnode`
- the voxel array itself

### 4.3 Phase Type

Each phase stores:

- `igas` as a logical phase flag in the current solver
- slip-system normals `dnca`
- slip directions `dbca`
- Schmid tensors in crystal basis `schca`
- hardening parameters `tau(:,1:3)` and `thet(:,0:1)`
- rate sensitivity `nrs`
- reference slip rates `gamd0`
- latent/self hardening matrix `hard`
- crystal elasticity tensor `C`

### 4.4 Process Type

Each loading process stores:

- masks for prescribed velocity-gradient entries: `iudot`
- masks for prescribed symmetric strain-rate entries: `idsim`
- masks for prescribed stress entries: `iscau`
- target tensors: `udot`, `dsim`, `tomtot`, `scauchy`
- timestep control data: `tdot`, `tdotref`, `tdotmin`, `devmmx`, `devmapp`, `ictrl`
- nonlinear solve controls: `error`, `erroral`, `itmax`
- update flags: `iupdate`, `iuphard`, `igamma`, `iwfields`, `iwstep`

### 4.5 Diagnostics Type

The diagnostics object stores one macrostep record:

- field residual norms `erre`, `errs`, `errsbc`
- target and achieved stresses `s33_target`, `s33`, `s33_sol`
- minimum `detF` values for solid and gas
- maximum volumetric, stretch, and rotation increments
- RMS and max per-step `L` changes
- BC-correction magnitude `dvelgradavg_mag`
- condition-number estimates for `c066mod` and `cgas`
- voxel locations where extreme minima occur

## 5. Input Parsing and Phase Initialization

### 5.1 Load Input

`lsevpfft/src/IO_functions.f90::load_input` is the authoritative parser for `fft.in`.

Its top-level duties are:

1. read `nph`, `npts1`, `npts2`, `npts3`;
2. initialize FFTW and allocate `micro`;
3. read voxel spacings `delt`;
4. read the microstructure path;
5. read one phase block per phase;
6. read restart settings and `complgas`;
7. read optional initial deformation-gradient data;
8. initialize per-voxel CRSS, trial CRSS, and Hall-Petch offsets;
9. read loading processes and enforce mask consistency;
10. compute reference solid/gas phase fractions.

### 5.2 Current igas Semantics in the Parser

This is the single most important phase-model reconciliation:

```fortran
read (11, *) ii
if (ii == 0) then
    props%phase(iph)%igas = .false.
else
    props%phase(iph)%igas = .true.
end if
```

So, in the current solver:

- `igas = 0` means solid
- `igas = 1` means gas
- `igas = 2` also means gas
- any other nonzero integer would also mean gas

Since `phase_type` only stores a logical `igas`, all nonzero modes collapse to one branch before the constitutive code is reached.

### 5.3 Consequence for gas material files

Still inside `load_input`, the crystal plasticity and elasticity files are only read when:

```fortran
if (.not. props%phase(iph)%igas) then
    call data_crystal(...)
    call data_crystal_elast(...)
end if
```

Therefore:

- solid phase files are read normally;
- gas phase files are ignored for any nonzero `igas`.

This means that the present Fortran solver never reads `cuel_gas.sx` or `cupl2_gas.sx` for a gas-like phase, even though the Python layer still writes them for `igas == 2`.

### 5.4 Crystal plasticity file reader

`IO_functions.f90::data_crystal` reads the `cupl2.sx` slip/hardening file.

For each active deformation mode it reads:

- `modex`
- `nsmx`
- `nrsx`
- `gamd0x`
- `twshx`
- `isectwx`
- `tau0xf`
- `tau0xb`
- `tau1x`
- `thet0x`
- `thet1x`
- `hpfacx`
- `hselfx`
- `hlatex`

The routine then:

1. expands mode-level values to system-level arrays;
2. builds normalized slip plane normals `dnca`;
3. builds normalized slip directions `dbca`;
4. constructs the symmetric Schmid tensor

$$
S^\alpha_{\text{crystal}} = \frac{1}{2}\left(n^\alpha \otimes b^\alpha + b^\alpha \otimes n^\alpha\right)
$$

and stores it in the five-component basis `schca`;
5. constructs the latent hardening matrix `hard` by filling every pair of systems in modes `(m,n)` with `hlatex(m,n)` and then replacing diagonal self terms with `hselfx(m)`.

### 5.5 Hall-Petch directional offset during initialization

After `data_grain` assigns orientations and phases, `load_input` initializes voxel CRSS values as:

$$
\tau_{\alpha}^{\pm} \leftarrow \tau_{\alpha,\text{file}}^{\pm} + h_{P,\alpha}\sqrt{|\cos \alpha|}
$$

where in code:

- the slip direction is rotated to sample axes as `dbsa = AG * dbca`;
- `cosalpha = abs(dbsa(2))`, i.e. the absolute projection onto the sample `y` direction.

This is an anisotropic Hall-Petch-like shift baked into the initialization step.

## 6. Kinematic and Constitutive Notation

The solver is a finite-strain FFT crystal plasticity code. The most useful notation for reading the implementation is:

- reference position: $X$
- current position: $x$
- deformation map: $x = \varphi(X,t)$
- deformation gradient:

$$
F = \nabla_X \varphi
$$

- Jacobian:

$$
J = \det F
$$

- velocity gradient:

$$
L = \nabla_x v
$$

- symmetric and antisymmetric parts:

$$
D = \operatorname{sym}(L) = \frac{1}{2}(L + L^T), \qquad
W = \operatorname{skw}(L) = \frac{1}{2}(L - L^T)
$$

- Cauchy stress: $\sigma$

The code uses two reduced tensor bases:

- a five-component deviatoric basis for Schmid tensors and resolved-shear calculations;
- a six-component symmetric basis for stresses, strain rates, tangents, and gas operators.

The basis/Voigt conversions are handled by `tensor_functions.f90::chg_basis` and `voigt_vpsc`.

### 6.1 Current-configuration averaging

The current voxel weight is updated as:

$$
w_i^{c} = \frac{w_i^{0} J_i}{\langle J \rangle_0}
$$

where `props%wgt` is the uniform reference weight and `detFavg` is the reference-weighted average of `J`.

These weights are then used for most macro averages.

### 6.2 PK1-like fluctuation field

Before the spectral correction, the solver forms a PK1-like stress fluctuation:

$$
P_i' = (\sigma_i - \langle \sigma \rangle) F_i^{-T} J_i
$$

implemented by `various_functions.f90::CauchyToPK1` after subtracting the current iteration's average stress.

This is the field convolved with the Green operator in Fourier space.

## 7. Top-Level Solver Algorithm

The executable starts in `lsevpfft/src/LS-EVPFFT.f90`.

The actual control flow is:

1. parse CLI options;
2. initialize FFTW output directories and clear stale VTK/`.out` files;
3. `load_input(micro, props, nthreads, fft_input_file)`;
4. `initial_guess(micro, props)`;
5. write initial VTK and initialize `results.csv` and `diagnostics.csv`;
6. loop over loading processes `ipr`;
7. inside each process, loop over macrosteps `imicro`;
8. inside each macrostep, run the equilibrium iteration until either convergence or accepted cutback;
9. after the step is accepted, update history variables and finite-strain state, write outputs, possibly ramp `tdot`, and rebuild the macro tangent.

### 7.1 Command-line controls

The solver accepts:

- `--in`
- `--out`
- `--nthreads`
- `--fftw_save`
- `--modal`
- `--substep`
- `--mix`

These correspond directly to the values passed by `utils/run.py`.

### 7.2 One macrostep in pseudocode

For one process and one macrostep, the current code does:

```text
take_snapshot
if first step or iupdate == 1:
    update_schmid
form_G
calc_c066mod_Goperr066mod

initialize iter = 0
while not converged and not accept_step:
    if iter == 0:
        calc_gas_compl
        dvelgradavg = 0
        mix = 1
        initialize erre, errs, errsbc

    compute iteration-average stress for PK1 fluctuation
    build sgPK1 = CauchyToPK1(sg - sgavg, Finv, detF)
    convolution_with_Goper
    correction_to_current
    adaptively mix velgrad
    solve_res
    calc_velgradavg_corr
    normalize errors
    update tdot if ictrl == 0
    test convergence
    if ictrl == -1:
        maybe rollback and split the step using kinematic_substep

accept step
accumulate ept, sgt, disgradsym, epvm
if iupdate == 1:
    update_plastic_defgrad
    update_grid_velgrad_node
    update_defgrad
    update_elstic_defgrad
    update_el_stiff
    update_tensor_config
if iuphard == 1:
    update_hardening
recompute averages
maybe write VTK
update_diagnostics(finalize=.true., scan_cond=.true.)
write diagnostics.csv and results.csv
maybe ramp tdot upward
calc_c0
```

The relevant implementation blocks are:

- `LS-EVPFFT.f90:169-509`
- `various_functions.f90`
- `types.f90::update_diagnostics`
- `IO_functions.f90::write_results`

## 8. Spectral Equilibrium Solve

### 8.1 Green operator construction

The Green operator is built in `various_functions.f90::form_G`.

At each Fourier mode, the code constructs a derivative tensor $D^{\text{dft}}$ and an acoustic tensor:

$$
A_{ik} = C^0_{ijkl} D^{\text{dft}}_{jl}
$$

then inverts $A$ and forms:

$$
G_{pqij} = -A^{-1}_{pi} D^{\text{dft}}_{qj}
$$

The macro reference stiffness `micro%c0` entering this expression is updated once per accepted macrostep by `calc_c0`.

#### `igamma = 0`: continuous derivative

For the continuous Fourier derivative:

$$
D^{\text{dft}}_{ij} = \xi_i \xi_j
$$

where $\xi$ is the spatial frequency vector returned by `spatial_freq`.

#### `igamma = 1`: discrete central-difference-like derivative

For the discrete form:

$$
D_{11} = 2(\cos \xi_x - 1), \quad
D_{22} = 2(\cos \xi_y - 1), \quad
D_{33} = 2(\cos \xi_z - 1)
$$

and the off-diagonals are built from products of sines, e.g.

$$
D_{21} = -\sin \xi_x \sin \xi_y
$$

with symmetry used to fill the remaining entries.

The code then performs up to four fixed-point smoothing iterations on the Green operator itself:

$$
G \leftarrow -G : C^0 : G_{\text{tmp}}
$$

with full mixing `xmix = 1`.

#### `igamma = 2`: modified tangent operator

For the modified operator, the code forms a complex wavevector:

$$
k^{\text{mod}}_i =
\frac{i}{4}\tan\left(\frac{\xi_i}{2}\right)
\prod_{m=1}^{3}\left(1 + e^{i\xi_m}\right)
$$

implemented exactly as the product-and-tangent form in `form_G`, then constructs

$$
D^{\text{dft}}_{ij} = \Re\left(k^{\text{mod}}_i \overline{k^{\text{mod}}_j}\right)
$$

after checking that the imaginary part is numerically negligible.

#### Real-space origin operator

After the Fourier-space tensor `micro%Goper` is built, the code inverse-transforms it and stores only the origin value:

$$
G^{r}_0 = \mathcal{F}^{-1}[G](0)
$$

as `micro%Goperr0`.

That origin term is later used in the modified augmented-Lagrangian correction.

### 8.2 Push-forwarded stiffness and modified AL operator

`various_functions.f90::calc_c066mod_Goperr066mod` updates the per-voxel mapped operators.

The macro reference stiffness is pushed to the current configuration as:

$$
c^{\text{mod}}_{ijkl}
=
\frac{1}{J}
\sum_{m,n}
C^0_{i m k n} F_{j m} F_{l n}
$$

then symmetrized in $(i,j)$ and $(k,l)$ and converted to 6x6 form as `c066mod`.

When modified AL is enabled, the code also pushes forward the real-space origin Green operator:

$$
G^{r,\text{mod}}_{ijkl}
=
J
\sum_{m,n}
G^r_{0, i m k n} F^{-1}_{m j} F^{-1}_{n l}
$$

symmetrizes it, converts it to 6x6 form, and stores

$$
C^{\text{mod}} G^{r,\text{mod}}
$$

as `c066modGoperr066mod`.

### 8.3 PK1 convolution and correction back to current configuration

`convolution_with_Goper` performs:

$$
\Delta L^{\text{ref}} = G * P'
$$

by FFTing `sgPK1`, multiplying with `Goper`, and inverse FFTing back to `dvelgradref`.

`correction_to_current` then maps that fluctuation correction into the current configuration:

$$
\Delta L_i^{\text{fluct}} = \Delta L_i^{\text{ref}} F_i^{-1}
$$

The fluctuation average is

$$
\langle \Delta L^{\text{fluct}} \rangle = \sum_i w_i^c \Delta L_i^{\text{fluct}}
$$

and the reference-average correction needed to enforce the desired macro average is:

$$
\left(\frac{\partial \dot u}{\partial X}\right)_{\text{avg}}
=
\left(\Delta \bar L - \langle \Delta L^{\text{fluct}} \rangle \right)
\left\langle F^{-1} \right\rangle^{-1}
$$

where `micro%dvelgradavg` carries the BC correction from the previous equilibrium iteration.

Finally each voxel is updated as:

$$
L_i \leftarrow L_i +
\left(
\Delta L_i^{\text{ref}} +
\left(\frac{\partial \dot u}{\partial X}\right)_{\text{avg}}
\right) F_i^{-1}
$$

This alternating structure is important:

- first `correction_to_current` applies the previous macro correction;
- then `solve_res` updates local stress;
- then `calc_velgradavg_corr` computes the new macro correction for the next iteration.

### 8.4 Adaptive Mixing of Velgrad

After `correction_to_current`, the code computes per-iteration relative changes:

$$
\delta L_i = \|L_i - L_i^{\text{old}}\|_F, \qquad
\delta L_i^{\text{rel}} = \frac{\delta L_i}{\max(\|\bar L\|_{\text{vm}}, \text{tiny})}
$$

and forms weighted RMS values separately over solid and gas voxels. The global control scalar is:

$$
dL_{\text{iter}} = \max(dL_{\text{iter,sol}}, dL_{\text{iter,gas}})
$$

The relaxation factor is then:

$$
\text{mix} =
\max\left(
\text{mix}_{\min},
\min\left(1,\frac{dL_{\text{target}}}{\max(dL_{\text{iter}}, dL_{\text{target}})}\right)
\right)
$$

with current constants:

- `mix_min = 1e-8`
- `dL_target = 5e-2`

and the actual update is:

$$
L_i \leftarrow \text{mix}\,L_i + (1-\text{mix})\,L_i^{\text{old}}
$$

when mixing is enabled.

## 9. Mixed Macroscopic Boundary Conditions

Mixed boundary conditions are enforced in `various_functions.f90::calc_velgradavg_corr` and `state_6x6_evpsc`.

### 9.1 What stresses are averaged

The routine computes two related stresses:

1. normalized solid average

$$
\langle \sigma \rangle_{\text{solid}}
=
\frac{1}{w_{\text{solid}}^c}
\sum_{i \in \text{solid}} w_i^c \sigma_i
$$

stored as `micro%sgavgsol`;

2. unnormalized solid contribution to the domain average

$$
\sigma_{\text{dom,solid contrib}}
=
\sum_{i \in \text{solid}} w_i^c \sigma_i
$$

stored as `micro%sgavg`.

This second quantity is not the exact all-voxel average if gas stresses are nonzero. It is closer to "domain average with gas stress omitted".

This subtlety matters because:

- `convert_load` in Python assumes the imposed `scauchy33` is a domain-average quantity;
- `calc_velgradavg_corr` compares `props%process(ipr)%scauchy` against `micro%sgavg`, not against `micro%sgavgsol`.

### 9.2 Stress boundary-condition mismatch

The stress mismatch tensor is built componentwise with the `iscau` mask:

$$
\Delta \sigma^{\text{BC}}_{ij}
=
\mathbb{I}^{\sigma}_{ij}\left(\sigma^{\text{target}}_{ij} - \sigma^{\text{avg}}_{ij}\right)
$$

and the mismatch norm is:

$$
\text{errsbc} = \|\Delta \sigma^{\text{BC}}\|_F
$$

### 9.3 Macro compliance relation

The current macro averages are mapped to six-component form:

- `dbar = basis(sym(micro%velgradavg))`
- `sbar = basis(micro%sgavg)`

and the macro compliance-like operator is the current `micro%de_ds`.

The code forms:

$$
d_{\text{vp},0} = d_{\text{bar}} - M_{\text{evp}} s_{\text{bar}}
$$

where `M_evp` is `micro%de_ds` in the VPSC six-component convention.

### 9.4 State 6x6 EVPSC

`state_6x6_evpsc` solves the mixed stress/strain-rate correction system in VPSC Voigt notation.

If we define:

- $P_d = \operatorname{diag}(\texttt{idsim})$
- $P_s = \operatorname{diag}(\texttt{iscau})$
- $M = M_{\text{evp}}$ in the VPSC convention
- $d_0 = d_{\text{vp},0}$

then the routine constructs a masked linear system of the form

$$
\left(P_s - P_d M\right)x
=
-P_d\left(d_{\text{bar}} - d_0\right) + M P_s s_{\text{bar}}
$$

up to the shear-scaling factors `profac` used by VPSC notation.

The solved vector `x` is then used to update both:

- the imposed macroscopic symmetric strain-rate `dsim`
- the imposed macroscopic stress tensor `scauchy`

so that the stress-controlled and strain-controlled components remain complementary.

### 9.5 Average-velocity-gradient correction

After the mixed system is solved:

$$
\Delta \bar D = D^{\text{target}} - D^{\text{current}}
$$

and the antisymmetric part is constructed directly from `iudot` and `udot`.

The final macro correction stored for the next global iteration is:

$$
\Delta \bar L = \Delta \bar D + \Delta \bar W
$$

as `micro%dvelgradavg`.

## 10. Solid Crystal Plasticity Model

The solid constitutive update lives primarily in:

- `various_functions.f90::update_schmid`
- `various_functions.f90::plastic_strain_rate`
- `various_functions.f90::elastic_strain_rate`
- `various_functions.f90::solve_res`
- `various_functions.f90::harden`
- `various_functions.f90::calc_c0`

### 10.1 Slip geometry and Schmid tensors

For each solid voxel and slip system, `update_schmid` builds the unsymmetrized Schmid tensor:

$$
M^\alpha = b^\alpha \otimes n^\alpha
$$

using the crystal-basis vectors from the phase definition.

It then maps that tensor to the current elastic/lattice configuration using:

$$
A = F^e R
$$

where in code `R` is `voxel%ag`, and computes:

$$
\widetilde M^\alpha = A M^\alpha A^{-1}
$$

The symmetric part

$$
S^\alpha = \operatorname{sym}(\widetilde M^\alpha)
$$

is then converted to the five-component basis and stored in `voxel%sch(:,alpha)`.

### 10.2 Resolved shear stress and sign-dependent resistance

The resolved shear stress on system $\alpha$ is:

$$
\tau^\alpha = S^\alpha : \sigma
$$

implemented in the five-component basis as a simple dot product of `sch(:,is)` with the deviatoric part of `sg`.

The code supports forward/backward asymmetry by selecting:

$$
\tau_{\text{trial}}^{\alpha,\pm}
=
\begin{cases}
\tau_{\text{trial}}^{\alpha,2}, & \tau^\alpha - \chi^\alpha < 0 \\
\tau_{\text{trial}}^{\alpha,1}, & \tau^\alpha - \chi^\alpha \ge 0
\end{cases}
$$

where `chi` is `voxel%kin(is)`.

The normalized driving stress is:

$$
r^\alpha = \frac{\tau^\alpha - \chi^\alpha}{\tau_{\text{trial}}^{\alpha,\pm}}
$$

### 10.3 Exact slip law used by the code

The current implementation is the power law:

$$
\dot\gamma^\alpha
=
\dot\gamma_0^\alpha
\left|r^\alpha\right|^{n^\alpha-1}
r^\alpha
$$

This is not a thresholded Macaulay-bracket law. For any nonzero driving stress, the slip rate is nonzero, although it may be extremely small.

In code notation:

$$
\texttt{dum} = \dot\gamma_0^\alpha |r^\alpha|^{n^\alpha - 1}
$$

$$
\dot\gamma^\alpha = r^\alpha \,\texttt{dum}
$$

### 10.4 Plastic strain-rate tensor

The plastic strain rate is assembled as:

$$
\dot\varepsilon^p = \sum_\alpha \dot\gamma^\alpha S^\alpha
$$

In the five-component basis this is implemented by summing `sch(jj,is) * rss2`.

### 10.5 Consistent derivative of the slip law

The code also forms the derivative of plastic strain rate with respect to stress:

$$
\frac{\partial \dot\varepsilon^p}{\partial \sigma}
=
\sum_\alpha
\left(
\frac{n^\alpha \dot\gamma_0^\alpha |r^\alpha|^{n^\alpha - 1}}
{\tau_{\text{trial}}^{\alpha,\pm}}
\right)
S^\alpha \otimes S^\alpha
$$

This is exactly the `dedotp_dsg` tensor built in `plastic_strain_rate`.

### 10.6 Elastic strain-rate contribution

The elastic contribution is treated as a rate form based on the stress increment from the previous accepted macrostep:

$$
\dot\varepsilon^e
=
S^{e}\frac{\sigma - \sigma_t}{\Delta t}
$$

where `S^e` is the voxel elastic compliance matrix `sg66`.

The consistent derivative is:

$$
\frac{\partial \dot\varepsilon^e}{\partial \sigma}
=
\frac{S^e}{\Delta t}
$$

### 10.7 Local nonlinear constitutive solve

The solid constitutive update in `solve_res` has:

- an outer fixed-point loop on slip resistance (`itmaxaltau = 10`);
- an inner Newton loop on stress (`itmaxal = 100`).

Let:

- $\lambda$ be the current global-iteration stress guess (`voxel%sg` at entry to the local solve, expressed in the 6-component basis);
- $\sigma$ be the local unknown in the same basis;
- $D$ be the symmetric part of the current voxel velocity gradient in 6-component form.

The total strain-rate response is:

$$
\dot\varepsilon(\sigma)
=
\dot\varepsilon^p(\sigma) + \dot\varepsilon^e(\sigma) + \dot\varepsilon^{th}
$$

where the thermal/eigenstrain term is added only for the first macrostep of a thermally active process.

The local residual is:

$$
R(\sigma)
=
\sigma - \lambda
+ C^{\text{mod}}\left(\dot\varepsilon(\sigma) - D\right)
- \mathcal{A}_{\text{modAL}} C^{\text{mod}} G_0^{r,\text{mod}} (\sigma - \lambda)
$$

where $\mathcal{A}_{\text{modAL}}$ is 1 when modified AL is enabled and 0 otherwise.

The Jacobian is:

$$
J(\sigma)
=
I + C^{\text{mod}}\frac{\partial \dot\varepsilon}{\partial \sigma}
- \mathcal{A}_{\text{modAL}} C^{\text{mod}} G_0^{r,\text{mod}}
$$

The Newton update solves:

$$
J \,\Delta \sigma = R
$$

and then applies:

$$
\sigma \leftarrow \sigma - \Delta \sigma
$$

The code also uses a simple line-search-like half-step backtracking when the residual norm grows between trial states.

### 10.8 Hardening law actually used by the source

The hardening update in `various_functions.f90::harden` is more specific than the generic "extended Voce" form described in the markdown notes.

First the accumulated slip increment is:

$$
\Delta \Gamma
=
\sum_\alpha |\dot\gamma^\alpha| \Delta t
$$

and the accumulated slip scalar updates as:

$$
\Gamma_{n+1} = \Gamma_n + \Delta \Gamma
$$

For each slip mode $m$, the code defines:

$$
f_m = \left|\frac{\theta_{0,m}}{\tau_{1,m}}\right|
$$

$$
e_0 = e^{-f_m \Gamma_n}, \qquad e_1 = e^{-f_m \Gamma_{n+1}}
$$

and then the mode-level hardening factor:

$$
V_m
=
\Delta \Gamma\, \theta_{1,m}
+ \tau_{1,m}(e_0 - e_1)
+ \theta_{1,m}\left(\Gamma_n e_0 - \Gamma_{n+1} e_1\right)
$$

For each slip system $\alpha$ in that mode, the interaction rate is:

$$
\Delta \tau_\alpha
=
\sum_\beta h_{\alpha\beta} |\dot\gamma^\beta|
$$

and the forward/backward slip resistances are updated as:

$$
\tau_{\alpha,n+1}^{+}
=
\tau_{\alpha,n}^{+}
+ \Delta \tau_\alpha \frac{V_m \Delta t}{\Delta \Gamma}
$$

$$
\tau_{\alpha,n+1}^{-}
=
\tau_{\alpha,n}^{-}
+ \Delta \tau_\alpha \frac{V_m \Delta t}{\Delta \Gamma}
$$

This is the exact discrete law used by the source. It is not identical to the simpler scalar law often written in the documentation as

$$
g(\Gamma)=\tau_0 + \tau_1\left(1-e^{-\theta_0 \Gamma/\tau_1}\right)+\theta_1 \Gamma
$$

although that scalar expression is still useful as intuition.

### 10.9 When hardening is applied

There are two hardening-related places in the solver:

1. inside `solve_res`, the outer fixed-point loop updates `trialtau` during the local constitutive iteration when `iuphard == 1`;
2. after an accepted macrostep, `update_hardening` calls `harden` again and commits:

   - `crss <- tautmp`
   - `trialtau <- tautmp`
   - `gacumgr <- gamacum`

Only solid voxels participate.

### 10.10 Macro tangent update

After each accepted macrostep, `calc_c0` rebuilds the macro tangent.

For each solid voxel it computes:

$$
\frac{\partial \dot\varepsilon}{\partial \sigma}
=
\frac{\partial \dot\varepsilon^e}{\partial \sigma} + \frac{\partial \dot\varepsilon^p}{\partial \sigma}
$$

and inverts it to obtain a consistent stress tangent:

$$
\frac{\partial \sigma}{\partial \varepsilon}
=
\left(
\frac{\partial \dot\varepsilon^e}{\partial \sigma} + \frac{\partial \dot\varepsilon^p}{\partial \sigma}
\right)^{-1}
$$

It also constructs a reference-timestep-scaled variant:

$$
\left(
\frac{\Delta t}{\Delta t_{\text{ref}}}\frac{\partial \dot\varepsilon^e}{\partial \sigma}
+ \frac{\partial \dot\varepsilon^p}{\partial \sigma}
\right)^{-1}
$$

to build a PK1 tangent:

$$
\left(\frac{\partial P}{\partial F}\right)_{ijkl}
=
\sum_{m,n}
\left(\frac{\partial \sigma}{\partial \varepsilon}\right)^{\text{ref}}_{imkn}
F^{-1}_{jm} F^{-1}_{ln} J
$$

The final macro objects are:

- `micro%de_ds`: average compliance used in mixed BC enforcement;
- `micro%c066`, `micro%c0`: macro stiffness used by the spectral solve;
- `micro%s0`: inverse macro stiffness.

Only solid voxels contribute.

## 11. Current Gas Model

### 11.1 What the current source actually implements

The current source implements exactly two physical branches:

1. solid voxels: full crystal plasticity with elasticity, slip kinetics, hardening, `Fp`, `Fe`, stiffness updates, and objective updates;
2. gas voxels: a special linear operator controlled by one scalar `props%complgas`.

There is no distinct `igas == 2` branch in the present Fortran source.

### 11.2 Gas compliance operator

`various_functions.f90::calc_gas_compl` builds a local 6x6 operator for every gas voxel:

$$
A_{\text{gas}}
=
I + \frac{C^{\text{mod}}}{\Delta t}\,\texttt{complgas}
- \mathcal{A}_{\text{modAL}}\, C^{\text{mod}} G_0^{r,\text{mod}}
$$

and stores:

$$
C_{\text{gas}} = A_{\text{gas}}^{-1}
$$

as `voxel%cgas`.

Solid voxels get `cgas = 0`.

### 11.3 Gas Stress Update in Solve Res

For a gas voxel, the code forms:

$$
r_{\text{gas}}
=
\lambda + C^{\text{mod}}\left(D + \frac{\sigma_t}{\Delta t}\,\texttt{complgas}\right)
- \mathcal{A}_{\text{modAL}}\, C^{\text{mod}} G_0^{r,\text{mod}} \lambda
$$

and updates stress by:

$$
\sigma = C_{\text{gas}} r_{\text{gas}}
$$

Then it enforces:

- `edotp = 0`
- `gamdot = 0`

for gas voxels.

So the gas phase is not an elastic crystal phase, not a plastic phase, and not a separate damper model in the current source. It is a linear operator with memory through the previous-step stress `sgt`.

### 11.4 Why gas material files do not matter today

The current Python layer writes gas files for `igas == 2`:

- `cuel_gas.sx`
- `cupl2_gas.sx`

but current Fortran never reads them because any nonzero `igas` is parsed as logical gas and skipped in `load_input`.

Therefore:

- `GasElastic` and `GasPlastic` in `utils/config.py` are presently documentation/plumbing for a not-currently-active model;
- `complgas` is the only gas constitutive parameter that currently affects the Fortran gas branch.

### 11.5 Historical and Planned igas 2 Behavior

Several docs and Python comments still describe a historically present or planned `igas == 2` "damper" behavior:

- isotropic elastic constants are read;
- slip kinetics provide dissipation;
- hardening evolution is skipped;
- `Fp`, `Fe`, and some bookkeeping updates are skipped.

This behavior is not implemented in the current Fortran solver, but it remains relevant as a reimplementation target.

To restore a true `igas == 2` branch, the minimum architectural change would be:

1. change `phase_type%igas` from logical back to an integer or enum-like mode;
2. keep the original integer value in `load_input`;
3. branch separately in:
   - input reading
   - `solve_res`
   - `update_plastic_defgrad`
   - `update_elstic_defgrad`
   - `update_el_stiff`
   - `update_tensor_config`
   - any place where `voxel%gas` currently skips constitutive or kinematic updates

As the code stands today, none of those distinctions exist.

## 12. Accepted-Step Finite-Strain Updates

Once a macrostep is accepted, the solver commits the converged increment.

### 12.1 Immediate state accumulation

In `LS-EVPFFT.f90` the solver first updates:

$$
\varepsilon^p_t \leftarrow \varepsilon^p_t + \dot\varepsilon^p \Delta t
$$

$$
\sigma_t \leftarrow \sigma
$$

$$
\nabla u^{\text{sym}} \leftarrow \nabla u^{\text{sym}} + \operatorname{sym}(L)\Delta t
$$

$$
\bar\varepsilon_{\text{vm}}^p \leftarrow \bar\varepsilon_{\text{vm}}^p + \varepsilon_{\text{vm}}(\dot\varepsilon^p)\Delta t
$$

and tracks the maximum solid von Mises strain-rate `micro%edotvmmx`.

### 12.2 Plastic deformation gradient

`update_plastic_defgrad` first computes the plastic velocity gradient in the intermediate configuration:

$$
L^p = \sum_\alpha \dot\gamma^\alpha \, b^\alpha_{\text{sample}} \otimes n^\alpha_{\text{sample}}
$$

where the sample-basis directions are obtained from `ag`.

Then it updates:

$$
F^p_{n+1} = \exp(L^p \Delta t)\, F^p_n
$$

Gas voxels are skipped.

### 12.3 Total deformation gradient

`update_defgrad` applies the finite-strain update:

$$
F_{n+1} = \exp(L \Delta t) F_n
$$

using `matrix_exp_adaptive`.

That adaptive exponential chooses:

- a cubic Taylor expansion when the matrix 1-norm is at most `1e-1`;
- a 13th-order Padé approximation with scaling and squaring otherwise.

Then it computes:

- `detF = det(F)`
- `F^{-1}` by LU inversion
- `defgradinc = exp(L dt)`

and updates the current-volume quantities:

$$
w_i^c = \frac{w_i^0 J_i}{\sum_j w_j^0 J_j}
$$

as well as:

- `defgradavg`
- `defgradinvavgc`
- `defgradinvavgcs`
- `defgradinvavgcg`
- `wphcsol`
- `wphcgas`

### 12.4 Elastic deformation gradient

`update_elstic_defgrad` updates:

$$
F^e = F (F^p F_{\text{ini}})^{-1}
$$

and the elastic increment:

$$
F^{e,\text{inc}} = F^e_{n+1} (F^e_n)^{-1}
$$

Gas voxels are skipped.

### 12.5 Current elastic stiffness

`update_el_stiff` push-forwards the crystal stiffness tensor with:

$$
A = F^e R
$$

and computes:

$$
c_{ijkl}
=
\frac{1}{\det F^e}
\sum_{I,J,K,L}
A_{iI} A_{jJ} C_{IJKL} A^T_{Kk} A^T_{Ll}
$$

This is converted to `cg66`, and its inverse gives `sg66`.

Again, gas voxels are skipped.

### 12.6 Objective rotation and stress push-forward

`update_tensor_config` performs the objective update after the configuration change.

For solids it:

1. computes the polar decomposition of `defgradinc`:

   $$
   F^{\text{inc}} = R^{\text{inc}} U^{\text{inc}}
   $$

2. rotates the stored symmetric tensors:

   $$
   \varepsilon^p \leftarrow R^{\text{inc}} \varepsilon^p (R^{\text{inc}})^T
   $$

   $$
   \nabla u^{\text{sym}} \leftarrow R^{\text{inc}} \nabla u^{\text{sym}} (R^{\text{inc}})^T
   $$

3. pushes the previous-step stress with the elastic increment:

   $$
   \sigma_t \leftarrow \frac{F^{e,\text{inc}} \sigma_t (F^{e,\text{inc}})^T}{\det F^{e,\text{inc}}}
   $$

For all voxels it then transforms the velocity gradient as:

$$
L \leftarrow L (F^{\text{inc}})^{-1}
$$

and applies an average correction so the macro average remains consistent after the objective update.

## 13. Timestep Control, Substepping, and Current Iteration Safeguards

The current code contains two active step-size controls.

### 13.1 Strain-controlled timestep rescaling

When `ictrl == 0`, the process input `tdot` is interpreted as a target equivalent strain increment rather than a physical time increment.

During the equilibrium loop the solver computes:

$$
\Delta \varepsilon_{\text{vm}} = \varepsilon_{\text{vm}}(\operatorname{sym}\bar L)
$$

and resets:

$$
\Delta t = \frac{\Delta \varepsilon_{\text{vm,target}}}{\Delta \varepsilon_{\text{vm}}}
$$

storing the updated value in both `tdot` and `tdotref`.

### 13.2 Kinematic substepping and rollback

When `ictrl == -1`, the solver can cut the step and roll back to the macrostep snapshot.

The routine `kinematic_substep` estimates the required number of substeps from the worst of:

- volumetric increment limit `volinc_lim = 5e-2`
- rotation angle limit `rotang_lim = 2.0`
- stretch increment limit `stretchinc_lim = 5e-2`

These are evaluated from:

$$
F^{\text{inc}} = \exp(L \Delta t)
$$

with:

- `volinc = |log(det(Finc))|`
- `rotang = acos((trace(Rinc)-1)/2)`
- `stretchinc = ||log Uinc||_F`

There is also a hard floor `detF_floor = 1e-3` in `LS-EVPFFT.f90`.

### 13.3 Step rejection path

If the equilibrium loop fails to converge or if the projected kinematic increments are too large, the code:

- computes a reduced `tdot`;
- increases the number of remaining macrosteps accordingly;
- resets `iter = 0` and `mix = 1` when needed;
- restores `velgrad`, `velgradold`, `sg`, `edotp`, `gamdot`, and `trialtau` from the snapshot via `rollback`.

## 14. Diagnostics and Output Files

### 14.1 Results CSV

`IO_functions.f90::write_results` writes one line per accepted macrostep with columns:

- run identifiers: `iter`, `macrostep`, `ipr`, `sim_time`, `wall_time`, `tdot`, `tdot_ref`
- averaged total strain tensor `eav..`
- averaged solid stress tensor `sav..`
- averaged plastic strain tensor `epav..`
- averaged elastic strain tensor `eelav..`
- average velocity gradient `vgrad..`
- average plastic strain rate `edotp..`
- scalars `evm`, `evmp`, `dvm`, `dvmp`, `svm`

Notable implementation detail:

- `sav..` is written from `micro%sgsolavg`, i.e. the normalized solid-only average stress after the step is accepted.

### 14.2 Diagnostics CSV

`types.f90::write_diagnostics` writes:

- `vm_max_inc`
- `tdotref`
- `erre`, `errs`, `errsbc`
- `s33_target`, `s33`, `s33_sol`
- min-`detF` statistics and relative changes
- max volumetric/stretch/rotation increment statistics and relative changes
- `dL_step_*`
- `dvelgradavg_mag`
- `cond_c066mod_max_sol`
- `cond_c066mod_max_gas`
- `cond_cgas_max`
- voxel indices of the worst solid/gas minima

### 14.3 What the diagnostics mean mathematically

`types.f90::update_diagnostics` recomputes per-step quantities from the converged field:

1. minimum Jacobian

   $$
   \min_i J_i
   $$

2. volumetric increment

   $$
   \max_i |\log \det(F_i^{\text{inc}})|
   $$

3. rotation increment

   $$
   \max_i \arccos\left(\frac{\operatorname{tr}(R_i^{\text{inc}})-1}{2}\right)
   $$

4. stretch increment

   $$
   \max_i \|\log U_i^{\text{inc}}\|_F
   $$

5. per-step velocity-gradient change relative to the macrostep snapshot

   $$
   \frac{\|L_i - L_i^{\text{snapshot}}\|_F}
   {\max(\|L_i^{\text{snapshot}}\|_F, L_{\text{floor}})}
   $$

6. 1-norm condition-number estimates for the worst `c066mod` and `cgas` matrices using LAPACK `dgetrf`/`dgecon`.

### 14.4 Python-side loaders

`utils/data_utils.py` defines:

- `SimResults.load(...)`
- `Diagnostics.load(...)`

The loaders:

- require fixed CSV column sets;
- cast integer and floating columns explicitly;
- optionally sweep VTK files to compute roughness and slip metrics;
- truncate the time history at the first nonfinite solver result.

## 15. Calibration, Loss Functions, and Postprocessing

### 15.1 What the docs recommend conceptually

The markdown notes on loss normalization and calibration argue for:

- deterministic response scaling rather than ad hoc weights;
- explicit early-transient weighting;
- optional slope or shape penalties;
- multiobjective combinations over stress and strain;
- logging slip activity and accumulated slip to determine whether early mismatch is actually plasticity-controlled.

Those recommendations are conceptually consistent with the constitutive model and remain good design guidance, but they are not the exact behavior of the current calibration code.

### 15.2 What the current calibration code actually computes

The current implementation is in `utils/calibrate.py` and `utils/calc_utils.py`.

The three basic aligned metrics are:

1. MAE

   $$
   \text{MAE}
   =
   \frac{1}{N}
   \sum_i
   \frac{|y_i^{\text{exp}} - y_{n(i)}^{\text{sim}}|}{\max(y^{\text{exp}})}
   $$

2. MSE

   $$
   \text{MSE}
   =
   \frac{1}{N}
   \sum_i
   \left(
   \frac{y_i^{\text{exp}} - y_{n(i)}^{\text{sim}}}{\max(y^{\text{exp}})}
   \right)^2
   $$

3. relative error

   $$
   \text{MRE}
   =
   \frac{1}{N}
   \sum_i
   \frac{|y_i^{\text{exp}} - y_{n(i)}^{\text{sim}}|}
   {\max(|y_i^{\text{exp}}|, |y_{n(i)}^{\text{sim}}|, \varepsilon)}
   $$

where `n(i)` is the nearest simulation sample to experimental abscissa `x_i`.

There is:

- no early weighting in the current code;
- no Huber penalty;
- no explicit slope/shape term;
- no deterministic scale parameter separate from `max(exp_y)`.

### 15.3 Creep objective

For creep calibration:

- the strain series is either total strain `eav33` or plastic strain `epav33`;
- the comparison is aligned on experimental time points;
- the stress objective compares simulated `sav33` against nominal applied load;
- the total objective is stress, strain, or their sum depending on `cfg.target`.

### 15.4 SRJ objective

For SRJ calibration:

- the code supports alignment on strain or time via `cfg.srj_align`;
- strain history is always compared as a function of time;
- the stress objective is intended to compare stress against either strain or time depending on alignment.

### 15.5 Combined loss

At the run level, `calibrate.py` combines the base constitutive/data misfit with an iteration penalty:

$$
\text{loss} = \text{base} + w_{\text{iters}} \,\overline{N}_{\text{iters}}
$$

and in multiobjective mode:

$$
\text{loss}_{\text{joint}}
=
w_s \,\text{loss}_{\text{srj}} + w_c \,\text{loss}_{\text{creep}}
$$

### 15.6 Postprocessing

`utils/postprocess.py` is intentionally thin:

1. load `SimResults`
2. load `Diagnostics`
3. call `plotting/mechanical.py`
4. call `plotting/diagnostics.py`

The plot layer can use VTK-derived fields to generate:

- surface-height measures
- roughness metrics `Sa`, `Sz`
- maximum and mean slip on free-surface-adjacent faces

when the VTK files exist and `skip_vtk` is false.

## 16. Material Parameters: Meaning, Priors, and Current Defaults

### 16.1 Parameter meanings

The solid plasticity parameter set used across the docs and source is:

- `tau0xf`: initial forward CRSS
- `tau0xb`: initial backward CRSS
- `tau1x`: saturation-like hardening scale
- `thet0`: early hardening-rate parameter
- `thet1`: late hardening-rate parameter
- `nrsx`: stress exponent
- `gamd0x`: reference slip rate
- `hselfx`: self-hardening coefficient
- `hlatex`: latent-hardening coefficient

The common rate-sensitivity notation mapping is:

$$
m = \frac{1}{n}
$$

with `nrsx = n`.

### 16.2 Documentation priors versus active code defaults

`docs/material_parameters.md` proposes materially different "v2" priors than the active imported defaults in `utils/config.py`.

The current `config.py` defines `SolidElastic` and `SolidPlastic` twice. Because Python keeps the later definitions, the second block is the active one imported by the rest of the repo.

That means the effective current defaults are approximately:

- `tau0xf = 200`
- `rb = 0.9`
- `tau0xb = 180`
- `tau1x = 10`
- `thet0 = 1000`
- `thet1 = 100`
- `nrsx = 100`
- `gamd0x = 1e-8`

By contrast, `docs/material_parameters.md` recommends much more aggressive solid priors for AM 316L, roughly:

- `tau0xf = 185 MPa`
- `tau1x = 200 MPa`
- `thet0 = 5000 MPa`
- `thet1 = 300 MPa`
- `nrsx = 20`
- `gamd0x = 1e-4 1/s`

So the documentation and the active Python defaults do not currently describe the same baseline material behavior.

### 16.3 Gas parameter defaults

`utils/config.py` still defines:

- `GasElastic(young, nu)`
- `GasPlastic(tau0xf, tau0xb, tau1x, thet0, thet1, nrsx, gamd0x, ...)`
- `COMPLGAS_BOUNDS`

but only `complgas` currently affects the Fortran gas branch.

The comment in `config.py` is accurate for the current codebase:

- gas plastic parameters are "written/read for `igas = 2` but never used" in the current active solver.

## 17. The 0D Bulk Surrogate Note

`docs/0d_bulk_surrogate_model.md` describes a conceptual material-point surrogate for uniaxial bulk response. It is not currently implemented as an executable module under `utils/`.

The surrogate note is still useful because it clarifies the intended constitutive interpretation of the LS-EVPFFT parameters:

1. Taylor-factor closure

   $$
   \tau \approx \frac{\sigma}{M}
   $$

2. effective single-channel slip law

   $$
   \dot\gamma
   =
   \dot\gamma_0
   \left(\frac{|\tau|}{g}\right)^n
   \operatorname{sign}(\tau)
   $$

3. accumulated slip

   $$
   \dot\Gamma = |\dot\gamma|
   $$

4. scalar hardening

   $$
   g(\Gamma)
   =
   \tau_0
   +
   \tau_1\left(1-e^{-\theta_0 \Gamma/\tau_1}\right)
   +
   \theta_1 \Gamma
   $$

5. uniaxial plastic strain estimate

   $$
   \dot\varepsilon_p \approx \frac{\dot\Gamma}{M}
   $$

This surrogate is best understood as a conceptual reduction of the same slip-law ingredients, not as a replacement for the full-field FFT solver.

## 18. Code Map by Algorithmic Responsibility

For practical navigation, the most important code-to-algorithm mapping is:

- main driver and macrostep control:
  `lsevpfft/src/LS-EVPFFT.f90`

- type definitions, diagnostics accumulation, diagnostics CSV:
  `lsevpfft/src/types.f90`

- input parsing, crystal file reading, elastic file reading, results CSV, VTK writing:
  `lsevpfft/src/IO_functions.f90`

- Schmid updates, Green operator, mapped stiffness, gas operator, FFT correction, constitutive solve, hardening, BC correction, finite-strain updates, macro tangent:
  `lsevpfft/src/various_functions.f90`

- FFTW plans and tensor FFT wrappers:
  `lsevpfft/src/fourier_functions.f90`

- basis transforms, tensor algebra, LU solves, polar decomposition, LAPACK helpers:
  `lsevpfft/src/tensor_functions.f90`

- Python run orchestration:
  `utils/run.py`

- Python `fft.in` dispatch:
  `utils/write_fft.py`

- Python config/defaults/bounds:
  `utils/config.py`

- calibration driver:
  `utils/calibrate.py`

- metrics and mechanical helper math:
  `utils/calc_utils.py`

- typed CSV/VTK loaders:
  `utils/data_utils.py`

- plotting/postprocessing:
  `utils/postprocess.py`, `utils/plotting/`

## 19. Baseline Summary

The current `mimosa` baseline can be summarized as follows.

At present, the repository is a large-strain FFT crystal plasticity solver plus a Python orchestration stack. The solid constitutive model is fully active and implemented in detail: slip-system geometry is read from `cupl2.sx`, Schmid tensors are updated in the current elastic configuration, slip rates follow a power law with forward/backward CRSS selection, elastic and plastic strain-rate contributions are combined in a local Newton solve, and hardening follows the exact discrete Voce-like update in `various_functions.f90::harden`.

The spectral solve is equally concrete: the code builds a Green operator from the current macro tangent, forms a PK1-like stress fluctuation, convolves it spectrally, maps the correction back to the current configuration, adaptively mixes the resulting velocity-gradient field, and then enforces mixed stress/strain-rate boundary conditions through a 6x6 masked solve.

The current gas model is much narrower than several docs and Python entry points imply. In the active Fortran source, any nonzero `igas` is just logical gas, and gas voxels are updated by a linear compliance operator governed by `complgas`; no distinct `igas == 2` damper branch exists anymore. Historical `igas == 2` discussion remains relevant because the Python layer and several docs still describe that design, and reimplementing it would require restoring an integer-valued gas mode throughout the Fortran solver.

The Python calibration and postprocessing stack is useful but simpler than some of the docs describe. The present calibration code uses aligned MAE/MSE/MRE objectives normalized by the maximum experimental response, not the richer weighted and shape-aware losses described in the notes, and the current `config.py` defaults do not match the newer material-prior recommendations in `docs/material_parameters.md`.
