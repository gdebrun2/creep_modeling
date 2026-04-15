program fft
!$  use omp_lib
    use types
    use IO_functions
    use fourier_functions
    use various_functions
    use tensor_functions

    implicit none

    type(micro_type) :: micro ! microstructure state (voxel fields, averages, operators)
    type(props_type) :: props ! simulation inputs (grid, phases, loading processes)
    type(diag_type) :: diag ! macrostep diagnostics (one record per step)

    double precision :: mix ! relaxation/mixing parameter for velgrad update
    double precision :: velgradavg_sym(3, 3) ! symmetric part of macroscopic avg velgrad
    double precision :: evmavg ! von Mises strain increment used for timestep control
    double precision :: edotvm ! per-voxel von Mises strain-rate (scratch)
    double precision :: velgrad_sym(3, 3) ! symmetric part of voxel velgrad (scratch)

    double precision :: aux33(3, 3) ! scratch 3x3 tensor
    double precision :: aux6(6) ! scratch 6-vector (Voigt/basis)
    double precision :: sgavg(3, 3) ! macroscopic average Cauchy stress accumulator

    integer :: start ! wall-clock start time (seconds)
    integer :: i, j, k ! generic loop indices
    integer :: nthreads = 0 ! requested number of threads (0 => default/max)
    integer :: nsteps_left ! steps remaining in current process
    integer :: nsteps_total_new ! replanned total number of steps in current process
    integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)
    integer :: is, nsyst_loc ! slip-system loop indices/counts for accepted-step averages

    double precision, parameter :: mix_min = 1.0d-8 ! lower bound on mixing/relaxation
    double precision, parameter :: dL_target = 5.0d-2 ! target rel. velgrad change for mix ~ 1
    double precision, parameter :: detF_floor = 1.0d-3 ! timestep-cut threshold for detF_min
    double precision :: Lscale ! scaling for relative dL_iter diagnostic
    double precision :: dL_iter_sol, dL_iter_gas, dL_iter ! per-iteration velgrad relative change diagnostics
    double precision :: sum_rel2_sol, sum_w_sol ! weighted RMS accumulators (solid)
    double precision :: sum_rel2_gas, sum_w_gas ! weighted RMS accumulators (gas)
    double precision :: dL, rel, wgtc ! per-voxel scratch for dL_iter

    logical :: skip_next = .false. ! command-line parsing: previous arg consumed next token
    character(len=256) :: arg ! command-line argument token buffer
    character(len=256) :: fft_input_file = 'fft.in' ! main input file path (overridable via CLI)
    character(len=256) :: output_dir = "evpfft_outputs/"
    character(len=256) :: fftw_save_path = ""
    integer :: nconverged = 0 ! count number of steps that have converged to try and ramp tdot up
    logical :: did_converge ! increment converged counter for tdot ramp up
    logical :: modal = .true. ! toggle on or off the modified AL
    integer :: tmp_int ! dummy
    logical :: kinematic_controlled ! to turn off velgrad mixing for srj case

    double precision :: sim_time = 0.0
    double precision :: wall_time = 0.0
    double precision :: tdot_used
    double precision :: tdot
    double precision :: tdot_plan
    double precision :: proc_time_target
    double precision :: proc_time_done
    double precision :: remaining_time

    integer :: pct = 0
    integer :: total_macrosteps = 0
    integer :: nsub
    ! double precision :: dt_sub
    logical :: do_mix = .true.
    ! logical :: overstep = .false.
    double precision :: max_volinc = 0.0
    double precision :: max_rotang = 0.0
    double precision :: max_stretchinc = 0.0
    double precision :: min_detF = huge(1.0d0)
    logical :: accept_step = .false.
    double precision :: gammaavg = 0.0d0
    double precision :: gavg = 0.0d0
    double precision :: geffavg = 0.0d0
    double precision :: sum_gamma = 0.0d0
    double precision :: sum_g = 0.0d0
    double precision :: sum_gdot = 0.0d0
    double precision :: sum_ggdot = 0.0d0
    double precision :: gsys = 0.0d0
    double precision :: nsolid = 0.0d0
    double precision :: nsys_solid = 0.0d0

    write (*, *)

    ! Read optional command line arguments
    do i = 1, command_argument_count()
        if (skip_next) then
            skip_next = .false.
            cycle
        end if
        call get_command_argument(i, arg)

        select case (arg)
        case ('--in')
            skip_next = .true.
            call get_command_argument(i + 1, arg)
            fft_input_file = trim(arg)
        case ('--out')
            skip_next = .true.
            call get_command_argument(i + 1, arg)
            output_dir = trim(arg)

        case ('--nthreads')
            skip_next = .true.
            call get_command_argument(i + 1, arg)
            read (arg, *) nthreads

        case ('--fftw_save')
            skip_next = .true.
            call get_command_argument(i + 1, arg)
            fftw_save_path = trim(arg)

        case ('--modal')
            skip_next = .true.
            call get_command_argument(i + 1, arg)
            read (arg, *) tmp_int
            modal = (tmp_int /= 0)



        case ('--mix')
            skip_next = .true.
            call get_command_argument(i + 1, arg)
            read (arg, *) tmp_int
            do_mix = (tmp_int /= 0)

        case ('-h', '--help')
            write (*, '(A)') 'Optional inputs:'
            write (*, '(A)') '--in FILE | name of input file (fft.in by default)'
            write (*, '(A)') '--out DIR | output directory (evpfft_outputs by default)'
            write (*, '(A)') '--nthreads n | number of threads to run with (# physical cpu cores by default)'
            write (*, '(A)') '--save FILE | path to FFTW wisdom file (default: <outdir>/fftw_wisdom.save)'
            write (*, "(A)") '--modal 0-1 | toggle modified AL (default 1)'
            stop

        case default
            write (*, '(3A)') 'Unknown command-line option: "', trim(arg), '" use -h or --help for argument list'
            stop

        end select
    end do

    call system('mkdir -p "'//trim(output_dir)//'"')
    call system('find "'//trim(output_dir)//'" -maxdepth 1 -type f -name "*.vtk" -delete')
    call system('find "'//trim(output_dir)//'" -maxdepth 1 -type f -name "*.out" -delete')

    call init_fftw_wisdom_path(fftw_save_path)
    call load_input(micro, props, nthreads, fft_input_file)
    call initial_guess(micro, props)
    call write_vtk_cell(micro, props, 0, output_dir)

    start = time()
    props%modAL = modal

    micro%c066 = 0.0
    do ip3 = 1, props%npts3
    do ip2 = 1, props%npts2
    do ip1 = 1, props%npts1
        if (.not. micro%voxel(ip1, ip2, ip3)%gas) micro%c066 = micro%c066 + &
                                                               micro%voxel(ip1, ip2, ip3)%cg66 * micro%voxel(ip1, ip2, ip3)%wgtc
    end do
    end do
    end do

    micro%c066 = micro%c066 / micro%wphcsol
    micro%c066 = micro%c066 * props%process(1)%tdot
    micro%s066 = micro%c066

    call lu_inverse(micro%s066, 6)
    call chg_basis(aux6, aux33, micro%c066, micro%c0, 3, 6)
    call chg_basis(aux6, aux33, micro%s066, micro%s0, 3, 6)
    micro%de_ds = micro%s0

    ! initialize results and diagnostics
    call diag%write_diagnostics(output_dir, 0, 0, 0, 0.0d0, 0.0d0, 0.0d0, init=.true.)
    call write_results(output_dir, micro, 0, 0, 0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                       0.0d0, 0.0d0, 0.0d0, init=.true.)

    do ipr = 1, props%nproc
        total_macrosteps = total_macrosteps + props%process(ipr)%nsteps
    end do

    ! Loop over processes
    imacroloop = 1
    ipr = 0
    write (*, *)
    do while (ipr < props%nproc)
        ipr = ipr + 1
        nconverged = 0
        kinematic_controlled = (props%process(ipr)%iscau(3) == 0)
        tdot_plan = props%process(ipr)%tdot
        proc_time_target = dble(props%process(ipr)%nsteps) * tdot_plan
        proc_time_done = 0.0d0

        ! Loop over number of steps
        imicro = 0
        do while (imicro < props%process(ipr)%nsteps)
            imicro = imicro + 1
            if (props%process(ipr)%ictrl == -1) then
                remaining_time = max(0.0d0, proc_time_target - proc_time_done)
                if (remaining_time > 0.0d0) then
                    props%process(ipr)%tdot = min(tdot_plan, remaining_time)
                else
                    props%process(ipr)%tdot = tdot_plan
                end if
            end if
            call take_snapshot(micro, props, ipr)

            did_converge = .false.
            accept_step = .false.

            if ((imicro == 1 .and. ipr == 1) .or. props%process(ipr)%iupdate == 1) call update_schmid(micro, props)
            call form_G(micro, props)
            call calc_c066mod_Goperr066mod(micro, props)

            iterNRavg_accum = 0.0
            iter = 0

            do while (.not. (did_converge .or. accept_step))

                if (iter == 0) then

                    call calc_gas_compl(micro, props)
                    micro%dvelgradavg = 0.0
                    micro%edotvmmx = 0
                    mix = 1.0
                    erre = 2.0 * props%process(ipr)%error
                    errs = 2.0 * props%process(ipr)%error
                    errsbc = 2.0 * props%process(ipr)%error
                end if

                iter = iter + 1

                pct = (100 * imacroloop) / max(1, total_macrosteps)
                call write_status(imacroloop, sim_time, wall_time, ipr, imicro, props%process(ipr)%tdot, iter, pct)

                ! PK1
                sgavg = 0.0
                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP PRIVATE(i, j, k) &
                !$OMP SHARED(props, micro) &
                !$OMP REDUCTION(+:sgavg)
                !$OMP DO COLLAPSE(3)
                do k = 1, props%npts3
                do j = 1, props%npts2
                do i = 1, props%npts1
                    sgavg = sgavg + micro%voxel(i, j, k)%sg * micro%voxel(i, j, k)%wgtc
                end do
                end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL
                micro%sgavg = sgavg

                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP PRIVATE(i, j, k) &
                !$OMP SHARED(props, micro)
                !$OMP DO COLLAPSE(3)
                do k = 1, props%npts3
                do j = 1, props%npts2
                do i = 1, props%npts1
                    micro%voxel(i, j, k)%sgPK1 = CauchyToPK1(micro%voxel(i, j, k)%sg - micro%sgavg, &
                                                             micro%voxel(i, j, k)%defgradinv, micro%voxel(i, j, k)%detF)
                end do
                end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL

                ! convolution
                call convolution_with_Goper(micro, props)

                ! to current configuration
                call correction_to_current(micro, props)

                ! Compute dL_iter for mixing (from raw velgrad update) and apply mixing.
                sum_rel2_sol = 0.0d0; sum_w_sol = 0.0d0
                sum_rel2_gas = 0.0d0; sum_w_gas = 0.0d0
                Lscale = max(micro%lvmavg, tiny(1.0d0))
                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP PRIVATE(i, j, k, dL, rel, wgtc) &
                !$OMP SHARED(props, micro, mix, Lscale, kinematic_controlled, do_mix) &
                !$OMP REDUCTION(+:sum_rel2_sol, sum_w_sol, sum_rel2_gas, sum_w_gas)
                !$OMP DO COLLAPSE(3)
                do k = 1, props%npts3
                do j = 1, props%npts2
                do i = 1, props%npts1
                    if (.not. kinematic_controlled) then
                        dL = sqrt(sum((micro%voxel(i, j, k)%velgrad - micro%voxel(i, j, k)%velgradold)**2))
                        rel = dL / Lscale
                        wgtc = micro%voxel(i, j, k)%wgtc
                        if (micro%voxel(i, j, k)%gas) then
                            sum_rel2_gas = sum_rel2_gas + (rel * rel) * wgtc
                            sum_w_gas = sum_w_gas + wgtc
                        else
                            sum_rel2_sol = sum_rel2_sol + (rel * rel) * wgtc
                            sum_w_sol = sum_w_sol + wgtc
                        end if
                    end if
                    if (do_mix) then
                        micro%voxel(i, j, k)%velgrad = mix * micro%voxel(i, j, k)%velgrad + (1.0 - mix) * micro%voxel(i, j, k)%velgradold
                    else
                        micro%voxel(i, j, k)%velgrad = mix * micro%voxel(i, j, k)%velgrad
                    end if

                end do
                end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL

                if (.not. kinematic_controlled) then
                    if (sum_w_sol > 0.0d0) then
                        dL_iter_sol = sqrt(sum_rel2_sol / sum_w_sol)
                    else
                        dL_iter_sol = 0.0d0
                    end if
                    if (sum_w_gas > 0.0d0) then
                        dL_iter_gas = sqrt(sum_rel2_gas / sum_w_gas)
                    else
                        dL_iter_gas = 0.0d0
                    end if

                    dL_iter = max(dL_iter_sol, dL_iter_gas)
                    mix = max(mix_min, min(1.0d0, dL_target / max(dL_iter, dL_target)))
                end if

                ! update stress field
                call solve_res(micro, props)
                ! modifies:
                ! micro%voxel(ip1, ip2, ip3)%edotp
                ! micro%voxel(ip1, ip2, ip3)%sg
                ! micro%voxel(ip1, ip2, ip3)%gamdot(:)
                ! micro%voxel(ip1, ip2, ip3)%trialtau(:, :)

                ! update average
                call calc_velgradavg_corr(micro, props)

                ! normalize errors
                erre = erre / micro%lvmavg
                errs = errs / micro%svmavg
                errsbc = errsbc / micro%svmavg

                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP PRIVATE(i, j, k) &
                !$OMP SHARED(props, micro)
                !$OMP DO COLLAPSE(3)
                do k = 1, props%npts3
                do j = 1, props%npts2
                do i = 1, props%npts1
                    ! modifies velgradold
                    micro%voxel(i, j, k)%velgradold = micro%voxel(i, j, k)%velgrad
                end do
                end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL

                ! readjust the time increment
                if (props%process(ipr)%ictrl == 0) then

                    velgradavg_sym = sym33(micro%velgradavg)
                    evmavg = vm_strain(velgradavg_sym) ! increment in von Mises strain
                    props%process(ipr)%tdot = props%process(ipr)%devmapp / evmavg
                    props%process(ipr)%tdotref = props%process(ipr)%tdot

                end if

                if ((errs <= props%process(ipr)%error .and. &
                     erre <= props%process(ipr)%error .and. &
                     errsbc <= props%process(ipr)%error)) then

                    did_converge = .true.

                end if

                if (props%process(ipr)%ictrl == -1) then
                    if (props%process(ipr)%tdot > props%process(ipr)%tdotmin) then
                        if (.not. did_converge .and. iter == props%process(ipr)%itmax) then
                            tdot = props%process(ipr)%tdot
                            nsub = kinematic_substep(micro, props, tdot, max_volinc, max_rotang, max_stretchinc, min_detF)
                            tdot = tdot / nsub
                            tdot = max(props%process(ipr)%tdotmin, min(tdot, 0.5d0 * props%process(ipr)%tdot))
                            nsub = ceiling(props%process(ipr)%tdot / tdot)
                            if (nsub > 1) then
                                remaining_time = max(0.0d0, proc_time_target - proc_time_done)
                                tdot_plan = tdot
                                nsteps_total_new = (imicro - 1) + ceiling(remaining_time / tdot_plan)
                                total_macrosteps = total_macrosteps + nsteps_total_new - props%process(ipr)%nsteps
                                props%process(ipr)%nsteps = nsteps_total_new
                                props%process(ipr)%tdot = min(tdot_plan, remaining_time)
                                nconverged = 0
                                iter = 0
                                call rollback(micro, props, ipr)
                            end if

                        else if (did_converge) then
                            tdot = props%process(ipr)%tdot
                            nsub = kinematic_substep(micro, props, tdot, max_volinc, max_rotang, max_stretchinc, min_detF)
                            tdot = tdot / nsub
                            if (min_detF <= detF_floor) then
                                tdot = max(props%process(ipr)%tdotmin, min(tdot, 0.5d0 * props%process(ipr)%tdot))
                            else
                                tdot = max(props%process(ipr)%tdotmin, tdot)
                            end if
                            nsub = ceiling(props%process(ipr)%tdot / tdot)
                            if (nsub > 1) then
                                remaining_time = max(0.0d0, proc_time_target - proc_time_done)
                                tdot_plan = tdot
                                nsteps_total_new = (imicro - 1) + ceiling(remaining_time / tdot_plan)
                                total_macrosteps = total_macrosteps + nsteps_total_new - props%process(ipr)%nsteps
                                props%process(ipr)%nsteps = nsteps_total_new
                                props%process(ipr)%tdot = min(tdot_plan, remaining_time)
                                did_converge = .false.
                                nconverged = 0
                                iter = 0
                                mix = 1.0d0
                                call rollback(micro, props, ipr)
                            end if
                        end if

                    else if (iter == props%process(ipr)%itmax) then
                        accept_step = .true.
                    end if

                end if

            end do

            ! update stress, plastic strain, displacement gradient, plastic von misses and max strain rate
            do ip3 = 1, props%npts3
            do ip2 = 1, props%npts2
            do ip1 = 1, props%npts1
                micro%voxel(ip1, ip2, ip3)%ept = micro%voxel(ip1, ip2, ip3)%ept + micro%voxel(ip1, ip2, ip3)%edotp * props%process(ipr)%tdot
                micro%voxel(ip1, ip2, ip3)%sgt = micro%voxel(ip1, ip2, ip3)%sg
                velgrad_sym = sym33(micro%voxel(ip1, ip2, ip3)%velgrad)
                micro%voxel(ip1, ip2, ip3)%disgradsym = micro%voxel(ip1, ip2, ip3)%disgradsym + &
                                                        velgrad_sym * props%process(ipr)%tdot
                micro%voxel(ip1, ip2, ip3)%epvm = micro%voxel(ip1, ip2, ip3)%epvm + &
                                                  vm_strain(micro%voxel(ip1, ip2, ip3)%edotp) * props%process(ipr)%tdot
                edotvm = vm_strain(velgrad_sym)
                if (edotvm > micro%edotvmmx .and. .not. micro%voxel(ip1, ip2, ip3)%gas) micro%edotvmmx = edotvm ! set the worst solid vm strain rate
            end do
            end do
            end do

            ! Updating
            if (props%process(ipr)%iupdate == 1) then
                call update_plastic_defgrad(micro, props)
                call update_grid_velgrad_node(micro, props)
                call update_defgrad(micro, props)
                call update_elstic_defgrad(micro, props)
                call update_el_stiff(micro, props)
                call update_tensor_config(micro, props)
            end if

            if (props%process(ipr)%iuphard == 1) call update_hardening(micro, props)

            ! update averages
            velgradavg_sym = sym33(micro%velgradavg)
            micro%lvmavg = vm_strain(micro%velgradavg)
            micro%dvmavg = vm_strain(velgradavg_sym)
            micro%evmavg = micro%evmavg + micro%dvmavg * props%process(ipr)%tdot
            micro%epavg = 0.0
            micro%eavg = 0.0
            micro%eelavg = 0.0
            micro%edotpavg = 0.0
            micro%sgsolavg = 0.0
            sum_gamma = 0.0d0
            sum_g = 0.0d0
            sum_gdot = 0.0d0
            sum_ggdot = 0.0d0
            nsolid = 0.0d0
            nsys_solid = 0.0d0
            do ip3 = 1, props%npts3
            do ip2 = 1, props%npts2
            do ip1 = 1, props%npts1
                if (.not. micro%voxel(ip1, ip2, ip3)%gas) then
                    micro%epavg = micro%epavg + micro%voxel(ip1, ip2, ip3)%ept * micro%voxel(ip1, ip2, ip3)%wgtc
                    micro%eavg = micro%eavg + micro%voxel(ip1, ip2, ip3)%disgradsym * micro%voxel(ip1, ip2, ip3)%wgtc
                    micro%eelavg = micro%eavg + (micro%voxel(ip1, ip2, ip3)%disgradsym - micro%voxel(ip1, ip2, ip3)%ept) &
                                   * micro%voxel(ip1, ip2, ip3)%wgtc
                    micro%edotpavg = micro%edotpavg + micro%voxel(ip1, ip2, ip3)%edotp * micro%voxel(ip1, ip2, ip3)%wgtc
                    micro%sgsolavg = micro%sgsolavg + micro%voxel(ip1, ip2, ip3)%sgt * micro%voxel(ip1, ip2, ip3)%wgtc
                    sum_gamma = sum_gamma + micro%voxel(ip1, ip2, ip3)%gacumgr
                    nsolid = nsolid + 1.0d0

                    nsyst_loc = props%phase(micro%voxel(ip1, ip2, ip3)%phase)%nsyst
                    do is = 1, nsyst_loc
                        gsys = 0.5d0 * (micro%voxel(ip1, ip2, ip3)%crss(is, 1) + &
                                        micro%voxel(ip1, ip2, ip3)%crss(is, 2))
                        sum_g = sum_g + gsys
                        nsys_solid = nsys_solid + 1.0d0
                        sum_gdot = sum_gdot + abs(micro%voxel(ip1, ip2, ip3)%gamdot(is))
                        sum_ggdot = sum_ggdot + abs(micro%voxel(ip1, ip2, ip3)%gamdot(is)) * gsys
                    end do
                end if
            end do
            end do
            end do
            micro%epavg = micro%epavg / micro%wphcsol
            micro%eavg = micro%eavg / micro%wphcsol
            micro%eelavg = micro%eelavg / micro%wphcsol
            micro%edotpavg = micro%edotpavg / micro%wphcsol
            micro%sgsolavg = micro%sgsolavg / micro%wphcsol
            if (nsolid > 0.0d0) then
                gammaavg = sum_gamma / nsolid
            else
                gammaavg = 0.0d0
            end if
            if (nsys_solid > 0.0d0) then
                gavg = sum_g / nsys_solid
            else
                gavg = 0.0d0
            end if
            if (sum_gdot > 0.0d0) then
                geffavg = sum_ggdot / sum_gdot
            else
                geffavg = gavg
            end if
            micro%evmpavg = vm_strain(micro%epavg)
            micro%dvmpavg = vm_strain(micro%edotpavg)

            ! Write vtk fields
            if (props%process(ipr)%iwfields == 1 .and. &
                (mod(imicro, abs(props%process(ipr)%iwstep)) == 0 .or. &
                 (props%process(ipr)%iwstep < 0 .and. imicro == props%process(ipr)%nsteps))) then

                call order_sys_per_activity(micro, props)

                ! Write to vtk
                call write_vtk_cell(micro, props, imacroloop, output_dir)

            end if

            call update_diagnostics(micro, diag, finalize=.true., scan_cond=.true.)
            diag%erre = erre
            diag%errs = errs
            diag%errsbc = errsbc
            diag%tdotref = props%process(ipr)%tdotref
            diag%s33_target = props%process(ipr)%scauchy(3, 3)
            diag%s33 = micro%sgavg(3, 3)
            diag%s33_sol = micro%sgsolavg(3, 3)
            tdot_used = props%process(ipr)%tdot
            if (props%process(ipr)%ictrl == -1) then
                proc_time_done = proc_time_done + tdot_used
                props%process(ipr)%tdot = tdot_plan
            end if
            diag%devm_inc_max = micro%edotvmmx * tdot_used
            sim_time = sim_time + tdot_used
            wall_time = dble(time() - start)
            call diag%write_diagnostics(output_dir, iter, imacroloop, ipr, sim_time, wall_time, tdot_used)
            call write_results(output_dir, micro, iter, imacroloop, ipr, sim_time, wall_time, tdot_used, &
                               diag%tdotref, gammaavg, gavg, geffavg)

            ! ramp tdot back up
            if (props%process(ipr)%ictrl == -1) then
                if (tdot_plan < props%process(ipr)%tdotref .and. nconverged >= 5) then
                    nconverged = 0
                    tdot_plan = min(2.0d0 * tdot_plan, props%process(ipr)%tdotref)
                    remaining_time = max(0.0d0, proc_time_target - proc_time_done)
                    nsteps_left = ceiling(remaining_time / tdot_plan)
                    nsteps_total_new = imicro + nsteps_left
                    total_macrosteps = total_macrosteps + nsteps_total_new - props%process(ipr)%nsteps
                    props%process(ipr)%nsteps = nsteps_total_new
                    if (remaining_time > 0.0d0) then
                        props%process(ipr)%tdot = min(tdot_plan, remaining_time)
                    else
                        props%process(ipr)%tdot = tdot_plan
                    end if
                end if
            else if (props%process(ipr)%tdot < props%process(ipr)%tdotref .and. nconverged >= 5) then
                nconverged = 0
                tdot = props%process(ipr)%tdot
                props%process(ipr)%tdot = min(2.0 * props%process(ipr)%tdot, props%process(ipr)%tdotref)
                nsteps_left = props%process(ipr)%nsteps - imicro
                nsteps_left = ceiling(dble(nsteps_left) * tdot / props%process(ipr)%tdot)
                total_macrosteps = total_macrosteps + nsteps_left - (props%process(ipr)%nsteps - imicro)
                props%process(ipr)%nsteps = imicro + nsteps_left
            end if

            if (did_converge) nconverged = nconverged + 1
            call calc_c0(micro, props)

            close (21)

            imacroloop = imacroloop + 1
        end do
    end do

    call write_status(imacroloop, sim_time, wall_time, ipr, imicro, props%process(ipr)%tdot, iter, pct, final=.true.)

    write (*, *)
    call exit

end program fft
