module IO_functions
    implicit none

contains

    subroutine load_input(micro, props, nthreads, fft_input_file)
        use types
        use tensor_functions
        use fourier_functions
        use various_functions
        implicit none

        type(micro_type), intent(out) :: micro ! microstructure state (voxel fields, averages)
        type(props_type), intent(out) :: props ! simulation inputs/material/process definitions
        integer, intent(inout) :: nthreads ! thread count for FFTW/OMP setup
        integer :: ii ! x-index (grid) / scratch
        integer :: jj ! y-index (grid) / scratch
        integer :: kk ! z-index (grid) / scratch
        integer :: ijv(6, 2) ! Voigt->tensor index map: (11,22,33,23,13,12)
        integer :: m ! implied-do index for ijv DATA statement
        integer :: n ! implied-do index for ijv DATA statement
        integer :: iph ! phase index
        integer :: i ! generic loop index
        integer :: j ! generic loop index
        integer :: k ! generic loop index
        integer :: i1 ! component loop index (1..3)
        integer :: ig ! grain index (thermo/eigenstrain input)
        integer :: j1 ! component loop index (1..3)
        integer :: ngrains ! number of grains in thermo file
        integer :: grain_nvox(ngrmx, nphmx) ! solid voxel counts per (grain, phase)
        data((ijv(n, m), m=1, 2), n=1, 6) / 1, 1, 2, 2, 3, 3, 2, 3, 1, 3, 1, 2/ ! (11,22,33,23,13,12)
        character(len=256) :: prosa ! scratch line (headers/comments) from input files
        character(len=256) :: filecryspl ! slip-system/hardening definition file
        character(len=256) :: filecrysel ! elastic stiffness definition file
        character(len=256) :: filetext ! texture/microstructure file
        character(len=256) :: filethermo ! per-grain eigenstrain/thermal strain file
        character(len=256) :: filedispl ! initial displacement field file
        double precision :: cosalpha ! |cos| between slip direction and sample y-axis
        double precision :: dvm ! von Mises strain-rate of imposed symmetric part
        double precision, allocatable :: eth(:, :, :) ! eigenstrain tensor per grain (3x3xngrains)
        double precision :: dum ! scratch scalar
        double precision :: dbsa(3) ! slip direction in sample axes (for directional HP effect)
        double precision :: grain_d ! equivalent grain size for current solid voxel
        double precision :: grain_dmin ! lower cap on equivalent grain size
        double precision :: grain_size(ngrmx, nphmx) ! equivalent grain size per (grain, phase)
        double precision :: hp_kgr_default ! startup Hall-Petch coefficient (stress*sqrt(length))
        double precision :: hp_shift ! startup Hall-Petch shift added to tau0
        double precision :: voxel_vol ! physical voxel volume from delt
        character(len=256), intent(in) :: fft_input_file ! input fft file

        open (unit=11, file=fft_input_file, status='old')
        read (11, *) props%nph
        read (11, *) props%npts1, props%npts2, props%npts3

        props%wgt = 1.0 / float(props%npts1 * props%npts2 * props%npts3)

!$      call setup_openmp(nthreads)
        call initialize_fftw(props%npts1, props%npts2, props%npts3, nthreads)

        call init_micro(micro, props)

        read (11, *) props%delt(1), props%delt(2), props%delt(3)

        props%delt(1) = 4
        props%delt(2) = 4
        props%delt(3) = 4

        read (11, '(a)') prosa
        read (11, '(a)') filetext

        allocate (props%phase(props%nph))

        ! loop over phases
        do iph = 1, props%nph
            read (11, '(a)') prosa

            ! read igas integer to logical
            read (11, *) ii
            if (ii == 0) then
                props%phase(iph)%igas = .false.
            else
                props%phase(iph)%igas = .true.
            end if

            read (11, '(a)') prosa
            read (11, '(a)') filecryspl
            read (11, '(a)') filecrysel

            ! reads slip and twinning modes for the phase
            if (.not. props%phase(iph)%igas) then
                open (unit=12, file=filecryspl, status='old')
                call data_crystal(props, iph)
                close (unit=12)

                open (unit=12, file=filecrysel, status='old')
                call data_crystal_elast(props, iph)
                close (unit=12)
            end if
        end do ! end of data input loop over all phases

        ! reads initial texture from filetext
        ! and calculates the local elastic stiffness
        open (unit=12, file=filetext, status='old')
        call data_grain(props, micro)
        close (unit=12)

        ! Startup grain-size Hall-Petch shift applied to tau0 only.
        ! Keep the default coefficient at zero unless you explicitly want
        ! grain-wise strengthening without changing the input file format.
        hp_kgr_default = 632 / 2.5 !200.0d0
        voxel_vol = props%delt(1) * props%delt(2) * props%delt(3)
        grain_dmin = 3.0d0 * voxel_vol**(1.0d0 / 3.0d0)
        grain_nvox = 0
        grain_size = 0.0d0

        do ii = 1, props%npts1
        do jj = 1, props%npts2
        do kk = 1, props%npts3
            iph = micro%voxel(ii, jj, kk)%phase
            if (.not. micro%voxel(ii, jj, kk)%gas) then
                ig = micro%voxel(ii, jj, kk)%grain
                if (ig < 1 .or. ig > ngrmx) then
                    write (*, '('' grain id out of bounds for Hall-Petch init: '',i8)') ig
                    stop
                end if
                grain_nvox(ig, iph) = grain_nvox(ig, iph) + 1
            end if
        end do
        end do
        end do

        do iph = 1, props%nph
        do ig = 1, ngrmx
            if (grain_nvox(ig, iph) > 0) then
                grain_size(ig, iph) = (dble(grain_nvox(ig, iph)) * voxel_vol)**(1.0d0 / 3.0d0)
            end if
        end do
        end do

        ! restart analysis
        read (11, '(a)') prosa
        read (11, *) i, props%imicrosave
        read (11, *) props%complgas
        props%res = .false.
        if (i == 1) props%res = .true.

        ! Initial displacement field (ITDISPINI & FILEDISPL)
        read (11, '(a)') prosa
        read (11, *) props%idispini
        if (props%idispini > 0) read (11, '(A)') filedispl

        if (props%idispini > 0) then
            call init_defgrad(filedispl, micro, props)
        end if

        do ii = 1, props%npts1
        do jj = 1, props%npts2
        do kk = 1, props%npts3
            micro%voxel(ii, jj, kk)%gacumgr = 0.0
            iph = micro%voxel(ii, jj, kk)%phase
            if (.not. micro%voxel(ii, jj, kk)%gas) then
                hp_shift = 0.0d0
                if (hp_kgr_default > 0.0d0) then
                    ig = micro%voxel(ii, jj, kk)%grain
                    grain_d = grain_size(ig, iph)
                    if (grain_d > 0.0d0) then
                        grain_d = max(grain_d, grain_dmin)
                        hp_shift = hp_kgr_default / sqrt(grain_d)
                    end if
                end if

                do i = 1, props%phase(iph)%nsyst

                    micro%voxel(ii, jj, kk)%crss(i, 1) = props%phase(iph)%tau(i, 1) + hp_shift
                    micro%voxel(ii, jj, kk)%crss(i, 2) = props%phase(iph)%tau(i, 2) + hp_shift

                    micro%voxel(ii, jj, kk)%trialtau(i, 1) = props%phase(iph)%tau(i, 1) + hp_shift
                    micro%voxel(ii, jj, kk)%trialtau(i, 2) = props%phase(iph)%tau(i, 2) + hp_shift

                    micro%voxel(ii, jj, kk)%kin(i) = 0.0

                    ! DIRECTIONAL HP EFFECT
                    do i1 = 1, 3
                        dum = 0.0
                        do j1 = 1, 3
                            dum = dum + micro%voxel(ii, jj, kk)%AG(i1, j1) * props%phase(iph)%dbca(j1, i)
                        end do
                        dbsa(i1) = dum
                    end do
                    cosalpha = abs(dbsa(2)) ! cos with y axis
                    micro%voxel(ii, jj, kk)%crss(i, 1) = micro%voxel(ii, jj, kk)%crss(i, 1) + props%phase(iph)%hpfac(i) * sqrt(cosalpha)
                    micro%voxel(ii, jj, kk)%crss(i, 2) = micro%voxel(ii, jj, kk)%crss(i, 2) + props%phase(iph)%hpfac(i) * sqrt(cosalpha)
                    micro%voxel(ii, jj, kk)%trialtau(i, 1) = micro%voxel(ii, jj, kk)%crss(i, 2)
                    micro%voxel(ii, jj, kk)%trialtau(i, 2) = micro%voxel(ii, jj, kk)%crss(i, 2)
                end do
            end if
        end do
        end do
        end do

        ! read number of processes and begin process loop
        read (11, *) props%nproc
        if (props%nproc < 1) stop 'cannot have less than one processes'

        allocate (props%process(props%nproc))

        do ipr = 1, props%nproc

            ! read boundary conditions on overall stress and strain-rate
            read (11, '(a)') prosa
            read (11, '(a)') prosa

            do i = 1, 3
                read (11, *) (props%process(ipr)%iudot(i, j), j=1, 3)
            end do

            if (props%process(ipr)%iudot(1, 1) + props%process(ipr)%iudot(2, 2) + props%process(ipr)%iudot(3, 3) == 2) then
                write (*, '(a)') 'check diagonal boundary conditions iudot'
                write (*, '(a)') 'cannot enforce only two deviatoric components'
                stop
            end if

            do i = 1, 3
            do j = 1, 3
            if (i /= j .and. props%process(ipr)%iudot(i, j) + props%process(ipr)%iudot(j, i) == 0) then
                write (*, '(a)') 'check off-diagonal boundary conditions iudot'
                stop
            end if
            end do
            end do

            read (11, *)
            do i = 1, 3
                read (11, *) (props%process(ipr)%udot(i, j), j=1, 3)
            end do

            ! symmetric strain-rate, antisymmetric rotation-rate tensors
            ! and indices of imposed components
            props%process(ipr)%dsim = sym33(props%process(ipr)%udot)
            props%process(ipr)%tomtot = antisym33(props%process(ipr)%udot)
            dvm = vm_strain(props%process(ipr)%dsim)

            do i = 1, 3
                props%process(ipr)%idsim(i) = props%process(ipr)%iudot(i, i)
            end do

            props%process(ipr)%idsim(4) = 0
            if (props%process(ipr)%iudot(2, 3) == 1 .and. props%process(ipr)%iudot(3, 2) == 1) props%process(ipr)%idsim(4) = 1
            props%process(ipr)%idsim(5) = 0
            if (props%process(ipr)%iudot(1, 3) == 1 .and. props%process(ipr)%iudot(3, 1) == 1) props%process(ipr)%idsim(5) = 1
            props%process(ipr)%idsim(6) = 0
            if (props%process(ipr)%iudot(1, 2) == 1 .and. props%process(ipr)%iudot(2, 1) == 1) props%process(ipr)%idsim(6) = 1

            read (11, *)
            read (11, *) props%process(ipr)%iscau(1), props%process(ipr)%iscau(6), props%process(ipr)%iscau(5)
            read (11, *) props%process(ipr)%iscau(2), props%process(ipr)%iscau(4)
            read (11, *) props%process(ipr)%iscau(3)

            do i = 1, 6
                if (props%process(ipr)%iscau(i) * props%process(ipr)%idsim(i) /= 0 .or. &
                    props%process(ipr)%iscau(i) + props%process(ipr)%idsim(i) /= 1) then
                    write (*, '(a)') ' check boundary condits on strain-rate and stress'
                    write (*, '('' idsim = '',6i3)') props%process(ipr)%idsim
                    write (*, '('' iscau = '',6i3)') props%process(ipr)%iscau
                    stop
                end if
            end do

            read (11, *)
            read (11, *) props%process(ipr)%scauchy(1, 1), props%process(ipr)%scauchy(1, 2), props%process(ipr)%scauchy(1, 3)
            read (11, *) props%process(ipr)%scauchy(2, 2), props%process(ipr)%scauchy(2, 3)
            read (11, *) props%process(ipr)%scauchy(3, 3)
            props%process(ipr)%scauchy(3, 2) = props%process(ipr)%scauchy(2, 3)
            props%process(ipr)%scauchy(3, 1) = props%process(ipr)%scauchy(1, 3)
            props%process(ipr)%scauchy(2, 1) = props%process(ipr)%scauchy(1, 2)

            read (11, *)
            read (11, *) props%process(ipr)%tdot, props%process(ipr)%tdotref, props%process(ipr)%tdotmin, props%process(ipr)%devmmx
            read (11, *) props%process(ipr)%ictrl

            if (props%process(ipr)%ictrl == 0) then
                props%process(ipr)%devmapp = props%process(ipr)%tdot
                props%process(ipr)%tdot = props%process(ipr)%tdot / dvm
                props%process(ipr)%tdotref = props%process(ipr)%tdot
            else if (props%process(ipr)%ictrl > 0) then
                write (*, '(a)') 'ictrl>0 not implemented'
                stop
            end if

            read (11, *)
            read (11, *) props%process(ipr)%nsteps
            read (11, *) props%process(ipr)%error, props%process(ipr)%erroral
            read (11, *) props%process(ipr)%itmax
            read (11, *) props%process(ipr)%irecover
            if (props%process(ipr)%irecover == 1) open (50, file='stress.in', status='old', access='sequential', form='unformatted')
            read (11, *) props%process(ipr)%isave
            read (11, *) props%process(ipr)%iupdate
            read (11, *) props%process(ipr)%iuphard
            read (11, *) props%process(ipr)%iwtex
            read (11, *) props%process(ipr)%iwfields, props%process(ipr)%iwstep
            read (11, *) props%process(ipr)%xc0, props%process(ipr)%igamma

            ! parameter requirements
            if (.not. (props%process(ipr)%nsteps > 0)) stop 'nsteps must be positive'
            if (.not. (props%process(ipr)%error > 0 .or. props%process(ipr)%erroral > 0)) stop 'error and erroral tolerances must be positive'
            if (.not. (props%process(ipr)%itmax > 0)) stop 'itmax must be positive'
            if (.not. (props%process(ipr)%irecover == 0 .or. props%process(ipr)%irecover == 1)) stop 'irecover must be 0 or 1'
            if (.not. (props%process(ipr)%isave == 0 .or. props%process(ipr)%isave == 1)) stop 'isave must be 0 or 1'
            if (.not. (props%process(ipr)%iupdate == 0 .or. props%process(ipr)%iupdate == 1)) stop 'iupdate must be 0 or 1'
            if (.not. (props%process(ipr)%iuphard == 0 .or. props%process(ipr)%iuphard == 1)) stop 'iuphard must be 0 or 1'
            if (.not. (props%process(ipr)%iwtex == 0 .or. props%process(ipr)%iwtex == 1)) stop 'iwtex must be 0 or 1'
            if (.not. (props%process(ipr)%iwfields == 0 .or. props%process(ipr)%iwfields == 1)) stop 'iwfields must be 0 or 1'
            if (.not. (props%process(ipr)%igamma == 0 .or. props%process(ipr)%igamma == 1 .or. props%process(ipr)%igamma == 2)) &
                stop 'igamma must be 0 (continuous), 1 (discrete), or 2 (modified)'

            ! initial eigenstrain field (ithermo & filethermo)
            read (11, *) props%process(ipr)%ithermo
            if (ipr == 1) allocate (eth(3, 3, ngrmx))
            if (props%process(ipr)%ithermo == 1) then

                read (11, '(a)') filethermo
                open (49, file=filethermo, status='old')
                read (49, *) ngrains

                if (ngrains > ngrmx) then
                    write (*, '(a, i0, a, i0)') 'more grains than max'
                    stop
                end if

                do ig = 1, ngrains
                    read (49, *) eth(1, 1, ig), eth(2, 2, ig), eth(3, 3, ig), eth(2, 3, ig), eth(3, 1, ig), eth(1, 2, ig)
                    eth(3, 2, ig) = eth(2, 3, ig)
                    eth(1, 3, ig) = eth(3, 1, ig)
                    eth(2, 1, ig) = eth(1, 2, ig)
                end do
            end if

            do k = 1, props%npts3
            do j = 1, props%npts2
            do i = 1, props%npts1
                micro%voxel(i, j, k)%eth(:, :) = eth(:, :, micro%voxel(i, j, k)%grain)
            end do
            end do
            end do

        end do ! nproc

        props%wphsol = 0.0
        props%wphgas = 0.0
        do k = 1, props%npts3
        do j = 1, props%npts2
        do i = 1, props%npts1
            if (micro%voxel(i, j, k)%gas) props%wphgas = props%wphgas + props%wgt
            if (.not. micro%voxel(i, j, k)%gas) props%wphsol = props%wphsol + props%wgt
        end do
        end do
        end do
        if (props%idispini == 0) then
            micro%wphcgas = props%wphgas
            micro%wphcsol = props%wphsol
        end if

        return

    end subroutine load_input

    subroutine data_crystal(props, iph)
        use types
        use tensor_functions
        implicit none

        type(props_type), intent(inout) :: props ! simulation inputs/material data
        integer, intent(in) :: iph ! phase index

        double precision :: sn(3) ! slip plane normal (unnormalized)
        double precision :: sb(3) ! slip direction / Burgers vector (unnormalized)
        double precision :: cdim(3) ! crystal lattice parameters from file
        double precision :: aux5(5) ! scratch: 5-comp deviatoric basis
        double precision :: aux33(3, 3) ! scratch: 3x3 tensor
        double precision :: aux55(5, 5) ! scratch: 5x5 matrix
        double precision :: aux3333(3, 3, 3, 3) ! scratch: 4th-order tensor
        double precision :: hselfx(10) ! self hardening per mode (file order)
        double precision :: hlatex(10, 10) ! latent hardening matrix (file order)
        double precision :: covera ! c/a ratio for HEX/TRI systems
        double precision :: gamd0x ! reference shear rate for current mode
        double precision :: hpfacx ! Hall-Petch factor for current mode
        double precision :: prod ! dot(n,b) orthogonality check
        double precision :: qnor ! |sb| norm
        double precision :: snor ! |sn| norm
        double precision :: tau0xb ! initial CRSS (backward direction)
        double precision :: tau0xf ! initial CRSS (forward direction)
        double precision :: tau1x ! saturation CRSS parameter (Voce tau1)
        double precision :: thet0x ! initial hardening rate parameter
        double precision :: thet1x ! final hardening rate parameter
        double precision :: twshx ! twinning shear (0 => slip mode)
        integer :: i ! generic loop index
        integer :: j ! generic loop index
        integer :: k ! generic loop index
        integer :: im ! active mode index (1..nmodes)
        integer :: is ! slip system index
        integer :: isectwx ! twinning section id from file
        integer :: iz ! skip counter for unused mode lines
        integer :: jm ! mode index (latent hardening)
        integer :: js ! system index within selected mode
        integer :: kount ! counter over selected (active) modes
        integer :: m ! component index (1..3)
        integer :: modex ! mode id read from crystal file
        integer :: nm ! loop index over modes in crystal file
        integer :: nmodesx ! number of modes defined in file
        double precision :: nrsx ! rate sensitivity exponent for current mode
        integer :: nsmx ! number of systems in current mode
        integer :: nsysx ! global slip system counter for the phase
        integer :: nmodes ! number of active modes selected
        integer :: isn(100, 4) ! slip plane indices from file (Miller/Miller-Bravais)
        integer :: isb(100, 4) ! slip direction indices from file (Miller/Miller-Bravais)
        integer :: mode(100) ! list of active mode ids
        integer :: nsyst ! total slip/twin system count for phase
        integer :: ntwmod ! number of twinning modes
        integer :: ntwsys ! number of twinning systems
        character(len=256) :: prosa ! scratch line (headers/comments)
        character(len=3) :: icryst ! crystal symmetry tag (HEX/CUB/...)

        read (12, '(a)') prosa
        read (12, '(a)') icryst
        read (12, *) (cdim(i), i=1, 3)
        covera = cdim(3) / cdim(1)
        read (12, *) nmodesx
        read (12, *) nmodes
        read (12, *) (mode(i), i=1, nmodes)
        props%phase(iph)%nmodes = nmodes

        if (nmodes > nmodmx) then
            write (*, '('' nmodes in phase'',i3,'' is'',i3)') iph, nmodes
            write (*, '('' change parameter nmodmx'')')
            stop
        end if

        ntwmod = 0
        ntwsys = 0
        nsyst = 0
        kount = 1

        ! start reading deformation modes from filecrys
        do nm = 1, nmodesx

            read (12, '(a)') prosa
            read (12, *) modex, nsmx, nrsx, gamd0x, twshx, isectwx
            read (12, *) tau0xf, tau0xb, tau1x, thet0x, thet1x, hpfacx
            read (12, *) hselfx(nm), (hlatex(nm, jm), jm=1, nmodesx)

            ! skips nsmx lines if the mode is not in the list.
            if (modex /= mode(kount)) then

                do iz = 1, nsmx
                    read (12, *)
                end do

            else

                if (thet0x < thet1x) then
                    write (*, '('' initial hardening lower than final hardening for mode'',i3,''  in phase'',i3)') kount, iph
                    stop
                end if

                ! case tau1=0 corresponds to linear hardening and is independent of tau0.
                ! avoid division by zero
                if (tau1x <= 1.e-6) then
                    tau1x = 1.e-6
                    thet0x = thet1x
                end if

                ! reorder hardening coefficients to account only for active modes
                hselfx(kount) = hselfx(nm)
                do i = 1, nmodes
                    hlatex(kount, i) = hlatex(nm, mode(i))
                end do

                ! systems given in four index notation: hexagonals and trigonals
                ! systems given in three index notation: cubic and orthorhombic
                if (icryst == 'HEX' .or. icryst == 'TRI' .or. icryst == 'hex' .or. icryst == 'tri') then
                    do j = 1, nsmx
                        read (12, *) (isn(j, k), k=1, 4), (isb(j, k), k=1, 4)
                    end do
                else if (icryst == 'CUB' .or. icryst == 'ORT' .or. icryst == 'cub' .or. icryst == 'ort') then
                    do j = 1, nsmx
                        read (12, *) (isn(j, k), k=1, 3), (isb(j, k), k=1, 3)
                    end do
                else
                    write (*, '('' cannot identify the crystal symmetry of phase '', i3)') iph
                    stop
                end if

                props%phase(iph)%nsm(kount) = nsmx
                if (twshx > 1.0d-15) ntwmod = ntwmod + 1

                if (ntwmod > ntwmmx) then
                    write (*, '('' ntwmod in phase'',i3,'' is'',i3)') iph, ntwmod
                    write (*, '('' change parameter ntwmmx in vpsc.dim'')')
                    stop
                end if

                do js = 1, props%phase(iph)%nsm(kount)

                    nsyst = nsyst + 1
                    nsysx = nsyst
                    if (twshx > 1.0d-15) ntwsys = ntwsys + 1

                    if (nsyst > nsysmx) then
                        write (*, '('' nsyst in phase'',i3,'' is'',i3)') iph, nsyst
                        write (*, '('' change parameter nsysmx in vpsc.dim'')')
                        write (*, *) twshx
                        stop
                    end if

                    ! defines rate sensitivity and crss for each system in the mode
                    props%phase(iph)%gamd0(nsysx) = gamd0x
                    props%phase(iph)%nrs(nsysx) = nrsx
                    props%phase(iph)%twsh(nsysx) = twshx
                    props%phase(iph)%tau(nsysx, 1) = tau0xf
                    props%phase(iph)%tau(nsysx, 2) = tau0xb
                    props%phase(iph)%tau(nsysx, 3) = tau1x
                    props%phase(iph)%thet(nsysx, 0) = thet0x
                    props%phase(iph)%thet(nsysx, 1) = thet1x
                    props%phase(iph)%hpfac(nsysx) = hpfacx

                    props%phase(iph)%isectw(nsysx) = isectwx

                    if (icryst == 'HEX' .or. icryst == 'TRI' .or. icryst == 'hex' .or. icryst == 'tri') then
                        sn(1) = float(isn(js, 1))
                        sn(2) = (float(isn(js, 1)) + 2.*float(isn(js, 2))) / sqrt(3.)
                        sn(3) = float(isn(js, 4)) / covera
                        sb(1) = 3./2.*float(isb(js, 1))
                        sb(2) = (float(isb(js, 1)) / 2.+float(isb(js, 2))) * sqrt(3.)
                        sb(3) = float(isb(js, 4)) * covera
                    else if (icryst == 'CUB' .or. icryst == 'ORT' .or. icryst == 'cub' .or. icryst == 'ort') then
                        do m = 1, 3
                            sn(m) = float(isn(js, m)) / cdim(m)
                            sb(m) = float(isb(js, m)) * cdim(m)
                        end do
                    end if

                    ! *** normalizes system vectors and checks normality
                    snor = sqrt(sn(1) * sn(1) + sn(2) * sn(2) + sn(3) * sn(3))
                    qnor = sqrt(sb(1) * sb(1) + sb(2) * sb(2) + sb(3) * sb(3))
                    prod = 0.
                    do j = 1, 3
                        props%phase(iph)%dnca(j, nsysx) = sn(j) / snor
                        props%phase(iph)%dbca(j, nsysx) = sb(j) / qnor
                        if (abs(props%phase(iph)%dnca(j, nsysx)) < 1.e-03) props%phase(iph)%dnca(j, nsysx) = 0.
                        if (abs(props%phase(iph)%dbca(j, nsysx)) < 1.e-03) props%phase(iph)%dbca(j, nsysx) = 0.
                        prod = prod + props%phase(iph)%dnca(j, nsysx) * props%phase(iph)%dbca(j, nsysx)
                    end do
                    if (prod >= 1.e-3) then
                        write (*, '(''system'',i4,''  in mode'',i4,'' in phase'',i4,''  is not orthogonal !!'')') js, nm, iph
                        stop
                    end if

                    ! define schmid vector in crystal axes for each system
                    do i = 1, 3
                    do j = 1, 3
                        aux33(i, j) = (props%phase(iph)%dnca(i, nsysx) * props%phase(iph)%dbca(j, nsysx) + &
                                       props%phase(iph)%dnca(j, nsysx) * props%phase(iph)%dbca(i, nsysx)) / 2.
                    end do
                    end do

                    call chg_basis(aux5, aux33, aux55, aux3333, 2, 5)

                    do i = 1, 5
                        props%phase(iph)%schca(i, nsysx) = aux5(i)
                    end do

                end do    ! end of loop over a given deformation mode
                props%phase(iph)%nsyst = nsyst
                kount = kount + 1

            end if

        end do    ! end of loop over all modes in a given phase

        ! initialize self & latent hardening coefs for each system of the phase.
        ! absolute units are accounted for by modulating factor in hardening law.

        i = 0
        do im = 1, nmodes
        do is = 1, props%phase(iph)%nsm(im)
            i = i + 1
            j = 0
            do jm = 1, nmodes
            do js = 1, props%phase(iph)%nsm(jm)
                j = j + 1
                props%phase(iph)%hard(i, j) = hlatex(im, jm)
            end do
            end do
            props%phase(iph)%hard(i, i) = hselfx(im)
        end do
        end do

        !c      do i=1,2
        !c      write(*,*) (hard(i,j,iph),j=1,2)
        !c      enddo
        !c      pause

        !     latent hardening of slip and twinning by twinning is based on the
        !     relative directions of the shear direction and the twin plane
        !     this approach is still being tested (30/4/99)

        !     nslsys=nsyst(iph)-ntwsys(iph)
        !     do is=1,nsyst(iph)
        !       do jt=nslsys+1,nsyst(iph)
        !         if(is.ne.jt) then
        !           cosa=dbca(1,is,iph)*dnca(1,jt,iph)+
        !    #           dbca(2,is,iph)*dnca(2,jt,iph)+
        !    #           dbca(3,is,iph)*dnca(3,jt,iph)
        !           cosa=abs(cosa)
        !           hard(is,jt,iph)=hard(is,jt,iph)*(0.5+1.0*cosa)
        !         endif
        !       enddo
        !     enddo
        !
        !     write(10,'(''  hardening matrix for phase'',i3)') iph
        !     do i=1,nsyst(iph)
        !       write(10,'(24f5.1)') (hard(i,j,iph),j=1,nsyst(iph))
        !     enddo

        !     verification of twinning data to be sure program will run properly

        if (nmodes > 1) then
            do i = 2, nsyst
                if (props%phase(iph)%twsh(i) < 1.0d-15 .and. props%phase(iph)%twsh(i - 1) > 1.0d-15) then
                    write (*, '(a)') ' warning! the twinning modes must follow the'
                    write (*, '(a)') ' slip modes   -->   reorder crystal file'
                    stop
                end if
            end do
        end if

        return
    end

    subroutine data_crystal_elast(props, iph)
        use types
        use tensor_functions
        implicit none

        type(props_type), intent(inout) :: props ! simulation inputs/material data
        integer, intent(in) :: iph ! phase index
        double precision :: dde(3, 3) ! 3x3 identity tensor
        double precision :: xid4(3, 3, 3, 3) ! symmetric 4th-order identity tensor
        double precision :: cc66v(6, 6) ! stiffness in Voigt notation (from file)
        double precision :: ccaux(3, 3, 3, 3) ! stiffness as 4th-order tensor
        double precision :: aux6(6) ! scratch: Voigt vector/tensor
        double precision :: aux33(3, 3) ! scratch: 3x3 tensor

        integer :: i ! loop index
        integer :: j ! loop index
        integer :: k ! loop index
        integer :: l ! loop index
        integer :: iso ! flag: 0=anisotropic stiffness, else isotropic
        double precision :: tla ! Lamé parameter lambda (isotropic)
        double precision :: tmu ! shear modulus mu (isotropic)
        double precision :: tnu ! Poisson ratio nu (isotropic)
        double precision :: young ! Young's modulus (isotropic)

        ! unitary tensors
        do i = 1, 3
        do j = 1, 3
            dde(i, j) = 0.d0
            if (i == j) dde(i, j) = 1.d0
        end do
        end do

        do i = 1, 3
        do j = 1, 3
        do k = 1, 3
        do l = 1, 3
            xid4(i, j, k, l) = (dde(i, k) * dde(j, l) + dde(i, l) * dde(j, k)) / 2.d0
        end do
        end do
        end do
        end do

        read (12, *) iso

        if (iso == 0) then

            do i = 1, 6
                read (12, *) (cc66v(i, j), j=1, 6)
            end do

            call voigt_vpsc(aux6, aux33, cc66v, ccaux, 3)
            do i = 1, 3
            do j = 1, 3
            do k = 1, 3
            do l = 1, 3
                props%phase(iph)%C(i, j, k, l) = ccaux(i, j, k, l)
            end do
            end do
            end do
            end do

        else

            read (12, *) young, tnu
            tmu = young / (2.d0 * (1.+tnu))
            tla = 2.d0 * tmu * tnu / (1.d0 - 2.d0 * tnu)

            do i = 1, 3
            do j = 1, 3
            do k = 1, 3
            do l = 1, 3
                props%phase(iph)%C(i, j, k, l) = tla * dde(i, j) * dde(k, l) + 2.d0 * tmu * xid4(i, j, k, l)
            end do
            end do
            end do
            end do

        end if

        return
    end

    subroutine data_grain(props, micro)
        use types
        use tensor_functions
        implicit none

        type(props_type), intent(inout) :: props ! simulation inputs/material data
        type(micro_type), intent(inout) :: micro ! microstructure grid (voxels)

        double precision :: aa(3, 3) ! sample->crystal rotation matrix
        double precision :: caux3333(3, 3, 3, 3) ! rotated stiffness tensor (4th order)
        double precision :: caux66(6, 6) ! rotated stiffness in Voigt notation
        double precision :: aux6(6) ! scratch: Voigt vector/tensor
        double precision :: aux33(3, 3) ! scratch: 3x3 tensor

        double precision :: dum ! scratch scalar
        double precision :: om ! Euler angle omega (deg)
        double precision :: ph ! Euler angle phi (deg)
        double precision :: th ! Euler angle theta (deg)
        integer :: ii ! voxel x-index from texture file
        integer :: jj ! voxel y-index from texture file
        integer :: kk ! voxel z-index from texture file
        integer :: i1 ! tensor index loop (1..3)
        integer :: i2 ! tensor index loop (1..3)
        integer :: j ! loop index
        integer :: j1 ! tensor index loop (1..3)
        integer :: j2 ! tensor index loop (1..3)
        integer :: jgr ! grain id from texture file
        integer :: jph ! phase id from texture file
        integer :: k ! loop index
        integer :: k1 ! tensor index loop (1..3)
        integer :: k2 ! tensor index loop (1..3)
        integer :: kkk ! voxel counter (1..npts1*npts2*npts3)
        integer :: l1 ! tensor index loop (1..3)
        integer :: l2 ! tensor index loop (1..3)
        integer :: nph1 ! number of voxels in phase 1 (diagnostic)

        nph1 = 0
        do kkk = 1, props%npts1 * props%npts2 * props%npts3

            read (12, *) ph, th, om, ii, jj, kk, jgr, jph

            if (jph == 1) nph1 = nph1 + 1

            micro%voxel(ii, jj, kk)%grain = jgr
            micro%voxel(ii, jj, kk)%phase = jph
            micro%voxel(ii, jj, kk)%gas = props%phase(jph)%igas

            if (props%phase(jph)%igas) cycle

            ! calculates the transformation matrix aa which transforms from
            ! sample to crystal. stores ag, which transforms from crystal to sample.
            call euler(2, ph * pi / 180.0, th * pi / 180.0, om * pi / 180.0, aa)

            do j = 1, 3
            do k = 1, 3
                micro%voxel(ii, jj, kk)%ag(j, k) = aa(k, j)
            end do
            end do

            do i1 = 1, 3
            do j1 = 1, 3
            do k1 = 1, 3
            do l1 = 1, 3
                dum = 0.
                do i2 = 1, 3
                do j2 = 1, 3
                do k2 = 1, 3
                do l2 = 1, 3
                    dum = dum + aa(i2, i1) * aa(j2, j1) * aa(k2, k1) * aa(l2, l1) * props%phase(jph)%C(i2, j2, k2, l2)
                end do
                end do
                end do
                end do
                caux3333(i1, j1, k1, l1) = dum
            end do
            end do
            end do
            end do

            call chg_basis(aux6, aux33, caux66, caux3333, 4, 6)

            micro%voxel(ii, jj, kk)%cg66(:, :) = caux66(:, :)

            call lu_inverse(caux66, 6)
            micro%voxel(ii, jj, kk)%sg66(:, :) = caux66(:, :)
            micro%voxel(ii, jj, kk)%sg66(:, :) = caux66(:, :)

        end do

        return
    end

!**************************************************
    subroutine init_defgrad(filedispl, micro, props)
        use types
        use tensor_functions
        use fourier_functions
        use various_functions
        implicit none

        character(len=256) :: filedispl ! initial displacement input file name
        type(props_type), intent(in) :: props ! simulation inputs/material data
        type(micro_type), intent(inout) :: micro ! microstructure state to update

        double precision :: dispini(3, props%npts1, props%npts2, props%npts3) ! initial displacement at cell centers
        double precision :: dispini_node(3, props%npts1 + 1, props%npts2 + 1, props%npts3 + 1) ! initial displacement at nodes
        double precision :: detFavg ! average det(F) for normalization
        double precision :: defgrad(3, 3, props%npts1, props%npts2, props%npts3) ! displacement gradient / defgrad fluctuation
        double precision :: defgradinv(3, 3) ! scratch: inverse deformation gradient
        integer :: i ! loop index
        integer :: j ! loop index
        integer :: k ! loop index
        integer :: l ! component loop index (1..3)
        integer :: ip1 ! x-index read from file
        integer :: ip2 ! y-index read from file
        integer :: ip3 ! z-index read from file
        integer :: ideriv ! derivative operator selection for gradient
        integer :: npts3node ! number of nodes along z (1 or npts3+1)

        write (*, '(A)') 'Initializing eigendisplacement field'

        if (props%npts3 == 1) then
            npts3node = 1
        else
            npts3node = props%npts3 + 1
        end if

        ! load initial displacement
        if (props%idispini == 1) then

            open (12, file=FILEDISPL, status='unknown')
            do k = 1, props%npts3
            do j = 1, props%npts2
            do i = 1, props%npts1
                read (12, *) ip1, ip2, ip3, (dispini(l, i, j, k), l=1, 3)
            end do
            end do
            end do
            close (12)

            call cell2node_tensor1(dispini, dispini_node, 3, props%npts1, props%npts2, props%npts3)

            ideriv = 1
            call gradient_tensor1(dispini, defgrad, ideriv, props%npts1, props%npts2, props%npts3)
        elseif (props%idispini == 2) then

            open (12, file=FILEDISPL, status='unknown')

            do k = 1, npts3node
            do j = 1, props%npts2 + 1
            do i = 1, props%npts1 + 1
                read (12, *) ip1, ip2, ip3, (dispini_node(l, i, j, k), l=1, 3)
            end do
            end do
            end do
            close (12)

            dispini = dispini_node(:, 1:props%npts1, 1:props%npts2, 1:props%npts3)

            ideriv = 2
            call gradient_tensor1(dispini, defgrad, ideriv, props%npts1, props%npts2, props%npts3)
        end if

        micro%xnode = micro%xnode + dispini_node(:, :, :, 1:npts3node)

        micro%defgradavg(1, 1) = props%delt(1)
        micro%defgradavg(2, 2) = props%delt(2)
        micro%defgradavg(3, 3) = props%delt(3)

        micro%defgradinvavgc = 0.0
        micro%defgradinvavgcs = 0.0
        micro%defgradinvavgcg = 0.0
        micro%wphcsol = 0.0
        micro%wphcgas = 0.0
        detFavg = 0.0
        do k = 1, props%npts3
        do j = 1, props%npts2
        do i = 1, props%npts1

            micro%voxel(i, j, k)%defgrad(:, :) = defgrad(:, :, i, j, k) + micro%defgradavg
            micro%voxel(i, j, k)%defgradini(:, :) = micro%voxel(i, j, k)%defgrad(:, :)
            defgradinv(:, :) = micro%voxel(i, j, k)%defgrad(:, :)
            micro%voxel(i, j, k)%detF = determinant33(defgradinv)

            call lu_inverse(defgradinv, 3)

            micro%voxel(i, j, k)%defgradinv(:, :) = defgradinv(:, :)
            micro%voxel(i, j, k)%defgradiniinv(:, :) = micro%voxel(i, j, k)%defgradinv(:, :)
            micro%defgradinvavgc(:, :) = micro%defgradinvavgc(:, :) + defgradinv(:, :) * props%wgt * micro%voxel(i, j, k)%detF
            if (.not. micro%voxel(i, j, k)%gas) micro%defgradinvavgcs(:, :) = micro%defgradinvavgcs(:, :) + &
                                                                              defgradinv(:, :) * props%wgt * micro%voxel(i, j, k)%detF
            if (.not. micro%voxel(i, j, k)%gas) micro%wphcsol = micro%wphcsol + props%wgt * micro%voxel(i, j, k)%detF
            if (micro%voxel(i, j, k)%gas) micro%defgradinvavgcg(:, :) = micro%defgradinvavgcg(:, :) + &
                                                                        defgradinv(:, :) * props%wgt * micro%voxel(i, j, k)%detF
            if (micro%voxel(i, j, k)%gas) micro%wphcgas = micro%wphcgas + props%wgt * micro%voxel(i, j, k)%detF
            micro%voxel(i, j, k)%wgtc = props%wgt * micro%voxel(i, j, k)%detF

            detFavg = detFavg + micro%voxel(i, j, k)%detF * props%wgt

        end do
        end do
        end do

        micro%defgradinvavgc = micro%defgradinvavgc / detFavg
        micro%wphcsol = micro%wphcsol / detFavg
        micro%wphcgas = micro%wphcgas / detFavg
        micro%defgradinvavgcs = micro%defgradinvavgcs / detFavg / micro%wphcsol
        if (micro%wphcgas > 0.0) then
            micro%defgradinvavgcg = micro%defgradinvavgcg / detFavg / micro%wphcgas
        else
            micro%defgradinvavgcg = id3
        end if
        micro%voxel%wgtc = micro%voxel%wgtc / detFavg

    end

    double precision function clamp_subnormal(a) result(b)
        implicit none
        double precision, intent(in) :: a
        b = a
        if (abs(b) < tiny(1.0d0)) b = 0.0d0
    end function clamp_subnormal

    subroutine write_vtk_cell(micro, props, istp, output_dir)
        use types
        use tensor_functions
        implicit none

        type(props_type), intent(in) :: props ! simulation inputs/material data
        type(micro_type), intent(in) :: micro ! microstructure state (fields, nodes)
        integer, intent(in) :: istp ! step index for output naming
        character(len=256), intent(in) ::output_dir
        double precision :: xtmp(3) ! scratch coordinate vector
        double precision :: uf(3, props%npts1, props%npts2, props%npts3) ! displacement fluctuation at cell centers
        integer :: id(props%npts1 + 1, props%npts2 + 1, props%npts3 + 1) ! node id lookup for VTK connectivity
        integer :: nsystmx ! max slip systems across phases (VTK array width)
        integer :: i ! generic loop index (components/systems)
        integer :: ic ! running node id counter
        integer :: ictp3d ! VTK 3D cell ordering/type selector (11 or 12)
        integer :: kx ! x-index (cell/node)
        integer :: ky ! y-index (cell/node)
        integer :: kz ! z-index (cell/node)
        integer :: npts1 ! local copy of grid points in x
        integer :: npts2 ! local copy of grid points in y
        integer :: npts3 ! local copy of grid points in z
        integer :: iph ! phase index (loop)
        integer :: ioerr ! iostat return for I/O error checking
        integer :: ncell ! number of cells in VTK grid
        integer :: np ! number of points (nodes) in VTK grid
        integer :: number_fields ! number of VTK field arrays written
        character(len=5) :: str ! zero-padded step string for filenames

        ictp3d = 12 ! 11-voxel, 12-hexahedron (different ordering)
        ! write (*, '(A)') 'Writing .vtk file'
        npts1 = props%npts1
        npts2 = props%npts2
        npts3 = props%npts3

        ! fluctuation of displacement at cell centers (it is periodic)
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1

            xtmp(1) = float(kx)
            xtmp(2) = float(ky)
            xtmp(3) = float(kz)

            do i = 1, 3
                uf(:, kx, ky, kz) = micro%voxel(kx, ky, kz)%xgr(:) - matmul(micro%defgradavg, xtmp)
            end do

        end do
        end do
        end do

        ! interpolated displacement fluctuation field at cell nodes
        if (npts3 == 1) then

            np = (npts1 + 1) * (npts2 + 1)

            write (str, '(I5.5)') istp
            open (123, file=trim(output_dir)//'//microstr_cell_fields_stp_' &
                  //str//'.vtk', status='replace', action='write')

            write (123, '(A)') '# vtk DataFile Version 3.0'
            write (123, '(A)') 'microstructure data'
            write (123, '(A)') 'ASCII'
            write (123, '(A)') 'DATASET UNSTRUCTURED_GRID'
            write (123, '(A,I10,A)') 'POINTS', np, ' double'

            ! write point coordinates and store point id
            ic = -1
            kz = 1
            do ky = 1, npts2 + 1
            do kx = 1, npts1 + 1
                ic = ic + 1
                id(kx, ky, kz) = ic
                write (123, '(3E18.8E3)') (micro%xnode(i, kx, ky, kz), i=1, 3)
            end do
            end do

            ! cells
            ncell = npts1 * npts2
            write (123, '(A,I10,I10)') 'CELLS', ncell, ncell * 5

            kz = 1
            do ky = 1, npts2
            do kx = 1, npts1

                write (123, '(5I10)') 4, id(kx, ky, kz), id(kx + 1, ky, kz), id(kx, ky + 1, kz), id(kx + 1, ky + 1, kz)

            end do
            end do

            write (123, '(A,I10)') 'CELL_TYPES', ncell
            kz = 1
            do ky = 1, npts2
            do kx = 1, npts1

                write (123, '(I5)') 8 ! 8-pixel, 9-quad (different ordering)

            end do
            end do

        else ! 3D

            np = (npts1 + 1) * (npts2 + 1) * (npts3 + 1)

            write (str, '(I5.5)') istp
            open (123, file=trim(output_dir)//'/microstr_cell_fields_stp_' &
                  //str//'.vtk', status='replace', action='write', iostat=ioerr)

            write (123, '(A)') '# vtk DataFile Version 3.0'
            write (123, '(A)') 'microstructure data'
            write (123, '(A)') 'ASCII'
            write (123, '(A)') 'DATASET UNSTRUCTURED_GRID'
            write (123, '(A,I10,A)') 'POINTS', np, ' double'

            ! write point coordinates and store point id
            ic = -1
            do kz = 1, npts3 + 1
            do ky = 1, npts2 + 1
            do kx = 1, npts1 + 1
                ic = ic + 1
                id(kx, ky, kz) = ic
                write (123, '(3E18.8E3)', iostat=ioerr) (micro%xnode(i, kx, ky, kz), i=1, 3)
                if (ioerr /= 0) stop 'ioerr in writevtk'
            end do
            end do
            end do

            ! cells
            ncell = npts1 * npts2 * npts3
            write (123, '(A,I10,I10)') 'CELLS', ncell, ncell * 9

            do kz = 1, npts3
            do ky = 1, npts2
            do kx = 1, npts1
                if (ictp3d == 11) then
                    write (123, '(9I10)', iostat=ioerr) 8, id(kx, ky, kz), id(kx + 1, ky, kz), &
                        id(kx, ky + 1, kz), id(kx + 1, ky + 1, kz), id(kx, ky, kz + 1), &
                        id(kx + 1, ky, kz + 1), id(kx, ky + 1, kz + 1), id(kx + 1, ky + 1, kz + 1)

                    if (ioerr /= 0) stop 'ioerr in writevtk'
                elseif (ictp3d == 12) then
                    write (123, '(9I10)', iostat=ioerr) 8, id(kx, ky, kz), id(kx + 1, ky, kz), &
                        id(kx + 1, ky + 1, kz), id(kx, ky + 1, kz), id(kx, ky, kz + 1), &
                        id(kx + 1, ky, kz + 1), id(kx + 1, ky + 1, kz + 1), id(kx, ky + 1, kz + 1)
                    if (ioerr /= 0) stop 'ioerr in writevtk'
                end if
            end do
            end do
            end do

            write (123, '(A,I10)') 'CELL_TYPES', ncell
            do kz = 1, npts3
            do ky = 1, npts2
            do kx = 1, npts1
                write (123, '(I5)') ictp3d
            end do
            end do
            end do

        end if

        ! Fields
        number_fields = 11
        write (123, '(A,I10)') 'CELL_DATA', ncell
        write (123, '(A,I10)') 'FIELD FieldData', number_fields

        nsystmx = 0
        do iph = 1, props%nph
            if (props%phase(iph)%nsyst > nsystmx) nsystmx = props%phase(iph)%nsyst
        end do

        write (123, '(A,I3,I10,A)') 'gamdot', nsystmx, ncell, ' double'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1

            write (123, '(100E18.8E3)') (clamp_subnormal(micro%voxel(kx, ky, kz)%gamdot(i)), i=1, nsystmx)

        end do
        end do
        end do

        write (123, '(A,I3,I10,A)') 'gamdot_acum', 1, ncell, ' double'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1

            write (123, '(100E18.8E3)') (clamp_subnormal(micro%voxel(kx, ky, kz)%gacumgr))

        end do
        end do
        end do

        write (123, '(A,I3,I10,A)') 'def_grad_p', 9, ncell, ' double'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1

            write (123, '(9E18.8E3)') clamp_subnormal(micro%voxel(kx, ky, kz)%defgradp(1, 1)), clamp_subnormal(micro%voxel(kx, ky, kz)%defgradp(2, 1)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%defgradp(3, 1)), clamp_subnormal(micro%voxel(kx, ky, kz)%defgradp(1, 2)), clamp_subnormal(micro%voxel(kx, ky, kz)%defgradp(2, 2)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%defgradp(3, 2)), clamp_subnormal(micro%voxel(kx, ky, kz)%defgradp(1, 3)), clamp_subnormal(micro%voxel(kx, ky, kz)%defgradp(2, 3)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%defgradp(3, 3))

        end do
        end do
        end do

        write (123, '(A,I3,I10,A)') 'stress', 6, ncell, ' double'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1
            write (123, '(6E18.8E3)') clamp_subnormal(micro%voxel(kx, ky, kz)%sg(1, 1)), clamp_subnormal(micro%voxel(kx, ky, kz)%sg(2, 2)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%sg(3, 3)), clamp_subnormal(micro%voxel(kx, ky, kz)%sg(1, 2)), clamp_subnormal(micro%voxel(kx, ky, kz)%sg(2, 3)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%sg(1, 3))
        end do
        end do
        end do

        write (123, '(A,I3,I10,A)') 'svm', 1, ncell, ' double'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1
            write (123, '(E18.8E3)') clamp_subnormal(vm_stress(micro%voxel(kx, ky, kz)%sg(:, :)))
        end do
        end do
        end do

        write (123, '(A,I3,I10,A)') 'def_grad', 9, ncell, ' double'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1
            write (123, '(9E18.8E3)') clamp_subnormal(micro%voxel(kx, ky, kz)%defgrad(1, 1)), clamp_subnormal(micro%voxel(kx, ky, kz)%defgrad(2, 1)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%defgrad(3, 1)), clamp_subnormal(micro%voxel(kx, ky, kz)%defgrad(1, 2)), clamp_subnormal(micro%voxel(kx, ky, kz)%defgrad(2, 2)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%defgrad(3, 2)), clamp_subnormal(micro%voxel(kx, ky, kz)%defgrad(1, 3)), clamp_subnormal(micro%voxel(kx, ky, kz)%defgrad(2, 3)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%defgrad(3, 3))
        end do
        end do
        end do

        write (123, '(A,I3,I10,A)') 'det_def_grad', 1, ncell, ' double'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1
            write (123, '(E18.8E3)') clamp_subnormal(micro%voxel(kx, ky, kz)%detF)
        end do
        end do
        end do

        write (123, '(A,I3,I10,A)') 'gr_id', 1, ncell, ' int'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1
            write (123, '(I10)') micro%voxel(kx, ky, kz)%grain
        end do
        end do
        end do

        write (123, '(A,I3,I10,A)') 'ph_id', 1, ncell, ' int'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1
            write (123, '(I10)') micro%voxel(kx, ky, kz)%phase
        end do
        end do
        end do

        write (123, '(A,I3,I10,A)') 'epvm', 1, ncell, ' double'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1
            write (123, '(E18.8E3)') clamp_subnormal(micro%voxel(kx, ky, kz)%epvm)
        end do
        end do
        end do

        write (123, '(A,I3,I10,A)') 'velgrad', 9, ncell, ' double'
        do kz = 1, npts3
        do ky = 1, npts2
        do kx = 1, npts1
            write (123, '(9E18.8E3)') clamp_subnormal(micro%voxel(kx, ky, kz)%velgrad(1, 1)), clamp_subnormal(micro%voxel(kx, ky, kz)%velgrad(2, 1)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%velgrad(3, 1)), clamp_subnormal(micro%voxel(kx, ky, kz)%velgrad(1, 2)), clamp_subnormal(micro%voxel(kx, ky, kz)%velgrad(2, 2)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%velgrad(3, 2)), clamp_subnormal(micro%voxel(kx, ky, kz)%velgrad(1, 3)), clamp_subnormal(micro%voxel(kx, ky, kz)%velgrad(2, 3)), &
                clamp_subnormal(micro%voxel(kx, ky, kz)%velgrad(3, 3))
        end do
        end do
        end do

        close (123)

    end

    subroutine write_results(output_dir, micro, iter, imacroloop, ipr, sim_time, wall_time, tdot, tdot_ref, gamma, g, geff, init)
        use types, only: micro_type
        implicit none
        type(micro_type), intent(in) :: micro
        character(len=*), intent(in) :: output_dir
        integer, intent(in) :: iter, imacroloop, ipr
        double precision, intent(in) :: sim_time, wall_time, tdot, tdot_ref, gamma, g, geff
        logical, intent(in), optional :: init

        character(len=512) :: csv_path
        character(len=2048) :: header
        integer :: u, ios
        logical :: do_init, exists
        integer :: nout

        nout = len_trim(output_dir)
        if (nout == 0) stop 'write_results: output_dir is empty'

        if (output_dir(nout:nout) == '/') then
            csv_path = trim(output_dir)//'results.csv'
        else
            csv_path = trim(output_dir)//'/results.csv'
        end if

        do_init = .false.
        if (present(init)) do_init = init

        header = 'iter,macrostep,ipr,sim_time,wall_time,tdot,tdot_ref,'// &
                 'eav11,eav22,eav33,eav23,eav13,eav12,'// &
                 'sav11,sav22,sav33,sav23,sav13,sav12,'// &
                 'epav11,epav22,epav33,epav23,epav13,epav12,'// &
                 'eelav11,eelav22,eelav33,eelav23,eelav13,eelav12,'// &
                 'vgrad11,vgrad22,vgrad33,vgrad23,vgrad13,vgrad12,'// &
                 'edotp11,edotp22,edotp33,edotp23,edotp13,edotp12,'// &
                 'evm,evmp,dvm,dvmp,svm,gamma,g,geff'

        if (do_init) then
            open (newunit=u, file=trim(csv_path), status='replace', action='write', iostat=ios)
            if (ios /= 0) stop 'results: failed to open results.csv for writing'
            write (u, '(A)') trim(header)
            close (u)
        end if

        inquire (file=trim(csv_path), exist=exists)
        if (.not. exists) then
            open (newunit=u, file=trim(csv_path), status='replace', action='write', iostat=ios)
            if (ios /= 0) stop 'write_results: failed to create results.csv'
            write (u, '(A)') trim(header)
            close (u)
        end if

        open (newunit=u, file=trim(csv_path), status='old', position='append', action='write', iostat=ios)
        if (ios /= 0) stop 'write_results: failed to open results.csv for append'

        write (u, '(*(g0,:,","))') iter, imacroloop, ipr, sim_time, wall_time, tdot, tdot_ref, &
            micro%eavg(1, 1), micro%eavg(2, 2), micro%eavg(3, 3), micro%eavg(2, 3), micro%eavg(1, 3), micro%eavg(1, 2), &
            micro%sgsolavg(1, 1), micro%sgsolavg(2, 2), micro%sgsolavg(3, 3), micro%sgsolavg(2, 3), micro%sgsolavg(1, 3), micro%sgsolavg(1, 2), &
            micro%epavg(1, 1), micro%epavg(2, 2), micro%epavg(3, 3), micro%epavg(2, 3), micro%epavg(1, 3), micro%epavg(1, 2), &
            micro%eelavg(1, 1), micro%eelavg(2, 2), micro%eelavg(3, 3), micro%eelavg(2, 3), micro%eelavg(1, 3), micro%eelavg(1, 2), &
            micro%velgradavg(1, 1), micro%velgradavg(2, 2), micro%velgradavg(3, 3), micro%velgradavg(2, 3), micro%velgradavg(1, 3), micro%velgradavg(1, 2), &
            micro%edotpavg(1, 1), micro%edotpavg(2, 2), micro%edotpavg(3, 3), micro%edotpavg(2, 3), micro%edotpavg(1, 3), micro%edotpavg(1, 2), &
            micro%evmavg, micro%evmpavg, micro%dvmavg, micro%dvmpavg, micro%svmavg, gamma, g, geff

        flush (u)
        close (u)

    end subroutine write_results

end module IO_functions
