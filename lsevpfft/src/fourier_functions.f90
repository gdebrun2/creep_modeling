module fourier_functions
!$  use omp_lib
    use FFTW3
    implicit none

    character(len=256) :: fftw_wisdom_file = 'fftw_wisdom.save'
!$  character(len=256) :: fftw_wisdom_file_omp = 'fftw_wisdom_omp.save'
    real(C_DOUBLE), parameter :: max_planning_time = 10.0

    ! FFTW PLANNING MODE: (http://www.fftw.org/fftw3_doc/Planner-Flags.html):
    integer(C_INT) :: FFTW_PLANNING_FLAG = FFTW_PATIENT
    integer(C_INT) :: FFTW_PLAN_COMBINED

    type(C_PTR) :: plan_advanced33, iplan_advanced33, plan_advanced3, iplan_advanced3
    type(C_PTR) :: plan_advanced333, iplan_advanced333, plan_advanced3333, iplan_advanced3333
    type(C_PTR) :: plan_advanced, iplan_advanced
    integer :: gridpoints

    integer :: fft_verbosity = 2
    complex(C_DOUBLE_COMPLEX), pointer :: fourgrid33(:, :, :, :, :)
    complex(C_DOUBLE_COMPLEX), pointer :: fourgrid3(:, :, :, :)
    complex(C_DOUBLE_COMPLEX), pointer :: fourgrid333(:, :, :, :, :, :)
    complex(C_DOUBLE_COMPLEX), pointer :: fourgrid3333(:, :, :, :, :, :, :)
    complex(C_DOUBLE_COMPLEX), pointer :: fourgrid(:, :, :)

    private :: gridpoints, fftw_wisdom_file, max_planning_time, FFTW_PLANNING_FLAG, FFTW_PLAN_COMBINED
!$  private :: fftw_wisdom_file_omp
    public :: plan_advanced33, iplan_advanced33, plan_advanced3, iplan_advanced3, fourgrid33, fourgrid3, fft_verbosity
    public :: plan_advanced333, iplan_advanced333, fourgrid333, plan_advanced3333, iplan_advanced3333, fourgrid3333
    public :: plan_advanced, iplan_advanced, fourgrid

! **************************************************
contains

    subroutine init_fftw_wisdom_path(path)
!$      use omp_lib
        character(len=*), intent(in) :: path
        logical :: openmp_enabled = .false.
        character(len=256) :: parent
        integer :: lastslash
        logical :: dir_exists

!$      openmp_enabled = .true.

        if (len_trim(path) > 0) then

            lastslash = index(trim(path), '/', back=.true.)
            if (lastslash == 0) then
                parent = '.'                      ! no slash => current directory
            else if (lastslash == 1) then
                parent = '/'                      ! parent is root
            else
                parent = path(:lastslash - 1)   ! substring up to before last '/'
            end if

            inquire (file=trim(parent), exist=dir_exists)
            if (.not. dir_exists) call system('mkdir -p "'//trim(parent)//'"')

            if (openmp_enabled) then
!$              fftw_wisdom_file_omp = trim(path)
            else
                fftw_wisdom_file = trim(path)
            end if
        end if

    end subroutine init_fftw_wisdom_path

    ! Setup for FFTW
    subroutine initialize_fftw(npts1, npts2, npts3, nthreads)
!$      use omp_lib
        use FFTW3
        implicit none

        integer, intent(in) :: npts1, npts2, npts3, nthreads
        logical :: openmp_enabled = .false.
        logical :: dir_e = .false.
        logical :: wisdom_read = .false.
        integer(C_INT) :: non_zero_success
!$      integer :: iret
        type(C_PTR) :: fftw_c_ptr33, fftw_c_ptr3, fftw_c_ptr333, fftw_c_ptr3333, fftw_c_ptr
        integer :: i
        integer(C_INT), dimension(3) :: gridshape

        write (*, '(A,I0,A,I0,A,I0,A)') 'Initializing FFTW on grid (x,y,z)=(', npts1, ',', npts2, ',', npts3, ')'
        gridshape = [npts1, npts2, npts3]
        gridpoints = npts3 * npts2 * npts1

!$      openmp_enabled = .true.
!$      iret = fftw_init_threads()
!$      if (iret == 0) stop '-> Error initializing FFTW threads'
!$      call fftw_plan_with_nthreads(nthreads)
        call fftw_set_timelimit(max_planning_time)

        ! Allocate fourier grids as a fully contigious array with row-major ordering (note the reverse indexing)
        fftw_c_ptr33 = fftw_alloc_complex(int(9 * gridpoints, C_SIZE_T))
        fftw_c_ptr3 = fftw_alloc_complex(int(3 * gridpoints, C_SIZE_T))
        fftw_c_ptr333 = fftw_alloc_complex(int(27 * gridpoints, C_SIZE_T))
        fftw_c_ptr3333 = fftw_alloc_complex(int(81 * gridpoints, C_SIZE_T))
        fftw_c_ptr = fftw_alloc_complex(int(1 * gridpoints, C_SIZE_T))
        call c_f_pointer(fftw_c_ptr33, fourgrid33, [npts3, npts2, npts1, 3, 3])
        call c_f_pointer(fftw_c_ptr3, fourgrid3, [npts3, npts2, npts1, 3])
        call c_f_pointer(fftw_c_ptr333, fourgrid333, [npts3, npts2, npts1, 3, 3, 3])
        call c_f_pointer(fftw_c_ptr3333, fourgrid3333, [npts3, npts2, npts1, 3, 3, 3, 3])
        call c_f_pointer(fftw_c_ptr, fourgrid, [npts3, npts2, npts1])

        ! Read previous fftw wisdom on best algorithm for user hardware (if exists)
        if (openmp_enabled) then
!$          inquire (file=trim(fftw_wisdom_file_omp), exist=dir_e)
!$          if (dir_e) then
!$              non_zero_success = fftw_import_wisdom_from_filename(trim(fftw_wisdom_file_omp)//c_null_char)
!$              if (non_zero_success == 0) then
!$                  write (*, '(A)') '-> ERROR IMPORTING FFTW WISDOM'
!$              else
!$                  wisdom_read = .true.
!$              end if
!$          end if
        else
            inquire (file=trim(fftw_wisdom_file), exist=dir_e)
            if (dir_e) then
                non_zero_success = fftw_import_wisdom_from_filename(trim(fftw_wisdom_file)//c_null_char)
                if (non_zero_success == 0) then
                    write (*, '(A)') '-> ERROR IMPORTING FFTW WISDOM'
                else
                    wisdom_read = .true.
                end if
            end if
        end if

        ! Loop 1 sets up plans assuming wisdom for the problem exists
        ! Loop 2 sets up plans only if loop 1 fails to set them up
        do i = 1, 2
            ! If wisdom found, do not create additional wisdom for problem
            FFTW_PLAN_COMBINED = FFTW_PLANNING_FLAG
            if (wisdom_read .and. i == 1) FFTW_PLAN_COMBINED = FFTW_PLAN_COMBINED + FFTW_WISDOM_ONLY
            if ((wisdom_read .and. i == 2) .or. (.not. wisdom_read .and. i == 1)) then
                write (*, '(A,F4.1,A)') '-> Planning FFTW algorithm. Max planning time: ', max_planning_time, 's'
            end if

            ! Build advanced plans for rank 1 & 2 tensor transforms
            ! This advanced plan is more efficient as it performs all transforms perfectly sequentially
            plan_advanced33 = fftw_plan_many_dft(3, gridshape, 9, &
                                                 fourgrid33, gridshape, 1, gridpoints, &
                                                 fourgrid33, gridshape, 1, gridpoints, FFTW_FORWARD, FFTW_PLAN_COMBINED)
            iplan_advanced33 = fftw_plan_many_dft(3, gridshape, 9, &
                                                  fourgrid33, gridshape, 1, gridpoints, &
                                                  fourgrid33, gridshape, 1, gridpoints, FFTW_BACKWARD, FFTW_PLAN_COMBINED)
            plan_advanced3 = fftw_plan_many_dft(3, gridshape, 3, &
                                                fourgrid3, gridshape, 1, gridpoints, &
                                                fourgrid3, gridshape, 1, gridpoints, FFTW_FORWARD, FFTW_PLAN_COMBINED)
            iplan_advanced3 = fftw_plan_many_dft(3, gridshape, 3, &
                                                 fourgrid3, gridshape, 1, gridpoints, &
                                                 fourgrid3, gridshape, 1, gridpoints, FFTW_BACKWARD, FFTW_PLAN_COMBINED)
            plan_advanced333 = fftw_plan_many_dft(3, gridshape, 27, &
                                                  fourgrid333, gridshape, 1, gridpoints, &
                                                  fourgrid333, gridshape, 1, gridpoints, FFTW_FORWARD, FFTW_PLAN_COMBINED)
            iplan_advanced333 = fftw_plan_many_dft(3, gridshape, 27, &
                                                   fourgrid333, gridshape, 1, gridpoints, &
                                                   fourgrid333, gridshape, 1, gridpoints, FFTW_BACKWARD, FFTW_PLAN_COMBINED)
            plan_advanced3333 = fftw_plan_many_dft(3, gridshape, 81, &
                                                   fourgrid3333, gridshape, 1, gridpoints, &
                                                   fourgrid3333, gridshape, 1, gridpoints, FFTW_FORWARD, FFTW_PLAN_COMBINED)
            iplan_advanced3333 = fftw_plan_many_dft(3, gridshape, 81, &
                                                    fourgrid3333, gridshape, 1, gridpoints, &
                                                    fourgrid3333, gridshape, 1, gridpoints, FFTW_BACKWARD, FFTW_PLAN_COMBINED)
            plan_advanced = fftw_plan_many_dft(3, gridshape, 1, &
                                               fourgrid, gridshape, 1, gridpoints, &
                                               fourgrid, gridshape, 1, gridpoints, FFTW_FORWARD, FFTW_PLAN_COMBINED)
            iplan_advanced = fftw_plan_many_dft(3, gridshape, 1, &
                                                fourgrid, gridshape, 1, gridpoints, &
                                                fourgrid, gridshape, 1, gridpoints, FFTW_BACKWARD, FFTW_PLAN_COMBINED)

            if (c_associated(plan_advanced33) .and. c_associated(iplan_advanced33) .and. &
                c_associated(plan_advanced3) .and. c_associated(iplan_advanced3) .and. &
                c_associated(plan_advanced333) .and. c_associated(iplan_advanced333) .and. &
                c_associated(plan_advanced3333) .and. c_associated(iplan_advanced3333) .and. &
                c_associated(plan_advanced) .and. c_associated(iplan_advanced)) then
                exit
            else if (i == 2) then
                stop '-> FFTW planning error, try deleting fftw_wisdom file(s)'
            end if
        end do

        ! Write wisdom

        if (openmp_enabled) then
!$          write (*, '(A)') 'FFTW wisdom file: '//trim(fftw_wisdom_file_omp)
!$          non_zero_success = fftw_export_wisdom_to_filename(trim(fftw_wisdom_file_omp)//c_null_char)
!$          if (non_zero_success == 0) write (*, '(A)') '-> ERROR EXPORTING FFTW WISDOM'
!$          if (fft_verbosity > 1) write (*, '(A, I0, A)') '-> Successfully setup parallel FFTW plans with ', nthreads, ' threads'
        else
            non_zero_success = fftw_export_wisdom_to_filename(trim(fftw_wisdom_file)//c_null_char)
            if (non_zero_success == 0) write (*, '(A)') '-> ERROR EXPORTING FFTW WISDOM'
            if (fft_verbosity > 1) write (*, '(A)') '-> Successfully setup serial FFTW plans'
        end if

    end subroutine initialize_fftw

    ! **************************************************
    subroutine fft_tensor33(plan, grid33)
        use FFTW3
        implicit none

        type(C_PTR), intent(in) :: plan
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: grid33(:, :, :, :, :)

        ! Call to FFTW
        call fftw_execute_dft(plan, grid33, grid33)

    end

    ! **************************************************
    subroutine ifft_tensor33(iplan, grid33)
        use FFTW3
        implicit none

        type(C_PTR), intent(in) :: iplan
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: grid33(:, :, :, :, :)

        ! Call to FFTW
        call fftw_execute_dft(iplan, grid33, grid33)

        grid33 = real(grid33, 8) / float(gridpoints)

    end

    ! **************************************************
    subroutine fft_tensor3(plan, grid3)
        use FFTW3
        implicit none

        type(C_PTR), intent(in) :: plan
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: grid3(:, :, :, :)

        ! Call to FFTW
        call fftw_execute_dft(plan, grid3, grid3)

    end

    ! **************************************************
    subroutine ifft_tensor3(iplan, grid3)
        use FFTW3
        implicit none

        type(C_PTR), intent(in) :: iplan
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: grid3(:, :, :, :)

        ! Call to FFTW
        call fftw_execute_dft(iplan, grid3, grid3)

        grid3 = real(grid3, 8) / float(gridpoints)

    end

    ! **************************************************
    subroutine fft_scalar(plan, grid)
        use FFTW3
        implicit none

        type(C_PTR), intent(in) :: plan
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: grid(:, :, :)

        ! Call to FFTW
        call fftw_execute_dft(plan, grid, grid)

    end

    ! **************************************************
    subroutine ifft_scalar(iplan, grid)
        use FFTW3
        implicit none

        type(C_PTR), intent(in) :: iplan
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: grid(:, :, :)

        ! Call to FFTW
        call fftw_execute_dft(iplan, grid, grid)

        grid = real(grid, 8) / float(gridpoints)

    end

    ! **************************************************
    subroutine fft_tensor333(plan, grid333)
        use FFTW3
        implicit none

        type(C_PTR), intent(in) :: plan
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: grid333(:, :, :, :, :, :)

        ! Call to FFTW
        call fftw_execute_dft(plan, grid333, grid333)

    end

    ! **************************************************
    subroutine ifft_tensor333(iplan, grid333)
        use FFTW3
        implicit none

        type(C_PTR), intent(in) :: iplan
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: grid333(:, :, :, :, :, :)

        ! Call to FFTW
        call fftw_execute_dft(iplan, grid333, grid333)

        grid333 = real(grid333, 8) / float(gridpoints)

    end

    ! **************************************************
    subroutine fft_tensor3333(plan, grid3333)
        use FFTW3
        implicit none

        type(C_PTR), intent(in) :: plan
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: grid3333(:, :, :, :, :, :, :)

        ! Call to FFTW
        call fftw_execute_dft(plan, grid3333, grid3333)

    end

    ! **************************************************
    subroutine ifft_tensor3333(iplan, grid3333)
        use FFTW3
        implicit none

        type(C_PTR), intent(in) :: iplan
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: grid3333(:, :, :, :, :, :, :)

        ! Call to FFTW
        call fftw_execute_dft(iplan, grid3333, grid3333)

        grid3333 = real(grid3333, 8) / float(gridpoints)

    end

    ! **************************************************
    ! Computes the spatial frequency xi1,xi2,xi3 for a given grid point x,y,z
    pure function spatial_freq(x, y, z, npts1, npts2, npts3) result(xi)
        double precision :: xi(3)
        integer, intent(in) :: x, y, z, npts1, npts2, npts3
        double precision, parameter :: twopi = 8.d0 * datan(1.d0)

        if (x <= npts1 / 2) then
            xi(1) = twopi * float(x - 1) / float(npts1)
        else if (x > npts1 / 2) then
            xi(1) = twopi * float(x - npts1 - 1) / float(npts1)
        else
            xi(1) = 0.
        end if

        if (y <= npts2 / 2) then
            xi(2) = twopi * float(y - 1) / float(npts2)
        else if (y > npts2 / 2) then
            xi(2) = twopi * float(y - npts2 - 1) / float(npts2)
        else
            xi(2) = 0.
        end if

        if (npts3 == 1) then
            xi(3) = 0.
        else if (z <= npts3 / 2) then
            xi(3) = twopi * float(z - 1) / float(npts3)
        else if (z > npts3 / 2) then
            xi(3) = twopi * float(z - npts3 - 1) / float(npts3)
        else
            xi(3) = 0.
        end if

    end function spatial_freq

    ! DEPRECATED SUBROUTINES
    ! USE THE SUBROUTINES ABOVE FOR FASTER FFT EXECUTION

    ! **************************************************
    subroutine fft_tensor2d(f, fhatr, fhati, idim1, idim2, npts1, npts2, npts3)
!$      use omp_lib
        use FFTW3
        implicit none
        integer :: idim1, idim2, ic1, ic2, npts1, npts2, npts3, kxx, kyy, kzz
        double precision :: f(idim1, idim2, npts1, npts2, npts3)
        double precision :: fhatr(idim1, idim2, npts1, npts2, npts3)
        double precision :: fhati(idim1, idim2, npts1, npts2, npts3)

        if (idim1 /= 3 .or. idim2 /= 3) stop 'fft_tensor2d only supports 3x3 tensor input due to FFTW plan setup'

        ! Column to row major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
        do ic2 = 1, idim2
            fourgrid33(kzz, kyy, kxx, ic2, ic1) = cmplx(f(ic1, ic2, kxx, kyy, kzz), 0.d0, 8)
        end do
        end do
        end do
        end do
        end do

        ! Call to FFTW
        call fftw_execute_dft(plan_advanced33, fourgrid33, fourgrid33)

        ! Row to column major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
        do ic2 = 1, idim2
            fhatr(ic1, ic2, kxx, kyy, kzz) = real(fourgrid33(kzz, kyy, kxx, ic2, ic1), 8)
            fhati(ic1, ic2, kxx, kyy, kzz) = -aimag(fourgrid33(kzz, kyy, kxx, ic2, ic1))
        end do
        end do
        end do
        end do
        end do

    end

    ! **************************************************
    subroutine ifft_tensor2d(fhatr, fhati, f, idim1, idim2, npts1, npts2, npts3)
!$      use omp_lib
        use FFTW3
        implicit none
        integer :: idim1, idim2, ic1, ic2, npts1, npts2, npts3, kxx, kyy, kzz
        double precision :: f(idim1, idim2, npts1, npts2, npts3)
        double precision :: fhatr(idim1, idim2, npts1, npts2, npts3)
        double precision :: fhati(idim1, idim2, npts1, npts2, npts3)

        if (idim1 /= 3 .or. idim2 /= 3) stop 'ifft_tensor2d only supports 3x3 tensor input due to FFTW plan setup'

        ! Column to row major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
        do ic2 = 1, idim2
            fourgrid33(kzz, kyy, kxx, ic2, ic1) = cmplx(fhatr(ic1, ic2, kxx, kyy, kzz), -fhati(ic1, ic2, kxx, kyy, kzz), 8)
        end do
        end do
        end do
        end do
        end do

        ! Call to FFTW
        call fftw_execute_dft(iplan_advanced33, fourgrid33, fourgrid33)

        ! Row to column major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
        do ic2 = 1, idim2
            f(ic1, ic2, kxx, kyy, kzz) = real(fourgrid33(kzz, kyy, kxx, ic2, ic1), 8) / gridpoints
        end do
        end do
        end do
        end do
        end do

    end

    ! **************************************************
    subroutine fft_tensor1d(f, fhatr, fhati, idim1, npts1, npts2, npts3)
!$      use omp_lib
        use FFTW3
        implicit none
        integer :: idim1, ic1, npts1, npts2, npts3, kxx, kyy, kzz
        double precision :: f(idim1, npts1, npts2, npts3)
        double precision :: fhatr(idim1, npts1, npts2, npts3)
        double precision :: fhati(idim1, npts1, npts2, npts3)

        if (idim1 /= 3) stop 'fft_tensor1d only supports 1x3 tensor input due to FFTW plan setup'

        ! Column to row major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
            fourgrid3(kzz, kyy, kxx, ic1) = cmplx(f(ic1, kxx, kyy, kzz), 0.d0, 8)
        end do
        end do
        end do
        end do

        ! Call to FFTW
        call fftw_execute_dft(plan_advanced3, fourgrid3, fourgrid3)

        ! Row to column major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
            fhatr(ic1, kxx, kyy, kzz) = real(fourgrid3(kzz, kyy, kxx, ic1), 8)
            fhati(ic1, kxx, kyy, kzz) = -aimag(fourgrid3(kzz, kyy, kxx, ic1))
        end do
        end do
        end do
        end do
    end

    ! **************************************************
    subroutine ifft_tensor1d(fhatr, fhati, f, idim1, npts1, npts2, npts3)
!$      use omp_lib
        use FFTW3
        implicit none
        integer :: idim1, ic1, npts1, npts2, npts3, kxx, kyy, kzz
        double precision :: f(idim1, npts1, npts2, npts3)
        double precision :: fhatr(idim1, npts1, npts2, npts3)
        double precision :: fhati(idim1, npts1, npts2, npts3)

        if (idim1 /= 3) stop 'ifft_tensor1d only supports 1x3 tensor input due to FFTW plan setup'

        ! Column to row major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
            fourgrid3(kzz, kyy, kxx, ic1) = cmplx(fhatr(ic1, kxx, kyy, kzz), -fhati(ic1, kxx, kyy, kzz), 8)
        end do
        end do
        end do
        end do

        ! Call to FFTW
        call fftw_execute_dft(iplan_advanced3, fourgrid3, fourgrid3)

        ! Row to column major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
            f(ic1, kxx, kyy, kzz) = real(fourgrid3(kzz, kyy, kxx, ic1), 8) / gridpoints
        end do
        end do
        end do
        end do
    end

    ! **************************************************
    subroutine fft_tensor3d(f, fhatr, fhati, idim1, idim2, idim3, npts1, npts2, npts3)
!$      use omp_lib
        use FFTW3
        implicit none
        integer :: idim1, idim2, idim3, ic1, ic2, ic3, npts1, npts2, npts3, kxx, kyy, kzz
        double precision :: f(idim1, idim2, idim3, npts1, npts2, npts3)
        double precision :: fhatr(idim1, idim2, idim3, npts1, npts2, npts3)
        double precision :: fhati(idim1, idim2, idim3, npts1, npts2, npts3)

        if (idim1 /= 3 .or. idim2 /= 3 .or. idim3 /= 3) stop 'fft_tensor3d only supports 3x3x3 tensor input due to FFTW plan setup'

        ! Column to row major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
        do ic2 = 1, idim2
        do ic3 = 1, idim3
            fourgrid333(kzz, kyy, kxx, ic3, ic2, ic1) = cmplx(f(ic1, ic2, ic3, kxx, kyy, kzz), 0.d0, 8)
        end do
        end do
        end do
        end do
        end do
        end do

        ! Call to FFTW
        call fftw_execute_dft(plan_advanced333, fourgrid333, fourgrid333)

        ! Row to column major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
        do ic2 = 1, idim2
        do ic3 = 1, idim3
            fhatr(ic1, ic2, ic3, kxx, kyy, kzz) = real(fourgrid333(kzz, kyy, kxx, ic3, ic2, ic1), 8)
            fhati(ic1, ic2, ic3, kxx, kyy, kzz) = -aimag(fourgrid333(kzz, kyy, kxx, ic3, ic2, ic1))
        end do
        end do
        end do
        end do
        end do
        end do

    end

    ! **************************************************
    subroutine ifft_tensor3d(fhatr, fhati, f, idim1, idim2, idim3, npts1, npts2, npts3)
!$      use omp_lib
        use FFTW3
        implicit none
        integer :: idim1, idim2, idim3, ic1, ic2, ic3, npts1, npts2, npts3, kxx, kyy, kzz
        double precision :: f(idim1, idim2, idim3, npts1, npts2, npts3)
        double precision :: fhatr(idim1, idim2, idim3, npts1, npts2, npts3)
        double precision :: fhati(idim1, idim2, idim3, npts1, npts2, npts3)

        if (idim1 /= 3 .or. idim2 /= 3 .or. idim3 /= 3) &
            stop 'ifft_tensor3d only supports 3x3x3 tensor input due to FFTW plan setup'

        ! Column to row major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
        do ic2 = 1, idim2
        do ic3 = 1, idim3
            fourgrid333(kzz, kyy, kxx, ic3, ic2, ic1) = cmplx(fhatr(ic1, ic2, ic3, kxx, kyy, kzz), -fhati(ic1, ic2, ic3, kxx, kyy, kzz), 8)
        end do
        end do
        end do
        end do
        end do
        end do

        ! Call to FFTW
        call fftw_execute_dft(iplan_advanced333, fourgrid333, fourgrid333)

        ! Row to column major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
        do ic2 = 1, idim2
        do ic3 = 1, idim3
            f(ic1, ic2, ic3, kxx, kyy, kzz) = real(fourgrid333(kzz, kyy, kxx, ic3, ic2, ic1), 8) / gridpoints
        end do
        end do
        end do
        end do
        end do
        end do

    end

    ! **************************************************
    subroutine ifft_tensor4d(fhatr, fhati, f, idim1, idim2, idim3, idim4, npts1, npts2, npts3)
!$      use omp_lib
        use FFTW3
        implicit none
        integer :: idim1, idim2, idim3, idim4, ic1, ic2, ic3, ic4, npts1, npts2, npts3, kxx, kyy, kzz
        double precision :: f(idim1, idim2, idim3, idim4, npts1, npts2, npts3)
        double precision :: fhatr(idim1, idim2, idim3, idim4, npts1, npts2, npts3)
        double precision :: fhati(idim1, idim2, idim3, idim4, npts1, npts2, npts3)

        if (idim1 /= 3 .or. idim2 /= 3 .or. idim3 /= 3 .or. idim4 /= &
            3) stop 'ifft_tensor4d only supports 3x3x3x3 tensor input due to FFTW plan setup'

        ! Column to row major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
        do ic2 = 1, idim2
        do ic3 = 1, idim3
        do ic4 = 1, idim4
            fourgrid3333(kzz, kyy, kxx, ic4, ic3, ic2, ic1) = cmplx(fhatr(ic1, ic2, ic3, ic4, kxx, kyy, kzz), -fhati(ic1, ic2, ic3, ic4, kxx, kyy, kzz), 8)
        end do
        end do
        end do
        end do
        end do
        end do
        end do

        ! Call to FFTW
        call fftw_execute_dft(iplan_advanced3333, fourgrid3333, fourgrid3333)

        ! Row to column major order
        do kxx = 1, npts1
        do kyy = 1, npts2
        do kzz = 1, npts3
        do ic1 = 1, idim1
        do ic2 = 1, idim2
        do ic3 = 1, idim3
        do ic4 = 1, idim4
            f(ic1, ic2, ic3, ic4, kxx, kyy, kzz) = real(fourgrid3333(kzz, kyy, kxx, ic4, ic3, ic2, ic1), 8) / gridpoints
        end do
        end do
        end do
        end do
        end do
        end do
        end do

    end

end module fourier_functions
