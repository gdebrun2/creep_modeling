module various_functions

    implicit none

contains

    ! Initialize voxel stress/velgrad fields for the first process (elastic predictor).
    subroutine initial_guess(micro, props)
!$      use omp_lib
        use types
        use tensor_functions
        ! use IO_functions
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (updated in-place)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: udot6(6) ! imposed udot in 6-component basis (Voigt-like)
        double precision :: aux66(6, 6) ! scratch 6x6 matrix (basis transforms)
        double precision :: aux3333(3, 3, 3, 3) ! scratch 4th-order tensor (basis transforms)
        double precision :: sg6(6) ! voxel stress in 6-component basis
        double precision :: sg33(3, 3) ! voxel stress as 3x3 tensor
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, sg6, aux66,aux3333, udot6, sg33) &
        !$OMP SHARED(props,micro)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            micro%voxel(ip1, ip2, ip3)%velgrad(:, :) = props%process(1)%udot

            if (.not. micro%voxel(ip1, ip2, ip3)%gas) then

                call chg_basis(udot6, props%process(1)%udot, aux66, aux3333, 2, 6)
                sg6 = matmul(micro%voxel(ip1, ip2, ip3)%cg66(:, :), udot6 * props%process(1)%tdot)
                call chg_basis(sg6, sg33, aux66, aux3333, 1, 6)
                micro%voxel(ip1, ip2, ip3)%sg(:, :) = sg33

            end if

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        micro%dvmavg = vm_strain(props%process(1)%dsim)
        micro%lvmavg = micro%dvmavg

    end

    ! Compute the spatial gradient of a 3-component vector field using FFT differentiation.
    subroutine gradient_tensor1(tensor1, grad_tensor1, ideriv, npts1, npts2, npts3)
        use fourier_functions
        implicit none

        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        integer, intent(in) :: ideriv ! derivative operator selector (spectral)
        double precision, intent(in) :: tensor1(3, npts1, npts2, npts3) ! input vector field (cell-centered)
        double precision, intent(out) :: grad_tensor1(3, 3, npts1, npts2, npts3) ! output gradient (3x3) per cell
        double precision :: xk(3) ! spatial frequency vector at (i,j,k)
        complex(kind(1.0d0)) :: kmod(3) ! complex wavevector modifier for derivative
        integer :: i, j, k, l, m ! loop indices

        call tensor1_to_fouriergrid(tensor1, npts1, npts2, npts3)
        call fft_tensor3(plan_advanced3, fourgrid3)

        do k = 1, npts3
        do j = 1, npts2
        do i = 1, npts1

            xk = spatial_freq(i, j, k, npts1, npts2, npts3)
            if (ideriv == 2) kmod = cmplx(0.0, 1.0) * 0.25 * &
                                    (1.0 + exp(cmplx(0.0, xk(1), kind=kind(0.0d0)))) * (1.0 + exp(cmplx(0.0, xk(2), kind=kind(0.0d0)))) * &
                                    (1.0 + exp(cmplx(0.0, xk(3), kind=kind(0.0d0))))

            do l = 1, 3
                if (ideriv == 0) then
                    kmod(l) = cmplx(0.0, xk(l), kind=kind(0.0d0))
                elseif (ideriv == 1) then
                    kmod(l) = cmplx(0.0, sin(xk(l)), kind=kind(0.0d0))
                elseif (ideriv == 2) then
                    kmod(l) = tan(xk(l) / 2.0) * kmod(l)
                end if
            end do

            do l = 1, 3
            do m = 1, 3
                fourgrid33(k, j, i, m, l) = fourgrid3(k, j, i, l) * kmod(m)
            end do
            end do

            if (((mod(npts1, 2) == 0 .and. i == npts1 / 2 + 1) .and. &
                 (mod(npts2, 2) == 0 .and. j == npts2 / 2 + 1) .and. &
                 ((mod(npts3, 2) == 0 .and. k == npts3 / 2 + 1) .and. npts3 > 1)) .and. ideriv == 2) then

                fourgrid33(k, j, i, :, :) = (0.0, 0.0)

            end if

        end do
        end do
        end do

        call ifft_tensor33(iplan_advanced33, fourgrid33)

        call fouriergrid_to_tensor2(grad_tensor1, npts1, npts2, npts3)

    end

    ! Pack a real-space vector field into the FFTW complex work array (fourgrid3).
    subroutine tensor1_to_fouriergrid(tensor1, npts1, npts2, npts3)
        use fourier_functions
        implicit none
        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        double precision, intent(in) :: tensor1(3, npts1, npts2, npts3) ! real-space vector field
        integer :: i, j, k, l ! loop indices

        do k = 1, npts3
        do j = 1, npts2
        do i = 1, npts1
        do l = 1, 3
            fourgrid3(k, j, i, l) = cmplx(tensor1(l, i, j, k), kind=kind(0.0d0))
        end do
        end do
        end do
        end do

    end

    ! Unpack the FFTW complex work array (fourgrid3) into a real-space vector field.
    subroutine fouriergrid_to_tensor1(tensor1, npts1, npts2, npts3)
        use fourier_functions
        implicit none
        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        double precision, intent(out) :: tensor1(3, npts1, npts2, npts3) ! real-space vector field (output)
        integer :: i, j, k, l ! loop indices

        do k = 1, npts3
        do j = 1, npts2
        do i = 1, npts1
        do l = 1, 3
            tensor1(l, i, j, k) = real(fourgrid3(k, j, i, l))
        end do
        end do
        end do
        end do

    end

    ! Pack a real-space 3x3 tensor field into the FFTW complex work array (fourgrid33).
    subroutine tensor2_to_fouriergrid(tensor2, npts1, npts2, npts3)
        use fourier_functions
        implicit none
        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        double precision, intent(in) :: tensor2(3, 3, npts1, npts2, npts3) ! real-space 2nd-order tensor field
        integer :: i, j, k, l, m ! loop indices

        do k = 1, npts3
        do j = 1, npts2
        do i = 1, npts1
        do l = 1, 3
        do m = 1, 3
            fourgrid33(k, j, i, m, l) = cmplx(tensor2(l, m, i, j, k), kind=kind(0.0d0))
        end do
        end do
        end do
        end do
        end do

    end

    ! Unpack the FFTW complex work array (fourgrid33) into a real-space 3x3 tensor field.
    subroutine fouriergrid_to_tensor2(tensor2, npts1, npts2, npts3)
        use fourier_functions
        implicit none
        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        double precision, intent(out) :: tensor2(3, 3, npts1, npts2, npts3) ! real-space 2nd-order tensor field (output)
        integer :: i, j, k, l, m ! loop indices

        do k = 1, npts3
        do j = 1, npts2
        do i = 1, npts1
        do l = 1, 3
        do m = 1, 3
            tensor2(l, m, i, j, k) = real(fourgrid33(k, j, i, m, l))
        end do
        end do
        end do
        end do
        end do

    end

    ! Pack a real-space 4th-order tensor field into the FFTW complex work array (fourgrid3333).
    subroutine tensor4_to_fouriergrid(tensor4, npts1, npts2, npts3)
        use fourier_functions
        implicit none
        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        double precision, intent(in) :: tensor4(3, 3, 3, 3, npts1, npts2, npts3) ! real-space 4th-order tensor field
        integer :: i, j, k, l, m, n, o ! loop indices

        do k = 1, npts3
        do j = 1, npts2
        do i = 1, npts1
        do l = 1, 3
        do m = 1, 3
        do n = 1, 3
        do o = 1, 3
            fourgrid3333(k, j, i, o, n, m, l) = cmplx(tensor4(l, m, n, o, i, j, k), kind=kind(0.0d0))
        end do
        end do
        end do
        end do
        end do
        end do
        end do

    end

    ! Unpack the FFTW complex work array (fourgrid3333) into a real-space 4th-order tensor field.
    subroutine fouriergrid_to_tensor4(tensor4, npts1, npts2, npts3)
        use fourier_functions
        implicit none
        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        double precision, intent(out) :: tensor4(3, 3, 3, 3, npts1, npts2, npts3) ! real-space 4th-order tensor field (output)
        integer :: i, j, k, l, m, n, o ! loop indices

        do k = 1, npts3
        do j = 1, npts2
        do i = 1, npts1
        do l = 1, 3
        do m = 1, 3
        do n = 1, 3
        do o = 1, 3
            tensor4(l, m, n, o, i, j, k) = real(fourgrid3333(k, j, i, o, n, m, l))
        end do
        end do
        end do
        end do
        end do
        end do
        end do

    end

    ! (Optional) Configure OpenMP threading (currently disabled via preprocessor directives).
    subroutine setup_openmp(nthreads)
!$      use omp_lib
        implicit none
        integer, intent(inout) :: nthreads ! requested/actual OpenMP thread count
!$      if (nthreads == 0) then
!$          nthreads = omp_get_max_threads()
!$      end if
!$      call omp_set_num_threads(nthreads)
!$      write (*, '(A, I0, A)') 'Running OMP parallel on ', nthreads, ' threads'
    end subroutine setup_openmp

    ! Write a 3-component vector field (cell-centered) to an ASCII file for debugging/postprocessing.
    subroutine write_tensor1_to_file(filename, tensor1, npts1, npts2, npts3)
        implicit none
        character(len=256), intent(in) :: filename ! output file path
        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        double precision, intent(in) :: tensor1(3, npts1, npts2, npts3) ! vector field to write
        integer :: i, j, k, l ! loop indices

        open (13, file=filename, status='replace')
        do i = 1, npts1
        do j = 1, npts2
        do k = 1, npts3
            write (13, '(3I10)') i, j, K
            write (13, '(3E16.8)') (tensor1(l, i, j, k), l=1, 3)
        end do
        end do
        end do
        close (13)

    end

    ! Write a 3x3 tensor field (cell-centered) to an ASCII file for debugging/postprocessing.
    subroutine write_tensor2_to_file(filename, tensor2, npts1, npts2, npts3)
        implicit none
        character(len=256), intent(in) :: filename ! output file path
        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        double precision, intent(in) :: tensor2(3, 3, npts1, npts2, npts3) ! 2nd-order tensor field to write
        integer :: i, j, k, l, m ! loop indices

        open (13, file=filename, status='replace')
        do i = 1, npts1
        do j = 1, npts2
        do k = 1, npts3
            write (13, '(3I10)') i, j, K
            write (13, '(3E16.8)') ((tensor2(l, m, i, j, k), m=1, 3), l=1, 3)
        end do
        end do
        end do
        close (13)

    end

    ! Write a generic idim1-component field (cell-centered) to an ASCII file.
    subroutine write_array1_to_file(filename, array1, idim1, npts1, npts2, npts3)
        implicit none
        character(len=256), intent(in) :: filename ! output file path
        integer, intent(in) :: idim1, npts1, npts2, npts3 ! component count and grid dimensions
        double precision, intent(in) :: array1(idim1, npts1, npts2, npts3) ! array field to write
        integer :: i, j, k, l ! loop indices
        character(len=3) :: str ! formatted component count string
        character(len=20) :: format ! Fortran format string for output

        write (str, '(I3.3)') idim1
        format = '('//str//'E24.16)'

        open (13, file=filename, status='replace')
        do i = 1, npts1
        do j = 1, npts2
        do k = 1, npts3
            write (13, '(3I10)') i, j, k
            write (13, format) (array1(l, i, j, k), l=1, idim1)
        end do
        end do
        end do
        close (13)

    end

    ! Write a generic idim1-by-idim2 field (cell-centered) to an ASCII file.
    subroutine write_array2_to_file(filename, array2, idim1, idim2, npts1, npts2, npts3)
        implicit none
        character(len=256), intent(in) :: filename ! output file path
        integer, intent(in) :: idim1, idim2, npts1, npts2, npts3 ! component counts and grid dimensions
        double precision, intent(in) :: array2(idim1, idim2, npts1, npts2, npts3) ! matrix-valued field to write
        integer :: i, j, k, l, m ! loop indices
        character(len=3) :: str ! formatted inner dimension string
        character(len=20) :: format ! Fortran format string for output

        write (str, '(I3.3)') idim2
        format = '('//str//'E24.16)'

        open (13, file=filename, status='replace')
        do i = 1, npts1
        do j = 1, npts2
        do k = 1, npts3
            write (13, '(3I10)') i, j, k
            do l = 1, idim1
                write (13, format) (array2(l, m, i, j, k), m=1, idim2)
            end do
        end do
        end do
        end do
        close (13)

    end

    ! Integrate a 3x3 gradient field to a 3-component vector field in Fourier space (inverse divergence).
    subroutine integrate_tensor2(grad_tensor1, tensor1, ideriv, npts1, npts2, npts3)
        use fourier_functions
        implicit none

        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        integer, intent(in) :: ideriv ! derivative operator selector (spectral)
        double precision, intent(in) :: grad_tensor1(3, 3, npts1, npts2, npts3) ! gradient field to integrate
        double precision, intent(out) :: tensor1(3, npts1, npts2, npts3) ! integrated vector field (up to constant)
        double precision :: xk(3) ! spatial frequency vector at (i,j,k)
        complex(kind(1.0d0)) :: kmod(3) ! complex wavevector modifier for derivative
        complex(kind(1.0d0)) :: dumcmplx ! scratch complex scalar
        integer :: i, j, k, l, m ! loop indices

        call tensor2_to_fouriergrid(grad_tensor1, npts1, npts2, npts3)
        call fft_tensor33(plan_advanced33, fourgrid33)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(i, j, k, l, m, xk, kmod, dumcmplx) &
        !$OMP SHARED(npts1, npts2, npts3, ideriv, fourgrid33, fourgrid3)
        !$OMP DO COLLAPSE(3)
        do k = 1, npts3
        do j = 1, npts2
        do i = 1, npts1

            xk = spatial_freq(i, j, k, npts1, npts2, npts3)
            if (ideriv == 2) kmod = cmplx(0.0, 1.0) * 0.25 * &
                                    (1.0 + exp(cmplx(0.0, xk(1), kind=kind(0.0d0)))) * (1.0 + exp(cmplx(0.0, xk(2), kind=kind(0.0d0)))) * &
                                    (1.0 + exp(cmplx(0.0, xk(3), kind=kind(0.0d0))))

            do l = 1, 3
                if (ideriv == 0) then
                    kmod(l) = cmplx(0.0, xk(l), kind=kind(0.0d0))
                elseif (ideriv == 1) then
                    kmod(l) = cmplx(0.0, sin(xk(l)), kind=kind(0.0d0))
                elseif (ideriv == 2) then
                    kmod(l) = tan(xk(l) / 2.0) * kmod(l)
                end if
            end do

            if (real(dot_product(kmod, kmod)) > 0.0) then

                do l = 1, 3
                    dumcmplx = (0.0, 0.0)
                    do m = 1, 3
                        dumcmplx = dumcmplx + conjg(kmod(m)) * fourgrid33(k, j, i, m, l) / &
                                   dot_product(kmod, kmod)
                    end do
                    fourgrid3(k, j, i, l) = dumcmplx
                end do

            else

                fourgrid3(k, j, i, :) = (0.0, 0.0)

            end if

            if (((mod(npts1, 2) == 0 .and. i == npts1 / 2 + 1) .and. &
                 (mod(npts2, 2) == 0 .and. j == npts2 / 2 + 1) .and. &
                 ((mod(npts3, 2) == 0 .and. k == npts3 / 2 + 1) .and. npts3 > 1)) .and. ideriv == 2) fourgrid3(k, j, i, :) = (0.0, 0.0)

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        call ifft_tensor3(iplan_advanced3, fourgrid3)

        call fouriergrid_to_tensor1(tensor1, npts1, npts2, npts3)

    end

    ! Update current voxel coordinates xgr by integrating the deformation gradient field.
    subroutine update_grid_defgrad(micro, props)
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (xgr updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: dx(3, props%npts1, props%npts2, props%npts3) ! integrated displacement fluctuation field
        double precision :: x(3) ! reference-space voxel center coordinate
        double precision :: defgrad(3, 3, props%npts1, props%npts2, props%npts3) ! defgrad field (cell-centered)
        integer :: ideriv ! derivative operator selector
        integer :: i, j, k ! loop indices

        ideriv = props%process(1)%igamma
        do k = 1, props%npts3
        do j = 1, props%npts2
        do i = 1, props%npts1
            defgrad(:, :, i, j, k) = micro%voxel(i, j, k)%defgrad
        end do
        end do
        end do
        call integrate_tensor2(defgrad, dx, ideriv, props%npts1, props%npts2, props%npts3)

        do k = 1, props%npts3
        do j = 1, props%npts2
        do i = 1, props%npts1

            x(1) = float(i)
            x(2) = float(j)
            x(3) = float(k)

            micro%voxel(i, j, k)%xgr(:) = matmul(micro%defgradavg, x) + dx(:, i, j, k)

        end do
        end do
        end do

    end

    ! Update per-voxel Schmid tensors based on current elastic defgrad and orientation.
    subroutine update_schmid(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (Schmid tensors updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: aux5(5) ! scratch: Schmid tensor in 5-component basis
        double precision :: aux33(3, 3) ! scratch: 3x3 Schmid tensor
        double precision :: aux55(5, 5) ! scratch: 5x5 matrix (basis transform)
        double precision :: aux3333(3, 3, 3, 3) ! scratch: 4th-order tensor (basis transform)
        double precision :: aux33r(3, 3) ! rotated Schmid tensor (current config)
        double precision :: aa(3, 3) ! crystal->current mapping (Fe*R)
        double precision :: aainv(3, 3) ! inverse of aa
        double precision :: aux33sym(3, 3) ! symmetric part of rotated Schmid tensor
        double precision :: aux33asym(3, 3) ! antisymmetric part (unused)
        integer :: i, j, k ! voxel indices
        integer :: ii, jj ! tensor index loops
        integer :: is ! slip system index
        integer :: jph ! phase index for voxel

        ! write (*, '(A)') 'Updating Schmid tensors'
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(i, j, k, ii, jj, jph, is, aux33, aa, aainv, aux33r, aux33sym, aux33asym) &
        !$OMP PRIVATE(aux5, aux55, aux3333) &
        !$OMP SHARED(props, micro)
        !$OMP DO COLLAPSE(3)
        do i = 1, props%npts1
        do j = 1, props%npts2
        do k = 1, props%npts3
            jph = micro%voxel(i, j, k)%phase
            if (.not. micro%voxel(i, j, k)%gas) then

                do is = 1, props%phase(jph)%nsyst

                    ! aux5=schca(:,is,jph)
                    ! call chg_basis(aux5,aux33,aux55,aux3333,1,5)

                    do ii = 1, 3
                    do jj = 1, 3
                        aux33(ii, jj) = props%phase(jph)%dbca(ii, is) * props%phase(jph)%dnca(jj, is) ! total Schmid tensor
                    end do
                    end do

                    ! faster calculation
                    aa = matmul(micro%voxel(i, j, k)%defgrade(:, :), micro%voxel(i, j, k)%ag(:, :))
                    aainv = aa
                    call lu_inverse(aainv, 3)

                    aux33r = matmul(matmul(aa, aux33), aainv)

                    aux33sym = sym33(aux33r)
                    aux33asym = antisym33(aux33r)

                    call chg_basis(aux5, aux33sym, aux55, aux3333, 2, 5)
                    micro%voxel(i, j, k)%sch(:, is) = aux5

                end do
            end if ! igas endif
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        return
    end

    ! Compute the weighted average of a tensor field over the full grid.
    subroutine average_tensor2(tensor2, wgtc, tensor2avg, dim1, dim2, npts1, npts2, npts3)
!$      use omp_lib
        implicit none

        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        integer, intent(in) :: dim1, dim2 ! leading tensor dimensions
        double precision, intent(in) :: tensor2(dim1, dim2, npts1, npts2, npts3) ! tensor field to average
        double precision, intent(in) :: wgtc(npts1, npts2, npts3) ! weights (typically voxel volume fractions)
        double precision, intent(out) :: tensor2avg(dim1, dim2) ! weighted average tensor
        integer :: i, j, k ! loop indices

        tensor2avg = 0.0
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(i,j,k) &
        !$OMP SHARED(npts1,npts2,npts3,tensor2,wgtc) &
        !$OMP REDUCTION(+:tensor2avg)
        !$OMP DO COLLAPSE(3)
        do i = 1, npts1
        do j = 1, npts2
        do k = 1, npts3
            tensor2avg = tensor2avg + tensor2(:, :, i, j, k) * wgtc(i, j, k)
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end

    ! Compute the weighted average of a tensor field over masked voxels only.
    subroutine average_tensor2_mask(tensor2, wgtc, mask, tensor2avg, dim1, dim2, npts1, npts2, npts3)
!$      use omp_lib
        implicit none

        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z
        integer, intent(in) :: dim1, dim2 ! leading tensor dimensions
        double precision, intent(in) :: tensor2(dim1, dim2, npts1, npts2, npts3) ! tensor field to average
        double precision, intent(in) :: wgtc(npts1, npts2, npts3) ! weights (typically voxel volume fractions)
        logical, intent(in) :: mask(npts1, npts2, npts3) ! mask selecting voxels to include
        double precision, intent(out) :: tensor2avg(dim1, dim2) ! masked weighted average tensor
        double precision :: wtot ! total weight of masked region
        integer :: i, j, k ! loop indices

        tensor2avg = 0.0
        wtot = 0.0
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(i,j,k) &
        !$OMP SHARED(npts1,npts2,npts3,tensor2,wgtc,mask) &
        !$OMP REDUCTION(+:tensor2avg, wtot)
        !$OMP DO COLLAPSE(3)
        do i = 1, npts1
        do j = 1, npts2
        do k = 1, npts3
            if (mask(i, j, k)) then
                tensor2avg = tensor2avg + tensor2(:, :, i, j, k) * wgtc(i, j, k)
                wtot = wtot + wgtc(i, j, k)
            end if
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        if (wtot > 0.0) tensor2avg = tensor2avg / wtot

    end

    ! Build the Fourier-space Green operator for the current reference stiffness and derivative scheme.
    subroutine form_G(micro, props)
        use tensor_functions
        use fourier_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (G operator written)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: xk(3) ! spatial frequency vector at (kxx,kyy,kzz)
        double precision :: g1(3, 3, 3, 3) ! Fourier-space Green operator at this k
        double precision :: D_dft(3, 3) ! discrete/continuous derivative operator (k_i k_j form)
        double precision :: a(3, 3) ! acoustic tensor inverse (in Fourier space)
        double precision :: g1tmp(3, 3, 3, 3) ! scratch Green operator (fixed-point smoothing)
        double precision :: aux3333(3, 3, 3, 3) ! scratch 4th-order tensor

        integer :: i, j, k, l, m, n, o, p, q ! tensor index loops
        integer :: ii, jj ! index loops (1..3)
        integer :: itt ! fixed-point iteration counter (igamma=1)
        integer :: kxx, kyy, kzz ! Fourier grid indices (x,y,z)
        double precision :: dum ! scratch scalar
        double precision :: r ! fixed-point residual metric (igamma=1)
        double precision :: rold ! previous residual metric (unused)
        double precision :: tol ! fixed-point tolerance (igamma=1)
        double precision :: xisum ! 2*pi*(kxx-1)/npts1 (discrete operator)
        double precision :: xjsum ! 2*pi*(kyy-1)/npts2 (discrete operator)
        double precision :: xksum ! 2*pi*(kzz-1)/npts3 (discrete operator)
        double precision :: xmix ! mixing factor for fixed-point smoothing
        complex(kind(1.0d0)) :: kmod(3) ! modified complex wavevector for igamma=2
        complex(kind(1.0d0)) :: D_dft_c(3, 3) ! complex D operator before taking real part

        do kzz = 1, props%npts3
        do kyy = 1, props%npts2
        do kxx = 1, props%npts1

            xk = spatial_freq(kxx, kyy, kzz, props%npts1, props%npts2, props%npts3)

            if (((mod(props%npts1, 2) == 0 .and. kxx == props%npts1 / 2 + 1) .or. &
                 (mod(props%npts2, 2) == 0 .and. kyy == props%npts2 / 2 + 1) .or. &
                 ((mod(props%npts3, 2) == 0 .and. kzz == props%npts3 / 2 + 1) .and. props%npts3 > 1)) .and. &
                props%process(ipr)%igamma == 0) then

                g1 = -micro%s0

            elseif ((((mod(props%npts1, 2) == 0 .and. kxx == props%npts1 / 2 + 1) .and. &
                      (mod(props%npts2, 2) == 0 .and. kyy == props%npts2 / 2 + 1)) .or. &
                     ((mod(props%npts1, 2) == 0 .and. kxx == props%npts1 / 2 + 1) .and. &
                      (mod(props%npts3, 2) == 0 .and. kzz == props%npts3 / 2 + 1)) .or. &
                     ((mod(props%npts2, 2) == 0 .and. kyy == props%npts2 / 2 + 1) .and. &
                      (mod(props%npts3, 2) == 0 .and. kzz == props%npts3 / 2 + 1))) .and. &
                    props%process(ipr)%igamma >= 2) then ! 3d

                g1 = 0.0

            else

                ! Build discrete and continuous D operator
                if (props%process(ipr)%igamma == 0) then

                    do ii = 1, 3
                    do jj = 1, 3
                        D_dft(ii, jj) = xk(ii) * xk(jj)
                    end do
                    end do

                elseif (props%process(ipr)%igamma == 1) then

                    xisum = 2.*pi * float(kxx - 1) / float(props%npts1)
                    xjsum = 2.*pi * float(kyy - 1) / float(props%npts2)
                    xksum = 2.*pi * float(kzz - 1) / float(props%npts3)

                    D_dft(1, 1) = 2.0 * (cos(xisum) - 1.0)
                    D_dft(2, 2) = 2.0 * (cos(xjsum) - 1.0)
                    D_dft(3, 3) = 2.0 * (cos(xksum) - 1.0)
                    D_dft(2, 1) = -sin(xisum) * sin(xjsum) ! 0.5*(cos(xisum+xjsum)-cos(xisum-xjsum))
                    D_dft(3, 1) = -sin(xisum) * sin(xksum) ! 0.5*(cos(xisum+xksum)-cos(xisum-xksum))
                    D_dft(3, 2) = -sin(xjsum) * sin(xksum) ! 0.5*(cos(xjsum+xksum)-cos(xjsum-xksum))
                    D_dft(1, 2) = D_dft(2, 1)
                    D_dft(1, 3) = D_dft(3, 1)
                    D_dft(2, 3) = D_dft(3, 2)

                elseif (props%process(ipr)%igamma == 2) then

                    kmod = cmplx(0.0, 1.0) * 0.25 * &
                           (1.0 + exp(cmplx(0.0, xk(1), kind=kind(0.0d0)))) * (1.0 + exp(cmplx(0.0, xk(2), kind=kind(0.0d0)))) * &
                           (1.0 + exp(cmplx(0.0, xk(3), kind=kind(0.0d0))))

                    do ii = 1, 3
                        kmod(ii) = tan(xk(ii) / 2.0) * kmod(ii)
                    end do

                    do ii = 1, 3
                    do jj = 1, 3
                        D_dft_c(ii, jj) = kmod(ii) * conjg(kmod(jj))
                    end do
                    end do
                    if (norm2(aimag(D_dft_c)) > 1.0e-7) then
                        write (*, *) 'norm2(aimag(D_dft_c))'
                        write (*, *) norm2(aimag(D_dft_c))
                        stop
                    end if
                    D_dft = real(D_dft_c)

                end if

                ! Compute acoustic tensor inverse
                if (kxx == 1 .and. kyy == 1 .and. kzz == 1) then
                    a = 0.
                else
                    do i = 1, 3
                    do k = 1, 3
                        a(i, k) = 0.
                        do j = 1, 3
                        do l = 1, 3
                            a(i, k) = a(i, k) + micro%c0(i, j, k, l) * D_dft(j, l)!*xk(j)*xk(l)
                        end do
                        end do
                    end do
                    end do

                    call lu_inverse(a, 3)
                end if

                do p = 1, 3
                do q = 1, 3
                do i = 1, 3
                do j = 1, 3
                    g1(p, q, i, j) = -a(p, i) * D_dft(q, j) !*xk(q)*xk(j)
                end do
                end do
                end do
                end do

            end if

            if (props%process(ipr)%igamma == 1) then

                g1tmp = g1
                tol = 1.0d-50
                r = 2.0 * tol
                rold = 2.0e2
                itt = 0
                xmix = 1.0
                do while (r > tol .and. itt < 4)
                    itt = itt + 1

                    do i = 1, 3
                    do j = 1, 3
                    do k = 1, 3
                    do l = 1, 3
                        dum = 0.0
                        do m = 1, 3
                        do n = 1, 3
                        do o = 1, 3
                        do p = 1, 3
                            dum = dum - g1(i, j, m, n) * micro%c0(m, n, o, p) * g1tmp(o, p, k, l)
                        end do
                        end do
                        end do
                        end do
                        aux3333(i, j, k, l) = dum
                    end do
                    end do
                    end do
                    end do

                    g1tmp = xmix * aux3333 + (1.0 - xmix) * g1tmp

                end do
                g1 = g1tmp

            end if

            micro%Goper(:, :, :, :, kxx, kyy, kzz) = g1

            do i = 1, 3
            do j = 1, 3
            do k = 1, 3
            do l = 1, 3
                fourgrid3333(kzz, kyy, kxx, l, k, j, i) = cmplx(micro%Goper(i, j, k, l, kxx, kyy, kzz), kind=kind(0.0d0))
            end do
            end do
            end do
            end do

        end do
        end do
        end do

        call ifft_tensor3333(iplan_advanced3333, fourgrid3333)
        do i = 1, 3
        do j = 1, 3
        do k = 1, 3
        do l = 1, 3
            micro%Goperr0(i, j, k, l) = real(fourgrid3333(1, 1, 1, l, k, j, i))
        end do
        end do
        end do
        end do

    end

    ! Update per-voxel mapped stiffness and (optionally) modified-AL stabilization operator terms.
    subroutine calc_c066mod_Goperr066mod(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (per-voxel operators updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: dum ! scratch scalar
        double precision :: c0mod(3, 3, 3, 3) ! push-forwarded macro stiffness to voxel config
        double precision :: c0modsym(3, 3, 3, 3) ! symmetrized c0mod
        double precision :: c066mod(6, 6) ! c0modsym in 6x6 basis
        double precision :: Goperr0mod(3, 3, 3, 3) ! modified real-space G at origin (current config)
        double precision :: Goperr0modsym(3, 3, 3, 3) ! symmetrized Goperr0mod
        double precision :: Goperr066mod(6, 6) ! Goperr0modsym in 6x6 basis
        double precision :: aux6(6) ! scratch 6-vector
        double precision :: aux33(3, 3) ! scratch 3x3 tensor
        integer :: i, j, k ! voxel indices
        integer :: ii, jj, kk, ll ! tensor index loops (1..3)
        integer :: m, n ! tensor contraction indices (1..3)

        ! Modify the reference stiffness
        do i = 1, props%npts1
        do j = 1, props%npts2
        do k = 1, props%npts3
            do ii = 1, 3
            do jj = 1, 3
            do kk = 1, 3
            do ll = 1, 3
                dum = 0.0
                do m = 1, 3
                do n = 1, 3
                    dum = dum + micro%c0(ii, m, kk, n) * micro%voxel(i, j, k)%defgrad(jj, m) * &
                          micro%voxel(i, j, k)%defgrad(ll, n)
                end do
                end do
                c0mod(ii, jj, kk, ll) = dum / micro%voxel(i, j, k)%detF
            end do
            end do
            end do
            end do
            do ii = 1, 3
            do jj = 1, 3
            do kk = 1, 3
            do ll = 1, 3
                c0modsym(ii, jj, kk, ll) = 0.25 * (c0mod(ii, jj, kk, ll) + c0mod(jj, ii, kk, ll) + &
                                                   c0mod(ii, jj, ll, kk) + c0mod(jj, ii, ll, kk)) ! symmetrize
            end do
            end do
            end do
            end do
            call chg_basis(aux6, aux33, c066mod, c0modsym, 4, 6)
            micro%voxel(i, j, k)%c066mod(:, :) = c066mod

            ! modify Goperr0
            if (props%modAL) then
                do ii = 1, 3
                do jj = 1, 3
                do kk = 1, 3
                do ll = 1, 3
                    dum = 0.0
                    do m = 1, 3
                    do n = 1, 3
                        dum = dum + micro%Goperr0(ii, m, kk, n) * &
                              micro%voxel(i, j, k)%defgradinv(m, jj) * &
                              micro%voxel(i, j, k)%defgradinv(n, ll)
                    end do
                    end do
                    Goperr0mod(ii, jj, kk, ll) = dum * micro%voxel(i, j, k)%detF
                end do
                end do
                end do
                end do
                do ii = 1, 3
                do jj = 1, 3
                do kk = 1, 3
                do ll = 1, 3
                    Goperr0modsym(ii, jj, kk, ll) = 0.25 * (Goperr0mod(ii, jj, kk, ll) + &
                                                            Goperr0mod(jj, ii, kk, ll) + &
                                                            Goperr0mod(ii, jj, ll, kk) + &
                                                            Goperr0mod(jj, ii, ll, kk)) ! symmetrize
                end do
                end do
                end do
                end do
                call chg_basis(aux6, aux33, Goperr066mod, Goperr0modsym, 4, 6)

                micro%voxel(i, j, k)%c066modGoperr066mod(:, :) = matmul(c066mod, Goperr066mod)
            end if

        end do
        end do
        end do

    end

    ! Build the gas-phase compliance operator used to solve stresses in dummy gas voxels.
    subroutine calc_gas_compl(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (gas compliance operator updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: aux66(6, 6) ! scratch 6x6 matrix (gas operator)
        integer :: i, j, k, ii ! loop indices

        do i = 1, props%npts1
        do j = 1, props%npts2
        do k = 1, props%npts3

            if (micro%voxel(i, j, k)%gas) then

                aux66 = 0.0
                do ii = 1, 6
                    aux66(ii, ii) = 1.0
                end do
                aux66 = aux66 + micro%voxel(i, j, k)%c066mod(:, :) / props%process(ipr)%tdot * &
                        props%complgas

                if (props%modAL) then
                    aux66 = aux66 - micro%voxel(i, j, k)%c066modGoperr066mod(:, :)
                end if

                call lu_inverse(aux66, 6)

                micro%voxel(i, j, k)%cgas(:, :) = aux66

            else

                micro%voxel(i, j, k)%cgas(:, :) = 0.0

            end if

        end do
        end do
        end do

    end

    ! Convert Cauchy stress to first Piola-Kirchhoff stress: P = det(F) * sigma * F^{-T}.
    function CauchyToPK1(Cauchy, Finv, detF) result(PK1)
        implicit none
        double precision, intent(in), dimension(3, 3) :: Cauchy ! Cauchy stress (3x3)
        double precision, intent(in), dimension(3, 3) :: Finv ! inverse deformation gradient F^{-1}
        double precision, intent(in) :: detF ! det(F) Jacobian
        double precision, dimension(3, 3) :: PK1 ! 1st Piola-Kirchhoff stress (3x3)

        PK1 = matmul(Cauchy, transpose(Finv)) * detF

    end function

    ! Apply the Green operator convolution in Fourier space to obtain velgrad corrections.
    subroutine convolution_with_Goper(micro, props)
        use tensor_functions
        use fourier_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (dvelgradref computed)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        complex(kind(1.0d0)) :: aux33cmplx(3, 3) ! scratch complex 3x3 (FFT space)
        complex(kind(1.0d0)) :: dumcmplx ! scratch complex scalar accumulator
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)
        integer :: i, j, k, l ! tensor index loops (1..3)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, i, j) &
        !$OMP SHARED(props, micro, fourgrid33)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1
            do i = 1, 3
            do j = 1, 3
                fourgrid33(ip3, ip2, ip1, j, i) = micro%voxel(ip1, ip2, ip3)%sgPK1(i, j)
            end do
            end do
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        call fft_tensor33(plan_advanced33, fourgrid33)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, i, j, k, l, aux33cmplx, dumcmplx) &
        !$OMP SHARED(props, micro, fourgrid33)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            if (ip1 == 1 .and. ip2 == 1 .and. ip3 == 1) then

                fourgrid33(ip3, ip2, ip1, :, :) = (0., 0.)

            else

                do j = 1, 3
                do i = 1, 3
                    aux33cmplx(i, j) = fourgrid33(ip3, ip2, ip1, j, i)
                end do
                end do

                do i = 1, 3
                do j = 1, 3
                    dumcmplx = (0., 0.)
                    do k = 1, 3
                    do l = 1, 3
                        dumcmplx = dumcmplx + micro%Goper(i, j, k, l, ip1, ip2, ip3) * aux33cmplx(k, l)
                    end do
                    end do
                    fourgrid33(ip3, ip2, ip1, j, i) = dumcmplx
                end do
                end do

            end if

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        call ifft_tensor33(iplan_advanced33, fourgrid33)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, i, j, k, l) &
        !$OMP SHARED(props, micro, fourgrid33)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1
            do i = 1, 3
            do j = 1, 3
                micro%voxel(ip1, ip2, ip3)%dvelgradref(i, j) = real(fourgrid33(ip3, ip2, ip1, j, i))
            end do
            end do
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end

    ! Map reference-config corrections back to current config and update per-voxel velgrad fields.
    subroutine correction_to_current(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (velgrad updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: velgrad(3, 3, props%npts1, props%npts2, props%npts3) ! velgrad fluctuations (current config)
        double precision :: velgradfluctavg(3, 3) ! weighted average velgrad fluctuation
        double precision :: velgradavg(3, 3) ! weighted average total velgrad
        double precision :: defgradinvavgc_inv(3, 3) ! inverse of avg(F^{-1}) (all voxels)
        double precision :: dudot_dX_avg(3, 3) ! BC correction term in reference gradient
        double precision :: velgradavggas(3, 3) ! average velgrad in gas voxels
        double precision :: velgradavgsol(3, 3) ! average velgrad in solid voxels
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)

        velgradfluctavg = 0.0
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3) &
        !$OMP SHARED(props, micro, velgrad) &
        !$OMP REDUCTION(+:velgradfluctavg)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1
            velgrad(:, :, ip1, ip2, ip3) = matmul(micro%voxel(ip1, ip2, ip3)%dvelgradref, &
                                                  micro%voxel(ip1, ip2, ip3)%defgradinv)
            velgradfluctavg = velgradfluctavg + micro%voxel(ip1, ip2, ip3)%wgtc * &
                              velgrad(:, :, ip1, ip2, ip3)
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        defgradinvavgc_inv = micro%defgradinvavgc
        call lu_inverse(defgradinvavgc_inv, 3)

        dudot_dX_avg = matmul(micro%dvelgradavg - velgradfluctavg, defgradinvavgc_inv)

        velgradavg = 0.0
        velgradavggas = 0.0
        velgradavgsol = 0.0
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3) &
        !$OMP SHARED(props, micro, dudot_dX_avg) &
        !$OMP REDUCTION(+:velgradavg, velgradavggas, velgradavgsol)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1
            micro%voxel(ip1, ip2, ip3)%velgrad = micro%voxel(ip1, ip2, ip3)%velgrad(:, :) + &
                                                 matmul(micro%voxel(ip1, ip2, ip3)%dvelgradref(:, :) + dudot_dX_avg, &
                                                        micro%voxel(ip1, ip2, ip3)%defgradinv(:, :))

            velgradavg = velgradavg + micro%voxel(ip1, ip2, ip3)%velgrad(:, :) * &
                         micro%voxel(ip1, ip2, ip3)%wgtc
            if (micro%voxel(ip1, ip2, ip3)%gas) then
                velgradavggas = velgradavggas + micro%voxel(ip1, ip2, ip3)%velgrad(:, :) * &
                                micro%voxel(ip1, ip2, ip3)%wgtc
            else
                velgradavgsol = velgradavgsol + micro%voxel(ip1, ip2, ip3)%velgrad(:, :) * &
                                micro%voxel(ip1, ip2, ip3)%wgtc
            end if
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        if (micro%wphcgas > 0.0) velgradavggas = velgradavggas / micro%wphcgas
        velgradavgsol = velgradavgsol / micro%wphcsol

        micro%velgradavg = velgradavg
        micro%velgradavggas = velgradavggas
        micro%velgradavgsol = velgradavgsol

    end

    ! Per-voxel constitutive solve: given velgrad, update stress and plastic strain-rate (AL/Newton).
    subroutine solve_res(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (stress/edotp updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions
        type(voxel_type) :: voxel_loc ! local voxel copy for constitutive solve
        integer, parameter :: itmaxal = 100 ! max Newton iterations (inner augmented-Lagrangian)
        integer, parameter :: itmaxaltau = 10 ! max hardening fixed-point iterations (outer loop)

        double precision :: lambda(6) ! stress-like multiplier in 6-component basis
        double precision :: sg(6) ! local stress unknown in 6-component basis
        double precision :: sgt(6) ! previous-step stress in 6-component basis
        double precision :: velgrad_sym33(3, 3) ! symmetric part of voxel velgrad (3x3)
        double precision :: velgrad_sym(6) ! symmetric velgrad in 6-component basis
        double precision :: edotth(6) ! eigenstrain/thermal strain-rate term in 6-basis
        double precision :: sgnorm ! norm of stress for normalization
        double precision :: dgnorm ! norm of strain-rate for normalization
        double precision :: taunrm ! norm of trial tau for normalization
        double precision :: crssnrm ! norm of CRSS for normalization (unused)
        double precision :: errtau ! hardening fixed-point error
        double precision :: itertauavg ! average outer-loop iterations (accumulator)
        double precision :: erral ! augmented-Lagrangian convergence measure
        double precision :: resoldn ! previous residual norm (line search)
        double precision :: iterNRavg ! average Newton iterations (accumulator)
        double precision :: edotp(6) ! plastic strain-rate in 6-component basis
        double precision :: dedotp_dsg(6, 6) ! d(edotp)/d(sg) in 6-basis
        double precision :: edote(6) ! elastic strain-rate in 6-component basis
        double precision :: dedote_dsg(6, 6) ! d(edote)/d(sg) in 6-basis
        double precision :: edot(6) ! total strain-rate in 6-component basis
        double precision :: dedot_dsg(6, 6) ! d(edot)/d(sg) in 6-basis
        double precision :: res(6) ! residual vector for local solve
        double precision :: jacobian(6, 6) ! Jacobian for Newton update
        double precision :: dsgnormNR ! Newton update norm
        double precision :: dsgnorm ! stress mismatch norm
        double precision :: ddgnorm ! strain-rate mismatch norm
        double precision :: errald ! strain-rate-based convergence measure
        double precision :: trialtau(nsysmx, 2) ! updated slip resistances (fwd/bwd)
        double precision :: trialtauold(nsysmx, 2) ! previous slip resistances (fwd/bwd)
        double precision :: errtauold ! previous hardening error (unused)
        double precision :: aux33(3, 3) ! scratch 3x3 tensor (basis transforms)
        double precision :: aux6(6) ! scratch 6-vector (basis transforms)
        double precision :: aux3333(3, 3, 3, 3) ! scratch 4th-order tensor (basis transforms)
        double precision :: aux66(6, 6) ! scratch 6x6 matrix (basis transforms)
        double precision :: resn ! residual norm (line search)
        double precision :: sgold(6) ! previous sg iterate (line search)
        double precision :: gamacum ! accumulated slip from hardening update
        integer :: itertau ! outer hardening iteration counter
        integer :: iterNR ! Newton iteration counter
        integer :: iunconv ! count of voxels not converged within itmaxal
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)
        integer :: icut ! line-search cut counter

        itertauavg = 0.0
        erre = 0.0
        errs = 0.0
        iterNRavg = 0.0
        itertauavg = 0.0
        iunconv = 0

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED(props, micro, imicro, ipr) &
        !$OMP REDUCTION(+:erre, errs, itertauavg, iterNRavg, iunconv) &
        !$OMP PRIVATE(ip1, ip2, ip3, voxel_loc, lambda, aux33, aux6, aux3333, aux66) &
        !$OMP PRIVATE(sg, sgt, velgrad_sym33, velgrad_sym, edotth, sgnorm, dgnorm, taunrm, crssnrm) &
        !$OMP PRIVATE(itertau, errtau, iterNR, erral, resoldn, icut, resn, edotp, dedotp_dsg) &
        !$OMP PRIVATE(edote, dedote_dsg, edot, dedot_dsg, res, sgold, jacobian, dsgnormNR) &
        !$OMP PRIVATE(dsgnorm, ddgnorm, errald, errtauold, trialtauold, gamacum, trialtau)
        !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            ! local array
            voxel_loc = micro%voxel(ip1, ip2, ip3)
            call chg_basis(lambda, voxel_loc%sg, aux66, aux3333, 2, 6)
            sg = lambda
            call chg_basis(sgt, voxel_loc%sgt, aux66, aux3333, 2, 6)
            velgrad_sym33 = sym33(voxel_loc%velgrad)
            call chg_basis(velgrad_sym, velgrad_sym33, aux66, aux3333, 2, 6)

            if (.not. voxel_loc%gas) then

                if (imicro == 1 .and. props%process(ipr)%ithermo == 1) then
                    call chg_basis(edotth, voxel_loc%eth(:, :), aux66, aux3333, 2, 6)
                end if

                ! norms
                sgnorm = norm2(sg)
                dgnorm = norm2(velgrad_sym)
                taunrm = norm2(voxel_loc%trialtau)
                crssnrm = norm2(voxel_loc%crss(:, 1))

                ! outer loop over slip resistances
                itertau = 0
                errtau = 2.0 * maxval([props%process(ipr)%error, props%process(ipr)%erroral])
                do while (itertau < itmaxaltau .and. errtau > &
                          maxval([props%process(ipr)%error, props%process(ipr)%erroral]))

                    itertau = itertau + 1
                    itertauavg = itertauavg + 1.0

                    iterNR = 0
                    erral = 2.0 * props%process(ipr)%erroral
                    resoldn = 1.0

                    ! inner loop over stress
                    do while (iterNR < itmaxal .and. erral > props%process(ipr)%erroral)

                        iterNR = iterNR + 1
                        iterNRavg = iterNRavg + 1.0

                        ! loop to find step length
                        icut = 0
                        resn = 2.0 * resoldn
                        do while (resn > resoldn .and. icut < 10)

                            ! constitutive response
                            call plastic_strain_rate(sg, voxel_loc, props, edotp, dedotp_dsg)
                            call elastic_strain_rate(sg, sgt, voxel_loc, props, edote, dedote_dsg)
                            edot = edotp + edote
                            dedot_dsg = dedotp_dsg + dedote_dsg
                            if (imicro == 1 .and. props%process(ipr)%ithermo == 1) edot = edot + edotth

                            ! residual
                            res = sg - lambda + matmul(voxel_loc%c066mod, edot - velgrad_sym)
                            if (props%modAL) res = res - matmul(voxel_loc%c066modGoperr066mod, sg - lambda)

                            resn = norm2(res)
                            if (iterNR == 1) resoldn = resn

                            if (resn > resoldn) sg = sgold + (sg - sgold) * 0.5

                            icut = icut + 1

                        end do ! end do while step length
                        resoldn = resn
                        sgold = sg

                        ! Build Jacobian
                        jacobian = id6 + matmul(voxel_loc%c066mod, dedot_dsg)
                        if (props%modAL) jacobian = jacobian - voxel_loc%c066modGoperr066mod

                        ! Calculate new stress by solving the system -[J][delt_sg] = [R]
                        call lu_eqsystem(jacobian, res, 6)
                        sg = sg - res

                        ! errors
                        dsgnormNR = sum((sg - sgold)**2)
                        dsgnorm = sum((sg - lambda)**2)
                        ddgnorm = sum((edot - velgrad_sym)**2)

                        erral = sqrt(dsgnormNR) / sgnorm
                        errald = sqrt(ddgnorm) / dgnorm

                    end do ! enddo for NR

                    if (iterNR >= itmaxal) iunconv = iunconv + 1

                    ! new guess for tau and error
                    errtauold = errtau
                    if (props%process(ipr)%iuphard == 1) then
                        trialtauold = voxel_loc%trialtau(:, :)
                        call harden(voxel_loc, props, trialtau, gamacum)

                        errtau = norm2(trialtauold - trialtau) / taunrm
                        voxel_loc%trialtau(:, :) = trialtau
                    else
                        errtau = 0.
                    end if

                end do ! enddo for fixed point

                call chg_basis(sg, voxel_loc%sg(:, :), aux66, aux3333, 1, 6)
                call chg_basis(edotp, voxel_loc%edotp(:, :), aux66, aux3333, 1, 6)
                micro%voxel(ip1, ip2, ip3)%edotp = voxel_loc%edotp
                micro%voxel(ip1, ip2, ip3)%sg = voxel_loc%sg
                micro%voxel(ip1, ip2, ip3)%gamdot(:) = voxel_loc%gamdot(:)
                micro%voxel(ip1, ip2, ip3)%trialtau(:, :) = voxel_loc%trialtau(:, :)

                errs = errs + dsgnorm * voxel_loc%wgtc
                erre = erre + ddgnorm * voxel_loc%wgtc

            else

                res = lambda + matmul(voxel_loc%c066mod, velgrad_sym + sgt / props%process(ipr)%tdot * props%complgas)
                if (props%modAL) res = res - matmul(voxel_loc%c066modGoperr066mod, lambda)
                sg = matmul(voxel_loc%cgas(:, :), res)

                ! verify
                ! res = sg + matmul(voxel_loc%c066mod(:,:), props%complgas*(sg-sgt)/props%process(ipr)%tdot) &
                !   - lambda - matmul(voxel_loc%c066mod(:,:), velgrad_sym)
                ! if (props%modAL) res = res-matmul(voxel_loc%c066modGoperr066mod(:,:), sg-xlambda)

                call chg_basis(sg, voxel_loc%sg(:, :), aux66, aux3333, 1, 6)

                voxel_loc%edotp(:, :) = 0.0
                voxel_loc%gamdot(:) = 0.0

                micro%voxel(ip1, ip2, ip3)%edotp = voxel_loc%edotp
                micro%voxel(ip1, ip2, ip3)%sg = voxel_loc%sg
                micro%voxel(ip1, ip2, ip3)%gamdot(:) = voxel_loc%gamdot(:)

            end if
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        iterNRavg = iterNRavg * props%wgt
        iterNRavg_accum = iterNRavg_accum + iterNRavg
        itertauavg = itertauavg * props%wgt

        if (erre < 0 .or. errs < 0) stop 'ERRE/ERRS negative -> convergence failed!'

        erre = sqrt(erre)
        errs = sqrt(errs)

        return

    end

    ! Crystal plasticity flow rule: compute plastic strain-rate and consistent tangent wrt stress.
    subroutine plastic_strain_rate(sg, voxel_loc, props, edotp, dedotp_dsg)
        use tensor_functions
        use types
        implicit none

        double precision, intent(in) :: sg(6) ! stress in 6-component basis
        type(voxel_type), intent(inout) :: voxel_loc ! voxel state (gamdot updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions
        double precision, intent(out) :: edotp(6), dedotp_dsg(6, 6) ! plastic rate and derivative wrt stress

        integer :: ii, is, isign, jj ! loop indices; slip system; forward/backward selector
        double precision :: rss, rss1, rss2 ! resolved shear stress and derivative terms
        double precision :: dum ! scratch scalar (rate prefactor)
        double precision :: nrs(nsysmx) ! rate sensitivity exponents for this phase

        nrs = props%phase(voxel_loc%phase)%nrs

        dedotp_dsg = 0.0
        edotp = 0.0

        do is = 1, props%phase(voxel_loc%phase)%nsyst

            ! Calculate resolved shear stress
            rss = voxel_loc%sch(1, is) * sg(1) + voxel_loc%sch(2, is) * sg(2) + &
                  voxel_loc%sch(3, is) * sg(3) + &
                  voxel_loc%sch(4, is) * sg(4) + voxel_loc%sch(5, is) * sg(5)

            if (rss - voxel_loc%kin(is) < 0.0) then
                isign = 2
            else
                isign = 1
            end if
            rss = (rss - voxel_loc%kin(is)) / voxel_loc%trialtau(is, isign)

            ! intermediaries
            dum = props%phase(voxel_loc%phase)%gamd0(is) * abs(rss)**(nrs(is) - 1) ! g*|rss|^(n-1) prefactor
            rss1 = nrs(is) * dum / voxel_loc%trialtau(is, isign)
            rss2 = rss * dum
            voxel_loc%gamdot(is) = rss2

            ! Calculate strain rate and strain rate derivative w.r.t stress
            do jj = 1, 5
                edotp(jj) = edotp(jj) + voxel_loc%sch(jj, is) * rss2
                do ii = 1, 5
                    dedotp_dsg(ii, jj) = dedotp_dsg(ii, jj) + voxel_loc%sch(ii, is) * &
                                         voxel_loc%sch(jj, is) * rss1
                end do
            end do

        end do

        return
    end

    ! Elastic update: compute elastic strain-rate from stress increment and compliance.
    subroutine elastic_strain_rate(sg, sgt, voxel_loc, props, edote, dedote_dsg)
        use tensor_functions
        use types
        implicit none

        double precision, intent(in) :: sg(6), sgt(6) ! current and previous-step stress (6-basis)
        type(voxel_type), intent(in) :: voxel_loc ! voxel elastic compliance (sg66)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions
        double precision, intent(out) :: edote(6), dedote_dsg(6, 6) ! elastic rate and derivative wrt stress

        edote = matmul(voxel_loc%sg66(:, :), (sg - sgt) / props%process(ipr)%tdot)
        dedote_dsg = voxel_loc%sg66(:, :) / props%process(ipr)%tdot

        return
    end

    ! Update slip resistances (CRSS) using a Voce-style hardening law and accumulated slip.
    subroutine harden(voxel_loc, props, tau, gamma1)
        use tensor_functions
        use types
        implicit none

        type(voxel_type), intent(in) :: voxel_loc ! voxel state (gamdot, CRSS, gacumgr)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions
        double precision, intent(out) :: tau(nsysmx, 2), gamma1 ! updated CRSS (fwd/bwd) and new accumulated slip

        double precision :: deltgam, gamma0, tau0, tau1, thet0, thet1, fact, exp0, exp1, voce, dtau ! hardening intermediates
        integer :: jph, is, islast, mode, ism, js ! phase and slip/mode indices

        jph = voxel_loc%phase

        ! Calculate gamma increment
        deltgam = 0.0
        do is = 1, props%phase(jph)%nsyst
            deltgam = deltgam + abs(voxel_loc%gamdot(is)) * props%process(ipr)%tdot
        end do
        gamma0 = voxel_loc%gacumgr
        gamma1 = gamma0 + deltgam

        is = 0 ! total slip system counter
        islast = 0
        tau = 0.0
        do mode = 1, props%phase(jph)%nmodes ! loop over slip modes

            islast = islast + props%phase(jph)%nsm(mode) ! get last slip system in the new mode
            tau0 = props%phase(jph)%tau(islast, 1)
            tau1 = props%phase(jph)%tau(islast, 3)
            thet0 = props%phase(jph)%thet(islast, 0)
            thet1 = props%phase(jph)%thet(islast, 1)

            fact = abs(thet0 / tau1)
            exp0 = exp(-fact * gamma0)
            exp1 = exp(-fact * gamma1)
            voce = deltgam * thet1 + tau1 * (exp0 - exp1) + thet1 * (gamma0 * exp0 - gamma1 * exp1)

            do ism = 1, props%phase(jph)%nsm(mode) ! loop over slip system within mode
                is = is + 1 ! total slip system counter

                dtau = 0.
                do js = 1, props%phase(jph)%nsyst
                    dtau = dtau + props%phase(jph)%hard(is, js) * abs(voxel_loc%gamdot(js))
                end do

                tau(is, 1) = voxel_loc%crss(is, 1) + dtau * voce * props%process(ipr)%tdot / deltgam
                tau(is, 2) = voxel_loc%crss(is, 2) + dtau * voce * props%process(ipr)%tdot / deltgam
            end do
        end do

        return
    end

    subroutine update_hardening(micro, props)
        use types
        implicit none

        type(micro_type), intent(inout) :: micro
        type(props_type), intent(in) :: props

        double precision :: tautmp(nsysmx, 2)
        double precision :: gamacum
        integer :: ip1, ip2, ip3

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, tautmp, gamacum) &
        !$OMP SHARED(micro, props)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1
            if (.not. micro%voxel(ip1, ip2, ip3)%gas) then
                call harden(micro%voxel(ip1, ip2, ip3), props, tautmp, gamacum)
                micro%voxel(ip1, ip2, ip3)%crss(:, :) = tautmp
                micro%voxel(ip1, ip2, ip3)%trialtau(:, :) = tautmp
                micro%voxel(ip1, ip2, ip3)%gacumgr = gamacum
            end if
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine update_hardening

    ! Enforce mixed macroscopic boundary conditions by computing the required avg velgrad correction.
    subroutine calc_velgradavg_corr(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (sgavg, dvelgradavg updated)
        type(props_type), intent(inout) :: props ! simulation inputs/material/process definitions

        double precision :: aux33(3, 3) ! stress BC mismatch tensor (3x3)
        double precision :: dbar(6) ! avg symmetric velgrad in 6-component basis
        double precision :: aux6(6) ! scratch 6-vector
        double precision :: aux66(6, 6) ! scratch 6x6 matrix
        double precision :: aux3333(3, 3, 3, 3) ! scratch 4th-order tensor
        double precision :: sbar(6) ! avg stress in 6-component basis
        double precision :: mevptg(6, 6) ! macro compliance-like operator in 6-basis
        double precision :: devp0(6) ! deviatoric strain-rate offset for mixed BCs
        double precision :: velgradavg_sym(3, 3) ! symmetric part of avg velgrad
        double precision :: velgradavg_asym(3, 3) ! antisymmetric part of avg velgrad
        double precision :: wasim(3, 3) ! target macroscopic spin (antisymmetric part)
        double precision :: dvelgradavg_sym(3, 3) ! correction to symmetric avg velgrad
        double precision :: dvelgradavg_asym(3, 3) ! correction to antisymmetric avg velgrad
        double precision :: sgavgsol(3, 3) ! average solid stress (3x3)
        double precision :: sgavg(3, 3) ! average overall stress (3x3)
        integer :: ijv(6, 2) ! Voigt->tensor index map: (11,22,33,23,13,12)
        integer :: i ! loop index
        integer :: ii, jj ! tensor index loops
        integer :: m, n ! implied-do indices for ijv DATA statement
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)
        data((ijv(n, m), m=1, 2), n=1, 6) / 1, 1, 2, 2, 3, 3, 2, 3, 1, 3, 1, 2/ ! (11,22,33,23,13,12)

        sgavgsol = 0.0
        sgavg = 0.0
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3) &
        !$OMP SHARED(props, micro) &
        !$OMP REDUCTION(+:sgavgsol, sgavg)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1
            if (.not. micro%voxel(ip1, ip2, ip3)%gas) then
                sgavgsol = sgavgsol + micro%voxel(ip1, ip2, ip3)%wgtc * micro%voxel(ip1, ip2, ip3)%sg
                sgavg = sgavg + micro%voxel(ip1, ip2, ip3)%wgtc * micro%voxel(ip1, ip2, ip3)%sg
            end if
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        sgavgsol = sgavgsol / micro%wphcsol
        micro%sgavg = sgavg
        micro%sgavgsol = sgavgsol
        micro%svmavg = vm_stress(micro%sgavg)

        ! error in stress b.c.
        do i = 1, 6
            ii = ijv(i, 1)
            jj = ijv(i, 2)

            aux33(ii, jj) = props%process(ipr)%iscau(i) * (props%process(ipr)%scauchy(ii, jj) - &
                                                           micro%sgavg(ii, jj))
            aux33(jj, ii) = props%process(ipr)%iscau(i) * (props%process(ipr)%scauchy(jj, ii) - &
                                                           micro%sgavg(jj, ii))
        end do
        errsbc = norm2(aux33)

        ! mixed boundary conditions
        velgradavg_sym = sym33(micro%velgradavg)
        velgradavg_asym = antisym33(micro%velgradavg)
        call chg_basis(dbar, velgradavg_sym, aux66, aux3333, 2, 6)
        call chg_basis(sbar, micro%sgavg, aux66, aux3333, 2, 6)

        call chg_basis(aux6, aux33, mevptg, micro%de_ds, 4, 6)
        devp0 = dbar - matmul(mevptg, sbar)

        call state_6x6_evpsc(props%process(ipr)%idsim, dbar, props%process(ipr)%iscau, sbar, mevptg, devp0, &
                             props%process(ipr)%scauchy, props%process(ipr)%dsim)

        ! correction in symmetric part of average velgrad
        dvelgradavg_sym = props%process(ipr)%dsim - velgradavg_sym

        ! anti-symmetric part of average velgrad
        do ii = 1, 3
        do jj = 1, 3
            if (ii == jj) then
                wasim(ii, jj) = 0.0
            else
                if (props%process(ipr)%iudot(ii, jj) == 1 .and. props%process(ipr)%iudot(jj, ii) == 0) then
                    wasim(ii, jj) = props%process(ipr)%udot(ii, jj) - props%process(ipr)%dsim(ii, jj)
                    wasim(jj, ii) = -wasim(ii, jj)
                elseif (props%process(ipr)%iudot(ii, jj) == 0 .and. props%process(ipr)%iudot(jj, ii) == 1) then
                    wasim(ii, jj) = wasim(ii, jj) ! already assigned
                elseif (props%process(ipr)%iudot(ii, jj) == 1 .and. props%process(ipr)%iudot(jj, ii) == 1) then
                    wasim(ii, jj) = 0.5 * (props%process(ipr)%udot(ii, jj) - props%process(ipr)%udot(jj, ii))
                else
                    wasim(jj, ii) = 0.0
                end if
            end if
        end do
        end do

        ! correction in anti-symmetric part of average velgrad
        dvelgradavg_asym = wasim - velgradavg_asym

        ! total and accumulated corrections average velgrad
        micro%dvelgradavg = dvelgradavg_sym + dvelgradavg_asym

        return

    end

    ! Solve the mixed stress/strain-rate boundary condition system in 6x6 (Voigt) form.
    subroutine state_6x6_evpsc(idsim, dbar, iscau, sbar, xmevptg, devpzero, scauchy, dsim)
        use tensor_functions
        implicit none

        integer, intent(in) :: idsim(6) ! mask for imposed strain-rate components (Voigt)
        integer, intent(in) :: iscau(6) ! mask for imposed stress components (Voigt)
        double precision, intent(in) :: xmevptg(6, 6) ! operator relating strain-rate and stress increments
        double precision, intent(in) :: devpzero(6) ! offset in strain-rate space for mixed BCs
        double precision, intent(inout) :: dsim(3, 3) ! symmetric strain-rate tensor (updated)
        double precision, intent(inout) :: dbar(6) ! strain-rate in 6-component basis (updated)
        double precision, intent(inout) :: scauchy(3, 3) ! macroscopic Cauchy stress tensor (updated)
        double precision, intent(inout) :: sbar(6) ! stress in 6-component basis (updated)
        double precision :: aux6(6) ! scratch 6-vector
        double precision :: aux66(6, 6) ! scratch 6x6 matrix
        double precision :: aux6a(6) ! scratch RHS vector for linear system
        double precision :: profac(6) ! Voigt shear scaling factors (engineering vs tensorial)
        double precision :: aux3333(3, 3, 3, 3) ! scratch 4th-order tensor
        double precision :: aux33(3, 3) ! scratch 3x3 tensor
        double precision :: sbarv(6) ! stress in VPSC Voigt convention
        double precision :: dbarv(6) ! strain-rate in VPSC Voigt convention
        double precision :: xmevptgv(6, 6) ! operator in VPSC Voigt convention
        double precision :: devpzerov(6) ! offset in VPSC Voigt convention
        double precision :: aux6b(6) ! scratch copy of RHS (unused)

        integer :: i, j ! loop indices

        do i = 1, 6
            profac(i) = 1.0 + (i / 4)
        end do

        call voigt_vpsc(dbarv, dsim, aux66, aux3333, 2)
        call voigt_vpsc(sbarv, scauchy, aux66, aux3333, 2)
        call chg_basis(aux6, aux33, xmevptg, aux3333, 3, 6)
        call voigt_vpsc(aux6, aux33, xmevptgv, aux3333, 4)
        call chg_basis(devpzero, aux33, aux66, aux3333, 1, 6)
        call voigt_vpsc(devpzerov, aux33, aux66, aux3333, 2)

        do i = 1, 6
            aux6(i) = dbarv(i) - devpzerov(i)
        end do

        do i = 1, 6
            aux6a(i) = -1.d0 * idsim(i) * aux6(i)
            do j = 1, 6
                aux6a(i) = aux6a(i) + xmevptgv(i, j) * iscau(j) * sbarv(j) * profac(j)
                aux66(i, j) = iscau(j) * (i / j) * (j / i) - idsim(j) * xmevptgv(i, j) * profac(j)
            end do
        end do

        aux6b = aux6a
        call lu_eqsystem(aux66, aux6a, 6)

        do i = 1, 6
            aux6(i) = idsim(i) * aux6(i) + iscau(i) * aux6a(i)
            sbarv(i) = iscau(i) * sbarv(i) + idsim(i) * aux6a(i)
        end do

        do i = 1, 6
            dbarv(i) = aux6(i) + devpzerov(i)
        end do

        call voigt_vpsc(dbarv, dsim, aux66, aux3333, 1)
        call chg_basis(dbar, dsim, aux66, aux3333, 2, 6)
        call voigt_vpsc(sbarv, scauchy, aux66, aux3333, 1)
        call chg_basis(sbar, scauchy, aux66, aux3333, 2, 6)

        return
    end

    ! Update the plastic deformation gradient Fp using exp(Lp*dt) per voxel.
    ! subroutine update_plastic_defgrad(micro, props)
    !     use tensor_functions
    !     use types
    !     implicit none

    !     type(micro_type), intent(inout) :: micro ! microstructure state (Fp updated)
    !     type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

    !     double precision :: aa(3, 3) ! scratch 3x3 matrix (unused)
    !     double precision :: Lpinter(3, 3) ! plastic velocity gradient in intermediate configuration
    !     double precision :: disgradincp(3, 3) ! incremental plastic displacement gradient (Lp*dt)
    !     double precision :: expdisgradincp(3, 3) ! exp(disgradincp) update for Fp
    !     integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)
    !     integer :: iph ! phase index (unused)

    !     do ip3 = 1, props%npts3
    !     do ip2 = 1, props%npts2
    !     do ip1 = 1, props%npts1

    !         if (.not. micro%voxel(ip1, ip2, ip3)%gas) then

    !             ! plastic velocity gradient in intermediate configuration
    !             call plastic_velgrad(micro%voxel(ip1, ip2, ip3), props, Lpinter)
    !             micro%voxel(ip1, ip2, ip3)%Lpinter(:, :) = Lpinter

    !             ! update of plastic deformation gradient
    !             disgradincp = Lpinter * props%process(ipr)%tdot
    !             expdisgradincp = matrix_exp_adaptive(disgradincp)

    !             micro%voxel(ip1, ip2, ip3)%defgradp(:, :) = matmul(expdisgradincp, &
    !                                                         micro%voxel(ip1, ip2, ip3)%defgradp(:, :))

    !         end if

    !     end do
    !     end do
    !     end do

    ! end

    subroutine update_plastic_defgrad(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (Fp updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: Lpinter(3, 3) ! plastic velocity gradient in intermediate configuration
        double precision :: disgradincp(3, 3) ! incremental plastic displacement gradient (Lp*dt)
        double precision :: expdisgradincp(3, 3) ! exp(disgradincp) update for Fp
        double precision :: tdot_loc ! cached timestep for thread-safe OpenMP region
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)

        tdot_loc = props%process(ipr)%tdot

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, Lpinter, disgradincp, expdisgradincp) &
        !$OMP SHARED(micro, props, tdot_loc)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            if (.not. micro%voxel(ip1, ip2, ip3)%gas) then

                ! plastic velocity gradient in intermediate configuration
                call plastic_velgrad(micro%voxel(ip1, ip2, ip3), props, Lpinter)
                micro%voxel(ip1, ip2, ip3)%Lpinter(:, :) = Lpinter

                ! update of plastic deformation gradient
                disgradincp = Lpinter * tdot_loc
                expdisgradincp = matrix_exp_adaptive(disgradincp)

                micro%voxel(ip1, ip2, ip3)%defgradp(:, :) = matmul(expdisgradincp, &
                                                                   micro%voxel(ip1, ip2, ip3)%defgradp(:, :))

            end if

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end
    ! Compute the plastic velocity gradient Lp from slip rates and system geometry.
    subroutine plastic_velgrad(voxel, props, Lp)
        use types
        implicit none

        type(voxel_type), intent(in) :: voxel ! voxel state (orientation, slip rates)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions
        double precision, intent(out) :: Lp(3, 3) ! plastic velocity gradient (sample axes)

        double precision :: aa(3, 3) ! crystal->sample rotation matrix
        double precision :: dnsa(3) ! slip plane normal in sample axes
        double precision :: dbsa(3) ! slip direction in sample axes
        integer :: i, j ! tensor index loops
        integer :: iph ! phase index
        integer :: isys ! slip system index

        aa = voxel%ag(:, :)
        iph = voxel%phase

        Lp = 0.0
        do isys = 1, props%phase(iph)%nsyst

            dnsa = matmul(aa, props%phase(iph)%dnca(:, isys))
            dbsa = matmul(aa, props%phase(iph)%dbca(:, isys))

            do i = 1, 3
                do j = 1, 3
                    Lp(i, j) = Lp(i, j) + dbsa(i) * dnsa(j) * voxel%gamdot(isys)
                end do
            end do

        end do

    end

    ! Approximate exp(mat) for a 3x3 matrix using a truncated Taylor expansion (up to cubic term).
    function matrix_exp(mat) result(mat_exp)
        use global, only: id3
        implicit none
        double precision, intent(in), dimension(3, 3) :: mat ! input 3x3 matrix
        double precision, dimension(3, 3) :: mat_exp ! approximate exp(mat)
        double precision, dimension(3, 3) :: mat2 ! mat^2
        double precision, dimension(3, 3) :: mat3 ! mat^3

        mat2 = matmul(mat, mat)
        mat3 = matmul(mat2, mat)
        mat_exp = id3 + mat + mat2 / 2.0 + mat3 / 6.0

    end function

    pure double precision function l1_norm33(mat) result(l1)
        implicit none
        double precision, intent(in) :: mat(3, 3)
        double precision :: colsum
        integer :: j

        l1 = 0.0d0
        do j = 1, 3
            colsum = abs(mat(1, j)) + abs(mat(2, j)) + abs(mat(3, j))
            if (colsum > l1) l1 = colsum
        end do

    end function l1_norm33

    function matrix_exp_pade(mat) result(mat_exp)
        ! 13th-order Padé with scaling & squaring for a 3x3 matrix.
        !
        ! Implements the classic Higham 2005 style [13/13] Padé:
        !   E ≈ (V - U)^{-1} (V + U)
        ! with scaling: A <- A / 2^s, then E <- E^(2^s) by repeated squaring.

        use global, only: id3
        ! use tensor_functions, only : lu_inverse
        implicit none
        double precision, intent(in), dimension(3, 3) :: mat ! input 3x3 matrix
        double precision, dimension(3, 3) :: mat_exp ! approximate exp(mat)
        double precision, dimension(3, 3) :: A, A2, A4, A6
        double precision, dimension(3, 3) :: U, V, P, Q
        double precision :: normA, theta13
        double precision :: b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13
        integer :: s, k
        double precision :: R(3, 3), X(3, 3)
        integer :: ipiv(3), info

        theta13 = 5.371920351148152D0
        b0 = 64764752532480000.D0
        b1 = 32382376266240000.D0
        b2 = 7771770303897600.D0
        b3 = 1187353796428800.D0
        b4 = 129060195264000.D0
        b5 = 10559470521600.D0
        b6 = 670442572800.D0
        b7 = 33522128640.D0
        b8 = 1323241920.D0
        b9 = 40840800.D0
        b10 = 960960.D0
        b11 = 16380.D0
        b12 = 182.D0
        b13 = 1.D0

        A = mat
        normA = l1_norm33(A)

        if (normA <= 0.0d0) then
            mat_exp = id3
            return
        end if

        ! s = max(0, ceil(log2(normA/theta13)))
        if (normA <= theta13) then
            s = 0
        else
            s = ceiling(log(normA / theta13) / log(2.0D0))
        end if
        A = A / (2.0D0**s)

        ! Powers
        A2 = matmul(A, A)
        A4 = matmul(A2, A2)
        A6 = matmul(A2, A4)

        ! Padé(13) coefficients from Higham (2005), arranged per standard expm implementations:
        !   U = A * (A6*(b13*A6 + b11*A4 + b9*A2) + b7*A6 + b5*A4 + b3*A2 + b1*I)
        !   V =      A6*(b12*A6 + b10*A4 + b8*A2) + b6*A6 + b4*A4 + b2*A2 + b0*I
        U = matmul(A, (matmul(A6, (b13 * A6 + b11 * A4 + b9 * A2)) + b7 * A6 + b5 * A4 + b3 * A2 + b1 * id3))
        V = matmul(A6, (b12 * A6 + b10 * A4 + b8 * A2)) + b6 * A6 + b4 * A4 + b2 * A2 + b0 * id3

        ! Compute (V - U)^{-1} (V + U)
        P = V + U
        Q = V - U
        R = Q

        X = P
        call dgetrf(3, 3, R, 3, ipiv, info)
        ! if (info /= 0) stop 'Singular (V-U)'

        if (info /= 0) then
            write (*, '(A,I0)') 'Singular (V-U), info=', info
            stop 1
        end if

        call dgetrs('N', 3, 3, R, 3, ipiv, X, 3, info)
        ! if (info /= 0) stop 'Solve failed'
        if (info /= 0) then
            write (*, '(A,I0)') 'Solve failed info=', info
            stop 1
        end if

        mat_exp = X

        ! R = Q
        ! call lu_inverse(R, 3)
        ! mat_exp = matmul(R, P)
        ! Squaring
        do k = 1, s
            mat_exp = matmul(mat_exp, mat_exp)

        end do

    end function

    function matrix_exp_adaptive(mat) result(mat_exp)
        implicit none
        double precision, intent(in), dimension(3, 3) :: mat ! input 3x3 matrix
        double precision, dimension(3, 3) :: mat_exp ! selected exp(mat) approximation
        double precision, parameter :: taylor_l1_threshold = 1.0d-1

        if (l1_norm33(mat) <= taylor_l1_threshold) then
            mat_exp = matrix_exp(mat)
        else
            mat_exp = matrix_exp_pade(mat)
        end if

    end function

    ! Update nodal coordinates by integrating the reference-config velgrad field to nodal displacements.
    subroutine update_grid_velgrad_node(micro, props)
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (xnode updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: velgradref(3, 3, props%npts1, props%npts2, props%npts3) ! velgrad mapped to reference config
        double precision :: velgradrefavg(3, 3) ! average velgradref
        double precision :: dv(3, props%npts1, props%npts2, props%npts3) ! integrated displacement fluctuation at cells
        double precision :: x(3) ! reference-space node coordinate
        double precision :: wgt(props%npts1, props%npts2, props%npts3) ! weights for averaging (uniform)
        double precision :: dv_node(3, props%npts1 + 1, props%npts2 + 1, props%npts3 + 1) ! nodal displacement fluctuation
        integer :: ideriv ! derivative operator selector
        integer :: ip1, ip2, ip3 ! voxel/node indices
        integer :: npts3node ! number of nodes along z (1 or npts3+1)
        double precision :: tdot_loc ! cached timestep for OpenMP-safe update

        tdot_loc = props%process(ipr)%tdot

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3) &
        !$OMP SHARED(micro, props, velgradref)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1
            velgradref(:, :, ip1, ip2, ip3) = matmul(micro%voxel(ip1, ip2, ip3)%velgrad(:, :), &
                                                     micro%voxel(ip1, ip2, ip3)%defgrad(:, :))
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        wgt = props%wgt

        call average_tensor2(velgradref, wgt, velgradrefavg, 3, 3, props%npts1, props%npts2, props%npts3)

        ideriv = props%process(ipr)%igamma

        if (ideriv /= 2) then

            call integrate_tensor2(velgradref, dv, ideriv, props%npts1, props%npts2, props%npts3)

            call cell2node_tensor1(dv, dv_node, 3, props%npts1, props%npts2, props%npts3)

        else

            call integrate_tensor2(velgradref, dv, ideriv, props%npts1, props%npts2, props%npts3)
            dv_node(:, 1:props%npts1, 1:props%npts2, 1:props%npts3) = dv
            dv_node(:, props%npts1 + 1, :, :) = dv_node(:, 1, :, :)
            dv_node(:, :, props%npts2 + 1, :) = dv_node(:, :, 1, :)
            dv_node(:, :, :, props%npts3 + 1) = dv_node(:, :, :, 1)

        end if

        if (props%npts3 == 1) then
            npts3node = 1
        else
            npts3node = props%npts3 + 1
        end if
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, x) &
        !$OMP SHARED(micro, props, npts3node, velgradrefavg, dv_node, tdot_loc)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, npts3node
        do ip2 = 1, props%npts2 + 1
        do ip1 = 1, props%npts1 + 1

            x(1) = float(ip1) - 0.5
            x(2) = float(ip2) - 0.5
            if (npts3node == 1) then
                x(3) = 0.0
            else
                x(3) = float(ip3) - 0.5
            end if

            micro%xnode(:, ip1, ip2, ip3) = micro%xnode(:, ip1, ip2, ip3) + &
                                            (matmul(velgradrefavg, x) + dv_node(:, ip1, ip2, ip3)) * tdot_loc

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end

    ! ! Update deformation gradient F, inverse F^{-1}, detF, and current voxel weights for the new step.
    ! subroutine update_defgrad(micro, props)
    !     use tensor_functions
    !     use types

    !     implicit none

    !     type(micro_type), intent(inout) :: micro ! microstructure state (F, F^{-1}, detF, weights updated)
    !     type(props_type), intent(in) :: props ! simulation inputs/material/process definitions
    !     double precision :: defgradold(3, 3) ! previous voxel deformation gradient
    !     double precision :: defgradinvold(3, 3) ! previous voxel inverse deformation gradient (unused)
    !     double precision :: disgradinc(3, 3) ! incremental displacement gradient (L*dt)
    !     double precision :: expdisgradinc(3, 3) ! exp(disgradinc) update for F
    !     double precision :: defgradinv(3, 3) ! scratch: inverse deformation gradient
    !     double precision :: detFavg ! average det(F) for normalization
    !     integer :: i, j, k ! loop indices (unused)
    !     integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)
    !     ! double precision :: Dinc(3,3) = 0
    !     ! double precision :: Winc(3,3) = 0
    !     ! double precision :: devDinc(3,3) = 0

    !     micro%defgradavg = 0.0
    !     micro%defgradinvavgc = 0.0
    !     micro%defgradinvavgcs = 0.0
    !     micro%defgradinvavgcg = 0.0
    !     micro%wphcsol = 0.0
    !     micro%wphcgas = 0.0
    !     detfavg = 0.0
    !     do ip3 = 1, props%npts3
    !     do ip2 = 1, props%npts2
    !     do ip1 = 1, props%npts1

    !         defgradold = micro%voxel(ip1, ip2, ip3)%defgrad(:, :)
    !         defgradinvold = micro%voxel(ip1, ip2, ip3)%defgradinv(:, :)

    !         disgradinc = micro%voxel(ip1, ip2, ip3)%velgrad(:, :) * props%process(ipr)%tdot
    !         ! vm_increment: micro%edotvmmx * props%process(ipr)%tdot
    !         ! Dinc = sym33(disgradinc) ! strain/stretch increment
    !         ! Winc = skew33(disgradin) ! rotation increment
    !         ! devDinc = Dinc - (tr33(Dinc)/3) * id3

    !         expdisgradinc = matrix_exp_pade(disgradinc)
    !         micro%voxel(ip1, ip2, ip3)%defgrad(:, :) = matmul(expdisgradinc, defgradold)
    !         micro%voxel(ip1, ip2, ip3)%detF = determinant33(micro%voxel(ip1, ip2, ip3)%defgrad(:, :))

    !         if (micro%voxel(ip1, ip2, ip3)%detF < 0) then
    !             write (*, '(A,E13.6,A,I0,A,I0,A,I0,A)') &
    !                 ' -> WARNING: detF=', micro%voxel(ip1, ip2, ip3)%detF, &
    !                 ' in voxel (', ip1, ',', ip2, ',', ip3, ')'
    !         endif

    !         defgradinv = micro%voxel(ip1, ip2, ip3)%defgrad
    !         call lu_inverse(defgradinv, 3)
    !         micro%voxel(ip1, ip2, ip3)%defgradinv(:, :) = defgradinv

    !         micro%defgradavg = micro%defgradavg + micro%voxel(ip1, ip2, ip3)%defgrad * props%wgt
    !         micro%defgradinvavgc = micro%defgradinvavgc + micro%voxel(ip1, ip2, ip3)%defgradinv * &
    !                                props%wgt * micro%voxel(ip1, ip2, ip3)%detF
    !         if (micro%voxel(ip1, ip2, ip3)%gas) then
    !             micro%defgradinvavgcg = micro%defgradinvavgcg + micro%voxel(ip1, ip2, ip3)%defgradinv * &
    !                                     props%wgt * micro%voxel(ip1, ip2, ip3)%detF
    !             micro%wphcgas = micro%wphcgas + props%wgt * micro%voxel(ip1, ip2, ip3)%detF
    !         else
    !             micro%defgradinvavgcs = micro%defgradinvavgcs + micro%voxel(ip1, ip2, ip3)%defgradinv * &
    !                                     props%wgt * micro%voxel(ip1, ip2, ip3)%detF
    !             micro%wphcsol = micro%wphcsol + props%wgt * micro%voxel(ip1, ip2, ip3)%detF
    !         end if
    !         micro%voxel(ip1, ip2, ip3)%defgradinc(:, :) = expdisgradinc
    !         micro%voxel(ip1, ip2, ip3)%wgtc = props%wgt * micro%voxel(ip1, ip2, ip3)%detF

    !         detFavg = detFavg + micro%voxel(ip1, ip2, ip3)%detF * props%wgt
    !     end do
    !     end do
    !     end do

    !     micro%defgradinvavgc = micro%defgradinvavgc / detFavg
    !     micro%wphcsol = micro%wphcsol / detFavg
    !     micro%wphcgas = micro%wphcgas / detFavg
    !     micro%defgradinvavgcs = micro%defgradinvavgcs / detFavg / micro%wphcsol
    !     if (micro%wphcgas > 0.0) then
    !         micro%defgradinvavgcg = micro%defgradinvavgcg / detFavg / micro%wphcgas
    !     else
    !         micro%defgradinvavgcg = id3
    !     end if
    !     do ip3 = 1, props%npts3
    !     do ip2 = 1, props%npts2
    !     do ip1 = 1, props%npts1
    !         micro%voxel(ip1, ip2, ip3)%wgtc = micro%voxel(ip1, ip2, ip3)%wgtc / detFavg
    !     end do
    !     end do
    !     end do

    ! end

    subroutine update_defgrad(micro, props)
        use global, only: id3
        use tensor_functions
        use types

        implicit none

        type(micro_type), intent(inout) :: micro
        type(props_type), intent(in) :: props

        double precision :: defgradold(3, 3)
        double precision :: disgradinc(3, 3)
        double precision :: expdisgradinc(3, 3)
        double precision :: defgradinv(3, 3)
        double precision :: detFavg
        double precision :: tdot_loc

        double precision :: defgradavg_acc(3, 3)
        double precision :: defgradinvavgc_acc(3, 3)
        double precision :: defgradinvavgcs_acc(3, 3)
        double precision :: defgradinvavgcg_acc(3, 3)
        double precision :: wphcsol_acc
        double precision :: wphcgas_acc
        double precision :: detFavg_acc

        integer :: ip1, ip2, ip3

        tdot_loc = props%process(ipr)%tdot

        micro%defgradavg = 0.0d0
        micro%defgradinvavgc = 0.0d0
        micro%defgradinvavgcs = 0.0d0
        micro%defgradinvavgcg = 0.0d0
        micro%wphcsol = 0.0d0
        micro%wphcgas = 0.0d0
        detFavg = 0.0d0

        defgradavg_acc = 0.0d0
        defgradinvavgc_acc = 0.0d0
        defgradinvavgcs_acc = 0.0d0
        defgradinvavgcg_acc = 0.0d0
        wphcsol_acc = 0.0d0
        wphcgas_acc = 0.0d0
        detFavg_acc = 0.0d0

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, defgradold, disgradinc, expdisgradinc, defgradinv) &
        !$OMP SHARED(micro, props, tdot_loc) &
        !$OMP REDUCTION(+:defgradavg_acc, defgradinvavgc_acc, defgradinvavgcs_acc, defgradinvavgcg_acc, &
        !$OMP              wphcsol_acc, wphcgas_acc, detFavg_acc)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            defgradold = micro%voxel(ip1, ip2, ip3)%defgrad(:, :)

            disgradinc = micro%voxel(ip1, ip2, ip3)%velgrad(:, :) * tdot_loc
            expdisgradinc = matrix_exp_adaptive(disgradinc)

            micro%voxel(ip1, ip2, ip3)%defgrad(:, :) = matmul(expdisgradinc, defgradold)
            micro%voxel(ip1, ip2, ip3)%detF = determinant33(micro%voxel(ip1, ip2, ip3)%defgrad(:, :))

            if (micro%voxel(ip1, ip2, ip3)%detF < 0.0d0) then
                write (*, '(A,E13.6,A,I0,A,I0,A,I0,A)') &
                    ' -> WARNING: detF=', micro%voxel(ip1, ip2, ip3)%detF, &
                    ' in voxel (', ip1, ',', ip2, ',', ip3, ')'
            end if

            defgradinv = micro%voxel(ip1, ip2, ip3)%defgrad
            call lu_inverse(defgradinv, 3)
            micro%voxel(ip1, ip2, ip3)%defgradinv(:, :) = defgradinv

            defgradavg_acc = defgradavg_acc + micro%voxel(ip1, ip2, ip3)%defgrad * props%wgt
            defgradinvavgc_acc = defgradinvavgc_acc + micro%voxel(ip1, ip2, ip3)%defgradinv * &
                                 props%wgt * micro%voxel(ip1, ip2, ip3)%detF

            if (micro%voxel(ip1, ip2, ip3)%gas) then
                defgradinvavgcg_acc = defgradinvavgcg_acc + micro%voxel(ip1, ip2, ip3)%defgradinv * &
                                      props%wgt * micro%voxel(ip1, ip2, ip3)%detF
                wphcgas_acc = wphcgas_acc + props%wgt * micro%voxel(ip1, ip2, ip3)%detF
            else
                defgradinvavgcs_acc = defgradinvavgcs_acc + micro%voxel(ip1, ip2, ip3)%defgradinv * &
                                      props%wgt * micro%voxel(ip1, ip2, ip3)%detF
                wphcsol_acc = wphcsol_acc + props%wgt * micro%voxel(ip1, ip2, ip3)%detF
            end if

            micro%voxel(ip1, ip2, ip3)%defgradinc(:, :) = expdisgradinc
            micro%voxel(ip1, ip2, ip3)%wgtc = props%wgt * micro%voxel(ip1, ip2, ip3)%detF

            detFavg_acc = detFavg_acc + micro%voxel(ip1, ip2, ip3)%detF * props%wgt

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        micro%defgradavg = defgradavg_acc
        micro%defgradinvavgc = defgradinvavgc_acc
        micro%defgradinvavgcs = defgradinvavgcs_acc
        micro%defgradinvavgcg = defgradinvavgcg_acc
        micro%wphcsol = wphcsol_acc
        micro%wphcgas = wphcgas_acc
        detFavg = detFavg_acc

        micro%defgradinvavgc = micro%defgradinvavgc / detFavg
        micro%wphcsol = micro%wphcsol / detFavg
        micro%wphcgas = micro%wphcgas / detFavg
        micro%defgradinvavgcs = micro%defgradinvavgcs / detFavg / micro%wphcsol
        if (micro%wphcgas > 0.0d0) then
            micro%defgradinvavgcg = micro%defgradinvavgcg / detFavg / micro%wphcgas
        else
            micro%defgradinvavgcg = id3
        end if

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3) &
        !$OMP SHARED(micro, props, detFavg)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1
            micro%voxel(ip1, ip2, ip3)%wgtc = micro%voxel(ip1, ip2, ip3)%wgtc / detFavg
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine update_defgrad

    ! Update elastic deformation gradient Fe from total F, plastic Fp, and initial F0.
    subroutine update_elstic_defgrad(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (Fe updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions
        double precision :: defgradeinvold(3, 3) ! inverse of previous elastic defgrad (Fe^{-1})
        double precision :: FpFiniinv(3, 3) ! inverse of (Fp * Fini)
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, defgradeinvold, FpFiniinv) &
        !$OMP SHARED(micro, props)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            if (.not. micro%voxel(ip1, ip2, ip3)%gas) then

                defgradeinvold = micro%voxel(ip1, ip2, ip3)%defgrade(:, :)
                call lu_inverse(defgradeinvold, 3)

                ! elastic deformation gradient from total, plastic and initial
                FpFiniinv = matmul(micro%voxel(ip1, ip2, ip3)%defgradp(:, :), micro%voxel(ip1, ip2, ip3)%defgradini(:, :))
                call lu_inverse(FpFiniinv, 3)

                micro%voxel(ip1, ip2, ip3)%defgrade(:, :) = matmul(micro%voxel(ip1, ip2, ip3)%defgrad(:, :), FpFiniinv)
                micro%voxel(ip1, ip2, ip3)%defgradeinc(:, :) = matmul(micro%voxel(ip1, ip2, ip3)%defgrade(:, :), defgradeinvold)

            end if

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end

    ! Push-forward crystal stiffness to current configuration and update voxel stiffness/compliance matrices.
    subroutine update_el_stiff(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (cg66/sg66 updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions
        double precision :: aa(3, 3) ! combined mapping (Fe * ag)
        double precision :: aat(3, 3) ! transpose of aa
        double precision :: c(3, 3, 3, 3) ! pushed-forward crystal stiffness tensor
        double precision :: dum ! scratch scalar accumulator
        double precision :: detFe ! det(Fe) Jacobian
        double precision :: c66(6, 6) ! stiffness in 6x6 basis (then inverted to compliance)
        double precision :: aux6(6) ! scratch 6-vector
        double precision :: aux33(3, 3) ! scratch 3x3 tensor
        integer :: i1, j1, k1, l1 ! tensor index loops (1..3)
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)
        integer :: jph ! phase index
        integer :: i2, j2, k2, l2 ! tensor contraction indices (1..3)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, jph, aa, aat, detFe, i1, j1, k1, l1, i2, j2, k2, l2, dum) &
        !$OMP PRIVATE(c, aux6, aux33, c66) &
        !$OMP SHARED(props, micro)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            if (.not. micro%voxel(ip1, ip2, ip3)%gas) then

                jph = micro%voxel(ip1, ip2, ip3)%phase

                ! elastic stiffness from initial to current configuration
                ! faster calculation
                aa = matmul(micro%voxel(ip1, ip2, ip3)%defgrade, micro%voxel(ip1, ip2, ip3)%ag)
                aat = aa
                aat = transpose(aat)
                detFe = determinant33(micro%voxel(ip1, ip2, ip3)%defgrade)

                do i1 = 1, 3
                do j1 = 1, 3
                do k1 = 1, 3
                do l1 = 1, 3
                    dum = 0.0
                    do i2 = 1, 3
                    do j2 = 1, 3
                    do k2 = 1, 3
                    do l2 = 1, 3
                        dum = dum + aa(i1, i2) * aa(j1, j2) * props%phase(jph)%c(i2, j2, k2, l2) * &
                              aat(k2, k1) * aat(l2, l1)
                    end do
                    end do
                    end do
                    end do
                    c(i1, j1, k1, l1) = dum / detFe
                end do
                end do
                end do
                end do

                ! Calculate 6x6 stiffness tensor and invert to get compliance tensor
                call chg_basis(aux6, aux33, c66, c, 4, 6)
                micro%voxel(ip1, ip2, ip3)%cg66 = c66
                call lu_inverse(c66, 6)
                micro%voxel(ip1, ip2, ip3)%sg66 = c66

            end if

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end

    ! Objective update: rotate stresses/strains and adjust velgrad to keep consistent averages after config update.
    subroutine update_tensor_config(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (objective updates applied)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: defgradinc(3, 3) ! incremental total defgrad update (exp(L*dt))
        double precision :: defgradeinc(3, 3) ! incremental elastic defgrad update (Fe * Fe_old^{-1})
        double precision :: detdefgradeinc ! det(defgradeinc) for push-forward of stress
        double precision :: drot(3, 3) ! rotation from polar decomposition
        double precision :: dstrech(3, 3) ! stretch from polar decomposition (unused)
        double precision :: defgradincinv(3, 3) ! inverse of defgradinc
        double precision :: velgradavg(3, 3) ! average velgrad after rotation update
        double precision :: velgradavg_acc(3, 3) ! reduced accumulator for avg velgrad
        double precision :: defgradinvavgcinv(3, 3) ! inverse of avg(F^{-1})
        double precision :: dudot_dXini_avg(3, 3) ! correction to maintain avg velgrad
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)

        velgradavg = 0.0d0
        velgradavg_acc = 0.0d0
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, defgradinc, defgradeinc, detdefgradeinc, drot, dstrech, defgradincinv) &
        !$OMP SHARED(micro, props) &
        !$OMP REDUCTION(+:velgradavg_acc)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            defgradinc = micro%voxel(ip1, ip2, ip3)%defgradinc
            defgradeinc = micro%voxel(ip1, ip2, ip3)%defgradeinc
            detdefgradeinc = determinant33(defgradeinc)
            defgradincinv = defgradinc
            call lu_inverse(defgradincinv, 3)

            if (.not. micro%voxel(ip1, ip2, ip3)%gas) then

                call polar_dcmp(defgradinc, drot, dstrech)
                ! call polar_dcmp(defgradeinctmp,drote,aux33)

                micro%voxel(ip1, ip2, ip3)%disgradsym = matmul(matmul(drot, &
                                                                      micro%voxel(ip1, ip2, ip3)%disgradsym), &
                                                               transpose(drot))
                micro%voxel(ip1, ip2, ip3)%ept = matmul(matmul(drot, micro%voxel(ip1, ip2, ip3)%ept), &
                                                        transpose(drot))
                micro%voxel(ip1, ip2, ip3)%sgt = matmul(matmul(defgradeinc, micro%voxel(ip1, ip2, ip3)%sgt), &
                                                        transpose(defgradeinc)) / detdefgradeinc

            end if

            micro%voxel(ip1, ip2, ip3)%velgrad = matmul(micro%voxel(ip1, ip2, ip3)%velgrad, defgradincinv)
            velgradavg_acc = velgradavg_acc + micro%voxel(ip1, ip2, ip3)%velgrad * micro%voxel(ip1, ip2, ip3)%wgtc

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        velgradavg = velgradavg_acc

        ! ensure same volume average
        defgradinvavgcinv = micro%defgradinvavgc
        call lu_inverse(defgradinvavgcinv, 3)

        dudot_dXini_avg = matmul(micro%velgradavg - velgradavg, defgradinvavgcinv)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3) &
        !$OMP SHARED(micro, props, dudot_dXini_avg)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            micro%voxel(ip1, ip2, ip3)%velgrad = micro%voxel(ip1, ip2, ip3)%velgrad + &
                                                 matmul(dudot_dXini_avg, &
                                                        micro%voxel(ip1, ip2, ip3)%defgradinv)

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end

    ! Record the indices of the top-3 most active slip systems per voxel (by |gamdot|).
    subroutine order_sys_per_activity(micro, props)
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (igamdotmx updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        double precision :: auxnsmx(nsysmx) ! scratch array of |gamdot| per system
        double precision :: mxval ! current maximum activity (unused)
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)
        integer :: i ! top-k counter (1..3)
        integer :: iph ! phase index
        integer :: isys ! slip system index
        integer :: mxind ! index of current max activity system

        ! ordering of system per activity
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1
            iph = micro%voxel(ip1, ip2, ip3)%phase
            if (.not. micro%voxel(ip1, ip2, ip3)%gas) then
                auxnsmx = 0.0
                do isys = 1, props%phase(iph)%nsyst
                    auxnsmx(isys) = abs(micro%voxel(ip1, ip2, ip3)%gamdot(isys))
                end do
                do i = 1, 3
                    mxval = maxval(auxnsmx)
                    mxind = maxloc(auxnsmx, 1)
                    auxnsmx(mxind) = -1.0
                    micro%voxel(ip1, ip2, ip3)%igamdotmx(i) = mxind
                end do
            else
                micro%voxel(ip1, ip2, ip3)%igamdotmx(:) = 1
            end if

        end do
        end do
        end do

    end

    ! Compute/update macroscopic tangent stiffness/compliance tensors from per-voxel consistent tangents.
    subroutine calc_c0(micro, props)
        use tensor_functions
        use types
        implicit none

        type(micro_type), intent(inout) :: micro ! microstructure state (macro tangents updated)
        type(props_type), intent(in) :: props ! simulation inputs/material/process definitions

        type(voxel_type) :: voxel_loc ! local voxel copy for tangent evaluation
        double precision :: dP6_dF6(6, 6) ! averaged PK1 tangent in 6x6 basis
        double precision :: ds_de(3, 3, 3, 3) ! averaged Cauchy tangent wrt strain (3x3x3x3)
        double precision :: detF ! voxel det(F)
        double precision :: Finv(3, 3) ! voxel F^{-1}
        double precision :: wgtc ! voxel current-config weight
        double precision :: sg66(6, 6) ! voxel compliance in 6x6 basis
        double precision :: sg(3, 3) ! voxel Cauchy stress (3x3)
        double precision :: sg6(6) ! voxel stress in 6-component basis
        double precision :: aux33(3, 3) ! scratch 3x3 tensor
        double precision :: aux6(6) ! scratch 6-vector
        double precision :: aux66(6, 6) ! scratch 6x6 matrix
        double precision :: aux3333(3, 3, 3, 3) ! scratch 4th-order tensor
        double precision :: edotp6(6) ! plastic strain-rate in 6-component basis
        double precision :: dedotp6_dsg6(6, 6) ! d(edotp6)/d(sg6)
        double precision :: dsg6_de6(6, 6) ! inverse of d(edot)/d(sg) (tangent)
        double precision :: dsg6_de6_tref(6, 6) ! tangent scaled to reference dt
        double precision :: dsg_de(3, 3, 3, 3) ! tangent in 3x3x3x3 form
        double precision :: dsg_de_tref(3, 3, 3, 3) ! tangent scaled to reference dt (used for PK1)
        double precision :: dPg6_dF6(6, 6) ! voxel PK1 tangent in 6x6 basis
        double precision :: dPg_dF(3, 3, 3, 3) ! voxel PK1 tangent in 4th-order tensor form
        double precision :: sgt(3, 3) ! voxel previous-step stress (3x3)
        double precision :: sgt6(6) ! voxel previous-step stress (6-basis)
        double precision :: edote6(6) ! elastic strain-rate in 6-basis
        double precision :: dedote6_dsg6(6, 6) ! d(edote6)/d(sg6)
        integer :: ip1, ip2, ip3 ! voxel indices (x,y,z)
        integer :: i, j, k, l, m, n ! tensor index loops

        dP6_dF6 = 0.0
        ds_de = 0.0
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            if (.not. micro%voxel(ip1, ip2, ip3)%gas) then

                ! Get values from grid arrays
                voxel_loc = micro%voxel(ip1, ip2, ip3)
                detF = voxel_loc%detF
                Finv = voxel_loc%defgradinv
                wgtc = voxel_loc%wgtc
                sg66 = voxel_loc%sg66
                sg = voxel_loc%sg
                sgt = voxel_loc%sgt

                ! Convert 3x3 tensors to 1x6 tensor basis
                call chg_basis(sg6, sg, aux66, aux3333, 2, 6)
                call chg_basis(sgt6, sgt, aux66, aux3333, 2, 6)

                ! get plastic strain rate and derivative
                call plastic_strain_rate(sg6, voxel_loc, props, edotp6, dedotp6_dsg6)
                call elastic_strain_rate(sg6, sgt6, voxel_loc, props, edote6, dedote6_dsg6)

                ! derivative
                dsg6_de6 = dedote6_dsg6 + dedotp6_dsg6
                dsg6_de6_tref = dedote6_dsg6 * props%process(ipr)%tdot / props%process(ipr)%tdotref + &
                                dedotp6_dsg6
                call lu_inverse(dsg6_de6, 6)
                call lu_inverse(dsg6_de6_tref, 6)
                call chg_basis(aux6, aux33, dsg6_de6_tref, dsg_de, 3, 6)
                call chg_basis(aux6, aux33, dsg6_de6_tref, dsg_de_tref, 3, 6)

                ! reference configuration
                do l = 1, 3
                do k = 1, 3
                do j = 1, 3
                do i = 1, 3
                    dPg_dF(i, j, k, l) = 0.0
                    do m = 1, 3
                    do n = 1, 3
                        dPg_dF(i, j, k, l) = dPg_dF(i, j, k, l) + dsg_de_tref(i, m, k, n) * &
                                             Finv(j, m) * Finv(l, n) * detF
                    end do
                    end do
                end do
                end do
                end do
                end do
                call chg_basis(aux6, aux33, dPg6_dF6, dPg_dF, 4, 6)

                ds_de = ds_de + dsg_de * wgtc
                dP6_dF6 = dP6_dF6 + dPg6_dF6 * props%wgt

            end if

        end do
        end do
        end do
        ds_de = ds_de / micro%wphcsol
        dP6_dF6 = dP6_dF6 / props%wphsol

        ! de_ds calculation (for b.c.)
        ! ds_de=ds_de*xc0
        call chg_basis(aux6, aux33, aux66, ds_de, 4, 6)
        call lu_inverse(aux66, 6)
        call chg_basis(aux6, aux33, aux66, micro%de_ds, 3, 6)

        ! c0 and s0
        micro%c066 = dP6_dF6 * props%process(ipr)%xc0
        aux66 = dP6_dF6 * props%process(ipr)%xc0
        call lu_inverse(aux66, 6)
        call chg_basis(aux6, aux33, micro%c066, micro%c0, 3, 6)
        call chg_basis(aux6, aux33, aux66, micro%s0, 3, 6)

    end

    ! Interpolate a cell-centered vector field to periodic nodes via 8-cell averaging.
    subroutine cell2node_tensor1(tensor1cell, tensor1node, idim, npts1, npts2, npts3)
        implicit none

        integer, intent(in) :: npts1, npts2, npts3 ! grid points in x,y,z (cells)
        integer, intent(in) :: idim ! vector dimension
        double precision, intent(in) :: tensor1cell(idim, npts1, npts2, npts3) ! cell-centered field
        double precision, intent(out) :: tensor1node(idim, npts1 + 1, npts2 + 1, npts3 + 1) ! node-interpolated field

        integer :: i, j, k ! node indices
        integer :: i1, j1, k1 ! wrapped lower-neighbor indices
        integer :: i2, j2, k2 ! wrapped upper-neighbor indices

        do k = 1, npts3 + 1
            k1 = k - 1
            k2 = k
            if (k1 == 0) k1 = npts3
            if (k2 == npts3 + 1) k2 = 1
            do j = 1, npts2 + 1
                j1 = j - 1
                j2 = j
                if (j1 == 0) j1 = npts2
                if (j2 == npts2 + 1) j2 = 1
                do i = 1, npts1 + 1
                    i1 = i - 1
                    i2 = i
                    if (i1 == 0) i1 = npts1
                    if (i2 == npts2 + 1) i2 = 1

                    tensor1node(:, i, j, k) = (tensor1cell(:, i2, j2, k2) + &
                                               tensor1cell(:, i1, j2, k2) + tensor1cell(:, i2, j1, k2) + &
                                               tensor1cell(:, i2, j2, k1) + &
                                               tensor1cell(:, i1, j1, k2) + tensor1cell(:, i1, j2, k1) + &
                                               tensor1cell(:, i2, j1, k1) + &
                                               tensor1cell(:, i1, j1, k1)) * 0.125

                end do
            end do
        end do

    end

    subroutine write_status(imacroloop, sim_time, wall_time, ipr, imicro, tdot, iter, pct, final)

        use iso_fortran_env, only: output_unit
        implicit none
        double precision, intent(in) :: sim_time, tdot
        double precision, intent(in) :: wall_time
        integer, intent(in) :: pct
        integer, intent(in) :: imacroloop, ipr, imicro, iter
        logical, intent(in), optional :: final

        character(len=*), parameter :: CR = achar(13)
        character(len=*), parameter :: ESC = achar(27)
        character(len=256) :: status

        logical :: is_final

        is_final = .false.
        if (present(final)) is_final = final

        if (.not. is_final) then

            if (sim_time < 3600) then
                write (status, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
                    'MACROSTEP=', imacroloop, ' SimTime=', nint(sim_time), 's Wall=', nint(wall_time), &
                    's | PROC=', ipr, ' STEP=', imicro, ' TDOT=', nint(tdot), 's ITER=', iter, ' | Completion: ', pct, '%'
            else

                write (status, '(A,I0,A,F5.1,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
                    'MACROSTEP=', imacroloop, ' SimTime=', sim_time / 3600, 'hrs Wall=', nint(wall_time), &
                    's | PROC=', ipr, ' STEP=', imicro, ' TDOT=', nint(tdot), 's ITER=', iter, ' | Completion: ', pct, '%'
            end if

        else

            if (sim_time < 3600) then
                write (status, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
                    'MACROSTEP=', imacroloop - 1, ' SimTime=', nint(sim_time), 's Wall=', nint(wall_time), &
                    's | PROC=', ipr, ' STEP=', imicro, ' TDOT=', nint(tdot), 's ITER=', iter, ' | Completion: ', 100, '%'
            else

                write (status, '(A,I0,A,F5.1,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
                    'MACROSTEP=', imacroloop - 1, ' SimTime=', sim_time / 3600, 'hrs Wall=', nint(wall_time), &
                    's | PROC=', ipr, ' STEP=', imicro, ' TDOT=', nint(tdot), 's ITER=', iter, ' | Completion: ', 100, '%'
            end if

        end if

        write (output_unit, '(A)', advance='no') CR//ESC//'[2K'//CR//trim(status)
        flush (output_unit)

    end

    integer function kinematic_substep(micro, props, tdot_step, max_volinc, max_rotang, max_stretchinc, min_detF) result(nsub)
        use types
        implicit none

        type(micro_type), intent(in) :: micro
        type(props_type), intent(in) :: props
        double precision, intent(in) :: tdot_step
        double precision, intent(out) :: max_volinc, max_rotang, max_stretchinc, min_detF

        double precision, parameter :: volinc_lim = 5.0d-2
        double precision, parameter :: rotang_lim = 2.0d0
        double precision, parameter :: stretchinc_lim = 5.0d-2
        double precision :: worst_ratio = 0.0d0

        call kinematic_increments(micro, props, tdot_step, max_volinc, max_rotang, max_stretchinc, min_detF)

        worst_ratio = max(1.0d0, max_volinc / volinc_lim, max_rotang / rotang_lim, max_stretchinc / stretchinc_lim)
        nsub = max(1, ceiling(worst_ratio))

    end function kinematic_substep

    subroutine kinematic_increments(micro, props, tdot_step, max_volinc, max_rotang, max_stretchinc, min_detF)
        use types
        use tensor_functions
        implicit none

        type(micro_type), intent(in) :: micro
        type(props_type), intent(in) :: props
        double precision, intent(in) :: tdot_step
        double precision, intent(out) :: max_volinc, max_rotang, max_stretchinc, min_detF

        double precision :: L(3, 3), Finc(3, 3), R(3, 3), U(3, 3)
        double precision :: Jcur, Jinc, Jnext, volinc, rotang, stretchinc
        double precision :: arg
        double precision :: Q(3, 3), w(3), work(64)
        double precision :: tmp(3, 3), logU(3, 3)
        integer :: ip1, ip2, ip3, ii, info, lwork

        max_volinc = 0.0d0
        max_rotang = 0.0d0
        max_stretchinc = 0.0d0
        min_detF = huge(1.0d0)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP PRIVATE(ip1, ip2, ip3, ii, info, lwork, L, Finc, R, U, Jcur, Jinc, Jnext, volinc, rotang, stretchinc, arg, Q, w, work, tmp, logU) &
        !$OMP SHARED(micro, props, tdot_step) &
        !$OMP REDUCTION(max:max_volinc, max_rotang, max_stretchinc) &
        !$OMP REDUCTION(min:min_detF)
        lwork = size(work)
        !$OMP DO COLLAPSE(3)
        do ip3 = 1, props%npts3
        do ip2 = 1, props%npts2
        do ip1 = 1, props%npts1

            L = micro%voxel(ip1, ip2, ip3)%velgrad
            Finc = matrix_exp_adaptive(L * tdot_step)

            Jcur = micro%voxel(ip1, ip2, ip3)%detF
            Jinc = determinant33(Finc)
            Jnext = Jcur * Jinc

            volinc = abs(log(max(Jinc, tiny(1.0d0))))

            call polar_dcmp(Finc, R, U)
            arg = 0.5d0 * (R(1, 1) + R(2, 2) + R(3, 3) - 1.0d0)
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

            max_volinc = max(max_volinc, volinc)
            max_rotang = max(max_rotang, rotang)
            max_stretchinc = max(max_stretchinc, stretchinc)
            min_detF = min(min_detF, Jnext)

        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine kinematic_increments

    subroutine take_snapshot(micro, props, ipr_local)
        use types
        implicit none

        type(micro_type), intent(inout) :: micro
        type(props_type), intent(inout) :: props
        integer, intent(in) :: ipr_local
        integer :: ip1_local, ip2_local, ip3_local

        props%process(ipr_local)%dsim_snapshot = props%process(ipr_local)%dsim
        props%process(ipr_local)%scauchy_snapshot = props%process(ipr_local)%scauchy

        !$OMP PARALLEL DO DEFAULT(NONE) SHARED(micro, props) PRIVATE(ip1_local, ip2_local, ip3_local) COLLAPSE(3) SCHEDULE(STATIC)
        do ip3_local = 1, props%npts3
        do ip2_local = 1, props%npts2
        do ip1_local = 1, props%npts1
            micro%voxel(ip1_local, ip2_local, ip3_local)%velgrad_snapshot = micro%voxel(ip1_local, ip2_local, ip3_local)%velgrad
            micro%voxel(ip1_local, ip2_local, ip3_local)%velgradold_snapshot = micro%voxel(ip1_local, ip2_local, ip3_local)%velgradold
            micro%voxel(ip1_local, ip2_local, ip3_local)%sg_snapshot = micro%voxel(ip1_local, ip2_local, ip3_local)%sg
            micro%voxel(ip1_local, ip2_local, ip3_local)%edotp_snapshot = micro%voxel(ip1_local, ip2_local, ip3_local)%edotp
            micro%voxel(ip1_local, ip2_local, ip3_local)%gamdot_snapshot = micro%voxel(ip1_local, ip2_local, ip3_local)%gamdot
            micro%voxel(ip1_local, ip2_local, ip3_local)%trialtau_snapshot = micro%voxel(ip1_local, ip2_local, ip3_local)%trialtau
        end do
        end do
        end do
        !$OMP END PARALLEL DO

    end subroutine take_snapshot

    subroutine rollback(micro, props, ipr_local)
        use types
        implicit none

        type(micro_type), intent(inout) :: micro
        type(props_type), intent(inout) :: props
        integer, intent(in) :: ipr_local
        integer :: ip1_local, ip2_local, ip3_local

        props%process(ipr_local)%dsim = props%process(ipr_local)%dsim_snapshot
        props%process(ipr_local)%scauchy = props%process(ipr_local)%scauchy_snapshot

        !$OMP PARALLEL DO DEFAULT(NONE) SHARED(micro, props) PRIVATE(ip1_local, ip2_local, ip3_local) COLLAPSE(3) SCHEDULE(STATIC)
        do ip3_local = 1, props%npts3
        do ip2_local = 1, props%npts2
        do ip1_local = 1, props%npts1
            micro%voxel(ip1_local, ip2_local, ip3_local)%velgrad = micro%voxel(ip1_local, ip2_local, ip3_local)%velgrad_snapshot
            micro%voxel(ip1_local, ip2_local, ip3_local)%velgradold = micro%voxel(ip1_local, ip2_local, ip3_local)%velgradold_snapshot
            micro%voxel(ip1_local, ip2_local, ip3_local)%sg = micro%voxel(ip1_local, ip2_local, ip3_local)%sg_snapshot
            micro%voxel(ip1_local, ip2_local, ip3_local)%edotp = micro%voxel(ip1_local, ip2_local, ip3_local)%edotp_snapshot
            micro%voxel(ip1_local, ip2_local, ip3_local)%gamdot = micro%voxel(ip1_local, ip2_local, ip3_local)%gamdot_snapshot
            micro%voxel(ip1_local, ip2_local, ip3_local)%trialtau = micro%voxel(ip1_local, ip2_local, ip3_local)%trialtau_snapshot
        end do
        end do
        end do
        !$OMP END PARALLEL DO

    end subroutine rollback

end module various_functions
