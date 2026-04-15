module global
    implicit none
    double precision, parameter  :: pi = 4.d0 * datan(1.d0) ! pi
    complex(kind(1.0d0)), parameter  :: ximag = (0.d0, 1.d0) ! imaginary i
    double precision, parameter, dimension(3, 3) :: id3 = reshape((/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/), shape(id3)) ! 3x3 identity matrix
    double precision, parameter, dimension(6, 6) :: id6 = reshape((/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &  ! 6x6 identity matrix (Voigt)
                                                                    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, &
                                                                    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, &
                                                                    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, &
                                                                    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, &
                                                                    0.0, 0.0, 0.0, 0.0, 0.0, 1.0/), shape(id6))
    double precision, parameter, dimension(3, 3, 3) :: e = reshape((/0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, & ! Levi-Civita permutation tensor
                                                                     0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0/), shape(e))
    integer, parameter :: nphmx = 3 ! max phases
    integer, parameter :: nmodmx = 3 ! max deformation modes per phase
    integer, parameter :: nsysmx = 24 ! max slip systems
    integer, parameter :: ntwmmx = 0 ! max twinning modes (0 => disabled)
    integer, parameter :: ngrmx = 1000 ! max grains
    integer :: ipr ! active load/process index
    integer :: imacroloop ! macrostep counter
    integer :: imicro ! active microstep index within a process
    integer :: iter ! global FFT equilibrium iteration index
    double precision :: erre ! strain/compatibility field error norm
    double precision :: errs ! stress/equilibrium field error norm
    double precision :: errsds ! legacy/additional error measure
    double precision :: errsbc ! boundary condition mismatch error norm
    double precision :: iterNRavg_accum ! accumulated avg Newton iterations

end module global
