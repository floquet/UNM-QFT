module mConstants

    use, intrinsic :: iso_fortran_env,  only : stdin => input_unit, stdout => output_unit

    use mSetPrecision,                  only : rp

    implicit none

    ! precision controlled numbers
    real ( rp ), parameter :: zero = 0.0_rp, one = 1.0_rp, two = 2.0_rp, three = 3.0_rp

    real ( rp ), parameter :: biggest = huge ( one )

    real ( rp ), parameter :: pi = acos ( -one )

end module mConstants
