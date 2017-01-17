module mConstants

    use, intrinsic :: iso_fortran_env,  only : stdin => input_unit, stdout => output_unit

    use mSetPrecision,                  only : ip, rp

    implicit none

    ! precision controlled numbers
    real ( rp ), parameter :: one = 1.0_rp, two = 2.0_rp, three = 3.0_rp

    real ( rp ), parameter :: pi = acos ( -one )

end module mConstants
