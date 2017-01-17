module mExtents

    use, intrinsic :: iso_fortran_env,  only : stdin => input_unit, stdout => output_unit

    use mSetPrecision,                  only : ip

    implicit none

    type :: extents
        integer ( ip ) :: Nphi, Ngphi, Ndphi
        integer ( ip ) :: Nsweeps, Ns, Nt
    end type extents

end module mExtents
