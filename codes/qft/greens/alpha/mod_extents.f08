module mExtents

    use mSetPrecision,                  only : ip, rp

    implicit none

    type :: extents
        integer ( ip ) :: Nphi, Ngphi, Ndphi
        integer ( ip ) :: Nsweeps, Ns, Nt
        integer ( ip ) :: volume_ip
        integer ( ip ) :: kount
        real ( rp ) :: volume_rp
        real ( rp ) :: as, at
        real ( rp ) :: avolume
    contains
        private
        procedure, public :: volume => volume_sub
    end type extents

    private :: volume_sub

contains

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine volume_sub ( me )

        class ( extents ), target :: me

            me % avolume   = me % as ** 3 * me % at
            me % volume_ip = me % Ns ** 3 * me % Nt
            me % volume_rp = real ( me % volume_ip, rp )

            me % kount = me % volume_ip * me % Nsweeps

    end subroutine volume_sub

end module mExtents
