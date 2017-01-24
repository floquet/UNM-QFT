! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
module mExtents

    use mSetPrecision,                  only : ip, rp

    implicit none

    type :: extents

        integer ( ip ) :: Nphi, Ngphi, Ndphi
        integer ( ip ) :: Nsweeps, Ns, Nt
        integer ( ip ) :: volume_ip
        integer ( ip ) :: kount

        real ( rp ) :: phistep, gphistep, dphistep
        real ( rp ) :: maxPhi, maxGPhi, maxDPhi, outoftable
        real ( rp ) :: as, at, df
        real ( rp ) :: volume_rp
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
