module mDomain  ! collection of cuboids

    use mSetPrecision,                  only : ip, rp

    implicit none

    type :: domain
        !type ( cuboid ), allocatable :: cuboids ( : )
        real ( rp )    :: DomainVolume
        integer ( ip ) :: nTotalCuboids
    end type domain

end module mDomain
