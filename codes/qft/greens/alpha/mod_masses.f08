module mMasses

    use mSetPrecision,                  only : ip, rp

    implicit none

    type :: masses
        real ( rp ) :: Mass
        real ( rp ) :: m
    end type masses

end module mMasses
