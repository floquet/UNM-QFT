module mCuboid

    use mSetPrecision,                  only : rp
    use mParameters,                    only : nDim
    use mPointInSpace,                  only : PointInSpace

    implicit none

    type :: cuboid
        type ( PointInSpace ) :: vertex ( 1 : 2 ** nDim )
        real ( rp ) :: volume
    end type cuboid

end module mCuboid
