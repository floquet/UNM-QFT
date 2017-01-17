module mPointInSpace

    use mParameters,                    only : nDim
    use mSetPrecision,                  only : rp

    implicit none

    type :: PointInSpace
        real ( rp ) :: addr ( 1 : nDim )
    end type PointInSpace

end module mPointInSpace
