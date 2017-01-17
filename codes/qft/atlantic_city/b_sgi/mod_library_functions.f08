module mLibraryFunctions

    use mSetPrecision,                  only : rp

    implicit none

    real ( rp ), parameter :: one = 1.0_rp, two = 2.0_rp, three = 3.0_rp

contains

    function H ( phi, gphi, pi, as, m, Mass ) result ( myH )
        ! inputs and outputs
        real ( rp ), intent ( in ) :: phi, gphi, pi, as, m, Mass
        real ( rp )                :: myH

            myH = two * pi ** 2 - X ( phi, gphi, as, m ) - Mass ** 4
            myH = myH / sqrt ( Y ( pi, phi, gphi, as, m, Mass ) )
            myH = myH + Mass ** 4

    end function H

    function D ( phi, gphi, pi, as, m, Mass ) result ( myD )
        ! inputs and outputs
        real ( rp ), intent ( in ) :: phi, gphi, pi, as, m, Mass
        real ( rp )                :: myD

            myD = one + ( two * pi ** 2 + X ( phi, gphi, as, m ) ) / Mass ** 4
            myD = myD * Y ( pi, phi, gphi, as, m, Mass ) ** ( -three / two )

    end function D

    function F ( phi, gphi, dphi, pi, as, m, Mass ) result ( myF )
        ! inputs and outputs
        real ( rp ), intent ( in ) :: phi, gphi, dphi, pi, as, m, Mass
        real ( rp )                :: myF

            myF = dphi * pi / sqrt ( one - Y ( pi, phi, gphi, as, m, Mass ) )

    end function F

    function X ( phi, gphi, as, m ) result ( myX )
        ! inputs and outputs
        real ( rp ), intent ( in ) :: phi, gphi, as, m
        real ( rp )                :: myX
        ! local variables
        real ( rp ) :: gphias, mphi

            gphias = gphi / as  ! rescale
            mphi = m * phi

            myX = dot_product ( [ gphias, mphi ], [ gphias, mphi ] )  ! L2 sum

    end function X

    function Y ( pi, phi, gphi, as, m, Mass ) result ( myY )
        ! inputs and outputs
        real ( rp ), intent ( in ) :: pi, phi, gphi, as, m, Mass
        real ( rp )                :: myY

            myY = ( pi ** 2 - X ( phi, gphi, as, m ) ) / Mass ** 4

    end function Y

end module mLibraryFunctions
