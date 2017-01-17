module mLibraryFunctions

    use mSetPrecision,                  only : rp

    implicit none

    ! Define working precision: Hansen and Tompkins, p. 22
    real ( rp ), parameter :: one = 1.0_rp, two = 2.0_rp

contains

    function H ( phi, gphi, pi, as, m, Mass ) result ( myH )
        ! inputs and outputs
        real ( rp ), intent ( in ) :: phi, gphi, pi, as, m, Mass
        real ( rp )                :: myH
        ! intermediates
        real ( rp ) :: gphias, mphi, mass4, x, y, z

            ! useful intermediates
            gphias = gphi / as
            mphi = m * phi
            mass4 = Mass ** 4
            x = dot_product ( [ gphias, mphi], [ gphias, mphi] )
            y = one - ( pi ** 2 - x ) / mass4
            z = sqrt ( y )

            myH = ( two * pi ** 2 - x - mass4 ) / y / z
            myH = myH + mass4

    end function H

    function D ( phi, gphi, pi, as, m, Mass ) result ( myD )
        ! inputs and outputs
        real ( rp ), intent ( in ) :: phi, gphi, pi, as, m, Mass
        real ( rp )                :: myD
        ! intermediates
        real ( rp ) :: gphias, mphi, mass4, x, y, z

            ! useful intermediates
            gphias = gphi / as
            mphi = m * phi
            mass4 = Mass ** 4
            x = dot_product ( [ gphias, mphi], [ gphias, mphi] )
            y = one - ( pi ** 2 - x ) / mass4
            z = sqrt ( y )

            myD = ( one + ( two * pi ** 2 + x ) / mass4 ) / y / z

    end function D

    function F ( phi, gphi, dphi, pi, as, m, Mass ) result ( myF )
        ! inputs and outputs
        real ( rp ), intent ( in ) :: phi, gphi, dphi, pi, as, m, Mass
        real ( rp ) :: myF
        ! intermediates
        real ( rp ) :: gphias, mphi, mass4, x, y, z

            ! useful intermediates
            gphias = gphi / as
            mphi = m * phi
            mass4 = Mass ** 4
            x = dot_product ( [ gphias, mphi], [ gphias, mphi] )
            y = one - ( pi ** 2 - x ) / mass4
            z = sqrt ( y )

            myF = dphi * pi / z

    end function F

end module mLibraryFunctions
