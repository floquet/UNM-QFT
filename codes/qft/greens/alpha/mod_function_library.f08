module mFunctionLibrary

    use mSetPrecision,                  only : ip, rp

    implicit none

contains
    function aA ( iu, ju, ku, lu, i, j, k, l ) result ( fcn_result )

        integer ( ip ), intent ( in ) :: iu, ju, ku, lu, i, j, k, l
        real ( rp ) :: fcn_result
        ! locals
        real ( rp ) :: phil, gphil, dphil, V

            phil = abs ( f ( i, j, k, l ) )
            gphil = sqrt( ( ( f ( iu,  j,  k, l ) - f ( i, j, k, l ) )**2 &
                         +  ( f (  i, ju,  k, l ) - f ( i, j, k, l ) )**2 &
                         +  ( f (  i,  j, ku, l ) - f ( i, j, k, l ) )**2 ) )
            dphil = abs ( f ( i, j, k, lu ) - f ( i, j, k, l ) )
            V = ( gphil / as )**2 + ( m * phil )**2
            fcn_result = at * sqrt ( V ) / ( as**3 * ( at**2 * V + dphil**2 ) )
     end function aA

end module mFunctionLibrary
