module mFunctionLibrary

    use mFields,                        only : fields
    use mMasses,                        only : masses
    use mSetPrecision,                  only : ip, rp

    implicit none

contains
    function aA ( iu, ju, ku, lu, i, j, k, l, myFields, myMasses ) result ( fcn_result )

        integer ( ip ),  intent ( in ) :: iu, ju, ku, lu, i, j, k, l
        type ( fields ), intent ( in ) :: myFields
        type ( masses ), intent ( in ) :: myMasses
        real ( rp ) :: fcn_result
        ! locals
        real ( rp ) :: phil, gphil, dphil, V

            phil = abs ( myFields % f ( i, j, k, l ) )
            gphil = sqrt( ( ( myFields % f ( iu,  j,  k,  l ) - myFields % f ( i, j, k, l ) )**2 &
                         +  ( myFields % f (  i, ju,  k,  l ) - myFields % f ( i, j, k, l ) )**2 &
                         +  ( myFields % f (  i,  j, ku,  l ) - myFields % f ( i, j, k, l ) )**2 ) )
            dphil = abs     ( myFields % f (  i,  j,  k, lu ) - myFields % f ( i, j, k, l ) )
            V = ( gphil / myFields % myExtents % as )**2 + ( myMasses % m * phil )**2
            fcn_result = myFields % myExtents % at * sqrt ( V ) / &
                       ( myFields % myExtents % as**3 * ( myFields % myExtents % at**2 * V + dphil**2 ) )
     end function aA

end module mFunctionLibrary
