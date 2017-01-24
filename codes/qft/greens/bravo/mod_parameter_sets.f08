module mParameterSets

    use mConstants,                     only : zero, one, stdout
    use mSetPrecision,                  only : rp

    implicit none

    integer, parameter :: nParameterSets = 2

    type :: ParameterSet
        real ( rp ) :: Mass, m
        real ( rp ) :: at, a
    end type ParameterSet

    type ( ParameterSet ), dimension ( 1 : nParameterSets ) :: ParameterCollection

contains

    function load_parameter_sets_fcn ( ) result ( success )

        integer :: success

        success = -1

        ParameterCollection ( 1 ) % Mass = one
        ParameterCollection ( 1 ) % m    = zero
        ParameterCollection ( 1 ) % at   = one
        ParameterCollection ( 1 ) % a    = one

        ParameterCollection ( 2 ) % Mass = zero
        ParameterCollection ( 2 ) % m    = one
        ParameterCollection ( 2 ) % at   = one
        ParameterCollection ( 2 ) % a    = one

        write ( stdout, '( /,"Parameter sets loaded.", / )' )

        success = 0

    end function load_parameter_sets_fcn

end module mParameterSets
