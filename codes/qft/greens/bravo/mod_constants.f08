! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
module mConstants

    use, intrinsic :: iso_fortran_env,  only : stdin => input_unit, stdout => output_unit

    use mSetPrecision,                  only : ip, rp

    implicit none

    ! random number seed
    integer :: mySeed ( 1 : 33 ) = [ 2048438402, -405862320, 899756676, -1809184714, 1085168783, -1646620091, 691402267, &
        -1912685341, 581541989, 951441836, 1078925186, 909269592, 925394843, 1315245928, -299535180, 1623079232, -924637127, &
        -1284492704, -1558884308, 1153552125, 786484735, -333531883, 771493458, 957931328, -1001042831, -921906651, 825589163, &
        1607494709, -1144785726, -1922710716, -1401228682, 2021422086, 696467395 ]

    ! precision controlled numbers
    real ( rp ), parameter :: zero = 0.0_rp, one = 1.0_rp, two = 2.0_rp, three = 3.0_rp
    real ( rp ), parameter :: half = one / two

    real ( rp ), parameter :: biggest = huge ( one ), mille = one / 1000.0_rp

    real ( rp ), parameter :: pi = acos ( -one )

    ! formats
    character ( len = * ), parameter :: fmt_generic = '( * ( g0 ) )'

end module mConstants
