! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
module mInputs

    use mSetPrecision,                  only : ip, rp

    type :: inputs
        integer ( ip ) :: index
        character ( len = 50 ) :: tablename, temp, farray, root, infile
    end type inputs

end module mInputs
