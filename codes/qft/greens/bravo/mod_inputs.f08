module mInputs

    use mSetPrecision,                  only : ip, rp

    type :: inputs
        integer ( ip ) :: index
        character ( len = 50 ) :: tablename, temp, farray, root, infile
    end type inputs

end module mInputs
