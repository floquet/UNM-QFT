! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
program integer_formats
    use, intrinsic :: iso_fortran_env, only : compiler_options, compiler_version
    implicit none
    integer :: myInteger = 0, myPower = 0

        do myPower = 0, 6
            myInteger = 10 ** myPower
            write ( *, 100 ) myInteger, myInteger, myInteger
        end do ! myPower

        write ( *, 200 ) 'Fortran compiler version: ', compiler_version ()
        write ( *, 200 ) 'Fortran compiler options: ', compiler_options ()

        stop 'successful run for integer_formats...'

    100 format ( 'padded: ', I7.7, '; unpadded: ', I7, '; g0: ', g0 )
    200 format ( g0, g0 )

    end program integer_formats

! dantopa@Muntz-Szasz:demos $ echo $gflags
! -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only
! dantopa@Muntz-Szasz:demos $
! dantopa@Muntz-Szasz:demos $ gfortran integer_formats.f08 -o integer_formats
! dantopa@Muntz-Szasz:demos $ ./integer_formats
! padded: 0000001; unpadded:       1; g0: 1
! padded: 0000010; unpadded:      10; g0: 10
! padded: 0000100; unpadded:     100; g0: 100
! padded: 0001000; unpadded:    1000; g0: 1000
! padded: 0010000; unpadded:   10000; g0: 10000
! padded: 0100000; unpadded:  100000; g0: 100000
! padded: 1000000; unpadded: 1000000; g0: 1000000
! Fortran compiler version: GCC version 7.0.0 20161120 (experimental)
! Fortran compiler options: -fPIC -mmacosx-version-min=10.12.1 -mtune=core2
! STOP successful run for integer_formats...
! dantopa@Muntz-Szasz:demos $
! dantopa@Muntz-Szasz:demos $ date
! Sun Nov 27 09:12:26 MST 2016
! dantopa@Muntz-Szasz:demos $ pwd
! /Users/dantopa/Documents/git desktop/UNM-QFT/demos
