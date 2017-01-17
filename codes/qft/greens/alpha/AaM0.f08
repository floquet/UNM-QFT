! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
program AaM0

    use, intrinsic :: iso_fortran_env,  only : compiler_version, compiler_options

    use mConstants,                     only : stdout, zero!, one
    use mFields,                        only : fields
    use mFileHandling,                  only : safeopen_readonly
    use mExtents,                       only : extents
    use mSetPrecision,                  only : ip, rp
    use mTimeStamp,                     only : timestamp

    implicit none

    ! rank 2
    ! rank 1
    real ( rp ) :: maxPhi = zero, maxGphi = zero, maxDphi = zero
    real ( rp ) :: as = zero, at = zero, Mass = zero, m = zero, df = zero
    real ( rp ) :: cpu_time_start = zero, cpu_time_stop = zero, cpu_time_elapsed = zero

    !integer ( ip ) :: Nphi = 0, Ngphi = 0, Ndphi = 0
    !integer ( ip ) :: Nsweeps = 0, Ns = 0, Nt = 0
    integer ( ip ) :: index = 0
    integer        :: io_in_run_parameters = 0

    ! derived types
    type ( fields ), target  :: myFields
    type ( extents), pointer :: extent => null ( )

    character ( len = 64 ) :: tablename = '', temp = '', farray = '', root = ''
    character ( len = *  ), parameter :: input_file = 'inM0m1a.1'

        call cpu_time ( cpu_time_start )

            io_in_run_parameters = safeopen_readonly ( input_file )
            write ( stdout, 100 ) 'Reading parameters in file ', input_file, '.'

            ! ./AaM0 < input_file
            extent => myFields % myExtents
            read ( io_in_run_parameters, * ); read ( io_in_run_parameters, * ) maxPhi, maxGphi, maxDphi
            read ( io_in_run_parameters, * ); read ( io_in_run_parameters, * ) extent % Nphi, &
                                                                               extent % Ngphi, &
                                                                               extent % Ndphi
            read ( io_in_run_parameters, * ); read ( io_in_run_parameters, * ) as, at, Mass, m, df
            read ( io_in_run_parameters, * ); read ( io_in_run_parameters, * ) tablename, temp, root, farray, index
            read ( io_in_run_parameters, * ); read ( io_in_run_parameters, * ) extent % Nsweeps, &
                                                                               extent % Ns, &
                                                                               extent % Nt
            close ( io_in_run_parameters )

            write ( stdout, 100 ) 'extent % Nphi = ', extent % Nphi
            write ( stdout, 100 ) 'myFields % myExtents % Ngphi = ', myFields % myExtents % Ngphi

            extent => null ( )

        call cpu_time ( cpu_time_stop  )
        cpu_time_elapsed = cpu_time_stop - cpu_time_start

        write ( stdout, 100 ) 'cpu seconds: ', cpu_time_elapsed
        write ( stdout, 100 ) 'timestamp: ', timestamp ( )

        write ( stdout, '( /, "Fortran compiler version:    ", g0    )' ) compiler_version ( )
        write ( stdout, '(    "Fortran compilation options: ", g0, / )' ) compiler_options ( )

        stop 'successful completion for AaM0.f08 . . .'

    100 format ( * ( g0 ) )

end program AaM0

! rditldmt@ITLDMT-MD-O2034:alpha $ date
! Mon Jan 16 18:20:52 CST 2017
! rditldmt@ITLDMT-MD-O2034:alpha $ pwd
! /Users/rditldmt/hpc/fortran/projects/qft/greens/alpha
! rditldmt@ITLDMT-MD-O2034:alpha $ make
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_set_precision.o mod_set_precision.f08
! make: Circular mod_constants.o <- mod_constants.o dependency dropped.
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_constants.o mod_constants.f08
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_extents.o mod_extents.f08
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_fields.o mod_fields.f08
! mod_fields.f08:52:37:
!
!      function allocate_rank_1_fcn ( me, array, length ) result  ( fcn_success )
!                                      1
! Warning: Unused dummy argument 'me' at (1) [-Wunused-dummy-argument]
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_file_handling.o mod_file_handling.f08
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_randoms.o mod_randoms.f08
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_time_stamp.o mod_time_stamp.f08
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o AaM0.o AaM0.f08
! gfortran -g -o AaM0 AaM0.o mod_constants.o mod_extents.o mod_fields.o mod_file_handling.o mod_randoms.o mod_set_precision.o mod_time_stamp.o
! rditldmt@ITLDMT-MD-O2034:alpha $ ./AaM0
! Reading parameters in file inM0m1a.1.
! extent % Nphi = 100
! myFields % myExtents % Ngphi = 100
! cpu seconds: 0.25299999999999975E-003
! timestamp: 2017-01-16  18:20:59  UCT-0600
!
! Fortran compiler version:    GCC version 7.0.0 20170115 (experimental)
! Fortran compilation options: -fPIC -feliminate-unused-debug-symbols -mmacosx-version-min=10.11.6 -mtune=core2 -auxbase-strip AaM0.o -g -Og -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Wpedantic -Wuse-without-only -ffpe-trap=denormal -fbacktrace -fcheck=bounds -fmax-errors=5
!
! STOP successful completion for AaM0.f08 . . .
