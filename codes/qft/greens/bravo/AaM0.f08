! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
program AaM0

    use, intrinsic :: iso_fortran_env,  only : compiler_version, compiler_options

    use mConstants,                     only : stdout, zero, mySeed, fmt_generic
    use mExtents,                       only : extents
    use mFields,                        only : fields
    use mFileHandling,                  only : safeopen_readonly
    use mInputs,                        only : inputs
    use mMasses,                        only : masses
    use mParameterSets,                 only : nParameterSets, ParameterCollection, load_parameter_sets_fcn
    use mRandoms,                       only : init_random_seed_sub!, SeedUsed
    use mSetPrecision,                  only : ip, rp
    use mTimeStamp,                     only : timestamp

    implicit none

    ! rank 1
    ! rank 0
    real ( rp ) :: cpu_time_start = zero, cpu_time_stop = zero, cpu_time_elapsed = zero
    real ( rp ) :: random = zero
    integer        :: io_in_run_parameters = 0, success = -1
    integer ( ip ) :: k = 0, kPS = 0

    ! derived types
    type ( fields ),  target  :: myFields
    type ( inputs ),  target  :: myInputs

    type ( extents ), pointer :: ex   => null ( )
    type ( masses ),  pointer :: mass => null ( )
    type ( inputs ),  pointer :: in   => null ( )

    character ( len = * ), parameter :: input_file = 'inM0m1a.1'

        call cpu_time ( cpu_time_start )

            ! call init_random_seed_sub ( FlagCheckOS = .true. )
            ! write ( stdout, 110 ) 'random number seed ', SeedUsed
            call init_random_seed_sub ( mySeed )
            write ( stdout, fmt_generic ) 'first 5 random numbers:'
            do k = 1, 5
                call random_number ( random )
                write ( stdout, fmt_generic ) k, '. ', random
            end do

            io_in_run_parameters = safeopen_readonly ( input_file ) ! file with run parameters
            write ( stdout, fmt_generic ) 'Reading parameters in file ', input_file, '.'

            ex     => myFields % myExtents ! use pointers to shorten code
            mass   => myFields % myMasses
            in     => myInputs
            ! ./AaM0 < input_file
            read ( io_in_run_parameters, * ); read ( io_in_run_parameters, * ) ex % maxPhi, ex % maxGphi, ex % maxDphi
            read ( io_in_run_parameters, * ); read ( io_in_run_parameters, * ) ex % Nphi, &
                                                                               ex % Ngphi, &
                                                                               ex % Ndphi
            read ( io_in_run_parameters, * ); read ( io_in_run_parameters, * ) ex % as, ex % at, &
                                                                               mass % Mass, mass % m, ex % df
            read ( io_in_run_parameters, * ); read ( io_in_run_parameters, * ) in % tablename, &
                                                                               in % temp, &
                                                                               in % root, &
                                                                               in % farray, &
                                                                               in % index
            read ( io_in_run_parameters, * ); read ( io_in_run_parameters, * ) ex % Nsweeps, &
                                                                               ex % Ns, &
                                                                               ex % Nt
            close ( io_in_run_parameters )

            call ex % volume ( ) ! build volume parameters

            call myFields % housekeeping ( )  ! allocate, initialize

            ! loop over parameter sets
            success = load_parameter_sets_fcn ( )
            if ( success /= 0 ) stop 'Fatal error: parameter sets failed to load.'

            do kPS = 1, nParameterSets
                mass % Mass = ParameterCollection ( kPS ) % Mass
                mass % m    = ParameterCollection ( kPS ) % m
                ex   % at   = ParameterCollection ( kPS ) % at
                ex   % as   = ParameterCollection ( kPS ) % a
                write ( stdout, fmt_generic ) kPS, ': Mass = ', mass % Mass
                write ( stdout, fmt_generic ) kPS, ': m    = ', mass % m

                ex % phistep  = ex % maxPhi  / real ( ex % Nphi,  rp )
                ex % gphistep = ex % maxGphi / real ( ex % Ngphi, rp )
                ex % dphistep = ex % maxDphi / real ( ex % Ndphi, rp )

                call myFields % thermalize ( temp = in % temp, farray = in % farray )
                write ( stdout, fmt_generic ) 'thermalized: farray = ', in % farray
                call myFields % update_f ( )
                write ( stdout, fmt_generic ) 'updated'
                call myFields % greens_two_point ( )
                call myFields % compute_sigma ( )
                call myFields % extrema ( )
                call myFields % writer ( io_output_handle = stdout, myInputs = myInputs )

            end do ! kPS
            ex   => null ( )
            mass => null ( )
            in   => null ( )

        call cpu_time ( cpu_time_stop  )
        cpu_time_elapsed = cpu_time_stop - cpu_time_start

        write ( stdout, * )
        write ( stdout, fmt_generic ) 'cpu seconds: ', cpu_time_elapsed
        write ( stdout, fmt_generic ) 'timestamp: ', timestamp ( )

        write ( stdout, * )
        write ( stdout, fmt_generic ) 'Fortran compiler version: ', compiler_version ( )
        write ( stdout, * )
        write ( stdout, fmt_generic ) 'Fortran compilation options: ', compiler_options ( )

        write ( stdout, * )
        stop 'successful completion for AaM0.f08 . . .'

end program AaM0

! rditldmt@ITLDMT-MD-O2034:bravo $ date
! Tue Jan 24 17:08:07 CST 2017

! rditldmt@ITLDMT-MD-O2034:bravo $ pwd
! /Users/rditldmt/Documents/GitHub_Desktop/UNM-QFT/codes/qft/greens/bravo

! rditldmt@ITLDMT-MD-O2034:bravo $ make
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_set_precision.o mod_set_precision.f08
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_constants.o mod_constants.f08
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_extents.o mod_extents.f08
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_masses.o mod_masses.f08
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_file_handling.o mod_file_handling.f08
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_inputs.o mod_inputs.f08
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_fields.o mod_fields.f08
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_fields_allocate.o mod_fields_allocate.f08
! mod_fields_allocate.f08:2:19:
!
!  submodule ( mFields ) smFieldsAllocate
!                    1
! Warning: USE statement at (1) has no ONLY qualifier [-Wuse-without-only]
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_fields_output.o mod_fields_output.f08
! mod_fields_output.f08:2:19:
!
!  submodule ( mFields ) smFieldsOutput
!                    1
! Warning: USE statement at (1) has no ONLY qualifier [-Wuse-without-only]
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_fields_utilities.o mod_fields_utilities.f08
! mod_fields_utilities.f08:2:19:
!
!  submodule ( mFields ) smFieldsUtilities
!                    1
! Warning: USE statement at (1) has no ONLY qualifier [-Wuse-without-only]
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_parameter_sets.o mod_parameter_sets.f08
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_randoms.o mod_randoms.f08
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_time_stamp.o mod_time_stamp.f08
! gfortran -c -g -ffpe-trap=denormal,invalid,zero -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o AaM0.o AaM0.f08
! gfortran -g -o AaM0 AaM0.o mod_constants.o mod_extents.o mod_fields.o mod_fields_allocate.o mod_fields_output.o mod_fields_utilities.o mod_file_handling.o mod_inputs.o mod_masses.o mod_parameter_sets.o mod_randoms.o mod_set_precision.o mod_time_stamp.o

! rditldmt@ITLDMT-MD-O2034:bravo $ ./AaM0
! first 5 random numbers:
! 1. 0.78157313515666638
! 2. 0.85276506658379059
! 3. 0.86172951474681547
! 4. 0.94357548881740327
! 5. 0.74924959946961245
! Reading parameters in file inM0m1a.1.
!
! Parameter sets loaded.
!
! 1: Mass = 1.0000000000000000
! 1: m    = 0.0000000000000000
! The temp is cold
! thermalized: farray = duh
! updated
! maxPhi,    maxGphi,    maxDphi
!        2.000       2.000       2.000
! Nphi,   Ngphi,    Ndphi
!        100       100       100
! as,   at,   Mass,   m,   df
!     1.000000    1.000000    1.000000    0.000000    0.010000
! tablename,     temp,    root,    farray,    index
! tableM0m1a.1 hot M0m1a.1  duh         2
! Nsweeps,   Ns,   Nt
!       1000        20        20
! Results from run 1
! E_0 = 0.10000000000000183E-003, sigma = 0.0000000000000000
! G( 0 ) = 0.50926127705580896E-006
! G( 1 ) = 0.30436188976563380E-006
! G( 2 ) = 0.26914277376277926E-006
! G( 3 ) = 0.26059225043890117E-006
! minimum of A = 0.0000000000000000
! maximum of A = 0.0000000000000000
! minimum of C = 0.0000000000000000
! maximum of C = 0.0000000000000000
! highphi  = 0.0000000000000000
! highgphi = 0.0000000000000000
! highdphi = 0.0000000000000000
! A ( highphi = 1, highgphi = 1, highdphi = 1 ) = 0.0000000000000000
! out of table = 0.0000000000000000
! minA = 0.17976931348623157E+309
! 2: Mass = 0.0000000000000000
! 2: m    = 1.0000000000000000
! The temp is cold
! thermalized: farray = duh
! updated
! maxPhi,    maxGphi,    maxDphi
!        2.000       2.000       2.000
! Nphi,   Ngphi,    Ndphi
!        100       100       100
! as,   at,   Mass,   m,   df
!     1.000000    1.000000    0.000000    1.000000    0.010000
! tablename,     temp,    root,    farray,    index
! tableM0m1a.1 hot M0m1a.1  duh         2
! Nsweeps,   Ns,   Nt
!       1000        20        20
! Results from run 1
! E_0 = 0.10000000000000183E-003, sigma = 0.0000000000000000
! G( 0 ) = 0.11730241732807630E-005
! G( 1 ) = 0.40461181515299013E-006
! G( 2 ) = 0.29223647172478797E-006
! G( 3 ) = 0.26947863812950206E-006
! minimum of A = 0.0000000000000000
! maximum of A = 0.0000000000000000
! minimum of C = 0.0000000000000000
! maximum of C = 0.0000000000000000
! highphi  = 0.0000000000000000
! highgphi = 0.0000000000000000
! highdphi = 0.0000000000000000
! A ( highphi = 1, highgphi = 1, highdphi = 1 ) = 0.0000000000000000
! out of table = 0.0000000000000000
! minA = 0.17976931348623157E+309
!
! cpu seconds: 323.54083099999997
! timestamp: 2017-01-24  17:13:57  UCT-0600
!
! Fortran compiler version: GCC version 7.0.1 20170122 (experimental)
!
! Fortran compilation options: -fPIC -feliminate-unused-debug-symbols -mmacosx-version-min=10.11.6 -mtune=core2 -auxbase-strip AaM0.o -g -Og -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Wpedantic -Wuse-without-only -ffpe-trap=denormal,invalid,zero -fbacktrace -fcheck=bounds -fmax-errors=5
!
! STOP successful completion for AaM0.f08 .
