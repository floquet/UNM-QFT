! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
program AaM0

    use, intrinsic :: iso_fortran_env,  only : compiler_version, compiler_options

    use mConstants,                     only : stdout, zero, one, mySeed, biggest, mille
    use mFields,                        only : housekeeping_sub, fields
    use mFileHandling,                  only : safeopen_readonly
    use mExtents,                       only : extents
    use mParameterSets,                 only : nParameterSets, ParameterCollection, load_parameter_sets_fcn
    use mRandoms,                       only : init_random_seed_sub!, SeedUsed
    use mSetPrecision,                  only : ip, rp
    use mTimeStamp,                     only : timestamp
    !use mFields,                        only : housekeeping_sub, fields

    implicit none

    ! rank 2
    ! rank 1
    real ( rp ) :: maxPhi = zero, maxGphi = zero, maxDphi = zero
    real ( rp ) :: as = zero, at = zero, Mass = zero, m = zero, df = zero
    real ( rp ) :: cpu_time_start = zero, cpu_time_stop = zero, cpu_time_elapsed = zero
    real ( rp ) :: random = zero
    real ( rp ) :: highphi = zero, highgphi = zero, highdphi = zero
    real ( rp ) :: minA = biggest, outoftable = zero
    real ( rp ) :: phistep = zero, gphistep = zero, dphistep = zero

    integer ( ip ) :: index = 0
    integer        :: io_in_run_parameters = 0, io_in_farray = 0, success = -1
    integer ( ip ) :: i = 0, j = 0, k = 0, l = 0, kPS = 0

    ! derived types
    type ( fields ), target   :: myFields
    type ( extents ), pointer :: extent => null ( )

    character ( len = 64 ) :: tablename = '', temp = '', farray = '', root = ''
    character ( len = *  ), parameter :: input_file = 'inM0m1a.1'

        call cpu_time ( cpu_time_start )

            ! call init_random_seed_sub ( FlagCheckOS = .true. )
            ! write ( stdout, 110 ) 'random number seed ', SeedUsed
            call init_random_seed_sub ( mySeed )
            write ( stdout, 100 ) 'first 10 random numbers:'
            do k = 1, 10
                call random_number ( random )
                write ( stdout, 100 ) k, '. ', random
            end do

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

            call housekeeping_sub ( myFields )

            ! loop over parameter sets
            success = load_parameter_sets_fcn ( )
            if ( success /= 0 ) stop 'Fatal error: parameter sets failed to load.'

            do kPS = 1, nParameterSets
                Mass = ParameterCollection ( kPS ) % Mass
                m    = ParameterCollection ( kPS ) % m
                at   = ParameterCollection ( kPS ) % at
                as   = ParameterCollection ( kPS ) % a
                write ( stdout, 100 ) kPS, ': Mass = ', Mass
                write ( stdout, 100 ) kPS, ': m    = ', m

                phistep  = maxPhi  / real ( extent % Nphi,  rp )
                gphistep = maxGphi / real ( extent % Ngphi, rp )
                dphistep = maxDphi / real ( extent % Ndphi, rp )

                if ( temp == 'hot' ) then
                    io_in_farray = safeopen_readonly ( farray )
                    read ( io_in_farray, * ) myFields % f
                    write ( stdout, 100 ) 'The temp is ', temp
                else if ( temp == 'cold') then
                    do i = 1, extent % Ns
                        do j = 1, extent % Ns
                            do k = 1, extent % Ns
                                do l = 1, extent % Nt
                                    call random_number ( random )
                                    myFields % f ( i, j, k, l ) = ( one - random ) * mille
                                end do ! l
                            end do ! k
                        end do ! j
                    end do ! i
                    write ( stdout, 100 ) 'The temp is ', temp
                    else
                        write ( stdout, 100 ) 'I need to know the temperature.'
                        stop 'Fatal error - no temperature.'
                    end if
            end do ! kPS
            extent => null ( )

        call cpu_time ( cpu_time_stop  )
        cpu_time_elapsed = cpu_time_stop - cpu_time_start

        write ( stdout, * )
        write ( stdout, 100 ) 'cpu seconds: ', cpu_time_elapsed
        write ( stdout, 100 ) 'timestamp: ', timestamp ( )

        write ( stdout, * )
        write ( stdout, 100 ) 'Fortran compiler version: ', compiler_version ( )
        write ( stdout, * )
        write ( stdout, 100 ) 'Fortran compilation options: ', compiler_options ( )

        write ( stdout, * )
        stop 'successful completion for AaM0.f08 . . .'

    100 format ( * ( g0 ) )
    !110 format ( * ( g0, ', ' ) )

end program AaM0

! rditldmt@ITLDMT-MD-O2034:alpha $ pwd
! /Users/rditldmt/Documents/GitHub Desktop/UNM-QFT/codes/qft/greens/alpha
! rditldmt@ITLDMT-MD-O2034:alpha $ date
! Tue Jan 17 16:17:40 CST 2017
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
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_parameter_sets.o mod_parameter_sets.f08
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_randoms.o mod_randoms.f08
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_time_stamp.o mod_time_stamp.f08
! gfortran -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o AaM0.o AaM0.f08
! AaM0.f08:94:7:
!
!      110 format ( * ( g0, ', ' ) )
!        1
! Warning: Label 110 at (1) defined but not used [-Wunused-label]
! gfortran -g -o AaM0 AaM0.o mod_constants.o mod_extents.o mod_fields.o mod_file_handling.o mod_parameter_sets.o mod_randoms.o mod_set_precision.o mod_time_stamp.o
! rditldmt@ITLDMT-MD-O2034:alpha $ ./AaM0
! first 10 random numbers:
! 1. 0.78157313515666638
! 2. 0.85276506658379059
! 3. 0.86172951474681547
! 4. 0.94357548881740327
! 5. 0.74924959946961245
! 6. 0.52760745048508229
! 7. 0.38451498498411230E-001
! 8. 0.64528496827780668
! 9. 0.87860794998854819
! 10. 0.45900666593867046E-001
! Reading parameters in file inM0m1a.1.
! extent % Nphi = 100
! myFields % myExtents % Ngphi = 100
!
! Parameter sets loaded.
!
! 1: Mass = 1.0000000000000000
! 1: m    = 0.0000000000000000
! 2: Mass = 0.0000000000000000
! 2: m    = 1.0000000000000000
!
! cpu seconds: 0.67799999999999978E-003
! timestamp: 2017-01-17  16:17:51  UCT-0600
!
! Fortran compiler version: GCC version 7.0.0 20170115 (experimental)
!
! Fortran compilation options: -fPIC -feliminate-unused-debug-symbols -mmacosx-version-min=10.11.6 -mtune=core2 -auxbase-strip AaM0.o -g -Og -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Wpedantic -Wuse-without-only -ffpe-trap=denormal -fbacktrace -fcheck=bounds -fmax-errors=5
!
! STOP successful completion for AaM0.f08 . . .
! rditldmt@ITLDMT-MD-O2034:alpha $ gcc --version
! Configured with: --prefix=/Applications/Xcode.app/Contents/Developer/usr --with-gxx-include-dir=/usr/include/c++/4.2.1
! Apple LLVM version 8.0.0 (clang-800.0.42.1)
! Target: x86_64-apple-darwin15.6.0
! Thread model: posix
! InstalledDir: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin
