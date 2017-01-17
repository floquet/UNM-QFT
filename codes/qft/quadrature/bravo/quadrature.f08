! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
program quadrature

    use, intrinsic :: iso_fortran_env,  only : compiler_version, compiler_options

    use mConstants,                     only : pi, stdout
    use mParameters,                    only : nDim
    use mPartition,                     only : partitions, GrandPartition
    use mSetPrecision,                  only : ip, rp
    use mTimeStamp,                     only : timestamp

    implicit none

    type ( GrandPartition ), target  :: allPartitions
    type ( partitions ),     pointer :: myPart => null ( )

    integer ( ip ) :: k = 0, kPartition = 0, n = 0
    integer        :: status = 0

    !character ( len = 256 ) :: alloc_msg = ''

        ! assemble partitions
        ! Partition I
        kPartition = 1
        myPart => allPartitions % partition ( kPartition )
        n = 10
        myPart % nCell = n
        status = myPart % allocate_mesh ( )
        myPart % mesh ( : ) = [ ( pi * real ( k, rp ) / real ( n, rp ), k = 0, n ) ]

        ! Partition II
        kPartition = 2
        myPart => allPartitions % partition ( kPartition )
        n = 4
        myPart % nCell = n
        status = myPart % allocate_mesh ( )
        myPart % mesh ( : ) = [ ( real ( k, rp ) / real ( n, rp ), k = 0, n ) ]

        ! Partition III
        kPartition = 3
        myPart => allPartitions % partition ( kPartition )
        n = 10
        myPart % nCell = n
        status = myPart % allocate_mesh ( )
        myPart % mesh ( : ) = [ ( real ( k, rp ) / real ( n, rp ), k = 0, n ) ]

        ! Partition IV
        kPartition = 4
        myPart => allPartitions % partition ( kPartition )
        n = 3
        myPart % nCell = n
        status = myPart % allocate_mesh ( )
        myPart % mesh ( : ) = [ ( real ( k, rp ), k = 0, n ) ]

        do kPartition = 1, nDim
            myPart => allPartitions % partition ( kPartition )
            write ( *, '( "" )' )
            write ( *, 900 ) 'length of partition ', kPartition, ' = ', myPart % nCell
            write ( *, 910 ) 'elements', myPart % mesh ( : )
        end do

        myPart => null ( )
        write ( stdout, 900 ) 'Partitions established.'

        ! count the cells
        allPartitions % nTotalCells = product ( [ ( allPartitions % partition ( kPartition ) % nCell, kPartition = 1, nDim ) ] )
        write ( stdout, 900 ) 'nTotalCells = product nCell ( 1 : nDim ) = ', allPartitions % nTotalCells

        write ( stdout, 900 ) 'time stamp: ', timestamp ( )

        write ( stdout, '( /, "Fortran compiler version:    ", g0    )' ) compiler_version ( )
        write ( stdout, '(    "Fortran compilation options: ", g0, / )' ) compiler_options ( )

    stop 'successful completion for quadrature . . .'

    900 format ( * ( g0 ) )
    910 format ( * ( g0, " - " ) )

end program quadrature

! rditldmt@ITLDMT-MD-O2034:bravo $ pwd
! /Users/rditldmt/hpc/fortran/projects/qft/quadrature/bravo
! rditldmt@ITLDMT-MD-O2034:bravo $ date
! Tue Jan 10 17:57:40 CST 2017
! rditldmt@ITLDMT-MD-O2034:bravo $ make
! mpif90 -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_set_precision.o mod_set_precision.f08
! mpif90 -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_constants.o mod_constants.f08
! mpif90 -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_mpif_sgi.o mod_mpif_sgi.f08
! mpif90 -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_parameters.o mod_parameters.f08
! mpif90 -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_partition.o mod_partition.f08
! mpif90 -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o mod_point_in_space.o mod_point_in_space.f08
! mpif90 -c -g -ffpe-trap=denormal -fbacktrace -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Og -pedantic -fcheck=bounds -fmax-errors=5 -Wuse-without-only -o quadrature.o quadrature.f08
! mpif90 -g -o quadrature mod_constants.o mod_mpif_sgi.o mod_parameters.o mod_partition.o mod_point_in_space.o mod_set_precision.o quadrature.o
! rditldmt@ITLDMT-MD-O2034:bravo $ ./quadrature
!
! length of partition 1 = 10
! elements - .0000000000000000 - .31415926535897931 - .62831853071795862 - .94247779607693793 - 1.2566370614359172 - 1.5707963267948966 - 1.8849555921538759 - 2.1991148575128552 - 2.5132741228718345 - 2.8274333882308138 - 3.1415926535897931 -
!
! length of partition 2 = 4
! elements - .0000000000000000 - .25000000000000000 - .50000000000000000 - .75000000000000000 - 1.0000000000000000 -
!
! length of partition 3 = 10
! elements - .0000000000000000 - .10000000000000001 - .20000000000000001 - .29999999999999999 - .40000000000000002 - .50000000000000000 - .59999999999999998 - .69999999999999996 - .80000000000000004 - .90000000000000002 - 1.0000000000000000 -
!
! length of partition 4 = 3
! elements - .0000000000000000 - 1.0000000000000000 - 2.0000000000000000 - 3.0000000000000000 -
!
! date and time: 20170110  175756.140  -0600
!
! Fortran compiler version:    GCC version 5.4.0
! Fortran compilation options: -I /opt/local/include/mpich-gcc5 -I /opt/local/include/mpich-gcc5 -fPIC -feliminate-unused-debug-symbols -mmacosx-version-min=10.11.6 -m64 -mtune=core2 -auxbase-strip quadrature.o -g -Og -Wall -Waliasing -Wconversion-extra -Wextra -Wsurprising -Wimplicit-procedure -Wintrinsics-std -Wpedantic -Wuse-without-only -ffpe-trap=denormal -fbacktrace -fcheck=bounds -fmax-errors=5
!
! STOP successful completion for quadrature . . .
