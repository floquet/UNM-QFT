! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
program quadrature

    use, intrinsic :: iso_fortran_env,  only : compiler_version, compiler_options

    use mConstants,                     only : pi, one, stdout
    ! use mParameters,                    only : nDim
    use mPartition,                     only : partitions, GrandPartition
    use mSetPrecision,                  only : ip

    implicit none

    type ( GrandPartition ), target :: allPartitions! ( 1 : nDim )
    type ( partitions ), pointer :: myPart => null ( )
    type ( partitions ) :: myPartition

    integer ( ip ) :: k = 0, kPartition = 0, n = 0
    integer        :: alloc_status = 0, status = 0

    character ( len = 256 ) :: alloc_msg = ''
    character ( len =   8 ) :: date      = ''  ! ccyymmdd
    character ( len =  10 ) :: time      = ''  ! hhmmss.sss
    character ( len =   5 ) :: zone      = ''  ! (+-)hhmm, difference from Coordinated Universal Time (UTC)


        ! assemble partitions
        ! Partition I
        kPartition = 1
        n = 3
        myPartition % nCell = n
        write ( *, 900 ) 'myPartition % nCell = ', myPartition % nCell
        allocate ( myPartition % mesh ( 0 : n ), stat = alloc_status, errmsg = alloc_msg )
        !write ( *, '( * ( g0, " - " ) )' ) 'allocate: myPartition % mesh = ', myPartition % mesh
        myPartition % mesh ( : ) = [ ( pi * k / n, k = 0, n ) ]
        !write ( *, '( * ( g0, " - " ) )' ) 'populate: myPartition % mesh = ', myPartition % mesh
        !allPartitions ( 1 ) % partition = myPartition
        allPartitions % partition ( 1 ) = myPartition
        write ( *, 900 ) 'allPartitions % partition ( 1 ) % nCell = ', allPartitions % partition ( 1 ) % nCell
        write ( *, 900 ) 'allPartitions % partition ( 1 ) % mesh  = ', allPartitions % partition ( 1 ) % mesh

        ! Partition II
        kPartition = 2
        myPart => allPartitions % partition ( kPartition )
        n = 4
        ! allPartitions % partition ( kPartition ) % nCell = n
        myPart % nCell = n
        allocate ( allPartitions % partition ( kPartition ) % mesh ( 0 : n ), stat = alloc_status, errmsg = alloc_msg )
        ! allPartitions % partition ( kPartition ) % mesh ( : ) = [ ( one * k / n, k = 0, n ) ]
        myPart % mesh ( : ) = [ ( one * k / n, k = 0, n ) ]
        write ( *, 900 ) 'allPartitions % partition ( 2 ) % nCell = ', allPartitions % partition ( kPartition ) % nCell
        write ( *, 900 ) 'allPartitions % partition ( 2 ) % mesh ( : )  = ', allPartitions % partition ( kPartition ) % mesh ( : )
        !write ( *, '( * ( g0, " - " ) )' ) 'output = ', allPartitions ( 1 ) % partition ( 1 ) % nCell
        !write ( *, '( * ( g0, " - " ) )' ) 'output = ', allPartitions ( 1 ) % partition ( 1 ) % mesh ( : )
        !myPartition % partition ( 1 ) % nCell = n
        !write ( *, '( * ( g0, " - " ) )' ) 'myPartition ( 1 ) % partition % nCell = ', myPartition ( 1 ) % partition % nCell
        !allocate ( myPartition ( 1 ) % partition % mesh ( 1 : n + 1 ), stat = alloc_status, errmsg = alloc_msg )
        !myPartition ( 1 ) % partition % mesh = [ ( k, k = 0, 3 ) ]

        ! Partition III
        kPartition = 3
        myPart => allPartitions % partition ( kPartition )
        n = 10
        myPart % nCell = n
        status = myPart % allocate_mesh ( )
        myPart % mesh ( : ) = [ ( one * k / n, k = 0, n ) ]
        write ( *, 900 ) 'allPartitions % partition ( 3 ) % nCell = ', allPartitions % partition ( kPartition ) % nCell
        write ( *, 900 ) 'allPartitions % partition ( 3 ) % mesh ( : )  = ', allPartitions % partition ( kPartition ) % mesh ( : )

        call date_and_time ( date = date, time = time, zone = zone )

        write ( stdout, '( /, "date and time: ", * ( g0, 2X ) )' ) date, time, zone
        write ( stdout, '( /, "Fortran compiler version:    ", g0    )' ) compiler_version ( )
        write ( stdout, '(    "Fortran compilation options: ", g0, / )' ) compiler_options ( )

    stop 'successful completion for quadrature . . .'

    900 format ( * ( g0, " - " ) )

end program quadrature
