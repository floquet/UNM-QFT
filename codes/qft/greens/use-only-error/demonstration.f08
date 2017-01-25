program demonstration
    use mModule, only : myType
    implicit none
    type ( myType ) :: example
        example % this = 1.0
        example % that = 2.0
        call example % adder ( )
        write ( *, * ) 'this + that = ', example % other
        write ( *, * ) 'expected value is 3'
    stop
end program demonstration
