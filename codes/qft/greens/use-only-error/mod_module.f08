module mModule
    implicit none
    type :: myType
        real :: this, that, other
    contains
        private
        procedure, public :: adder => adder_sub
    end type myType

    private :: adder_sub

    interface
        module subroutine adder_sub ( me )
            class ( myType ), target :: me
        end subroutine adder_sub
    end interface

end module mModule
