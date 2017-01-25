submodule ( mModule ) mSubModule
    implicit none
contains
    module subroutine adder_sub ( me )
        class ( myType ), target :: me
        me % other = me % this + me % that
    end subroutine adder_sub
end submodule mSubModule
