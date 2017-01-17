module mReplicate
    implicit none
contains

    subroutine duplicate ( anchors, from, to, length, shift )
        ! passed
        integer, intent ( inout ), dimension ( 1 : 1000 ) :: anchors
        integer, intent ( in )                            :: from, to, length, shift
        ! local
        integer :: k

            do k = 0, length - 1
                anchors ( to + k ) = anchors ( from + k ) + shift
            end do

    end subroutine duplicate

    subroutine replicate ( anchors, length, position, copies, shift )
        ! passed
        integer, intent ( inout ), dimension ( 1 : 1000 ) :: anchors
        integer, intent ( inout ) :: position
        integer, intent ( in )    :: length, copies, shift
        ! local
        integer :: replication, l

            do replication = 1, copies ! replications
                do l = 1, length ! sweep subset
                    position = position + 1
                    anchors ( position ) = anchors ( position - length ) + shift
                end do
            end do

    end subroutine replicate

end module mReplicate
