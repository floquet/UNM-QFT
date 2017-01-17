module mStride
    implicit none

    integer, parameter :: nDim = 3

contains

    subroutine stride_up ( shape, strides, nUnits )
        ! passed: rank 1
        integer, dimension ( 1 : nDim ), intent ( in )  :: shape
        ! passed: rank 0
        integer, dimension ( 1 : ndim ), intent ( out ) :: strides
        integer,                         intent ( out ) :: nUnits
        ! local
        integer :: j = 0, k = 0

            nUnits = product ( shape )
            strides ( : ) = nUnits

            do k = 1, nDim
                strides ( k ) = nUnits / product ( [ ( shape ( j ), j = k, nDim ) ] )
            end do

    end subroutine stride_up

    subroutine stride_dn ( shape, strides, nUnits )
        ! passed: rank 1
        integer, dimension ( 1 : nDim ), intent ( in )  :: shape
        ! passed: rank 0
        integer, dimension ( 1 : ndim ), intent ( out ) :: strides
        integer,                         intent ( out ) :: nUnits
        ! local
        integer :: j = 0, k = 0

            nUnits = product ( shape )
            strides ( : ) = nUnits

            do k = 1, nDim
                strides ( k ) = nUnits / product ( [ ( shape ( j ), j = nDim - k + 1, nDim ) ] )
            end do

    end subroutine stride_dn

end module mStride
