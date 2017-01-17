program demo
    implicit none
    real :: x = 2.0
        write ( * , 100 ) x, x ** ( 3 / 2 )
 100    format( g0, ' ** ( 3 / 2 ) = ', g0 )
    stop
end program demo
