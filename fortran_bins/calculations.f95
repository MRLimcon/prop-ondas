module calc
    implicit none
    
contains

    function free_wave_equation_2d(ic_lenx, ic_leny, sol_len, c, ic, dt) result(array)
        implicit none

        integer, intent(in) :: ic_lenx, ic_leny, sol_len
        real, intent(in) :: ic(ic_lenx, ic_leny), c(ic_lenx, ic_leny), dt
        real :: array(ic_lenx, ic_leny, sol_len), acceleration(ic_lenx, ic_leny), velocity(ic_lenx, ic_leny)
        integer :: i

        array(:, :, 1) = ic
        velocity = 0

        do i = 2, sol_len, 1
            acceleration = 0

            acceleration(1, :) = ( array(2, :, i-1) - array(1, :, i-1) )
            acceleration(ic_lenx, :) = (array(ic_lenx-1, :, i-1) - array(ic_lenx, :, i-1))
            acceleration(2:ic_lenx-1, :) = array(3:ic_lenx, :, i-1) + array(1:ic_lenx-2, :, i-1) &
                - (2*array(2:ic_lenx-1, :, i-1))

            acceleration(:, 1) = acceleration(:, 1) + ( array(:, 2, i-1) - array(:, 1, i-1) )
            acceleration(:, ic_leny) = acceleration(:, ic_leny) & 
                + (array(:, ic_leny-1, i-1) - array(:, ic_leny, i-1))
            acceleration(:, 2:ic_leny-1) = acceleration(:, 2:ic_leny-1) + array(:, 3:ic_leny, i-1) &
                + array(:, 1:ic_leny-2, i-1) - (2*array(:, 2:ic_leny-1, i-1))

            acceleration = acceleration*c

            velocity = velocity + (acceleration*dt)
            array(:, :, i) = array(:, :, i-1) + (velocity*dt)

        end do

    end function free_wave_equation_2d

end module calc