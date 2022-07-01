module calc
    implicit none
    
contains

    function free_wave_equation_2d(ic_lenx, ic_leny, delta_x, delta_t, sol_len, c, ic, vr) result(array)
        implicit none

        integer, intent(in) :: ic_lenx, ic_leny, sol_len
        real, intent(in) :: ic(ic_lenx, ic_leny), delta_x, delta_t, c, vr
        real :: array(ic_lenx, ic_leny, sol_len), acceleration(ic_lenx, ic_leny)
        real :: constant
        integer :: i, half_ic_len

        constant = ((delta_t**2) * c**2)/(delta_x**2)
        array(:, :, 1) = ic
        half_ic_len = ic_lenx/7

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

            acceleration(1:3*half_ic_len-1, :) = acceleration(1:3*half_ic_len-1, :)*constant
            acceleration(3*half_ic_len:4*half_ic_len, :) = acceleration(3*half_ic_len:4*half_ic_len, :)*vr*constant
            acceleration(4*half_ic_len+1:ic_lenx, :) = acceleration(4*half_ic_len+1:ic_lenx, :)*constant

            if ( i == 2 ) then
                array(:, :, i) = array(:, :, i-1) + (acceleration)
            else
                array(:, :, i) = (2*array(:, :, i-1)) - (array(:, :, i-2)) + (acceleration)
            end if

        end do

    end function free_wave_equation_2d

end module calc