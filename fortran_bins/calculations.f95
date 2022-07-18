module calc
    implicit none
    
contains

    function free_wave_equation_2d(ic_lenx, ic_leny, sol_len, c, ew, ew_len, dt) result(array)
        implicit none

        integer, intent(in) :: ic_lenx, ic_leny, sol_len, ew_len
        real, intent(in) :: c(ic_lenx, ic_leny), dt, ew(ew_len)
        real :: array(ic_lenx, ic_leny, sol_len), acceleration(ic_lenx, ic_leny), velocity(ic_lenx, ic_leny)
        integer :: i, half_lenx, half_leny

        array = 0
        velocity = 0
        half_lenx = ic_lenx/2
        half_leny = ic_leny/2
        array(half_lenx, half_leny, 1) = ew(1)

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

            if ( i <= ew_len ) then
                array(half_lenx, half_leny, i) = ew(i)
            end if
            
        end do

    end function free_wave_equation_2d

    function isotropic_elastodynamic_2d(ic_lenx, ic_leny, sol_len, mu, lambda, rho, ic, dt) result(array)
        implicit none

        integer, intent(in) :: ic_lenx, ic_leny, sol_len
        real, intent(in) :: ic(ic_lenx, ic_leny), mu(ic_lenx, ic_leny), dt
        real, intent(in) :: lambda(ic_lenx, ic_leny), rho(ic_lenx, ic_leny)
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

            acceleration = (acceleration*mu) + (acceleration*(mu + lambda))
            acceleration = acceleration/rho

            velocity = velocity + (acceleration*dt)
            array(:, :, i) = array(:, :, i-1) + (velocity*dt)

        end do

    end function isotropic_elastodynamic_2d

    function anisotropic_elastodynamic_2d(ic_lenx, ic_leny, sol_len, c, ic, rho, dt) result(array)
        implicit none

        integer, intent(in) :: ic_lenx, ic_leny, sol_len
        real, intent(in) :: ic(ic_lenx, ic_leny), c(ic_lenx, ic_leny, 2, 2), rho(ic_lenx, ic_leny), dt
        real :: array(ic_lenx, ic_leny, sol_len), local_acceleration(ic_lenx, ic_leny, 2)
        real :: acceleration(ic_lenx, ic_leny), final_acceleration(ic_lenx, ic_leny), velocity(ic_lenx, ic_leny)
        integer :: i, j, k

        array(:, :, 1) = ic
        velocity = 0

        do i = 2, sol_len, 1
            acceleration = 0.
            final_acceleration = 0.
            local_acceleration = 0.

            local_acceleration(1, :, 1) = ( array(2, :, i-1) - array(1, :, i-1) )
            local_acceleration(ic_lenx, :, 1) = -(array(ic_lenx-1, :, i-1) - array(ic_lenx, :, i-1))
            local_acceleration(2:ic_lenx-1, :, 1) = (array(3:ic_lenx, :, i-1) - array(1:ic_lenx-2, :, i-1))/2 

            local_acceleration(:, 1, 2) = ( array(:, 2, i-1) - array(:, 1, i-1) )
            local_acceleration(:, ic_leny, 2) = -(array(:, ic_leny-1, i-1) - array(:, ic_leny, i-1))
            local_acceleration(:, 2:ic_leny-1, 2) = (array(:, 3:ic_leny, i-1) - array(:, 1:ic_leny-2, i-1))/2

            do concurrent (j = 1: ic_lenx)
                do k = 1, ic_leny
                    acceleration(j, k) = sum(local_acceleration(j, k, 1)*c(j, k, 1, :)) + &
                        sum(local_acceleration(j, k, 2)*c(j, k, 2, :))
                end do
            end do

            final_acceleration(1, :) = ( acceleration(2, :) - acceleration(1, :) )
            final_acceleration(ic_lenx, :) = -(acceleration(ic_lenx-1, :) - acceleration(ic_lenx, :))
            final_acceleration(2:ic_lenx-1, :) = (acceleration(3:ic_lenx, :) - acceleration(1:ic_lenx-2, :))/2 

            final_acceleration(:, 1) = final_acceleration(:, 1) + &
                ( acceleration(:, 2) - acceleration(:, 1) )
            final_acceleration(:, ic_leny) = final_acceleration(:, ic_leny) - &
                (acceleration(:, ic_leny-1) - acceleration(:, ic_leny))
            final_acceleration(:, 2:ic_leny-1) = final_acceleration(:, 2:ic_leny-1) + &
                (acceleration(:, 3:ic_leny) - acceleration(:, 1:ic_leny-2))/2

            velocity = velocity + (final_acceleration*dt/rho)
            array(:, :, i) = array(:, :, i-1) + (velocity*dt)

        end do

    end function anisotropic_elastodynamic_2d

end module calc