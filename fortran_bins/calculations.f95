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

    function elastodynamic_2d(ic_lenx, ic_leny, sol_len, ew, ew_len, mu, lambda, rho, ic, dt) result(array)
        implicit none

        integer, intent(in) :: ic_lenx, ic_leny, sol_len, ew_len
        real, intent(in) :: ic(ic_lenx, ic_leny, 2), mu(ic_lenx, ic_leny), dt, ew(ew_len, 2)
        real, intent(in) :: lambda(ic_lenx, ic_leny), rho(ic_lenx, ic_leny)
        real :: grad(ic_lenx, ic_leny), laplacian(ic_lenx, ic_leny), derivative(ic_lenx, ic_leny, 2)
        real :: array(sol_len, ic_lenx, ic_leny, 2), acceleration(ic_lenx, ic_leny, 2), velocity(ic_lenx, ic_leny, 2)
        integer :: i, half_lenx, half_leny

        array(1, :, :, :) = ic
        velocity = 0
        half_lenx = ic_lenx/2
        half_leny = ic_leny/2
        array(1, half_lenx, half_leny, :) = ew(1, :)

        do i = 2, sol_len, 1
            acceleration = 0
            laplacian = 0
            grad = 0
            derivative = 0

            !taking the laplacian in u
            laplacian(1, :) = ( array(i-1, 2, :, 1) - array(i-1, 1, :, 1) )
            laplacian(ic_lenx, :) = (array(i-1, ic_lenx-1, :, 1) - array(i-1, ic_lenx, :, 1))
            laplacian(2:ic_lenx-1, :) = array(i-1, 3:ic_lenx, :, 1) + array(i-1, 1:ic_lenx-2, :, 1) &
                - (2*array(i-1, 2:ic_lenx-1, :, 1))

            laplacian(:, 1) = laplacian(:, 1) + ( array(i-1, :, 2, 1) - array(i-1, :, 1, 1) )
            laplacian(:, ic_leny) = laplacian(:, ic_leny) & 
                + (array(i-1, :, ic_leny-1, 1) - array(i-1, :, ic_leny, 1))
            laplacian(:, 2:ic_leny-1) = laplacian(:, 2:ic_leny-1) + array(i-1, :, 3:ic_leny, 1) &
                + array(i-1, :, 1:ic_leny-2, 1) - (2*array(i-1, :, 2:ic_leny-1, 1))

            acceleration(:, :, 1) = mu*laplacian
            laplacian = 0

            !taking the laplacian in v
            laplacian(1, :) = ( array(i-1, 2, :, 2) - array(i-1, 1, :, 2) )
            laplacian(ic_lenx, :) = (array(i-1, ic_lenx-1, :, 2) - array(i-1, ic_lenx, :, 2))
            laplacian(2:ic_lenx-1, :) = array(i-1, 3:ic_lenx, :, 2) + array(i-1, 1:ic_lenx-2, :, 2) &
                - (2*array(i-1, 2:ic_lenx-1, :, 1))

            laplacian(:, 1) = laplacian(:, 1) + ( array(i-1, :, 2, 2) - array(i-1, :, 1, 2) )
            laplacian(:, ic_leny) = laplacian(:, ic_leny) & 
                + (array(i-1, :, ic_leny-1, 2) - array(i-1, :, ic_leny, 2))
            laplacian(:, 2:ic_leny-1) = laplacian(:, 2:ic_leny-1) + array(i-1, :, 3:ic_leny, 2) &
                + array(i-1, :, 1:ic_leny-2, 2) - (2*array(i-1, :, 2:ic_leny-1, 2))

            acceleration(:, :, 2) = mu*laplacian

            ! taking the divergence
            grad(1, :) = ( array(i-1, 2, :, 1) - array(i-1, 1, :, 1) )
            grad(ic_lenx, :) = -(array(i-1, ic_lenx-1, :, 1) - array(i-1, ic_lenx, :, 1))
            grad(2:ic_lenx-1, :) = (array(i-1, 3:ic_lenx, :, 1) - array(i-1, 1:ic_lenx-2, :, 1))/2

            grad(:, 1) = grad(:, 1) + ( array(i-1, :, 2, 2) - array(i-1, :, 1, 2) )
            grad(:, ic_leny) = grad(:, ic_leny) & 
                - (array(i-1, :, ic_leny-1, 2) - array(i-1, :, ic_leny, 2))
            grad(:, 2:ic_leny-1) = grad(:, 2:ic_leny-1) + (array(i-1, :, 3:ic_leny, 2) &
                - array(i-1, :, 1:ic_leny-2, 2))/2
            
            ! differentiating in x and y, for u and v
            derivative(1, :, 1) = ( grad(2, :) - grad(1, :) )
            derivative(ic_lenx, :, 1) = -(grad(ic_lenx-1, :) - grad(ic_lenx, :))
            derivative(2:ic_lenx-1, :, 1) = (grad(3:ic_lenx, :) - grad(1:ic_lenx-2, :))/2

            derivative(:, 1, 2) = ( grad(:, 2) - grad(:, 1) )
            derivative(:, ic_leny, 2) = -(grad(:, ic_leny-1) - grad(:, ic_leny))
            derivative(:, 2:ic_leny-1, 2) = (grad(:, 3:ic_leny) - grad(:, 1:ic_leny-2))/2

            ! integrating
            acceleration(:, :, 1) = acceleration(:, :, 1) + (derivative(:, :, 1)*(mu + lambda))
            acceleration(:, :, 2) = acceleration(:, :, 2) + (derivative(:, :, 2)*(mu + lambda))

            velocity = velocity + (acceleration*dt/rho)
            array(i, :, :, :) = array(i-1, :, :, :) + (velocity*dt)

            if ( i <= ew_len ) then
                array(i, half_lenx, half_leny, :) = ew(i, :)
            end if

        end do

    end function elastodynamic_2d

end module calc