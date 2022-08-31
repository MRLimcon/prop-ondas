module calc
    implicit none
    
contains

    function free_wave_equation_2d(lenx, leny, sol_len, c, ew, ew_len, dt) result(array)
        implicit none

        integer, intent(in) :: lenx, leny, sol_len, ew_len
        real, intent(in) :: c(lenx, leny), dt, ew(ew_len)
        real :: array(lenx, leny, sol_len), acceleration(lenx, leny), velocity(lenx, leny)
        integer :: i, half_lenx, half_leny

        array = 0
        velocity = 0
        half_lenx = lenx/2
        half_leny = leny/2
        array(half_lenx, half_leny, 1) = ew(1)

        do i = 2, sol_len, 1
            acceleration = 0

            acceleration(1, :) = ( array(2, :, i-1) - array(1, :, i-1) )
            acceleration(lenx, :) = (array(lenx-1, :, i-1) - array(lenx, :, i-1))
            acceleration(2:lenx-1, :) = array(3:lenx, :, i-1) + array(1:lenx-2, :, i-1) &
                - (2*array(2:lenx-1, :, i-1))

            acceleration(:, 1) = acceleration(:, 1) + ( array(:, 2, i-1) - array(:, 1, i-1) )
            acceleration(:, leny) = acceleration(:, leny) & 
                + (array(:, leny-1, i-1) - array(:, leny, i-1))
            acceleration(:, 2:leny-1) = acceleration(:, 2:leny-1) + array(:, 3:leny, i-1) &
                + array(:, 1:leny-2, i-1) - (2*array(:, 2:leny-1, i-1))

            acceleration = acceleration*c

            velocity = velocity + (acceleration*dt)
            array(:, :, i) = array(:, :, i-1) + (velocity*dt)

            if ( i <= ew_len ) then
                array(half_lenx, half_leny, i) = ew(i)
            end if
            
        end do

    end function free_wave_equation_2d

    function get_wave_acceleration(lenx, leny, c, dx, array) result(acceleration)
        implicit none

        integer, intent(in) :: lenx, leny
        real, intent(in) :: c(lenx, leny-2), array(lenx, leny), dx
        real ::  acceleration(lenx, leny-2)

        acceleration = 0

        acceleration(2:lenx-1, :) = array(3:lenx, 2:leny-1) + array(1:lenx-2, 2:leny-1) &
                - (2*array(2:lenx-1, 2:leny-1))

        acceleration(:, :) = acceleration(:, :) + array(:, 3:leny) &
                + array(:, 1:leny-2) - (2*array(:, 2:leny-1))

        acceleration = acceleration*c/(dx**2)

    end function get_wave_acceleration

    function get_solid_acceleration(lenx, leny, mu, l, dx, solution) result(acceleration)
        implicit none

        integer, intent(in) :: lenx, leny
        real, intent(in) :: mu(lenx, leny),  dx, l(lenx, leny), solution(lenx, leny, 2)
        real :: grad(lenx, leny), laplacian(lenx, leny, 2), derivative(lenx, leny, 2)
        real ::  acceleration(lenx, leny, 2)

        acceleration = 0
        laplacian = 0
        grad = 0
        derivative = 0

        ! u_tt = (l + mu)*(divergent . (u, v))_x + (mu)*(laplacian u)
        ! v_tt = (l + mu)*(divergent . (u, v))_v + (mu)*(laplacian v)
        !taking the laplacian in u

        !laplacian(1, :, 1) = ( solution(2, :, 1) - solution(1, :, 1) )
        !laplacian(lenx, :, 1) = (solution(lenx-1, :, 1) - solution(lenx, :, 1))
        laplacian(2:lenx-1, :, 1) = solution(3:lenx, :, 1) + solution(1:lenx-2, :, 1) &
            - (2*solution(2:lenx-1, :, 1))

        !laplacian(:, 1, 1) = laplacian(:, 1, 1) + ( solution(:, 2, 1) - solution(:, 1, 1) )
        !laplacian(:, leny, 1) = laplacian(:, leny, 1) & 
        !    + (solution(:, leny-1, 1) - solution(:, leny, 1))
        laplacian(:, 2:leny-1, 1) = laplacian(:, 2:leny-1, 1) + solution(:, 3:leny, 1) &
            + solution(:, 1:leny-2, 1) - (2*solution(:, 2:leny-1, 1))

        !taking the laplacian in v

        !laplacian(1, :, 2) = ( solution(2, :, 2) - solution(1, :, 2) )
        !laplacian(lenx, :, 2) = (solution(lenx-1, :, 2) - solution(lenx, :, 2))
        laplacian(2:lenx-1, :, 2) = solution(3:lenx, :, 2) + solution(1:lenx-2, :, 2) &
            - (2*solution(2:lenx-1, :, 2))

        !laplacian(:, 1, 2) = laplacian(:, 1, 2) + ( solution(:, 2, 2) - solution(:, 1, 2) )
        !laplacian(:, leny, 2) = laplacian(:, leny, 2) & 
        !    + (solution(:, leny-1, 2) - solution(:, leny, 2))
        laplacian(:, 2:leny-1, 2) = laplacian(:, 2:leny-1, 2) + solution(:, 3:leny, 2) &
            + solution(:, 1:leny-2, 2) - (2*solution(:, 2:leny-1, 2))

        ! taking the divergence

        !grad(1, :) = ( solution(2, :, 1) - solution(1, :, 1) )
        !grad(lenx, :) = -(solution(lenx-1, :, 1) - solution(lenx, :, 1))
        grad(2:lenx-1, :) = (solution(3:lenx, :, 1) - solution(1:lenx-2, :, 1))/2

        !grad(:, 1) = grad(:, 1) + ( solution(:, 2, 2) - solution(:, 1, 2) )
        !grad(:, leny) = grad(:, leny) & 
        !    - (solution(:, leny-1, 2) - solution(:, leny, 2))
        grad(:, 2:leny-1) = grad(:, 2:leny-1) + (solution(:, 3:leny, 2) &
            - solution(:, 1:leny-2, 2))/2
        
        ! differentiating in x and y, for u and v

        !derivative(1, :, 1) = ( grad(2, :) - grad(1, :) )
        !derivative(lenx, :, 1) = -(grad(lenx-1, :) - grad(lenx, :))
        derivative(2:lenx-1, :, 1) = (grad(3:lenx, :) - grad(1:lenx-2, :))/2

        !derivative(:, 1, 2) = ( grad(:, 2) - grad(:, 1) )
        !derivative(:, leny, 2) = -(grad(:, leny-1) - grad(:, leny))
        derivative(:, 2:leny-1, 2) = (grad(:, 3:leny) - grad(:, 1:leny-2))/2

        ! integrating
        acceleration(:, :, 1) = ((laplacian(:, :, 1)/(dx**2))*mu) + ((derivative(:, :, 1)/(dx))*(mu + l))
        acceleration(:, :, 2) = ((laplacian(:, :, 2)/(dx**2))*mu) + ((derivative(:, :, 2)/(dx))*(mu + l))

    end function get_solid_acceleration

    function elastodynamic_2d(lenx, leny, sol_len, ar_len, ar_steps, ew, ew_len, mu, l, rho, dt, dx, lower, upper) result(array)
        implicit none

        integer, intent(in) :: lenx, leny, sol_len, ew_len, ar_len, ar_steps, lower, upper
        real, intent(in) :: mu(lenx, leny), dt, ew(ew_len, 2), dx
        real, intent(in) :: l(lenx, leny), rho(lenx, leny)
        real :: array(ar_len, lenx, leny, 2), acceleration(lenx, leny, 2), last(lenx, leny, 2)
        real :: solution(lenx, leny, 2), k
        integer :: i, j, half_lenx, half_leny

        solution = 0
        last = 0
        half_lenx = lenx/2
        half_leny = leny/2
        solution(half_lenx, half_leny, :) = ew(1, :)
        array(1, :, :, :) = solution
        j = 2

        do i = 2, sol_len, 1
            acceleration(:, 1:lower, :) = get_solid_acceleration(lenx, lower, mu(:, 1:lower), l(:, 1:lower), &
                    dx, solution(:, 1:lower, :))
            acceleration(:, upper:lenx, :) = get_solid_acceleration(lenx, leny-upper+1, mu(:, upper:lenx), &
                    l(:, upper:lenx), dx, solution(:, upper:lenx, :))
            acceleration(:, lower:upper, 1) = get_wave_acceleration(lenx, upper-lower+3, l(:, lower:upper), &
                    dx, solution(:, lower-1:upper+1, 1))
            acceleration(:, lower:upper, 2) = get_wave_acceleration(lenx, upper-lower+3, l(:, lower:upper), &
                    dx, solution(:, lower-1:upper+1, 2))

            solution = (2*solution) - (last) + ((dt**2)*acceleration)
            last = solution

            if ( i <= ew_len ) then
                solution(half_lenx, half_leny, :) = ew(i, :)
            end if

            if ( mod(i, ar_steps) == 0 .and. j <= ar_len ) then
                k = (100.0*j/ar_len)
                write(*,*) j, "/", ar_len, "-", k, "%"
                array(j, :, :, :) = solution
                j = j + 1
            else if ( j > ar_len ) then
                exit
            end if

        end do

    end function elastodynamic_2d

end module calc