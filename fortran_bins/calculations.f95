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
            array(:, :, i) = array(:, :, i-1) + (velocity*dt) + (acceleration*(dt**2)/2)

            if ( i <= ew_len ) then
                array(half_lenx, half_leny, i) = ew(i)
            end if
            
        end do

    end function free_wave_equation_2d

end module calc

module elastodynamic
    implicit none

    real, allocatable :: wave_accel(:, :)
    real, allocatable :: left_grad(:, :), left_laplacian(:, :, :), left_derivative(:, :, :), left_accel(:, :, :)
    real, allocatable :: right_grad(:, :), right_laplacian(:, :, :), right_derivative(:, :, :), right_accel(:, :, :)

    real, allocatable :: v(:, :, :, :), s(:, :, :, :)
    real, allocatable :: k_runge(:, :, :, :), runge_accel(:, :, :)

contains

    subroutine get_wave_acceleration(w_lenx, w_leny, c, dx, array)
        implicit none

        integer, intent(in) :: w_lenx, w_leny
        real, intent(in) :: c(w_lenx, w_leny-2), array(w_lenx, w_leny), dx

        wave_accel(2:w_lenx-1, :) = array(3:w_lenx, 2:w_leny-1) + array(1:w_lenx-2, 2:w_leny-1) &
                - (2*array(2:w_lenx-1, 2:w_leny-1))

        wave_accel(:, :) = wave_accel(:, :) + array(:, 3:w_leny) &
                + array(:, 1:w_leny-2) - (2*array(:, 2:w_leny-1))

        wave_accel = wave_accel*c/(dx**2)

    end subroutine get_wave_acceleration

    subroutine get_solid_acceleration(lenx, leny, mu, l, dx, solution, sign)
        implicit none

        integer, intent(in) :: lenx, leny
        logical, intent(in) :: sign
        real, intent(in) :: mu(lenx, leny),  dx, l(lenx, leny), solution(lenx, leny, 2)

        ! u_tt = (l + mu)*(divergent . (u, v))_x + (mu)*(laplacian u)
        ! v_tt = (l + mu)*(divergent . (u, v))_v + (mu)*(laplacian v)

        if ( sign ) then
            !taking the laplacian in u
            left_laplacian(2:lenx-1, :, 1) = solution(3:lenx, :, 1) + solution(1:lenx-2, :, 1) &
            - (2*solution(2:lenx-1, :, 1))
            left_laplacian(:, 2:leny-1, 1) = left_laplacian(:, 2:leny-1, 1) + solution(:, 3:leny, 1) &
                + solution(:, 1:leny-2, 1) - (2*solution(:, 2:leny-1, 1))

            !taking the laplacian in v
            left_laplacian(2:lenx-1, :, 2) = solution(3:lenx, :, 2) + solution(1:lenx-2, :, 2) &
                - (2*solution(2:lenx-1, :, 2))
            left_laplacian(:, 2:leny-1, 2) = left_laplacian(:, 2:leny-1, 2) + solution(:, 3:leny, 2) &
                + solution(:, 1:leny-2, 2) - (2*solution(:, 2:leny-1, 2))

            ! taking the divergence
            left_grad(2:lenx-1, :) = (solution(3:lenx, :, 1) - solution(1:lenx-2, :, 1))/2
            left_grad(:, 2:leny-1) = left_grad(:, 2:leny-1) + (solution(:, 3:leny, 2) &
                - solution(:, 1:leny-2, 2))/2
            
            ! differentiating in x and y, for u and v
            left_derivative(2:lenx-1, :, 1) = (left_grad(3:lenx, :) - left_grad(1:lenx-2, :))/2
            left_derivative(:, 2:leny-1, 2) = (left_grad(:, 3:leny) - left_grad(:, 1:leny-2))/2

            ! integrating
            left_accel(:, :, 1) = ((left_laplacian(:, :, 1)/(dx**2))*mu) + ((left_derivative(:, :, 1)/(dx))*(mu + l))
            left_accel(:, :, 2) = ((left_laplacian(:, :, 2)/(dx**2))*mu) + ((left_derivative(:, :, 2)/(dx))*(mu + l))
        else 
            !taking the laplacian in u
            right_laplacian(2:lenx-1, :, 1) = solution(3:lenx, :, 1) + solution(1:lenx-2, :, 1) &
            - (2*solution(2:lenx-1, :, 1))
            right_laplacian(:, 2:leny-1, 1) = right_laplacian(:, 2:leny-1, 1) + solution(:, 3:leny, 1) &
                + solution(:, 1:leny-2, 1) - (2*solution(:, 2:leny-1, 1))

            !taking the laplacian in v
            right_laplacian(2:lenx-1, :, 2) = solution(3:lenx, :, 2) + solution(1:lenx-2, :, 2) &
                - (2*solution(2:lenx-1, :, 2))
            right_laplacian(:, 2:leny-1, 2) = right_laplacian(:, 2:leny-1, 2) + solution(:, 3:leny, 2) &
                + solution(:, 1:leny-2, 2) - (2*solution(:, 2:leny-1, 2))

            ! taking the divergence
            right_grad(2:lenx-1, :) = (solution(3:lenx, :, 1) - solution(1:lenx-2, :, 1))/2
            right_grad(:, 2:leny-1) = right_grad(:, 2:leny-1) + (solution(:, 3:leny, 2) &
                - solution(:, 1:leny-2, 2))/2
            
            ! differentiating in x and y, for u and v
            right_derivative(2:lenx-1, :, 1) = (right_grad(3:lenx, :) - right_grad(1:lenx-2, :))/2
            right_derivative(:, 2:leny-1, 2) = (right_grad(:, 3:leny) - right_grad(:, 1:leny-2))/2

            ! integrating
            right_accel(:, :, 1) = ((right_laplacian(:, :, 1)/(dx**2))*mu) + ((right_derivative(:, :, 1)/(dx))*(mu + l))
            right_accel(:, :, 2) = ((right_laplacian(:, :, 2)/(dx**2))*mu) + ((right_derivative(:, :, 2)/(dx))*(mu + l))
        end if

    end subroutine get_solid_acceleration

    subroutine runge(lenx, leny, mu, l, dt, dx, lower, upper, solution, velocity)
        implicit none

        real, intent(in) :: solution(lenx, leny, 2), l(lenx, leny), mu(lenx, leny), dt, dx
        real, intent(in) :: velocity(lenx, leny, 2)
        real :: effective_dt
        integer, intent(in) :: lenx, leny, lower, upper
        integer :: i, j

        s(1, :, :, :) = solution

        do i = 1, 4

            if ( i == 1 ) then
                j = 1
            else 
                j = i-1
            end if

            if ( i == 2 .or. i == 3 ) then
                effective_dt = dt/2
            else 
                effective_dt = dt
            end if
            
            call get_solid_acceleration(lenx, lower, mu(:, 1:lower), l(:, 1:lower), &
                    dx, s(j, :, 1:lower, :), .True.)
            call get_solid_acceleration(lenx, leny-upper+1, mu(:, upper:lenx), &
                    l(:, upper:lenx), dx, s(j, :, upper:lenx, :), .False.)

            k_runge(i, :, 1:lower, :) = left_accel
            k_runge(i, :, upper:lenx, :) = right_accel

            call get_wave_acceleration(lenx, upper-lower+3, l(:, lower:upper), &
                    dx, s(j, :, lower-1:upper+1, 1))
            k_runge(i, :, lower:upper, 1) = wave_accel

            call get_wave_acceleration(lenx, upper-lower+3, l(:, lower:upper), &
                    dx, s(j, :, lower-1:upper+1, 2))
            k_runge(i, :, lower:upper, 2) = wave_accel

            v(i, :, :, :) = velocity + (k_runge(i, :, :, :)*effective_dt)
            s(i, :, :, :) = solution + (v(i, :, :, :)*effective_dt) + (k_runge(i, :, :, :)*(effective_dt**2)/2)
            k_runge(i, :, :, :) = k_runge(i, :, :, :)*effective_dt
        end do

        runge_accel = (k_runge(1, :, :, :) + (2*k_runge(2, :, :, :)) &
                + (2*k_runge(3, :, :, :)) + k_runge(4, :, :, :))/6

    end subroutine runge

    function elastodynamic_2d(lenx, leny, sol_len, ar_len, ar_steps, ew, ew_len, mu, l, rho, dt, dx, lower, upper) result(array)
        implicit none

        integer, intent(in) :: lenx, leny, sol_len, ew_len, ar_len, ar_steps, lower, upper
        real, intent(in) :: mu(lenx, leny), dt, ew(ew_len, 2), dx
        real, intent(in) :: l(lenx, leny), rho(lenx, leny)
        real :: array(ar_len, lenx, leny), acceleration(lenx, leny, 2), velocity(lenx, leny, 2)
        real :: solution(lenx, leny, 2), k
        integer :: i, j, half_lenx, half_leny

        solution = 0
        velocity = 0
        half_lenx = lenx/2
        half_leny = leny/2
        solution(half_lenx, half_leny, :) = ew(1, :)
        j = 1

        allocate(wave_accel(lenx, upper-lower+1))

        allocate(left_grad(lenx, lower))
        allocate(left_laplacian(lenx, lower, 2))
        allocate(left_derivative(lenx, lower, 2))
        allocate(left_accel(lenx, lower, 2))

        allocate(right_grad(lenx, leny-upper+1))
        allocate(right_laplacian(lenx, leny-upper+1, 2))
        allocate(right_derivative(lenx, leny-upper+1, 2))
        allocate(right_accel(lenx, leny-upper+1, 2))

        allocate(v(4, lenx, leny, 2))
        allocate(s(4, lenx, leny, 2))
        allocate(k_runge(4, lenx, leny, 2))
        allocate(runge_accel(lenx, leny, 2))

        wave_accel = 0

        left_grad = 0
        left_laplacian = 0
        left_derivative = 0
        left_accel = 0

        right_grad = 0
        right_laplacian = 0
        right_derivative = 0
        right_accel = 0

        k_runge = 0
        runge_accel = 0
        s = 0
        v = 0

        do i = 2, sol_len, 1

            call runge(lenx, leny, mu, l, dt, dx, lower, upper, solution, velocity)

            acceleration = runge_accel
            velocity = velocity + (acceleration)
            solution = solution + (velocity*dt) + (acceleration*dt/2)

            if ( i <= ew_len ) then
                solution(half_lenx, half_leny, :) = ew(i, :)
            end if

            if ( mod(i, ar_steps) == 0 .and. j <= ar_len ) then
                k = (100.0*j/ar_len)
                write(*,*) j, "/", ar_len, "-", k, "%"

                array(j, :, :) = acceleration(:, :, 2) - acceleration(:, :, 1)
                j = j + 1
            else if ( j > ar_len ) then
                exit
            end if

        end do

        deallocate(wave_accel)

        deallocate(left_grad)
        deallocate(left_laplacian)
        deallocate(left_derivative)
        deallocate(left_accel)

        deallocate(right_grad)
        deallocate(right_laplacian)
        deallocate(right_derivative)
        deallocate(right_accel)
        
        deallocate(v)
        deallocate(s)
        deallocate(k_runge)
        deallocate(runge_accel)

    end function elastodynamic_2d
    
end module elastodynamic

module electromagnetic
    implicit none

    real, allocatable :: kb_runge(:, :, :, :, :), ke_runge(:, :, :, :, :), curl_obj(:, :, :, :)
    real, allocatable :: sb_runge(:, :, :, :, :), se_runge(:, :, :, :, :)
    real, allocatable :: b_runge(:, :, :, :), e_runge(:, :, :, :)
    real, parameter :: pi = 3.1415926535
    
contains

    subroutine curl(lenx, leny, lenz, obj, dx)
        implicit none

        integer, intent(in) :: lenx, leny, lenz
        real, intent(in) :: obj(lenx, leny, lenz, 3), dx

        curl_obj = 0

        curl_obj(:, 2:leny-1, :, 1) = obj(:, 3:leny, :, 3) - obj(:, 1:leny-2, :, 3)
        curl_obj(:, :, 2:lenz-1, 1) = curl_obj(:, :, 2:lenz-1, 1) - &
                (obj(:, :, 3:lenz, 2) - obj(:, :, 1:lenz-2, 2))

        curl_obj(:, :, 2:lenz-1, 2) = (obj(:, :, 3:lenz, 1) - obj(:, :, 1:lenz-2, 1))
        curl_obj(2:lenx-1, :, :, 2) = curl_obj(2:lenx-1, :, :, 2) - &
                (obj(3:lenx, :, :, 3) - obj(1:lenx-2, :, :, 3))

        curl_obj(2:lenx-1, :, :, 3) = (obj(3:lenx, :, :, 2) - obj(1:lenx-2, :, :, 2))
        curl_obj(:, 2:leny-1, :, 3) = curl_obj(:, 2:leny-1, :, 3) - &
                (obj(:, 3:leny, :, 1) - obj(:, 1:leny-2, :, 1))

        curl_obj = curl_obj/(2*dx)

    end subroutine

    subroutine runge(lenx, leny, lenz, B, E, conductivity, permi, mu, c, dx, dt)
        implicit none

        integer, intent(in) :: lenx, leny, lenz
        real, intent(in) :: B(lenx, leny, lenz, 3), E(lenx, leny, lenz, 3), dx, dt, permi(lenx, leny, lenz)
        real, intent(in) :: conductivity(lenx, leny, lenz, 3), c(lenx, leny, lenz, 3), mu(lenx, leny, lenz, 3)
        real :: effective_dt
        integer :: i, j

        se_runge(1, :, :, :, :) = E
        sb_runge(1, :, :, :, :) = B 

        do i = 1, 4

            if ( i == 1 ) then
                j = 1
            else 
                j = i-1
            end if

            if ( i == 2 .or. i == 3 ) then
                effective_dt = dt/2
            else 
                effective_dt = dt
            end if

            call curl(lenx, leny, lenz, se_runge(j, :, :, :, :), dx)
            kb_runge(i, :, :, :, :) = -c*(curl_obj)
            call curl(lenx, leny, lenz, sb_runge(j, :, :, :, :), dx)
            ke_runge(i, :, :, :, :) = (curl_obj/c) - (4*pi*conductivity*se_runge(j, :, :, :, :)/c)

            se_runge(i, :, :, :, :) = E + (effective_dt*ke_runge(i, :, :, :, :))
            sb_runge(i, :, :, :, :) = B + (effective_dt*kb_runge(i, :, :, :, :))

        end do

        e_runge = (ke_runge(1, :, :, :, :) + (2*ke_runge(2, :, :, :, :)) &
            + (2*ke_runge(3, :, :, :, :)) + ke_runge(4, :, :, :, :))/6
        b_runge = (kb_runge(1, :, :, :, :) + (2*kb_runge(2, :, :, :, :)) &
            + (2*kb_runge(3, :, :, :, :)) + kb_runge(4, :, :, :, :))/6

    end subroutine runge

    function electromagnetic_3d(lenx, leny, lenz, sol_len, ew_len, ew, ew_format, &
                ar_len, ar_steps, conductivity, mu, permi, c, dx, dt) result(array)
        implicit none

        integer, intent(in) :: lenx, leny, lenz, ar_len, ar_steps, ew_len, sol_len
        logical, intent(in) :: ew_format(lenx, leny, lenz, 3)
        real, intent(in) :: ew(ew_len, lenx, leny, lenz, 3), conductivity(lenx, leny, lenz), c(lenx, leny, lenz), dx, dt
        real, intent(in) :: mu(lenx, leny, lenz), permi(lenx, leny, lenz)
        real :: array(ar_len, lenx, leny, lenz, 3), E(lenx, leny, lenz, 3), B(lenx, leny, lenz, 3)
        integer :: i, j, k

        E = 0
        B = 0
        j = 1
        where(ew_format) E = ew(1, :, :, :, :)

        allocate(kb_runge(4, lenx, leny, lenz, 3))
        allocate(ke_runge(4, lenx, leny, lenz, 3))
        allocate(sb_runge(4, lenx, leny, lenz, 3))
        allocate(se_runge(4, lenx, leny, lenz, 3))
        allocate(b_runge(lenx, leny, lenz, 3))
        allocate(e_runge(lenx, leny, lenz, 3))
        allocate(curl_obj(lenx, leny, lenz, 3))

        do i = 2, sol_len
            call runge(lenx, leny, lenz, B, E, conductivity, permi, mu, c, dx, dt)

            B = B + (b_runge*dt)
            E = E + (e_runge*dt)

            if ( i <= ew_len ) then
                where(ew_format) E = ew(i, :, :, :, :)
            end if

            if ( mod(i, ar_steps) == 0 .and. j <= ar_len ) then
                k = (100.0*j/ar_len)
                write(*,*) j, "/", ar_len, "-", k, "%"

                array(j, :, :, :, :) = E
                j = j + 1
            else if ( j > ar_len ) then
                exit
            end if
        end do

        deallocate(kb_runge)
        deallocate(ke_runge)
        deallocate(sb_runge)
        deallocate(se_runge)
        deallocate(b_runge)
        deallocate(e_runge)
        deallocate(curl_obj)

    
    end function electromagnetic_3d
    
end module electromagnetic