module utils
    implicit none

    real, parameter :: pi = 3.1415926535, shift(2) = (/0.25, -0.25/)
    integer, parameter :: turn_const = 50
    
contains
    

    function fix_ew(ew_len, ew) result(array)
        integer, intent(in) :: ew_len
        real, intent(in) :: ew(ew_len)
        real :: max_val
        integer :: array(2)
        integer :: i, init, end_int

        max_val = maxval(abs(ew))*0.02

        do i = 1, ew_len
            if ( all(abs(ew(1:i)) <= max_val) ) then
                init = i
            else if ( all(abs(ew(i:ew_len)) <= max_val) ) then
                end_int = i
                exit
            end if
        end do

        array(1) = init 
        array(2) = end_int

    end function fix_ew

    function make_permeable(lenx, leny, matrix_const, liquid_const, percent) result(array)
        integer, intent(in) :: lenx, leny
        real, intent(in) :: matrix_const, liquid_const, percent
        real :: r(lenx, leny), array(lenx, leny)
        logical :: check(lenx, leny)

        call random_number(r)

        check = r >= percent
        where (check) array = matrix_const
        where (.not. check) array = liquid_const

    end function make_permeable

    function make_3d_permeable(lenx, leny, lenz, matrix_const, liquid_const, percent) result(array)
        integer, intent(in) :: lenx, leny, lenz
        real, intent(in) :: matrix_const, liquid_const, percent
        real :: r(lenx, leny, lenz), array(lenx, leny, lenz)
        logical :: check(lenx, leny, lenz)

        call random_number(r)

        check = r >= percent
        where (check) array = matrix_const
        where (.not. check) array = liquid_const

    end function make_3d_permeable

    function make_stripes(lenx, leny, matrix_const, liquid_const, height, layers) result(array)
        integer, intent(in) :: lenx, leny, layers, height
        real, intent(in) :: matrix_const, liquid_const
        real :: array(lenx, 0:leny-1)
        integer :: i

        do concurrent (i = 0: layers-1) 
            block
                integer :: min, max

                min = i*height
                max = (i+1)*height

                if (max > leny-1) then
                    max = leny-1
                end if 

                if (mod(i, 2) == 0) then
                    array(:, min: max) = matrix_const
                else 
                    array(:, min: max) = liquid_const
                end if

            end block

        end do

    end function make_stripes

    function make_3d_stripes(lenx, leny, lenz, matrix_const, liquid_const, height, layers) result(array)
        integer, intent(in) :: lenx, leny, lenz, layers, height
        real, intent(in) :: matrix_const, liquid_const
        real :: array(lenx, leny, 0:lenz-1)
        integer :: i

        do concurrent (i = 0: layers-1) 
            block
                integer :: min, max

                min = i*height
                max = (i+1)*height

                if (max > lenz-1) then
                    max = lenz-1
                end if 

                if (mod(i, 2) == 0) then
                    array(:, :, min: max) = matrix_const
                else 
                    array(:, :, min: max) = liquid_const
                end if

            end block

        end do

    end function make_3d_stripes

    function make_circle(lenx, leny, X, Y, center, radius) result(result_array)
        integer, intent(in) :: lenx, leny
        real, intent(in) :: X(lenx, leny), Y(lenx, leny), radius, center(2)
        real :: distances(lenx, leny)
        logical :: result_array(lenx, leny)

        distances = sqrt(((X-center(1))**2) + ((Y-center(2))**2))
        result_array = distances <= radius

    end function make_circle

    function make_sphere(lenx, leny, lenz, X, Y, Z, center, radius) result(result_array)
        integer, intent(in) :: lenx, leny, lenz
        real, intent(in) :: X(lenx, leny, lenz), Y(lenx, leny, lenz), Z(lenx, leny, lenz), radius, center(3)
        real :: distances(lenx, leny, lenz)
        logical :: result_array(lenx, leny, lenz)

        distances = sqrt(((X-center(1))**2) + ((Y-center(2))**2) + ((Z-center(3))**2))
        result_array = distances <= radius

    end function make_sphere

    function make_cilinder(lenx, leny, lenz, X, Y, Z, height, center, radius) result(result_array)
        integer, intent(in) :: lenx, leny, lenz
        real, intent(in) :: X(lenx, leny, lenz), Y(lenx, leny, lenz), Z(lenx, leny, lenz), radius, center(3)
        real, intent(in) :: height
        real :: distances(lenx, leny, lenz)
        logical :: result_array(lenx, leny, lenz)

        distances = sqrt(((X-center(1))**2) + ((Y-center(2))**2))
        result_array = (distances <= radius) .and. (abs(Z - center(3)) <= (height/2))

    end function make_cilinder

    function ring(radius, declination, t) result(points)
        real, intent(in) :: radius, declination, t
        real :: points(3)

        points(1) = radius*cos(t)
        points(2) = radius*sin(t)*cos(declination)
        points(3) = radius*sin(t)*sin(declination)

    end function ring

    function ring_derivative(radius, declination, t) result(points)
        real, intent(in) :: radius, declination, t
        real :: points(3)

        points(1) = -radius*sin(t)
        points(2) = radius*cos(t)*cos(declination)
        points(3) = radius*cos(t)*sin(declination)
    end function ring_derivative

    subroutine make_ring_coil(lenx, leny, lenz, X, Y, Z, dx, center, radius, radius_b, &
                declination, coil_format, coil_derivative)
        integer, intent(in) :: lenx, leny, lenz
        real, intent(in) :: X(lenx, leny, lenz), Y(lenx, leny, lenz), Z(lenx, leny, lenz)
        real, intent(in) :: dx, radius, radius_b, center(3), declination
        real, intent(inout) :: coil_derivative(lenx, leny, lenz, 3)
        real :: distances(lenx, leny, lenz), f_distances(lenx, leny, lenz), l_distances(lenx, leny, lenz)
        real :: point(3), f_point(3), l_point(3), derivative(3), vec_size, vec_sizes(lenx, leny, lenz)
        integer :: steps, i, j, k, l
        logical :: local_points(lenx, leny, lenz)
        logical, intent(inout) :: coil_format(lenx, leny, lenz, 3)

        steps = (2*pi)/dx
        coil_format = .False.
        coil_derivative = 0

        do j = 1, 2
            do k = 1, 2
                do l = 1, 2

                    do i = 0, steps

                        if ( i == 0 ) then
                            point = ring(radius, declination, i*dx) + center
                            f_point = ring(radius, declination, (i*dx)+dx) + center
                            l_point = ring(radius, declination, (i*dx)-dx) + center
            
                            distances = sqrt(((X-point(1) + (shift(j)*dx))**2) + ((Y-point(2) + (shift(k)*dx))**2) &
                                + ((Z-point(3) + (shift(l)*dx))**2))
                            f_distances = sqrt(((X-f_point(1) + (shift(j)*dx))**2) + ((Y-f_point(2) + (shift(k)*dx))**2) &
                                + ((Z-f_point(3) + (shift(l)*dx))**2))
                            l_distances = sqrt(((X-l_point(1) + (shift(j)*dx))**2) + ((Y-l_point(2) + (shift(k)*dx))**2) &
                                + ((Z-l_point(3) + (shift(l)*dx))**2))
                        else
                            f_point = ring(radius, declination, (i*dx)+dx) + center
                            f_distances = sqrt(((X-f_point(1) + (shift(j)*dx))**2) + ((Y-f_point(2) + (shift(k)*dx))**2) &
                                + ((Z-f_point(3) + (shift(l)*dx))**2))
                        end if
            
                        local_points = distances <= radius_b .and. distances < f_distances .and. distances < l_distances
            
                        derivative = ring_derivative(radius, declination, i*dx)
                        vec_size = sqrt((derivative(1)**2) + (derivative(2)**2) + (derivative(3)**2))
                        derivative = derivative/vec_size
            
                        where(local_points) coil_derivative(:, :, :, 1) = coil_derivative(:, :, :, 1) + derivative(1)
                        where(local_points) coil_derivative(:, :, :, 2) = coil_derivative(:, :, :, 2) + derivative(2)
                        where(local_points) coil_derivative(:, :, :, 3) = coil_derivative(:, :, :, 3) + derivative(3)
            
                        coil_format(:, :, :, 1) = coil_format(:, :, :, 1) .or. local_points
                        coil_format(:, :, :, 2) = coil_format(:, :, :, 2) .or. local_points
                        coil_format(:, :, :, 3) = coil_format(:, :, :, 3) .or. local_points
            
                        l_point = point
                        point = f_point
                        l_distances = distances
                        distances = f_distances

                    end do

                end do
            end do
        end do

        vec_sizes = sqrt((coil_derivative(:, :, :, 1)**2) + (coil_derivative(:, :, :, 2)**2) &
            + (coil_derivative(:, :, :, 3)**2))

        where(vec_sizes > 8) coil_derivative(:, :, :, 1) = coil_derivative(:, :, :, 1)*8/vec_sizes
        where(vec_sizes > 8) coil_derivative(:, :, :, 2) = coil_derivative(:, :, :, 2)*8/vec_sizes
        where(vec_sizes > 8) coil_derivative(:, :, :, 3) = coil_derivative(:, :, :, 3)*8/vec_sizes

        coil_derivative = coil_derivative*(1.0/8.0)

    end subroutine make_ring_coil

    function make_ring(lenx, leny, lenz, X, Y, Z, dx, center, radius, radius_b, &
        declination) result(ring_format)
        integer, intent(in) :: lenx, leny, lenz
        real, intent(in) :: X(lenx, leny, lenz), Y(lenx, leny, lenz), Z(lenx, leny, lenz)
        real, intent(in) :: dx, radius, radius_b, center(3), declination
        integer :: steps, i
        logical :: ring_format(lenx, leny, lenz)

        steps = (2*pi)/dx
        ring_format = .False.

        do concurrent (i = 0: steps)
            block
                logical :: local_points(lenx, leny, lenz)
                real :: point(3), distances(lenx, leny, lenz)

                point = ring(radius, declination, i*dx) + center

                distances = sqrt(((X-point(1))**2) + ((Y-point(2))**2) + ((Z-point(3))**2))

                local_points = distances <= radius_b

                ring_format = ring_format .or. local_points
            end block
        end do

    end function make_ring

    function make_array(lenx, leny, lenz) result(array)
        integer, intent(in) :: lenx, leny, lenz
        real :: array(lenx, leny, lenz, 3)

        array = 0
    end function make_array

    function make_logical_array(lenx, leny, lenz) result(array)
        integer, intent(in) :: lenx, leny, lenz
        logical :: array(lenx, leny, lenz, 3)

        array = .False.
    end function make_logical_array

    function get_coil_response(lenx, leny, lenz, steps, derivative, format, array) result(result_array)
        integer, intent(in) :: lenx, leny, lenz, steps
        real, intent(in) :: derivative(lenx, leny, lenz, 3), array(steps, lenx, leny, lenz, 3)
        logical, intent(in) :: format(lenx, leny, lenz, 3)
        real :: result_array(steps)
        integer :: i

        do concurrent (i = 1: steps)
            result_array(i) = sum(pack(derivative*array(i, :, :, :, :), format))           
        end do

    end function get_coil_response

end module utils
