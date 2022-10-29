module utils
    implicit none

    real, parameter :: pi = 3.1415926535
    
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

    function coil(const_h, const_theta, radius, t) result(points)
        real, intent(in) :: const_h, const_theta, radius, t
        real :: points(3)

        points(1) = radius*cos(const_theta*t)
        points(2) = radius*sin(const_theta*t)
        points(3) = const_h*t
    end function coil

    function make_coil_format(lenx, leny, lenz, X, Y, Z, dx, height, num, center, radius, const) result(result_array)
        integer, intent(in) :: lenx, leny, lenz, num
        real, intent(in) :: X(lenx, leny, lenz), Y(lenx, leny, lenz), Z(lenx, leny, lenz), radius, center(3)
        real, intent(in) :: height, const, dx
        real :: distances(lenx, leny, lenz), min, max
        integer :: steps, i
        logical :: result_array(lenx, leny, lenz)

        min = center(3) - (height/2)
        max = min + height
        steps = height/(2*pi*dx)

        do i = 1, steps
            
        end do

        distances = sqrt(((X-center(1))**2) + ((Y-center(2))**2))
        result_array = (distances <= radius) .and. (abs(Z - center(3)) <= (height/2))

    end function make_coil_format

end module utils
