module utils
    implicit none
    
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
        integer :: i, j

        call random_number(r)

        do concurrent (i = 1: lenx) 
            do j = 1, leny
                if ( r(i, j) >= percent ) then
                    array(i, j) = matrix_const
                else 
                    array(i, j) = liquid_const
                end if
            end do            
        end do

    end function make_permeable

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

    function make_circle(lenx, leny, X, Y, center, radius) result(result_array)
        integer, intent(in) :: lenx, leny
        real, intent(in) :: X(lenx, leny), Y(lenx, leny), radius, center(2)
        real :: distances(lenx, leny)
        logical :: result_array(lenx, leny)

        distances = sqrt(((X-center(1))**2) + ((Y-center(2))**2))
        result_array = distances <= radius

    end function make_circle

end module utils
