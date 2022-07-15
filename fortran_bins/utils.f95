module utils
    implicit none
    
contains
    
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


    function calc_semblance(data, lenx, leny, T, s, dt) result(result_array)
        integer, intent(in) :: lenx, leny
        real, intent(in) :: data(lenx, leny), T, s, dt
        integer :: m, i, reso
        real :: a, p, time, temp, step, result_array, numerador, denominador


    end function calc_semblance

end module utils

module perlin
    implicit none
    INTEGER, DIMENSION(0:511), PARAMETER :: p = (/151,160,137,91,90,15,    &
    131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23, &
    190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33, &
    88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166, &
    77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244, &
    102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196, &
    135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123, &
    5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42, &
    223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9, &
    129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228, &
    251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107, &
    49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254, &
    138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180, &
    151,160,137,91,90,15,    &
    131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23, &
    190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33, &
    88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166, &
    77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244, &
    102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196, &
    135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123, &
    5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42, &
    223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9, &
    129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228, &
    251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107, &
    49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254, &
    138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180/)
        
CONTAINS

    function generate_perlin(lenx, leny) result(array)
        integer, intent(in) :: lenx, leny
        integer :: j, k
        REAL :: r(lenx, leny, 3), array(lenx, leny)

        call random_number(r)

        do concurrent (k = 1: lenx)
            do j = 1, leny
                block
                    array(k, j) = NOISE(r(k, j, 1), r(k, j, 2), r(k, j, 3))
                end block
            end do
        end do
    end function generate_perlin
        
    FUNCTION NOISE(x_in, y_in, z_in) result(random)
        REAL, INTENT(IN) :: x_in, y_in, z_in
        REAL :: x, y, z
        INTEGER :: xx, yy, zz, a, aa, ab, b, ba, bb
        REAL :: u, v, w
        REAL :: random
        
        x = x_in
        y = y_in
        z = z_in
        
        xx = IAND(FLOOR(x), 255)
        yy = IAND(FLOOR(y), 255)
        zz = IAND(FLOOR(z), 255)
        
        x = x - FLOOR(x)
        y = y - FLOOR(y)
        z = z - FLOOR(z)
        
        u = FADE(x)
        v = FADE(y)
        w = FADE(z)
        
        a  = p(xx)   + yy
        aa = p(a)    + zz
        ab = p(a+1)  + zz
        b  = p(xx+1) + yy
        ba = p(b)    + zz
        bb = p(b+1)  + zz
        
        random = LERP(w, LERP(v, LERP(u, GRAD(p(aa),   x,   y,   z),     &
                                        GRAD(p(ba),   x-1, y,   z)),    &
                                LERP(u, GRAD(p(ab),   x,   y-1, z),     &
                                        GRAD(p(bb),   x-1, y-1, z))),   &
                        LERP(v, LERP(u, GRAD(p(aa+1), x,   y,   z-1),   &
                                        GRAD(p(ba+1), x-1, y,   z-1)),  &
                                LERP(u, GRAD(p(ab+1), x,   y-1, z-1),   &
                                        GRAD(p(bb+1), x-1, y-1, z-1))))
    END FUNCTION
        
        
    REAL FUNCTION FADE(t)
        REAL, INTENT(IN) :: t
        
        FADE = t ** 3 * (t * ( t * 6 - 15) + 10)
    END FUNCTION
        
        
    REAL FUNCTION LERP(t, a, b)
        REAL, INTENT(IN) :: t, a, b
        
        LERP = a + t * (b - a)
    END FUNCTION
        
        
    REAL FUNCTION GRAD(hash, x, y, z)
        INTEGER, INTENT(IN) :: hash
        REAL, INTENT(IN) :: x, y, z
        INTEGER :: h
        REAL :: u, v
        
        h = IAND(hash, 15)
        
        u = MERGE(x, y, h < 8)
        
        v = MERGE(y, MERGE(x, z, h == 12 .OR. h == 14), h < 4)
        
        GRAD = MERGE(u, -u, IAND(h, 1) == 0) + MERGE(v, -v, IAND(h, 2) == 0)
    END FUNCTION
    
end module perlin