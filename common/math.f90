module math_mod
    implicit none
    real, parameter :: pi = 3.14159265358979323846264338327950288419716939937510
    real, parameter :: pi2 = pi*0.5
    real, parameter :: inf = huge(0.)
contains

    function cross_product(v1, v2) result(cross)
        implicit none
        
        real, dimension(3), intent(in) :: v1, v2
        real, dimension(3) :: cross

        ! Calculate Cross Product
        cross(1) = v1(2) * v2(3) - v1(3) * v2(2)
        cross(2) = v1(3) * v2(1) - v1(1) * v2(3)
        cross(3) = v1(1) * v2(2) - v1(2) * v2(1)
    end function cross_product

    
end module math_mod