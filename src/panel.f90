! Subroutines for working with an individual panel
module panel_mod

    use base_geom_mod
    implicit none

    
    type panel

        real, dimension(3) :: normal, dC_f, centr, centr_pro
        integer, dimension(3) :: vert_ind
        real :: c_p=0, norm_mag, area, m_surf, theta, phi, shadowed=0
        real :: x_min_pro, x_max_pro, y_min_pro, y_max_pro, z_min_pro, z_max_pro ! Outline of projected panel
        real, dimension(3) :: s12_pro, s23_pro, s31_pro ! Projected side vectors
        integer :: N=3
        type(vertex_pointer), dimension(:), allocatable :: vertices
        real, dimension(:), allocatable :: panel_shadowing
        logical :: seperated=.false.

        contains

        procedure :: calc_norm => panel_calc_norm
        procedure :: calc_centroid => panel_calc_centroid
        procedure :: calc_projected_outline => panel_calc_projected_outline
        procedure :: calc_pressure_newton => panel_calc_pressure_newton
        procedure :: calc_pressure_prandtl => panel_calc_pressure_prandtl
        procedure :: get_vertex_loc => panel_get_vertex_loc



    end type panel

contains

    subroutine panel_calc_norm(this, vertices)
        ! Calculates the normal vector of a panel
        implicit none
        class(panel), intent(inout) :: this
        type(vertex), dimension(:), allocatable, intent(in) :: vertices
        
        real, dimension(3) :: v1, v2, cross
        
        ! Calculate panel vectors
        v1 = vertices(this%vert_ind(1))%location - vertices(this%vert_ind(2))%location
        v2 = vertices(this%vert_ind(1))%location - vertices(this%vert_ind(3))%location
        
        ! Calculate Cross Product
        cross(1) = v1(2) * v2(3) - v1(3) * v2(2)
        cross(2) = v1(3) * v2(1) - v1(1) * v2(3)
        cross(3) = v1(1) * v2(2) - v1(2) * v2(1)
        
        ! Calculate the magnitude of the normal
        this%norm_mag = sqrt((cross(1)**2 + cross(2)**2 + cross(3)**2))

        this%normal = cross / this%norm_mag

        

    end subroutine panel_calc_norm

    subroutine panel_calc_centroid(this, vertices)
    
        implicit none
        
        class(panel), intent(inout) :: this
        type(vertex), dimension(:), allocatable, intent(in) :: vertices
        real, dimension(3) :: sum
        integer :: i

        ! Get average of corner points
        sum = 0.
        do i=1,this%N
            sum = sum + vertices(this%vert_ind(i))%location
        end do

        ! Set centroid
        this%centr = sum/this%N

    end subroutine panel_calc_centroid

    subroutine panel_calc_pressure_newton(this, m, c_pmax)
        ! Calculates the pressure on a panel
        real, intent(in) :: m, c_pmax
        class(panel), intent(inout) :: this

               ! Prandtl-meyer expansion on shadow region
                !if (.not. this%seperated) then 
                    !this%c_p = - ((gamma + 1)/2 * this%theta**2 * (sqrt(1 + ( 4 / ((gamma + 1) * m * this%theta))**2) - 1))
           
        ! Newtonian Method on impact region
        this%c_p = c_pmax*sin(this%theta)**2
        this%m_surf = m * cos(this%theta)
    
    end subroutine panel_calc_pressure_newton

    subroutine panel_calc_pressure_prandtl(this, gamma, m)
        
        class(panel), intent(inout) :: this
        real, intent(in) :: gamma, m

        ! Ensure that 
        this%c_p = - ((gamma + 1)/2 * this%theta**2 * (sqrt(1 + ( 4 / ((gamma + 1) * m * this%theta))**2) - 1))

    end subroutine panel_calc_pressure_prandtl

    subroutine panel_calc_projected_outline(this, projected_verts)

        implicit none
        
        class(panel), intent(inout) :: this
        type(vertex), dimension(:), allocatable, intent(in) :: projected_verts

        real, dimension(3) :: x_s, y_s, z_s

        x_s = [projected_verts(this%vert_ind(1))%location(1), &
                projected_verts(this%vert_ind(2))%location(1), &
                projected_verts(this%vert_ind(3))%location(1)]

        this%x_min_pro = min(x_s(1),x_s(2),x_s(3))
        this%x_max_pro = max(x_s(1),x_s(2),x_s(3))

        y_s = [projected_verts(this%vert_ind(1))%location(2), &
                projected_verts(this%vert_ind(2))%location(2), &
                projected_verts(this%vert_ind(3))%location(2)]

        this%y_min_pro = min(y_s(1),y_s(2),y_s(3))
        this%y_max_pro = max(y_s(1),y_s(2),y_s(3))

        z_s = [projected_verts(this%vert_ind(1))%location(3), &
                projected_verts(this%vert_ind(2))%location(3), &
                projected_verts(this%vert_ind(3))%location(3)]

        this%z_min_pro = min(z_s(1),z_s(2),z_s(3))
        this%z_max_pro = max(z_s(1),z_s(2),z_s(3))

    end subroutine panel_calc_projected_outline

    function panel_get_vertex_loc(this, i) result(loc)

        implicit none
        class(panel), intent(in) :: this
        integer, intent(in) :: i
        real, dimension(3) :: loc

        loc = this%vertices(i)%ptr%location

    end function panel_get_vertex_loc


end module panel_mod
