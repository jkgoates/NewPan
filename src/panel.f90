! Subroutines for working with an individual panel
module panel_mod

    use base_geom_mod
    use math_mod
    use linalg_mod
    implicit none

    
    type panel

        real, dimension(3) :: normal, dC_f, centr, centr_pro, u
        real :: c_p=0, area, m_surf, theta, phi, shadowed=0
        real :: x_min_pro, x_max_pro, y_min_pro, y_max_pro, z_min_pro, z_max_pro ! Outline of projected panel
        real, dimension(3) :: s12_pro, s23_pro, s31_pro ! Projected side vectors
        integer :: index ! Panel index in mesh array
        integer :: N ! Number of sides
        type(vertex_pointer), dimension(:), allocatable :: vertices
        integer, dimension(3) :: abutting_panels ! Indices of panels abutting this one
        integer, dimension(3) :: edges ! Indices of the edges of this panel
        real, dimension(:), allocatable :: panel_shadowing
        logical :: seperated=.false.
        real :: a_1, a_2, a_3, a_4=1.0 ! Surface Equation Coefficients F(x,y,z) = a_1*x + a_2*y + a_3*z + a_4 = 0

        contains

        ! Panel Initialization
        procedure :: init => panel_init
        procedure :: calc_derived_geom => panel_calc_derived_geom

        ! Geometry calculations
        procedure :: calc_norm => panel_calc_norm
        procedure :: calc_centroid => panel_calc_centroid
        procedure :: calc_area => panel_calc_area
        procedure :: calc_velocity_vector => panel_calc_velocity_vector
        procedure :: calc_projected_outline => panel_calc_projected_outline

        ! Pressure calculations
        procedure :: calc_pressure_newton => panel_calc_pressure_newton
        procedure :: calc_pressure_prandtl => panel_calc_pressure_prandtl

        ! Getters
        procedure :: get_vertex_loc => panel_get_vertex_loc
        procedure :: get_vertex_index => panel_get_vertex_index


        ! Property equations
        procedure :: calc_surface_equation => panel_calc_surface_equation



    end type panel

contains

    subroutine panel_init(this, v1, v2, v3, index)
        implicit none
        
        class(panel), intent(inout) :: this
        type(vertex), intent(inout), target :: v1, v2, v3
        integer, intent(in) :: index

        integer :: i

        ! Set number of sides
        this%N = 3

        ! Allocate vertex arrays
        allocate(this%vertices(this%N))

        ! Assign vertex pointers
        this%vertices(1)%ptr => v1
        this%vertices(2)%ptr => v2
        this%vertices(3)%ptr => v3
        
        ! Store the index of the panel
        this%index = index

        ! Store that this panel is attached to its vertices
        do i = 1, this%N
            call this%vertices(i)%ptr%panels%append(this%index)
        end do

        ! Calculate panel geometries only dependent on vertex locations
        call this%calc_derived_geom()

        ! Initialize a few other things
        this%abutting_panels = 0

    end subroutine panel_init

    subroutine panel_calc_derived_geom(this)
        ! Initializes geometry based on the location of vertices

        implicit none
        
        class(panel), intent(inout) :: this

        ! Calculate Normal Vector
        call this%calc_norm()

        ! Calculate area
        call this%calc_area()

        ! Calculate centroid
        call this%calc_centroid()

    end subroutine panel_calc_derived_geom

    subroutine panel_calc_norm(this)
        ! Calculates the normal vector of a panel
        implicit none
        class(panel), intent(inout) :: this
        
        real, dimension(3) :: d1, d2, cross
        
        ! Calculate panel vectors
        d1 = this%get_vertex_loc(2) - this%get_vertex_loc(1)
        d2 = this%get_vertex_loc(3) - this%get_vertex_loc(2)
        
        ! Calculate Cross Product
        cross = cross_product(d1, d2)
        this%normal = cross / norm2(cross)

    end subroutine panel_calc_norm

    subroutine panel_calc_area(this)

        implicit none
        
        class(panel), intent(inout) :: this
        real,dimension(3) :: d1, d2

        ! Get side vectors
        d1 = this%get_vertex_loc(2) - this%get_vertex_loc(1)
        d2 = this%get_vertex_loc(3) - this%get_vertex_loc(2)

        ! Calculate area from  cross product
        this%area = 0.5*norm2(cross_product(d1, d2))

        ! Check for zero area
        if (this%area < 1.e-12) then
            write(*,*) "!!! Panel", this%index, "has zero area. Quitting..."
            write(*,*) "!!! Vertex 1 (", this%get_vertex_index(1), "):", this%get_vertex_loc(1)
            write(*,*) "!!! Vertex 2 (", this%get_vertex_index(2), "):", this%get_vertex_loc(2)
            write(*,*) "!!! Vertex 3 (", this%get_vertex_index(3), "):", this%get_vertex_loc(3)
            stop
        end if

    end subroutine panel_calc_area

    subroutine panel_calc_centroid(this)
    
        implicit none
        
        class(panel), intent(inout) :: this
        real, dimension(3) :: sum
        integer :: i

        ! Get average of corner points
        sum = 0.
        do i=1,this%N
            sum = sum + this%get_vertex_loc(i)
        end do

        ! Set centroid
        this%centr = sum/this%N

    end subroutine panel_calc_centroid

    subroutine panel_calc_velocity_vector(this, v_inf)
        ! Calculates the panel velocity direction vector based on panel geometry
        implicit none
        
        class(panel), intent(inout) :: this
        real, dimension(3), intent(inout) :: v_inf

        this%u = cross_product(cross_product(v_inf, this%normal), this%normal)

        this%u = this%u / sqrt(this%u(1)**2 + this%u(2)**2 + this%u(3)**2)

    end subroutine panel_calc_velocity_vector

    subroutine panel_calc_pressure_newton(this, m, c_pmax)
        implicit none
        
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

        this%c_p = - ((gamma + 1)/2 * this%theta**2 * (sqrt(1 + ( 4 / ((gamma + 1) * m * this%theta))**2) - 1))

    end subroutine panel_calc_pressure_prandtl

    subroutine panel_calc_projected_outline(this, projected_verts)

        implicit none
        
        class(panel), intent(inout) :: this
        type(vertex), dimension(:), allocatable, intent(in) :: projected_verts

        real, dimension(3) :: x_s, y_s, z_s

        x_s = [projected_verts(this%get_vertex_index(1))%location(1), &
                projected_verts(this%get_vertex_index(2))%location(1), &
                projected_verts(this%get_vertex_index(3))%location(1)]

        this%x_min_pro = min(x_s(1),x_s(2),x_s(3))
        this%x_max_pro = max(x_s(1),x_s(2),x_s(3))

        y_s = [projected_verts(this%get_vertex_index(1))%location(2), &
                projected_verts(this%get_vertex_index(2))%location(2), &
                projected_verts(this%get_vertex_index(3))%location(2)]

        this%y_min_pro = min(y_s(1),y_s(2),y_s(3))
        this%y_max_pro = max(y_s(1),y_s(2),y_s(3))

        z_s = [projected_verts(this%get_vertex_index(1))%location(3), &
                projected_verts(this%get_vertex_index(2))%location(3), &
                projected_verts(this%get_vertex_index(3))%location(3)]

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

    function panel_get_vertex_index(this, i) result(index)

        implicit none
        
        class(panel), intent(in) :: this
        integer, intent(in) :: i
        integer :: index

        index = this%vertices(i)%ptr%ind

    end function panel_get_vertex_index

    subroutine panel_calc_surface_equation(this, vertices)

        implicit none
        
        class(panel), intent(inout) :: this
        type(vertex), dimension(:), allocatable, intent(in) :: vertices

        real, dimension(this%N, this%N) :: A
        real, dimension(:), allocatable :: x
        real, dimension(this%N) :: b
        real, dimension(this%N) :: loc
        integer :: i, j

        ! Fill array A
        do i = 1, this%N
            do j = 1, this%N
                loc = this%get_vertex_loc(j)
                ! Avoid singular matrices
                if (loc(i) /= 0) then
                    A(i,j) = loc(i) 
                else
                    A(i,j) = loc(i) + 1.0e-12
                end if
            end do
            ! If vert(i) is (0,0,0) a_4 = 0.0
            if (all(this%get_vertex_loc(i) == 0.0)) then
                    this%a_4 = 0.0
            end if
        end do

        ! Fill array b
        do i = 1, this%N
            b(i) = -this%a_4
        end do

        ! Solve equation
        call lu_solve(this%N, A, b, x)

        ! Get equation coefficients
        this%a_1 = x(1) 
        this%a_2 = x(2) 
        this%a_3 = x(3)

    
    end subroutine panel_calc_surface_equation

end module panel_mod
