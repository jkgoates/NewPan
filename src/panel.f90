! Subroutines for working with an individual panel
module panel_mod

    use base_geom_mod
    use math_mod
    use linalg_mod
    use flow_mod
    implicit none

    
    type panel

        real, dimension(3) :: normal, dC_f, centr, centr_pro,velocity 
        real :: c_p=0, area, M, V, theta, phi, shadowed=0
        real :: x_min_pro, x_max_pro, y_min_pro, y_max_pro, z_min_pro, z_max_pro ! Outline of projected panel
        real, dimension(3) :: s12_pro, s23_pro, s31_pro ! Projected side vectors
        integer :: index ! Panel index in mesh array
        integer :: N ! Number of sides
        type(vertex_pointer), dimension(:), allocatable :: vertices
        integer, dimension(3) :: abutting_panels ! Indices of panels abutting this one
        integer, dimension(3) :: edges ! Indices of the edges of this panel
        logical :: seperated=.false.
        integer :: order

        !! Panel equations
        real :: a1, a2, a3, a4=1.0 ! Surface Equation Coefficients F(x,y,z) = a1*x + a2*y + a3*z + a4 = 0
        ! Velocity Equations
        real :: b1, b2, b3, b4, b5, b6 ! u component coefficients u = b1*indep(1) + b2*indep(2) + b3
        real :: c1, c2, c3, c4, c5, c6 ! v component coefficients v = c1*indep(1) + c2*indep(2) + c3
        real :: d1, d2, d3, d4, d5, d6 ! w component coefficients w = d1*indep(1) + d2*indep(2) + c3

        ! Variable dependencies
        integer, dimension(2) :: indep=(/2,3/) ! Independent Variables for the panel. x=1,y=2,z=3
        integer :: dep=1 ! Dependent Variable for the panel.

        contains

        ! Panel Initialization
        procedure :: init => panel_init
        procedure :: calc_derived_geom => panel_calc_derived_geom
        procedure :: get_characteristic_length => panel_get_characteristic_length

        ! Geometry calculations
        procedure :: calc_norm => panel_calc_norm
        procedure :: calc_centroid => panel_calc_centroid
        procedure :: calc_area => panel_calc_area
        procedure :: calc_projected_outline => panel_calc_projected_outline
        procedure :: define_dependent_direction => panel_define_dependent_direction

        ! Pressure calculations
        procedure :: calc_pressure_newton => panel_calc_pressure_newton
        procedure :: calc_pressure_prandtl => panel_calc_pressure_prandtl

        ! Getters
        procedure :: get_vertex_loc => panel_get_vertex_loc
        procedure :: get_vertex_index => panel_get_vertex_index
        procedure :: get_vertex_velocity => panel_get_vertex_velocity
        procedure :: get_panel_velocity => panel_get_panel_velocity

        ! Panel Centroid Properties
        procedure :: calc_centroid_mach => panel_calc_centroid_mach
        procedure :: calc_centroid_velocity => panel_calc_centroid_velocity
        procedure :: calc_velocity_vector => panel_calc_velocity_vector

        ! Property distributions
        procedure :: set_distribution => panel_set_distribution
        procedure :: calc_surface_equation => panel_calc_surface_equation
        procedure :: calc_velocity_equations => panel_calc_velocity_equations

        ! Streamline functions
        procedure :: rk4 => panel_rk4
        procedure :: calc_streamline => panel_calc_streamline
        procedure :: point_on_panel => panel_point_on_panel
        procedure :: restrict_point_to_panel => panel_restrict_point_to_panel
        procedure :: find_edge_intersection => panel_find_edge_intersection


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


    subroutine panel_define_dependent_direction(this,freestream)
        ! Defines which directions are independent and dependent on the panel for streamline calculations
        implicit none
        
        class(panel), intent(inout) :: this
        type(flow), intent(in) :: freestream

        integer :: dep

        if (abs(this%normal(maxloc(abs(this%normal),dim=1))) > abs(2*this%normal(maxloc(abs(freestream%v_inf),dim=1)))) then
            dep = maxloc(abs(this%normal),dim=1)
        else
            dep = maxloc(abs(freestream%v_inf), dim=1)
        end if


        select case(dep)
        case(1)
            continue
        case(2)
            this%dep = dep
            this%indep = (/1,3/)
        case(3)
            this%dep = dep
            this%indep = (/1,2/)
        end select

    end subroutine panel_define_dependent_direction

    
    function panel_get_characteristic_length(this) result (l)
        ! Returns the square root of the panel area

        implicit none
        
        class(panel), intent(in) :: this

        real :: l

        l = sqrt(this%area)

    end function panel_get_characteristic_length


    subroutine panel_calc_pressure_newton(this, c_pmax)
        implicit none
        
        ! Calculates the pressure on a panel
        real, intent(in) :: c_pmax
        class(panel), intent(inout) :: this

        ! Newtonian Method on impact region
        this%c_p = c_pmax*sin(this%theta)**2
    
    end subroutine panel_calc_pressure_newton


    subroutine panel_calc_pressure_prandtl(this, gamma, M_inf, M_1)
        
        class(panel), intent(inout) :: this
        real, intent(in) :: gamma, M_inf, M_1

        real :: p2_p1

        !this%c_p = - ((gamma + 1)/2 * this%theta**2 * (sqrt(1 + ( 4 / ((gamma + 1) * M_inf * this%theta))**2) - 1))

        p2_p1 = (1 - (gamma - 1)/2 * M_1 * -this%theta)**(2*gamma/(gamma-1))
        this%c_p = 2/(gamma * M_inf**2) * (p2_p1 - 1)

    end subroutine panel_calc_pressure_prandtl


    subroutine panel_calc_centroid_mach(this, P, M_inf, gamma)

        class(panel), intent(inout) :: this
        real, intent(in) :: P, M_inf, gamma

        real :: P_stag_to_panel

        ! Avoid problems
        if (this%c_p < 0.0) then
            P_stag_to_panel = P*((0.0*gamma*M_inf**2)/2 + 1)
        else 
            P_stag_to_panel = P*((this%c_p*gamma*M_inf**2)/2 + 1)
        end if

        this%M = sqrt((P_stag_to_panel**(-(gamma - 1)/gamma) - 1) * (2 / (gamma - 1)))

    end subroutine panel_calc_centroid_mach


    subroutine panel_calc_centroid_velocity(this, M_inf, gamma)

        implicit none
        
        class(panel), intent(inout) :: this
        real, intent(in) :: M_inf, gamma

        this%V = (this%M/M_inf)*sqrt((1 + (gamma - 1)/2 * M_inf**2)/(1 + (gamma - 1)/2 * this%M**2))


    end subroutine panel_calc_centroid_velocity


    subroutine panel_calc_velocity_vector(this, v_inf)
        ! Calculates the panel velocity direction vector based on panel geometry
        implicit none
        
        class(panel), intent(inout) :: this
        real, dimension(3), intent(inout) :: v_inf

        this%velocity = cross_product(cross_product(this%normal, v_inf), this%normal)

        this%velocity = this%V * this%velocity / sqrt(this%velocity(1)**2 + this%velocity(2)**2 + this%velocity(3)**2)

    end subroutine panel_calc_velocity_vector


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


    function panel_get_vertex_velocity(this,i) result(velocity)

        implicit none
        
        class(panel), intent(in) :: this
        integer, intent(in) :: i
        real, dimension(3) :: velocity

        velocity = this%vertices(i)%ptr%V

    end function panel_get_vertex_velocity


    function panel_get_panel_velocity(this, point) result(velocity)

        implicit none
        
        class(panel), intent(in) :: this
        real, dimension(3), intent(in) :: point

        real, dimension(3) :: velocity

        ! Linear distribution
        velocity(1) = this%b1 + this%b2*point(this%indep(1)) + this%b3*point(this%indep(2))
        velocity(2) = this%c1 + this%c2*point(this%indep(1)) + this%c3*point(this%indep(2))
        velocity(3) = this%d1 + this%d2*point(this%indep(1)) + this%d3*point(this%indep(2))

        ! Second order
        if (this%order == 2) then
            velocity(1) = velocity(1) + this%b4*point(this%indep(1))**2 + this%b5*point(this%indep(2))**2 &
                          + this%b6*point(this%indep(1))*point(this%indep(2))
            velocity(2) = velocity(2) + this%c4*point(this%indep(1))**2 + this%c5*point(this%indep(2))**2 &
                          + this%c6*point(this%indep(1))*point(this%indep(2))
            velocity(3) = velocity(3) + this%d4*point(this%indep(1))**2 + this%d5*point(this%indep(2))**2 &
                          + this%d6*point(this%indep(1))*point(this%indep(2))
        end if

    end function panel_get_panel_velocity


    subroutine panel_set_distribution(this, order)

        implicit none
        
        class(panel), intent(inout) :: this
        integer, intent(in) :: order
        
        if (any(this%abutting_panels == 0)) then
            this%order = 1
        else
            this%order = order
        end if


    end subroutine


    subroutine panel_calc_surface_equation(this)

        ! Defines the panel by a surface equation F(x,y,z) = Ax + By + Cz + D = 0
        implicit none
        
        class(panel), intent(inout) :: this

        real, dimension(this%N, this%N) :: A
        real, dimension(:), allocatable :: x
        real, dimension(this%N) :: b
        real, dimension(this%N) :: loc
        integer :: i, j

        this%a1 = this%normal(1)
        this%a2 = this%normal(2)
        this%a3 = this%normal(3)
        this%a4 = - (this%a1*this%centr(1) + this%a2*this%centr(2) + this%a3*this%centr(3))

    end subroutine panel_calc_surface_equation

    
    subroutine panel_calc_velocity_equations(this, neighbors)
        ! Calculate the linear velocity distribution over the panel
        implicit none
        
        class(panel), intent(inout) :: this
        type(panel), dimension(3), intent(inout) :: neighbors

        real, dimension(this%N*this%order, this%N*this%order) :: A, A_STORED
        real, dimension(:), allocatable :: xu, xv, xw
        real, dimension(this%N*this%order) :: bu, bv, bw
        real, dimension(this%N) :: loc, v1, v2, v3, v4, v5, v6
        integer :: i, j, k, neighbor_vert
        integer, dimension(3) :: extra_verts, this_verts
        type(panel) :: neighbor

        if (this%order == 2) then
            this_verts = (/this%get_vertex_index(1), this%get_vertex_index(2), this%get_vertex_index(3)/)

            ! Find the unique vertex for the neighboring vertices
            neighbors_loop: do i = 1, 3
                neighbor = neighbors(i)
            
                neighbor_verts: do j = 1, 3
                    neighbor_vert = neighbor%get_vertex_index(j)

                    if (any(this_verts == neighbor_vert)) then
                        cycle neighbor_verts
                    else
                        extra_verts(i) = j
                        cycle neighbors_loop
                    end if

                end do neighbor_verts

            end do neighbors_loop
        end if

        ! Retrieve point velocities
        v1 = this%get_vertex_velocity(1)
        v2 = this%get_vertex_velocity(2)
        v3 = this%get_vertex_velocity(3)
        if (this%order == 2) then
            v4 = neighbors(1)%get_vertex_velocity(extra_verts(1))
            v5 = neighbors(2)%get_vertex_velocity(extra_verts(2))
            v6 = neighbors(3)%get_vertex_velocity(extra_verts(3))
        end if

        ! Fill A matrix

        do j = 1, this%N*this%order
            if (j <= this%N) loc = this%get_vertex_loc(j)
            if (this%order == 2) then
                if (j > this%N) loc = neighbors(j-this%N)%get_vertex_loc(extra_verts(j-this%N))
            end if

            if (any(loc == 0.0)) loc = loc + 1e-12

            A(j,1) = 1
            A(j,2) = loc(this%indep(1))
            A(j,3) = loc(this%indep(2))
            if (this%order ==2) then 
                A(j,4) = loc(this%indep(1))**2
                A(j,5) = loc(this%indep(2))**2
                A(j,6) = loc(this%indep(1))*loc(this%indep(2))
            end if
        end do


        ! Fill b vectors
        if (this%order == 2) then
            bu = (/v1(1),v2(1),v3(1),v4(1),v5(1),v6(1)/)
            bv = (/v1(2),v2(2),v3(2),v4(2),v5(2),v6(2)/)
            bw = (/v1(3),v2(3),v3(3),v4(3),v5(3),v6(3)/)
        else
            bu = (/v1(1),v2(1),v3(1)/)
            bv = (/v1(2),v2(2),v3(2)/)
            bw = (/v1(3),v2(3),v3(3)/)
        end if


        ! Solve
        A_STORED = A
        call lu_solve(this%N*this%order, A, bu, xu)
        A = A_STORED
        call lu_solve(this%N*this%order, A, bv, xv)
        A = A_STORED
        call lu_solve(this%N*this%order, A, bw, xw)

        ! Get coefficients
        this%b1 = xu(1)
        this%b2 = xu(2)
        this%b3 = xu(3)
        if (this%order == 2) then
            this%b4 = xu(4)
            this%b5 = xu(5)
            this%b6 = xu(6)
        end if

        this%c1 = xv(1)
        this%c2 = xv(2)
        this%c3 = xv(3)
        if (this%order == 2) then
            this%c4 = xv(4)
            this%c5 = xv(5)
            this%c6 = xv(6)
        end if

        this%d1 = xw(1)
        this%d2 = xw(2)
        this%d3 = xw(3)
        if (this%order == 2) then
            this%d4 = xw(4)
            this%d5 = xw(5)
            this%d6 = xw(6)
        end if

    end subroutine panel_calc_velocity_equations


    function panel_rk4(this, point, delta_s) result(new_point)
        ! Calculates the next point along a streamline on a panel surface

        implicit none
        
        class(panel), intent(in) :: this
        real, dimension(3), intent(in) :: point
        real, intent(in) :: delta_s

        real, dimension(3) :: k1,k2,k3,k4, new_point

        k1 = this%get_panel_velocity(point)/magnitude(this%get_panel_velocity(point))
        k2 = this%get_panel_velocity(point + 0.5*delta_s*k1)/magnitude(this%get_panel_velocity(point + 0.5*delta_s*k1))
        k3 = this%get_panel_velocity(point + 0.5*delta_s*k3)/magnitude(this%get_panel_velocity(point + 0.5*delta_s*k2))
        k4 = this%get_panel_velocity(point + delta_s*k3)/magnitude(this%get_panel_velocity(point + delta_s*k3))
        new_point = point + delta_s*onesixth*(k1 + 2*k2 + 2*k3 + k4)

    end function


    subroutine panel_calc_streamline(this, start_point, delta_s, end_point, panel_edge)
        ! Calculates a streamline across a panel starting at a point and cut off at the panel edge

        implicit none
        
        class(panel), intent(in) :: this
        real, dimension(3), intent(inout) :: start_point
        real, intent(in) :: delta_s
        integer, intent(inout) :: panel_edge
        real, dimension(3), intent(inout) :: end_point

        real, dimension(3) :: next_point, previous_point, intersection_point
        logical :: complete

        integer :: i, stat, point_stat

        !write(*,*) "Panel: ", this%index

        call this%restrict_point_to_panel(start_point)
        ! Ensure start point is in the panel boundary
        point_stat = this%point_on_panel(start_point)
        if (point_stat == 2) then
            
            write(*,*) "!!! Streamline Calculation Failed", point_stat
            write(*,*) "!!! Point at location ", start_point, " is not on panel ", this%index, ". Quitting..."
            stop
        end if


        ! Restrict Starting point to the panel surface
        call this%restrict_point_to_panel(start_point)

        complete = .false.
        previous_point = start_point
        ! Integrate point until it leaves the panel
        do while (.not. complete)
            ! Integrate
            next_point = this%rk4(previous_point, delta_s)

            ! Restrict to panel plane
            point_stat = this%point_on_panel(next_point)
            if(point_stat == 2) then
                call this%restrict_point_to_panel(next_point)
            end if

            ! Check if point has left the panel
            point_stat = this%point_on_panel(next_point)
            if (point_stat == 2) then
                write(*,*) "!!! Point to panel restriction failed. Point ", next_point, &
                            " on panel ", this%index, ". Quitting..."
                stop
            else if (point_stat == 0) then
                !write(*,*) "Point: ", next_point
                previous_point = next_point
                cycle
            else
                complete = .true. 
            end if
            
        end do
        
        ! Find edge intersection
        call this%find_edge_intersection(previous_point, next_point, end_point, panel_edge, stat)

        if (stat == 1) then
            write(*,*) "!!! Streamline Calculation Failed"
            write(*,*) "!!! Couldn't Calculate Streamline Edge intersection on panel ", this%index, ". Quitting..."
            stop
        end if
        

    end subroutine panel_calc_streamline


    function panel_point_on_panel(this, point) result(stat)
        ! Checks whether a point is on this panel

        ! STAT
        ! 0: on panel
        ! 1: outside panel
        ! 2: not on panel plane
        implicit none
        
        class(panel),intent(in) :: this
        real, dimension(3), intent(in) :: point

        integer :: stat, i
        real, dimension(this%N) :: delta

        delta = 0


        ! Check that the point fits the surface equation
        if (abs(this%a1*point(1)+this%a2*point(2)+this%a3*point(3)+this%a4) > 1.0e-12) then
            stat = 2
        else
            ! Loop through edges and check if point is inside
            do i = 1, this%N
                if (i==this%N) then
                    delta(i) = dot_product(cross_product(this%get_vertex_loc(i)-point, &
                               this%get_vertex_loc(i-2)-point),this%normal)
                else
                    delta(i) = dot_product(cross_product(this%get_vertex_loc(i)-point, &
                               this%get_vertex_loc(i+1)-point),this%normal)
                end if
            end do
            !write(*,*) "Velocity: ", this%get_panel_velocity(point)
            !write(*,*) "Delta: ", delta
            !write(*,*)

            if (delta(1) >= 0.0 .and. delta(2) >= 0.0 .and. delta(3) >= 0.0) then
                stat = 0
            else
                stat = 1
            end if
        end if

    end function panel_point_on_panel


    subroutine panel_restrict_point_to_panel(this, point)
        ! Restricts a point to the plane of the panel

        implicit none
        
        class(panel), intent(in) :: this
        real, dimension(3), intent(inout) :: point

        select case(this%dep)
        case(1)
            point(1) = - (this%a2*point(2) + this%a3*point(3) + this%a4) / this%a1
        case(2)
            point(2) = - (this%a1*point(1) + this%a3*point(3) + this%a4) / this%a2
        case(3)
            point(3) = - (this%a1*point(1) + this%a2*point(2) + this%a4) / this%a3
        end select

    end subroutine


    subroutine panel_find_edge_intersection(this, p1, p2, intersection, panel_edge, stat)
        ! Finds the intersection between two points and one of the panel's edges
        ! stat:
        !   0 = found intersection
        !   1 = didn't find interesection

        implicit none
        
        class(panel), intent(in) :: this
        real, dimension(3), intent(in) :: p1, p2
        real, dimension(3), intent(out) :: intersection
        integer, intent(inout) :: stat, panel_edge

        integer :: i, j
        real, dimension(2,2) :: A
        real, dimension(2) :: b
        real, dimension(:), allocatable :: x
        real, dimension(3) :: vert1, vert2, check

        ! Loop through the edges
        do i = 1,this%N

            vert1 = this%get_vertex_loc(i)
            if ( i == 3) then
                vert2 = this%get_vertex_loc(i-2)
            else
                vert2 = this%get_vertex_loc(i+1)
            end if


            ! Fill A Matrix
            do j = 1, 2
                A(j, 1) = vert2(this%indep(j)) - vert1(this%indep(j))
                if (A(j,1) == 0.0) A(j,1) = A(j,1) + 1e-12
                A(j, 2) = p1(this%indep(j)) - p2(this%indep(j))
                if (A(j,2) == 0.0) A(j,2) = A(j,2) + 1e-12
            end do

            ! b vector
            b = (/p1(this%indep(1)) - vert1(this%indep(1)), &
                p1(this%indep(2)) - vert1(this%indep(2))/)

            ! Solve
            call lu_solve(2, A, b, x)
            
            ! Check if it worked...
            check = x(1) * (vert2 - vert1) + x(2) * (p1-p2) - (p1 - vert1)

            !!write(*,*) "Check: ", check
            !write(*,*) x
            !write(*,*)
            
            if (all(abs(check) < 1e-12)) then
                ! Make sure the intersection is on the edge and vector
                if (all(x <= 1) .and. x(1) >= 0.0 .and. x(2) >= -1e-12) then
                    panel_edge = this%edges(i)
                    intersection = p1 + (x(2)+1e-4) * (p2-p1)
                    stat = 0
                    exit
                else 
                    stat = 1
                end if
            else
                stat = 1
            end if
        end do
        !write(*,*)
        !write(*,*)
        !write(*,*)




    end subroutine


end module panel_mod
