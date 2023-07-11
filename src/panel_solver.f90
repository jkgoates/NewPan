module panel_solver_mod
    use panel_mod
    use base_geom_mod
    use flow_mod
    use surface_mesh_mod
    use math_mod
    use linked_list_mod

    implicit none
    
    type panel_solver

        type(flow) :: freestream
        character(len=:), allocatable :: windward_method, leeward_method
        real :: c_pmax


        contains

            ! Initialization
            procedure :: init => panel_solver_init
            procedure :: parse_solver_settings => panel_solver_parse_solver_settings

            ! Pressure calculations
            procedure :: newton_max_pressure => panel_solver_newton_max_pressure
            procedure :: calc_pressures_kaufman => panel_solver_calc_pressures_kaufman
            procedure :: calc_pressures_str_newton => panel_solver_calc_pressures_str_newton
            procedure :: calc_pressures_mod_newton => panel_solver_calc_pressures_mod_newton
            procedure :: calc_pressures_free_prandtl => panel_solver_calc_pressures_free_prandtl

            ! Panel calculations
            procedure :: calc_angles => panel_solver_calc_angles
            procedure :: calc_forces => panel_solver_calc_forces

            ! Solve
            procedure :: solve => panel_solver_solve

            ! Surface property calculations
            procedure :: surface_velocities => panel_solver_surface_velocities
            procedure :: point_velocities => panel_solver_point_velocities
            procedure :: velocity_equations => panel_solver_velocity_equations

    end type panel_solver

contains

    subroutine panel_solver_init(this, solver_settings, body_mesh, freestream_flow)

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(json_value), pointer,intent(in) :: solver_settings
        type(surface_mesh), intent(inout) :: body_mesh
        type(flow), intent(inout) :: freestream_flow
        
        ! Store
        this%freestream = freestream_flow

        ! Get solver settings
        call this%parse_solver_settings(solver_settings)

        ! Calculate max pressure coefficent
        call this%newton_max_pressure()

        ! Calculate panel angles
        call this%calc_angles(body_mesh)

        ! Initialize matching point calculations
        if (this%windward_method == 'kaufman') then
            call this%freestream%init_kaufman()
            call this%freestream%get_prandtl_max_turning_angle(this%freestream%m_sep)
        end if

    end subroutine panel_solver_init

    subroutine panel_solver_parse_solver_settings(this, solver_settings)

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(json_value), pointer, intent(in) :: solver_settings

        ! Get solver methods
        call json_xtnsn_get(solver_settings, 'windward_method', this%windward_method, 'modified-newtonian')
        call json_xtnsn_get(solver_settings, 'leeward_method', this%leeward_method, 'none')
        
    end subroutine panel_solver_parse_solver_settings

    subroutine panel_solver_solve(this, body_mesh)

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body_mesh

        write(*,*)
        write(*,'(a)', advance='no') " Solving for inviscid pressures... "

        ! Calculate impact pressures
        select case (this%windward_method)
        case ('kaufman')
            call this%calc_pressures_kaufman(body_mesh)
        case ('straight-newtonian')
            call this%calc_pressures_str_newton(body_mesh)
        case ('modified-newtonian')
            call this%calc_pressures_mod_newton(body_mesh)
        case default
            write(*,*) "!!! Invalid windward method selected. Quitting..."
            stop
        end select

        ! Calculate shadow pressures
        select case (this%leeward_method)
        case ('none')
            continue
        case ('kaufman')
            continue
        case ('prandtl-meyer')
            if (this%windward_method /= 'kaufman') then
                call this%calc_pressures_free_prandtl(body_mesh)
            end if
        case default
            write(*,*) "!!! Invalid leeward method selected. Quitting..."
            stop
        end select

        write(*,*) "Done."


        ! Calculate pressure dependent properties
        write(*,'(a)', advance='no') " Calculating surface properties... "

        call this%surface_velocities(body_mesh)
        call this%point_velocities(body_mesh)
        call this%velocity_equations(body_mesh)

        write(*,*) "Done."


    end subroutine panel_solver_solve

    subroutine panel_solver_calc_forces(this, body_mesh)
        ! Calculates the force coefficients on the model
        implicit none

        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body_mesh

        integer :: i
        
        ! Loop through panels and sum up discrete force coefficients
        do i = 1, body_mesh%N_panels
    
            body_mesh%panels(i)%dC_f = body_mesh%panels(i)%c_p * body_mesh%panels(i)%normal * body_mesh%panels(i)%area
            body_mesh%C_f = body_mesh%C_f + body_mesh%panels(i)%dC_f
    
        end do
    
        ! Make C_f non-dimensional
        body_mesh%C_f = body_mesh%C_f / body_mesh%S_ref

    end subroutine panel_solver_calc_forces

    subroutine panel_solver_newton_max_pressure(this)
        ! Calculate max pressure coefficent for use in newtonian based methods
        implicit none

        class(panel_solver), intent(inout) :: this

        real :: gamma, M_inf
        gamma = this%freestream%gamma
        M_inf = this%freestream%M_inf

        if (this%windward_method == 'modified-newtonian' .or. this%windward_method == 'kaufman') then
            this%c_pmax = (2/(gamma*M_inf**2))*( (( (((gamma + 1)**2) * M_inf**2 )/ &
                ((4*gamma*M_inf**2)-2*(gamma-1)))**(gamma/(gamma-1))) * &
                ((1-gamma+(2*gamma*M_inf**2))/(gamma+1)) -1)
        else if (this%windward_method == 'straight-newtonian') then
            this%c_pmax = 2
        else
            write(*,*) "!!! Invalid solver method. Quitting..."
            stop
        end if

    end subroutine panel_solver_newton_max_pressure

    subroutine panel_solver_calc_angles(this, body)
        implicit none

        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body
        real, dimension(3) :: norm
        integer :: i

        do i = 1, body%N_panels

            norm = body%panels(i)%normal

            ! Calculate inclination of the panel normal vector with respect to the freestream
            body%panels(i)%phi = acos(dot_product(this%freestream%v_inf,norm))

            ! Calculate panel inclination
            body%panels(i)%theta = body%panels(i)%phi -pi/2
        end do

    end subroutine panel_solver_calc_angles
    
    
    subroutine panel_solver_calc_pressures_str_newton(this, body_mesh)

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body_mesh
        
        integer :: i

        ! Loop through panels
        do i = 1, body_mesh%N_panels
            ! Check if it faces the freestream

            if (body_mesh%panels(i)%theta > 0) then
                call body_mesh%panels(i)%calc_pressure_newton(this%c_pmax)
            end if

        end do
        
    end subroutine panel_solver_calc_pressures_str_newton

    subroutine panel_solver_calc_pressures_mod_newton(this, body_mesh)

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body_mesh
        
        integer :: i

        ! Loop through panels
        do i = 1, body_mesh%N_panels
            ! Check if it faces the freestream
            if (body_mesh%panels(i)%theta > 0) then
                call body_mesh%panels(i)%calc_pressure_newton(this%c_pmax)
            end if

        end do

    end subroutine panel_solver_calc_pressures_mod_newton

    subroutine panel_solver_calc_pressures_kaufman(this, body_mesh)

        implicit none

        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body_mesh
        
        integer :: i
        real :: m_i, m, gamma, P_infty_ratio

        

        ! Declutter
        gamma = this%freestream%gamma
        m = this%freestream%M_inf

        do i = 1, body_mesh%N_panels
            ! Modified Newtonian for blunt panels
            if (body_mesh%panels(i)%theta > this%freestream%delta_q) then
                body_mesh%panels(i)%c_p = this%c_pmax * sin(body_mesh%panels(i)%theta)**2
            ! Prandtl meyer expansion from matching point to assumed flow seperation
            else if (body_mesh%panels(i)%theta > this%freestream%theta_sep) then
                m_i = this%freestream%solve_prandtl_meyer_mach(body_mesh%panels(i)%theta)
                body_mesh%panels(i)%c_p = 2 / (gamma * m**2) * ((1/this%freestream%P_free_to_stag * &
                            ((2 / (2 + (gamma - 1) * m_i**2))**(gamma / (gamma - 1)))) - 1)
            ! Continued expansion
            else if (body_mesh%panels(i)%theta > this%freestream%theta_min) then
                call body_mesh%panels(i)%calc_pressure_prandtl(this%freestream%gamma, this%freestream%M_inf, &
                                                               this%freestream%m_sep)
            else
                body_mesh%panels(i)%c_p = 0
                body_mesh%panels(i)%seperated = .true.
            end if
        end do

    end subroutine panel_solver_calc_pressures_kaufman

    subroutine panel_solver_calc_pressures_free_prandtl(this, body_mesh)

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body_mesh

        integer :: i

        call this%freestream%get_prandtl_max_turning_angle(this%freestream%M_inf)

        do i = 1, body_mesh%N_panels
            ! Check panel faces away from freestream
            if (body_mesh%panels(i)%theta < 0 .and. body_mesh%panels(i)%theta > this%freestream%theta_min) then
                call body_mesh%panels(i)%calc_pressure_prandtl(this%freestream%gamma, this%freestream%M_inf, &
                                                               this%freestream%M_inf)
            else if (body_mesh%panels(i)%theta <= this%freestream%theta_min) then
                body_mesh%panels(i)%c_p = 0
                body_mesh%panels(i)%seperated = .true.
            end if
        end do

    end subroutine panel_solver_calc_pressures_free_prandtl

    subroutine panel_solver_surface_velocities(this, body_mesh)
    ! Calculates the surface velocity ratio based on surface pressure coefficients
        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body_mesh
        
        integer :: i
        real, dimension(3) :: v_inf

        ! Initialize
        v_inf = reshape(this%freestream%v_inf, [3])
        ! Loop through panels
        do i = 1, body_mesh%N_panels
            
            ! Calculate the surface mach number
            call body_mesh%panels(i)%calc_centroid_mach(this%freestream%P_free_to_stag, &
                                                this%freestream%M_inf, this%freestream%gamma)

            ! Calculate the surface velocity ratio
            call body_mesh%panels(i)%calc_centroid_velocity(this%freestream%M_inf, this%freestream%gamma)

            ! Calculate velocity vectors
            call body_mesh%panels(i)%calc_velocity_vector(v_inf)

        end do

    end subroutine panel_solver_surface_velocities

    subroutine panel_solver_point_velocities(this, body_mesh)

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body_mesh

        type(list) :: adj_panels

        real, dimension(3) :: sum_num, sum_den

        integer :: i, j, pan

        ! Loop through vertices and calculate velocities
        do i = 1, body_mesh%N_verts

            adj_panels = body_mesh%vertices(i)%get_panels()

            sum_num = 0
            sum_den = 0

            ! Loop through panels and sum up vertices
            do j = 1, adj_panels%len()

                call adj_panels%get(j,pan)
                sum_num = sum_num + body_mesh%panels(pan)%velocity
                sum_den = sum_den + 1

            end do

            body_mesh%vertices(i)%V = sum_num/sum_den
        end do

    end subroutine panel_solver_point_velocities


    subroutine panel_solver_velocity_equations(this, body_mesh)

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(surface_mesh), intent(inout) :: body_mesh

        type(panel), dimension(3) :: neighbors
        integer :: i



        ! Loop through panels and calculate the velocity equations
        do i = 1, body_mesh%N_panels
            neighbors = (/body_mesh%panels(body_mesh%panels(i)%abutting_panels(1)), &
                          body_mesh%panels(body_mesh%panels(i)%abutting_panels(2)), &
                          body_mesh%panels(body_mesh%panels(i)%abutting_panels(3))/)
            call body_mesh%panels(i)%calc_velocity_equations(neighbors)
        end do

    end subroutine panel_solver_velocity_equations


end module panel_solver_mod