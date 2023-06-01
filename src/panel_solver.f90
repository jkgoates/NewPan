module panel_solver_mod
    use panel_mod
    use base_geom_mod

    contains


    subroutine panel_solver_calc_forces(C_f, N_panels, panels, ref_area)
        ! Calculates the force coefficients on the model
        implicit none

        real, dimension(3), intent(inout) :: C_f
        integer, intent(in) :: N_panels
        type(panel), dimension(:), allocatable, intent(inout) :: panels
        real, intent(in) :: ref_area

        integer :: i
        
        ! Loop through panels and sum up discrete force coefficients
        do i = 1, N_panels
    
            panels(i)%dC_f = panels(i)%c_p * panels(i)%normal * panels(i)%area
            C_f = C_f + panels(i)%dC_f
    
        end do
    
        ! Make C_f non-dimensional
        C_f = C_f / ref_area

    end subroutine panel_solver_calc_forces

    subroutine panel_solver_max_pressure(windward_method, gamma, m, c_pmax)
        ! Calculate max pressure coefficent for use in newtonian based methods
        implicit none

        character(len=:), allocatable, intent(in) :: windward_method
        real, intent(in) :: gamma, m
        real, intent(out) :: c_pmax


        if (windward_method == 'modified-newtonian') then
            write(*,*) "Solving for pressure coefficients using Modified Newtonian Method..."
            c_pmax = (2/(gamma*m**2))*( (( (((gamma + 1)**2) * m**2 )/((4*gamma*m**2)-2*(gamma-1)))**(gamma/(gamma-1))) * &
                ((1-gamma+(2*gamma*m**2))/(gamma+1)) -1)
        else if (windward_method == 'straight-newtonian') then
            c_pmax = 2
        else
            write(*,*) "!!! Invalid solver method. Quitting..."
            stop
        end if

    end subroutine panel_solver_max_pressure

    subroutine panel_solver_calc_angles(panels, freestream, pi, N_panels)
        implicit none

        real, dimension(:), allocatable, intent(in) :: freestream
        real, intent(in) :: pi
        type(panel), dimension(:), allocatable, intent(inout) :: panels
        integer, intent(in) :: N_panels

        real, dimension(3) :: norm
        integer :: i

        do i = 1, N_panels

            norm = panels(i)%normal

            ! Calculate inclination of the panel normal with respect to the freestream
            panels(i)%phi = acos(dot_product(freestream,norm))

            ! Calculate panel inclination
            panels(i)%theta = panels(i)%phi -pi/2
        end do

    end subroutine panel_solver_calc_angles
    
    
    subroutine panel_solver_calc_pressures(freestream, gamma, leeward_method, N_panels, panels, m, pi, c_pmax)
        ! Calculates pressure coefficients on each panel
        implicit none
        
        real, dimension(:), allocatable, intent(in) :: freestream
        real, intent(in) :: gamma, m, pi, c_pmax
        character(len=:), allocatable, intent(in) :: leeward_method
        integer, intent(in) :: N_panels
        type(panel), dimension(:), allocatable, intent(inout) :: panels

        real, dimension(3) :: norm
        integer :: i

        ! Loop through and calculate pressure for all panels
        do i = 1, N_panels
            

            
            ! Calculate pressure
            call panel_calc_pressure(panels(i), leeward_method, gamma, pi, m, c_pmax)

        end do

    end subroutine panel_solver_calc_pressures

    function panel_solver_prandtl_meyer_max_angle(m, gamma) result(omega)
        ! Identify the maximum turning angle for a prandtl meyer expansion
        implicit none

        real, intent(in) :: m, gamma

        real :: omega

        omega = sqrt((gamma + 1)/(gamma -1)) * atan(sqrt((gamma - 1)/(gamma + 1) * (m**2 -1))) - atan(sqrt(m**2 -1))
        write(*,*) omega

    end function panel_solver_prandtl_meyer_max_angle

    subroutine panel_solver_calc_seperation(leeward_method, N_panels, panels, m, gamma)
        implicit none

        character(len=:), allocatable, intent(in) :: leeward_method
        real, intent(in) :: gamma, m
        integer, intent(in) :: N_panels
        type(panel), dimension(:), allocatable, intent(inout) :: panels

        real :: omega
        integer :: i

        if (leeward_method == "prandtl-meyer") then
            omega = panel_solver_prandtl_meyer_max_angle(m, gamma)

            do i = 1, N_panels
                call panel_prandtl_meyer_seperation(panels(i), omega)

            end do
        end if

    end subroutine panel_solver_calc_seperation

end module panel_solver_mod