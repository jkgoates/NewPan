module flow_mod
    implicit none
    
    type flow
    
        real, dimension(3) :: freestream=[-1,0,0]
        real :: m=20, gamma=1.4, m_q=0.0, q, delta_q, pi=3.14159268, theta_min
        real :: P_free_to_stag, m_sep

        contains
        
        procedure :: get_prandtl_meyer_matching_angle => flow_get_prandtl_meyer_matching_angle
        procedure :: solve_prandtl_meyer_mach => flow_solve_prandtl_meyer_mach
        procedure :: get_seperation_mach => flow_get_seperation_mach
        procedure :: init_kaufman => flow_init_kaufman

    end type flow
    
contains

    subroutine flow_init_kaufman(this)

        implicit none
        
        class(flow), intent(inout) :: this

        call this%get_prandtl_meyer_matching_angle()
        call this%get_seperation_mach()



    end subroutine flow_init_kaufman

    subroutine flow_get_prandtl_meyer_matching_angle(this)
        implicit none

        class(flow), intent(inout) :: this
        real :: P_free_to_stag, M_guess=1.1
        real, dimension(50) :: M_q, P_q, Qs
        integer, dimension(2) :: ind
        integer :: i, j

        ! Calculate freestream to stagnation pressure ratio
        P_free_to_stag = (2/((this%gamma + 1)*this%m**2))**(this%gamma/(this%gamma - 1)) &
                        * ((2 * this%gamma * this%m**2 - (this%gamma - 1))/(this%gamma + 1))**(1/(this%gamma - 1))
        ! Make an initial guess at mach number
        M_q(1) = M_guess

        ! Calculate matching point pressure ratio
        Qs(1) = (2/(2 + (this%gamma - 1) * M_q(1)**2))**(this%gamma/(this%gamma - 1))

        P_q(1) = Qs(1) * (1 - (this%gamma**2 * M_q(1)**4 * Qs(1))/(4 * (M_q(1)**2 - 1) * (1 - Qs(1))))

        
        ! Make a second guess at mach number
        M_q(2) = M_guess + .2

        ! Interpolate and repeat to find matching point
        do i = 2, 50
            Qs(i) = (2/(2 + (this%gamma - 1) * M_q(i)**2))**(this%gamma/(this%gamma - 1))

            P_q(i) = Qs(i) * (1 - (this%gamma**2 * M_q(i)**4 * Qs(i))/(4 * (M_q(i)**2 - 1) * (1 - Qs(i))))
            ! Loop through guesses and find the closest guesses
            ind = [1,2]
            if(i>=3) then
                do j = 1, i
                    if (j == ind(1) .or. j == ind(2)) then
                        cycle
                    end if
                    if ( abs(P_free_to_stag - P_q(j)) < abs(P_free_to_stag - P_q(ind(2))) ) then
                        if (abs(P_free_to_stag - P_q(ind(2))) < abs(P_free_to_stag - P_q(ind(1)))) then
                            ind(1) = ind(2)
                        end if
                        ind(2) = j
                        cycle
                    end if
                    if ( abs(P_free_to_stag - P_q(j)) < abs(P_free_to_stag - P_q(ind(1))) ) then
                        ind(1) = j
                        cycle
                    end if
                end do
            end if
            
            ! Check for convergence
            if (abs(P_q(i)-P_free_to_stag) < 1e-6) then
                this%M_q = M_q(i)
                this%Q = Qs(i)
                exit
            end if

            ! Interpolate if not at end
            if (i /= 20) then
                M_q(i+1) = (M_q(ind(2)) - M_q(ind(1))) / (P_q(ind(2)) - P_q(ind(1))) &
                            * (P_free_to_stag - P_q(ind(1))) + M_q(ind(1))

            end if


        end do
        
        if (this%m_q < 1 .or. this%m_q > 2) then 
            write(*,*) "!!! Solution did not converge or reported an invalid matching mach number, quitting..."
            stop
        else
            write(*,'(a i2 a)') "Solution converged in ", i, " iterations."
            write(*,*) "Matching Mach Number: ", this%m_q
            this%P_free_to_stag = P_free_to_stag
        end if
        ! Calculate angle of matching point
        this%delta_q = asin(sqrt((this%q - P_free_to_stag)/(1 - P_free_to_stag)))
        write(*,*) "Matching inclination: ", this%delta_q


    end subroutine flow_get_prandtl_meyer_matching_angle

    function flow_solve_prandtl_meyer_mach(this, theta) result(m_2)
        implicit none

        class(flow), intent(inout) :: this
        real, intent(in) :: theta

        real :: m_2, nu_1, nu_2, nu_max, max_defl
        real, dimension(50) :: nu, Ms
        integer, dimension(2) :: ind
        integer :: i, j

        ! Evaluate the prandtl meyer equation at mach one
        nu_1 = sqrt((this%gamma + 1) / (this%gamma - 1)) * atan2(sqrt((this%gamma - 1) / &
                (this%gamma + 1) * (this%m_q**2 -1)),1.) - atan2(sqrt(this%m_q**2 - 1),1.)

        ! Calculate max turning angle
        nu_max = sqrt((this%gamma + 1) / (this%gamma - 1)) * atan2(sqrt((this%gamma - 1) / &
                    (this%gamma + 1) * (this%m_sep**2 -1)),1.) - atan2(sqrt(this%m_sep**2 - 1),1.)
        max_defl = nu_max - nu_1
        this%theta_min = this%delta_q - max_defl

        ! Initialize interpolation arrays
        nu(1) = nu_1
        nu(2) = nu_max
        Ms(1) = this%m_q
        Ms(2) = this%m_sep

        ! Make sure theta does not excede theta_max
        if (theta >= this%theta_min) then

            ! Calculate nu_2
            nu_2 = nu_1 + this%delta_q - theta

            

            ind = [1,2]


            do i = 2, 50
                ! Recalculate prandtl-meyer function
                nu(i) = sqrt((this%gamma + 1) / (this%gamma - 1)) * atan2(sqrt((this%gamma - 1) / &
                        (this%gamma + 1) * (Ms(i)**2 -1)),1.) - atan2(sqrt(Ms(i)**2 - 1),1.)
                

                ! Loop through guesses and find the closest guesses

                do j = 1, i
                    if (j == ind(1) .or. j == ind(2)) then
                        cycle
                    end if
                    if ( abs(nu_2 - nu(j)) < abs(nu_2 - nu(ind(2))) ) then
                        if (abs(nu_2 - nu(ind(2))) < abs(nu_2 - nu(ind(1)))) then
                            ind(1) = ind(2)
                        end if
                        ind(2) = j
                        cycle
                    end if
                    if ( abs(nu_2 - nu(j)) < abs(nu_2 - nu(ind(1))) ) then
                        ind(1) = j
                        cycle
                    end if
                end do

                ! Check for convergence
                if (abs(nu(i) - nu_2) < 1e-5) then
                    m_2 = Ms(i)
                    exit
                end if

                ! Interpolate if not at end
                if (i /= 50) then
                    Ms(i+1) = (Ms(ind(2)) - Ms(ind(1))) / (nu(ind(2)) - nu(ind(1))) &
                                * (nu_2 - nu(ind(1))) + Ms(ind(1))

                end if


            end do

            ! Make sure calculated mach has not exceeded seperation mach
            if (Ms(i) < this%m_sep) then
                m_2 = Ms(i)
            end if
        
        else 
            write(*,*) "!!! Theta beyond seperation angle. Quitting..."
            stop
        end if



    end function flow_solve_prandtl_meyer_mach

    subroutine flow_get_seperation_mach(this)
        implicit none

        class(flow), intent(inout) :: this

        this%m_sep = sqrt((2 / (this%P_free_to_stag**((this%gamma-1)/this%gamma)) - 2)/(this%gamma -1))

        write(*,*) "Separation Mach Number: ", this%m_sep

    end subroutine

    
end module flow_mod