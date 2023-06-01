module flow_mod
    implicit none
    
    type flow
    
        real, dimension(3) :: freestream=[1,0,0]
        real :: m=5, gamma=1.4, m_q=0.0, q, delta_q

        contains
        
        procedure :: get_prandtl_meyer_matching_angle => flow_get_prandtl_meyer_matching_angle

    end type flow
    
contains

    subroutine flow_get_prandtl_meyer_matching_angle(this)
        implicit none

        class(flow), intent(inout) :: this
        real :: P_free_to_stag, M_guess=1.1
        real, dimension(20) :: M_q, P_q, Qs
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
        M_q(2) = M_guess + .4

        ! Interpolate and repeat to find matching point
        do i = 2, 20
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
            if (P_q(ind(1)) == P_q(ind(2))) then
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
        end if
        ! Calculate angle of matching point
        this%delta_q = asin(sqrt((this%q - P_free_to_stag)/(1 - P_free_to_stag)))
        write(*,*) "Matching inclination: ", this%delta_q


    end subroutine
    
end module flow_mod