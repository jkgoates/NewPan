program main
    implicit none
    
    ! Initialize variables
    real :: gamma=1.4, m=100, c_pmax, c_p, pi=3.1415926, theta
    real, dimension(3) :: freestream, norm, p_1, p_2, p_3

    freestream = [100,0,0]

    p_1 = [0,0,0]
    p_2 = [1,1,1]
    p_3 = [0,1,0]

    norm = [25,0,0]
    
    ! Normalize freestream and normal vectors
    freestream = freestream/(sqrt(freestream(1)**2 + freestream(2)**2 + freestream(3)**2))
    
    norm = norm/(sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2))

    ! Calculate max pressure coefficent    
    c_pmax = (2/(gamma*m**2))*( (( (((gamma + 1)**2) * m**2 )/((4*gamma*m**2)-2*(gamma-1)))**(gamma/(gamma-1))) * &
        ((1-gamma+(2*gamma*m**2))/(gamma+1)) -1)

    theta = acos(dot_product(freestream,norm)) - pi/2

    write(*,*) theta

    ! Check if panel faces freestream
    if ( (theta <= pi/2) .and. (theta >= -pi/2) ) then
        c_p = 0
        write(*,*) 'Panel does not face the freestream'
    else
        c_p = c_pmax*sin(theta)**2
    end if

    ! Print pressure coefficent
    write(*,*) c_pmax
    write(*,*) c_p

end program main