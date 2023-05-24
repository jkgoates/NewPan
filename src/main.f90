program main
   
    implicit none
    
    ! Initialize variables
    real :: gamma=1.4, m=10, c_pmax, c_p, pi=3.1415926, theta
    real, dimension(3) :: freestream, norm

    ! Panel Type
    type :: panel
        real, dimension(3) :: normal
        integer, dimension(3) :: vert_ind
        real :: cp
        integer :: N=3
    
    
        
    
    contains

    end type panel

    ! Vertex type
    type :: vertex
        integer :: ind
        real, dimension(3) :: location
        
    
    
    contains
        
    
    end type vertex

    ! Initialize variables
    integer :: i, j, N, i1, i2, i3, unit, N_panels, N_verts
    type(panel), dimension(:), allocatable :: panels
    type(vertex), dimension(:), allocatable :: vertices
    character(len=:), allocatable :: mesh_file
    character(len=200) :: dummy_read, line

    real, dimension(:,:), allocatable :: vertex_locs
    
    mesh_file = "meshes/ehlers_spindle_medium.vtk"
    ! Open vtk
    open(newunit=unit, file=mesh_file)

        ! Determine number of vertices
        read(unit,*) ! Header
        read(unit,*) ! Header
        read(unit,*) ! Header
        read(unit,*) ! Header
        read(unit,*) dummy_read, N_verts, dummy_read

        ! Allocate vertex arrays
        allocate(vertex_locs(3,N_verts))
        allocate(vertices(N_verts))

        ! Read One from each line
        do i=1,N_verts
            read(unit,*) vertex_locs(1,i), vertex_locs(2,i), vertex_locs(3,i)
        end do

        ! Store array of vertices
        do i = 1, N_verts
            do j = 1,3

                ! Store Locations
                vertices(i)%location(j) = vertex_locs(j,i)
            end do

            ! Store index
            vertices(i)%ind = i
        end do

        ! Get to start of polygons
        read(unit, '(a)') line
        do while (index(line, 'POLYGONS') == 0)
            read(unit,'(a)') line
        end do

        ! Deternine number of panels
        read(line,*) dummy_read, N_panels, dummy_read

        ! Allocate panel array
        allocate(panels(N_panels))

        ! Read in panels
        do i = 1, N_panels
            
            ! Get vertex indices
            read(unit,'(a)') dummy_read

            read(dummy_read,*) N, i1, i2, i3

            panels(i)%vert_ind = [i1+1,i2+1,i3+1]

            ! Calculate Normal
            call calc_norm(panels(i))

        end do

    ! Close vtk
    close(unit)


    freestream = [100,0,0]

    ! Normalize freestream vector
    freestream = freestream/(sqrt(freestream(1)**2 + freestream(2)**2 + freestream(3)**2))
    

    ! Calculate max pressure coefficent    
    c_pmax = (2/(gamma*m**2))*( (( (((gamma + 1)**2) * m**2 )/((4*gamma*m**2)-2*(gamma-1)))**(gamma/(gamma-1))) * &
        ((1-gamma+(2*gamma*m**2))/(gamma+1)) -1)
    

    ! Loop through and calculate pressure for all panels
    do i = 1, N_panels
        
        norm = panels(i)%normal
        write(*,*) "Panel Normal", norm
        ! Normalize normal vector
        norm = norm/(sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2))
        write(*,*) "Normalized", norm
        ! Calculate Panel inclination with respect to the freestream
        theta = acos(dot_product(freestream,norm))
        write(*,*) theta
        ! Check if panel faces freestream
        if ( (theta <= pi/2) .and. (theta >= -pi/2) ) then
            c_p = 0
            write(*,*) 'Panel does not face the freestream'
        else
            c_p = c_pmax*sin(theta - pi/2)**2
        end if

        panels(i)%cp = c_p
        write(*,*) panels(i)%cp
        write(*,*)
    end do



    ! Write out new vtk


contains

    ! Calculates the normal vector of a panel
    subroutine calc_norm(pan)
        type(panel), intent(inout) :: pan
        real, dimension(3) :: v1, v2, cross
        !write(*,*) pan%vert_ind
        !write(*,*) "1.", vertices(pan%vert_ind(1))%location
        ! Calculate panel vectors
        v1 = vertices(pan%vert_ind(3))%location - vertices(pan%vert_ind(1))%location
        v2 = vertices(pan%vert_ind(2))%location - vertices(pan%vert_ind(3))%location
        !write(*,*) "Vector 1", v1
        !write(*,*) "Vector 2", v2
        

        ! Calculate Cross Product
        cross(1) = v1(2) * v2(3) - v1(3) * v2(2)
        cross(2) = v1(3) * v2(1) - v1(1) * v2(3)
        cross(3) = v1(1) * v2(2) - v1(2) * v2(1)

        pan%normal = cross

        

    end subroutine calc_norm

end program main



