program main
   
    implicit none

    ! Panel Type
    type :: panel
        real, dimension(3) :: normal, dC_f
        integer, dimension(3) :: vert_ind
        real :: cp, norm_mag, area
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
    real :: gamma=1, m=100, c_pmax, c_p, pi=3.1415926, theta
    real, dimension(3) :: freestream, norm, C_f
    integer :: i, j, N=3, i1, i2, i3, unit, N_panels, N_verts
    type(panel), dimension(:), allocatable :: panels
    type(vertex), dimension(:), allocatable :: vertices
    character(len=:), allocatable :: mesh_file, result_file
    character(len=200) :: dummy_read, line
    real :: ref_area
    real, dimension(:,:), allocatable :: vertex_locs

    
    mesh_file = "meshes/random_sphere_ultra_fine_sample_8.vtk"
    ref_area = pi*1**2

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

            ! Calculate Magnitude of normal
            panels(i)%norm_mag = sqrt(sqrt(panels(i)%normal(1)**2 + panels(i)%normal(2)**2 + panels(i)%normal(3)**2))
            

            ! Calculate Panel Area
            panels(i)%area = panels(i)%norm_mag/2

        end do

    ! Close vtk
    close(unit)


    freestream = [-100,0,0]

    ! Normalize freestream vector
    freestream = freestream/(sqrt(freestream(1)**2 + freestream(2)**2 + freestream(3)**2))
    

    ! Calculate max pressure coefficent    
    c_pmax = (2/(gamma*m**2))*( (( (((gamma + 1)**2) * m**2 )/((4*gamma*m**2)-2*(gamma-1)))**(gamma/(gamma-1))) * &
        ((1-gamma+(2*gamma*m**2))/(gamma+1)) -1)
    

    ! Loop through and calculate pressure for all panels
    do i = 1, N_panels
        
        norm = panels(i)%normal
        
        ! Normalize normal vector
        norm = norm/(sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2))
        
        ! Calculate Panel inclination with respect to the freestream
        theta = acos(dot_product(freestream,norm))
        
        ! Check if panel faces freestream
        if ( (theta <= pi/2) .and. (theta >= -pi/2) ) then
            c_p = 0
        else
            c_p = c_pmax*sin(theta - pi/2)**2
        end if

        panels(i)%cp = c_p
    end do

    ! Loop through panels and calculate discrete force coefficients
    C_f = [0,0,0]
    do i = 1, N_panels

        panels(i)%dC_f = panels(i)%cp * panels(i)%normal * panels(i)%area/ panels(i)%norm_mag
        C_f = C_f + panels(i)%dC_f

    end do

    ! Make C_f non-dimensional
    C_f = C_f / ref_area
    write(*,*) "Force Coefficients", C_f


    ! Write out new vtk
    result_file = "results/random_sphere_ultra_fine_MN_mach10.vtk"

    open(newunit=unit, file=result_file)

    ! File header
    write(unit,'(a)') "# vtk DataFile Version 3.0"
    write(unit,'(a)') "NewPan results file."
    write(unit,'(a)') "ASCII"

    ! Vertex header
    write(unit,'(a)') "DATASET POLYDATA"
    write(unit, '(a i20 a)') "POINTS", N_verts, " float"

    ! Write out vertices
    do i=1,N_verts
        write(unit, '(e20.12, e20.12, e20.12)') vertices(i)%location(1), vertices(i)%location(2), vertices(i)%location(3)
    end do

    ! Panel Header
    write(unit, '(a i20 i20)') "POLYGONS", N_panels, N_panels*4

    ! Write out Panels
    do i=1,N_panels
        write(unit, '(i20)', advance='no') 3

        do j=1,panels(i)%N
            write(unit, '(i20)', advance='no') panels(i)%vert_ind(j) -1
        end do
        
        write(unit,*)
    end do

    ! Cell Data header
    write(unit, '(a i20)') "CELL_DATA", N_panels

    ! Write out cell normals
    write(unit, '(a)') "NORMALS normals float"
    do i=1,N_panels

        write(unit, '(e20.12, e20.12, e20.12)') panels(i)%normal(1), panels(i)%normal(2), panels(i)%normal(3)

    end do

    ! Write C_p
    write(unit, '(a, a, a)') "SCALARS ", "C_p", " float 1"
    write(unit, '(a)') "LOOKUP_TABLE default"
    do i = 1, N_panels
        write(unit, '(e20.12)') panels(i)%cp
    end do

    ! Write dC_f
    write(unit, '(a, a, a)') "VECTORS ", "dC_f", " float"
    do i = 1, N_panels
        write(unit, '(e20.12, e20.12, e20.12)') panels(i)%dC_f(1), panels(i)%dC_f(2), panels(i)%dC_f(3)
    end do
    

    close(unit)





contains

    ! Calculates the normal vector of a panel
    subroutine calc_norm(pan)
        implicit none
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



