module vtk_mod


    use panel_mod
    use base_geom_mod

    contains

    subroutine load_surface_vtk(mesh_file, N_verts, N_panels, vertices, panels)
        ! Loads surface mesh from a vtk file.
        implicit none

        character(len=:), allocatable, intent(in) :: mesh_file
        integer, intent(out) :: N_verts, N_panels
        type(vertex), dimension(:), allocatable, intent(out) :: vertices
        type(panel), dimension(:), allocatable, intent(out) :: panels

        character(len=200) :: line
        integer :: ind, ver, unit

        ! Determine version
        open(newunit=unit, file=mesh_file)
        read(unit, '(a)') line
        ind = index(line, 'Version')
        read(line(ind+8:ind+8), *) ver
        close(unit)

        ! Load Based on version
        select case (ver)

        case (3)
            call load_surface_vtk_ver_3(mesh_file, N_verts, N_panels, vertices, panels)

        case default
            write(*,*) "!!! VTK file version ", ver, " not recognized. Quitting..."
            stop
            
        end select
        
    end subroutine load_surface_vtk

    subroutine load_surface_vtk_ver_3(mesh_file, N_verts, N_panels, vertices, panels)
        ! Loads surface mesh from vtk ver 3
        implicit none

        character(len=:), allocatable, intent(in) :: mesh_file
        integer, intent(out) :: N_verts, N_panels
        type(vertex), dimension(:), allocatable, intent(out) :: vertices
        type(panel), dimension(:), allocatable, intent(out) :: panels

        character(len=200) :: dummy_read, line
        real,dimension(:,:), allocatable :: vertex_locs
        integer :: i, j, N=3, i1, i2, i3, unit

        
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
            call panel_calc_norm(panels(i), vertices)

            ! Calculate Panel Area
            panels(i)%area = panels(i)%norm_mag/2

        end do

        ! Close vtk
        close(unit)


    end subroutine load_surface_vtk_ver_3

    subroutine vtk_out_write(body_file, N_verts, N_panels, vertices, panels)

        implicit none

        character(len=:), allocatable, intent(in) :: body_file
        integer, intent(in) :: N_verts, N_panels
        type(vertex), dimension(:), allocatable, intent(in) :: vertices
        type(panel), dimension(:), allocatable, intent(in) :: panels
        

        integer :: unit, i, j
        real :: sep
        
        open(newunit=unit, file=body_file)

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
            write(unit, '(e20.12)') panels(i)%c_p
        end do

        ! Write surface mach number
        write(unit, '(a, a, a)') "SCALARS ", "Mach_surface", " float 1"
        write(unit, '(a)') "LOOKUP_TABLE default"
        do i = 1, N_panels
            write(unit, '(e20.12)') panels(i)%m_surf
        end do

        ! Write surface mach number
        write(unit, '(a, a, a)') "SCALARS ", "Seperation", " float 1"
        write(unit, '(a)') "LOOKUP_TABLE default"
        do i = 1, N_panels
            if (panels(i)%seperated) then
                sep = 1.
                write(unit, '(e20.12)') sep
            else
                sep = 0
                write(unit, '(e20.12)') sep
            end if
        end do

        ! Write dC_f
        write(unit, '(a, a, a)') "VECTORS ", "dC_f", " float"
        do i = 1, N_panels
            write(unit, '(e20.12, e20.12, e20.12)') panels(i)%dC_f(1), panels(i)%dC_f(2), panels(i)%dC_f(3)
        end do

        ! Write inclination angle
        write(unit, '(a, a, a)') "SCALARS ", "theta", " float 1"
        write(unit, '(a)') "LOOKUP_TABLE default"
        do i = 1, N_panels
            write(unit, '(e20.12)') panels(i)%theta
        end do

        close(unit)

    end subroutine vtk_out_write



end module vtk_mod