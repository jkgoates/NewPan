module surface_mesh_mod
    
    use json_mod
    use json_xtnsn_mod
    use base_geom_mod
    use panel_mod
    use vtk_mod
    use mesh_mod
    use flow_mod
    
    implicit none
    

    type, extends(mesh) :: surface_mesh

        real, dimension(:), allocatable :: C_p_str, C_p_mod, C_p_kau ! Pressure Coefficients
        !real, dimension(:), allocatable :: dC_f ! Forces
        real, dimension(3) :: C_f=[0,0,0] ! Forces
        real :: S_ref ! Reference parameters

        contains

        ! Initialization
        procedure :: init => surface_mesh_init
        procedure :: load_mesh_file => surface_mesh_load_mesh_file
        
        ! Panel Shadowing
        procedure :: find_shadowed_panels => surface_mesh_find_shadowed_panels
        procedure :: project => surface_mesh_project
        procedure :: possible_panel_shadowing => surface_mesh_possible_panel_shadowing
        procedure :: panels_determine_overlap => surface_mesh_panels_determine_overlap
        procedure :: calc_point_overlap => surface_mesh_calc_point_overlap



    end type surface_mesh

contains
    

    subroutine surface_mesh_init(this, settings)

        implicit none
        
        class(surface_mesh), intent(inout) :: this
        type(json_value), pointer, intent(in) :: settings

        character(len=:), allocatable :: mesh_file
        logical :: found

        ! Load mesh from file
        call json_get(settings, 'file', mesh_file)
        mesh_file = trim(mesh_file)
        call this%load_mesh_file(mesh_file)

        ! Store settings
        call json_xtnsn_get(settings, 'reference.area', this%S_ref, 1.)

    end subroutine surface_mesh_init


    subroutine surface_mesh_load_mesh_file(this, mesh_file)
        ! Loads the mesh from file

        implicit none
        
        class(surface_mesh), intent(inout) :: this
        character(len=:), allocatable, intent(in) :: mesh_file

        logical :: file_exists
        integer :: loc
        character(len=:), allocatable :: extension

        ! Check if mesh file exists
        inquire(file=mesh_file, exist=file_exists)
        if ( .not. file_exists ) then
            write(*,*) "!!! Mesh file ", mesh_file, " does not exist. Quitting..."
            stop
        end if

        ! Determine the type of mesh file
        loc = index(mesh_file, '.')
        extension = mesh_file(loc:len(mesh_file))

        ! Load mesh file
        write(*,'(a, a, a)', advance='no') "      Reading surface mesh in from file ", mesh_file, "..."

        select case (extension)

        case ('.vtk')
            call load_surface_vtk(mesh_file, this%N_verts, this%N_panels, this%vertices, this%panels)

        case default
            write(*,*) "Machline cannot read ", extension, " type mesh files. Quitting..."
            stop
        
        end select
        write(*,*) "Done."

        ! Display mesh info
        write(*,'(a, i7, a, i7, a)') "      Surface mesh has ", this%N_verts, " vertices and ", &
                                                    this%N_panels, " panels."
    end subroutine surface_mesh_load_mesh_file

    subroutine surface_mesh_find_shadowed_panels(this, freestream)

        implicit none
        
        class(surface_mesh), intent(inout) :: this
        type(flow), intent(in) :: freestream
        write(*,*)
        write(*,'(a)',advance='no') " Finding shadowed panels..."

        ! Project the vertices and centroids onto the freestream
        call this%project(freestream)

        ! Find possible intersections to simplify calculations
        call this%possible_panel_shadowing(freestream)

        ! Determine overlap
        call this%panels_determine_overlap()

        write(*,*) "Done."

    end subroutine surface_mesh_find_shadowed_panels

    subroutine surface_mesh_project(this, freestream)

        implicit none
        
        class(surface_mesh), intent(inout) :: this
        type(flow), intent(in) :: freestream

        real, dimension(:,:), allocatable :: v_inf_T, v_inf, v_inf_m, proj_mat, &
                                            temp_vert_1, temp_vert_2, temp_centr_1, temp_centr_2
        real, dimension(3,3) :: iden
        integer :: i, j


        v_inf = reshape(freestream%v_inf, [3,1])

        ! Take the transpose of the freestream and multiply
        v_inf_T = transpose(v_inf)
        v_inf_m = matmul(v_inf, v_inf_T)

        ! Initialize identity matrix
        do i = 1,3
            do j = 1, 3
                if (i == j) then
                    iden(i,j) = 1
                else
                    iden(i,j) = 0
                end if
            end do
        end do

        ! Projection matrix
        proj_mat = iden - v_inf_m

        ! Allocate projected vertices
        allocate(this%projected_verts(this%N_verts))

        ! Loop through vertices and project onto the freestream plane

        do i = 1, this%N_verts

            temp_vert_1 = reshape(this%vertices(i)%location, [3,1])
            temp_vert_2 = matmul(proj_mat, temp_vert_1)

            do j = 1, 3
                this%projected_verts(i)%location(j) = temp_vert_2(j,1)
            end do 

            this%projected_verts(i)%ind = this%vertices(i)%ind
        end do

        ! Loop through panels and project centroid onto the freestream plane
        do i = 1, this%N_panels

            temp_centr_1 = reshape(this%panels(i)%centr, [3,1])
            temp_centr_2 = matmul(proj_mat, temp_centr_1)

            do j = 1, 3
                this%panels(i)%centr_pro(j) = temp_centr_2(j,1)
            end do

        end do

    end subroutine surface_mesh_project

    subroutine surface_mesh_possible_panel_shadowing(this, freestream)
        
        implicit none
        
        class(surface_mesh), intent(inout) :: this
        type(flow), intent(in) :: freestream
        integer :: i,j,k,l
        logical :: check_1, check_2, check_3
        
        ! Calculate panel limits
        !$OMP parallel do schedule(static)
        do i = 1, this%N_panels
            call this%panels(i)%calc_projected_outline(this%projected_verts)
        end do

        !$OMP parallel do private(j,k,l,check_1,check_2,check_3) &
        !$OMP & schedule(dynamic) shared(this, freestream)
        do i = 1, this%N_panels

            ! Allocate shadowing array
            allocate(this%panels(i)%panel_shadowing(this%N_panels))
            do j=1,this%N_panels
                this%panels(i)%panel_shadowing(j) = 0.0
            end do
            
            ! Only panels facing the freestream can have a shadowing effect
            if (this%panels(i)%theta <= 0) cycle

            ! Loop through all other panels and find possible intersections
            secondary: do j = 1, this%N_panels

                ! Skip if the same panel
                if (j==i) cycle secondary

                ! Skip if already shadowed
                !if (this%panels(j)%shadowed) cycle secondary

                ! Panels that share a vertex cannot shadow each other
                do k = 1,3
                    do l = 1,3
                        if (this%panels(i)%vert_ind(k) == this%panels(j)%vert_ind(l)) then
                            cycle secondary
                        end if
                    end do
                end do

                ! Check if outlines intersect
                check_1 = this%panels(i)%x_min_pro > this%panels(j)%x_max_pro .or. &
                          this%panels(j)%x_min_pro > this%panels(i)%x_max_pro

                if (check_1) cycle secondary

                check_2 = this%panels(i)%y_min_pro > this%panels(j)%y_max_pro .or. &
                          this%panels(j)%y_min_pro > this%panels(i)%y_max_pro

                if (check_2) cycle secondary
                
                check_3 = this%panels(i)%z_min_pro > this%panels(j)%z_max_pro .or. &
                          this%panels(j)%z_min_pro > this%panels(i)%z_max_pro

                if (check_3) cycle secondary

                ! Panels can only shadow other panels downstream of themselves
                if (dot_product(this%panels(i)%centr,-freestream%v_inf) < &
                    dot_product(this%panels(j)%centr,-freestream%v_inf)) then
                    cycle secondary
                end if


                ! If the loop has gotten to this point, there's a chance that the panel i shadows panel j
                !$OMP critical
                this%panels(i)%panel_shadowing(j) = 1.0
                !this%panels(j)%shadowed = .true.
                !$OMP end critical



            end do secondary

        end do

    end subroutine surface_mesh_possible_panel_shadowing

    subroutine surface_mesh_panels_determine_overlap(this)

        implicit none
        
        class(surface_mesh), intent(inout) :: this

        integer ::i,j
        logical :: overlapped

        ! Calculate side vectors for all panels
        !$OMP parallel do schedule(static)
        do i = 1, this%N_panels
            this%panels(i)%s12_pro = this%projected_verts(this%panels(i)%vert_ind(2))%location &
                                    -this%projected_verts(this%panels(i)%vert_ind(1))%location

            this%panels(i)%s23_pro = this%projected_verts(this%panels(i)%vert_ind(3))%location &
                                    -this%projected_verts(this%panels(i)%vert_ind(2))%location

            this%panels(i)%s31_pro = this%projected_verts(this%panels(i)%vert_ind(1))%location &
                                    -this%projected_verts(this%panels(i)%vert_ind(3))%location

        end do

        !$OMP parallel do private(j,overlapped) &
        !$OMP & schedule(dynamic) shared(this)
        ! Loop through panels and calculate overlap with each possible panel
        do i = 1, this%N_panels
            ! Downstream facing panels cannot shadow
            if (this%panels(i)%theta <= 0.0) cycle

            secondary: do j = 1, size(this%panels(i)%panel_shadowing)
                ! Skip if the same panel 
                if (i == j) cycle secondary

                ! Skip if panel i does not possibly shadow panel j 
                if (this%panels(i)%panel_shadowing(j) == 0) cycle secondary

                ! Otherwise run the panel overlap calculations
                call this%calc_point_overlap(this%panels(i), this%panels(j), overlapped)


                !$OMP critical
                ! Check if panel is already overlapped
                if (overlapped) then
                    this%panels(j)%shadowed = .true.
                    if (j==5360) write(*,*) i
                end if
                !$OMP end critical


            end do secondary
        end do


    end subroutine surface_mesh_panels_determine_overlap

    subroutine surface_mesh_calc_point_overlap(this, pan_1, pan_2, overlapped)

        implicit none
        
        class(surface_mesh), intent(inout) :: this
        type(panel), intent(inout) :: pan_1, pan_2
        logical, intent(inout) :: overlapped

        real, dimension(:), allocatable :: side, base_point, b, proj, point, centr, proj_point, err, c
        real, dimension(:,:), allocatable :: P, side_vec, side_vec_T, b_vec, proj_vec
        type(panel) :: main_panel, test_panel
        integer :: i, j, k
        real :: left_right
        integer, dimension(3) :: lefts = 0
        overlapped = .false.

        ! Evaluate for each panel
        panels: do i = 1, 2
            ! Identify panels
            if (i == 1) then
                main_panel = pan_1
                test_panel = pan_2
            else
                main_panel = pan_2
                test_panel = pan_1
            end if

            lefts = 0

            ! Get the centroid of the main panel
            centr = main_panel%centr_pro

            ! Each side of the panel
            sides: do j = 1, 3

                ! Calculate the projection matrix for the side
                select case(j)
                case (1)
                side = main_panel%s12_pro
                case (2)
                side = main_panel%s23_pro
                case (3)
                side = main_panel%s31_pro
                end select

                side_vec = reshape(side, [3,1])
                side_vec_T = transpose(side_vec)

                P = matmul(side_vec, side_vec_T)/dot_product(side, side)

                base_point = this%projected_verts(main_panel%vert_ind(j))%location

                ! Each point on the other panel
                points: do k = 1,3

                    point = this%projected_verts(test_panel%vert_ind(k))%location
                    b = point - base_point
                    b_vec = reshape(b, [3,1])
                    
                    ! Project point onto side
                    proj_vec = matmul(P,b_vec)
                    proj = reshape(proj_vec, [3])
                    proj_point = proj + base_point

                    ! Determine error and centroid vector
                    err = point - proj_point
                    c = centr - proj_point

                    ! Do the left right check
                    left_right = dot_product(err,c)

                    ! Update the lefts array
                    if (left_right > 0) lefts(k) = lefts(k) + 1


                end do points

            end do sides

            if(lefts(1) == 3) overlapped = .true.
            if(lefts(2) == 3) overlapped = .true.
            if(lefts(3) == 3) overlapped = .true.

        end do panels

    end subroutine surface_mesh_calc_point_overlap



end module surface_mesh_mod