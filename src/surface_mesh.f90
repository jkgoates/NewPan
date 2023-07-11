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
        integer :: N_edges
        type(edge), allocatable, dimension(:) :: edges
        !real, dimension(:), allocatable :: dC_f ! Forces
        real, dimension(3) :: C_f=[0,0,0] ! Forces
        real :: S_ref ! Reference parameters

        contains

        ! Initialization
        procedure :: init => surface_mesh_init
        procedure :: load_mesh_file => surface_mesh_load_mesh_file
        procedure :: locate_adjacent_panels => surface_mesh_locate_adjacent_panels
        procedure :: store_adjacent_vertices => surface_mesh_store_adjacent_vertices
        procedure :: check_panels_adjacent => surface_mesh_check_panels_adjacent
        procedure :: surface_fit => surface_mesh_surface_fit

        ! Flow dependent Initialization
        procedure :: define_dependent_direction => surface_mesh_define_dependent_direction

        ! Getters
        procedure :: get_avg_characteristic_panel_length => surface_mesh_get_avg_characteristic_panel_length
        
        ! Panel Shadowing
        procedure :: find_shadowed_panels => surface_mesh_find_shadowed_panels
        procedure :: project => surface_mesh_project
        procedure :: shadowed_vertices => surface_mesh_shadowed_vertices
        procedure :: panel_over_point => surface_mesh_panel_over_point

        ! Initialization Based on flow
        procedure :: init_with_flow => surface_mesh_init_with_flow



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

        ! Locate adjacent panels
        call this%locate_adjacent_panels()

        ! Determine surface element equations
        call this%surface_fit()


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

    subroutine surface_mesh_locate_adjacent_panels(this)
        ! Loops through panels to determine which are adjacent

        implicit none

        class(surface_mesh), intent(inout), target :: this

        integer :: i, j, i_edge, edge_index_i, edge_index_j
        integer,dimension(2) :: i_endpoints
        integer, dimension(this%N_panels*3) :: panel1, panel2, vertex1, vertex2, edge_index1, edge_index2

        write(*,'(a)',advance='no') "      Locating adjacent panels..."

        ! Loop through each panel
        this%N_edges = 0

        !$OMP parallel do schedule(dynamic) private(j, i_endpoints, i_edge, edge_index_i, edge_index_j)
        do i = 1, this%N_panels

            ! Loop through each potential neighbor
            neighbor_loop: do j = i+1, this%N_panels

                ! Check if we've found all neighbors for this panel
                if (all(this%panels(i)%abutting_panels /= 0)) exit neighbor_loop

                ! Check if these are abutting
                if (this%check_panels_adjacent(i, j, i_endpoints, edge_index_i, edge_index_j)) then

                    !$OMP critical
                    ! Update number of edges
                    this%N_edges = this%N_edges + 1
                    i_edge = this%N_edges

                    ! Store vertices being adjacent to one another
                    call this%store_adjacent_vertices(i_endpoints, i_edge)

                    ! Store adjacent panels and panel edges
                    edge_index1(i_edge) = edge_index_i
                    edge_index2(i_edge) = edge_index_j

                    ! Store information in arrays
                    panel1(i_edge) = i
                    panel2(i_edge) = j
                    vertex1(i_edge) = i_endpoints(1)
                    vertex2(i_edge) = i_endpoints(2)

                    !$OMP end critical

                end if

            end do neighbor_loop

        end do

        ! Check for panels abutting empty space and add those edges
        do i=1, this%N_panels

            ! Check for an edge with no abutting panel
            do j = 1, this%panels(i)%N
                if (this%panels(i)%abutting_panels(j) == 0) then

                    ! Get endpoint indices
                    i_endpoints(1) = this%panels(i)%get_vertex_index(j)
                    i_endpoints(2) = this%panels(i)%get_vertex_index(mod(j, this%panels(i)%N) + 1)

                    ! Set up an edge
                    this%N_edges = this%N_edges + 1
                    i_edge = this%N_edges
                    panel1(i_edge) = i
                    panel2(i_edge) = 0 ! Placeholder
                    vertex1(i_edge) = i_endpoints(1)
                    vertex2(i_edge) = i_endpoints(2)
                    edge_index1(i_edge) = j
                    edge_index2(i_edge) = 0

                    ! Store adjacent vertices
                    call this%store_adjacent_vertices(i_endpoints, i_edge)

                end if
            end do
        end do

        ! Allocate edge storage
        allocate(this%edges(this%N_edges))

        ! Initialize edges
        do i = 1, this%N_edges

            ! Initialize
            call this%edges(i)%init(vertex1(i), vertex2(i), panel1(i), panel2(i))

            ! Store more information
            this%edges(i)%edge_index_for_panel(1) = edge_index1(i)
            this%edges(i)%edge_index_for_panel(2) = edge_index2(i)
            this%panels(panel1(i))%edges(this%edges(i)%edge_index_for_panel(1)) = i
            if (panel2(i) <= this%N_panels .and. panel2(i) > 0) then
                this%panels(panel2(i))%edges(this%edges(i)%edge_index_for_panel(2)) = i
            end if

        end do

        write(*,'(a, i7, a)') "Done. Found ", this%N_edges, " edges."


    end subroutine surface_mesh_locate_adjacent_panels

    subroutine surface_mesh_store_adjacent_vertices(this, i_verts, i_edge)
        ! Stores that the given two vertices are adjacent

        implicit none
        
        class(surface_mesh), intent(inout) :: this
        integer, dimension(2), intent(in) :: i_verts
        integer, intent(in) :: i_edge

        ! Store that the vertices are adjacent
        if (.not. this%vertices(i_verts(1))%adjacent_vertices%is_in(i_verts(2))) then
            call this%vertices(i_verts(1))%adjacent_vertices%append(i_verts(2))
        end if
        if (.not. this%vertices(i_verts(2))%adjacent_vertices%is_in(i_verts(1))) then
            call this%vertices(i_verts(2))%adjacent_vertices%append(i_verts(1))
        end if

        ! Store that they touch this edge
        call this%vertices(i_verts(1))%adjacent_edges%append(i_edge)
        call this%vertices(i_verts(2))%adjacent_edges%append(i_edge)

    end subroutine surface_mesh_store_adjacent_vertices

    function surface_mesh_check_panels_adjacent(this, i, j, i_endpoints, edge_index_i, edge_index_j) result(adjacent)
        ! Checks whether  panels i and j are adjacent

        class(surface_mesh), intent(inout) :: this
        integer, intent(in) :: i, j
        integer, dimension(2), intent(out) :: i_endpoints
        integer, intent(out) :: edge_index_i, edge_index_j

        logical :: adjacent

        logical :: already_found_shared
        integer :: m, n, m1, n1, temp

        ! Initialize
        adjacent = .false.

        ! Initialize for this panel pair
        already_found_shared = .false.

        ! Check if the panels are abutting
        abutting_loop: do m = 1, this%panels(i)%N
            do n = 1, this%panels(j)%N

                ! Check if they have the same index
                if (this%panels(i)%get_vertex_index(m) == this%panels(j)%get_vertex_index(n)) then

                    ! Previously found a shared vertex, so we have abutting panels
                    if (already_found_shared) then

                        adjacent = .true.

                        ! Store the second shared vertex
                        i_endpoints(2) = this%panels(i)%get_vertex_index(m)

                        ! Check order; edge should proceed counterclockwise about panel i
                        if (m1 == 1 .and. m == this%panels(i)%N) then
                            temp = i_endpoints(1)
                            i_endpoints(1) = i_endpoints(2)
                            i_endpoints(2) = temp
                        end if

                        ! Store adjacent panels and panel edges
                        ! This stores the adjacent panels and edges according to the index of that edge
                        ! for the current panel

                        ! Store that i is adjacent to j
                        if ( (n1 == 1 .and. n == this%panels(j)%N) .or. (n == 1 .and. n1 == this%panels(j)%N) ) then
                            this%panels(j)%abutting_panels(this%panels(j)%N) = i
                            edge_index_j = this%panels(j)%N
                        else
                            n1 = min(n,n1)
                            this%panels(j)%abutting_panels(n1) = i
                            edge_index_j = n1
                        end if

                        ! Store that j is adjacent to i
                        if (m1 == 1 .and. m== this%panels(i)%N) then ! Nth edge
                            this%panels(i)%abutting_panels(m) = j
                            edge_index_i = m
                        else ! 1st or 2nd edge
                            this%panels(i)%abutting_panels(m1) = j
                            edge_index_i = m1
                        end if

                        return

                    ! First shared vertex
                    else

                        already_found_shared = .true.
                        i_endpoints(1) = this%panels(i)%get_vertex_index(m)
                        m1 = m
                        n1 = n

                    end if
                end if
            
            end do

        end do abutting_loop

    end function surface_mesh_check_panels_adjacent



    subroutine surface_mesh_init_with_flow(this, freestream)

        implicit none
        
        class(surface_mesh), intent(inout) :: this
        type(flow), intent(in) :: freestream

        real, dimension(3) :: v_inf
        integer :: i

        call this%define_dependent_direction(freestream)

    end subroutine surface_mesh_init_with_flow



    subroutine surface_mesh_find_shadowed_panels(this, freestream)

        implicit none
        
        class(surface_mesh), intent(inout) :: this
        type(flow), intent(in) :: freestream
        write(*,*)
        write(*,'(a)',advance='no') " Finding shadowed panels..."

        ! Project the vertices and centroids onto the freestream
        call this%project(freestream)

        ! Find the shadowed vertices
        call this%shadowed_vertices(freestream)

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



    subroutine surface_mesh_shadowed_vertices(this, freestream)

        implicit none
        
        class(surface_mesh), intent(inout) :: this
        type(flow), intent(in) :: freestream

        integer :: i, j, k, l
        logical :: check_x, check_y, check_z
        
        ! Calculate side vectors for all panels
        do i = 1, this%N_panels
            this%panels(i)%s12_pro = this%projected_verts(this%panels(i)%get_vertex_index(2))%location &
                                    -this%projected_verts(this%panels(i)%get_vertex_index(1))%location

            this%panels(i)%s23_pro = this%projected_verts(this%panels(i)%get_vertex_index(3))%location &
                                    -this%projected_verts(this%panels(i)%get_vertex_index(2))%location

            this%panels(i)%s31_pro = this%projected_verts(this%panels(i)%get_vertex_index(1))%location &
                                    -this%projected_verts(this%panels(i)%get_vertex_index(3))%location

        end do
        
        panel_loop: do i = 1, this%N_panels
            ! Calculate panel limits
            call this%panels(i)%calc_projected_outline(this%projected_verts)

            ! Only panels facing the freestream can shadow vertices
            if (this%panels(i)%theta <= 0) cycle panel_loop

            ! Loop through vertices
            vertex_loop: do j = 1, this%N_verts

                ! Skip if the vertex is already shadowed
                if (this%projected_verts(j)%shadowed) cycle vertex_loop

                ! A panel cannot shadow it's own vertices
                !!! This check needs to be refined
                if (this%panels(i)%get_vertex_index(1) == j) cycle vertex_loop
                if (this%panels(i)%get_vertex_index(2) == j) cycle vertex_loop
                if (this%panels(i)%get_vertex_index(3) == j) cycle vertex_loop

                ! Check if outlines intersect
                check_x = (this%panels(i)%x_min_pro-0.01 > this%projected_verts(j)%location(1)) .or. &
                          (this%panels(i)%x_max_pro+0.01 < this%projected_verts(j)%location(1))

                if (check_x) cycle vertex_loop

                check_y = (this%panels(i)%y_min_pro-0.01 > this%projected_verts(j)%location(2)) .or. &
                          (this%panels(i)%y_max_pro+0.01 < this%projected_verts(j)%location(2))

                if (check_y) cycle vertex_loop

                check_z = (this%panels(i)%z_min_pro-0.01 > this%projected_verts(j)%location(3)) .or. &
                          (this%panels(i)%z_max_pro+0.01 < this%projected_verts(j)%location(3))

                if (check_z) cycle vertex_loop

                ! Panels can only shadow vertices downstream of themselves
                !!! This check needs to be refined
                if (dot_product(this%panels(i)%centr,freestream%v_inf) > &
                    dot_product(this%vertices(j)%location,freestream%v_inf)) then
                    cycle vertex_loop
                end if

                ! If to this point, see if the panel shadows the vertex

                this%projected_verts(j)%shadowed = this%panel_over_point(this%panels(i), this%projected_verts(j))

            end do vertex_loop

        end do panel_loop


        ! Loop through panels and see how many of their vertices are shadowed
        do i = 1, this%N_panels

            ! Loop through the vertices
            do j = 1,3
                if(this%projected_verts(this%panels(i)%get_vertex_index(j))%shadowed) then
                    this%panels(i)%shadowed = this%panels(i)%shadowed + 1
                end if

            end do
            
        end do

    end subroutine surface_mesh_shadowed_vertices


    function surface_mesh_panel_over_point(this, pan, vert) result(shadowed)
        
        implicit none
        
        class(surface_mesh), intent(inout) :: this
        type(panel), intent(inout) :: pan
        type(vertex), intent(inout) :: vert

        logical :: shadowed
        real, dimension(:), allocatable :: side, base_point, b, proj, point, centr, proj_point, err, c
        real, dimension(:,:), allocatable :: P, side_vec, side_vec_T, b_vec, proj_vec
        integer :: i, j, k
        real :: left_right
        integer :: lefts


        ! Get the centroid of the panel
        centr = pan%centr_pro

        lefts = 0
        ! Each side of the panel
        sides: do j = 1, 3

            ! Calculate the projection matrix for the side
            select case(j)
            case (1)
            side = pan%s12_pro
            case (2)
            side = pan%s23_pro
            case (3)
            side = pan%s31_pro
            end select

            side_vec = reshape(side, [3,1])
            side_vec_T = transpose(side_vec)

            P = matmul(side_vec, side_vec_T)/dot_product(side, side)

            base_point = this%projected_verts(pan%get_vertex_index(j))%location

            point = vert%location
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
            if (left_right >= 0) lefts = lefts + 1


        end do sides

        if (lefts == 3) then
            shadowed = .true.
        else if (lefts /= 3) then
            shadowed = .false.
        end if

    end function surface_mesh_panel_over_point


    subroutine surface_mesh_surface_fit(this)
        implicit none
        
        class(surface_mesh), intent(inout) :: this

        integer :: i

        ! Loop through panels and calc surface fit
        do i = 1, this%N_panels
            call this%panels(i)%calc_surface_equation()

        end do

    end subroutine surface_mesh_surface_fit


    subroutine surface_mesh_define_dependent_direction(this,freestream)
        implicit none

        class(surface_mesh),intent(inout) :: this
        type(flow), intent(in) :: freestream

        integer :: i

        ! Loop through panels and define the independent and dependent directions of the panel
        do i = 1, this%N_panels
            call this%panels(i)%define_dependent_direction(freestream)

        end do
    end subroutine surface_mesh_define_dependent_direction

    function surface_mesh_get_avg_characteristic_panel_length(this) result(l_avg)
        ! Returns the average characteristic length of the panels

        implicit none
        
        class(surface_mesh), intent(in):: this

        real :: l_avg

        integer :: i

        ! Loop through panels
        l_avg = 0.
        do i = 1, this%N_panels
            l_avg = l_avg + this%panels(i)%get_characteristic_length()
        end do
        l_avg = l_avg/this%N_panels

    end function surface_mesh_get_avg_characteristic_panel_length

end module surface_mesh_mod