module base_geom_mod

    use math_mod
    use linked_list_mod

    implicit none

    type vertex
        integer :: ind
        real, dimension(3) :: location
        logical :: shadowed = .false.
        type(list) :: adjacent_vertices ! List of indices of vertices for the vertices which share an edge with this vertex
        type(list) :: adjacent_edges ! List of indices for the edges which touch this vertex
        type(list) :: panels ! List of indices for the panels which connect to this vertex

    end type vertex

    type vertex_pointer
        type(vertex), pointer :: ptr
    end type vertex_pointer

    type edge
        ! A mesh edge

        integer, dimension(2) :: top_verts ! Indices of the end vertices in the mesh vertex array belonging to the top panel
        integer, dimension(2) :: bot_verts ! Indices of the end vertices in the mesh vertex array belonging to the top panel
        integer, dimension(2) :: panels ! Indices of the top and bottom panels
        integer, dimension(2) :: edge_index_for_panel

        contains

        procedure :: init => edge_init

    end type edge

contains

    subroutine edge_init(this, i1, i2, top_panel, bottom_panel)

        implicit none
        
        class(edge), intent(inout) :: this
        integer, intent(in) :: i1, i2
        integer, intent(in) :: top_panel, bottom_panel

        ! Store indices
        this%top_verts(1) = i1
        this%top_verts(2) = i2

        ! Store panels
        this%panels(1) = top_panel
        this%panels(2) = bottom_panel

        ! Set defaults
        this%bot_verts = this%top_verts

    end subroutine edge_init

end module base_geom_mod