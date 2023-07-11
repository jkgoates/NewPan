module base_geom_mod

    use math_mod
    use linked_list_mod

    implicit none

    type vertex
        integer :: ind
        real, dimension(3) :: location
        real, dimension(3) :: V
        logical :: shadowed = .false.
        type(list) :: adjacent_vertices ! List of indices of vertices for the vertices which share an edge with this vertex
        type(list) :: adjacent_edges ! List of indices for the edges which touch this vertex
        type(list) :: panels ! List of indices for the panels which connect to this vertex

        contains

        procedure :: get_panels => vertex_get_panels

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
        procedure :: get_opposing_panel => edge_get_opposing_panel

    end type edge

contains

    function vertex_get_panels(this) result(panels)

        implicit none
        
        class(vertex), intent(inout) :: this
        type(list) :: panels

        panels = this%panels
        


    end function vertex_get_panels


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


    function edge_get_opposing_panel(this, i_panel) result(i_oppose)
        ! Returns the index of the panel opposite this one on the edge

        implicit none
        
        class(edge), intent(in) :: this
        integer, intent(in) :: i_panel

        integer :: i_oppose

        if(i_panel == this%panels(1)) then
            i_oppose = this%panels(2)
        else if (i_panel == this%panels(2)) then
            i_oppose = this%panels(1)
        else
            i_oppose = 0
        end if

    end function  edge_get_opposing_panel


end module base_geom_mod