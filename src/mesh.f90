module mesh_mod

    use base_geom_mod
    use panel_mod
    use flow_mod

    implicit none
    
    type mesh
        integer :: N_verts, N_panels = 0
        type(vertex), allocatable, dimension(:) :: vertices
        type(panel), allocatable, dimension(:) :: panels
        
    
    contains
        
    
    end type mesh
contains
    
end module mesh_mod

