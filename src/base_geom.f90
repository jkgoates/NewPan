module base_geom_mod

    implicit none

    type vertex
        integer :: ind
        real, dimension(3) :: location
        logical :: shadowed = .false.
    end type vertex

    type vertex_pointer
        type(vertex), pointer :: ptr
    end type vertex_pointer

end module base_geom_mod