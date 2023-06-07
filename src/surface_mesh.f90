module surface_mesh_mod
    
    use json_mod
    use json_xtnsn_mod
    use base_geom_mod
    use panel_mod
    use vtk_mod
    use mesh_mod
    
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

end module surface_mesh_mod