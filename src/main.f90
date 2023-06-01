program main

    use json_mod
    use json_xtnsn_mod
    use panel_mod
    use base_geom_mod
    use vtk_mod
    use panel_solver_mod
    use flow_mod

    implicit none


    ! Initialize variables
    real :: gamma, m, c_pmax, pi=3.1415926, ref_area, count_rate, runtime
    real, dimension(3) :: C_f=[0,0,0]
    real, dimension(:), allocatable :: freestream
    integer :: i, j, unit, N_panels, N_verts, start_count, end_count, i_unit

    character(len=:), allocatable :: mesh_file, body_file, result_file, &
                                    report_file, windward_method, leeward_method
    character(100) :: input_file
    type(flow) :: flows
    type(panel), dimension(:), allocatable :: panels
    type(vertex), dimension(:), allocatable :: vertices
    type(json_file) :: input_json
    type(json_value), pointer :: flow_settings, &
                                 geom_settings, &
                                 solver_settings, &
                                 output_settings, &
                                 report_json, p_parent

    logical :: exists, found


    ! Initialize json
    call json_initialize()

    ! Get input file from command line
    call getarg(1, input_file)
    input_file = trim(input_file)

    ! Get input file from user
    if (input_file == '') then
        write(*,*) "Please specity an input file:"
        read(*,'(a)') input_file
        input_file = trim(input_file)
    end if

    ! Check it exists
    inquire(file=input_file, exist=exists)
    if (.not. exists) then
        write(*,*) "!!! the file ", input_file, " does not exist. Quitting..."
        stop
    end if
    
    ! Start Timer
    call system_clock(start_count, count_rate)

    ! Load settings from input file
    call input_json%load_file(filename=input_file)
    call json_check()
    call input_json%get('flow', flow_settings, found)
    call input_json%get('geometry', geom_settings, found)
    call input_json%get('solver', solver_settings, found)
    call input_json%get('output', output_settings, found)

    ! Welcome message
    write(*,*) ""
    write(*,*) "      ------------------------------"
    write(*,*) "     -----/  ___  ------------------------"
    write(*,*) "    -----/  /   \         " 
    write(*,*) "   -----|  |      \                  "
    write(*,*) "  ------|  |   _    \     NewPan (c) 2023 USU Aerolab  "
    write(*,*) " -------|  |  |_|    )               v1.0"   
    write(*,*) "  ------|  |        /               "
    write(*,*) "   -----|  |      /               "
    write(*,*) "    -----\  \___/        " 
    write(*,*) "     -----\       ------------------------"
    write(*,*) "      -----------------------------"
    write(*,*) ""
    write(*,*) "Got input file: ", input_file
    write(*,*) ""
    write(*,*) ""

    ! Initialize report JSON
    call json_value_create(report_json)
    call to_object(report_json, 'report')
    call json_value_create(p_parent)
    call to_object(p_parent, 'info')
    call json_value_add(report_json, p_parent)
    call json_value_add(p_parent, 'generated_by', 'NewPan (c) 2023 USU Aerolab')
    call json_value_add(p_parent, 'executed', fdate())
    nullify(p_parent)
    
    ! Obtain input settings
    call json_xtnsn_get(output_settings, 'body_file', body_file, 'none')
    call json_xtnsn_get(geom_settings, 'file', mesh_file, 'none')
    call json_xtnsn_get(flow_settings, 'gamma', gamma, 1.)
    call json_xtnsn_get(flow_settings, 'mach_number', m, 1.)
    call json_get(flow_settings, 'freestream_direction', freestream, found)
    call json_xtnsn_get(geom_settings, 'reference.area', ref_area, 1.)
    call json_xtnsn_get(solver_settings, 'windward_method', windward_method, 'modified-newtonian')
    call json_xtnsn_get(solver_settings, 'leeward_method', leeward_method, 'none')

    ! Read Surface mesh
    call load_surface_vtk(mesh_file, N_verts, N_panels, vertices, panels)


    ! Normalize freestream vector
    freestream = freestream/(sqrt(freestream(1)**2 + freestream(2)**2 + freestream(3)**2))

    ! Calculate max pressure coefficient
    call panel_solver_max_pressure(windward_method, gamma, m, c_pmax)

    ! Calculate panel angles
    call panel_solver_calc_angles(panels, freestream, pi, N_panels)
    
    ! Calculate Prandtl meyer seperation
    call panel_solver_calc_seperation(leeward_method, N_panels, panels, m, gamma)

    ! Calculate pressures
    call panel_solver_calc_pressures(freestream, gamma, leeward_method, N_panels, panels, m, pi, c_pmax)
    
    ! Calculate force coefficients
    call panel_solver_calc_forces(C_f, N_panels, panels, ref_area)
    write(*,'(a20)') "Force Coefficients:"
    write(*, '(a20, e14.6)') "  C_x:", C_f(1)
    write(*, '(a20, e14.6)') "  C_y:", C_f(2)
    write(*, '(a20, e14.6)') "  C_z:", C_f(3)
    write(*,*)

    ! Write Mesh info to json
    call json_value_create(p_parent)
    call to_object(p_parent, 'mesh_info')
    call json_value_add(report_json, p_parent)
    call json_value_add(p_parent, 'N_body_panels', N_panels)
    call json_value_add(p_parent, 'N_body_vertices', N_verts)
    nullify(p_parent)

    ! Write solver results
    call json_value_create(p_parent)
    call to_object(p_parent, 'total_forces')
    call json_value_add(report_json, p_parent)
    call json_value_add(p_parent, 'C_x', C_f(1))
    call json_value_add(p_parent, 'C_y', C_f(2))
    call json_value_add(p_parent, 'C_z', C_f(3))
    nullify(p_parent)

    ! Find Time
    call system_clock(end_count)
    runtime = real(end_count - start_count)/count_rate

    ! Write timing to json
    call json_value_create(p_parent)
    call to_object(p_parent, 'timing')
    call json_value_add(report_json,p_parent)
    call json_value_add(p_parent, 'runtime', runtime)
    nullify(p_parent)

    ! Write input to json
    call json_value_create(p_parent)
    call to_object(p_parent, 'input')
    call json_value_add(report_json, p_parent)
    call json_value_add(p_parent, flow_settings)
    nullify(p_parent)

    ! Write out new vtk
    call vtk_out_write(body_file, N_verts, N_panels, vertices, panels)

    ! Write report
    call json_xtnsn_get(output_settings, 'report_file', report_file, 'none')
    if (report_file /= 'none') then
        open(newunit=i_unit, file=report_file, status='REPLACE')
        call json_print(report_json, i_unit)
        close(i_unit)
        write(*,'(a20 a)') "    Report: ", report_file
        write(*,*)
    end if

    ! Goodbye
    write(*,'(a, f10.4, a)') " New Pan Execution Time", runtime, " s"

    write(*,*)

end program main