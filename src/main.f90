program main

    use json_mod
    use json_xtnsn_mod
    use panel_mod
    use base_geom_mod
    use vtk_mod
    use panel_solver_mod
    use surface_mesh_mod
    use flow_mod
    use math_mod

    implicit none


    ! Initialize variables
    real :: count_rate, runtime
    real, dimension(3) :: C_f=[0,0,0]
    integer :: unit, start_count, end_count, i_unit

    character(len=:), allocatable :: mesh_file, body_file, result_file, &
                                    report_file, windward_method, leeward_method
    character(100) :: input_file
    type(flow) :: freestream_flow
    type(surface_mesh) :: body_mesh
    type(panel_solver) :: solver
    type(json_file) :: input_json
    type(json_value), pointer :: flow_settings, &
                                 geom_settings, &
                                 solver_settings, &
                                 output_settings, &
                                 report_json, p_parent

    logical :: exists, found, verbose
    integer :: solver_stat


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

    ! Get verbose toggle
    call json_xtnsn_get(output_settings, 'verbose', verbose, default_value=.true.)

    if (verbose) then
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
        write(*,*) "Reading and analyzing surface mesh"
    end if

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
    call json_xtnsn_get(solver_settings, 'windward_method', windward_method, 'modified-newtonian')
    call json_xtnsn_get(solver_settings, 'leeward_method', leeward_method, 'none')

    ! Read Surface mesh
    call body_mesh%init(geom_settings)

    ! Initialize flow
    call freestream_flow%init(flow_settings)

    if (verbose) then
        write(*,*)
        write(*,*) "Initializing based on flow properties"
    end if

    ! Get result files
    call json_xtnsn_get(output_settings, 'body_file', body_file, 'none')

    ! Initialize panel solver
    call solver%init(solver_settings, body_mesh, freestream_flow)

    ! Calculate shadowing
    call body_mesh%find_shadowed_panels(freestream_flow)

    ! Solve
    call solver%solve(body_mesh)

    ! Calculate force coefficients
    call solver%calc_forces(body_mesh)
    write(*,'(a20)') "Force Coefficients:"
    write(*, '(a20, e14.6)') "  C_x:", body_mesh%C_f(1)
    write(*, '(a20, e14.6)') "  C_y:", body_mesh%C_f(2)
    write(*, '(a20, e14.6)') "  C_z:", body_mesh%C_f(3)
    write(*,*)

    ! Write Mesh info to json
    call json_value_create(p_parent)
    call to_object(p_parent, 'mesh_info')
    call json_value_add(report_json, p_parent)
    call json_value_add(p_parent, 'N_body_panels', body_mesh%N_panels)
    call json_value_add(p_parent, 'N_body_vertices', body_mesh%N_verts)
    nullify(p_parent)

    ! Write solver results
    call json_value_create(p_parent)
    call to_object(p_parent, 'total_forces')
    call json_value_add(report_json, p_parent)
    call json_value_add(p_parent, 'C_x', body_mesh%C_f(1))
    call json_value_add(p_parent, 'C_y', body_mesh%C_f(2))
    call json_value_add(p_parent, 'C_z', body_mesh%C_f(3))
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
    call vtk_out_write(body_file, body_mesh%N_verts, body_mesh%N_panels, body_mesh%vertices, body_mesh%panels)

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