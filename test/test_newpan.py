import os
import json
import shutil
import subprocess as sp


class NewPanError(Exception):
    pass


def run_newpan(input_file, remove_input=False, remove_results=True):
    # Runs NewPan and delivers the output; cleans up output files

    # Create results directory
    if not os.path.exists("test/results"):
        os.mkdir("test/results")

    # Run NewPan
    result = sp.run(["./newpan.exe", input_file], capture_output=True, text=True)

    # Remove input
    if remove_input:
        os.remove(input_file)

    # Read in report
    report_file = "test/results/report.json"
    if os.path.exists(report_file):
        with open(report_file, 'r') as report_handle:
            report = json.load(report_handle)

        # Check if  NewPan Thinks it was successful
        success = "NewPan Execution Time" in result.stdout

    else:
        
        # If no report file was generated, NewPan was not successful
        success = False

    # Clean up results directory
    if remove_results:
        shutil.rmtree("test/results")

    if not success:
        print(result.stdout)
        print(result.stderr)
        raise NewPanError
    
    else:
        Cx = float(report['total_forces']['C_x'])
        Cy = float(report['total_forces']['C_y'])
        Cz = float(report['total_forces']['C_z'])
        return Cx, Cy, Cz
    

def test_sphere_modified_newtonian_prandtl_meyer_from_free():
    # Tests the sphere case with the modified newtonian on the windward side and a prandtl-meyer expansion from the freestream

    Cx, Cy, Cz = run_newpan("test/input_files/sphere_input.json")

    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(Cx - 0.939448189087236) < 1e-12)
    assert(abs(Cy) < 1e-12)
    assert(abs(Cz) < 1e-12)

def test_ogive_kaufman():
    # Tests the ogive case using kaufman's method

    Cx, Cy, Cz = run_newpan("test/input_files/ogive_input.json")

    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(Cx - 0.12226021493975) < 1e-12)
    assert(abs(Cy - -1.82718119920202e-09) < 1e-12)
    assert(abs(Cz - -1.82669564475116e-09) < 1e-12)

def test_biconvex_airfoil_straight_newton():
    # Tests the biconvex airfoil case with the straight newtonian method and the prandtl-meyer expansion

    Cx, Cy, Cz = run_newpan("test/input_files/biconvex_airfoil_input.json")

    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(Cx - 0.029979455043952) < 1e-12)
    assert(abs(Cy - -1.2035980047744e-07) < 1e-12)
    assert(abs(Cz - 6.75945685934791e-08) < 1e-12)
    