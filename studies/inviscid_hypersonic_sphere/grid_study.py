import numpy as np
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, './studies')

from case_running_functions import run_newpan, write_input_file

RERUN_NEWPAN = True
study_dir = "studies/inviscid_hypersonic_sphere/"
plot_dir = study_dir + "plots/"

def run_study(M, grid, method):
    """Runs a case for given Mach number and mesh density"""

    # Storage locations
    case_name = "M_{0}_{1}_{2}".format(round(M, 3), grid, method)
    mesh_file = "meshes/random_sphere_{0}.vtk".format(grid)
    results = study_dir + "results/"
    results_file = study_dir + "results/"+case_name+".vtk"
    report_file = study_dir + "reports/"+case_name+".json"



    # Specify methods
    if method == "mod-newton":
        windward_method = "modified-newtonian"
        leeward_method = "none"
    elif method == "mod-newton-pm":
        windward_method = "modified-newtonian"
        leeward_method = "prandtl-meyer"
    elif method == "kaufman":
        windward_method = "kaufman"
        leeward_method = "none"

    # Write out input file

    input_dict = {
        "flow": {
        "freestream_direction": [-1.0,0.0,0.0],
        "mach_number" : M,
        "gamma" : 1.4
        },
        "geometry" : {
            "file": "meshes/random_sphere_ultra_fine_sample_8.vtk",
            "reference": {
                "area": 3.14159268
            }
        },
        "solver": {
            "windward_method" : windward_method,
            "leeward_method" : leeward_method
        },
        "output": {
            "body_file": results_file,
            "report_file": report_file
        }
        
    }

    # Dump
    input_file = study_dir + "input.json"
    write_input_file(input_dict, input_file)

    # Run case
    report = run_newpan(input_file, run=RERUN_NEWPAN)

    # Pull out forces
    C_F = np.zeros(3)
    try:
        C_F[0] = report["total_forces"]["C_x"]
        C_F[1] = report["total_forces"]["C_y"]
        C_F[2] = report["total_forces"]["C_z"]
    except KeyError:
        C_F[0] = np.nan
        C_F[1] = np.nan
        C_F[2] = np.nan
    except TypeError:
        C_F[0] = np.nan
        C_F[1] = np.nan
        C_F[2] = np.nan

    return C_F

if __name__ == "__main__":
    # Retrieve experimental data
    exp_data = np.genfromtxt(study_dir + "/data/exp_data_3.csv", delimiter=',')
    M_exp = exp_data[:,0]
    CD_exp = exp_data[:,1]


    # Parameters
    grids = ['coarse', 'medium', 'fine']
    methods = ['mod-newton', 'mod-newton-pm', 'kaufman']
    Ms = list(M_exp[1:])
    CD = np.zeros((len(grids), len(Ms), 3))

    # Loop
    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):
            for k, method in enumerate(methods):

                C_F = run_study(M,grid,method)
                CD[i,j,k] = C_F[0]
    
    # Plot
    plt.figure()
    for i, grid in enumerate(grids):
        if i == 1:
            plt.plot(Ms, CD[i,:,0], 'kH', mfc='none', markersize=10-3*i, label="Newton")
            plt.plot(Ms, CD[i,:,1], 'ko', mfc='none', markersize=10-3*i, label="Newton-Prandtl")
            plt.plot(Ms, CD[i,:,2], 'k^', mfc='none', markersize=10-3*i, label="Kaufman")
        else:
            plt.plot(Ms, CD[i,:,0], 'kH', mfc='none', markersize=10-3*i)
            plt.plot(Ms, CD[i,:,1], 'ko', mfc='none', markersize=10-3*i)
            plt.plot(Ms, CD[i,:,2], 'k^', mfc='none', markersize=10-3*i)

    plt.plot(M_exp, CD_exp, 'k-', label='Bailey et al. (1971)', markersize=3)
    plt.xlabel('$M_\infty$')
    plt.ylabel('$C_d$')
    plt.legend()
    plt.savefig(plot_dir+"c_d_over_M.pdf")
    plt.savefig(plot_dir+"c_d_over_M.svg")
    plt.close()

