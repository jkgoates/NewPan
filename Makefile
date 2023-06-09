
# make for NewPan. A Modified-Newtonian Panel Method 

# Deactivate implicit rules
.SUFFIXES:

# Directories
SRC_DIR = src
COM_DIR = common
BIN_DIR = bin

# List common files (ordered based on dependency)
COMMON_FILES = json.f95 json_xtnsn.f95 math.f90
COMMON_PATHS = $(addprefix $(COM_DIR)/, $(COMMON_FILES))

# List source files (ordered based on dependency)
SRC_FILES = flow.f90 base_geom.f90 panel.f90 mesh.f90 vtk.f90 surface_mesh.f90 panel_solver.f90
SRC_PATHS = $(addprefix $(SRC_DIR)/, $(SRC_FILES))

# Main
MAIN_PATH = src/main.f90

# Compiler
COMPILER = gfortran

# Flags
DEBUG_FLAGS = -fbounds-check -fbacktrace -g
OMP_FLAG = -fopenmp
FLAGS = -O2 -fdefault-real-8
#DISLIN_FLAGS = -luser32 -lgdi32 -lopengl32

# DISLIN Command 
#DISLIN_PATH = .\lib\dislin_win\dismg.a
#DISLIN_LINUX = -I ./lib/dislin_linux/gf -L ./lib/dislin_linux/ -ldislin

# Program name
PROGRAM = newpan.exe

# Default make
default:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH)

# Linux make
#linux:
#	$(COMPILER) $(FLAGS) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH) $(DISLIN_LINUX)
	
# Debug option
debug:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) $(DEBUG_FLAGS) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH) 

# Debug with all warnings
wall:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) $(DEBUG_FLAGS) -Wall -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH) 

# Debug with all warnings
debug-serial:
	$(COMPILER) $(FLAGS) $(DEBUG_FLAGS) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH)

# Serial compilation (without OpenMP)
serial:
	$(COMPILER) $(FLAGS) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH) 

# Cleanup
clean:
	rm -rf *.mod *.exe $(SRC_DIR)/*.mod $(COM_DIR)/*.mod