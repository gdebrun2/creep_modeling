# ******* Edit these lines *******
# Fortran compiler
FC := gfortran

# LAPACK location (Example: "C:/msys64/mingw64/bin"). This line may be left blank on some systems.
LAPACK_Windows := /Users/gtdebru/miniconda3/envs/LSEVPFFT/lib
LAPACK_Darwin := /Users/gtdebru/miniconda3/envs/LSEVPFFT/lib
LAPACK_Linux := /home/gtdebru/miniconda3/envs/lsevpfft/lib

# FFTW location (Example: "D:/LosAlamos/NewVersion/fftw339")
FFTW3_Windows := /Users/gtdebru/miniconda3/envs/LSEVPFFT/lib
FFTW3_Darwin := /Users/gtdebru/miniconda3/envs/lsevpfft
FFTW3_Linux := /home/gtdebru/miniconda3/envs/LSEVPFFT/lib
# ********************************

# executable base threads
PROGRAM := LS-EVPFFT
SRC_DIR := src
BUILD_DIR := build

help:
	@echo
	@echo   Usage: make [program] [FLAGS]
	@echo
	@echo   Available programs:
	@echo     $(PROGRAM)
	@echo     $(PROGRAM)-debug
	@echo     $(PROGRAM)-test
	@echo
	@echo   make option USE_OPENMP=1 compiles with OpenMP parallelization
	@echo   See README.md for guidance.

# Detect operating system
ifeq ($(OS),Windows_NT)
  detected_OS := Windows
else
  detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

OS:
	@echo $(detected_OS)

# Detect operating system
ifeq ($(detected_OS),Windows) 
  LAPACK := $(LAPACK_Windows)
  FFTW3 := $(FFTW3_Windows)
endif
ifeq ($(detected_OS),Darwin) 
  LAPACK := $(LAPACK_Darwin)
  FFTW3 := $(FFTW3_Darwin)
endif
ifeq ($(detected_OS),Linux) 
  LAPACK := $(LAPACK_Linux)
  FFTW3 := $(FFTW3_Linux)
endif

# Check for fftw folder
ifeq ("$(wildcard $(FFTW3))","")
  $(error Directory "$(FFTW3)" does not exist, check FFTW3 flag in makefile)
endif

# Check for fftw3.f03
ifeq ("$(wildcard $(SRC_DIR)/fftw3.f03)","")
  $(error Missing file "fftw3.f03" - copy the file to this directory from your fftw download (/fftw-3.x.x/api/fftw3.f03))
endif

# Change settings depending on OS
ifeq ($(detected_OS),Windows)
  EXE := .exe
else
  EXE :=
endif
# Remove files
# NOTE:
#   We intentionally keep *.mod files after compilation so that editors/
#   language servers (VS Code + Modern Fortran/fortls) can resolve `use xyz`
#   statements. Removing module files makes the code compile fine, but breaks
#   IntelliSense/linting.
ifeq ($(detected_OS),Windows)
  RF := del /q $(BUILD_DIR)\\*.o 2>nul
else
  RF := rm -f $(BUILD_DIR)/*.o
endif

# If you really want a fully clean build directory (including *.mod), use:
#   make clean_mod

# openmp make option
ifdef USE_OPENMP
  OMPFLAGS := -fopenmp 
  ifeq ($(detected_OS),Windows)
	FFTWOMP :=
  else
	FFTWOMP := -lfftw3_threads -lm
  endif
  PARALLEL_MSG := OpenMP parallelized 
else
  OMPFLAGS :=
  FFTWOMP :=
  PARALLEL_MSG := serial 
endif

# Setup flags
ifeq ($(strip $(LAPACK)),)
  LAPACKFLAGS := -llapack
else
  LAPACKFLAGS := -L$(LAPACK) -llapack
endif
ifeq ($(detected_OS),Windows)
  FFTWFLAGS := -L${FFTW3} -lfftw3 -lm $(FFTWOMP)
else
  FFTWFLAGS := -L${FFTW3}/lib -lfftw3 -lm $(FFTWOMP)
endif

# Build options
GENERALFLAGS := -m64 -fdefault-real-8 -fdefault-double-8 -fno-align-commons
MODFLAGS := -J$(BUILD_DIR) -I$(BUILD_DIR)
LIBRARYFLAGS := $(LAPACKFLAGS) $(FFTWFLAGS) $(OMPFLAGS)
DEBUGFLAGS := -Wall -Wline-truncation -Wextra -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface \
            -Wunused-parameter -Wconversion -fwhole-file -pedantic -fbacktrace -fbounds-check -fcheck=all -ftrapv -O0 -C -g -fimplicit-none \
			-ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan -ffree-line-length-0
SPEEDFLAGS := -O3 -ffast-math -fno-protect-parens -fno-signed-zeros -fvariable-expansion-in-unroller -funroll-loops

# dependencies.  Note: order is important.
modules = \
$(SRC_DIR)/fftw3.f90 \
$(SRC_DIR)/fourier_functions.f90 \
$(SRC_DIR)/global.f90\
$(SRC_DIR)/types.f90 \
$(SRC_DIR)/tensor_functions.f90 \
$(SRC_DIR)/various_functions.f90 \
$(SRC_DIR)/IO_functions.f90 

all : $(PROGRAM) $(PROGRAM)-debug $(PROGRAM)-test

## Compile the main version
$(PROGRAM) : $(modules) $(SRC_DIR)/LS-EVPFFT.f90 
	@echo Compiling $(PARALLEL_MSG)$@:
	@mkdir -p $(BUILD_DIR)
	@echo $(FC) $^ -o $@$(EXE) $(GENERALFLAGS) $(SPEEDFLAGS) $(MODFLAGS) $(LIBRARYFLAGS)
	@$(FC) $^ -o $@$(EXE) $(GENERALFLAGS) $(SPEEDFLAGS) $(MODFLAGS) $(LIBRARYFLAGS)
	@$(RF)

## Compile the debug version
$(PROGRAM)-debug : $(modules) $(SRC_DIR)/LS-EVPFFT.f90 
	@echo Compiling $(PARALLEL_MSG)$@:
	@mkdir -p $(BUILD_DIR)
	@echo $(FC) $^ -o $@$(EXE) $(GENERALFLAGS) $(MODFLAGS) $(LIBRARYFLAGS) $(DEBUGFLAGS)
	@$(FC) $^ -o $@$(EXE) $(GENERALFLAGS) $(MODFLAGS) $(LIBRARYFLAGS) $(DEBUGFLAGS)
	@$(RF)

## Unit test
$(PROGRAM)-test : $(modules) $(SRC_DIR)/unit_tests.f90 
	@echo Compiling $(PARALLEL_MSG)$@:
	@mkdir -p $(BUILD_DIR)
	@echo $(FC) $^ -o $@$(EXE) $(GENERALFLAGS) $(MODFLAGS) $(LIBRARYFLAGS)
	@$(FC) $^ -o $@$(EXE) $(GENERALFLAGS) $(MODFLAGS) $(LIBRARYFLAGS)
	@$(RF)

# Cleanup
clean :
	@$(RF)

clean_mod :
	@rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.mod 2>/dev/null || true
