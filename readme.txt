1) Install miniconda or anaconda. 


2) From anaconda prompt create environment and install packages 
   For osx and linux:

    conda create -n LSEVPFFT
    conda activate LSEVPFFT
    conda install gfortran -c conda-forge  
    conda install lapack -c conda-forge
    conda install fftw -c conda-forge
    conda install make -c conda-forge
    conda install cxx-compiler -c conda-forge

    Last step is not really needed but helps set the environment variables properly.

   For Win:

    conda create -n LSEVPFFT
    conda activate LSEVPFFT
    conda install make -c conda-forge
    conda install lapack -c conda-forge
    conda install fftw -c conda-forge

    Appropriate version of gfortran can be installed as follows.
    a) Download archive (tested working version)
         GCC 8.5.0 + MinGW-w64 9.0.0 (MSVCRT) - release 1
           - Win64: 7-Zip archive* | Zip archive
       from https://winlibs.com/.
    b) Unzip to desired folder.
    c) Add path folder\mingw64\bin to User Enviroment Variable named Path.
   
3) Adjust the FFTW3 variable in makefile to be directory of your enviroment in conda. 
   LAPACK variable in makefile may be left blank or may need to be specified as well (depending on the system).

4) Compile LS-EVPFFT (serial, serial debug, openmp parallelized) using following commands
    make -B LS-EVPFFT
    make -B LS-EVPFFT-debug
    make -B LS-EVPFFT USE_OPENMP=1

5) Run serial
   osx and linux:
    ./LS-EVPFFT
   Win:
    LS-EVPFFT.exe

  Run openmp version
   osx and linux:
    ./LS-EVPFFT --nthreads n
   Win:
    LS-EVPFFT.exe --nthreads n
  where n is number of threads.
  
A thourough codebase reference can be found in docs/codebase_reference.md
