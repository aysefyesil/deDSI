# deDSI

# Description: 
deDSI is a highly parallel hybrid C++/Python library that creates a dynamical system interface between dedalus (http://dedalus-project.org/) and the NSolver of Channelflow (https://www.channelflow.ch/). 

This library aims to find the solutions ( equilibria, traveling waves and periodic orbits) of non-linear partial differential equations, continue such solutions in parameter space and calculate the eigenvalues of such solutions using Arnoldi iterations. 

For the equation definition and the time integration, deDSI utilizes Dedalus library. For aforementioned dynamical system calculations deDSI utilizes the highly optimize NSolver library.

# Release Notes:

In the current version, Python part of the library is aimed to be a template for the user to easily modify it according to their own equations at hand. The equations provided with this version are Busse's Annulus equations. Also, in the current version 1-D and 2-D equations can be solved, however following the example of the library, functions that handle the third dimension can be implemented easily.

# Installation: 

Operating system: Linux, Unix

#### Following packages should be installed before installing deDSI

 dedalus 

 Channelflow [Nsolver comes with channelflow]

 Boost 3.5  

#### Additional to the installing and linking the above packages, following packages has to be explicitly linked:

 Python 3.5 (or the same python version with dedalus)

 Numpy 3.5

 Eigen 

a Makefile is provided for the guidance. 

# Starting & Contents:

Using this library starts with the implementation (dedalus) of listed functions in the interface. Implementation should be done according to the equation that user needs to solve. Guidelines are provided as follows:

**1.** Python functions that are mandatory are listed in the dedalus interface: ded/de_module_ns_dsi.py

**2.** User should use the same name of the functions listed in the interface and implement their own functions for their equations.

**3.** An example implementation of those functions for Busse's Annulus model can be found under the folder: ded/busse_ann/ 

**4.** Makefile should be run

**5.** Executables will be created under the folder: ./execs/

**6.** User will find the following executables:

 * execs/findsolnx      -> for finding the solutions (roots)

 * execs/continuationx  -> for continuing the solutions in parameter space

 * execs/findeigenvalsx -> for finding the eigenvalues of given solution
 
**7.** Some post-processing functions are also provided under ded/utils/

# Usage:

## Help

For all the executables above, user can use the help option as:

./findsolnx -h

Help provides all possible algorithmic, solver and optimization parameters and their default values.

## Examples:
All following examples are for the Busse's Annulus system. For more detailed mathematical explanations check channelflow webpage, for the explanations of the flags please refer the help flag.

### 1. findsolnx

Commandline arguments for searching for solutions (roots) of the given equation.

#### a. equilibrium solution
execs/findsolnx -eqb -fname solution_candidate.h5 -valRa 4200 -valBe 3.2E5 -valCc 0 -valPr 1 -T 0.01 -es 1e-13 -d 100 -dmax 1000

#### b. travelling wave
travelling wave solution can be searched through the periodic dimension, here x (-xrel -ax), starting values can be written :

execs/findsolnx -eqb -fname solution_candidate.h5 -valRa 1.675e7 -valBe 3.16E5 -valCc 0 -valPr 1 -valLx 0.6283185307179586 -Nx 64 -Ny 64 -T 0.0203 -es 1e-12 -d 1 -tstep 1e-5 -xrel -ax ax_candidate.asc

#### c. periodic orbit

execs/findsolnx -orb -fname solution_candidate.h5 -valRa 7.2e4 -valBe 3.8e5 -valCc 0 -valPr 1 -valLx 1.0471975511965976 -Nx 256 -Ny 64 -T T_candidate.asc -es 1e-9 -d 1e-8 -edt 1e-10 -Dtmax 1e4 -Dtmin 1e-4 -xrel -ax ax_candidate.asc -dmax 1e4 -Nn 100

### 2. continuationx
Commandline arguments for continuing a solution:

execs/continuationx -eqb -fname solution.h5 -Sn 0 -Nx 64 -Ny 64 -cont Ra -valRa 5000 -valBe 0 -valCc 0 -valPr 1 -Bwd -T 0.01 -es 1e-8 -dmu 1e1 -dsmax 1e1 -dsmin 1e-12 -errmin 5e-7 -errmax 2e1 -Nn 60 -Ng 1000 -d 1e3 -dmax 1e5 -dmin 1e-12 -dg 17 

### 3. findeigenvalsx
Commandline arguments for finding the eigenvalues of a solution:

execs/findeigenvalsx -fname solution.h5 -N 100 -valRa 3100 -valBe 0 -valCc 0 -valPr 1 -T 0.01

# Author: 

Ayse Yesil 

# Contributors:

 Mirko Farano -- contributed to the implementation of eigenvalue calculator of the library

# License: 

deDSI is released under the [GNU GPL version 2](./LICENSE)


