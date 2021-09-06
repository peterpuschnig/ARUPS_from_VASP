# ARUPS_from_VASP
This repository contains a fortran code which simulates angle-resolved ultra-violet photoemission spectroscopy (ARUPS) intensity maps based on the one-step-model of photoemission. As initial states it uses Kohn-Sham orbitals from a VASP calculation (WAVECAR-file) and as final state, either a plane wave or a damped plane wave can be used. The program simulates ARUPS intensity maps either as "band-maps" I=I(k,Eb) for a given emission plane, as constant binding state "k-maps" I=I(kx,ky), or as constant intial state scans "CIS-scan" I=I(k,E_photon) for varying photon energies. More information about the method an application results can be found in the following publications:

[1] D. Lüftner, P. Hurdax, G. Koller, M. G. Ramsey, S. Weiß, X. Yang, S. Soubatch, F.S. Tautz, V. Feyer, A. Gottwald and P. Puschnig,
Understanding the photoemission distribution of strongly interacting two-dimensional overlayers,
Phys. Rev. B 96, 125402 (2017).
https://doi.org/10.1103/PhysRevB.96.125402

[2] P. Puschnig, S. Berkebile, A.J. Fleming, G. Koller, K. Emtsev, T. Seyller, J.D. Riley, C. Ambrosch-Draxl, F.P. Netzer, and M.G. Ramsey,
Reconstruction of Molecular Orbital Densities from Photoemission Data
Science 326, 702-706 (2009).
https://doi.org/10.1126/science.1176105

Please cite these publications when using results from this code in your publications.

## Authors
- Peter Puschnig (peter.puschnig@uni-graz.at)
- Daniel Lüftner (daniel.lueftner@hotmail.com)
- Christoph Dösinger (christoph.doesinger@unileoben.ac.at)

## Quick-Start

Compilation of program:

    cd src
    gfortran arups_spinpol.f90
    mv a.out ../bin/arups_spinpol.x

Usage:

    cd ../example_bandmap/
    ../bin/arups_spinpol.x

The program runs about 2 seconds for this example data (single-layer graphene) and produces the text outputfile `ARUPS.out` containing the ARUPS intensity data and the log-file `ARUPS.log`. 
Compare your output with the sample outputfile provided in the `example_bandmap/output/` folder which also contains the `ARUPS.png` displaying the data containd in `ARUPS.out`.

## Example-data

### VASP data

The `VASP` folder contains all necessary input files for running VASP for the example strucure of a single-layer graphene. 
In particular, it already contains the `WAVECAR` file containing the Kohn-Sham energies and orbitals, which is the only input rquired to run the ARUPS program.

### Band-map

In order to simulate an ARUPS band map go to the respective folder,

    cd ../example_bandmap/
    
If necessary, edit the text file `ARUPS.in` which provides the input data for running the ARUPS program and 
make sure that the `WAVECAR` is located in the same folder or create a soft link to it, e.g.: `ln -s ../VASP/WAVECAR .`. Execute the program  

    ../bin/arups_spinpol.x    
    
The resulting text output file `ARUPS.out` contains the computed ARUPS intensity. Lines 1-11 repeat the most important input data.
Line 12 gives the size of the two-dimensional intensity array, e.g. `251  176` means 251 points along the binding energy axis and 176 points along the k-axis.
Line 13 and 14 provide the minimum , maximum and step-size values for the energy and k axes, respectively. The remaining lines contain the 251x76 intensity values for the 2D-array.
You may use the files in `example_bandmap/output` as a reference.

### k-map

In order to simulate an ARUPS k-map for a constant binding energy go to the respective folder,

    cd ../example_kmap/
    
If necessary, edit the text file `ARUPS.in` which provides the input data for running the ARUPS program and 
make sure that the `WAVECAR` is located in the same folder or create a soft link to it, e.g.: `ln -s ../VASP/WAVECAR .`. Execute the program  

    ../bin/arups_spinpol.x    
    
The resulting text output file `ARUPS.out` contains the computed ARUPS intensity. Lines 1-11 repeat the most important input data.
Line 12 gives the size of the two-dimensional intensity array, e.g. `126  126` means 126 points along kx-axis and 126 points along the ky-axis.
Line 13 and 14 provide the minimum , maximum and step-size values for both k axes, respectively. The remaining lines contain the 126x126 intensity values for the 2D-array.
You may use the files in `example_kmap/output` as a reference.


### CIS scan

In order to simulate an ARUPS constant initial state (CIS) scan for a constant binding energy and chosen emission plane go to the respective folder,

    cd ../example_CIS/
    
If necessary, edit the text file `ARUPS.in` which provides the input data for running the ARUPS program and 
make sure that the `WAVECAR` is located in the same folder or create a soft link to it, e.g.: `ln -s ../VASP/WAVECAR .`. Execute the program  

    ../bin/arups_spinpol.x    
    
The resulting text output file `ARUPS.out` contains the computed ARUPS intensity. Lines 1-11 repeat the most important input data.
Line 12 gives the size of the two-dimensional intensity array, e.g. `131  151` means 131 points for the varying photon energy and 151 points along the k-axis.
Line 13 and 14 provide the minimum , maximum and step-size values for the photon energy and k-axis, respectively. The remaining lines contain the 131x151 intensity values for the 2D-array.
You may use the files in `example_CIS/output` as a reference.

## Generating VASP files for your own systems

In addition to a standard scf-run from VASP, the ARUPS program in most cases requires an additional non-scf run to produce the necessary `WAVECAR` file.
To this end, use the following additional settings in the `INCAR` file (refer to https://www.vasp.at/wiki/index.php/The_VASP_Manual for more information)

    ICHARG = 11      # this switches the non-self-consiostent calculation on
    ISYM = 0         # switches off symmetry since teh ARUPS program cannot deal with reduced k-point lists
    LDIPOL = .TRUE.  # turn on the dipole correction in your slab calculation
    IDIPOL = 3       # the surface normal of your slab must point in z-direction
    
In most cases, it is also necessary to increase the number of k-points in the non-scf run to achieve a desired resoluition of the resulting ARUPS intensity maps. 
On the exmaple of the single-layer layer graphene with lattice parameters of a=2.47 Angström, b= 2.47 Angström and c=20 Angström, the example `KPOINT` file

    Gamma-centered k-grid                           
    0
    Gamma
    30 30 1 
    0 0 0

leads to a spacing of points in reciprocal space of Delta-kx = 2*pi/(30*a) = 0.085 1/Angström, Delta-ky = 2*pi/(30*b) = 0.085 1/Angström and Delta-kz = 2*pi/(1*c) = 0.31 1/Angström. 
Taking into account the max(Delta-kx, Delta-ky, Delta-kz) = 0.31 1/Angström, the k-space broadening settings in line 9 of `ARUPS.in` should be set to a similar value. 
Thus, if the desired k-resolution is smaller, e.g., Delta-k = 0.05 1/Angström, the non-scf VASP calculation should be run with a k-grid of N_x = [2*pi/a*Delta-k] = 51,  
N_y = [2*pi/b*Delta-k] = 51 and N_z = [2*pi/c*Delta-k] = 6 leading to the following `KPOINTS` file:

    Gamma-centered k-grid                           
    0
    Gamma
    51 51 6 
    0 0 0

For more details, please refer to Ref. [1] mentioned above: https://doi.org/10.1103/PhysRevB.96.125402 

## Compilation of the ARUPS program

### Serial version

The serial version of the ARUPS program is located in `src` and is called `rups_spinpol.f90`. It can be compiled, for instance, with `gfortran` in a simple manner:

    gfortran arups_spinpol.f90 -o arups_spinpol.x
    
### MPI version

The MPI-parallel version of the ARUPS program is located in `src` and is called `arups_spinpol_mpi_v2.f90`. It can be compiled, for instance, with `mpifort` in a simple manner:

    mpifort arups_spinpol_mpi_v2.f90 -o arups_spinpol_mpi_v2.f90.x
    
and executed, e.g. with using 4 parallel processes, as follows:

    mpirun -np 4 ../bin/arups_spinpol_mpi_v2.x
    
### Spinor version

The program `arups_spinor.f90` in `src` contains an experimental, serial version of the program which can be applied to spionor wave functions which are, for instance, 
required in the case spin-orbital coupling. Please refer to the Master's thesis of Christoph Dösinger for further information: https://homepage.uni-graz.at/de/peter.puschnig/theses/ 

------
September 6th, 2021
