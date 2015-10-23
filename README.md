c12sim
======

A Geant4 simulation for clas12.

Overview
--------

Unlike gemc, this project is meant to focus only on CLAS12 related detectors.  
It uses pure geant4 everywhere possible. The only new external dependency is a 
project called [clasdigi](https://github.com/whit2333/clasdigi) which provides 
a few handy libraries.  These libraries can be easily extended and used in 
further analysis, e.g., with ROOT.

One library provides the ROOT io classes to easily extract detector hit/event 
information from the simulation. Each detector has their own implementations, 
unlike gemc, which uses one catch-all hit to gather copious and unnecessary hit 
information. Instead, with the event and hit information well defined before 
hand, the goal is clear and the information to be extracted is cleanly defined, 
thus avoiding unnecessary computation and IO.


Build and Install
-----------------

Standard cmake build assuming you have geant4 installed.
Note that you need a compiler that supports c++14, otherwise you will have 
problems. (UPDATE YOUR COMPILER!)

### First install ClasDigi

    git clone https://github.com/whit2333/clasdigi.git
    mkdir clasdigi_build && cd clasdigi_build
    cmake ../clasdigi/. -DCMAKE_INSTALL_PREFIX=$HOME -DPROJECT_USE_ROOT6=1
    make install

    git clone git@gitlab.com:whit2333/c12sim.git
    mkdir c12sim_build
    cd c12sim_build
    make ../c12sim/. -DCMAKE_INSTALL_PREFIX=../c12run
    make install

Run the simulation

    cd ../c12run
    ./bin/c12sim share/C12SIM/examples/vis.mac 

Using git 
--------- 

First you will want to do the following if you have never used git:

    git config --global user.name  "Your Name"
    git config --global user.email "youremail@wherever"

See the [git configuration 
documentation](https://git-scm.com/book/en/v2/Customizing-Git-Git-Configuration) 
for more details.

Running options
---------------

The default  will run 1000 events  and save histograms to run file "0".
The run number can be specified with the "--run" flag

    ./bin/c12sim  --run=100

This creates the file EBL_sim_output_100.root.




