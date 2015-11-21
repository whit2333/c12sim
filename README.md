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

Standard cmake build assuming you have geant4 and root installed.
Also note that you need a compiler that supports c++14, otherwise you will have 
problems. (UPDATE YOUR COMPILER!)

The following sections show how to build clasdigi and c12sim. The installation 
of c12sim is selected to be a separate director but this can be any standard 
location as it only installs the binary and a shared data directory. The shared 
data (installed in share/c12sim) contains various default macros used for 
visualization along with other useful examples.

### First install ClasDigi

    git clone https://github.com/whit2333/clasdigi.git
    mkdir clasdigi_build && cd clasdigi_build
    cmake ../clasdigi/. -DCMAKE_INSTALL_PREFIX=$HOME
    make install

### Then build and install c12sim

    git clone git@gitlab.com:whit2333/c12sim.git
    mkdir c12sim_build
    cd c12sim_build
    make ../c12sim/. -DCMAKE_INSTALL_PREFIX=../c12run
    make install

### Run the simulation

    cd ../c12run
    ./bin/c12sim

Running options
---------------

### Basics

Running the simulation without any arguments will open a Qt gui along with a 
visualization of the detectors. 

### Run number

The run number can be specified with the "--run" flag

    c12sim  --run=100

This sets the run number to 100. The default value is 0.
When the standard geant4 command to generate events is invoked, ie,

    /run/beamOn 5000

This creates the output file <code>data/rootfiles/clas12sim100.root</code>. If 
<code> /run/beamOn</code> is invoked again, the run number is incremented and 
<code>data/rootfiles/clas12sim101.root</code> is created. Note that the event 
numbers of this second file start from where the previous run left off.

### Output

The output file (base)name and directory can be provided as arguments. For 
example,

    c12sim --output=MYOUTPUT --dir=some/directory --run=99

causes the output file of the first run to be 
<code>some/directory/MYOUTPUT99.root</code>




Using git
--------- 

First you will want to do the following if you have never used git:

    git config --global user.name  "Your Name"
    git config --global user.email "youremail@wherever"

See the [git configuration 
documentation](https://git-scm.com/book/en/v2/Customizing-Git-Git-Configuration) 
for more details.


