c12sim
======

A Geant4 simulation for clas12.

Build and Install
-----------------

Standard cmake build assuming you have geant4 installed.

    git clone https://github.com/whit2333/c12sim 
    mkdir c12sim_build
    cd c12sim_buildc
    make ../EBL_sim/. -DCMAKE_INSTALL_PREFIX=../ebl_run
    make install

Run the simulation

    cd ../ebl_run
    ./bin/ebl1  

Running options
---------------

The default  will run 1000 events  and save histograms to run file "0".
The run number can be specified with the "--run" flag

    ./bin/ebl1  --run=100

This creates the file EBL_sim_output_100.root.




