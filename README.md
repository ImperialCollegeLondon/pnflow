##  pnflow - classical pore-network flow simulation

This repository included the classical network flow simulation code called pnflow. 
[pnextract](https://github.com/aliraeini/pnextract) network extraction 
code is also included for convenience.

These pnflow code is essentially a cleaned up version of the [poreflow code] by
[Valvatne and Blunt (2004)].  The poreflow code was restructured with some 
post-processing features added in preparation for the generalized network 
model [Raeini, Bijeljic and Blunt 2018], sponsored by [TOTAL]. 
A recent validation of the pnflow and [pnextract] codes is published by [Bultreys et al. 2018].



### Instructions
Instructions for extracting a network from a micro-CT image are given in
the doc folder; see also the [pnextract] github page.

To run a network flow simulation, copy the sample input file src/doc/input_pnflow.dat
and the extracted networks (*_link1.dat, _link1.dat, _node1.dat, 
_node2.dat) into a new folder.  Then edit the input_pnflow.dat file and set the NETWORK
keyword and other flow parameters. Finally run, in a Windows Command Prompt: 
  
   bin\pnflow.exe  input_pnflow.dat

* You may need to modify the command above and provide the full path for the bin\pnflow.exe,
  if it exists in a different directory than your command prompt working directory.

* To open a command-prompt in Windows, hold the *Shift* key and *right-click*
  into the folder where the input_pnflow.dat is copied and click on the *Open Command Window Here* menu.

###  Build instructions:
The code is already compiled to bin/pnextract.exe and bin/pnflow.exe, Win64 
executables (extract the bin.7z to see these files) using mingw compilers.

The compilation can be done in Linux by running './AllMake' bash script.

The './AllMakeMinGW' bash script compiles the code for Windows machines.
Run './AllClean' beforhand, to avoid mixing the intermidiate Linux and Windows object files.

###  Dependencies
The pnflow code depends on the [Hypre] library for solving linear equations. 
This library, along with [pnextract] dependencies, is included in the thirdparty/ folder. 

###  Licence

The code is release as a free, using a zlib-style licence, see the file 
src/pnflow/netsim.h for a copy of the licence.

For contact and further information see [Imperial College - pore-scale consortium] website,
or send me an email:   a.qaseminejad-raeini09@imperial.ac.uk




[Imperial College - pore-scale consortium]: http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling
[poreflow code]: http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling/software/two-phase-network-modelling-code
[Valvatne and Blunt (2004)]:  https://doi.org/10.1029/2003WR002627
[Bultreys et al. 2018]: https://doi.org/10.1103/PhysRevE.97.053104
[Raeini, Bijeljic and Blunt 2018]: https://doi.org/10.1103/PhysRevE.97.023308
[Hypre]: https://github.com/LLNL/hypre
[TOTAL]: https://www.total.com
[pnextract]:  (https://github.com/aliraeini/pnextract)
