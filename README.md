##  pnflow - classical pore-network (extraction and) flow simulation

This repository included the classical network flow simulation code called pnflow. 
[pnextract network extraction 
code](https://github.com/aliraeini/pnextract) is also included [here](src/pnm/pnextract) for convenience.

These pnflow code is essentially a cleaned up version of the [poreflow code] by
[Valvatne and Blunt (2004)].  The poreflow code was restructured with some 
post-processing features added in preparation for the generalized network 
model [Raeini, Bijeljic and Blunt 2018], sponsored by [TOTAL]. 
A recent validation of the pnflow and [pnextract] codes is published by [Bultreys et al. 2018].



### Notice:
* USE ``;`` (or double empty lines) to mark the end of input_pnflow.dat keyword data. 
  ``#`` can not be used for this purpose any more. 

* Previous code has been moved to branch [bu20190607] (discontinued) 

### Instructions for Windows:

Instructions for extracting a network from a micro-CT image are given in
the doc folder; see also the [pnextract] README.

To run a network flow simulation, copy the sample input file src/doc/input_pnflow.dat
and the extracted networks (with suffixes  _link1.dat, _link1.dat, _node1.dat and 
_node2.dat) into a new folder.  Then edit the input_pnflow.dat file and set the NETWORK
keyword and other flow parameters. Finally run, in a Windows Command Prompt: 
  
    PATH\TO\bin\pnflow.exe  input_pnflow.dat

* You may need to modify the command above and, instead of `` PATH\TO\bin\ ``, provide the full path to the pnflow.exe, if it exists in a different directory than your command prompt working directory.

* To open a command-prompt in Windows, hold the *Shift* key and *right-click*
  into the folder where the input_pnflow.dat is copied and click on the *Open Command Window Here* menu.

###  Compiling:
The code is already compiled to bin/pnextract.exe and bin/pnflow.exe, Win64 
executables (extract the bin.7z to see these files) using MinGW compilers.

The compilation can be done in Linux by running in the top-level directory (where this file is):    

    make -j

"-j" flag is to compile in parallel and can be omitted. 

Running ``make mgw`` instead cross-compiles the code for Windows machines, if MinGW compilers are ailable in your system and set in [src/AllMakeMinGW](src/AllMakeMinGW) script.
Run ``make clean`` before switching between Windows and Linux compilations, to avoid mixing the intermidiate Linux and Windows object files.

### Installation

In Linux, you can source the src/bashrc file to set te paths to the compiled binaries:     

     source PATH/TO/src/bashrc

In Windows, you just need to know how to run standalone exe files from Command Prompt, as discussed above.

###  Dependencies:
The pnflow code depends on the [Hypre] library which along with other [pnextract] dependencies are included in 
the [thirdparty](thirdparty) folder. 


Note, may need to install jpg and lzma libraries, prerequisites of libtif and oxelImage/libvoxel libraries (to check...), for the make command above to succeed. In Ubuntu this can be installed by running the following commands:      

    sudo apt install libjpeg-dev liblzma-dev
      
###  Licence:

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
[pnextract]:  src/pnm/pnextract
[bu20190607]:  https://github.com/aliraeini/pnflow/tree/bu20190607
