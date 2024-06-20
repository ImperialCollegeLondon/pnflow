##  pnflow - classical pore-network (extraction and) flow simulation

This repository includes the classical network flow simulation code called pnflow. 
[pnextract network extraction code](https://github.com/aliraeini/pnextract) is also 
included [here in pnextract folder](src/pnm/pnextract) for convenience. These two 
codes are kept compatible and refered to as pore-network model (PNM) together.

The pnflow code is essentially a cleaned up version of the [poreflow code] by
[Valvatne and Blunt (2004)].  The poreflow code was restructured with some 
post-processing features added in preparation for the generalized network 
model [Raeini, Bijeljic and Blunt 2018], sponsored by [TOTAL]. 
A recent validation of the pnflow and [pnextract] codes is published by [Bultreys et al. 2018], [Foroughi et al. 2020], and [Foroughi et al. 2021].


### Release notes

- In 2020-2021, the code has been through further restructuring and clean up as part of a sponsorship from [Wintershall Dea].
The common components of PNM codes are kept in sync with privately developed codes and the closed-source generalized network model to reduce development effort and to allow future collaborations. 
However, this code is still well behind my local branch. If you like to work on this code development, please drop me an email for updates. 
The physical assumptions used in the model are not affected, so if you are a user of the code, you can continue using its older versions.  The previous version of the code can be downloaded from [pnm2019 branch](https://github.com/aliraeini/pnflow/tree/pnm2019).

- pnflow input keywords are changed significantly (for compatibility with the generalized network model) see the new [input_pnflow.dat](https://github.com/aliraeini/pnflow/blob/master/doc/input_pnflow.dat) file in doc folder for more information.

* The file [pnflow/ChangeLog](pnflow/ChangeLog) is used to document the main changes or extensions of the pnflow code.

----------------------------------------

### Instructions for MS Windows

Instructions for extracting a network from a micro-CT image are given in
the doc folder; see also the [pnextract] README.md file.

To run a network flow simulation, copy the sample input file src/doc/input_pnflow.dat
and the extracted networks (with suffixes  _link1.dat, _link1.dat, _node1.dat and 
_node2.dat) into a new folder.  Then edit the input_pnflow.dat file and set the NETWORK
keyword and other flow parameters. Finally run, in a Windows Command Prompt: 
  
    PATH\TO\bin\pnflow.exe  input_pnflow.dat

* Replace `PATH\TO\bin\`` with the path to the extracted pnflow.exe.

* To open a command-prompt in Windows, hold the `Shift` key and `right-click`
  into the folder where the input_pnflow.dat is copied and click the `Open Command Window Here` menu.


The instructions are similar in Linux. Additionally, you can source the src/script/bashrc file to set the paths to the compiled binaries:     

     source PATH/TO/src/script/bashrc

###  Build instructions

Download and extract [bin.7z](../../bin.7z) for pnextract.exe and pnflow.exe, 
which are Win64 executables compiled using MinGW compilers.  

See [../script/README.md](../script/README.md) for build instructions.
The pnflow code depends on the [Hypre] library which along with other [pnextract] dependencies are included in 
the [../../pkgs](pkgs) folder. 


###  Licence

The code is release as a free, using a zlib-style licence, see the file 
src/pnflow/netsim.h for a copy of the licence.

For contact and further information see [Imperial College - Pore-scale Consortium] website,
raise an `issue` on github mentioning @ForoughiSajjad @aliraeini in your message,
or send me an email: s.foroughi@imperial.ac.uk or a.q.raeini@imperial.ac.uk




[Imperial College - Pore-scale Consortium]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling
[poreflow code]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/software/two-phase-network-modelling-code
[Valvatne and Blunt (2004)]:  https://doi.org/10.1029/2003WR002627
[Bultreys et al. 2018]: https://doi.org/10.1103/PhysRevE.97.053104
[Raeini, Bijeljic and Blunt 2018]: https://doi.org/10.1103/PhysRevE.97.023308
[Foroughi et al. 2020]: https://doi.org/10.1103/PhysRevE.102.023302
[Foroughi et al. 2021]: https://doi.org/10.1007/s11242-021-01609-y
[Hypre]: https://github.com/LLNL/hypre
[TOTAL]: https://www.total.com
[pnextract]:  src/pnm/pnextract
[bu20190607]:  https://github.com/aliraeini/pnflow/tree/bu20190607
[Wintershall Dea]: https://wintershalldea.com
