


 ----------------------------------------------------------------
* Please see the apps/doc directory for description of the codes *
 ----------------------------------------------------------------




###### Compiling ######

To compile, open a terminal in the upper most directory and run:

 make 

once everything compiled successfully, to clean the temporary files, type:

 make clean

The above command can be run inside most of the subfolders, wherever a 
makefile or Makefile is present.  The libraries, those with a makefile,
should be compiled before the apps that contain "Makefile"s.

Compilation requires gnu (Linux) make, cmake, a c++ compiler with -std=c++11
support and an MPI. The compilation is tested using g++ (version 5+) (default)
and using intel-2018 compilers.


###### Contact and References ######

For contacts and references please see: 
http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling
or contact Ali Q. Raeini, email: a.qaseminejad-raeini09@imperial.ac.uk

More details are given in the apps/doc directory.

