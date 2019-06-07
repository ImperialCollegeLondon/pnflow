#ifdef WIN32
#pragma warning(disable:4786)
#endif


#include "netsim.h"
//#include "mpi.h"

#ifndef __DATE__
 #define __DATE__  "2016"
#endif 

#include "globals.cpp"
int main(int argc, char *argv[])
{

 try
 {
	//int myid(0), num_procs(1);
	//MPI_Init(&argc, &argv);
	//MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	//MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    string inputFileName;
    cout<<"\n Network Model Code version 2 alpha, built: "<<__DATE__ << endl;

    if (argc > 1)
        inputFileName = argv[1];
    else
    {
        cout<< "Please input data file : ";
        cin >> inputFileName;
    }

    InputData input(inputFileName);


    Netsim netsim(input);


    bool propOut(false), swOut(false);
    input.output(propOut, swOut);

    bool calcKr, calcI, EntreL, EntreR, ExitL, ExitR;
    double sw(1.0), pc(0.0); ///. Sw initial condition
    double requestedFinalSw, requestedFinalPc, deltaSw, deltaPc, deltaPcIncFactor;

    while(input.satTarget(requestedFinalSw, requestedFinalPc, deltaSw, deltaPc, deltaPcIncFactor, calcKr, calcI, EntreL, EntreR, ExitL, ExitR))
    {
        if(requestedFinalSw > 1.0)
        {
            cerr << "\nError: Saturations to be given as fractions \n" << endl;       exit(-1);
        }

        if(requestedFinalSw < sw)
        {
            netsim.Drainage(input, sw, pc, requestedFinalSw, requestedFinalPc, deltaSw, deltaPc, deltaPcIncFactor,
               calcKr, calcI, EntreL, EntreR, ExitL, ExitR, swOut);
        }
        else
        {
            netsim.Imbibition(input, sw, pc, requestedFinalSw, requestedFinalPc, deltaSw, deltaPc, deltaPcIncFactor,
                calcKr, calcI, EntreL, EntreR, ExitL, ExitR, swOut);
        }

    }
    
   //MPI_Finalize();

 }
 catch (std::exception &exc)
 {
	std::cerr << "\n\n Exception on processing: \n" << exc.what() << "Aborting! \n"  << endl;
	return 1;
 }
 catch (...)
 {
	std::cerr << "\n\n Unknown exception! \n Aborting! \n" << endl;
	return 1;
 }

 return 0;
}




// f2c, used for in old bu/AMG solver, on linux needs a reference to the entry point MAIN__ 
// otherwise it will not link. This function is never called
//extern "C" { void MAIN__() { cout<<"This function is never called"<<endl;  exit(-1); } }


