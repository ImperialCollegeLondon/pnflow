


class Solver
{

public:
 
        virtual ~Solver() {};

    virtual double flowrate(double inPrs, double outPrs, const Fluid& fluid, double& flowErr, double& elap,
        double satWat = 1.0, bool writeVelocity = false, bool writeMat = false, bool resistivitySolve = false) = 0;


};



//#include "netlib_amg_1995_prototypes.h"

//#include "netlib_amg_solver.h"


#include "hypre.h"
