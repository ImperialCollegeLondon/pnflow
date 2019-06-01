

#ifndef hypreSolve_H
#define hypreSolve_H


///. to define 'integer'
#include <vector>
#include <cmath>
#include "_hypre_utilities.h"
#include "HYPRE.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"

using namespace std;
class Element;

class hypreSolver : public Solver
{

public:
 
	hypreSolver(const vector<Element*>& network, const vector<Element*>& inlet, const vector<Element*>& outlet,
	  int outletIdx, int debugMode, string matFileName, bool matlabFormat);


	double flowrate(double inPrs, double outPrs, const Fluid& fluid, double& flowErr, double& elap, double a=1.0, bool b=false, bool c=false, bool resistivitySolve=false);
	//static void inithypreSolve(double eps, int scaleFact, int slvrOutput, bool verboseSlvr, bool useGrav);


private:

	double getFlowRate(HYPRE_IJVector xxx, const Fluid*, double&, int&, bool) const;
	void fillMatrixHypre(HYPRE_IJMatrix AAA, HYPRE_IJVector bbb, HYPRE_IJVector xxx, double, double, const Fluid *, bool writeVelocity, bool resistivitySolve);

	double getHypreSP(HYPRE_IJVector xxx, const Fluid* fluid, double& maxError, int& nErrors) const;
	void fillMatrixHypreSP(HYPRE_IJMatrix AAA, HYPRE_IJVector bbb, HYPRE_IJVector xxx, double inletPrs, double outletPrs, const Fluid *fluid, bool writeVelocity, bool resistivitySolve);

    static bool                                             MPIINITIALISED;
    static bool                                             INITIALISED;
    static bool                                             USE_GRAVITY;
    static double                                           TOLERANCE;
    static const double                                     SCALE_FACTOR;


	const int                                     			m_debugMode;

	const vector<Element*>&                            m_netelems;
	vector< FourSome<int, double, double, double> >         m_networkOutlets;
	vector< FourSome<int, double, double, double> >         m_networkInlets;
	const string                                       m_matrixFileName;
	const vector<Element*>&                              m_inPores;
	const vector<Element*>&                              m_outPores;
	const int                                               m_numPoresp1;




	vector< pair<const Element*, double> >		m_throatConductances;
	vector< int>                						m_rowedPores;
	vector< int>                						m_poreiRows;


};


#endif
