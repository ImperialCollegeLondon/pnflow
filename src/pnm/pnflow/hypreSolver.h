

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
class Elem;

/// A little storage class for four elemnts
template<typename TOne, typename TTwo, typename TThree, typename TFour>
class FourSome
{
public:

	FourSome() {}
	FourSome(TOne one, TTwo two, TThree three, TFour four) : first_(one), second_(two), third_(three), fourth_(four) {}

	const TOne& first() const {return first_;}
	const TTwo& second() const {return second_;}
	const TThree& third() const {return third_;}
	const TFour& fourth() const {return fourth_;}
	void first(TOne entry) {first_ = entry;}
	void second(TTwo entry) {second_ = entry;}
	void third(TThree entry) {third_ = entry;}
	void fourth(TFour entry) {fourth_ = entry;}

private:

	TOne                first_;
	TTwo                second_;
	TThree              third_;
	TFour               fourth_;

};

class hypreSolver
{

public:
 
	hypreSolver(const vector<Elem*>& network, const vector<Elem*>& inlet, const vector<Elem*>& outlet,
	  int, int, int debugMode, string matFileName, bool matlabFormat);


	double flowrate(double inPrs, double outPrs, const Fluid& fluid, double& flowErr, double& elap, double a=1., bool b=false, bool c=false);
	//static void inithypreSolve(double eps, int scaleFact, int slvrOutput, bool verboseSlvr, bool useGrav);


private:

	double getFlowRate(HYPRE_IJVector xxx, const Fluid*, double&, int&) const;
	void fillMatrixHypre(HYPRE_IJMatrix AAA, HYPRE_IJVector bbb, HYPRE_IJVector xxx, double, double, const Fluid&, bool writeVelocity);

	double getHypreSP(HYPRE_IJVector xxx, const Fluid* fluid, double& maxError, int& nErrors) const;
	void fillMatrixHypreSP(HYPRE_IJMatrix AAA, HYPRE_IJVector bbb, HYPRE_IJVector xxx, double inletPrs, double outletPrs, const Fluid& fluid, bool writeVelocity);

	static bool                                        MPIINITIALISED;
	static bool                                        INITIALISED;
	static bool                                        USE_GRAVITY;
	static double                                      TOLERANCE;
	static const double                                SCALE_FACTOR;


	const int                                          debugMode_;

	const vector<Elem*>&                            elemans_;
	vector< FourSome<int, double, double, double> >    networkOutlets_;
	vector< FourSome<int, double, double, double> >    networkInlets_;
	const vector<Elem*>&                            inPors_;
	const vector<Elem*>&                            outPores_;
	const int                                          nBSs_;
	const int                                          nBpPors_;




	vector< pair<const Elem*, double> >          throatConductances_;
	vector< int>                                    rowedPores_;
	vector< int>                                    poreiRows_;

	const string                                       matrixFileName_;

};


#endif
