#ifndef INPUTDATA_H
#define INPUTDATA_H

#include "InputFile.h"


typedef struct
{
	int         index;
	double      x;
	double      y;
	double      z;
	int         connNum;
	double      volume;
	double      radius;
	double      shapeFact;
	double      clayVol;
} PoreStruct;

typedef struct
{
	int     index;
	int     poreOne;
	int     poreTwo;
	double  radius;
	double  shapeFact;
	double  lenPoreOne;
	double  lenPoreTwo;
	double  lenThroat;
	double  lenTot;
	double  volume;
	double  clayVol;
} ThroatStruct;

/// A little storage class for three elements
template<typename TOne, typename TTwo, typename TThree>
class ThreeSome
{
public:

	ThreeSome() {}
	ThreeSome(TOne one, TTwo two, TThree three) : first_(one), second_(two), third_(three) {}

	const TOne& first() const {return first_;}
	const TTwo& second() const {return second_;}
	const TThree& third() const {return third_;}
	void first(TOne entry) {first_ = entry;}
	void second(TTwo entry) {second_ = entry;}
	void third(TThree entry) {third_ = entry;}
	bool operator==(ThreeSome b)
	{return first_ == b.first() && second_ == b.second() && third_ == b.third() ;}
	bool operator!=(ThreeSome b)
	{return first_ != b.first() || second_ != b.second() || third_ != b.third() ;}

private:

	TOne                first_;
	TTwo                second_;
	TThree              third_;

};

class InputData : public InputFile
{
public:

	InputData(const InputFile&);

	ststr title() const;
	int randSeed() const;
	void prsBdrs(bool& usePrsBdr, bool& reportPrsBdr, int& numPlanes) const;
	void solverTune(double& eps, double& scaleFact, int& slvrOutput, bool& verbose, double& condCutOff) const;
	void poreFillWgt(stvec< double >& weights)            			const;
	void poreFillAlg(ststr& algorithm)                     	 		const;
	void relPermCompression(bool& useComp, double& krThres, 			//
			double& deltaSw, bool& wettPhase, bool& nonWettPhase)  	const;
	void satConvergence(int& minNumFillings, double& initStepSize, 	//
			   double& cutBack, double& maxIncrFact, bool& stable) 	const;
	void fillingList(bool& drainage, bool& imbibition, bool& location) const;
	void output(bool& propOut, bool& swOut) 							const;
	void resFormat(bool& matlabFormat, bool& excelFormat, bool& mcpFormat) 				const;
	void fluid(double& intfacTen, double& watVisc, double& oilVisc, //
			   double& watRes, double& oilRes, double& watDens, double& oilDens) const;
	void aCloseShave(double& shaveSetting)                        		const;
	void calcBox(double& inletBdr, double& outletBdr)					const;
	void prsDiff(double& inletPrs, double& outletPrs, bool& useGravInKr) const;

	void initConAng(int& wettClass, double& min, double& max, double& delta, double& eta, ststr& mod, double& sep) const;
	void writeNetwork(bool& writeNet, bool& binary, ststr& netName) const;
	void solverDebug(bool& watMat, bool& oilMat, bool& resMat, bool& watVel, bool& oilVel, bool& resVel, bool& matlab, bool& initOnly) const;
	void getModifyRadDistOptions(int& throatModel, int& poreModel, ststr& throatOpt, ststr& poreOpt,
		bool& maintainLtoR, bool& toFile, int& numPts) const;
	void getModifyPoro(double& netPoroTrgt, double& clayPoroTrgt) const;
	void modifyConnNum(double& targetConnNum, ststr& model) const;
	void getModifyGDist(int& throatModel, int& poreModel, ststr& throatOpt, ststr& poreOpt,
		bool& toFile, int& numPts) const;
	void getModifyModelSize(double& scaleFactor) const;
	void clearNetworkData(); //finishedLoadingNetwork
	void poreLocation(int idx, double& xPos) const;

	///. non constant public member functions:

	void network(int& numPores, int& numThroats, double& xDim, double& yDim, double& zDim);//. TODO document
	void poreData(int idx, double& x, double& y, double& z, int& n, stvec< int >& throats,
		stvec<int>& pores, double& c, double& vcl, double& e, double& g)	;        		//something wierd hapenning, but probably this is a const function;
	void throatData(int idx, int& p1, int& p2, double& v, double& vcl, double& r, double& g, double& lp1,
		double& lp2, double& lt, double& lTot); //* similar to poreData


private:


	typedef std::multimap< int, ThreeSome<int, int, double> >::value_type MultiConnValType;
	typedef std::multimap< int, ThreeSome< int, int, double > >::iterator MapItr;


	void getRadDist(std::istream& data, int& model, ststr& options) const;
	void getOutletData(MapItr itr, int numNetsInFront, int& throatIdx, int& thatIdx) const;
	void getInletData(MapItr itr, int numNetsInFront, int& throatIdx, int& thatIdx) const;
	int findClosestPore(const stvec< ThreeSome< PoreStruct*, int*, int* > > pores, double xPos,
		double yPos, double zPos, double& totalLen) const;



	inline bool outletThroat(int throatIdx) const;
	inline bool inletThroat(int throatIdx) const;
	void findClosestPoreForPBC(const stvec< ThreeSome< double, double, double > >& pos,
		stvec< int >& idx, stvec< double >& p2BLen) const;


	///.   non constant private member functions:

	inline void appendPoreData(int thisPoreIdx, int throatIdx, int thatPoreIdx);
	void loadPoreData();
	void loadThroatData();
	void findBoundaryPores();


	static const int                                        DUMMY_INDEX;

	std::ifstream                                                poreConn_;
	std::ifstream                                                poreProp_;
	std::ifstream                                                throatConn_;
	std::ifstream                                                throatProp_;

	int                                                     connectionsRemoved_;
	int                                                     origNumPores_;
	int                                                     origNumThroats_;
	int                                                     numNetInSeries_;
	bool                                                    binaryFiles_;
	bool                                                    useAvrXOverThroatLen_;
	bool                                                    useAvrPbcThroatLen_;
	bool                                                    addPeriodicBC_;
	stvec< ThreeSome< PoreStruct*, int*, int* > >          poreData_;
	stvec< ThroatStruct* >                                 throatData_;
	stvec< ThreeSome< PoreStruct*, int*, int* > >          inletPores_;
	stvec< ThreeSome< PoreStruct*, int*, int* > >          outletPores_;
	std::multimap< int, ThreeSome< int, int, double > >    outletConnections_;
	std::multimap< int, ThreeSome< int, int, double > >    inletConnections_;
	std::set< std::pair< int, int > >                      xyPbcConn_;
	std::set< std::pair< int, int > >                      xzPbcConn_;
	stvec< ThreeSome< int, int, double > >                 pbcData_;
	stvec< int >                                           throatHash_;
	stvec< int >                                           reverseThroatHash_;
	double                                                  origXDim_;
	double                                                  origYDim_;
	double                                                  origZDim_;
	double                                                  averageThroatLength_;
	double                                                  averagePoreHalfLength_;
	double                                                  networkSeparation_;
};

inline bool InputData::outletThroat(int throatIdx) const
{
	ThroatStruct *throat = throatData_[throatIdx-1];
	return (throat->poreOne == 0 || throat->poreTwo == 0);
}

inline bool InputData::inletThroat(int throatIdx) const
{
	ThroatStruct *throat = throatData_[throatIdx-1];
	return (throat->poreOne == -1 || throat->poreTwo == -1);
}





inline void InputData::appendPoreData(int thisPoreIdx, int throatIdx, int thatPoreIdx)
{
	int oldConnNum = poreData_[thisPoreIdx-1].first()->connNum;
	int *newThroatConnList = new int[oldConnNum+1];
	int *newPoreConnList = new int[oldConnNum+1];
	for(int i = 0; i < oldConnNum; ++i)
	{
		newPoreConnList[i] = poreData_[thisPoreIdx-1].second()[i];
		newThroatConnList[i] = poreData_[thisPoreIdx-1].third()[i];
	}
	newPoreConnList[oldConnNum] = thatPoreIdx;
	newThroatConnList[oldConnNum] = throatIdx;
	poreData_[thisPoreIdx-1].first()->connNum++;
	delete[] poreData_[thisPoreIdx-1].second();
	delete[] poreData_[thisPoreIdx-1].third();
	poreData_[thisPoreIdx-1].second(newPoreConnList);
	poreData_[thisPoreIdx-1].third(newThroatConnList);
}


#endif

