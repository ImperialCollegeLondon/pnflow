#ifndef INPUTDATA_H
#define INPUTDATA_H

#include "inputFile.h"
#include "threeSome.h"


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

class InputData : public InputFile
{
public:

    InputData(const string&);

    std::string title() const;
    int randSeed();
    void prsBdrs(bool& usePrsBdr, bool& reportPrsBdr, int& numPlanes) const;
    void solverTune(double& eps, double& scaleFact, int& slvrOutput, bool& verbose, double& condCutOff) const;
    void poreFillWgt(vector< double >& weights)            			const;
    void poreFillAlg(string& algorithm)                     	 		const;
    void relPermDef(string& flowRef/*, bool& strictTrpCond*/)   	 		const;
    void relPermCompression(bool& useComp, double& krThres, 			//
            double& deltaSw, bool& wettPhase, bool& nonWettPhase)  	const;
    void satConvergence(int& minNumFillings, double& initStepSize, 	//
               double& cutBack, double& maxIncrFact, bool& stable) 	const;
    //void trapping(bool& drainEnds, double& gMult) const;
    void fillingList(bool& drainage, bool& imbibition, bool& location) const;
    void output(bool& propOut, bool& swOut) 							const;
    void resFormat(bool& matlabFormat, bool& excelFormat, bool& mcpFormat) 				const;
    void fluid(double& intfacTen, double& watVisc, double& oilVisc, //
               double& watRes, double& oilRes, double& watDens, double& oilDens) const;
    double clayEdit()                           		const;
    void aCloseShave(double& shaveSetting)                        		const;
    void apexPrs(bool& doApexAnalysis)									const;
    void gravityConst(double& gravX, double& gravY, double& gravZ)	const;
    void calcBox(double& inletBdr, double& outletBdr)					const;
    void prsDiff(double& inletPrs, double& outletPrs, bool& useGravInKr) const;

    //void fracWetting(bool& applyFracWet, double& whatFrac, bool& volBased, double& min,
        //double& max, double& delta, double& eta, string& fracModel, int& clustDiam, bool& oilInWat) const;
    //void initConAng(double& min, double& max, double& delta, double& eta) const;
    void initConAng(int& wettClass, double& min, double& max, double& delta, double& eta, string& mod, double& sep) const;
    //void equilConAng(int& wettClass, double& min, double& max, double& delta, double& eta, string& mod, double& sep) const;
    void writeNetwork(bool& writeNet, bool& binary, string& netName) const;
    void solverDebug(bool& watMat, bool& oilMat, bool& resMat, bool& watVel, bool& oilVel, bool& resVel,
        string& fileName, bool& matlab, bool& initOnly) const;
    void getModifyRadDistOptions(int& throatModel, int& poreModel, string& throatOpt, string& poreOpt,
        bool& maintainLtoR, bool& toFile, int& numPts) const;
    void getModifyPoro(double& netPoroTrgt, double& clayPoroTrgt) const;
    void modifyConnNum(double& targetConnNum, string& model) const;
    bool sourceNode(int& sourceNode) const;
    void matBal(bool& reportMatBal) const;
    void getModifyGDist(int& throatModel, int& poreModel, string& throatOpt, string& poreOpt,
        bool& toFile, int& numPts) const;
    void getModifyModelSize(double& scaleFactor) const;
    void clearNetworkData(); //finishedLoadingNetwork
    void poreLocation(int idx, double& xPos) const;

    ///. non constant public member functions:

    bool satTarget(double& sw, double& pc, double& dSw, double& dPc,  double& dPcIncFact, 
                   bool& kr, bool& I, bool& entreL, bool& entreR, bool& exitL, bool& exitR) ; //increment m_workingSatEntry

    void network(int& numPores, int& numThroats, double& xDim, double& yDim, double& zDim);//. TODO document
    void poreData(int idx, double& x, double& y, double& z, int& n, vector< int >& throats,
        vector< int >& pores, double& c, double& vcl, double& e, double& g)	;        		//something wierd hapenning, but probably this is a const function;
    void throatData(int idx, int& p1, int& p2, double& v, double& vcl, double& r, double& g, double& lp1,
        double& lp2, double& lt, double& lTot); //* similar to poreData
	const string& netFileBaseName() const {return m_netFileBase;};


private:

    string                                                  m_netFileBase;

    typedef multimap< int, ThreeSome<int, int, double> >::value_type MultiConnValType;
    typedef multimap< int, ThreeSome< int, int, double > >::iterator MapItr;


    void getRadDist(istream& data, int& model, string& options) const;
    void getOutletData(MapItr itr, int numNetsInFront, int& throatIdx, int& thatIdx) const;
    void getInletData(MapItr itr, int numNetsInFront, int& throatIdx, int& thatIdx) const;
    int findClosestPore(const vector< ThreeSome< PoreStruct*, int*, int* > > pores, double xPos,
        double yPos, double zPos, double& totalLen) const;



    inline bool outletThroat(int throatIdx) const;
    inline bool inletThroat(int throatIdx) const;
    void findClosestPoreForPBC(const vector< ThreeSome< double, double, double > >& pos,
        vector< int >& idx, vector< double >& p2BLen) const;


    ///.   non constant private member functions:

    inline void appendPoreData(int thisPoreIdx, int throatIdx, int thatPoreIdx);
    void loadPoreData();
    void loadThroatData();
    void findBoundaryPores();


    static const int                                        DUMMY_INDEX;

    ifstream                                                m_poreConn;
    ifstream                                                m_poreProp;
    ifstream                                                m_throatConn;
    ifstream                                                m_throatProp;

    int                                                     m_connectionsRemoved;
    int                                                     m_origNumPores;
    int                                                     m_origNumThroats;
    int                                                     m_numNetInSeries;
    int                                                     m_workingSatEntry;
    int                                                     m_numInletThroats;
    bool                                                    m_binaryFiles;
    bool                                                    m_useAvrXOverThroatLen;
    bool                                                    m_useAvrPbcThroatLen;
    bool                                                    m_addPeriodicBC;
    vector< ThreeSome< PoreStruct*, int*, int* > >          m_poreData;
    vector< ThroatStruct* >                                 m_throatData;
    vector< ThreeSome< PoreStruct*, int*, int* > >          m_inletPores;
    vector< ThreeSome< PoreStruct*, int*, int* > >          m_outletPores;
    multimap< int, ThreeSome< int, int, double > >          m_outletConnections;
    multimap< int, ThreeSome< int, int, double > >          m_inletConnections;
    set< pair< int, int > >                                 m_xyPbcConn;
    set< pair< int, int > >                                 m_xzPbcConn;
    vector< ThreeSome< int, int, double > >                 m_pbcData;
    vector< int >                                           m_throatHash;
    vector< int >                                           m_reverseThroatHash;
    double                                                  m_origXDim;
    double                                                  m_origYDim;
    double                                                  m_origZDim;
    double                                                  m_averageThroatLength;
    double                                                  m_averagePoreHalfLength;
    double                                                  m_networkSeparation;
};

inline bool InputData::outletThroat(int throatIdx) const
{
    ThroatStruct *throat = m_throatData[throatIdx-1];
    return (throat->poreOne == 0 || throat->poreTwo == 0);
}

inline bool InputData::inletThroat(int throatIdx) const
{
    ThroatStruct *throat = m_throatData[throatIdx-1];
    return (throat->poreOne == -1 || throat->poreTwo == -1);
}

/**
// Echos the keywords and associated data to the given output stream
*/
/*inline void InputData::echoKeywords(ostream& out) const
{
    for(unsigned i = 0; i < m_parsedData.size(); ++i)
        out << m_parsedData[i].first    << endl
            << m_parsedData[i].second   << endl
            << endl;
}*/

/**
// Retrieves data sting based on supplied keyword, and connects a
// string stream to it. Data is removed from storage after having been
// retrived
*/
/*inline bool InputFile::getData(istringstream& data, const string& keyword) const
{
    for(unsigned i = 0; i < m_parsedData.size(); ++i)
        if(m_parsedData[i].first == keyword)
        {
            data.str(m_parsedData[i].second);
            return true;
        }
    return false;
}*/





inline void InputData::appendPoreData(int thisPoreIdx, int throatIdx, int thatPoreIdx)
{
    int oldConnNum = m_poreData[thisPoreIdx-1].first()->connNum;
    int *newThroatConnList = new int[oldConnNum+1];
    int *newPoreConnList = new int[oldConnNum+1];
    for(int i = 0; i < oldConnNum; ++i)
    {
        newPoreConnList[i] = m_poreData[thisPoreIdx-1].second()[i];
        newThroatConnList[i] = m_poreData[thisPoreIdx-1].third()[i];
    }
    newPoreConnList[oldConnNum] = thatPoreIdx;
    newThroatConnList[oldConnNum] = throatIdx;
    m_poreData[thisPoreIdx-1].first()->connNum++;
    delete[] m_poreData[thisPoreIdx-1].second();
    delete[] m_poreData[thisPoreIdx-1].third();
    m_poreData[thisPoreIdx-1].second(newPoreConnList);
    m_poreData[thisPoreIdx-1].third(newThroatConnList);
}

#endif

