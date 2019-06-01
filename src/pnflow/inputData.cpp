#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>
#include <algorithm>
#include <map>
#include <set>
using namespace std;

#include "threeSome.h"
#include "inputData.h"

/**
//
//  Available keywords are:
//
//X  TITLE               Title of run
//X  RAND_SEED           Seed for random number generator
//X  SAT_TARGET          Specifying saturation targets
//X  RELPERM_DEF         USe flowrate at single or residual saturation to calculate kr
//X  PRS_BDRS            What type of pressure boundaries should be used
//X  SAT_COMPRESS        Reduce saturation interval once kr drops below given threshold
//X  MAT_BAL             Calculates material balace of the fluids present
//X  APEX_PRS            Record tme minimum (drainage) pc observed in the last step
//X  POINT_SOURCE        Rather than inject across ine inlet face, inject from a single pore index
//X  SOLVER_TUNE         Tuning options for rel perm solver
//X  PORE_FILL_WGT       Weights for pore body filling mechnism
//X  PORE_FILL_ALG       Pore body filling mechnism algorithm
//X  MODIFY_RAD_DIST     Modify inscribed radii distribution
//X  MODIFY_G_DIST       Modify shape factor distribution
//X  MODIFY_PORO         Modify porosity
//X  MODIFY_CONN_NUM     Reduce connection number
//X  MODIFY_MOD_SIZE     Modify absolute model size
//X  SAT_COVERGENCE      Saturation covergence tuning
//X  TRAPPING            Trapping options
//X  FILLING_LIST        Create list for movie post processing
//X  OUTPUT              Create water sat map
//X  RES_FORMAT          Which output format to use
//X  GRAV_CONST          Vectorized gravity constant
//X  FLUID               Fluid properties
//X  CLAY_EDIT           Edit clay content
//X  A_CLOSE_SHAVE       Shave off boundaries to remove end-effects in net generation
//X  CALC_BOX            Where to set pressure boundaries for rel perm solver
//X  PRS_DIFF            Define the pressure differential across network for rel perm calculations
//X  NETWORK             The network to be used
//  NET_SERIES          Add networks in series
//  PERIODIC_PBC        Add periodic bc
//X  FRAC_CON_ANG        Fractional wetting options
//X  INIT_CON_ANG        Initial contact angle
//X  EQUIL_CON_ANG       Equilibration contact angle
//X  WRITE_NET           Write new network to file
//  SOLVER_DBG          Debugging options for solver
//
*/














const int       InputData::DUMMY_INDEX = -99;

InputData::InputData(const string& inputFileName)
:InputFile(inputFileName)
{
    m_workingSatEntry = 0;
    m_numInletThroats = 0;
    m_averageThroatLength = 0.0;
    m_networkSeparation = 0.0;
    m_connectionsRemoved = 0;
    m_useAvrXOverThroatLen = false;
    m_addPeriodicBC = false;
    m_useAvrPbcThroatLen = false;
}



/**
// All basic keywords are retived from the parsed storage. If not present the
// default values are used.
*/
string InputData::title() const
{
    istringstream data;
    
	
    if(getData(data, "TITLE"))
    {
        cout<< "Reading TITLE: "; cout.flush();
		string baseFileName;
        data >> baseFileName;
        cout<< baseFileName << endl;
        if(!data) errorMsg("TITLE");
        errorInDataCheck(data, "TITLE");
        return baseFileName;
    }
    else
        return m_baseFileName;
}

int InputData::randSeed()
{	int seedNum;
    istringstream data;
    string keyword("RAND_SEED");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> seedNum;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
        seedNum = (unsigned)time( NULL );
   
    return seedNum;
}


bool InputData::satTarget(double& requestedFinalSw, double& requestedFinalPc, double& deltaSw, double& deltaPc, double& deltaPcIncFactor,
                         bool& calcKr, bool& calcI, bool& entreL, bool& entreR, bool& exitL, bool& exitR)
{
    istringstream data;
    string keyword("SAT_CONTROL");
    char relPerm, resI,  entreLC, entreRC, exitLC, exitRC;

    if(getData(data, keyword))
    {
        for(int i = 0; i < m_workingSatEntry; ++i)
        {
            data >> requestedFinalSw;
            data >> requestedFinalPc >> deltaSw >> deltaPc >> deltaPcIncFactor >> relPerm >> resI>> entreLC>> entreRC>> exitLC>> exitRC;
        }

        if(data >> requestedFinalSw)
        {
            if (verbose) cout<< "Reading " << keyword << ",    "; cout.flush();
            data >> requestedFinalPc >> deltaSw >> deltaPc >> deltaPcIncFactor >> relPerm >> resI>> entreLC>> entreRC>> exitLC>> exitRC;
            calcKr = (relPerm == 'T' || relPerm == 't');
            calcI = (resI == 'T' || resI == 't');
            entreL = (entreLC == 'T' || entreLC == 't');
            entreR = (entreRC == 'T' || entreRC == 't');
            exitL = (exitLC == 'T' || exitLC == 't');
            exitR = (exitRC == 'T' || exitRC == 't');
            cout<<"Target Sw: "<<requestedFinalSw<<", Entry BC:" <<entreL << entreR<<", Exit BC: " <<exitL << exitR<<endl;
            if(!data) errorMsg(keyword);
            ++m_workingSatEntry;
            return true;
        }
    }
    else
        missingDataErr(keyword);

    return false;
}

void InputData::prsBdrs(bool& usePrsBdr, bool& reportPrsBdr, int& numPlanes) const
{
    istringstream data;
    string keyword("PRS_BDRS");
    char usePrs, numPl;

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> usePrs >> numPl >> numPlanes;
        usePrsBdr = (usePrs == 'T' || usePrs == 't');
        reportPrsBdr = (numPl == 'T' || numPl == 't');

        if(!reportPrsBdr) numPlanes = 0;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        usePrsBdr = false;
        reportPrsBdr = false;
        numPlanes = 0;
    }
}

void InputData::matBal(bool& reportMatBal) const
{
    istringstream data;
    string keyword("MAT_BAL");
    char mb;

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> mb;
        reportMatBal = (mb == 'T' || mb == 't');

        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        reportMatBal = false;
    }
}

void InputData::apexPrs(bool& doApexAnalysis)const
{
    istringstream data;
    string keyword("APEX_PRS");
    char apex;

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> apex;
        doApexAnalysis = (apex == 'T' || apex == 't');

        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        doApexAnalysis = false;
    }
}

bool InputData::sourceNode(int& sourceNode) const
{
    istringstream data;
    string keyword("POINT_SOURCE");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> sourceNode;

        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
        return true;

    }
    else
    {
        sourceNode = 0;
        return false;
    }
}

void InputData::solverTune(double& eps, double& scaleFact, int& slvrOutput, bool& verbose, double& condCutOff) const
{
    istringstream data;
    char verb('F');
    string keyword("SOLVER_TUNE");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> eps >> scaleFact >> slvrOutput >> verb >> condCutOff;
        verbose = (verb == 'T' || verb == 't');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        eps = 1.0E-15;
        scaleFact = 5;
        slvrOutput = 0;
        condCutOff = 0.0;
        verbose = false;
    }
}

/**
// Filling weights:
// Oren1/2 = 0.0, 0.5, 1.0, 5.0, 10.0, 50.0
// Blunt2 = 0.0, 15000, 15000, 15000, 15000, 15000
// Blunt1 = 0.0, 50E-6, 50E-6, 100E-6, 200E-6, 500E-6
*/
void InputData::poreFillWgt(vector< double >& weights) const
{
    istringstream data;
    string keyword("PORE_FILL_WGT");
    weights.resize(6);

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword ; cout.flush();
        data >> weights[0] >> weights[1] >> weights[2] >> weights[3] >> weights[4] >> weights[5];
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        weights[0] = 0.0;
        weights[1] = 15000.0;
        weights[2] = 15000.0;
        weights[3] = 15000.0;
        weights[4] = 15000.0;
        weights[5] = 15000.0;
        cout<< "Using default pore filling weights" ; cout.flush();

    }
   cout<< " " << endl;

}

void InputData::poreFillAlg(string& algorithm) const
{
    istringstream data;
    string keyword("PORE_FILL_ALG");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> algorithm;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        algorithm = "blunt2";
    }
}

void InputData::relPermDef(string& flowRef/*, bool& strictTrpCond*/) const
{
    istringstream data;
    string keyword("RELPERM_DEF");
    //char cond;

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> flowRef ;
        //strictTrpCond = (cond == 't' || cond == 'T');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        flowRef = "single";
        //strictTrpCond = true;
    }
}

void InputData::getModifyRadDistOptions(int& throatModel, int& poreModel, string& throatOptions, string& poreOptions,
                              bool& maintainLtoR, bool& writeDistToFile, int& numPtsRDist)const
{
    istringstream data;
    char toFile, maintainAR;
    string keyword("MODIFY_RAD_DIST");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        getRadDist(data, throatModel, throatOptions);
        getRadDist(data, poreModel, poreOptions);
        data >> maintainAR >> toFile >> numPtsRDist;
        writeDistToFile = (toFile == 't' || toFile == 'T');
        maintainLtoR = (maintainAR == 't' || maintainAR == 'T');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        throatModel = 0;
        poreModel = 0;
        writeDistToFile = false;
        maintainLtoR = false;
    }
}

void InputData::getModifyGDist(int& throatModel, int& poreModel, string& throatOptions, string& poreOptions,
                            bool& writeDistToFile, int& numPts) const
{
    istringstream data;
    char toFile;
    string keyword("MODIFY_G_DIST");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        getRadDist(data, throatModel, throatOptions);
        getRadDist(data, poreModel, poreOptions);
        data >> toFile >> numPts;
        writeDistToFile = (toFile == 't' || toFile == 'T');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
   }
    else
    {
        throatModel = 0;
        poreModel = 0;
        writeDistToFile = false;
    }
}

void InputData::getModifyPoro(double& netPoroTrgt, double& clayPoroTrgt) const
{
    istringstream data;
    string keyword("MODIFY_PORO");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> netPoroTrgt >> clayPoroTrgt;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    //else
    //{
        //netPoroTrgt = -1.0;
        //clayPoroTrgt = -1.0;
    //}
}

void InputData::modifyConnNum(double& targetConnNum, string& model)const
{
    istringstream data;
    string keyword("MODIFY_CONN_NUM");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> targetConnNum >> model;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        targetConnNum = -1.0;
    }
}

void InputData::getModifyModelSize(double& scaleFactor) const
{
    istringstream data;
    string keyword("MODIFY_MOD_SIZE");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> scaleFactor;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
   }
    else
    {
        scaleFactor = -1.0;
    }
}


void InputData::getRadDist(istream& data, int& model, string& options) const
{
    ostringstream optionStr;
    data >> model;

    //if(model == 1)
    //{
        //string fileName;
        //double lowCutOff, highCutOff;
        //data >> fileName >> lowCutOff >> highCutOff;
        //optionStr << fileName << "  " << lowCutOff << "  " << highCutOff;
    //}
    //else if(model == 2 || model == 5)
    //{
        //double min, max, delta, eta;
        //data >> min >> max >> delta >> eta;
        //optionStr << min << "  " << max << "  " << delta << "  " << eta;
    //}
    //else if(model > 2 && model < 4)
    //{
        //double tmp;
        //data >> tmp;
        //optionStr << tmp;
    //}
    //else if(model == 4)
    //{
        //double aConst;
        //char above, below;
        //data >> aConst >> above >> below;
        //optionStr << aConst <<" "<< above << below;
        ////cout<< aConst << above << below;
    //}
    //else
    {
		data.get (*(optionStr.rdbuf()));
	}

    options = optionStr.str();
}


void InputData::satConvergence(int& minNumFillings, double& initStepSize,
      double& cutBack, double& maxIncrFact, bool& stable) const
{
    istringstream data;
    string keyword("SAT_COVERGENCE");
    char stab('F');

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> minNumFillings >> initStepSize >> cutBack >> maxIncrFact >> stab;
        stable = (stab == 'T' || stab == 't');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        stable = false;
        minNumFillings = 10;
        initStepSize = 0.02;
        cutBack = 0.8;
        maxIncrFact = 2.0;
    }
}

/* //void InputData::trapping(bool& drainEnds, double& KrwatcornAtSw0) const
//{
    //istringstream data;
    //string keyword("TRAPPING");
    //char drnEnds;//fromEntry, fromExit, 
//
    //if(getData(data, keyword))
    //{
        //cout<< "Reading " << keyword << endl;
        //data >> drnEnds >> KrwatcornAtSw0;
        injEntry = (fromEntry == 'T' || fromEntry == 't');
        injExit = (fromExit == 'T' || fromExit == 't');
        //drainEnds = (drnEnds == 'T' || drnEnds == 't');
        //if(!data) errorMsg(keyword);
        //errorInDataCheck(data, keyword);
    //}
    //else
    //{
        injEntry = true;
        injExit = false;
        //drainEnds = true;
        //KrwatcornAtSw0 = 0.0;
    //}
//}*/

void InputData::fillingList(bool& isDrainage, bool& imbibition, bool& location) const
{
    istringstream data;
    string keyword("FILLING_LIST");
    char drain, imb, loc;

    if(getData(data, keyword))
    {
		        cout<< "Reading " << keyword<<"\n ERROR disabled\n  ERROR \n  ERROR \n  ERROR \n  ERROR \n  ERROR \n " << endl;

        if (verbose) cout<< "Reading " << keyword << endl;
        data >> drain >> imb >> loc;
        isDrainage = (drain == 'T' || drain == 't');
        imbibition = (imb == 'T' || imb == 't');
        location = (loc == 'T' || loc == 't');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
   }
    else
    {
       isDrainage = false;
       imbibition = false;
       location = false;
    }
}

void InputData::output(bool& propOut, bool& swOut) const
{
    istringstream data;
    string keyword("OUTPUT");
    char prop, sw;

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> prop >> sw;
        propOut = (prop == 'T' || prop == 't');
        swOut = (sw == 'T' || sw == 't');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
       prop = false;
       sw = false;
    }
}

void InputData::resFormat(bool& matlabFormat, bool& excelFormat, bool& mcpFormat) const
{
    istringstream data;
    string keyword("RES_FORMAT");
    string resForm;

    if(getData(data, keyword))
    {
        if (verbose)  cout<< "Reading " << keyword << endl;
		data >> resForm;
        matlabFormat = (resForm == "MATLAB" || resForm == "matlab" || resForm == "Matlab");
        excelFormat = (resForm == "EXCEL" || resForm == "excel" || resForm == "Excel");
        if (!excelFormat) excelFormat = (resForm == "excelAndMicroPorosity" || resForm == "EXCELANDMICROPOROSITY" || resForm == "ExcelAndMicroPorosity");
        if (excelFormat) 
			mcpFormat = (resForm == "excelAndMicroPorosity" || resForm == "EXCELANDMICROPOROSITY" || resForm == "ExcelAndMicroPorosity");
        else
			mcpFormat = (resForm == "upscaling" || resForm == "UPSCALING");
        
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        matlabFormat = false;
        excelFormat = false;
    }
}

void InputData::gravityConst(double& gravX, double& gravY, double& gravZ)const
{
    istringstream data;
    string keyword("GRAV_CONST");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> gravX >> gravY >> gravZ;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        gravX = 0.0;
        gravY = 0.0;
        gravZ = -9.81;
    }
}

void InputData::fluid(double& intfacTen, double& watVisc, double& oilVisc, double& watResist,
                      double& oilResist, double& watDens, double& oilDens) const
{
    istringstream data;
    string keyword("FLUID");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> intfacTen >> watVisc >> oilVisc >> watResist >> oilResist >> watDens >> oilDens;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);

        intfacTen *= 1.0E-3;                       // Get it into N/m
        watVisc *= 1.0E-3;                          // Get it into Pa.s
        oilVisc *= 1.0E-3;                          // Get it into Pa.s
    }
    else
    {
        intfacTen = 30.0E-3;
        watVisc = 1.0E-3;
        oilVisc = 1.0E-3;
        watResist = 1.0;
        oilResist = 1000.0;
        watDens = 1000.0;
        oilDens = 1000.0;
    }
}

void InputData::aCloseShave(double& shaveSetting) const
{
    istringstream data;
    string keyword("A_CLOSE_SHAVE");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> shaveSetting;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        shaveSetting = 1.0;
    }
}

void InputData::relPermCompression(bool& useComp, double& krThres,
         double& deltaSw, bool& wettPhase, bool& nonWettPhase) const
{
    istringstream data;
    string keyword("SAT_COMPRESS");
    char wettP('F'), nonWettP('F');

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        useComp = true;
        data >> krThres >> deltaSw >> wettP >> nonWettP;
        wettPhase = (wettP == 'T' || wettP == 't');
        nonWettPhase = (nonWettP == 'T' || nonWettP == 't');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        useComp = false;
        krThres = 0.0;
        deltaSw = 0.1;
    }
}

double InputData::clayEdit() const
{
	double clayEditFact;
    istringstream data;
    string keyword("CLAY_EDIT");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> clayEditFact;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
	}
    else
    {
        clayEditFact = 0.0;
    }
    return clayEditFact;
}

void InputData::calcBox(double& inletBdr, double& outletBdr) const
{
    istringstream data;
    string keyword("CALC_BOX");

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> inletBdr >> outletBdr;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
   }
    else
    {
        inletBdr = 0.5;
        outletBdr = 1.0;
    }
}

void InputData::prsDiff(double& inletPrs, double& outletPrs, bool& useGravInKr) const
{
    istringstream data;
    string keyword("PRS_DIFF");
    char grav;

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> inletPrs >> outletPrs >> grav;
        useGravInKr = (grav == 'T' || grav == 't');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        inletPrs = 1.0;
        outletPrs = 0.0;
        useGravInKr = false;
    }
}


/**
// The nework data is supplied back as the pores and throats are created. The input
// streams are intially just opened and the headers read.
*/
void InputData::network(int& numPores, int& numThroats, double& xDim, double& yDim, double& zDim)
{
    istringstream data, dataNet, dataPbc;
    //string m_m_netFileBase;
    string keyword("NETWORK");
    string keywordNet("NET_SERIES");
    string keywordPbc("PERIODIC_BC");

    if(getData(dataPbc, keywordPbc))
    {
        if (verbose) cout<< "Reading " << keywordPbc << endl;
        char usePbc, avrLen;
        dataPbc >> usePbc >> avrLen;
        m_useAvrPbcThroatLen = (avrLen == 'T' || avrLen == 't');
        m_addPeriodicBC = (usePbc == 'T' || usePbc == 't');
        if(!dataPbc) errorMsg(keywordPbc);
        errorInDataCheck(dataPbc, keywordPbc);
    }
    else
    {
        m_addPeriodicBC = false;
        m_useAvrPbcThroatLen = false;
    }

    if(getData(dataNet, keywordNet))
    {
        if (verbose) cout<< "Reading " << keywordNet << endl;
        char avrLen;
        dataNet >> m_numNetInSeries >> avrLen >> m_networkSeparation;
        m_useAvrXOverThroatLen = (avrLen == 'T' || avrLen == 't');
        if(!dataNet) errorMsg(keywordNet);
        errorInDataCheck(dataNet, keywordNet);
    }
    else
    {
        m_numNetInSeries = 1;
        m_useAvrXOverThroatLen = false;
    }

    char binFile;
    if(getData(data, keyword))
    {
		std::string line;
		std::getline(data, line);
	    std::istringstream iss(line);
        if (verbose) cout<< "Reading " << keyword << endl;
        iss >> binFile >> m_netFileBase;
        m_binaryFiles = (binFile == 'T' || binFile == 't');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
        missingDataErr(keyword);

    if(m_binaryFiles)
    {
        string porePropFile(m_netFileBase + "_node.bin");
        m_poreProp.open(porePropFile.c_str(), ios::binary);

        string throatPropFile(m_netFileBase + "_link.bin");
        m_throatProp.open(throatPropFile.c_str(), ios::binary);

        m_poreProp.read((char *)(&numPores), sizeof(int));
        m_poreProp.read((char *)(&xDim), sizeof(double));
        m_poreProp.read((char *)(&yDim), sizeof(double));
        m_poreProp.read((char *)(&zDim), sizeof(double));
        m_throatProp.read((char *)(&numThroats), sizeof(int));
    }
    else
    {
        string poreConnFile(m_netFileBase + "_node1.dat");          // Open file containing pore connection data
        m_poreConn.open(poreConnFile.c_str());

        string porePropFile(m_netFileBase + "_node2.dat");          // Open file containing pore geometry data
        m_poreProp.open(porePropFile.c_str());

        string throatConnFile(m_netFileBase + "_link1.dat");        // Open file containing throat connection data
        m_throatConn.open(throatConnFile.c_str());

        string throatPropFile(m_netFileBase + "_link2.dat");        // Open file containing throat geometry data
        m_throatProp.open(throatPropFile.c_str());


        m_poreConn >> numPores >> xDim >> yDim >> zDim;
        m_throatConn >> numThroats;
    }
		if (!m_poreProp || !m_throatProp || (!m_binaryFiles && (!m_poreConn || !m_throatConn)))
		{
			std::cout<< "======================================== " << std::endl
				<<"Error: Unable to open network data files " << m_netFileBase << std::endl
				<<" failed on !"<< bool(m_poreProp) <<"  || !"<< bool(m_throatProp)  <<" || (!"<< bool(m_binaryFiles) <<" && (!"<<bool(m_poreConn)<<" || !"<<bool(m_throatConn)<<"))" << std::endl
				<<"            "<< !m_poreProp <<"  ||  "<< !m_throatProp  <<" || "<< (!m_binaryFiles && (!m_poreConn || !m_throatConn)) << std::endl
				<< "======================================== " << endl;
			exit( -1 );
		}
		outD<<"network:" << m_netFileBase <<(binFile?" binary ":" ascii ")<<" Nt:"<<numThroats<<" Np:"<<numPores<< std::endl;
    m_origNumPores = numPores;
    m_origNumThroats = numThroats;
    m_origXDim = xDim;
    m_origYDim = yDim;
    m_origZDim = zDim;



    if (m_numNetInSeries < 1)
    {
        cerr << "==================================================" << endl
            << "Error: There need to be at least one net in series" << endl
            << "==================================================" << endl;
        exit( -1 );
    }

    loadPoreData();
    loadThroatData();

    numPores *= m_numNetInSeries;
    numThroats = m_origNumThroats*m_numNetInSeries - m_connectionsRemoved*(m_numNetInSeries-1);
    xDim = xDim*static_cast< double >(m_numNetInSeries) + static_cast< double >(m_numNetInSeries-1)*m_networkSeparation;
}

/**
// As the pores are created the data is read from file and supplied back together with contact angles that are
// stored in vectors.
//
// The format of pore network files are:
// *_node1.dat:
// index, x_pos, y_pos, z_pos, connection num, connecting nodes..., at inlet?, at outlet?, connecting links...
//
// *_node2.dat:
// index, volume, radius, shape factor, clay volume
*/
void InputData::loadPoreData()
{
    m_poreData.resize(m_origNumPores);
    for(int i = 0; i < m_origNumPores; ++i)
    {
        if ((!m_binaryFiles && !m_poreConn) || !m_poreProp)
        {
            cerr << "=========================" << endl
                 << "Error while reading pores." << endl
                 << i << endl
                 << (!m_binaryFiles && !m_poreConn) << " || "<< !m_poreProp << endl
                 << "==========================" << endl;
            exit( -1 );
        }

        PoreStruct *poreProp = new PoreStruct;
        int *connThroats, *connPores;

        if(m_binaryFiles)
        {
            m_poreProp.read((char *)(poreProp), sizeof(*poreProp));

            connThroats = new int[poreProp->connNum];
            connPores = new int[poreProp->connNum];
            m_poreProp.read((char *)(connPores), poreProp->connNum*sizeof(int));
            m_poreProp.read((char *)(connThroats), poreProp->connNum*sizeof(int));
        }
        else
        {
            int idx;
            bool isAtInletRes, isAtOutletRes;

            m_poreProp >> poreProp->index
                >> poreProp->volume
                >> poreProp->radius
                >> poreProp->shapeFact
                >> poreProp->clayVol;

            m_poreConn >> idx
                >> poreProp->x
                >> poreProp->y
                >> poreProp->z
                >> poreProp->connNum;

            ensure(idx == poreProp->index);

            connThroats = new int[poreProp->connNum];
            connPores = new int[poreProp->connNum];

            for(int k = 0; k < poreProp->connNum; ++k)
                m_poreConn >> connPores[k];

            m_poreConn >> isAtInletRes >> isAtOutletRes;

            for(int j = 0; j < poreProp->connNum; ++j)
                m_poreConn >> connThroats[j];
        }
        ThreeSome< PoreStruct*, int*, int* > elem(poreProp, connPores, connThroats);
        m_poreData[i] = elem;

        for(int conn = 0; conn < poreProp->connNum; ++conn)
        {
            if(connPores[conn] == -1)
                m_inletPores.push_back(elem);
            else if(connPores[conn] == 0)
                m_outletPores.push_back(elem);
        }
    }

    if(m_addPeriodicBC) findBoundaryPores();
    m_poreProp.close();
    m_poreConn.close();
}

void InputData::findBoundaryPores()
{
    int nDir = static_cast< int>(pow(m_origNumPores, 1.0/3.0))+1;
    double xStep(m_origXDim/nDir), yStep(m_origYDim/nDir), zStep(m_origZDim/nDir);
    for(int i = 0; i < nDir; ++i)
    {
        double xPos = xStep/2.0 + i*xStep;
        ensure(xPos < m_origXDim);
        for(int j = 0; j < nDir; ++j)
        {
            double yPos = yStep/2.0 + j*yStep;
            double zPos = zStep/2.0 + j*zStep;
            ensure(xPos < m_origXDim && yPos < m_origYDim && zPos < m_origZDim);
            vector< ThreeSome< double, double, double > > positonData;
            ThreeSome< double, double, double > xyMinus(xPos, yPos, 0.0);
            positonData.push_back(xyMinus);
            ThreeSome< double, double, double > xyPluss(xPos, yPos, m_origZDim);
            positonData.push_back(xyPluss);
            ThreeSome< double, double, double > xzMinus(xPos, 0.0, zPos);
            positonData.push_back(xzMinus);
            ThreeSome< double, double, double > xzPluss(xPos, m_origYDim, zPos);
            positonData.push_back(xzPluss);

            vector< int > poreIndecies;
            vector< double > p2BdrLength;
            findClosestPoreForPBC(positonData, poreIndecies, p2BdrLength);
            pair< int, int > xyConn(poreIndecies[0], poreIndecies[1]);
            pair< int, int > xzConn(poreIndecies[2], poreIndecies[3]);

            if(xyConn.first != xyConn.second && m_xyPbcConn.count(xyConn) == 0)  // Detect 2D networks and duplicates
            {
                ThreeSome<int, int, double> pbcConn(xyConn.first, xyConn.second, p2BdrLength[0]+p2BdrLength[1]);
                m_pbcData.push_back(pbcConn);
            }

            if(xzConn.first != xzConn.second && m_xzPbcConn.count(xzConn) == 0)
            {
                ThreeSome<int, int, double> pbcConn(xzConn.first, xzConn.second, p2BdrLength[2]+p2BdrLength[3]);
                m_pbcData.push_back(pbcConn);
            }

            m_xyPbcConn.insert(xyConn);
            m_xzPbcConn.insert(xzConn);
        }
    }
}

void InputData::findClosestPoreForPBC(const vector< ThreeSome< double, double, double > >& positonData,
                                      vector< int >& poreIndecies, vector< double >& p2BdrLength) const
{
    poreIndecies.resize(positonData.size());
    p2BdrLength.resize(positonData.size(), 1.0E21);
    for(size_t i = 0; i < m_poreData.size(); ++i)
    {
        for(size_t j = 0; j < positonData.size(); ++j)
        {
            double len = sqrt(pow(positonData[j].first()-m_poreData[i].first()->x, 2.0)+
                pow(positonData[j].second()-m_poreData[i].first()->y, 2.0)+
                pow(positonData[j].third()-m_poreData[i].first()->z, 2.0));
            if(len < p2BdrLength[j])
            {
                p2BdrLength[j] = len;
                poreIndecies[j] = static_cast< int >(i)+1;
            }
        }
    }
}

void InputData::getOutletData(MapItr itr, int numNetsInFront, int& throatIdx, int& thatIdx) const
{
    ThreeSome<int,int,double> entry = (*itr).second;
    thatIdx = entry.second() + m_origNumPores*(numNetsInFront+1);
    bool throatToOutlet = outletThroat(entry.first());
    bool throatToInlet = inletThroat(entry.first());
    int hashedIdx(entry.first());
    if(numNetsInFront > 0 && throatToOutlet)
    {
        hashedIdx = m_reverseThroatHash[hashedIdx-1] +
            numNetsInFront*m_origNumThroats - (numNetsInFront-1)*m_connectionsRemoved;
    }
    else if(throatToInlet)
    {
        ensure(m_reverseThroatHash[hashedIdx-1] != DUMMY_INDEX);
        hashedIdx = m_reverseThroatHash[hashedIdx-1] +
            (numNetsInFront+1)*m_origNumThroats - numNetsInFront*m_connectionsRemoved;
    }
    throatIdx = hashedIdx;
}

void InputData::getInletData(MapItr itr, int numNetsInFront, int& throatIdx, int& thatIdx) const
{
    ThreeSome<int,int,double> entry = (*itr).second;
    thatIdx = entry.second() + m_origNumPores*(numNetsInFront-1);
    bool throatToOutlet = outletThroat(entry.first());
    bool throatToInlet = inletThroat(entry.first());
    int hashedIdx(entry.first());
    if(numNetsInFront > 1 && throatToOutlet)
    {
        hashedIdx = m_reverseThroatHash[hashedIdx-1] +
            (numNetsInFront-1)*m_origNumThroats - (numNetsInFront-2)*m_connectionsRemoved;
    }
    else if(throatToInlet)
    {
        hashedIdx = m_reverseThroatHash[hashedIdx-1] +
            numNetsInFront*m_origNumThroats - (numNetsInFront-1)*m_connectionsRemoved;
    }
    throatIdx = hashedIdx;
}

void InputData::poreLocation(int idx, double& xPos) const
{
    int stdIdx(idx), numNetsInFront(0);
    if(idx > m_origNumPores)
    {
		cout<<"Error(?)";
        numNetsInFront = (idx-1)/m_origNumPores;
        stdIdx = idx - numNetsInFront*m_origNumPores;
    }

    PoreStruct *poreProp = m_poreData[stdIdx-1].first();
    xPos = poreProp->x + numNetsInFront*(m_origXDim+m_networkSeparation);
}

/**
// There is a lot of memory assosiated with storing all network data.
// Best to clean up after ourself before proceeding
*/
void InputData::clearNetworkData()
{
    for(size_t i = 0; i < m_poreData.size(); ++i)
    {
        delete[] m_poreData[i].second();
        delete[] m_poreData[i].third();
        delete m_poreData[i].first();
    }

    for(size_t j = 0; j < m_throatData.size(); ++j)
    {
        delete m_throatData[j];
    }

    m_poreData.clear();
    m_throatData.clear();
    m_inletPores.clear();
    m_outletPores.clear();
    m_outletConnections.clear();
    m_inletConnections.clear();
    m_throatHash.clear();
    m_reverseThroatHash.clear();
}

/**
// As the throats are created the data is read from file and supplied back together with contact angles
// that are stored in vectors.
//
// The format of throat network files are:
// *_link1.dat:
// index, pore 1 index, pore 2 index, radius, shape factor, total length (pore center to pore center)
//
// *_link2.dat:
// index, pore 1 index, pore 2 index, length pore 1, length pore 2, length throat, volume, clay volume
*/
void InputData::loadThroatData()
{
    m_throatData.resize(m_origNumThroats);
    double lenSumPore(0.0), lenSumThroat(0.0);
    for(int i = 0; i < m_origNumThroats; ++i)
    {
        if ((!m_binaryFiles && !m_throatConn) || !m_throatProp)
        {
            cerr << "=========================== " << endl
                << "Error while reading throats." << endl
                << "throat index: "<<i << endl
                << bool(m_throatConn)<<bool(m_throatProp) << endl
                << "============================ " << endl;
            exit( -1 );
        }

        ThroatStruct *throatProp = new ThroatStruct;

        if(m_binaryFiles)
        {
            m_throatProp.read((char *)(throatProp), sizeof(*throatProp));
        }
        else
        {
            int idx, tmp;

            m_throatConn >> throatProp->index
                >> throatProp->poreOne
                >> throatProp->poreTwo
                >> throatProp->radius
                >> throatProp->shapeFact
                >> throatProp->lenTot;

            m_throatProp >> idx
                >> tmp
                >> tmp
                >> throatProp->lenPoreOne
                >> throatProp->lenPoreTwo
                >> throatProp->lenThroat
                >> throatProp->volume
                >> throatProp->clayVol;

            ensure(idx == throatProp->index);
        }
        m_throatData[i] = throatProp;
        lenSumThroat += throatProp->lenThroat;
        lenSumPore += (throatProp->lenPoreOne+throatProp->lenPoreTwo)/2.0;
        if(throatProp->poreOne == -1 || throatProp->poreTwo == -1) ++m_numInletThroats;
    }
    m_throatConn.close();
    m_throatProp.close();
    m_averageThroatLength = lenSumThroat/m_origNumThroats;
    m_averagePoreHalfLength = lenSumPore/m_origNumThroats;

    if(m_addPeriodicBC)
    {
        int index = m_origNumThroats+1;
        for(size_t conn = 0; conn < m_pbcData.size(); ++conn)
        {
            ThroatStruct *pbcThroat = new ThroatStruct;
            double randNum = double(rand()) / double(RAND_MAX);
            //double randNum = 0.5;   // delete me
            int randThroatIdx = (randNum*m_origNumThroats);
            if(randThroatIdx >= m_origNumThroats) randThroatIdx = m_origNumThroats-1;
            pbcThroat->clayVol = 0.0;
            pbcThroat->index = index++;
            pbcThroat->poreOne = m_pbcData[conn].first();
            pbcThroat->poreTwo = m_pbcData[conn].second();
            pbcThroat->radius = m_throatData[randThroatIdx]->radius;
            pbcThroat->shapeFact = m_throatData[randThroatIdx]->shapeFact;
            pbcThroat->volume = 0.0;
            if(m_useAvrPbcThroatLen || m_pbcData[conn].third() == 0.0)
            {
                pbcThroat->lenPoreOne = m_averagePoreHalfLength;
                pbcThroat->lenPoreTwo = m_averagePoreHalfLength;
                pbcThroat->lenThroat = m_averageThroatLength;
                pbcThroat->lenTot = m_averageThroatLength+2.0*m_averagePoreHalfLength;
            }
            else
            {
                pbcThroat->lenTot = m_pbcData[conn].third();
                pbcThroat->lenThroat = m_pbcData[conn].third() *
                    (m_averageThroatLength/(m_averageThroatLength+2.0*m_averagePoreHalfLength));
                pbcThroat->lenPoreOne = (pbcThroat->lenTot-pbcThroat->lenThroat)/2.0;
                pbcThroat->lenPoreTwo = (pbcThroat->lenTot-pbcThroat->lenThroat)/2.0;
            }
            m_throatData.push_back(pbcThroat);
            appendPoreData(pbcThroat->poreOne, pbcThroat->index, pbcThroat->poreTwo);
            appendPoreData(pbcThroat->poreTwo, pbcThroat->index, pbcThroat->poreOne);
        }
        m_origNumThroats += m_pbcData.size();
    }

    set<int> addedInletThroats;

    if(m_numNetInSeries > 1)
    {
        for(int outT = 0; outT < m_origNumThroats; ++outT)
        {
            if(m_throatData[outT]->poreOne == 0 || m_throatData[outT]->poreTwo == 0)
            {
                int outFacePore = m_throatData[outT]->poreOne == 0 ? m_throatData[outT]->poreTwo: m_throatData[outT]->poreOne;
                double xPos = m_poreData[outFacePore-1].first()->x - m_origXDim - m_networkSeparation;
                double yPos = m_poreData[outFacePore-1].first()->y;
                double zPos = m_poreData[outFacePore-1].first()->z;
                double p2pLength(0.0);
                int inFacePore = findClosestPore(m_inletPores, xPos, yPos, zPos, p2pLength);
                ThreeSome<int, int, double> outEntry(outT+1, inFacePore, p2pLength);
                ThreeSome<int, int, double> inEntry(outT+1, outFacePore, p2pLength);
                m_inletConnections.insert(MultiConnValType(inFacePore, inEntry));
                m_outletConnections.insert(MultiConnValType(outFacePore, outEntry));
            }
        }

        for(int inT = 0; inT < m_origNumThroats; ++inT)
        {
            if(m_throatData[inT]->poreOne == -1 || m_throatData[inT]->poreTwo == -1)
            {
                int inFacePore = m_throatData[inT]->poreOne == -1 ? m_throatData[inT]->poreTwo: m_throatData[inT]->poreOne;
                if(m_inletConnections.count(inFacePore) == 0)
                {
                    addedInletThroats.insert(inT+1);
                    double xPos = m_poreData[inFacePore-1].first()->x + m_origXDim + m_networkSeparation;
                    double yPos = m_poreData[inFacePore-1].first()->y;
                    double zPos = m_poreData[inFacePore-1].first()->z;
                    double p2pLength(0.0);
                    int outFacePore = findClosestPore(m_outletPores, xPos, yPos, zPos, p2pLength);
                    ThreeSome<int, int, double> outEntry(inT+1, inFacePore, p2pLength);
                    ThreeSome<int, int, double> inEntry(inT+1, outFacePore, p2pLength);
                    m_inletConnections.insert(MultiConnValType(inFacePore, inEntry));
                    m_outletConnections.insert(MultiConnValType(outFacePore, outEntry));
                }
                else
                    ++m_connectionsRemoved;
            }
        }
    }

    int runningIndex(-1);
    m_throatHash.resize(m_origNumThroats-m_connectionsRemoved);      // Since inlet throats do not exist in serial
    for(int j = 1; j <=  m_origNumThroats; ++j)                       // nets, we end up having to take into account
    {                                                               // an offset for the indecies
        if((m_throatData[j-1]->poreOne != -1 && m_throatData[j-1]->poreTwo != -1) || addedInletThroats.count(j) > 0)
            m_throatHash[++runningIndex] = j;
    }

    runningIndex = 1;
    m_reverseThroatHash.resize(m_origNumThroats);
    for(int k = 0; k < m_origNumThroats; ++k)
    {
        int hashedIdx(DUMMY_INDEX);
        if((m_throatData[k]->poreOne != -1 && m_throatData[k]->poreTwo != -1) || addedInletThroats.count(k+1) > 0)
            hashedIdx = runningIndex++;

        m_reverseThroatHash[k] = hashedIdx;
    }
}

int InputData::findClosestPore(const vector< ThreeSome< PoreStruct*, int*, int* > > pores, double xPos,
                               double yPos, double zPos, double& totalLen) const
{
    totalLen = 1.0E21;
    int index(DUMMY_INDEX);
    for(size_t i = 0; i < pores.size(); ++i)
    {
        double len = sqrt(pow(xPos-pores[i].first()->x, 2.0)+
            pow(yPos-pores[i].first()->y, 2.0)+
            pow(zPos-pores[i].first()->z, 2.0));
        if(len < totalLen)
        {
            totalLen = len;
            index = pores[i].first()->index;
        }
    }
    ensure(index != DUMMY_INDEX);
    return index;
}



void InputData::poreData(int idx, double& xPos, double& yPos, double& zPos, int& connNum, vector< int >& connThroats,
                         vector< int >& connPores, double& vol, double& volCl, double& rad, double& shapeFact)
{
    int stdIdx(idx), numNetsInFront(0);
    if(idx > m_origNumPores)
    {
        numNetsInFront = (idx-1)/m_origNumPores;
        stdIdx = idx - numNetsInFront*m_origNumPores;
    }

    const PoreStruct *poreProp = m_poreData[stdIdx-1].first();
    xPos = poreProp->x + numNetsInFront*(m_origXDim+m_networkSeparation);
    yPos = poreProp->y;
    zPos = poreProp->z;
    connNum = poreProp->connNum;
    vol = poreProp->volume;
    volCl = poreProp->clayVol;
    rad = poreProp->radius;
    shapeFact = poreProp->shapeFact;
    connThroats.resize(connNum);
    connPores.resize(connNum);
    bool outletPore(false), inletPore(false);
    for(int i = 0; i < connNum; ++i)
    {
        connPores[i] = m_poreData[stdIdx-1].second()[i];
        connThroats[i] = m_poreData[stdIdx-1].third()[i];

        if(connPores[i] != 0 && connPores[i] != -1 && numNetsInFront > 0)
        {
            connPores[i] += numNetsInFront*m_origNumPores;
            int throatStdIdx = m_reverseThroatHash[connThroats[i]-1];
            connThroats[i] = throatStdIdx + numNetsInFront*m_origNumThroats - (numNetsInFront-1)*m_connectionsRemoved;
        }
        else if(connPores[i] == 0)
            outletPore = true;
        else if(connPores[i] == -1)
            inletPore = true;
    }

    if(m_numNetInSeries > 1 && (outletPore || inletPore))
    {
        //int  numHookedUp(0);//numInst(0),
        pair<MapItr,MapItr> mapConns;

        if(outletPore)
        {
            //numInst = (m_outletConnections.count(stdIdx));
            mapConns = m_outletConnections.equal_range(stdIdx);
        }
        else
        {
            //numInst = (m_inletConnections.count(stdIdx));
            mapConns = m_inletConnections.equal_range(stdIdx);
        }
        MapItr itr = mapConns.first;

        for(int j = 0; j < connNum; ++j)
        {
            if(numNetsInFront < m_numNetInSeries-1 && connPores[j] == 0)
            {
                ensure(itr != mapConns.second);
                getOutletData(itr, numNetsInFront, connThroats[j], connPores[j]);

                //++numHookedUp;
                ++itr;
            }
            else if(numNetsInFront > 0 && connPores[j] == -1)
            {
                ensure(itr != mapConns.second);
                getInletData(itr, numNetsInFront, connThroats[j], connPores[j]);
                //++numHookedUp;
                ++itr;
            }
            else if(numNetsInFront > 0 && connPores[j] == 0)
            {
                int throatStdIdx = m_reverseThroatHash[connThroats[j]-1];
                connThroats[j] = throatStdIdx + numNetsInFront*m_origNumThroats - (numNetsInFront-1)*m_connectionsRemoved;
            }
        }

        while(itr != mapConns.second)
        {
            int poreIdx(DUMMY_INDEX), throatIdx(DUMMY_INDEX);

            if(numNetsInFront < m_numNetInSeries-1 && outletPore)
                getOutletData(itr, numNetsInFront, throatIdx, poreIdx);
            else if(numNetsInFront > 0 && inletPore)
                getInletData(itr, numNetsInFront, throatIdx, poreIdx);

            if(poreIdx != DUMMY_INDEX)
            {
                connPores.push_back(poreIdx);
                connThroats.push_back(throatIdx);
                ++connNum;
            }
            ++itr;
            //++numHookedUp;
        }
    }
}


void InputData::throatData(int idx, int& pore1, int& pore2, double& vol, double& volCl, double& rad, double& shapeFact,
                           double& lenPore1, double& lenPore2, double& lenThroat, double& lenTot)
{
    int stdIdx(idx), numNetsInFront(0);
    if(idx > m_origNumThroats)
    {
        numNetsInFront = 1+(idx-m_origNumThroats-1)/(m_origNumThroats-m_connectionsRemoved);
        int tmpIdx = idx - numNetsInFront*m_origNumThroats + (numNetsInFront-1)*m_connectionsRemoved;
        stdIdx = m_throatHash[tmpIdx-1];
    }

    pore1 = m_throatData[stdIdx-1]->poreOne;
    pore2 = m_throatData[stdIdx-1]->poreTwo;
    vol = m_throatData[stdIdx-1]->volume;
    volCl = m_throatData[stdIdx-1]->clayVol;
    rad = m_throatData[stdIdx-1]->radius;
    shapeFact = m_throatData[stdIdx-1]->shapeFact;
    lenPore1 = m_throatData[stdIdx-1]->lenPoreOne;
    lenPore2 = m_throatData[stdIdx-1]->lenPoreTwo;
    lenThroat = m_throatData[stdIdx-1]->lenThroat;
    lenTot = m_throatData[stdIdx-1]->lenTot;

    if(idx > m_origNumThroats)                      // Point to correct pores in the subsequent networks
    {
        if(pore1 != 0 && pore1 != -1)pore1 += numNetsInFront*m_origNumPores;
        if(pore2 != 0 && pore2 != -1)pore2 += numNetsInFront*m_origNumPores;
    }

    double totalLength(0.0);
    pair<MapItr,MapItr> mapPos;
    bool inlet(false), good2go(false);
    if(numNetsInFront != m_numNetInSeries-1 && m_throatData[stdIdx-1]->poreOne == 0)
    {
        good2go = true;
        mapPos = m_outletConnections.equal_range(m_throatData[stdIdx-1]->poreTwo);
    }
    else if(numNetsInFront != m_numNetInSeries-1 && m_throatData[stdIdx-1]->poreTwo == 0)
    {
        good2go = true;
        mapPos = m_outletConnections.equal_range(m_throatData[stdIdx-1]->poreOne);
    }
    else if(numNetsInFront > 0 && m_throatData[stdIdx-1]->poreOne == -1)
    {
        good2go = true;
        inlet = true;
        mapPos = m_inletConnections.equal_range(m_throatData[stdIdx-1]->poreTwo);
    }
    else if(numNetsInFront > 0 && m_throatData[stdIdx-1]->poreTwo == -1)
    {
        good2go = true;
        inlet = true;
        mapPos = m_inletConnections.equal_range(m_throatData[stdIdx-1]->poreOne);
    }

    if(good2go)
    {
        MapItr itr = mapPos.first;
        while(itr != mapPos.second)
        {
            if((*itr).second.first() == stdIdx)
            {
                int poreIdx = (*itr).second.second();

                if(inlet)
                    poreIdx += (numNetsInFront-1)*m_origNumPores;
                else
                    poreIdx += (numNetsInFront+1)*m_origNumPores;

                totalLength = (*itr).second.third();
                if(m_throatData[stdIdx-1]->poreOne == 0 || m_throatData[stdIdx-1]->poreOne == -1)
                {
                    pore1 = poreIdx;
                    lenPore1 = lenPore2;
                }
                else
                {
                    pore2 = poreIdx;
                    lenPore2 = lenPore1;
                }
            }
            ++itr;
        }
        if(m_useAvrXOverThroatLen)
        {
            lenThroat = m_averageThroatLength;
            lenTot = 2.0*lenPore2 + lenThroat;
        }
        else if(totalLength > 2.0*lenPore2)
        {
            lenTot = totalLength;
            lenThroat = totalLength - 2.0*lenPore2;
        }
        else
        {
            lenPore1 = 1.0E-8;
            lenPore2 = 1.0E-8;
            lenThroat = 1.0E-8;
            lenTot = 3.0E-8;
        }
        double xSectArea(pow(rad, 2.0) / (4.0*shapeFact));  // Oren networks have 0 volume for in/outlet throats
        vol = lenThroat*xSectArea;                          // Assign volume based on a tube
    }
}


//void InputData::fracWetting(bool& applyFracWet, double& whatFrac, bool& volBased, double& min,
    //double& max, double& delta, double& eta, string& fracModel, int& clustDiam, bool& oilInWat) const
//{
    //istringstream data;
    //ostringstream optionStr;
    //string keyword("FRAC_CON_ANG");
    //char volModel, oInW('T');
//
    //if(getData(data, keyword))
    //{
        //cout<< "Reading " << keyword << endl;
        //data >> whatFrac >> volModel  >> min >> max >> delta >> eta >> fracModel;
        //min *= acos(-1.0) / 180.0;
        //max *= acos(-1.0) / 180.0;
        //volBased = (volModel == 'T' || volModel == 't');
        //if (!(volModel == 'T' || volModel == 't' || volModel == 'F' || volModel == 'f'))
			//cout<<"      Error Wrong Model, expected T/t or F/f"<<endl;
        //applyFracWet = true;
        //if(!data) errorMsg(keyword);
//
        //data >> clustDiam >> oInW;
        //oilInWat = (oInW == 'T' || oInW == 't');
        //errorInDataCheck(data, keyword);
    //}
    //else
    //{
        //applyFracWet = false;
    //}
    //cout<<"FRAC_CON_ANG read"<<endl;
    //cout<< " " << whatFrac << " " << volModel  << " " << min << " " << max << " " << delta << " " << eta << " " << fracModel<<endl;
    //cout<< " " << clustDiam << " " << oInW<<endl;
//
//
//}

//void InputData::initConAng(double& min, double& max, double& delta, double& eta) const
//{
    //istringstream data;
    //string keyword("INIT_CON_ANG");
//
    //if(getData(data, keyword))
    //{
        //cout<< "Reading " << keyword << endl;
        //data >> min >> max >> delta >> eta;
//
        //min *= acos(-1.0) / 180.0;
        //max *= acos(-1.0) / 180.0;
        //if(!data) errorMsg(keyword);
        //errorInDataCheck(data, keyword);
    //}
    //else
    //{
        //min = 0.0;
        //max = 0.0;
        //delta = 0.2;
        //eta = 3.0;
    //}
    //cout<<"InitCA: min = "<<min<<", max = "<<max<<endl;
//}

void InputData::initConAng(int& wettClass, double& min, double& max, double& delta, double& eta,
                            string& distModel, double& sepAng) const
{
    istringstream data;
    string keywordIni("INIT_CONT_ANG");
    string keywordEqi("EQUIL_CON_ANG");

    if(getData(data, keywordIni))
    {
        if (verbose) cout<< "Reading " << keywordIni << endl;
        data >> wettClass >> min >> max >> delta >> eta;
        if(!data) errorMsg(keywordIni);

        data >> distModel;                  // Maintain backward compatibility
        if(!data) distModel = "rand";       // by allowing both model and
                                            // separation angle not to be present
        data >> sepAng;
        if(!data) sepAng = 25.2;

        errorInDataCheck(data, keywordIni);

        min *= acos(-1.0) / 180.0;
        max *= acos(-1.0) / 180.0;
        sepAng *= acos(-1.0) / 180.0;
    }
    else if(getData(data, keywordEqi))
    {
        cout<< "Reading equilibrium contact angles from keyword " << keywordEqi << endl;
        data >> wettClass >> min >> max >> delta >> eta;
        if(!data) errorMsg(keywordEqi);

        data >> distModel;                  // Maintain backward compatibility
        if(!data) distModel = "rand";       // by allowing both model and
                                            // separation angle not to be present
        data >> sepAng;
        if(!data) sepAng = 25.2;

        errorInDataCheck(data, keywordEqi);

        min *= acos(-1.0) / 180.0;
        max *= acos(-1.0) / 180.0;
        sepAng *= acos(-1.0) / 180.0;
    }
    else missingDataErr(keywordIni);
    //{
        //missingDataErr("neither" + keywordIni+" nor "+keywordEqi);
    //}
}

/*
//void InputData::equilConAng(int& wettClass, double& min, double& max, double& delta, double& eta,
                            //string& distModel, double& sepAng) const
//{
    //istringstream data;
    //string keywordIni("INIT_CONT_ANG");
    //string keywordEqi("EQUIL_CON_ANG");
//
     //if(getData(data, keywordEqi))
    //{
        //cout<< "Reading " << keywordEqi << endl;
        //data >> wettClass >> min >> max >> delta >> eta;
        //if(!data) errorMsg(keywordEqi);
//
        //data >> distModel;                  // Maintain backward compatibility
        //if(!data) distModel = "rand";       // by allowing both model and
                                            //// separation angle not to be present
        //data >> sepAng;
        //if(!data) sepAng = 25.2;
//
        //errorInDataCheck(data, keywordEqi);
//
        //min *= acos(-1.0) / 180.0;
        //max *= acos(-1.0) / 180.0;
        //sepAng *= acos(-1.0) / 180.0;
    //}
    //else if(getData(data, keywordIni))
    //{
        //cout<< "============================================================="<< endl;
        //cout<< "Missing keyword: " << keywordEqi << endl;
        //cout<< "Reading equilibrium contact angles from keyword: " << keywordIni << endl << endl;
        //cout<< "============================================================="<< endl;
        //data >> wettClass >> min >> max >> delta >> eta;
        //if(!data) errorMsg(keywordIni);
//
        //data >> distModel;                  // Maintain backward compatibility
        //if(!data) distModel = "rand";       // by allowing both model and
                                            //// separation angle not to be present
        //data >> sepAng;
        //if(!data) sepAng = 25.2;
//
        //errorInDataCheck(data, keywordIni);
//
        //min *= acos(-1.0) / 180.0;
        //max *= acos(-1.0) / 180.0;
        //sepAng *= acos(-1.0) / 180.0;
    //}
    //else 
    //{
        //missingDataErr("neither" + keywordIni+" nor "+keywordEqi);
    //}
//}



//void InputData::equilConAng(int& wettClass, double& min, double& max, double& delta, double& eta,
                            //string& distModel, double& sepAng) const
//{
    //istringstream data;
    //string keyword("EQUIL_CON_ANG");
//
    //if(getData(data, keyword))
    //{
        //cout<< "Reading " << keyword << endl;
        //data >> wettClass >> min >> max >> delta >> eta;
        //if(!data) errorMsg(keyword);
//
        //data >> distModel;                  // Maintain backward compatibility
        //if(!data) distModel = "rand";       // by allowing both model and
                                            //// separation angle not to be present
        //data >> sepAng;
        //if(!data) sepAng = 25.2;
//
        //errorInDataCheck(data, keyword);
//
        //min *= acos(-1.0) / 180.0;
        //max *= acos(-1.0) / 180.0;
        //sepAng *= acos(-1.0) / 180.0;
    //}
    //else
    //{
        //missingDataErr(keyword);
    //}
//}*/

void InputData::writeNetwork(bool& writeNet, bool& writeBin, string& netName) const
{
    istringstream data;
    string keyword("WRITE_NET");
    char inBinary;

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> inBinary >> netName;
        writeBin = (inBinary == 'T' || inBinary == 't');
        writeNet = true;
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
  }
    else
    {
        writeNet = false;
        writeBin = false;
        netName = "nothing";
    }
}

void InputData::solverDebug(bool& watMat, bool& oilMat, bool& resMat, bool& watVel, bool& oilVel, bool& resVel,
							string& fileName, bool& matlab, bool& initOnly)const
{
    istringstream data;
    string keyword("SOLVER_DBG");
    char writeWat, writeOil, writeRes, writeWatVel, writeOilVel, writeResVel, mat, forInit;

    if(getData(data, keyword))
    {
        if (verbose) cout<< "Reading " << keyword << endl;
        data >> writeWat >> writeOil >> writeRes >> writeWatVel >> writeOilVel >> writeResVel >> fileName >> mat >> forInit;
        watMat = (writeWat == 'T' || writeWat == 't');
        oilMat = (writeOil == 'T' || writeOil == 't');
        resMat = (writeRes == 'T' || writeRes == 't');
        watVel = (writeWatVel == 'T' || writeWatVel == 't');
        oilVel = (writeOilVel == 'T' || writeOilVel == 't');
        resVel = (writeResVel == 'T' || writeResVel == 't');
        matlab = (mat == 'T' || mat == 't');
        initOnly = (forInit == 'T' || forInit == 't');
        if(!data) errorMsg(keyword);
        errorInDataCheck(data, keyword);
    }
    else
    {
        watMat = false;
        oilMat = false;
        resMat = false;
        watVel = false;
        oilVel = false;
        resVel = false;
        matlab = false;
        initOnly = false;
        fileName = "nothing";
    }
}




