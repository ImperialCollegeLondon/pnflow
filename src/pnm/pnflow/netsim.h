/*!

Copyright (c):  2004, Imperial College Pore-Scale Consortium
                2004, Per Valvatne and Martin J Blunt,  original code
                2015, Ali Q Raeini: clean-up and post-processing


This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgement in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be marked as such, and must not be
   misrepresented as being the original software.
3. This copyright notice may not be altered or removed from any source
   distribution.


--------------------------------------------------------------------
See the Imperial College Consortium on Pore-Scale Modelling and Imaging
for contact details:
https://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling/

*/


//!  Netsim, main network driver,
//!  implements quasi-static network flow simulation routines



#ifndef NETSIM_H
#define NETSIM_H

#ifdef WIN32
#pragma warning(disable:4786)
#endif




class Solver;
class Fluid;
class CommonData;
class ElemModel;

#include "Element.h"
#include "elem_Model.h"
#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"
#include "compareFuncs.h"
#include "sortedEvents.h"
#include "vtuWriter.h"
#include "utility.h"
#include "netsim_data.h"


///. class to write both to std::cout, and .prt file
class mstream
{
	public:
	ofstream prtFile;
	mstream(string fileName){ if (!fileName.empty()) prtFile.open(fileName.c_str()); };
	~mstream(void){};
	mstream& operator<<(ostream& (*pfun)(ostream&)) {pfun(prtFile); pfun(std::cout); return *this;}
	ofstream& fileStream() {return prtFile;};
};

template <class T>
mstream& operator<< (mstream& st, T val)
{
  (st.prtFile) << val;
  std::cout<< val;
  return st;
}

///. class to write .dbg  file for debugging
//class dbgstream
//{
	//public:
	//ofstream prtFile;
	//dbgstream(string fileName,int dbgMode) : debugMode(dbgMode) { if (dbgMode && !fileName.empty()) prtFile.open(fileName.c_str()); };
	//~dbgstream(void){};
	//dbgstream& operator<<(ostream& (*pfun)(ostream&)) {pfun(prtFile); return *this;}
	//void flush() {prtFile.flush();};
	//int debugMode;
//};

//template <class T>
//dbgstream& operator<< (dbgstream& st, T val)
//{
	//if (st.debugMode)
		//(st.prtFile) << val;
  //return st;
//}
 //#define dbgstream OnDemandStream
 //#define dbgFile OnDemandStream::dbgFile











#define m_layerEventsCh m_eventsCh


class Netsim
{
public:

    Netsim(InputData & input);
    //Netsim(const Netsim& net);
    ~Netsim();

    void init(InputData& input);

    void initializeDrainage(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, bool wantRelPerm, bool wantResIdx, bool entreL, bool entreR, bool exitL, bool exitR);
    int  floodingCycle() const {return m_comn.floodingCycle();}
	void Drainage(InputData& input, double& Sw, double& Pc, double requestedFinalSw, double requestedFinalPc,
              double deltaSw, double deltaPc, double deltaPcIncFactor, bool calcKr, bool calcI, bool entreL, bool entreR, bool exitL, bool exitR, bool swOut);
    void finaliseDrainage();
   
    bool oilInjection() const {return m_comn.injectant() == &m_oil;}
    //int floodingCycle() const {return m_comn.floodingCycle();}
    void initializeImbibition(SortedEvents<Apex*,PceImbCmp>& m_eventsCh, bool wantRelPerm, bool wantResIdx, bool entreL, bool entreR, bool exitL, bool exitR, InputData& input);
    void Imbibition(InputData& input, double& sw, double& pc, double requestedFinalSw, double requestedFinalPc,
              double deltaSw, double deltaPc, double deltaPcIncFactor, bool calcKr, bool calcI, bool entreL, bool entreR, bool exitL, bool exitR, bool swOut);
    void finaliseImbibition();


    const Oil& oil() const {return  m_oil;}
    const Water& water() const {return  m_water;}


	int  dispCycle() const {return m_comn.floodingCycle();}
	bool isDrainage() const {return m_comn.isDrainageCycle();}


private:

    void formatResults(bool, bool, bool);
    void initNetwork(InputData& input);
    int setupInsidePoreHashings(InputData& input, vector< pair< int, double > >& insidePoreHashs) const;
    void modifyNetwork(InputData& input);
    void applyFWettabilityChange(InputData& input);
    void SolveSinglePhaseFlow(InputData& input);
    int readAndCreateThroats(InputData& input, vector< pair<int, int> >& pores,
        vector<Element*>& inT, vector<Element*>& outT, const vector< pair< int, double > >& insidePoreHashs,
        vector< int >& insideThroatHashs, int newNumPores);
    void writeNetworkToFile(const InputData& input) const;
    void readAndCreatePores(InputData& input, /*vector< pair<int, int> >& connectingPores,
        const vector< pair< int, double > >& insidePoreHashs,*/ const vector< int >& insideThroatHashs, int newNumPores);
    void createInAndOutletPore(int index, double xPos, double yPos, double zPos, vector<Element*>& connThroats);
    int reIndex(int index) const;
    void modifyConnNum_removeConnectionsFromNetwork(InputData& input);

    void createMatlabLocationData() const;
    //void recordPrsProfiles(const Fluid* fluid);
    void recordRes(bool wantRelPerm, bool wantResIdx);
    //void writePrsProfileData(const Fluid* fluid, ostream& resStream);
    void writeResultData(bool relPermIncluded, bool resIdxIncluded);
    string calcUSBMindex() const;
    inline void recordUSBMData(bool isDrainage);
    inline void writePrtData(ostringstream& out);
    inline double cpuTimeElapsed(clock_t start);

    void updateSatAndConductances(double Pc);
    void solve_forRelPermResIndex(bool wantRelPerm, bool wantResIdx);
	void solveforRelPerm();
	void solveforResIndex();
	bool prsOrVoltDrop(const Fluid* fluid, int resistSolve, double& prsDrop) const;
    bool avrPrsOrVolt(const Fluid* fluid, int resist, int pPlane, double& res, double& stdDev, int& numVal) const;

	//template<class compType> double nextCentrInjPc(const SortedEvents<Apex*,compType>&    m_eventsCh) const /// slow function
	//{	
		//double nextEventPc = 1.0e64*m_dd;
		//if (!m_eventsCh.empty())
			//nextEventPc = m_eventsCh.peek()->gravCorrectedEntryPress();
		//if(!m_layerEventsCh.empty())
		//{
			//double nextLayerPc = m_layerEventsCh.peek()->gravCorrectedEntryPress() ;
			//nextEventPc =
			//m_dd*min( m_dd*nextEventPc,
				 //m_dd*nextLayerPc
			   //);
		  //}
		//return nextEventPc;
	//}

	template<class compType> double nextCentrInjPc(const SortedEvents<Apex*,compType>&    m_eventsCh) const /// slow function
	{	
		if (!m_eventsCh.empty())	return m_eventsCh.peek()->gravCorrectedEntryPress();
		else return  1.0e64*m_dd;
	}



    void singleDrainStep(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, double SwTarget, double PcTarget, bool& residualSat);
    void updateCentreInj_Oil(Element* elem, bool& insideBox);
    void update_layerInjs_Oil(Element* elem);

    void popUpdateOilInj(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, bool& insideBox, double & pc, double localPcTarget);
    void popUpdateCentreInj_Oil(Element *currElemCh, SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, double & pc);
    void popUpdate_layerInjs_Oil(Apex* apex, SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, double & pc);
    void findMarkStoreTrappedWater_reCalcOlPc(SortedEvents<Apex*,PceDrainCmp>& m_layerEventsCh, Element* elem, FluidBlob startPt, double localPc);
    void untrap_OilGanglia(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, Element* elem);
    inline bool insertReCalcDrainEntryPrs(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, Element* elem, double cappPrs);
	template<class compType> inline void reCalcOilLayerPc_markTrappedFilms(SortedEvents<Apex*,compType>& m_layerEventsCh, Element* elem, double cappPrs);
	
    void checkUntrapWaterIfUnstableConfigsDrain(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh);///. rarely does anything
    void addElemTo_layerDrainEvents(SortedEvents<Apex*,PceDrainCmp>& m_layerEventsCh, Element* elem);
    inline void clearTrappedWatFromEvents(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh);   

    void untrap_WaterGanglia(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, Element* elem, FluidBlob blob);////. Imb+Drain


    void singleImbibeStep(SortedEvents<Apex*,PceImbCmp>& m_eventsCh, double SwTarget, double& PcTarget, bool& residualSat);
    void updateCentreInj_Water(Element* elem, bool& insideBox);
    void popUpdateWaterInj(SortedEvents<Apex*,PceImbCmp>& m_eventsCh, bool& insideBox, double & pc, double localPcTarget);
    void popUpdateCentreInj_Water(Element *currElemCh, SortedEvents<Apex*,PceImbCmp>& m_eventsCh, double & pc);
    void popUpdate_layerInj_Water(Apex* apex, SortedEvents<Apex*,PceImbCmp>& m_eventsCh, double & pc);
    void findMarkStoreTrappedOil(SortedEvents<Apex*,PceImbCmp>& m_layerEventsCh, Element* elem, double localPc);
    void untrap_WaterGanglia(SortedEvents<Apex*,PceImbCmp>& m_eventsCh, Element* elem, FluidBlob blob);////. Imb+Drain
    inline void insertReCalcImbibeEntryPrs(SortedEvents<Apex*,PceImbCmp>& m_eventsCh, Element* elem, double cappPrs);
    void checkUntrapOilIfUnstableConfigsImb(SortedEvents<Apex*,PceImbCmp>& m_eventsCh);///. rarely does anything
    void addElemTo_layerImbibeEvents(SortedEvents<Apex*,PceImbCmp>& m_layerEventsCh, Element* elem);
    inline void clearTrappedOilFromEvents(SortedEvents<Apex*,PceImbCmp>&    m_eventsCh);


    static const double                             MAX_FLOW_ERR;
    static const int                                DUMMY_IDX;

    const InputData&                                m_input;
    CommonData                                      m_comn;
   



    const Oil&                              m_oil;// TODELETE
    const Water&                              m_water;// TODELETE

    const Fluid*                                    invadingFluid;//pointers to above variables
    const Fluid*                                    retreatingFluid;//pointers to above variables
	
	int 											m_dd; ///. ? oil inj or wat inj ?

    vector<Element*>                              m_rockLattice;
    vector<Element*>                              m_krInletBoundary;
    vector<Element*>                              m_krOutletBoundary;
    int                                             m_sourceNode;
    vector<const Element*>                        m_elemans;  // duplicate, for merging

    vector< double >                                m_pressurePlanesLoc;
    vector< vector<Element*> >                    m_pressurePlanes;


    int                                             m_maxNonZeros;
    int                                             m_numPores;
    int                                             m_numThroats;
    int                                             m_numIsolatedElems;
    int                                             m_numPressurePlanes;
    int                                             m_minNumFillings;
    int                                             m_totNumFillings;


    bool                                            m_useAvrPrsAsBdr;
    bool                                            m_wantRelPerm;
    bool                                            m_wantResIdx;   
    bool                                            m_includeGravityInRelPerm;
    double									        m_inletSolverPrs;
    double									        m_outletSolverPrs;
    double                                          m_deltaPo;
    double                                          m_deltaPw;
    double                                          m_deltaV;
   
    double                                          m_singlePhaseWaterQ;
    double                                          m_singlePhaseOilQ;
    double                                          m_singlePhaseDprs;
    double                                          m_singlePhaseCurrent;
    double                                          m_singlePhaseDvolt;
   
    double                                          m_oilFlowRate;
    double                                          m_watFlowRate;
	double                                          m_current;

    mutable double                                  m_maxOilFlowErr;
    mutable double                                  m_maxWatFlowErr;
    mutable double                                  m_maxResIdxErr;
   
   
    double                                          m_keepFraction;
    double                                          m_satBoxStart;
    double                                          m_satBoxEnd;
    double                                          m_solverBoxStart;
    double                                          m_solverBoxEnd;
    double                                          m_totalVoidVolume;
    double                                          m_totalFlowVolume;
    double                                          m_totalClayVolume;
    double                                          m_satBoxVolume;
    double                                          m_xSize;
    double                                          m_ySize;
    double                                          m_zSize;
   

    double                                          m_satWater;
    double                                          m_cappPressCh;
    const double &                                  m_cappPress;

    double                                          m_maxCycleCappPress;
    double                                          m_minCycleCappPress;



    double                                          m_initStepSize;
    double                                          m_extrapCutBack;
    double                                          m_maxFillIncrease;
    bool                                            m_StableFilling;   
    bool                                            m_injAtLeftRes;
    bool                                            m_injAtRightRes;
    TrappingCriteria                                m_trappingCriteria;

   
    Solver*                                         m_solver;


    string  										m_baseFileName;
    string                                          m_matrixFileName;

   
    mstream 										m_out;
    ofstream                                        m_drainListOut;
    ofstream                                        m_imbListOut;
    //dbgstream                                        dbgFile;   


    bool                                            m_apexPrsReported;
    bool                                            m_reportMaterialBal;
    bool                                            m_prtPressureProfile;
    bool                                            m_writeWatMatrix;
    bool                                            m_writeOilMatrix;
    bool                                            m_writeResMatrix;
    bool                                            m_writeWatVelocity;
    bool                                            m_writeOilVelocity;
    bool                                            m_writeResVelocity;
    bool                                            m_writeSlvMatrixAsMatlab;
    bool                                            m_writeDrainList;
    bool                                            m_writeImbList;
   
   
    vector< string >                                m_results;
    //vector< vector< string > >                      m_watPrsProfiles;
    //vector< vector< string > >                      m_oilPrsProfiles;
    vector< pair< double, double > >                m_usbmDataDrainage;
    vector< pair< double, double > >                m_usbmDataImbibition;
    vector< double >                                m_amottDataDrainage;
    vector< double >                                m_amottDataImbibition;
    double                                          m_amottOilIdx;
    double                                          m_amottWaterIdx;
    vector< pair<double, double> >                  m_resultWaterFlowRate;
    vector< pair<double, double> >                  m_resultOilFlowRate;
    vector< double >                                m_resultWaterSat;
    vector< double >                                m_resultCappPress;
    vector< double >                                m_resultResistivityIdx;
    vector< double >                                m_resultWaterMass;
    vector< double >                                m_resultOilMass;

   
    double                                          m_cpuTimeTotal;
    double                                          m_cpuTimeCoal;
    double                                          m_cpuTimeKrw;
    double                                          m_cpuTimeKro;
	double                                          m_cpuTimeResIdx;
    double                                          m_cpuTimeTrapping;

    vtuWriter                                       m_vtkWriter;
    bool 											useHypre;				
};


inline void Netsim::writePrtData(ostringstream& out)
{
    m_out << out.str() << endl;
}

inline double Netsim::cpuTimeElapsed(clock_t start)
{
    return (static_cast< double >(clock() - start) / CLOCKS_PER_SEC);
}


inline void Netsim::recordUSBMData(bool isDrainage)
{
    pair<double, double> entry(m_cappPress, m_satWater);
    if(isDrainage && m_cappPress > 0.0)
        m_usbmDataDrainage.push_back(entry);
    else if(!isDrainage && m_cappPress < 0.0)
        m_usbmDataImbibition.push_back(entry);
}




#endif
