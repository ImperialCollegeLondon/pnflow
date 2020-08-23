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


//!  FlowDomain, main network driver,
//!  implements quasi-static network flow simulation routines



#ifndef NETSIM_H
#define NETSIM_H





#include "Element.h"
#include "elem_Model.h"
#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"
#include "compareFuncs.h"
#include "sortedEvents.h"
#include "results3D.h"
#include "FlowData.h"



void  setElemProps(const InputFile& input, const stvec<Element*>& elemans, size_t n6pPors, mstream& out_, const CommonData& comn);

class hypreSolver;
class Fluid;
class CommonData;
class ElemModel;












class FlowDomain
{
public:

	FlowDomain(InputFile& input, unsigned int randSeed);
	~FlowDomain();


	//void initializeDrainage(SortedEvents<Apex*,PceDrainCmp>& evnts, bool wantRelPerm, bool wantResIdx, bool entreL, bool entreR, bool exitL, bool exitR);
	int  dispCycle() const {return comn_.dispCycle();}
	void Drainage(double requestedFinalSw, double requestedFinalPc, double deltaSw, bool calcKr, bool calcI, bool entreL, bool entreR, bool exitL, bool exitR);
	//void finaliseDrainage();

	//void initializeImbibition(SortedEvents<Apex*,PceImbCmp>& evnts, bool wantRelPerm, bool wantResIdx, bool entreL, bool entreR, bool exitL, bool exitR, InputData& input);
	void Imbibition(double requestedFinalSw, double requestedFinalPc, double deltaSw, bool calcKr, bool calcI, bool entreL, bool entreR, bool exitL, bool exitR);
	//void finaliseImbibition();




	InputData&   input()  { return input_;}


private:

	void formatResults(bool, bool, bool);
	void initNetworkOld(InputData& input);
	void modifyNetwork(InputData& input);
	void applyFWettabilityChange(InputData& input);
	void SolveSinglePhaseFlow(InputData& input);
	void writeNetworkToFile(const InputData& input) const;
	void modifyConnNum_removeConnectionsFromNetwork(InputData& input);

	void createMatlabLocationData() const;
	void writeResultData(bool relPermIncluded, bool resIdxIncluded);
	ststr calcUSBMindex() const;
	inline void recordUSBMData(bool isDrainage);
	inline double cpuTimeElapsed(clock_t start);

	void updateSatAndConductances(double Pc);
	void solve_forRelPermResIndex(bool wantRelPerm, bool wantResIdx);
	bool prsOrVoltDrop(const Fluid& fluid, double& prsDrop) const;
	bool avrPrsOrVolt(const Fluid& fluid, int pPlane, double& res, double& stdDev, int& numVal) const;


	template<class compType> double nextCentrInjPc(const SortedEvents<Apex*,compType>&    evnts) const /// slow function
	{
		if (!evnts.empty())	return evnts.peek()->gravCorrectedEntryPress();
		else return  1.0e64*dd_;
	}



	void singleDrainStep(SortedEvents<Apex*,PceDrainCmp>& evnts, double SwTarget, double PcTarget, bool& residualSat);
	void updateCentreInj_Oil(Element* elem, bool& insideBox);
	void update_layerInjs_Oil(Element* elem);

	void popUpdateOilInj(SortedEvents<Apex*,PceDrainCmp>& evnts, bool& insideBox, double & pc, double localPcTarget);
	void popUpdateCentreInj_Oil(Element *currElemCh, SortedEvents<Apex*,PceDrainCmp>& evnts, double & pc);
	void popUpdate_layerInjs_Oil(Apex* apex, SortedEvents<Apex*,PceDrainCmp>& evnts, double & pc);
	void findMarkStoreTrappedWater_reCalcOlPc(SortedEvents<Apex*,PceDrainCmp>& evnts, Element* elem, FluidBlob startPt, double localPc);
	void untrap_OilGanglia(SortedEvents<Apex*,PceDrainCmp>& evnts, Element* elem);
	inline bool insertReCalcDrainEntryPrs(SortedEvents<Apex*,PceDrainCmp>& evnts, Element* elem, double cappPrs);
	template<class compType> inline void reCalcOilLayerPc_markTrappedFilms(SortedEvents<Apex*,compType>& evnts, Element* elem, double cappPrs);

	void checkUntrapWaterIfUnstableConfigsDrain(SortedEvents<Apex*,PceDrainCmp>& evnts);///. rarely does anything
	void addElemTo_layerDrainEvents(SortedEvents<Apex*,PceDrainCmp>& evnts, Element* elem);
	inline void clearTrappedWatFromEvents(SortedEvents<Apex*,PceDrainCmp>& evnts);

	void untrap_WaterGanglia(SortedEvents<Apex*,PceDrainCmp>& evnts, Element* elem, FluidBlob blob);////. Imb+Drain


	void singleImbibeStep(SortedEvents<Apex*,PceImbCmp>& evnts, double SwTarget, double& PcTarget, bool& residualSat);
	void updateCentreInj_Water(Element* elem, bool& insideBox);
	void popUpdateWaterInj(SortedEvents<Apex*,PceImbCmp>& evnts, bool& insideBox, double & pc, double localPcTarget);
	void popUpdateCentreInj_Water(Element *currElemCh, SortedEvents<Apex*,PceImbCmp>& evnts, double & pc);
	void popUpdate_layerInj_Water(Apex* apex, SortedEvents<Apex*,PceImbCmp>& evnts, double & pc);
	void findMarkStoreTrappedOil(SortedEvents<Apex*,PceImbCmp>& evnts, Element* elem, double localPc);
	void untrap_WaterGanglia(SortedEvents<Apex*,PceImbCmp>& evnts, Element* elem, FluidBlob blob);////. Imb+Drain
	inline void insertReCalcImbibeEntryPrs(SortedEvents<Apex*,PceImbCmp>& evnts, Element* elem, double cappPrs);
	void checkUntrapOilIfUnstableConfigsImb(SortedEvents<Apex*,PceImbCmp>& evnts);///. rarely does anything
	void addElemTo_layerImbibeEvents(SortedEvents<Apex*,PceImbCmp>& evnts, Element* elem);
	inline void clearTrappedOilFromEvents(SortedEvents<Apex*,PceImbCmp>&    evnts);


	InputData                              input_;
	CommonData                             comn_;


	const Fluid&                           oil_;// TODELETE
	const Fluid&                           water_;// TODELETE

	const Fluid*                           invadingFluid;//pointers to above variables
	const Fluid*                           retreatingFluid;//pointers to above variables

	int                                    dd_; ///. ? oil inj or wat inj ?

	stvec<Element*>                        elemans_;
	stvec<Element*>                        krInletBoundary_;
	stvec<Element*>                        krOutletBoundary_;
	int                                    sourceNode_;

	stvec< double >                        pressurePlanesLoc_;
	stvec< stvec<Element*> >               pressurePlanes_;


	int                                     nPors_;
	int                                     nBpPors_; ///<  number of boundary and internal nodes TODO: clean up
	int                                     nBSs_; ///<  number of boundary nodes
	int                                     nTrots_;
	int                                     numPressurePlanes_;


	bool                                    useAvrPrsAsBdr_;
	bool                                    wantRelPerm_;
	bool                                    wantResIdx_;
	bool                                    includeGravityInRelPerm_;
	double                                  inletSolverPrs_;
	double                                  outletSolverPrs_;
	double                                  deltaPo_;
	double                                  deltaPw_;
	double                                  deltaV_;

	double                                  singlePhaseWaterQ_;
	double                                  singlePhaseOilQ_;
	double                                  singlePhaseDprs_;
	double                                  singlePhaseCurrent_;
	double                                  singlePhaseDvolt_;

	double                                  oilFlowRate_;
	double                                  watFlowRate_;
	double                                  current_;

	mutable double                          maxOilFlowErr_;
	mutable double                          maxWatFlowErr_;
	mutable double                          maxResIdxErr_;



	double                                  solverBoxStart_;
	double                                  solverBoxEnd_;
	double                                  flowVolume_; //!< void volume of elements inside CALC_BOX
	double                                  clayVolume_; //!< clay volume of elements inside CALC_BOX
	double                                  satBoxVolume_;
	dbl3                                    box_, X0_;

	double coreFrac_;  //!< cross-sectional area fraction, excluding core sample sleeve etc

	double                                  Sw_;
	double                                  cappPressCh_;
	const double &                          Pc_;

	double                                  maxCyclePc_;
	double                                  minCyclePc_;
	int                                     totNumFillings_;



	double                                  initStepSize_;
	double                                  extrapCutBack_;
	double                                  maxFillIncrease_;
	bool                                    stableFilling_;

	TrappingCriteria                        trappingCriteria_;
	double                                  KcIncr_; // this is PcIncr_
	double                                  KcIncrFactor_; // this is PcIncrFactor_
	int                                     minNumFillings_;	//double                                  minVolFillings_;
	double      minPcStep_; // to sort out this


	hypreSolver*                            solver_;


	ststr                                   title_;

	mstream                                 out_;
	std::ofstream                           drainListOut_;
	std::ofstream                           imbListOut_;
	//dbgstream                               dbgFile;


	bool                                    reportMaterialBal_;
	bool                                    prtPressureProfile_;
	bool                                    writeWatMatrix_;
	bool                                    writeOilMatrix_;
	bool                                    writeResMatrix_;
	bool                                    writeWatVelocity_;
	bool                                    writeOilVelocity_;
	bool                                    writeResVelocity_;
	bool                                    writeSlvMatrixAsMatlab_;
	bool                                    writeDrainList_;
	bool                                    writeImbList_;



	stvec< std::pair< double, double > >    usbmDataDrainage_;
	stvec< std::pair< double, double > >    usbmDataImbibition_;
	stvec< double >                         amottDataDrainage_;
	stvec< double >                         amottDataImbibition_;
	double                                  amottOilIdx_;
	double                                  amottWaterIdx_;
	stvec< std::pair<double, double> >      resultWaterFlowRate_;
	stvec< std::pair<double, double> >      resultOilFlowRate_;
	stvec< double >                         resultWaterSat_;
	stvec< double >                         resultCappPress_;
	stvec< double >                         resultResistivityIdx_;
	stvec< double >                         resultWaterMass_;
	stvec< double >                         resultOilMass_;

	double                                  cpuTimeTotal_;
	double                                  cpuTimeCoal_;
	double                                  cpuTimeKrw_;
	double                                  cpuTimeKro_;
	double                                  cpuTimeResIdx_;
	double                                  cpuTimeTrapping_;




 public:
	results3D                               results3D_; //: order is important
};



inline double FlowDomain::cpuTimeElapsed(clock_t start)
{
	return (static_cast< double >(clock() - start) / CLOCKS_PER_SEC);
}


inline void FlowDomain::recordUSBMData(bool isDrainage)
{
	std::pair<double, double> entry(Pc_, Sw_);
	if(isDrainage && Pc_ > 0.0)        usbmDataDrainage_.push_back(entry);
	else if(!isDrainage && Pc_ < 0.0)  usbmDataImbibition_.push_back(entry);
}




#endif
