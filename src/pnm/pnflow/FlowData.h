#ifndef FLOWDATA_H
#define FLOWDATA_H


#include "InputFile.h"
#include "fluid.h"
#include <algorithm>
#include <random>
#include "inputData.h"
#include "typses.h"



class Element;
class InputData;
enum TrappingCriteria {escapeToInlet = 0, escapeToOutlet, escapeToEither, escapeToBoth};
enum FluidBlob {filmBlob = 0, bulkBlob};



class RockType
{
 public:
	RockType(const InputFile& input, int index) : index_(index)
	{
		//#ifdef EXPERIMENTAL
			//maxRcRandnessF_ = input.getOr(0.0, "maxPcRandness");
			//maxRcRandness_  = input.getOr(0.0, "maxRcRandness");
		//#else
			maxRcRandnessF_ = input.getOr(0.01, "maxPcRandness");
			maxRcRandness_ = input.getOr(0.50, "maxRcRandness"); //between 0 and 1
		//#endif
	}
	virtual double Sw(double Pc, int icycle) const = 0;
	virtual double Pc(double Sw, int icycle) const = 0;
	virtual ~RockType() {};
	int index() const { return index_;}
//protected :
	const int index_;
	double maxRcRandnessF_;
	double maxRcRandness_;

};

class SolidWall : public RockType
{
 public:
	SolidWall(const InputFile& input, int index) : RockType(input, index) {}
	double Sw(double Pc, int icycle) const {return 1.0;};
	double Pc(double Sw, int icycle) const {return 1.0e29*(2*(icycle%2)-1.0);};
};


//template <typename T>
struct WeibulParam
{
	WeibulParam() : minV(0.0), maxV(0.0), deltaExp(-0.2), etaExp(-3.0), correlation("rand") {}
	WeibulParam(std::istream& ins) : minV(0.0), maxV(0.0), deltaExp(-0.2), etaExp(-3.0), correlation("rand") { read(ins); }

	inline void read(std::istream& ins)
	{
		ins>>minV;
		if(ins.good()) ins >>maxV;    else maxV=minV;

		if(ins.good()) ins >>deltaExp;
		if(ins.good()) ins >>etaExp;
		if(ins.good()) ins >>correlation;
		std::cout<<"WeibulParam: "<<minV<<"  "<<maxV<<"  "<<deltaExp<<"  "<<etaExp<<"  "<<correlation<<"  "<<std::endl;
	}
 
	double minV, maxV;
	double deltaExp, etaExp;
	std::string correlation;
};


class CommonData
{
 public:

	CommonData(const InputData& input, unsigned int randSeed)
	:	input_(input), 
		watoil_("WaterOilInterface",input,1),
		water_(WTR,input, "Water"),
		oil_(OIL,input, "Oil"), 
		elec_(ELEC,input, "Electricity"),
		clay_(CLAY,input, "Clay"), 
			maxEverCappPress_(0.0),
			minEverCappPress_(0.0),   
		   trappedRegionsOil_(trappedRegionsOil__),
		   trappedRegionsWat_(trappedRegionsWat__),
		   randomGenerator_(randSeed), uniformRand01_(0.0,1.0),
		   informative(input.informative)
	{
		std::cout<<"creating flowData"<<std::endl;
		rockTypes_.push_back( new SolidWall(input,0) );
		floodingCycle_ = 0;
		circWatCondMultFact_ = input.getOr(0.0,"SURFACE_FILM_COND_FACT");
		input.lookupOr("GRAV_CONST",dbl3(0.,0.,-9.81));
		input.poreFillWgt(poreFillWeights_);
		input.poreFillAlg(poreFillAlg_);

		maxPcLastDrainCycle_ = 0.0;
		minPcLastImbCycle_ = 0.0;
		numSquares_ = 0;
		numTriangles_ = 0;
		numCircles_ = 0;
		injectant_ = 0;
		guessPc_ = 0;
		drainagePhase_ = false;

		//double interfacTen, watDens(water_.density()), oilDens(oil_.density()) , watVisc, oilVisc, watResist, oilResist;
		//input.fluid(interfacTen, watVisc, oilVisc, watResist, oilResist, watDens, oilDens);
		//double clayResistivity(watResist), surfaceResistivity(1.0e12);
		//input.getVar(clayResistivity, "clayResistivity");
		//input.getVar(surfaceResistivity, "surfaceResistivity");

		//double epsilonW(80.1*8.854e-12), zetaW(0.061), zetaW_OW(0.061), epsilonO(80.1*8.854e-12*0.0001), zetaOW(-0.061), zetaO(0.0);


		//water_.setFluidProps(watVisc, interfacTen, watResist, watDens, epsilonW, zetaW);
		//clay_.setFluidProps(watVisc*1.0e12, 0, clayResistivity, watDens, epsilonW, zetaW);
		//rockWaterSurface_.setFluidProps(watVisc*1.0e12, 0, surfaceResistivity, watDens, epsilonW, zetaW);
		//rockWaterSurface__OW.setFluidProps(watVisc*1.0e12, 0, surfaceResistivity, watDens, epsilonW, zetaW_OW);
		//watoil_.setFluidProps(watVisc*1.0e12, interfacTen, surfaceResistivity, watDens, epsilonW, zetaOW); ///. only zetaOW is used
		//oil_.setFluidProps(oilVisc, interfacTen, oilResist, oilDens, epsilonO, zetaO); ///. zetaO is not used
		rockWaterSurface_.resistivity_=1.0e300;

		std::cout<<"WO.IFT:"<<watoil_.sigma()<<std::endl;
		std::cout<<"W.viscosity:"<<water_.viscosity()<<std::endl;
		std::cout<<"W.resistivity:"<<water_.resistivity()<<std::endl;
		std::cout<<"W.density:"<<water_.density()<<std::endl;
		std::cout<<"O.viscosity:"<<oil_.viscosity()<<std::endl;
		std::cout<<"O.resistivity:"<<oil_.resistivity()<<std::endl;
		std::cout<<"O.density:"<<oil_.density()<<std::endl; 

	}


	~CommonData(){}


	double poreFillWeights(int i) const {return poreFillWeights_[i];}
	const std::string& poreFillAlg() const {return poreFillAlg_;}



	int dispCycle() const {return floodingCycle_;}
	void incrementFloodCycle() {++floodingCycle_;}
	bool isDrainage() const {return floodingCycle_%2==1;}


	double maxEverPc() const {return maxEverCappPress_;}
	double maxPcLastDrainCycle() const {return maxPcLastDrainCycle_;}
	inline void setMaxPcLastDrainCycle(double pc);
	double minEverCappPress() const {return minEverCappPress_;}
	double minPcLastImbCycle() const {return minPcLastImbCycle_;}
	inline void setMinPcLastCycle(double pc);
	void GuessCappPress(double pc) { guessPc_ = pc; }
	double GuessCappPress() const  { return guessPc_; }
	const dbl3& gravConst() const      { return gravConst_; }

	double KrwatcornAtSw0() const {return circWatCondMultFact_;}

	void injectant(const Fluid* injFluid) {injectant_ = injFluid;}
	const Fluid* injectant() const {return injectant_;}
	const Fluid& noFluid() const {return nof_;}
	const Fluid& oil() const {return oil_;}
	const Fluid& water() const { return water_; }
	const Fluid& elec() const { return elec_; }
	const Fluid& clay() const { return clay_; }
	const Fluid& surface() const { return rockWaterSurface_; }
	const Fluid& surfaceOW() const { return rockWaterSurface__OW; }
	const Interface& oilWater() const { return watoil_; }
	double sigmaOW() const { return watoil_.sigma(); }

	size_t newOilTrappingIndex() const {return trappedRegionsOil_.size();}
	size_t newWatTrappingIndex() const {return trappedRegionsWat_.size();}

	inline void removeTrappedOilElem(int idx, Element* elem);
	inline void removeTrappeWatElem(int idx, std::pair<Element*, FluidBlob> elem);
	void addTrappeWatElem(int idx, std::pair<Element*, FluidBlob> elem) {    trappedRegionsWat__[idx].push_back(elem); }

	void addTrappedRegionOil(const std::vector<Element*>& elems) {trappedRegionsOil__.push_back(elems);}
	void addTrappedRegionWat(const std::vector< std::pair<Element*, FluidBlob> >& elems) {trappedRegionsWat__.push_back(elems);}
	void removeTrappedRegionOil(int region) {trappedRegionsOil__[region].clear();}
	void removeTrappedRegionWat(int region) {trappedRegionsWat__[region].clear();}

	const std::vector<Element*>& trappedRegionsOil(int region) const {return trappedRegionsOil_[region];}
	const std::vector< std::pair<Element*,FluidBlob> >& trappedRegionsWat(int region) const {return trappedRegionsWat_[region];}

	/// for writing statistics
	const std::vector< std::vector<Element*> >	&						trappedOilRegions() const {return trappedRegionsOil_;}
	const std::vector< std::vector< std::pair<Element*,FluidBlob> > > &    trappedWatRegions() const {return trappedRegionsWat_;}



	int numPores() const {return nPors_;}
	int numThroats() const {return nTrots_;}
	void setNumElem(int numP, int numT) {nPors_ = numP; nTrots_ = numT;}

	size_t numTrappedOilRegions() const;
	size_t numTrappedOilElems() const;
	size_t numTrappedWatRegions() const;
	size_t numTrappedWatElems() const;    

	void countCircle()   const {++numCircles_;}
	void countTriangle() const {++numTriangles_;}
	void countSquare()   const {++numSquares_;}
	void countPorous()   const {++numPorous_;}
	int numCircles() const {return numCircles_;}
	int numTriangles() const {return numTriangles_;}
	int numSquares() const {return numSquares_;}
	int numPorous() const {return numPorous_;}    

	std::ofstream                        dbgOut;

	std::default_random_engine&  randomGenerator() const  { return randomGenerator_; }
	double  rand01() const   { return uniformRand01_(randomGenerator_); }
	inline double weibull(double minv, double maxv, double delta, double eta) const // delta = pow(beta,eta)
	{	if(delta < 0.0 || eta < 0.0)    return minv + (maxv-minv)*rand01(); // Uniform Distribution
		else  return minv+(maxv-minv) * pow(-delta*log(1.0-rand01()*(1.0-exp(-1.0/delta))), 1.0/eta);  /// Weibull truncated up to 1.0, scaled between min and max
		  /// wiki:  eta -> k ,  delta -> lambda^k
	}
	const InputFile&	input() const { return input_; }
	const stvec<RockType*>& rockTypes() const {return rockTypes_;}
	stvec<RockType*>&	rockTypesCh() { return rockTypes_; }

	void addKcSwKrsQs(const stvec<double>& KcSwKrsQs)
	{
		int icy = dispCycle();
		while(int(KcSwKrsQsss_.size())<=icy) KcSwKrsQsss_.push_back(stvec<stvec<double> >());
			if(KcSwKrsQsss_[icy].size()) {
				if(isDrainage())
				{
					ensure(KcSwKrsQsss_[icy].back()[2]+1.0e-8>=KcSwKrsQs[2]);
					ensure(KcSwKrsQsss_[icy].back()[3]<=KcSwKrsQs[3]+1.0e-8);
				}
				else
				{
					ensure(KcSwKrsQsss_[icy].back()[2]<=KcSwKrsQs[2]+1.0e-8);
					ensure(KcSwKrsQsss_[icy].back()[3]+1.0e-8>=KcSwKrsQs[3]);
				}
			}

		KcSwKrsQsss_[icy].push_back(KcSwKrsQs);
	}

	void writeResultData(bool wantRelPerm, bool wantResIdx, const std::string& title) const;

	stvec<stvec<stvec<double> > >   KcSwKrsQsss_; //  used in svg only 

 private:
	const InputFile&               input_;


	Interface                                   watoil_;
	Fluid                                       nof_; ///. TODO make const
	Fluid                                       water_;
	Fluid                                       oil_; ///. TODO make const
	Fluid                                       elec_;
	Fluid                                       clay_;
	Fluid                                       rockWaterSurface_;   
	Fluid                                       rockWaterSurface__OW;

	const Fluid*                                  injectant_;

	double										  maxEverCappPress_;
	double                                        maxPcLastDrainCycle_;
	double										  minEverCappPress_;
	double                                        minPcLastImbCycle_;
	double                                        guessPc_;
	double										  circWatCondMultFact_;
	dbl3                                        gravConst_;
	int													nPors_;
	int													nTrots_;
	bool                                                drainagePhase_;

	int                                                 floodingCycle_;

	stvec<RockType*>               rockTypes_;

	std::vector< double >									poreFillWeights_;
	std::string												poreFillAlg_;
	std::vector< std::vector<Element*> >							trappedRegionsOil__;
	std::vector< std::vector< std::pair<Element*,FluidBlob> > >		trappedRegionsWat__;
	const std::vector< std::vector<Element*> >	&						trappedRegionsOil_;
	const std::vector< std::vector< std::pair<Element*,FluidBlob> > >	&    trappedRegionsWat_;


	mutable std::default_random_engine              randomGenerator_;
	mutable std::uniform_real_distribution<double>  uniformRand01_;

	///.  not used data
	mutable int										numSquares_;
	mutable int										numTriangles_;
	mutable int										numCircles_;
	mutable int										numPorous_;

 public:
	int										informative;
};

inline size_t CommonData::numTrappedOilElems() const
{
	size_t num(0);
	for(const auto& reg:trappedRegionsOil_)  num += reg.size();
	return num;
}

inline size_t CommonData::numTrappedWatElems() const
{
	size_t num(0);
	for(const auto& reg:trappedRegionsWat_)  num += reg.size();
	return num;
}

inline size_t CommonData::numTrappedOilRegions() const
{
	size_t num(0);
	for(const auto& reg:trappedRegionsOil_)  num += !reg.empty();
	return num;
}

inline size_t CommonData::numTrappedWatRegions() const
{
	size_t num(0);
	for(const auto& reg:trappedRegionsWat_)  num += !reg.empty();
	return num;
}

inline void CommonData::setMaxPcLastDrainCycle(double pc)
{
	maxEverCappPress_ = std::max(pc, maxEverCappPress_);
	maxPcLastDrainCycle_ = pc;
}

inline void CommonData::setMinPcLastCycle(double pc)
{
	minEverCappPress_ = std::min ( pc, minEverCappPress_);
	minPcLastImbCycle_ = pc;
}

///. deletes elem from trappedRegionsOil //called from checkUntrapOilIfUnstableConfigsImb
inline void CommonData::removeTrappedOilElem(int idx, Element* elem)
{
	trappedRegionsOil__[idx].erase(
		std::remove(trappedRegionsOil__[idx].begin(),trappedRegionsOil__[idx].end(),elem),
		trappedRegionsOil__[idx].end());
}

///. deletes elem blob from trappedRegionsWat //called from checkUntrapWaterIfUnstableConfigsDrain
inline void CommonData::removeTrappeWatElem(int idx, std::pair<Element*, FluidBlob> elem)
{
	trappedRegionsWat__[idx].erase(
		std::remove(trappedRegionsWat__[idx].begin(),trappedRegionsWat__[idx].end(),elem),
		trappedRegionsWat__[idx].end());
}




#endif  //FLOWDATA_H

