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
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

*/


//!  FlowDomain, main network driver,
//!  implements quasi-static network flow simulation routines



#ifndef CNMDATA_H
#define CNMDATA_H



#include "inputData.h"

#include "FlowData.h"







class CommonData : public GNMData
{
 public:

	CommonData(const InputData& input)
	:	GNMData(input,""),
			maxEverCappPress_(0.),
			minEverCappPress_(0.),   
		   trappedRegionsOil_(trappedRegionsOil__),
		   trappedRegionsWat_(trappedRegionsWat__),
		   informative(input.informative)
	{ 
		circWatCondMultFact_ = input.getOr("SURFACE_FILM_COND_FACT", 0.);
		input.getOr("GRAV_CONST",dbl3(0.,0.,-9.81));
		input.poreFillWgt(poreFillWeights_);
		input.poreFillAlg(poreFillAlg_);

		maxPcLastDrainCycle_ = 0.;
		minPcLastImbCycle_ = 0.;
		numSquares_ = 0;
		numTriangles_ = 0;
		numCircles_ = 0;
		injectant_ = 0;
		guessPc_ = 0;
		drainagePhase_ = false;

		//double interfacTen, watDens(water_.density()), oilDens(oil_.density()) , watVisc, oilVisc, watResist, oilResist;
		//input.fluid(interfacTen, watVisc, oilVisc, watResist, oilResist, watDens, oilDens);
		//double clayResistivity(watResist), surfaceResistivity(1e12);
		//input.giv("clayResistivity", clayResistivity);
		//input.giv("surfaceResistivity", surfaceResistivity);

		//double epsilonW(80.1*8.854e-12), zetaW(0.061), zetaW_OW(0.061), epsilonO(80.1*8.854e-12*0.0001), zetaOW(-0.061), zetaO(0.);


		//water_.setFluidProps(watVisc, interfacTen, watResist, watDens, epsilonW, zetaW);
		//clay_.setFluidProps(watVisc*1e12, 0, clayResistivity, watDens, epsilonW, zetaW);
		//rockWaterSurface_.setFluidProps(watVisc*1e12, 0, surfaceResistivity, watDens, epsilonW, zetaW);
		//rockWaterSurface__OW.setFluidProps(watVisc*1e12, 0, surfaceResistivity, watDens, epsilonW, zetaW_OW);
		//watoil_.setFluidProps(watVisc*1e12, interfacTen, surfaceResistivity, watDens, epsilonW, zetaOW); ///. only zetaOW is used
		//oil_.setFluidProps(oilVisc, interfacTen, oilResist, oilDens, epsilonO, zetaO); ///. zetaO is not used
		rockWaterSurface_.resistivity_=1e300;

		std::cout<<"WO.IFT:"<<watoil_.sigma()<<std::endl;
		std::cout<<"W.viscosity:"<<water_.viscosity()<<std::endl;
		std::cout<<"W.resistivity:"<<water_.resistivity()<<std::endl;
		std::cout<<"W.density:"<<water_.density()<<std::endl;
		std::cout<<"O.viscosity:"<<oil_.viscosity()<<std::endl;
		std::cout<<"O.resistivity:"<<oil_.resistivity()<<std::endl;
		std::cout<<"O.density:"<<oil_.density()<<std::endl; 

	}


	const Fluid& surface() const { return rockWaterSurface_; }
	const Fluid& surfaceOW() const { return rockWaterSurface__OW; }
	const Interface& oilWater() const { return watoil_; }




	double poreFillWeights(int i) const {return poreFillWeights_[i];}
	const std::string& poreFillAlg() const {return poreFillAlg_;}


	void setMaxPcLastDrainCycle(double pc)
	{	maxEverCappPress_ = std::max(pc, maxEverCappPress_);  	maxPcLastDrainCycle_ = pc;	}
	double maxEverPc() const {return maxEverCappPress_;}
	double maxPcLastDrainCycle() const {return maxPcLastDrainCycle_;}
	double minEverCappPress() const {return minEverCappPress_;}
	double minPcLastImbCycle() const {return minPcLastImbCycle_;}
	void setMinPcLastCycle(double pc) {	minEverCappPress_= std::min(pc, minEverCappPress_);	minPcLastImbCycle_= pc;  }

	void GuessCappPress(double pc) { guessPc_ = pc; }
	double GuessCappPress() const  { return guessPc_; }
	const dbl3& gravConst() const      { return gravConst_; }
	double rhogh(double rho, dbl3 dh) const {return rho * (gravConst_ & dh);}

	double KrwatcornAtSw0() const {return circWatCondMultFact_;}

	void injectant(const Fluid* injFluid) {injectant_ = injFluid;}
	const Fluid* injectant() const {return injectant_;}
	const Fluid& noFluid() const {return nof_;}




	size_t newOilTrappingIndex() const {  return trappedRegionsOil_.size();  }
	size_t newWatTrappingIndex() const {  return trappedRegionsWat_.size();  }

	size_t numTrappedOilElems() const  {
		size_t num(0); 	for(const auto& reg:trappedRegionsOil_)  { num += reg.size(); }  return num;	 }

	size_t numTrappedWatElems() const  {
		size_t num(0);	for(const auto& reg:trappedRegionsWat_) { num += reg.size(); }  return num;  }

	size_t numTrappedOilRegions() const  {
		size_t num(0);	for(const auto& reg:trappedRegionsOil_)  { num += !reg.empty(); }  return num;  }

	size_t numTrappedWatRegions() const  {
		size_t num(0);	for(const auto& reg:trappedRegionsWat_)  { num += !reg.empty(); }  	return num;  }


	///. deletes elem from trappedRegionsOil //called from checkUntrapOilIfUnstableConfigsImb
	void removeTrappedOilElem(int idx, Elem* elem)  {
		trappedRegionsOil__[idx].erase(
			std::remove(trappedRegionsOil__[idx].begin(),trappedRegionsOil__[idx].end(),elem),
			trappedRegionsOil__[idx].end());
	}

	///. deletes elem blob from trappedRegionsWat //called from checkUntrapWaterIfUnstableConfigsDrain
	void removeTrappeWatElem(int idx, std::pair<Elem*, FluidBlob> elem)  {
		trappedRegionsWat__[idx].erase(
			std::remove(trappedRegionsWat__[idx].begin(),trappedRegionsWat__[idx].end(),elem),
			trappedRegionsWat__[idx].end());
	}




	void addTrappeWatElem(int idx, std::pair<Elem*, FluidBlob> elem) {    trappedRegionsWat__[idx].push_back(elem); }

	void addTrappedRegionOil(const std::vector<Elem*>& elems) {trappedRegionsOil__.push_back(elems);}
	void addTrappedRegionWat(const std::vector< std::pair<Elem*, FluidBlob> >& elems) {trappedRegionsWat__.push_back(elems);}
	void removeTrappedRegionOil(int region) {trappedRegionsOil__[region].clear();}
	void removeTrappedRegionWat(int region) {trappedRegionsWat__[region].clear();}

	const std::vector<Elem*>& trappedRegionsOil(int region) const {return trappedRegionsOil_[region];}
	const std::vector< std::pair<Elem*,FluidBlob> >& trappedRegionsWat(int region) const {return trappedRegionsWat_[region];}

	/// for writing statistics
	const std::vector< std::vector<Elem*> >	&						trappedOilRegions() const {return trappedRegionsOil_;}
	const std::vector< std::vector< std::pair<Elem*,FluidBlob> > > &    trappedWatRegions() const {return trappedRegionsWat_;}



	int numPores() const {return nPors_;}
	int numThroats() const {return nTrots_;}
	void setNumElem(int numP, int numT) {nPors_ = numP; nTrots_ = numT;}

	void countCircle()   const {++numCircles_;}
	void countTriangle() const {++numTriangles_;}
	void countSquare()   const {++numSquares_;}
	void countPorous()   const {++numPorous_;}
	int numCircles() const {return numCircles_;}
	int numTriangles() const {return numTriangles_;}
	int numSquares() const {return numSquares_;}
	int numPorous() const {return numPorous_;}    

	std::ofstream                        dbgOut;

 private:
	Fluid                                       nof_; ///. TODO make const
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
	int   												nPors_;
	int   												nTrots_;
	bool                                                drainagePhase_;


	std::vector<double>   								poreFillWeights_;
	std::string												poreFillAlg_;
	std::vector< std::vector<Elem*> >							trappedRegionsOil__;
	std::vector< std::vector< std::pair<Elem*,FluidBlob> > >		trappedRegionsWat__;
	const std::vector< std::vector<Elem*> >	&						trappedRegionsOil_;
	const std::vector< std::vector< std::pair<Elem*,FluidBlob> > >	&    trappedRegionsWat_;


	///.  not used data
	mutable int										numSquares_;
	mutable int										numTriangles_;
	mutable int										numCircles_;
	mutable int										numPorous_;

 public:
	int										informative;
};





#endif
